"""Tests of ``mellor.llm.tools``, ``mellor.llm.prompts`` and ``mellor.llm.agent``.

The pipeline modules (data, library, runner, analysis, figures, report, sources)
are replaced by fakes through ``monkeypatch``; the model is a ``ScriptedBackend``.
``get_partial`` is tested on the toolkit's recorded water trajectory with the real
``epsrvis``.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from mellor import analysis as analysis_module
from mellor import data as data_module
from mellor import figures as figures_module
from mellor import library as library_module
from mellor import report as report_module
from mellor import runner as runner_module
from mellor import sources as sources_module
from mellor.llm import prompts
from mellor.llm.agent import DISCUSSION_NAME, TRANSCRIPT_NAME, Mellor, PipelineOutcome
from mellor.llm.backend import AssistantMessage, ScriptedBackend, ToolCall
from mellor.llm.tools import (
    Session,
    SessionTools,
    Tool,
    ToolError,
    ToolRegistry,
    compact_analysis,
    default_tools,
    jsonable,
    resolve_data_files,
    spec_summary,
    trajectory_folder,
)
from mellor.runner import RunError, RunResult
from mellor.spec import SystemSpec

REPO = Path(__file__).resolve().parents[2]
WATER_DATA = REPO / "experimental_data" / "water" / "data"
WATER_FILES = sorted(WATER_DATA.glob("*.mint01"))
EXAMPLE_TRAJECTORY = (
    REPO
    / "visualisation"
    / "python_visualisation_toolkit"
    / "example_H2O"
    / "data"
    / "trajectory"
)

pytestmark = pytest.mark.skipif(not WATER_FILES, reason="water .mint01 files not found")


# --------------------------------------------------------------------------- helpers
def _label(path: Path) -> str:
    stem = path.stem
    return stem.split("-", 1)[1] if "-" in stem else stem


def water_spec_dict(files: list[Path]) -> dict:
    """A valid water specification for the given data files."""
    datasets = []
    for path in files:
        label = _label(path)
        if "D2O" in label:
            isotopes = {"HW": [[2, 1.0]]}
        elif "HDO" in label:
            isotopes = {"HW": [[1, 0.5], [2, 0.5]]}
        else:
            isotopes = {"HW": [[0, 1.0]]}
        datasets.append(
            {
                "file": str(path),
                "label": label,
                "isotopes": isotopes,
                "exchangeable": ["HW"],
            }
        )
    return {
        "name": "waterbox",
        "title": "Liquid water at 25 C (SANDALS)",
        "kind": "molecular",
        "atom_types": [
            {
                "label": "OW",
                "element": "O",
                "epsilon": 0.6503,
                "sigma": 3.165492,
                "mass": 15.9994,
                "charge": -0.82,
            },
            {
                "label": "HW",
                "element": "H",
                "epsilon": 0.0,
                "sigma": 0.0,
                "mass": 2.0,
                "charge": 0.41,
            },
        ],
        "molecules": [
            {
                "name": "Water",
                "count": 1000,
                "atoms": [
                    {"label": "OW", "x": 0.0, "y": 0.0, "z": 0.0},
                    {"label": "HW", "x": 0.976, "y": 0.0, "z": 0.0},
                    {"label": "HW", "x": -0.2875, "y": 0.9327, "z": 0.0},
                ],
                "bonds": [[1, 2], [1, 3]],
            }
        ],
        "density": 0.1,
        "temperature": 300.0,
        "datasets": datasets,
        "schedule": {"ereq": 3.0},
        "notes": "test specification",
    }


def _manifest(spec: SystemSpec, quick: bool) -> dict:
    runs = []
    for index, (stage, niter, ereq) in enumerate(
        spec.schedule.iterations(quick), start=1
    ):
        runs.append(
            {
                "run": index,
                "stage": stage,
                "niter": niter,
                "ereq": ereq,
                "duration_s": 0.5 * index,
                "summary": {
                    "r_factor": 0.05 / index,
                    "energy": -10.0 - index,
                    "pressure": 1.0,
                    "ep_amplitude": ereq,
                },
            }
        )
    return {
        "name": spec.name,
        "title": spec.title,
        "kind": spec.kind,
        "quick": quick,
        "epsr_folder": "epsr",
        "trajectory": "trajectory",
        "runs": runs,
        "final": dict(runs[-1]["summary"]),
    }


def _fake_analysis(result: RunResult) -> dict:
    labels = [dataset.label for dataset in result.spec.datasets]
    return {
        "system": {
            "name": result.spec.name,
            "n_atoms": result.spec.n_atoms,
            "box_side": np.float64(31.07),
        },
        "schedule": {"quick": result.quick},
        "runs": [
            {
                "run": run["run"],
                "stage": run["stage"],
                "niter": run["niter"],
                "ereq": run["ereq"],
                "r_factor": run["summary"]["r_factor"],
                "r_factor_per_dataset": dict.fromkeys(
                    labels, run["summary"]["r_factor"]
                ),
                "energy": run["summary"]["energy"],
                "pressure": run["summary"]["pressure"],
                "ep_amplitude": run["summary"]["ep_amplitude"],
                "duration_s": run["duration_s"],
                "chisq": np.nan,
            }
            for run in result.runs
        ],
        "fit": {
            label: {
                "initial_r_factor": 0.05,
                "final_r_factor": 0.0123,
                "r2": np.float32(0.98),
            }
            for label in labels
        },
        "convergence": {"verdict": "converged", "ep_at_limit": False},
        "shells": {"OW-OW": {"r_peak": 2.8, "n_coord": np.float64(4.4)}},
        "flags": ["quick schedule"] if result.quick else [],
        "figures": {},
    }


class Fakes:
    """The fake pipeline modules and what they recorded."""

    def __init__(self, monkeypatch, work_dir: Path, *, run_fails: bool = False):
        self.work_dir = work_dir
        self.run_fails = run_fails
        self.run_calls: list[dict] = []
        self.pdf_calls: list[dict] = []
        self.figure_calls = 0
        self.downloads: list[str] = []
        self.searches: list[dict] = []
        m = monkeypatch

        # data
        def describe_all(paths):
            return [
                {
                    "file": Path(p).name,
                    "source": f"sls{18498 + i}.mdcs01",
                    "title": f"{_label(Path(p))} in can : 25C",
                    "instrument": "SANDALS (ISIS)",
                    "run_number": str(18498 + i),
                    "n_points": 1001,
                    "q_min": 0.1,
                    "q_max": 50.0,
                    "f_high_q": 0.0,
                    "sample": {
                        "hydrogen": "D" if "D2O" in Path(p).name else "H",
                        "candidates": ["water"],
                    },
                }
                for i, p in enumerate(paths)
            ]

        m.setattr(data_module, "describe_all", describe_all)
        m.setattr(
            data_module,
            "read_mint01",
            lambda path: SimpleNamespace(path=Path(path), title=Path(path).stem),
        )

        # library
        def make_spec(key, datasets, *, name=None, **adjust):
            if key != "water":
                raise KeyError(key)
            spec = SystemSpec.from_dict(water_spec_dict([d.path for d in datasets]))
            if name:
                spec.name = name
            for field_name, value in adjust.items():
                if field_name == "ereq":
                    spec.schedule.ereq = float(value)
                else:
                    setattr(spec, field_name, value)
            spec.source = "library:water"
            return spec

        def auto_spec(datasets, *, name=None):
            spec = make_spec("water", datasets, name=name)
            return spec, "the titles mention H2O, D2O and HDO"

        m.setattr(
            library_module,
            "describe_library",
            lambda: [
                {
                    "key": "water",
                    "title": "Liquid water",
                    "description": "SPC/E-like water",
                    "atom_types": ["OW", "HW"],
                    "keywords": ["water", "h2o", "d2o"],
                },
                {
                    "key": "silica",
                    "title": "Amorphous silica",
                    "description": "SiO2 glass",
                    "atom_types": ["Si", "O"],
                    "keywords": ["silica"],
                },
            ],
        )
        m.setattr(library_module, "make_spec", make_spec)
        m.setattr(library_module, "auto_spec", auto_spec)
        m.setattr(library_module, "KNOWN_SYSTEMS", {"water": None, "silica": None})

        # runner
        def run_spec(
            spec,
            folder,
            *,
            quick=False,
            record=True,
            overwrite=False,
            timeout_s=None,
            progress=None,
        ):
            self.run_calls.append(
                {
                    "spec": spec,
                    "folder": Path(folder),
                    "quick": quick,
                    "overwrite": overwrite,
                }
            )
            spec.validate(check_files=True)
            folder = Path(folder)
            folder.mkdir(parents=True, exist_ok=True)
            if progress is not None:
                progress("setup: makemole Water")
            if self.run_fails:
                raise RunError(
                    "the script failed (exit 1)",
                    log="\n".join(f"line {i}" for i in range(60)),
                )
            spec.save(folder / "spec.json")
            manifest = _manifest(spec, quick)
            (folder / "manifest.json").write_text(json.dumps(manifest))
            (folder / "run.log").write_text("ok\n")
            return RunResult(
                folder=folder, spec=spec, manifest=manifest, duration_s=1.5, quick=quick
            )

        def load_result(folder):
            folder = Path(folder)
            spec = SystemSpec.load(folder / "spec.json")
            manifest = json.loads((folder / "manifest.json").read_text())
            return RunResult(
                folder=folder,
                spec=spec,
                manifest=manifest,
                quick=bool(manifest.get("quick")),
            )

        m.setattr(runner_module, "run_spec", run_spec)
        m.setattr(runner_module, "load_result", load_result)

        # analysis
        m.setattr(analysis_module, "analyse", _fake_analysis)
        m.setattr(
            analysis_module,
            "save_analysis",
            lambda analysis, path: Path(path).write_text(json.dumps(jsonable(analysis)))
            and Path(path),
        )
        m.setattr(
            analysis_module,
            "load_analysis",
            lambda path: json.loads(Path(path).read_text()),
        )
        m.setattr(
            analysis_module,
            "analysis_markdown",
            lambda analysis: "## Fit\n\nfinal R-factor 0.0123 for every data set\n",
        )

        # figures
        def make_figures(result, out_dir=None, *, dpi=150):
            self.figure_calls += 1
            out = Path(out_dir) if out_dir else result.folder / "figures"
            out.mkdir(parents=True, exist_ok=True)
            paths = {}
            for name in ("01_convergence", "11_fit_final"):
                path = out / f"{name}.png"
                path.write_bytes(b"")
                paths[name] = path
            return paths

        m.setattr(figures_module, "make_figures", make_figures)

        # report
        def build_pdf(
            result,
            analysis,
            figures,
            discussion_markdown,
            out_path=None,
            *,
            model_name="",
            data_descriptions=None,
        ):
            out = Path(out_path) if out_path else result.folder / "report.pdf"
            out.write_bytes(b"%PDF-1.4 fake\n")
            self.pdf_calls.append(
                {
                    "discussion": discussion_markdown,
                    "model_name": model_name,
                    "figures": figures,
                    "analysis": analysis,
                    "data_descriptions": data_descriptions,
                    "out": out,
                }
            )
            return out

        m.setattr(report_module, "build_pdf", build_pdf)

        # sources
        def find_similar_data(
            query, *, sources=("local", "dissolve", "zenodo"), limit=10
        ):
            self.searches.append({"query": query, "sources": sources, "limit": limit})
            return {
                "query": query,
                "candidates": [
                    {
                        "source": "local",
                        "title": "water",
                        "url": "https://example.org/x.mint01",
                    }
                ],
                "notes": [],
            }

        def download(url, folder, *, timeout=60.0):
            self.downloads.append(url)
            path = Path(folder) / url.rsplit("/", 1)[-1]
            path.write_text("# fake\n")
            return path

        m.setattr(sources_module, "find_similar_data", find_similar_data)
        m.setattr(sources_module, "download", download)


@pytest.fixture
def work_dir(tmp_path: Path) -> Path:
    return tmp_path / "work" / "water"


@pytest.fixture
def fakes(monkeypatch, work_dir: Path) -> Fakes:
    return Fakes(monkeypatch, work_dir)


@pytest.fixture
def session(work_dir: Path) -> Session:
    return Session(work_dir=work_dir, data_files=list(WATER_FILES), quick=True)


# --------------------------------------------------------------------------- jsonable / registry
def test_jsonable_converts_numpy_paths_and_nan():
    value = {
        "a": np.float64(1.5),
        "b": np.array([1, 2]),
        "c": Path("/x"),
        "d": (1, 2),
        "e": float("nan"),
        "f": np.int32(3),
        "g": {1: np.bool_(True)},
    }
    out = jsonable(value)
    assert out == {
        "a": 1.5,
        "b": [1, 2],
        "c": "/x",
        "d": [1, 2],
        "e": None,
        "f": 3,
        "g": {"1": True},
    }
    json.dumps(out)


def _registry() -> ToolRegistry:
    registry = ToolRegistry()

    def add(a: int, b: int = 1):
        return {"sum": a + b}

    def boom():
        raise ValueError("kaboom")

    def picky():
        raise ToolError("not like this", "try that")

    def big():
        return {"text": "x" * 20_000}

    def numpy_result():
        return {"value": np.float32(2.5), "array": np.arange(3)}

    registry.register(
        Tool(
            "add",
            "add",
            {
                "type": "object",
                "properties": {"a": {"type": "integer"}, "b": {"type": "integer"}},
                "required": ["a"],
            },
            add,
        )
    )
    registry.register(Tool("boom", "fails", {"type": "object", "properties": {}}, boom))
    registry.register(
        Tool("picky", "tool error", {"type": "object", "properties": {}}, picky)
    )
    registry.register(Tool("big", "big", {"type": "object", "properties": {}}, big))
    registry.register(
        Tool(
            "numpy_result", "numpy", {"type": "object", "properties": {}}, numpy_result
        )
    )
    return registry


def test_registry_schemas_and_names():
    registry = _registry()
    assert registry.names() == ["add", "boom", "picky", "big", "numpy_result"]
    assert "add" in registry and len(registry) == 5
    schemas = registry.schemas(["big", "add"])
    assert [s["function"]["name"] for s in schemas] == ["add", "big"]
    assert schemas[0]["type"] == "function" and "parameters" in schemas[0]["function"]
    with pytest.raises(KeyError):
        registry.schemas(["nope"])


def test_registry_error_capture():
    registry = _registry()
    events: list = []
    registry.events = events
    assert registry.call("add", {"a": 2, "b": 3}) == {"sum": 5}
    unknown = registry.call("nope", {})
    assert "unknown tool" in unknown["error"] and "add" in unknown["hint"]
    bad = registry.call("add", {"c": 1})
    assert bad["error"].startswith("bad arguments for add") and "a, b" in bad["hint"]
    assert registry.call("add", '{"a": 1}') == {"sum": 2}
    assert "not valid JSON" in registry.call("add", "{bad")["error"]
    assert "must be a JSON object" in registry.call("add", [1])["error"]
    failed = registry.call("boom", {})
    assert (
        failed["error"].startswith("ValueError: kaboom") and "boom" in failed["error"]
    )
    assert failed["hint"]
    picky = registry.call("picky", {})
    assert picky == {"error": "not like this", "hint": "try that"}
    assert registry.call("numpy_result", {}) == {"value": 2.5, "array": [0, 1, 2]}
    assert [e["tool"] for e in events] == [
        "add",
        "nope",
        "add",
        "add",
        "add",
        "add",
        "boom",
        "picky",
        "numpy_result",
    ]
    assert (
        events[0]["ok"] is True
        and events[1]["ok"] is False
        and events[6]["error"].startswith("ValueError")
    )


def test_registry_call_json_truncation():
    registry = _registry()
    text = registry.call_json("big", {}, max_chars=500)
    assert len(text) < 700
    assert "[truncated:" in text and "more characters" in text
    short = registry.call_json("add", {"a": 1})
    assert json.loads(short) == {"sum": 2}
    assert ToolRegistry.serialise({"nan": float("nan")}) == '{"nan": null}'


# --------------------------------------------------------------------------- session
def test_session_load_existing(fakes: Fakes, work_dir: Path):
    spec = SystemSpec.from_dict(water_spec_dict(WATER_FILES))
    spec.source = "library:water"
    result = runner_module.run_spec(spec, work_dir, quick=True, overwrite=True)
    analysis_module.save_analysis(_fake_analysis(result), work_dir / "analysis.json")
    figures_module.make_figures(result)
    (work_dir / "figures" / "summary.md").write_text("# summary\n")
    (work_dir / "report.pdf").write_bytes(b"%PDF")
    (work_dir / "data").mkdir()
    (work_dir / "data" / "SLS18502-D2O.mint01").write_text("# copy\n")

    session = Session(work_dir=str(work_dir))
    session.load_existing()
    assert session.result is not None and session.result.spec.name == "waterbox"
    assert session.quick is True and session.spec_source == "library:water"
    assert session.analysis["fit"]["D2O"]["final_r_factor"] == 0.0123
    assert set(session.figures) == {"01_convergence", "11_fit_final", "summary.md"}
    assert session.report_path == work_dir / "report.pdf"
    assert session.data_files == [work_dir / "data" / "SLS18502-D2O.mint01"]
    assert session.events[-1]["type"] == "load_existing"


def test_session_load_existing_errors(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        Session(work_dir=tmp_path / "missing").load_existing()
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError, match="no finished run"):
        Session(work_dir=empty).load_existing()


# --------------------------------------------------------------------------- the tools
def test_default_tools_registers_every_documented_tool(session: Session):
    registry = default_tools(session)
    assert registry.names() == [
        "inspect_data",
        "list_known_systems",
        "library_spec",
        "validate_spec",
        "build_and_run",
        "get_results",
        "get_partial",
        "list_files",
        "read_text_file",
        "write_report",
        "find_similar_data",
        "download_data",
        "finish",
    ]
    assert registry.events is session.events
    schemas = registry.schemas()
    json.dumps(schemas)
    build = next(s for s in schemas if s["function"]["name"] == "build_and_run")
    assert build["function"]["parameters"]["required"] == ["spec"]
    assert (
        "atom_types"
        in build["function"]["parameters"]["properties"]["spec"]["properties"]
    )
    assert len(json.dumps(schemas)) < 20_000


def test_inspect_data_and_library(fakes: Fakes, session: Session):
    tools = SessionTools(session)
    described = tools.inspect_data()
    assert described["n_files"] == len(WATER_FILES)
    assert described["files"][0]["instrument"] == "SANDALS (ISIS)"
    assert session.data_descriptions == described["files"]
    with pytest.raises(ToolError, match="not found"):
        tools.inspect_data(files=["/nowhere/x.mint01"])
    assert Session(work_dir=session.work_dir).data_files == []
    with pytest.raises(ToolError, match="no data files"):
        SessionTools(Session(work_dir=session.work_dir)).inspect_data()

    systems = tools.list_known_systems()["systems"]
    assert [s["key"] for s in systems] == ["water", "silica"]

    result = tools.library_spec("Water", adjust={"ereq": 5})
    assert result["system"] == "water"
    assert result["spec"]["schedule"]["ereq"] == 5.0 and result["warnings"] == []
    assert (
        result["summary"]["n_atoms"] == 3000
        and result["summary"]["datasets"][0]["label"]
    )
    assert (
        session.draft_spec is not None and session.draft_spec.source == "library:water"
    )
    with pytest.raises(ToolError, match="unknown system"):
        tools.library_spec("plutonium")
    with pytest.raises(ToolError, match="adjust"):
        tools.library_spec("water", adjust=[1])


def test_validate_spec_and_file_resolution(fakes: Fakes, session: Session):
    tools = SessionTools(session)
    spec = water_spec_dict(WATER_FILES)
    for dataset in spec["datasets"]:
        dataset["file"] = Path(
            dataset["file"]
        ).name  # bare names, as a model writes them
    result = tools.validate_spec(spec)
    assert result["valid"] is True and result["errors"] == []
    assert session.draft_spec.datasets[0].file == str(WATER_FILES[0])
    assert result["summary"]["box_side"] == pytest.approx(31.07, abs=0.01)

    broken = water_spec_dict(WATER_FILES)
    del broken["schedule"]
    result = tools.validate_spec(broken)
    assert result["valid"] is False and "schedule.ereq" in result["errors"][0]

    bad_file = water_spec_dict(WATER_FILES)
    bad_file["datasets"][0]["file"] = "missing.mint01"
    result = tools.validate_spec(bad_file)
    assert result["valid"] is False and any("not found" in e for e in result["errors"])

    assert tools.validate_spec(json.dumps(spec))["valid"] is True
    assert tools.validate_spec("{bad")["valid"] is False
    assert tools.validate_spec({"spec": spec})["valid"] is True
    assert tools.validate_spec([1])["valid"] is False


def test_build_and_run_success(fakes: Fakes, session: Session):
    tools = SessionTools(session)
    lines: list[str] = []
    session.progress = lines.append
    result = tools.build_and_run(water_spec_dict(WATER_FILES))
    assert result["ok"] is True
    assert result["spec_source"] == "model" and result["quick"] is True
    assert [run["stage"] for run in result["runs"]] == [
        "scratch",
        "reference potential",
        "EP refinement",
        "EP refinement",
    ]
    assert result["fit"]["D2O"]["final_r_factor"] == 0.0123
    assert result["convergence"]["verdict"] == "converged"
    assert result["shells"]["OW-OW"]["n_coord"] == 4.4
    assert result["flags"] == ["quick schedule"]
    assert result["figures"] == ["01_convergence", "11_fit_final"]
    json.dumps(result)
    assert lines == ["setup: makemole Water"]
    assert (
        fakes.run_calls[0]["folder"] == session.work_dir
        and fakes.run_calls[0]["overwrite"] is True
    )
    assert session.result is not None and session.spec_source == "model"
    assert session.analysis["figures"] == {
        name: str(path) for name, path in session.figures.items()
    }
    assert json.loads((session.work_dir / "analysis.json").read_text())["fit"]["D2O"][
        "r2"
    ] == pytest.approx(0.98, abs=1e-3)
    assert [e["type"] for e in session.events] == ["run"]

    # a second successful run goes into a sub-folder, keeping the first
    again = tools.build_and_run(water_spec_dict(WATER_FILES), quick=False)
    assert again["ok"] is True and again["quick"] is False
    assert fakes.run_calls[1]["folder"] == session.work_dir / "rerun_1"
    assert session.result.folder == session.work_dir / "rerun_1"
    assert (session.work_dir / "manifest.json").exists()


def test_build_and_run_invalid_spec(fakes: Fakes, session: Session):
    tools = SessionTools(session)
    broken = water_spec_dict(WATER_FILES)
    broken["density"] = -1
    result = tools.build_and_run(broken)
    assert result["ok"] is False and any("density" in e for e in result["errors"])
    assert fakes.run_calls == [] and session.result is None


def test_build_and_run_run_error(monkeypatch, work_dir: Path, session: Session):
    fakes = Fakes(monkeypatch, work_dir, run_fails=True)
    tools = SessionTools(session)
    result = tools.build_and_run(water_spec_dict(WATER_FILES))
    assert result["ok"] is False and "failed" in result["error"]
    assert len(result["log_tail"]) == 40 and result["log_tail"][-1] == "line 59"
    assert session.result is None
    assert session.events[-1]["type"] == "run_failed"
    assert fakes.run_calls


def test_build_and_run_analysis_failure_is_reported(
    fakes: Fakes, monkeypatch, session: Session
):
    def analyse(result):
        raise ZeroDivisionError("no trajectory")

    monkeypatch.setattr(analysis_module, "analyse", analyse)
    result = SessionTools(session).build_and_run(water_spec_dict(WATER_FILES))
    assert result["ok"] is True
    assert result["analysis_error"].startswith("ZeroDivisionError")
    assert result["runs"] and result["figures"] == []
    assert session.analysis is None and session.result is not None


def test_get_results_sections(fakes: Fakes, session: Session):
    tools = SessionTools(session)
    with pytest.raises(ToolError, match="no results"):
        tools.get_results()
    tools.build_and_run(water_spec_dict(WATER_FILES))
    assert set(tools.get_results()) >= {"fit", "runs", "shells", "figures"}
    assert tools.get_results("fit") == {"fit": session.analysis["fit"]}
    assert tools.get_results("spec")["spec"]["name"] == "waterbox"
    assert tools.get_results("manifest")["manifest"]["runs"]
    assert set(tools.get_results("figures")["figures"]) == {
        "01_convergence",
        "11_fit_final",
    }
    with pytest.raises(ToolError, match="unknown section"):
        tools.get_results("mystery")


@pytest.mark.skipif(
    not (EXAMPLE_TRAJECTORY / "meta.json").is_file(),
    reason="example trajectory missing",
)
def test_get_partial_on_example_trajectory(fakes: Fakes, session: Session):
    tools = SessionTools(session)
    with pytest.raises(ToolError, match="no run"):
        tools.get_partial("OW-HW")
    tools.build_and_run(water_spec_dict(WATER_FILES))
    assert trajectory_folder(session.result) is None
    with pytest.raises(ToolError, match="no recorded trajectory"):
        tools.get_partial("OW-HW")
    os.symlink(EXAMPLE_TRAJECTORY, session.result.folder / "trajectory")
    assert trajectory_folder(session.result) == session.result.folder / "trajectory"

    result = tools.get_partial("hw-ow")
    assert result["pair"] == "OW-HW" and result["kind"] == "gr_inter"
    assert 40 <= len(result["x"]) <= 80 and len(result["x"]) == len(result["y"])
    assert max(result["x"]) <= 8.0 and result["x_label"] == "r (Å)"
    shell = result["first_shell"]
    assert 1.6 < shell["r_peak"] < 2.1  # the O...H hydrogen bond
    assert 1.0 < shell["coordination_number"] < 3.0 and shell["g_peak"] > 0.5
    json.dumps(result)

    sq = tools.get_partial("OW-OW", kind="sq_inter")
    assert (
        sq["x_label"].startswith("Q") and "first_shell" not in sq and len(sq["x"]) >= 40
    )
    total = tools.get_partial("OW-OW", kind="gr_total", r_max=4.0)
    assert max(total["x"]) <= 4.0 and "first_shell" in total
    assert tools.get_partial("HW-HW", kind="ep")["y_label"] == "U_EP (kJ/mol)"
    with pytest.raises(ToolError, match="unknown pair"):
        tools.get_partial("Si-O")
    with pytest.raises(ToolError, match="unknown kind"):
        tools.get_partial("OW-OW", kind="banana")


def test_list_and_read_files(fakes: Fakes, session: Session):
    tools = SessionTools(session)
    tools.build_and_run(water_spec_dict(WATER_FILES))
    listing = tools.list_files()
    assert listing["folder"] == str(session.work_dir)
    assert any(
        entry.startswith("figures/01_convergence.png") for entry in listing["files"]
    )
    assert "spec.json" in " ".join(listing["files"])
    figures = tools.list_files("figures")
    assert any(entry.startswith("11_fit_final.png") for entry in figures["files"])

    text = tools.read_text_file("spec.json", max_chars=300)
    assert text["truncated"] is True and text["text"].endswith("[truncated]")
    assert (
        json.loads(tools.read_text_file("manifest.json")["text"])["name"] == "waterbox"
    )
    data = tools.read_text_file(str(WATER_FILES[0]), max_chars=500)
    assert data["text"].startswith("#") and data["chars"] > 500
    with pytest.raises(ToolError, match="binary"):
        tools.read_text_file("figures/01_convergence.png")
    with pytest.raises(ToolError, match="not found"):
        tools.read_text_file("nothing.txt")
    with pytest.raises(ToolError, match="outside"):
        tools.read_text_file("/etc/hostname")
    with pytest.raises(ToolError, match="is a folder"):
        tools.read_text_file("figures")


def test_write_report_find_and_download(fakes: Fakes, session: Session):
    tools = SessionTools(session)
    with pytest.raises(ToolError, match="no run"):
        tools.write_report("text")
    tools.inspect_data()
    tools.build_and_run(water_spec_dict(WATER_FILES))
    session.model_name = "test-model"
    with pytest.raises(ToolError, match="empty"):
        tools.write_report("  ")
    written = tools.write_report("## Discussion\n\nGood fit.")
    assert written["ok"] is True and Path(written["report"]).is_file()
    assert session.report_path == session.work_dir / "report.pdf"
    call = fakes.pdf_calls[-1]
    assert (
        call["discussion"].startswith("## Discussion")
        and call["model_name"] == "test-model"
    )
    assert call["data_descriptions"] == session.data_descriptions
    assert session.events[-1]["type"] == "report"

    found = tools.find_similar_data("D2O water", sources="local, dissolve", limit=3)
    assert found["candidates"][0]["source"] == "local"
    assert fakes.searches[-1] == {
        "query": "D2O water",
        "sources": ("local", "dissolve"),
        "limit": 3,
    }
    with pytest.raises(ToolError, match="empty"):
        tools.find_similar_data("  ")

    downloaded = tools.download_data("https://example.org/data/SLS99999-H2O.mint01")
    assert downloaded["added_to_session"] is True
    assert session.data_files[-1] == session.data_folder / "SLS99999-H2O.mint01"
    assert fakes.downloads == ["https://example.org/data/SLS99999-H2O.mint01"]
    assert tools.finish("bye") == {"finished": True, "message": "bye"}


def test_spec_summary_and_compact_analysis():
    spec = SystemSpec.from_dict(water_spec_dict(WATER_FILES))
    summary = spec_summary(spec)
    assert summary["molecules"] == [{"name": "Water", "count": 1000, "atoms": 3}]
    assert (
        summary["schedule"]["total_iterations"] == 300
        and summary["schedule"]["quick_total_iterations"] == 31
    )
    assert summary["mass_density_g_cm3"] == pytest.approx(1.107, abs=0.01)  # HW mass 2
    json.dumps(summary)
    compact = compact_analysis(
        {
            "runs": [
                {"run": 1, "stage": "scratch", "extra": 1, "r_factor": np.float64(0.1)}
            ],
            "fit": {},
            "flags": ("a",),
        }
    )
    assert compact["runs"] == [
        {"run": 1, "stage": "scratch", "r_factor": 0.1}
    ] and compact["flags"] == ["a"]
    assert compact["convergence"] is None


def test_resolve_data_files_keeps_existing_paths(session: Session):
    spec = SystemSpec.from_dict(water_spec_dict(WATER_FILES))
    spec.datasets[0].file = "unknown-name.mint01"
    resolve_data_files(spec, session)
    assert spec.datasets[0].file == "unknown-name.mint01"
    assert spec.datasets[1].file == str(WATER_FILES[1])


# --------------------------------------------------------------------------- prompts
def test_system_prompt_is_short_and_complete():
    words = prompts.SYSTEM_PROMPT.split()
    assert 300 < len(words) < 900
    for needle in (
        "Mellor",
        "EPSR",
        "Empirical Potential Structure Refinement",
        "ereq",
        "build_and_run",
        "library",
        "isotopes",
        "one tool call",
    ):
        assert needle in prompts.SYSTEM_PROMPT, needle


def test_format_setup_task(fakes: Fakes):
    descriptions = data_module.describe_all(WATER_FILES)
    for description, path in zip(descriptions, WATER_FILES, strict=True):
        description["path"] = str(path)
    text = prompts.format_setup_task(
        descriptions, library_module.describe_library(), "use ereq 5", "/tmp/work", True
    )
    assert f"Data files ({len(WATER_FILES)})" in text
    assert str(WATER_FILES[0]) in text and "SANDALS (ISIS)" in text
    assert "- water: Liquid water" in text and "keywords: water, h2o, d2o" in text
    assert "Work folder: /tmp/work" in text and "Quick preview run" in text
    assert "User request: use ereq 5" in text
    full = prompts.format_setup_task([], [], None, "w", False)
    assert (
        "Full schedule" in full
        and "(none)" in full
        and "(the library is empty)" in full
    )
    assert "{" not in full.replace('{"', "")  # no unfilled placeholders


def test_format_narrative_and_chat():
    text = prompts.format_narrative_task("## Fit\n\nR 0.01", "why {this} model")
    assert "## Fit" in text and "why {this} model" in text and "300-600 words" in text
    assert "{analysis_markdown}" not in text
    assert "(none)" in prompts.format_narrative_task("", "")
    chat = prompts.format_chat_context(
        "summary text", ["spec.json", "figures/"], "/w/report.pdf"
    )
    assert (
        "summary text" in chat
        and "spec.json, figures/" in chat
        and "/w/report.pdf" in chat
    )
    chat = prompts.format_chat_context(None, [f"f{i}" for i in range(50)], None)
    assert "10 more" in chat and "not written yet" in chat and "get_results" in chat


# --------------------------------------------------------------------------- the agent
NARRATIVE = (
    "## What was modelled\n\nLiquid water, 1000 molecules.\n\n## Quality of the fit\n\n"
    "The final R-factor is 0.0123 for every data set.\n\n## What the structure shows\n\n"
    "A hydrogen-bonded network.\n\n## Limitations\n\nQuick schedule.\n\n## Next steps\n\nRun the full schedule.\n"
)


def scripted_model(messages: list[dict]) -> AssistantMessage:
    """A model that follows the setup procedure, then writes the narrative, then chats."""
    last = messages[-1]
    if last["role"] == "tool":
        content = json.loads(last["content"])
        if last["name"] == "inspect_data":
            return AssistantMessage(
                "Water data; using the library.",
                [ToolCall("library_spec", {"system": "water"})],
            )
        if last["name"] == "library_spec":
            spec = content["spec"]
            spec["title"] = "Liquid water at 25 C, SANDALS (model-adapted)"
            return AssistantMessage("", [ToolCall("build_and_run", {"spec": spec})])
        if last["name"] == "build_and_run":
            return AssistantMessage("Done." if content.get("ok") else "Failed.")
        if last["name"] == "get_results":
            fit = content["fit"]
            return AssistantMessage(
                f"The final R-factor of D2O is {fit['D2O']['final_r_factor']}."
            )
    if last["role"] == "user":
        if "Set up and run" in last["content"]:
            return AssistantMessage(
                "Let me look at the data.", [ToolCall("inspect_data", {})]
            )
        if "Discussion" in last["content"]:
            return AssistantMessage(NARRATIVE)
        if "R-factor" in last["content"]:
            return AssistantMessage("", [ToolCall("get_results", {"section": "fit"})])
    return AssistantMessage("I do not know.")


def _read_transcript(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def test_run_pipeline_with_scripted_model(fakes: Fakes, session: Session):
    backend = ScriptedBackend(scripted_model)
    events: list[dict] = []
    agent = Mellor(backend, session, on_event=events.append)
    assert (
        agent.messages[0]["role"] == "system"
        and "Mellor" in agent.messages[0]["content"]
    )
    assert session.model_name == "scripted"

    outcome = agent.run_pipeline("Refine the water data.")
    assert isinstance(outcome, PipelineOutcome)
    assert outcome.spec_source == "library:water"  # a library spec adapted by the model
    assert session.spec_source == "library:water"
    assert outcome.steps == 4 and agent.spec_failures == 0
    assert (
        outcome.report_path == session.work_dir / "report.pdf"
        and outcome.report_path.is_file()
    )
    assert (
        "Report:" in outcome.final_message and "quick schedule" in outcome.final_message
    )

    # the session
    assert session.result is not None and session.result.spec.title.endswith(
        "(model-adapted)"
    )
    assert session.analysis["fit"]["D2O"]["final_r_factor"] == 0.0123
    assert set(session.figures) == {"01_convergence", "11_fit_final"}
    assert session.data_descriptions and session.data_descriptions[0]["path"] == str(
        WATER_FILES[0]
    )
    assert [e["type"] for e in session.events] == [
        "tool_call",
        "tool_call",
        "run",
        "tool_call",
        "narrative",
        "report",
    ]
    assert session.events[-2]["source"] == "model"

    # the narrative reached the report, and discussion.md was written
    assert fakes.pdf_calls[-1]["discussion"] == NARRATIVE.strip()
    assert fakes.pdf_calls[-1]["model_name"] == "scripted"
    assert fakes.pdf_calls[-1]["data_descriptions"] == session.data_descriptions
    assert fakes.pdf_calls[-1]["figures"] == session.figures
    assert (session.work_dir / DISCUSSION_NAME).read_text() == NARRATIVE
    assert agent.narrative_source == "model"

    # the conversation and the transcript
    roles = [m["role"] for m in agent.messages]
    assert roles == [
        "system",
        "user",
        "assistant",
        "tool",
        "assistant",
        "tool",
        "assistant",
        "tool",
        "user",
        "assistant",
    ]
    assert (
        agent.messages[3]["name"] == "inspect_data"
        and agent.messages[3]["tool_call_id"] == "call_1"
    )
    assert (
        "Discussion" in agent.messages[8]["content"]
        and "final R-factor 0.0123" in agent.messages[8]["content"]
    )
    narrative_request = backend.requests[-1]
    assert (
        narrative_request[-1]["role"] == "user"
        and "Specification notes: test specification"
        in narrative_request[-1]["content"]
    )
    transcript = _read_transcript(session.work_dir / TRANSCRIPT_NAME)
    assert transcript[0]["role"] == "pipeline" and transcript[-1]["role"] == "pipeline"
    bookkeeping = [t for t in transcript if t["role"] == "narrative"]
    assert len(bookkeeping) == 1 and bookkeeping[0]["source"] == "model"
    assert bookkeeping[0]["model"] == "scripted"
    assert [t["role"] for t in transcript[1:-1] if t["role"] != "narrative"] == roles[
        1:
    ]
    assert all("time" in entry for entry in transcript)
    assert transcript[2]["raw"] == "Let me look at the data."

    # the events for the display
    kinds = [e["type"] for e in events]
    assert kinds[:3] == ["phase", "assistant", "tool_call"]
    assert (
        kinds.count("tool_result") == 3
        and "progress" in kinds
        and kinds[-1] == "report"
    )
    assert any(
        e["type"] == "tool_result" and e["name"] == "build_and_run" and e["ok"]
        for e in events
    )

    # chat afterwards: a question answered through get_results
    answer = agent.chat("What is the final R-factor?")
    assert answer == "The final R-factor of D2O is 0.0123."
    assert (
        agent.messages[-2]["role"] == "tool"
        and agent.messages[-2]["name"] == "get_results"
    )
    assert (
        "follow-up questions" in agent.messages[-4]["content"]
        and "User: What is the final R-factor?" in agent.messages[-4]["content"]
    )
    assert outcome.steps + 2 == agent.steps
    second = agent.chat("Anything else?")
    assert (
        second == "I do not know."
        and "follow-up questions" not in agent.messages[-2]["content"]
    )


def test_pipeline_fallback_after_invalid_specs(fakes: Fakes, session: Session):
    replies = []
    # distinct broken specs: identical calls would trip the repeated-call guard first
    for attempt in range(3):
        broken = water_spec_dict(WATER_FILES)
        broken["atom_types"][0]["mass"] = 0.0
        broken["title"] = f"attempt {attempt}"
        replies.append(
            AssistantMessage("Trying.", [ToolCall("build_and_run", {"spec": broken})])
        )
    replies.append(AssistantMessage(""))  # no narrative
    backend = ScriptedBackend(replies)
    agent = Mellor(backend, session, max_spec_attempts=3)
    outcome = agent.run_pipeline()

    assert agent.spec_failures == 3
    assert outcome.steps == 5  # 3 attempts, then an empty narrative reply and a nudge
    assert outcome.spec_source == "fallback:library:water"
    assert session.result.spec.source == "fallback:library:water"
    assert (
        "fell back to the built-in library system 'water'" in session.result.spec.notes
    )
    assert "3 invalid specification(s)" in session.result.spec.notes
    assert fakes.run_calls[-1]["folder"] == session.work_dir
    assert outcome.report_path.is_file()
    assert any(
        e["type"] == "fallback" and e["system"] == "water" for e in session.events
    )
    assert any("stopped: 3 invalid" in note for note in agent.notes)

    # narrative fallback: the discussion is the factual summary, clearly labelled
    discussion = fakes.pdf_calls[-1]["discussion"]
    assert discussion.startswith("_The language model did not produce a discussion")
    assert (
        "final R-factor 0.0123" in discussion and agent.narrative_source == "fallback"
    )
    assert "automatic summary" in outcome.final_message
    transcript = _read_transcript(session.work_dir / TRANSCRIPT_NAME)
    assert any(t["role"] == "fallback" for t in transcript)
    assert any(t["role"] == "note" and "fell back" in t["content"] for t in transcript)
    fallback_note = next(
        m
        for m in agent.messages
        if m["role"] == "user" and m["content"].startswith("[Mellor]")
    )
    assert "library fallback ('water')" in fallback_note["content"]


def test_pipeline_fallback_when_model_never_runs(fakes: Fakes, session: Session):
    backend = ScriptedBackend(["I cannot set this up.", NARRATIVE])
    agent = Mellor(backend, session)
    outcome = agent.run_pipeline()
    assert outcome.spec_source == "fallback:library:water"
    assert "stopped without a successful build_and_run" in session.result.spec.notes
    assert (
        fakes.pdf_calls[-1]["discussion"] == NARRATIVE.strip()
        and agent.narrative_source == "model"
    )
    assert outcome.steps == 2


def test_pipeline_fallback_run_failure_raises(
    monkeypatch, work_dir: Path, session: Session
):
    Fakes(monkeypatch, work_dir, run_fails=True)
    agent = Mellor(ScriptedBackend(["No."]), session)
    with pytest.raises(RuntimeError, match="fallback run failed"):
        agent.run_pipeline()


def test_pipeline_fallback_without_library_match_raises(
    fakes: Fakes, monkeypatch, session: Session
):
    def auto_spec(datasets, *, name=None):
        raise LookupError("no template matches")

    monkeypatch.setattr(library_module, "auto_spec", auto_spec)
    agent = Mellor(ScriptedBackend(["No."]), session)
    with pytest.raises(RuntimeError, match="no library system matches"):
        agent.run_pipeline()


def test_pipeline_rejects_missing_data(work_dir: Path):
    with pytest.raises(ValueError, match="no data files"):
        Mellor(ScriptedBackend([]), Session(work_dir=work_dir)).run_pipeline()
    with pytest.raises(FileNotFoundError):
        Mellor(
            ScriptedBackend([]),
            Session(work_dir=work_dir, data_files=[work_dir / "x.mint01"]),
        ).run_pipeline()


def test_repeated_identical_call_guard(fakes: Fakes, session: Session):
    backend = ScriptedBackend(
        [AssistantMessage("", [ToolCall("list_known_systems", {})]) for _ in range(6)]
    )
    agent = Mellor(backend, session)
    answer = agent.chat("Which systems do you know?")
    assert agent.steps == 3
    assert answer.startswith("Mellor gave no answer")
    assert any("repeated three times" in note for note in agent.notes)
    tool_messages = [m for m in agent.messages if m["role"] == "tool"]
    assert (
        len(tool_messages) == 3 and "already made twice" in tool_messages[-1]["content"]
    )
    assert len([e for e in session.events if e["type"] == "tool_call"]) == 2


def test_varying_calls_are_not_repeats(fakes: Fakes, session: Session):
    replies = [
        AssistantMessage("", [ToolCall("get_results", {"section": s})])
        for s in ("fit", "runs", "fit", "runs")
    ]
    replies.append(AssistantMessage("Summary given."))
    agent = Mellor(ScriptedBackend(replies), session)
    SessionTools(session).build_and_run(water_spec_dict(WATER_FILES))
    assert agent.chat("Tell me everything.") == "Summary given."
    assert agent.steps == 5


def test_empty_reply_guard_and_nudge(fakes: Fakes, session: Session):
    backend = ScriptedBackend(["", "<tool_call>{broken", ""])
    agent = Mellor(backend, session)
    answer = agent.chat("Hello?")
    assert agent.steps == 2 and answer.startswith("Mellor gave no answer")
    assert any("two empty or malformed replies" in note for note in agent.notes)
    assert (
        agent.messages[-2]["role"] == "user"
        and "empty" in agent.messages[-2]["content"]
    )

    backend = ScriptedBackend(["<tool_call>{broken", AssistantMessage("Recovered.")])
    agent = Mellor(backend, session)
    assert agent.chat("Again?") == "Recovered."
    nudges = [
        m
        for m in agent.messages
        if m["role"] == "user" and "could not be parsed" in m["content"]
    ]
    assert len(nudges) == 1


def test_max_steps_guard(fakes: Fakes, session: Session):
    def always_call(messages):
        return AssistantMessage(
            "", [ToolCall("get_results", {"section": str(len(messages))})]
        )

    agent = Mellor(ScriptedBackend(always_call), session, max_steps=4)
    agent.chat("loop")
    assert agent.steps == 4 and any("limit of 4" in note for note in agent.notes)


def test_finish_tool_and_phase_restriction(fakes: Fakes, session: Session):
    replies = [
        AssistantMessage("", [ToolCall("find_similar_data", {"query": "water"})]),
        AssistantMessage("", [ToolCall("finish", {"message": "All done."})]),
    ]
    agent = Mellor(ScriptedBackend(replies), session)
    assert agent._loop(allowed_tools=["get_results", "finish"]) == "All done."
    refused = [m for m in agent.messages if m["role"] == "tool"][0]
    assert "not available in this phase" in refused["content"]
    assert fakes.searches == []


def test_chat_on_loaded_run_gets_context(
    fakes: Fakes, session: Session, work_dir: Path
):
    SessionTools(session).build_and_run(water_spec_dict(WATER_FILES))
    loaded = Session(work_dir=work_dir)
    loaded.load_existing()
    replies = [AssistantMessage("The run converged.")]
    agent = Mellor(ScriptedBackend(replies), loaded)
    assert agent.chat("Did it converge?") == "The run converged."
    context = agent.messages[1]["content"]
    assert "Summary of the run" in context and "final R-factor 0.0123" in context
    assert "figures/01_convergence.png" in context and "spec.json" in context


def test_transcript_survives_folder_recreation(
    fakes: Fakes, session: Session, monkeypatch
):
    original = runner_module.run_spec

    def wiping_run_spec(spec, folder, **kwargs):
        for path in Path(folder).iterdir():
            if path.is_file():
                path.unlink()
        return original(spec, folder, **kwargs)

    monkeypatch.setattr(runner_module, "run_spec", wiping_run_spec)
    agent = Mellor(ScriptedBackend(scripted_model), session)
    agent.run_pipeline()
    transcript = _read_transcript(session.work_dir / TRANSCRIPT_NAME)
    assert transcript[0]["role"] == "pipeline" and transcript[0]["content"] == "start"
    assert len(transcript) == len(agent._transcript)


def test_on_event_failure_does_not_stop_the_agent(fakes: Fakes, session: Session):
    def broken(event):
        raise RuntimeError("display broke")

    agent = Mellor(ScriptedBackend(["Fine."]), session, on_event=broken)
    assert agent.chat("hi") == "Fine."
    transcript = _read_transcript(session.work_dir / TRANSCRIPT_NAME)
    assert any("on_event callback failed" in t.get("content", "") for t in transcript)
