"""``mellor.cli``: the fast paths (spec, search, info, usage and error handling, a
run and a report through fakes of the EPSR-bound modules, the chat loop) and,
marked slow, the real quick water run through ``mellor run --system water``
followed by ``mellor report`` and ``mellor chat`` on the same folder.
"""

from __future__ import annotations

import builtins
import contextlib
import io
import json
import sys
from pathlib import Path
from typing import Any

import pytest
from example_run import WATER_FILES
from pypdf import PdfReader

from mellor import analysis as analysis_module
from mellor import cli
from mellor import figures as figures_module
from mellor import report as report_module
from mellor import runner as runner_module
from mellor import sources as sources_module
from mellor.llm import backend as backend_module
from mellor.llm.agent import DISCUSSION_NAME, NO_MODEL_NOTE, TRANSCRIPT_NAME
from mellor.llm.backend import AssistantMessage, ScriptedBackend, ToolCall
from mellor.runner import RunResult
from mellor.spec import SystemSpec

FILES = [str(path) for path in WATER_FILES]
NARRATIVE = (
    "## Discussion\n\nThe refinement of liquid water finished; the empirical potential "
    "reached its requested amplitude and the hydrogen-bonded network is reproduced.\n"
)

pytestmark = pytest.mark.skipif(
    not all(path.is_file() for path in WATER_FILES),
    reason="the water .mint01 files are not in the repository",
)


# --------------------------------------------------------------------------- fakes
class Fakes:
    """Stand-ins for the EPSR-bound modules; records the report calls."""

    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        self.pdf_calls: list[dict[str, Any]] = []
        self.run_calls: list[dict[str, Any]] = []
        monkeypatch.setattr(runner_module, "run_spec", self.run_spec)
        monkeypatch.setattr(analysis_module, "analyse", self.analyse)
        monkeypatch.setattr(figures_module, "make_figures", self.make_figures)
        monkeypatch.setattr(report_module, "build_pdf", self.build_pdf)

    def run_spec(
        self,
        spec: SystemSpec,
        folder: Path,
        *,
        quick: bool = False,
        record: bool = True,
        overwrite: bool = False,
        timeout_s: float | None = None,
        progress: Any = None,
    ) -> RunResult:
        spec.validate(check_files=True)
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        self.run_calls.append({"folder": folder, "quick": quick, "source": spec.source})
        if progress is not None:
            progress("$ fake run")
            progress("run 1/4: scratch")
        spec.save(folder / "spec.json")
        runs = []
        for index, (stage, niter, ereq) in enumerate(
            spec.schedule.iterations(quick), start=1
        ):
            summary = {
                "r_factor": 0.02 / index,
                "r_factor_per_dataset": [0.02 / index] * len(spec.datasets),
                "energy": -10.0 * index,
                "pressure": 2.0,
                "ep_amplitude": ereq,
            }
            runs.append(
                {
                    "run": index,
                    "stage": stage,
                    "niter": niter,
                    "ereq": ereq,
                    "duration_s": 0.1,
                    "summary": summary,
                }
            )
        manifest = {
            "name": spec.name,
            "title": spec.title,
            "kind": spec.kind,
            "quick": quick,
            "complete": True,
            "epsr_folder": "epsr",
            "trajectory": None,
            "runs": runs,
            "final": dict(runs[-1]["summary"]),
        }
        (folder / "manifest.json").write_text(json.dumps(manifest))
        (folder / "run.log").write_text("fake run\n")
        return RunResult(
            folder=folder, spec=spec, manifest=manifest, duration_s=0.2, quick=quick
        )

    @staticmethod
    def analyse(result: RunResult) -> dict[str, Any]:
        runs = [
            {
                "run": run["run"],
                "stage": run["stage"],
                "niter": run["niter"],
                "ereq": run["ereq"],
                "r_factor": run["summary"]["r_factor"],
                "energy": run["summary"]["energy"],
                "pressure": run["summary"]["pressure"],
                "ep_amplitude": run["summary"]["ep_amplitude"],
                "duration_s": run["duration_s"],
            }
            for run in result.runs
        ]
        return {
            "system": {
                "name": result.spec.name,
                "title": result.spec.title,
                "n_atoms": 3000,
            },
            "schedule": {"quick": result.quick, "runs": [], "total_iterations": 31},
            "runs": runs,
            "fit": {},
            "convergence": {
                "verdict": "a fake verdict.",
                "r_factor": runs[-1]["r_factor"],
            },
            "shells": {},
            "flags": ["fake run"],
            "figures": {},
        }

    @staticmethod
    def make_figures(
        result: RunResult, out_dir: Any = None, *, dpi: int = 150
    ) -> dict[str, Path]:
        out = Path(out_dir) if out_dir else result.folder / "figures"
        out.mkdir(parents=True, exist_ok=True)
        path = out / "01_convergence.png"
        path.write_bytes(b"")
        return {"01_convergence": path}

    def build_pdf(
        self,
        result: RunResult,
        analysis: dict[str, Any],
        figures: dict[str, Path],
        discussion_markdown: str,
        out_path: Any = None,
        *,
        model_name: str = "",
        data_descriptions: Any = None,
    ) -> Path:
        out = Path(out_path) if out_path else result.folder / "report.pdf"
        out.write_bytes(b"%PDF-1.4 fake\n")
        self.pdf_calls.append(
            {
                "discussion": discussion_markdown,
                "model_name": model_name,
                "figures": dict(figures),
                "analysis": analysis,
                "out": out,
            }
        )
        return out


class FakeModel:
    """A backend that is not a ``ScriptedBackend``, so the CLI treats it as a model."""

    name = "fake-model"

    def __init__(self, script: Any) -> None:
        self._inner = ScriptedBackend(script)
        self.requests = self._inner.requests

    def chat(
        self, messages: list[dict[str, Any]], tools: Any = None
    ) -> AssistantMessage:
        return self._inner.chat(messages, tools)


def model_script(messages: list[dict[str, Any]]) -> AssistantMessage:
    """Setup through the library, the narrative, then short chat answers."""
    last = messages[-1]
    if last["role"] == "tool":
        content = json.loads(last["content"])
        if last["name"] == "inspect_data":
            return AssistantMessage(
                "Water.", [ToolCall("library_spec", {"system": "water"})]
            )
        if last["name"] == "library_spec":
            return AssistantMessage(
                "", [ToolCall("build_and_run", {"spec": content["spec"]})]
            )
        if last["name"] == "build_and_run":
            return AssistantMessage("Done.")
        if last["name"] == "get_results":
            return AssistantMessage(f"The verdict: {content['convergence']['verdict']}")
    text = last.get("content") or ""
    if text.startswith("Set up and run an EPSR refinement"):
        return AssistantMessage("Looking.", [ToolCall("inspect_data", {})])
    if 'Write the "Discussion" section' in text:
        return AssistantMessage(NARRATIVE)
    if "User: verdict?" in text or text == "verdict?":
        return AssistantMessage(
            "", [ToolCall("get_results", {"section": "convergence"})]
        )
    if "crash" in text:
        raise RuntimeError("the model crashed")
    return AssistantMessage("It is water.")


@pytest.fixture
def fakes(monkeypatch: pytest.MonkeyPatch) -> Fakes:
    return Fakes(monkeypatch)


@pytest.fixture
def fake_model(monkeypatch: pytest.MonkeyPatch) -> FakeModel:
    model = FakeModel(model_script)
    monkeypatch.setattr(backend_module, "make_backend", lambda *args, **kwargs: model)
    return model


def _read_transcript(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _run_args(folder: Path, *extra: str) -> list[str]:
    return [
        "run",
        *FILES,
        "--out",
        str(folder),
        "--quick",
        "--no-chat",
        "--quiet",
        *extra,
    ]


# --------------------------------------------------------------------------- usage
def test_usage_errors(capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.main([]) == 2
    assert cli.main(["--help"]) == 0
    assert "run" in capsys.readouterr().out
    assert cli.main(["run", "--system", "water", "--no-chat"]) == 2
    assert "no data files" in capsys.readouterr().err
    assert cli.main(["run", *FILES, "--spec", "x.json", "--system", "water"]) == 2
    assert "not allowed with" in capsys.readouterr().err
    assert cli.main(["search"]) == 2
    assert cli.main(["search", "--download", "http://example.org/a.mint01"]) == 2
    assert "--to" in capsys.readouterr().err
    assert cli.main(["spec", FILES[0]]) == 2  # --system is required
    assert cli.main(["nonsense"]) == 2


def test_handled_errors_and_interrupt(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def boom(args: Any) -> int:
        raise RuntimeError("boom")

    monkeypatch.setattr(cli, "info_command", boom)
    assert cli.main(["info"]) == 1
    assert capsys.readouterr().err.strip() == "error: boom"

    def interrupt(args: Any) -> int:
        raise KeyboardInterrupt

    monkeypatch.setattr(cli, "info_command", interrupt)
    assert cli.main(["info"]) == 130
    assert "interrupted" in capsys.readouterr().err


# --------------------------------------------------------------------------- spec
def test_spec_writes_a_valid_specification(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    target = tmp_path / "specs" / "water.json"
    assert cli.main(["spec", *FILES, "--system", "water", "-o", str(target)]) == 0
    captured = capsys.readouterr()
    assert captured.out == "" and f"wrote {target}" in captured.err
    spec = SystemSpec.load(target)
    assert spec.validate(check_files=True) == []
    assert spec.name == "waterbox" and spec.source == "library:water"
    assert sorted(dataset.label for dataset in spec.datasets) == ["D2O", "H2O", "HDO"]
    assert all(Path(dataset.file).is_file() for dataset in spec.datasets)

    # to stdout, with a run name
    assert cli.main(["spec", *FILES, "--system", "water", "--name", "mywater"]) == 0
    printed = SystemSpec.from_json(capsys.readouterr().out)
    assert (
        printed.name == "mywater"
        and printed.to_dict()["datasets"] == spec.to_dict()["datasets"]
    )


def test_spec_errors(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.main(["spec", *FILES, "--system", "argon"]) == 1
    assert "unknown system 'argon'" in capsys.readouterr().err
    assert cli.main(["spec", *FILES, "--system", "silica"]) == 1
    assert "looks like water" in capsys.readouterr().err
    assert (
        cli.main(["spec", str(tmp_path / "missing.mint01"), "--system", "water"]) == 1
    )
    assert "not found" in capsys.readouterr().err


# --------------------------------------------------------------------------- search
def test_search_local(capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.main(["search", "water SANDALS", "--sources", "local"]) == 0
    out = capsys.readouterr().out
    assert "[local]" in out and "candidate(s) for 'water SANDALS'" in out
    for path in WATER_FILES:
        assert path.name in out and str(path) in out


def test_search_download_and_notes(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    calls: list[Any] = []

    def download(url: str, folder: Path, *, timeout: float = 60.0) -> Path:
        calls.append((url, Path(folder)))
        path = Path(folder) / "a.mint01"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# fake\n")
        return path

    def find_similar_data(
        query: str, *, sources: Any = None, limit: int = 10
    ) -> dict[str, Any]:
        calls.append((query, sources, limit))
        return {
            "query": query,
            "candidates": [
                {
                    "source": "dissolve",
                    "title": "Dissolve example: water",
                    "description": "Liquid water",
                    "url": "https://github.com/x/examples/water",
                    "files": ["H2O.mint01"],
                    "score": 1.5,
                    "extra": {
                        "downloads": {
                            "H2O.mint01": "https://raw.example.org/H2O.mint01"
                        }
                    },
                }
            ],
            "notes": ["zenodo: Zenodo unreachable: offline"],
        }

    monkeypatch.setattr(sources_module, "download", download)
    monkeypatch.setattr(sources_module, "find_similar_data", find_similar_data)
    target = tmp_path / "downloads"
    assert cli.main(
        ["search", "water", "--sources", "dissolve,zenodo", "--limit", "3",
         "--download", "https://raw.example.org/H2O.mint01", "--to", str(target)]
    ) == 0  # fmt: skip
    out = capsys.readouterr().out
    assert calls[0] == ("water", "dissolve,zenodo", 3)
    assert calls[1] == ("https://raw.example.org/H2O.mint01", target)
    assert "[dissolve] score 1.5  Dissolve example: water" in out
    assert "download: https://raw.example.org/H2O.mint01" in out
    assert "note: zenodo: Zenodo unreachable: offline" in out
    assert f"downloaded {target / 'a.mint01'}" in out


# --------------------------------------------------------------------------- info
def test_info(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(cli, "_installed", lambda name: False)  # no torch import
    monkeypatch.delenv("HF_TOKEN", raising=False)
    monkeypatch.delenv("MELLOR_BACKEND", raising=False)
    monkeypatch.delenv("MELLOR_MODEL", raising=False)
    assert cli.main(["info"]) == 0
    out = capsys.readouterr().out
    assert out.startswith(f"mellor {cli.__version__}")
    assert "transformers: not installed" in out and "hf-api: HF_TOKEN not set" in out
    assert "default backend: auto -> none" in out
    assert f"default model: {backend_module.DEFAULT_MODEL}" in out
    assert "library systems: water, silica, dmso_water, benzene, cu_btc" in out
    assert "experimental_data" in out

    monkeypatch.setenv("HF_TOKEN", "hf_dummy")
    monkeypatch.setenv("MELLOR_MODEL", "org/model")
    assert cli.main(["info"]) == 0
    out = capsys.readouterr().out
    assert "hf-api: HF_TOKEN set" in out and "auto -> hf-api" in out
    assert "default model: org/model" in out


# --------------------------------------------------------------------------- progress
def test_event_printer() -> None:
    stream = io.StringIO()
    show = cli.event_printer(stream)
    show({"type": "phase", "phase": "setup"})
    show({"type": "assistant", "content": "Let me look.", "tool_calls": []})
    show(
        {
            "type": "assistant",
            "content": "",
            "tool_calls": [{"name": "inspect_data", "arguments": {}}],
        }
    )
    show({"type": "assistant", "content": "", "tool_calls": []})
    show(
        {
            "type": "tool_call",
            "name": "build_and_run",
            "arguments": {
                "spec": {"name": "waterbox", "kind": "molecular"},
                "quick": True,
            },
        }
    )
    show(
        {
            "type": "tool_result",
            "name": "build_and_run",
            "ok": True,
            "preview": '{"ok": true}',
        }
    )
    show(
        {
            "type": "tool_result",
            "name": "validate_spec",
            "ok": False,
            "preview": '{"error": "x"}',
        }
    )
    show({"type": "progress", "line": "run 1/4: scratch"})
    show({"type": "note", "content": "a note"})
    show({"type": "fallback", "system": "water", "reason": "the model stopped"})
    show({"type": "report", "path": "/tmp/report.pdf"})
    show({"type": "mystery", "value": 1})
    lines = stream.getvalue().splitlines()
    assert lines == [
        "== setup ==",
        "[model] Let me look.",
        "[model] calls inspect_data",
        "[model] (empty reply)",
        '[tool] build_and_run({"spec": "<spec waterbox>", "quick": true})',
        '[tool] build_and_run -> ok: {"ok": true}',
        '[tool] validate_spec -> error: {"error": "x"}',
        "  run 1/4: scratch",
        "[note] a note",
        "[fallback] library system 'water': the model stopped",
        "[report] /tmp/report.pdf",
        '[mystery] {"type": "mystery", "value": 1}',
    ]
    assert len(cli._short("x" * 500)) == cli.PREVIEW_CHARS


# --------------------------------------------------------------------------- run (fakes)
def test_run_system_without_model(
    fakes: Fakes, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    folder = tmp_path / "run"
    assert cli.main(["run", *FILES, "--out", str(folder), "--system", "water", "--quick", "--backend", "scripted", "--no-chat"]) == 0  # fmt: skip
    captured = capsys.readouterr()
    assert "run  stage" in captured.out and "EP refinement" in captured.out
    assert f"Report: {folder / 'report.pdf'}" in captured.out
    assert f"Run folder: {folder}" in captured.out
    assert "Specification: library:water." in captured.out
    assert "no language model was used" in captured.out
    assert "== run ==" in captured.err and "[tool] build_and_run(" in captured.err
    assert (
        "[tool] build_and_run -> ok" in captured.err
        and "  run 1/4: scratch" in captured.err
    )
    assert "== narrative ==" in captured.err and "[report]" in captured.err

    assert fakes.run_calls == [
        {"folder": folder, "quick": True, "source": "library:water"}
    ]
    assert SystemSpec.load(folder / "spec.json").source == "library:water"
    assert (folder / "report.pdf").is_file() and (folder / "analysis.json").is_file()
    discussion = (folder / DISCUSSION_NAME).read_text()
    assert discussion.startswith(f"_{NO_MODEL_NOTE}_") and "## Runs" in discussion
    assert fakes.pdf_calls[-1]["discussion"] + "\n" == discussion
    assert fakes.pdf_calls[-1]["model_name"].startswith("none:")
    assert fakes.pdf_calls[-1]["figures"] == {
        "01_convergence": folder / "figures" / "01_convergence.png"
    }
    transcript = _read_transcript(folder / TRANSCRIPT_NAME)
    assert (
        transcript[0]["spec_given"] is True
        and transcript[0]["narrative_by_model"] is False
    )
    assert (
        transcript[1]["role"] == "spec" and transcript[1]["source"] == "library:water"
    )
    narrative = [entry for entry in transcript if entry["role"] == "narrative"]
    assert narrative == [{**narrative[0], "source": "summary", "model": "scripted"}]


def test_run_spec_file_and_quiet(
    fakes: Fakes, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    spec_path = tmp_path / "water.json"
    assert cli.main(["spec", *FILES, "--system", "water", "-o", str(spec_path)]) == 0
    capsys.readouterr()
    folder = tmp_path / "from_spec"
    # no data files on the command line: they come from the specification
    assert cli.main(["run", "--spec", str(spec_path), "--out", str(folder), "--name", "mine", "--backend", "scripted", "--no-chat", "--quiet"]) == 0  # fmt: skip
    captured = capsys.readouterr()
    assert captured.err == ""
    assert "Specification: user." in captured.out
    spec = SystemSpec.load(folder / "spec.json")
    assert spec.source == "user" and spec.name == "mine"
    assert fakes.run_calls[-1]["source"] == "user"
    assert (folder / "report.pdf").is_file()

    # relative data-file names are resolved next to the specification file
    data = json.loads(spec_path.read_text())
    for dataset in data["datasets"]:
        dataset["file"] = Path(dataset["file"]).name
    (tmp_path / "relative.json").write_text(json.dumps(data))
    for path in WATER_FILES:
        (tmp_path / path.name).write_bytes(path.read_bytes())
    assert cli.main(["run", "--spec", str(tmp_path / "relative.json"), "--out", str(tmp_path / "rel"), "--backend", "scripted", "--no-chat", "--quiet"]) == 0  # fmt: skip
    loaded = SystemSpec.load(tmp_path / "rel" / "spec.json")
    assert [Path(d.file).parent for d in loaded.datasets] == [tmp_path] * 3


def test_run_invalid_spec_file(
    fakes: Fakes, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"name": "x", "schedule": {"ereq": 3}}))
    assert cli.main(["run", *FILES, "--spec", str(bad), "--out", str(tmp_path / "r"), "--backend", "scripted", "--no-chat", "--quiet"]) == 1  # fmt: skip
    err = capsys.readouterr().err
    assert "error: invalid system specification" in err and "no atom types" in err
    assert fakes.run_calls == []
    assert (
        cli.main(["run", *FILES, "--spec", str(tmp_path / "none.json"), "--no-chat"])
        == 1
    )
    assert "specification file not found" in capsys.readouterr().err
    (tmp_path / "notjson.json").write_text("{")
    assert (
        cli.main(["run", *FILES, "--spec", str(tmp_path / "notjson.json"), "--no-chat"])
        == 1
    )
    assert "not valid JSON" in capsys.readouterr().err


def test_run_auto_spec_without_model(
    fakes: Fakes, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    folder = tmp_path / "auto"
    assert cli.main(["run", *FILES, "--out", str(folder), "--backend", "scripted", "--no-chat", "--quiet"]) == 0  # fmt: skip
    captured = capsys.readouterr()
    assert "note: no language model; Chose the library system 'water'" in captured.err
    assert (
        fakes.run_calls[-1]["source"] == "library:water"
        and fakes.run_calls[-1]["quick"] is False
    )


def test_run_default_work_dir_and_overwrite(
    fakes: Fakes,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.chdir(tmp_path)
    assert cli.main(["run", *FILES, "--system", "water", "--backend", "scripted", "--no-chat", "--quiet"]) == 0  # fmt: skip
    folder = tmp_path / "work" / "waterbox"
    assert fakes.run_calls[-1]["folder"].resolve() == folder.resolve()
    assert (folder / "manifest.json").is_file()
    capsys.readouterr()
    # a second run into the same folder is refused without --overwrite
    assert cli.main(["run", *FILES, "--system", "water", "--backend", "scripted", "--no-chat", "--quiet"]) == 1  # fmt: skip
    assert (
        "already holds a run" in capsys.readouterr().err and len(fakes.run_calls) == 1
    )
    assert cli.main(["run", *FILES, "--system", "water", "--backend", "scripted", "--no-chat", "--quiet", "--overwrite"]) == 0  # fmt: skip
    assert len(fakes.run_calls) == 2
    assert cli.main(["run", *FILES, "--system", "water", "--name", "second", "--backend", "scripted", "--no-chat", "--quiet"]) == 0  # fmt: skip
    assert fakes.run_calls[-1]["folder"] == tmp_path / "work" / "second"


def test_run_auto_spec_needs_a_recognisable_title(
    fakes: Fakes, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    unknown = tmp_path / "mystery.mint01"
    unknown.write_text(
        "# run1.mdcs01\n# a sample nobody knows\n0.1 0.5 0.01\n0.2 0.4 0.01\n"
    )
    assert cli.main(["run", str(unknown), "--backend", "scripted", "--no-chat", "--quiet"]) == 1  # fmt: skip
    err = capsys.readouterr().err
    assert "no built-in system matches" in err and "--system KEY" in err
    assert fakes.run_calls == []


def test_run_missing_data_and_unknown_system(
    fakes: Fakes, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert (
        cli.main(
            ["run", str(tmp_path / "gone.mint01"), "--system", "water", "--no-chat"]
        )
        == 1
    )
    assert "data file(s) not found" in capsys.readouterr().err
    assert (
        cli.main(
            ["run", *FILES, "--system", "argon", "--no-chat", "--backend", "scripted"]
        )
        == 1
    )
    assert "unknown system 'argon'" in capsys.readouterr().err
    assert fakes.run_calls == []


def test_run_without_any_backend(
    fakes: Fakes,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def no_backend(*args: Any, **kwargs: Any) -> Any:
        raise RuntimeError("no language-model backend is available")

    monkeypatch.setattr(backend_module, "make_backend", no_backend)
    folder = tmp_path / "nomodel"
    assert (
        cli.main(
            ["run", *FILES, "--out", str(folder), "--system", "water", "--no-chat"]
        )
        == 0
    )
    captured = capsys.readouterr()
    assert (
        "note: no language model (no language-model backend is available)"
        in captured.err
    )
    assert fakes.pdf_calls[-1]["model_name"].startswith("none:")
    # without a specification the model is needed
    assert cli.main(["run", *FILES, "--out", str(tmp_path / "x"), "--no-chat"]) == 1
    assert "error: no language-model backend is available" in capsys.readouterr().err
    assert len(fakes.run_calls) == 1


def test_run_with_a_model(
    fakes: Fakes,
    fake_model: FakeModel,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    folder = tmp_path / "model"
    assert cli.main(["run", *FILES, "--out", str(folder), "--quick", "--no-chat", "--prompt", "Water at 25 C.", "--name", "named"]) == 0  # fmt: skip
    captured = capsys.readouterr()
    assert "== setup ==" in captured.err and "[model] Looking." in captured.err
    assert (
        "[tool] inspect_data({})" in captured.err
        and "[tool] library_spec(" in captured.err
    )
    assert "Specification: library:water." in captured.out
    assert fakes.pdf_calls[-1]["model_name"] == "fake-model"
    assert fakes.pdf_calls[-1]["discussion"] == NARRATIVE.strip()
    assert (folder / DISCUSSION_NAME).read_text() == NARRATIVE.strip() + "\n"
    setup_task = fake_model.requests[0][-1]["content"]
    assert (
        "User request: Water at 25 C." in setup_task
        and "Use 'named' as the run name" in setup_task
    )
    transcript = _read_transcript(folder / TRANSCRIPT_NAME)
    assert [entry for entry in transcript if entry["role"] == "narrative"][0][
        "source"
    ] == "model"

    # --system with a model: the model only writes the discussion
    folder2 = tmp_path / "model_system"
    assert cli.main(["run", *FILES, "--out", str(folder2), "--system", "water", "--no-chat", "--quiet"]) == 0  # fmt: skip
    assert fakes.pdf_calls[-1]["model_name"] == "fake-model"
    assert fakes.pdf_calls[-1]["discussion"] == NARRATIVE.strip()
    request = fake_model.requests[-1]
    assert (
        request[1]["role"] == "user"
        and "The specification was given (library:water)" in request[1]["content"]
    )
    assert 'Write the "Discussion" section' in request[-1]["content"]


def test_run_narrative_survives_a_model_failure(
    fakes: Fakes,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def crashing(messages: list[dict[str, Any]]) -> AssistantMessage:
        raise RuntimeError("the model crashed")

    model = FakeModel(crashing)
    monkeypatch.setattr(backend_module, "make_backend", lambda *args, **kwargs: model)
    folder = tmp_path / "crash"
    assert cli.main(["run", *FILES, "--out", str(folder), "--system", "water", "--no-chat", "--quiet"]) == 0  # fmt: skip
    out = capsys.readouterr().out
    assert (
        "the model failed while writing the discussion (RuntimeError: the model crashed)"
        in out
    )
    assert (folder / "report.pdf").is_file()
    assert fakes.pdf_calls[-1]["discussion"].startswith(
        "_The language model did not produce a discussion"
    )
    assert fakes.pdf_calls[-1]["model_name"] == "fake-model"


# --------------------------------------------------------------------------- chat
def _fake_run_folder(fakes: Fakes, tmp_path: Path) -> Path:
    folder = tmp_path / "run"
    assert cli.main(["run", *FILES, "--out", str(folder), "--system", "water", "--backend", "scripted", "--no-chat", "--quiet"]) == 0  # fmt: skip
    return folder


def test_chat_from_a_pipe(
    fakes: Fakes,
    fake_model: FakeModel,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    folder = _fake_run_folder(fakes, tmp_path)
    capsys.readouterr()
    monkeypatch.setattr(
        sys, "stdin", io.StringIO("verdict?\n\nwhat is it?\n/quit\nignored\n")
    )
    seen = len(fake_model.requests)
    assert cli.main(["chat", str(folder), "--quiet"]) == 0
    out = capsys.readouterr().out
    assert "you> verdict?\nMellor> The verdict: a fake verdict.\n" in out
    assert "you> what is it?\nMellor> It is water.\n" in out
    assert "ignored" not in out
    first = fake_model.requests[seen][-1]["content"]
    assert "Summary of the run" in first and "User: verdict?" in first
    assert cli.main(["chat", str(tmp_path / "nowhere"), "--quiet"]) == 1
    assert "does not exist" in capsys.readouterr().err


def test_chat_interactive_loop(
    fakes: Fakes,
    fake_model: FakeModel,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    folder = _fake_run_folder(fakes, tmp_path)
    capsys.readouterr()
    answers = iter(["", "what is it?", "crash", "verdict?", "/quit", "never"])
    monkeypatch.setattr(builtins, "input", lambda prompt="": next(answers))
    monkeypatch.setattr(cli, "_is_terminal", lambda stream: True)
    assert cli.main(["chat", str(folder), "--quiet"]) == 0
    captured = capsys.readouterr()
    assert "Ask Mellor about the run (/quit to leave)." in captured.out
    assert (
        "Mellor> It is water." in captured.out
        and "Mellor> The verdict: a fake verdict." in captured.out
    )
    assert "error: RuntimeError: the model crashed" in captured.err
    assert "never" not in captured.out

    # end of input leaves the loop too
    monkeypatch.setattr(
        builtins, "input", lambda prompt="": (_ for _ in ()).throw(EOFError())
    )
    assert cli.main(["chat", str(folder), "--quiet"]) == 0


def test_run_drops_into_the_chat_on_a_terminal(
    fakes: Fakes,
    fake_model: FakeModel,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    answers = iter(["what is it?", "/quit"])
    monkeypatch.setattr(builtins, "input", lambda prompt="": next(answers))
    monkeypatch.setattr(cli, "_is_terminal", lambda stream: True)
    assert cli.main(["run", *FILES, "--out", str(tmp_path / "chatty"), "--system", "water", "--quiet"]) == 0  # fmt: skip
    out = capsys.readouterr().out
    assert "Mellor> It is water." in out
    # --no-chat, and no model, both skip the chat
    with contextlib.suppress(StopIteration):
        assert cli.main(["run", *FILES, "--out", str(tmp_path / "quiet"), "--system", "water", "--quiet", "--no-chat"]) == 0  # fmt: skip
    assert "Mellor>" not in capsys.readouterr().out


# --------------------------------------------------------------------------- report (fakes)
def test_report_reuses_the_discussion(
    fakes: Fakes,
    fake_model: FakeModel,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    folder = tmp_path / "model_run"
    assert cli.main(["run", *FILES, "--out", str(folder), "--system", "water", "--no-chat", "--quiet"]) == 0  # fmt: skip
    (folder / "report.pdf").unlink()
    (folder / "analysis.json").unlink()
    capsys.readouterr()
    assert cli.main(["report", str(folder)]) == 0
    out = capsys.readouterr().out
    assert f"Report: {folder / 'report.pdf'}" in out and "run  stage" in out
    assert (folder / "report.pdf").is_file() and (folder / "analysis.json").is_file()
    call = fakes.pdf_calls[-1]
    assert (
        call["discussion"] == NARRATIVE.strip() + "\n"
        and call["model_name"] == "fake-model"
    )
    assert call["figures"] == {
        "01_convergence": folder / "figures" / "01_convergence.png"
    }
    assert json.loads((folder / "analysis.json").read_text())["figures"] == {
        "01_convergence": str(folder / "figures" / "01_convergence.png")
    }

    # a discussion given on the command line replaces it
    given = tmp_path / "mine.md"
    given.write_text("## My discussion\n\nWritten by hand.\n")
    assert cli.main(["report", str(folder), "--discussion", str(given)]) == 0
    call = fakes.pdf_calls[-1]
    assert call["discussion"] == given.read_text() and call["model_name"].startswith(
        "none:"
    )
    assert (folder / DISCUSSION_NAME).read_text() == given.read_text()
    assert (
        cli.main(["report", str(folder), "--discussion", str(tmp_path / "no.md")]) == 1
    )
    assert "discussion file not found" in capsys.readouterr().err


def test_report_without_a_discussion(
    fakes: Fakes, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    folder = _fake_run_folder(fakes, tmp_path)
    (folder / DISCUSSION_NAME).unlink()
    (folder / TRANSCRIPT_NAME).unlink()
    assert cli.main(["report", str(folder)]) == 0
    call = fakes.pdf_calls[-1]
    assert (
        call["discussion"].startswith(f"_{NO_MODEL_NOTE}_")
        and "## Runs" in call["discussion"]
    )
    assert call["model_name"].startswith("none:")
    assert (folder / DISCUSSION_NAME).read_text().startswith(f"_{NO_MODEL_NOTE}_")
    # a no-model run keeps its label when the transcript says so
    folder2 = tmp_path / "summary_run"
    assert cli.main(["run", *FILES, "--out", str(folder2), "--system", "water", "--backend", "scripted", "--no-chat", "--quiet"]) == 0  # fmt: skip
    assert cli.main(["report", str(folder2)]) == 0
    assert fakes.pdf_calls[-1]["model_name"].startswith("none:")
    assert cli.main(["report", str(tmp_path / "missing")]) == 1
    assert "does not exist" in capsys.readouterr().err


# --------------------------------------------------------------------------- slow: the real thing
@pytest.fixture(scope="module")
def cli_run(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """``mellor run --system water --quick --backend scripted --no-chat`` for real."""
    folder = tmp_path_factory.mktemp("cli") / "water"
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        status = cli.main(
            ["run", *FILES, "--out", str(folder), "--system", "water", "--quick",
             "--backend", "scripted", "--no-chat"]
        )  # fmt: skip
    return {
        "folder": folder,
        "status": status,
        "out": out.getvalue(),
        "err": err.getvalue(),
    }


def _pdf_text(path: Path) -> str:
    return "\n".join(page.extract_text() for page in PdfReader(str(path)).pages)


@pytest.mark.slow
def test_cli_run_water_quick(cli_run: dict[str, Any]) -> None:
    folder = cli_run["folder"]
    assert cli_run["status"] == 0, cli_run["err"][-2000:]
    out = cli_run["out"]
    assert (
        "run  stage" in out
        and "EP refinement" in out
        and f"Report: {folder / 'report.pdf'}" in out
    )
    assert "Specification: library:water." in out
    err = cli_run["err"]
    assert "[tool] build_and_run(" in err and "run 4/4" in err and "[report]" in err
    for name in ("report.pdf", "analysis.json", "manifest.json", "run_waterbox.py", "spec.json",
                 "discussion.md", TRANSCRIPT_NAME, "figures/01_convergence.png", "trajectory/meta.json"):  # fmt: skip
        assert (folder / name).is_file(), name
    text = _pdf_text(folder / "report.pdf")
    assert (
        "Not written by a language model" in text and "without a language model" in text
    )
    for heading in report_module.SECTION_HEADINGS:
        assert heading in text, heading
    assert "quick schedule" in text


@pytest.mark.slow
def test_cli_report_and_chat_on_the_real_run(
    cli_run: dict[str, Any],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    folder = cli_run["folder"]
    assert cli_run["status"] == 0
    before = (folder / "report.pdf").stat().st_size
    assert cli.main(["report", str(folder)]) == 0
    out = capsys.readouterr().out
    assert f"Report: {folder / 'report.pdf'}" in out and "Figures:" in out
    assert (folder / "report.pdf").stat().st_size > 0.5 * before
    assert "Not written by a language model" in _pdf_text(folder / "report.pdf")

    given = tmp_path / "discussion.md"
    given.write_text(
        "## Discussion\n\nA handwritten discussion of the quick water run.\n"
    )
    assert cli.main(["report", str(folder), "--discussion", str(given)]) == 0
    capsys.readouterr()
    assert "handwritten" in _pdf_text(folder / "report.pdf")

    # chat on the loaded run: the model reads the real analysis through get_results
    def answers(messages: list[dict[str, Any]]) -> AssistantMessage:
        last = messages[-1]
        if last["role"] == "tool":
            shells = json.loads(last["content"])["shells"]
            return AssistantMessage(
                f"The O-H first peak is at {shells['OW-HW']['r_peak']:.2f} A."
            )
        return AssistantMessage("", [ToolCall("get_results", {"section": "shells"})])

    model = FakeModel(answers)
    monkeypatch.setattr(backend_module, "make_backend", lambda *args, **kwargs: model)
    monkeypatch.setattr(sys, "stdin", io.StringIO("Where is the O-H peak?\n"))
    assert cli.main(["chat", str(folder), "--quiet"]) == 0
    out = capsys.readouterr().out
    analysis = analysis_module.load_analysis(folder / "analysis.json")
    expected = analysis["shells"]["OW-HW"]["r_peak"]
    assert f"Mellor> The O-H first peak is at {expected:.2f} A." in out
