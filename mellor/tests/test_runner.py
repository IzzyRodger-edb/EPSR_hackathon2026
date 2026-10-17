"""``mellor.runner``: the failure paths with a stand-in script (fast) and the real
water run with the quick schedule (slow, about a minute)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
from example_run import WATER_FILES
from specs import water_spec

from mellor import analysis, figures, library, script
from mellor.data import read_mint01
from mellor.runner import (
    RUN_PRODUCTS,
    RunError,
    RunResult,
    format_runs_table,
    load_result,
    run_spec,
)
from mellor.spec import SpecError

FAILING = """\
import sys
print("hello from the stand-in", flush=True)
print("arguments:", sys.argv[1:], flush=True)
print("boom", file=sys.stderr, flush=True)
sys.exit(3)
"""

SUCCEEDING = """\
import json
import sys
from pathlib import Path

print("stand-in run", flush=True)
summary = {
    "r_factor": 0.02, "r_factor_per_dataset": [0.03, 0.02, 0.01], "energy": -1.5,
    "pressure": 2.0, "ep_amplitude": 0.0, "energy_inter": -1.5, "energy_intra": 0.1,
    "chisq": 0.5,
}
manifest = {
    "name": "waterbox", "quick": "--quick" in sys.argv, "record": "--no-record" not in sys.argv,
    "complete": True, "duration_s": 0.25, "epsr_folder": "epsr",
    "trajectory": None if "--no-record" in sys.argv else "trajectory",
    "runs": [{"run": 1, "stage": "scratch", "niter": 1, "ereq": 0.0, "duration_s": 0.1, "summary": summary}],
    "final": summary,
}
Path("manifest.json").write_text(json.dumps(manifest))
print("done", flush=True)
"""

SLEEPING = """\
import time
print("started", flush=True)
time.sleep(30)
"""

NO_MANIFEST = """\
print("nothing written")
"""


def stand_in(body: str):
    """A ``write_script`` replacement that writes ``body`` as the run script."""

    def write(spec, inputs, folder):
        path = Path(folder) / script.script_name(spec)
        path.write_text(body)
        return path

    return write


# --------------------------------------------------------------------------- fast
def test_run_spec_rejects_a_missing_data_file(tmp_path: Path) -> None:
    spec = water_spec()
    spec.datasets[0].file = str(tmp_path / "missing.mint01")
    with pytest.raises(SpecError) as info:
        run_spec(spec, tmp_path / "run")
    assert "missing.mint01" in str(info.value)
    assert not (tmp_path / "run").exists()


def test_run_spec_failing_script(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(script, "write_script", stand_in(FAILING))
    spec = water_spec()
    folder = tmp_path / "run"
    lines: list[str] = []
    with pytest.raises(RunError) as info:
        run_spec(spec, folder, quick=True, progress=lines.append)
    error = info.value
    assert "run_waterbox.py failed with exit status 3" in str(error)
    assert "hello from the stand-in" in error.log and "boom" in error.log
    assert "arguments: ['--out', '.', '--quick']" in error.log
    assert (folder / "run.log").read_text() == error.log + "\n"
    assert lines[0].startswith("$ ") and lines[0].endswith(f"(in {folder})")
    assert sys.executable in lines[0] and " -I " in lines[0]
    assert lines[1:] == error.log.splitlines()[1:]
    # the folder holds everything the script needed
    assert (folder / "spec.json").is_file() and (folder / "run_waterbox.py").is_file()
    assert sorted(path.name for path in (folder / "data").iterdir()) == sorted(
        path.name for path in WATER_FILES
    )
    for original in WATER_FILES:
        assert (folder / "data" / original.name).read_bytes() == original.read_bytes()
    assert not (folder / "manifest.json").exists()


def test_run_spec_success_overwrite_and_load(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(script, "write_script", stand_in(SUCCEEDING))
    spec = water_spec()
    folder = tmp_path / "run"
    folder.mkdir()
    # leftovers of an earlier run, and files that must survive
    (folder / "epsr").mkdir()
    (folder / "epsr" / "stale.txt").write_text("old")
    (folder / "figures").mkdir()
    (folder / "figures" / "01_convergence.png").write_bytes(b"old")
    for name in ("manifest.json", "report.pdf", "analysis.json", "run.log"):
        (folder / name).write_text("old")
    (folder / "notes.txt").write_text("mine")
    (folder / "data").mkdir()
    (folder / "data" / "extra.mint01").write_text("# extra\n")

    with pytest.raises(FileExistsError, match="not empty"):
        run_spec(spec, folder)

    result = run_spec(spec, folder, quick=True, overwrite=True)
    assert result.folder == folder.resolve() and result.quick is True
    assert result.duration_s > 0
    assert result.manifest["quick"] is True and result.runs[0]["stage"] == "scratch"
    assert result.final_summary["r_factor"] == 0.02
    assert result.trajectory_path == folder / "trajectory" and not result.has_trajectory
    assert result.epsr_folder == folder / "epsr" and result.complete
    for name in RUN_PRODUCTS:
        assert (
            name in ("manifest.json", "run.log") or not (folder / name).exists()
        ), name
    assert (folder / "notes.txt").read_text() == "mine"
    assert (folder / "data" / "extra.mint01").read_text() == "# extra\n"
    assert "stand-in run" in (folder / "run.log").read_text()
    assert result.script_path == folder / "run_waterbox.py"

    loaded = load_result(folder)
    assert loaded.manifest == result.manifest
    assert loaded.spec.to_dict() == spec.to_dict()
    assert loaded.quick is True and loaded.duration_s == 0.25
    assert loaded.folder == folder.resolve()

    # a specification that points at the copies in data/ runs again (no self-copy)
    for dataset in spec.datasets:
        dataset.file = str(folder / "data" / dataset.file_name)
    again = run_spec(spec, folder, record=False, overwrite=True)
    assert again.manifest["record"] is False and again.trajectory_path is None
    assert "--no-record" in (folder / "run.log").read_text().splitlines()[0]


def test_run_spec_timeout(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(script, "write_script", stand_in(SLEEPING))
    with pytest.raises(RunError, match="did not finish within 1.0 s") as info:
        run_spec(water_spec(), tmp_path / "run", timeout_s=1.0)
    assert "started" in info.value.log and "killed after 1.0 s" in info.value.log


def test_run_spec_without_manifest(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(script, "write_script", stand_in(NO_MANIFEST))
    with pytest.raises(RunError, match="wrote no manifest.json") as info:
        run_spec(water_spec(), tmp_path / "run")
    assert "nothing written" in info.value.log


def test_progress_callback_failure_does_not_stop_the_run(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(script, "write_script", stand_in(SUCCEEDING))
    calls: list[str] = []

    def broken(line: str) -> None:
        calls.append(line)
        if len(calls) == 2:
            raise ValueError("listener gone")

    result = run_spec(water_spec(), tmp_path / "run", progress=broken)
    assert result.runs
    assert len(calls) == 2
    log = (tmp_path / "run" / "run.log").read_text()
    assert "progress callback failed (ValueError: listener gone)" in log
    assert "done" in log.splitlines()[-1]


def test_load_result_errors(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="does not exist"):
        load_result(tmp_path / "nowhere")
    folder = tmp_path / "run"
    folder.mkdir()
    with pytest.raises(FileNotFoundError, match="spec.json is missing"):
        load_result(folder)
    water_spec().save(folder / "spec.json")
    with pytest.raises(FileNotFoundError, match="manifest.json is missing"):
        load_result(folder)
    (folder / "manifest.json").write_text("{broken")
    with pytest.raises(ValueError, match="not valid JSON"):
        load_result(folder)
    (folder / "manifest.json").write_text("[]")
    with pytest.raises(ValueError, match="JSON object"):
        load_result(folder)
    (folder / "manifest.json").write_text("{}")
    result = load_result(folder)
    assert result.runs == [] and result.quick is False and result.duration_s == 0.0
    assert result.trajectory_path == folder / "trajectory"


def test_run_result_properties_and_table(tmp_path: Path) -> None:
    manifest = {
        "epsr_folder": "epsr",
        "trajectory": None,
        "complete": True,
        "runs": [
            {
                "run": 1,
                "stage": "scratch",
                "niter": 1,
                "ereq": 0.0,
                "duration_s": 1.25,
                "summary": {
                    "r_factor": 0.021,
                    "energy": 1.36e9,
                    "pressure": 3.56e9,
                    "ep_amplitude": 0.0,
                },
            },
            {
                "run": 2,
                "stage": "EP refinement",
                "niter": 100,
                "ereq": 3.0,
                "duration_s": 150.0,
                "summary": {
                    "r_factor": 2.22e-4,
                    "energy": -44.89,
                    "pressure": 3.87,
                    "ep_amplitude": 2.9999,
                },
            },
        ],
        "final": {"r_factor": 2.22e-4},
    }
    result = RunResult(
        folder=str(tmp_path), spec=water_spec(), manifest=manifest, quick=True
    )
    assert result.folder == tmp_path
    assert result.spec_path == tmp_path / "spec.json"
    assert result.script_path == tmp_path / "run_waterbox.py"
    assert result.epsr_folder == tmp_path / "epsr"
    assert result.trajectory_path is None and not result.has_trajectory
    assert result.log_path == tmp_path / "run.log"
    assert result.data_folder == tmp_path / "data"
    assert len(result.runs) == 2 and result.final_summary == {"r_factor": 2.22e-4}
    assert result.complete

    table = result.table()
    lines = table.splitlines()
    assert lines[0].split() == [
        "run",
        "stage",
        "niter",
        "ereq",
        "R-factor",
        "energy",
        "pressure",
        "EP",
        "ampl",
        "time",
    ]
    assert lines[1].split() == [
        "1",
        "scratch",
        "1",
        "0",
        "2.100e-02",
        "1.36e+09",
        "3.56e+09",
        "0.000",
        "1.2s",
    ]
    assert lines[2].split() == [
        "2",
        "EP",
        "refinement",
        "100",
        "3",
        "2.220e-04",
        "-44.89",
        "3.87",
        "3.000",
        "150.0s",
    ]
    assert format_runs_table([]) == lines[0] + "\n(no runs)"
    sparse = format_runs_table([{"run": 1, "stage": "scratch"}])
    assert sparse.splitlines()[1].split() == [
        "1",
        "scratch",
        "-",
        "-",
        "-",
        "-",
        "-",
        "-",
        "-s",
    ]


# --------------------------------------------------------------------------- slow: the real thing
@pytest.fixture(scope="module")
def water_quick_run(tmp_path_factory: pytest.TempPathFactory):
    """``run_spec`` of the library's water system with the quick schedule (31 iterations)."""
    spec = library.make_spec("water", [read_mint01(path) for path in WATER_FILES])
    folder = tmp_path_factory.mktemp("water_quick") / "run"
    lines: list[str] = []
    result = run_spec(spec, folder, quick=True, progress=lines.append)
    return result, lines


@pytest.mark.slow
def test_run_spec_water_quick(water_quick_run) -> None:
    result, lines = water_quick_run
    folder = result.folder
    assert result.quick is True and result.duration_s > 0 and result.complete
    for name in (
        "spec.json",
        "run_waterbox.py",
        "data",
        "epsr",
        "trajectory",
        "manifest.json",
        "run.log",
    ):
        assert (folder / name).exists(), name
    assert sorted(path.name for path in (folder / "data").iterdir()) == sorted(
        path.name for path in WATER_FILES
    )
    manifest = result.manifest
    assert manifest["name"] == "waterbox" and manifest["kind"] == "molecular"
    assert manifest["quick"] is True and manifest["complete"] is True
    assert manifest["epsrvis_case"] == "water"
    assert [run["stage"] for run in result.runs] == [
        "scratch",
        "reference potential",
        "EP refinement",
        "EP refinement",
    ]
    assert [run["niter"] for run in result.runs] == [1, 10, 10, 10]
    assert result.final_summary == result.runs[-1]["summary"]
    assert 0 < result.final_summary["r_factor"] < result.runs[0]["summary"]["r_factor"]
    assert result.has_trajectory and result.trajectory_path == folder / "trajectory"
    assert (folder / "epsr" / "waterbox.EPSR.u01").is_file()
    assert (folder / "epsr" / "waterbox.pcof").is_file()
    assert json.loads((folder / "spec.json").read_text()) == result.spec.to_dict()

    log = (folder / "run.log").read_text()
    assert log.startswith("$ ") and "--quick" in log.splitlines()[0]
    assert "run 4/4" in log and "done in" in log
    assert lines[0] == log.splitlines()[0]
    assert any("run 4/4" in line for line in lines)
    assert "\n".join(lines) + "\n" == log

    loaded = load_result(folder)
    assert loaded.manifest == manifest
    assert loaded.spec.to_dict() == result.spec.to_dict()
    assert loaded.quick is True
    assert loaded.duration_s == pytest.approx(manifest["duration_s"])
    table = result.table()
    assert "EP refinement" in table and table.splitlines()[0].startswith("run")
    assert len(table.splitlines()) == 5


@pytest.mark.slow
def test_analysis_of_the_quick_run(water_quick_run) -> None:
    result, _ = water_quick_run
    out = analysis.analyse(result)
    assert (
        out["schedule"]["total_iterations"] == 31 and out["schedule"]["quick"] is True
    )
    assert out["schedule"]["recorded_iterations"] == 31
    assert [run["last_iteration"] for run in out["runs"]] == [1, 11, 21, 31]
    assert out["runs"][-1]["r_factor"] == pytest.approx(
        result.final_summary["r_factor"], rel=1e-6
    )
    assert out["runs"][-1]["r_factor_trajectory"] == pytest.approx(
        out["runs"][-1]["r_factor"], rel=1e-4
    )
    assert set(out["fit"]) == {"D2O", "HDO", "H2O"}
    assert out["fit"]["D2O"]["r2"] is not None and out["fit"]["D2O"]["n_points"] > 100
    assert set(out["shells"]) == {"OW-OW", "OW-HW", "HW-HW"}
    assert out["shells"]["OW-HW"]["n_coord_reverse"] == pytest.approx(
        out["shells"]["OW-HW"]["n_coord"] * 0.5
    )
    assert any("quick schedule" in flag for flag in out["flags"])
    assert out["convergence"]["verdict"]
    text = analysis.analysis_markdown(out)
    assert "## Runs" in text and "| 4 | EP refinement | 10 | 3 |" in text
    path = analysis.save_analysis(out, result.folder / "analysis.json")
    assert analysis.load_analysis(path) == out


@pytest.mark.slow
def test_figures_of_the_quick_run(water_quick_run) -> None:
    result, _ = water_quick_run
    written = figures.make_figures(result)
    assert {
        "00_data",
        "01_convergence",
        "02_fq_fit",
        "05_partials",
        "07_first_shell",
        "08_local_structure",
        "10_final_outputs",
        "11_fit_final",
        "summary.csv",
        "summary.md",
    } <= set(written)
    for name, path in written.items():
        assert path.is_file() and path.stat().st_size > 0, name
