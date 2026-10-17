"""``mellor.analysis`` on the toolkit's example water run (300 iterations, 4 runs)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from example_run import EXAMPLE, build_example_result, example_available

from mellor import analysis
from mellor.runner import RunResult

pytestmark = pytest.mark.skipif(
    not example_available(), reason="toolkit example data or water data not found"
)


@pytest.fixture(scope="module")
def example_result(tmp_path_factory: pytest.TempPathFactory) -> RunResult:
    return build_example_result(tmp_path_factory.mktemp("example_run"))


@pytest.fixture(scope="module")
def example_analysis(example_result: RunResult) -> dict:
    return analysis.analyse(example_result)


@pytest.fixture(scope="module")
def trajectory():
    from epsrvis import Trajectory

    return Trajectory.load(EXAMPLE / "trajectory")


# --------------------------------------------------------------------------- sections
def test_analysis_sections(example_analysis: dict) -> None:
    assert list(example_analysis) == [
        "system",
        "schedule",
        "runs",
        "fit",
        "convergence",
        "shells",
        "flags",
        "figures",
    ]
    assert example_analysis["figures"] == {}
    system = example_analysis["system"]
    assert system["name"] == "waterbox" and system["kind"] == "molecular"
    assert system["source"] == "library:water"
    assert system["n_atoms"] == 3000 and system["n_molecules"] == 1000
    assert system["box_side"] == pytest.approx(31.072, abs=1e-3)
    assert system["density"] == pytest.approx(0.1, rel=1e-5)
    assert 0.9 < system["mass_density_g_cm3"] < 1.2
    assert system["temperature"] == 300.0
    assert [t["label"] for t in system["atom_types"]] == ["OW", "HW"]
    assert [t["count"] for t in system["atom_types"]] == [1000, 2000]
    assert system["atom_types"][0]["fraction"] == pytest.approx(1 / 3)
    assert system["atom_types"][0]["epsilon"] == 0.6503
    assert system["pairs"] == ["OW-OW", "OW-HW", "HW-HW"]
    assert [d["label"] for d in system["datasets"]] == ["D2O", "HDO", "H2O"]
    assert system["datasets"][0]["isotopes"] == {"HW": [[2, 1.0]]}
    assert system["epsrvis_case"] == "water"
    assert system["molecules"] == [{"name": "Water", "count": 1000, "atoms": 3}]

    schedule = example_analysis["schedule"]
    assert schedule["total_iterations"] == 300 and schedule["quick"] is False
    assert schedule["recorded"] is True and schedule["recorded_iterations"] == 300
    assert schedule["complete"] is True and schedule["ereq"] == 3.0
    assert [r["niter"] for r in schedule["runs"]] == [1, 99, 100, 100]
    assert schedule["setup_steps"][0] == "makemole Water"


def test_runs_match_the_manifest_and_the_trajectory(
    example_analysis: dict, trajectory
) -> None:
    runs = example_analysis["runs"]
    assert [r["run"] for r in runs] == [1, 2, 3, 4]
    assert [r["stage"] for r in runs] == [
        "scratch",
        "reference potential",
        "EP refinement",
        "EP refinement",
    ]
    assert [r["ereq"] for r in runs] == [0.0, 0.0, 3.0, 3.0]
    assert [r["last_iteration"] for r in runs] == [1, 100, 200, 300]
    last = runs[-1]
    assert last["r_factor"] == pytest.approx(
        float(trajectory.history["rfactot"][-1]), rel=1e-6
    )
    assert 2.0e-4 < last["r_factor"] < 2.5e-4
    assert last["r_factor_trajectory"] == pytest.approx(last["r_factor"], rel=1e-4)
    assert set(last["r_factor_per_dataset"]) == {"D2O", "HDO", "H2O"}
    assert last["r_factor_per_dataset"]["D2O"] == pytest.approx(2.127e-4, rel=1e-2)
    assert last["r2"]["D2O"] == pytest.approx(0.973, abs=0.005)
    assert last["r2"]["HDO"] == pytest.approx(0.846, abs=0.005)
    assert last["chi2"]["H2O"] > 0
    assert last["ep_amplitude"] == pytest.approx(3.0, abs=0.01)
    assert last["pressure"] == pytest.approx(3.87, abs=0.01)
    assert last["energy"] == pytest.approx(-44.89, abs=0.01)
    assert last["energy_intra"] == pytest.approx(3.776, abs=0.01)
    assert last["duration_s"] == 150.0
    assert runs[0]["r_factor"] == pytest.approx(0.021, abs=1e-3)
    assert runs[1]["r_factor"] == pytest.approx(6.14e-4, rel=1e-2)
    assert runs[0]["r_factor"] > runs[1]["r_factor"] > last["r_factor"]


def test_fit_section(example_analysis: dict) -> None:
    fit = example_analysis["fit"]
    assert list(fit) == ["D2O", "HDO", "H2O"]
    d2o = fit["D2O"]
    assert d2o["file"] == "SLS18502-D2O.mint01"
    assert d2o["trajectory_dataset"] == "SLS18502-D2O"
    assert d2o["initial_r_factor"] == pytest.approx(5.81e-2, rel=1e-2)
    assert d2o["reference_r_factor"] == pytest.approx(1.21e-3, rel=1e-2)
    assert d2o["final_r_factor"] == pytest.approx(2.13e-4, rel=1e-2)
    assert d2o["improvement_factor"] == pytest.approx(
        d2o["initial_r_factor"] / d2o["final_r_factor"]
    )
    assert d2o["ep_improvement_factor"] == pytest.approx(
        d2o["reference_r_factor"] / d2o["final_r_factor"]
    )
    assert d2o["r2"] == pytest.approx(0.973, abs=0.005)
    assert d2o["r_factor_trajectory"] == pytest.approx(d2o["final_r_factor"], rel=1e-4)
    assert d2o["q_min"] == pytest.approx(0.2, abs=1e-6)
    assert d2o["q_max"] == pytest.approx(29.9, abs=1e-3)
    assert d2o["n_points"] == 595
    assert d2o["chi2"] > 0
    assert fit["H2O"]["trajectory_dataset"] == "SLS18498-H2O"


def test_shells(example_analysis: dict) -> None:
    shells = example_analysis["shells"]
    assert list(shells) == ["OW-OW", "OW-HW", "HW-HW"]
    oh = shells["OW-HW"]
    assert oh["centre"] == "OW" and oh["neighbour"] == "HW"
    assert oh["r_peak"] == pytest.approx(1.83, abs=0.1)
    assert oh["g_peak"] == pytest.approx(0.97, abs=0.1)
    assert oh["r_min"] == pytest.approx(2.4, abs=0.1)
    assert oh["n_coord"] == pytest.approx(1.76, abs=0.1)
    # N(OW around HW) = N(HW around OW) * c_OW / c_HW, c_OW / c_HW = 1/2
    assert oh["n_coord_reverse"] == pytest.approx(oh["n_coord"] * 0.5, rel=1e-6)
    assert oh["partner_density"] == pytest.approx(0.1 * 2 / 3, rel=1e-4)
    assert "HW around each OW" in oh["description"]
    assert "OW around each HW" in oh["description"]
    oo = shells["OW-OW"]
    assert oo["r_peak"] == pytest.approx(2.8, abs=0.1)
    assert oo["r_min"] == pytest.approx(3.45, abs=0.1)
    assert oo["n_coord"] == pytest.approx(5.1, abs=0.3)
    assert oo["n_coord_reverse"] is None
    assert shells["HW-HW"]["g_peak"] > 1.0


def test_convergence(example_analysis: dict) -> None:
    conv = example_analysis["convergence"]
    assert conv["last_run"] == 4 and conv["ereq"] == 3.0
    assert conv["ep_at_limit"] is True
    assert conv["ep_fraction_of_ereq"] == pytest.approx(1.0, abs=0.01)
    assert conv["ep_on_from_iteration"] == 101
    assert 101 <= conv["ep_limit_reached_at_iteration"] <= 200
    assert conv["ep_limit_reached_after"] == conv["ep_limit_reached_at_iteration"] - 100
    assert conv["last_run_iterations"] == 100 and conv["quarter_iterations"] == 25
    assert conv["last_run_first_iteration"] == 201
    trends = conv["trends"]
    assert set(trends) == {"energy", "pressure", "ep_amplitude", "r_factor"}
    assert trends["energy"]["start"] < 0 and abs(trends["energy"]["change_percent"]) < 5
    assert trends["ep_amplitude"]["end"] == pytest.approx(3.0, abs=0.01)
    assert trends["r_factor"]["end"] == pytest.approx(2.2e-4, rel=0.1)
    assert conv["r_factor_initial"] == pytest.approx(0.021, abs=1e-3)
    assert conv["r_factor_before_last_run"] == pytest.approx(2.10e-4, rel=1e-2)
    assert conv["r_factor_change_percent"] == pytest.approx(5.5, abs=1.0)
    assert isinstance(conv["r_factor_rising"], bool)
    assert conv["energy_drift_percent"] == trends["energy"]["change_percent"]
    verdict = conv["verdict"]
    assert "the EP amplitude reached the ereq limit of 3 kJ/mol after" in verdict
    assert (
        "the R-factor rose from 2.1e-04 to 2.2e-04 over the last run (+6 % between run ends, "
        "-0.1 % between its first- and last-quarter means)" in verdict
    )
    assert "from 2.1e-02 at the start" in verdict
    assert "the energy is -44.9 kJ/mol per molecule" in verdict
    assert "the pressure is 3.9 kbar" in verdict
    assert verdict.endswith(".")


def test_flags(example_analysis: dict) -> None:
    flags = example_analysis["flags"]
    assert any("at the ereq limit" in flag for flag in flags)
    assert not any("quick" in flag for flag in flags)
    assert not any("trajectory" in flag for flag in flags)
    assert not any("poor fit" in flag for flag in flags)
    assert not any("high pressure" in flag for flag in flags)
    assert not any("language model" in flag for flag in flags)


def test_markdown(example_analysis: dict) -> None:
    text = analysis.analysis_markdown(example_analysis)
    for heading in (
        "# Liquid water (1000 molecules)",
        "## System",
        "## Schedule",
        "## Runs",
        "## Fit quality per data set",
        "## Convergence",
        "## First coordination shells",
        "## Flags",
    ):
        assert heading in text, heading
    assert "| 4 | EP refinement | 100 | 3 | 2.22e-04 |" in text
    assert "| D2O | `SLS18502-D2O.mint01` | 5.81e-02 | 1.21e-03 | 2.13e-04 | ×" in text
    assert "| R² | χ² | Q range (Å⁻¹) |" in text
    assert "| OW-HW | 1.83 | 0.97 | 2.40 | 1.76 | 0.88 |" in text
    assert "| OW-OW |" in text and "| same |" in text
    assert "ereq limit of 3 kJ/mol" in text
    assert "| first quarter | last quarter | change |" in text
    assert "- the EP amplitude is at the ereq limit" in text
    assert "1000 × Water (3 atoms)" in text and "31.072 Å" in text


def test_save_and_load(example_analysis: dict, tmp_path: Path) -> None:
    path = analysis.save_analysis(example_analysis, tmp_path / "sub" / "analysis.json")
    assert path == tmp_path / "sub" / "analysis.json"
    text = path.read_text()
    assert "NaN" not in text and "Infinity" not in text
    assert json.loads(text) == example_analysis
    assert analysis.load_analysis(path) == example_analysis
    (tmp_path / "bad.json").write_text("[1, 2]")
    with pytest.raises(ValueError, match="analysis object"):
        analysis.load_analysis(tmp_path / "bad.json")
    (tmp_path / "broken.json").write_text("{")
    with pytest.raises(ValueError, match="not valid JSON"):
        analysis.load_analysis(tmp_path / "broken.json")


def test_plain() -> None:
    value = {
        "a": np.float32(1.5),
        "b": np.array([1, 2]),
        "c": float("nan"),
        "d": (Path("x"), np.int64(3), np.bool_(True)),
        "e": float("inf"),
        "f": None,
        "g": "s",
        "h": True,
        3: np.float64("nan"),
    }
    assert analysis.plain(value) == {
        "a": 1.5,
        "b": [1, 2],
        "c": None,
        "d": ["x", 3, True],
        "e": None,
        "f": None,
        "g": "s",
        "h": True,
        "3": None,
    }
    json.dumps(analysis.plain(value), allow_nan=False)


# --------------------------------------------------------------------------- degraded runs
def test_without_trajectory(tmp_path: Path) -> None:
    result = build_example_result(tmp_path / "run", with_trajectory=False)
    assert result.trajectory_path is None and not result.has_trajectory
    out = analysis.analyse(result)
    assert out["schedule"]["recorded"] is False
    assert out["schedule"]["recorded_iterations"] is None
    assert any("not recorded" in flag for flag in out["flags"])
    runs = out["runs"]
    assert len(runs) == 4 and "r2" not in runs[-1] and "last_iteration" not in runs[-1]
    assert runs[-1]["r_factor"] == pytest.approx(2.22e-4, rel=1e-2)
    assert "r2" not in out["fit"]["D2O"]
    assert out["fit"]["D2O"]["final_r_factor"] == pytest.approx(2.13e-4, rel=1e-2)
    conv = out["convergence"]
    assert conv["ep_at_limit"] is True and conv["trends"] == {}
    assert conv["r_factor_change_percent"] == pytest.approx(5.5, abs=1.0)
    assert "ereq limit of 3 kJ/mol" in conv["verdict"]
    assert "after" not in conv["verdict"].split(";")[0]
    assert "(+6 % between run ends)" in conv["verdict"]
    # the first shells come from EPSR's own .g01 / .z01 files
    shells = out["shells"]
    assert list(shells) == ["OW-OW", "OW-HW", "HW-HW"]
    assert shells["OW-HW"]["r_peak"] == pytest.approx(1.83, abs=0.1)
    assert shells["OW-HW"]["r_min"] == pytest.approx(2.4, abs=0.1)
    assert shells["OW-HW"]["n_coord"] == pytest.approx(1.76, abs=0.1)
    assert shells["OW-HW"]["n_coord_reverse"] == pytest.approx(0.88, abs=0.06)
    assert shells["OW-OW"]["n_coord"] == pytest.approx(5.1, abs=0.3)
    assert shells["OW-OW"]["n_coord_reverse"] is None
    assert shells["OW-HW"]["partner_density"] == pytest.approx(0.1 * 2 / 3, rel=1e-4)
    text = analysis.analysis_markdown(out)
    assert "trajectory not recorded" in text
    assert "| R² |" not in text
    assert "| OW-HW | 1.83 |" in text
    json.dumps(out, allow_nan=False)


def test_manifest_only(tmp_path: Path) -> None:
    result = build_example_result(
        tmp_path / "run", with_trajectory=False, with_epsr=False, with_data=False
    )
    out = analysis.analyse(result)
    assert out["shells"] == {}
    assert out["runs"][-1]["r_factor"] == pytest.approx(2.22e-4, rel=1e-2)
    assert out["system"]["atom_types"][1]["count"] == 2000  # from the specification
    text = analysis.analysis_markdown(out)
    assert "No first-shell numbers" in text
    json.dumps(out, allow_nan=False)


def test_runs_derived_from_the_trajectory_without_manifest_runs(
    example_result: RunResult,
) -> None:
    manifest = {
        key: value
        for key, value in example_result.manifest.items()
        if key not in ("runs", "final")
    }
    result = RunResult(
        folder=example_result.folder, spec=example_result.spec, manifest=manifest
    )
    out = analysis.analyse(result)
    assert [r["run"] for r in out["runs"]] == [1, 2, 3, 4]
    assert [r["niter"] for r in out["runs"]] == [1, 99, 100, 100]
    assert [r["stage"] for r in out["runs"]] == [
        "scratch",
        "reference potential",
        "EP refinement",
        "EP refinement",
    ]
    assert [r["ereq"] for r in out["runs"]] == [0.0, 0.0, 3.0, 3.0]
    assert out["runs"][-1]["r_factor"] == pytest.approx(2.22e-4, rel=1e-2)
    assert out["runs"][-1]["r_factor_per_dataset"]["HDO"] == pytest.approx(
        2.30e-4, rel=1e-2
    )
    assert out["fit"]["D2O"]["initial_r_factor"] == pytest.approx(5.81e-2, rel=1e-2)


def test_quick_and_model_flags(example_result: RunResult) -> None:
    manifest = dict(example_result.manifest)
    manifest["quick"] = True
    manifest["complete"] = False
    spec = example_result.spec
    spec.source = "model"
    result = RunResult(
        folder=example_result.folder, spec=spec, manifest=manifest, quick=True
    )
    try:
        out = analysis.analyse(result)
    finally:
        spec.source = "library:water"
    flags = out["flags"]
    assert any(flag.startswith("quick schedule (300 iterations") for flag in flags)
    assert any("did not complete: 4 of 4" in flag for flag in flags)
    assert any("proposed by the language model" in flag for flag in flags)
    assert out["schedule"]["quick"] is True
