"""``mellor.script``: the standalone script, rendered, run, and compared with the wrapper workflow."""

from __future__ import annotations

import ast
import json
import shutil
import subprocess
import sys

import pytest
from specs import REPOSITORY, cu_btc_spec, silica_spec, water_spec

from mellor import inputs, script
from mellor.spec import Schedule

WRAPPER_TESTS = REPOSITORY / "python" / "tests"


# --------------------------------------------------------------------------- rendering
def rendered(spec):
    return script.render_script(spec, inputs.generate(spec))


def test_script_name():
    assert script.script_name(water_spec()) == "run_waterbox.py"


def test_epsrvis_case():
    assert script.epsrvis_case(water_spec()) == "water"
    assert script.epsrvis_case(silica_spec()) == "silica_bulk"
    assert script.epsrvis_case(cu_btc_spec()) == "CuBTC"


def test_render_script_is_self_contained_python():
    spec = water_spec()
    generated = inputs.generate(spec)
    text = script.render_script(spec, generated)
    tree = ast.parse(text)
    imported = {
        (node.module if isinstance(node, ast.ImportFrom) else alias.name).split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    }
    assert "mellor" not in imported
    assert {"epsr", "epsrvis"} <= imported
    # the input texts are embedded verbatim and the configuration is complete
    for name, content in generated.files().items():
        assert content in text, name
    for name, content in generated.answer_files().items():
        assert content in text, name
    for constant in (
        "NAME = 'waterbox'",
        "KIND = 'molecular'",
        "MIXATO_ANSWERS = ['1', 'Water', '1000', '0.1', 'waterbox']",
        "FMOLE_ITERATIONS = 100",
        "SCHEDULE = [",
        "QUICK_SCHEDULE = [",
        "MELLOR_VERSION = ",
    ):
        assert constant in text, constant
    assert (
        "('EP refinement', 100, 3.0)" in text and "('EP refinement', 10, 3.0)" in text
    )
    assert "class RecorderAdapter" in text


def test_render_script_embeds_the_schedule():
    spec = water_spec()
    spec.schedule = Schedule(
        ereq=5.0,
        equilibration_iterations=20,
        refinement_runs=3,
        refinement_iterations=30,
        quick_equilibration_iterations=2,
        quick_refinement_iterations=3,
    )
    text = rendered(spec)
    namespace = _configuration_namespace(text)
    assert (
        namespace["SCHEDULE"]
        == [("scratch", 1, 0.0), ("reference potential", 20, 0.0)]
        + [("EP refinement", 30, 5.0)] * 3
    )
    assert (
        namespace["QUICK_SCHEDULE"]
        == [("scratch", 1, 0.0), ("reference potential", 2, 0.0)]
        + [("EP refinement", 3, 5.0)] * 3
    )
    assert namespace["DATASETS"] == [
        ("D2O", "SLS18502-D2O.mint01", "D2O.NWTStot.wts"),
        ("HDO", "SLS18500-HDO5050.mint01", "HDO.NWTStot.wts"),
        ("H2O", "SLS18498-H2O.mint01", "H2O.NWTStot.wts"),
    ]
    assert (
        namespace["INPUT_FILES"]["Water.mol"]
        == inputs.generate(spec).mol_files["Water.mol"]
    )
    assert namespace["PDB_FILE"] is None


def test_render_script_pdb_and_atomic():
    cu_btc = _configuration_namespace(rendered(cu_btc_spec()))
    assert cu_btc["KIND"] == "pdb" and cu_btc["PDB_FILE"] == "CuBTC.pdb"
    assert cu_btc["READPDB_ANSWERS"] == ["0.68", "n"]
    assert cu_btc["CHANGEATO_COMMANDS"][-3:] == ["save", "y", "e"]
    assert cu_btc["MOLECULES"] == [] and cu_btc["MIXATO_ANSWERS"] == []
    silica = _configuration_namespace(rendered(silica_spec()))
    assert silica["KIND"] == "atomic" and silica["MOLECULES"] == [
        ("Si", 1500),
        ("O", 3000),
    ]
    assert silica["EPSRVIS_CASE"] == "silica_bulk"


def test_render_script_rejects_foreign_inputs():
    with pytest.raises(ValueError, match="generated for"):
        script.render_script(water_spec(), inputs.generate(silica_spec()))


def test_write_script(tmp_path):
    spec = water_spec()
    path = script.write_script(spec, inputs.generate(spec), tmp_path / "run")
    assert path == tmp_path / "run" / "run_waterbox.py"
    assert path.read_text().startswith("#!/usr/bin/env python3\n")


def test_script_help_runs_without_mellor(tmp_path):
    """The script imports and parses its arguments outside the test process."""
    spec = water_spec()
    path = script.write_script(spec, inputs.generate(spec), tmp_path)
    completed = subprocess.run(
        [sys.executable, "-I", str(path), "--help"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert completed.returncode == 0, completed.stderr
    assert "--quick" in completed.stdout and "--no-record" in completed.stdout


def _configuration_namespace(text: str) -> dict:
    """Execute only the configuration block of a rendered script."""
    marker = "# The calculation. Nothing below needs editing."
    configuration = text.split(marker)[0]
    configuration = configuration.replace("import epsr\nimport epsrvis\n", "")
    namespace: dict = {}
    exec(compile(configuration, "configuration", "exec"), namespace)
    return namespace


# --------------------------------------------------------------------------- running it
@pytest.fixture(scope="session")
def water_quick_run(tmp_path_factory):
    """The generated water script run once with ``--quick`` (31 EPSR iterations)."""
    spec = water_spec()
    generated = inputs.generate(spec)
    folder = tmp_path_factory.mktemp("water_script")
    path = script.write_script(spec, generated, folder)
    (folder / "data").mkdir()
    for data_file in generated.data_files:
        shutil.copy2(data_file, folder / "data" / data_file.name)
    completed = subprocess.run(
        [sys.executable, path.name, "--out", ".", "--quick"],
        cwd=folder,
        capture_output=True,
        text=True,
    )
    return folder, completed


@pytest.mark.slow
def test_script_quick_run(water_quick_run):
    folder, completed = water_quick_run
    assert completed.returncode == 0, completed.stdout + completed.stderr
    epsr_folder = folder / "epsr"
    for name in (
        "waterbox.EPSR.g01",
        "waterbox.EPSR.u01",
        "waterbox.ato",
        "waterbox.pcof",
        "D2O.NWTStot.wts",
        "mixato.stdin",
    ):
        assert (epsr_folder / name).is_file(), name
    assert not (epsr_folder / "Water.mol").exists()

    manifest = json.loads((folder / "manifest.json").read_text())
    assert manifest["name"] == "waterbox" and manifest["kind"] == "molecular"
    assert manifest["quick"] is True and manifest["complete"] is True
    assert manifest["epsr_folder"] == "epsr" and manifest["trajectory"] == "trajectory"
    assert manifest["inp_name"] == "waterbox.EPSR.inp"
    assert manifest["n_atoms"] == 3000 and manifest["n_molecules"] == 1000
    assert manifest["box_side"] == pytest.approx(31.0723, abs=1e-3)
    assert manifest["atom_types"] == ["OW", "HW"]
    assert manifest["pairs"] == ["OW-OW", "OW-HW", "HW-HW"]
    assert [d["label"] for d in manifest["datasets"]] == ["D2O", "HDO", "H2O"]
    assert manifest["datasets"][0] == {
        "label": "D2O",
        "file": "SLS18502-D2O.mint01",
        "wts": "D2O.NWTStot.wts",
    }
    assert manifest["versions"]["mellor"] == script.MELLOR_VERSION
    assert set(manifest["versions"]) == {"epsr", "epsrvis", "mellor"}
    assert manifest["setup_steps"][0] == "makemole Water"

    runs = manifest["runs"]
    assert [run["stage"] for run in runs] == [
        "scratch",
        "reference potential",
        "EP refinement",
        "EP refinement",
    ]
    assert [run["niter"] for run in runs] == [1, 10, 10, 10]
    assert [run["ereq"] for run in runs] == [0.0, 0.0, 3.0, 3.0]
    assert [run["overrides"] for run in runs] == [[], ["niter 10"], ["ereq 3.0"], []]
    assert all(run["duration_s"] > 0 for run in runs)
    r_factors = [run["summary"]["r_factor"] for run in runs]
    assert r_factors[-1] < r_factors[0]
    assert len(runs[0]["summary"]["r_factor_per_dataset"]) == 3
    assert runs[-1]["summary"]["ep_amplitude"] > 0
    assert manifest["final"] == runs[-1]["summary"]

    import epsrvis

    trajectory = epsrvis.Trajectory.load(folder / "trajectory")
    assert trajectory.nsteps == 31
    assert trajectory.case == "water"
    assert trajectory.pairs == ["OW-OW", "OW-HW", "HW-HW"]
    assert trajectory.datasets == ["SLS18502-D2O", "SLS18500-HDO5050", "SLS18498-H2O"]
    assert list(trajectory.history["run"][:12]) == [1] + [2] * 10 + [3]
    assert trajectory.history["ereq"][-1] == 3.0
    assert trajectory.history["rfactot"][-1] == pytest.approx(r_factors[-1], rel=1e-6)
    # the Fortran-named fields resolve through the adapter (no NaN columns)
    for column in ("foldsv", "foldsb", "sentb", "epcup", "epccoeff", "chisq"):
        assert not any(map(lambda v: v != v, trajectory.history[column])), column
    assert trajectory.positions.shape == (31, 3000, 3)
    assert "run 4/4" in completed.stdout and "R-factor" in completed.stdout


@pytest.mark.slow
def test_script_reproduces_the_wrapper_workflow(water_quick_run, tmp_path):
    """The script's EPSR folder equals the one of the wrapper test workflow
    (``validation/water/steps.txt`` shortened to 10 + 10 + 10 iterations)."""
    folder, completed = water_quick_run
    assert completed.returncode == 0, completed.stdout + completed.stderr
    sys.path.insert(0, str(WRAPPER_TESTS))
    from epsr_cases import compare_dirs
    from epsr_workflow import (
        WATER,
        load_steps,
        prepare_case_folder,
        run_steps_with_wrapper,
        shorten_schedule,
    )

    steps = shorten_schedule(load_steps(WATER), 10, 10)
    reference = prepare_case_folder(WATER, tmp_path / "wrapper")
    for _result in run_steps_with_wrapper(steps, reference, WATER.case_folder):
        pass

    comparison = compare_dirs(reference, folder / "epsr")
    assert comparison.passed, comparison.report()
    assert comparison.all_identical, comparison.report()
    reverse = compare_dirs(folder / "epsr", reference)
    assert reverse.passed, reverse.report()
