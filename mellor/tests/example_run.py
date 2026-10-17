"""A run folder built by hand on the toolkit's example data (water, 300 iterations).

Used by ``test_analysis.py`` and ``test_figures.py``: the folder looks like what
``mellor.runner.run_spec`` leaves behind (``data/``, ``spec.json``, ``epsr/``,
``trajectory/``, ``manifest.json``, ``run.log``), with the trajectory and the EPSR
output files of ``visualisation/python_visualisation_toolkit/example_H2O/data`` and
a manifest whose run summaries are taken from the trajectory's history at the end
of each of the four runs (1 + 99 + 100 + 100 iterations, ereq 0, 0, 3, 3).
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path
from typing import Any

import numpy as np

from mellor import library
from mellor.data import read_mint01
from mellor.runner import RunResult
from mellor.spec import SystemSpec

REPOSITORY = Path(__file__).resolve().parents[2]
EXAMPLE = (
    REPOSITORY
    / "visualisation"
    / "python_visualisation_toolkit"
    / "example_H2O"
    / "data"
)
WATER_DATA = REPOSITORY / "experimental_data" / "water" / "data"
#: In the order of the example's ``waterbox.EPSR.inp`` (D2O, HDO, H2O).
WATER_FILES = [
    WATER_DATA / name
    for name in (
        "SLS18502-D2O.mint01",
        "SLS18500-HDO5050.mint01",
        "SLS18498-H2O.mint01",
    )
]

#: The example's schedule (``example_H2O/README.md``): ``(stage, niter, ereq)``.
RUNS = [
    ("scratch", 1, 0.0),
    ("reference potential", 99, 0.0),
    ("EP refinement", 100, 3.0),
    ("EP refinement", 100, 3.0),
]


def example_available() -> bool:
    return (
        (EXAMPLE / "trajectory" / "meta.json").is_file()
        and (EXAMPLE / "epsr_output" / "waterbox.EPSR.u01").is_file()
        and all(path.is_file() for path in WATER_FILES)
    )


def water_spec() -> SystemSpec:
    """The library's water specification for the three SANDALS files."""
    spec = library.make_spec("water", [read_mint01(path) for path in WATER_FILES])
    spec.validate()
    return spec


def example_manifest(spec: SystemSpec, trajectory: Path) -> dict[str, Any]:
    """A manifest (docs/architecture.md §6) whose summaries come from the trajectory."""
    from epsrvis import Trajectory

    traj = Trajectory.load(trajectory)
    history = traj.history
    runs = []
    end = -1
    for number, (stage, niter, ereq) in enumerate(RUNS, start=1):
        end += niter
        summary = {
            "iteration": int(history["iter"][end]),
            "energy_inter": float(history["foldsv"][end]),
            "energy_intra": float(history["foldsb"][end]),
            "energy": float(history["energd"][end]),
            "pressure": float(history["pressd"][end]),
            "ep_amplitude": float(history["energabs"][end]),
            "r_factor": float(history["rfactot"][end]),
            "r_factor_raw": float(history["rfactot"][end]),
            "r_factor_per_dataset": [
                float(value) for value in np.atleast_1d(history["rfacd"][end])
            ],
            "ep_up": float(history["epcup"][end]),
            "ep_coefficient": float(history["epccoeff"][end]),
            "gr_roughness": float(history["sentb"][end]),
            "chisq": float(history["chisq"][end]),
        }
        runs.append(
            {
                "run": number,
                "stage": stage,
                "niter": niter,
                "ereq": ereq,
                "overrides": [],
                "duration_s": 1.5 * niter,
                "summary": summary,
            }
        )
    assert end == traj.nsteps - 1, "the schedule does not match the trajectory"
    return {
        "name": spec.name,
        "title": spec.title,
        "kind": spec.kind,
        "quick": False,
        "record": True,
        "created": "2026-10-17T12:00:00",
        "complete": True,
        "duration_s": sum(run["duration_s"] for run in runs),
        "python": "3.12.3",
        "versions": {"epsr": "0.1.0", "epsrvis": "0.1.0", "mellor": "0.1.0"},
        "script": f"run_{spec.name}.py",
        "epsr_folder": "epsr",
        "trajectory": "trajectory",
        "inp_name": spec.inp_name,
        "epsrvis_case": "water",
        "temperature": spec.temperature,
        "box_side": float(history["box"][-1]),
        "n_atoms": int(traj.meta["natom"]),
        "n_molecules": int(traj.meta["nmol"]),
        "density": float(history["rho"][-1]),
        "atom_types": list(traj.atom_types),
        "pairs": list(traj.pairs),
        "datasets": [
            {"label": d.label, "file": d.file_name, "wts": d.weights_file}
            for d in spec.datasets
        ],
        "setup_steps": [
            "makemole Water",
            "mixato < mixato.stdin",
            "randomise waterbox",
            "delete Water.mol",
            "fmole waterbox 100",
            "nwts D2O.NWTS.dat",
            "nwts HDO.NWTS.dat",
            "nwts H2O.NWTS.dat",
        ],
        "schedule": [{"stage": s, "niter": n, "ereq": e} for s, n, e in RUNS],
        "runs": runs,
        "final": dict(runs[-1]["summary"]),
    }


def _link_or_copy(source: Path, target: Path) -> None:
    try:
        os.symlink(source, target, target_is_directory=True)
    except OSError:  # no symlinks (Windows without privileges): copy instead
        shutil.copytree(source, target)


def build_example_result(
    folder: Path,
    *,
    with_trajectory: bool = True,
    with_epsr: bool = True,
    with_data: bool = True,
) -> RunResult:
    """Build the run folder in ``folder`` and return its :class:`RunResult`."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    spec = water_spec()
    if with_trajectory:
        _link_or_copy(EXAMPLE / "trajectory", folder / "trajectory")
    if with_epsr:
        _link_or_copy(EXAMPLE / "epsr_output", folder / "epsr")
    if with_data:
        (folder / "data").mkdir(exist_ok=True)
        for path in WATER_FILES:
            shutil.copy2(path, folder / "data" / path.name)
    manifest = example_manifest(spec, EXAMPLE / "trajectory")
    if not with_trajectory:
        manifest["trajectory"] = None
        manifest["record"] = False
    spec.save(folder / "spec.json")
    (folder / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (folder / "run.log").write_text("run folder built by tests/example_run.py\n")
    return RunResult(
        folder=folder,
        spec=spec,
        manifest=manifest,
        duration_s=float(manifest["duration_s"]),
        quick=False,
    )
