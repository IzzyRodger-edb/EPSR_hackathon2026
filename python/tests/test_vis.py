"""``epsr.vis``: recording EPSR iteration by iteration and the analysis figures."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epsr

from epsr_cases import compare_dirs

pytest.importorskip("matplotlib")
pytest.importorskip("scipy")
import matplotlib  # noqa: E402

matplotlib.use("Agg")  # draw without a display
from epsr import vis  # noqa: E402

#: Short runs: two iterations of one Monte Carlo cycle each.
SHORT_RUN = {"niter": 2, "ntimes": 1}


def test_recorded_run_writes_the_same_files(template: Path, tmp_path: Path) -> None:
    """Recording only reads results: the folder ends up as after ``run()``."""
    plain = epsr.Project.from_template(template, tmp_path / "plain")
    plain_simulation = plain.simulation(SHORT_RUN)
    plain_simulation.run()

    recorded = epsr.Project.from_template(template, tmp_path / "recorded")
    recorder = vis.TrajectoryRecorder()
    summary = vis.run_recorded(recorded.simulation(SHORT_RUN), recorder)

    comparison = compare_dirs(plain.folder, recorded.folder)
    assert comparison.all_identical, comparison.report()
    assert summary == plain_simulation.summary

    trajectory = recorder.trajectory()
    assert trajectory.nsteps == 2
    assert list(trajectory.history["iter"]) == [1, 2]
    # the last snapshot is the state the run ended in
    np.testing.assert_array_equal(trajectory.gr_inter[-1], plain_simulation.gr_inter())
    np.testing.assert_array_equal(
        trajectory.ep[-1], plain_simulation.potential.empirical()
    )
    for dataset in range(plain_simulation.n_datasets):
        np.testing.assert_array_equal(
            trajectory.fq_sim[-1, dataset], plain_simulation.simulated_fq(dataset)
        )
    np.testing.assert_array_equal(trajectory.positions[-1], plain_simulation.positions)


def test_per_iteration_metrics_match_epsr(template: Path, tmp_path: Path) -> None:
    """The recomputed R-factors equal EPSR's own values on the last iteration.

    ``rfactor_inst`` uses the current-configuration F(Q) that the recorder builds
    from ``Simulation.current_sq()``; matching EPSR's raw R-factor checks that it
    tells the inter- and intramolecular halves apart correctly.
    """
    project = epsr.Project.from_template(template, tmp_path / "work")
    recorder = vis.TrajectoryRecorder()
    summary = vis.run_recorded(project.simulation(SHORT_RUN), recorder)
    fit = vis.metrics.fit_metrics(recorder.trajectory())
    np.testing.assert_allclose(
        fit["rfactor"][-1], summary.r_factor_per_dataset, rtol=1e-4
    )
    np.testing.assert_allclose(
        fit["rfactor_inst_tot"][-1], summary.r_factor_raw, rtol=1e-3
    )


def test_refine_recorded_matches_refine(template: Path, tmp_path: Path) -> None:
    arguments = {
        "n_runs": 2,
        "overrides": SHORT_RUN,
        "initial_overrides": {"ereq": 2.0},
    }
    plain = epsr.Project.from_template(template, tmp_path / "plain")
    plain_history = epsr.refine(plain, progress=False, **arguments)

    recorded = epsr.Project.from_template(template, tmp_path / "recorded")
    recorder = vis.TrajectoryRecorder()
    history = vis.refine_recorded(recorded, recorder=recorder, **arguments)

    comparison = compare_dirs(plain.folder, recorded.folder)
    assert comparison.all_identical, comparison.report()
    np.testing.assert_array_equal(history.r_factor, plain_history.r_factor)
    trajectory = recorder.trajectory()
    assert list(trajectory.history["run"]) == [1, 1, 2, 2]
    assert list(trajectory.history["ereq"]) == [2.0] * 4  # read from the potential


def test_report_and_round_trip(template: Path, tmp_path: Path) -> None:
    project = epsr.Project.from_template(template, tmp_path / "work")
    recorder = vis.TrajectoryRecorder()
    vis.refine_recorded(project, 2, recorder, overrides=SHORT_RUN)
    trajectory = vis.Trajectory.load(recorder.save(tmp_path / "trajectory"))
    np.testing.assert_array_equal(trajectory.sq_inter, recorder.trajectory().sq_inter)

    written = vis.make_report(trajectory, tmp_path / "figures", verbose=False)
    figures = [name for name in written if name[:2].isdigit()]
    assert len(figures) >= 8  # 08_local_structure needs a known system
    assert all(path.stat().st_size > 0 for path in written.values())
