"""``epsr.refine`` (the repeated-run loop) and ``epsr.History``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epsr
from epsr.io import read_erg

from conftest import CompletedRun

#: Value of ``nsumt`` (iterations already accumulated) in the template's input file.
TEMPLATE_ACCUMULATED_ITERATIONS = 102


def test_refine_three_runs(project: epsr.Project) -> None:
    erg_file = project.output("erg")
    erg_lines_before = len(read_erg(erg_file))
    callback_calls = []

    def record_call(run_number: int, simulation: epsr.Simulation) -> None:
        callback_calls.append((run_number, simulation.summary.iteration))

    history = epsr.refine(project, 3, callback=record_call)

    assert len(history) == 3
    assert callback_calls == [(1, 1), (2, 1), (3, 1)]
    assert len(read_erg(erg_file)) == erg_lines_before + 3
    np.testing.assert_array_equal(history.run, [1, 2, 3])
    assert history.r_factor.shape == (3,)
    assert history.r_factor_per_dataset.shape == (3, 3)
    assert all(isinstance(summary, epsr.Summary) for summary in history.summaries)
    accumulated_iterations = int(project.input_file().values["nsumt"])
    assert accumulated_iterations == TEMPLATE_ACCUMULATED_ITERATIONS + 3


def test_history_matches_erg(project: epsr.Project) -> None:
    history = epsr.refine(project, 2)
    erg_history = project.history()
    assert "parab_a" in erg_history.columns
    for column in ("energy", "ep_amplitude", "pressure"):
        np.testing.assert_allclose(
            erg_history.column(column)[-2:], history.column(column), rtol=1e-5
        )


def test_refine_stops_early(project: epsr.Project) -> None:
    history = epsr.refine(
        project, 10, callback=lambda run_number, simulation: run_number == 2
    )
    assert len(history) == 2
    history = epsr.refine(
        project, 10, until=lambda current_history: len(current_history) >= 1
    )
    assert len(history) == 1


def test_refine_overrides_and_snapshots(project: epsr.Project) -> None:
    history = epsr.refine(
        project,
        2,
        overrides={"ntimes": 2},
        initial_overrides={"ereq": 12.0},
        snapshot_every=1,
    )
    assert len(history) == 2
    assert project.simulation().control.ntimes == 2
    assert project.simulation().potential.params.ereq == pytest.approx(12.0)
    snapshots_folder = project.folder.parent / f"{project.folder.name}.snapshots"
    assert sorted(path.name for path in snapshots_folder.iterdir()) == [
        "run0001",
        "run0002",
    ]


def test_refine_accepts_path(project: epsr.Project) -> None:
    assert len(epsr.refine(project.folder, 1)) == 1


def test_history_export(project: epsr.Project, tmp_path: Path) -> None:
    history = epsr.refine(project, 2)
    csv_file = history.to_csv(tmp_path / "history.csv")
    header = csv_file.read_text().splitlines()[0].split(",")
    assert header[:3] == ["run", "iteration", "energy"]
    assert "r_factor_per_dataset_2" in header
    csv_data = np.loadtxt(csv_file, delimiter=",", skiprows=1)
    assert csv_data.shape == (2, len(header))

    assert history[0]["run"] == 1
    assert len(list(history)) == 2
    with pytest.raises(AttributeError):
        _ = history.no_such_column

    pandas = pytest.importorskip("pandas")
    dataframe = history.to_dataframe()
    assert isinstance(dataframe, pandas.DataFrame) and len(dataframe) == 2


def test_history_from_summaries(completed_run: CompletedRun) -> None:
    summary = completed_run.simulation.summary
    history = epsr.History.from_summaries([summary, summary])
    assert len(history) == 2
    assert history.r_factor[0] == pytest.approx(summary.r_factor)
