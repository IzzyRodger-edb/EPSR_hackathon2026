"""The matplotlib helpers in ``epsr.plot`` draw without errors."""

from __future__ import annotations

from pathlib import Path

import pytest

import epsr

from conftest import CompletedRun

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")  # draw without a display
pyplot = pytest.importorskip("matplotlib.pyplot")


def test_plots(completed_run: CompletedRun, tmp_path: Path) -> None:
    simulation = completed_run.simulation

    fits_axes = epsr.plot.fits(simulation)
    assert len(fits_axes.lines) == 2 * simulation.n_datasets  # data and fit per dataset

    partials_axes = epsr.plot.partials(simulation, [("O1", "H1"), ("O1", "O1")])
    assert len(partials_axes.lines) == 2
    all_pairs_axes = epsr.plot.partials(simulation, kind="sq_intra")
    assert len(all_pairs_axes.lines) == simulation.n_pairs

    for which in ("empirical", "reference", "total"):
        epsr.plot.potentials(simulation, [("O1", "H1")], which=which)

    history = epsr.History.from_summaries([simulation.summary, simulation.summary])
    history_axes = epsr.plot.history(history, ["r_factor", "energy", "not_a_column"])
    assert len(history_axes) == 2  # unknown columns are skipped

    image_file = tmp_path / "history.png"
    history_axes[0].figure.savefig(image_file)
    assert image_file.stat().st_size > 0
    pyplot.close("all")
