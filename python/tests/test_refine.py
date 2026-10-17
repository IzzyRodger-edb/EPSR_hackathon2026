"""``epsr.refine`` (the repeated-run loop) and ``epsr.History``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epsr
from epsr.io import read_erg


def test_refine_three_runs(project: epsr.Project) -> None:
    erg = project.output("erg")
    n0 = len(read_erg(erg))
    seen = []
    hist = epsr.refine(project, 3, callback=lambda k, sim: seen.append((k, sim.summary.iteration)))
    assert len(hist) == 3
    assert seen == [(1, 1), (2, 1), (3, 1)]
    assert len(read_erg(erg)) == n0 + 3
    np.testing.assert_array_equal(hist.run, [1, 2, 3])
    assert hist.r_factor.shape == (3,)
    assert hist.r_factor_per_dataset.shape == (3, 3)
    assert all(isinstance(s, epsr.Summary) for s in hist.summaries)
    # nsumt (number of accumulated iterations) advances with every run
    assert int(project.input_file().values["nsumt"]) == 102 + 3


def test_history_matches_erg(project: epsr.Project) -> None:
    hist = epsr.refine(project, 2)
    erg = project.history()
    assert "parab_a" in erg.columns
    np.testing.assert_allclose(erg.energy[-2:], hist.energy, rtol=1e-5)
    np.testing.assert_allclose(erg.ep_amplitude[-2:], hist.ep_amplitude, rtol=1e-5)
    np.testing.assert_allclose(erg.pressure[-2:], hist.pressure, rtol=1e-5)


def test_refine_stops_early(project: epsr.Project) -> None:
    hist = epsr.refine(project, 10, callback=lambda k, sim: k == 2)
    assert len(hist) == 2
    hist = epsr.refine(project, 10, until=lambda h: len(h) >= 1)
    assert len(hist) == 1


def test_refine_overrides_and_snapshots(project: epsr.Project) -> None:
    hist = epsr.refine(
        project,
        2,
        overrides={"ntimes": 2},
        initial_overrides={"ereq": 12.0},
        snapshot_every=1,
    )
    assert len(hist) == 2
    assert project.simulation().control.ntimes == 2
    assert project.simulation().potential.params.ereq == pytest.approx(12.0)
    snaps = project.folder.parent / f"{project.folder.name}.snapshots"
    assert sorted(p.name for p in snaps.iterdir()) == ["run0001", "run0002"]


def test_refine_accepts_path(project: epsr.Project) -> None:
    assert len(epsr.refine(project.folder, 1)) == 1


def test_history_export(project: epsr.Project, tmp_path: Path) -> None:
    hist = epsr.refine(project, 2)
    csv = hist.to_csv(tmp_path / "h.csv")
    header = csv.read_text().splitlines()[0].split(",")
    assert header[:3] == ["run", "iteration", "energy"]
    assert "r_factor_per_dataset_2" in header
    data = np.loadtxt(csv, delimiter=",", skiprows=1)
    assert data.shape == (2, len(header))
    assert hist[0]["run"] == 1 and len(list(hist)) == 2
    with pytest.raises(AttributeError):
        _ = hist.no_such_column
    pd = pytest.importorskip("pandas")
    df = hist.to_dataframe()
    assert isinstance(df, pd.DataFrame) and len(df) == 2


def test_history_from_summaries(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    _, sim = ran
    h = epsr.History.from_summaries([sim.summary, sim.summary])
    assert len(h) == 2
    assert h.r_factor[0] == pytest.approx(sim.summary.r_factor)
