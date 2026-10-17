"""Tests of the toolkit on the two example data sets."""

import numpy as np
import matplotlib

matplotlib.use("Agg")

from epsrvis import OutputFolder, Trajectory, TrajectoryRecorder, metrics, plots, record_run
from epsrvis.report import make_report

from conftest import H2O, SIO2


def test_trajectories_load(water, silica):
    assert water.nsteps == 300 and silica.nsteps == 300
    assert water.pairs == ["OW-OW", "OW-HW", "HW-HW"]
    assert silica.pairs == ["Si-Si", "Si-O", "O-O"]
    assert water.fq_sim.shape == (300, 3, 600)
    assert silica.positions.shape == (300, 4500, 3)
    assert [r for _, r, _ in water.run_starts] == [1, 2, 3, 4]


def test_rfactor_matches_epsr(water, silica):
    """The per-step R-factor equals EPSR's rfacd wherever EPSR computes it (end of each run)."""
    for t in (water, silica):
        fit = metrics.fit_metrics(t)
        for s in plots.key_steps(t):
            np.testing.assert_allclose(fit["rfactor"][s], t.history["rfacd"][s], rtol=2e-5)


def test_final_outputs_match_trajectory(water):
    """The last snapshot equals the EPSR output files written on that iteration."""
    out = OutputFolder.find(H2O / "epsr_output")
    g = out.table("g")
    np.testing.assert_allclose(g.values, water.gr_inter[-1][:, : g.values.shape[1]], atol=1e-6)
    p = out.table("p")
    np.testing.assert_allclose(p.values, water.ep[-1][:, : p.values.shape[1]], atol=1e-5)
    erg = out.erg()
    np.testing.assert_allclose(erg["energy"][-1], water.history["energd"][-1], rtol=1e-5)


def test_recorder_replays_wrapper_api(water, replay, tmp_path):
    """record_run on a wrapper-like Simulation reproduces the recorded trajectory."""
    rec = TrajectoryRecorder(case="water")
    pcof = H2O / "epsr_output" / "waterbox.pcof"
    for start, niter in ((1, 1), (2, 99), (101, 100)):   # first three runs
        sim = replay(water, start, niter)
        record_run(sim, rec, ereq=water.history["ereq"][start - 1])
        assert sim.finished
    t = rec.trajectory()
    assert t.nsteps == 200
    for name in ("fq_data", "fq_sim", "fq_inst", "gr_inter", "sq_inter", "ep", "pos"):
        np.testing.assert_array_equal(getattr(t, name), getattr(water, name)[:200], err_msg=name)
    for key in ("run", "iter", "ereq", "energd", "pressd", "energabs", "rot_rej", "rho"):
        np.testing.assert_allclose(t.history[key], water.history[key][:200], rtol=1e-6, err_msg=key)
    # save -> load round trip, then the analysis runs on it
    back = Trajectory.load(rec.save(tmp_path / "traj"))
    np.testing.assert_array_equal(back.gr_inter, t.gr_inter)
    assert back.pairs == water.pairs and back.datasets == water.datasets
    fit = metrics.fit_metrics(back)
    np.testing.assert_allclose(fit["rfactor"][99], water.history["rfacd"][99], rtol=2e-5)
    # ereq read from a .pcof
    rec2 = TrajectoryRecorder()
    rec2.start_run(replay(water, 1, 1), pcof=pcof)
    assert rec2._ereq == 3.0


def test_recorder_without_optional_quantities(water, replay):
    """A minimal wrapper (no errors, no current-configuration F(Q), no MC counters) still works."""

    class Minimal(replay):
        def __getattribute__(self, name):
            if name in ("data_errors", "simulated_fq_current", "mc_counters"):
                raise AttributeError(name)
            return super().__getattribute__(name)
    rec = TrajectoryRecorder(case="water", record_positions=False)
    record_run(Minimal(water, 2, 20), rec)
    t = rec.trajectory()
    fit = metrics.fit_metrics(t)
    assert np.all(np.isnan(fit["rfactor_inst"])) and np.all(np.isfinite(fit["rfactor"]))
    assert metrics.acceptance(t) == {}
    fig = plots.convergence(t)  # must not fail on the missing columns
    assert fig is not None


def test_structure_metrics(water, silica):
    s = silica.nsteps - 1
    n, frac = metrics.coordination_distribution(silica.positions[s], silica.box(s), silica.atom_labels, "Si", "O", 2.0)
    assert 0.75 < frac[4] < 0.9  # mostly 4-fold Si
    q = metrics.tetrahedral_order(silica.positions[s], silica.box(s), silica.atom_labels, "Si", "O")
    assert 0.8 < q.mean() < 1.0
    th, hist, _ = metrics.bond_angle_distribution(water.positions[-1], water.box(-1), water.atom_labels,
                                                  "OW", "HW", 1.25, intramolecular_only=True,
                                                  molecule=water.molecule_of_atom)
    assert 100 < th[np.argmax(hist)] < 115  # H-O-H near 107 degrees
    sh = metrics.shell_history(water, every=50)
    assert 4.5 < sh["OW-OW"]["n_coord"][-1] < 6


def test_report(silica, tmp_path):
    written = make_report(silica.subset(range(0, 300, 10)), tmp_path, verbose=False)
    assert len([k for k in written if k.endswith(".png") or k[:2].isdigit()]) == 9
    assert (tmp_path / "summary.md").read_text().startswith("| quantity |")
