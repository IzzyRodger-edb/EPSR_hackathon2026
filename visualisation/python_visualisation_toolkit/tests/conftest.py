"""Shared fixtures: the example trajectories, and a stand-in for the wrapper ``Simulation``."""

from pathlib import Path

import numpy as np
import pytest

from epsrvis import Trajectory

ROOT = Path(__file__).resolve().parents[1]
H2O = ROOT / "example_H2O" / "data"
SIO2 = ROOT / "example_SiO2" / "data"


@pytest.fixture(scope="session")
def water():
    return Trajectory.load(H2O / "trajectory")


@pytest.fixture(scope="session")
def silica():
    return Trajectory.load(SIO2 / "trajectory")


class ReplaySimulation:
    """Mimics the planned PyO3 ``epsr.Simulation`` (rust/docs/PYTHON_WRAPPER_NOTES.md section 5).

    ``iterate(i)`` moves to the next recorded step of a trajectory, and the
    accessors return that step's arrays, so the recorder can be tested end to
    end without the wrapper.  Uses the *suggested Python* names (properties
    ``q``/``r``/``positions``, all-pairs ``gr_inter()``, summary dict with
    friendly keys) to exercise the recorder's name mapping.
    """

    def __init__(self, traj: Trajectory, start: int, niter: int):
        self.t, self.k, self.niter = traj, start - 2, niter  # iterate() advances to step start-1 (0-based)

    # grids and labels
    q = property(lambda self: self.t.q.astype(np.float32))
    r = property(lambda self: self.t.r.astype(np.float32))
    pair_labels = property(lambda self: [tuple(p) for p in self.t.meta["pairs"]])

    def atom_types(self):
        return self.t.atom_types

    def n_datasets(self):
        return len(self.t.datasets)

    def data_q_range(self, d):
        return tuple(self.t.meta["q_range"][d])

    def get_parameter(self, name, dataset=None):
        return self.t.meta["datasets"][dataset] if name == "datafile" else None

    def atom_labels(self):
        return list(self.t.atom_labels)

    def molecule_atoms(self):
        mol = self.t.molecule_of_atom
        edges = np.flatnonzero(np.diff(np.r_[-1, mol, -2]))
        return [range(a, b) for a, b in zip(edges[:-1], edges[1:])]

    # the simulation
    def iterate(self, i):
        self.k += 1

    def write_fort35(self):
        self.finished = True

    def box_size(self):
        return self.t.box(self.k)

    @property
    def positions(self):
        return np.asarray(self.t.positions[self.k])

    def data_fq(self, d):
        return self.t.fq_data[self.k, d]

    def data_errors(self, d):
        return self.t.fq_err[self.k, d]

    def simulated_fq(self, d):
        return self.t.fq_sim[self.k, d]

    def simulated_fq_current(self, d):
        return self.t.fq_inst[self.k, d]

    def gr_inter(self):
        return np.asarray(self.t.gr_inter[self.k])

    def gr_intra(self):
        return np.asarray(self.t.gr_intra[self.k])

    def sq_inter(self, ip):          # per-pair form, as the Rust accessor
        return self.t.sq_inter[self.k, ip]

    def sq_intra(self, ip):
        return self.t.sq_intra[self.k, ip]

    def empirical_potential(self):
        return np.asarray(self.t.ep[self.k])

    def reference_potential(self):
        return np.asarray(self.t.ref_pot[self.k])

    def mc_counters(self):
        h = self.t.history
        return [[h[f"{m}_att"][self.k] for m in ("intra", "irot", "rot", "trans")],
                [h[f"{m}_rej"][self.k] for m in ("intra", "irot", "rot", "trans")]]

    @property
    def summary(self):
        h, k = self.t.history, self.k
        return {"energy": h["energd"][k], "pressure": h["pressd"][k], "ep_amplitude": h["energabs"][k],
                "r_factor": h["rfactot"][k], "r_factor_per_dataset": list(h["rfacd"][k]),
                "foldsv": h["foldsv"][k], "foldsb": h["foldsb"][k], "chisq": h["chisq"][k],
                "sentb": h["sentb"][k], "epcup": h["epcup"][k], "epccoeff": h["epccoeff"][k]}


@pytest.fixture
def replay():
    return ReplaySimulation
