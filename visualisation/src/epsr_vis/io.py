"""Readers for EPSR results.

Two sources are supported:

* a **trajectory** recorded by ``visualisation/recorder`` (one snapshot after
  every EPSR iteration: ``meta.json``, ``history.tsv`` and float32 arrays), and
* the ordinary **EPSR output files** of a run folder (``.EPSR.u01``, ``.g01``,
  ``.erg``, ...), which only hold the state after the last iteration.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import numpy as np

# Columns of history.tsv that are not plain floats.
_INT_COLS = ("step", "run", "iter", "intra_att", "intra_rej", "irot_att", "irot_rej",
             "rot_att", "rot_rej", "trans_att", "trans_rej")


@dataclass
class Trajectory:
    """Per-iteration record of an EPSR refinement (see recorder/src/main.rs).

    Array attributes are memory-mapped with a leading *step* axis, e.g.
    ``fq_sim[step, dataset, q]`` or ``gr_inter[step, pair, r]``.  Step ``k``
    (0-based) is EPSR iteration ``history["step"][k]`` (1-based, counted over all
    runs).
    """

    path: Path
    meta: dict = field(repr=False)
    history: dict = field(repr=False)

    @classmethod
    def load(cls, path: str | Path) -> "Trajectory":
        path = Path(path)
        meta = json.loads((path / "meta.json").read_text())
        history = _read_history(path / "history.tsv")
        return cls(path=path, meta=meta, history=history)

    # ---- sizes and labels
    @property
    def case(self) -> str:
        return self.meta["case"]

    @property
    def nsteps(self) -> int:
        return len(self.history["step"])

    @property
    def q(self) -> np.ndarray:
        return np.asarray(self.meta["q"], dtype=float)

    @property
    def r(self) -> np.ndarray:
        return np.asarray(self.meta["r"], dtype=float)

    @property
    def pairs(self) -> list[str]:
        return [f"{a}-{b}" for a, b in self.meta["pairs"]]

    @property
    def datasets(self) -> list[str]:
        """Data set names, from the data file names (e.g. ``SLS18502-D2O``)."""
        return [Path(d).stem for d in self.meta["datasets"]]

    @property
    def atom_types(self) -> list[str]:
        return self.meta["atom_types"]

    @property
    def atom_labels(self) -> np.ndarray:
        return np.asarray(self.meta["atom_labels"])

    @property
    def molecule_of_atom(self) -> np.ndarray:
        return np.asarray(self.meta["molecule_of_atom"], dtype=int)

    def q_slice(self, dataset: int) -> slice:
        """Q-grid indices that EPSR fits for ``dataset`` (``data_q_range``)."""
        lo, hi = self.meta["q_range"][dataset]
        return slice(lo, hi + 1)

    @property
    def steps(self) -> np.ndarray:
        return self.history["step"]

    @property
    def run_starts(self) -> list[tuple[int, int, float]]:
        """(first step index, run number, ereq) of each EPSR invocation."""
        out = []
        run = self.history["run"]
        for k in range(self.nsteps):
            if k == 0 or run[k] != run[k - 1]:
                out.append((k, int(run[k]), float(self.history["ereq"][k])))
        return out

    # ---- arrays (memory-mapped)
    def _array(self, name: str) -> np.ndarray:
        a, b = self.meta["arrays"][name]
        mm = np.memmap(self.path / f"{name}.f32", dtype="<f4", mode="r")
        n = min(mm.size // (a * b), self.nsteps)  # a trajectory still being recorded can be ahead of history.tsv
        return mm[: n * a * b].reshape(n, a, b)

    @cached_property
    def fq_data(self):
        return self._array("fq_data")

    @cached_property
    def fq_err(self):
        return self._array("fq_err")

    @cached_property
    def fq_sim(self):
        return self._array("fq_sim")

    @cached_property
    def fq_inst(self):
        return self._array("fq_inst")

    @cached_property
    def gr_inter(self):
        return self._array("gr_inter")

    @cached_property
    def gr_intra(self):
        return self._array("gr_intra")

    @cached_property
    def sq_inter(self):
        return self._array("sq_inter")

    @cached_property
    def sq_intra(self):
        return self._array("sq_intra")

    @cached_property
    def ep(self):
        return self._array("ep")

    @cached_property
    def ref_pot(self):
        return self._array("ref_pot")

    @cached_property
    def positions(self):
        return self._array("pos")

    @cached_property
    def initial_positions(self) -> np.ndarray:
        return np.fromfile(self.path / "pos_initial.f32", dtype="<f4").reshape(-1, 3)

    def box(self, step: int) -> float:
        """Cubic box side (Å) at step index ``step``; the box spans [-L/2, L/2)."""
        return float(self.history["box"][step])

    def gr_total(self, step: int) -> np.ndarray:
        """Partial g(r) including intramolecular contributions, [pair, r]."""
        return np.asarray(self.gr_inter[step]) + np.asarray(self.gr_intra[step])

    def concentrations(self) -> dict[str, float]:
        """Atom fraction of each atom type."""
        labels = self.atom_labels
        return {t: float(np.mean(labels == t)) for t in self.atom_types}


def _read_history(path: Path) -> dict[str, np.ndarray]:
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:] if ln.strip()]
    rows = [r for r in rows if len(r) == len(head)]  # last line may be incomplete while recording
    out: dict[str, np.ndarray] = {}
    for j, name in enumerate(head):
        col = [r[j] for r in rows]
        if name == "rfacd":
            out[name] = np.array([[float(x) for x in c.split(",")] for c in col])
        elif name in _INT_COLS:
            out[name] = np.array([int(x) for x in col])
        else:
            out[name] = np.array([float(x) for x in col])
    return out


# ---------------------------------------------------------------- EPSR output files

@dataclass
class EpsrTable:
    """One ``<run>.EPSR.?01`` file: x column, then (value, error) per column label."""

    x: np.ndarray
    labels: list[str]
    values: np.ndarray  # [column, x]
    errors: np.ndarray  # [column, x]

    def __getitem__(self, label: str) -> np.ndarray:
        return self.values[self.labels.index(label)]


def read_epsr_table(path: str | Path) -> EpsrTable:
    """Read an EPSR column file such as ``waterbox.EPSR.g01``.

    The first line is ``#  r  OW-OW  OW-HW ...``; data rows are x followed by
    (value, error) pairs.  ``.z01`` has one column per *ordered* pair and no errors.
    """
    path = Path(path)
    with open(path, encoding="latin-1") as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path, comments="#", ndmin=2)
    labels = header[1:]
    ncol = data.shape[1] - 1
    if ncol == 2 * len(labels):
        values, errors = data[:, 1::2].T, data[:, 2::2].T
    else:  # no error columns
        values, errors = data[:, 1:].T, np.zeros_like(data[:, 1:].T)
        labels = labels[:ncol]
    return EpsrTable(x=data[:, 0], labels=labels, values=values, errors=errors)


# Column meaning of the .EPSR.erg file (one line per iteration, written in epsr.rs).
ERG_COLUMNS = ("energy", "pressure", "rfactor_raw", "ep_amplitude", "number_density",
               "r_space_coeff", "quality_factor", "parabola_a", "parabola_b", "parabola_c")


def read_erg(path: str | Path) -> dict[str, np.ndarray]:
    """Read ``<run>.EPSR.erg``: the convergence history EPSR itself keeps.

    Note: ``rfactor_raw`` (and the r-space/quality columns) are only computed on
    the last iteration of each ``epsr`` invocation; on other lines they repeat the
    previous value.  Use a recorded :class:`Trajectory` for per-iteration R-factors.
    """
    a = np.loadtxt(path, ndmin=2)
    return {name: a[:, j] for j, name in enumerate(ERG_COLUMNS[: a.shape[1]])}


@dataclass
class OutputFolder:
    """The final EPSR output files of one run folder (``prefix`` e.g. ``waterbox.EPSR``)."""

    path: Path
    prefix: str

    @classmethod
    def find(cls, path: str | Path) -> "OutputFolder":
        path = Path(path)
        u = sorted(path.glob("*.EPSR.u01"))
        if not u:
            raise FileNotFoundError(f"no *.EPSR.u01 in {path}")
        return cls(path=path, prefix=u[0].name[: -len(".u01")])

    def table(self, ext: str) -> EpsrTable:
        """``ext`` is the two-letter code without the run number, e.g. ``"g"`` -> ``.g01``."""
        return read_epsr_table(self.path / f"{self.prefix}.{ext}01")

    def erg(self) -> dict[str, np.ndarray]:
        return read_erg(self.path / f"{self.prefix}.erg")


# What each EPSR output table holds (see validation/PYTHON_WRAPPER_WORKFLOW.md §4.1).
TABLES = {
    "t": "measured F(Q) of each data set (intramolecular part added back)",
    "u": "simulated F(Q) of each data set",
    "v": "difference data - simulation, per data set",
    "w": "Fourier transform of the data, G(r) per data set",
    "x": "Fourier transform of the simulation, G(r) per data set",
    "f": "simulated intermolecular partial S(Q)",
    "s": "simulated intramolecular partial S(Q)",
    "q": "estimated partial S(Q) from the data",
    "d": "difference estimated - simulated partial S(Q)",
    "g": "simulated intermolecular partial g(r)",
    "y": "simulated intramolecular partial g(r)",
    "r": "partial g(r) estimated from the data (Fourier transform of .q01)",
    "p": "empirical potential per pair (kJ/mol)",
    "o": "reference potential per pair (kJ/mol)",
    "l": "running integral 4*pi*rho*int r^2 (g_inter(r) - 1) dr per pair",
    "z": "running coordination number N_ab(r) from g_inter, per ordered pair (a-b, b-a)",
}
