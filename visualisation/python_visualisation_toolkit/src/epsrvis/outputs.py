"""Readers for the ordinary EPSR output files of a simulation folder.

These hold the state after the *last* iteration of the last ``epsr`` call
only.  They need no recording and work on any EPSR (Fortran or Rust) folder.
The Python wrapper's planned ``epsr.read_table`` helper
(rust/docs/PYTHON_WRAPPER_NOTES.md section 5) is :func:`read_epsr_table`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


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
