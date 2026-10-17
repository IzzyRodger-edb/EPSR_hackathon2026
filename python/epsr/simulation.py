"""One EPSR run: :class:`Simulation`, its :class:`Summary` and the input dataclasses.

A :class:`Simulation` is one invocation of the Fortran ``epsr`` program: it reads
every input file in a folder (``Epsr::setup``), runs ``niter`` iterations
(Monte Carlo cycles followed by an empirical-potential update) and, at the end
of the last iteration, **rewrites the input files in place** (``.ato``,
``.pcof``, ``.EPSR.inp``, ``.EPSR.inpa``), appends to ``.erg`` and writes every
result file. That is how EPSR carries its state from one run to the next. Work
on a copy of your folder (:meth:`epsr.Project.from_template`).
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import asdict, dataclass, field, fields
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Any, Union

import numpy as np
import numpy.typing as npt

from . import _core

if TYPE_CHECKING:
    from .potential import Potential

__all__ = [
    "Simulation",
    "Summary",
    "Control",
    "Dataset",
    "Fit",
    "Overrides",
    "PairKey",
    "render_overrides",
    "find_inp",
]

F32 = npt.NDArray[np.float32]
#: A pair is selected by its index in ``pair_labels`` or by its two atom-type labels.
PairKey = Union[int, tuple[str, str]]
#: Overrides: raw ``.EPSR.inp`` lines, or a mapping rendered by :func:`render_overrides`.
Overrides = Union[Mapping[str, Any], Sequence[str], None]


@dataclass(frozen=True)
class Summary:
    """Results of the most recent iteration (what EPSR writes to ``.out``/``.erg``).

    Each field gives the Fortran variable name in its comment; the raw dict with
    Fortran names is in :attr:`fortran`.
    """

    iteration: int  # iter
    energy_inter: float  # foldsv: running intermolecular energy per molecule (kJ/mol)
    energy_intra: float  # foldsb: running intramolecular energy per molecule (kJ/mol)
    energy: float  # energd: intermolecular energy per molecule from `press` (kJ/mol)
    pressure: float  # pressd (kbar)
    ep_amplitude: float  # energabs: absolute energy of the empirical potential
    r_factor: float  # rfactot: overall R-factor
    r_factor_raw: float  # rfacraw
    r_factor_per_dataset: tuple[float, ...]  # rfacd
    ep_up: float  # epcup
    ep_coefficient: float  # epccoeff: r-space quality coefficient
    gr_roughness: float  # sentb: roughness of the averaged g(r)
    chisq: float  # chisq
    fortran: Mapping[str, Any] = field(repr=False, compare=False, default_factory=dict)

    @classmethod
    def from_fortran(cls, d: Mapping[str, Any]) -> Summary:
        return cls(
            iteration=int(d["iter"]),
            energy_inter=float(d["foldsv"]),
            energy_intra=float(d["foldsb"]),
            energy=float(d["energd"]),
            pressure=float(d["pressd"]),
            ep_amplitude=float(d["energabs"]),
            r_factor=float(d["rfactot"]),
            r_factor_raw=float(d["rfacraw"]),
            r_factor_per_dataset=tuple(float(x) for x in d["rfacd"]),
            ep_up=float(d["epcup"]),
            ep_coefficient=float(d["epccoeff"]),
            gr_roughness=float(d["sentb"]),
            chisq=float(d["chisq"]),
            fortran=dict(d),
        )

    def as_dict(self) -> dict[str, Any]:
        """The descriptive fields as a plain dict (without :attr:`fortran`)."""
        d = asdict(self)
        d.pop("fortran")
        return d


@dataclass(frozen=True)
class Dataset:
    """One ``data n`` block of the ``.EPSR.inp`` file (Fortran names)."""

    infili: str  # datafile as given in the input
    fnamedata: str  # data file path as resolved by the engine
    fnamewts: str  # weights file path
    nrtype: int  # data format
    rshmin: float  # minimum r for background subtraction (Å)
    szeros: float  # Q=0 limit (0 = use first data point)
    tweak: float  # scale factor on the data
    efilereq: float  # requested EP amplitude for this dataset
    porodadd: float
    nqmind: int
    nqmaxd: int
    porodamp: float

    @property
    def name(self) -> str:
        """The data file name without its extension (the label used in result files)."""
        return Path(self.fnamedata).stem

    @classmethod
    def from_fortran(cls, d: Mapping[str, Any]) -> Dataset:
        return cls(**{f.name: d[f.name] for f in fields(cls)})


@dataclass(frozen=True)
class Control:
    """Run-control values read from the input (Fortran names, ``Control`` in ``epsr.rs``)."""

    nq: int
    qstep: float
    ireset: int
    iinit: int
    ntimes: int
    niter: int
    intra: int
    rotfreq: int
    inter: int
    nisumm: int
    ngrsamples: int
    fwhm: float
    fwhmq: float
    nsmoop: int
    pressreq: float
    ereqmin: float
    ereqmax: float
    ereqstepread: float
    epcthresh: float
    ereqbias: float
    navmaxparab: int
    nparabstartfactor: int
    epcrmin: float
    epcfrac: float
    num_threads: int
    nmolcellin: int
    thresholde: float
    stepmisave: float
    steprisave: float
    rmaxgrsave: float
    irefsizefactor: int
    porodq: float
    porodpwr: float
    nporodq: int
    pressure: Mapping[str, Any]

    @classmethod
    def from_fortran(cls, d: Mapping[str, Any]) -> Control:
        return cls(**{f.name: d[f.name] for f in fields(cls)})


@dataclass(frozen=True)
class Fit:
    """Data and simulated F(Q) of one dataset over its fitted Q range."""

    name: str
    q: F32
    data: F32
    errors: F32
    simulated: F32
    q_range: tuple[int, int]  # 0-based inclusive indices into the full Q grid

    @property
    def residual(self) -> F32:
        return self.data - self.simulated


def _fmt_value(v: Any) -> str:
    if isinstance(v, bool):
        return "T" if v else "F"
    if isinstance(v, (list, tuple, np.ndarray)):
        return " ".join(_fmt_value(x) for x in v)
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    return str(v)


def render_overrides(overrides: Overrides) -> list[str]:
    """Turn overrides into ``.EPSR.inp`` lines for ``Epsr::setup_with_lines``.

    * A sequence of strings is passed through unchanged
      (``["ntimes 10", "data 2", "tweak 1.05"]``).
    * A mapping renders each item as ``"key value"``. Lists and tuples become
      space-separated values (``{"potfac": [1.0, 1.0]}`` gives
      ``"potfac 1.0 1.0"``), booleans become ``T``/``F``.
    * The key ``"data"`` takes a mapping from **EPSR dataset numbers (1-based,
      as in the input file)** to per-dataset keywords:
      ``{"data": {2: {"tweak": 1.05}}}`` gives ``["data 2", "tweak 1.05"]``.
      Dataset blocks are emitted after all global keywords.

    Keywords are those of the ``.EPSR.inp`` file and the ``.pcof`` header
    (``ereq``, ``rminpt``, ...). Overrides are written back to those files at
    the end of the run, so they persist into the next run.
    """
    if overrides is None:
        return []
    if isinstance(overrides, str):
        return [overrides]
    if not isinstance(overrides, Mapping):
        lines = list(overrides)
        if not all(isinstance(x, str) for x in lines):
            raise TypeError("override lines must be strings")
        return lines
    lines = []
    data = None
    for key, value in overrides.items():
        if key == "data":
            data = value
            continue
        if isinstance(value, Mapping):
            raise TypeError(f"override {key!r}: nested mappings are only allowed under 'data'")
        lines.append(f"{key} {_fmt_value(value)}")
    if data is not None:
        if not isinstance(data, Mapping):
            raise TypeError("override 'data' must map dataset numbers (1-based) to keyword mappings")
        for n in sorted(data):
            if not isinstance(n, int) or n < 1:
                raise ValueError(f"dataset numbers in overrides are 1-based integers, got {n!r}")
            lines.append(f"data {n}")
            for key, value in data[n].items():
                lines.append(f"{key} {_fmt_value(value)}")
    return lines


def find_inp(folder: str | PathLike[str]) -> str:
    """Name of the single ``*.EPSR.inp`` file in ``folder``."""
    folder = Path(folder)
    found = sorted(p.name for p in folder.glob("*.EPSR.inp"))
    if len(found) != 1:
        what = "no" if not found else f"{len(found)}"
        raise FileNotFoundError(f"{what} *.EPSR.inp files in {folder}; pass inp_name explicitly")
    return found[0]


class Simulation:
    """One EPSR run on a simulation folder.

    Parameters
    ----------
    folder:
        Folder holding the ``.EPSR.inp``, ``.ato``, ``.pcof``, weights and data files.
    inp_name:
        Name of the ``.EPSR.inp`` file in ``folder``. Detected automatically when
        the folder contains exactly one.
    overrides:
        Parameter changes applied after reading the input file, see
        :func:`render_overrides`. They are saved back to the input files.
    verbose:
        Print the engine's progress messages. They are written by Rust directly
        to the process's standard output, not through ``sys.stdout``.

    Notes
    -----
    * **Files are modified in place** when the last iteration finishes. Use
      :meth:`epsr.Project.from_template` to work on a copy.
    * One object is one EPSR invocation: call :meth:`run` once, or
      :meth:`iterate` for ``1..niter`` in order. For the next run create a new
      ``Simulation`` (or use :func:`epsr.refine`).
    * ``run`` and ``iterate`` release the GIL. Simulations on *different*
      folders may run in parallel threads; never run two on the same folder.
    * Arrays are float32 copies (EPSR computes in single precision).
    """

    def __init__(
        self,
        folder: str | PathLike[str],
        inp_name: str | None = None,
        *,
        overrides: Overrides = None,
        verbose: bool = False,
    ) -> None:
        self.folder = Path(folder)
        self.inp_name = inp_name if inp_name is not None else find_inp(self.folder)
        self.overrides = render_overrides(overrides)
        self._sim = _core.Simulation(str(self.folder), self.inp_name, self.overrides, verbose)
        self._pair_labels: list[tuple[str, str]] = [tuple(p) for p in self._sim.pair_labels()]  # type: ignore[misc]
        self._pair_lookup = {}
        for i, (a, b) in enumerate(self._pair_labels):
            self._pair_lookup.setdefault((a, b), i)
            self._pair_lookup.setdefault((b, a), i)

    def __repr__(self) -> str:
        state = "completed" if self.completed else f"iteration {self.last_iteration}/{self.control.niter}"
        return f"<epsr.Simulation {self.folder / self.inp_name} ({state})>"

    # ---------------------------------------------------------------- running
    def run(self) -> Summary:
        """Run all ``niter`` iterations and write every output file (``Epsr::run``)."""
        self._sim.run()
        return self.summary

    def iterate(self, iteration: int) -> Summary:
        """Run iteration ``iteration`` (``1..niter``, in order); see :meth:`iterations`."""
        return Summary.from_fortran(self._sim.iterate(iteration))

    def iterations(self) -> Iterator[Summary]:
        """Run the remaining iterations one by one, yielding each :class:`Summary`.

        Equivalent to :meth:`run`: the output files are written after the last one.
        """
        for i in range(self.last_iteration + 1, self.control.niter + 1):
            yield self.iterate(i)

    @property
    def last_iteration(self) -> int:
        """Number of the last completed iteration (0 before the first)."""
        return self._sim.last_iter

    @property
    def completed(self) -> bool:
        """True once all ``niter`` iterations have run and the files are written."""
        return self._sim.last_iter >= self._sim.control["niter"]

    # ---------------------------------------------------------------- metadata
    @property
    def run_name(self) -> str:
        """The run name, e.g. ``200DMSO400Water.EPSR`` (prefix of every output file)."""
        return self.inp_name[: -len(".inp")] if self.inp_name.endswith(".inp") else self.inp_name

    def output_path(self, ext: str) -> Path:
        """Path of an output file, e.g. ``output_path("g01")``."""
        return self.folder / f"{self.run_name}.{ext.lstrip('.')}"

    @property
    def summary(self) -> Summary:
        """Results of the most recent iteration."""
        return Summary.from_fortran(self._sim.summary)

    @property
    def control(self) -> Control:
        return Control.from_fortran(self._sim.control)

    @property
    def datasets(self) -> list[Dataset]:
        return [Dataset.from_fortran(d) for d in self._sim.datasets]

    @property
    def n_datasets(self) -> int:
        return self._sim.n_datasets

    @property
    def n_pairs(self) -> int:
        return self._sim.n_pairs

    @property
    def nq(self) -> int:
        return self._sim.nq

    @property
    def nr(self) -> int:
        return self._sim.nr

    @property
    def n_atoms(self) -> int:
        return self._sim.n_atoms

    @property
    def n_molecules(self) -> int:
        return self._sim.n_molecules

    @property
    def density(self) -> float:
        """Atomic number density (atoms/Å³, ``rho``)."""
        return self._sim.rho

    @property
    def box_size(self) -> float:
        """Box side (Å, ``dmax``)."""
        return self._sim.box_size()

    @property
    def pair_labels(self) -> list[tuple[str, str]]:
        """Atom-type labels of each pair, in the order of every pair-resolved array."""
        return list(self._pair_labels)

    @property
    def atom_types(self) -> list[str]:
        return self._sim.atom_types()

    @property
    def atom_labels(self) -> list[str]:
        return self._sim.atom_labels()

    @property
    def molecule_atoms(self) -> list[range]:
        """For each molecule, the range of its atom indices (0-based)."""
        return [range(a, b) for a, b in self._sim.molecule_atoms()]

    @property
    def weights_files(self) -> list[str]:
        """Paths of the weights file of each dataset."""
        return self._sim.weights_files()

    @property
    def positions(self) -> F32:
        """``(n_atoms, 3)`` atom positions in Å."""
        return self._sim.atom_positions()

    def get_parameter(self, name: str, dataset: int | None = None) -> str | None:
        """Current value of an input keyword; ``dataset`` is a 0-based index."""
        return self._sim.get_parameter(name, dataset)

    def pair_index(self, a: str | PairKey, b: str | None = None) -> int:
        """Index of the pair ``(a, b)`` in :attr:`pair_labels` (order-insensitive).

        Accepts ``pair_index("O", "H")``, ``pair_index(("O", "H"))`` or an int.
        """
        if isinstance(a, (int, np.integer)) and b is None:
            i = int(a)
            if not -self.n_pairs <= i < self.n_pairs:
                raise IndexError(f"pair index {i} out of range (n_pairs = {self.n_pairs})")
            return i % self.n_pairs
        if isinstance(a, tuple) and b is None:
            a, b = a
        if not isinstance(a, str) or not isinstance(b, str):
            raise TypeError("select a pair by index, by two labels, or by a (label, label) tuple")
        try:
            return self._pair_lookup[(a, b)]
        except KeyError:
            raise KeyError(f"no atom-type pair {a}-{b}; pairs are {self._pair_labels}") from None

    # ---------------------------------------------------------------- grids
    @property
    def q(self) -> F32:
        """Q grid (Å⁻¹)."""
        return self._sim.q_grid()

    @property
    def r(self) -> F32:
        """r grid (Å)."""
        return self._sim.r_grid()

    # ---------------------------------------------------------------- partials
    def _select(self, arr: F32, pair: PairKey | None) -> F32:
        return arr if pair is None else arr[self.pair_index(pair)]

    def gr_inter(self, pair: PairKey | None = None) -> F32:
        """Running-average intermolecular partial g(r): ``(n_pairs, nr)``, or ``(nr,)`` for one pair."""
        return self._select(self._sim.gr_inter(), pair)

    def gr_intra(self, pair: PairKey | None = None) -> F32:
        """Running-average intramolecular partial g(r)."""
        return self._select(self._sim.gr_intra(), pair)

    def sq_inter(self, pair: PairKey | None = None) -> F32:
        """Running-average intermolecular partial S(Q): ``(n_pairs, nq)``."""
        return self._select(self._sim.sq_inter(), pair)

    def sq_intra(self, pair: PairKey | None = None) -> F32:
        """Running-average intramolecular partial S(Q)."""
        return self._select(self._sim.sq_intra(), pair)

    def current_gr(self) -> F32:
        """``(2, n_pairs, nr)`` histograms of the current configuration (0 = inter, 1 = intra)."""
        return self._sim.current_gr()

    def current_sq(self) -> F32:
        """``(2, n_pairs, nq)`` partial S(Q) of the current configuration (0 = inter, 1 = intra)."""
        return self._sim.current_sq()

    # ---------------------------------------------------------------- data
    def _dataset_index(self, i: int) -> int:
        if not -self.n_datasets <= i < self.n_datasets:
            raise IndexError(f"dataset index {i} out of range (n_datasets = {self.n_datasets})")
        return i % self.n_datasets

    def data_fq(self, i: int) -> F32:
        """Measured F(Q) of dataset ``i`` (0-based) on the Q grid; zero outside the data."""
        return self._sim.data_fq(self._dataset_index(i))

    def data_errors(self, i: int) -> F32:
        """Errors of the measured F(Q) of dataset ``i`` on the Q grid."""
        return self._sim.data_errors(self._dataset_index(i))

    def simulated_fq(self, i: int) -> F32:
        """Simulated total F(Q) of dataset ``i``."""
        return self._sim.simulated_fq(self._dataset_index(i))

    def data_q_range(self, i: int) -> tuple[int, int]:
        """First and last Q index (0-based, inclusive) fitted for dataset ``i``."""
        return self._sim.data_q_range(self._dataset_index(i))

    def fit(self, i: int) -> Fit:
        """Data, errors and simulation of dataset ``i`` over its fitted Q range."""
        i = self._dataset_index(i)
        lo, hi = self.data_q_range(i)
        s = slice(lo, hi + 1)
        return Fit(
            name=self.datasets[i].name,
            q=self.q[s],
            data=self.data_fq(i)[s],
            errors=self.data_errors(i)[s],
            simulated=self.simulated_fq(i)[s],
            q_range=(lo, hi),
        )

    # ---------------------------------------------------------------- potentials
    @property
    def potential(self) -> Potential:
        """Handle for reading and perturbing the potentials (:class:`epsr.Potential`)."""
        from .potential import Potential

        return Potential(self)
