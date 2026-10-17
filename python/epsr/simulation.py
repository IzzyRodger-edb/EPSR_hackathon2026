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
from typing import TYPE_CHECKING, Any

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

Float32Array = npt.NDArray[np.float32]
#: A pair is selected by its index in ``pair_labels`` or by its two atom-type labels.
PairKey = int | tuple[str, str]
#: Overrides: raw ``.EPSR.inp`` lines, or a mapping rendered by :func:`render_overrides`.
Overrides = Mapping[str, Any] | Sequence[str] | None


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
    def from_fortran(cls, fortran_values: Mapping[str, Any]) -> Summary:
        """Build from the engine's summary dict (Fortran variable names)."""
        return cls(
            iteration=int(fortran_values["iter"]),
            energy_inter=float(fortran_values["foldsv"]),
            energy_intra=float(fortran_values["foldsb"]),
            energy=float(fortran_values["energd"]),
            pressure=float(fortran_values["pressd"]),
            ep_amplitude=float(fortran_values["energabs"]),
            r_factor=float(fortran_values["rfactot"]),
            r_factor_raw=float(fortran_values["rfacraw"]),
            r_factor_per_dataset=tuple(
                float(value) for value in fortran_values["rfacd"]
            ),
            ep_up=float(fortran_values["epcup"]),
            ep_coefficient=float(fortran_values["epccoeff"]),
            gr_roughness=float(fortran_values["sentb"]),
            chisq=float(fortran_values["chisq"]),
            fortran=dict(fortran_values),
        )

    def as_dict(self) -> dict[str, Any]:
        """The descriptive fields as a plain dict (without :attr:`fortran`)."""
        values = asdict(self)
        values.pop("fortran")
        return values


def _from_fortran_dict(cls: Any, fortran_values: Mapping[str, Any]) -> Any:
    """Build a dataclass whose field names equal the engine's dict keys."""
    return cls(
        **{
            dataclass_field.name: fortran_values[dataclass_field.name]
            for dataclass_field in fields(cls)
        }
    )


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
    def from_fortran(cls, fortran_values: Mapping[str, Any]) -> Dataset:
        return _from_fortran_dict(cls, fortran_values)


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
    def from_fortran(cls, fortran_values: Mapping[str, Any]) -> Control:
        return _from_fortran_dict(cls, fortran_values)


@dataclass(frozen=True)
class Fit:
    """Data and simulated F(Q) of one dataset over its fitted Q range."""

    name: str
    q: Float32Array
    data: Float32Array
    errors: Float32Array
    simulated: Float32Array
    q_range: tuple[int, int]  # 0-based inclusive indices into the full Q grid

    @property
    def residual(self) -> Float32Array:
        return self.data - self.simulated


def _format_value(value: Any) -> str:
    """Write a Python value the way EPSR input files expect it."""
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, (list, tuple, np.ndarray)):
        return " ".join(_format_value(item) for item in value)
    if isinstance(value, (float, np.floating)):
        return repr(float(value))
    return str(value)


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
        if not all(isinstance(line, str) for line in lines):
            raise TypeError("override lines must be strings")
        return lines

    lines = []
    dataset_overrides = None
    for keyword, value in overrides.items():
        if keyword == "data":
            dataset_overrides = value
            continue
        if isinstance(value, Mapping):
            raise TypeError(
                f"override {keyword!r}: nested mappings are only allowed under 'data'"
            )
        lines.append(f"{keyword} {_format_value(value)}")

    if dataset_overrides is not None:
        if not isinstance(dataset_overrides, Mapping):
            raise TypeError(
                "override 'data' must map dataset numbers (1-based) to keyword mappings"
            )
        for dataset_number in sorted(dataset_overrides):
            if not isinstance(dataset_number, int) or dataset_number < 1:
                raise ValueError(
                    f"dataset numbers in overrides are 1-based integers, got {dataset_number!r}"
                )
            lines.append(f"data {dataset_number}")
            for keyword, value in dataset_overrides[dataset_number].items():
                lines.append(f"{keyword} {_format_value(value)}")
    return lines


def find_inp(folder: str | PathLike[str]) -> str:
    """Name of the single ``*.EPSR.inp`` file in ``folder``."""
    folder = Path(folder)
    input_files = sorted(path.name for path in folder.glob("*.EPSR.inp"))
    if len(input_files) != 1:
        how_many = "no" if not input_files else str(len(input_files))
        raise FileNotFoundError(
            f"{how_many} *.EPSR.inp files in {folder}; pass inp_name explicitly"
        )
    return input_files[0]


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
    * **Files are modified in place** when the last iteration finishes, and
      creating the object already rewrites the ``.qdr`` trend file. Use
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
        #: The Rust engine object (``epsr._core.Simulation``), with Fortran names.
        self._engine = _core.Simulation(
            str(self.folder), self.inp_name, self.overrides, verbose
        )

        self._pair_labels: list[tuple[str, str]] = [
            (first_type, second_type)
            for first_type, second_type in self._engine.pair_labels()
        ]
        # Both orders of a pair map to the same index.
        self._pair_index_by_labels: dict[tuple[str, str], int] = {}
        for index, (first_type, second_type) in enumerate(self._pair_labels):
            self._pair_index_by_labels.setdefault((first_type, second_type), index)
            self._pair_index_by_labels.setdefault((second_type, first_type), index)

    def __repr__(self) -> str:
        if self.completed:
            state = "completed"
        else:
            state = f"iteration {self.last_iteration}/{self.control.niter}"
        return f"<epsr.Simulation {self.folder / self.inp_name} ({state})>"

    # ---------------------------------------------------------------- running
    def run(self) -> Summary:
        """Run all ``niter`` iterations and write every output file (``Epsr::run``)."""
        self._engine.run()
        return self.summary

    def iterate(self, iteration: int) -> Summary:
        """Run iteration ``iteration`` (``1..niter``, in order); see :meth:`iterations`."""
        return Summary.from_fortran(self._engine.iterate(iteration))

    def iterations(self) -> Iterator[Summary]:
        """Run the remaining iterations one by one, yielding each :class:`Summary`.

        Equivalent to :meth:`run`: the output files are written after the last one.
        """
        for iteration in range(self.last_iteration + 1, self.control.niter + 1):
            yield self.iterate(iteration)

    @property
    def last_iteration(self) -> int:
        """Number of the last completed iteration (0 before the first)."""
        return self._engine.last_iter

    @property
    def completed(self) -> bool:
        """True once all ``niter`` iterations have run and the files are written."""
        return self._engine.last_iter >= self._engine.control["niter"]

    # ---------------------------------------------------------------- metadata
    @property
    def run_name(self) -> str:
        """The run name, e.g. ``200DMSO400Water.EPSR`` (prefix of every output file)."""
        return self.inp_name.removesuffix(".inp")

    def output_path(self, extension: str) -> Path:
        """Path of an output file, e.g. ``output_path("g01")``."""
        return self.folder / f"{self.run_name}.{extension.lstrip('.')}"

    @property
    def summary(self) -> Summary:
        """Results of the most recent iteration."""
        return Summary.from_fortran(self._engine.summary)

    @property
    def control(self) -> Control:
        return Control.from_fortran(self._engine.control)

    @property
    def datasets(self) -> list[Dataset]:
        return [Dataset.from_fortran(values) for values in self._engine.datasets]

    @property
    def n_datasets(self) -> int:
        return self._engine.n_datasets

    @property
    def n_pairs(self) -> int:
        return self._engine.n_pairs

    @property
    def nq(self) -> int:
        """Number of Q points."""
        return self._engine.nq

    @property
    def nr(self) -> int:
        """Number of r points."""
        return self._engine.nr

    @property
    def n_atoms(self) -> int:
        return self._engine.n_atoms

    @property
    def n_molecules(self) -> int:
        return self._engine.n_molecules

    @property
    def density(self) -> float:
        """Atomic number density (atoms/Å³, ``rho``)."""
        return self._engine.rho

    @property
    def box_size(self) -> float:
        """Box side (Å, ``dmax``)."""
        return self._engine.box_size()

    @property
    def pair_labels(self) -> list[tuple[str, str]]:
        """Atom-type labels of each pair, in the order of every pair-resolved array."""
        return list(self._pair_labels)

    @property
    def atom_types(self) -> list[str]:
        return self._engine.atom_types()

    @property
    def atom_labels(self) -> list[str]:
        return self._engine.atom_labels()

    @property
    def molecule_atoms(self) -> list[range]:
        """For each molecule, the range of its atom indices (0-based)."""
        return [range(start, stop) for start, stop in self._engine.molecule_atoms()]

    @property
    def weights_files(self) -> list[str]:
        """Paths of the weights file of each dataset."""
        return self._engine.weights_files()

    @property
    def positions(self) -> Float32Array:
        """``(n_atoms, 3)`` atom positions in Å."""
        return self._engine.atom_positions()

    def get_parameter(self, name: str, dataset: int | None = None) -> str | None:
        """Current value of an input keyword; ``dataset`` is a 0-based index."""
        return self._engine.get_parameter(name, dataset)

    def pair_index(
        self, pair_or_type: str | PairKey, other_type: str | None = None
    ) -> int:
        """Index of a pair in :attr:`pair_labels` (the order of the two types does not matter).

        Accepts ``pair_index("O1", "H1")``, ``pair_index(("O1", "H1"))`` or an int.
        """
        if isinstance(pair_or_type, (int, np.integer)) and other_type is None:
            index = int(pair_or_type)
            if not -self.n_pairs <= index < self.n_pairs:
                raise IndexError(
                    f"pair index {index} out of range (n_pairs = {self.n_pairs})"
                )
            return index % self.n_pairs

        labels: tuple[object, object]
        if isinstance(pair_or_type, tuple) and other_type is None:
            labels = pair_or_type
        else:
            labels = (pair_or_type, other_type)
        first_type, second_type = labels
        if not isinstance(first_type, str) or not isinstance(second_type, str):
            raise TypeError(
                "select a pair by index, by two labels, or by a (label, label) tuple"
            )
        try:
            return self._pair_index_by_labels[(first_type, second_type)]
        except KeyError:
            raise KeyError(
                f"no atom-type pair {first_type}-{second_type}; pairs are {self._pair_labels}"
            ) from None

    # ---------------------------------------------------------------- grids
    @property
    def q(self) -> Float32Array:
        """Q grid (Å⁻¹)."""
        return self._engine.q_grid()

    @property
    def r(self) -> Float32Array:
        """r grid (Å)."""
        return self._engine.r_grid()

    # ---------------------------------------------------------------- partials
    def _select_pair(
        self,
        all_pairs: Float32Array,
        pair: str | PairKey | None,
        other_type: str | None,
    ) -> Float32Array:
        """All rows, or the row of one pair."""
        if pair is None:
            return all_pairs
        return all_pairs[self.pair_index(pair, other_type)]

    def gr_inter(
        self, pair: str | PairKey | None = None, other_type: str | None = None
    ) -> Float32Array:
        """Running-average intermolecular partial g(r).

        ``(n_pairs, nr)``, or ``(nr,)`` for one pair given as an index, two labels
        (``gr_inter("O1", "H1")``) or a label tuple.
        """
        return self._select_pair(self._engine.gr_inter(), pair, other_type)

    def gr_intra(
        self, pair: str | PairKey | None = None, other_type: str | None = None
    ) -> Float32Array:
        """Running-average intramolecular partial g(r)."""
        return self._select_pair(self._engine.gr_intra(), pair, other_type)

    def sq_inter(
        self, pair: str | PairKey | None = None, other_type: str | None = None
    ) -> Float32Array:
        """Running-average intermolecular partial S(Q): ``(n_pairs, nq)``."""
        return self._select_pair(self._engine.sq_inter(), pair, other_type)

    def sq_intra(
        self, pair: str | PairKey | None = None, other_type: str | None = None
    ) -> Float32Array:
        """Running-average intramolecular partial S(Q)."""
        return self._select_pair(self._engine.sq_intra(), pair, other_type)

    def current_gr(self) -> Float32Array:
        """``(2, n_pairs, nr)`` histograms of the current configuration (0 = inter, 1 = intra)."""
        return self._engine.current_gr()

    def current_sq(self) -> Float32Array:
        """``(2, n_pairs, nq)`` partial S(Q) of the current configuration (0 = inter, 1 = intra)."""
        return self._engine.current_sq()

    # ---------------------------------------------------------------- data
    def _dataset_index(self, dataset: int) -> int:
        """Check a 0-based dataset index (negative values count from the end)."""
        if not -self.n_datasets <= dataset < self.n_datasets:
            raise IndexError(
                f"dataset index {dataset} out of range (n_datasets = {self.n_datasets})"
            )
        return dataset % self.n_datasets

    def data_fq(self, dataset: int) -> Float32Array:
        """Measured F(Q) of a dataset (0-based) on the Q grid; zero outside the data."""
        return self._engine.data_fq(self._dataset_index(dataset))

    def data_errors(self, dataset: int) -> Float32Array:
        """Errors of the measured F(Q) of a dataset on the Q grid."""
        return self._engine.data_errors(self._dataset_index(dataset))

    def simulated_fq(self, dataset: int) -> Float32Array:
        """Simulated total F(Q) of a dataset."""
        return self._engine.simulated_fq(self._dataset_index(dataset))

    def data_q_range(self, dataset: int) -> tuple[int, int]:
        """First and last Q index (0-based, inclusive) fitted for a dataset."""
        return self._engine.data_q_range(self._dataset_index(dataset))

    def fit(self, dataset: int) -> Fit:
        """Data, errors and simulation of a dataset over its fitted Q range."""
        dataset = self._dataset_index(dataset)
        first_index, last_index = self.data_q_range(dataset)
        fitted_range = slice(first_index, last_index + 1)
        return Fit(
            name=self.datasets[dataset].name,
            q=self.q[fitted_range],
            data=self.data_fq(dataset)[fitted_range],
            errors=self.data_errors(dataset)[fitted_range],
            simulated=self.simulated_fq(dataset)[fitted_range],
            q_range=(first_index, last_index),
        )

    # ---------------------------------------------------------------- potentials
    @property
    def potential(self) -> Potential:
        """Handle for reading and perturbing the potentials (:class:`epsr.Potential`)."""
        from .potential import Potential

        return Potential(self)
