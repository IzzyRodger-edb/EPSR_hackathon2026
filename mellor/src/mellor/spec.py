"""The system specification: everything Mellor needs to set up one EPSR refinement.

A :class:`SystemSpec` is the contract between the language model and the code. The
model (or the built-in library, or a user) produces one as JSON; ``mellor.inputs``
turns it into EPSR input files and ``mellor.script`` into a standalone Python script.
It is deliberately explicit: every number that ends up in an EPSR file is here.

Units: lengths in Å, energies in kJ/mol, charges in e, masses in amu, densities in
atoms/Å³, temperatures in K. Atom indices inside a molecule are 1-based (as in the
EPSR ``.mol`` files). Mass number 0 means the natural isotopic mixture.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any

__all__ = [
    "LABEL_MAX_LENGTH",
    "KINDS",
    "EPSR_KEYWORDS",
    "SpecError",
    "AtomType",
    "Atom",
    "Molecule",
    "Isotope",
    "Dataset",
    "Schedule",
    "PdbSetup",
    "SystemSpec",
    "SPEC_JSON_SCHEMA",
    "spec_from_dict",
    "spec_to_dict",
]

#: ``makemole`` truncates atom-type labels to three characters.
LABEL_MAX_LENGTH = 3
#: How the box is built.
KINDS = ("molecular", "atomic", "pdb")
#: Global keywords of the ``.EPSR.inp`` template that a specification may override
#: (the values of ``validation/water/waterbox.EPSR.inp`` are the defaults).
EPSR_KEYWORDS = (
    "feedback",
    "potfac",
    "rspcrmin",
    "rspcfrac",
    "porodpwr",
    "nmolcell",
    "control-p",
    "ref_intra",
    "sizefactor",
    "nq",
    "qstep",
    "qinscalefc",
    "ntimes",
    "intra",
    "rotfreq",
    "inter",
    "cellst",
    "rmaxgr",
    "ngrsamples",
    "fwhm",
    "fwhmq",
    "nsmoop",
    "revlorch",
    "qwidthqmax",
    "mplicities",
    "hklqmin",
    "diffuse",
    "nlatrefine",
    "latrefine",
    "rejrate",
    "qmin",
)
#: Per-dataset keywords of a ``data n`` block (besides ``datafile`` and ``wtsfile``).
DATASET_KEYWORDS = ("nrtype", "rshmin", "szeros", "tweak", "efilereq", "porodadd")

_LABEL_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9]{0,2}$")
_NAME_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_-]{0,23}$")
_ELEMENT_PATTERN = re.compile(r"^[A-Z][a-z]?$")


class SpecError(ValueError):
    """A specification that cannot be run. ``errors`` lists every problem found."""

    def __init__(self, errors: list[str]):
        self.errors = list(errors)
        super().__init__("invalid system specification:\n  - " + "\n  - ".join(errors))


# --------------------------------------------------------------------------- parts
@dataclass(frozen=True)
class AtomType:
    """An EPSR atom type with its reference-potential parameters.

    ``nsph``, ``rmin``, ``rmax`` and ``cg`` are the four trailing values of a
    ``potential`` line in a ``.mol`` file; ``makemole`` defaults them to
    ``1, 0, 0.5*sigma, 1`` when they are absent (``rmax`` None means that default).
    """

    label: str
    element: str
    epsilon: float = 0.0
    sigma: float = 0.0
    mass: float = 0.0
    charge: float = 0.0
    nsph: int = 1
    rmin: float = 0.0
    rmax: float | None = None
    cg: float = 1.0

    @property
    def rmax_value(self) -> float:
        return 0.5 * self.sigma if self.rmax is None else self.rmax


@dataclass(frozen=True)
class Atom:
    """One atom of a molecule template: its type label and coordinates (Å)."""

    label: str
    x: float
    y: float
    z: float

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass
class Molecule:
    """A molecule template (one ``.mol`` file) and how many copies go into the box.

    ``bonds`` are 1-based atom-index pairs. ``bond_lengths`` and ``angles`` may give
    the values written to the ``.mol`` file explicitly (keyed ``"i-j"`` and
    ``"i-j-k"``, 1-based, ``j`` the vertex); otherwise they are computed from the
    coordinates. ``rotations`` are ``rot`` lines: bonds about which internal
    rotation moves are made (1-based pairs).
    """

    name: str
    atoms: list[Atom]
    count: int
    bonds: list[tuple[int, int]] = field(default_factory=list)
    bond_lengths: dict[str, float] = field(default_factory=dict)
    angles: dict[str, float] = field(default_factory=dict)
    rotations: list[tuple[int, int]] = field(default_factory=list)

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def labels(self) -> list[str]:
        return [atom.label for atom in self.atoms]

    def neighbours(self) -> dict[int, list[int]]:
        """Bonded neighbours of every atom (1-based), in bond order."""
        result: dict[int, list[int]] = {
            index: [] for index in range(1, self.n_atoms + 1)
        }
        for first, second in self.bonds:
            if second not in result[first]:
                result[first].append(second)
            if first not in result[second]:
                result[second].append(first)
        return result

    def distance(self, first: int, second: int) -> float:
        a, b = self.atoms[first - 1], self.atoms[second - 1]
        return math.dist(a.xyz, b.xyz)

    def bond_length(self, first: int, second: int) -> float:
        """Length written to the ``.mol`` file for the bond ``first-second``."""
        for key in (f"{first}-{second}", f"{second}-{first}"):
            if key in self.bond_lengths:
                return float(self.bond_lengths[key])
        return self.distance(first, second)

    def angle_triples(self) -> list[tuple[int, int, int]]:
        """Every ``(i, j, k)`` with ``i`` and ``k`` bonded to the vertex ``j`` (``i < k``).

        Ordered by vertex, then by the neighbour order of the bond list; this is the
        set of ``angle`` lines ``makemole`` needs (it prompts for any missing one).
        """
        triples = []
        neighbours = self.neighbours()
        for vertex in range(1, self.n_atoms + 1):
            around = neighbours[vertex]
            for position, first in enumerate(around):
                for second in around[position + 1 :]:
                    triples.append((first, vertex, second))
        return triples

    def angle_value(self, first: int, vertex: int, second: int) -> float:
        """Angle in degrees written to the ``.mol`` file for ``first-vertex-second``."""
        for key in (f"{first}-{vertex}-{second}", f"{second}-{vertex}-{first}"):
            if key in self.angles:
                return float(self.angles[key])
        a = self.atoms[first - 1].xyz
        b = self.atoms[vertex - 1].xyz
        c = self.atoms[second - 1].xyz
        u = (a[0] - b[0], a[1] - b[1], a[2] - b[2])
        v = (c[0] - b[0], c[1] - b[1], c[2] - b[2])
        norm = math.hypot(*u) * math.hypot(*v)
        if norm == 0.0:
            raise SpecError(
                [
                    f"molecule {self.name}: coincident atoms in angle {first}-{vertex}-{second}"
                ]
            )
        cosine = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / norm
        return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))

    def is_connected(self) -> bool:
        if self.n_atoms <= 1:
            return True
        neighbours = self.neighbours()
        seen = {1}
        stack = [1]
        while stack:
            current = stack.pop()
            for other in neighbours[current]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        return len(seen) == self.n_atoms


@dataclass(frozen=True)
class Isotope:
    """One isotope of an atom type in a sample: mass number (0 = natural) and fraction."""

    mass_number: int
    fraction: float


@dataclass
class Dataset:
    """One measured data set (a ``.mint01`` file) and how the sample was labelled.

    ``isotopes`` maps atom-type labels to the isotopic composition of that type in
    this sample; types that are missing are natural. ``exchangeable`` lists the
    types whose hydrogen exchanges between molecules (``iexchange 1`` in ``nwts``).
    The remaining fields are the per-dataset keywords of the ``.EPSR.inp``;
    ``rshmin`` None means "choose from the system" (see ``inputs.auto_rshmin``).
    """

    file: str
    label: str
    isotopes: dict[str, list[Isotope]] = field(default_factory=dict)
    exchangeable: list[str] = field(default_factory=list)
    nrtype: int = 5
    rshmin: float | None = None
    szeros: float = 0.0
    tweak: float = 1.0
    efilereq: float = 1.0
    porodadd: float = 0.0

    @property
    def file_name(self) -> str:
        return Path(self.file).name

    @property
    def weights_file(self) -> str:
        """The ``.NWTStot.wts`` file ``nwts`` writes for this sample."""
        return f"{self.label}.NWTStot.wts"

    @property
    def nwts_file(self) -> str:
        return f"{self.label}.NWTS.dat"

    def composition(self, label: str) -> list[Isotope]:
        return list(self.isotopes.get(label) or [Isotope(0, 1.0)])


@dataclass
class Schedule:
    """The sequence of EPSR runs.

    1. scratch start: one iteration without a ``.pcof`` (creates it, ``ereq`` 0);
    2. ``equilibration_iterations`` on the reference potential only (``ereq`` 0);
    3. ``refinement_runs`` runs of ``refinement_iterations`` with ``ereq`` set.

    ``quick_*`` are the iteration counts used with ``--quick`` (tests, previews).
    """

    ereq: float
    equilibration_iterations: int = 99
    refinement_runs: int = 2
    refinement_iterations: int = 100
    ntimes: int = 5
    fmole_iterations: int = 100
    randomise: bool = True
    quick_equilibration_iterations: int = 10
    quick_refinement_iterations: int = 10

    def iterations(self, quick: bool = False) -> list[tuple[str, int, float]]:
        """``(stage, niter, ereq)`` of every EPSR run in order."""
        equilibration = (
            self.quick_equilibration_iterations
            if quick
            else self.equilibration_iterations
        )
        refinement = (
            self.quick_refinement_iterations if quick else self.refinement_iterations
        )
        runs = [("scratch", 1, 0.0), ("reference potential", equilibration, 0.0)]
        runs += [("EP refinement", refinement, self.ereq)] * self.refinement_runs
        return runs

    def total_iterations(self, quick: bool = False) -> int:
        return sum(niter for _, niter, _ in self.iterations(quick))


@dataclass
class PdbSetup:
    """How a ``pdb`` system is built: ``readpdb`` then ``changeato``.

    ``file`` is the PDB file (atom names are the EPSR atom-type labels, residues the
    molecules). ``changeato`` receives the temperature, the Monte Carlo step sizes and
    the ``atom_types`` parameters in PDB order; ``extra_commands`` are appended before
    ``save``.
    """

    file: str
    vdw_factor: float = 0.68
    dihedrals: bool = False
    tolerance: float = 0.0
    step_sizes: dict[str, float] = field(
        default_factory=lambda: {
            "stepmi": 0.1,
            "stepri": 0.3,
            "stepr": 0.3,
            "stepm": 1.0,
        }
    )
    extra_commands: list[str] = field(default_factory=list)

    @property
    def file_name(self) -> str:
        return Path(self.file).name


# --------------------------------------------------------------------------- the spec
@dataclass
class SystemSpec:
    """A complete description of one EPSR refinement. See the module docstring."""

    name: str
    title: str
    kind: str
    atom_types: list[AtomType]
    molecules: list[Molecule]
    density: float
    temperature: float
    datasets: list[Dataset]
    schedule: Schedule
    epsr_keywords: dict[str, Any] = field(default_factory=dict)
    pdb: PdbSetup | None = None
    notes: str = ""
    source: str = ""

    # ---- derived quantities
    @property
    def n_molecules(self) -> int:
        return sum(molecule.count for molecule in self.molecules)

    @property
    def n_atoms(self) -> int:
        return sum(molecule.count * molecule.n_atoms for molecule in self.molecules)

    @property
    def box_side(self) -> float:
        """Side of the cubic box ``mixato`` builds for this density (Å)."""
        return (self.n_atoms / self.density) ** (1.0 / 3.0)

    @property
    def labels(self) -> list[str]:
        return [atom_type.label for atom_type in self.atom_types]

    def atom_type(self, label: str) -> AtomType:
        for atom_type in self.atom_types:
            if atom_type.label == label:
                return atom_type
        raise KeyError(label)

    def type_counts(self) -> dict[str, int]:
        """Number of atoms of each type in the box."""
        counts = dict.fromkeys(self.labels, 0)
        for molecule in self.molecules:
            for atom in molecule.atoms:
                counts[atom.label] = counts.get(atom.label, 0) + molecule.count
        return counts

    def type_fractions(self) -> dict[str, float]:
        counts = self.type_counts()
        total = sum(counts.values()) or 1
        return {label: count / total for label, count in counts.items()}

    def mass_density(self) -> float:
        """Mass density in g/cm³ from the atomic number density and the masses."""
        fractions = self.type_fractions()
        mean_mass = sum(
            self.atom_type(label).mass * fraction
            for label, fraction in fractions.items()
        )
        return self.density * mean_mass * 1.66053907  # amu/Å³ -> g/cm³

    @property
    def inp_name(self) -> str:
        return f"{self.name}.EPSR.inp"

    @property
    def box_name(self) -> str:
        """Name of the box: ``<name>.ato`` and ``<name>.pcof``."""
        return self.name

    # ---- validation
    def problems(self, *, check_files: bool = False) -> tuple[list[str], list[str]]:
        """``(errors, warnings)``. Errors make the specification unusable."""
        errors: list[str] = []
        warnings: list[str] = []

        if not _NAME_PATTERN.match(self.name or ""):
            errors.append(
                "name must be 1-24 letters, digits, '-' or '_' and start with a letter"
            )
        if self.kind not in KINDS:
            errors.append(f"kind must be one of {KINDS}, not {self.kind!r}")
        if not self.title:
            warnings.append("title is empty")

        seen: set[str] = set()
        for atom_type in self.atom_types:
            label = atom_type.label
            if not _LABEL_PATTERN.match(label or ""):
                errors.append(
                    f"atom type label {label!r} must be 1-{LABEL_MAX_LENGTH} letters/digits starting with a letter"
                )
            if label in seen:
                errors.append(f"duplicate atom type label {label!r}")
            seen.add(label)
            if not _ELEMENT_PATTERN.match(atom_type.element or ""):
                errors.append(
                    f"atom type {label}: element {atom_type.element!r} is not a chemical symbol"
                )
            if atom_type.epsilon < 0 or atom_type.sigma < 0:
                errors.append(f"atom type {label}: epsilon and sigma must be >= 0")
            if atom_type.mass <= 0:
                errors.append(f"atom type {label}: mass must be > 0")
            if abs(atom_type.charge) > 4:
                warnings.append(
                    f"atom type {label}: unusual charge {atom_type.charge} e"
                )
        if not self.atom_types:
            errors.append("no atom types")

        if not self.molecules and self.kind != "pdb":
            errors.append("no molecules")
        used: set[str] = set()
        molecule_names: set[str] = set()
        for molecule in self.molecules:
            if not _NAME_PATTERN.match(molecule.name or ""):
                errors.append(
                    f"molecule name {molecule.name!r} must be letters, digits, '-' or '_'"
                )
            if molecule.name in molecule_names:
                errors.append(f"duplicate molecule name {molecule.name!r}")
            molecule_names.add(molecule.name)
            if molecule.count < 1:
                errors.append(f"molecule {molecule.name}: count must be >= 1")
            if not molecule.atoms:
                errors.append(f"molecule {molecule.name}: no atoms")
            for atom in molecule.atoms:
                used.add(atom.label)
                if atom.label not in seen:
                    errors.append(
                        f"molecule {molecule.name}: atom type {atom.label!r} is not defined"
                    )
            n = molecule.n_atoms
            for first, second in molecule.bonds:
                if not (1 <= first <= n and 1 <= second <= n) or first == second:
                    errors.append(
                        f"molecule {molecule.name}: bad bond {first}-{second} (atoms are 1..{n})"
                    )
                elif molecule.distance(first, second) > 3.0:
                    warnings.append(
                        f"molecule {molecule.name}: bond {first}-{second} is {molecule.distance(first, second):.2f} Å long"
                    )
            for first, second in molecule.rotations:
                if (first, second) not in molecule.bonds and (
                    second,
                    first,
                ) not in molecule.bonds:
                    errors.append(
                        f"molecule {molecule.name}: rotation {first}-{second} is not a bond"
                    )
            if self.kind == "atomic" and n != 1:
                errors.append(
                    f"kind 'atomic' needs single-atom molecules; {molecule.name} has {n}"
                )
            if self.kind == "molecular" and n > 1 and not molecule.bonds:
                errors.append(f"molecule {molecule.name}: {n} atoms but no bonds")
            if n > 1 and molecule.bonds and not molecule.is_connected():
                errors.append(
                    f"molecule {molecule.name}: the bonds do not connect all atoms"
                )
            for position, atom in enumerate(molecule.atoms, start=1):
                for other_position, other in enumerate(
                    molecule.atoms[: position - 1], start=1
                ):
                    if math.dist(atom.xyz, other.xyz) < 0.3:
                        errors.append(
                            f"molecule {molecule.name}: atoms {other_position} and {position} are closer than 0.3 Å"
                        )
        if self.kind != "pdb":
            for label in self.labels:
                if label not in used:
                    errors.append(
                        f"atom type {label} is defined but used by no molecule"
                    )

        if self.kind == "pdb":
            if self.pdb is None:
                errors.append("kind 'pdb' needs the 'pdb' setup")
            elif check_files and not Path(self.pdb.file).is_file():
                errors.append(f"PDB file {self.pdb.file} not found")
        elif self.pdb is not None:
            warnings.append("'pdb' setup is ignored for kind " + self.kind)

        if not (self.density > 0):
            errors.append("density must be > 0 atoms/Å³")
        elif not (0.005 <= self.density <= 0.2):
            warnings.append(
                f"density {self.density} atoms/Å³ is unusual (water is 0.1, silica 0.066)"
            )
        if not (self.temperature > 0):
            errors.append("temperature must be > 0 K")

        if not self.datasets:
            errors.append("no datasets")
        labels_seen: set[str] = set()
        for dataset in self.datasets:
            if not _NAME_PATTERN.match(dataset.label or ""):
                errors.append(
                    f"dataset label {dataset.label!r} must be letters, digits, '-' or '_'"
                )
            if dataset.label in labels_seen:
                errors.append(f"duplicate dataset label {dataset.label!r}")
            labels_seen.add(dataset.label)
            if not dataset.file:
                errors.append(f"dataset {dataset.label}: no file")
            elif check_files and not Path(dataset.file).is_file():
                errors.append(f"dataset {dataset.label}: file {dataset.file} not found")
            for label, isotopes in dataset.isotopes.items():
                if label not in seen:
                    errors.append(
                        f"dataset {dataset.label}: isotopes given for unknown atom type {label!r}"
                    )
                total = sum(isotope.fraction for isotope in isotopes)
                if not isotopes or abs(total - 1.0) > 1e-6:
                    errors.append(
                        f"dataset {dataset.label}: isotope fractions of {label} sum to {total:g}, not 1"
                    )
                for isotope in isotopes:
                    if isotope.mass_number < 0 or isotope.fraction < 0:
                        errors.append(
                            f"dataset {dataset.label}: negative mass number or fraction for {label}"
                        )
            for label in dataset.exchangeable:
                if label not in seen:
                    errors.append(
                        f"dataset {dataset.label}: exchangeable type {label!r} is not defined"
                    )
            if dataset.nrtype not in (1, 2, 3, 4, 5):
                errors.append(f"dataset {dataset.label}: nrtype must be 1-5")
            if dataset.rshmin is not None and dataset.rshmin < 0:
                errors.append(f"dataset {dataset.label}: rshmin must be >= 0")

        schedule = self.schedule
        if schedule.ereq < 0:
            errors.append("schedule.ereq must be >= 0")
        if schedule.equilibration_iterations < 1 or schedule.refinement_iterations < 1:
            errors.append("schedule iterations must be >= 1")
        if schedule.refinement_runs < 0:
            errors.append("schedule.refinement_runs must be >= 0")
        if schedule.refinement_runs > 0 and schedule.ereq == 0:
            warnings.append(
                "ereq is 0: the empirical potential stays off during the refinement runs"
            )
        if schedule.ntimes < 1 or schedule.fmole_iterations < 0:
            errors.append("schedule.ntimes must be >= 1 and fmole_iterations >= 0")

        for keyword in self.epsr_keywords:
            if keyword not in EPSR_KEYWORDS:
                errors.append(
                    f"unknown .EPSR.inp keyword {keyword!r}; allowed: {', '.join(EPSR_KEYWORDS)}"
                )
        return errors, warnings

    def validate(self, *, check_files: bool = False) -> list[str]:
        """Raise :class:`SpecError` on errors; return the warnings."""
        errors, warnings = self.problems(check_files=check_files)
        if errors:
            raise SpecError(errors)
        return warnings

    # ---- JSON
    def to_dict(self) -> dict[str, Any]:
        return spec_to_dict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.write_text(self.to_json() + "\n")
        return path

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> SystemSpec:
        return spec_from_dict(data)

    @classmethod
    def from_json(cls, text: str) -> SystemSpec:
        return spec_from_dict(json.loads(text))

    @classmethod
    def load(cls, path: str | Path) -> SystemSpec:
        return cls.from_json(Path(path).read_text())


# --------------------------------------------------------------------------- JSON conversion
def _plain(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return {f.name: _plain(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value


def spec_to_dict(spec: SystemSpec) -> dict[str, Any]:
    """A JSON-serialisable dict (isotopes as ``[mass_number, fraction]`` pairs)."""
    data = _plain(spec)
    for dataset in data["datasets"]:
        dataset["isotopes"] = {
            label: [
                [isotope["mass_number"], isotope["fraction"]] for isotope in isotopes
            ]
            for label, isotopes in dataset["isotopes"].items()
        }
    return data


def _pairs(value: Iterable[Any] | None, what: str) -> list[tuple[int, int]]:
    result = []
    for item in value or []:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            raise SpecError([f"{what}: expected [i, j] pairs, got {item!r}"])
        result.append((int(item[0]), int(item[1])))
    return result


def _isotopes(value: Any, where: str) -> list[Isotope]:
    result = []
    for item in value or []:
        if isinstance(item, Mapping):
            result.append(Isotope(int(item["mass_number"]), float(item["fraction"])))
        elif isinstance(item, (list, tuple)) and len(item) == 2:
            result.append(Isotope(int(item[0]), float(item[1])))
        else:
            raise SpecError(
                [
                    f"{where}: isotopes must be [mass_number, fraction] pairs, got {item!r}"
                ]
            )
    return result


def _only(data: Mapping[str, Any], cls: type, where: str) -> dict[str, Any]:
    allowed = {f.name for f in fields(cls)}
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise SpecError(
            [f"{where}: unknown field(s) {unknown}; allowed: {sorted(allowed)}"]
        )
    return dict(data)


def spec_from_dict(data: Mapping[str, Any]) -> SystemSpec:
    """Build a :class:`SystemSpec` from JSON data (tolerant of missing optional fields).

    Raises :class:`SpecError` for unknown fields or malformed values; call
    :meth:`SystemSpec.validate` afterwards for the physics checks.
    """
    if not isinstance(data, Mapping):
        raise SpecError(["the specification must be a JSON object"])
    try:
        top = _only(data, SystemSpec, "spec")
        atom_types = [
            AtomType(**_only(item, AtomType, f"atom_types[{i}]"))
            for i, item in enumerate(top.get("atom_types") or [])
        ]
        molecules = []
        for i, item in enumerate(top.get("molecules") or []):
            item = _only(item, Molecule, f"molecules[{i}]")
            atoms = []
            for j, atom in enumerate(item.get("atoms") or []):
                if isinstance(atom, Mapping):
                    atoms.append(
                        Atom(**_only(atom, Atom, f"molecules[{i}].atoms[{j}]"))
                    )
                elif isinstance(atom, (list, tuple)) and len(atom) == 4:
                    atoms.append(
                        Atom(
                            str(atom[0]), float(atom[1]), float(atom[2]), float(atom[3])
                        )
                    )
                else:
                    raise SpecError(
                        [f"molecules[{i}].atoms[{j}]: expected {{label, x, y, z}}"]
                    )
            molecules.append(
                Molecule(
                    name=str(item.get("name", "")),
                    atoms=atoms,
                    count=int(item.get("count", 0)),
                    bonds=_pairs(item.get("bonds"), f"molecules[{i}].bonds"),
                    bond_lengths={
                        str(k): float(v)
                        for k, v in (item.get("bond_lengths") or {}).items()
                    },
                    angles={
                        str(k): float(v) for k, v in (item.get("angles") or {}).items()
                    },
                    rotations=_pairs(
                        item.get("rotations"), f"molecules[{i}].rotations"
                    ),
                )
            )
        datasets = []
        for i, item in enumerate(top.get("datasets") or []):
            item = _only(item, Dataset, f"datasets[{i}]")
            isotopes = {
                str(label): _isotopes(value, f"datasets[{i}].isotopes[{label}]")
                for label, value in (item.get("isotopes") or {}).items()
            }
            datasets.append(
                Dataset(
                    file=str(item.get("file", "")),
                    label=str(item.get("label", "")),
                    isotopes=isotopes,
                    exchangeable=[str(v) for v in item.get("exchangeable") or []],
                    nrtype=int(item.get("nrtype", 5)),
                    rshmin=(
                        None if item.get("rshmin") is None else float(item["rshmin"])
                    ),
                    szeros=float(item.get("szeros", 0.0)),
                    tweak=float(item.get("tweak", 1.0)),
                    efilereq=float(item.get("efilereq", 1.0)),
                    porodadd=float(item.get("porodadd", 0.0)),
                )
            )
        schedule_data = top.get("schedule")
        if not isinstance(schedule_data, Mapping) or "ereq" not in schedule_data:
            raise SpecError(
                [
                    "schedule.ereq (requested empirical-potential amplitude, kJ/mol) is required"
                ]
            )
        schedule = Schedule(**_only(schedule_data, Schedule, "schedule"))
        pdb_data = top.get("pdb")
        pdb = (
            PdbSetup(**_only(pdb_data, PdbSetup, "pdb"))
            if isinstance(pdb_data, Mapping)
            else None
        )
        return SystemSpec(
            name=str(top.get("name", "")),
            title=str(top.get("title", "")),
            kind=str(top.get("kind", "molecular")),
            atom_types=atom_types,
            molecules=molecules,
            density=float(top.get("density", 0.0)),
            temperature=float(top.get("temperature", 0.0)),
            datasets=datasets,
            schedule=schedule,
            epsr_keywords=dict(top.get("epsr_keywords") or {}),
            pdb=pdb,
            notes=str(top.get("notes", "")),
            source=str(top.get("source", "")),
        )
    except SpecError:
        raise
    except (TypeError, ValueError, KeyError) as error:
        raise SpecError([f"malformed specification: {error}"]) from error


# --------------------------------------------------------------------------- JSON schema (for the model)
SPEC_JSON_SCHEMA: dict[str, Any] = {
    "type": "object",
    "description": "A complete EPSR system specification (units: Å, kJ/mol, e, amu, atoms/Å³, K).",
    "required": [
        "name",
        "title",
        "kind",
        "atom_types",
        "molecules",
        "density",
        "temperature",
        "datasets",
        "schedule",
    ],
    "properties": {
        "name": {
            "type": "string",
            "description": "Run name, e.g. 'waterbox' (letters, digits, - or _).",
        },
        "title": {
            "type": "string",
            "description": "One-line description of the system and the data.",
        },
        "kind": {
            "type": "string",
            "enum": list(KINDS),
            "description": "molecular: makemole/mixato/randomise/fmole; atomic: single-atom molecules (no fmole); pdb: readpdb + changeato on a PDB file.",
        },
        "atom_types": {
            "type": "array",
            "description": "EPSR atom types with reference-potential parameters.",
            "items": {
                "type": "object",
                "required": ["label", "element", "epsilon", "sigma", "mass", "charge"],
                "properties": {
                    "label": {
                        "type": "string",
                        "maxLength": LABEL_MAX_LENGTH,
                        "description": "At most 3 characters, e.g. OW, HW, Si.",
                    },
                    "element": {"type": "string"},
                    "epsilon": {
                        "type": "number",
                        "description": "Lennard-Jones well depth, kJ/mol (0 for hydrogen).",
                    },
                    "sigma": {
                        "type": "number",
                        "description": "Lennard-Jones diameter, Å.",
                    },
                    "mass": {
                        "type": "number",
                        "description": "amu (2.0 for H so EPSR can treat it as H or D).",
                    },
                    "charge": {"type": "number", "description": "partial charge, e."},
                },
            },
        },
        "molecules": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["name", "atoms", "count"],
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "Template name (file <name>.mol).",
                    },
                    "count": {
                        "type": "integer",
                        "description": "Number of molecules in the box.",
                    },
                    "atoms": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "required": ["label", "x", "y", "z"],
                            "properties": {
                                "label": {"type": "string"},
                                "x": {"type": "number"},
                                "y": {"type": "number"},
                                "z": {"type": "number"},
                            },
                        },
                    },
                    "bonds": {
                        "type": "array",
                        "items": {
                            "type": "array",
                            "items": {"type": "integer"},
                            "minItems": 2,
                            "maxItems": 2,
                        },
                        "description": "1-based atom index pairs.",
                    },
                    "rotations": {
                        "type": "array",
                        "items": {"type": "array", "items": {"type": "integer"}},
                        "description": "Bonds about which internal rotations are made (optional).",
                    },
                },
            },
        },
        "density": {
            "type": "number",
            "description": "Atomic number density, atoms/Å³ (water 0.1, silica 0.0664).",
        },
        "temperature": {"type": "number", "description": "K"},
        "datasets": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["file", "label", "isotopes"],
                "properties": {
                    "file": {
                        "type": "string",
                        "description": "Path of the .mint01 file.",
                    },
                    "label": {
                        "type": "string",
                        "description": "Short sample name, e.g. D2O (used for the weights files).",
                    },
                    "isotopes": {
                        "type": "object",
                        "description": 'Atom-type label -> list of [mass_number, fraction]; mass number 0 = natural, e.g. {"HW": [[2, 1.0]]} for D2O, {"HW": [[1, 0.5], [2, 0.5]]} for a 50:50 H/D mixture. Omitted types are natural.',
                        "additionalProperties": {
                            "type": "array",
                            "items": {
                                "type": "array",
                                "items": {"type": "number"},
                                "minItems": 2,
                                "maxItems": 2,
                            },
                        },
                    },
                    "exchangeable": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Atom types whose H/D exchanges between molecules (e.g. water HW).",
                    },
                    "nrtype": {
                        "type": "integer",
                        "description": "5 for Gudrun .mint01 files.",
                    },
                    "rshmin": {
                        "type": "number",
                        "description": "Minimum r (Å) for background subtraction; below the shortest intramolecular distance. Leave out to let Mellor choose.",
                    },
                    "tweak": {"type": "number"},
                    "efilereq": {"type": "number"},
                },
            },
        },
        "schedule": {
            "type": "object",
            "required": ["ereq"],
            "properties": {
                "ereq": {
                    "type": "number",
                    "description": "Requested empirical-potential amplitude, kJ/mol (water 3, benzene 20, silica 50).",
                },
                "equilibration_iterations": {
                    "type": "integer",
                    "description": "Iterations on the reference potential (default 99).",
                },
                "refinement_runs": {
                    "type": "integer",
                    "description": "Number of EP refinement runs (default 2).",
                },
                "refinement_iterations": {
                    "type": "integer",
                    "description": "Iterations per refinement run (default 100).",
                },
                "ntimes": {
                    "type": "integer",
                    "description": "MC cycles between potential refinements (default 5).",
                },
                "fmole_iterations": {
                    "type": "integer",
                    "description": "Geometry relaxation iterations (default 100).",
                },
            },
        },
        "epsr_keywords": {
            "type": "object",
            "description": 'Overrides of .EPSR.inp keywords, e.g. {"qmin": [0.1, 0.0], "fwhmq": 0.02, "qwidthqmax": [0.2, 5.0]}.',
        },
        "pdb": {
            "type": "object",
            "description": "Only for kind 'pdb'.",
            "properties": {
                "file": {"type": "string"},
                "vdw_factor": {"type": "number"},
                "dihedrals": {"type": "boolean"},
            },
        },
        "notes": {
            "type": "string",
            "description": "Why this model was chosen (shown in the report).",
        },
    },
}
