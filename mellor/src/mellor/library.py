"""Built-in systems: validated EPSR models the model can start from.

Each entry knows how to turn a set of ``.mint01`` files into a :class:`SystemSpec`,
assigning the isotopic composition of every data set from its title (``H2O``,
``D2O``, ``HDO``...). The water, silica and Cu-BTC entries reproduce the inputs of
``validation/`` exactly; DMSO/water comes from ``EPSR-test``; benzene from the
Dissolve benzene example (``examples/benzene/input.txt`` of the Dissolve
repository). No argon entry: the repository holds no validated argon inputs, and
the library copies numbers, it does not invent them.

The entries are the physics the model must not invent: geometries, Lennard-Jones
parameters, charges, densities, temperatures and the ``ereq`` values. Every number
below is copied from the file named next to it; ``tests/test_library.py`` parses
those files and checks the copies.
"""

from __future__ import annotations

import copy
import functools
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .data import Mint01, guess_sample
from .spec import (
    Atom,
    AtomType,
    Dataset,
    Isotope,
    Molecule,
    PdbSetup,
    Schedule,
    SystemSpec,
)

__all__ = [
    "SystemTemplate",
    "KNOWN_SYSTEMS",
    "ADJUSTABLE",
    "MolFile",
    "describe_library",
    "make_spec",
    "auto_spec",
    "read_mol",
    "repository_root",
]


@dataclass(frozen=True)
class SystemTemplate:
    """A known system and the function that builds its specification."""

    key: str  # e.g. "water"
    title: str
    description: str  # what the model is, where the parameters come from
    keywords: tuple[str, ...]  # words in a data title that suggest this system
    atom_types: tuple[str, ...]
    build: Callable[
        ..., SystemSpec
    ]  # build(datasets, name=None, **adjust) -> SystemSpec
    needs_pdb: bool = False


#: The ``adjust`` keywords :func:`make_spec` accepts.
ADJUSTABLE = (
    "density",
    "temperature",
    "ereq",
    "counts",
    "count",
    "title",
    "equilibration_iterations",
    "refinement_runs",
    "refinement_iterations",
)

#: Composition of "null" water (zero mean coherent scattering length) as written in
#: ``EPSR-test/DMSODWaterNull.NWTS.dat``: ``0 0.64 2 0.36``.
NULL_WATER = (Isotope(0, 0.64), Isotope(2, 0.36))

_NAME_CHARS = re.compile(r"[^A-Za-z0-9_-]")


def repository_root() -> Path:
    """The EPSR_hackathon repository (parent of ``mellor/``)."""
    return Path(__file__).resolve().parents[3]


# --------------------------------------------------------------------------- the physics
@dataclass(frozen=True)
class _System:
    """Everything a template needs; the public :class:`SystemTemplate` wraps it."""

    key: str
    title: str
    description: str
    origin: str  # the files the numbers were copied from
    keywords: tuple[str, ...]
    kind: str
    default_name: str
    atom_types: tuple[AtomType, ...]
    molecules: tuple[Molecule, ...]  # with their counts in the box
    density: float
    temperature: float
    ereq: float
    hydrogen_types: tuple[
        str, ...
    ] = ()  # types whose isotopes follow the sample's H/D labelling
    exchangeable: tuple[str, ...] = ()  # types with ``iexchange 1`` in every sample
    fixed_isotopes: dict[str, tuple[Isotope, ...]] = field(
        default_factory=dict
    )  # the same in every sample
    epsr_keywords: dict[str, Any] = field(default_factory=dict)
    rshmin: float | None = None
    pdb: PdbSetup | None = None


def _single_atom(label: str, count: int) -> Molecule:
    return Molecule(name=label, atoms=[Atom(label, 0.0, 0.0, 0.0)], count=count)


# ---- water: validation/water/Water.mol, mixato.stdin, *.NWTS.dat, waterbox.EPSR.inp, steps.txt
_WATER = _System(
    key="water",
    title="Liquid water (1000 molecules)",
    description=(
        "Rigid-geometry water (O-H 0.976 Å, H-O-H 107.134°) with the Dissolve water "
        "example force field (OW ε 0.6503 kJ/mol, σ 3.165492 Å, q -0.82 e; HW q +0.41 e, "
        "mass 2 so H or D), 1000 molecules at 0.1 atoms/Å³ and 300 K, ereq 3 kJ/mol. "
        "Fits H2O, D2O and HDO neutron data; HW follows the sample labelling and is "
        "exchangeable. Validated against the Fortran EPSR in validation/water."
    ),
    origin="validation/water (Water.mol, mixato.stdin, H2O/D2O/HDO.NWTS.dat, waterbox.EPSR.inp, steps.txt)",
    keywords=("water", "h2o", "d2o", "hdo"),
    kind="molecular",
    default_name="waterbox",
    atom_types=(
        AtomType(
            "OW",
            "O",
            epsilon=0.6503,
            sigma=3.165492,
            mass=15.9994,
            charge=-0.82,
            rmax=1.55,
        ),
        AtomType("HW", "H", epsilon=0.0, sigma=0.0, mass=2.0, charge=0.41, rmax=0.0),
    ),
    molecules=(
        Molecule(
            name="Water",
            atoms=[
                Atom("OW", 0.0, 0.0, 0.0),
                Atom("HW", 0.976000011, 0.0, 0.0),
                Atom("HW", -0.287536860, 0.932683527, 0.0),
            ],
            count=1000,
            bonds=[(1, 2), (1, 3)],
            bond_lengths={"1-2": 0.976, "1-3": 0.976},
            angles={"2-1-3": 107.134},
        ),
    ),
    density=0.1,
    temperature=300.0,
    ereq=3.0,
    hydrogen_types=("HW",),
    exchangeable=("HW",),
    epsr_keywords={"qmin": [0.05, 0.0]},
    rshmin=0.7,
)

# ---- silica: validation/silica_bulk/Si.mol, O.mol, mixato.stdin, SiO2.NWTS.dat, silica.EPSR.inp, steps.txt
_SILICA = _System(
    key="silica",
    title="Amorphous silica (1500 Si + 3000 O)",
    description=(
        "Atomic SiO2 glass: Si (ε 0.175 kJ/mol, σ 1.03 Å, +2 e) and O (ε 0.165, σ 3.5, "
        "-1 e) from the Dissolve silica example, 1500 Si + 3000 O at 0.0664 atoms/Å³ "
        "and 298 K, ereq 50 kJ/mol, qmin 0.1 Å⁻¹, rshmin 1.2 Å. Natural isotopes, one "
        "NIMROD data set. Validated against the Fortran EPSR in validation/silica_bulk."
    ),
    origin="validation/silica_bulk (Si.mol, O.mol, mixato.stdin, SiO2.NWTS.dat, silica.EPSR.inp, steps.txt)",
    keywords=("silica", "sio2", "glass"),
    kind="atomic",
    default_name="silica",
    atom_types=(
        AtomType(
            "Si", "Si", epsilon=0.175, sigma=1.03, mass=28.085501, charge=2.0, rmax=1.11
        ),
        AtomType(
            "O", "O", epsilon=0.165, sigma=3.5, mass=15.9994, charge=-1.0, rmax=1.55
        ),
    ),
    molecules=(_single_atom("Si", 1500), _single_atom("O", 3000)),
    density=0.0664,
    temperature=298.0,
    ereq=50.0,
    epsr_keywords={"qmin": [0.1, 0.0]},
    rshmin=1.2,
)

# ---- DMSO/water: EPSR-test/DMSO.mol, WaterMolecule.mol, 200DMSO400Water.{ato,EPSR.inp,pcof}, DMSODWater*.NWTS.dat
_DMSO = Molecule(
    name="DMSO",
    atoms=[
        Atom("S", 0.0, 0.0, 0.0),
        Atom("C", 1.81980002, 0.0, 0.0),
        Atom("C", -0.279000014, 1.79830003, 0.0),
        Atom("O", -0.383199990, -0.447299987, 1.44169998),
        Atom("H", 2.22189999, -1.01440001, -0.105300002),
        Atom("H", 2.22189999, 0.404799998, 0.936100006),
        Atom("H", 2.22189999, 0.597500026, -0.826499999),
        Atom("H", 0.249799997, 2.28730011, -0.826499999),
        Atom("H", 5.93999997e-02, 2.25769997, 0.936100006),
        Atom("H", -1.34300005, 2.04020000, -0.105400003),
    ],
    count=200,
    bonds=[(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (2, 7), (3, 8), (3, 9), (3, 10)],
    bond_lengths={
        "1-2": 1.81980,
        "1-3": 1.81981,
        "1-4": 1.55738,
        "2-5": 1.09626,
        "2-6": 1.09628,
        "2-7": 1.09626,
        "3-8": 1.09629,
        "3-9": 1.09629,
        "3-10": 1.09623,
    },
    angles={
        "2-1-3": 98.81892,
        "2-1-4": 104.24419,
        "3-1-4": 104.24656,
        "1-2-5": 111.51799,
        "1-2-6": 111.51751,
        "1-2-7": 111.51787,
        "5-2-6": 106.80777,
        "5-2-7": 107.30040,
        "6-2-7": 107.93754,
        "1-3-8": 111.51997,
        "1-3-9": 111.51669,
        "1-3-10": 111.52218,
        "8-3-9": 107.93495,
        "8-3-10": 107.29300,
        "9-3-10": 106.81195,
    },
    rotations=[(1, 2), (1, 3)],
)
_DMSO_WATER_MOLECULE = Molecule(
    name="WaterMolecule",
    atoms=[
        Atom("O1", 0.0, 0.0, 0.0),
        Atom("H1", 0.950999975, 0.0, 0.0),
        Atom("H1", -0.289000005, 0.906000018, 0.0),
    ],
    count=400,
    bonds=[(1, 2), (1, 3)],
    bond_lengths={"1-2": 0.976, "1-3": 0.976},
    angles={"2-1-3": 104.5},
)
_DMSO_WATER = _System(
    key="dmso_water",
    title="DMSO in water, 1:2 (200 DMSO + 400 water)",
    description=(
        "The EPSR test case: 200 DMSO (S, C, O, H with the EPSR-test DMSO.mol geometry, "
        "rotations about both S-C bonds) and 400 water (O1, H1) at 0.08807 atoms/Å³ and "
        "300 K, ereq 10 kJ/mol. DMSO is fully deuterated in every sample; the water "
        "hydrogen H1 follows the sample labelling (D2O, H2O, null water 0.64 H/0.36 D) "
        "and is exchangeable. SANDALS data D6DMSO_D2O/H2O/HDO."
    ),
    origin="EPSR-test (DMSO.mol, WaterMolecule.mol, 200DMSO400Water.EPSR.inp, 200DMSO400Water.pcof, DMSODWaterD/H/Null.NWTS.dat)",
    keywords=("dmso", "dimethyl sulfoxide"),
    kind="molecular",
    default_name="dmso_water",
    atom_types=(  # in the order of the .ato file (first appearance), as the nwts components
        AtomType(
            "S",
            "S",
            epsilon=0.845170,
            sigma=3.56,
            mass=32.066002,
            charge=0.54,
            rmax=1.6,
        ),
        AtomType(
            "C", "C", epsilon=0.7605, sigma=3.88, mass=12.0111, charge=0.0, rmax=1.85
        ),
        AtomType(
            "O", "O", epsilon=0.27614, sigma=2.94, mass=15.9994, charge=-0.54, rmax=1.55
        ),
        AtomType("H", "H", epsilon=0.0, sigma=0.0, mass=2.0, charge=0.0, rmax=0.0),
        AtomType(
            "O1",
            "O",
            epsilon=0.65,
            sigma=3.166,
            mass=15.9994,
            charge=-0.8476,
            rmax=1.55,
        ),
        AtomType("H1", "H", epsilon=0.0, sigma=0.0, mass=2.0, charge=0.4238, rmax=0.0),
    ),
    molecules=(_DMSO, _DMSO_WATER_MOLECULE),
    density=8.80699679e-02,  # rho of 200DMSO400Water.EPSR.inp (box 33.1213 Å, 3200 atoms)
    temperature=300.0,
    ereq=10.0,  # ereq 0.10000E+02 of 200DMSO400Water.pcof
    hydrogen_types=("H1",),
    exchangeable=("H1",),
    fixed_isotopes={"H": (Isotope(2, 1.0),)},
    epsr_keywords={"qmin": [0.05, 0.0]},
    rshmin=0.7,
)

# ---- benzene: Dissolve examples/benzene/input.txt (not in this repository)
_AVOGADRO = 6.02214076e23
_BENZENE_MOLAR_MASS = 6 * 12.011 + 6 * 1.008  # g/mol with the real masses
#: 0.876 g/cm³ -> atoms/Å³ (12 atoms per molecule): about 0.0811.
_BENZENE_DENSITY = 0.876 / _BENZENE_MOLAR_MASS * _AVOGADRO * 1e-24 * 12
_BENZENE = _System(
    key="benzene",
    title="Liquid benzene (200 molecules)",
    description=(
        "Benzene from the Dissolve benzene example: planar C6H6 (CA ε 0.29288 kJ/mol, "
        "σ 3.55 Å, -0.115 e; HA ε 0.12552, σ 2.42, +0.115 e, mass 2 so H or D), 200 "
        "molecules at 0.876 g/cm³ (0.0811 atoms/Å³) and 300 K, ereq 20 kJ/mol. HA "
        "follows the sample labelling (C6H6, C6D6, 50:50); nothing is exchangeable."
    ),
    origin="Dissolve examples/benzene/input.txt (github.com/disorderedmaterials/dissolve)",
    keywords=("benzene", "c6h6", "c6d6"),
    kind="molecular",
    default_name="benzene",
    atom_types=(
        AtomType("CA", "C", epsilon=0.29288, sigma=3.55, mass=12.011, charge=-0.115),
        AtomType("HA", "H", epsilon=0.12552, sigma=2.42, mass=2.0, charge=0.115),
    ),
    molecules=(
        Molecule(
            name="Benzene",
            atoms=[
                Atom("CA", -1.399, 0.160, 0.0),
                Atom("CA", -0.561, 1.293, 0.0),
                Atom("CA", 0.839, 1.132, 0.0),
                Atom("CA", 1.399, -0.160, 0.0),
                Atom("CA", 0.560, -1.293, 0.0),
                Atom("CA", -0.839, -1.132, 0.0),
                Atom("HA", 1.483, 2.001, 0.0),
                Atom("HA", 2.472, -0.284, 0.0),
                Atom("HA", 0.991, -2.284, 0.0),
                Atom("HA", -1.483, -2.000, 0.0),
                Atom("HA", -2.472, 0.282, 0.0),
                Atom("HA", -0.990, 2.284, 0.0),
            ],
            count=200,
            bonds=[
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 6),
                (6, 1),
                (7, 3),
                (4, 8),
                (5, 9),
                (6, 10),
                (1, 11),
                (2, 12),
            ],
        ),
    ),
    density=_BENZENE_DENSITY,
    temperature=300.0,
    ereq=20.0,
    hydrogen_types=("HA",),
)

# ---- Cu-BTC: validation/cu-btc (CuBTC.pdb, readpdb.stdin, changeato.stdin, CuBTC.NWTS.dat, CuBTC.EPSR.inp, README.md)
_CU_BTC_CELL = 26.3336  # Å, CRYST1 of CuBTC.pdb
_CU_BTC_ATOMS = 624  # 48 Cu, 192 O, 288 C, 96 H
_CU_BTC = _System(
    key="cu_btc",
    title="Cu-BTC (HKUST-1) metal-organic framework, one unit cell",
    description=(
        "Crystalline Cu-BTC read from validation/cu-btc/CuBTC.pdb (one cubic cell, "
        "a 26.3336 Å, 624 atoms: 24 Cu2 paddlewheels + 32 BTC linkers) with readpdb "
        "(van der Waals factor 0.68, no dihedrals) and the Dissolve Cu-BTC force field "
        "applied with changeato (CU, O, C1, C2, C3, H); 300 K, ereq 10 kJ/mol, Bragg "
        "peaks to 5 Å⁻¹ (qwidthqmax 0.2 5.0, hklqmin 0.0 10.0 0.0), rshmin 0.7 Å. "
        "Natural isotopes, one NIMROD data set of the empty framework. The crystal "
        "positions are kept; Cu-O coordination bonds are intermolecular (see README)."
    ),
    origin="validation/cu-btc (CuBTC.pdb, readpdb.stdin, changeato.stdin, CuBTC.NWTS.dat, CuBTC.EPSR.inp, README.md)",
    keywords=("cubtc", "cu-btc", "hkust-1", "mof"),
    kind="pdb",
    default_name="CuBTC",
    atom_types=(  # changeato.stdin atoms 1..6; masses and sphere radii as EPSR wrote them in CuBTC.ato
        AtomType(
            "CU",
            "Cu",
            epsilon=0.02092,
            sigma=3.144,
            mass=63.546,
            charge=1.098,
            rmax=1.6,
        ),
        AtomType(
            "O",
            "O",
            epsilon=0.401664,
            sigma=3.033,
            mass=15.9994,
            charge=-0.665,
            rmax=1.55,
        ),
        AtomType(
            "C1",
            "C",
            epsilon=0.39748,
            sigma=3.473,
            mass=12.0111,
            charge=0.778,
            rmax=1.85,
        ),
        AtomType(
            "C2",
            "C",
            epsilon=0.39748,
            sigma=3.473,
            mass=12.0111,
            charge=-0.092,
            rmax=1.85,
        ),
        AtomType(
            "C3",
            "C",
            epsilon=0.39748,
            sigma=3.473,
            mass=12.0111,
            charge=-0.014,
            rmax=1.85,
        ),
        AtomType("H", "H", epsilon=0.06276, sigma=2.846, mass=1.008, charge=0.109),
    ),
    molecules=(),
    density=_CU_BTC_ATOMS
    / _CU_BTC_CELL**3,  # informational: readpdb takes the box from the PDB
    temperature=300.0,
    ereq=10.0,
    epsr_keywords={
        "qwidthqmax": [0.2, 5.0],
        "hklqmin": [0.0, 10.0, 0.0],
        "qmin": [0.05, 0.0],
    },
    rshmin=0.7,
    pdb=PdbSetup(
        file=str(repository_root() / "validation" / "cu-btc" / "CuBTC.pdb"),
        vdw_factor=0.68,
        dihedrals=False,
        tolerance=0.0,
        step_sizes={"stepmi": 0.1, "stepri": 0.3, "stepr": 0.3, "stepm": 1.0},
    ),
)

_SYSTEMS: tuple[_System, ...] = (_WATER, _SILICA, _DMSO_WATER, _BENZENE, _CU_BTC)


# --------------------------------------------------------------------------- building specs
def _safe_label(stem: str) -> str:
    label = _NAME_CHARS.sub("_", stem).strip("_") or "data"
    if not label[0].isalpha():
        label = "d" + label
    return label[:24]


def _unique(label: str, used: set[str]) -> str:
    candidate = label
    number = 2
    while candidate in used:
        suffix = f"_{number}"
        candidate = label[: 24 - len(suffix)] + suffix
        number += 1
    used.add(candidate)
    return candidate


def _hydrogen_composition(guess: dict[str, Any]) -> tuple[list[Isotope], str]:
    """The isotopes of a hydrogen type that follows the sample, and a description."""
    hydrogen = guess.get("hydrogen")
    if hydrogen is None:
        return [Isotope(0, 1.0)], "no H/D labelling recognised, treated as natural"
    if hydrogen == "H":
        return [Isotope(0, 1.0)], "natural hydrogen"
    if hydrogen == "D":
        return [Isotope(2, 1.0)], "deuterium"
    if guess.get("sample") == "Null":
        return list(NULL_WATER), "null water (0.64 H / 0.36 D)"
    fraction = float(guess.get("hd_fraction") or 0.5)
    composition = [Isotope(1, round(1.0 - fraction, 6)), Isotope(2, round(fraction, 6))]
    return composition, f"H/D mixture, {fraction:.0%} D"


def _describe_composition(label: str, isotopes: list[Isotope]) -> str:
    pairs = " ".join(
        f"{isotope.mass_number} {isotope.fraction:g}" for isotope in isotopes
    )
    return f"{label} [{pairs}]"


def _dataset(
    system: _System, data: Mint01, used_labels: set[str]
) -> tuple[Dataset, str]:
    """The :class:`Dataset` of one file under ``system`` and the note explaining it."""
    guess = guess_sample(data.title, data.name)
    candidates = guess["candidates"]
    if candidates and system.key not in candidates:
        raise ValueError(
            f"{data.name}: the title {data.title!r} looks like {' or '.join(candidates)}, "
            f"not {system.key}; use the '{candidates[0]}' system or check the data"
        )
    isotopes: dict[str, list[Isotope]] = {
        label: list(fixed) for label, fixed in system.fixed_isotopes.items()
    }
    parts: list[str] = []
    if system.hydrogen_types:
        composition, text = _hydrogen_composition(guess)
        for label in system.hydrogen_types:
            isotopes[label] = list(composition)
        parts.append(f"{'/'.join(system.hydrogen_types)}: {text}")
    elif guess["hydrogen"] is not None:
        parts.append(
            f"the {guess['sample']} labelling is ignored: {system.key} has no sample-dependent hydrogen"
        )
    else:
        parts.append("natural isotopes")
    for label, fixed in system.fixed_isotopes.items():
        parts.append(f"{label} always {_describe_composition('', list(fixed)).strip()}")
    if system.exchangeable:
        parts.append(f"exchangeable: {', '.join(system.exchangeable)}")

    label = _unique(guess["sample"] or _safe_label(data.path.stem), used_labels)
    dataset = Dataset(
        file=str(Path(data.path).resolve()),
        label=label,
        isotopes=isotopes,
        exchangeable=list(system.exchangeable),
        rshmin=system.rshmin,
    )
    note = f"{data.name} ({data.title!r}) -> data set {label}: " + "; ".join(parts)
    return dataset, note


def _build(
    system: _System, datasets: list[Mint01], name: str | None = None, **adjust: Any
) -> SystemSpec:
    unknown = sorted(set(adjust) - set(ADJUSTABLE))
    if unknown:
        raise ValueError(
            f"unknown adjust keyword(s) {unknown}; allowed: {', '.join(ADJUSTABLE)}"
        )
    if not datasets:
        raise ValueError(f"{system.key}: at least one data set is needed")

    molecules = copy.deepcopy(list(system.molecules))
    counts = dict(adjust.get("counts") or {})
    if "count" in adjust:
        if not molecules:
            raise ValueError(
                f"{system.key}: 'count' cannot be adjusted, the system has no molecule templates"
            )
        counts.setdefault(molecules[0].name, adjust["count"])
    names = [molecule.name for molecule in molecules]
    for molecule_name, count in counts.items():
        if molecule_name not in names:
            raise ValueError(
                f"{system.key}: no molecule {molecule_name!r} to count; molecules: {names}"
            )
        molecules[names.index(molecule_name)].count = int(count)

    used_labels: set[str] = set()
    spec_datasets: list[Dataset] = []
    notes = [
        f"Built from the Mellor library system '{system.key}' ({system.title}); parameters copied from {system.origin}."
    ]
    for data in datasets:
        dataset, note = _dataset(system, data, used_labels)
        spec_datasets.append(dataset)
        notes.append(note)
    for keyword in (
        "density",
        "temperature",
        "ereq",
        "equilibration_iterations",
        "refinement_runs",
        "refinement_iterations",
        "title",
    ):
        if keyword in adjust:
            notes.append(f"adjusted {keyword} = {adjust[keyword]!r}")
    if counts:
        notes.append(f"adjusted molecule counts: {counts}")

    schedule = Schedule(
        ereq=float(adjust.get("ereq", system.ereq)),
        equilibration_iterations=int(
            adjust.get("equilibration_iterations", Schedule.equilibration_iterations)
        ),
        refinement_runs=int(adjust.get("refinement_runs", Schedule.refinement_runs)),
        refinement_iterations=int(
            adjust.get("refinement_iterations", Schedule.refinement_iterations)
        ),
    )
    return SystemSpec(
        name=name or system.default_name,
        title=str(adjust.get("title", system.title)),
        kind=system.kind,
        atom_types=list(system.atom_types),
        molecules=molecules,
        density=float(adjust.get("density", system.density)),
        temperature=float(adjust.get("temperature", system.temperature)),
        datasets=spec_datasets,
        schedule=schedule,
        epsr_keywords=copy.deepcopy(system.epsr_keywords),
        pdb=copy.deepcopy(system.pdb),
        notes="\n".join(notes),
        source=f"library:{system.key}",
    )


def _template(system: _System) -> SystemTemplate:
    return SystemTemplate(
        key=system.key,
        title=system.title,
        description=system.description,
        keywords=system.keywords,
        atom_types=tuple(atom_type.label for atom_type in system.atom_types),
        build=functools.partial(_build, system),
        needs_pdb=system.pdb is not None,
    )


#: key -> template. Keys: water, silica, dmso_water, benzene, cu_btc.
KNOWN_SYSTEMS: dict[str, SystemTemplate] = {
    system.key: _template(system) for system in _SYSTEMS
}
_PHYSICS: dict[str, _System] = {system.key: system for system in _SYSTEMS}


def _lookup(key: str) -> _System:
    if key not in _PHYSICS:
        raise LookupError(
            f"unknown system {key!r}; known systems: {', '.join(KNOWN_SYSTEMS)}"
        )
    return _PHYSICS[key]


def describe_library() -> list[dict[str, Any]]:
    """JSON-serialisable listing for the model: key, title, description, atom types, keywords."""
    listing = []
    for system in _SYSTEMS:
        listing.append(
            {
                "key": system.key,
                "title": system.title,
                "description": system.description,
                "kind": system.kind,
                "atom_types": [atom_type.label for atom_type in system.atom_types],
                "keywords": list(system.keywords),
                "molecules": [
                    {
                        "name": molecule.name,
                        "count": molecule.count,
                        "atoms": molecule.n_atoms,
                    }
                    for molecule in system.molecules
                ],
                "density": round(system.density, 6),
                "temperature": system.temperature,
                "ereq": system.ereq,
                "hydrogen_types": list(system.hydrogen_types),
                "needs_pdb": system.pdb is not None,
                "default_name": system.default_name,
            }
        )
    return listing


def make_spec(
    key: str, datasets: list[Mint01], *, name: str | None = None, **adjust: Any
) -> SystemSpec:
    """The specification of a known system for the given data files.

    ``adjust`` may override top-level numbers (``density``, ``temperature``,
    ``ereq``, ``counts`` as a dict of molecule name -> count, ``count`` of the first
    molecule, ``title``, ``equilibration_iterations``, ``refinement_runs``,
    ``refinement_iterations``; see :data:`ADJUSTABLE`). Every data set's isotopes
    come from ``data.guess_sample`` on its title; a data set with no recognisable
    hydrogen labelling is treated as natural (said in ``spec.notes``). The data set
    labels are the recognised sample names (``H2O``, ``D2O``, ``HDO``, ``Null``...)
    or the file stems. ``Dataset.file`` is the absolute path of each file.

    Raises ``LookupError`` for an unknown key and ``ValueError`` when a title names
    a different system (a water file with the silica template) or ``adjust`` is
    wrong. The result is not validated: call ``spec.validate()``.
    """
    _lookup(key)
    return KNOWN_SYSTEMS[key].build(datasets, name=name, **adjust)


def auto_spec(
    datasets: list[Mint01], *, name: str | None = None
) -> tuple[SystemSpec, str]:
    """The best library match for the data titles, and the reason as text.

    A system must be a candidate (``data.guess_sample``) of every data set whose
    title names any system; data sets that name nothing are accepted as natural.
    Raises ``LookupError`` when no template matches.
    """
    if not datasets:
        raise LookupError("auto_spec: no data sets given")
    guesses = [guess_sample(data.title, data.name) for data in datasets]
    named = [guess["candidates"] for guess in guesses if guess["candidates"]]
    summary = "; ".join(
        f"{data.name}: {data.title!r} -> {', '.join(guess['candidates']) or 'no known system'}"
        for data, guess in zip(datasets, guesses, strict=True)
    )
    if not named:
        raise LookupError(
            f"no built-in system matches the data titles ({summary}); known systems: {', '.join(KNOWN_SYSTEMS)}"
        )
    common = [
        key for key in KNOWN_SYSTEMS if all(key in candidates for candidates in named)
    ]
    if not common:
        raise LookupError(
            f"the data titles do not all point to one built-in system ({summary})"
        )
    best = min(
        common, key=lambda key: sum(candidates.index(key) for candidates in named)
    )
    spec = make_spec(best, datasets, name=name)
    reason = (
        f"Chose the library system '{best}' ({KNOWN_SYSTEMS[best].title}): it matches "
        f"{len(named)} of {len(datasets)} data title(s) ({summary})."
    )
    return spec, reason


# --------------------------------------------------------------------------- .mol reading
@dataclass
class MolFile:
    """A parsed ``.mol`` (gmol) template: its molecule, atom types and temperature."""

    name: str
    molecule: Molecule
    atom_types: list[AtomType]
    temperature: float


def read_mol(path: str | Path) -> MolFile:
    """Parse an EPSR ``.mol`` file (``makemole`` input) into spec objects.

    Reads ``atom`` (index, label, x, y, z, neighbours), ``bond``, ``angle``, ``rot``,
    ``changelabel``, ``potential`` (label, epsilon, sigma, mass, charge, element,
    nsph, rmin, rmax, cg) and ``temperature`` lines. Used by the tests to check the
    library against ``validation/`` and ``EPSR-test``; handy for user templates.
    """
    path = Path(path)
    atoms: list[tuple[str, float, float, float, list[int]]] = []
    labels: dict[str, str] = {}
    bonds: list[tuple[int, int]] = []
    bond_lengths: dict[str, float] = {}
    angles: dict[str, float] = {}
    rotations: list[tuple[int, int]] = []
    atom_types: list[AtomType] = []
    temperature = 300.0
    for raw in path.read_text().splitlines():
        words = raw.split()
        if not words:
            continue
        keyword = words[0].lower()
        try:
            if keyword == "atom":
                n_neighbours = int(words[6]) if len(words) > 6 else 0
                neighbours = [int(word) for word in words[7 : 7 + n_neighbours]]
                atoms.append(
                    (
                        words[2],
                        float(words[3]),
                        float(words[4]),
                        float(words[5]),
                        neighbours,
                    )
                )
            elif keyword == "bond":
                first, second = int(words[1]), int(words[2])
                bonds.append((first, second))
                if len(words) > 3:
                    bond_lengths[f"{first}-{second}"] = float(words[3])
            elif keyword == "angle":
                angles[f"{int(words[1])}-{int(words[2])}-{int(words[3])}"] = float(
                    words[4]
                )
            elif keyword == "rot":
                rotations.append((int(words[1]), int(words[2])))
            elif keyword == "changelabel":
                labels[words[1]] = words[2]
            elif keyword == "potential":
                values = [float(word) for word in words[2:6]]
                extra = [float(word) for word in words[7:11]]
                atom_types.append(
                    AtomType(
                        label=words[1],
                        element=words[6],
                        epsilon=values[0],
                        sigma=values[1],
                        mass=values[2],
                        charge=values[3],
                        nsph=int(extra[0]) if len(extra) > 0 else 1,
                        rmin=extra[1] if len(extra) > 1 else 0.0,
                        rmax=extra[2] if len(extra) > 2 else None,
                        cg=extra[3] if len(extra) > 3 else 1.0,
                    )
                )
            elif keyword == "temperature":
                temperature = float(words[1])
        except (IndexError, ValueError) as error:
            raise ValueError(
                f"{path.name}: cannot parse line {raw!r}: {error}"
            ) from error
    if not atoms:
        raise ValueError(f"{path.name}: no atom lines")
    if not bonds:  # no explicit bond lines: take them from the neighbour lists
        for index, (_, _, _, _, neighbours) in enumerate(atoms, start=1):
            bonds.extend((index, other) for other in neighbours if index < other)
    molecule = Molecule(
        name=path.stem,
        atoms=[Atom(labels.get(label, label), x, y, z) for label, x, y, z, _ in atoms],
        count=1,
        bonds=bonds,
        bond_lengths=bond_lengths,
        angles=angles,
        rotations=rotations,
    )
    return MolFile(
        name=path.stem,
        molecule=molecule,
        atom_types=atom_types,
        temperature=temperature,
    )
