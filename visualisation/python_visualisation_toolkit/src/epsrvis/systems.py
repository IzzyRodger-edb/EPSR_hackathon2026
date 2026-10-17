"""Per-system settings: atom appearance, bonds and the structural analyses to run.

Add an entry to ``CASES`` (keyed by the case folder name stored in
``meta.json``) to support a new system.  Unknown cases fall back to
``default_case``, which guesses elements from the atom-type labels.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# Conventional (CPK / Jmol) element colours and display radii (Å, ~0.4 x van der Waals).
ELEMENTS = {
    "H": dict(color="#f0f0f0", radius=0.30, mass=1.008),
    "D": dict(color="#f0f0f0", radius=0.30, mass=2.014),
    "O": dict(color="#e0302a", radius=0.62, mass=15.9994),
    "Si": dict(color="#d9b48f", radius=0.75, mass=28.0855),
    "C": dict(color="#606060", radius=0.65, mass=12.011),
    "N": dict(color="#3050f8", radius=0.62, mass=14.007),
    "S": dict(color="#e6c619", radius=0.80, mass=32.06),
}


@dataclass
class Case:
    title: str
    element: dict[str, str]  # atom type -> element
    # bonds drawn in the 3D view: (type a, type b, cutoff Å, intramolecular only?)
    bonds: list[tuple[str, str, float, bool]] = field(default_factory=list)
    # coordination-number distributions: (centre, neighbour, cutoff Å, intermolecular only?, description)
    coordination: list[tuple[str, str, float, bool, str]] = field(default_factory=list)
    # bond-angle distributions: (outer, centre, cutoff Å, mode, description); mode "intra" | "inter" | "all"
    angles: list[tuple[str, str, float, str, str]] = field(default_factory=list)
    stages: dict[int, str] = field(default_factory=dict)  # run number -> stage label

    def mass(self, atom_type: str) -> float:
        return ELEMENTS[self.element[atom_type]]["mass"]

    def color(self, atom_type: str) -> str:
        return ELEMENTS[self.element[atom_type]]["color"]

    def radius(self, atom_type: str) -> float:
        return ELEMENTS[self.element[atom_type]]["radius"]


CASES = {
    "water": Case(
        title="Water: 1000 H$_2$O, H$_2$O/D$_2$O/HDO (SANDALS)",
        element={"OW": "O", "HW": "H"},
        bonds=[("OW", "HW", 1.25, True)],
        coordination=[
            ("OW", "OW", 3.35, True, "O neighbours of O (r < 3.35 Å)"),
            ("OW", "HW", 2.45, True, "H-bonds per O (O···H < 2.45 Å)"),
        ],
        angles=[
            ("HW", "OW", 1.25, "intra", "H–O–H (intramolecular)"),
            ("OW", "OW", 3.35, "inter", "O···O···O (tetrahedrality)"),
        ],
        stages={1: "scratch", 2: "reference potential", 3: "EP refinement", 4: "EP refinement"},
    ),
    "silica_bulk": Case(
        title="Amorphous SiO$_2$: 1500 Si + 3000 O (NIMROD)",
        element={"Si": "Si", "O": "O"},
        bonds=[("Si", "O", 2.0, False)],
        coordination=[
            ("Si", "O", 2.0, False, "O around Si (r < 2.0 Å)"),
            ("O", "Si", 2.0, False, "Si around O (r < 2.0 Å)"),
        ],
        angles=[
            ("O", "Si", 2.0, "all", "O–Si–O"),
            ("Si", "O", 2.0, "all", "Si–O–Si"),
        ],
        stages={1: "scratch", 2: "reference potential", 3: "EP refinement", 4: "EP refinement"},
    ),
}


def default_case(atom_types: list[str]) -> Case:
    def guess(t: str) -> str:
        for n in (2, 1):
            e = t[:n].capitalize()
            if e in ELEMENTS:
                return e
        return "C"

    return Case(title="EPSR simulation", element={t: guess(t) for t in atom_types})


def get_case(name: str, atom_types: list[str]) -> Case:
    return CASES.get(name) or default_case(atom_types)
