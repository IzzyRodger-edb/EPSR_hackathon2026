"""Test specifications written as plain dicts, with the numbers of ``validation/``.

The library module builds the same systems for the agent; the tests here do not
depend on it. Every number is copied from ``validation/water/Water.mol``,
``validation/silica_bulk/Si.mol`` and ``O.mol``, the ``.EPSR.inp`` files and
``validation/cu-btc/changeato.stdin``.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

from mellor.spec import SystemSpec, spec_from_dict

REPOSITORY = Path(__file__).resolve().parents[2]
EXPERIMENTAL_DATA = REPOSITORY / "experimental_data"
VALIDATION = REPOSITORY / "validation"
WATER_DATA = EXPERIMENTAL_DATA / "water" / "data"
SILICA_DATA = EXPERIMENTAL_DATA / "silica_bulk" / "data"
CU_BTC_DATA = EXPERIMENTAL_DATA / "cu-btc" / "data"

#: Water: 1000 molecules, Dissolve geometry and force field, three SANDALS samples in
#: the order of validation/water/waterbox.EPSR.inp (D2O, HDO, H2O).
WATER: dict[str, Any] = {
    "name": "waterbox",
    "title": "Water, 1000 H2O at 0.1 atoms/A^3 and 300 K: H2O, D2O and HDO (SANDALS)",
    "kind": "molecular",
    "atom_types": [
        {
            "label": "OW",
            "element": "O",
            "epsilon": 0.6503,
            "sigma": 3.165492,
            "mass": 15.9994,
            "charge": -0.82,
            "rmax": 1.55,
        },
        {
            "label": "HW",
            "element": "H",
            "epsilon": 0.0,
            "sigma": 0.0,
            "mass": 2.0,
            "charge": 0.41,
            "rmax": 0.0,
        },
    ],
    "molecules": [
        {
            "name": "Water",
            "count": 1000,
            "atoms": [
                ["OW", 0.0, 0.0, 0.0],
                ["HW", 0.976, 0.0, 0.0],
                ["HW", -0.28753686, 0.93268352, 0.0],
            ],
            "bonds": [[1, 2], [1, 3]],
            "bond_lengths": {"1-2": 0.976, "1-3": 0.976},
            "angles": {"2-1-3": 107.134},
        }
    ],
    "density": 0.1,
    "temperature": 300.0,
    "datasets": [
        {
            "file": str(WATER_DATA / "SLS18502-D2O.mint01"),
            "label": "D2O",
            "isotopes": {"HW": [[2, 1.0]]},
            "exchangeable": ["HW"],
        },
        {
            "file": str(WATER_DATA / "SLS18500-HDO5050.mint01"),
            "label": "HDO",
            "isotopes": {"HW": [[1, 0.5], [2, 0.5]]},
            "exchangeable": ["HW"],
        },
        {
            "file": str(WATER_DATA / "SLS18498-H2O.mint01"),
            "label": "H2O",
            "isotopes": {},
            "exchangeable": ["HW"],
        },
    ],
    "schedule": {"ereq": 3.0},
    "epsr_keywords": {"qmin": [0.05, 0.0]},
    "notes": "Test specification copied from validation/water.",
}

#: Amorphous silica: 1500 Si + 3000 O single-atom molecules, one NIMROD data set.
SILICA: dict[str, Any] = {
    "name": "silica",
    "title": "Amorphous silica, 1500 Si + 3000 O at 0.0664 atoms/A^3 and 298 K (NIMROD)",
    "kind": "atomic",
    "atom_types": [
        {
            "label": "Si",
            "element": "Si",
            "epsilon": 0.175,
            "sigma": 1.03,
            "mass": 28.0855,
            "charge": 2.0,
            "rmax": 1.11,
        },
        {
            "label": "O",
            "element": "O",
            "epsilon": 0.165,
            "sigma": 3.5,
            "mass": 15.9994,
            "charge": -1.0,
            "rmax": 1.55,
        },
    ],
    "molecules": [
        {"name": "Si", "count": 1500, "atoms": [["Si", 0.0, 0.0, 0.0]]},
        {"name": "O", "count": 3000, "atoms": [["O", 0.0, 0.0, 0.0]]},
    ],
    "density": 0.0664,
    "temperature": 298.0,
    "datasets": [
        {
            "file": str(SILICA_DATA / "SiO2_NIMROD.mint01"),
            "label": "SiO2",
            "isotopes": {},
            "rshmin": 1.2,
        }
    ],
    "schedule": {"ereq": 50.0},
    "epsr_keywords": {"qmin": [0.1, 0.0]},
    "notes": "Test specification copied from validation/silica_bulk.",
}

#: Cu-BTC (HKUST-1): one unit cell read from a PDB file, Dissolve force field.
CU_BTC: dict[str, Any] = {
    "name": "CuBTC",
    "title": "Cu-BTC (HKUST-1), one unit cell (624 atoms) at 300 K, empty framework (NIMROD)",
    "kind": "pdb",
    "atom_types": [
        {
            "label": "CU",
            "element": "Cu",
            "epsilon": 0.02092,
            "sigma": 3.144,
            "mass": 63.546,
            "charge": 1.098,
        },
        {
            "label": "O",
            "element": "O",
            "epsilon": 0.401664,
            "sigma": 3.033,
            "mass": 15.9994,
            "charge": -0.665,
        },
        {
            "label": "C1",
            "element": "C",
            "epsilon": 0.39748,
            "sigma": 3.473,
            "mass": 12.0111,
            "charge": 0.778,
        },
        {
            "label": "C2",
            "element": "C",
            "epsilon": 0.39748,
            "sigma": 3.473,
            "mass": 12.0111,
            "charge": -0.092,
        },
        {
            "label": "C3",
            "element": "C",
            "epsilon": 0.39748,
            "sigma": 3.473,
            "mass": 12.0111,
            "charge": -0.014,
        },
        {
            "label": "H",
            "element": "H",
            "epsilon": 0.06276,
            "sigma": 2.846,
            "mass": 1.008,
            "charge": 0.109,
        },
    ],
    "molecules": [],
    "density": 624 / 26.3336**3,
    "temperature": 300.0,
    "datasets": [
        {
            "file": str(CU_BTC_DATA / "Empty-CuBTC.mint01"),
            "label": "CuBTC",
            "isotopes": {},
            "rshmin": 0.7,
        }
    ],
    "schedule": {"ereq": 10.0},
    "epsr_keywords": {
        "qmin": [0.05, 0.0],
        "qwidthqmax": [0.2, 5.0],
        "hklqmin": [0.0, 10.0, 0.0],
    },
    "pdb": {"file": str(VALIDATION / "cu-btc" / "CuBTC.pdb"), "vdw_factor": 0.68},
    "notes": "Test specification copied from validation/cu-btc.",
}


def make(data: dict[str, Any], **changes: Any) -> SystemSpec:
    """A :class:`SystemSpec` from one of the dicts above, with top-level changes."""
    merged = copy.deepcopy(data)
    merged.update(changes)
    return spec_from_dict(merged)


def water_spec(**changes: Any) -> SystemSpec:
    return make(WATER, **changes)


def silica_spec(**changes: Any) -> SystemSpec:
    return make(SILICA, **changes)


def cu_btc_spec(**changes: Any) -> SystemSpec:
    return make(CU_BTC, **changes)
