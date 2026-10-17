"""Tests of ``mellor.library``: every template is checked against the files its numbers come from."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from mellor import library
from mellor.data import Mint01, read_mint01
from mellor.library import (
    KNOWN_SYSTEMS,
    auto_spec,
    describe_library,
    make_spec,
    read_mol,
)
from mellor.spec import Isotope, SpecError, SystemSpec

ROOT = library.repository_root()
VALIDATION = ROOT / "validation"
EPSR_TEST = ROOT / "EPSR-test"
WATER_FILES = sorted((ROOT / "experimental_data" / "water" / "data").glob("*.mint01"))
SILICA_FILE = ROOT / "experimental_data" / "silica_bulk" / "data" / "SiO2_NIMROD.mint01"
CU_BTC_FILE = ROOT / "experimental_data" / "cu-btc" / "data" / "Empty-CuBTC.mint01"
DMSO_FILES = sorted(EPSR_TEST.glob("D6DMSO_*.mint01"))

HEADER = (
    "# {source}\n# {title}\n#     1  3\n#    10\n# Tester\n# 01-JAN-2020 12:00:00\n"
)
ROWS = "  0.1  -0.2  0.01\n  0.2  -0.1  0.01\n  0.3  0.1  0.02\n  0.4  0.0  0.0\n"


def write_mint01(path: Path, title: str, source: str = "SLS00001.mint01") -> Mint01:
    path.write_text(HEADER.format(source=source, title=title) + ROWS)
    return read_mint01(path)


# --------------------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def water_data() -> list[Mint01]:
    return [read_mint01(path) for path in WATER_FILES]


@pytest.fixture(scope="module")
def silica_data() -> list[Mint01]:
    return [read_mint01(SILICA_FILE)]


@pytest.fixture(scope="module")
def cu_btc_data() -> list[Mint01]:
    return [read_mint01(CU_BTC_FILE)]


@pytest.fixture(scope="module")
def dmso_data() -> list[Mint01]:
    return [read_mint01(path) for path in DMSO_FILES]


@pytest.fixture
def benzene_data(tmp_path: Path) -> list[Mint01]:
    titles = {
        "C6H6": "C6H6 benzene in 1mm can",
        "C6D6": "C6D6 benzene in 1mm can",
        "5050": "50:50 C6H6/C6D6 benzene",
    }
    return [
        write_mint01(tmp_path / f"{stem}.mint01", title)
        for stem, title in titles.items()
    ]


# --------------------------------------------------------------------------- parsers of the reference files
def composition(dataset, label: str) -> list[tuple[int, float]]:
    return [
        (isotope.mass_number, isotope.fraction)
        for isotope in dataset.composition(label)
    ]


def numbers(tokens: list[str]) -> list[float]:
    """Leading numeric tokens (the value of an EPSR keyword line, before its description)."""
    values: list[float] = []
    for token in tokens:
        try:
            values.append(float(token))
        except ValueError:
            break
    return values


def nwts_components(path: Path) -> list[tuple[str, int, list[tuple[int, float]]]]:
    """``(label, iexchange, [(mass number, abundance), ...])`` per component of a ``.NWTS.dat``."""
    components: list[tuple[str, int, list[tuple[int, float]]]] = []
    label, exchange, abundances = "", 0, []
    for line in path.read_text().splitlines():
        words = line.split()
        if not words:
            continue
        if words[0] == "atom":
            label = words[1]
        elif words[0] == "iexchange":
            exchange = int(words[1])
        elif words[0] == "abundances":
            values = numbers(words[1:])
            abundances = [
                (int(values[i]), values[i + 1]) for i in range(0, len(values), 2)
            ]
            components.append((label, exchange, abundances))
    return components


def inp_values(path: Path) -> dict[str, list[list[float]]]:
    """keyword -> the numeric values of every line with that keyword, in file order."""
    values: dict[str, list[list[float]]] = {}
    for line in path.read_text().splitlines():
        words = line.split()
        if words:
            values.setdefault(words[0], []).append(numbers(words[1:]))
    return values


def inp_data_blocks(path: Path) -> list[tuple[str, str, float]]:
    """``(datafile, wtsfile, rshmin)`` of every ``data n`` block of a ``.EPSR.inp``."""
    blocks: list[tuple[str, str, float]] = []
    datafile = wtsfile = ""
    for line in path.read_text().splitlines():
        words = line.split()
        if not words:
            continue
        if words[0] == "datafile":
            datafile = words[1]
        elif words[0] == "wtsfile":
            wtsfile = words[1]
        elif words[0] == "rshmin":
            blocks.append((datafile, wtsfile, float(words[1])))
    return blocks


def mixato_answers(path: Path) -> tuple[list[tuple[str, int]], float, str]:
    lines = path.read_text().split()
    n_types = int(lines[0])
    types = [(lines[1 + 2 * i], int(lines[2 + 2 * i])) for i in range(n_types)]
    return types, float(lines[1 + 2 * n_types]), lines[2 + 2 * n_types]


def steps_ereq(path: Path) -> float:
    match = re.search(
        r"^set\s+\S+\.pcof\s+ereq\s+([\d.]+)", path.read_text(), re.MULTILINE
    )
    assert match, path
    return float(match.group(1))


def check_against_epsr_files(spec: SystemSpec, inp_path: Path, folder: Path) -> None:
    """The data blocks, weights files and ``qmin`` of a validated ``.EPSR.inp`` match ``spec``."""
    blocks = inp_data_blocks(inp_path)
    assert len(blocks) == len(spec.datasets)
    by_file = {dataset.file_name: dataset for dataset in spec.datasets}
    for datafile, wtsfile, rshmin in blocks:
        dataset = by_file[datafile]
        assert dataset.rshmin == rshmin
        assert dataset.nrtype == 5
        components = nwts_components(
            folder / wtsfile.replace(".NWTStot.wts", ".NWTS.dat")
        )
        assert [label for label, _, _ in components] == spec.labels
        for label, exchange, abundances in components:
            assert composition(dataset, label) == abundances, (datafile, label)
            assert (label in dataset.exchangeable) == (exchange == 1), (datafile, label)
    values = inp_values(inp_path)
    assert spec.epsr_keywords["qmin"] == values["qmin"][0]
    for keyword, value in spec.epsr_keywords.items():
        assert values[keyword][0] == value, keyword


def check_molecule_against_mol(
    spec: SystemSpec, name: str, mol_path: Path, count: int
) -> None:
    """The template molecule ``name`` and its atom types equal the parsed ``.mol`` file."""
    mol = read_mol(mol_path)
    molecule = next(molecule for molecule in spec.molecules if molecule.name == name)
    assert molecule.count == count
    assert molecule.labels == mol.molecule.labels
    assert [atom.xyz for atom in molecule.atoms] == [
        atom.xyz for atom in mol.molecule.atoms
    ]
    assert molecule.bonds == mol.molecule.bonds
    assert molecule.bond_lengths == mol.molecule.bond_lengths
    assert molecule.angles == mol.molecule.angles
    assert molecule.rotations == mol.molecule.rotations
    assert {(i, j, k) for i, j, k in molecule.angle_triples()} == {
        tuple(int(index) for index in key.split("-")) for key in mol.molecule.angles
    }
    for atom_type in mol.atom_types:
        assert spec.atom_type(atom_type.label) == atom_type, atom_type.label
    assert spec.temperature == mol.temperature


# --------------------------------------------------------------------------- the listing
def test_known_systems() -> None:
    assert list(KNOWN_SYSTEMS) == ["water", "silica", "dmso_water", "benzene", "cu_btc"]
    for key, template in KNOWN_SYSTEMS.items():
        assert template.key == key
        assert template.title and template.description and template.keywords
        assert template.atom_types and callable(template.build)
        assert template.needs_pdb == (key == "cu_btc")
    assert KNOWN_SYSTEMS["water"].atom_types == ("OW", "HW")
    assert KNOWN_SYSTEMS["cu_btc"].atom_types == ("CU", "O", "C1", "C2", "C3", "H")


def test_describe_library() -> None:
    listing = describe_library()
    json.dumps(listing)
    assert [item["key"] for item in listing] == list(KNOWN_SYSTEMS)
    for item in listing:
        assert {
            "key",
            "title",
            "description",
            "atom_types",
            "keywords",
            "molecules",
            "density",
            "temperature",
            "ereq",
            "kind",
        } <= set(item)
    water = listing[0]
    assert water["molecules"] == [{"name": "Water", "count": 1000, "atoms": 3}]
    assert (
        water["density"] == 0.1
        and water["temperature"] == 300.0
        and water["ereq"] == 3.0
    )
    assert listing[-1]["needs_pdb"] is True and listing[-1]["molecules"] == []


# --------------------------------------------------------------------------- water
def test_water_spec(water_data: list[Mint01]) -> None:
    spec = make_spec("water", water_data)
    assert spec.validate(check_files=True) == []
    assert spec.name == "waterbox" and spec.kind == "molecular"
    assert spec.source == "library:water"
    assert round(spec.box_side, 4) == 31.0723
    assert spec.n_atoms == 3000 and spec.n_molecules == 1000
    assert spec.labels == ["OW", "HW"]
    assert spec.density == 0.1 and spec.temperature == 300.0
    assert spec.schedule.ereq == 3.0
    assert spec.schedule.iterations() == [
        ("scratch", 1, 0.0),
        ("reference potential", 99, 0.0),
        ("EP refinement", 100, 3.0),
        ("EP refinement", 100, 3.0),
    ]
    assert spec.epsr_keywords == {"qmin": [0.05, 0.0]}
    assert [dataset.label for dataset in spec.datasets] == ["H2O", "HDO", "D2O"]
    h2o, hdo, d2o = spec.datasets
    assert composition(h2o, "HW") == [(0, 1.0)]
    assert composition(hdo, "HW") == [(1, 0.5), (2, 0.5)]
    assert composition(d2o, "HW") == [(2, 1.0)]
    for dataset, data in zip(spec.datasets, water_data, strict=True):
        assert (
            dataset.file == str(data.path.resolve())
            and Path(dataset.file).is_absolute()
        )
        assert dataset.exchangeable == ["HW"]
        assert dataset.rshmin == 0.7
        assert composition(dataset, "OW") == [(0, 1.0)]
        assert data.name in spec.notes
    assert spec.mass_density() == pytest.approx(0.1 * (15.9994 + 4.0) / 3 * 1.66053907)


def test_water_matches_validation_files(water_data: list[Mint01]) -> None:
    spec = make_spec("water", water_data)
    check_molecule_against_mol(spec, "Water", VALIDATION / "water" / "Water.mol", 1000)
    types, density, box = mixato_answers(VALIDATION / "water" / "mixato.stdin")
    assert [(molecule.name, molecule.count) for molecule in spec.molecules] == types
    assert spec.density == density and spec.name == box
    assert spec.schedule.ereq == steps_ereq(VALIDATION / "water" / "steps.txt")
    check_against_epsr_files(
        spec, VALIDATION / "water" / "waterbox.EPSR.inp", VALIDATION / "water"
    )


# --------------------------------------------------------------------------- silica
def test_silica_spec(silica_data: list[Mint01]) -> None:
    spec = make_spec("silica", silica_data)
    assert spec.validate(check_files=True) == []
    assert spec.kind == "atomic" and spec.name == "silica"
    assert round(spec.box_side, 4) == 40.7707
    assert spec.n_atoms == 4500 and spec.type_counts() == {"Si": 1500, "O": 3000}
    assert spec.temperature == 298.0 and spec.schedule.ereq == 50.0
    assert spec.epsr_keywords == {"qmin": [0.1, 0.0]}
    (dataset,) = spec.datasets
    assert dataset.label == "SiO2_NIMROD" and dataset.rshmin == 1.2
    assert dataset.isotopes == {} and dataset.exchangeable == []
    assert "natural" in spec.notes


def test_silica_matches_validation_files(silica_data: list[Mint01]) -> None:
    spec = make_spec("silica", silica_data)
    check_molecule_against_mol(spec, "Si", VALIDATION / "silica_bulk" / "Si.mol", 1500)
    check_molecule_against_mol(spec, "O", VALIDATION / "silica_bulk" / "O.mol", 3000)
    types, density, box = mixato_answers(VALIDATION / "silica_bulk" / "mixato.stdin")
    assert [(molecule.name, molecule.count) for molecule in spec.molecules] == types
    assert spec.density == density and spec.name == box
    assert spec.schedule.ereq == steps_ereq(VALIDATION / "silica_bulk" / "steps.txt")
    check_against_epsr_files(
        spec, VALIDATION / "silica_bulk" / "silica.EPSR.inp", VALIDATION / "silica_bulk"
    )


# --------------------------------------------------------------------------- DMSO / water
def test_dmso_water_spec(dmso_data: list[Mint01]) -> None:
    spec = make_spec("dmso_water", dmso_data)
    assert spec.validate(check_files=True) == []
    assert spec.name == "dmso_water" and spec.kind == "molecular"
    assert spec.labels == ["S", "C", "O", "H", "O1", "H1"]
    assert [(molecule.name, molecule.count) for molecule in spec.molecules] == [
        ("DMSO", 200),
        ("WaterMolecule", 400),
    ]
    assert spec.n_atoms == 3200
    assert (
        round(spec.box_side, 4) == 33.1213
    )  # box side of EPSR-test/200DMSO400Water.ato
    assert spec.temperature == 300.0 and spec.schedule.ereq == 10.0
    assert [dataset.label for dataset in spec.datasets] == ["D2O", "H2O", "Null"]
    d2o, h2o, null = spec.datasets
    for dataset in spec.datasets:
        assert composition(dataset, "H") == [
            (2, 1.0)
        ]  # DMSO fully deuterated in every sample
        assert dataset.exchangeable == ["H1"] and dataset.rshmin == 0.7
    assert composition(d2o, "H1") == [(2, 1.0)]
    assert composition(h2o, "H1") == [(0, 1.0)]
    assert composition(null, "H1") == [(0, 0.64), (2, 0.36)]
    assert spec.molecules[0].rotations == [(1, 2), (1, 3)]
    assert len(spec.molecules[0].angles) == 15 == len(spec.molecules[0].angle_triples())


def test_dmso_water_matches_epsr_test_files(dmso_data: list[Mint01]) -> None:
    spec = make_spec("dmso_water", dmso_data)
    check_molecule_against_mol(spec, "DMSO", EPSR_TEST / "DMSO.mol", 200)
    check_molecule_against_mol(
        spec, "WaterMolecule", EPSR_TEST / "WaterMolecule.mol", 400
    )
    inp = EPSR_TEST / "200DMSO400Water.EPSR.inp"
    assert spec.density == inp_values(inp)["rho"][0][0]
    pcof = inp_values(EPSR_TEST / "200DMSO400Water.pcof")
    assert spec.schedule.ereq == pcof["ereq"][0][0]
    check_against_epsr_files(spec, inp, EPSR_TEST)


# --------------------------------------------------------------------------- benzene
def test_benzene_spec(benzene_data: list[Mint01]) -> None:
    spec = make_spec("benzene", benzene_data)
    assert spec.validate(check_files=True) == []
    assert spec.density == pytest.approx(0.0811, abs=2e-4)
    assert spec.mass_density() == pytest.approx(
        0.876 * (6 * 12.011 + 6 * 2.0) / (6 * 12.011 + 6 * 1.008), rel=1e-3
    )
    assert spec.labels == ["CA", "HA"] and spec.n_atoms == 2400
    assert spec.temperature == 300.0 and spec.schedule.ereq == 20.0
    (benzene,) = spec.molecules
    assert benzene.n_atoms == 12 and len(benzene.bonds) == 12 and benzene.is_connected()
    assert len(benzene.angle_triples()) == 18
    assert all(1.0 < benzene.bond_length(i, j) < 1.5 for i, j in benzene.bonds)
    assert benzene.angle_value(1, 2, 3) == pytest.approx(120.0, abs=1.0)
    assert spec.atom_type("CA").rmax_value == pytest.approx(0.5 * 3.55)
    assert [dataset.label for dataset in spec.datasets] == ["C6H6", "C6D6", "HD50"]
    h6, d6, mixture = spec.datasets
    assert composition(h6, "HA") == [(0, 1.0)]
    assert composition(d6, "HA") == [(2, 1.0)]
    assert composition(mixture, "HA") == [(1, 0.5), (2, 0.5)]
    assert all(
        dataset.exchangeable == [] and dataset.rshmin is None
        for dataset in spec.datasets
    )


# --------------------------------------------------------------------------- Cu-BTC
def test_cu_btc_spec(cu_btc_data: list[Mint01]) -> None:
    spec = make_spec("cu_btc", cu_btc_data)
    assert spec.validate(check_files=True) == []
    assert spec.kind == "pdb" and spec.name == "CuBTC" and spec.molecules == []
    assert spec.pdb is not None and Path(spec.pdb.file).is_file()
    assert spec.pdb.file_name == "CuBTC.pdb"
    assert (
        spec.pdb.vdw_factor == 0.68
        and spec.pdb.dihedrals is False
        and spec.pdb.tolerance == 0.0
    )
    assert spec.pdb.step_sizes == {
        "stepmi": 0.1,
        "stepri": 0.3,
        "stepr": 0.3,
        "stepm": 1.0,
    }
    assert spec.labels == ["CU", "O", "C1", "C2", "C3", "H"]
    assert spec.temperature == 300.0 and spec.schedule.ereq == 10.0
    assert spec.epsr_keywords == {
        "qwidthqmax": [0.2, 5.0],
        "hklqmin": [0.0, 10.0, 0.0],
        "qmin": [0.05, 0.0],
    }
    assert spec.density == pytest.approx(624 / 26.3336**3)
    (dataset,) = spec.datasets
    assert dataset.label == "Empty-CuBTC" and dataset.rshmin == 0.7
    assert dataset.isotopes == {} and dataset.exchangeable == []
    assert spec.atom_type("H").mass == 1.008 and spec.atom_type("CU").element == "Cu"


def test_cu_btc_matches_validation_files(cu_btc_data: list[Mint01]) -> None:
    spec = make_spec("cu_btc", cu_btc_data)
    folder = VALIDATION / "cu-btc"
    pdb_lines = (folder / "CuBTC.pdb").read_text().splitlines()
    cell = float(
        next(line for line in pdb_lines if line.startswith("CRYST1")).split()[1]
    )
    n_atoms = sum(line.startswith("ATOM") for line in pdb_lines)
    assert (cell, n_atoms) == (26.3336, 624)
    assert spec.density == pytest.approx(n_atoms / cell**3)
    answers = (folder / "readpdb.stdin").read_text().split()
    assert spec.pdb is not None
    assert spec.pdb.vdw_factor == float(answers[0]) and spec.pdb.dihedrals == (
        answers[1] == "y"
    )
    # changeato.stdin: temp, tol, step sizes, then per atom epsilon / sigma / charge (atmass for H)
    atoms: list[dict[str, float]] = []
    settings: dict[str, float] = {}
    for line in (folder / "changeato.stdin").read_text().splitlines():
        words = line.split()
        if words[0] == "atom":
            atoms.append({})
        elif len(words) == 2 and words[0] not in ("save",):
            (atoms[-1] if atoms else settings)[words[0]] = float(words[1])
    assert (
        settings["temp"] == spec.temperature and settings["tol"] == spec.pdb.tolerance
    )
    assert {key: settings[key] for key in spec.pdb.step_sizes} == spec.pdb.step_sizes
    assert len(atoms) == len(spec.atom_types)
    for values, atom_type in zip(atoms, spec.atom_types, strict=True):
        assert (values["epsilon"], values["sigma"], values["charge"]) == (
            atom_type.epsilon,
            atom_type.sigma,
            atom_type.charge,
        )
        if "atmass" in values:
            assert atom_type.mass == values["atmass"]
    assert spec.schedule.ereq == steps_ereq(folder / "steps.txt")
    check_against_epsr_files(spec, folder / "CuBTC.EPSR.inp", folder)


# --------------------------------------------------------------------------- auto_spec and errors
def test_auto_spec_water(water_data: list[Mint01]) -> None:
    spec, reason = auto_spec(water_data)
    assert spec.source == "library:water"
    assert "water" in reason and all(data.name in reason for data in water_data)
    assert {dataset.label: composition(dataset, "HW") for dataset in spec.datasets} == {
        "H2O": [(0, 1.0)],
        "HDO": [(1, 0.5), (2, 0.5)],
        "D2O": [(2, 1.0)],
    }
    assert spec.validate(check_files=True) == []


def test_auto_spec_other_systems(silica_data, dmso_data, cu_btc_data) -> None:
    assert auto_spec(silica_data)[0].source == "library:silica"
    assert auto_spec(dmso_data)[0].source == "library:dmso_water"
    assert auto_spec(cu_btc_data, name="cell")[0].name == "cell"


def test_auto_spec_failures(water_data, silica_data, tmp_path: Path) -> None:
    with pytest.raises(LookupError, match="do not all point to one"):
        auto_spec(water_data + silica_data)
    unknown = write_mint01(tmp_path / "x.mint01", "mystery liquid")
    with pytest.raises(LookupError, match="no built-in system"):
        auto_spec([unknown])
    with pytest.raises(LookupError):
        auto_spec([])
    # a data set that names nothing is accepted beside ones that name the system
    spec, _ = auto_spec(water_data[:1] + [unknown])
    assert spec.source == "library:water" and composition(spec.datasets[1], "HW") == [
        (0, 1.0)
    ]


def test_wrong_mixture_is_caught(water_data, silica_data, dmso_data) -> None:
    with pytest.raises(ValueError, match="looks like water, not silica"):
        make_spec("silica", water_data)
    with pytest.raises(ValueError, match="looks like silica, not water"):
        make_spec("water", silica_data)
    with pytest.raises(ValueError, match="looks like dmso_water, not water"):
        make_spec("water", dmso_data)
    with pytest.raises(ValueError, match="looks like water, not benzene"):
        make_spec("benzene", water_data)


def test_validate_catches_a_bad_isotope_type(water_data: list[Mint01]) -> None:
    spec = make_spec("water", water_data)
    spec.datasets[0].isotopes = {"H1": [Isotope(2, 1.0)]}
    with pytest.raises(SpecError, match="unknown atom type 'H1'"):
        spec.validate()


def test_unknown_key_and_adjust(water_data: list[Mint01]) -> None:
    with pytest.raises(LookupError, match="unknown system 'argon'"):
        make_spec("argon", water_data)
    with pytest.raises(ValueError, match="unknown adjust keyword"):
        make_spec("water", water_data, pressure=1.0)
    with pytest.raises(ValueError, match="no molecule 'Ice'"):
        make_spec("water", water_data, counts={"Ice": 10})
    with pytest.raises(ValueError, match="at least one data set"):
        make_spec("water", [])


def test_adjust(water_data: list[Mint01]) -> None:
    spec = make_spec(
        "water",
        water_data,
        name="small",
        density=0.0999,
        temperature=298.0,
        ereq=2.5,
        counts={"Water": 500},
        title="Small water box",
        equilibration_iterations=5,
        refinement_runs=1,
        refinement_iterations=7,
    )
    assert spec.validate(check_files=True) == []
    assert spec.name == "small" and spec.title == "Small water box"
    assert spec.density == 0.0999 and spec.temperature == 298.0
    assert spec.molecules[0].count == 500 and spec.n_atoms == 1500
    assert spec.schedule.ereq == 2.5
    assert spec.schedule.iterations() == [
        ("scratch", 1, 0.0),
        ("reference potential", 5, 0.0),
        ("EP refinement", 7, 2.5),
    ]
    assert "adjusted density = 0.0999" in spec.notes
    assert make_spec("water", water_data, count=250).molecules[0].count == 250
    assert (
        KNOWN_SYSTEMS["water"].build(water_data).molecules[0].count == 1000
    )  # the template is untouched


def test_unknown_labelling_is_natural_and_labels_are_unique(tmp_path: Path) -> None:
    mystery = write_mint01(tmp_path / "run1.mint01", "mystery liquid")
    first = write_mint01(tmp_path / "a.mint01", "H2O at 25C")
    second = write_mint01(tmp_path / "b.mint01", "H2O at 50C")
    spec = make_spec("water", [mystery, first, second])
    assert [dataset.label for dataset in spec.datasets] == ["run1", "H2O", "H2O_2"]
    assert composition(spec.datasets[0], "HW") == [(0, 1.0)]
    assert "no H/D labelling recognised" in spec.notes
    assert spec.validate(check_files=True) == []


@pytest.mark.parametrize("key", list(KNOWN_SYSTEMS))
def test_json_round_trip(
    key: str, water_data, silica_data, dmso_data, cu_btc_data, benzene_data
) -> None:
    datasets = {
        "water": water_data,
        "silica": silica_data,
        "dmso_water": dmso_data,
        "cu_btc": cu_btc_data,
        "benzene": benzene_data,
    }[key]
    spec = make_spec(key, datasets)
    as_dict = spec.to_dict()
    json.dumps(as_dict)
    again = SystemSpec.from_json(spec.to_json())
    assert again.to_dict() == as_dict
    assert again.validate(check_files=True) == []


# --------------------------------------------------------------------------- read_mol
def test_read_mol_water() -> None:
    mol = read_mol(VALIDATION / "water" / "Water.mol")
    assert mol.name == "Water" and mol.temperature == 300.0
    assert mol.molecule.labels == ["OW", "HW", "HW"]
    assert mol.molecule.bonds == [(1, 2), (1, 3)]
    assert mol.molecule.bond_lengths == {"1-2": 0.976, "1-3": 0.976}
    assert mol.molecule.angles == {"2-1-3": 107.134}
    assert [atom_type.label for atom_type in mol.atom_types] == ["OW", "HW"]
    ow = mol.atom_types[0]
    assert (
        ow.element,
        ow.epsilon,
        ow.sigma,
        ow.mass,
        ow.charge,
        ow.nsph,
        ow.rmin,
        ow.rmax,
        ow.cg,
    ) == ("O", 0.6503, 3.165492, 15.9994, -0.82, 1, 0.0, 1.55, 1.0)


def test_read_mol_bonds_from_neighbours_and_errors(tmp_path: Path) -> None:
    path = tmp_path / "Pair.mol"
    path.write_text(
        " .gmol 0\n atom 1 001 0.0 0.0 0.0 1 2\n atom 2 002 1.0 0.0 0.0 1 1\nchangelabel 001 A\nchangelabel 002 B\npotential A 0.1 3.0 12.0 0.0 C\n"
    )
    mol = read_mol(path)
    assert mol.molecule.bonds == [(1, 2)] and mol.molecule.labels == ["A", "B"]
    assert mol.atom_types[0].rmax is None and mol.atom_types[0].element == "C"
    empty = tmp_path / "Empty.mol"
    empty.write_text(" .gmol 0\n")
    with pytest.raises(ValueError, match="no atom lines"):
        read_mol(empty)
    bad = tmp_path / "Bad.mol"
    bad.write_text(" atom 1 001 x y z 0\n")
    with pytest.raises(ValueError, match="cannot parse"):
        read_mol(bad)
