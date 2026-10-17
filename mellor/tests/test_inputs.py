"""``mellor.inputs``: the generated files reproduce the validated inputs of ``validation/``."""

from __future__ import annotations

import math
from pathlib import Path

import epsr
import pytest
from specs import VALIDATION, cu_btc_spec, silica_spec, water_spec

from mellor import inputs
from mellor.spec import EPSR_KEYWORDS, Atom, AtomType, SpecError

WATER_CASE = VALIDATION / "water"
SILICA_CASE = VALIDATION / "silica_bulk"
CU_BTC_CASE = VALIDATION / "cu-btc"


# --------------------------------------------------------------------------- helpers
def words_equal(first: str, second: str) -> bool:
    """List-directed comparison: word by word, numbers as floats."""
    first_words, second_words = first.split(), second.split()
    if len(first_words) != len(second_words):
        return False
    for one, other in zip(first_words, second_words, strict=True):
        try:
            value, reference = float(one), float(other)
        except ValueError:
            if one != other:
                return False
        else:
            if not math.isclose(value, reference, rel_tol=1e-6, abs_tol=1e-5):
                return False
    return True


def assert_lines_equivalent(text: str, reference_path: Path) -> None:
    """Every non-blank line of ``text`` equals the reference line, list-directed."""
    lines = [line for line in text.splitlines() if line.strip()]
    reference = [
        line for line in reference_path.read_text().splitlines() if line.strip()
    ]
    assert len(lines) == len(
        reference
    ), f"{len(lines)} lines, expected {len(reference)}"
    for line, reference_line in zip(lines, reference, strict=True):
        assert words_equal(line, reference_line), f"{line!r} != {reference_line!r}"


def parse_keyword_file(text: str, name: str, folder: Path) -> epsr.io.KeywordFile:
    path = folder / name
    path.write_text(text, encoding="latin-1")
    return epsr.read_keyword_file(path)


def assert_keyword_files_equal(
    generated: epsr.io.KeywordFile, reference: epsr.io.KeywordFile
) -> None:
    assert generated.title == reference.title
    assert list(generated.values) == list(reference.values)
    for keyword, value in reference.values.items():
        assert words_equal(generated.values[keyword], value), (
            keyword,
            generated.values[keyword],
            value,
        )
    assert generated.block_keyword == reference.block_keyword
    assert len(generated.blocks) == len(reference.blocks)
    for number, (block, reference_block) in enumerate(
        zip(generated.blocks, reference.blocks, strict=True), start=1
    ):
        assert list(block) == list(reference_block), f"block {number}"
        for keyword, value in reference_block.items():
            assert words_equal(block[keyword], value), (
                number,
                keyword,
                block[keyword],
                value,
            )


# --------------------------------------------------------------------------- .mol
def test_render_mol_water_equals_validation():
    spec = water_spec()
    text = inputs.render_mol(spec.molecules[0], spec)
    assert_lines_equivalent(text, WATER_CASE / "Water.mol")


def test_render_mol_water_is_byte_identical():
    """The layout mirrors makemole's own rewrite format, down to the byte."""
    spec = water_spec()
    text = inputs.render_mol(spec.molecules[0], spec)
    assert text == (WATER_CASE / "Water.mol").read_text()


@pytest.mark.parametrize("index, file_name", [(0, "Si.mol"), (1, "O.mol")])
def test_render_mol_silica_equals_validation(index, file_name):
    spec = silica_spec()
    text = inputs.render_mol(spec.molecules[index], spec)
    assert_lines_equivalent(text, SILICA_CASE / file_name)
    assert text == (SILICA_CASE / file_name).read_text()


def test_render_mol_structure_follows_the_bond_list():
    spec = water_spec()
    lines = inputs.render_mol(spec.molecules[0], spec).splitlines()
    atom_lines = [line.split() for line in lines if line.startswith(" atom")]
    assert [words[2] for words in atom_lines] == ["001", "002", "003"]
    # neighbour count and neighbour list of the oxygen, then of each hydrogen
    assert atom_lines[0][6:] == ["2", "2", "3"]
    assert atom_lines[1][6:] == ["1", "1"]
    assert atom_lines[2][6:] == ["1", "1"]
    assert [line for line in lines if line.startswith("bond")] == [
        "bond 001 002    0.97600",
        "bond 001 003    0.97600",
    ]
    assert [line for line in lines if line.startswith("angle")] == [
        "angle 002 001 003  107.13400"
    ]
    potential = [line.split() for line in lines if line.startswith("potential")]
    assert potential[0][1:] == [
        "OW", "0.650300", "3.165492", "15.999400", "-0.820000", "O",
        "0.100000E+01", "0.000000", "1.550000", "0.100000E+01",
    ]  # fmt: skip
    assert potential[1][1] == "HW" and potential[1][9] == "0.000000"


def test_render_mol_default_rmax_and_rotations():
    spec = water_spec()
    molecule = spec.molecules[0]
    molecule.rotations = [(1, 2)]
    spec.atom_types[0] = AtomType("OW", "O", 0.6503, 3.165492, 15.9994, -0.82)
    lines = inputs.render_mol(molecule, spec).splitlines()
    assert "rot   1   2" in lines
    potential = next(line.split() for line in lines if line.startswith("potential OW"))
    assert float(potential[9]) == pytest.approx(0.5 * 3.165492, abs=1e-6)


def test_render_mol_refuses_too_many_atoms():
    spec = water_spec()
    molecule = spec.molecules[0]
    molecule.atoms = [Atom("OW", 0.0, 0.0, 0.4 * i) for i in range(1000)]
    with pytest.raises(SpecError, match="three digits"):
        inputs.render_mol(molecule, spec)


# --------------------------------------------------------------------------- .NWTS.dat
@pytest.mark.parametrize(
    "index, file_name", [(0, "D2O.NWTS.dat"), (1, "HDO.NWTS.dat"), (2, "H2O.NWTS.dat")]
)
def test_render_nwts_water_equals_validation(index, file_name, tmp_path):
    spec = water_spec()
    text = inputs.render_nwts(spec.datasets[index], spec)
    generated = parse_keyword_file(text, file_name, tmp_path)
    reference = epsr.read_keyword_file(WATER_CASE / file_name)
    assert_keyword_files_equal(generated, reference)
    assert generated.blocks[1]["iexchange"] == "1"


def test_render_nwts_silica_equals_validation(tmp_path):
    spec = silica_spec()
    text = inputs.render_nwts(spec.datasets[0], spec)
    generated = parse_keyword_file(text, "SiO2.NWTS.dat", tmp_path)
    reference = epsr.read_keyword_file(SILICA_CASE / "SiO2.NWTS.dat")
    assert_keyword_files_equal(generated, reference)


def test_render_nwts_layout():
    spec = water_spec()
    lines = inputs.render_nwts(spec.datasets[1], spec).splitlines()
    assert lines[0].startswith("HDO.NWTS" + " " * 15 + "Title of this file")
    assert lines[1].startswith(
        "fnameato    waterbox.ato" + " " * 15 + "Name of .ato file"
    )
    assert "component   1" in lines and "component   2" in lines
    assert lines[-1] == "q"
    assert "abundances  1 0.5 2 0.5" in " ".join(lines)


# --------------------------------------------------------------------------- .EPSR.inp
def test_render_inp_water_equals_validation(tmp_path):
    spec = water_spec()
    generated = parse_keyword_file(inputs.render_inp(spec), spec.inp_name, tmp_path)
    reference = epsr.read_keyword_file(WATER_CASE / "waterbox.EPSR.inp")
    # rho is derived from the .ato file by EPSR; the generated file carries the real density
    assert generated.values.pop("rho") == "0.1"
    reference.values.pop("rho")
    assert_keyword_files_equal(generated, reference)
    assert [block["datafile"] for block in generated.blocks] == [
        "SLS18502-D2O.mint01", "SLS18500-HDO5050.mint01", "SLS18498-H2O.mint01",
    ]  # fmt: skip


def test_render_inp_silica_equals_validation(tmp_path):
    spec = silica_spec()
    generated = parse_keyword_file(inputs.render_inp(spec), spec.inp_name, tmp_path)
    reference = epsr.read_keyword_file(SILICA_CASE / "silica.EPSR.inp")
    generated.values.pop("rho")
    reference.values.pop("rho")
    assert_keyword_files_equal(generated, reference)
    assert generated.values["qmin"] == "0.1 0.0"
    assert generated.blocks[0]["rshmin"] == "1.2"


def test_render_inp_cu_btc_equals_validation(tmp_path):
    spec = cu_btc_spec()
    generated = parse_keyword_file(inputs.render_inp(spec), spec.inp_name, tmp_path)
    reference = epsr.read_keyword_file(CU_BTC_CASE / "CuBTC.EPSR.inp")
    generated.values.pop("rho")
    reference.values.pop("rho")
    assert_keyword_files_equal(generated, reference)


def test_render_inp_scratch_start_and_overrides():
    spec = water_spec(
        epsr_keywords={
            "qmin": [0.05, 0.0],
            "fwhmq": 0.03,
            "nq": 500,
            "potfac": [1.0, 0.5],
        }
    )
    spec.schedule.ntimes = 7
    text = inputs.render_inp(spec)
    values = {
        line.split()[0]: line.split()[1:]
        for line in text.splitlines()
        if line.strip() and line[0].isalpha()
    }
    assert text.startswith("waterbox.EPSR" + " " * 15 + "Title of this file\n")
    assert values["ireset"][0] == "1" and values["iinit"][0] == "1"
    assert values["nsumt"][0] == "-1" and values["niter"][0] == "1"
    assert values["ntimes"][0] == "7"
    assert values["fwhmq"][0] == "0.03" and values["nq"][0] == "500"
    assert values["potfac"][:2] == ["1.0", "0.5"]
    assert (
        values["fnameato"][0] == "waterbox.ato"
        and values["fnamepcof"][0] == "waterbox.pcof"
    )
    assert values["ndata"][0] == "3"
    assert text.rstrip().endswith("\nq")
    assert inputs.render_inp(spec, quick=True) == text


def test_every_epsr_keyword_is_in_the_template():
    _, global_entries, block_entries = inputs._template_entries()
    keywords = {keyword for keyword, _, _ in global_entries}
    assert set(EPSR_KEYWORDS) <= keywords
    assert [keyword for keyword, _, _ in block_entries] == [
        "datafile", "wtsfile", "nrtype", "rshmin", "szeros", "tweak", "efilereq", "porodadd",
    ]  # fmt: skip


def test_format_value():
    assert inputs.format_value(0.1) == "0.1"
    assert inputs.format_value(300.0) == "300.0"
    assert inputs.format_value([0.05, 0.0]) == "0.05 0.0"
    assert inputs.format_value(True) == "T" and inputs.format_value(False) == "F"
    assert inputs.format_value(5) == "5"
    assert inputs.format_value(1e-05) == "1E-05"
    assert inputs.format_value("none 0 0") == "none 0 0"


# --------------------------------------------------------------------------- answers
def test_mixato_answers():
    assert inputs.mixato_answers(water_spec()) == [
        "1",
        "Water",
        "1000",
        "0.1",
        "waterbox",
    ]
    assert inputs.mixato_answers(silica_spec()) == [
        "2",
        "Si",
        "1500",
        "O",
        "3000",
        "0.0664",
        "silica",
    ]
    assert (
        inputs.mixato_answers(water_spec())
        == (WATER_CASE / "mixato.stdin").read_text().split()
    )
    assert (
        inputs.mixato_answers(silica_spec())
        == (SILICA_CASE / "mixato.stdin").read_text().split()
    )
    assert inputs.mixato_answers(cu_btc_spec()) == []


def test_readpdb_answers():
    assert (
        inputs.readpdb_answers(cu_btc_spec())
        == (CU_BTC_CASE / "readpdb.stdin").read_text().splitlines()
    )
    spec = cu_btc_spec()
    spec.pdb.dihedrals = True
    assert inputs.readpdb_answers(spec) == ["0.68", "y"]
    assert inputs.readpdb_answers(water_spec()) == []


def test_changeato_commands():
    assert (
        inputs.changeato_commands(cu_btc_spec())
        == (CU_BTC_CASE / "changeato.stdin").read_text().splitlines()
    )


def test_changeato_atmass_only_when_the_mass_is_not_epsrs_element_mass():
    spec = cu_btc_spec()
    commands = inputs.changeato_commands(spec)
    assert (
        commands.count("atmass 1.008") == 1
        and sum(c.startswith("atmass") for c in commands) == 1
    )
    spec.atom_types[0] = AtomType("CU", "Cu", 0.02092, 3.144, 63.55, 1.098)
    spec.pdb.extra_commands = ["qradius 1 0.5"]
    commands = inputs.changeato_commands(spec)
    assert "atmass 63.55" in commands
    assert commands[-4:] == ["qradius 1 0.5", "save", "y", "e"]
    assert inputs.EPSR_ELEMENT_MASSES["H"] == 1.00797


def test_pdb_atom_labels_and_order_check():
    labels = inputs.pdb_atom_labels(CU_BTC_CASE / "CuBTC.pdb")
    assert labels == ["CU", "O", "C1", "C2", "C3", "H"]
    spec = cu_btc_spec()
    spec.atom_types = [spec.atom_types[1], spec.atom_types[0]] + spec.atom_types[2:]
    with pytest.raises(SpecError, match="PDB atom names"):
        inputs.generate(spec)


# --------------------------------------------------------------------------- rshmin
def test_auto_rshmin():
    assert inputs.auto_rshmin(water_spec()) == 0.7
    assert inputs.auto_rshmin(silica_spec()) == 1.2
    assert inputs.auto_rshmin(cu_btc_spec()) == 0.7
    spec = water_spec()
    spec.molecules[0].bond_lengths = {"1-2": 0.5, "1-3": 0.976}
    assert inputs.auto_rshmin(spec) == pytest.approx(0.3)


def test_dataset_rshmin_in_the_inp():
    spec = water_spec()
    spec.datasets[0].rshmin = 0.9
    blocks = [
        line.split()
        for line in inputs.render_inp(spec).splitlines()
        if line.startswith("rshmin")
    ]
    assert [words[1] for words in blocks] == ["0.9", "0.7", "0.7"]


# --------------------------------------------------------------------------- generate / write
def test_generate_water():
    spec = water_spec()
    generated = inputs.generate(spec)
    assert set(generated.mol_files) == {"Water.mol"}
    assert list(generated.nwts_files) == [
        "D2O.NWTS.dat",
        "HDO.NWTS.dat",
        "H2O.NWTS.dat",
    ]
    assert generated.inp_name == "waterbox.EPSR.inp"
    assert generated.inp_text.startswith("waterbox.EPSR")
    assert all(path.is_absolute() for path in generated.data_files)
    assert [path.name for path in generated.data_files] == [
        "SLS18502-D2O.mint01", "SLS18500-HDO5050.mint01", "SLS18498-H2O.mint01",
    ]  # fmt: skip
    assert generated.pdb_file is None
    assert generated.setup_steps == [
        "makemole Water", "mixato < mixato.stdin", "randomise waterbox", "delete Water.mol",
        "fmole waterbox 100", "nwts D2O.NWTS.dat", "nwts HDO.NWTS.dat", "nwts H2O.NWTS.dat",
    ]  # fmt: skip
    assert set(generated.files()) == {
        "Water.mol",
        "D2O.NWTS.dat",
        "HDO.NWTS.dat",
        "H2O.NWTS.dat",
        "waterbox.EPSR.inp",
    }
    assert generated.answer_files() == {
        "mixato.stdin": "1\nWater\n1000\n0.1\nwaterbox\n"
    }


def test_generate_silica_and_cu_btc():
    silica = inputs.generate(silica_spec())
    assert set(silica.mol_files) == {"Si.mol", "O.mol"}
    assert silica.setup_steps == [
        "makemole Si", "makemole O", "mixato < mixato.stdin", "randomise silica",
        "delete Si.mol", "delete O.mol", "nwts SiO2.NWTS.dat",
    ]  # fmt: skip
    cu_btc = inputs.generate(cu_btc_spec())
    assert cu_btc.mol_files == {} and cu_btc.mixato_answers == []
    assert cu_btc.pdb_file == CU_BTC_CASE / "CuBTC.pdb"
    assert cu_btc.pdb_name == "CuBTC.pdb"
    assert cu_btc.setup_steps == [
        "readpdb CuBTC.pdb < readpdb.stdin", "changeato CuBTC.ato < changeato.stdin", "nwts CuBTC.NWTS.dat",
    ]  # fmt: skip
    assert set(cu_btc.answer_files()) == {"readpdb.stdin", "changeato.stdin"}
    assert cu_btc.answer_files()["readpdb.stdin"] == "0.68\nn\n"
    assert set(cu_btc.files()) == {
        "CuBTC.NWTS.dat",
        "CuBTC.EPSR.inp",
        "vanderWaalsRadii.txt",
    }
    assert "vanderWaalsRadii.txt" not in inputs.generate(silica_spec()).files()


def test_van_der_waals_radii_equal_validation():
    """``readpdb`` reads the radii from the folder; the text is EPSR's table."""
    text = inputs.render_vdw_radii()
    reference = (CU_BTC_CASE / "vanderWaalsRadii.txt").read_text()
    rows = [line.split() for line in text.splitlines()]
    reference_rows = [line.split() for line in reference.splitlines() if line.split()]
    assert len(rows) == len(reference_rows) == 103
    for row, reference_row in zip(rows, reference_rows, strict=True):
        assert words_equal(" ".join(row), " ".join(reference_row)), (row, reference_row)
    assert rows[28] == ["29", "Cu", "2"]


def test_generate_validates_first():
    spec = water_spec(density=-1.0)
    with pytest.raises(SpecError, match="density"):
        inputs.generate(spec)
    spec = water_spec()
    spec.molecules[0].name = "waterbox"
    with pytest.raises(SpecError, match="run name"):
        inputs.generate(spec)
    spec = water_spec()
    spec.datasets[1].file = spec.datasets[0].file
    with pytest.raises(SpecError, match="same file name"):
        inputs.generate(spec)


def test_write_water(tmp_path):
    spec = water_spec()
    generated = inputs.generate(spec)
    written = inputs.write(generated, tmp_path / "epsr")
    names = {path.name for path in written}
    assert names == {
        "Water.mol", "D2O.NWTS.dat", "HDO.NWTS.dat", "H2O.NWTS.dat", "waterbox.EPSR.inp", "mixato.stdin",
        "SLS18502-D2O.mint01", "SLS18500-HDO5050.mint01", "SLS18498-H2O.mint01",
    }  # fmt: skip
    assert all(path.is_file() for path in written)
    assert (tmp_path / "epsr" / "Water.mol").read_text() == generated.mol_files[
        "Water.mol"
    ]
    assert (
        tmp_path / "epsr" / "SLS18498-H2O.mint01"
    ).read_bytes() == generated.data_files[2].read_bytes()


def test_write_cu_btc_copies_the_pdb(tmp_path):
    generated = inputs.generate(cu_btc_spec())
    written = inputs.write(generated, tmp_path)
    names = {path.name for path in written}
    assert {
        "CuBTC.pdb",
        "readpdb.stdin",
        "changeato.stdin",
        "CuBTC.NWTS.dat",
        "CuBTC.EPSR.inp",
        "Empty-CuBTC.mint01",
        "vanderWaalsRadii.txt",
    } <= names
    assert (tmp_path / "changeato.stdin").read_text() == (
        CU_BTC_CASE / "changeato.stdin"
    ).read_text()


def test_write_missing_data_file(tmp_path):
    spec = water_spec()
    spec.datasets[0].file = str(tmp_path / "missing.mint01")
    generated = inputs.generate(spec)
    with pytest.raises(FileNotFoundError, match="missing.mint01"):
        inputs.write(generated, tmp_path / "epsr")
