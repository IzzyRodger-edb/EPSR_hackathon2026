"""The auxiliary programs through ``epsr.tools``."""

from __future__ import annotations

from pathlib import Path

import pytest

import epsr
from epsr import tools

from epsr_cases import CASE_FILES_TO_SKIP, TEMPLATE, TEST_PROBLEMS, copy_files, load_case

NWTS_SAMPLES = ["DMSODWaterD", "DMSODWaterH", "DMSODWaterNull"]
CASES_FOLDER = TEST_PROBLEMS / "200DMSO400Water"


def _case_folder(case_name: str) -> Path:
    case_folder = CASES_FOLDER / case_name
    if not case_folder.exists():
        pytest.skip("test_problems not available")
    return case_folder


def _file_contents(folder: Path) -> dict[str, bytes]:
    return {path.name: path.read_bytes() for path in folder.iterdir() if path.is_file()}


def test_commands() -> None:
    command_names = tools.commands()
    assert len(command_names) == 35
    expected_subset = {"nwts", "xwts", "partials", "makeato", "mixato", "fmole", "dockato"}
    assert expected_subset <= set(command_names)


def test_unknown_tool(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unknown EPSR tool"):
        tools.run(tmp_path, "no_such_tool")
    with pytest.raises(FileNotFoundError):
        tools.run(tmp_path / "missing", "nwts", ["x.NWTS.dat"])


@pytest.mark.parametrize("sample", NWTS_SAMPLES)
def test_nwts_reproduces_shipped_weights(template: Path, tmp_path: Path, sample: str) -> None:
    """EPSR-test-ex20: nwts rebuilds the .wts files shipped with EPSR-test byte for byte."""
    work_folder = copy_files(template, tmp_path / "work")
    weights_file = f"{sample}.NWTStot.wts"
    (work_folder / weights_file).unlink()
    tools.nwts(work_folder, f"{sample}.NWTS.dat")
    assert (work_folder / weights_file).read_bytes() == (template / weights_file).read_bytes()


def test_keywords_table() -> None:
    nwts_table = tools.keywords("nwts")
    assert nwts_table.extinp == ".NWTS.dat"
    assert nwts_table.dataprompt == "component"
    assert nwts_table.global_keywords[:2] == ["fnameato", "output"]
    assert nwts_table.block_keywords == [
        "atom",
        "atsymbol",
        "iexchange",
        "abundances",
        "isubs",
        "abund2",
    ]
    assert nwts_table.defaults()["output"] == "1"
    assert nwts_table.symbols[0].name == "fnameato"
    assert "partials" in epsr._core.tool_tables()
    with pytest.raises(ValueError):
        tools.keywords("no_such_table")


def test_write_keyword_file_round_trip(template: Path, tmp_path: Path) -> None:
    """Rebuild DMSODWaterD.NWTS.dat from its parsed contents; nwts gives the same weights."""
    work_folder = copy_files(template, tmp_path / "work")
    original = epsr.read_keyword_file(template / "DMSODWaterD.NWTS.dat")
    (work_folder / "DMSODWaterD.NWTS.dat").unlink()

    rebuilt_path = tools.write_keyword_file(
        work_folder / "Rebuilt.NWTS.dat", "nwts", original.values, original.blocks
    )
    rebuilt = epsr.read_keyword_file(rebuilt_path)
    assert rebuilt.title == "Rebuilt"
    assert rebuilt.values == original.values
    assert rebuilt.blocks == original.blocks

    tools.nwts(work_folder, "Rebuilt.NWTS.dat")
    rebuilt_weights = (work_folder / "Rebuilt.NWTStot.wts").read_bytes()
    assert rebuilt_weights == (template / "DMSODWaterD.NWTStot.wts").read_bytes()


def test_write_keyword_file_defaults_and_errors(tmp_path: Path) -> None:
    written_path = tools.write_keyword_file(
        tmp_path / "test.PARTIALS.dat", "partials", {"fnameato": "x.ato"}
    )
    written = epsr.read_keyword_file(written_path)
    assert written.values["fnameato"] == "x.ato"
    for keyword, default in tools.keywords("partials").defaults().items():
        if keyword != "fnameato" and default.strip():
            assert written.values[keyword] == default.strip()
    with pytest.raises(ValueError, match="unknown"):
        tools.write_keyword_file(tmp_path / "bad.dat", "partials", {"bogus": 1})


def test_interactive_tool_with_answers(tmp_path: Path) -> None:
    """makemole (EPSR-test-ex28) driven by answers given as a list of lines."""
    case_folder = _case_folder("EPSR-test-ex28")
    work_folder = copy_files(case_folder, tmp_path / "work", CASE_FILES_TO_SKIP)
    command = next(
        command for command in load_case(case_folder).commands if command.program == "makemole"
    )
    answers = []
    if command.stdin_file:
        answers = (case_folder / command.stdin_file).read_text().splitlines()

    contents_before = _file_contents(work_folder)
    tools.run(work_folder, command.program, command.args, answers)
    assert _file_contents(work_folder) != contents_before  # the tool wrote or updated files


def test_missing_keyword_file(tmp_path: Path) -> None:
    copy_files(TEMPLATE, tmp_path)
    tools.run(tmp_path, "nwts", ["missing.NWTS.dat"])  # faithful to the Fortran: a silent no-op
    assert not (tmp_path / "missing.NWTS.dat").exists()
    with pytest.raises(FileNotFoundError):
        tools.nwts(tmp_path, "missing.NWTS.dat")


def test_tool_errors_are_epsr_errors(tmp_path: Path) -> None:
    """EPSR-test-ex30: makeato reproduces the Fortran crash as EpsrUnsupportedError."""
    case_folder = _case_folder("EPSR-test-ex30")
    answers = (case_folder / "makeato_Ow.stdin").read_text()
    with pytest.raises(epsr.EpsrUnsupportedError, match="SIGSEGV"):
        tools.run(tmp_path, "makeato", ["Ow"], answers)
    assert (tmp_path / "Ow.ato").stat().st_size == 0  # an empty file, as the Fortran leaves it
