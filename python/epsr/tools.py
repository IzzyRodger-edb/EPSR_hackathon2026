"""The auxiliary EPSR programs (setup and analysis tools).

All 35 programs of the EPSR suite are available through :func:`run`, with the
same arguments as on the command line (``epsr-rs <folder> <command> <args>``).

* **Keyword-file programs** (``nwts``, ``xwts``, ``partials``, ``coord``,
  ``sharm``, ...) take the name of their ``.dat`` file as argument and rewrite
  it with defaults filled in. :func:`write_keyword_file` creates such a file
  from a dict.
* **Interactive programs** (``makeato``, ``mixato``, ``fmole``, ``changeato``,
  ``readxyz``...) read their answers from standard input; pass them as
  ``answers``. The ``*.stdin`` files of the ``test_problems/*/EPSR-test-ex*``
  cases show the expected order of answers for each program.

Outputs are files in the folder; read the column files with :func:`epsr.io.read_table`.

**Working directory.** Like the Fortran, several programs (the box builders
``makeato``, ``mixato``, ``fmole``, ``addato``, ``randomise``, ``makelattice``,
``makelatticeato``) read ``system_commands.txt`` and write ``dirlist.txt`` in the
*current working directory*, because the EPSR shell runs every program from
inside the simulation folder. :func:`run` therefore changes the working directory
to ``folder`` for the duration of the call and passes ``./`` as the folder, exactly
as ``epsr-rs`` is used by the EPSR verification scripts. Because the working
directory is process-wide, tool runs are serialised by a lock (one tool at a time
per process), and other threads should not rely on relative paths while a tool
runs. :class:`epsr.Simulation` does not depend on the working directory and is
not affected.
"""

from __future__ import annotations

import os
import threading
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any

from . import _core

__all__ = [
    "commands",
    "run",
    "keywords",
    "ToolKeyword",
    "ToolKeywords",
    "write_keyword_file",
    "nwts",
    "xwts",
    "partials",
    "coord",
    "triangles",
    "torangles",
    "writexyz",
]

#: Serialises tool runs: they change the process-wide working directory.
_WORKING_DIRECTORY_LOCK = threading.RLock()

#: Column layout of the keyword files EPSR writes (Fortran format ``(a11,1x,a,15x,a)``).
KEYWORD_COLUMN_WIDTH = 11
DESCRIPTION_GAP = " " * 15


def commands() -> list[str]:
    """Names of the 35 auxiliary programs."""
    return list(_core.tool_commands())


@contextmanager
def _working_directory(folder: Path) -> Iterator[None]:
    """Temporarily make ``folder`` the working directory (one thread at a time)."""
    with _WORKING_DIRECTORY_LOCK:
        previous_directory = os.getcwd()
        os.chdir(folder)
        try:
            yield
        finally:
            os.chdir(previous_directory)


def _standard_input_text(answers: str | Sequence[str]) -> str:
    """The answers as the text a program reads from standard input (newline-terminated)."""
    if isinstance(answers, str):
        if answers and not answers.endswith("\n"):
            return answers + "\n"
        return answers
    return "".join(f"{answer}\n" for answer in answers)


def run(
    folder: str | PathLike[str],
    command: str,
    args: Sequence[str] = (),
    answers: str | Sequence[str] = "",
    *,
    verbose: bool = False,
) -> None:
    """Run the auxiliary program ``command`` in ``folder``.

    Parameters
    ----------
    folder:
        Working folder; files are read from and written to it in place.
    command:
        Program name, see :func:`commands`.
    args:
        Words after the command name on the EPSR command line (normally the
        input file name), e.g. ``["DMSODWaterD.NWTS.dat"]``.
    answers:
        Standard input for interactive programs: a string, or a sequence of
        lines (joined with newlines).

    Runs with the working directory set to ``folder`` (see the module notes).
    Raises :class:`epsr.EpsrError` subclasses on failure. Releases the GIL.

    Like the Fortran programs, keyword-file tools do nothing when their input
    file is missing; the typed helpers (:func:`nwts`, :func:`partials`...)
    raise ``FileNotFoundError`` instead.
    """
    folder = Path(folder).resolve()
    if not folder.is_dir():
        raise FileNotFoundError(f"folder {folder} does not exist")
    arguments = [str(argument) for argument in args]
    with _working_directory(folder):
        _core.run_tool("./", command, arguments, _standard_input_text(answers), verbose)


@dataclass(frozen=True)
class ToolKeyword:
    """One entry of a program's keyword table."""

    name: str
    value_type: str  # "c" text, "i" integer, "r" real
    flag: int  # Fortran ``iflag``
    default: str
    description: str


@dataclass(frozen=True)
class ToolKeywords:
    """The keyword table of a keyword-file program (from ``symtab_tools.rs``)."""

    name: str
    extinp: str  # input-file extension, e.g. ".NWTS.dat"
    dataprompt: str  # name of the repeated block ("component", "distribution"...) or ""
    dataprompttext: str
    symbols: tuple[ToolKeyword, ...]
    #: Fortran markers, e.g. ``nsymboluptodata``: the number of keywords before the first block.
    marks: Mapping[str, int]

    @property
    def n_global(self) -> int:
        """Number of keywords before the first repeated block."""
        if self.dataprompt and "nsymboluptodata" in self.marks:
            return self.marks["nsymboluptodata"]
        return len(self.symbols)

    @property
    def global_keywords(self) -> list[str]:
        return [keyword.name for keyword in self.symbols[: self.n_global]]

    @property
    def block_keywords(self) -> list[str]:
        end_of_block = self.marks.get("nsymboldata", len(self.symbols))
        return [keyword.name for keyword in self.symbols[self.n_global : end_of_block]]

    def defaults(self) -> dict[str, str]:
        return {keyword.name: keyword.default for keyword in self.symbols}

    def descriptions(self) -> dict[str, str]:
        return {keyword.name: keyword.description for keyword in self.symbols}


def keywords(command: str) -> ToolKeywords:
    """The keyword table of ``command`` (``nwts``, ``partials``, ``coord``...)."""
    table = _core.tool_keywords(command)
    symbols = tuple(
        ToolKeyword(
            name=str(name),
            value_type=str(value_type),
            flag=int(flag),
            default=str(default),
            description=str(description),
        )
        for name, value_type, flag, default, description in table["symbols"]
    )
    return ToolKeywords(
        name=table["name"],
        extinp=table["extinp"],
        dataprompt=table["dataprompt"],
        dataprompttext=table["dataprompttext"],
        symbols=symbols,
        marks=dict(table["marks"]),
    )


def _format_value(value: Any) -> str:
    """Write a Python value the way EPSR keyword files expect it."""
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, (list, tuple)):
        return " ".join(_format_value(item) for item in value)
    return str(value)


def _check_keywords(given: Mapping[str, Any], allowed: list[str], what: str) -> None:
    unknown = set(given) - set(allowed)
    if unknown:
        raise ValueError(f"unknown {what} {sorted(unknown)}; valid: {allowed}")


def write_keyword_file(
    path: str | PathLike[str],
    command: str,
    values: Mapping[str, Any],
    blocks: Sequence[Mapping[str, Any]] | None = None,
    *,
    title: str | None = None,
) -> Path:
    """Write a keyword input file for a keyword-file program.

    ``values`` sets global keywords; keywords not given take the table
    defaults. ``blocks`` gives one mapping per repeated block (``component``
    for ``nwts``, ...); keywords missing in a block take the defaults. Without
    ``blocks`` no block section is written, and programs that derive their
    blocks from another file (``nwts`` from the ``.ato`` file) fill them in.
    The file uses the layout EPSR writes (keyword, value, description, ``q``);
    the program rewrites it on exit.
    """
    path = Path(path)
    table = keywords(command)
    blocks = list(blocks or [])
    _check_keywords(values, table.global_keywords, f"{command} keywords")
    if blocks and not table.dataprompt:
        raise ValueError(f"{command} input files have no repeated blocks")
    for block in blocks:
        _check_keywords(block, table.block_keywords, f"{command} block keywords")

    if title is None:
        title = (
            path.name.removesuffix(table.extinp) if path.name.endswith(table.extinp) else path.stem
        )
    defaults = table.defaults()
    descriptions = table.descriptions()

    def keyword_line(keyword: str, value: Any) -> str:
        return (
            f"{keyword:<{KEYWORD_COLUMN_WIDTH}} {_format_value(value)}"
            f"{DESCRIPTION_GAP}{descriptions[keyword]}"
        )

    lines = [f"{title}{DESCRIPTION_GAP}Title of this file"]
    for keyword in table.global_keywords:
        lines.append(keyword_line(keyword, values.get(keyword, defaults[keyword])))
    for block_number, block in enumerate(blocks, start=1):
        lines += ["", f"{table.dataprompt} {block_number:>3}", ""]
        for keyword in table.block_keywords:
            lines.append(keyword_line(keyword, block.get(keyword, defaults[keyword])))
    lines.append("q")
    path.write_text("\n".join(lines) + "\n", encoding="latin-1")
    return path


def _run_keyword_tool(
    command: str, folder: str | PathLike[str], dat_file: str, verbose: bool
) -> None:
    """Run a keyword-file program, refusing a missing input file.

    The Fortran program (and :func:`run`) silently does nothing in that case.
    """
    input_file = Path(folder) / dat_file
    if not input_file.is_file():
        raise FileNotFoundError(f"{command} input file {input_file} does not exist")
    run(folder, command, [dat_file], verbose=verbose)


def nwts(folder: str | PathLike[str], dat_file: str, *, verbose: bool = False) -> None:
    """Neutron weights from a ``.NWTS.dat`` file. Raises FileNotFoundError if it is missing."""
    _run_keyword_tool("nwts", folder, dat_file, verbose)


def xwts(folder: str | PathLike[str], dat_file: str, *, verbose: bool = False) -> None:
    """X-ray weights from a ``.XWTS.dat`` file. Raises FileNotFoundError if it is missing."""
    _run_keyword_tool("xwts", folder, dat_file, verbose)


def partials(folder: str | PathLike[str], dat_file: str, *, verbose: bool = False) -> None:
    """Partial g(r) from a ``.PARTIALS.dat`` file. Raises FileNotFoundError if it is missing."""
    _run_keyword_tool("partials", folder, dat_file, verbose)


def coord(folder: str | PathLike[str], dat_file: str, *, verbose: bool = False) -> None:
    """Coordination numbers from a ``.COORD.dat`` file.

    Raises FileNotFoundError if the file is missing.
    """
    _run_keyword_tool("coord", folder, dat_file, verbose)


def triangles(folder: str | PathLike[str], dat_file: str, *, verbose: bool = False) -> None:
    """Bond angles from a ``.TRIANGLES.dat`` file. Raises FileNotFoundError if it is missing."""
    _run_keyword_tool("triangles", folder, dat_file, verbose)


def torangles(folder: str | PathLike[str], dat_file: str, *, verbose: bool = False) -> None:
    """Torsion angles from a ``.TORANGLES.dat`` file. Raises FileNotFoundError if it is missing."""
    _run_keyword_tool("torangles", folder, dat_file, verbose)


def writexyz(folder: str | PathLike[str], ato_file: str, answers: str | Sequence[str] = "") -> None:
    """Run ``writexyz`` on ``ato_file`` (answers as in the interactive program)."""
    run(folder, "writexyz", [ato_file], answers)
