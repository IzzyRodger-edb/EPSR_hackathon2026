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
Tools hold no global state, so different folders may be processed in parallel
threads; never run two tools on the same folder at once.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any

from . import _core

__all__ = [
    "commands",
    "run",
    "keywords",
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


def commands() -> list[str]:
    """Names of the 35 auxiliary programs."""
    return list(_core.tool_commands())


def _answers(answers: str | Sequence[str]) -> str:
    if isinstance(answers, str):
        return answers if not answers or answers.endswith("\n") else answers + "\n"
    return "".join(f"{a}\n" for a in answers)


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

    Raises :class:`epsr.EpsrError` subclasses on failure. Releases the GIL.
    """
    folder = Path(folder)
    if not folder.is_dir():
        raise FileNotFoundError(f"folder {folder} does not exist")
    _core.run_tool(str(folder), command, [str(a) for a in args], _answers(answers), verbose)


@dataclass(frozen=True)
class ToolKeywords:
    """The keyword table of a keyword-file program (from ``symtab_tools.rs``)."""

    name: str
    extinp: str  # input-file extension, e.g. ".NWTS.dat"
    dataprompt: str  # name of the repeated block ("component", "distribution"...) or ""
    dataprompttext: str
    #: ``(keyword, type, iflag, default, description)``; type is ``c``/``i``/``r``.
    symbols: tuple[tuple[str, str, int, str, str], ...]
    marks: Mapping[str, int]

    @property
    def n_global(self) -> int:
        """Number of keywords before the first repeated block."""
        if self.dataprompt and "nsymboluptodata" in self.marks:
            return self.marks["nsymboluptodata"]
        return len(self.symbols)

    @property
    def global_keywords(self) -> list[str]:
        return [s[0] for s in self.symbols[: self.n_global]]

    @property
    def block_keywords(self) -> list[str]:
        end = self.marks.get("nsymboldata", len(self.symbols))
        return [s[0] for s in self.symbols[self.n_global : end]]

    def defaults(self) -> dict[str, str]:
        return {s[0]: s[3] for s in self.symbols}

    def descriptions(self) -> dict[str, str]:
        return {s[0]: s[4] for s in self.symbols}


def keywords(command: str) -> ToolKeywords:
    """The keyword table of ``command`` (``nwts``, ``partials``, ``coord``...)."""
    d = _core.tool_keywords(command)
    return ToolKeywords(
        name=d["name"],
        extinp=d["extinp"],
        dataprompt=d["dataprompt"],
        dataprompttext=d["dataprompttext"],
        symbols=tuple(tuple(s) for s in d["symbols"]),  # type: ignore[misc]
        marks=dict(d["marks"]),
    )


def _fmt(v: Any) -> str:
    if isinstance(v, bool):
        return "T" if v else "F"
    if isinstance(v, (list, tuple)):
        return " ".join(_fmt(x) for x in v)
    return str(v)


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
    kw = keywords(command)
    allowed = set(kw.global_keywords)
    unknown = set(values) - allowed
    if unknown:
        raise ValueError(f"unknown {command} keywords {sorted(unknown)}; valid: {kw.global_keywords}")
    block_allowed = set(kw.block_keywords)
    blocks = list(blocks or [])
    if blocks and not kw.dataprompt:
        raise ValueError(f"{command} input files have no repeated blocks")
    for b in blocks:
        bad = set(b) - block_allowed
        if bad:
            raise ValueError(f"unknown {command} block keywords {sorted(bad)}; valid: {kw.block_keywords}")
    if title is None:
        title = path.name[: -len(kw.extinp)] if path.name.endswith(kw.extinp) else path.stem
    defaults, descr = kw.defaults(), kw.descriptions()

    def line(k: str, v: Any) -> str:
        return f"{k:<11} {_fmt(v)}{' ' * 15}{descr[k]}"

    out = [f"{title}{' ' * 15}Title of this file"]
    for k in kw.global_keywords:
        out.append(line(k, values.get(k, defaults[k])))
    for n, b in enumerate(blocks, start=1):
        out += ["", f"{kw.dataprompt} {n:>3}", ""]
        for k in kw.block_keywords:
            out.append(line(k, b.get(k, defaults[k])))
    out.append("q")
    path.write_text("\n".join(out) + "\n", encoding="latin-1")
    return path


def _keyword_tool(command: str):  # type: ignore[no-untyped-def]
    def tool(folder: str | PathLike[str], dat_file: str, *, verbose: bool = False) -> None:
        run(folder, command, [dat_file], verbose=verbose)

    tool.__name__ = command
    tool.__doc__ = f"Run ``{command}`` on the keyword file ``dat_file`` in ``folder``."
    return tool


nwts = _keyword_tool("nwts")
xwts = _keyword_tool("xwts")
partials = _keyword_tool("partials")
coord = _keyword_tool("coord")
triangles = _keyword_tool("triangles")
torangles = _keyword_tool("torangles")


def writexyz(folder: str | PathLike[str], ato_file: str, answers: str | Sequence[str] = "") -> None:
    """Run ``writexyz`` on ``ato_file`` (answers as in the interactive program)."""
    run(folder, "writexyz", [ato_file], answers)
