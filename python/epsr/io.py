"""Read-only parsers for the text files EPSR writes.

The engine is the only writer of EPSR files; these helpers only read them, using
numpy. Nothing here changes a simulation folder.
"""

from __future__ import annotations

import importlib
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

__all__ = [
    "Table",
    "KeywordFile",
    "read_table",
    "read_table_series",
    "read_erg",
    "read_keyword_file",
    "ERG_COLUMNS",
]

FloatArray = npt.NDArray[np.float64]

#: Columns of the ``.erg`` file, one line per iteration (written by ``epsr.rs``,
#: Fortran format ``(10(1x,e13.6))``).
ERG_COLUMNS = (
    "energy",  # energd: intermolecular energy per molecule (kJ/mol)
    "pressure",  # pressd (kbar)
    "r_factor_raw",  # rfacraw
    "ep_amplitude",  # energabs: absolute energy of the empirical potential
    "density",  # rho (atoms/Å³)
    "ep_up",  # epcup
    "ep_coefficient",  # epccoeff: r-space coefficient
    "parab_a",  # aparab: trend-line fit
    "parab_b",  # bparab
    "parab_c",  # cparab
)

#: In keyword files the value is separated from its description by a long run of blanks.
VALUE_DESCRIPTION_SEPARATOR = re.compile(r"\s{8,}")

#: Marker EPSR writes after the run name on the first line of a keyword file.
TITLE_MARKER = "Title of this file"


@dataclass
class Table:
    """One EPSR column file (``.g01``, ``.f01``, ``.u01``, ...).

    ``columns[label]`` is the value column of each header label and
    ``errors[label]`` the matching error/spread column, when the file has one.
    Labels that repeat in the header (``.z01`` lists each pair twice) get a
    ``#2`` suffix on the second occurrence.
    """

    path: Path
    x_label: str
    x: FloatArray
    columns: dict[str, FloatArray]
    errors: dict[str, FloatArray] = field(default_factory=dict)

    @property
    def labels(self) -> list[str]:
        return list(self.columns)

    def __getitem__(self, label: str) -> FloatArray:
        return self.columns[label]

    def __len__(self) -> int:
        return len(self.x)

    def to_dataframe(self) -> Any:
        """Return a ``pandas.DataFrame`` (requires the ``pandas`` extra)."""
        pandas = import_optional("pandas", extra="pandas")
        data: dict[str, Any] = {self.x_label: self.x}
        for label, values in self.columns.items():
            data[label] = values
            if label in self.errors:
                data[f"{label} err"] = self.errors[label]
        return pandas.DataFrame(data)


def _make_labels_unique(labels: Iterable[str]) -> list[str]:
    """Append ``#2``, ``#3``... to labels that occur more than once."""
    times_seen: dict[str, int] = {}
    unique_labels = []
    for label in labels:
        count = times_seen.get(label, 0) + 1
        times_seen[label] = count
        unique_labels.append(label if count == 1 else f"{label}#{count}")
    return unique_labels


def read_table(path: str | PathLike[str]) -> Table:
    """Read one EPSR column file.

    The first line is a ``#`` header with the x label and one label per
    column group; the numeric rows are ``x y1 [e1] y2 [e2] ...``.
    """
    path = Path(path)
    with path.open(encoding="latin-1") as file:
        header = file.readline()
    if not header.startswith("#"):
        raise ValueError(f"{path} is not an EPSR column file (no '#' header line)")

    header_words = header[1:].split()
    if header_words:
        x_label, labels = header_words[0], _make_labels_unique(header_words[1:])
    else:
        x_label, labels = "x", []

    data = np.loadtxt(path, comments="#", ndmin=2, encoding="latin-1")
    n_value_columns = data.shape[1] - 1
    columns: dict[str, FloatArray] = {}
    errors: dict[str, FloatArray] = {}

    if labels and n_value_columns == 2 * len(labels):
        # value and error column for every label
        for index, label in enumerate(labels):
            columns[label] = data[:, 1 + 2 * index]
            errors[label] = data[:, 2 + 2 * index]
    elif labels and n_value_columns == len(labels):
        # one value column per label, no errors
        for index, label in enumerate(labels):
            columns[label] = data[:, 1 + index]
    else:
        # header does not match the data: number the columns
        for index in range(n_value_columns):
            columns[f"col{index + 1}"] = data[:, 1 + index]

    return Table(path=path, x_label=x_label, x=data[:, 0], columns=columns, errors=errors)


def read_table_series(first_file: str | PathLike[str]) -> Table:
    """Read ``name.EXT01`` and any continuation files ``EXT02``, ``EXT03``...

    EPSR writes at most 100 columns per file; wide outputs continue in the
    next numbered file. The tables are merged into one.
    """
    first_file = Path(first_file)
    number_match = re.search(r"(\d\d)$", first_file.name)
    if not number_match:
        return read_table(first_file)

    name_without_number = first_file.name[: number_match.start()]
    table = read_table(first_file)
    file_number = int(number_match.group(1)) + 1
    while (next_file := first_file.with_name(f"{name_without_number}{file_number:02d}")).exists():
        continuation = read_table(next_file)
        for label, values in continuation.columns.items():
            key = label if label not in table.columns else f"{label}#{file_number}"
            table.columns[key] = values
            if label in continuation.errors:
                table.errors[key] = continuation.errors[label]
        file_number += 1
    return table


def read_erg(path: str | PathLike[str]) -> FloatArray:
    """Read the ``.erg`` energy history: an ``(n_iterations, 10)`` array.

    Column meanings are listed in :data:`ERG_COLUMNS`.
    """
    data = np.loadtxt(Path(path), ndmin=2)
    if data.size == 0:
        return np.empty((0, len(ERG_COLUMNS)))
    return data.reshape(-1, len(ERG_COLUMNS))


@dataclass
class KeywordFile:
    """A parsed keyword file (``.EPSR.inp``, ``.pcof`` header, ``*.NWTS.dat``...).

    ``values`` holds the keywords before the first block; ``blocks`` holds one
    dict per repeated block (``data n`` in ``.EPSR.inp``, ``component n`` in
    ``.NWTS.dat``...). Each value is the raw text between the keyword and its
    description.
    """

    path: Path
    title: str
    values: dict[str, str]
    block_keyword: str | None = None
    blocks: list[dict[str, str]] = field(default_factory=list)
    descriptions: dict[str, str] = field(default_factory=dict)
    #: Lines after ``q`` (the pair section of a ``.pcof`` file).
    trailer: list[str] = field(default_factory=list)

    def get(self, key: str, default: str | None = None) -> str | None:
        return self.values.get(key, default)

    def words(self, key: str) -> list[str]:
        """The value of ``key`` split into words (empty if absent)."""
        value = self.values.get(key)
        return value.split() if value else []


def _block_keywords() -> set[str]:
    """Words that start a repeated block: ``data`` and every tool's ``dataprompt``."""
    block_keywords = {"data"}
    try:
        from . import _core

        for table_name in _core.tool_tables():
            dataprompt = _core.tool_keywords(table_name)["dataprompt"]
            if dataprompt:
                block_keywords.add(dataprompt)
    except ImportError:  # extension not built: fall back to the common ones
        block_keywords.update({"component", "distribution", "centre"})
    return block_keywords


def _has_title_line(first_line: str) -> bool:
    """EPSR writes ``<name>   Title of this file`` first; ``.pcof`` files have no title."""
    return first_line.rstrip().endswith(TITLE_MARKER) or len(first_line.split()) == 1


def read_keyword_file(path: str | PathLike[str]) -> KeywordFile:
    """Parse an EPSR keyword file.

    The optional first line is the title (``<name>  Title of this file``).
    Lines are ``keyword  value(s)  description``; values are separated from
    the description by a run of at least 8 blanks (EPSR writes 15).
    A line ``q`` ends the keyword section.
    """
    path = Path(path)
    lines = path.read_text(encoding="latin-1").splitlines()
    has_title = bool(lines) and _has_title_line(lines[0])
    title = lines[0].split()[0] if has_title else ""
    first_keyword_line = 1 if has_title else 0

    block_keywords = _block_keywords()
    values: dict[str, str] = {}
    descriptions: dict[str, str] = {}
    blocks: list[dict[str, str]] = []
    block_keyword: str | None = None
    current_section = values
    trailer: list[str] = []

    for line_number, line in enumerate(lines[first_keyword_line:], start=first_keyword_line):
        words = line.split()
        if not words:
            continue
        if words == ["q"]:
            trailer = lines[line_number + 1 :]
            break
        is_block_header = words[0] in block_keywords and len(words) == 2 and words[1].isdigit()
        if is_block_header:
            block_keyword = words[0]
            current_section = {}
            blocks.append(current_section)
            continue
        keyword = words[0]
        text_after_keyword = line.lstrip()[len(keyword) :].lstrip()
        value_and_description = VALUE_DESCRIPTION_SEPARATOR.split(text_after_keyword, maxsplit=1)
        current_section[keyword] = value_and_description[0].strip()
        if len(value_and_description) > 1 and keyword not in descriptions:
            descriptions[keyword] = value_and_description[1].strip()

    return KeywordFile(
        path=path,
        title=title,
        values=values,
        block_keyword=block_keyword,
        blocks=blocks,
        descriptions=descriptions,
        trailer=trailer,
    )


def import_optional(module_name: str, extra: str) -> Any:
    """Import an optional dependency, with an installation hint if it is missing."""
    try:
        return importlib.import_module(module_name)
    except ImportError as error:  # pragma: no cover - depends on environment
        raise ImportError(
            f"{module_name} is required for this feature; install it with "
            f"`uv sync --extra {extra}` (or `uv add epsr-py[{extra}]`)"
        ) from error
