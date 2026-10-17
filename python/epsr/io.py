"""Read-only parsers for the text files EPSR writes.

The engine is the only writer of EPSR files; these helpers only read them, using
numpy. Nothing here changes a simulation folder.
"""

from __future__ import annotations

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

_SPLIT_VALUE = re.compile(r"\s{8,}")


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
    x: npt.NDArray[np.float64]
    columns: dict[str, npt.NDArray[np.float64]]
    errors: dict[str, npt.NDArray[np.float64]] = field(default_factory=dict)

    @property
    def labels(self) -> list[str]:
        return list(self.columns)

    def __getitem__(self, label: str) -> npt.NDArray[np.float64]:
        return self.columns[label]

    def __len__(self) -> int:
        return len(self.x)

    def to_dataframe(self) -> Any:
        """Return a ``pandas.DataFrame`` (requires the ``pandas`` extra)."""
        pd = _import_optional("pandas", "pandas")
        data: dict[str, Any] = {self.x_label: self.x}
        for k, v in self.columns.items():
            data[k] = v
            if k in self.errors:
                data[f"{k} err"] = self.errors[k]
        return pd.DataFrame(data)


def _unique(labels: Iterable[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for lab in labels:
        n = seen.get(lab, 0) + 1
        seen[lab] = n
        out.append(lab if n == 1 else f"{lab}#{n}")
    return out


def read_table(path: str | PathLike[str]) -> Table:
    """Read one EPSR column file.

    The first line is a ``#`` header with the x label and one label per
    column group; the numeric rows are ``x y1 [e1] y2 [e2] ...``.
    """
    path = Path(path)
    with path.open(encoding="latin-1") as f:
        header = f.readline()
    if not header.startswith("#"):
        raise ValueError(f"{path} is not an EPSR column file (no '#' header line)")
    words = header[1:].split()
    x_label, labels = (words[0], words[1:]) if words else ("x", [])
    data = np.loadtxt(path, comments="#", ndmin=2, encoding="latin-1")
    ncol = data.shape[1] - 1
    labels = _unique(labels)
    columns: dict[str, npt.NDArray[np.float64]] = {}
    errors: dict[str, npt.NDArray[np.float64]] = {}
    if labels and ncol == 2 * len(labels):
        for i, lab in enumerate(labels):
            columns[lab] = data[:, 1 + 2 * i]
            errors[lab] = data[:, 2 + 2 * i]
    elif labels and ncol == len(labels):
        for i, lab in enumerate(labels):
            columns[lab] = data[:, 1 + i]
    else:
        for i in range(ncol):
            columns[f"col{i + 1}"] = data[:, 1 + i]
    return Table(path=path, x_label=x_label, x=data[:, 0], columns=columns, errors=errors)


def read_table_series(first: str | PathLike[str]) -> Table:
    """Read ``name.EXT01`` and any continuation files ``EXT02``, ``EXT03``...

    EPSR writes at most 100 columns per file; wide outputs continue in the
    next numbered file. The tables are merged into one.
    """
    first = Path(first)
    m = re.search(r"(\d\d)$", first.name)
    if not m:
        return read_table(first)
    stem = first.name[: m.start()]
    table = read_table(first)
    n = int(m.group(1)) + 1
    while (nxt := first.with_name(f"{stem}{n:02d}")).exists():
        more = read_table(nxt)
        for lab, col in more.columns.items():
            key = lab if lab not in table.columns else f"{lab}#{n}"
            table.columns[key] = col
            if lab in more.errors:
                table.errors[key] = more.errors[lab]
        n += 1
    return table


def read_erg(path: str | PathLike[str]) -> npt.NDArray[np.float64]:
    """Read the ``.erg`` energy history: an ``(n_iterations, 10)`` array.

    Column meanings are listed in :data:`ERG_COLUMNS`.
    """
    data = np.loadtxt(Path(path), ndmin=2)
    return data.reshape(-1, len(ERG_COLUMNS)) if data.size else np.empty((0, len(ERG_COLUMNS)))


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
        v = self.values.get(key)
        return v.split() if v else []


def _block_keywords() -> set[str]:
    kws = {"data"}
    try:
        from . import _core

        for name in _core.tool_tables():
            dp = _core.tool_keywords(name)["dataprompt"]
            if dp:
                kws.add(dp)
    except ImportError:  # extension not built: fall back to the common ones
        kws.update({"component", "distribution", "centre"})
    return kws


def read_keyword_file(path: str | PathLike[str]) -> KeywordFile:
    """Parse an EPSR keyword file.

    Lines are ``keyword  value(s)  description``; values are separated from the
    description by a run of at least 8 blanks (the files EPSR writes use 15).
    A line ``q`` ends the keyword section.
    """
    path = Path(path)
    lines = path.read_text(encoding="latin-1").splitlines()
    title = lines[0].split()[0] if lines and lines[0].split() else ""
    blockkw = _block_keywords()
    values: dict[str, str] = {}
    descr: dict[str, str] = {}
    blocks: list[dict[str, str]] = []
    block_keyword: str | None = None
    current = values
    trailer: list[str] = []
    for i, line in enumerate(lines[1:], start=1):
        w = line.split()
        if not w:
            continue
        if w[0] == "q" and len(w) == 1:
            trailer = lines[i + 1 :]
            break
        if w[0] in blockkw and len(w) == 2 and w[1].isdigit():
            block_keyword = w[0]
            current = {}
            blocks.append(current)
            continue
        kw = w[0]
        rest = line.lstrip()[len(kw) :].lstrip()
        parts = _SPLIT_VALUE.split(rest, maxsplit=1)
        current[kw] = parts[0].strip()
        if len(parts) > 1 and kw not in descr:
            descr[kw] = parts[1].strip()
    return KeywordFile(
        path=path,
        title=title,
        values=values,
        block_keyword=block_keyword,
        blocks=blocks,
        descriptions=descr,
        trailer=trailer,
    )


def _import_optional(module: str, extra: str) -> Any:
    import importlib

    try:
        return importlib.import_module(module)
    except ImportError as e:  # pragma: no cover - depends on environment
        raise ImportError(
            f"{module} is required for this feature; install it with `uv sync --extra {extra}` "
            f"(or `uv add epsr-py[{extra}]`)"
        ) from e
