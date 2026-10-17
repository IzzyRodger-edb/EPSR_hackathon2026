"""Run folders, upload validation and safe path handling.

Two rules hold everywhere below.

*The run folder is the trust boundary.* The model's tools can read any file
inside the folder it is given (integration.md section 6), so every run gets its
own folder under a per-user root, and no request is ever allowed to address a
path outside the caller's root.

*Data file names must be blank-free.* EPSR builds the data path as
`directory + name` and truncates it at the first blank, so a file called
`my data.mint01` silently becomes `my`.
"""

from __future__ import annotations

import json
import re
import secrets
import unicodedata
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

#: Mellor's own figure names, in display order (visualisation toolkit README).
#: Used as the fallback when `mellor.figures.FIGURE_CAPTIONS` is unavailable.
FALLBACK_CAPTIONS: dict[str, str] = {
    "01_convergence": "Convergence: energies, pressure, R-factor and R-squared per data "
    "set, empirical-potential amplitude against ereq, density and Monte Carlo acceptance.",
    "02_fq_fit": "Measured against simulated F(Q) for each data set, with the difference "
    "below, at the end of each run.",
    "03_fq_evolution": "Simulated F(Q) through the refinement (light to dark).",
    "04_residual_map": "Data minus simulation over Q and iteration.",
    "05_partials": "Partial pair distribution functions g(r) and partial structure "
    "factors S(Q).",
    "06_potentials": "The empirical potential per atom pair through the refinement, "
    "against the reference potential.",
    "07_first_shell": "First-shell peak position, peak height and coordination number "
    "against iteration.",
    "08_local_structure": "Coordination-number and bond-angle distributions taken from "
    "the atomic coordinates.",
    "09_displacement": "Root-mean-square atom displacement since the first iteration.",
}

#: Figures this app draws itself with seaborn, from the tidy CSV tables
#: (integration.md section 5.1). Captions double as LaTeX captions on upload.
SEABORN_CAPTIONS: dict[str, str] = {
    "sb_fit_overlay": "Measured and simulated F(Q) for each isotopic sample at the end "
    "of the refinement.",
    "sb_residuals": "Residual F(Q) (data minus simulation) for each isotopic sample.",
    "sb_convergence": "R-factor against EPSR iteration, one line per refinement stage.",
    "sb_partials": "Total partial pair distribution functions g(r) for each atom pair.",
    "sb_potentials": "Refined empirical potential against the reference potential, per "
    "atom pair.",
}

_SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def new_run_id() -> str:
    """Sortable, unguessable, filesystem-safe."""
    return datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(3)


def slugify(text: str, *, default: str = "run", limit: int = 40) -> str:
    cleaned = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    cleaned = _SAFE_NAME.sub("-", cleaned).strip("-._")
    return (cleaned[:limit] or default).lower()


def safe_data_filename(filename: str) -> str:
    """A blank-free, traversal-free name for an uploaded data file.

    Keeps the extension the user gave (EPSR never checks it; `.mint01` is only a
    convention) but guarantees the name EPSR will see contains no blank.
    """
    base = Path(filename.replace("\\", "/")).name
    stem, dot, ext = base.rpartition(".")
    if not dot:
        stem, ext = base, "mint01"
    stem = _SAFE_NAME.sub("-", stem).strip("-._")[:60]
    ext = _SAFE_NAME.sub("", ext)[:12] or "mint01"
    return f"{stem or 'data'}.{ext}"


def resolve_within(base: Path, relative: str) -> Path:
    """Join `relative` onto `base`, refusing anything that escapes `base`.

    Raises ValueError rather than returning a path, so a caller cannot forget to
    check the result.
    """
    base = base.resolve()
    if relative.startswith(("/", "\\")) or ".." in Path(relative.replace("\\", "/")).parts:
        raise ValueError("path escapes the run folder")
    candidate = (base / relative).resolve()
    if candidate != base and base not in candidate.parents:
        raise ValueError("path escapes the run folder")
    return candidate


@dataclass(frozen=True)
class DataFileInfo:
    """What we can tell the user about an uploaded file before anything runs."""

    name: str
    title: str
    points: int
    q_min: float
    q_max: float
    warnings: tuple[str, ...]


def validate_mint01(text: str, *, name: str = "data") -> DataFileInfo:
    """Check an uploaded file against the Gudrun `.mint01` contract.

    Mirrors `mellor.data.read_mint01`, which raises `ValueError` on anything that
    is not a three-column Gudrun table, so the user finds out in the upload form
    instead of ten minutes into a job. Header lines are those whose *first
    whitespace token* contains `#`; the second header line is the sample title.
    """
    headers: list[str] = []
    rows: list[tuple[float, float, float]] = []
    bad_lines: list[int] = []

    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if "#" in line.split()[0]:
            headers.append(line.lstrip("#").strip())
            continue
        fields = line.split()
        if len(fields) < 3:
            bad_lines.append(lineno)
            continue
        try:
            row = (float(fields[0]), float(fields[1]), float(fields[2]))
        except ValueError:
            bad_lines.append(lineno)
            continue
        rows.append(row)

    if bad_lines:
        shown = ", ".join(str(n) for n in bad_lines[:5])
        raise ValueError(
            f"{name}: line(s) {shown} are not three numeric columns. A .mint01 file holds "
            "Q, F(Q) and the error on F(Q), with header lines starting with '#'."
        )
    if len(rows) < 10:
        raise ValueError(f"{name}: only {len(rows)} data rows found; this is not a .mint01 file.")
    if len(rows) > 5000:
        raise ValueError(f"{name}: {len(rows)} points exceeds EPSR's limit of 5000 (MQD).")

    q_values = [r[0] for r in rows]
    warnings: list[str] = []
    if any(b <= a for a, b in zip(q_values, q_values[1:])):
        warnings.append("Q is not strictly increasing; EPSR expects ascending histogram edges.")
    # EPSR interpolates onto q[i] = qstep*(i-1), i = 1..nq; with the defaults
    # (nq 600, qstep 0.05) everything above 30 A^-1 is read and then discarded.
    if q_values[-1] > 30.0:
        warnings.append(
            f"data extends to {q_values[-1]:.1f} A^-1; with EPSR's defaults (nq 600, "
            "qstep 0.05) everything above 30 A^-1 is read and discarded."
        )
    if q_values[0] > 1.0:
        warnings.append(f"the first point is at Q = {q_values[0]:.2f} A^-1; low-Q data is missing.")
    if rows[-1][1] != 0.0 or rows[-1][2] != 0.0:
        warnings.append(
            "the last row is not a zero sentinel; EPSR treats Q as a bin edge and drops the "
            "final point, so one real data point will be lost."
        )

    title = headers[1] if len(headers) > 1 else (headers[0] if headers else name)
    return DataFileInfo(
        name=name,
        title=title.strip() or name,
        points=len(rows),
        q_min=q_values[0],
        q_max=q_values[-1],
        warnings=tuple(warnings),
    )


def dir_size(path: Path) -> int:
    total = 0
    if not path.exists():
        return 0
    for item in path.rglob("*"):
        if item.is_file():
            try:
                total += item.stat().st_size
            except OSError:
                pass
    return total


def human_size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ("B", "kB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} GB"


def read_json(path: Path) -> dict | None:
    try:
        with path.open("r", encoding="utf-8") as handle:
            loaded = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    return loaded if isinstance(loaded, dict) else None


def tail_text(path: Path, *, limit: int = 20000) -> str:
    """The last `limit` bytes of a text file, for the log pane."""
    try:
        size = path.stat().st_size
        with path.open("rb") as handle:
            if size > limit:
                handle.seek(size - limit)
                handle.readline()
            return handle.read().decode("utf-8", "replace")
    except OSError:
        return ""
