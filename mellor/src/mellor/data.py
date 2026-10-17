"""Reading and describing Gudrun ``.mint01`` structure-factor files.

A ``.mint01`` file is a text table with header lines starting with ``#``:

    # sls18498.mdcs01                 <- source file (instrument run)
    #  H2O in 1mm can 5 in position 5 : 25C     <- sample title
    #     1  1001                     <- (record count)
    #    10
    # Bychkov                         <- user
    # 28-FEB-2001 16:15:03            <- date
    # ...spec.bad / groups_def.dat ...
    Q [Å⁻¹]   F(Q) [barn/sr/atom]   error     <- data columns (3)

Trailing points with zero value and zero error are padding. The run number is in
the source name (``sls18498`` -> SANDALS run 18498, ``NIMROD00001122`` -> NIMROD).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "Mint01",
    "read_mint01",
    "describe",
    "describe_all",
    "plot_data",
    "guess_sample",
]

#: Instrument prefixes of the source file name -> instrument name.
INSTRUMENT_PREFIXES = {
    "SLS": "SANDALS (ISIS)",
    "NIMROD": "NIMROD (ISIS)",
    "GEM": "GEM (ISIS)",
    "POLARIS": "POLARIS (ISIS)",
    "D4": "D4 (ILL)",
}

_DATE_PATTERN = re.compile(r"\b\d{1,2}-[A-Za-z]{3}-\d{4}\b")
_RUN_PATTERN = re.compile(r"[-_]?0*(\d+)")
_GENERIC_SOURCE_PATTERN = re.compile(r"^[A-Za-z]+[-_]?0*(\d+)")
_NUMERIC_LINE_PATTERN = re.compile(r"^[-+\d.\sEeDd]+$")

# The dataviz reference palette: categorical slots in their validated order, a
# blue ramp for more than eight series, and the chart ink.
_SERIES_COLOURS = (
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#6250d6",
    "#e34948",
)
_SEQUENTIAL_BLUES = (
    "#86b6ef",
    "#6da7ec",
    "#5598e7",
    "#3987e5",
    "#2a78d6",
    "#256abf",
    "#1c5cab",
    "#184f95",
    "#104281",
    "#0d366b",
)
_INK = {
    "primary": "#0b0b0b",
    "secondary": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
    "surface": "#fcfcfb",
}


@dataclass
class Mint01:
    """One ``.mint01`` file: its header and the data columns."""

    path: Path
    source: str  # first header line without '#', e.g. 'sls18498.mdcs01'
    title: str  # second header line, the sample description
    header: list[str]  # every header line, without the leading '#'
    q: np.ndarray  # Å⁻¹, all rows
    f: np.ndarray  # F(Q) (differential cross section), all rows
    err: np.ndarray  # error, all rows
    user: str = ""
    date: str = ""

    @property
    def name(self) -> str:
        return self.path.name

    @property
    def instrument(self) -> str:
        """Instrument guessed from the source name ('' if unknown)."""
        return _instrument_and_run(self.source)[0]

    @property
    def run_number(self) -> str:
        """Run number from the source name ('' if none)."""
        return _instrument_and_run(self.source)[1]

    @property
    def valid(self) -> np.ndarray:
        """Boolean mask of the rows that carry data (not the zero padding at the end)."""
        finite = np.isfinite(self.q) & np.isfinite(self.f) & np.isfinite(self.err)
        return finite & ~((self.f == 0.0) & (self.err == 0.0))

    @property
    def q_range(self) -> tuple[float, float]:
        """Lowest and highest Q with data."""
        q = self.q[self.valid]
        if q.size == 0:
            raise ValueError(
                f"{self.name}: no rows with data (every value and error is zero)"
            )
        return float(q.min()), float(q.max())

    @property
    def q_step(self) -> float:
        """Median Q spacing."""
        q = self.q[self.valid]
        if q.size < 2:
            return 0.0
        return float(np.median(np.diff(q)))

    @property
    def n_points(self) -> int:
        return int(np.count_nonzero(self.valid))


def _instrument_and_run(source: str) -> tuple[str, str]:
    """``(instrument, run number)`` from a source name such as ``sls18498.mdcs01``."""
    stem = source.strip().split(".", 1)[0]
    upper = stem.upper()
    for prefix in sorted(INSTRUMENT_PREFIXES, key=len, reverse=True):
        if upper.startswith(prefix):
            match = _RUN_PATTERN.match(stem[len(prefix) :])
            return INSTRUMENT_PREFIXES[prefix], match.group(1) if match else ""
    match = _GENERIC_SOURCE_PATTERN.match(stem)
    return "", match.group(1) if match else ""


def _header_fields(header: list[str]) -> tuple[str, str, str, str]:
    """``(source, title, user, date)`` from the header lines (Gudrun's fixed layout)."""
    source = header[0] if header else ""
    title = header[1] if len(header) > 1 else ""
    user = header[4] if len(header) > 4 else ""
    date = header[5] if len(header) > 5 else ""
    if not _DATE_PATTERN.search(date):
        date = next((line for line in header if _DATE_PATTERN.search(line)), "")
    if _NUMERIC_LINE_PATTERN.match(user) or "\\" in user or "/" in user:
        user = ""
    return source, title, user, date


def _number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def read_mint01(path: str | Path) -> Mint01:
    """Read a ``.mint01`` file. Raises ``ValueError`` if it has no three-column data."""
    path = Path(path)
    try:
        text = path.read_text(encoding="latin-1")
    except OSError as error:
        raise ValueError(f"cannot read {path}: {error}") from error

    header: list[str] = []
    rows: list[list[float]] = []
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if not rows:
                header.append(stripped[1:].strip())
            continue
        tokens = stripped.split()
        try:
            values = [_number(token) for token in tokens]
        except ValueError:
            raise ValueError(
                f"{path.name}, line {number}: not a numeric data row: {stripped!r}"
            ) from None
        if len(values) != 3:
            raise ValueError(
                f"{path.name}, line {number}: expected 3 columns (Q, F(Q), error), found {len(values)}"
            )
        rows.append(values)
    if not rows:
        raise ValueError(
            f"{path.name}: no three-column data rows (only header or empty lines)"
        )

    table = np.array(rows, dtype=float)
    source, title, user, date = _header_fields(header)
    data = Mint01(
        path=path,
        source=source,
        title=title,
        header=header,
        q=table[:, 0],
        f=table[:, 1],
        err=table[:, 2],
        user=user,
        date=date,
    )
    if data.n_points == 0:
        raise ValueError(
            f"{path.name}: every row has zero value and zero error; no data to fit"
        )
    return data


def describe(data: Mint01) -> dict[str, Any]:
    """A JSON-serialisable description for the model and the report.

    Keys: ``file``, ``source``, ``title``, ``instrument``, ``run_number``, ``user``,
    ``date``, ``n_points``, ``q_min``, ``q_max``, ``q_step``, ``f_min``, ``f_max``,
    ``f_high_q`` (mean F over the last 10 % of the Q range, the high-Q level),
    ``mean_error``, ``sample`` (from :func:`guess_sample`).
    """
    valid = data.valid
    q, f, err = data.q[valid], data.f[valid], data.err[valid]
    q_min, q_max = data.q_range
    high_q = q >= q_max - 0.1 * (q_max - q_min)
    return {
        "file": data.name,
        "path": str(data.path.resolve()),
        "source": data.source,
        "title": data.title,
        "instrument": data.instrument,
        "run_number": data.run_number,
        "user": data.user,
        "date": data.date,
        "n_points": data.n_points,
        "q_min": round(q_min, 6),
        "q_max": round(q_max, 6),
        "q_step": round(data.q_step, 6),
        "f_min": round(float(f.min()), 6),
        "f_max": round(float(f.max()), 6),
        "f_high_q": round(float(f[high_q].mean()), 6),
        "mean_error": round(float(err.mean()), 6),
        "sample": guess_sample(data.title, data.name),
    }


def describe_all(paths: list[str | Path]) -> list[dict[str, Any]]:
    return [describe(read_mint01(path)) for path in paths]


# --------------------------------------------------------------------------- sample hints
#: Fraction of deuterium in "null" water (zero mean coherent scattering length):
#: 0.64 H / 0.36 D, the composition of ``EPSR-test/DMSODWaterNull.NWTS.dat``.
NULL_WATER_D_FRACTION = 0.36

#: Library key -> (canonical keyword, patterns in a lower-cased title or file name).
#: Water is only a candidate when no other system is mentioned (it is the solvent
#: of most samples), so the order matters: solutes first.
_SYSTEM_PATTERNS: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    ("dmso_water", "dmso", (r"dmso", r"dimethyl\s*sulph?oxide")),
    ("benzene", "benzene", (r"benzene", r"c6h6", r"c6d6")),
    ("cu_btc", "cubtc", (r"cu[\s_-]?btc", r"hkust", r"(?<![a-z])mof(?![a-z])")),
    ("silica", "silica", (r"sio2", r"silica")),
    ("argon", "argon", (r"argon",)),
    ("water", "water", (r"water", r"h2o", r"d2o", r"hdo")),
)
_WATER_KEY = "water"
_PERCENT_PAIR = re.compile(r"(?<![\d.])(\d{1,3})\s*[/:]\s*(\d{1,3})(?![\d.])")


def _find(pattern: str, text: str) -> re.Match[str] | None:
    """Search ``pattern`` as a word of letters and digits, so ``2h2o`` and ``h2o/d2o`` match."""
    return re.search(rf"(?<![a-z]){pattern}(?![a-z])", text)


def _hydrogen_labelling(
    text: str,
) -> tuple[str | None, float | None, str | None, list[str]]:
    """``(hydrogen, hd_fraction, sample name, keywords)`` from a lower-cased text."""
    keywords: list[str] = []
    h2o = _find(r"h2o", text) or _find(r"light\s+water", text)
    d2o = _find(r"d2o", text) or _find(r"heavy\s+water", text)
    hdo = _find(r"hdo", text)
    null = _find(r"null", text)
    fifty = re.search(
        r"50\s*[/:]\s*50|5050|(?<![a-z])h\s*/\s*d(?![a-z])|equimolar", text
    )
    for found, word in ((h2o, "h2o"), (d2o, "d2o"), (hdo, "hdo"), (null, "null")):
        if found:
            keywords.append(word)
    if fifty:
        keywords.append("50:50")

    if null:
        return "HD", NULL_WATER_D_FRACTION, "Null", keywords
    if hdo or fifty or (h2o and d2o):
        fraction = 0.5
        if h2o and d2o:
            pair = _PERCENT_PAIR.search(text)
            if pair and int(pair.group(1)) + int(pair.group(2)) == 100:
                first, second = int(pair.group(1)), int(pair.group(2))
                fraction = (second if h2o.start() < d2o.start() else first) / 100.0
        water = h2o or d2o or hdo or _find(r"waters?", text)
        sample = "HDO" if water else f"HD{round(fraction * 100)}"
        return "HD", fraction, sample, keywords
    if d2o:
        return "D", 1.0, "D2O", keywords
    if h2o:
        return "H", 0.0, "H2O", keywords
    if _find(r"c6d6", text) or re.search(r"benzene-d6|d6-benzene", text):
        return "D", 1.0, "C6D6", keywords + ["c6d6"]
    if _find(r"c6h6", text):
        return "H", 0.0, "C6H6", keywords + ["c6h6"]
    if re.search(r"deuterated|perdeuterated", text):
        return "D", 1.0, "D", keywords + ["deuterated"]
    if re.search(r"protiated|hydrogenated", text):
        return "H", 0.0, "H", keywords + ["protiated"]
    return None, None, None, keywords


def guess_sample(title: str, file_name: str = "") -> dict[str, Any]:
    """Hints about the sample from the title and file name.

    Returns ``{"keywords": [...], "hydrogen": "H" | "D" | "HD" | None, "hd_fraction":
    float | None, "candidates": [library keys, best first]}``. Recognises e.g.
    ``H2O``, ``D2O``, ``HDO``/``5050``/``50:50``, ``SiO2``/``silica``, ``CuBTC``,
    ``DMSO``, ``benzene``/``C6H6``/``C6D6``, ``argon``. ``hd_fraction`` is the
    deuterium fraction of the exchangeable hydrogen (1.0 for D2O, 0.5 for HDO 50:50).

    Also returns ``"sample"``: the short sample name used as the data set label
    (``H2O``, ``D2O``, ``HDO``, ``Null`` for null water, ``C6H6``, ``C6D6``) or None.
    The title is read first; the file name only when the title says nothing.
    """
    title_text = (title or "").lower()
    file_text = Path(file_name).stem.lower() if file_name else ""
    hydrogen, fraction, sample, keywords = _hydrogen_labelling(title_text)
    if hydrogen is None and file_text:
        hydrogen, fraction, sample, keywords = _hydrogen_labelling(file_text)

    text = f"{title_text} {file_text}"
    candidates: list[str] = []
    for key, keyword, patterns in _SYSTEM_PATTERNS:
        if any(re.search(pattern, text) for pattern in patterns):
            if key == _WATER_KEY and candidates:
                continue  # water is the solvent when another system is named
            candidates.append(key)
            if keyword not in keywords:
                keywords.append(keyword)
    if (
        not candidates
        and hydrogen is not None
        and sample in ("H2O", "D2O", "HDO", "Null")
    ):
        candidates.append(_WATER_KEY)
        keywords.append("water")
    return {
        "keywords": keywords,
        "hydrogen": hydrogen,
        "hd_fraction": fraction,
        "candidates": candidates,
        "sample": sample,
    }


# --------------------------------------------------------------------------- figure
def plot_data(datasets: list[Mint01], path: str | Path, *, offset: float = 0.5) -> Path:
    """Overlay of the measured F(Q) of every file (vertically offset), saved as PNG."""
    from matplotlib.figure import Figure  # no pyplot: works without a display

    if not datasets:
        raise ValueError("plot_data: no data sets to plot")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_series = len(datasets)
    if n_series <= len(_SERIES_COLOURS):
        colours = _SERIES_COLOURS[:n_series]
    else:  # beyond the categorical slots: one blue ramp, light to dark, in file order
        colours = tuple(
            _SEQUENTIAL_BLUES[round(i * (len(_SEQUENTIAL_BLUES) - 1) / (n_series - 1))]
            for i in range(n_series)
        )

    figure = Figure(figsize=(8.0, 4.8), dpi=150, facecolor=_INK["surface"])
    axes = figure.add_subplot(111)
    axes.set_facecolor(_INK["surface"])
    for index, (data, colour) in enumerate(zip(datasets, colours, strict=True)):
        valid = data.valid
        q, f = data.q[valid], data.f[valid] + index * offset
        label = f"{data.name}: {data.title}" if data.title else data.name
        axes.plot(q, f, color=colour, linewidth=1.2, label=label)
        if index:
            axes.axhline(index * offset, color=_INK["grid"], linewidth=0.6, zorder=0)
        if n_series <= 4 and q.size:
            axes.annotate(
                data.path.stem,
                xy=(float(q[-1]), float(f[-1])),
                xytext=(4, 0),
                textcoords="offset points",
                fontsize=7,
                color=_INK["secondary"],
                va="center",
            )
    axes.set_xlabel("Q (Å⁻¹)", color=_INK["secondary"])
    axes.set_ylabel(
        "F(Q) (barn sr⁻¹ atom⁻¹)" + (f", offset {offset:g}" if n_series > 1 else ""),
        color=_INK["secondary"],
    )
    axes.set_title(
        "Measured structure factors", fontsize=10, color=_INK["primary"], loc="left"
    )
    axes.grid(True, color=_INK["grid"], linewidth=0.6)
    axes.tick_params(colors=_INK["muted"], labelsize=8)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color(_INK["axis"])
    if n_series > 1:
        axes.legend(fontsize=7, frameon=False, loc="best", labelcolor=_INK["secondary"])
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    return path
