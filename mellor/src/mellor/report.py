"""The PDF report (reportlab).

Sections (the headings are :data:`SECTION_HEADINGS`):

* Title page: title, "prepared by Mellor" with the model name, date, the data files.
1. The data: one table row per ``.mint01`` (file, title, instrument, run, Q range,
   points) and the data figure.
2. The model: composition (molecules and counts), box, density, temperature; atom
   types with ε, σ, charge, mass; isotopic composition per data set; the
   specification's notes (and whether it came from the model, the library or the user).
3. The calculation: the setup steps, the EPSR schedule table (run, stage, iterations,
   ereq), the ``.EPSR.inp`` keywords that differ from the defaults, software versions.
4. Results: run table (R-factor total and per data set, energy, pressure, EP
   amplitude), fit-quality table (R-factor before/after, R²), first-shell table, the
   figures with captions.
5. Discussion: the model's narrative, headed "Written by the language model
   (<model>); every number in this report was computed by EPSR".
6. Reproducibility: how to re-run (``python run_<name>.py --out <folder>``), the run
   folder layout, SHA-256 hashes of the data files, the specification and the script.
Appendix A: the generated input files (``.mol``, ``.NWTS.dat``, ``.EPSR.inp``);
Appendix B: the standalone Python script.

Every number comes from the run folder (the manifest, the analysis dict, the
specification). Missing analysis keys and missing figures are tolerated: the
corresponding rows or figures are left out, never invented.

Fonts: DejaVu Sans / Sans Mono (from matplotlib's ``mpl-data/fonts/ttf``) so that
Å, ², ⁻¹ and Greek letters render.
"""

from __future__ import annotations

import contextlib
import datetime as dt
import hashlib
import math
import re
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_RIGHT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, StyleSheet1
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Flowable,
    HRFlowable,
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

from . import __version__
from .figures import FIGURE_CAPTIONS
from .inputs import EpsrInputs, dataset_keywords, generate
from .runner import RunResult
from .spec import Dataset, Molecule, SystemSpec

__all__ = [
    "SECTION_HEADINGS",
    "NO_MODEL_PREFIX",
    "build_pdf",
    "markdown_to_flowables",
    "make_styles",
    "register_fonts",
    "escape_text",
    "format_r_factor",
    "format_fixed",
    "format_general",
    "format_int",
]

# --------------------------------------------------------------------------- layout
PAGE_SIZE = A4
MARGIN_SIDE = 18 * mm
MARGIN_TOP = 20 * mm
MARGIN_BOTTOM = 20 * mm
FRAME_WIDTH = PAGE_SIZE[0] - 2 * MARGIN_SIDE
FRAME_HEIGHT = PAGE_SIZE[1] - MARGIN_TOP - MARGIN_BOTTOM
#: Tallest figure, so that the image and its caption always fit on one page.
MAX_FIGURE_HEIGHT = 0.72 * FRAME_HEIGHT
#: Shown for a value that is missing.
DASH = "–"

#: The section headings, in order (the title page carries no heading).
SECTION_HEADINGS: tuple[str, ...] = (
    "1. The data",
    "2. The model",
    "3. The calculation",
    "4. Results",
    "5. Discussion",
    "6. Reproducibility",
    "Appendix A: generated input files",
    "Appendix B: the Python script",
)

#: The run folder (docs/architecture.md §3): path -> content.
FOLDER_LAYOUT: tuple[tuple[str, str], ...] = (
    ("data/", "copies of the input .mint01 files"),
    ("spec.json", "the system specification that was run"),
    (
        "run_<name>.py",
        "the standalone Python script: the input file to check and re-run",
    ),
    (
        "epsr/",
        "the EPSR simulation folder: .mol, .ato, .NWTS.dat, .wts, .EPSR.inp, .pcof and every output file",
    ),
    ("trajectory/", "per-iteration record in the epsrvis folder format"),
    ("manifest.json", "per-run summaries, timings and versions"),
    ("run.log", "console output of the script"),
    ("analysis.json", "the analysis dict (the numbers of section 4)"),
    ("figures/", "PNG figures, summary.csv and summary.md"),
    ("report.pdf", "this report"),
    ("mellor_transcript.jsonl", "every model message and tool call"),
)

INK = colors.HexColor("#1f2430")
GREY = colors.HexColor("#5b6472")
LIGHT_GREY = colors.HexColor("#8a93a0")
GRID = colors.HexColor("#c9ced6")
HEAD_BG = colors.HexColor("#eef1f5")
NOTE_BG = colors.HexColor("#f4f6f9")
LINK = colors.HexColor("#2a5d9f")

# --------------------------------------------------------------------------- fonts
#: role -> (registered name, TTF file, built-in fallback).
_FONT_FILES: dict[str, tuple[str, str, str]] = {
    "sans": ("DejaVuSans", "DejaVuSans.ttf", "Helvetica"),
    "bold": ("DejaVuSans-Bold", "DejaVuSans-Bold.ttf", "Helvetica-Bold"),
    "italic": ("DejaVuSans-Oblique", "DejaVuSans-Oblique.ttf", "Helvetica-Oblique"),
    "bold_italic": (
        "DejaVuSans-BoldOblique",
        "DejaVuSans-BoldOblique.ttf",
        "Helvetica-BoldOblique",
    ),
    "mono": ("DejaVuSansMono", "DejaVuSansMono.ttf", "Courier"),
    "mono_bold": ("DejaVuSansMono-Bold", "DejaVuSansMono-Bold.ttf", "Courier-Bold"),
}
_FONTS: dict[str, str] = {}


def _font_directories() -> list[Path]:
    """Where the DejaVu TTF files are looked for: matplotlib's copy first."""
    directories: list[Path] = []
    with contextlib.suppress(ImportError):
        import matplotlib

        directories.append(Path(matplotlib.get_data_path()) / "fonts" / "ttf")
    directories += [
        Path("/usr/share/fonts/truetype/dejavu"),
        Path("/usr/share/fonts/TTF"),
        Path("/usr/share/fonts/dejavu"),
        Path("/usr/local/share/fonts"),
    ]
    return directories


def register_fonts() -> dict[str, str]:
    """Register the DejaVu fonts with reportlab (once) and return role -> font name.

    Roles: ``sans``, ``bold``, ``italic``, ``bold_italic``, ``mono``, ``mono_bold``.
    When a TTF file cannot be found the built-in Helvetica/Courier is used for that
    role and a ``RuntimeWarning`` says so (those fonts lack ⁻¹ and Greek letters).
    """
    if _FONTS:
        return dict(_FONTS)
    directories = _font_directories()
    registered = set(pdfmetrics.getRegisteredFontNames())
    missing: list[str] = []
    for role, (name, file_name, fallback) in _FONT_FILES.items():
        if name in registered:
            _FONTS[role] = name
            continue
        path = next(
            (
                directory / file_name
                for directory in directories
                if (directory / file_name).is_file()
            ),
            None,
        )
        if path is None:
            _FONTS[role] = fallback
            missing.append(file_name)
            continue
        pdfmetrics.registerFont(TTFont(name, str(path)))
        _FONTS[role] = name
    pdfmetrics.registerFontFamily(
        _FONTS["sans"],
        normal=_FONTS["sans"],
        bold=_FONTS["bold"],
        italic=_FONTS["italic"],
        boldItalic=_FONTS["bold_italic"],
    )
    pdfmetrics.registerFontFamily(
        _FONTS["mono"],
        normal=_FONTS["mono"],
        bold=_FONTS["mono_bold"],
        italic=_FONTS["mono"],
        boldItalic=_FONTS["mono_bold"],
    )
    if missing:
        warnings.warn(
            "DejaVu font files not found ("
            + ", ".join(missing)
            + "): the report uses Helvetica/Courier, so Å⁻¹ and Greek letters may not render",
            RuntimeWarning,
            stacklevel=2,
        )
    return dict(_FONTS)


# --------------------------------------------------------------------------- styles
def make_styles() -> StyleSheet1:
    """The paragraph styles of the report (``Body``, ``H1``..``H3``, ``Bullet``,
    ``Code``, ``Cell``, ...), all in the DejaVu fonts."""
    fonts = register_fonts()
    sans, bold, italic, mono = (
        fonts["sans"],
        fonts["bold"],
        fonts["italic"],
        fonts["mono"],
    )
    styles = StyleSheet1()

    def add(name: str, **settings: Any) -> None:
        styles.add(ParagraphStyle(name, **settings))

    add("Body", fontName=sans, fontSize=9.5, leading=13, spaceAfter=5, textColor=INK)
    add("Small", parent=styles["Body"], fontSize=8, leading=10.5, textColor=GREY)
    add("Meta", parent=styles["Body"], fontSize=10, leading=14)
    add("MetaKey", parent=styles["Meta"], fontName=bold)
    add("Title", fontName=bold, fontSize=22, leading=27, spaceAfter=8, textColor=INK)
    add(
        "Subtitle",
        fontName=sans,
        fontSize=12.5,
        leading=16,
        textColor=GREY,
        spaceAfter=14,
    )
    add(
        "H1",
        fontName=bold,
        fontSize=15,
        leading=19,
        spaceBefore=16,
        spaceAfter=7,
        textColor=INK,
    )
    add(
        "H2",
        fontName=bold,
        fontSize=11.5,
        leading=15,
        spaceBefore=11,
        spaceAfter=4,
        textColor=INK,
    )
    add(
        "H3",
        fontName=bold,
        fontSize=10,
        leading=13,
        spaceBefore=8,
        spaceAfter=3,
        textColor=INK,
    )
    add("MdH1", parent=styles["H2"])
    add("MdH2", parent=styles["H3"])
    add("MdH3", parent=styles["H3"], fontSize=9.5, leading=12.5)
    add(
        "Bullet",
        parent=styles["Body"],
        leftIndent=14,
        bulletIndent=3,
        spaceAfter=2,
        bulletFontName=sans,
        bulletFontSize=9.5,
    )
    add("Quote", parent=styles["Body"], fontName=italic, leftIndent=14, textColor=GREY)
    add(
        "Caption",
        fontName=italic,
        fontSize=8,
        leading=10.5,
        textColor=GREY,
        spaceBefore=3,
        spaceAfter=12,
        alignment=TA_CENTER,
    )
    add(
        "Note",
        parent=styles["Body"],
        backColor=NOTE_BG,
        borderColor=GRID,
        borderWidth=0.5,
        borderPadding=5,
        spaceBefore=8,
        spaceAfter=12,
    )
    add("Code", fontName=mono, fontSize=6.8, leading=8.4, spaceAfter=6, textColor=INK)
    add("CodeBlock", parent=styles["Code"], fontSize=7.5, leading=9.5, leftIndent=6)
    add("Cell", fontName=sans, fontSize=7.8, leading=9.8, textColor=INK)
    add("CellRight", parent=styles["Cell"], alignment=TA_RIGHT)
    add("CellHead", parent=styles["Cell"], fontName=bold)
    add("CellMono", parent=styles["Cell"], fontName=mono, fontSize=6.8, leading=8.6)
    return styles


#: Names of the sample stylesheet accepted in place of the report's own.
_STYLE_FALLBACKS: dict[str, tuple[str, ...]] = {
    "Body": ("BodyText", "Normal"),
    "MdH1": ("Heading2", "Heading1"),
    "MdH2": ("Heading3", "Heading2"),
    "MdH3": ("Heading4", "Heading3"),
    "Bullet": ("BodyText", "Normal"),
    "Quote": ("Italic", "BodyText", "Normal"),
    "Code": ("Code",),
    "CodeBlock": ("Code",),
    "Cell": ("BodyText", "Normal"),
    "CellHead": ("BodyText", "Normal"),
    "CellRight": ("BodyText", "Normal"),
}


def _style(styles: Any, name: str) -> ParagraphStyle:
    """``styles[name]`` from the report stylesheet, a dict or the sample stylesheet."""
    if styles is not None:
        for candidate in (name, *_STYLE_FALLBACKS.get(name, ())):
            try:
                return styles[candidate]
            except (KeyError, TypeError):
                continue
    return make_styles()[name]


# --------------------------------------------------------------------------- text and numbers
def escape_text(text: Any) -> str:
    """``&``, ``<`` and ``>`` escaped for reportlab's paragraph markup."""
    if text is None:
        return ""
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _as_float(value: Any) -> float | None:
    """A finite float, or None for anything that is not a number."""
    if isinstance(value, bool) or value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def format_r_factor(value: Any, significant: int = 3) -> str:
    """``2.27e-4``: scientific notation with ``significant`` digits and a bare exponent."""
    number = _as_float(value)
    if number is None:
        return DASH
    if number == 0.0:
        return "0"
    mantissa, exponent = f"{number:.{max(significant - 1, 0)}e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def format_fixed(value: Any, decimals: int) -> str:
    """Fixed-point with ``decimals`` decimals (energies 1, distances 2)."""
    number = _as_float(value)
    return DASH if number is None else f"{number:.{decimals}f}"


def format_general(value: Any, significant: int = 5) -> str:
    """Up to ``significant`` significant digits, no trailing zeros (force-field numbers)."""
    number = _as_float(value)
    return DASH if number is None else f"{number:.{significant}g}"


def format_int(value: Any) -> str:
    number = _as_float(value)
    return DASH if number is None else f"{int(round(number))}"


def _format_any(value: Any) -> str:
    """A readable rendering of an analysis value of unknown type."""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if value is None:
        return DASH
    if isinstance(value, (int, float)):
        return format_general(value, 4)
    if isinstance(value, Mapping):
        return "; ".join(f"{_pretty_key(k)} {_format_any(v)}" for k, v in value.items())
    if isinstance(value, (list, tuple)):
        return ", ".join(_format_any(item) for item in value)
    return str(value)


def _trend_rows(trends: Mapping[str, Any]) -> list[list[str]]:
    """Rows of the convergence-trend table (empty when there is no trend)."""
    names: dict[str, tuple[str, int | None]] = {
        "energy": ("Energy (kJ/mol)", 1),
        "pressure": ("Pressure (kbar)", 2),
        "ep_amplitude": ("EP amplitude (kJ/mol)", 2),
        "r_factor": ("R-factor", None),
    }
    rows = [["Quantity", "First quarter", "Last quarter", "Change", "Change (%)"]]
    for key, value in trends.items():
        trend = _mapping(value)
        if not trend:
            continue
        label, decimals = names.get(key, (_pretty_key(key), None))
        change = _as_float(trend.get("change_percent"))
        rows.append(
            [
                label,
                _format_trend(trend.get("start"), decimals),
                _format_trend(trend.get("end"), decimals),
                _format_trend(trend.get("change"), decimals),
                f"{change:+.1f}" if change is not None else DASH,
            ]
        )
    return rows if len(rows) > 1 else []


def _format_trend(value: Any, decimals: int | None) -> str:
    return format_r_factor(value) if decimals is None else format_fixed(value, decimals)


def _pretty_key(key: Any) -> str:
    return str(key).replace("_", " ")


def _shorten(text: str, limit: int = 72) -> str:
    return text if len(text) <= limit else text[: limit - 1].rstrip() + "…"


def _mapping(value: Any) -> dict[str, Any]:
    return dict(value) if isinstance(value, Mapping) else {}


def _records(value: Any, key_name: str) -> list[dict[str, Any]]:
    """A list of dicts from a list of dicts or a mapping name -> dict (or number)."""
    if isinstance(value, Mapping):
        records = []
        for key, item in value.items():
            if isinstance(item, Mapping):
                record = dict(item)
                record.setdefault(key_name, key)
            else:
                record = {key_name: key, "value": item}
            records.append(record)
        return records
    if isinstance(value, (list, tuple)):
        return [dict(item) for item in value if isinstance(item, Mapping)]
    return []


def _first(record: Mapping[str, Any], *keys: str, default: Any = None) -> Any:
    """The first of ``keys`` present in ``record`` with a value that is not None."""
    for key in keys:
        if record.get(key) is not None:
            return record[key]
    return default


def _per_dataset(value: Any, labels: Sequence[str]) -> dict[str, Any]:
    """A per-data-set mapping from a dict (label -> value) or a list in label order."""
    if isinstance(value, Mapping):
        return {str(key): item for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        names = list(labels) + [
            f"data set {i + 1}" for i in range(len(labels), len(value))
        ]
        return dict(zip(names[: len(value)], value, strict=True))
    return {}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


#: A ``model_name`` that is this word, or starts with it followed by ``:``, states
#: that no language model wrote the discussion (``"none: <what it is instead>"``).
NO_MODEL_PREFIX = "none"


def _no_model(model_name: str) -> bool:
    lowered = (model_name or "").strip().lower()
    return lowered == NO_MODEL_PREFIX or lowered.startswith(NO_MODEL_PREFIX + ":")


def _describe_source(source: str) -> str:
    """Plain words for ``SystemSpec.source``."""
    source = (source or "").strip()
    if not source:
        return "not recorded"
    if source == "model":
        return "proposed by the language model"
    if source == "user":
        return "given by the user"
    if source.startswith("fallback:library:"):
        key = source.rsplit(":", 1)[-1]
        return (
            f"the built-in library system '{key}', used because the language model "
            "did not produce a valid specification"
        )
    if source.startswith("library:"):
        return f"the built-in library system '{source.split(':', 1)[1]}'"
    return source


def _formula(molecule: Molecule) -> str:
    """``OW, HW×2``: the atom types of a template with their multiplicities."""
    counts: dict[str, int] = {}
    for atom in molecule.atoms:
        counts[atom.label] = counts.get(atom.label, 0) + 1
    return ", ".join(f"{label}×{n}" if n > 1 else label for label, n in counts.items())


# --------------------------------------------------------------------------- flowable helpers
def _paragraph(
    markup: str, style: ParagraphStyle, plain: str, **kwargs: Any
) -> Paragraph:
    """A Paragraph from markup, or from the escaped plain text if the markup is malformed."""
    try:
        return Paragraph(markup, style, **kwargs)
    except Exception:  # reportlab's parser: show the text as written rather than fail
        return Paragraph(escape_text(plain), style, **kwargs)


def _cell(value: Any, styles: Any, style_name: str) -> Flowable:
    if isinstance(value, Flowable):
        return value
    text = DASH if value is None else str(value)
    return Paragraph(escape_text(text), _style(styles, style_name))


def _cell_text(value: Any) -> str:
    if isinstance(value, Paragraph):
        return value.getPlainText()
    if isinstance(value, Flowable):
        return ""
    return "" if value is None else str(value)


#: Cell padding (left + right) plus a little slack added to every measured width,
#: so that a word measured at exactly the cell width does not wrap on rounding.
_CELL_PADDING = 9.0
#: No column is wider than this fraction of the frame; longer text wraps.
_MAX_COLUMN_FRACTION = 0.45
#: The widest columns of an over-wide table are capped no lower than this.
_SHRINK_FLOOR = 60.0


def _text_width(text: str, style: ParagraphStyle) -> float:
    return pdfmetrics.stringWidth(text, style.fontName, style.fontSize)


def _auto_widths(
    cells: list[list[Flowable]], total: float, header: bool, fallback: ParagraphStyle
) -> list[float]:
    """Column widths from the font metrics: body cells get the width of their text
    (so numbers never wrap), a header cell the width of its longest word (headers
    wrap by word), long prose is capped at :data:`_MAX_COLUMN_FRACTION` of the
    frame; spare width goes to the headers first, then to every column; an
    over-wide table caps its widest (prose) columns first (:func:`_cap_widest`)."""
    n_columns = max(len(row) for row in cells)
    cap = _MAX_COLUMN_FRACTION * total
    needs: list[float] = []
    wants: list[float] = []
    for column in range(n_columns):
        need = want = 0.0
        for index, row in enumerate(cells):
            if column >= len(row):
                continue
            value = row[column]
            style = value.style if isinstance(value, Paragraph) else fallback
            content = _cell_text(value)
            if header and index == 0:
                words = content.split() or [""]
                need = max(need, max(_text_width(word, style) for word in words))
                want = max(want, _text_width(content, style))
            else:
                need = max(need, _text_width(content, style))
        needs.append(min(need + _CELL_PADDING, cap))
        wants.append(min(max(want + _CELL_PADDING, needs[-1]), cap))
    total_need, total_want = sum(needs), sum(wants)
    if total_want <= total:
        extra = total - total_want
        return [want + extra * want / total_want for want in wants]
    if total_need <= total:
        extra = total - total_need
        growth = [want - need for want, need in zip(wants, needs, strict=True)]
        pool = sum(growth) or 1.0
        return [
            need + extra * grow / pool for need, grow in zip(needs, growth, strict=True)
        ]
    return _cap_widest(needs, total)


def _cap_widest(needs: list[float], total: float) -> list[float]:
    """Widths summing to ``total``: the widest columns are capped at one common
    level (their text wraps) and the narrower ones keep their width; when even a
    cap of :data:`_SHRINK_FLOOR` is not enough, everything shrinks proportionally."""
    order = sorted(needs, reverse=True)
    for k in range(1, len(order) + 1):
        level = (total - sum(order[k:])) / k
        next_need = order[k] if k < len(order) else 0.0
        if level >= next_need:
            if level < _SHRINK_FLOOR:
                break
            return [min(need, level) for need in needs]
    scale = total / (sum(needs) or 1.0)
    return [need * scale for need in needs]


def _table(
    rows: list[list[Any]],
    styles: Any,
    *,
    header: bool = True,
    widths: Sequence[float] | None = None,
    align: Sequence[str] | None = None,
) -> Table:
    """A table with a light grid. Cells are strings (escaped) or ready flowables;
    ``widths`` are fractions of the frame width; ``align`` gives ``"right"`` per column.
    """
    if not rows:
        raise ValueError("a table needs at least one row")
    n_columns = max(len(row) for row in rows)
    padded = [list(row) + [""] * (n_columns - len(row)) for row in rows]
    data: list[list[Flowable]] = []
    for index, row in enumerate(padded):
        cells = []
        for column, value in enumerate(row):
            if header and index == 0:
                style_name = "CellHead"
            elif align and column < len(align) and align[column] == "right":
                style_name = "CellRight"
            else:
                style_name = "Cell"
            cells.append(_cell(value, styles, style_name))
        data.append(cells)
    if widths:
        scale = FRAME_WIDTH / sum(widths)
        column_widths = [width * scale for width in widths]
    else:
        column_widths = _auto_widths(data, FRAME_WIDTH, header, _style(styles, "Cell"))
    table = Table(
        data, colWidths=column_widths, repeatRows=1 if header else 0, hAlign="LEFT"
    )
    commands: list[tuple[Any, ...]] = [
        ("GRID", (0, 0), (-1, -1), 0.4, GRID),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ]
    if header:
        commands.append(("BACKGROUND", (0, 0), (-1, 0), HEAD_BG))
    table.setStyle(TableStyle(commands))
    return table


def _code_columns(style: ParagraphStyle) -> int:
    """How many characters of ``style``'s mono font fit on one line of the frame."""
    width = pdfmetrics.stringWidth("M", style.fontName, style.fontSize) or 1.0
    usable = (
        FRAME_WIDTH - getattr(style, "leftIndent", 0) - getattr(style, "rightIndent", 0)
    )
    return max(40, int(usable / width) - 1)


def _preformatted(text: str, style: ParagraphStyle) -> Preformatted:
    """Verbatim text in the mono font; long lines are wrapped (continuations start with ↳)."""
    text = text.replace("\r\n", "\n").replace("\r", "\n").expandtabs(4)
    return Preformatted(
        text, style, maxLineLength=_code_columns(style), newLineChars="↳ "
    )


def _scaled_image(path: Path | None) -> Image | None:
    """The PNG/JPEG at ``path`` scaled to the frame width (aspect ratio kept, height
    capped), or None when the file is missing or not an image."""
    if path is None or not path.is_file():
        return None
    if path.suffix.lower() not in {".png", ".jpg", ".jpeg"}:
        return None
    try:
        width, height = ImageReader(str(path)).getSize()
    except Exception:  # unreadable or truncated image: the figure is skipped
        return None
    if width <= 0 or height <= 0:
        return None
    scale = min(FRAME_WIDTH / width, MAX_FIGURE_HEIGHT / height)
    return Image(str(path), width=width * scale, height=height * scale)


class _Footer:
    """Draws "Mellor report — <title> — page N" at the foot of every page."""

    def __init__(self, title: str, font: str) -> None:
        self.title = title
        self.font = font

    def __call__(self, canvas: Any, document: Any) -> None:
        canvas.saveState()
        canvas.setFont(self.font, 7.5)
        canvas.setFillColor(GREY)
        canvas.drawCentredString(
            PAGE_SIZE[0] / 2,
            11 * mm,
            f"Mellor report — {self.title} — page {document.page}",
        )
        canvas.restoreState()


# --------------------------------------------------------------------------- the report
class _ReportBuilder:
    """Builds the story (list of flowables) section by section."""

    def __init__(
        self,
        result: RunResult,
        analysis: Mapping[str, Any] | None,
        figures: Mapping[str, Any] | None,
        discussion: str,
        model_name: str,
        data_descriptions: Sequence[Mapping[str, Any]] | None,
    ) -> None:
        self.result = result
        self.spec: SystemSpec = result.spec
        self.manifest = _mapping(result.manifest)
        self.analysis = _mapping(analysis)
        self.figures: dict[str, Path] = {
            str(name): Path(path) for name, path in (figures or {}).items() if path
        }
        self.discussion = discussion or ""
        self.model_name = (model_name or "").strip()
        self.data_descriptions = [
            dict(item)
            for item in (data_descriptions or [])
            if isinstance(item, Mapping)
        ]
        self.fonts = register_fonts()
        self.styles = make_styles()
        self.story: list[Flowable] = []
        self.figure_number = 0
        self.title = self.spec.title or self.spec.name
        self.quick = bool(result.quick or self.manifest.get("quick"))
        self.generated: EpsrInputs | None = None
        self.generate_error = ""
        try:
            self.generated = generate(self.spec)
        except Exception as error:  # the report is still written; the error is shown
            self.generate_error = f"{type(error).__name__}: {error}"

    # ---- small helpers
    def heading(self, text: str) -> None:
        self.story.append(Paragraph(escape_text(text), self.styles["H1"]))

    def subheading(self, text: str, style: str = "H2") -> None:
        self.story.append(Paragraph(escape_text(text), self.styles[style]))

    def plain(self, text: str, style: str = "Body") -> None:
        self.story.append(Paragraph(escape_text(text), self.styles[style]))

    def markup(self, text: str, style: str = "Body") -> None:
        self.story.append(_paragraph(text, self.styles[style], text))

    def bullets(self, items: Sequence[str], numbered: bool = False) -> None:
        for number, item in enumerate(items, start=1):
            bullet = f"{number}." if numbered else "•"
            self.story.append(
                Paragraph(escape_text(item), self.styles["Bullet"], bulletText=bullet)
            )

    def table(self, rows: list[list[Any]], **settings: Any) -> None:
        self.story.append(_table(rows, self.styles, **settings))

    def code(self, text: str, style: str = "Code") -> None:
        self.story.append(_preformatted(text, self.styles[style]))

    def figure(self, name: str) -> bool:
        """Add the figure ``name`` (image + caption); False when it is missing."""
        path = self.figures.get(name)
        image = _scaled_image(path)
        if image is None or path is None:
            return False
        self.figure_number += 1
        caption = FIGURE_CAPTIONS.get(name) or name.replace("_", " ")
        text = (
            f"<b>Figure {self.figure_number}.</b> {escape_text(caption)} "
            f"<font color='#8a93a0'>[{escape_text(path.name)}]</font>"
        )
        self.story.append(
            KeepTogether([Spacer(1, 4), image, Paragraph(text, self.styles["Caption"])])
        )
        return True

    def dataset_labels(self) -> list[str]:
        labels = [dataset.label for dataset in self.spec.datasets]
        if labels:
            return labels
        return [
            str(item.get("label"))
            for item in self.manifest.get("datasets") or []
            if isinstance(item, Mapping) and item.get("label")
        ]

    def locate_data_file(self, dataset: Dataset) -> Path | None:
        candidates = [
            self.result.data_folder / dataset.file_name,
            Path(dataset.file),
            Path(self.result.folder) / dataset.file,
        ]
        return next((path for path in candidates if path.is_file()), None)

    def run_records(self) -> list[dict[str, Any]]:
        """One flat dict per EPSR run: the analysis ``runs`` entries, or the manifest
        runs with their ``summary`` merged in."""
        runs = self.analysis.get("runs")
        source = runs if isinstance(runs, (list, tuple)) and runs else self.result.runs
        records = []
        for index, run in enumerate(source, start=1):
            if not isinstance(run, Mapping):
                continue
            summary = run.get("summary")
            record = dict(summary) if isinstance(summary, Mapping) else {}
            record.update(
                {key: value for key, value in run.items() if key != "summary"}
            )
            record.setdefault("run", index)
            records.append(record)
        return records

    def schedule_rows(self) -> list[dict[str, Any]]:
        schedule = _mapping(self.analysis.get("schedule"))
        for candidate in (schedule.get("runs"), self.manifest.get("schedule")):
            rows = [dict(item) for item in candidate or [] if isinstance(item, Mapping)]
            if rows:
                return rows
        return [
            {"stage": stage, "niter": niter, "ereq": ereq}
            for stage, niter, ereq in self.spec.schedule.iterations(self.quick)
        ]

    # ---- sections
    def build(self) -> list[Flowable]:
        self.title_page()
        self.data_section()
        self.model_section()
        self.calculation_section()
        self.results_section()
        self.discussion_section()
        self.reproducibility_section()
        self.appendix_inputs()
        self.appendix_script()
        return self.story

    def title_page(self) -> None:
        styles = self.styles
        story = self.story
        story.append(Spacer(1, 28 * mm))
        story.append(Paragraph(escape_text(self.title), styles["Title"]))
        story.append(Paragraph("EPSR refinement report", styles["Subtitle"]))
        if _no_model(self.model_name):
            prepared = "Mellor, without a language model"
        elif self.model_name:
            prepared = f"Mellor, with the language model {self.model_name}"
        else:
            prepared = "Mellor (no language model named)"
        today = dt.date.today().strftime("%d %B %Y")
        rows: list[list[Any]] = [
            ["Prepared by", prepared],
            ["Date", today],
            ["Run name", self.spec.name],
            ["Run folder", str(Path(self.result.folder).resolve())],
            ["Specification", _describe_source(self.spec.source)],
        ]
        created = self.manifest.get("created")
        if created:
            rows.append(["Run started", str(created)])
        files = self.data_file_lines()
        rows.append(
            [
                "Data files",
                Paragraph(
                    "<br/>".join(escape_text(line) for line in files), styles["Meta"]
                ),
            ]
        )
        cells = [
            [
                Paragraph(escape_text(key), styles["MetaKey"]),
                (
                    value
                    if isinstance(value, Flowable)
                    else Paragraph(escape_text(value), styles["Meta"])
                ),
            ]
            for key, value in rows
        ]
        meta = Table(
            cells, colWidths=[0.2 * FRAME_WIDTH, 0.8 * FRAME_WIDTH], hAlign="LEFT"
        )
        meta.setStyle(
            TableStyle(
                [
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("TOPPADDING", (0, 0), (-1, -1), 3),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                    ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ]
            )
        )
        story.append(meta)
        story.append(Spacer(1, 8 * mm))
        key_rows = self.key_results()
        if key_rows:
            story.append(Paragraph("Key results", styles["H2"]))
            self.table([["Quantity", "Value"]] + key_rows, widths=[0.45, 0.55])
            story.append(Spacer(1, 6 * mm))
        story.append(
            Paragraph(
                "<b>How to read this report.</b> Every table, number and figure was computed "
                "by EPSR and the analysis code from the files in the run folder. The "
                "language model chose the model system and wrote the discussion (section 5), "
                "which is labelled as its own text. The Python script in Appendix B is the "
                "complete input: running it reproduces the calculation without Mellor.",
                styles["Note"],
            )
        )
        story.append(PageBreak())

    def data_file_lines(self) -> list[str]:
        lines = []
        rows = self.data_rows()
        for row in rows:
            line = row["file"]
            if row["label"] != DASH:
                line = f"{row['label']}: {line}"
            if row["title"] != DASH:
                line += f" — {row['title']}"
            lines.append(line)
        return lines or ["none"]

    def key_results(self) -> list[list[str]]:
        records = self.run_records()
        if not records:
            return []
        first, last = records[0], records[-1]
        rows = []
        if (
            _as_float(first.get("r_factor")) is not None
            and _as_float(last.get("r_factor")) is not None
        ):
            rows.append(
                [
                    "R-factor, first run → last run",
                    f"{format_r_factor(first.get('r_factor'))} → {format_r_factor(last.get('r_factor'))}",
                ]
            )
        total = sum(int(_as_float(record.get("niter")) or 0) for record in records)
        rows.append(["EPSR runs / iterations", f"{len(records)} / {total}"])
        if _as_float(last.get("energy")) is not None:
            rows.append(["Final energy (kJ/mol)", format_fixed(last.get("energy"), 1)])
        if _as_float(last.get("pressure")) is not None:
            rows.append(
                ["Final pressure (kbar)", format_fixed(last.get("pressure"), 2)]
            )
        if _as_float(last.get("ep_amplitude")) is not None:
            text = format_fixed(last.get("ep_amplitude"), 2)
            if _as_float(last.get("ereq")) is not None:
                text += f" (requested {format_fixed(last.get('ereq'), 2)})"
            rows.append(["Final EP amplitude (kJ/mol)", text])
        return rows

    # ---- 1. the data
    def data_rows(self) -> list[dict[str, str]]:
        descriptions = self.data_descriptions or self.describe_data_files()
        labels = {dataset.file_name: dataset.label for dataset in self.spec.datasets}
        rows = []
        for item in descriptions:
            file_name = Path(str(item.get("file") or item.get("path") or "")).name
            q_min, q_max = _as_float(item.get("q_min")), _as_float(item.get("q_max"))
            q_range = (
                f"{q_min:.2f} – {q_max:.2f}"
                if q_min is not None and q_max is not None
                else DASH
            )
            rows.append(
                {
                    "label": labels.get(file_name, DASH),
                    "file": file_name or DASH,
                    "title": str(item.get("title") or DASH),
                    "instrument": str(item.get("instrument") or DASH),
                    "run": str(item.get("run_number") or DASH),
                    "q_range": q_range,
                    "points": format_int(item.get("n_points")),
                }
            )
        return rows

    def describe_data_files(self) -> list[dict[str, Any]]:
        """``data.describe`` of every data file found in the run folder."""
        from .data import describe, read_mint01

        descriptions: list[dict[str, Any]] = []
        for dataset in self.spec.datasets:
            path = self.locate_data_file(dataset)
            if path is None:
                descriptions.append(
                    {"file": dataset.file_name, "title": "(file not found)"}
                )
                continue
            try:
                descriptions.append(describe(read_mint01(path)))
            except (ValueError, OSError) as error:
                descriptions.append(
                    {"file": dataset.file_name, "title": f"(not readable: {error})"}
                )
        return descriptions

    def data_section(self) -> None:
        self.heading(SECTION_HEADINGS[0])
        rows = self.data_rows()
        self.plain(
            "The measured structure factors (Gudrun .mint01 files: Q, F(Q), error) given "
            "to Mellor. The data set label is the name used for the sample in the rest of "
            "the report and in the EPSR files."
        )
        if rows:
            header = [
                "Data set",
                "File",
                "Title",
                "Instrument",
                "Run",
                "Q range (Å⁻¹)",
                "Points",
            ]
            body = [
                [
                    r["label"],
                    r["file"],
                    r["title"],
                    r["instrument"],
                    r["run"],
                    r["q_range"],
                    r["points"],
                ]
                for r in rows
            ]
            self.table(
                [header] + body,
                align=["left", "left", "left", "left", "right", "right", "right"],
            )
        else:
            self.plain("No data files are recorded for this run.", "Small")
        self.figure("00_data")

    # ---- 2. the model
    def model_section(self) -> None:
        spec = self.spec
        styles = self.styles
        self.heading(SECTION_HEADINGS[1])
        self.subheading("Composition")
        if spec.molecules:
            rows: list[list[Any]] = [
                [
                    "Molecule",
                    "Atoms",
                    "Atom types",
                    "Molecules in the box",
                    "Atoms in the box",
                ]
            ]
            for molecule in spec.molecules:
                rows.append(
                    [
                        molecule.name,
                        str(molecule.n_atoms),
                        _formula(molecule),
                        str(molecule.count),
                        str(molecule.count * molecule.n_atoms),
                    ]
                )
            self.table(rows, align=["left", "right", "left", "right", "right"])
        if spec.kind == "pdb" and spec.pdb is not None:
            self.plain(
                f"The box is read from the PDB file {spec.pdb.file_name} (readpdb, van der "
                f"Waals factor {spec.pdb.vdw_factor:g}); the crystal positions are kept."
            )
        elif not spec.molecules:
            self.plain("The specification lists no molecules.", "Small")

        self.subheading("Box")
        system = _mapping(self.analysis.get("system"))
        engine: dict[str, Any] = dict(self.manifest)
        engine.update(
            {key: value for key, value in system.items() if value is not None}
        )
        spec_box = (
            format_fixed(spec.box_side, 2)
            if spec.density > 0 and spec.n_atoms > 0
            else DASH
        )
        try:
            spec_mass = format_fixed(spec.mass_density(), 3)
        except (KeyError, ZeroDivisionError):
            spec_mass = DASH
        rows = [
            ["Quantity", "Specification", "As read by the engine"],
            ["Kind", spec.kind, str(engine.get("kind") or DASH)],
            ["Molecules", str(spec.n_molecules), format_int(engine.get("n_molecules"))],
            ["Atoms", str(spec.n_atoms), format_int(engine.get("n_atoms"))],
            ["Box side (Å)", spec_box, format_fixed(engine.get("box_side"), 2)],
            [
                "Atomic number density (atoms/Å³)",
                format_general(spec.density, 4),
                format_general(engine.get("density"), 4),
            ],
            [
                "Mass density (g/cm³)",
                spec_mass,
                format_fixed(engine.get("mass_density_g_cm3"), 3),
            ],
            [
                "Temperature (K)",
                format_fixed(spec.temperature, 1),
                format_fixed(engine.get("temperature"), 1),
            ],
        ]
        self.table(rows, align=["left", "right", "right"])

        self.subheading("Atom types and reference potential")
        counts = spec.type_counts()
        rows = [
            [
                "Type",
                "Element",
                "ε (kJ/mol)",
                "σ (Å)",
                "Charge (e)",
                "Mass (amu)",
                "Atoms in the box",
            ]
        ]
        for atom_type in spec.atom_types:
            rows.append(
                [
                    atom_type.label,
                    atom_type.element,
                    format_general(atom_type.epsilon),
                    format_general(atom_type.sigma),
                    format_general(atom_type.charge),
                    format_general(atom_type.mass),
                    str(counts.get(atom_type.label, 0)),
                ]
            )
        self.table(
            rows, align=["left", "left", "right", "right", "right", "right", "right"]
        )
        self.plain(
            "Lennard-Jones ε and σ with the partial charges define the reference "
            "potential; EPSR refines an empirical potential on top of it.",
            "Small",
        )

        self.subheading("Isotopic composition of the samples")
        rows = [["Data set", "File", "Composition", "Exchangeable H"]]
        for dataset in spec.datasets:
            rows.append(
                [
                    dataset.label,
                    dataset.file_name,
                    _paragraph(
                        self.composition_markup(dataset), styles["Cell"], dataset.label
                    ),
                    ", ".join(dataset.exchangeable) or "none",
                ]
            )
        if len(rows) > 1:
            self.table(rows)
        else:
            self.plain("No data sets are recorded in the specification.", "Small")

        self.subheading("Provenance and notes")
        self.plain(f"Specification source: {_describe_source(spec.source)}.")
        for line in spec.notes.splitlines():
            if line.strip():
                self.plain(line.strip())
        _errors, spec_warnings = spec.problems()
        if spec_warnings:
            self.plain("Specification warnings:", "Small")
            self.bullets(spec_warnings)

    def composition_markup(self, dataset: Dataset) -> str:
        spec = self.spec
        parts = []
        natural = []
        for label in spec.labels:
            isotopes = dataset.composition(label)
            try:
                element = spec.atom_type(label).element
            except KeyError:
                element = "?"
            if len(isotopes) == 1 and isotopes[0].mass_number == 0:
                natural.append(label)
                continue
            terms = []
            for isotope in isotopes:
                percent = f"{isotope.fraction * 100:g} %"
                if isotope.mass_number:
                    terms.append(
                        f"<super>{isotope.mass_number}</super>{escape_text(element)} {percent}"
                    )
                else:
                    terms.append(f"natural {escape_text(element)} {percent}")
            parts.append(f"<b>{escape_text(label)}</b>: " + " + ".join(terms))
        if not parts:
            return "natural abundance for every atom type"
        if natural:
            parts.append(
                "other types natural ("
                + ", ".join(escape_text(x) for x in natural)
                + ")"
            )
        return "; ".join(parts)

    # ---- 3. the calculation
    def calculation_section(self) -> None:
        spec = self.spec
        self.heading(SECTION_HEADINGS[2])
        self.subheading("Setup")
        schedule_info = _mapping(self.analysis.get("schedule"))
        steps = [
            str(step)
            for step in (
                self.manifest.get("setup_steps")
                or schedule_info.get("setup_steps")
                or []
            )
        ]
        if not steps and self.generated is not None:
            steps = list(self.generated.setup_steps)
        self.plain(
            "The simulation box was built with the EPSR setup programs, in this order "
            "(the script runs each of them through the epsr wrapper):"
        )
        if steps:
            self.bullets(steps, numbered=True)
        else:
            self.plain("The setup steps are not recorded.", "Small")

        self.subheading("EPSR schedule")
        schedule_rows = self.schedule_rows()
        rows: list[list[Any]] = [["Run", "Stage", "Iterations", "ereq (kJ/mol)"]]
        total = 0
        for index, item in enumerate(schedule_rows, start=1):
            niter = _as_float(item.get("niter"))
            total += int(niter or 0)
            rows.append(
                [
                    str(item.get("run") or index),
                    str(item.get("stage") or DASH),
                    format_int(item.get("niter")),
                    format_fixed(item.get("ereq"), 2),
                ]
            )
        self.table(rows, align=["right", "left", "right", "right"])
        schedule = spec.schedule
        summary = (
            f"{total} EPSR iterations in total; ntimes {schedule.ntimes} Monte Carlo cycles "
            f"per iteration; requested empirical-potential amplitude ereq "
            f"{format_fixed(schedule.ereq, 2)} kJ/mol"
        )
        if spec.kind == "molecular":
            summary += (
                f"; fmole {schedule.fmole_iterations} geometry-relaxation iterations"
            )
        self.plain(summary + ".")
        if self.quick:
            self.plain(
                "The quick schedule (--quick) was used: fewer iterations than the "
                "specification's full schedule, for a preview rather than a converged result.",
                "Note",
            )

        self.subheading(f"{spec.inp_name} keywords")
        if spec.epsr_keywords:
            rows = [["Keyword", "Value"]]
            for keyword, value in spec.epsr_keywords.items():
                rows.append([str(keyword), _format_any(value)])
            self.table(rows, widths=[0.3, 0.7])
            self.plain(
                "Every other keyword keeps the template default "
                "(validation/water/waterbox.EPSR.inp, the EPSR test-case values).",
                "Small",
            )
        else:
            self.plain(
                "No global keyword differs from the template defaults "
                "(validation/water/waterbox.EPSR.inp).",
                "Small",
            )
        dataset_rows = self.dataset_keyword_rows()
        if dataset_rows:
            self.plain("Per-data-set keywords of the data blocks:", "Small")
            self.table(
                dataset_rows, align=["left"] + ["right"] * (len(dataset_rows[0]) - 1)
            )

        self.subheading("Software")
        versions = _mapping(
            self.manifest.get("versions")
            or _mapping(self.analysis.get("system")).get("versions")
        )
        rows = [["Component", "Version"]]
        rows.append(
            [
                "epsr (Python wrapper of the Rust EPSR engine)",
                str(versions.get("epsr") or self.installed_version("epsr")),
            ]
        )
        rows.append(
            [
                "epsrvis (visualisation toolkit)",
                str(versions.get("epsrvis") or self.installed_version("epsrvis")),
            ]
        )
        rows.append(["mellor", str(versions.get("mellor") or __version__)])
        rows.append(["Python", str(self.manifest.get("python") or DASH)])
        if self.manifest.get("script"):
            rows.append(["Script", str(self.manifest.get("script"))])
        self.table(rows, widths=[0.6, 0.4])

    def dataset_keyword_rows(self) -> list[list[str]]:
        keyword_sets: list[tuple[str, dict[str, Any]]] = []
        for dataset in self.spec.datasets:
            try:
                keywords = dict(dataset_keywords(dataset, self.spec))
            except (
                Exception
            ):  # fall back to the dataclass fields as written in the spec
                keywords = {
                    "nrtype": dataset.nrtype,
                    "rshmin": dataset.rshmin,
                    "tweak": dataset.tweak,
                    "efilereq": dataset.efilereq,
                    "szeros": dataset.szeros,
                    "porodadd": dataset.porodadd,
                }
            keyword_sets.append((dataset.label, keywords))
        if not keyword_sets:
            return []
        columns: list[str] = []
        for _label, keywords in keyword_sets:
            for key in keywords:
                if key not in columns:
                    columns.append(key)
        rows = [["Data set"] + columns]
        for label, keywords in keyword_sets:
            rows.append(
                [label] + [_format_any(keywords.get(column)) for column in columns]
            )
        return rows

    @staticmethod
    def installed_version(package: str) -> str:
        from importlib import metadata

        try:
            return metadata.version(package)
        except metadata.PackageNotFoundError:
            return "unknown"

    # ---- 4. results
    def results_section(self) -> None:
        self.heading(SECTION_HEADINGS[3])
        self.runs_table()
        self.fit_table()
        self.convergence_summary()
        self.shells_table()
        flags = self.analysis.get("flags")
        if isinstance(flags, (list, tuple)) and flags:
            self.subheading("Warnings")
            self.bullets([str(flag) for flag in flags])
        names = [name for name in FIGURE_CAPTIONS if name != "00_data"]
        names += sorted(name for name in self.figures if name not in FIGURE_CAPTIONS)
        start = len(self.story)
        added = [self.figure(name) for name in names]
        if any(added):
            self.story.insert(start, Paragraph("Figures", self.styles["H2"]))
        else:
            self.plain("No figures are available for this run.", "Small")

    def runs_table(self) -> None:
        records = self.run_records()
        if not records:
            self.plain("No EPSR runs are recorded for this folder.", "Note")
            return
        self.subheading("EPSR runs")
        labels = self.dataset_labels()
        per_dataset = [
            _per_dataset(r.get("r_factor_per_dataset"), labels) for r in records
        ]
        columns = [label for label in labels if any(label in p for p in per_dataset)]
        for mapping in per_dataset:
            columns += [key for key in mapping if key not in columns]
        header = ["Run", "Stage", "Iter.", "ereq", "R-factor"]
        header += [f"R ({column})" for column in columns]
        header += [
            "Energy (kJ/mol)",
            "Pressure (kbar)",
            "EP ampl. (kJ/mol)",
            "Time (s)",
        ]
        rows: list[list[Any]] = [header]
        for record, mapping in zip(records, per_dataset, strict=True):
            row = [
                str(record.get("run", DASH)),
                str(record.get("stage") or DASH),
                format_int(record.get("niter")),
                format_fixed(record.get("ereq"), 2),
                format_r_factor(record.get("r_factor")),
            ]
            row += [format_r_factor(mapping.get(column)) for column in columns]
            row += [
                format_fixed(record.get("energy"), 1),
                format_fixed(record.get("pressure"), 2),
                format_fixed(record.get("ep_amplitude"), 2),
                format_fixed(record.get("duration_s"), 1),
            ]
            rows.append(row)
        self.table(rows, align=["right", "left"] + ["right"] * (len(header) - 2))
        self.plain(
            "Values at the last iteration of each run, from the engine's summary "
            "(manifest.json). The R-factor is EPSR's measure of the fit to the data "
            "(lower is better); the EP amplitude is the magnitude of the empirical "
            "potential against the requested ereq.",
            "Small",
        )

    def fit_table(self) -> None:
        fit = _records(self.analysis.get("fit"), "dataset")
        if not fit:
            return
        self.subheading("Fit quality per data set")
        has_reference = any(
            _as_float(record.get("reference_r_factor")) is not None for record in fit
        )
        header = ["Data set", "Initial R-factor"]
        if has_reference:
            header.append("After the reference potential")
        header += ["Final R-factor", "Improvement", "Final R²", "Fitted Q range (Å⁻¹)"]
        rows: list[list[Any]] = [header]
        for record in fit:
            initial = _as_float(
                _first(record, "initial_r_factor", "r_factor_initial", "initial")
            )
            final = _as_float(
                _first(record, "final_r_factor", "r_factor_final", "final")
            )
            improvement = _as_float(_first(record, "improvement_factor", "improvement"))
            if improvement is None and initial and final:
                improvement = initial / final
            q_min, q_max = _as_float(record.get("q_min")), _as_float(
                record.get("q_max")
            )
            q_range = _first(record, "q_range", "fitted_q_range")
            if isinstance(q_range, (list, tuple)) and len(q_range) == 2:
                q_min, q_max = _as_float(q_range[0]), _as_float(q_range[1])
            if q_min is not None and q_max is not None:
                q_text = f"{format_fixed(q_min, 2)} – {format_fixed(q_max, 2)}"
            else:
                q_text = DASH
            row = [
                str(_first(record, "dataset", "label", "name", default=DASH)),
                format_r_factor(initial),
            ]
            if has_reference:
                row.append(format_r_factor(record.get("reference_r_factor")))
            row += [
                format_r_factor(final),
                f"{improvement:.1f}×" if improvement is not None else DASH,
                format_fixed(_first(record, "r2", "r2_final", "final_r2"), 3),
                q_text,
            ]
            rows.append(row)
        self.table(rows, align=["left"] + ["right"] * (len(header) - 1))
        note = "R-factor of each data set at the end of the first run (initial)"
        if has_reference:
            note += ", after the equilibration on the reference potential alone"
        note += (
            " and at the end of the last run; R² of the final fit over the fitted "
            "Q range (from the recorded trajectory)."
        )
        self.plain(note, "Small")

    def convergence_summary(self) -> None:
        convergence = _mapping(self.analysis.get("convergence"))
        if not convergence:
            return
        self.subheading("Convergence of the last run")
        verdict = _first(convergence, "verdict", "summary", "text")
        if verdict:
            self.plain(str(verdict))
        rows: list[list[Any]] = [["Quantity", "Value"]]
        shown = {"verdict", "summary", "text", "trends"}
        last_run = convergence.get("last_run")
        if last_run is not None:
            stage = convergence.get("last_run_stage")
            rows.append(["Last run", f"{last_run}" + (f" ({stage})" if stage else "")])
            shown |= {"last_run", "last_run_stage"}
        ep = _as_float(convergence.get("ep_amplitude"))
        ereq = _as_float(convergence.get("ereq"))
        if ep is not None or ereq is not None:
            text = (
                f"{format_fixed(ep, 2)} of the requested {format_fixed(ereq, 2)} kJ/mol"
            )
            fraction = _as_float(convergence.get("ep_fraction_of_ereq"))
            if fraction is not None:
                text += f" ({fraction:.0%})"
            rows.append(["EP amplitude", text])
            shown |= {"ep_amplitude", "ereq", "ep_fraction_of_ereq"}
        if "ep_at_limit" in convergence:
            text = "yes" if convergence.get("ep_at_limit") else "no"
            after = _as_float(convergence.get("ep_limit_reached_after"))
            if after is not None:
                text += f", reached after {int(after)} iteration(s) with the EP on"
            rows.append(["EP amplitude at the ereq limit", text])
            shown |= {
                "ep_at_limit",
                "ep_limit_reached_after",
                "ep_limit_reached_at_iteration",
                "ep_on_from_iteration",
            }
        r_parts = []
        for key, name in (
            ("r_factor_initial", "initial"),
            ("r_factor_before_last_run", "before the last run"),
            ("r_factor", "final"),
        ):
            if _as_float(convergence.get(key)) is not None:
                r_parts.append(f"{format_r_factor(convergence.get(key))} ({name})")
        if r_parts:
            text = " → ".join(r_parts)
            change = _as_float(convergence.get("r_factor_change_percent"))
            if change is not None:
                text += f"; change over the last run {change:+.1f} %"
            if convergence.get("r_factor_rising"):
                text += "; rising"
            rows.append(["R-factor", text])
            shown |= {
                "r_factor_initial",
                "r_factor_before_last_run",
                "r_factor",
                "r_factor_change_percent",
                "r_factor_rising",
            }
        if _as_float(convergence.get("energy")) is not None:
            text = f"{format_fixed(convergence.get('energy'), 1)} kJ/mol"
            drift = _as_float(convergence.get("energy_drift_percent"))
            if drift is not None:
                text += f"; drift over the last run {drift:+.1f} %"
            rows.append(["Energy", text])
            shown |= {"energy", "energy_drift_percent"}
        if _as_float(convergence.get("pressure")) is not None:
            rows.append(
                ["Pressure", f"{format_fixed(convergence.get('pressure'), 2)} kbar"]
            )
            shown.add("pressure")
        iterations = _as_float(convergence.get("last_run_iterations"))
        if iterations is not None:
            text = f"{int(iterations)}"
            quarter = _as_float(convergence.get("quarter_iterations"))
            if quarter is not None:
                text += f" (the trends compare the first and the last {int(quarter)})"
            rows.append(["Iterations in the last run", text])
            shown |= {
                "last_run_iterations",
                "quarter_iterations",
                "last_run_first_iteration",
            }
        for key, value in convergence.items():
            if key not in shown:
                rows.append([_pretty_key(key), _format_any(value)])
        if len(rows) > 1:
            self.table(rows, widths=[0.35, 0.65])
        trend_rows = _trend_rows(_mapping(convergence.get("trends")))
        if trend_rows:
            self.plain(
                "Trends over the last run: the mean of its first quarter against the "
                "mean of its last quarter.",
                "Small",
            )
            self.table(trend_rows, align=["left"] + ["right"] * 4)

    def shells_table(self) -> None:
        shells = _records(self.analysis.get("shells"), "pair")
        if not shells:
            return
        self.subheading("First coordination shells")
        has_reverse = any(
            _as_float(record.get("n_coord_reverse")) is not None for record in shells
        )
        header = [
            "Pair",
            "First peak r (Å)",
            "g(r) at the peak",
            "First minimum (Å)",
            "Coordination number",
        ]
        if has_reverse:
            header.append("Reverse coordination")
        header.append("Partner density (atoms/Å³)")
        rows: list[list[Any]] = [header]
        for record in shells:
            pair = _first(record, "pair", "label", "name")
            if pair is None and record.get("centre") and record.get("neighbour"):
                pair = f"{record['centre']}-{record['neighbour']}"
            row = [
                str(pair or DASH),
                format_fixed(
                    _first(record, "r_peak", "peak_position", "first_peak", "r_max"), 2
                ),
                format_fixed(_first(record, "g_peak", "peak_height", "g_max"), 2),
                format_fixed(_first(record, "r_min", "first_minimum", "minimum"), 2),
                format_fixed(
                    _first(record, "n_coord", "coordination_number", "coordination"), 2
                ),
            ]
            if has_reverse:
                row.append(format_fixed(record.get("n_coord_reverse"), 2))
            row.append(format_general(_first(record, "partner_density", "density"), 4))
            rows.append(row)
        self.table(rows, align=["left"] + ["right"] * (len(header) - 1))
        self.plain(
            "From the final intermolecular partial g(r) of each pair A-B: the "
            "coordination number is the number of B atoms around an A atom within the "
            "first minimum (reverse: A atoms around a B atom); the partner density is "
            "the number density of B.",
            "Small",
        )

    # ---- 5. discussion
    def discussion_section(self) -> None:
        self.heading(SECTION_HEADINGS[4])
        if _no_model(self.model_name):
            why = self.model_name.partition(":")[2].strip()
            opening = "<b>Not written by a language model" + (
                f" ({escape_text(why)})." if why else "."
            )
            opening += "</b> Every number in this report was computed by EPSR and the analysis code from the run folder."
        else:
            model = self.model_name or "model not named"
            opening = (
                f"<b>Written by the language model ({escape_text(model)}).</b> Every "
                "number in this report was computed by EPSR and the analysis code from "
                "the run folder; the text below interprets them and produced none of them."
            )
        self.markup(opening, "Note")
        if self.discussion.strip():
            self.story.extend(markdown_to_flowables(self.discussion, self.styles))
        else:
            self.plain("No discussion was written for this run.", "Small")

    # ---- 6. reproducibility
    def reproducibility_section(self) -> None:
        self.heading(SECTION_HEADINGS[5])
        folder = Path(self.result.folder).resolve()
        script_name = self.result.script_path.name
        command = f"python {script_name} --out {folder}"
        if self.quick:
            command += " --quick"
        self.plain(
            f"The script {script_name} is the complete input of this calculation: it "
            "embeds every EPSR input file (Appendix A), builds the box with the setup "
            "programs, runs the schedule of section 3 through the epsr wrapper and writes "
            "manifest.json and the trajectory. To repeat the calculation into the run "
            "folder (from the folder that holds the script):"
        )
        self.code(command, "CodeBlock")
        self.plain(
            "Requirements: Python 3.10 or newer with the epsr wrapper (it builds the Rust "
            "EPSR engine) and the epsrvis toolkit from this repository, e.g. 'uv sync' "
            "in mellor/. Add --quick for the short preview schedule and --no-record to "
            "skip the trajectory. Mellor itself is not needed."
        )
        self.subheading("Run folder")
        rows: list[list[Any]] = [["Path", "Content", "Present"]]
        for layout_path, content in FOLDER_LAYOUT:
            actual = script_name if layout_path == "run_<name>.py" else layout_path
            if (folder / actual.rstrip("/")).exists():
                present = "yes"
            elif actual == "report.pdf":
                present = "being written"
            else:
                present = "no"
            rows.append([actual, content, present])
        self.table(rows)

        self.subheading("Checksums (SHA-256)")
        self.plain(
            "Hashes of the inputs, so that a re-run can be checked against exactly these files.",
            "Small",
        )
        rows = [["File", "Role", "Bytes", "SHA-256"]]
        entries: list[tuple[str, Path]] = []
        for dataset in self.spec.datasets:
            found = self.locate_data_file(dataset)
            entries.append((f"data: {dataset.label}", found or Path(dataset.file)))
        if self.spec.pdb is not None:
            candidates = [folder / self.spec.pdb.file_name, Path(self.spec.pdb.file)]
            entries.append(
                (
                    "PDB file",
                    next((p for p in candidates if p.is_file()), candidates[-1]),
                )
            )
        entries.append(("specification", self.result.spec_path))
        entries.append(("script", self.result.script_path))
        for role, file_path in entries:
            if not file_path.is_file():
                rows.append([file_path.name, role, DASH, "file not found"])
                continue
            try:
                rows.append(
                    [
                        file_path.name,
                        role,
                        str(file_path.stat().st_size),
                        _sha256(file_path),
                    ]
                )
            except OSError as error:
                rows.append([file_path.name, role, DASH, f"not readable: {error}"])
        mono = self.styles["CellMono"]
        for row in rows[1:]:
            row[3] = Paragraph(escape_text(row[3]), mono)
        self.table(
            rows,
            widths=[0.24, 0.12, 0.08, 0.56],
            align=["left", "left", "right", "left"],
        )

    # ---- appendices
    def appendix_inputs(self) -> None:
        self.story.append(PageBreak())
        self.heading(SECTION_HEADINGS[6])
        if self.generated is not None:
            files = dict(self.generated.files())
            files.update(self.generated.answer_files())
            self.plain(
                "Regenerated from the specification (spec.json) by mellor.inputs.generate: "
                "the same texts are embedded in the script and written into the epsr folder. "
                "The .stdin files are the answers given to the interactive setup programs."
            )
        else:
            files = self.input_files_from_disk()
            self.plain(
                "The input files could not be regenerated from the specification "
                f"({self.generate_error}); the copies found in the epsr folder are shown instead.",
                "Note",
            )
        if not files:
            self.plain("No input files are available.", "Small")
        for name, text in files.items():
            self.subheading(name, "H3")
            self.code(text)

    def input_files_from_disk(self) -> dict[str, str]:
        folder = self.result.epsr_folder
        if not folder.is_dir():
            return {}
        files: dict[str, str] = {}
        patterns = [self.spec.inp_name, "*.NWTS.dat", "*.mol", "*.stdin"]
        for pattern in patterns:
            for path in sorted(folder.glob(pattern)):
                with contextlib.suppress(OSError, UnicodeDecodeError):
                    files[path.name] = path.read_text(encoding="utf-8")
        return files

    def appendix_script(self) -> None:
        self.story.append(PageBreak())
        self.heading(SECTION_HEADINGS[7])
        path = self.result.script_path
        text: str | None = None
        if path.is_file():
            try:
                text = path.read_text(encoding="utf-8")
                self.plain(
                    f"{path.name}, as found in the run folder ({len(text.splitlines())} lines)."
                )
            except (OSError, UnicodeDecodeError) as error:
                self.plain(f"{path.name} could not be read: {error}", "Note")
        if text is None and self.generated is not None:
            from .script import render_script

            text = render_script(self.spec, self.generated)
            self.plain(
                f"{path.name} was not found in the run folder; the script below was "
                "rendered again from the specification.",
                "Note",
            )
        if text is None:
            self.plain("The script is not available.", "Small")
            return
        self.code(text)


# --------------------------------------------------------------------------- public API
def build_pdf(
    result: RunResult,
    analysis: dict[str, Any],
    figures: dict[str, Path],
    discussion_markdown: str,
    out_path: str | Path | None = None,
    *,
    model_name: str = "",
    data_descriptions: list[dict[str, Any]] | None = None,
) -> Path:
    """Write the report (default ``<run>/report.pdf``) and return its path.

    ``analysis`` is the dict of ``mellor.analysis.analyse`` (any subset of its keys),
    ``figures`` maps figure names (``mellor.figures.FIGURE_CAPTIONS``) to PNG paths
    (missing files are skipped), ``discussion_markdown`` the model's narrative and
    ``data_descriptions`` the output of ``mellor.data.describe_all`` (read from the
    run folder's data files when not given). ``model_name`` names the language model
    that wrote the discussion; ``"none"`` or ``"none: <what the text is instead>"``
    states that no model wrote it (:data:`NO_MODEL_PREFIX`).
    """
    out_path = (
        Path(out_path) if out_path is not None else Path(result.folder) / "report.pdf"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    builder = _ReportBuilder(
        result, analysis, figures, discussion_markdown, model_name, data_descriptions
    )
    story = builder.build()
    document = SimpleDocTemplate(
        str(out_path),
        pagesize=PAGE_SIZE,
        leftMargin=MARGIN_SIDE,
        rightMargin=MARGIN_SIDE,
        topMargin=MARGIN_TOP,
        bottomMargin=MARGIN_BOTTOM,
        title=f"Mellor report: {builder.title}",
        author="Mellor",
        subject="EPSR refinement report",
        creator=f"Mellor {__version__} (reportlab)",
    )
    footer = _Footer(_shorten(builder.title), builder.fonts["sans"])
    document.build(story, onFirstPage=footer, onLaterPages=footer)
    return out_path


# --------------------------------------------------------------------------- markdown
_FENCE = re.compile(r"^\s*(```|~~~)\s*([\w+-]*)\s*$")
_HEADING = re.compile(r"^(#{1,6})\s+(.+?)\s*#*\s*$")
_RULE = re.compile(r"^\s*([-*_])(\s*\1){2,}\s*$")
_BULLET = re.compile(r"^(\s*)[-*+]\s+(.*)$")
_NUMBERED = re.compile(r"^(\s*)(\d{1,3})[.)]\s+(.*)$")
_QUOTE = re.compile(r"^\s*>\s?(.*)$")
_TABLE_ROW = re.compile(r"^\s*\|.*\|\s*$")
_TABLE_SEPARATOR = re.compile(r"^\s*\|?\s*:?-{1,}:?\s*(\|\s*:?-{1,}:?\s*)*\|?\s*$")
_CODE_SPAN = re.compile(r"(`[^`\n]+`)")
_BOLD_ITALIC = re.compile(r"\*\*\*(?=\S)(.+?)(?<=\S)\*\*\*")
_BOLD = re.compile(r"\*\*(?=\S)(.+?)(?<=\S)\*\*")
_BOLD_UNDERSCORE = re.compile(r"(?<!\w)__(?=\S)(.+?)(?<=\S)__(?!\w)")
_ITALIC = re.compile(r"(?<![\w*])\*(?=[^\s*])(.+?)(?<=[^\s*])\*(?![\w*])")
_ITALIC_UNDERSCORE = re.compile(r"(?<!\w)_(?=[^\s_])(.+?)(?<=[^\s_])_(?!\w)")
_LINK = re.compile(r"\[([^\]]+)\]\((https?://[^\s)]+)\)")
#: Marks a hard line break inside a paragraph until the markup is assembled.
HARD_BREAK = "\x00"


def _inline(text: str, mono: str) -> str:
    """Escaped paragraph markup for one line of Markdown (code spans, bold, italic, links)."""
    pieces = []
    for part in _CODE_SPAN.split(text):
        if len(part) >= 2 and part.startswith("`") and part.endswith("`"):
            pieces.append(f'<font face="{mono}">{escape_text(part[1:-1])}</font>')
            continue
        chunk = escape_text(part)
        chunk = _BOLD_ITALIC.sub(r"<b><i>\1</i></b>", chunk)
        chunk = _BOLD.sub(r"<b>\1</b>", chunk)
        chunk = _BOLD_UNDERSCORE.sub(r"<b>\1</b>", chunk)
        chunk = _ITALIC.sub(r"<i>\1</i>", chunk)
        chunk = _ITALIC_UNDERSCORE.sub(r"<i>\1</i>", chunk)
        chunk = _LINK.sub(r'<a href="\2" color="#2a5d9f"><u>\1</u></a>', chunk)
        pieces.append(chunk)
    return "".join(pieces)


def _split_table_row(line: str) -> list[str]:
    inner = line.strip()
    if inner.startswith("|"):
        inner = inner[1:]
    if inner.endswith("|") and not inner.endswith("\\|"):
        inner = inner[:-1]
    cells = inner.replace("\\|", "\x00").split("|")
    return [cell.replace("\x00", "|").strip() for cell in cells]


def _table_alignments(separator: str) -> list[str]:
    alignments = []
    for cell in _split_table_row(separator):
        if cell.endswith(":") and not cell.startswith(":"):
            alignments.append("right")
        else:
            alignments.append("left")
    return alignments


def _bullet_style(base: ParagraphStyle, level: int) -> ParagraphStyle:
    if level <= 0:
        return base
    return ParagraphStyle(
        f"{base.name}-level{level}",
        parent=base,
        leftIndent=base.leftIndent + 12 * level,
        bulletIndent=base.bulletIndent + 12 * level,
    )


def markdown_to_flowables(text: str, styles: Any = None) -> list[Any]:
    """A small Markdown subset (headings, paragraphs, bullet lists, **bold**, *italic*,
    `code`, fenced code blocks, simple pipe tables) as reportlab flowables.

    Also numbered lists, ``>`` quotes, horizontal rules, links and hard line breaks
    (two trailing spaces). Anything else is rendered as a plain paragraph. ``styles``
    is the sheet of :func:`make_styles` (default), a dict of ParagraphStyles or the
    reportlab sample stylesheet.
    """
    if styles is None:
        styles = make_styles()
    mono = register_fonts()["mono"]
    body = _style(styles, "Body")
    bullet_base = _style(styles, "Bullet")
    lines = (text or "").replace("\r\n", "\n").replace("\r", "\n").split("\n")
    flowables: list[Any] = []
    paragraph: list[str] = []

    def flush() -> None:
        if not paragraph:
            return
        # a line ending in two spaces is a hard break; the marker survives escaping
        joined = "".join(
            line.strip() + (HARD_BREAK if line.endswith("  ") else " ")
            for line in paragraph
        ).strip(" " + HARD_BREAK)
        markup = _inline(joined, mono).replace(HARD_BREAK, "<br/>")
        flowables.append(_paragraph(markup, body, joined.replace(HARD_BREAK, "\n")))
        paragraph.clear()

    index = 0
    count = len(lines)
    while index < count:
        line = lines[index]
        stripped = line.strip()
        if not stripped:
            flush()
            index += 1
            continue
        fence = _FENCE.match(line)
        if fence:
            flush()
            marker = fence.group(1)
            block: list[str] = []
            index += 1
            while index < count and lines[index].strip() != marker:
                block.append(lines[index])
                index += 1
            index += 1  # the closing fence
            flowables.append(
                _preformatted("\n".join(block), _style(styles, "CodeBlock"))
            )
            continue
        heading = _HEADING.match(stripped)
        if heading:
            flush()
            level = len(heading.group(1))
            style_name = "MdH1" if level == 1 else "MdH2" if level == 2 else "MdH3"
            title = heading.group(2)
            flowables.append(
                _paragraph(_inline(title, mono), _style(styles, style_name), title)
            )
            index += 1
            continue
        if _RULE.match(stripped):
            flush()
            flowables.append(
                HRFlowable(
                    width="100%", thickness=0.5, color=GRID, spaceBefore=4, spaceAfter=8
                )
            )
            index += 1
            continue
        if _TABLE_ROW.match(line):
            flush()
            block = []
            while index < count and _TABLE_ROW.match(lines[index]):
                block.append(lines[index])
                index += 1
            flowables.append(_markdown_table(block, styles, mono))
            continue
        bullet = _BULLET.match(line)
        numbered = _NUMBERED.match(line)
        if bullet or numbered:
            flush()
            if bullet:
                indent, item, marker = bullet.group(1), bullet.group(2), "•"
            else:
                assert numbered is not None
                indent, item = numbered.group(1), numbered.group(3)
                marker = f"{numbered.group(2)}."
            level = min(len(indent.expandtabs(4)) // 2, 3)
            index += 1
            while (
                index < count
                and lines[index].strip()
                and lines[index][0] in " \t"
                and not (_BULLET.match(lines[index]) or _NUMBERED.match(lines[index]))
                and not _FENCE.match(lines[index])
            ):
                item += " " + lines[index].strip()
                index += 1
            style = _bullet_style(bullet_base, level)
            flowables.append(
                _paragraph(_inline(item, mono), style, item, bulletText=marker)
            )
            continue
        quote = _QUOTE.match(line)
        if quote:
            flush()
            quoted = []
            while index < count and _QUOTE.match(lines[index]):
                match = _QUOTE.match(lines[index])
                assert match is not None
                quoted.append(match.group(1).strip())
                index += 1
            joined = " ".join(part for part in quoted if part)
            flowables.append(
                _paragraph(_inline(joined, mono), _style(styles, "Quote"), joined)
            )
            continue
        paragraph.append(line)
        index += 1
    flush()
    return flowables


def _markdown_table(block: list[str], styles: Any, mono: str) -> Flowable:
    """A pipe table (header row when the second line is a ``|---|`` separator)."""
    header = len(block) >= 2 and bool(_TABLE_SEPARATOR.match(block[1]))
    alignments = _table_alignments(block[1]) if header else []
    source_rows = [block[0]] + block[2:] if header else list(block)
    rows: list[list[Any]] = []
    for row_index, line in enumerate(source_rows):
        cells = []
        for column, cell in enumerate(_split_table_row(line)):
            if header and row_index == 0:
                style_name = "CellHead"
            elif column < len(alignments) and alignments[column] == "right":
                style_name = "CellRight"
            else:
                style_name = "Cell"
            cells.append(
                _paragraph(_inline(cell, mono), _style(styles, style_name), cell)
            )
        rows.append(cells)
    if not rows:
        return Spacer(1, 0)
    return _table(rows, styles, header=header, align=alignments)
