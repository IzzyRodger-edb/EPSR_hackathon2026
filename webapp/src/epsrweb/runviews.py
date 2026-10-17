"""Turning a run folder into what the results page needs.

Kept apart from the routes because it is the part most likely to need changing:
it reads files Mellor writes, and Mellor's `analysis.json` is a rich document
whose exact shape may grow. Everything here therefore tolerates missing keys and
falls back to showing the raw JSON rather than failing.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import events, mellor_adapter, plots, storage

#: Subfolders of a run that may be served to the owner. `trajectory/` is
#: deliberately absent: it is tens to hundreds of megabytes and is for Mellor's
#: own analysis, not for download (integration.md section 5).
SERVABLE_DIRS = ("figures", plots.WEB_FIGURES_DIR, "tables", "data", "uploads")
SERVABLE_FILES = (
    "report.pdf", "analysis.json", "manifest.json", "spec.json", "discussion.md",
    "run.log", "worker.log", "mellor_transcript.jsonl",
)


@dataclass
class Figure:
    key: str
    caption: str
    source: str          # "mellor" or "seaborn"
    png: str | None      # path relative to the run folder
    pdf: str | None

    @property
    def best(self) -> str:
        return self.pdf or self.png or ""


def list_figures(run_dir: Path) -> list[Figure]:
    """Mellor's figures first, then the app's seaborn ones."""
    captions = mellor_adapter.figure_captions()
    found: list[Figure] = []

    figures_dir = run_dir / "figures"
    if figures_dir.is_dir():
        for png in sorted(figures_dir.glob("*.png")):
            pdf = figures_dir / f"{png.stem}.pdf"
            found.append(Figure(
                key=png.stem,
                caption=captions.get(png.stem, png.stem.replace("_", " ")),
                source="mellor",
                png=f"figures/{png.name}",
                pdf=f"figures/{pdf.name}" if pdf.is_file() else None,
            ))

    for key, files in plots.existing_seaborn_figures(run_dir).items():
        found.append(Figure(
            key=key,
            caption=captions.get(key, key.replace("sb_", "").replace("_", " ")),
            source="seaborn",
            png=f"{plots.WEB_FIGURES_DIR}/{files['png']}" if files.get("png") else None,
            pdf=f"{plots.WEB_FIGURES_DIR}/{files['pdf']}" if files.get("pdf") else None,
        ))
    return found


def figure_by_key(run_dir: Path, key: str) -> Figure | None:
    for figure in list_figures(run_dir):
        if figure.key == key:
            return figure
    return None


def _as_rows(value: Any) -> list[dict[str, Any]]:
    """A list of flat dicts, whatever shape the key happened to have."""
    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]
    if isinstance(value, dict):
        rows = []
        for name, item in value.items():
            if isinstance(item, dict):
                rows.append({"name": name, **item})
            else:
                rows.append({"name": name, "value": item})
        return rows
    return []


def analysis_view(run_dir: Path) -> dict[str, Any]:
    """The pieces of `analysis.json` the results page shows.

    Mellor documents the contents as "runs, fit quality per data set,
    convergence verdict, first-shell table, flags". Anything else in the file is
    still available through the raw JSON download.
    """
    analysis = storage.read_json(run_dir / "analysis.json")
    if analysis is None:
        return {"present": False}
    fit = _as_rows(analysis.get("fit") or analysis.get("datasets") or analysis.get("fit_quality"))
    runs = _as_rows(analysis.get("runs"))
    shells = _as_rows(analysis.get("shells") or analysis.get("first_shell"))
    convergence = analysis.get("convergence")
    flags = analysis.get("flags") or []
    return {
        "present": True,
        "fit": fit,
        "fit_columns": _columns(fit),
        "runs": runs,
        "run_columns": _columns(runs),
        "shells": shells,
        "shell_columns": _columns(shells),
        "convergence": convergence if isinstance(convergence, dict) else None,
        "verdict": (convergence or {}).get("verdict") if isinstance(convergence, dict) else None,
        "flags": [str(flag) for flag in flags] if isinstance(flags, list) else [str(flags)],
        "demo": bool(analysis.get("demo")),
    }


def _columns(rows: list[dict[str, Any]]) -> list[str]:
    """Union of the keys, in first-seen order, so a table can be rendered."""
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    return columns


def downloads(run_dir: Path) -> list[dict[str, str]]:
    """Every file worth offering, with a label and a size."""
    offers: list[dict[str, str]] = []
    labels = {
        "report.pdf": "PDF report",
        "analysis.json": "analysis.json (every computed number)",
        "spec.json": "spec.json (the model system that was run)",
        "manifest.json": "manifest.json (per-run summaries)",
        "discussion.md": "discussion.md (the written discussion)",
        "mellor_transcript.jsonl": "mellor_transcript.jsonl (every model message and tool call)",
        "run.log": "run.log (console output of the generated script)",
        "worker.log": "worker.log (this app's own log for the job)",
    }
    for name, label in labels.items():
        path = run_dir / name
        if path.is_file():
            offers.append({"name": name, "label": label,
                           "size": storage.human_size(path.stat().st_size)})
    for script in sorted(run_dir.glob("run_*.py")):
        offers.append({"name": script.name,
                       "label": f"{script.name} (standalone script: reproduce this run)",
                       "size": storage.human_size(script.stat().st_size)})
    for table in sorted((run_dir / "tables").glob("*.csv")):
        offers.append({"name": f"tables/{table.name}",
                       "label": f"tables/{table.name} (long-format, for pandas)",
                       "size": storage.human_size(table.stat().st_size)})
    return offers


def log_lines(run_dir: Path, limit: int = 400) -> list[str]:
    """The event stream rendered as text, for the log pane on first load."""
    lines = [events.summarise(event) for event in events.read_events(run_dir)]
    lines = [line for line in lines if line.strip()]
    return lines[-limit:]


def run_summary(row: sqlite3.Row, run_dir: Path) -> dict[str, Any]:
    """One row for the dashboard."""
    return {
        "id": row["id"],
        "title": row["title"],
        "status": row["status"],
        "created_at": row["created_at"],
        "finished_at": row["finished_at"],
        "backend": row["backend"],
        "quick": bool(row["quick"]),
        "demo": bool(row["demo"]),
        "spec_source": row["spec_source"],
        "error": row["error"],
        "size": storage.human_size(storage.dir_size(run_dir)),
        "has_report": (run_dir / "report.pdf").is_file(),
    }
