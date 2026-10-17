"""The app's own seaborn figures, drawn from Mellor's tidy CSV tables.

Mellor's figures are matplotlib PNGs from the `epsrvis` toolkit. These are the
app's: integration.md section 5.1 says not to parse EPSR files but to read the
long-format tables in `<run>/tables/`, which is what happens here.

Every figure is written twice -- PNG for the gallery and PDF for LaTeX -- because
the figures are meant to end up in a paper, and raster axes in a vector document
are the first thing a referee notices.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # no display on a server; must precede pyplot

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

LOG = logging.getLogger(__name__)

#: Where the app writes its own figures inside a run folder, kept apart from
#: `figures/` so that regenerating Mellor's report never clobbers them.
WEB_FIGURES_DIR = "webfigures"

_THEME = {
    "context": "paper",
    "style": "ticks",
    "palette": "colorblind",
    "font_scale": 1.05,
    "rc": {
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "axes.grid": True,
        "grid.alpha": 0.25,
        "legend.frameon": False,
    },
}


def _read(run_dir: Path, name: str) -> pd.DataFrame | None:
    path = run_dir / "tables" / f"{name}.csv"
    if not path.is_file():
        return None
    try:
        frame = pd.read_csv(path)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        LOG.warning("could not read %s: %s", path, exc)
        return None
    return frame if not frame.empty else None


def _save(fig: plt.Figure, out_dir: Path, name: str) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, str] = {}
    for suffix in ("png", "pdf"):
        target = out_dir / f"{name}.{suffix}"
        fig.savefig(target)
        written[suffix] = target.name
    plt.close(fig)
    return written


def _fit_overlay(fits: pd.DataFrame) -> plt.Figure:
    long = fits.melt(
        id_vars=["dataset", "q"], value_vars=["data", "simulated"],
        var_name="curve", value_name="fq",
    )
    grid = sns.relplot(
        data=long, x="q", y="fq", hue="curve", col="dataset", kind="line",
        col_wrap=min(3, long["dataset"].nunique()), height=3.1, aspect=1.25,
        facet_kws={"sharey": False}, linewidth=1.2,
    )
    grid.set_axis_labels(r"$Q\ (\mathrm{\AA}^{-1})$", r"$F(Q)$")
    grid.set_titles("{col_name}")
    grid.legend.set_title("")
    grid.figure.suptitle("Measured and simulated structure factor", y=1.02)
    return grid.figure


def _residuals(fits: pd.DataFrame) -> plt.Figure:
    grid = sns.relplot(
        data=fits, x="q", y="residual", col="dataset", kind="line",
        col_wrap=min(3, fits["dataset"].nunique()), height=3.1, aspect=1.25,
        facet_kws={"sharey": False}, linewidth=1.0, color="#b4432f",
    )
    for ax in grid.axes.flat:
        ax.axhline(0.0, color="0.4", linewidth=0.8, zorder=0)
    grid.set_axis_labels(r"$Q\ (\mathrm{\AA}^{-1})$", r"data $-$ simulation")
    grid.set_titles("{col_name}")
    grid.figure.suptitle("Residual structure factor", y=1.02)
    return grid.figure


def _convergence(history: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    hue = "run" if "run" in history.columns and history["run"].nunique() > 1 else None
    sns.lineplot(data=history, x="step", y="r_factor", hue=hue, ax=ax, linewidth=1.4)
    ax.set(xlabel="EPSR iteration", ylabel="R-factor", yscale="log")
    ax.set_title("Fit quality against iteration")
    if hue:
        ax.legend(title="run", ncols=2)
    return fig


def _partials(partials: pd.DataFrame) -> plt.Figure:
    column = "g_total" if "g_total" in partials.columns else "g_inter"
    near = partials[partials["r"] <= 8.0]
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    sns.lineplot(data=near, x="r", y=column, hue="pair", ax=ax, linewidth=1.3)
    ax.axhline(1.0, color="0.6", linewidth=0.8, linestyle=":", zorder=0)
    ax.set(xlabel=r"$r\ (\mathrm{\AA})$", ylabel=r"$g(r)$")
    ax.set_title("Partial pair distribution functions")
    ax.legend(title="pair", ncols=2)
    return fig


def _potentials(potentials: pd.DataFrame) -> plt.Figure:
    long = potentials.melt(
        id_vars=["pair", "r"], value_vars=[
            c for c in ("empirical", "reference") if c in potentials.columns
        ], var_name="potential", value_name="energy",
    )
    grid = sns.relplot(
        data=long[long["r"] <= 8.0], x="r", y="energy", hue="potential", col="pair",
        kind="line", col_wrap=min(3, long["pair"].nunique()), height=2.9, aspect=1.2,
        facet_kws={"sharey": False}, linewidth=1.2,
    )
    for ax in grid.axes.flat:
        ax.axhline(0.0, color="0.5", linewidth=0.8, zorder=0)
    grid.set_axis_labels(r"$r\ (\mathrm{\AA})$", "energy (kJ/mol)")
    grid.set_titles("{col_name}")
    grid.legend.set_title("")
    grid.figure.suptitle("Empirical and reference potential", y=1.02)
    return grid.figure


def make_seaborn_figures(run_dir: Path) -> dict[str, dict[str, str]]:
    """Draw every figure the run's tables support.

    Returns `{figure name: {"png": filename, "pdf": filename}}`. A table that is
    missing (for instance after a `--no-record` run, which has no per-iteration
    history) simply yields fewer figures; it is never an error.
    """
    out_dir = run_dir / WEB_FIGURES_DIR
    fits, history = _read(run_dir, "fits"), _read(run_dir, "history")
    partials, potentials = _read(run_dir, "partials"), _read(run_dir, "potentials")

    jobs: list[tuple[str, object, pd.DataFrame | None]] = [
        ("sb_fit_overlay", _fit_overlay, fits),
        ("sb_residuals", _residuals, fits),
        ("sb_convergence", _convergence, history),
        ("sb_partials", _partials, partials),
        ("sb_potentials", _potentials, potentials),
    ]

    written: dict[str, dict[str, str]] = {}
    with sns.axes_style(_THEME["style"]):
        sns.set_theme(**_THEME)
        for name, draw, frame in jobs:
            if frame is None:
                continue
            try:
                written[name] = _save(draw(frame), out_dir, name)  # type: ignore[operator]
            except Exception as exc:  # one bad table must not lose the others
                LOG.warning("figure %s failed: %s", name, exc)
                plt.close("all")
    return written


def existing_seaborn_figures(run_dir: Path) -> dict[str, dict[str, str]]:
    """What is already on disk, for serving a finished run without redrawing."""
    out_dir = run_dir / WEB_FIGURES_DIR
    found: dict[str, dict[str, str]] = {}
    if not out_dir.is_dir():
        return found
    for png in sorted(out_dir.glob("sb_*.png")):
        entry = {"png": png.name}
        if (out_dir / f"{png.stem}.pdf").is_file():
            entry["pdf"] = f"{png.stem}.pdf"
        found[png.stem] = entry
    return found
