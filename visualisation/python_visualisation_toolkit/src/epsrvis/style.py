"""Shared plotting style.

Series colours are the first slots of a colour-vision-deficiency-validated
categorical palette, assigned in fixed order (pair / data set index), never
cycled.  Measured data are always drawn in ink (black), simulations in colour.
"Progress through the refinement" uses a single-hue light->dark blue ramp.
"""

from __future__ import annotations

import matplotlib as mpl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"

# one-hue sequential ramp (light -> dark) for "iteration number"
_RAMP = ["#b7d3f6", "#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]
PROGRESS = LinearSegmentedColormap.from_list("epsr_progress", _RAMP)
# diverging pair with a neutral midpoint, for residual maps
DIVERGING = LinearSegmentedColormap.from_list(
    "epsr_diverging", ["#1c5cab", "#86b6ef", "#f2f1ee", "#f0a07c", "#b8461b"])


def series(i: int) -> str:
    return SERIES[i % len(SERIES)]


def progress_colors(n: int) -> list:
    return [PROGRESS(x) for x in np.linspace(0.0, 1.0, n)]


def apply() -> None:
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK_2,
        "axes.titlecolor": INK,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": 9,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": INK_2,
        "ytick.color": INK_2,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "lines.linewidth": 1.6,
        "font.size": 9,
        "figure.dpi": 110,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })
