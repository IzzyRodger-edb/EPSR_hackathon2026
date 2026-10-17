#!/usr/bin/env python3
"""Plot the partial g(r) of a validation case as computed by the three codes.

usage: plot_gr_agreement.py CASE FORTRAN_DIR RUST_DIR WRAPPER_DIR OUTPUT.png

CASE is water, silica_bulk or cu-btc; it selects the run name and the atom pairs
shown. The three folders hold the results of the same calculation (steps.txt)
run by the original Fortran EPSR, by the Rust CLI (epsr-rs) and by the Python
wrapper. The intermolecular partial g(r) is read from the <run>.g01 file of each
folder, the file the engine writes and the wrapper's Simulation.gr_inter returns.

For every pair the top panel overlays the three curves; the bottom panel shows the
differences of the Rust and Python results from the Fortran, drawn with the same
encodings. The largest absolute difference of each pair over the plotted range is
printed and written on the figure. The number of EPSR iterations in the subtitle
is read from the Fortran folder's <run>.erg file (one line per iteration).

Run with the wrapper's environment, from the python/ folder:

    uv run python ../validation/plot_gr_agreement.py cu-btc work/fortran work/rust work/wrapper cu-btc.png
"""

import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np

import epsr

matplotlib.use("Agg")  # write the figure to a file, no display needed
import matplotlib.pyplot as plt  # noqa: E402  (after the backend is chosen)


@dataclass(frozen=True)
class CaseLayout:
    """What to show for one validation case."""

    run_name: str  # prefix of the output files
    pairs: tuple[str, ...]  # column labels of the .g01 file
    r_max: float  # Å, the range shown
    title: str


CASES = {
    "water": CaseLayout(
        run_name="waterbox.EPSR",
        pairs=("OW-OW", "OW-HW", "HW-HW"),
        r_max=10.0,
        title="Water: 1000 H$_2$O fitted to SANDALS H$_2$O / D$_2$O / HDO data",
    ),
    "silica_bulk": CaseLayout(
        run_name="silica.EPSR",
        pairs=("Si-O", "O-O", "Si-Si"),
        r_max=10.0,
        title="Amorphous silica: 1500 Si + 3000 O fitted to NIMROD data",
    ),
    "cu-btc": CaseLayout(
        run_name="CuBTC.EPSR",
        pairs=("CU-O", "CU-CU", "O-O", "C3-H"),
        r_max=12.0,
        title="Cu-BTC (HKUST-1): one unit cell fitted to NIMROD data, with Bragg peaks",
    ),
}


@dataclass(frozen=True)
class CodeStyle:
    """How one code's curves are drawn: colour, legend text and matplotlib line options."""

    colour: str
    legend_label: str
    line_options: dict


SURFACE = "#fcfcfb"
TEXT_PRIMARY, TEXT_SECONDARY, GRIDLINE = "#0b0b0b", "#52514e", "#e1e0d9"
MARKER_EVERY = (
    12  # one marker per this many r points, so the lines beneath stay visible
)

# The three curves normally coincide, so each code gets its own encoding (a solid
# line, a dashed line, markers on a subset of points) in both panels. The colours
# are slots 1-3 of the data-visualisation reference palette, validated
# colour-blind safe for three series on a light surface.
FORTRAN, RUST, WRAPPER = "Fortran EPSR", "Rust CLI (epsr-rs)", "Python wrapper"
CODES = (FORTRAN, RUST, WRAPPER)  # the order of the folder arguments
STYLES = {
    FORTRAN: CodeStyle(
        colour="#2a78d6",
        legend_label="Fortran EPSR (solid line)",
        line_options={"linestyle": "-", "linewidth": 2.4},
    ),
    RUST: CodeStyle(
        colour="#eb6834",
        legend_label="Rust CLI epsr-rs (dashed line)",
        line_options={"linestyle": "--", "linewidth": 1.6},
    ),
    WRAPPER: CodeStyle(
        colour="#1baf7a",
        legend_label=f"Python wrapper (every {MARKER_EVERY}th point marked)",
        line_options={
            "linestyle": "none",
            "marker": "o",
            "markersize": 4.5,
            "markevery": MARKER_EVERY,
            "markeredgecolor": SURFACE,  # a surface ring keeps markers legible on the lines
            "markeredgewidth": 0.8,
        },
    ),
}


def read_gr(folder: Path, run_name: str) -> epsr.io.Table:
    return epsr.read_table(folder / f"{run_name}.g01")


def count_iterations(folder: Path, run_name: str) -> int:
    """EPSR iterations run in a folder: one line per iteration in the .erg file."""
    return len(epsr.read_erg(folder / f"{run_name}.erg"))


def style_axes(axes: plt.Axes) -> None:
    axes.set_facecolor(SURFACE)
    axes.grid(True, color=GRIDLINE, linewidth=1, linestyle="-")
    axes.set_axisbelow(True)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color(GRIDLINE)
    axes.tick_params(colors=TEXT_SECONDARY, labelsize=9)


def plot_pair(
    gr_axes: plt.Axes,
    difference_axes: plt.Axes,
    pair: str,
    tables: dict[str, epsr.io.Table],
    r_max: float,
) -> float:
    """Draw one pair's g(r) from the three codes and their differences; return max |Δ|."""
    fortran_table = tables[FORTRAN]
    shown = fortran_table.x <= r_max
    r = fortran_table.x[shown]
    fortran_gr = fortran_table[pair][shown]

    largest_difference = 0.0
    for code in CODES:
        style = STYLES[code]
        gr = tables[code][pair][shown]
        gr_axes.plot(
            r, gr, color=style.colour, label=style.legend_label, **style.line_options
        )
        if code == FORTRAN:
            continue
        difference = gr - fortran_gr
        largest_difference = max(largest_difference, float(np.max(np.abs(difference))))
        difference_axes.plot(
            r,
            difference,
            color=style.colour,
            label=f"{style.legend_label} − Fortran",
            **style.line_options,
        )

    first_type, second_type = pair.split("-")
    gr_axes.set_title(
        f"{first_type}–{second_type}", fontsize=11, color=TEXT_PRIMARY, loc="left"
    )
    gr_axes.set_xlim(0, r_max)
    gr_axes.set_ylim(bottom=min(0.0, float(fortran_gr.min())))
    # A fine scale: identical results give a flat line even where 0.1 would show.
    difference_half_range = max(0.1, 2.0 * largest_difference)
    difference_axes.set_xlim(0, r_max)
    difference_axes.set_ylim(-difference_half_range, difference_half_range)
    difference_axes.axhline(
        0.0, color=GRIDLINE, linewidth=1, zorder=0
    )  # under the data
    difference_axes.text(
        0.02,
        0.80,
        f"max |Δ| = {largest_difference:.1e}" if largest_difference else "max |Δ| = 0",
        transform=difference_axes.transAxes,
        fontsize=9,
        color=TEXT_SECONDARY,
    )
    difference_axes.set_xlabel("r (Å)", color=TEXT_SECONDARY)
    return largest_difference


def subtitle(n_iterations: int, all_identical: bool) -> str:
    text = f"Partial g(r) after {n_iterations} EPSR iterations from the three codes. "
    if all_identical:
        return text + "The curves coincide, so each code is drawn differently."
    return text + "The curves are overlaid; the lower panels show the differences."


def main() -> int:
    if len(sys.argv) != 6:
        print(__doc__)
        return 2
    case_name, *folder_names, output_name = sys.argv[1:]
    if case_name not in CASES:
        print(f"unknown case {case_name!r}; choose from {', '.join(CASES)}")
        return 2
    layout = CASES[case_name]
    folders = dict(zip(CODES, (Path(name) for name in folder_names), strict=True))
    tables = {
        code: read_gr(folder, layout.run_name) for code, folder in folders.items()
    }
    n_iterations = count_iterations(folders[FORTRAN], layout.run_name)

    n_pairs = len(layout.pairs)
    figure, axes_grid = plt.subplots(
        2,
        n_pairs,
        figsize=(3.6 * n_pairs + 0.6, 5.4),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.12, "wspace": 0.28},
        sharex="col",
        squeeze=False,
    )
    figure.patch.set_facecolor(SURFACE)
    for axes in axes_grid.flat:
        style_axes(axes)

    largest_differences = {}
    for column, pair in enumerate(layout.pairs):
        gr_axes, difference_axes = axes_grid[0, column], axes_grid[1, column]
        largest_differences[pair] = plot_pair(
            gr_axes, difference_axes, pair, tables, layout.r_max
        )
        if column == 0:
            gr_axes.set_ylabel("intermolecular g(r)", color=TEXT_SECONDARY)
            difference_axes.set_ylabel("difference from Fortran", color=TEXT_SECONDARY)

    all_identical = all(
        difference == 0.0 for difference in largest_differences.values()
    )
    figure.suptitle(layout.title, x=0.01, ha="left", fontsize=12, color=TEXT_PRIMARY)
    figure.text(
        0.01,
        0.925,
        subtitle(n_iterations, all_identical),
        fontsize=9,
        color=TEXT_SECONDARY,
    )
    # The legend has its own row under the subtitle, so long labels never meet the title.
    handles, labels = axes_grid[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.005, 0.905),
        ncol=3,
        frameon=False,
        fontsize=9,
        labelcolor=TEXT_PRIMARY,
    )
    figure.subplots_adjust(left=0.07, right=0.99, top=0.80, bottom=0.11)
    figure.savefig(output_name, dpi=150, facecolor=SURFACE)

    for pair, difference in largest_differences.items():
        print(
            f"{case_name:12s} {pair:7s} max |g(r) difference| from Fortran "
            f"(r <= {layout.r_max:g} Å): {difference:.3e}"
        )
    print(f"wrote {output_name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
