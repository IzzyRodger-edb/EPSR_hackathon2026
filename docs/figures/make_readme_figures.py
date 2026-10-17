"""Make the plain-language figures in the top-level README.

Run from the repository root with any Python that has numpy and matplotlib, e.g.

    visualisation/python_visualisation_toolkit/.venv/bin/python docs/figures/make_readme_figures.py

Inputs (all in the repository):
  visualisation/data/water/trajectory/   every iteration of the 300-iteration water refinement
  experimental_data/water/reference_partials/   Soper (2000) water partials
  docs/figures/data/                      EPSR energy logs (.erg) of the Fortran/Rust runs, see README there
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
TRAJ = ROOT / "visualisation/data/water/trajectory"
DATA = OUT / "data"

# reference palette (dataviz skill), light mode
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8984"
GRID = "#e6e5e1"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.size": 11,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK2,
        "axes.titlecolor": INK,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "legend.frameon": False,
        "lines.linewidth": 2,
    }
)


def headline(fig, title: str, subtitle: str) -> None:
    fig.text(0.012, 0.975, title, fontsize=15, fontweight="bold", color=INK, va="top")
    fig.text(0.012, 0.905, subtitle, fontsize=11, color=INK2, va="top")


# ------------------------------------------------------------------ trajectory
meta = json.loads((TRAJ / "meta.json").read_text())
q = np.asarray(meta["q"])
r = np.asarray(meta["r"])
nd, nq, npair, nr = meta["ndata"], meta["nq"], meta["npairs"], meta["nr"]


def load(name: str, *shape: int) -> np.ndarray:
    return np.fromfile(TRAJ / f"{name}.f32", dtype="<f4").reshape(-1, *shape).astype(float)


fq_data, fq_sim = load("fq_data", nd, nq), load("fq_sim", nd, nq)
gr = load("gr_inter", npair, nr) + load("gr_intra", npair, nr)
nsteps = fq_sim.shape[0]


def r_factor(step: int, d: int) -> float:
    lo, hi = meta["q_range"][d]
    diff = fq_data[step, d, lo : hi + 1] - fq_sim[step, d, lo : hi + 1]
    return float(np.mean(diff**2))


rf = np.array([[r_factor(s, d) for d in range(nd)] for s in range(nsteps)])
SAMPLES = {  # dataset index -> plain name, colour
    0: ("heavy water (D₂O)", BLUE),
    1: ("50:50 mixture (HDO)", ORANGE),
    2: ("ordinary water (H₂O)", AQUA),
}


def sci(x: float) -> str:
    m, e = f"{x:.1e}".split("e")
    return f"{m} × 10$^{{{int(e)}}}$"


# ------------------------------------------- 1. the model learns to match the data
def fig_fit() -> None:
    d = 0  # D2O: the sample that starts furthest from the data
    stages = [
        (0, "Start", "molecules placed at random"),
        (99, "After 100 steps", "physics model only"),
        (nsteps - 1, "After 300 steps", "refined against the data"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), sharey=True)
    lo, hi = meta["q_range"][d]
    keep = slice(lo, int(np.searchsorted(q, 16.0)))
    for ax, (step, title, sub) in zip(axes, stages):
        ax.plot(q[keep], fq_data[step, d, keep], "o", ms=3.2, color=INK, alpha=0.55, label="measured")
        ax.plot(q[keep], fq_sim[step, d, keep], color=BLUE, label="simulated")
        ax.set_title(f"{title}: {sub}", fontsize=11.5)
        ax.text(0.97, 0.06, f"mismatch {sci(rf[step, d])}", transform=ax.transAxes,
                ha="right", color=INK, fontsize=11,
                bbox=dict(boxstyle="round,pad=0.3", fc=SURFACE, ec=GRID))
        ax.set_xlabel("Q (Å⁻¹)")
        ax.set_xlim(0, 16)
    axes[0].set_ylabel("scattering signal F(Q)")
    axes[0].legend(loc="upper right")
    headline(fig, "The simulation learns to reproduce a real experiment",
             "Heavy water measured on the SANDALS instrument (ISIS) against EPSR's simulated signal. "
             "Black dots are the experiment; the blue line is the model.")
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.14, top=0.76, wspace=0.08)
    fig.savefig(OUT / "fit_before_after.png", dpi=150)
    plt.close(fig)


# ------------------------------------------------ 2. mismatch over the refinement
def fig_mismatch() -> None:
    steps = np.arange(1, nsteps + 1)
    fig, ax = plt.subplots(figsize=(13, 4.1))
    ax.axvspan(1, 100, color=GRID, alpha=0.45, lw=0)
    for x, label in [(50, "physics model only"), (200, "data-driven correction switched on")]:
        ax.text(x, 0.075, label, ha="center", color=INK2, fontsize=10.5)
    for d, (name, colour) in SAMPLES.items():
        ax.plot(steps, rf[:, d], color=colour, label=name)
        ax.annotate(name, (steps[-1], rf[-1, d]), xytext=(8, (d - 1) * 14),
                    textcoords="offset points", va="center", color=INK, fontsize=10.5)
    ax.set_yscale("log")
    ax.set_ylim(1.2e-4, 1.1e-1)
    ax.set_xlim(1, nsteps)
    ax.set_xlabel("refinement step")
    ax.set_ylabel("mismatch with experiment\n(R-factor, log scale)")
    first, last = rf[0].mean(), rf[-1].mean()
    headline(fig, f"The mismatch with experiment falls about {first / last:.0f}-fold",
             "All three isotopic samples of water are fitted by the same model at the same time.")
    fig.subplots_adjust(left=0.08, right=0.86, bottom=0.15, top=0.78)
    fig.savefig(OUT / "mismatch_history.png", dpi=150)
    plt.close(fig)


# ------------------------------- 3. structure agrees with the published reference
def read_isis(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = []
    for line in path.read_text(encoding="latin-1").splitlines():
        try:
            vals = [float(x) for x in line.split()[:4]]
        except ValueError:
            continue
        if len(vals) == 4:
            rows.append(vals)
    a = np.array(rows)
    a = a[(a[:, 1] > 0) | (np.arange(len(a)) == 0)]  # drop the all-zero terminator row
    return a[:, 1], a[:, 2]


def fig_structure() -> None:
    ref = ROOT / "experimental_data/water/reference_partials"
    panels = [
        (0, "h2o298aoogor.txt", "oxygen ↔ oxygen", [(2.79, "neighbouring molecules\n2.8 Å apart")]),
        (1, "h2o298aohgor.txt", "oxygen ↔ hydrogen", [(1.80, "hydrogen bond\n1.8 Å")]),
        (2, "h2o298ahhgor.txt", "hydrogen ↔ hydrogen", [(1.55, "the two H atoms\nof one molecule")]),
    ]
    g = gr[-1]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
    for ax, (p, fname, title, notes) in zip(axes, panels):
        rr, gg = read_isis(ref / fname)
        ax.plot(rr, gg, color=INK, lw=2.4, label="published (Soper 2000)")
        ax.plot(r, g[p], color=ORANGE, lw=1.6, label="this work")
        ax.set_title(title)
        ax.set_xlim(0, 8)
        ax.set_ylim(0, 3.0)
        ax.set_xlabel("distance r (Å)")
        for x, txt in notes:
            ax.annotate(txt, (x, np.interp(x, r, g[p])),
                        xytext=(x + 1.3, 2.3), color=INK2, fontsize=10,
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    axes[0].set_ylabel("how likely a neighbour is\nat distance r  (g(r))")
    axes[1].text(1.12, 2.93, "O–H bond within a molecule (off scale)", color=INK2, fontsize=9, va="top")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper right", bbox_to_anchor=(0.99, 0.99), ncol=2)
    headline(fig, "The refined model recovers the known structure of water",
             "Pair distribution functions after 300 steps against published reference partials (Soper, Chem. "
             "Phys. 258, 121 (2000)).\nA peak marks a favoured distance between two atoms; 1 means no preference.")
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.14, top=0.72, wspace=0.18)
    fig.savefig(OUT / "structure_vs_published.png", dpi=150)
    plt.close(fig)


# ---------------------------------------- 4. Rust == Fortran, and why that is hard
def read_erg(path: Path) -> np.ndarray:
    """Column 1 of an EPSR .erg file: energy per molecule after each iteration (kJ/mol)."""
    return np.array([float(line.split()[0]) for line in path.read_text().splitlines() if line.strip()])


def fig_identical() -> None:
    fortran = read_erg(DATA / "water_fortran.erg")
    rust = read_erg(DATA / "water_rust.erg")
    c_f = read_erg(DATA / "chaos_fortran.erg")
    c_r = read_erg(DATA / "chaos_rust.erg")
    c_n = read_erg(DATA / "chaos_nudged.erg")
    assert np.array_equal(fortran, rust) and np.array_equal(c_f, c_r)

    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 4.4), gridspec_kw={"width_ratios": [1.15, 1]})
    steps = np.arange(1, len(fortran) + 1)
    show = steps >= 20  # the first iterations start from a random, overlapping box
    a.plot(steps[show], fortran[show], color=BLUE, lw=4, label="original Fortran")
    a.plot(steps[show], rust[show], color=ORANGE, lw=1.6, ls=(0, (4, 3)), label="our Rust translation")
    a.set_title("Same refinement, two programs")
    a.set_xlabel("refinement step")
    a.set_ylabel("energy per molecule (kJ/mol)")
    a.legend(loc="upper right")
    a.text(0.03, 0.07, f"all {len(fortran)} values identical\nto every printed digit",
           transform=a.transAxes, color=INK, fontsize=11,
           bbox=dict(boxstyle="round,pad=0.35", fc=SURFACE, ec=GRID))
    a.set_xlim(0, len(fortran) + 2)

    k = np.arange(1, len(c_f) + 1)
    b.plot(k, np.abs(c_n - c_f), color=INK2, marker="o", ms=4, lw=1.6,
           label="Rust, one molecule moved by 0.000001 Å")
    b.plot(k, np.abs(c_r - c_f), color=ORANGE, marker="o", ms=4, lw=2, label="Rust, unchanged input")
    b.set_title("Why that is a demanding test")
    b.set_xlabel("further refinement steps")
    b.set_ylabel("difference from Fortran (kJ/mol)")
    b.legend(loc="upper left")
    b.set_xlim(0, len(c_f) + 1)
    b.set_ylim(bottom=-0.02 * max(np.abs(c_n - c_f).max(), 1e-3))

    headline(fig, "The Rust translation reproduces the original program exactly",
             "Left: the 300-step water refinement run by both programs. Right: the simulation is chaotic, so moving one "
             "molecule by a millionth\nof an ångström changes the result after 13 steps. The translation never "
             "differs, because it performs exactly the same arithmetic.")
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.14, top=0.72, wspace=0.2)
    fig.savefig(OUT / "rust_matches_fortran.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    print("R-factor check, step 1:", rf[0])
    fig_fit()
    fig_mismatch()
    fig_structure()
    fig_identical()
    print("written to", OUT)
