"""Write the full set of diagnostic figures and summary tables for a trajectory."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

from . import plots
from .systems import get_case
from .trajectory import Trajectory

#: Figure name -> (function, what it shows).  ``make_report`` writes ``<name>.png``.
FIGURES = {
    "01_convergence": (
        plots.convergence,
        "energies, pressure, R-factor, R², EP amplitude, density, MC acceptance",
    ),
    "02_fq_fit": (
        plots.fq_fit,
        "measured vs simulated F(Q) and their difference at the end of each run",
    ),
    "03_fq_evolution": (plots.fq_evolution, "simulated F(Q) through the refinement"),
    "04_residual_map": (plots.residual_map, "data − simulation over Q and iteration"),
    "05_partials": (
        plots.partials,
        "partial g(r) and S(Q) at the middle and end of each run",
    ),
    "06_potentials": (plots.potentials, "empirical potential, and its effect on g(r)"),
    "07_first_shell": (
        plots.shells,
        "first-shell peak position, height and coordination number vs iteration",
    ),
    "08_local_structure": (
        plots.structure_from_coordinates,
        "coordination and bond-angle distributions",
    ),
    "09_displacement": (
        plots.displacement,
        "RMS atom displacement since the first iteration",
    ),
}


def masses_for(traj: Trajectory) -> dict[str, float] | None:
    """Atomic masses of the atom types, from the system settings (None if unknown)."""
    case = get_case(traj.case, traj.atom_types)
    try:
        return {t: case.mass(t) for t in traj.atom_types}
    except KeyError:
        return None


def summary_markdown(rows: list[dict]) -> str:
    """The rows of :func:`epsrvis.plots.summary_rows` as a Markdown table (one column per row)."""
    keys = list(rows[0])
    out = [
        "| quantity | " + " | ".join(f"iter {r['iteration']}" for r in rows) + " |",
        "|---|" + "---|" * len(rows),
    ]
    for k in keys[1:]:
        fmt = (
            (lambda v: f"{v:.3e}")
            if k.startswith(("rfactor", "chi2"))
            else (lambda v: f"{v:.4g}")
        )
        out.append(f"| {k} | " + " | ".join(fmt(r[k]) for r in rows) + " |")
    return "\n".join(out) + "\n"


def make_report(
    traj: Trajectory,
    out_dir: str | Path,
    figures=None,
    gif: bool = False,
    gif_every: int = 3,
    dpi: int = 150,
    verbose: bool = True,
) -> dict[str, Path]:
    """Write the figures of :data:`FIGURES` (or the subset ``figures``), ``summary.csv`` and ``summary.md``.

    Returns {name: path} of everything written.  With ``gif=True`` also writes
    ``fq_fit.gif`` (F(Q) fit and partial g(r), animated over the iterations).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    masses = masses_for(traj)
    written: dict[str, Path] = {}
    for name in figures or FIGURES:
        func, _ = FIGURES[name]
        fig = func(traj, masses) if func is plots.convergence else func(traj)
        if fig is None:  # e.g. no structural analyses defined for this system
            continue
        path = out / f"{name}.png"
        fig.savefig(path, dpi=dpi)
        plt.close(fig)
        written[name] = path
        if verbose:
            print(f"wrote {path}")
    rows = plots.summary_rows(traj, masses)
    with open(out / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (out / "summary.md").write_text(summary_markdown(rows))
    written["summary.csv"] = out / "summary.csv"
    written["summary.md"] = out / "summary.md"
    if gif:
        path = out / "fq_fit.gif"
        plots.fq_animation(traj, path, every=gif_every)
        written["fq_fit.gif"] = path
        if verbose:
            print(f"wrote {path}")
    return written
