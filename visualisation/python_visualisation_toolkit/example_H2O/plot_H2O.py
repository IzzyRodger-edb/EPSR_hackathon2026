"""Example: analyse and plot the EPSR refinement of liquid water (H2O / D2O / HDO).

Run from python_visualisation_toolkit/:

    uv run python example_H2O/plot_H2O.py                 # figures -> example_H2O/figures/
    uv run python example_H2O/plot_H2O.py --interactive   # also open the 3D viewer
    uv run python example_H2O/plot_H2O.py --movies        # also write the animated GIFs (slower)

Data (in example_H2O/data/):
    trajectory/   one snapshot per EPSR iteration (300 iterations, 4 EPSR runs)
    epsr_output/  the normal EPSR output files after the last iteration

See example_H2O/README.md for what each output shows.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from epsrvis import Trajectory, make_report, metrics, plots, style

HERE = Path(__file__).resolve().parent


def print_overview(traj: Trajectory) -> None:
    print(f"System: {traj.case}, {traj.meta['natom']} atoms in {traj.meta['nmol']} molecules, "
          f"box {traj.box(-1):.3f} Å, {traj.history['rho'][-1]:.4f} atoms/Å³")
    print(f"Data sets: {', '.join(traj.datasets)};  pairs: {', '.join(traj.pairs)}")
    print(f"{traj.nsteps} EPSR iterations in {len(traj.run_starts)} runs:")
    for k, run, ereq in traj.run_starts:
        print(f"  run {run}: from iteration {int(traj.steps[k])}, ereq = {ereq:g} kJ/mol")


def print_fit_table(traj: Trajectory) -> None:
    """R-factor and R² per data set at the end of each run."""
    fit = metrics.fit_metrics(traj)
    h = traj.history
    print("\nEnd of each run:")
    print(f"{'iter':>5} {'E_inter':>9} {'P(kbar)':>8} {'EP amp':>7}  " +
          "  ".join(f"{d:>21}" for d in traj.datasets))
    print(f"{'':>5} {'kJ/mol':>9} {'':>8} {'kJ/mol':>7}  " + "  ".join(f"{'R':>10} {'R²':>10}" for _ in traj.datasets))
    for s in plots.key_steps(traj):
        cols = "  ".join(f"{fit['rfactor'][s, d]:10.3e} {fit['r2'][s, d]:10.4f}" for d in range(len(traj.datasets)))
        print(f"{int(h['step'][s]):5d} {h['energd'][s]:9.4g} {h['pressd'][s]:8.4g} {h['energabs'][s]:7.3f}  {cols}")


def hydrogen_bond_figure(traj: Trajectory, every: int = 10):
    """Water-specific analysis built from the toolkit's coordinate functions.

    Left: mean number of hydrogen bonds per molecule and mean tetrahedral order
    parameter q of the O atoms, against EPSR iteration.  Right: distribution of q
    at the end of the reference-potential stage and at the end of the refinement.
    """
    labels, mol = traj.atom_labels, traj.molecule_of_atom
    steps = np.arange(0, traj.nsteps, every)
    steps = np.unique(np.append(steps, traj.nsteps - 1))
    n_hb, q_mean = [], []
    for s in steps:
        pos, box = np.asarray(traj.positions[s], dtype=float), traj.box(s)
        # H-bond: O...H closer than 2.45 Å, H on another molecule (first minimum of the O-H g(r))
        i, j, _ = metrics.neighbour_pairs(pos, box, labels, "OW", "HW", 2.45)
        n_bonds = np.sum(mol[i] != mol[j])
        n_hb.append(2.0 * n_bonds / np.sum(labels == "OW"))  # each bond is shared by two molecules
        q_mean.append(metrics.tetrahedral_order(pos, box, labels, "OW", "OW").mean())

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    x = traj.steps[steps]
    axes[0].plot(x, n_hb, color=style.series(0), marker="o", ms=3)
    axes[0].set_title("Hydrogen bonds per molecule")
    axes[1].plot(x, q_mean, color=style.series(1), marker="o", ms=3)
    axes[1].set_title("Mean tetrahedral order q (O)")
    for a in axes[:2]:
        a.set_xlabel("EPSR iteration")
        for k, _, _ in traj.run_starts[1:]:
            a.axvline(traj.steps[k] - 0.5, color=style.MUTED, lw=0.8, ls=(0, (3, 3)))
    ref_end = max(k for k in range(traj.nsteps) if traj.history["ereq"][k] == 0)
    for c, s, lab in ((style.INK_2, ref_end, "reference potential only"), (style.series(0), traj.nsteps - 1, "with EP")):
        qv = metrics.tetrahedral_order(traj.positions[s], traj.box(s), labels, "OW", "OW")
        hist, edges = np.histogram(qv, bins=40, range=(-0.5, 1.0), density=True)
        axes[2].plot(0.5 * (edges[1:] + edges[:-1]), hist, color=c, label=f"{lab} (iter {int(traj.steps[s])})")
    axes[2].set_title("Distribution of q")
    axes[2].set_xlabel("q")
    axes[2].legend(loc="upper left", fontsize=7)
    fig.suptitle("Water: hydrogen bonding and tetrahedrality", x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig, dict(steps=x, n_hbonds=np.array(n_hb), q_mean=np.array(q_mean))


def isotope_contrast_figure(traj: Trajectory):
    """How the fit to each isotopic sample improves: R-factor normalised to its first value."""
    fit = metrics.fit_metrics(traj)
    fig, a = plt.subplots(figsize=(6.5, 3.6))
    for d, name in enumerate(traj.datasets):
        a.plot(traj.steps, fit["rfactor"][:, d] / fit["rfactor"][0, d], color=style.series(d), label=name)
    a.set_yscale("log")
    a.set_xlabel("EPSR iteration")
    a.set_ylabel("R / R(iteration 1)")
    a.set_title("Relative improvement of the fit per isotopic sample")
    a.legend(fontsize=8)
    fig.tight_layout()
    return fig


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", default=HERE / "data", type=Path)
    p.add_argument("--out", default=HERE / "figures", type=Path)
    p.add_argument("--interactive", action="store_true", help="open the pyvista viewer at the end")
    p.add_argument("--movies", action="store_true", help="write the animated GIFs (a few minutes)")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    # 1. load the per-iteration record (memory-mapped)
    traj = Trajectory.load(a.data / "trajectory")
    print_overview(traj)
    print_fit_table(traj)

    # 2. the standard set of figures + summary tables
    print()
    make_report(traj, a.out, gif=a.movies)

    # 3. water-specific analyses, composed from the toolkit functions
    fig, hb = hydrogen_bond_figure(traj)
    fig.savefig(a.out / "10_hydrogen_bonds.png")
    plt.close(fig)
    print(f"wrote {a.out / '10_hydrogen_bonds.png'}  "
          f"(final: {hb['n_hbonds'][-1]:.2f} H-bonds/molecule, <q> = {hb['q_mean'][-1]:.3f})")
    fig = isotope_contrast_figure(traj)
    fig.savefig(a.out / "11_isotope_fit.png")
    plt.close(fig)
    print(f"wrote {a.out / '11_isotope_fit.png'}")

    # 4. overview from the plain EPSR output files (works for any EPSR folder, no recording needed)
    fig = plots.final_outputs(a.data / "epsr_output")
    fig.savefig(a.out / "12_final_epsr_output.png")
    plt.close(fig)
    print(f"wrote {a.out / '12_final_epsr_output.png'}")

    # 5. 3D configuration next to the partial g(r) (pyvista, off-screen); slab of 10 Å through the box centre
    from epsrvis.viewer import ConfigurationViewer

    for it in (1, 100, 300):
        v = ConfigurationViewer(traj, slab=10.0, off_screen=True)
        v.screenshot(a.out / f"config_iter{it:03d}.png", step=it - 1)
        v.pl.close()
        print(f"wrote {a.out / f'config_iter{it:03d}.png'}")
    if a.movies:
        ConfigurationViewer(traj, slab=10.0, off_screen=True).movie(a.out / "config_movie.gif", every=4)
        print(f"wrote {a.out / 'config_movie.gif'}")
    if a.interactive:
        ConfigurationViewer(traj, slab=10.0).show()


if __name__ == "__main__":
    main()
