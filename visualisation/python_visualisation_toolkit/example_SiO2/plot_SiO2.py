"""Example: analyse and plot the EPSR refinement of amorphous silica (SiO2).

Run from python_visualisation_toolkit/:

    uv run python example_SiO2/plot_SiO2.py                 # figures -> example_SiO2/figures/
    uv run python example_SiO2/plot_SiO2.py --interactive   # also open the 3D viewer
    uv run python example_SiO2/plot_SiO2.py --movies        # also write the animated GIFs (slower)

Data (in example_SiO2/data/):
    trajectory/   one snapshot per EPSR iteration (300 iterations, 4 EPSR runs)
    epsr_output/  the normal EPSR output files after the last iteration

See example_SiO2/README.md for what each output shows.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from epsrvis import Trajectory, make_report, metrics, plots, style

HERE = Path(__file__).resolve().parent
SI_O_CUTOFF = 2.0  # Å, first minimum of the Si-O g(r)


def print_overview(traj: Trajectory) -> None:
    print(f"System: {traj.case}, {traj.meta['natom']} atoms "
          f"({', '.join(f'{t}: {np.sum(traj.atom_labels == t)}' for t in traj.atom_types)}), "
          f"box {traj.box(-1):.3f} Å, {traj.history['rho'][-1]:.4f} atoms/Å³")
    print(f"Data set: {', '.join(traj.datasets)};  pairs: {', '.join(traj.pairs)}")
    print(f"{traj.nsteps} EPSR iterations in {len(traj.run_starts)} runs:")
    for k, run, ereq in traj.run_starts:
        print(f"  run {run}: from iteration {int(traj.steps[k])}, ereq = {ereq:g} kJ/mol")


def print_fit_table(traj: Trajectory) -> None:
    fit = metrics.fit_metrics(traj)
    h = traj.history
    print("\nEnd of each run:")
    print(f"{'iter':>5} {'E_inter':>9} {'P(kbar)':>8} {'EP amp':>7} {'ereq':>6} {'R':>10} {'R²':>7}")
    for s in plots.key_steps(traj):
        print(f"{int(h['step'][s]):5d} {h['energd'][s]:9.4g} {h['pressd'][s]:8.4g} {h['energabs'][s]:7.2f} "
              f"{h['ereq'][s]:6.1f} {fit['rfactor'][s, 0]:10.3e} {fit['r2'][s, 0]:7.4f}")


def convergence_check(traj: Trajectory, last: int = 50) -> None:
    """Is the refinement finished?  Trends over the last ``last`` iterations."""
    h = traj.history
    fit = metrics.fit_metrics(traj)
    x = traj.steps[-last:].astype(float)
    print(f"\nConvergence over the last {last} iterations (linear trend per 100 iterations):")
    for name, y in (("R-factor", fit["rfactor_tot"][-last:]), ("pressure (kbar)", h["pressd"][-last:]),
                    ("energy (kJ/mol)", h["energd"][-last:]), ("EP amplitude", h["energabs"][-last:])):
        slope = np.polyfit(x, y, 1)[0] * 100
        print(f"  {name:16s} {y[-1]:11.4g}   trend {slope:+.3g} per 100 iterations")
    if h["energabs"][-1] < 0.95 * h["ereq"][-1]:
        print(f"  -> not converged: the EP amplitude ({h['energabs'][-1]:.1f}) has not reached "
              f"ereq ({h['ereq'][-1]:.0f}); run more iterations.")


def network_figure(traj: Trajectory, every: int = 10):
    """Silica-specific analysis: how the tetrahedral network forms during the refinement.

    (a) fraction of Si with 3, 4, 5 O neighbours; (b) fraction of bridging O
    (2 Si) and non-bridging O (1 Si); (c) mean tetrahedral order q of the Si.
    """
    labels = traj.atom_labels
    steps = np.unique(np.append(np.arange(0, traj.nsteps, every), traj.nsteps - 1))
    si_n = {3: [], 4: [], 5: []}
    o_n = {1: [], 2: []}
    q_si = []
    for s in steps:
        pos, box = np.asarray(traj.positions[s], dtype=float), traj.box(s)
        n, frac = metrics.coordination_distribution(pos, box, labels, "Si", "O", SI_O_CUTOFF)
        for k in si_n:
            si_n[k].append(frac[k] if k < len(frac) else 0.0)
        n, frac = metrics.coordination_distribution(pos, box, labels, "O", "Si", SI_O_CUTOFF)
        for k in o_n:
            o_n[k].append(frac[k] if k < len(frac) else 0.0)
        q_si.append(metrics.tetrahedral_order(pos, box, labels, "Si", "O").mean())
    x = traj.steps[steps]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    for i, (k, v) in enumerate(si_n.items()):
        axes[0].plot(x, v, color=style.series(i), label=f"{k}-fold Si", marker="o", ms=3)
    axes[0].set_title("Si coordination (O within 2.0 Å)")
    axes[0].set_ylabel("fraction of Si")
    for i, (k, v) in enumerate(o_n.items()):
        axes[1].plot(x, v, color=style.series(i), label=("non-bridging O (1 Si)" if k == 1 else "bridging O (2 Si)"),
                     marker="o", ms=3)
    axes[1].set_title("Oxygen bonding")
    axes[1].set_ylabel("fraction of O")
    axes[2].plot(x, q_si, color=style.series(0), marker="o", ms=3)
    axes[2].set_title("Mean tetrahedral order q of Si (4 nearest O)")
    for a in axes:
        a.set_xlabel("EPSR iteration")
        for k, _, _ in traj.run_starts[1:]:
            a.axvline(traj.steps[k] - 0.5, color=style.MUTED, lw=0.8, ls=(0, (3, 3)))
    axes[0].legend(fontsize=7)
    axes[1].legend(fontsize=7)
    fig.suptitle("Silica: formation of the tetrahedral network", x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig, dict(steps=x, si4=np.array(si_n[4]), bridging=np.array(o_n[2]), q_si=np.array(q_si))


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
    convergence_check(traj)

    # 2. the standard set of figures + summary tables
    print()
    make_report(traj, a.out, gif=a.movies)

    # 3. silica-specific analysis, composed from the toolkit functions
    fig, net = network_figure(traj)
    fig.savefig(a.out / "10_network.png")
    plt.close(fig)
    print(f"wrote {a.out / '10_network.png'}  (final: {100 * net['si4'][-1]:.1f} % 4-fold Si, "
          f"{100 * net['bridging'][-1]:.1f} % bridging O, <q_Si> = {net['q_si'][-1]:.3f})")

    # 4. overview from the plain EPSR output files
    fig = plots.final_outputs(a.data / "epsr_output")
    fig.savefig(a.out / "11_final_epsr_output.png")
    plt.close(fig)
    print(f"wrote {a.out / '11_final_epsr_output.png'}")

    # 5. 3D configuration + partial g(r); Si coloured by number of O neighbours; 12 Å slab
    from epsrvis.viewer import ConfigurationViewer

    for it, color_by in ((1, "type"), (100, "coordination"), (300, "coordination"), (300, "type")):
        v = ConfigurationViewer(traj, slab=12.0, color_by=color_by, off_screen=True)
        name = f"config_iter{it:03d}{'_coordination' if color_by == 'coordination' else ''}.png"
        v.screenshot(a.out / name, step=it - 1)
        v.pl.close()
        print(f"wrote {a.out / name}")
    if a.movies:
        ConfigurationViewer(traj, slab=12.0, color_by="coordination", off_screen=True).movie(
            a.out / "config_movie.gif", every=4)
        print(f"wrote {a.out / 'config_movie.gif'}")
    if a.interactive:
        ConfigurationViewer(traj, slab=12.0, color_by="coordination").show()


if __name__ == "__main__":
    main()
