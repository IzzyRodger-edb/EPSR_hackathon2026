"""Command line: ``uv run epsr-vis <command> ...`` (see README.md)."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _masses(traj):
    from .cases import get_case
    case = get_case(traj.case, traj.atom_types)
    try:
        return {t: case.mass(t) for t in traj.atom_types}
    except KeyError:
        return None


def cmd_report(a):
    import matplotlib.pyplot as plt
    from . import plots
    from .io import Trajectory

    traj = Trajectory.load(a.trajectory)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    masses = _masses(traj)
    figs = {
        "01_convergence": lambda: plots.convergence(traj, masses),
        "02_fq_fit": lambda: plots.fq_fit(traj),
        "03_fq_evolution": lambda: plots.fq_evolution(traj),
        "04_residual_map": lambda: plots.residual_map(traj),
        "05_partials": lambda: plots.partials(traj),
        "06_potentials": lambda: plots.potentials(traj),
        "07_first_shell": lambda: plots.shells(traj),
        "08_local_structure": lambda: plots.structure_from_coordinates(traj),
        "09_displacement": lambda: plots.displacement(traj),
    }
    for name, make in figs.items():
        fig = make()
        if fig is None:
            continue
        fig.savefig(out / f"{name}.png")
        plt.close(fig)
        print(f"wrote {out / name}.png")
    rows = plots.summary_rows(traj, masses)
    with open(out / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(out / "summary.md", "w") as fh:
        keys = list(rows[0])
        fh.write("| quantity | " + " | ".join(f"iter {r['iteration']}" for r in rows) + " |\n")
        fh.write("|---|" + "---|" * len(rows) + "\n")
        for k in keys[1:]:
            vals = [r[k] for r in rows]
            fmt = (lambda v: f"{v:.3e}") if k.startswith(("rfactor", "chi2")) else (lambda v: f"{v:.4g}")
            fh.write(f"| {k} | " + " | ".join(fmt(v) for v in vals) + " |\n")
    print(f"wrote {out / 'summary.csv'}, {out / 'summary.md'}")
    if a.gif:
        plots.fq_animation(traj, out / "fq_fit.gif", every=a.every)
        print(f"wrote {out / 'fq_fit.gif'}")


def _viewer(a, off_screen):
    from .io import Trajectory
    from .viewer import ConfigurationViewer

    return ConfigurationViewer(Trajectory.load(a.trajectory), gr=a.gr, slab=a.slab, bonds=not a.no_bonds,
                               color_by=a.color_by, off_screen=off_screen, atom_scale=a.atom_scale)


def cmd_view(a):
    v = _viewer(a, off_screen=False)
    v.set_step(a.step - 1 if a.step > 0 else a.step)
    v.show()


def cmd_snapshot(a):
    v = _viewer(a, off_screen=True)
    v.screenshot(a.out, step=a.step - 1 if a.step > 0 else a.step)
    print(f"wrote {a.out}")


def cmd_movie(a):
    v = _viewer(a, off_screen=True)
    v.movie(a.out, every=a.every, fps=a.fps)
    print(f"wrote {a.out}")


def cmd_final(a):
    import matplotlib.pyplot as plt
    from . import plots

    fig = plots.final_outputs(a.folder)
    fig.savefig(a.out)
    plt.close(fig)
    print(f"wrote {a.out}")


def main(argv=None):
    p = argparse.ArgumentParser(prog="epsr-vis", description="Analysis and visualisation of EPSR refinements")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("report", help="all diagnostic figures + summary table for a recorded trajectory")
    r.add_argument("trajectory", help="trajectory folder written by the recorder")
    r.add_argument("--out", required=True, help="output folder for the figures")
    r.add_argument("--gif", action="store_true", help="also write an animated F(Q)/g(r) GIF")
    r.add_argument("--every", type=int, default=2, help="GIF: one frame every N iterations")
    r.set_defaults(func=cmd_report)

    def viewer_args(s):
        s.add_argument("trajectory")
        s.add_argument("--step", type=int, default=-1, help="EPSR iteration to show (1-based; -1 = last)")
        s.add_argument("--gr", choices=["inter", "total"], default="inter",
                       help="intermolecular g(r) only (as .g01), or including intramolecular")
        s.add_argument("--slab", type=float, default=None, help="only atoms with |z| < SLAB/2 (Å)")
        s.add_argument("--no-bonds", action="store_true")
        s.add_argument("--color-by", choices=["type", "coordination"], default="type")
        s.add_argument("--atom-scale", type=float, default=0.6, help="sphere size (0.6 ball-and-stick, 1 large)")

    v = sub.add_parser("view", help="interactive 3D configuration + partial g(r), with an iteration slider")
    viewer_args(v)
    v.set_defaults(func=cmd_view)

    s = sub.add_parser("snapshot", help="off-screen PNG of the 3D view + g(r) at one iteration")
    viewer_args(s)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_snapshot)

    m = sub.add_parser("movie", help="off-screen GIF/MP4 of the configuration + g(r) through the refinement")
    viewer_args(m)
    m.add_argument("--out", required=True)
    m.add_argument("--every", type=int, default=5)
    m.add_argument("--fps", type=int, default=10)
    m.set_defaults(func=cmd_movie)

    f = sub.add_parser("final", help="overview figure from the final EPSR output files of a run folder")
    f.add_argument("folder")
    f.add_argument("--out", required=True)
    f.set_defaults(func=cmd_final)

    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
