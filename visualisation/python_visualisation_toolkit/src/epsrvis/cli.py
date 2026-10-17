"""Command line: ``epsrvis <command> ...`` (see README.md)."""

from __future__ import annotations

import argparse


def cmd_report(a):
    from .report import make_report
    from .trajectory import Trajectory

    make_report(Trajectory.load(a.trajectory), a.out, gif=a.gif, gif_every=a.every)


def _viewer(a, off_screen):
    from .trajectory import Trajectory
    from .viewer import ConfigurationViewer

    return ConfigurationViewer(
        Trajectory.load(a.trajectory),
        gr=a.gr,
        slab=a.slab,
        bonds=not a.no_bonds,
        color_by=a.color_by,
        off_screen=off_screen,
        atom_scale=a.atom_scale,
    )


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
    p = argparse.ArgumentParser(
        prog="epsrvis", description="Analysis and visualisation of EPSR refinements"
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser(
        "report",
        help="all diagnostic figures + summary table for a recorded trajectory",
    )
    r.add_argument("trajectory", help="trajectory folder written by the recorder")
    r.add_argument("--out", required=True, help="output folder for the figures")
    r.add_argument(
        "--gif", action="store_true", help="also write an animated F(Q)/g(r) GIF"
    )
    r.add_argument(
        "--every", type=int, default=2, help="GIF: one frame every N iterations"
    )
    r.set_defaults(func=cmd_report)

    def viewer_args(s):
        s.add_argument("trajectory")
        s.add_argument(
            "--step",
            type=int,
            default=-1,
            help="EPSR iteration to show (1-based; -1 = last)",
        )
        s.add_argument(
            "--gr",
            choices=["inter", "total"],
            default="inter",
            help="intermolecular g(r) only (as .g01), or including intramolecular",
        )
        s.add_argument(
            "--slab", type=float, default=None, help="only atoms with |z| < SLAB/2 (Å)"
        )
        s.add_argument("--no-bonds", action="store_true")
        s.add_argument("--color-by", choices=["type", "coordination"], default="type")
        s.add_argument(
            "--atom-scale",
            type=float,
            default=0.6,
            help="sphere size (0.6 ball-and-stick, 1 large)",
        )

    v = sub.add_parser(
        "view",
        help="interactive 3D configuration + partial g(r), with an iteration slider",
    )
    viewer_args(v)
    v.set_defaults(func=cmd_view)

    s = sub.add_parser(
        "snapshot", help="off-screen PNG of the 3D view + g(r) at one iteration"
    )
    viewer_args(s)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_snapshot)

    m = sub.add_parser(
        "movie",
        help="off-screen GIF/MP4 of the configuration + g(r) through the refinement",
    )
    viewer_args(m)
    m.add_argument("--out", required=True)
    m.add_argument("--every", type=int, default=5)
    m.add_argument("--fps", type=int, default=10)
    m.set_defaults(func=cmd_movie)

    f = sub.add_parser(
        "final", help="overview figure from the final EPSR output files of a run folder"
    )
    f.add_argument("folder")
    f.add_argument("--out", required=True)
    f.set_defaults(func=cmd_final)

    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
