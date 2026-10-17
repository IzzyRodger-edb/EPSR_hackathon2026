"""Water movie with the experimental partials fading out and back in -> H2O_config_movie_2.mp4.

Run from python_visualisation_toolkit/:

    uv run python presentation_assets/make_H2O_movie_2.py

Same content as H2O_config_movie.mp4, but the experimental g(r) fades out over the
first FADE_OUT_ITERS iterations, so the refinement is seen alone; while it is gone
the EPSR curves are drawn opaque.  The final distribution is held for HOLD seconds,
then the experiment fades back in (EPSR back to translucent) and is held again.
The R^2 of each EPSR partial against the experiment (over the plotted 0-10 A, on the
experimental r grid) is shown at the right of the g(r) panel for every frame.
"""

from __future__ import annotations

import argparse

import numpy as np

from epsrvis import Trajectory, style
from epsrvis.viewer import ConfigurationViewer

from make_presentation_assets import CALC_OPACITY, HERE, ROOT, SIZE, water_reference

FADE_OUT_ITERS = 100   # experiment gone by the end of the reference-potential stage
HOLD = 3.0             # s, final distribution without experiment
FADE_IN = 1.5          # s
HOLD_END = 3.0         # s, final distribution with experiment
R2_RMAX = 10.0         # Å, R^2 over the plotted range


def r_squared(traj, ref, step):
    """R^2 = 1 - SS_res/SS_tot of each EPSR partial (inter + intra) against the reference."""
    g = traj.gr_total(step)
    out = {}
    for ip, pair in enumerate(traj.pairs):
        rr, gg = (np.asarray(x, float) for x in ref[pair])
        m = rr <= R2_RMAX
        calc = np.interp(rr[m], traj.r, g[ip])
        exp = gg[m]
        out[pair] = 1.0 - np.sum((calc - exp) ** 2) / np.sum((exp - exp.mean()) ** 2)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--every", type=int, default=2, help="one frame every N EPSR iterations")
    p.add_argument("--fps", type=int, default=12)
    a = p.parse_args()

    traj = Trajectory.load(ROOT / "example_H2O" / "data" / "trajectory")
    ref = water_reference()
    v = ConfigurationViewer(traj, off_screen=True, window_size=SIZE, gr="total", slab=10.0,
                            reference=ref, calc_opacity=CALC_OPACITY)
    exp_lines = [ln for ln in v.chart.plots("line") if ln.label.endswith("experiment")]
    calc_lines = [ln for ln in v.chart.plots("line") if ln.label.endswith("EPSR")]
    title = v.chart.title

    def set_fade(f):
        """f = 1: experiment solid, EPSR translucent; f = 0: experiment hidden, EPSR opaque."""
        v.chart.title = title if f >= 0.5 else title.split(" (")[0]
        for ln in exp_lines:
            ln.pen.SetOpacity(int(round(255 * f)))
        for ln in calc_lines:
            ln.pen.SetOpacity(int(round(255 * (CALC_OPACITY * f + (1 - f)))))

    # R^2 readout at the right of the g(r) panel, below the legend
    fs = int(round(10 * SIZE[1] / 820))
    v.pl.subplot(0, 1)
    v.pl.add_text("R² vs experiment", position=(0.76, 0.80), viewport=True, font_size=fs, color="black")
    r2_text = {
        pair: v.pl.add_text("", position=(0.76, 0.77 - 0.028 * ip), viewport=True, font_size=fs,
                            color=style.series(ip))
        for ip, pair in enumerate(traj.pairs)
    }
    v.pl.subplot(0, 0)

    def set_r2(step):
        for pair, val in r_squared(traj, ref, step).items():
            r2_text[pair].SetInput(f"{pair}:  {val:.4f}")

    steps = list(range(0, traj.nsteps, a.every))
    if steps[-1] != traj.nsteps - 1:
        steps.append(traj.nsteps - 1)
    iters = traj.history["step"]

    v.pl.open_movie(str(HERE / "H2O_config_movie_2.mp4"), framerate=a.fps, quality=9, macro_block_size=8)
    v.pl.subplot(0, 0)
    for s in steps:
        v.set_step(s)
        set_r2(s)
        set_fade(max(0.0, 1.0 - (int(iters[s]) - 1) / (FADE_OUT_ITERS - 1)))
        v.pl.camera.Azimuth(180.0 / len(steps))  # half a turn over the iterations
        v.pl.write_frame()
    for _ in range(int(HOLD * a.fps)):
        v.pl.write_frame()
    n = int(FADE_IN * a.fps)
    for f in np.linspace(0.0, 1.0, n + 1)[1:]:
        set_fade(0.5 - 0.5 * np.cos(np.pi * f))  # ease in/out
        v.pl.write_frame()
    for _ in range(int(HOLD_END * a.fps)):
        v.pl.write_frame()
    v.pl.close()
    print("wrote H2O_config_movie_2.mp4")


if __name__ == "__main__":
    main()
