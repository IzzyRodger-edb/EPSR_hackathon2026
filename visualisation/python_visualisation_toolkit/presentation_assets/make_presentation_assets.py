"""Build the presentation assets: MP4 movies and still frames of the configuration + partial g(r).

Run from python_visualisation_toolkit/:

    uv run python presentation_assets/make_presentation_assets.py            # everything (~10 min)
    uv run python presentation_assets/make_presentation_assets.py --stills   # only the still images (~1 min)

Water is compared with the experimental partials of Soper (2000), from the ISIS
Disordered Materials Database (see experimental_data/SOURCE.md): experiment solid,
EPSR in the same colour but translucent. The database has no partials for
silica, so the silica movie shows the calculated g(r) only.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from epsrvis import Trajectory, load_reference_partials, plots
from epsrvis.viewer import ConfigurationViewer

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
EXP = HERE / "experimental_data"

SIZE = (1920, 1080)       # full-HD frames
CALC_OPACITY = 0.4        # translucency of the calculated g(r) where experiment is shown


def water_reference():
    # Soper's O/H labels -> the EPSR atom types of the water example
    return load_reference_partials({
        "OW-OW": EXP / "h2o298aoogor.txt",
        "OW-HW": EXP / "h2o298aohgor.txt",
        "HW-HW": EXP / "h2o298ahhgor.txt",
    })


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stills", action="store_true", help="only the still images, no movies")
    p.add_argument("--every", type=int, default=2, help="movie: one frame every N EPSR iterations")
    p.add_argument("--fps", type=int, default=12)
    a = p.parse_args()

    water = Trajectory.load(ROOT / "example_H2O" / "data" / "trajectory")
    silica = Trajectory.load(ROOT / "example_SiO2" / "data" / "trajectory")
    ref = water_reference()

    # static comparison figure (matplotlib)
    fig = plots.partials_vs_reference(water, ref, step=-1, total=True, calc_alpha=CALC_OPACITY,
                                      source="experiment (Soper 2000, 298 K)")
    fig.savefig(HERE / "H2O_partials_vs_experiment.png", dpi=200)
    plt.close(fig)
    print("wrote H2O_partials_vs_experiment.png")

    # still frames, same layout as the movies
    def viewer(traj, **kw):
        return ConfigurationViewer(traj, off_screen=True, window_size=SIZE, **kw)

    water_kw = dict(gr="total", slab=10.0, reference=ref, calc_opacity=CALC_OPACITY)
    silica_kw = dict(slab=12.0, color_by="coordination")
    for it in (100, 300):
        v = viewer(water, **water_kw)
        v.screenshot(HERE / f"H2O_config_iter{it:03d}.png", step=it - 1)
        v.pl.close()
        v = viewer(silica, **silica_kw)
        v.screenshot(HERE / f"SiO2_config_iter{it:03d}.png", step=it - 1)
        v.pl.close()
        print(f"wrote H2O_config_iter{it:03d}.png, SiO2_config_iter{it:03d}.png")
    if a.stills:
        return

    # movies (H.264 MP4); same content as the examples' config_movie.gif, at full HD
    viewer(water, **water_kw).movie(HERE / "H2O_config_movie.mp4", every=a.every, fps=a.fps)
    print("wrote H2O_config_movie.mp4")
    viewer(silica, **silica_kw).movie(HERE / "SiO2_config_movie.mp4", every=a.every, fps=a.fps)
    print("wrote SiO2_config_movie.mp4")


if __name__ == "__main__":
    main()
