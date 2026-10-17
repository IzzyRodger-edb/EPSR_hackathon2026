"""Build the explanatory fluid assets: an example fluid configuration and two cartoons.

Run from python_visualisation_toolkit/:

    uv run python presentation_assets/make_fluid_cartoons.py            # re-render (reuses the simulation, ~15 s)
    uv run python presentation_assets/make_fluid_cartoons.py --rerun    # rerun the simulation too (~1 min more)

Writes, next to this script:

    fluid_configuration.xyz      atom positions of a 2D Lennard-Jones fluid (Metropolis Monte Carlo)
    fluid_configuration_gr.dat   its simulated g(r)
    fluid_configuration.png      a plain rendering of that configuration
    fluid_gr_cartoon.png         cartoon of a typical fluid pair distribution function g(r)
    fluid_config_cartoon.png     cartoon of the configuration around one atom, with its
                                 coordination shells drawn at the g(r) peak positions
    fluid_jiggle.gif             the configuration's atoms jiggling (the Monte Carlo continued,
                                 one sweep per frame)

The fluid is 2D so that every atom is visible on a slide. Units are reduced
Lennard-Jones units (sigma = 1); the cartoon g(r) uses the peak positions
measured from the simulation, so the two cartoons line up.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Wedge

HERE = Path(__file__).resolve().parent

N = 400            # atoms
RHO = 0.70         # reduced number density (atoms / sigma^2)
T = 1.0            # reduced temperature: a dense 2D LJ liquid
RCUT = 3.0
SWEEPS_EQ = 1500
SWEEPS_PROD = 500
SEED = 7

ATOM = "#2a78d6"       # categorical slot 1 (blue)
CENTRE = "#eb6834"     # slot 2 (orange): the reference atom
SHELL = "#1baf7a"      # slot 3 (aqua): first-shell neighbours
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#d9d8d4"
SURFACE = "#fcfcfb"


# ---------------------------------------------------------------- simulation

def lj_energy(dr2):
    inv6 = (1.0 / dr2) ** 3
    return np.where(dr2 < RCUT ** 2, 4.0 * (inv6 * inv6 - inv6), 0.0)


def atom_energy(pos, i, ri, L):
    d = pos - ri
    d -= L * np.round(d / L)
    dr2 = np.einsum("ij,ij->i", d, d)
    dr2[i] = np.inf
    return lj_energy(dr2).sum()


def mc_sweep(pos, L, step, rng):
    """One Metropolis sweep (N single-atom trial moves), in place; returns the number accepted."""
    acc = 0
    for i in rng.permutation(len(pos)):
        new = (pos[i] + rng.uniform(-step, step, 2)) % L
        dE = atom_energy(pos, i, new, L) - atom_energy(pos, i, pos[i], L)
        if dE < 0 or rng.random() < np.exp(-dE / T):
            pos[i] = new
            acc += 1
    return acc


def simulate(rng):
    L = np.sqrt(N / RHO)
    n = int(np.ceil(np.sqrt(N)))
    g = (np.arange(n) + 0.5) * L / n
    pos = np.array([(x, y) for x in g for y in g])[:N].copy()
    step = 0.15
    hist_r = np.linspace(0, L / 2, 301)
    hist = np.zeros(len(hist_r) - 1)
    nsamp = 0
    for sweep in range(SWEEPS_EQ + SWEEPS_PROD):
        acc = mc_sweep(pos, L, step, rng)
        if sweep < SWEEPS_EQ:                       # tune to ~40 % acceptance
            step *= 1.05 if acc / N > 0.4 else 0.95
        elif sweep % 5 == 0:
            d = pos[:, None, :] - pos[None, :, :]
            d -= L * np.round(d / L)
            r = np.sqrt((d ** 2).sum(-1))[np.triu_indices(N, 1)]
            hist += np.histogram(r, hist_r)[0]
            nsamp += 1
    rc = 0.5 * (hist_r[1:] + hist_r[:-1])
    shell = np.pi * (hist_r[1:] ** 2 - hist_r[:-1] ** 2)
    gr = 2 * hist / (nsamp * N * RHO * shell)
    return pos, L, rc, gr


def peaks(r, g, n=3):
    gs = np.convolve(g, np.ones(5) / 5, mode="same")
    idx = [i for i in range(2, len(gs) - 2) if gs[i] > gs[i - 1] and gs[i] >= gs[i + 1] and gs[i] > 1.02]
    out = []
    for i in idx:                                   # drop noise wiggles within 0.6 sigma of a kept peak
        if not out or r[i] - out[-1][0] > 0.6:
            out.append((r[i], gs[i]))
        elif gs[i] > out[-1][1]:
            out[-1] = (r[i], gs[i])
    return [p[0] for p in out[:n]]


# ---------------------------------------------------------------- outputs

def write_xyz(pos, L, path):
    with open(path, "w") as f:
        f.write(f"{len(pos)}\n")
        f.write(f'Lattice="{L:.6f} 0 0 0 {L:.6f} 0 0 0 1" Properties=species:S:1:pos:R:3 '
                f"2D Lennard-Jones fluid, rho*={RHO}, T*={T}, sigma=1\n")
        for x, y in pos:
            f.write(f"Ar {x:12.6f} {y:12.6f} {0.0:12.6f}\n")


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ax.spines.values():
        s.set_visible(False)


def plot_configuration(pos, L, path):
    fig, ax = plt.subplots(figsize=(8, 8), facecolor=SURFACE)
    style(ax)
    box = plt.Rectangle((0, 0), L, L, fill=False, ec=INK2, lw=1.2)
    ax.add_patch(box)
    for x, y in pos:
        for dx in (-L, 0, L):                       # periodic images fill the edges, clipped to the box
            for dy in (-L, 0, L):
                ax.add_patch(Circle((x + dx, y + dy), 0.5, fc=ATOM, ec=SURFACE, lw=1.0)).set_clip_path(box)
    ax.set_xlim(-0.6, L + 0.6)
    ax.set_ylim(-0.6, L + 0.6)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"A fluid: {len(pos)} atoms, no long-range order", color=INK, fontsize=15, pad=12)
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)


def cartoon_gr(r, r1, d, amp=2.0, decay=0.5, width=0.09):
    """Idealised fluid g(r): zero inside the core, a sharp first peak, decaying oscillations to 1."""
    rising = (1 + amp) * np.exp(-((r - r1) / width) ** 2)
    osc = 1 + amp * np.exp(-(r - r1) / decay) * np.cos(2 * np.pi * (r - r1) / d)
    return np.where(r < r1, rising, osc)


def plot_gr_cartoon(shells, path):
    r1 = shells[0]
    d = shells[1] - shells[0] if len(shells) > 1 else 0.95 * r1
    r = np.linspace(0, 5.0, 2000)
    g = cartoon_gr(r, r1, d)

    fig, ax = plt.subplots(figsize=(9, 5.2), facecolor=SURFACE)
    style(ax)
    ax.axhline(1, color=GRID, lw=1.2, zorder=0)
    ax.axvspan(0, 0.9, color=GRID, alpha=0.45, lw=0, zorder=0)
    ax.plot(r, g, color=ATOM, lw=3, solid_capstyle="round")

    kw = dict(color=INK, fontsize=12, ha="center")
    ax.text(0.45, 1.55, "excluded\ncore\n(atoms can't\noverlap)", color=INK2, fontsize=11, ha="center", va="center")
    ax.annotate("1st shell:\nnearest neighbours", (r1, g[np.argmin(abs(r - r1))]), (r1 + 0.95, 2.75),
                arrowprops=dict(arrowstyle="-", color=INK2, lw=1), va="center", ha="left", fontsize=12, color=INK)
    if len(shells) > 1:
        r2 = shells[1]
        ax.annotate("2nd shell", (r2, g[np.argmin(abs(r - r2))]), (r2 + 0.35, 1.7),
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=1), va="center", ha="left", fontsize=12, color=INK)
    ax.text(4.55, 1.12, "g(r) → 1\nno order at long range", va="bottom", **kw)

    ax.set_xlim(0, 5.0)
    ax.set_ylim(0, 3.2)
    ax.set_xlabel("r  (atomic diameters)", color=INK2, fontsize=13)
    ax.set_ylabel("g(r)", color=INK2, fontsize=13)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_xticks(range(6))
    ax.tick_params(colors=INK2, length=0, labelsize=11)
    for s in ("left", "bottom"):
        ax.spines[s].set_visible(True)
        ax.spines[s].set_color(INK2)
    ax.set_title("Pair distribution function of a fluid", color=INK, fontsize=15, pad=12)
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)


def plot_config_cartoon(pos, L, shells, path):
    # the atom nearest the box centre is the reference atom; positions relative to it (minimum image)
    i0 = np.argmin(((pos - L / 2) ** 2).sum(1))
    d = pos - pos[i0]
    d -= L * np.round(d / L)
    rr = np.sqrt((d ** 2).sum(1))
    r1, r2 = shells[0], shells[1]
    first_min = 0.5 * (r1 + r2)
    R = 3.4                                         # radius of the circular window
    dr = 0.35

    fig, ax = plt.subplots(figsize=(11, 7.6), facecolor=SURFACE)
    style(ax)
    window = Circle((0, 0), R, fc="none", ec=GRID, lw=1.5, zorder=8)
    ax.add_patch(window)
    ax.add_patch(Wedge((0, 0), r2 + dr / 2, 0, 360, width=dr, fc=GRID, ec="none", alpha=0.7, zorder=1))
    for (x, y), r in zip(d, rr):
        if r > R + 0.6:
            continue
        col = CENTRE if r == 0 else SHELL if r < first_min else ATOM
        ax.add_patch(Circle((x, y), 0.46, fc=col, ec=SURFACE, lw=2, zorder=3, clip_path=window))
    for rs in (r1, r2):
        ax.add_patch(Circle((0, 0), rs, fill=False, ec=INK2, lw=1.6, ls=(0, (5, 4)), zorder=4))

    ang = np.radians(-35)
    ax.annotate("", (r2 * np.cos(ang), r2 * np.sin(ang)), (0, 0), zorder=7,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=2, mutation_scale=18))
    ax.text(0.5 * r2 * np.cos(ang) + 0.3 * np.sin(-ang), 0.5 * r2 * np.sin(ang) + 0.3 * np.cos(ang), "r",
            fontsize=20, style="italic", color=INK, zorder=7, ha="center", va="center",
            bbox=dict(fc=SURFACE, ec="none", pad=1))

    # key on the right
    x0, y = R + 1.0, 2.3
    items = [
        ("disc", CENTRE, "reference atom"),
        ("disc", SHELL, f"1st shell: nearest neighbours\n(the first g(r) peak, r ≈ {r1:.1f} diameters)"),
        ("ring", INK2, f"2nd shell (second peak, r ≈ {r2:.1f})"),
        ("band", GRID, "thin shell r → r + dr:\ng(r) counts the atoms in it"),
    ]
    for kind, col, label in items:
        if kind == "disc":
            ax.add_patch(Circle((x0, y), 0.3, fc=col, ec="none"))
        elif kind == "ring":
            ax.plot([x0 - 0.35, x0 + 0.35], [y, y], color=col, lw=1.6, ls=(0, (5, 4)))
        else:
            ax.add_patch(plt.Rectangle((x0 - 0.35, y - 0.15), 0.7, 0.3, fc=col, ec="none", alpha=0.7))
        ax.text(x0 + 0.65, y, label, va="center", ha="left", fontsize=13, color=INK, linespacing=1.4)
        y -= 1.35

    ax.text(-R, -R - 0.55,
            "g(r): how many atoms sit at distance r from a typical atom,\n"
            "relative to a structureless fluid of the same density",
            ha="left", va="top", color=INK2, fontsize=12)
    ax.set_xlim(-R - 0.3, R + 8.6)
    ax.set_ylim(-R - 1.7, R + 0.3)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("A fluid seen from one atom", color=INK, fontsize=16, pad=10, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)


def make_jiggle_gif(pos, L, path, frames=120, fps=15, step=0.2, rng=None):
    """Continue the Monte Carlo from `pos`, one sweep per frame, and save the atoms jiggling as a GIF."""
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.collections import EllipseCollection

    rng = rng or np.random.default_rng(SEED + 1)
    pos = pos.copy()
    images = np.array([(dx, dy) for dx in (-L, 0, L) for dy in (-L, 0, L)])

    def tiled(p):                                   # periodic images fill the edges, clipped to the box
        return (p[None, :, :] + images[:, None, :]).reshape(-1, 2)

    fig, ax = plt.subplots(figsize=(6.4, 6.4), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for sp in ax.spines.values():                   # the axes frame is the box; atoms clip to it
        sp.set_color(INK2)
        sp.set_linewidth(1.2)
    atoms = EllipseCollection(1.0, 1.0, 0.0, units="xy", offsets=tiled(pos), offset_transform=ax.transData,
                              facecolors=ATOM, edgecolors=SURFACE, linewidths=0.8)
    ax.add_collection(atoms)
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("A fluid: atoms jiggle and wander, with no long-range order", color=INK, fontsize=12, pad=8)
    fig.tight_layout()

    def update(_):
        mc_sweep(pos, L, step, rng)
        atoms.set_offsets(tiled(pos))
        return (atoms,)

    FuncAnimation(fig, update, frames=frames, blit=False).save(path, writer=PillowWriter(fps=fps), dpi=100,
                                                              savefig_kwargs={"facecolor": SURFACE})
    plt.close(fig)


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    L = float(lines[1].split('"')[1].split()[0])
    pos = np.array([[float(v) for v in ln.split()[1:3]] for ln in lines[2:2 + int(lines[0])]])
    return pos, L


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rerun", action="store_true", help="rerun the simulation even if its outputs exist")
    a = p.parse_args()

    xyz, dat = HERE / "fluid_configuration.xyz", HERE / "fluid_configuration_gr.dat"
    if a.rerun or not (xyz.exists() and dat.exists()):
        pos, L, r, g = simulate(np.random.default_rng(SEED))
    else:                                           # re-render only (~seconds)
        pos, L = read_xyz(xyz)
        r, g = np.loadtxt(dat, unpack=True)
    shells = peaks(r, g)
    print(f"box {L:.2f} sigma, g(r) peaks at " + ", ".join(f"{p:.2f}" for p in shells))

    write_xyz(pos, L, xyz)
    np.savetxt(dat, np.c_[r, g], fmt="%.5f",
               header="r (sigma)  g(r)   simulated, averaged over the production sweeps")
    plot_configuration(pos, L, HERE / "fluid_configuration.png")
    plot_gr_cartoon(shells, HERE / "fluid_gr_cartoon.png")
    plot_config_cartoon(pos, L, shells, HERE / "fluid_config_cartoon.png")
    make_jiggle_gif(pos, L, HERE / "fluid_jiggle.gif")
    for f in ("fluid_configuration.xyz", "fluid_configuration_gr.dat", "fluid_configuration.png",
              "fluid_gr_cartoon.png", "fluid_config_cartoon.png", "fluid_jiggle.gif"):
        print("wrote", f)


if __name__ == "__main__":
    main()
