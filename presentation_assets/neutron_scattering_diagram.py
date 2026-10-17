"""Vertical explainer graphic: how a neutron scattering experiment works.

Renders neutron_scattering_diagram.{png,svg} next to this script, sized to sit
in the left-hand gap of the "Motivation" slide.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle, Wedge

plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "svg.fonttype": "none"})

NAVY = "#0d2233"
BLUE = "#2f7ed8"
NEUTRON = "#e8772e"
GREY = "#8a96a3"

fig = plt.figure(figsize=(3.3, 4.8))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 10)
ax.set_ylim(2.6, 17)
ax.set_aspect("equal")
ax.axis("off")

CX = 2.6  # centre line of the icon column
TX = 5.4  # left edge of the caption column


def caption(y, n, title, body):
    ax.text(TX, y + 0.35, f"{n}. {title}", fontsize=11, fontweight="bold", color=NAVY, va="bottom")
    ax.text(TX, y + 0.2, body, fontsize=9, color=NAVY, va="top", linespacing=1.15)


def down_arrow(y0, y1):
    ax.add_patch(FancyArrowPatch((CX, y0), (CX, y1), arrowstyle="-|>", mutation_scale=14, lw=2, color=NAVY))


# 1. Source: accelerator ring + target
ax.add_patch(Circle((CX, 15.5), 1.0, fill=False, lw=3, ec=GREY))
ax.add_patch(Circle((CX + 1.0, 15.5), 0.13, color="#c0392b", zorder=3))  # proton bunch
ax.annotate("", xy=(CX, 14.05), xytext=(CX + 1.0, 15.35),
            arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.5, connectionstyle="arc3,rad=-0.4"))
ax.add_patch(Rectangle((CX - 0.55, 13.55), 1.1, 0.45, fc="#5b6672", ec=NAVY, lw=1.2, zorder=3))
caption(15.2, 1, "Neutron source", "Protons from an\naccelerator smash into\na metal target,\nknocking out neutrons")

# neutron beam
for y in np.linspace(12.85, 11.65, 4):
    ax.add_patch(Circle((CX, y), 0.14, color=NEUTRON, zorder=3))
down_arrow(13.4, 10.75)

# 2. Sample with atoms, scattered neutrons, detector arc
rng = np.random.default_rng(4)
ax.add_patch(FancyBboxPatch((CX - 0.75, 9.0), 1.5, 1.6, boxstyle="round,pad=0.05,rounding_size=0.15",
                            fc="#eaf2fb", ec=NAVY, lw=1.5, zorder=2))
# liquid-like atoms: a jittered grid, so neighbours never overlap
gx, gy = np.meshgrid(np.linspace(CX - 0.48, CX + 0.48, 3), np.linspace(9.25, 10.35, 3))
pts = np.column_stack([gx.ravel(), gy.ravel()]) + rng.uniform(-0.09, 0.09, (9, 2))
for p in pts:
    ax.add_patch(Circle(p, 0.16, color=BLUE, zorder=3))

R = 2.15
for ang_deg in (-165, -135, -110, -70, -45, -15):
    a = np.deg2rad(ang_deg)
    x1, y1 = CX + 0.9 * np.cos(a), 9.8 + 0.9 * np.sin(a)
    x2, y2 = CX + (R - 0.25) * np.cos(a), 9.8 + (R - 0.25) * np.sin(a)
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=8,
                                 lw=1.2, ls=(0, (3, 2)), color=NEUTRON, zorder=1))
for start in range(-175, -5, 12):  # segmented detector bank
    ax.add_patch(Wedge((CX, 9.8), R + 0.3, start, start + 10, width=0.3, fc=GREY, ec=NAVY, lw=0.6))
caption(10.4, 2, "Sample", "Neutrons bounce off\nthe atoms' nuclei")
caption(8.05, 3, "Detectors", "Count how many\nneutrons scatter\nat each angle")

down_arrow(7.3, 6.05)

# 4. Measured pattern
ax.add_patch(FancyBboxPatch((CX - 1.9, 3.35), 3.8, 2.5, boxstyle="round,pad=0.05,rounding_size=0.2",
                            fc="white", ec=NAVY, lw=1.5))
q = np.linspace(0.3, 12, 300)
s = 1 - 0.85 * np.exp(-((q / 1.8) ** 2)) + 0.9 * np.sin(1.9 * q - 3.2) / (0.9 * q + 0.2) * (1 - np.exp(-((q / 2.2) ** 4)))
xs = CX - 1.6 + (q - q.min()) / (q.max() - q.min()) * 3.2
ys = 3.65 + (s - s.min()) / (s.max() - s.min()) * 1.9
ax.plot(xs, ys, color=BLUE, lw=1.8)
ax.plot([CX - 1.6, CX + 1.6], [3.6, 3.6], color=NAVY, lw=0.8)
ax.text(CX, 3.35 - 0.2, "scattering angle", fontsize=8, color=NAVY, ha="center", va="top", style="italic")
caption(5.1, 4, "Scattering pattern", "Its ripples encode\nhow far apart the\natoms sit")

out = Path(__file__).with_name("neutron_scattering_diagram")
for ext, kw in ((".png", {"dpi": 300}), (".svg", {})):
    fig.savefig(out.with_suffix(ext), transparent=True, bbox_inches="tight", pad_inches=0.05, **kw)
print(out)
