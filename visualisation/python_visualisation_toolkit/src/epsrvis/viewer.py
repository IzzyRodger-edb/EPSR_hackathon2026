"""3D view of the EPSR configuration (pyvista) side by side with the partial g(r).

Left: atoms as spheres in CPK colours, bonds, the simulation box.
Right: a 2D chart of the partial pair distribution functions of the *same*
iteration.  A slider steps through the recorded EPSR iterations.

    from epsrvis import Trajectory
    from epsrvis.viewer import ConfigurationViewer
    v = ConfigurationViewer(Trajectory.load("example_H2O/data/trajectory"))
    v.show()                                  # interactive window with a slider
    v.screenshot("water.png", step=-1)        # off-screen
    v.movie("water.gif", every=5)             # off-screen animation
"""

from __future__ import annotations

import numpy as np
import pyvista as pv

from . import metrics as M
from . import style
from .systems import get_case
from .trajectory import Trajectory

BACKGROUND = "#262624"
CHART_BG = "#fcfcfb"


class ConfigurationViewer:
    """Linked 3D configuration + partial g(r) view of a recorded trajectory.

    Parameters
    ----------
    gr : "inter" | "total"
        Intermolecular g(r) only (as in ``.EPSR.g01``) or including the
        intramolecular part (adds e.g. the O-H bond peak of water).
    slab : float or None
        Only draw atoms with |z| < slab/2 (Å), to see inside the box.
    atom_scale : float
        Sphere radius as a fraction of the element's display radius (0.6 gives
        ball-and-stick, 1.0 space-filling-like).
    color_by : "type" | "coordination"
        "coordination" colours the centre atoms of the case's first coordination
        analysis (e.g. Si by number of O neighbours) and greys out the others.
    """

    def __init__(self, traj: Trajectory, gr: str = "inter", slab: float | None = None,
                 bonds: bool = True, color_by: str = "type", off_screen: bool = False,
                 atom_scale: float = 0.6,
                 window_size=(1600, 820)):
        self.traj = traj
        self.case = get_case(traj.case, traj.atom_types)
        self.gr_mode = gr
        self.slab = slab
        self.draw_bonds = bonds
        self.color_by = color_by
        self.labels = traj.atom_labels
        self.mol = traj.molecule_of_atom
        self.step = traj.nsteps - 1
        self.pl = pv.Plotter(shape=(1, 2), off_screen=off_screen, window_size=window_size,
                             border=False)
        self._sphere = {t: pv.Sphere(radius=atom_scale * self.case.radius(t), theta_resolution=14, phi_resolution=10)
                        for t in traj.atom_types}
        self._atom_meshes: dict[str, pv.PolyData] = {}
        self._bond_mesh: pv.PolyData | None = None
        self._lines = []
        self._build()

    # ---------------------------------------------------------------- geometry
    def _visible(self, pos: np.ndarray) -> np.ndarray:
        if self.slab is None:
            return np.ones(len(pos), dtype=bool)
        return np.abs(pos[:, 2]) < 0.5 * self.slab

    def _atoms(self, pos: np.ndarray, vis: np.ndarray, t: str) -> pv.PolyData:
        sel = (self.labels == t) & vis
        pts = pv.PolyData(pos[sel] if sel.any() else np.zeros((0, 3)))
        if self.color_by == "coordination" and self.case.coordination:
            a, b, cut, inter, _ = self.case.coordination[0]
            if t == a and sel.any():
                pts["neighbours"] = self._coordination(pos, a, b, cut, inter)[sel]
        if pts.n_points == 0:
            return pv.PolyData()
        return pts.glyph(geom=self._sphere[t], scale=False, orient=False)

    def _coordination(self, pos, a, b, cut, inter):
        i, j, _ = M.neighbour_pairs(pos, self.traj.box(self.step), self.labels, a, b, cut)
        if inter:
            keep = self.mol[i] != self.mol[j]
            i = i[keep]
        return np.bincount(i, minlength=len(pos)).astype(float)

    def _bonds(self, pos: np.ndarray, vis: np.ndarray) -> pv.PolyData:
        box = self.traj.box(self.step)
        segs = []
        for a, b, cut, intra in self.case.bonds:
            i, j, _ = M.neighbour_pairs(pos, box, self.labels, a, b, cut)
            keep = vis[i] & vis[j]
            if intra:
                keep &= self.mol[i] == self.mol[j]
            # only bonds that do not cross the periodic boundary (direct distance < cutoff)
            direct = np.linalg.norm(pos[i] - pos[j], axis=1)
            keep &= direct < cut
            if a == b:
                keep &= i < j
            segs.append(np.stack([i[keep], j[keep]], axis=1))
        if not segs or sum(len(s) for s in segs) == 0:
            return pv.PolyData()
        seg = np.concatenate(segs)
        lines = np.column_stack([np.full(len(seg), 2), seg]).ravel()
        return pv.PolyData(pos.astype(float), lines=lines)

    def _gr(self, step: int) -> np.ndarray:
        return self.traj.gr_total(step) if self.gr_mode == "total" else np.asarray(self.traj.gr_inter[step])

    # ---------------------------------------------------------------- scene
    def _build(self):
        pl, traj = self.pl, self.traj
        pos = np.asarray(traj.positions[self.step], dtype=float)
        vis = self._visible(pos)
        pl.subplot(0, 0)
        pl.set_background(BACKGROUND, all_renderers=False)
        for t in traj.atom_types:
            mesh = self._atoms(pos, vis, t)
            self._atom_meshes[t] = mesh
            if "neighbours" in mesh.array_names:
                pl.add_mesh(mesh, scalars="neighbours", cmap="viridis", clim=[-0.5, 6.5], n_colors=7,
                            annotations={float(i): str(i) for i in range(7)},
                            smooth_shading=True,
                            scalar_bar_args=dict(title=f"{self.case.coordination[0][1]} around {t}", color="white",
                                                 n_labels=0, vertical=True, position_x=0.9,
                                                 position_y=0.25, width=0.04, height=0.45))
            else:
                op = 0.35 if (self.color_by == "coordination" and self.case.coordination
                              and t != self.case.coordination[0][0]) else 1.0
                pl.add_mesh(mesh, color=self.case.color(t), smooth_shading=True, opacity=op,
                            specular=0.3, name=f"atoms_{t}")
        if self.draw_bonds and self.case.bonds:
            self._bond_mesh = self._bonds(pos, vis)
            pl.add_mesh(self._bond_mesh, color="#bdbcb6", line_width=2.5, render_lines_as_tubes=True,
                        name="bonds")
        L = traj.box(self.step)
        pl.add_mesh(pv.Box(bounds=(-L / 2, L / 2) * 3).outline(), color="#c3c2b7", line_width=1.5)
        legend = [[f"{t} atoms" if t == self.case.element[t] else f"{t} ({self.case.element[t]})", self.case.color(t)]
                  for t in traj.atom_types]
        pl.add_legend(legend, bcolor=BACKGROUND, border=False, size=(0.13, 0.045 * len(legend)), loc="upper right",
                      face="circle")
        self._title = pl.add_text(self._title_text(), position="lower_left", font_size=10, color="white")
        pl.camera_position = "iso"
        pl.camera.zoom(1.15)

        pl.subplot(0, 1)
        pl.set_background(CHART_BG, all_renderers=False)
        chart = pv.Chart2D()
        g = self._gr(self.step)
        r = traj.r
        for ip, pair in enumerate(traj.pairs):
            self._lines.append(chart.line(r, g[ip], color=style.series(ip), width=2.5, label=pair))
        chart.line([r[0], r[-1]], [1, 1], color="#8a8985", width=1.0, style="--")
        chart.x_label = "r (Å)"
        chart.y_label = "g(r)" + (" (inter + intra)" if self.gr_mode == "total" else " (intermolecular)")
        chart.title = "Partial pair distribution functions"
        chart.x_range = [0, min(10.0, float(r[-1]))]
        gmax = max(float(np.nanmax(self._gr(s)[:, r > 0.5])) for s in (0, traj.nsteps - 1))
        chart.y_range = [-0.1, min(gmax * 1.1, 6.0 if self.gr_mode == "total" else gmax * 1.1)]
        self.chart = chart
        pl.add_chart(chart)
        pl.subplot(0, 0)

    def _title_text(self) -> str:
        h = self.traj.history
        s = self.step
        fit = getattr(self, "_fit", None)
        if fit is None:
            self._fit = fit = M.fit_metrics(self.traj)
        stage = self.case.stages.get(int(h["run"][s]), f"run {int(h['run'][s])}")
        return (f"{self.traj.case}: EPSR iteration {int(h['step'][s])} ({stage}, ereq {h['ereq'][s]:g})\n"
                f"E_inter {h['energd'][s]:.2f} kJ/mol, P {h['pressd'][s]:.2f} kbar\n"
                f"R-factor {fit['rfactor_tot'][s]:.2e}, EP amplitude {h['energabs'][s]:.2f} kJ/mol")

    def set_step(self, step: int):
        """Show recorded step index ``step`` (negative counts from the end)."""
        traj = self.traj
        step = int(step) % traj.nsteps
        self.step = step
        pos = np.asarray(traj.positions[step], dtype=float)
        vis = self._visible(pos)
        for t, mesh in self._atom_meshes.items():
            new = self._atoms(pos, vis, t)
            mesh.copy_from(new)
        if self._bond_mesh is not None:
            self._bond_mesh.copy_from(self._bonds(pos, vis))
        g = self._gr(step)
        for ip, line in enumerate(self._lines):
            line.update(traj.r, g[ip])
        self._title.SetText(0, self._title_text())

    # ---------------------------------------------------------------- output
    def show(self):
        """Interactive window; the slider selects the EPSR iteration."""
        self.pl.subplot(0, 0)
        self.pl.add_slider_widget(
            lambda v: self.set_step(int(round(v)) - 1), rng=[1, self.traj.nsteps], value=self.step + 1,
            title="EPSR iteration", pointa=(0.1, 0.93), pointb=(0.6, 0.93), style="modern", fmt="%.0f",
            interaction_event="always", color="white", title_color="white")
        self.pl.show()

    def screenshot(self, path, step: int = -1):
        self.set_step(step)
        self.pl.render()
        self.pl.screenshot(str(path))

    def movie(self, path, every: int = 5, fps: int = 10, steps=None):
        """GIF (or .mp4 with imageio-ffmpeg) of the configuration and g(r) through the refinement."""
        steps = list(steps) if steps is not None else list(range(0, self.traj.nsteps, every))
        if steps[-1] != self.traj.nsteps - 1:
            steps.append(self.traj.nsteps - 1)
        path = str(path)
        if path.endswith(".gif"):
            self.pl.open_gif(path, fps=fps)
        else:
            self.pl.open_movie(path, framerate=fps)
        self.pl.subplot(0, 0)
        for s in steps:
            self.set_step(s)
            self.pl.camera.Azimuth(180.0 / len(steps))  # half a turn over the movie
            self.pl.write_frame()
        self.pl.close()
