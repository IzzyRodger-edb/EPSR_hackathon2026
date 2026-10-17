"""Live plots while the Python wrapper runs EPSR.

Pass a :class:`LiveMonitor` as the ``callback`` of
:func:`epsrvis.recorder.record_run`::

    mon = LiveMonitor(every=5)
    record_run(sim, rec, pcof=..., callback=mon)

Every ``every`` iterations it redraws the R-factor history and the current
F(Q) fit of each data set (in a Jupyter notebook use ``%matplotlib widget``
or ``%matplotlib inline`` with ``clear_output=True``).
"""

from __future__ import annotations

import matplotlib.pyplot as plt

from . import style


class LiveMonitor:
    def __init__(self, every: int = 5, clear_output: bool = False):
        self.every = every
        self.clear_output = clear_output
        self.fig = None

    def __call__(self, sim, recorder) -> None:
        if recorder.nsteps % self.every:
            return
        traj = recorder.trajectory()
        self.draw(traj)

    def draw(self, traj) -> None:
        from .metrics import fit_metrics

        style.apply()
        nd = len(traj.datasets)
        if self.fig is None:
            self.fig, self.axes = plt.subplots(
                1, nd + 1, figsize=(4.2 * (nd + 1), 3.4), squeeze=False
            )
        axes = self.axes[0]
        for a in axes:
            a.cla()
        fit = fit_metrics(traj)
        x = traj.steps
        for d, name in enumerate(traj.datasets):
            axes[0].plot(x, fit["rfactor"][:, d], color=style.series(d), label=name)
            sl = traj.q_slice(d)
            a = axes[d + 1]
            a.plot(
                traj.q[sl],
                traj.fq_data[-1, d, sl],
                "o",
                ms=1.8,
                color=style.INK,
                label="data",
            )
            a.plot(
                traj.q[sl], traj.fq_sim[-1, d, sl], color=style.series(d), label="EPSR"
            )
            a.set_title(
                f"{name}: R = {fit['rfactor'][-1, d]:.2e}, R² = {fit['r2'][-1, d]:.3f}",
                fontsize=9,
            )
            a.set_xlabel("Q (Å⁻¹)")
        axes[0].set_yscale("log")
        axes[0].set_title(f"R-factor, iteration {int(x[-1])}")
        axes[0].set_xlabel("EPSR iteration")
        axes[0].legend(fontsize=7)
        self.fig.tight_layout()
        if self.clear_output:
            from IPython.display import clear_output, display

            clear_output(wait=True)
            display(self.fig)
        else:
            self.fig.canvas.draw_idle()
            plt.pause(0.001)
