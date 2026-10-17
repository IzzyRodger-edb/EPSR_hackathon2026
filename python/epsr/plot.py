"""Matplotlib helpers (optional extra ``plot``).

These replace the PGPLOT/gnuplot front ends of the EPSR suite, which are not ported.
Each function draws on ``ax`` (a new figure when omitted) and returns the axes.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from .io import _import_optional

if TYPE_CHECKING:
    from .refine import History
    from .simulation import PairKey, Simulation

__all__ = ["fits", "partials", "potentials", "history"]


def _axes(ax: Any) -> Any:
    if ax is not None:
        return ax
    plt = _import_optional("matplotlib.pyplot", "plot")
    _, ax = plt.subplots()
    return ax


def fits(sim: Simulation, datasets: Sequence[int] | None = None, *, offset: float = 1.0, ax: Any = None) -> Any:
    """Data (points) and simulated F(Q) (lines) of each dataset, offset vertically."""
    ax = _axes(ax)
    idx = range(sim.n_datasets) if datasets is None else datasets
    for k, i in enumerate(idx):
        f = sim.fit(i)
        line = ax.plot(f.q, f.simulated + k * offset, label=f.name)[0]
        ax.plot(f.q, f.data + k * offset, ".", ms=2, color=line.get_color())
    ax.set_xlabel("Q (Å$^{-1}$)")
    ax.set_ylabel("F(Q)")
    ax.legend()
    return ax


def partials(
    sim: Simulation,
    pairs: Sequence[PairKey] | None = None,
    *,
    kind: str = "gr_inter",
    offset: float = 1.0,
    ax: Any = None,
) -> Any:
    """Partial g(r) or S(Q) (``kind`` = ``gr_inter``, ``gr_intra``, ``sq_inter``, ``sq_intra``)."""
    ax = _axes(ax)
    data = getattr(sim, kind)()
    x = sim.r if kind.startswith("gr") else sim.q
    sel = range(sim.n_pairs) if pairs is None else [sim.pair_index(p) for p in pairs]
    for k, ip in enumerate(sel):
        a, b = sim.pair_labels[ip]
        ax.plot(x, data[ip] + k * offset, label=f"{a}-{b}")
    ax.set_xlabel("r (Å)" if kind.startswith("gr") else "Q (Å$^{-1}$)")
    ax.set_ylabel(kind.replace("_", " "))
    ax.legend(fontsize="small", ncol=2)
    return ax


def potentials(
    sim: Simulation,
    pairs: Sequence[PairKey] | None = None,
    *,
    which: str = "empirical",
    ax: Any = None,
) -> Any:
    """Empirical, reference or total potential of each pair (kJ/mol)."""
    ax = _axes(ax)
    pot = sim.potential
    data = getattr(pot, which)()
    sel = range(sim.n_pairs) if pairs is None else [sim.pair_index(p) for p in pairs]
    for ip in sel:
        a, b = sim.pair_labels[ip]
        ax.plot(pot.r, data[ip], label=f"{a}-{b}")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlabel("r (Å)")
    ax.set_ylabel(f"{which} potential (kJ/mol)")
    if which != "empirical":
        finite = data[np.isfinite(data)]
        if finite.size:
            ax.set_ylim(min(-1.0, float(np.percentile(finite, 1))), 10.0)
    ax.legend(fontsize="small", ncol=2)
    return ax


def history(hist: History, columns: Sequence[str] = ("r_factor", "energy", "ep_amplitude"), *, axes: Any = None) -> Any:
    """One panel per column of a :class:`~epsr.History` against run number."""
    plt = _import_optional("matplotlib.pyplot", "plot")
    cols = [c for c in columns if c in hist.columns]
    if axes is None:
        _, axes = plt.subplots(len(cols), 1, sharex=True, squeeze=False)
        axes = axes[:, 0]
    for ax, c in zip(axes, cols):
        ax.plot(hist.run, hist.column(c), "-o", ms=3)
        ax.set_ylabel(c.replace("_", " "))
    axes[-1].set_xlabel("run")
    return axes
