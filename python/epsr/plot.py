"""Matplotlib helpers (optional extra ``plot``).

These replace the PGPLOT/gnuplot front ends of the EPSR suite, which are not ported.
Each function draws on ``ax`` (a new figure when omitted) and returns the axes.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from .io import import_optional

if TYPE_CHECKING:
    from .refine import History
    from .simulation import PairKey, Simulation

__all__ = ["fits", "partials", "potentials", "history"]

Q_AXIS_LABEL = "Q (Å$^{-1}$)"
R_AXIS_LABEL = "r (Å)"
#: Upper limit of the y axis for reference and total potentials, which are huge at small r.
POTENTIAL_PLOT_CEILING = 10.0


def _pyplot() -> Any:
    return import_optional("matplotlib.pyplot", extra="plot")


def _axes_or_new(ax: Any) -> Any:
    """The given axes, or the axes of a new figure."""
    if ax is not None:
        return ax
    _, new_axes = _pyplot().subplots()
    return new_axes


def _selected_pairs(simulation: Simulation, pairs: Sequence[PairKey] | None) -> list[int]:
    """Indices of the requested pairs, or of all pairs."""
    if pairs is None:
        return list(range(simulation.n_pairs))
    return [simulation.pair_index(pair) for pair in pairs]


def _pair_label(simulation: Simulation, pair_index: int) -> str:
    first_type, second_type = simulation.pair_labels[pair_index]
    return f"{first_type}-{second_type}"


def fits(
    simulation: Simulation,
    datasets: Sequence[int] | None = None,
    *,
    offset: float = 1.0,
    ax: Any = None,
) -> Any:
    """Data (points) and simulated F(Q) (lines) of each dataset, offset vertically."""
    ax = _axes_or_new(ax)
    dataset_indices = range(simulation.n_datasets) if datasets is None else datasets
    for position, dataset in enumerate(dataset_indices):
        fit = simulation.fit(dataset)
        shift = position * offset
        (simulated_line,) = ax.plot(fit.q, fit.simulated + shift, label=fit.name)
        ax.plot(fit.q, fit.data + shift, ".", ms=2, color=simulated_line.get_color())
    ax.set_xlabel(Q_AXIS_LABEL)
    ax.set_ylabel("F(Q)")
    ax.legend()
    return ax


def partials(
    simulation: Simulation,
    pairs: Sequence[PairKey] | None = None,
    *,
    kind: str = "gr_inter",
    offset: float = 1.0,
    ax: Any = None,
) -> Any:
    """Partial g(r) or S(Q) (``kind`` = ``gr_inter``, ``gr_intra``, ``sq_inter``, ``sq_intra``)."""
    ax = _axes_or_new(ax)
    all_partials = getattr(simulation, kind)()
    is_real_space = kind.startswith("gr")
    x_values = simulation.r if is_real_space else simulation.q
    for position, pair_index in enumerate(_selected_pairs(simulation, pairs)):
        ax.plot(
            x_values,
            all_partials[pair_index] + position * offset,
            label=_pair_label(simulation, pair_index),
        )
    ax.set_xlabel(R_AXIS_LABEL if is_real_space else Q_AXIS_LABEL)
    ax.set_ylabel(kind.replace("_", " "))
    ax.legend(fontsize="small", ncol=2)
    return ax


def potentials(
    simulation: Simulation,
    pairs: Sequence[PairKey] | None = None,
    *,
    which: str = "empirical",
    ax: Any = None,
) -> Any:
    """Empirical, reference or total potential of each pair (kJ/mol)."""
    ax = _axes_or_new(ax)
    potential = simulation.potential
    all_potentials = getattr(potential, which)()
    for pair_index in _selected_pairs(simulation, pairs):
        ax.plot(potential.r, all_potentials[pair_index], label=_pair_label(simulation, pair_index))
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlabel(R_AXIS_LABEL)
    ax.set_ylabel(f"{which} potential (kJ/mol)")
    if which != "empirical":
        finite_values = all_potentials[np.isfinite(all_potentials)]
        if finite_values.size:
            lowest_shown = min(-1.0, float(np.percentile(finite_values, 1)))
            ax.set_ylim(lowest_shown, POTENTIAL_PLOT_CEILING)
    ax.legend(fontsize="small", ncol=2)
    return ax


def history(
    run_history: History,
    columns: Sequence[str] = ("r_factor", "energy", "ep_amplitude"),
    *,
    axes: Any = None,
) -> Any:
    """One panel per column of a :class:`~epsr.History` against run number."""
    columns_to_plot = [name for name in columns if name in run_history.columns]
    if axes is None:
        _, axes_grid = _pyplot().subplots(len(columns_to_plot), 1, sharex=True, squeeze=False)
        axes = axes_grid[:, 0]
    for ax, column_name in zip(axes, columns_to_plot, strict=True):
        ax.plot(run_history.run, run_history.column(column_name), "-o", ms=3)
        ax.set_ylabel(column_name.replace("_", " "))
    axes[-1].set_xlabel("run")
    return axes
