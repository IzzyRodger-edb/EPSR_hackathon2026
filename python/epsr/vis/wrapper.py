"""Recorded EPSR runs with the wrapper: drop-in replacements for ``run`` and ``refine``.

``Simulation.run()`` and :func:`epsr.refine` run whole EPSR invocations in Rust,
so only the state after the last iteration is visible.  These functions run the
same iterations one at a time (``Simulation.iterations()``) and record a
snapshot after each one.  The simulation files are exactly those of
``run()``/``refine()``: recording only reads results.

    from epsr import vis
    recorder = vis.TrajectoryRecorder(case="water")
    vis.run_recorded(project.simulation(), recorder)
    history = vis.refine_recorded(project, 2, recorder, overrides={"niter": 100},
                                  initial_overrides={"ereq": 3.0})
    vis.make_report(recorder.trajectory(), "figures/")
"""

from __future__ import annotations

from collections.abc import Callable
from os import PathLike

from ..project import Project
from ..refine import History, _combine_overrides
from ..simulation import Overrides, Simulation, Summary
from .recorder import TrajectoryRecorder, record_run

__all__ = ["run_recorded", "refine_recorded"]


def run_recorded(
    simulation: Simulation,
    recorder: TrajectoryRecorder,
    *,
    iteration_callback: Callable[[Simulation, TrajectoryRecorder], None] | None = None,
) -> Summary:
    """``simulation.run()``, recording every iteration; returns the final summary.

    ``iteration_callback(simulation, recorder)`` runs after each iteration
    (e.g. :class:`epsr.vis.live.LiveMonitor`).  The requested EP amplitude
    ``ereq`` is taken from ``simulation.potential.params``.
    """
    record_run(simulation, recorder, callback=iteration_callback)
    return simulation.summary


def refine_recorded(
    project: Project | str | PathLike[str],
    n_runs: int,
    recorder: TrajectoryRecorder,
    *,
    overrides: Overrides = None,
    initial_overrides: Overrides = None,
    callback: Callable[[int, Simulation], bool | None] | None = None,
    iteration_callback: Callable[[Simulation, TrajectoryRecorder], None] | None = None,
    verbose: bool = False,
) -> History:
    """:func:`epsr.refine` with every iteration recorded into ``recorder``.

    Same arguments and the same :class:`epsr.History` result (one row per run);
    ``callback(run_number, simulation)`` after each run may return True to stop.
    """
    if not isinstance(project, Project):
        project = Project(project)
    history = History()
    for run_number in range(1, n_runs + 1):
        first_run_overrides = initial_overrides if run_number == 1 else None
        lines = _combine_overrides(overrides, first_run_overrides)
        simulation = project.simulation(lines, verbose=verbose)
        run_recorded(simulation, recorder, iteration_callback=iteration_callback)
        history.append(simulation.summary)
        if callback is not None and callback(run_number, simulation):
            break
    return history
