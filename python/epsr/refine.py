"""Repeated EPSR runs: :func:`refine`, :func:`refine_many` and :class:`History`."""

from __future__ import annotations

import sys
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from os import PathLike
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

from .io import ERG_COLUMNS, import_optional, read_erg
from .project import Project
from .simulation import Overrides, Simulation, Summary, render_overrides

__all__ = ["History", "refine", "refine_many"]

#: Summary fields recorded for every run (in addition to ``run`` and ``r_factor_per_dataset``).
SUMMARY_COLUMNS = (
    "iteration",
    "energy",
    "energy_inter",
    "energy_intra",
    "pressure",
    "ep_amplitude",
    "r_factor",
    "r_factor_raw",
    "ep_up",
    "ep_coefficient",
    "gr_roughness",
    "chisq",
)


class History:
    """Column-oriented record of EPSR runs.

    Each column is available as an attribute returning a numpy array, e.g.
    ``history.r_factor``. A history built from :func:`refine` has the columns of
    :class:`~epsr.Summary` plus ``run`` and ``r_factor_per_dataset`` (2-D); one
    read from an ``.erg`` file (:meth:`from_erg`) has :data:`epsr.io.ERG_COLUMNS`.
    """

    def __init__(self, columns: dict[str, list[Any]] | None = None) -> None:
        self._columns: dict[str, list[Any]] = columns if columns is not None else {}
        self.summaries: list[Summary] = []

    @classmethod
    def from_summaries(cls, summaries: Sequence[Summary]) -> History:
        history = cls()
        for summary in summaries:
            history.append(summary)
        return history

    @classmethod
    def from_erg(cls, path: str | PathLike[str]) -> History:
        """Read the ``.erg`` file, one row per EPSR iteration ever run in the folder."""
        erg_data = read_erg(path)
        columns: dict[str, list[Any]] = {
            name: list(erg_data[:, column_index]) for column_index, name in enumerate(ERG_COLUMNS)
        }
        columns["run"] = list(range(1, len(erg_data) + 1))
        return cls(columns)

    def append(self, summary: Summary) -> None:
        """Add one run's summary as a new row."""
        run_number = len(self) + 1
        self._columns.setdefault("run", []).append(run_number)
        for name in SUMMARY_COLUMNS:
            self._columns.setdefault(name, []).append(getattr(summary, name))
        self._columns.setdefault("r_factor_per_dataset", []).append(summary.r_factor_per_dataset)
        self.summaries.append(summary)

    @property
    def columns(self) -> list[str]:
        return list(self._columns)

    def column(self, name: str) -> npt.NDArray[Any]:
        try:
            return np.asarray(self._columns[name])
        except KeyError:
            raise KeyError(f"no column {name!r}; columns are {self.columns}") from None

    def __getattr__(self, name: str) -> npt.NDArray[Any]:
        # Called only for names that are not normal attributes: treat them as columns.
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            return self.column(name)
        except KeyError as error:
            raise AttributeError(str(error)) from None

    def __len__(self) -> int:
        return len(self._columns.get("run", []))

    def __getitem__(self, row: int) -> dict[str, Any]:
        """One row as a dict of column values."""
        return {name: values[row] for name, values in self._columns.items()}

    def __iter__(self) -> Iterator[dict[str, Any]]:
        return (self[row] for row in range(len(self)))

    def __repr__(self) -> str:
        return f"<epsr.History {len(self)} rows, columns {self.columns}>"

    def _flat_columns(self) -> dict[str, npt.NDArray[Any]]:
        """All columns as 1-D arrays; a 2-D column ``name`` becomes ``name_0``, ``name_1``..."""
        flat: dict[str, npt.NDArray[Any]] = {}
        for name in self._columns:
            values = self.column(name)
            if values.ndim == 2:
                for sub_index in range(values.shape[1]):
                    flat[f"{name}_{sub_index}"] = values[:, sub_index]
            else:
                flat[name] = values
        return flat

    def to_csv(self, path: str | PathLike[str]) -> Path:
        """Write the history as CSV (2-D columns are split into ``name_0, name_1...``)."""
        flat = self._flat_columns()
        path = Path(path)
        table = np.column_stack(list(flat.values())) if flat else np.empty((0, 0))
        np.savetxt(path, table, delimiter=",", header=",".join(flat), comments="")
        return path

    def to_dataframe(self) -> Any:
        """Return a ``pandas.DataFrame`` (requires the ``pandas`` extra)."""
        pandas = import_optional("pandas", extra="pandas")
        return pandas.DataFrame(self._flat_columns())


def _with_progress_bar(run_numbers: Iterable[int], enabled: bool | None, total: int) -> Any:
    """Wrap ``run_numbers`` in a tqdm progress bar when enabled and tqdm is installed."""
    if enabled is None:
        enabled = sys.stdout.isatty()
    if not enabled:
        return run_numbers
    try:
        from tqdm.auto import tqdm
    except ImportError:
        return run_numbers
    return tqdm(run_numbers, total=total, desc="EPSR", unit="run")


def refine(
    project: Project | str | PathLike[str],
    n_runs: int,
    *,
    overrides: Overrides = None,
    initial_overrides: Overrides = None,
    callback: Callable[[int, Simulation], bool | None] | None = None,
    until: Callable[[History], bool] | None = None,
    snapshot_every: int | None = None,
    verbose: bool = False,
    progress: bool | None = None,
) -> History:
    """Run EPSR ``n_runs`` times in a row on one folder, as the EPSR shell loop does.

    Each run reads the folder afresh into a new :class:`Simulation`, runs its
    ``niter`` iterations and saves the state back to the folder.

    Parameters
    ----------
    project:
        A :class:`Project` or a folder path. The folder is modified in place.
    overrides:
        Applied to every run (see :func:`epsr.render_overrides`).
    initial_overrides:
        Applied to the first run only, after ``overrides`` (e.g. ``{"ireset": 1}``,
        which must not be repeated).
    callback:
        ``callback(run_number, simulation)`` after each run (1-based); return
        ``True`` to stop.
    until:
        ``until(history)`` after each run; stop when it returns true.
    snapshot_every:
        Copy the folder to ``<folder>.snapshots/runNNNN`` every this many runs.
    progress:
        Show a ``tqdm`` progress bar (if installed). Default: when stdout is a terminal.
    """
    if not isinstance(project, Project):
        project = Project(project)
    history = History()

    for run_number in _with_progress_bar(range(1, n_runs + 1), progress, n_runs):
        first_run_overrides = initial_overrides if run_number == 1 else None
        override_lines = _combine_overrides(overrides, first_run_overrides)
        simulation = project.simulation(override_lines, verbose=verbose)
        simulation.run()
        history.append(simulation.summary)

        if snapshot_every and run_number % snapshot_every == 0:
            project.snapshot(f"run{run_number:04d}", overwrite=True)

        callback_asks_to_stop = callback is not None and bool(callback(run_number, simulation))
        target_reached = until is not None and until(history)
        if callback_asks_to_stop or target_reached:
            break
    return history


def _combine_overrides(overrides: Overrides, extra_overrides: Overrides) -> list[str]:
    """Combine two override sets into one list of input lines.

    Two mappings are merged (``extra_overrides`` wins; ``data`` blocks merged per
    dataset) so that all global keywords precede the dataset blocks. Otherwise
    the line lists are concatenated.
    """
    both_are_mappings = isinstance(overrides, Mapping) and isinstance(extra_overrides, Mapping)
    if not both_are_mappings:
        return render_overrides(overrides) + render_overrides(extra_overrides)

    assert isinstance(overrides, Mapping) and isinstance(extra_overrides, Mapping)
    merged: dict[str, Any] = {**overrides, **extra_overrides}
    if "data" in overrides and "data" in extra_overrides:
        merged_datasets = {
            dataset_number: dict(keywords) for dataset_number, keywords in overrides["data"].items()
        }
        for dataset_number, keywords in extra_overrides["data"].items():
            merged_datasets.setdefault(dataset_number, {}).update(keywords)
        merged["data"] = merged_datasets
    return render_overrides(merged)


def refine_many(
    projects: Sequence[Project | str | PathLike[str]],
    n_runs: int,
    *,
    workers: int | None = None,
    **refine_options: Any,
) -> list[History]:
    """Refine several *different* folders in parallel threads.

    The engine releases the GIL while it runs, so the folders really run
    concurrently. Keyword arguments are passed to :func:`refine` (``progress``
    defaults to off). Never pass the same folder twice.
    """
    folders = [
        Path(project.folder if isinstance(project, Project) else project).resolve()
        for project in projects
    ]
    if len(set(folders)) != len(folders):
        raise ValueError("refine_many: the same folder appears more than once")
    refine_options.setdefault("progress", False)
    with ThreadPoolExecutor(max_workers=workers or len(projects)) as pool:
        futures = [pool.submit(refine, project, n_runs, **refine_options) for project in projects]
        return [future.result() for future in futures]
