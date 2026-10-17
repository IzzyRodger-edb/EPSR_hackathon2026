"""Repeated EPSR runs: :func:`refine`, :func:`refine_many` and :class:`History`."""

from __future__ import annotations

import sys
from collections.abc import Callable, Iterator, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from os import PathLike
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

from .io import ERG_COLUMNS, _import_optional, read_erg
from .project import Project
from .simulation import Overrides, Simulation, Summary, render_overrides

__all__ = ["History", "refine", "refine_many"]

_SUMMARY_COLUMNS = (
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
    ``hist.r_factor``. A history built from :func:`refine` has the columns of
    :class:`~epsr.Summary` plus ``run`` and ``r_factor_per_dataset`` (2-D); one
    read from an ``.erg`` file (:meth:`from_erg`) has :data:`epsr.io.ERG_COLUMNS`.
    """

    def __init__(self, columns: dict[str, list[Any]] | None = None) -> None:
        self._cols: dict[str, list[Any]] = columns if columns is not None else {}
        self.summaries: list[Summary] = []

    @classmethod
    def from_summaries(cls, summaries: Sequence[Summary]) -> History:
        h = cls()
        for s in summaries:
            h.append(s)
        return h

    @classmethod
    def from_erg(cls, path: str | PathLike[str]) -> History:
        """Read the ``.erg`` file, one row per EPSR iteration ever run in the folder."""
        data = read_erg(path)
        cols: dict[str, list[Any]] = {name: list(data[:, i]) for i, name in enumerate(ERG_COLUMNS)}
        cols["run"] = list(range(1, len(data) + 1))
        return cls(cols)

    def append(self, summary: Summary) -> None:
        n = len(self) + 1
        self._cols.setdefault("run", []).append(n)
        for name in _SUMMARY_COLUMNS:
            self._cols.setdefault(name, []).append(getattr(summary, name))
        self._cols.setdefault("r_factor_per_dataset", []).append(summary.r_factor_per_dataset)
        self.summaries.append(summary)

    @property
    def columns(self) -> list[str]:
        return list(self._cols)

    def column(self, name: str) -> npt.NDArray[Any]:
        try:
            return np.asarray(self._cols[name])
        except KeyError:
            raise KeyError(f"no column {name!r}; columns are {self.columns}") from None

    def __getattr__(self, name: str) -> npt.NDArray[Any]:
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            return self.column(name)
        except KeyError as e:
            raise AttributeError(str(e)) from None

    def __len__(self) -> int:
        return len(self._cols.get("run", []))

    def __getitem__(self, i: int) -> dict[str, Any]:
        return {k: v[i] for k, v in self._cols.items()}

    def __iter__(self) -> Iterator[dict[str, Any]]:
        return (self[i] for i in range(len(self)))

    def __repr__(self) -> str:
        return f"<epsr.History {len(self)} rows, columns {self.columns}>"

    def _flat(self) -> dict[str, npt.NDArray[Any]]:
        out: dict[str, npt.NDArray[Any]] = {}
        for k in self._cols:
            a = self.column(k)
            if a.ndim == 2:
                for j in range(a.shape[1]):
                    out[f"{k}_{j}"] = a[:, j]
            else:
                out[k] = a
        return out

    def to_csv(self, path: str | PathLike[str]) -> Path:
        """Write the history as CSV (2-D columns are split into ``name_0, name_1...``)."""
        flat = self._flat()
        path = Path(path)
        np.savetxt(
            path,
            np.column_stack(list(flat.values())) if flat else np.empty((0, 0)),
            delimiter=",",
            header=",".join(flat),
            comments="",
        )
        return path

    def to_dataframe(self) -> Any:
        """Return a ``pandas.DataFrame`` (requires the ``pandas`` extra)."""
        pd = _import_optional("pandas", "pandas")
        return pd.DataFrame(self._flat())


def _progress(iterable: Sequence[int], enabled: bool | None, total: int) -> Any:
    if enabled is None:
        enabled = sys.stdout.isatty()
    if not enabled:
        return iterable
    try:
        from tqdm.auto import tqdm
    except ImportError:
        return iterable
    return tqdm(iterable, total=total, desc="EPSR", unit="run")


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
        ``callback(k, sim)`` after run ``k`` (1-based); return ``True`` to stop.
    until:
        ``until(history)`` after each run; stop when it returns true.
    snapshot_every:
        Copy the folder to ``<folder>.snapshots/runNNNN`` every this many runs.
    progress:
        Show a ``tqdm`` progress bar (if installed). Default: when stdout is a terminal.
    """
    proj = project if isinstance(project, Project) else Project(project)
    hist = History()
    for k in _progress(range(1, n_runs + 1), progress, n_runs):
        lines = _combine_overrides(overrides, initial_overrides if k == 1 else None)
        sim = proj.simulation(lines, verbose=verbose)
        sim.run()
        hist.append(sim.summary)
        if snapshot_every and k % snapshot_every == 0:
            proj.snapshot(f"run{k:04d}", overwrite=True)
        stop = bool(callback(k, sim)) if callback is not None else False
        if stop or (until is not None and until(hist)):
            break
    return hist


def _combine_overrides(overrides: Overrides, extra: Overrides) -> list[str]:
    """Combine two override sets into one list of input lines.

    Two mappings are merged (``extra`` wins; ``data`` blocks merged per dataset)
    so that all global keywords precede the dataset blocks. Otherwise the line
    lists are concatenated.
    """
    if isinstance(overrides, Mapping) and isinstance(extra, Mapping):
        merged: dict[str, Any] = {**overrides, **extra}
        if "data" in overrides and "data" in extra:
            data = {n: dict(v) for n, v in overrides["data"].items()}
            for n, v in extra["data"].items():
                data.setdefault(n, {}).update(v)
            merged["data"] = data
        return render_overrides(merged)
    return render_overrides(overrides) + render_overrides(extra)


def refine_many(
    projects: Sequence[Project | str | PathLike[str]],
    n_runs: int,
    *,
    workers: int | None = None,
    **kwargs: Any,
) -> list[History]:
    """Refine several *different* folders in parallel threads.

    The engine releases the GIL while it runs, so the folders really run
    concurrently. Keyword arguments are passed to :func:`refine` (``progress``
    defaults to off). Never pass the same folder twice.
    """
    folders = [Path(p.folder if isinstance(p, Project) else p).resolve() for p in projects]
    if len(set(folders)) != len(folders):
        raise ValueError("refine_many: the same folder appears more than once")
    kwargs.setdefault("progress", False)
    with ThreadPoolExecutor(max_workers=workers or len(projects)) as pool:
        futures = [pool.submit(refine, p, n_runs, **kwargs) for p in projects]
        return [f.result() for f in futures]
