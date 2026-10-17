"""Tidy (long-format) tables of a run, for pandas, seaborn and web front ends.

The figures Mellor makes come from ``epsrvis`` (matplotlib). An application that
wants its own plots should not parse EPSR files: :func:`tidy_tables` returns the
final fits, partial g(r), potentials and the convergence history as lists of flat
row dicts, and :func:`write_tables` writes them as CSV files (``<run>/tables/``)
that ``pandas.read_csv`` and ``seaborn`` consume directly::

    import pandas as pd, seaborn as sns
    fits = pd.read_csv(run_dir / "tables" / "fits.csv")
    sns.relplot(data=fits, x="q", y="simulated", hue="dataset", kind="line")

The recorded trajectory (``<run>/trajectory``) is the preferred source; without it
the EPSR output files in ``<run>/epsr`` are used (no errors, no per-iteration
R-factors, as EPSR only keeps those at run ends).
"""

from __future__ import annotations

import contextlib
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "TABLES",
    "tidy_tables",
    "write_tables",
    "fits_rows",
    "partials_rows",
    "potentials_rows",
    "history_rows",
    "runs_rows",
]

#: table name -> columns and meaning.
TABLES = {
    "fits": "final fit of each data set over its fitted Q range: dataset, q, data, error, simulated, residual",
    "partials": "final partial g(r) of each pair: pair, r, g_inter, g_intra, g_total",
    "potentials": "final potentials of each pair (kJ/mol): pair, r, empirical, reference",
    "history": "one row per EPSR iteration: step, run, iteration, ereq, energy, pressure, ep_amplitude, r_factor, r_factor[<dataset>]...",
    "runs": "one row per EPSR run (from manifest.json): run, stage, niter, ereq, r_factor, energy, pressure, ep_amplitude, duration_s",
}


def _float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def _trajectory(folder: Path) -> Any | None:
    path = folder / "trajectory"
    if not (path / "meta.json").is_file():
        return None
    from epsrvis import Trajectory

    return Trajectory.load(path)


def _outputs(folder: Path) -> Any | None:
    manifest = folder / "manifest.json"
    epsr_folder = folder / "epsr"
    if manifest.is_file():
        with contextlib.suppress(OSError, ValueError):
            epsr_folder = folder / json.loads(manifest.read_text()).get(
                "epsr_folder", "epsr"
            )
    if not epsr_folder.is_dir():
        return None
    from epsrvis.outputs import OutputFolder

    try:
        return OutputFolder.find(epsr_folder)
    except (FileNotFoundError, ValueError, StopIteration):
        return None


# --------------------------------------------------------------------------- tables
def fits_rows(folder: str | Path) -> list[dict[str, Any]]:
    """Measured and simulated F(Q) of every data set at the end of the run."""
    folder = Path(folder)
    rows: list[dict[str, Any]] = []
    traj = _trajectory(folder)
    if traj is not None:
        last = traj.nsteps - 1
        q = traj.q
        for d, name in enumerate(traj.datasets):
            sl = traj.q_slice(d)
            data = np.asarray(traj.fq_data[last, d], dtype=float)
            sim = np.asarray(traj.fq_sim[last, d], dtype=float)
            err = (
                np.asarray(traj.fq_err[last, d], dtype=float)
                if traj.fq_err is not None
                else None
            )
            for i in range(sl.start, sl.stop):
                rows.append(
                    {
                        "dataset": name,
                        "q": _float(q[i]),
                        "data": _float(data[i]),
                        "error": _float(err[i]) if err is not None else None,
                        "simulated": _float(sim[i]),
                        "residual": _float(data[i] - sim[i]),
                    }
                )
        return rows
    out = _outputs(folder)
    if out is None:
        return rows
    data_table, sim_table = out.table("t"), out.table("u")
    for d, name in enumerate(data_table.labels):
        mask = data_table.values[d] != 0
        x = data_table.x[mask]
        data = data_table.values[d][mask]
        sim = np.interp(x, sim_table.x, sim_table.values[d])
        for q_value, data_value, sim_value in zip(x, data, sim, strict=True):
            rows.append(
                {
                    "dataset": name,
                    "q": _float(q_value),
                    "data": _float(data_value),
                    "error": None,
                    "simulated": _float(sim_value),
                    "residual": _float(data_value - sim_value),
                }
            )
    return rows


def partials_rows(folder: str | Path) -> list[dict[str, Any]]:
    """Final partial g(r) (intermolecular, intramolecular and their sum) per pair."""
    folder = Path(folder)
    rows: list[dict[str, Any]] = []
    traj = _trajectory(folder)
    if traj is not None:
        last = traj.nsteps - 1
        inter = np.asarray(traj.gr_inter[last], dtype=float)
        intra = (
            np.asarray(traj.gr_intra[last], dtype=float)
            if traj.gr_intra is not None
            else np.zeros_like(inter)
        )
        for ip, pair in enumerate(traj.pairs):
            for r_value, g_inter, g_intra in zip(
                traj.r, inter[ip], intra[ip], strict=True
            ):
                rows.append(
                    {
                        "pair": pair,
                        "r": _float(r_value),
                        "g_inter": _float(g_inter),
                        "g_intra": _float(g_intra),
                        "g_total": _float(g_inter + g_intra),
                    }
                )
        return rows
    out = _outputs(folder)
    if out is None:
        return rows
    g, y = out.table("g"), out.table("y")
    for ip, pair in enumerate(g.labels):
        intra = y.values[ip] if ip < len(y.labels) else np.zeros_like(g.values[ip])
        for r_value, g_inter, g_intra in zip(g.x, g.values[ip], intra, strict=True):
            rows.append(
                {
                    "pair": pair,
                    "r": _float(r_value),
                    "g_inter": _float(g_inter),
                    "g_intra": _float(g_intra),
                    "g_total": _float(g_inter + g_intra),
                }
            )
    return rows


def potentials_rows(folder: str | Path) -> list[dict[str, Any]]:
    """Final empirical and reference potentials per pair (kJ/mol)."""
    folder = Path(folder)
    rows: list[dict[str, Any]] = []
    traj = _trajectory(folder)
    if traj is not None:
        last = traj.nsteps - 1
        ep = np.asarray(traj.ep[last], dtype=float)
        ref = (
            np.asarray(traj.ref_pot[last], dtype=float)
            if traj.ref_pot is not None
            else np.full_like(ep, np.nan)
        )
        for ip, pair in enumerate(traj.pairs):
            for r_value, e_value, ref_value in zip(
                traj.r, ep[ip], ref[ip], strict=True
            ):
                rows.append(
                    {
                        "pair": pair,
                        "r": _float(r_value),
                        "empirical": _float(e_value),
                        "reference": _float(ref_value),
                    }
                )
        return rows
    out = _outputs(folder)
    if out is None:
        return rows
    p = out.table("p")
    try:
        o = out.table("o")
    except (FileNotFoundError, OSError):
        o = None
    for ip, pair in enumerate(p.labels):
        ref = (
            np.interp(p.x, o.x, o.values[ip])
            if o is not None and ip < len(o.labels)
            else np.full_like(p.x, np.nan)
        )
        for r_value, e_value, ref_value in zip(p.x, p.values[ip], ref, strict=True):
            rows.append(
                {
                    "pair": pair,
                    "r": _float(r_value),
                    "empirical": _float(e_value),
                    "reference": _float(ref_value),
                }
            )
    return rows


def history_rows(folder: str | Path) -> list[dict[str, Any]]:
    """One row per EPSR iteration (from the trajectory; from ``.erg`` without one)."""
    folder = Path(folder)
    rows: list[dict[str, Any]] = []
    traj = _trajectory(folder)
    if traj is not None:
        from epsrvis import metrics

        fit = metrics.fit_metrics(traj)
        h = traj.history
        for k in range(traj.nsteps):
            row = {
                "step": int(h["step"][k]),
                "run": int(h["run"][k]),
                "iteration": int(h["iter"][k]),
                "ereq": _float(h.get("ereq", np.zeros(traj.nsteps))[k]),
                "energy": _float(h["energd"][k]),
                "pressure": _float(h["pressd"][k]),
                "ep_amplitude": _float(h["energabs"][k]),
                "r_factor": _float(fit["rfactor_tot"][k]),
            }
            for d, name in enumerate(traj.datasets):
                row[f"r_factor[{name}]"] = _float(fit["rfactor"][k, d])
                row[f"r2[{name}]"] = _float(fit["r2"][k, d])
            rows.append(row)
        return rows
    out = _outputs(folder)
    if out is None:
        return rows
    erg = out.erg()
    n = len(erg["energy"])
    for k in range(n):
        rows.append(
            {
                "step": k + 1,
                "run": None,
                "iteration": None,
                "ereq": None,
                "energy": _float(erg["energy"][k]),
                "pressure": _float(erg["pressure"][k]),
                "ep_amplitude": _float(erg["ep_amplitude"][k]),
                "r_factor": _float(erg["rfactor_raw"][k]),
            }
        )
    return rows


def runs_rows(folder: str | Path) -> list[dict[str, Any]]:
    """One row per EPSR run from ``manifest.json`` (empty without a manifest)."""
    manifest = Path(folder) / "manifest.json"
    if not manifest.is_file():
        return []
    try:
        runs = json.loads(manifest.read_text()).get("runs", [])
    except (OSError, ValueError):
        return []
    rows = []
    for item in runs:
        summary = item.get("summary", {}) or {}
        rows.append(
            {
                "run": item.get("run"),
                "stage": item.get("stage"),
                "niter": item.get("niter"),
                "ereq": _float(item.get("ereq")),
                "r_factor": _float(summary.get("r_factor")),
                "energy": _float(summary.get("energy")),
                "pressure": _float(summary.get("pressure")),
                "ep_amplitude": _float(summary.get("ep_amplitude")),
                "duration_s": _float(item.get("duration_s")),
            }
        )
    return rows


def tidy_tables(folder: str | Path) -> dict[str, list[dict[str, Any]]]:
    """Every table of :data:`TABLES` for a run folder (empty lists where no source exists)."""
    return {
        "fits": fits_rows(folder),
        "partials": partials_rows(folder),
        "potentials": potentials_rows(folder),
        "history": history_rows(folder),
        "runs": runs_rows(folder),
    }


def write_tables(
    folder: str | Path, out_dir: str | Path | None = None
) -> dict[str, Path]:
    """Write the tables as CSV files (default ``<folder>/tables``); returns name -> path."""
    folder = Path(folder)
    out = Path(out_dir) if out_dir is not None else folder / "tables"
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for name, rows in tidy_tables(folder).items():
        if not rows:
            continue
        columns: list[str] = []
        for row in rows:
            for key in row:
                if key not in columns:
                    columns.append(key)
        path = out / f"{name}.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: _cell(row.get(key)) for key in columns})
        written[name] = path
    return written


def _cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.7g}"
    return value
