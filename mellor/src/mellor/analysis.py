"""Numbers from a finished run: fit quality, convergence, structure.

Everything here is computed from the recorded trajectory (``epsrvis.Trajectory``,
``epsrvis.metrics``) and the EPSR output files; nothing comes from the model.
The result is one JSON-serialisable dict that the report, the narrative prompt and
the chat tools all use.

Without a recorded trajectory (``--no-record``) the per-iteration quantities (R²,
χ², trends over the last run) are left out and a flag says so; the per-run numbers
come from ``manifest.json`` and the first shells from EPSR's own ``.g01``/``.z01``
files when they exist.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .runner import RunResult

__all__ = [
    "analyse",
    "analysis_markdown",
    "load_analysis",
    "save_analysis",
    "plain",
    "EP_LIMIT_TOLERANCE",
    "HIGH_PRESSURE_KBAR",
    "POOR_FIT_R_FACTOR",
    "ENERGY_DRIFT_PERCENT",
    "R_FACTOR_RISE_PERCENT",
]

#: The EP amplitude counts as "at the ereq limit" within this fraction of ereq.
EP_LIMIT_TOLERANCE = 0.02
#: |pressure| above this (kbar) is flagged: a liquid at its measured density should
#: sit at a few kbar in EPSR (water: 5.4; the unconverged silica case: 109).
HIGH_PRESSURE_KBAR = 10.0
#: A final R-factor above this is flagged as a poor fit (liquids on good neutron
#: data reach a few 1e-4; the Cu-BTC crystal with intermolecular bonds only: 1.6e-2).
POOR_FIT_R_FACTOR = 1e-2
#: Energy change between the first and last quarter of the last run (per cent of the
#: first-quarter mean) above which the run counts as still drifting.
ENERGY_DRIFT_PERCENT = 2.0
#: R-factor rise over the last run (per cent) above which it is flagged.
R_FACTOR_RISE_PERCENT = 5.0


# --------------------------------------------------------------------------- plain values
def plain(value: Any) -> Any:
    """``value`` with numpy scalars and arrays turned into Python numbers and lists,
    NaN and infinities into None, paths into strings (recursively)."""
    if isinstance(value, dict):
        return {str(key): plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    if isinstance(value, np.ndarray):
        return [plain(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        value = value.item()
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def _number(value: Any) -> float | None:
    """``value`` as a finite float, or None."""
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _ratio(numerator: Any, denominator: Any) -> float | None:
    top, bottom = _number(numerator), _number(denominator)
    if top is None or bottom is None or bottom == 0.0:
        return None
    return top / bottom


def _percent_change(start: Any, end: Any) -> float | None:
    first, last = _number(start), _number(end)
    if first is None or last is None or first == 0.0:
        return None
    return (last - first) / abs(first) * 100.0


# --------------------------------------------------------------------------- the analysis
def analyse(result: RunResult) -> dict[str, Any]:
    """The analysis dict of a run.

    Keys (all plain Python types):

    * ``system``: name, title, kind, n_atoms, n_molecules, box_side, density,
      mass_density_g_cm3 (when masses are known), temperature, atom types with their
      counts, pairs, data set labels and files;
    * ``schedule``: the runs (stage, niter, ereq) and total iterations, quick flag;
    * ``runs``: per EPSR run, at its last iteration: r_factor, r_factor_per_dataset
      (dict label -> value), energy, energy_intra, pressure, ep_amplitude, chisq,
      duration_s; plus from the trajectory ``r2`` and ``chi2`` per data set
      (``epsrvis.metrics.fit_metrics``);
    * ``fit``: per data set: initial and final R-factor, improvement factor, final R²,
      fitted Q range;
    * ``convergence``: trend of energy/pressure/EP amplitude over the last run (mean of
      the first and last quarter), ``ep_at_limit`` (EP amplitude within 2 % of ereq),
      R-factor change over the last run in per cent, and a short verdict text;
    * ``shells``: per pair: first-peak position and height, first minimum, coordination
      number (``epsrvis.metrics.first_shell`` on the final intermolecular g(r) with
      ``running_coordination``), with the pair's partner density;
    * ``flags``: list of plain-English warnings (high pressure, EP at limit, R-factor
      rising, unrecorded trajectory, quick schedule, library fallback...);
    * ``figures``: empty here (filled by ``figures.make_figures``).

    Without a trajectory (``--no-record``) the trajectory-based entries are omitted
    and a flag says so; the manifest values are still reported.
    """
    traj, trajectory_problem = _load_trajectory(result)
    datasets = _datasets(result)
    labels = [label for label, _file in datasets]
    traj_index = _trajectory_indices(traj, datasets)
    runs = _runs_section(result, traj, labels, traj_index)
    analysis: dict[str, Any] = {
        "system": _system_section(result, traj),
        "schedule": _schedule_section(result, traj),
        "runs": runs,
        "fit": _fit_section(result, traj, runs, datasets, traj_index),
        "convergence": _convergence_section(result, traj, runs),
        "shells": _shells_section(result, traj),
        "flags": [],
        "figures": {},
    }
    analysis["flags"] = _flags(result, analysis, trajectory_problem)
    return plain(analysis)


# --------------------------------------------------------------------------- inputs
def _load_trajectory(result: RunResult) -> tuple[Any | None, str | None]:
    """``(Trajectory, None)`` or ``(None, why there is none)``."""
    path = result.trajectory_path
    if path is None:
        return None, "the trajectory was not recorded (--no-record)"
    if not (path / "meta.json").is_file():
        return None, f"no recorded trajectory at {path}"
    from epsrvis import Trajectory

    try:
        traj = Trajectory.load(path)
    except Exception as error:  # a corrupt recording must not hide the manifest numbers
        return (
            None,
            f"the trajectory at {path} could not be read ({type(error).__name__}: {error})",
        )
    if traj.nsteps == 0:
        return None, f"the trajectory at {path} holds no iterations"
    return traj, None


def _datasets(result: RunResult) -> list[tuple[str, str]]:
    """``(label, file name)`` of every data set in ``.EPSR.inp`` order."""
    items = result.manifest.get("datasets") or []
    if items:
        return [
            (str(item.get("label", f"dataset{i + 1}")), str(item.get("file", "")))
            for i, item in enumerate(items)
        ]
    return [(dataset.label, dataset.file_name) for dataset in result.spec.datasets]


def _trajectory_indices(traj: Any, datasets: list[tuple[str, str]]) -> dict[str, int]:
    """label -> index of that data set along the trajectory's data-set axis.

    The trajectory names data sets by file stem; they are matched by stem, with the
    position in the ``.EPSR.inp`` as the fallback.
    """
    if traj is None:
        return {}
    names = list(traj.datasets)
    by_stem = {name: index for index, name in enumerate(names)}
    indices: dict[str, int] = {}
    for position, (label, file) in enumerate(datasets):
        stem = Path(file).stem if file else ""
        if stem in by_stem:
            indices[label] = by_stem[stem]
        elif position < len(names):
            indices[label] = position
    return indices


def _run_end_steps(traj: Any) -> dict[int, int]:
    """run number -> index of the last recorded step of that run."""
    starts = traj.run_starts
    ends: dict[int, int] = {}
    for position, (_first, run, _ereq) in enumerate(starts):
        last = (
            starts[position + 1][0] - 1
            if position + 1 < len(starts)
            else traj.nsteps - 1
        )
        ends[int(run)] = int(last)
    return ends


def _fit_metrics(traj: Any) -> dict[str, np.ndarray]:
    from epsrvis import metrics

    return metrics.fit_metrics(traj)


# --------------------------------------------------------------------------- sections
def _system_section(result: RunResult, traj: Any) -> dict[str, Any]:
    spec, manifest = result.spec, result.manifest
    n_atoms = int(manifest.get("n_atoms") or spec.n_atoms)
    if traj is not None:
        fractions = traj.concentrations()
        counts = {
            label: int(round(fraction * n_atoms))
            for label, fraction in fractions.items()
        }
    else:
        counts = spec.type_counts()
        total = sum(counts.values())
        fractions = (
            {label: count / total for label, count in counts.items()} if total else {}
        )
    density = _number(manifest.get("density")) or spec.density
    section: dict[str, Any] = {
        "name": spec.name,
        "title": spec.title,
        "kind": spec.kind,
        "source": spec.source,
        "n_atoms": n_atoms,
        "n_molecules": int(manifest.get("n_molecules") or spec.n_molecules),
        "molecules": [
            {"name": molecule.name, "count": molecule.count, "atoms": molecule.n_atoms}
            for molecule in spec.molecules
        ],
        "box_side": _number(manifest.get("box_side"))
        or (spec.box_side if spec.density > 0 and spec.n_atoms else None),
        "density": density,
        "temperature": spec.temperature,
        "atom_types": [
            {
                "label": atom_type.label,
                "element": atom_type.element,
                "count": counts.get(atom_type.label),
                "fraction": fractions.get(atom_type.label),
                "epsilon": atom_type.epsilon,
                "sigma": atom_type.sigma,
                "charge": atom_type.charge,
                "mass": atom_type.mass,
            }
            for atom_type in spec.atom_types
        ],
        "pairs": list(manifest.get("pairs") or _pairs(spec.labels)),
        "datasets": [
            {
                "label": dataset.label,
                "file": dataset.file_name,
                "isotopes": {
                    label: [
                        [isotope.mass_number, isotope.fraction] for isotope in isotopes
                    ]
                    for label, isotopes in dataset.isotopes.items()
                },
                "exchangeable": list(dataset.exchangeable),
            }
            for dataset in spec.datasets
        ],
        "epsrvis_case": manifest.get("epsrvis_case"),
        "versions": manifest.get("versions"),
        "created": manifest.get("created"),
    }
    mass_density = _mass_density(spec, fractions, density)
    if mass_density is not None:
        section["mass_density_g_cm3"] = mass_density
    return section


def _pairs(labels: list[str]) -> list[str]:
    """Pair labels in EPSR order (``A-A, A-B, ..., B-B, ...``)."""
    return [
        f"{first}-{second}" for i, first in enumerate(labels) for second in labels[i:]
    ]


def _mass_density(
    spec: Any, fractions: dict[str, float], density: float
) -> float | None:
    """g/cm³ from the atomic number density and the atom-type masses (None if unknown)."""
    if not fractions or not density:
        return None
    try:
        mean_mass = sum(
            spec.atom_type(label).mass * fraction
            for label, fraction in fractions.items()
        )
    except KeyError:
        return None
    if mean_mass <= 0:
        return None
    return density * mean_mass * 1.66053907


def _schedule_section(result: RunResult, traj: Any) -> dict[str, Any]:
    manifest, schedule = result.manifest, result.spec.schedule
    runs = manifest.get("schedule") or [
        {"stage": stage, "niter": niter, "ereq": ereq}
        for stage, niter, ereq in schedule.iterations(result.quick)
    ]
    return {
        "runs": [dict(run) for run in runs],
        "total_iterations": sum(int(run.get("niter") or 0) for run in runs),
        "quick": bool(result.quick),
        "ereq": schedule.ereq,
        "ntimes": schedule.ntimes,
        "complete": result.complete,
        "recorded": traj is not None,
        "recorded_iterations": int(traj.nsteps) if traj is not None else None,
        "duration_s": _number(manifest.get("duration_s")) or result.duration_s,
        "setup_steps": list(manifest.get("setup_steps") or []),
    }


def _runs_section(
    result: RunResult, traj: Any, labels: list[str], traj_index: dict[str, int]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    manifest_runs = result.runs
    if not manifest_runs and traj is not None:
        manifest_runs = _runs_from_trajectory(traj)
    fit = _fit_metrics(traj) if traj is not None else None
    ends = _run_end_steps(traj) if traj is not None else {}
    for position, run in enumerate(manifest_runs):
        summary = run.get("summary") or {}
        per_dataset = list(summary.get("r_factor_per_dataset") or [])
        row: dict[str, Any] = {
            "run": int(run.get("run") or position + 1),
            "stage": run.get("stage"),
            "niter": run.get("niter"),
            "ereq": _number(run.get("ereq")),
            "r_factor": _number(summary.get("r_factor")),
            "r_factor_per_dataset": {
                label: _number(per_dataset[i]) if i < len(per_dataset) else None
                for i, label in enumerate(labels)
            },
            "energy": _number(summary.get("energy")),
            "energy_inter": _number(summary.get("energy_inter")),
            "energy_intra": _number(summary.get("energy_intra")),
            "pressure": _number(summary.get("pressure")),
            "ep_amplitude": _number(summary.get("ep_amplitude")),
            "chisq": _number(summary.get("chisq")),
            "duration_s": _number(run.get("duration_s")),
            "overrides": list(run.get("overrides") or []),
        }
        step = ends.get(row["run"])
        if traj is not None and fit is not None and step is not None:
            row["last_iteration"] = int(traj.steps[step])
            row["r_factor_trajectory"] = float(fit["rfactor_tot"][step])
            row["r2"] = {
                label: float(fit["r2"][step, d]) for label, d in traj_index.items()
            }
            row["chi2"] = {
                label: float(fit["chi2"][step, d]) for label, d in traj_index.items()
            }
        rows.append(row)
    return rows


def _runs_from_trajectory(traj: Any) -> list[dict[str, Any]]:
    """Run summaries in the manifest's layout, from a trajectory without a manifest."""
    from epsrvis.systems import get_case

    stages = get_case(traj.case, traj.atom_types).stages
    history = traj.history
    rows = []
    for run, end in _run_end_steps(traj).items():
        start = min(k for k, r, _ in traj.run_starts if int(r) == run)
        rows.append(
            {
                "run": run,
                "stage": stages.get(run, f"run {run}"),
                "niter": int(end - start + 1),
                "ereq": float(history["ereq"][end]) if "ereq" in history else 0.0,
                "summary": {
                    "r_factor": float(history["rfactot"][end]),
                    "r_factor_per_dataset": [
                        float(v) for v in np.atleast_1d(history["rfacd"][end])
                    ],
                    "energy": float(history["energd"][end]),
                    "energy_inter": float(history["foldsv"][end]),
                    "energy_intra": float(history["foldsb"][end]),
                    "pressure": float(history["pressd"][end]),
                    "ep_amplitude": float(history["energabs"][end]),
                    "chisq": float(history["chisq"][end]),
                },
            }
        )
    return rows


def _fit_section(
    result: RunResult,
    traj: Any,
    runs: list[dict[str, Any]],
    datasets: list[tuple[str, str]],
    traj_index: dict[str, int],
) -> dict[str, Any]:
    first = runs[0] if runs else {}
    last = runs[-1] if runs else {}
    reference = next(
        (run for run in reversed(runs) if not run.get("ereq")), None
    )  # the last run on the reference potential alone
    fit = _fit_metrics(traj) if traj is not None else None
    section: dict[str, Any] = {}
    for label, file in datasets:
        initial = (first.get("r_factor_per_dataset") or {}).get(label)
        final = (last.get("r_factor_per_dataset") or {}).get(label)
        after_reference = (
            (reference.get("r_factor_per_dataset") or {}).get(label)
            if reference
            else None
        )
        entry: dict[str, Any] = {
            "file": file,
            "initial_r_factor": initial,
            "reference_r_factor": after_reference,
            "final_r_factor": final,
            "improvement_factor": _ratio(initial, final),
            "ep_improvement_factor": (
                _ratio(after_reference, final) if reference is not last else None
            ),
        }
        d = traj_index.get(label)
        if traj is not None and fit is not None and d is not None:
            q = traj.q[traj.q_slice(d)]
            entry.update(
                {
                    "r2": float(fit["r2"][-1, d]),
                    "chi2": float(fit["chi2"][-1, d]),
                    "r_factor_trajectory": float(fit["rfactor"][-1, d]),
                    "q_min": float(q[0]) if q.size else None,
                    "q_max": float(q[-1]) if q.size else None,
                    "n_points": int(q.size),
                    "trajectory_dataset": traj.datasets[d],
                }
            )
        section[label] = entry
    return section


def _convergence_section(
    result: RunResult, traj: Any, runs: list[dict[str, Any]]
) -> dict[str, Any]:
    last = runs[-1] if runs else {}
    previous = runs[-2] if len(runs) > 1 else None
    ereq = _number(last.get("ereq")) or 0.0
    ep = _number(last.get("ep_amplitude"))
    section: dict[str, Any] = {
        "last_run": last.get("run"),
        "last_run_stage": last.get("stage"),
        "ereq": ereq,
        "ep_amplitude": ep,
        "ep_fraction_of_ereq": _ratio(ep, ereq) if ereq else None,
        "ep_at_limit": bool(
            ereq and ep is not None and ep >= (1.0 - EP_LIMIT_TOLERANCE) * ereq
        ),
        "energy": _number(last.get("energy")),
        "pressure": _number(last.get("pressure")),
        "r_factor": _number(last.get("r_factor")),
        "r_factor_initial": _number(runs[0].get("r_factor")) if runs else None,
        "r_factor_before_last_run": (
            _number(previous.get("r_factor")) if previous else None
        ),
        "r_factor_change_percent": (
            _percent_change(previous.get("r_factor"), last.get("r_factor"))
            if previous
            else None
        ),
        "trends": {},
        "r_factor_rising": False,
        "energy_drift_percent": None,
    }
    if traj is not None and traj.nsteps > 0:
        section.update(_trajectory_convergence(traj))
    change = section["r_factor_change_percent"]
    if section["trends"].get("r_factor"):
        change = section["trends"]["r_factor"].get("change_percent")
    section["r_factor_rising"] = bool(
        change is not None and change > R_FACTOR_RISE_PERCENT
    )
    section["verdict"] = _verdict(section)
    return section


def _trajectory_convergence(traj: Any) -> dict[str, Any]:
    """Quarter means over the last run and when the EP reached its limit."""
    history = traj.history
    fit = _fit_metrics(traj)
    first_step = traj.run_starts[-1][0]
    n = traj.nsteps - first_step
    quarter = max(1, n // 4)
    series = {
        "energy": np.asarray(history["energd"], dtype=float),
        "pressure": np.asarray(history["pressd"], dtype=float),
        "ep_amplitude": np.asarray(history["energabs"], dtype=float),
        "r_factor": np.asarray(fit["rfactor_tot"], dtype=float),
    }
    trends: dict[str, Any] = {}
    for name, values in series.items():
        segment = values[first_step:]
        start = float(np.nanmean(segment[:quarter]))
        end = float(np.nanmean(segment[-quarter:]))
        trends[name] = {
            "start": start,
            "end": end,
            "change": end - start,
            "change_percent": _percent_change(start, end),
        }
    out: dict[str, Any] = {
        "last_run_iterations": int(n),
        "last_run_first_iteration": int(traj.steps[first_step]),
        "quarter_iterations": int(quarter),
        "trends": trends,
        "energy_drift_percent": trends["energy"]["change_percent"],
    }
    ereq_history = np.asarray(history.get("ereq", np.zeros(traj.nsteps)), dtype=float)
    ep_on = np.where(ereq_history > 0)[0]
    if ep_on.size:
        first_on = int(ep_on[0])
        out["ep_on_from_iteration"] = int(traj.steps[first_on])
        amplitude = series["ep_amplitude"]
        reached = [
            k
            for k in range(first_on, traj.nsteps)
            if amplitude[k] >= (1.0 - EP_LIMIT_TOLERANCE) * ereq_history[k]
        ]
        if reached:
            out["ep_limit_reached_at_iteration"] = int(traj.steps[reached[0]])
            out["ep_limit_reached_after"] = int(reached[0] - first_on + 1)
        else:
            out["ep_limit_reached_at_iteration"] = None
            out["ep_limit_reached_after"] = None
    return out


def _verdict(section: dict[str, Any]) -> str:
    """A factual sentence about the EP amplitude, the R-factor, the energy and the pressure."""
    parts: list[str] = []
    ereq, ep = section.get("ereq") or 0.0, section.get("ep_amplitude")
    if not ereq:
        parts.append("the empirical potential was off (ereq 0)")
    elif ep is None:
        parts.append(f"ereq was {ereq:g} kJ/mol but the EP amplitude was not reported")
    elif section.get("ep_at_limit"):
        text = f"the EP amplitude reached the ereq limit of {ereq:g} kJ/mol"
        after = section.get("ep_limit_reached_after")
        if after is not None:
            text += (
                f" after {after} iteration{'s' if after != 1 else ''} with the EP on"
            )
        parts.append(text)
    else:
        text = f"the EP amplitude is {ep:.2f} kJ/mol, {ep / ereq:.0%} of the requested {ereq:g} kJ/mol"
        trend = (section.get("trends") or {}).get("ep_amplitude")
        if trend and trend.get("change") is not None:
            text += (
                " and still rising"
                if trend["change"] > 0.02 * ereq
                else " and no longer rising"
            )
        parts.append(text)

    r_final = section.get("r_factor")
    r_before = section.get("r_factor_before_last_run")
    change = section.get(
        "r_factor_change_percent"
    )  # end of the previous run -> end of this one
    trend = (section.get("trends") or {}).get("r_factor", {}).get("change_percent")
    if r_final is not None and r_before is not None:
        verb = (
            "fell" if r_final < r_before else "rose" if r_final > r_before else "stayed"
        )
        text = f"the R-factor {verb} from {r_before:.1e} to {r_final:.1e} over the last run"
        details = [f"{change:+.0f} % between run ends"] if change is not None else []
        if trend is not None:
            details.append(f"{trend:+.1f} % between its first- and last-quarter means")
        if details:
            text += f" ({', '.join(details)})"
        initial = section.get("r_factor_initial")
        if initial is not None and initial != r_before:
            text += f", from {initial:.1e} at the start"
        parts.append(text)
    elif r_final is not None:
        parts.append(f"the final R-factor is {r_final:.1e}")

    energy = section.get("energy")
    drift = section.get("energy_drift_percent")
    if energy is not None:
        text = f"the energy is {energy:.1f} kJ/mol per molecule"
        if drift is not None:
            text += (
                f", drifting by {drift:+.1f} % over the last run"
                if abs(drift) > ENERGY_DRIFT_PERCENT
                else f" and flat over the last run ({drift:+.1f} %)"
            )
        parts.append(text)
    pressure = section.get("pressure")
    if pressure is not None:
        text = f"the pressure is {pressure:.1f} kbar"
        if abs(pressure) > HIGH_PRESSURE_KBAR:
            text += ", high for a condensed phase at its measured density"
        parts.append(text)
    return "; ".join(parts) + "." if parts else "no run summaries were found."


def _shells_section(result: RunResult, traj: Any) -> dict[str, Any]:
    if traj is not None:
        return _shells_from_trajectory(traj)
    return _shells_from_outputs(result)


def _shell_entry(
    centre: str,
    neighbour: str,
    shell: dict[str, float],
    n_reverse: float | None,
    partner_density: float | None,
) -> dict[str, Any]:
    r_peak, g_peak = _number(shell.get("r_peak")), _number(shell.get("g_peak"))
    r_min, n_coord = _number(shell.get("r_min")), _number(shell.get("n_coord"))
    if centre == neighbour:
        n_reverse = None
    entry: dict[str, Any] = {
        "centre": centre,
        "neighbour": neighbour,
        "r_peak": r_peak,
        "g_peak": g_peak,
        "r_min": r_min,
        "n_coord": n_coord,
        "n_coord_reverse": n_reverse,
        "partner_density": partner_density,
    }
    if r_peak is None:
        entry["description"] = (
            f"no intermolecular first peak found for {centre}-{neighbour}"
        )
    else:
        text = (
            f"first peak at {r_peak:.2f} Å (g = {g_peak:.2f}), minimum at {r_min:.2f} Å"
        )
        if n_coord is not None:
            text += (
                f"; {n_coord:.2f} {neighbour} around each {centre} within {r_min:.2f} Å"
            )
        if n_reverse is not None:
            text += f"; {n_reverse:.2f} {centre} around each {neighbour}"
        entry["description"] = text
    return entry


def _shells_from_trajectory(traj: Any) -> dict[str, Any]:
    from epsrvis import metrics

    last = traj.nsteps - 1
    g = np.asarray(traj.gr_inter[last], dtype=float)
    coordination = metrics.running_coordination(traj, last, total=False)
    rho = float(traj.history["rho"][last])
    fractions = traj.concentrations()
    shells: dict[str, Any] = {}
    for ip, (centre, neighbour) in enumerate(traj.meta["pairs"]):
        shell = metrics.first_shell(traj.r, g[ip], coordination[ip])
        n_reverse = None
        n_coord = _number(shell.get("n_coord"))
        if centre != neighbour and n_coord is not None and fractions.get(neighbour):
            n_reverse = n_coord * fractions[centre] / fractions[neighbour]
        shells[f"{centre}-{neighbour}"] = _shell_entry(
            centre, neighbour, shell, n_reverse, rho * fractions.get(neighbour, 0.0)
        )
    return shells


def _shells_from_outputs(result: RunResult) -> dict[str, Any]:
    """First shells from EPSR's own files: ``.g01`` (intermolecular g(r)) and ``.z01``
    (running coordination numbers, two columns per pair: N_ab then N_ba)."""
    from epsrvis import metrics
    from epsrvis.outputs import OutputFolder

    try:
        outputs = OutputFolder.find(result.epsr_folder)
        g_table = outputs.table("g")
        z_table = outputs.table("z")
    except (FileNotFoundError, OSError, ValueError):
        return {}
    density = _number(result.manifest.get("density")) or result.spec.density
    fractions = result.spec.type_fractions() if result.spec.molecules else {}
    shells: dict[str, Any] = {}
    for ip, label in enumerate(g_table.labels):
        centre, _, neighbour = label.partition("-")
        forward = reverse = None
        if 2 * ip + 1 < len(z_table.values):
            forward = np.interp(g_table.x, z_table.x, z_table.values[2 * ip])
            reverse = np.interp(g_table.x, z_table.x, z_table.values[2 * ip + 1])
        coordination = forward if forward is not None else np.zeros_like(g_table.x)
        shell = metrics.first_shell(g_table.x, g_table.values[ip], coordination)
        if forward is None:
            shell["n_coord"] = math.nan
        n_reverse = None
        r_min = _number(shell.get("r_min"))
        if reverse is not None and r_min is not None and centre != neighbour:
            n_reverse = float(reverse[int(np.searchsorted(g_table.x, r_min))])
        partner = density * fractions[neighbour] if neighbour in fractions else None
        shells[label] = _shell_entry(centre, neighbour, shell, n_reverse, partner)
    return shells


def _flags(
    result: RunResult, analysis: dict[str, Any], trajectory_problem: str | None
) -> list[str]:
    flags: list[str] = []
    schedule, convergence = analysis["schedule"], analysis["convergence"]
    if result.quick:
        flags.append(
            f"quick schedule ({schedule['total_iterations']} iterations in total): a preview, not a converged refinement"
        )
    if not result.complete:
        done, planned = len(analysis["runs"]), len(schedule["runs"])
        flags.append(
            f"the run did not complete: {done} of {planned} EPSR runs finished"
        )
    if trajectory_problem:
        flags.append(
            f"{trajectory_problem}: R², χ², the trends over the last run and the per-iteration "
            "figures are not available; the numbers come from manifest.json and the EPSR output files"
        )
    source = result.spec.source or ""
    if source == "model":
        flags.append(
            "the model system was proposed by the language model, not taken from the validated "
            "library: check the force field, the density and the isotopic compositions"
        )
    elif source.startswith("fallback"):
        flags.append(
            f"library fallback ({source}): the language model produced no valid specification, "
            "so the built-in system was used"
        )
    if convergence.get("ep_at_limit"):
        flags.append(
            f"the EP amplitude is at the ereq limit ({convergence['ep_amplitude']:.2f} of "
            f"{convergence['ereq']:g} kJ/mol): the fit may improve with a larger ereq, but a large "
            "empirical potential means the reference potential is poor"
        )
    elif convergence.get("ereq") and convergence.get("ep_amplitude") is not None:
        trend = (convergence.get("trends") or {}).get("ep_amplitude") or {}
        if (
            trend.get("change") is not None
            and trend["change"] > 0.02 * convergence["ereq"]
        ):
            flags.append(
                f"the EP amplitude ({convergence['ep_amplitude']:.2f} kJ/mol) was still rising towards "
                f"ereq {convergence['ereq']:g} kJ/mol at the end: the refinement is not converged"
            )
    pressure = convergence.get("pressure")
    if pressure is not None and abs(pressure) > HIGH_PRESSURE_KBAR:
        flags.append(
            f"high pressure ({pressure:.1f} kbar): the reference potential or the density may be off"
        )
    if convergence.get("r_factor_rising"):
        trend = (
            (convergence.get("trends") or {}).get("r_factor", {}).get("change_percent")
        )
        if trend is not None:
            flags.append(
                f"the R-factor rose by {trend:.0f} % between the first and last quarter of the last run"
            )
        else:
            flags.append(
                f"the R-factor rose by {convergence['r_factor_change_percent']:.0f} % over the last run"
            )
    drift = convergence.get("energy_drift_percent")
    if drift is not None and abs(drift) > ENERGY_DRIFT_PERCENT:
        flags.append(
            f"the energy changed by {drift:+.1f} % between the first and last quarter of the last "
            "run: not yet equilibrated"
        )
    for label, entry in analysis["fit"].items():
        final = entry.get("final_r_factor")
        if final is not None and final > POOR_FIT_R_FACTOR:
            flags.append(
                f"poor fit of {label} (R-factor {final:.1e}, above {POOR_FIT_R_FACTOR:.0e}): "
                "the model may be missing something the data see"
            )
    return flags


# --------------------------------------------------------------------------- markdown
def _sci(value: Any) -> str:
    number = _number(value)
    return "–" if number is None else f"{number:.2e}"


def _fixed(value: Any, decimals: int = 2) -> str:
    number = _number(value)
    return "–" if number is None else f"{number:.{decimals}f}"


def _general(value: Any, digits: int = 4) -> str:
    number = _number(value)
    return "–" if number is None else f"{number:.{digits}g}"


def _table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def analysis_markdown(analysis: dict[str, Any]) -> str:
    """A factual summary in Markdown (used when the model writes no narrative, and as
    the model's input): system, schedule, fit-quality table, convergence verdict,
    first-shell table, flags."""
    system = analysis.get("system") or {}
    schedule = analysis.get("schedule") or {}
    runs = analysis.get("runs") or []
    fit = analysis.get("fit") or {}
    convergence = analysis.get("convergence") or {}
    shells = analysis.get("shells") or {}
    flags = analysis.get("flags") or []
    out: list[str] = []

    title = system.get("title") or system.get("name") or "EPSR refinement"
    out.append(f"# {title}")
    out.append("")
    out.append("## System")
    out.append("")
    molecules = ", ".join(
        f"{m['count']} × {m['name']} ({m['atoms']} atoms)"
        for m in system.get("molecules") or []
    )
    composition = molecules or f"{system.get('n_atoms')} atoms"
    box = _fixed(system.get("box_side"), 3)
    line = (
        f"- {system.get('kind', '?')} system `{system.get('name', '?')}`: {composition}; "
        f"{system.get('n_atoms')} atoms in a cubic box of {box} Å; "
        f"density {_general(system.get('density'), 4)} atoms/Å³"
    )
    if system.get("mass_density_g_cm3") is not None:
        line += f" ({_fixed(system['mass_density_g_cm3'], 3)} g/cm³)"
    line += f"; temperature {_general(system.get('temperature'), 4)} K"
    out.append(line)
    types = "; ".join(
        f"{t['label']} ({t['element']}"
        + (f", {t['count']} atoms" if t.get("count") is not None else "")
        + f", ε {_general(t.get('epsilon'), 4)} kJ/mol, σ {_general(t.get('sigma'), 4)} Å, q {_general(t.get('charge'), 3)} e)"
        for t in system.get("atom_types") or []
    )
    out.append(f"- atom types: {types}")
    datasets = "; ".join(
        f"{d['label']}: `{d['file']}`"
        + (
            " ("
            + ", ".join(
                f"{label} " + "/".join(f"{int(m)}:{f:g}" for m, f in pairs)
                for label, pairs in (d.get("isotopes") or {}).items()
            )
            + ")"
            if d.get("isotopes")
            else ""
        )
        for d in system.get("datasets") or []
    )
    out.append(f"- data sets: {datasets}")
    if system.get("source"):
        out.append(f"- specification source: {system['source']}")
    out.append("")

    out.append("## Schedule")
    out.append("")
    out.append(
        _table(
            ["run", "stage", "iterations", "ereq (kJ/mol)"],
            [
                [
                    str(i + 1),
                    str(run.get("stage")),
                    str(run.get("niter")),
                    _general(run.get("ereq"), 3),
                ]
                for i, run in enumerate(schedule.get("runs") or [])
            ],
        )
    )
    out.append("")
    notes = [f"{schedule.get('total_iterations')} iterations in total"]
    if schedule.get("quick"):
        notes.append("quick schedule")
    notes.append("complete" if schedule.get("complete") else "incomplete")
    notes.append(
        f"{schedule.get('recorded_iterations')} iterations recorded"
        if schedule.get("recorded")
        else "trajectory not recorded"
    )
    if schedule.get("duration_s"):
        notes.append(f"{schedule['duration_s']:.0f} s")
    out.append("- " + "; ".join(notes))
    out.append("")

    out.append("## Runs")
    out.append("")
    out.append(
        _table(
            [
                "run",
                "stage",
                "iterations",
                "ereq",
                "R-factor",
                "energy (kJ/mol)",
                "pressure (kbar)",
                "EP amplitude (kJ/mol)",
                "time (s)",
            ],
            [
                [
                    str(run.get("run")),
                    str(run.get("stage")),
                    str(run.get("niter")),
                    _general(run.get("ereq"), 3),
                    _sci(run.get("r_factor")),
                    _general(run.get("energy"), 4),
                    _general(run.get("pressure"), 4),
                    _fixed(run.get("ep_amplitude"), 3),
                    _fixed(run.get("duration_s"), 1),
                ]
                for run in runs
            ],
        )
    )
    out.append("")
    out.append(
        "The R-factor is the mean squared difference between the measured and the simulated F(Q) over the fitted Q range, averaged over the data sets."
    )
    out.append("")

    out.append("## Fit quality per data set")
    out.append("")
    has_r2 = any(entry.get("r2") is not None for entry in fit.values())
    header = [
        "data set",
        "file",
        "R initial",
        "R after reference",
        "R final",
        "improvement",
    ]
    if has_r2:
        header += ["R²", "χ²", "Q range (Å⁻¹)"]
    rows = []
    for label, entry in fit.items():
        row = [
            label,
            f"`{entry.get('file', '')}`",
            _sci(entry.get("initial_r_factor")),
            _sci(entry.get("reference_r_factor")),
            _sci(entry.get("final_r_factor")),
            (
                f"×{entry['improvement_factor']:.0f}"
                if entry.get("improvement_factor")
                else "–"
            ),
        ]
        if has_r2:
            q_range = (
                f"{_fixed(entry.get('q_min'), 2)}–{_fixed(entry.get('q_max'), 2)}"
                if entry.get("q_min") is not None
                else "–"
            )
            row += [_fixed(entry.get("r2"), 3), _general(entry.get("chi2"), 3), q_range]
        rows.append(row)
    out.append(_table(header, rows))
    out.append("")
    out.append(
        "R initial: after the scratch start; R after reference: at the end of the last run without the empirical potential; R final: at the end of the last run."
    )
    out.append("")

    out.append("## Convergence")
    out.append("")
    out.append(convergence.get("verdict") or "no convergence information.")
    trends = convergence.get("trends") or {}
    if trends:
        out.append("")
        quarter = convergence.get("quarter_iterations")
        out.append(
            f"Over the last run (run {convergence.get('last_run')}, {convergence.get('last_run_iterations')} iterations; "
            f"means of the first and last {quarter} iterations):"
        )
        out.append("")
        units = {
            "energy": "kJ/mol",
            "pressure": "kbar",
            "ep_amplitude": "kJ/mol",
            "r_factor": "",
        }
        out.append(
            _table(
                ["quantity", "first quarter", "last quarter", "change"],
                [
                    [
                        name.replace("_", " ")
                        + (f" ({units[name]})" if units.get(name) else ""),
                        (
                            _sci(t.get("start"))
                            if name == "r_factor"
                            else _general(t.get("start"), 4)
                        ),
                        (
                            _sci(t.get("end"))
                            if name == "r_factor"
                            else _general(t.get("end"), 4)
                        ),
                        (
                            f"{t['change_percent']:+.1f} %"
                            if t.get("change_percent") is not None
                            else "–"
                        ),
                    ]
                    for name, t in trends.items()
                ],
            )
        )
    out.append("")

    out.append("## First coordination shells (intermolecular g(r), final iteration)")
    out.append("")
    if shells:
        out.append(
            _table(
                [
                    "pair A-B",
                    "r_peak (Å)",
                    "g(r_peak)",
                    "r_min (Å)",
                    "N(B around A)",
                    "N(A around B)",
                ],
                [
                    [
                        pair,
                        _fixed(entry.get("r_peak"), 2),
                        _fixed(entry.get("g_peak"), 2),
                        _fixed(entry.get("r_min"), 2),
                        _fixed(entry.get("n_coord"), 2),
                        (
                            _fixed(entry.get("n_coord_reverse"), 2)
                            if entry.get("n_coord_reverse") is not None
                            else "same"
                        ),
                    ]
                    for pair, entry in shells.items()
                ],
            )
        )
        out.append("")
        out.append(
            "Coordination numbers integrate the intermolecular g(r) to the first minimum (EPSR's .z01 convention)."
        )
    else:
        out.append(
            "No first-shell numbers: neither a trajectory nor the EPSR .g01/.z01 files were available."
        )
    out.append("")

    out.append("## Flags")
    out.append("")
    out += [f"- {flag}" for flag in flags] or ["- none"]
    out.append("")
    return "\n".join(out)


# --------------------------------------------------------------------------- files
def save_analysis(analysis: dict[str, Any], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(plain(analysis), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


def load_analysis(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except ValueError as error:
        raise ValueError(f"{path} is not valid JSON: {error}") from error
    if not isinstance(data, dict):
        raise ValueError(f"{path} does not hold an analysis object")
    return data
