"""Record a running EPSR simulation, iteration by iteration, from Python.

This is the bridge to the Python (PyO3) wrapper of ``epsr-core``.  EPSR only
writes its result files on the last iteration of each ``epsr`` call, so to
follow a refinement step by step the wrapper has to be driven one iteration
at a time and its in-memory results copied after each one::

    import epsr                                   # the PyO3 wrapper
    from epsrvis import TrajectoryRecorder, record_run

    rec = TrajectoryRecorder(case="water")
    for run in range(4):
        sim = epsr.Simulation("work/water", "waterbox.EPSR.inp")   # Epsr::setup
        record_run(sim, rec, pcof="work/water/waterbox.pcof")      # == sim.run(), recorded
    traj = rec.trajectory()        # analyse / plot in memory
    rec.save("work/water/trajectory")

:func:`record_run` does exactly what ``Epsr::run`` does (``iterate(1..=niter)``
then ``write_fort35``), so the simulation files are the same as with
``sim.run()``.  This is what the Rust ``visualisation/recorder`` does, and its
output was checked to be identical to the validated runs.

The wrapper API is not fixed yet, so the recorder is duck-typed (see
:class:`SimulationLike`).  Each quantity is looked up under the names
proposed in ``rust/docs/PYTHON_WRAPPER_NOTES.md`` and under the Rust
accessor names, as a property or a method, per pair/data set or as a 2D
array.  Optional quantities (current-configuration F(Q), Monte Carlo
counters, step sizes) are recorded when available and skipped otherwise.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Protocol, runtime_checkable

import numpy as np

from .trajectory import HISTORY_COLUMNS, Trajectory


@runtime_checkable
class SimulationLike(Protocol):
    """What the recorder needs from a wrapper ``Simulation`` (Rust name in brackets).

    Required
        ``iterate(i)`` [``Epsr::iterate``], ``niter`` [``ctl.niter``],
        Q and r grids [``q_grid()``, ``r_grid()``],
        ``data_fq(id)``, ``simulated_fq(id)``, ``data_q_range(id)``,
        ``gr_inter``, ``gr_intra``, ``sq_inter``, ``sq_intra``,
        ``empirical_potential``, ``reference_potential`` (per pair ``(ip)`` or all pairs ``()``),
        ``pair_labels``, ``atom_types``, ``n_datasets``,
        ``atom_positions()``/``positions`` [Å], ``atom_labels``, ``molecule_atoms``, ``box_size``,
        ``summary`` (dict or object with the ``Summary`` fields).

    Optional
        ``write_fort35()``; ``data_errors(id)``; ``simulated_fq_current(id)``
        (F(Q) of the current configuration); ``mc_counters()`` (2x4 attempts/rejects
        as in the .out file); ``step_sizes()``; ``get_parameter(name, dataset)``.
    """

    def iterate(self, iteration: int) -> Any: ...


# Summary field -> accepted names (Rust field first, then the names suggested for the wrapper).
SUMMARY_ALIASES = {
    "foldsv": ("foldsv", "energy_running", "inter_energy_running"),
    "foldsb": ("foldsb", "energy_intra", "intra_energy"),
    "energd": ("energd", "energy", "inter_energy"),
    "pressd": ("pressd", "pressure"),
    "energabs": ("energabs", "ep_amplitude"),
    "rfactot": ("rfactot", "r_factor"),
    "chisq": ("chisq", "chi_square", "chi2"),
    "sentb": ("sentb", "gr_noise"),
    "epcup": ("epcup", "r_space_coefficient"),
    "epccoeff": ("epccoeff", "quality_factor"),
    "rfacd": ("rfacd", "r_factor_per_dataset"),
}


# ---------------------------------------------------------------- duck-typed access

def _value(obj, *names, args=()):
    """First of ``names`` found on ``obj``: called with ``args`` if callable, else returned."""
    for n in names:
        if hasattr(obj, n):
            v = getattr(obj, n)
            return v(*args) if callable(v) else v
    raise AttributeError(f"{type(obj).__name__} has none of {names}")


def _optional(obj, *names, args=()):
    try:
        return _value(obj, *names, args=args)
    except (AttributeError, TypeError, NotImplementedError):
        return None


def _per_index(obj, names, n):
    """[n, len] array from an accessor that takes an index, or returns all rows at once."""
    for name in names:
        if not hasattr(obj, name):
            continue
        attr = getattr(obj, name)
        if not callable(attr):
            a = np.asarray(attr, dtype=np.float32)
            if a.ndim == 2:
                return a
            continue
        try:
            a = np.asarray(attr(), dtype=np.float32)  # all rows at once: (n, len)
            if a.ndim == 2:
                return a
        except TypeError:
            pass
        return np.stack([np.asarray(attr(i), dtype=np.float32) for i in range(n)])
    raise AttributeError(f"{type(obj).__name__} has none of {names}")


def _summary_get(summary, name):
    for alias in SUMMARY_ALIASES[name]:
        if isinstance(summary, dict) and alias in summary:
            return summary[alias]
        if not isinstance(summary, dict) and hasattr(summary, alias):
            return getattr(summary, alias)
    return None


def read_keyword(path: str | Path, keyword: str) -> str | None:
    """First value of ``keyword`` in an EPSR keyword file (``.pcof``, ``.EPSR.inp``)."""
    try:
        text = Path(path).read_bytes().decode("latin-1")
    except OSError:
        return None
    for line in text.splitlines():
        if line.strip() == "q":
            break
        if line[:12].strip() == keyword:
            parts = line[12:].split()
            return parts[0] if parts else None
    return None


# ---------------------------------------------------------------- recorder

class TrajectoryRecorder:
    """Collects one snapshot per EPSR iteration; see the module docstring.

    Parameters
    ----------
    case : str
        System name stored in the trajectory; selects colours, bonds and the
        structural analyses in :mod:`epsrvis.systems` (``"water"``,
        ``"silica_bulk"``, or anything else for generic defaults).
    record_positions : bool
        Store atom positions every step (needed by the 3D viewer and the
        coordinate-based analyses; ~12 bytes per atom per step).
    """

    def __init__(self, case: str = "", record_positions: bool = True):
        self.case = case
        self.record_positions = record_positions
        self.meta: dict | None = None
        self._hist: dict[str, list] = {}
        self._arr: dict[str, list] = {}
        self._initial_positions = None
        self._run = 0
        self._iter = 0
        self._ereq = 0.0

    # ---- metadata, once
    def _init_meta(self, sim) -> None:
        q = np.asarray(_value(sim, "q", "q_grid"), dtype=float)
        r = np.asarray(_value(sim, "r", "r_grid"), dtype=float)
        nd = int(_value(sim, "n_datasets"))
        pairs = [list(p) for p in _value(sim, "pair_labels")]
        labels = [str(x).strip() for x in _value(sim, "atom_labels")]
        mol = np.zeros(len(labels), dtype=int)
        for im, rng in enumerate(_value(sim, "molecule_atoms")):
            lo, hi = (rng.start, rng.stop) if isinstance(rng, range) else (rng[0], rng[1])
            mol[lo:hi] = im
        names = []
        for d in range(nd):
            n = _optional(sim, "get_parameter", args=("datafile", d))
            names.append(str(n).split()[0] if n else f"dataset{d + 1}")
        self.meta = {
            "case": self.case,
            "nq": len(q), "nr": len(r), "natom": len(labels), "nmol": int(mol.max()) + 1 if len(mol) else 0,
            "npairs": len(pairs), "ndata": nd,
            "q": q.tolist(), "r": r.tolist(),
            "atom_types": [str(t).strip() for t in _value(sim, "atom_types")],
            "pairs": pairs,
            "datasets": names,
            "q_range": [list(map(int, _value(sim, "data_q_range", args=(d,)))) for d in range(nd)],
            "atom_labels": labels,
            "molecule_of_atom": mol.tolist(),
        }
        if self.record_positions:
            self._initial_positions = self._positions(sim)

    def _positions(self, sim):
        return np.asarray(_value(sim, "positions", "atom_positions"), dtype=np.float32).reshape(-1, 3)

    # ---- per run / per iteration
    def start_run(self, sim, ereq: float | None = None, pcof: str | Path | None = None) -> None:
        """Call after creating the ``Simulation`` of a new ``epsr`` invocation.

        ``ereq`` (requested EP amplitude) lives in the ``.pcof``, not in the
        simulation object: pass it, or the ``.pcof`` path to read it from.
        """
        if self.meta is None:
            self._init_meta(sim)
        else:
            nq, nr = len(_value(sim, "q", "q_grid")), len(_value(sim, "r", "r_grid"))
            if (nq, nr) != (self.meta["nq"], self.meta["nr"]):
                raise ValueError("Q or r grid changed between runs; start a new recorder")
        self._run += 1
        self._iter = 0
        if ereq is None and pcof is not None:
            v = read_keyword(pcof, "ereq")
            ereq = float(v) if v is not None else 0.0
        self._ereq = float(ereq or 0.0)

    def record(self, sim) -> None:
        """Snapshot the simulation; call after every ``sim.iterate(i)``."""
        if self.meta is None:
            self.start_run(sim)
        self._iter += 1
        m = self.meta
        nd, npairs = m["ndata"], m["npairs"]
        a = self._arr
        a.setdefault("fq_data", []).append(_per_index(sim, ("data_fq",), nd))
        a.setdefault("fq_err", []).append(_per_index(sim, ("data_errors",), nd) if hasattr(sim, "data_errors")
                                          else np.zeros((nd, m["nq"]), np.float32))
        a.setdefault("fq_sim", []).append(_per_index(sim, ("simulated_fq",), nd))
        inst = None
        if hasattr(sim, "simulated_fq_current"):
            inst = _per_index(sim, ("simulated_fq_current",), nd)
        a.setdefault("fq_inst", []).append(inst if inst is not None else np.full((nd, m["nq"]), np.nan, np.float32))
        a.setdefault("gr_inter", []).append(_per_index(sim, ("gr_inter",), npairs))
        a.setdefault("gr_intra", []).append(_per_index(sim, ("gr_intra",), npairs))
        a.setdefault("sq_inter", []).append(_per_index(sim, ("sq_inter",), npairs))
        a.setdefault("sq_intra", []).append(_per_index(sim, ("sq_intra",), npairs))
        a.setdefault("ep", []).append(_per_index(sim, ("empirical_potential",), npairs))
        a.setdefault("ref_pot", []).append(_per_index(sim, ("reference_potential",), npairs))
        if self.record_positions:
            a.setdefault("pos", []).append(self._positions(sim))

        summary = _value(sim, "summary")
        box = float(_value(sim, "box_size"))
        row = {"step": len(a["fq_sim"]), "run": self._run, "iter": self._iter, "ereq": self._ereq}
        for key in SUMMARY_ALIASES:
            v = _summary_get(summary, key)
            if key == "rfacd":
                row[key] = np.asarray(v if v is not None else np.zeros(nd), dtype=float)
            else:
                row[key] = float(v) if v is not None else np.nan
        row["box"] = box
        row["rho"] = m["natom"] / box**3
        counters = _optional(sim, "mc_counters")
        if counters is not None:
            c = np.asarray(counters).reshape(2, 4)  # rows: attempts, rejects; cols: intra, internal rot, rotation, translation
            for j, mv in enumerate(("intra", "irot", "rot", "trans")):
                row[f"{mv}_att"], row[f"{mv}_rej"] = int(c[0, j]), int(c[1, j])
        steps = _optional(sim, "step_sizes")
        if isinstance(steps, dict):
            for k in ("sizefactor", "stepm", "stepr", "stepmi", "stepri"):
                if k in steps:
                    row[k] = float(steps[k])
        for k, v in row.items():
            self._hist.setdefault(k, []).append(v)

    # ---- results
    @property
    def nsteps(self) -> int:
        return len(self._hist.get("step", []))

    def trajectory(self) -> Trajectory:
        """The recorded steps as an in-memory :class:`~epsrvis.trajectory.Trajectory`."""
        if self.meta is None:
            raise RuntimeError("nothing recorded yet")
        hist = {}
        for k, v in self._hist.items():
            hist[k] = np.stack(v) if k == "rfacd" else np.asarray(v)
        arrays = {k: np.stack(v) for k, v in self._arr.items()}
        meta = dict(self.meta)
        meta["arrays"] = {k: list(v.shape[1:]) for k, v in arrays.items()}
        return Trajectory(meta=meta, history=hist, arrays=arrays, initial_positions=self._initial_positions)

    def save(self, path: str | Path) -> Path:
        """Write the recording in the trajectory folder format."""
        return self.trajectory().save(path)


def record_run(sim, recorder: TrajectoryRecorder, ereq: float | None = None,
               pcof: str | Path | None = None, callback=None) -> TrajectoryRecorder:
    """One full ``epsr`` invocation, recorded: ``iterate(1..=niter)`` + ``write_fort35()``.

    Equivalent to ``sim.run()`` for the simulation files.  ``callback(sim, recorder)``
    is called after each iteration, e.g. to update a live plot.
    """
    recorder.start_run(sim, ereq=ereq, pcof=pcof)
    niter = int(_value(sim, "niter"))
    for i in range(1, niter + 1):
        sim.iterate(i)
        recorder.record(sim)
        if callback is not None:
            callback(sim, recorder)
    fin = getattr(sim, "write_fort35", None)
    if callable(fin):
        fin()
    return recorder


__all__ = ["SimulationLike", "TrajectoryRecorder", "record_run", "read_keyword", "SUMMARY_ALIASES",
           "HISTORY_COLUMNS"]
