"""The per-iteration record of an EPSR refinement.

A :class:`Trajectory` holds, for every recorded EPSR iteration ("step"),
the measured and simulated F(Q), partial g(r)/S(Q), potentials, atom
positions and the scalar convergence numbers.  It can come from

* a :class:`~epsrvis.recorder.TrajectoryRecorder` fed by the Python wrapper
  (arrays in memory), or
* a folder on disk written by :meth:`Trajectory.save` or by the Rust
  ``epsr-recorder`` (arrays memory-mapped, so large runs are cheap to open).

Disk format (one folder)::

    meta.json      grids, labels, sizes, fitted Q range per data set,
                   molecule index of every atom, and the shape of each array
    history.tsv    one line per step (tab separated, header line first)
    <name>.f32     little-endian float32, one [a, b] block per step
    pos_initial.f32  [natom, 3], optional: configuration before step 1

Array names and shapes per step (see :data:`ARRAYS`)::

    fq_data, fq_err, fq_sim, fq_inst   [n_datasets, nq]
    gr_inter, gr_intra, ep, ref_pot    [n_pairs, nr]
    sq_inter, sq_intra                 [n_pairs, nq]
    pos                                [n_atoms, 3]
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

#: Per-step arrays: name -> (first axis kind, second axis kind, description).
ARRAYS = {
    "fq_data": ("dataset", "q", "measured F(Q) on the simulation Q grid"),
    "fq_err": ("dataset", "q", "error of the measured F(Q)"),
    "fq_sim": ("dataset", "q", "simulated F(Q), running average (what EPSR fits)"),
    "fq_inst": ("dataset", "q", "simulated F(Q) of the current configuration only"),
    "gr_inter": ("pair", "r", "intermolecular partial g(r), running average"),
    "gr_intra": ("pair", "r", "intramolecular partial g(r), running average"),
    "sq_inter": ("pair", "q", "intermolecular partial S(Q), running average"),
    "sq_intra": ("pair", "q", "intramolecular partial S(Q), running average"),
    "ep": ("pair", "r", "empirical potential (kJ/mol)"),
    "ref_pot": ("pair", "r", "reference potential, Lennard-Jones + Coulomb (kJ/mol)"),
    "pos": (
        "atom",
        "xyz",
        "atom positions (Å); box centred on the origin, molecules whole",
    ),
}

#: history.tsv columns.  Names follow the Rust ``Summary`` fields.
HISTORY_COLUMNS = {
    "step": "EPSR iteration, counted over all runs (1-based)",
    "run": "index of the epsr invocation (setup + run) the step belongs to (1-based)",
    "iter": "iteration within that run (1-based)",
    "ereq": "requested EP amplitude (kJ/mol), from the .pcof",
    "foldsv": "intermolecular energy per molecule, running (kJ/mol)",
    "foldsb": "intramolecular energy per molecule (kJ/mol)",
    "energd": "intermolecular energy per molecule from the pressure routine (kJ/mol)",
    "pressd": "intermolecular pressure (kbar)",
    "energabs": "EP amplitude actually reached (kJ/mol)",
    "rfactot": "EPSR R-factor; only computed on the last iteration of a run, else 0",
    "chisq": "EPSR chi-square of the partial S(Q); last iteration of a run only",
    "sentb": "noise of the averaged g(r)",
    "epcup": "r-space quality coefficient (last iteration of a run only)",
    "epccoeff": "combined quality factor (last iteration of a run only)",
    "box": "cubic box side (Å)",
    "rho": "atomic number density (atoms/Å^3)",
    "sizefactor": "EPSR sizefactor",
    "stepm": "molecule translation step (Å)",
    "stepr": "molecule rotation step",
    "stepmi": "intramolecular translation step (Å)",
    "stepri": "internal rotation step",
    "intra_att": "intramolecular moves attempted (last MC cycle)",
    "intra_rej": "intramolecular moves rejected",
    "irot_att": "internal rotations attempted",
    "irot_rej": "internal rotations rejected",
    "rot_att": "whole-molecule rotations attempted",
    "rot_rej": "whole-molecule rotations rejected",
    "trans_att": "whole-molecule translations attempted",
    "trans_rej": "whole-molecule translations rejected",
    "rfacd": "EPSR R-factor per data set (comma separated in the file)",
}

_INT_COLS = {
    "step",
    "run",
    "iter",
    "intra_att",
    "intra_rej",
    "irot_att",
    "irot_rej",
    "rot_att",
    "rot_rej",
    "trans_att",
    "trans_rej",
}


@dataclass
class Trajectory:
    """Per-iteration record of an EPSR refinement.

    Arrays have a leading *step* axis, e.g. ``fq_sim[step, dataset, q]`` or
    ``gr_inter[step, pair, r]``.  Step index ``k`` (0-based) is EPSR iteration
    ``history["step"][k]``.  Missing arrays are ``None``; missing history
    columns are simply absent from :attr:`history`.
    """

    meta: dict = field(repr=False)
    history: dict[str, np.ndarray] = field(repr=False)
    arrays: dict[str, np.ndarray] = field(repr=False)
    path: Path | None = None
    initial_positions: np.ndarray | None = field(default=None, repr=False)

    # ---------------------------------------------------------------- I/O
    @classmethod
    def load(cls, path: str | Path) -> Trajectory:
        """Open a trajectory folder (arrays are memory-mapped, read-only)."""
        path = Path(path)
        meta = json.loads((path / "meta.json").read_text())
        history = read_history(path / "history.tsv")
        nsteps = len(history["step"])
        arrays = {}
        for name, (a, b) in meta.get("arrays", {}).items():
            f = path / f"{name}.f32"
            if not f.exists():
                continue
            mm = np.memmap(f, dtype="<f4", mode="r")
            n = min(
                mm.size // (a * b), nsteps
            )  # a run still being recorded can be ahead of history.tsv
            arrays[name] = mm[: n * a * b].reshape(n, a, b)
        init = path / "pos_initial.f32"
        pos0 = np.fromfile(init, dtype="<f4").reshape(-1, 3) if init.exists() else None
        return cls(
            meta=meta, history=history, arrays=arrays, path=path, initial_positions=pos0
        )

    def save(self, path: str | Path) -> Path:
        """Write the trajectory in the folder format (readable by :meth:`load`)."""
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        meta = dict(self.meta)
        meta["arrays"] = {k: list(np.shape(v)[1:]) for k, v in self.arrays.items()}
        (path / "meta.json").write_text(json.dumps(meta, indent=1))
        write_history(path / "history.tsv", self.history)
        for name, a in self.arrays.items():
            np.ascontiguousarray(a, dtype="<f4").tofile(path / f"{name}.f32")
        if self.initial_positions is not None:
            np.ascontiguousarray(self.initial_positions, dtype="<f4").tofile(
                path / "pos_initial.f32"
            )
        return path

    # ---------------------------------------------------------------- labels and sizes
    @property
    def case(self) -> str:
        """System name (selects the settings in :mod:`epsrvis.systems`)."""
        return self.meta.get("case", "")

    @property
    def nsteps(self) -> int:
        return len(self.history["step"])

    @property
    def steps(self) -> np.ndarray:
        """EPSR iteration number of each step."""
        return self.history["step"]

    @property
    def q(self) -> np.ndarray:
        return np.asarray(self.meta["q"], dtype=float)

    @property
    def r(self) -> np.ndarray:
        return np.asarray(self.meta["r"], dtype=float)

    @property
    def pairs(self) -> list[str]:
        """Pair labels such as ``"OW-HW"``, in the order of the pair axis."""
        return [f"{a}-{b}" for a, b in self.meta["pairs"]]

    @property
    def datasets(self) -> list[str]:
        """Data set names (data file name without extension)."""
        return [Path(d).stem for d in self.meta["datasets"]]

    @property
    def atom_types(self) -> list[str]:
        return list(self.meta["atom_types"])

    @property
    def atom_labels(self) -> np.ndarray:
        return np.asarray(self.meta["atom_labels"])

    @property
    def molecule_of_atom(self) -> np.ndarray:
        return np.asarray(self.meta["molecule_of_atom"], dtype=int)

    def q_slice(self, dataset: int) -> slice:
        """Q-grid indices EPSR fits for ``dataset``."""
        lo, hi = self.meta["q_range"][dataset]
        return slice(lo, hi + 1)

    @property
    def run_starts(self) -> list[tuple[int, int, float]]:
        """(first step index, run number, ereq) of each EPSR invocation."""
        run = self.history["run"]
        ereq = self.history.get("ereq", np.zeros(self.nsteps))
        return [
            (k, int(run[k]), float(ereq[k]))
            for k in range(self.nsteps)
            if k == 0 or run[k] != run[k - 1]
        ]

    def concentrations(self) -> dict[str, float]:
        """Atom fraction of each atom type."""
        labels = self.atom_labels
        return {t: float(np.mean(labels == t)) for t in self.atom_types}

    # ---------------------------------------------------------------- arrays
    def __getattr__(self, name):
        # fq_sim, gr_inter, ... ; avoids recursion before the dataclass fields exist
        if name in ARRAYS:
            return self.__dict__.get("arrays", {}).get(name)
        raise AttributeError(name)

    @property
    def positions(self) -> np.ndarray | None:
        return self.arrays.get("pos")

    def has(self, name: str) -> bool:
        return name in self.arrays and len(self.arrays[name]) == self.nsteps

    def box(self, step: int) -> float:
        """Cubic box side (Å) at step index ``step``; the box spans [-L/2, L/2)."""
        return float(self.history["box"][step])

    def gr_total(self, step: int) -> np.ndarray:
        """Partial g(r) including the intramolecular part, [pair, r]."""
        g = np.asarray(self.arrays["gr_inter"][step], dtype=float)
        if "gr_intra" in self.arrays:
            g = g + np.asarray(self.arrays["gr_intra"][step], dtype=float)
        return g

    def subset(self, steps) -> Trajectory:
        """A new in-memory trajectory with only the given step indices."""
        idx = np.asarray(steps)
        hist = {k: np.asarray(v)[idx] for k, v in self.history.items()}
        arrs = {k: np.asarray(v)[idx] for k, v in self.arrays.items()}
        return Trajectory(
            meta=self.meta,
            history=hist,
            arrays=arrs,
            initial_positions=self.initial_positions,
        )


# ---------------------------------------------------------------- history.tsv


def read_history(path: str | Path) -> dict[str, np.ndarray]:
    lines = Path(path).read_text().splitlines()
    head = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:] if ln.strip()]
    rows = [
        r for r in rows if len(r) == len(head)
    ]  # the last line may be incomplete while recording
    out: dict[str, np.ndarray] = {}
    for j, name in enumerate(head):
        col = [r[j] for r in rows]
        if name == "rfacd":
            out[name] = (
                np.array([[float(x) for x in c.split(",")] for c in col])
                if col
                else np.zeros((0, 0))
            )
        elif name in _INT_COLS:
            out[name] = np.array([int(float(x)) for x in col], dtype=int)
        else:
            out[name] = np.array([float(x) for x in col])
    return out


def write_history(path: str | Path, history: dict[str, np.ndarray]) -> None:
    cols = [c for c in HISTORY_COLUMNS if c in history] + [
        c for c in history if c not in HISTORY_COLUMNS
    ]
    n = len(history["step"])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for k in range(n):
            vals = []
            for c in cols:
                v = history[c][k]
                if c == "rfacd":
                    vals.append(",".join(f"{x:e}" for x in np.atleast_1d(v)))
                elif c in _INT_COLS:
                    vals.append(str(int(v)))
                else:
                    vals.append(f"{float(v):e}")
            fh.write("\t".join(vals) + "\n")
