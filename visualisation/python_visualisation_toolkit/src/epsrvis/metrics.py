"""Quality-of-fit and structural metrics computed from a recorded trajectory.

EPSR itself only evaluates the R-factor on the last iteration of each ``epsr``
invocation (it is 0 or stale in between); everything here is evaluated at every
recorded step, with the same definitions where EPSR has one.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.spatial import cKDTree

from .trajectory import Trajectory


# ---------------------------------------------------------------- fit quality

def fit_metrics(traj: Trajectory) -> dict[str, np.ndarray]:
    """Per-step goodness of fit of the simulated F(Q) to each data set.

    Returns arrays of shape [step, dataset] (``*_tot`` keys: [step]):

    * ``rfactor`` - EPSR R-factor: mean of (data - sim)^2 over the fitted Q points,
      using the running-average simulation (identical to ``rfacd`` in the .out file)
    * ``rfactor_inst`` - the same for the current configuration only (EPSR's
      "raw" R-factor): the gap to ``rfactor`` is the configurational noise
    * ``r2`` - coefficient of determination 1 - SS_res / SS_tot over the fitted range
    * ``chi2`` - mean ((data - sim) / error)^2 over points with a non-zero error
    """
    nsteps, nd = traj.nsteps, len(traj.datasets)
    out = {k: np.zeros((nsteps, nd)) for k in ("rfactor", "rfactor_inst", "r2", "chi2")}
    for d in range(nd):
        sl = traj.q_slice(d)
        data = np.asarray(traj.fq_data[:, d, sl], dtype=float)
        sim = np.asarray(traj.fq_sim[:, d, sl], dtype=float)
        inst = (np.asarray(traj.fq_inst[:, d, sl], dtype=float) if traj.fq_inst is not None
                else np.full_like(data, np.nan))
        err = (np.asarray(traj.fq_err[:, d, sl], dtype=float) if traj.fq_err is not None
               else np.zeros_like(data))
        res = data - sim
        out["rfactor"][:, d] = np.mean(res**2, axis=1)
        out["rfactor_inst"][:, d] = np.mean((data - inst) ** 2, axis=1)
        ss_tot = np.sum((data - data.mean(axis=1, keepdims=True)) ** 2, axis=1)
        out["r2"][:, d] = 1.0 - np.sum(res**2, axis=1) / ss_tot
        with np.errstate(divide="ignore", invalid="ignore"):
            z = np.where(err > 0, res / err, np.nan)
        with warnings.catch_warnings():  # no errors recorded -> NaN without a warning
            warnings.simplefilter("ignore", RuntimeWarning)
            out["chi2"][:, d] = np.nanmean(z**2, axis=1)
    out["rfactor_tot"] = out["rfactor"].mean(axis=1)
    out["rfactor_inst_tot"] = out["rfactor_inst"].mean(axis=1)
    return out


def acceptance(traj: Trajectory) -> dict[str, np.ndarray]:
    """Monte Carlo acceptance ratio per move type (last MC cycle of each iteration).

    NaN where the move type was not attempted in that cycle.
    """
    h = traj.history
    out = {}
    if "rot_att" not in h:
        return out
    for move, name in (("rot", "whole-molecule rotation"), ("trans", "whole-molecule translation"),
                       ("intra", "intramolecular shake"), ("irot", "internal rotation")):
        att, rej = h[f"{move}_att"].astype(float), h[f"{move}_rej"].astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            out[name] = np.where(att > 0, 1.0 - rej / att, np.nan)
    return out


def mass_density(traj: Trajectory, masses: dict[str, float]) -> np.ndarray:
    """Mass density (g/cm^3) per step, given the mass (g/mol) of each atom type."""
    m = sum(masses[lab] for lab in traj.atom_labels)
    vol = traj.history["box"] ** 3  # Å^3
    return m / 0.602214076 / vol


# ---------------------------------------------------------------- g(r) features

def shell_volumes(r: np.ndarray) -> np.ndarray:
    """4/3 pi (r_hi^3 - r_lo^3) for each r bin, as EPSR's ``uniform`` (first bin half width)."""
    dr = r[1] - r[0]
    edges = np.concatenate([[0.0], r[1:] - 0.5 * dr, [r[-1] + 0.5 * dr]])
    v = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    v[0] = 0.0
    return v


def running_coordination(traj: Trajectory, step: int, total: bool = False) -> np.ndarray:
    """N_ab(r) for each pair [pair, r]: average number of b atoms within r of an a atom.

    With ``total=False`` (EPSR's ``.z01`` convention) only intermolecular
    neighbours are counted; ``total=True`` adds the intramolecular ones.
    """
    g = traj.gr_total(step) if total else np.asarray(traj.gr_inter[step], dtype=float)
    rho = float(traj.history["rho"][step])
    c = traj.concentrations()
    vol = shell_volumes(traj.r)
    out = np.zeros_like(g)
    for ip, (_, b) in enumerate(traj.meta["pairs"]):
        out[ip] = np.cumsum(g[ip] * vol * rho * c[b])
    return out


def first_shell(r: np.ndarray, g: np.ndarray, coord: np.ndarray, rmin: float = 0.5,
                threshold: float = 0.3, window: float = 1.2) -> dict:
    """Position/height of the first g(r) peak, the following minimum and N(r_min).

    The minimum is the lowest point of the smoothed g(r) within ``window`` Å after the peak.
    """
    gs = _smooth(g)
    above = np.where((r >= rmin) & (gs > threshold))[0]
    if above.size == 0:
        return dict(r_peak=np.nan, g_peak=np.nan, r_min=np.nan, n_coord=np.nan)
    ip = above[0]
    # climb to the local maximum, then take the lowest point within ``window`` Å after it
    while ip + 1 < len(g) and gs[ip + 1] >= gs[ip]:
        ip += 1
    after = np.where((r > r[ip]) & (r <= r[ip] + window))[0]
    im = int(after[np.argmin(gs[after])]) if after.size else ip
    return dict(r_peak=float(r[ip]), g_peak=float(g[ip]), r_min=float(r[im]), n_coord=float(coord[im]))


def _smooth(y: np.ndarray, w: int = 5) -> np.ndarray:
    k = np.ones(w) / w
    return np.convolve(y, k, mode="same")


def shell_history(traj: Trajectory, total: bool = False, every: int = 1) -> dict[str, dict[str, np.ndarray]]:
    """First-shell features of each pair over the run: {pair: {r_peak, g_peak, r_min, n_coord}}.

    ``n_coord`` is integrated to the first minimum of the *last* step for every
    step, so its changes reflect the structure rather than a moving minimum.
    """
    steps = np.arange(0, traj.nsteps, every)
    last = traj.nsteps - 1
    g_last = traj.gr_total(last) if total else np.asarray(traj.gr_inter[last], dtype=float)
    n_last = running_coordination(traj, last, total=total)
    rmin_fixed = [first_shell(traj.r, g_last[ip], n_last[ip])["r_min"] for ip in range(len(traj.pairs))]
    out = {p: {k: np.zeros(len(steps)) for k in ("r_peak", "g_peak", "r_min", "n_coord")} for p in traj.pairs}
    for j, s in enumerate(steps):
        g = traj.gr_total(s) if total else np.asarray(traj.gr_inter[s], dtype=float)
        n = running_coordination(traj, s, total=total)
        for ip, p in enumerate(traj.pairs):
            for k, v in first_shell(traj.r, g[ip], n[ip]).items():
                out[p][k][j] = v
            if np.isfinite(rmin_fixed[ip]):
                out[p]["n_coord"][j] = n[ip][np.searchsorted(traj.r, rmin_fixed[ip])]
    out["_steps"] = steps
    out["_r_min_fixed"] = rmin_fixed
    return out


# ---------------------------------------------------------------- from coordinates

def wrapped(pos: np.ndarray, box: float) -> np.ndarray:
    """Positions folded into [0, L) for periodic neighbour searches."""
    return np.mod(np.asarray(pos, dtype=float) + 0.5 * box, box)


def neighbour_pairs(pos: np.ndarray, box: float, labels: np.ndarray, a: str, b: str, cutoff: float):
    """All (i, j, distance) with label(i)=a, label(j)=b, minimum-image distance < cutoff."""
    w = wrapped(pos, box)
    w = np.minimum(w, np.nextafter(box, 0))  # cKDTree needs data strictly < boxsize
    ia = np.where(labels == a)[0]
    ib = np.where(labels == b)[0]
    tb = cKDTree(w[ib], boxsize=box)
    ta = cKDTree(w[ia], boxsize=box)
    sm = ta.sparse_distance_matrix(tb, cutoff, output_type="coo_matrix")
    i, j, d = ia[sm.row], ib[sm.col], sm.data
    keep = (i != j) & (d > 0)
    return i[keep], j[keep], d[keep]


def coordination_distribution(pos, box, labels, a, b, cutoff, exclude_molecule=None):
    """Histogram of the number of b neighbours (within ``cutoff``) of each a atom.

    ``exclude_molecule``: molecule index per atom; neighbours in the same molecule
    are then not counted (e.g. O...H hydrogen bonds in water).
    Returns (counts n, fraction of a atoms with n neighbours).
    """
    i, j, _ = neighbour_pairs(pos, box, labels, a, b, cutoff)
    if exclude_molecule is not None:
        keep = exclude_molecule[i] != exclude_molecule[j]
        i = i[keep]
    na = np.sum(labels == a)
    per_atom = np.bincount(i, minlength=len(labels))[labels == a]
    n = np.arange(per_atom.max() + 2 if per_atom.size else 1)
    frac = np.bincount(per_atom, minlength=len(n)) / max(na, 1)
    return n, frac[: len(n)]


def bond_angle_distribution(pos, box, labels, centre, outer, cutoff, bins=90, exclude_molecule=None,
                            intramolecular_only=False, molecule=None):
    """Distribution of outer-centre-outer angles (degrees) for neighbours within ``cutoff``.

    Use ``intramolecular_only=True`` with ``molecule`` to restrict to atoms of the
    same molecule (e.g. the H-O-H angle of water).
    """
    pos = np.asarray(pos, dtype=float)
    i, j, _ = neighbour_pairs(pos, box, labels, centre, outer, cutoff)
    if intramolecular_only and molecule is not None:
        keep = molecule[i] == molecule[j]
        i, j = i[keep], j[keep]
    elif exclude_molecule is not None:
        keep = exclude_molecule[i] != exclude_molecule[j]
        i, j = i[keep], j[keep]
    order = np.argsort(i, kind="stable")
    i, j = i[order], j[order]
    angles = []
    starts = np.searchsorted(i, np.unique(i))
    ends = np.append(starts[1:], len(i))
    for s, e in zip(starts, ends):
        if e - s < 2:
            continue
        c = pos[i[s]]
        v = pos[j[s:e]] - c
        v -= box * np.round(v / box)  # minimum image
        v /= np.linalg.norm(v, axis=1, keepdims=True)
        cosm = v @ v.T
        iu = np.triu_indices(e - s, 1)
        angles.append(np.degrees(np.arccos(np.clip(cosm[iu], -1, 1))))
    angles = np.concatenate(angles) if angles else np.zeros(0)
    hist, edges = np.histogram(angles, bins=bins, range=(0, 180), density=True)
    return 0.5 * (edges[1:] + edges[:-1]), hist, angles


def mean_squared_displacement(traj: Trajectory, every: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """MSD (Å^2) of the atoms relative to the first EPSR step, unwrapped across the box.

    EPSR is Monte Carlo, so this is not a dynamical quantity: it shows how far the
    configuration has moved (and so how decorrelated it is) during refinement.
    """
    steps = np.arange(0, traj.nsteps, every)
    ref = np.asarray(traj.positions[0], dtype=float)
    prev = ref.copy()
    disp = np.zeros_like(ref)
    out = np.zeros(len(steps))
    for k, s in enumerate(steps):
        cur = np.asarray(traj.positions[s], dtype=float)
        d = cur - prev
        L = traj.box(s)
        d -= L * np.round(d / L)
        disp += d
        prev = cur
        out[k] = np.mean(np.sum(disp**2, axis=1))
    return steps, out


def tetrahedral_order(pos, box, labels, centre, neighbour, k: int = 4) -> np.ndarray:
    """Tetrahedral order parameter q of every ``centre`` atom (Errington & Debenedetti 2001).

    q = 1 - 3/8 * sum over the 6 pairs (j<l) of the 4 nearest ``neighbour`` atoms
    of (cos psi_jl + 1/3)^2.  q = 1 for a perfect tetrahedron, 0 on average for
    an ideal gas.  Uses the minimum-image convention.
    """
    pos = np.asarray(pos, dtype=float)
    w = np.minimum(wrapped(pos, box), np.nextafter(box, 0))
    ic = np.where(labels == centre)[0]
    inb = np.where(labels == neighbour)[0]
    tree = cKDTree(w[inb], boxsize=box)
    same = centre == neighbour
    _, idx = tree.query(w[ic], k=k + 1 if same else k)
    if same:
        idx = idx[:, 1:]  # drop the atom itself
    v = pos[inb][idx] - pos[ic][:, None, :]
    v -= box * np.round(v / box)
    v /= np.linalg.norm(v, axis=2, keepdims=True)
    cos = np.einsum("nik,njk->nij", v, v)
    iu = np.triu_indices(k, 1)
    return 1.0 - 3.0 / 8.0 * np.sum((cos[:, iu[0], iu[1]] + 1.0 / 3.0) ** 2, axis=1)
