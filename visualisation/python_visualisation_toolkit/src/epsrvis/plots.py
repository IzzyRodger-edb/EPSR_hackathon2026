"""Matplotlib figures describing an EPSR refinement, step by step.

Every function takes a :class:`~epsrvis.trajectory.Trajectory` and returns a Figure.
The x axis "EPSR iteration" counts iterations over all ``epsr`` invocations;
the vertical lines mark where a new invocation (stage) starts.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from . import metrics as M
from . import style
from .systems import Case, get_case
from .trajectory import Trajectory

style.apply()


# ---------------------------------------------------------------- helpers

def case_of(traj: Trajectory) -> Case:
    return get_case(traj.case, traj.atom_types)


def key_steps(traj: Trajectory) -> list[int]:
    """Step indices worth comparing: first step and the last step of every run."""
    starts = [k for k, _, _ in traj.run_starts]
    ends = [s - 1 for s in starts[1:]] + [traj.nsteps - 1]
    out = sorted({0, *ends})
    return out


def spread_steps(traj: Trajectory, n: int = 8) -> list[int]:
    """~n step indices, denser at the start where things change fastest."""
    s = np.unique(np.round(np.geomspace(1, traj.nsteps, n)).astype(int) - 1)
    return [int(x) for x in s]


def stage_steps(traj: Trajectory) -> list[int]:
    """Middle and last step of every run of 10 or more iterations (skips the noisy random start)."""
    starts = [k for k, _, _ in traj.run_starts] + [traj.nsteps]
    out = []
    for a, b in zip(starts[:-1], starts[1:]):
        if b - a >= 10:
            out += [a + (b - a) // 2 - 1, b - 1]
    return out or [traj.nsteps - 1]


def step_label(traj: Trajectory, k: int) -> str:
    h = traj.history
    return f"iteration {int(h['step'][k])} (run {int(h['run'][k])}, ereq {h['ereq'][k]:g})"


def _mark_stages(ax, traj: Trajectory, case: Case, label: bool) -> None:
    """Dashed line where each ``epsr`` invocation starts; optionally its stage name above the axes."""
    x = traj.steps
    starts = traj.run_starts
    for n, (k, run, ereq) in enumerate(starts):
        if k > 0:
            ax.axvline(x[k] - 0.5, color=style.MUTED, lw=0.8, ls=(0, (3, 3)), zorder=0)
        end = starts[n + 1][0] if n + 1 < len(starts) else traj.nsteps
        if label and (end - k) >= 0.05 * traj.nsteps:  # skip labels of very short runs (scratch start)
            name = case.stages.get(run, f"run {run}")
            if ereq > 0:
                name += f", ereq {ereq:g}"
            ax.text(x[k] + 1, 1.01, name, transform=ax.get_xaxis_transform(), fontsize=7,
                    color=style.INK_2, va="bottom")
        if label:
            ax.set_title(ax.get_title(loc="left"), loc="left", pad=16)


def _robust_ylim(ax, ys: list[np.ndarray], skip_frac: float = 0.08, pad: float = 0.12) -> bool:
    """Scale the y axis to the values after the first ``skip_frac`` of the steps.

    The scratch start has overlapping atoms (energies ~1e9) that would flatten
    everything else; returns True when early values are off-scale.
    """
    vals = np.concatenate([np.asarray(y)[int(len(y) * skip_frac):] for y in ys])
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return False
    lo, hi = vals.min(), vals.max()
    span = hi - lo or abs(hi) or 1.0
    ax.set_ylim(lo - pad * span, hi + pad * span)
    allv = np.concatenate([np.asarray(y) for y in ys])
    return bool(np.any((allv > hi + pad * span) | (allv < lo - pad * span)))


def _offscale_note(ax, text="earlier values off-scale (random start)"):
    ax.text(0.99, 0.96, text, transform=ax.transAxes, ha="right", va="top", fontsize=7, color=style.MUTED)


# ---------------------------------------------------------------- convergence dashboard

def convergence(traj: Trajectory, masses: dict[str, float] | None = None) -> Figure:
    """Small multiples of the key metrics against EPSR iteration."""
    case = case_of(traj)
    h = traj.history
    x = traj.steps
    fit = M.fit_metrics(traj)
    acc = M.acceptance(traj)
    fig, axes = plt.subplots(4, 2, figsize=(11, 11), sharex=True)
    ax = axes.ravel()

    # 1 intermolecular energy
    a = ax[0]
    a.plot(x, h["energd"], color=style.series(0))
    if _robust_ylim(a, [h["energd"]]):
        _offscale_note(a)
    a.set_title("Intermolecular energy")
    a.set_ylabel("kJ/mol per molecule")

    # 2 intramolecular energy
    a = ax[1]
    a.plot(x, h["foldsb"], color=style.series(0))
    if np.all(h["foldsb"] == 0):
        a.text(0.5, 0.6, "no intramolecular terms (atomic system)", transform=a.transAxes, ha="center",
               color=style.MUTED, fontsize=8)
    else:
        _robust_ylim(a, [h["foldsb"]])
    a.set_title("Intramolecular energy")
    a.set_ylabel("kJ/mol per molecule")

    # 3 pressure
    a = ax[2]
    a.plot(x, h["pressd"], color=style.series(0))
    a.axhline(0, color=style.MUTED, lw=0.8)
    if _robust_ylim(a, [h["pressd"]]):
        _offscale_note(a)
    a.set_title("Intermolecular pressure")
    a.set_ylabel("kbar")

    # 4 R-factor per data set
    a = ax[3]
    for d, name in enumerate(traj.datasets):
        a.plot(x, fit["rfactor"][:, d], color=style.series(d), label=name)
    if np.any(np.isfinite(fit["rfactor_inst_tot"])):
        a.plot(x, fit["rfactor_inst_tot"], color=style.INK_2, lw=1.0, ls=(0, (2, 2)),
               label="current configuration (mean)")
    ends = [k for k in key_steps(traj) if np.nan_to_num(h["rfactot"][k]) > 0]
    a.plot(x[ends], h["rfactot"][ends], "D", ms=5, color=style.INK, mfc="none", label="EPSR total (.out)")
    a.set_yscale("log")
    a.set_title("R-factor  ⟨(data − sim)²⟩")
    a.legend(loc="upper right")

    # 5 R^2 per data set
    a = ax[4]
    for d, name in enumerate(traj.datasets):
        a.plot(x, fit["r2"][:, d], color=style.series(d), label=name)
    a.set_ylim(max(-0.05, min(0.0, np.nanmin(fit["r2"][int(0.08 * len(x)):]) - 0.05)), 1.0)
    a.set_title("R² of the F(Q) fit")
    a.legend(loc="lower right")

    # 6 EP amplitude vs requested
    a = ax[5]
    a.plot(x, h["energabs"], color=style.series(0), label="EP amplitude (energabs)")
    a.step(x, h["ereq"], where="mid", color=style.INK_2, lw=1.0, ls=(0, (4, 2)), label="ereq (requested)")
    a.set_title("Empirical potential amplitude")
    a.set_ylabel("kJ/mol")
    a.legend(loc="upper left")

    # 7 density
    a = ax[6]
    if masses:
        a.plot(x, M.mass_density(traj, masses), color=style.series(0))
        a.set_ylabel("g/cm³")
        a.set_title("Density")
    else:
        a.plot(x, h["rho"], color=style.series(0))
        a.set_ylabel("atoms/Å³")
        a.set_title("Atomic number density")
    rho = h["rho"]
    a.text(0.01, 0.04, f"ρ = {rho[-1]:.5f} atoms/Å³, box {h['box'][-1]:.3f} Å", transform=a.transAxes,
           fontsize=7, color=style.INK_2)
    y = a.lines[0].get_ydata()
    c = float(np.mean(y))
    a.set_ylim(c * (1 - 0.02), c * (1 + 0.02))  # constant-volume MC: a flat line, shown at ±2 %

    # 8 MC acceptance
    a = ax[7]
    rot = acc.get("whole-molecule rotation", np.full(len(x), np.nan))
    trans = acc.get("whole-molecule translation", np.full(len(x), np.nan))
    if acc and np.allclose(np.nan_to_num(rot, nan=-1), np.nan_to_num(trans, nan=-1)):
        # rotation and translation are one combined move, accepted or rejected together
        acc = {"whole-molecule move (rotation + translation)": trans,
               **{k: v for k, v in acc.items() if k not in ("whole-molecule rotation", "whole-molecule translation")}}
    i = 0
    for name, v in acc.items():
        if np.any(np.isfinite(v)):
            a.plot(x, v, color=style.series(i), label=name, lw=1.2)
            i += 1
    if not a.lines:
        a.text(0.5, 0.5, "Monte Carlo counters not recorded", transform=a.transAxes, ha="center",
               color=style.MUTED, fontsize=8)
    a.set_ylim(0, 1)
    a.set_title("Monte Carlo acceptance ratio")
    if a.lines:
        a.legend(loc="lower right")

    for k, a_ in enumerate(ax):
        _mark_stages(a_, traj, case, label=k < 2)
    for a_ in axes[-1]:
        a_.set_xlabel("EPSR iteration")
    fig.suptitle(f"{case.title} — convergence", x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


# ---------------------------------------------------------------- F(Q)

def fq_fit(traj: Trajectory, steps: list[int] | None = None) -> Figure:
    """Measured vs simulated F(Q) for each data set (rows) at selected steps (columns).

    The difference data − simulation is drawn below each curve, offset.
    """
    steps = steps or key_steps(traj)
    nd = len(traj.datasets)
    fit = M.fit_metrics(traj)
    q = traj.q
    fig, axes = plt.subplots(nd, len(steps), figsize=(3.6 * len(steps), 2.9 * nd), squeeze=False,
                             sharex=True, sharey="row")
    for d in range(nd):
        sl = traj.q_slice(d)
        data = traj.fq_data[-1, d, sl]
        off = 1.2 * (np.nanmax(data) - np.nanmin(data)) / 2 + 0.1
        for j, s in enumerate(steps):
            a = axes[d, j]
            sim = traj.fq_sim[s, d, sl]
            a.plot(q[sl], traj.fq_data[s, d, sl], "o", ms=2.2, color=style.INK, label="data", zorder=3)
            a.plot(q[sl], sim, color=style.series(d), label="EPSR", zorder=4)
            a.plot(q[sl], traj.fq_data[s, d, sl] - sim - off, color=style.INK_2, lw=1.0, label="data − EPSR")
            a.axhline(-off, color=style.MUTED, lw=0.6)
            a.text(0.98, 0.95, f"R = {fit['rfactor'][s, d]:.2e}\nR² = {fit['r2'][s, d]:.3f}",
                   transform=a.transAxes, ha="right", va="top", fontsize=7.5, color=style.INK_2)
            if d == 0:
                a.set_title(step_label(traj, s), fontsize=8.5)
            if j == 0:
                a.set_ylabel(f"{traj.datasets[d]}\nF(Q)")
            if d == nd - 1:
                a.set_xlabel("Q (Å⁻¹)")
    axes[0, 0].legend(loc="lower right", fontsize=7)
    fig.suptitle(f"{case_of(traj).title} — fit to the measured structure factors", x=0.01, ha="left",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


def fq_evolution(traj: Trajectory, steps: list[int] | None = None) -> Figure:
    """Simulated F(Q) at successive steps (light -> dark) against the data."""
    steps = steps or spread_steps(traj)
    nd = len(traj.datasets)
    q = traj.q
    cols = style.progress_colors(len(steps))
    fig, axes = plt.subplots(1, nd, figsize=(4.6 * nd, 3.8), squeeze=False)
    for d in range(nd):
        a = axes[0, d]
        sl = traj.q_slice(d)
        for c, s in zip(cols, steps):
            a.plot(q[sl], traj.fq_sim[s, d, sl], color=c, lw=1.2, label=f"iter {int(traj.steps[s])}")
        a.plot(q[sl], traj.fq_data[-1, d, sl], "o", ms=1.8, color=style.INK, label="data", zorder=5)
        a.set_title(traj.datasets[d])
        a.set_xlabel("Q (Å⁻¹)")
        if d == 0:
            a.set_ylabel("F(Q)")
    axes[0, -1].legend(loc="upper right", fontsize=7)
    fig.suptitle("Simulated F(Q) through the refinement", x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def residual_map(traj: Trajectory) -> Figure:
    """Heat map of data − simulation over Q (x) and EPSR iteration (y), per data set."""
    nd = len(traj.datasets)
    fig, axes = plt.subplots(1, nd, figsize=(4.6 * nd, 4.6), squeeze=False)
    q = traj.q
    for d in range(nd):
        a = axes[0, d]
        sl = traj.q_slice(d)
        res = np.asarray(traj.fq_data[:, d, sl]) - np.asarray(traj.fq_sim[:, d, sl])
        lim = np.nanpercentile(np.abs(res[int(0.1 * len(res)):]), 99) or 1e-3
        im = a.pcolormesh(q[sl], traj.steps, res, cmap=style.DIVERGING, vmin=-lim, vmax=lim, shading="nearest",
                          rasterized=True)
        for k, _, _ in traj.run_starts[1:]:
            a.axhline(traj.steps[k] - 0.5, color=style.INK_2, lw=0.6, ls=(0, (3, 3)))
        a.set_title(traj.datasets[d])
        a.set_xlabel("Q (Å⁻¹)")
        a.grid(False)
        if d == 0:
            a.set_ylabel("EPSR iteration")
        fig.colorbar(im, ax=a, label="data − EPSR", shrink=0.85)
    fig.suptitle("Where in Q the fit is still wrong", x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


# ---------------------------------------------------------------- partials and potentials

def partials(traj: Trajectory, steps: list[int] | None = None) -> Figure:
    """Partial g(r) (inter + intra, top) and intermolecular partial S(Q) (bottom) per pair."""
    steps = steps or stage_steps(traj)
    cols = style.progress_colors(len(steps))
    np_ = len(traj.pairs)
    fig, axes = plt.subplots(2, np_, figsize=(4.4 * np_, 6.4), squeeze=False)
    r, q = traj.r, traj.q
    for ip, pair in enumerate(traj.pairs):
        a = axes[0, ip]
        for c, s in zip(cols, steps):
            a.plot(r, traj.gr_total(s)[ip], color=c, lw=1.2, label=f"iter {int(traj.steps[s])}")
        a.axhline(1, color=style.MUTED, lw=0.6)
        a.set_xlim(0, min(12, r[-1]))
        g = traj.gr_total(steps[-1])[ip]
        a.set_ylim(-0.1, min(np.nanmax(g[r > 0.5]) * 1.15, 8) if np.any(g > 0) else 2)
        a.set_title(f"g(r)  {pair}")
        a.set_xlabel("r (Å)")
        b = axes[1, ip]
        for c, s in zip(cols, steps):
            b.plot(q, traj.sq_inter[s, ip], color=c, lw=1.2)
        b.set_xlim(0, min(20, q[-1]))
        b.set_title(f"S(Q) intermolecular  {pair}")
        b.set_xlabel("Q (Å⁻¹)")
    axes[0, -1].legend(loc="upper right", fontsize=7)
    fig.suptitle("Partial pair distribution functions and structure factors", x=0.01, ha="left",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def potentials(traj: Trajectory, steps: list[int] | None = None) -> Figure:
    """Empirical potential per pair through the refinement (top) and its effect on g(r) (bottom)."""
    steps = steps or [k for k in stage_steps(traj) if traj.history["ereq"][k] > 0] or stage_steps(traj)
    cols = style.progress_colors(len(steps))
    np_ = len(traj.pairs)
    r = traj.r
    fig, axes = plt.subplots(2, np_, figsize=(4.4 * np_, 6.0), squeeze=False)
    for ip, pair in enumerate(traj.pairs):
        a = axes[0, ip]
        for c, s in zip(cols, steps):
            a.plot(r, traj.ep[s, ip], color=c, lw=1.2, label=f"iter {int(traj.steps[s])}")
        a.axhline(0, color=style.MUTED, lw=0.6)
        a.set_title(f"Empirical potential  {pair}")
        if ip == 0:
            a.set_ylabel("kJ/mol")
        # effect of the EP on the structure: g(r) before the EP was switched on vs final
        b = axes[1, ip]
        r_ref = [k for k in range(traj.nsteps) if traj.history["ereq"][k] == 0]
        if r_ref:
            g0 = traj.gr_total(r_ref[-1])[ip]
            b.plot(r, g0, color=style.INK_2, lw=1.0, ls=(0, (4, 2)),
                   label=f"reference only (iter {int(traj.steps[r_ref[-1]])})")
        g1 = traj.gr_total(traj.nsteps - 1)[ip]
        b.plot(r, g1, color=style.series(ip), label=f"with EP (iter {int(traj.steps[-1])})")
        b.axhline(1, color=style.MUTED, lw=0.6)
        b.set_xlim(0, min(10, r[-1]))
        b.set_ylim(-0.1, min(np.nanmax(g1[r > 0.5]) * 1.15, 6))
        b.set_title(f"Effect of the EP on g(r)  {pair}", fontsize=9)
        b.set_xlabel("r (Å)")
        if ip == 0:
            b.set_ylabel("g(r)")
        b.legend(loc="upper right", fontsize=7)
    axes[0, -1].legend(loc="upper right", fontsize=7)
    for a in axes[0]:
        a.set_xlim(0, min(12, r[-1]))
    fig.suptitle("Empirical potential", x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


# ---------------------------------------------------------------- structure

def shells(traj: Trajectory, every: int = 2) -> Figure:
    """First-shell peak position, height and coordination number of each pair vs iteration."""
    case = case_of(traj)
    sh = M.shell_history(traj, total=False, every=every)
    x = traj.steps[sh["_steps"]]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), sharex=True)
    for ip, pair in enumerate(traj.pairs):
        d = sh[pair]
        rc = sh["_r_min_fixed"][ip]
        for a, key in zip(axes, ("r_peak", "g_peak", "n_coord")):
            a.plot(x, d[key], color=style.series(ip), lw=1.3,
                   label=f"{pair} (r < {rc:.2f} Å)" if key == "n_coord" else pair)
    axes[0].set_title("First intermolecular peak position (Å)")
    axes[1].set_title("First intermolecular peak height g(r)")
    axes[2].set_title("Coordination number (to final first minimum)")
    for a in axes:
        a.set_xlabel("EPSR iteration")
        _mark_stages(a, traj, case, label=False)
        vals = np.concatenate([l.get_ydata()[len(x) // 10:] for l in a.lines if len(l.get_ydata())])
        vals = vals[np.isfinite(vals)]
        if vals.size:
            span = vals.max() - vals.min() or 1
            a.set_ylim(vals.min() - 0.1 * span, vals.max() + 0.1 * span)
    axes[2].legend(loc="best")
    fig.suptitle("First coordination shell (intermolecular g(r))", x=0.01, ha="left", fontsize=12,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def structure_from_coordinates(traj: Trajectory, steps: list[int] | None = None) -> Figure | None:
    """Coordination-number and bond-angle distributions from the atom positions."""
    case = case_of(traj)
    if not (case.coordination or case.angles):
        return None
    steps = steps or [k for k in key_steps(traj) if k > 0][-3:]
    cols = style.progress_colors(len(steps) + 1)[1:]
    labels = traj.atom_labels
    mol = traj.molecule_of_atom
    panels = len(case.coordination) + len(case.angles)
    fig, axes = plt.subplots(1, panels, figsize=(4.0 * panels, 3.6), squeeze=False)
    axes = axes[0]
    for c, s in zip(cols, steps):
        pos = np.asarray(traj.positions[s], dtype=float)
        box = traj.box(s)
        for k, (a_, b_, cut, inter, desc) in enumerate(case.coordination):
            n, frac = M.coordination_distribution(pos, box, labels, a_, b_, cut,
                                                  exclude_molecule=mol if inter else None)
            w = 0.8 / len(steps)
            off = (steps.index(s) - (len(steps) - 1) / 2) * w
            axes[k].bar(n + off, frac, width=w * 0.9, color=c, label=f"iter {int(traj.steps[s])}")
            axes[k].set_title(desc, fontsize=9)
            axes[k].set_xlabel("number of neighbours")
            axes[k].set_xticks(n)
        for k, (outer, centre, cut, mode, desc) in enumerate(case.angles):
            a = axes[len(case.coordination) + k]
            th, hist, _ = M.bond_angle_distribution(
                pos, box, labels, centre, outer, cut,
                exclude_molecule=mol if mode == "inter" else None,
                intramolecular_only=mode == "intra", molecule=mol, bins=60)
            a.plot(th, hist, color=c, lw=1.4, label=f"iter {int(traj.steps[s])}")
            a.set_title(desc, fontsize=9)
            a.set_xlabel("angle (°)")
    for a in axes:
        a.grid(axis="x", visible=False)
    axes[0].set_ylabel("fraction")
    axes[-1].legend(loc="upper left", fontsize=7)
    fig.suptitle("Local structure from the atomic configuration", x=0.01, ha="left", fontsize=12,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def displacement(traj: Trajectory) -> Figure:
    """How far the atoms have moved since the first iteration (configuration decorrelation)."""
    case = case_of(traj)
    steps, msd = M.mean_squared_displacement(traj)
    fig, a = plt.subplots(figsize=(6, 3.4))
    a.plot(traj.steps[steps], np.sqrt(msd), color=style.series(0))
    a.set_xlabel("EPSR iteration")
    a.set_ylabel("RMS displacement (Å)")
    a.set_title("Atom displacement since the first iteration")
    _mark_stages(a, traj, case, label=False)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- animation

def fq_animation(traj: Trajectory, path, every: int = 2, fps: int = 12):
    """GIF of the F(Q) fit and the partial g(r), one frame per ``every`` iterations."""
    from matplotlib.animation import FuncAnimation, PillowWriter

    nd, np_ = len(traj.datasets), len(traj.pairs)
    q, r = traj.q, traj.r
    fit = M.fit_metrics(traj)
    fig, axes = plt.subplots(2, max(nd, np_), figsize=(4.2 * max(nd, np_), 6.2), squeeze=False)
    sim_lines, txt, gr_lines = [], [], []
    for d in range(nd):
        a = axes[0, d]
        sl = traj.q_slice(d)
        a.plot(q[sl], traj.fq_data[-1, d, sl], "o", ms=1.8, color=style.INK, label="data")
        (ln,) = a.plot(q[sl], traj.fq_sim[0, d, sl], color=style.series(d), label="EPSR")
        sim_lines.append((ln, sl))
        lo, hi = np.nanmin(traj.fq_data[-1, d, sl]), np.nanmax(traj.fq_data[-1, d, sl])
        a.set_ylim(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo))
        a.set_title(traj.datasets[d])
        a.set_xlabel("Q (Å⁻¹)")
        txt.append(a.text(0.98, 0.95, "", transform=a.transAxes, ha="right", va="top", fontsize=8))
    for ip, pair in enumerate(traj.pairs):
        a = axes[1, ip]
        (ln,) = a.plot(r, traj.gr_total(0)[ip], color=style.series(ip))
        gr_lines.append(ln)
        a.axhline(1, color=style.MUTED, lw=0.6)
        g = traj.gr_total(traj.nsteps - 1)[ip]
        a.set_ylim(-0.1, min(np.nanmax(g[r > 0.5]) * 1.2, 8))
        a.set_xlim(0, min(10, r[-1]))
        a.set_title(f"g(r)  {pair}")
        a.set_xlabel("r (Å)")
    for a in axes.ravel()[nd:max(nd, np_)]:
        a.set_visible(False)
    for a in axes[1, np_:]:
        a.set_visible(False)
    title = fig.suptitle("", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    frames = list(range(0, traj.nsteps, every))
    if frames[-1] != traj.nsteps - 1:
        frames.append(traj.nsteps - 1)

    def update(s):
        for d, (ln, sl) in enumerate(sim_lines):
            ln.set_ydata(traj.fq_sim[s, d, sl])
            txt[d].set_text(f"R = {fit['rfactor'][s, d]:.2e}   R² = {fit['r2'][s, d]:.3f}")
        g = traj.gr_total(s)
        for ip, ln in enumerate(gr_lines):
            ln.set_ydata(g[ip])
        title.set_text(step_label(traj, s))
        return []

    anim = FuncAnimation(fig, update, frames=frames, blit=False)
    anim.save(path, writer=PillowWriter(fps=fps))
    plt.close(fig)


# ---------------------------------------------------------------- summary table

def summary_rows(traj: Trajectory, masses: dict[str, float] | None = None) -> list[dict]:
    """One row of key numbers at the end of each run."""
    fit = M.fit_metrics(traj)
    h = traj.history
    rows = []
    for s in key_steps(traj):
        n = M.running_coordination(traj, s, total=False)
        g = np.asarray(traj.gr_inter[s], dtype=float)
        row = {
            "iteration": int(h["step"][s]), "run": int(h["run"][s]), "ereq": float(h["ereq"][s]),
            "energy_inter": float(h["energd"][s]), "energy_intra": float(h["foldsb"][s]),
            "pressure": float(h["pressd"][s]), "ep_amplitude": float(h["energabs"][s]),
            "rho_atoms": float(h["rho"][s]),
            "rfactor_total": float(fit["rfactor_tot"][s]),
        }
        if masses:
            row["density_g_cm3"] = float(M.mass_density(traj, masses)[s])
        for d, name in enumerate(traj.datasets):
            row[f"rfactor[{name}]"] = float(fit["rfactor"][s, d])
            row[f"R2[{name}]"] = float(fit["r2"][s, d])
            row[f"chi2[{name}]"] = float(fit["chi2"][s, d])
        for ip, pair in enumerate(traj.pairs):
            fs = M.first_shell(traj.r, g[ip], n[ip])
            row[f"r_peak[{pair}]"] = fs["r_peak"]
            row[f"N_coord[{pair}]"] = fs["n_coord"]
            row[f"r_min[{pair}]"] = fs["r_min"]
        rows.append(row)
    return rows


# ---------------------------------------------------------------- final EPSR output files only

def final_outputs(folder) -> Figure:
    """Overview of a finished EPSR run from its own output files (no trajectory needed).

    Top: F(Q) data (``.t01``) vs simulation (``.u01``) per data set.
    Bottom: partial g(r) (``.g01`` + ``.y01``), empirical potential (``.p01``),
    and the ``.erg`` history (energy, R-factor at run ends).
    """
    from .outputs import OutputFolder

    f = OutputFolder.find(folder)
    t, u = f.table("t"), f.table("u")
    g, y, p = f.table("g"), f.table("y"), f.table("p")
    erg = f.erg()
    nd = len(t.labels)
    fig = plt.figure(figsize=(4.4 * max(nd, 3), 7.4))
    gs = fig.add_gridspec(2, max(nd, 3))
    for d, name in enumerate(t.labels):
        a = fig.add_subplot(gs[0, d])
        m = t.values[d] != 0
        a.plot(t.x[m], t.values[d][m], "o", ms=2, color=style.INK, label="data (.t01)")
        ui = np.interp(t.x[m], u.x, u.values[d])
        a.plot(t.x[m], ui, color=style.series(d), label="EPSR (.u01)")
        a.set_title(name)
        a.set_xlabel("Q (Å⁻¹)")
        if d == 0:
            a.set_ylabel("F(Q)")
            a.legend(loc="lower right")
    a = fig.add_subplot(gs[1, 0])
    for ip, pair in enumerate(g.labels):
        a.plot(g.x, g.values[ip] + y.values[ip], color=style.series(ip), label=pair)
    a.axhline(1, color=style.MUTED, lw=0.6)
    a.set_xlim(0, 10)
    a.set_ylim(-0.1, min(6, 1.1 * np.max((g.values + y.values)[:, g.x > 0.5])))
    a.set_title("Partial g(r) (.g01 + .y01)")
    a.set_xlabel("r (Å)")
    a.legend()
    a = fig.add_subplot(gs[1, 1])
    for ip, pair in enumerate(p.labels):
        a.plot(p.x, p.values[ip], color=style.series(ip), label=pair)
    a.axhline(0, color=style.MUTED, lw=0.6)
    a.set_xlim(0, 12)
    a.set_title("Empirical potential (.p01)")
    a.set_xlabel("r (Å)")
    a.set_ylabel("kJ/mol")
    a = fig.add_subplot(gs[1, 2])
    it = np.arange(1, len(erg["energy"]) + 1)
    a.plot(it, erg["energy"], color=style.series(0))
    if _robust_ylim(a, [erg["energy"]]):
        _offscale_note(a)
    a.set_title("Energy per molecule (.erg)")
    a.set_xlabel("line of .erg (iteration)")
    a.set_ylabel("kJ/mol")
    fig.suptitle(f"{f.prefix}: final EPSR output", x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig
