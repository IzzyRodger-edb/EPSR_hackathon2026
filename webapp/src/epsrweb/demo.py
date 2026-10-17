"""Demonstration mode: a run folder with the right shape and no physics.

The Mellor package is a sibling of this app in the same repository and needs a
Rust build. When it is not importable the app still has to be usable -- for the
UI, for the Overleaf path, for tests -- so this module writes a run folder that
obeys the same contract as a real one (integration.md section 5): the same file
names, the same CSV columns, the same `analysis.json` keys.

The numbers in it are not a refinement. The measured F(Q) columns are the user's
own uploaded data; everything else is a smooth curve chosen to look like the
right kind of plot. Every artefact it writes is stamped DEMONSTRATION, the run is
flagged `demo` in the database, and the UI says so on the results page. Nothing
here should ever be read as a result.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pandas as pd

BANNER = "DEMONSTRATION DATA - Mellor is not installed; these numbers are not a refinement."

#: Guessed atom pairs, so the partials plot carries plausible labels.
_PAIRS_BY_KEYWORD = {
    "sio2": ("Si-Si", "Si-O", "O-O"),
    "silica": ("Si-Si", "Si-O", "O-O"),
    "benzene": ("C-C", "C-H", "H-H"),
    "dmso": ("C-C", "C-O", "O-H"),
    "cubtc": ("Cu-Cu", "Cu-O", "O-O"),
}
_DEFAULT_PAIRS = ("O-O", "O-H", "H-H")


def _pairs_for(titles: Sequence[str]) -> tuple[str, ...]:
    joined = " ".join(titles).lower().replace(" ", "").replace("-", "")
    for keyword, pairs in _PAIRS_BY_KEYWORD.items():
        if keyword in joined:
            return pairs
    return _DEFAULT_PAIRS


def _read_mint01(path: Path) -> tuple[str, np.ndarray, np.ndarray, np.ndarray]:
    headers: list[str] = []
    rows: list[tuple[float, float, float]] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if "#" in line.split()[0]:
            headers.append(line.lstrip("#").strip())
            continue
        fields = line.split()
        if len(fields) >= 3:
            try:
                rows.append((float(fields[0]), float(fields[1]), float(fields[2])))
            except ValueError:
                continue
    title = headers[1] if len(headers) > 1 else path.stem
    array = np.asarray(rows, dtype=float)
    return title.strip() or path.stem, array[:, 0], array[:, 1], array[:, 2]


def _pseudo_fit(q: np.ndarray, data: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """A curve that follows the data the way a converged simulation would.

    A short rolling mean (EPSR's running average smooths configurational noise)
    plus a small Q-dependent error, so that the residual plot is not identically
    zero and the fit plot is not a perfect overlay.
    """
    smooth = pd.Series(data).rolling(7, center=True, min_periods=1).mean().to_numpy()
    drift = 0.012 * np.exp(-q / 6.0) * np.sin(1.7 * q + 0.4)
    return smooth + drift + rng.normal(0.0, 0.0015, size=q.size)


def _write_tables(run_dir: Path, datasets: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]],
                  pairs: tuple[str, ...], steps: int) -> None:
    rng = np.random.default_rng(20261017)
    tables = run_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    fits: list[dict] = []
    for title, q, data, err in datasets:
        simulated = _pseudo_fit(q, data, rng)
        keep = q <= 30.0  # above nq*qstep EPSR discards the data
        for qi, di, ei, si in zip(q[keep], data[keep], err[keep], simulated[keep]):
            fits.append({"dataset": title, "q": qi, "data": di, "error": ei,
                         "simulated": si, "residual": di - si})
    pd.DataFrame(fits).to_csv(tables / "fits.csv", index=False)

    r = np.arange(0.0, 20.0, 0.02)
    first_peak = {"O-O": 2.80, "O-H": 1.83, "H-H": 2.35, "Si-Si": 3.09, "Si-O": 1.59,
                  "C-C": 1.40, "C-H": 1.09, "C-O": 1.43, "Cu-Cu": 2.63, "Cu-O": 1.95}
    partial_rows: list[dict] = []
    potential_rows: list[dict] = []
    for pair in pairs:
        r0 = first_peak.get(pair, 2.6)
        inter = 1.0 - np.exp(-((r / (0.8 * r0)) ** 8)) + 1.1 * np.exp(-((r - r0) / 0.18) ** 2) \
            + 0.35 * np.exp(-((r - 1.75 * r0) / 0.5) ** 2)
        inter[r < 0.6 * r0] = 0.0
        intra = 2.4 * np.exp(-((r - r0) / 0.06) ** 2) if pair.count("-") and r0 < 2.0 else r * 0.0
        for ri, gi, gj in zip(r, inter, intra):
            partial_rows.append({"pair": pair, "r": ri, "g_inter": gi, "g_intra": gj,
                                 "g_total": gi + gj})
        sigma, epsilon = 0.89 * r0, 0.65
        safe = np.maximum(r, 0.4)
        reference = 4 * epsilon * ((sigma / safe) ** 12 - (sigma / safe) ** 6)
        empirical = 1.4 * np.exp(-((r - r0) / 0.25) ** 2) - 0.9 * np.exp(-((r - 1.4 * r0) / 0.4) ** 2)
        for ri, ei, refi in zip(r, empirical, np.clip(reference, -5.0, 20.0)):
            potential_rows.append({"pair": pair, "r": ri, "empirical": ei, "reference": refi})
    pd.DataFrame(partial_rows).to_csv(tables / "partials.csv", index=False)
    pd.DataFrame(potential_rows).to_csv(tables / "potentials.csv", index=False)

    stages = [("scratch", 0.0), ("equilibrate", 0.0), ("refine", 3.0)]
    per_run = max(5, steps // len(stages))
    history: list[dict] = []
    runs: list[dict] = []
    step = 0
    r_factor = 1.1e-3
    for index, (stage, ereq) in enumerate(stages, start=1):
        for iteration in range(1, per_run + 1):
            step += 1
            decay = np.exp(-iteration / (per_run / 2.2))
            r_factor = (2.2e-4 + (1.1e-3 - 2.2e-4) * decay) if ereq else (6.1e-4 + 2e-4 * decay)
            amplitude = ereq * (1 - decay) if ereq else 0.0
            row = {"step": step, "run": f"{index}_{stage}", "iteration": iteration,
                   "ereq": ereq, "energy": -35.0 - 10.0 * (1 - decay),
                   "pressure": 5.5 - 2.0 * (1 - decay), "ep_amplitude": amplitude,
                   "r_factor": r_factor}
            for title, *_ in datasets:
                row[f"r_factor[{title}]"] = r_factor * (0.8 + 0.4 * rng.random())
                row[f"r2[{title}]"] = 1.0 - 25.0 * r_factor
            history.append(row)
        runs.append({"run": f"{index}_{stage}", "stage": stage, "niter": per_run, "ereq": ereq,
                     "r_factor": r_factor, "energy": -35.0 - 10.0 * index,
                     "pressure": 5.5 - 0.6 * index, "ep_amplitude": ereq,
                     "duration_s": 12.0 * index})
    pd.DataFrame(history).to_csv(tables / "history.csv", index=False)
    pd.DataFrame(runs).to_csv(tables / "runs.csv", index=False)


def _write_figures(run_dir: Path) -> list[str]:
    """Four of Mellor's nine figure names, drawn from the tables just written."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = run_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    tables = run_dir / "tables"
    fits = pd.read_csv(tables / "fits.csv")
    history = pd.read_csv(tables / "history.csv")
    partials = pd.read_csv(tables / "partials.csv")
    potentials = pd.read_csv(tables / "potentials.csv")
    names: list[str] = []

    fig, axes = plt.subplots(2, 2, figsize=(9, 6))
    axes[0, 0].semilogy(history["step"], history["r_factor"])
    axes[0, 0].set(xlabel="iteration", ylabel="R-factor")
    axes[0, 1].plot(history["step"], history["energy"])
    axes[0, 1].set(xlabel="iteration", ylabel="energy (kJ/mol per molecule)")
    axes[1, 0].plot(history["step"], history["pressure"])
    axes[1, 0].set(xlabel="iteration", ylabel="pressure (kbar)")
    axes[1, 1].plot(history["step"], history["ep_amplitude"])
    axes[1, 1].set(xlabel="iteration", ylabel="EP amplitude (kJ/mol)")
    fig.suptitle(BANNER, fontsize=8, color="#b4432f")
    fig.tight_layout()
    fig.savefig(figures / "01_convergence.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    names.append("01_convergence")

    datasets = list(dict.fromkeys(fits["dataset"]))
    fig, axes = plt.subplots(len(datasets), 1, figsize=(7, 2.6 * len(datasets)), squeeze=False)
    for ax, title in zip(axes[:, 0], datasets):
        part = fits[fits["dataset"] == title]
        ax.plot(part["q"], part["data"], color="k", linewidth=0.8, label="data")
        ax.plot(part["q"], part["simulated"], color="#2f6fb4", linewidth=1.0, label="EPSR")
        ax.plot(part["q"], part["residual"] - part["data"].max() * 0.6, color="#b4432f",
                linewidth=0.7, label="difference")
        ax.set(ylabel="F(Q)", title=title)
        ax.legend(fontsize=7, frameon=False)
    axes[-1, 0].set_xlabel("Q (1/A)")
    fig.suptitle(BANNER, fontsize=8, color="#b4432f")
    fig.tight_layout()
    fig.savefig(figures / "02_fq_fit.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    names.append("02_fq_fit")

    fig, (left, right) = plt.subplots(1, 2, figsize=(9, 3.4))
    for pair, group in partials.groupby("pair"):
        near = group[group["r"] <= 8]
        left.plot(near["r"], near["g_total"], label=pair)
    left.axhline(1.0, color="0.6", linestyle=":", linewidth=0.8)
    left.set(xlabel="r (A)", ylabel="g(r)")
    left.legend(fontsize=7, frameon=False)
    for pair, group in potentials.groupby("pair"):
        near = group[group["r"] <= 8]
        right.plot(near["r"], near["empirical"], label=pair)
    right.axhline(0.0, color="0.6", linewidth=0.8)
    right.set(xlabel="r (A)", ylabel="empirical potential (kJ/mol)")
    right.legend(fontsize=7, frameon=False)
    fig.suptitle(BANNER, fontsize=8, color="#b4432f")
    fig.tight_layout()
    fig.savefig(figures / "05_partials.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    names.append("05_partials")

    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    for pair, group in potentials.groupby("pair"):
        near = group[group["r"] <= 8]
        ax.plot(near["r"], near["reference"], linestyle="--", linewidth=0.9,
                label=f"{pair} reference")
        ax.plot(near["r"], near["empirical"], linewidth=1.2, label=f"{pair} EP")
    ax.set(xlabel="r (A)", ylabel="energy (kJ/mol)", ylim=(-6, 8))
    ax.legend(fontsize=6, frameon=False, ncols=2)
    fig.suptitle(BANNER, fontsize=8, color="#b4432f")
    fig.tight_layout()
    fig.savefig(figures / "06_potentials.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    names.append("06_potentials")

    runs = pd.read_csv(tables / "runs.csv")
    runs.to_csv(figures / "summary.csv", index=False)
    (figures / "summary.md").write_text(
        f"<!-- {BANNER} -->\n\n" + _markdown_table(runs), encoding="utf-8"
    )
    return names


def _markdown_table(frame: pd.DataFrame) -> str:
    """A GitHub-flavoured table without pulling in `tabulate`."""
    header = "| " + " | ".join(str(c) for c in frame.columns) + " |"
    rule = "|" + "|".join("---" for _ in frame.columns) + "|"
    rows = [
        "| " + " | ".join(
            f"{v:.4g}" if isinstance(v, float) else str(v) for v in record
        ) + " |"
        for record in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, rule, *rows])


def _write_analysis(run_dir: Path, datasets: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]],
                    pairs: tuple[str, ...], system: str) -> dict:
    runs = pd.read_csv(run_dir / "tables" / "runs.csv")
    history = pd.read_csv(run_dir / "tables" / "history.csv")
    final = history.iloc[-1]
    first_peak = {"O-O": (2.80, 5.1), "O-H": (1.83, 1.76), "H-H": (2.35, 4.2),
                  "Si-Si": (3.09, 4.0), "Si-O": (1.59, 3.88), "C-C": (1.40, 2.0),
                  "C-H": (1.09, 1.0), "C-O": (1.43, 1.0), "Cu-Cu": (2.63, 1.0),
                  "Cu-O": (1.95, 4.0)}
    analysis = {
        "demo": True,
        "banner": BANNER,
        "system": system,
        "runs": runs.to_dict(orient="records"),
        "fit": [
            {"dataset": title,
             "r_factor": float(final.get(f"r_factor[{title}]", final["r_factor"])),
             "r2": float(final.get(f"r2[{title}]", 0.95)),
             "q_min": float(q.min()), "q_max": float(min(q.max(), 30.0)),
             "points": int(q.size)}
            for title, q, _data, _err in datasets
        ],
        "convergence": {
            "verdict": "demonstration output; no refinement was performed",
            "r_factor_final": float(final["r_factor"]),
            "ep_amplitude": float(final["ep_amplitude"]),
            "ereq": float(final["ereq"]),
            "at_ereq": bool(final["ep_amplitude"] >= 0.98 * max(final["ereq"], 1e-9)),
        },
        "shells": [
            {"pair": pair, "position": first_peak.get(pair, (2.6, 4.0))[0],
             "coordination": first_peak.get(pair, (2.6, 4.0))[1], "cutoff": 3.45}
            for pair in pairs
        ],
        "flags": ["demo_mode: Mellor is not installed; every number here is synthetic"],
    }
    (run_dir / "analysis.json").write_text(json.dumps(analysis, indent=2), encoding="utf-8")
    return analysis


def _write_report(run_dir: Path, analysis: dict, figure_names: list[str], title: str) -> Path:
    """A PDF with the same role as Mellor's `report.pdf`, built with matplotlib."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    path = run_dir / "report.pdf"
    runs = pd.DataFrame(analysis["runs"])
    with PdfPages(path) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.08, 0.92, "EPSR refinement report", fontsize=22)
        fig.text(0.08, 0.885, title, fontsize=12, color="0.3")
        fig.text(0.08, 0.85, BANNER, fontsize=9, color="#b4432f", wrap=True)
        lines = [
            "",
            "Fit quality per data set",
            *[f"   {row['dataset']}: R-factor {row['r_factor']:.3g}, R2 {row['r2']:.3f}, "
              f"Q {row['q_min']:.2f} to {row['q_max']:.1f} 1/A"
              for row in analysis["fit"]],
            "",
            "Refinement stages",
            *[f"   {row['run']}: {row['niter']} iterations, ereq {row['ereq']:.0f}, "
              f"R-factor {row['r_factor']:.3g}" for row in runs.to_dict(orient="records")],
            "",
            "First shell",
            *[f"   {row['pair']}: first peak {row['position']:.2f} A, "
              f"coordination {row['coordination']:.2f}" for row in analysis["shells"]],
            "",
            f"Verdict: {analysis['convergence']['verdict']}",
        ]
        fig.text(0.08, 0.80, "\n".join(lines), fontsize=9.5, va="top", family="monospace")
        pdf.savefig(fig)
        plt.close(fig)

        for name in figure_names:
            image_path = run_dir / "figures" / f"{name}.png"
            if not image_path.is_file():
                continue
            image = mpimg.imread(image_path)
            height, width = image.shape[0], image.shape[1]
            fig = plt.figure(figsize=(8.27, min(11.0, 8.27 * height / width + 0.8)))
            ax = fig.add_axes((0.04, 0.04, 0.92, 0.86))
            ax.imshow(image)
            ax.axis("off")
            fig.text(0.04, 0.95, name, fontsize=12)
            pdf.savefig(fig)
            plt.close(fig)
    return path


def run_demo(
    run_dir: Path,
    data_files: Sequence[Path],
    *,
    prompt: str | None,
    system: str | None,
    quick: bool,
    emit: Callable[[dict], None],
) -> dict:
    """Write a complete, clearly-labelled demonstration run folder."""
    emit({"type": "note", "content": BANNER})
    emit({"type": "phase", "phase": "setup"})

    datasets = []
    for path in data_files:
        title, q, data, err = _read_mint01(path)
        datasets.append((title, q, data, err))
        emit({"type": "progress",
              "line": f"read {path.name}: '{title}', {q.size} points, "
                      f"Q = {q.min():.3f} to {q.max():.2f} 1/A"})
    if not datasets:
        raise ValueError("no readable data files")

    pairs = _pairs_for([title for title, *_ in datasets])
    chosen = system or "unidentified (demo)"
    emit({"type": "progress", "line": f"model system: {chosen}; atom pairs {', '.join(pairs)}"})

    steps = 30 if quick else 300
    stages = [("scratch start", 0.0), ("equilibration", 0.0), ("refinement", 3.0)]
    for index, (stage, ereq) in enumerate(stages, start=1):
        emit({"type": "progress", "line": f"run {index}: {stage}, ereq {ereq:.0f}"})
        for iteration in range(1, 4):
            emit({"type": "progress",
                  "line": f"   iteration {iteration * steps // 3 // len(stages)} of "
                          f"{steps // len(stages)}  (simulated)"})
            time.sleep(0.25)

    _write_tables(run_dir, datasets, pairs, steps)
    emit({"type": "progress", "line": "wrote tables/{fits,partials,potentials,history,runs}.csv"})
    figure_names = _write_figures(run_dir)
    emit({"type": "progress", "line": f"wrote figures: {', '.join(figure_names)}"})
    analysis = _write_analysis(run_dir, datasets, pairs, chosen)

    spec = {"demo": True, "system": chosen, "note": BANNER,
            "datasets": [{"file": p.name, "title": t} for p, (t, *_) in zip(data_files, datasets)],
            "user_prompt": prompt}
    (run_dir / "spec.json").write_text(json.dumps(spec, indent=2), encoding="utf-8")

    discussion = f"""> **{BANNER}**

The uploaded data was read and plotted, and a run folder with the structure of a
real Mellor run was written, but no Monte Carlo simulation took place: the EPSR
engine (`epsr`, a Rust extension) and the agent (`mellor`) are not installed in
this environment.

What is real on this page: the measured F(Q) curves, their Q ranges and point
counts, and the warnings raised while validating the files. What is synthetic:
the simulated F(Q), every partial g(r), both potentials, the whole convergence
history, the R-factors and the first-shell table.

To get real results, install the agent next to this app and restart the server:

```bash
cd ../mellor && uv sync        # builds the Rust engine; a few minutes
```

With `mellor` importable, this page will show a refinement of
{', '.join(t for t, *_ in datasets)}, run by the engine that was verified byte for
byte against the original Fortran EPSR.
"""
    (run_dir / "discussion.md").write_text(discussion, encoding="utf-8")

    manifest = {"demo": True, "banner": BANNER, "runs": analysis["runs"],
                "data_files": [p.name for p in data_files], "quick": quick,
                "versions": {"mellor": None, "epsr": None, "epsrweb": "0.1.0"}}
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    script = run_dir / f"run_{chosen.split()[0].lower()}.py"
    script.write_text(
        '"""Placeholder. A real Mellor run writes a standalone script here that\n'
        "reproduces the refinement with only `epsr` and `epsrvis` installed.\n"
        f'\n{BANNER}\n"""\n\nraise SystemExit("demonstration run: nothing to reproduce")\n',
        encoding="utf-8",
    )

    emit({"type": "phase", "phase": "report"})
    report = _write_report(run_dir, analysis, figure_names, chosen)
    emit({"type": "progress", "line": f"wrote {report.name}"})
    return {
        "report": str(report),
        "spec_source": "demo",
        "message": BANNER,
        "demo": True,
    }
