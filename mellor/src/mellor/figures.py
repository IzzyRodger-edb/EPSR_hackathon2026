"""Figures of a run, made with the repository's visualisation toolkit (``epsrvis``).

``make_figures`` writes into ``<run>/figures``:

* ``00_data.png``: the measured F(Q) of every input file (``mellor.data.plot_data``);
* ``01_convergence`` ... ``09_displacement`` (whichever apply), ``summary.csv`` and
  ``summary.md`` from ``epsrvis.make_report`` on the recorded trajectory;
* ``10_final_outputs.png`` from ``epsrvis.plots.final_outputs`` on the EPSR folder
  (works without a trajectory);
* ``11_fit_final.png``: data vs simulation per data set at the last iteration, with
  the residual below (from the trajectory arrays; ``epsrvis.plots.fq_fit`` on the
  last step).

The system settings of ``epsrvis.systems`` are used when the case name is known
(``water``, ``silica_bulk``); otherwise the toolkit's generic defaults apply.

Figures are drawn with matplotlib's Agg backend (no display needed). A figure that
cannot be made is skipped and the reason is logged (``logging``, logger
``mellor.figures``); only when the trajectory exists and not one of its figures could
be made is an error raised.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib

from . import data as data_module
from .runner import RunResult
from .script import epsrvis_case

matplotlib.use("Agg")

__all__ = ["make_figures", "FIGURE_CAPTIONS"]

log = logging.getLogger(__name__)

#: figure name -> caption used in the report.
FIGURE_CAPTIONS: dict[str, str] = {
    "00_data": "Measured structure factors F(Q) of the input data files (offset vertically).",
    "01_convergence": "Convergence: energies, pressure, R-factor per data set, R², EP amplitude against the requested ereq, density and Monte Carlo acceptance, iteration by iteration.",
    "02_fq_fit": "Measured (points) and simulated (lines) F(Q) and their difference at the end of each EPSR run.",
    "03_fq_evolution": "Simulated F(Q) through the refinement (light to dark) against the data.",
    "04_residual_map": "Data minus simulation over Q (horizontal) and iteration (vertical).",
    "05_partials": "Partial g(r) (top) and intermolecular partial S(Q) (bottom) per atom pair.",
    "06_potentials": "Empirical potential per pair through the refinement and its effect on g(r).",
    "07_first_shell": "First-peak position, height and coordination number of each pair against iteration.",
    "08_local_structure": "Coordination-number and bond-angle distributions from the atom coordinates.",
    "09_displacement": "Root-mean-square atom displacement since the first iteration.",
    "10_final_outputs": "Overview from the final EPSR output files: F(Q) fits, partial g(r), empirical potential and the energy history.",
    "11_fit_final": "Final fit of every data set with the residual.",
}


def make_figures(
    result: RunResult, out_dir: str | Path | None = None, *, dpi: int = 150
) -> dict[str, Path]:
    """Write the figures and return ``{name: path}`` (names as in :data:`FIGURE_CAPTIONS`,
    plus ``summary.csv``/``summary.md``). Figures that cannot be made (no trajectory,
    no structural analysis for the system) are skipped, not errors."""
    import epsrvis
    import matplotlib.pyplot as plt
    from epsrvis import plots
    from epsrvis.report import FIGURES

    out = Path(out_dir) if out_dir is not None else result.folder / "figures"
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    problems: dict[str, str] = {}

    # 00: the measured data
    try:
        path = _data_figure(result, out / "00_data.png")
    except Exception as error:  # a bad data file must not cost the other figures
        problems["00_data"] = _describe(error)
    else:
        if path is not None:
            written["00_data"] = path

    # 01-09 and the summary tables, then 11: from the recorded trajectory
    traj = _trajectory(result)
    if traj is not None:
        for name in FIGURES:
            # one call per figure, so that a figure which fails for an unusual
            # system does not take the other eight (and the summary tables) with it
            try:
                written.update(
                    epsrvis.make_report(
                        traj, out, figures=[name], dpi=dpi, verbose=False
                    )
                )
            except Exception as error:
                problems[name] = _describe(error)
        try:
            figure = plots.fq_fit(traj, steps=[traj.nsteps - 1])
            written["11_fit_final"] = _save(figure, out / "11_fit_final.png", dpi)
        except Exception as error:
            problems["11_fit_final"] = _describe(error)
        if not any(name in written for name in FIGURES):
            plt.close("all")
            raise RuntimeError(
                "no figure could be made from the trajectory: "
                + "; ".join(f"{name}: {why}" for name, why in problems.items())
            )

    # 10: from the EPSR output files (no trajectory needed)
    try:
        figure = plots.final_outputs(result.epsr_folder)
    except FileNotFoundError as error:
        log.info("10_final_outputs skipped: %s", error)
    except Exception as error:
        problems["10_final_outputs"] = _describe(error)
    else:
        written["10_final_outputs"] = _save(figure, out / "10_final_outputs.png", dpi)

    for name, why in problems.items():
        log.warning("figure %s could not be made: %s", name, why)
    plt.close("all")
    return written


# --------------------------------------------------------------------------- helpers
def _describe(error: BaseException) -> str:
    return f"{type(error).__name__}: {error}"


def _save(figure: Any, path: Path, dpi: int) -> Path:
    import matplotlib.pyplot as plt

    figure.savefig(path, dpi=dpi)
    plt.close(figure)
    return path


def _trajectory(result: RunResult) -> Any | None:
    """The recorded trajectory with its ``epsrvis`` case name set, or None."""
    path = result.trajectory_path
    if path is None or not (path / "meta.json").is_file():
        return None
    from epsrvis import Trajectory

    try:
        traj = Trajectory.load(path)
    except Exception as error:
        log.warning(
            "the trajectory at %s could not be read: %s", path, _describe(error)
        )
        return None
    if traj.nsteps == 0:
        log.warning("the trajectory at %s holds no iterations", path)
        return None
    if not traj.case:
        # the script stores the case name in the trajectory and the manifest; a
        # trajectory recorded elsewhere gets the name the script would have used
        traj.meta["case"] = result.manifest.get("epsrvis_case") or epsrvis_case(
            result.spec
        )
    return traj


def _data_files(result: RunResult) -> list[Path]:
    """The measured data files: the copies in ``<run>/data``, else the originals."""
    originals = {
        dataset.file_name: Path(dataset.file) for dataset in result.spec.datasets
    }
    names = [
        str(item["file"])
        for item in result.manifest.get("datasets") or []
        if item.get("file")
    ] or list(originals)
    files: list[Path] = []
    for name in names:
        copy = result.data_folder / name
        if copy.is_file():
            files.append(copy)
        elif name in originals and originals[name].is_file():
            files.append(originals[name])
    return files


def _data_figure(result: RunResult, path: Path) -> Path | None:
    files = _data_files(result)
    if not files:
        log.warning(
            "00_data skipped: no data files found in %s nor at the specification's paths",
            result.data_folder,
        )
        return None
    datasets = [data_module.read_mint01(file) for file in files]
    return data_module.plot_data(datasets, path)
