"""Repeated EPSR runs with convergence monitoring, as the EPSR shell loop does.

Run from the ``python/`` folder:  uv run python examples/02_refinement_loop.py [n_runs]
"""

import os
import sys
from pathlib import Path

import epsr

EXAMPLES_FOLDER = Path(__file__).resolve().parent
TEMPLATE = EXAMPLES_FOLDER.parents[1] / "EPSR-test"
OUTPUT_ROOT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLES_FOLDER.parent / "work"))
WORK_FOLDER = OUTPUT_ROOT / "02_refine"
N_RUNS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
TARGET_R_FACTOR = 1.0e-4

project = epsr.Project.from_template(TEMPLATE, WORK_FOLDER, overwrite=True)


def report_progress(run_number: int, simulation: epsr.Simulation) -> None:
    summary = simulation.summary
    print(
        f"run {run_number:3d}  R {summary.r_factor:.4e}  E {summary.energy:8.3f}  "
        f"P {summary.pressure:6.3f}  EP {summary.ep_amplitude:.2f}"
    )


def target_reached(history: epsr.History) -> bool:
    """Stop early once the fit is good enough."""
    return history.r_factor[-1] < TARGET_R_FACTOR


history = epsr.refine(
    project,
    n_runs=N_RUNS,
    overrides={"ntimes": 5},
    callback=report_progress,
    until=target_reached,
    snapshot_every=5,
    progress=False,
)

print(f"{len(history)} runs; R-factor {history.r_factor[0]:.4e} -> {history.r_factor[-1]:.4e}")
csv_file = history.to_csv(OUTPUT_ROOT / "02_history.csv")
print("history written to", csv_file)

try:
    import matplotlib

    matplotlib.use("Agg")  # write the figure to a file, no display needed
    axes = epsr.plot.history(history)
    plot_file = OUTPUT_ROOT / "02_history.png"
    axes[0].figure.savefig(plot_file, dpi=120)
    print("plot written to", plot_file)
except ImportError:
    print("(install the 'plot' extra for a convergence plot)")
