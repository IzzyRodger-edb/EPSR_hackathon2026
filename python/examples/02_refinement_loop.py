"""Repeated EPSR runs with convergence monitoring, as the EPSR shell loop does.

Run from the ``python/`` folder:  uv run python examples/02_refinement_loop.py [n_runs]
"""

import sys
from pathlib import Path

import epsr

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE.parents[1] / "EPSR-test"
WORK = HERE.parent / "work" / "02_refine"
N_RUNS = int(sys.argv[1]) if len(sys.argv) > 1 else 10

proj = epsr.Project.from_template(TEMPLATE, WORK, overwrite=True)


def report(k: int, sim: epsr.Simulation) -> None:
    s = sim.summary
    print(f"run {k:3d}  R {s.r_factor:.4e}  E {s.energy:8.3f}  P {s.pressure:6.3f}  EP {s.ep_amplitude:.2f}")


hist = epsr.refine(
    proj,
    n_runs=N_RUNS,
    overrides={"ntimes": 5},
    callback=report,
    until=lambda h: h.r_factor[-1] < 1.0e-4,  # stop early once the fit is good enough
    snapshot_every=5,
    progress=False,
)

print(f"{len(hist)} runs; R-factor {hist.r_factor[0]:.4e} -> {hist.r_factor[-1]:.4e}")
csv = hist.to_csv(WORK.parent / "02_history.csv")
print("history written to", csv)

try:
    import matplotlib

    matplotlib.use("Agg")
    axes = epsr.plot.history(hist)
    axes[0].figure.savefig(WORK.parent / "02_history.png", dpi=120)
    print("plot written to", WORK.parent / "02_history.png")
except ImportError:
    print("(install the 'plot' extra for a convergence plot)")
