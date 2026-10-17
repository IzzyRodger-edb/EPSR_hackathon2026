"""One EPSR run on a copy of the DMSO/water test case, then inspect the results.

Run from the ``python/`` folder:  uv run python examples/01_single_run.py
"""

from pathlib import Path

import numpy as np

import epsr

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE.parents[1] / "EPSR-test"  # 200 DMSO + 400 water, three neutron datasets
WORK = HERE.parent / "work" / "01_single_run"

# Copy the folder first: every EPSR run rewrites its input files in place.
proj = epsr.Project.from_template(TEMPLATE, WORK, overwrite=True)

# Read the inputs, overriding two run parameters (saved back to the .EPSR.inp file).
sim = proj.simulation(overrides={"ntimes": 5, "niter": 1})
print(sim)

summary = sim.run()  # releases the GIL while the Rust engine runs
print(f"R-factor {summary.r_factor:.4e}   energy {summary.energy:.3f} kJ/mol   pressure {summary.pressure:.3f} kbar")

# Grids and per-dataset fits
q, r = sim.q, sim.r
print(f"Q grid: {q.size} points up to {q[-1]:.2f} 1/Å;  r grid: {r.size} points up to {r[-1]:.2f} Å")
for i, ds in enumerate(sim.datasets):
    fit = sim.fit(i)
    rms = float(np.sqrt(np.mean(fit.residual**2)))
    print(f"  dataset {i} {ds.name:12s} Q {fit.q[0]:.2f}..{fit.q[-1]:.2f}  rms(data - fit) = {rms:.4f}")

# Partial g(r) and the empirical potential, (n_pairs, nr) float32 arrays
g = sim.gr_inter()
u = sim.potential.empirical()
print("pairs:", sim.pair_labels[:6], "...")
# Atom types: S, C, O, H belong to DMSO; O1, H1 to water.
short = r < 4.0
for a, b, what in (("O1", "H1", "water-water hydrogen bond"), ("O", "H1", "DMSO O ... water H")):
    peak = r[short][np.argmax(sim.gr_inter(a, b)[short])]
    print(f"{a}-{b} intermolecular g(r): first peak at {peak:.2f} Å ({what})")
print(f"EP range: {u.min():.3f} .. {u.max():.3f} kJ/mol")
print("positions:", sim.positions.shape, "Å, box", round(sim.box_size, 3), "Å")
