"""One EPSR run on a copy of the DMSO/water test case, then inspect the results.

Run from the ``python/`` folder:  uv run python examples/01_single_run.py
"""

import os
from pathlib import Path

import numpy as np

import epsr

EXAMPLES_FOLDER = Path(__file__).resolve().parent
TEMPLATE = EXAMPLES_FOLDER.parents[1] / "EPSR-test"  # 200 DMSO + 400 water, three neutron datasets
OUTPUT_ROOT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLES_FOLDER.parent / "work"))
WORK_FOLDER = OUTPUT_ROOT / "01_single_run"

# Copy the folder first: every EPSR run rewrites its input files in place.
project = epsr.Project.from_template(TEMPLATE, WORK_FOLDER, overwrite=True)

# Read the inputs, overriding two run parameters (saved back to the .EPSR.inp file).
simulation = project.simulation(overrides={"ntimes": 5, "niter": 1})
print(simulation)

summary = simulation.run()  # releases the GIL while the Rust engine runs
print(f"R-factor {summary.r_factor:.4e}")
print(f"energy   {summary.energy:.3f} kJ/mol per molecule, pressure {summary.pressure:.3f} kbar")

# Grids and per-dataset fits
q_grid, r_grid = simulation.q, simulation.r
print(f"Q grid: {q_grid.size} points up to {q_grid[-1]:.2f} 1/Å")
print(f"r grid: {r_grid.size} points up to {r_grid[-1]:.2f} Å")
for dataset_index, dataset in enumerate(simulation.datasets):
    fit = simulation.fit(dataset_index)
    root_mean_square_residual = float(np.sqrt(np.mean(fit.residual**2)))
    print(
        f"  dataset {dataset_index} {dataset.name:12s} Q {fit.q[0]:.2f}..{fit.q[-1]:.2f}  "
        f"rms(data - fit) = {root_mean_square_residual:.4f}"
    )

# Partial g(r) and the empirical potential are (n_pairs, nr) float32 arrays.
empirical_potential = simulation.potential.empirical()
print("pairs:", simulation.pair_labels[:6], "...")

# Atom types: S, C, O, H belong to DMSO; O1, H1 to water.
short_distances = r_grid < 4.0
pairs_to_show = (
    ("O1", "H1", "water-water hydrogen bond"),
    ("O", "H1", "DMSO O ... water H"),
)
for first_type, second_type, meaning in pairs_to_show:
    gr = simulation.gr_inter(first_type, second_type)
    first_peak = r_grid[short_distances][np.argmax(gr[short_distances])]
    pair_name = f"{first_type}-{second_type}"
    print(f"{pair_name} intermolecular g(r): first peak at {first_peak:.2f} Å ({meaning})")

print(f"EP range: {empirical_potential.min():.3f} .. {empirical_potential.max():.3f} kJ/mol")
print("positions:", simulation.positions.shape, "Å, box", round(simulation.box_size, 3), "Å")
