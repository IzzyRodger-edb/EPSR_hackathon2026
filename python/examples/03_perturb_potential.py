"""Perturbing the EPSR potential from Python.

Run from the ``python/`` folder:  uv run python examples/03_perturb_potential.py

Three ways to perturb the potential:
(a) the standard EPSR control keywords, as overrides;
(b) direct, in-memory edits of the empirical-potential (EP) coefficients;
(c) a fixed EP with refinement switched off (Monte Carlo only).
"""

import os
from pathlib import Path

import numpy as np

import epsr

EXAMPLES_FOLDER = Path(__file__).resolve().parent
TEMPLATE = EXAMPLES_FOLDER.parents[1] / "EPSR-test"
OUTPUT_ROOT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLES_FOLDER.parent / "work"))
WORK_FOLDER = OUTPUT_ROOT / "03_perturb"
WATER_HYDROGEN_BOND = ("O1", "H1")  # water oxygen ... water hydrogen

project = epsr.Project.from_template(TEMPLATE, WORK_FOLDER, overwrite=True)

# (a) Standard knobs: EP amplitude limit, feedback factor, refinement factor.
simulation = project.simulation(overrides={"ereq": 15.0, "feedback": 0.8, "potfac": [1.0, 1.0]})
potential = simulation.potential
print(potential)
print("parameters:", potential.params)

# (b) Direct manipulation, in memory, before running.
original_coefficients = potential.ep_coefficients()  # (n_pairs, ncoeffp - 1)
hydrogen_bond_pair = potential.pair_index(*WATER_HYDROGEN_BOND)
print(f"O1-H1 EP minimum before: {potential.empirical(*WATER_HYDROGEN_BOND).min():.3f} kJ/mol")

damped_coefficients = original_coefficients.copy()
damped_coefficients[hydrogen_bond_pair] *= 0.5  # damp the water hydrogen-bond EP
potential.set_ep_coefficients(damped_coefficients)  # writes pcoefp and rebuilds the EP
print(
    f"O1-H1 EP minimum after damping: {potential.empirical(*WATER_HYDROGEN_BOND).min():.3f} kJ/mol"
)

potential.perturb_ep(WATER_HYDROGEN_BOND, 0.002)  # add a constant to every coefficient of one pair
potential.scale_ep(0.9)  # scale the whole EP (never beyond ereq)
sum_matches = np.allclose(
    potential.total(), potential.reference() + potential.empirical(), atol=1e-3
)
print("total = reference + empirical:", sum_matches)

summary = simulation.run()  # the modified EP is refined further and saved to the .pcof file
print(f"after run: R-factor {summary.r_factor:.4e}, EP amplitude {summary.ep_amplitude:.3f}")

# (c) Keep the EP fixed and only run Monte Carlo: potfac 0 disables refinement.
fixed_ep_simulation = project.simulation(overrides={"potfac": [0.0, 0.0]})
coefficients_before_run = fixed_ep_simulation.potential.ep_coefficients()
fixed_ep_simulation.run()
coefficients_after_run = project.simulation().potential.ep_coefficients()
print(
    "EP unchanged by a run with potfac 0:",
    np.allclose(coefficients_after_run, coefficients_before_run),
)

# Reference-potential parameters are read-only (they come from the .ato file).
for atom_type in potential.type_params:
    print(
        f"  {atom_type.label:3s} epsilon {atom_type.epsilon:.4f} kJ/mol  "
        f"sigma {atom_type.sigma:.3f} Å  charge {atom_type.charge:+.4f} e"
    )
