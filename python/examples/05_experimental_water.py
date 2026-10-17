"""Refine liquid water against real neutron data (SANDALS: H2O, D2O and HDO 50:50).

Run from the ``python/`` folder:

    uv run python examples/05_experimental_water.py            # full schedule, ~3 minutes
    uv run python examples/05_experimental_water.py --quick    # 10 iterations per run, ~20 s

The system (1000 water molecules at 0.1 atoms/Å³, 300 K) is built from scratch with
the EPSR setup programs, equilibrated on the reference potential and then refined
with the empirical potential, following ``validation/water/steps.txt``. The setup
inputs (molecule template, neutron-weight files, EPSR input file) come from
``validation/water/``; the measured data come from ``experimental_data/water/data/``.
"""

import os
import shutil
import sys
from pathlib import Path

import epsr
from epsr import tools

EXAMPLES_FOLDER = Path(__file__).resolve().parent
REPOSITORY = EXAMPLES_FOLDER.parents[1]
SETUP_INPUTS = REPOSITORY / "validation" / "water"
NEUTRON_DATA = REPOSITORY / "experimental_data" / "water" / "data"
OUTPUT_ROOT = Path(
    os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLES_FOLDER.parent / "work")
)
WORK_FOLDER = OUTPUT_ROOT / "05_experimental_water"

QUICK = "--quick" in sys.argv
EQUILIBRATION_ITERATIONS = 10 if QUICK else 99
REFINEMENT_ITERATIONS = 10 if QUICK else 100
REQUESTED_EP_AMPLITUDE = 3.0  # kJ/mol, EReq of the Dissolve water example

# 1. Gather the inputs in a fresh folder: setup files plus the measured data.
epsr.copy_folder(SETUP_INPUTS, WORK_FOLDER, overwrite=True)
for data_file in NEUTRON_DATA.glob("*.mint01"):
    shutil.copy2(data_file, WORK_FOLDER / data_file.name)

# 2. Build the box with the EPSR setup programs.
tools.run(WORK_FOLDER, "makemole", ["Water"])  # Water.mol -> Water.ato template
tools.run(WORK_FOLDER, "mixato", answers=["1", "Water", "1000", "0.1", "waterbox"])
tools.run(WORK_FOLDER, "randomise", ["waterbox"])  # random positions and orientations
(WORK_FOLDER / "Water.mol").unlink()  # fmole refuses to run while the .mol is present
tools.run(WORK_FOLDER, "fmole", ["waterbox", "100"])  # relax the molecule geometry
for sample in ("H2O", "D2O", "HDO"):
    tools.nwts(WORK_FOLDER, f"{sample}.NWTS.dat")  # neutron weights per isotopic sample

# 3. EPSR: scratch start, then equilibrate on the reference potential (ereq = 0).
project = epsr.Project(WORK_FOLDER, "waterbox.EPSR.inp")
project.simulation().run()
equilibrated = project.simulation(overrides={"niter": EQUILIBRATION_ITERATIONS}).run()

# 4. Switch on the empirical potential and refine twice.
history = epsr.refine(
    project,
    n_runs=2,
    overrides={"niter": REFINEMENT_ITERATIONS},
    initial_overrides={"ereq": REQUESTED_EP_AMPLITUDE},
    progress=False,
)

# 5. Results.
simulation = project.simulation()  # reads the accumulated results back from the folder
print(f"R-factor: {equilibrated.r_factor:.3e} on the reference potential")
print(
    f"          {history.r_factor[-1]:.3e} after refinement with the empirical potential"
)
for dataset_index, dataset in enumerate(simulation.datasets):
    fit = simulation.fit(dataset_index)
    print(f"  {dataset.name:18s} fitted over Q = {fit.q[0]:.2f} to {fit.q[-1]:.2f} 1/Å")

r = simulation.r
for first_type, second_type, low, high, meaning in (
    ("OW", "OW", 2.0, 4.0, "nearest-neighbour waters"),
    ("OW", "HW", 1.5, 2.3, "hydrogen bond"),
):
    gr = simulation.gr_inter(first_type, second_type)
    window = (r > low) & (r < high)
    peak = r[window][gr[window].argmax()]
    print(f"{first_type}-{second_type} first peak at {peak:.2f} Å ({meaning})")

try:
    import matplotlib

    matplotlib.use("Agg")  # write the figure to a file, no display needed
    axes = epsr.plot.fits(simulation)
    figure_file = OUTPUT_ROOT / "05_experimental_water_fits.png"
    axes.figure.savefig(figure_file, dpi=120)
    print("fits to the three datasets written to", figure_file)
except ImportError:
    print("(install the 'plot' extra for a figure of the fits)")
