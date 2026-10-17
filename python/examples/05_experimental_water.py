"""Refine liquid water against real neutron data (SANDALS: H2O, D2O and HDO 50:50).

Run from the ``python/`` folder:

    uv run python examples/05_experimental_water.py            # full schedule, ~4 minutes
    uv run python examples/05_experimental_water.py --quick    # 10 iterations per run, ~30 s
    uv run python examples/05_experimental_water.py --no-3d    # skip the 3D (pyvista) image

The system (1000 water molecules at 0.1 atoms/Å³, 300 K) is built from scratch with
the EPSR setup programs, equilibrated on the reference potential and then refined
with the empirical potential, following ``validation/water/steps.txt``. The setup
inputs (molecule template, neutron-weight files, EPSR input file) come from
``validation/water/``; the measured data come from ``experimental_data/water/data/``.

With the ``plot`` extra installed, every EPSR iteration is recorded (``epsr.vis``)
and the script ends with figures of the whole refinement in
``work/05_experimental_water_figures/``: convergence, F(Q) fits, partials, empirical
potential, local structure, a comparison of the partial g(r) with the experimental
partials of Soper (2000), and (``viewer`` extra) the 3D configuration. Recording only
reads results: the numbers are the same as without it.
"""

import os
import shutil
import sys
from pathlib import Path

import epsr
from epsr import tools

try:
    from epsr import vis
except ImportError:  # the 'plot' extra (matplotlib, scipy) is not installed
    vis = None

EXAMPLES_FOLDER = Path(__file__).resolve().parent
REPOSITORY = EXAMPLES_FOLDER.parents[1]
SETUP_INPUTS = REPOSITORY / "validation" / "water"
NEUTRON_DATA = REPOSITORY / "experimental_data" / "water" / "data"
REFERENCE_PARTIALS = REPOSITORY / "experimental_data" / "water" / "reference_partials"
OUTPUT_ROOT = Path(
    os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLES_FOLDER.parent / "work")
)
WORK_FOLDER = OUTPUT_ROOT / "05_experimental_water"
FIGURE_FOLDER = OUTPUT_ROOT / "05_experimental_water_figures"

QUICK = "--quick" in sys.argv
DRAW_3D = "--no-3d" not in sys.argv
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

# Record every EPSR iteration for the figures (step 6) when epsr.vis is available.
recorder = vis.TrajectoryRecorder(case="water") if vis is not None else None
project = epsr.Project(WORK_FOLDER, "waterbox.EPSR.inp")


def run(simulation: epsr.Simulation) -> epsr.Summary:
    """One EPSR run, recorded iteration by iteration when epsr.vis is available."""
    if recorder is None:
        return simulation.run()
    return vis.run_recorded(simulation, recorder)


# 3. EPSR: scratch start, then equilibrate on the reference potential (ereq = 0).
run(project.simulation())
equilibrated = run(project.simulation(overrides={"niter": EQUILIBRATION_ITERATIONS}))

# 4. Switch on the empirical potential and refine twice.
refinement = {
    "n_runs": 2,
    "overrides": {"niter": REFINEMENT_ITERATIONS},
    "initial_overrides": {"ereq": REQUESTED_EP_AMPLITUDE},
}
if recorder is None:
    history = epsr.refine(project, progress=False, **refinement)
else:
    history = vis.refine_recorded(project, recorder=recorder, **refinement)

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

# 6. Figures of the whole refinement, iteration by iteration (epsr.vis).
if recorder is None:
    print("(install the 'plot' extra for figures of the refinement)")
    sys.exit()
import matplotlib.pyplot as plt  # noqa: E402  (available: epsr.vis imported)

trajectory = recorder.trajectory()
recorder.save(WORK_FOLDER / "trajectory")  # reload later: vis.Trajectory.load(...)
print(f"\n{trajectory.nsteps} EPSR iterations recorded; figures in {FIGURE_FOLDER}")
vis.make_report(trajectory, FIGURE_FOLDER, verbose=False)  # 01-09 + summary tables
fit_quality = vis.metrics.fit_metrics(trajectory)
for dataset_index, name in enumerate(trajectory.datasets):
    print(
        f"  {name:18s} R² of the fit: {fit_quality['r2'][0, dataset_index]:7.3f}"
        f" (first iteration) -> {fit_quality['r2'][-1, dataset_index]:.3f} (last)"
    )

# The partial g(r) against the experimental partials of Soper (2000), 298 K:
# experiment solid, EPSR in the same colour but translucent.
reference_partials = vis.load_reference_partials(
    {
        "OW-OW": REFERENCE_PARTIALS / "h2o298aoogor.txt",
        "OW-HW": REFERENCE_PARTIALS / "h2o298aohgor.txt",
        "HW-HW": REFERENCE_PARTIALS / "h2o298ahhgor.txt",
    }
)
figure = vis.plots.partials_vs_reference(
    trajectory, reference_partials, source="experiment (Soper 2000, 298 K)"
)
figure.savefig(FIGURE_FOLDER / "10_partials_vs_experiment.png", dpi=150)
plt.close(figure)

# Hydrogen bonds per molecule in the final configuration (O···H < 2.45 Å, other molecule).
labels, molecule = trajectory.atom_labels, trajectory.molecule_of_atom
oxygen, hydrogen, _ = vis.metrics.neighbour_pairs(
    trajectory.positions[-1], trajectory.box(-1), labels, "OW", "HW", 2.45
)
hydrogen_bonds = (molecule[oxygen] != molecule[hydrogen]).sum()
print(
    f"  {2 * hydrogen_bonds / trajectory.meta['nmol']:.2f} hydrogen bonds per molecule"
)

if DRAW_3D:
    try:
        from epsr.vis.viewer import ConfigurationViewer

        viewer = ConfigurationViewer(
            trajectory,
            gr="total",
            slab=10.0,
            reference=reference_partials,
            off_screen=True,
        )
        viewer.screenshot(FIGURE_FOLDER / "11_configuration.png", step=-1)
        print("  3D configuration + partial g(r): 11_configuration.png")
    except ImportError:
        print("  (install the 'viewer' extra for the 3D configuration image)")
