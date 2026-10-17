"""Refine amorphous silica (SiO2) against real neutron data (NIMROD).

Run from the ``python/`` folder:

    uv run python examples/06_experimental_SiO2.py            # full schedule, ~12 minutes
    uv run python examples/06_experimental_SiO2.py --quick    # 10 iterations per run, ~1 minute
    uv run python examples/06_experimental_SiO2.py --no-3d    # skip the 3D (pyvista) image

The system (1500 Si + 3000 O atoms at 0.0664 atoms/Å³, 298 K; each atom is a
one-atom "molecule") is built from scratch with the EPSR setup programs,
equilibrated on the reference potential and then refined with the empirical
potential, following ``validation/silica_bulk/steps.txt``. That calculation was run
by the original Fortran EPSR, the Rust CLI and this wrapper, with byte-identical
output files. The setup inputs (atom templates, neutron-weight file, EPSR input
file) come from ``validation/silica_bulk/``; the measured data come from
``experimental_data/silica_bulk/data/``.

With the ``plot`` extra installed, every EPSR iteration is recorded (``epsr.vis``)
and the script ends with figures of the whole refinement in
``work/06_experimental_SiO2_figures/``: convergence, F(Q) fit, partials, empirical
potential, local structure, the formation of the tetrahedral network and (``viewer``
extra) the 3D configuration with Si coloured by its number of O neighbours.
Recording only reads results: the numbers are the same as without it.
"""

import os
import shutil
import sys
from pathlib import Path

import numpy as np

import epsr
from epsr import tools

try:
    from epsr import vis
except ImportError:  # the 'plot' extra (matplotlib, scipy) is not installed
    vis = None

EXAMPLES_FOLDER = Path(__file__).resolve().parent
REPOSITORY = EXAMPLES_FOLDER.parents[1]
SETUP_INPUTS = REPOSITORY / "validation" / "silica_bulk"
NEUTRON_DATA = REPOSITORY / "experimental_data" / "silica_bulk" / "data"
OUTPUT_ROOT = Path(
    os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLES_FOLDER.parent / "work")
)
WORK_FOLDER = OUTPUT_ROOT / "06_experimental_SiO2"
FIGURE_FOLDER = OUTPUT_ROOT / "06_experimental_SiO2_figures"

QUICK = "--quick" in sys.argv
DRAW_3D = "--no-3d" not in sys.argv
EQUILIBRATION_ITERATIONS = 10 if QUICK else 99
REFINEMENT_ITERATIONS = 10 if QUICK else 100
REQUESTED_EP_AMPLITUDE = 50.0  # kJ/mol, EReq of the Dissolve silica example
SI_O_BOND_CUTOFF = 2.0  # Å, first minimum of the Si-O g(r)

# 1. Gather the inputs in a fresh folder: setup files plus the measured data.
epsr.copy_folder(SETUP_INPUTS, WORK_FOLDER, overwrite=True)
for data_file in NEUTRON_DATA.glob("*.mint01"):
    shutil.copy2(data_file, WORK_FOLDER / data_file.name)

# 2. Build the box with the EPSR setup programs.
tools.run(WORK_FOLDER, "makemole", ["Si"])  # Si.mol -> Si.ato (a one-atom "molecule")
tools.run(WORK_FOLDER, "makemole", ["O"])
# mixato asks: number of types, then each template and count, density, box name
mixato_answers = ["2", "Si", "1500", "O", "3000", "0.0664", "silica"]
tools.run(WORK_FOLDER, "mixato", answers=mixato_answers)  # -> silica.ato
tools.run(WORK_FOLDER, "randomise", ["silica"])  # random positions
tools.nwts(WORK_FOLDER, "SiO2.NWTS.dat")  # neutron weights (natural isotopes)

# Record every EPSR iteration for the figures (step 6) when epsr.vis is available.
recorder = vis.TrajectoryRecorder(case="silica_bulk") if vis is not None else None
project = epsr.Project(WORK_FOLDER, "silica.EPSR.inp")


def run(simulation: epsr.Simulation) -> epsr.Summary:
    """One EPSR run, recorded iteration by iteration when epsr.vis is available."""
    if recorder is None:
        return simulation.run()
    return vis.run_recorded(simulation, recorder)


# 3. EPSR: scratch start (creates the .pcof), then equilibrate on the reference
#    potential only (ereq = 0).
run(project.simulation())
equilibrated = run(project.simulation(overrides={"niter": EQUILIBRATION_ITERATIONS}))

# 4. Switch on the empirical potential (ereq 50 kJ/mol) and refine twice; the second
#    run restarts from the files written by the first.
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
final = history.summaries[-1]
print(f"R-factor: {equilibrated.r_factor:.3e} on the reference potential")
print(f"          {final.r_factor:.3e} after refinement with the empirical potential")
print(
    f"Energy {final.energy:.1f} kJ/mol per atom, pressure {final.pressure:.1f} kbar, "
    f"EP amplitude {final.ep_amplitude:.1f} of {REQUESTED_EP_AMPLITUDE:.0f} kJ/mol"
)
fit = simulation.fit(0)
print(f"  {fit.name:18s} fitted over Q = {fit.q[0]:.2f} to {fit.q[-1]:.2f} 1/Å")

r = simulation.r
for first_type, second_type, low, high in (
    ("Si", "O", 1.2, 2.0),
    ("O", "O", 2.2, 3.0),
    ("Si", "Si", 2.6, 3.6),
):
    gr = simulation.gr_inter(first_type, second_type)
    window = (r > low) & (r < high)
    peak = r[window][gr[window].argmax()]
    print(f"{first_type}-{second_type} first peak at {peak:.2f} Å")

try:
    import matplotlib

    matplotlib.use("Agg")  # write the figure to a file, no display needed
    axes = epsr.plot.fits(simulation)
    figure_file = OUTPUT_ROOT / "06_experimental_SiO2_fits.png"
    axes.figure.savefig(figure_file, dpi=120)
    print("fit to the dataset written to", figure_file)
except ImportError:
    print("(install the 'plot' extra for a figure of the fit)")

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
print(
    f"  R² of the fit: {fit_quality['r2'][0, 0]:.3f} (first iteration) -> "
    f"{fit_quality['r2'][-1, 0]:.3f} (last)"
)

# Formation of the tetrahedral network: Si coordination, bridging O, order parameter.
labels = trajectory.atom_labels
sampled = np.unique(np.r_[np.arange(0, trajectory.nsteps, 5), trajectory.nsteps - 1])
fourfold_si, bridging_o, order_si = [], [], []
for step in sampled:
    positions, box = trajectory.positions[step], trajectory.box(step)
    _, fraction = vis.metrics.coordination_distribution(
        positions, box, labels, "Si", "O", SI_O_BOND_CUTOFF
    )
    fourfold_si.append(fraction[4] if len(fraction) > 4 else 0.0)
    _, fraction = vis.metrics.coordination_distribution(
        positions, box, labels, "O", "Si", SI_O_BOND_CUTOFF
    )
    bridging_o.append(fraction[2] if len(fraction) > 2 else 0.0)
    order_si.append(
        vis.metrics.tetrahedral_order(positions, box, labels, "Si", "O").mean()
    )
figure, axes = plt.subplots(1, 3, figsize=(13, 3.6))
iterations = trajectory.steps[sampled]
for axis, values, title in (
    (axes[0], fourfold_si, "Fraction of 4-coordinated Si (O within 2.0 Å)"),
    (axes[1], bridging_o, "Fraction of bridging O (2 Si within 2.0 Å)"),
    (axes[2], order_si, "Mean tetrahedral order q of Si"),
):
    axis.plot(iterations, values, marker="o", markersize=3, color=vis.style.series(0))
    for first_step, _, _ in trajectory.run_starts[1:]:
        axis.axvline(trajectory.steps[first_step] - 0.5, color="#8a8985", ls=":")
    axis.set_title(title, fontsize=9)
    axis.set_xlabel("EPSR iteration")
figure.suptitle("Silica: formation of the tetrahedral network", x=0.01, ha="left")
figure.tight_layout()
figure.savefig(FIGURE_FOLDER / "10_network.png", dpi=150)
plt.close(figure)
print(
    f"  4-coordinated Si {100 * fourfold_si[-1]:.1f} %, bridging O "
    f"{100 * bridging_o[-1]:.1f} %, tetrahedral order <q> {order_si[-1]:.3f}"
)

# Is the refinement finished?  The EP amplitude should have reached ereq.
if final.ep_amplitude < 0.95 * REQUESTED_EP_AMPLITUDE:
    print(
        f"  not converged: the EP amplitude ({final.ep_amplitude:.1f} kJ/mol) is still "
        f"below ereq ({REQUESTED_EP_AMPLITUDE:.0f}); run more iterations"
    )

if DRAW_3D:
    try:
        from epsr.vis.viewer import ConfigurationViewer

        viewer = ConfigurationViewer(
            trajectory, slab=12.0, color_by="coordination", off_screen=True
        )
        viewer.screenshot(FIGURE_FOLDER / "11_configuration.png", step=-1)
        print("  3D configuration (Si coloured by O neighbours): 11_configuration.png")
    except ImportError:
        print("  (install the 'viewer' extra for the 3D configuration image)")
