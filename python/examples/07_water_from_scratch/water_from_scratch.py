"""Refine liquid water against neutron data, with every input defined in Python.

Run from the ``python/`` folder:

    uv run python examples/07_water_from_scratch/water_from_scratch.py            # ~4 minutes
    uv run python examples/07_water_from_scratch/water_from_scratch.py --quick    # ~30 s
    uv run python examples/07_water_from_scratch/water_from_scratch.py --no-3d    # no 3D image

Unlike ``05_experimental_water.py``, which copies the setup files from
``validation/water/``, this script starts from an empty folder holding only the
three measured data files (SANDALS: H2O, D2O, HDO 50:50, ``validation/water/*.mint01``).
Everything else is defined below and written by the script:

* the molecule: geometry, intramolecular restraints and the reference potential
  (Lennard-Jones epsilon and sigma, mass and charge of each atom type), written as
  the ``Water.mol`` text that ``makemole`` reads;
* the box (``mixato``, ``randomise``, ``fmole``);
* the neutron weights of each isotopic sample (``.NWTS.dat`` files from Python dicts);
* a minimal EPSR control file (``.EPSR.inp``): file names and per-dataset settings
  only. EPSR fills in every other keyword with its default.

The schedule is that of ``validation/water/steps.txt``: a scratch start, 99
iterations on the reference potential, then two runs of 100 iterations with the
empirical potential (EReq 3 kJ/mol). Everything is written to ``output/`` next to
this script: the EPSR run folder (``output/epsr_run/``), the recorded trajectory, the
figures of the visualisation toolkit (``output/figures/``) and ``output/summary.txt``.
"""

import math
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

EXAMPLE_FOLDER = Path(__file__).resolve().parent
REPOSITORY = EXAMPLE_FOLDER.parents[2]
MEASURED_DATA = REPOSITORY / "validation" / "water"  # only the *.mint01 files are used
# Experimental partial g(r) of Soper (2000), used only to compare with the result.
REFERENCE_PARTIALS = REPOSITORY / "experimental_data" / "water" / "reference_partials"
OUTPUT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLE_FOLDER)) / "output"
RUN_FOLDER = OUTPUT / "epsr_run"
FIGURE_FOLDER = OUTPUT / "figures"

QUICK = "--quick" in sys.argv
DRAW_3D = "--no-3d" not in sys.argv
EQUILIBRATION_ITERATIONS = 10 if QUICK else 99
REFINEMENT_ITERATIONS = 10 if QUICK else 100
REQUESTED_EP_AMPLITUDE = 3.0  # kJ/mol (EReq)

# --- The system ---------------------------------------------------------------------

#: Reference potential of each atom type: Lennard-Jones epsilon (kJ/mol) and sigma (Å),
#: mass (amu), charge (e) and element. EPSR has no built-in values: every type needs
#: them. Pairs are mixed by Lorentz-Berthelot: eps_ij = sqrt(eps_i eps_j),
#: sigma_ij = (sigma_i + sigma_j)/2, plus Coulomb q_i q_j / r.
ATOM_TYPES = {
    "OW": {
        "epsilon": 0.6503,
        "sigma": 3.165492,
        "mass": 15.9994,
        "charge": -0.82,
        "element": "O",
    },
    "HW": {"epsilon": 0.0, "sigma": 0.0, "mass": 2.0, "charge": 0.41, "element": "H"},
}
#: Sphere radius for each type's entry on the potential line (only used by EPSR's
#: coarse-graining options; given so that the .ato matches the validated one).
RADMAXSPH = {"OW": 1.55, "HW": 0.0}

OH_BOND = 0.976  # Å
HOH_ANGLE = 107.134  # degrees
#: Atoms: label and position (Å). The molecule lies in the xy plane.
ATOMS = [
    ("OW", (0.0, 0.0, 0.0)),
    ("HW", (OH_BOND, 0.0, 0.0)),
    (
        "HW",
        (
            OH_BOND * math.cos(math.radians(HOH_ANGLE)),
            OH_BOND * math.sin(math.radians(HOH_ANGLE)),
            0.0,
        ),
    ),
]
BONDS = [(1, 2), (1, 3)]  # 1-based atom numbers
ANGLES = [(2, 1, 3)]
#: Restraint stiffness (EPSR "temperatures") and the temperature of the simulation.
SETTINGS = {
    "temperature": 300,
    "vibtemp": 65,
    "angtemp": 1,
    "dihtemp": 10,
    "ecoredcore": "0 0",
}

N_MOLECULES = 1000
ATOMIC_DENSITY = 0.1  # atoms/Å³

#: Measured datasets and the deuterium content of the hydrogen in each sample:
#: ``abundances`` lists mass numbers and fractions (0 = natural isotope mix).
DATASETS = [
    ("SLS18502-D2O.mint01", "D2O", "2 1"),
    ("SLS18500-HDO5050.mint01", "HDO", "1 0.5 2 0.5"),
    ("SLS18498-H2O.mint01", "H2O", "0 1"),
]
DATA_TYPE = 5  # nrtype: the .mint01 files are interference differential cross sections
R_MIN_BACKGROUND = (
    0.7  # rshmin (Å): below this, g(r) structure is treated as background
)


def mol_file_text() -> str:
    """The ``.mol`` text that ``makemole`` turns into the molecule template ``.ato``."""

    def number(value: float) -> str:
        return repr(float(np.float32(value)))  # the single precision EPSR reads

    neighbours = {atom: [] for atom in range(1, len(ATOMS) + 1)}
    for first, second in BONDS:
        neighbours[first].append(second)
        neighbours[second].append(first)
    lines = [" .gmol 0"]
    for atom, (_, position) in enumerate(ATOMS, start=1):
        coordinates = " ".join(number(x) for x in position)
        bonded = " ".join(str(n) for n in neighbours[atom])
        lines.append(
            f" atom {atom} {atom:03d} {coordinates} {len(neighbours[atom])} {bonded}"
        )
    for first, second in BONDS:
        lines.append(f"bond {first:03d} {second:03d} {OH_BOND}")
    for first, middle, last in ANGLES:
        lines.append(f"angle {first:03d} {middle:03d} {last:03d} {HOH_ANGLE}")
    for atom, (label, _) in enumerate(ATOMS, start=1):
        lines.append(f"changelabel {atom:03d} {label}")
    lines.append("pseudofac 1.0")
    for label, p in ATOM_TYPES.items():
        lines.append(
            f"potential {label} {p['epsilon']} {p['sigma']} {p['mass']} {p['charge']}"
            f" {p['element']} 1 0 {RADMAXSPH[label]} 1"
        )
    lines += [f"{keyword} {value}" for keyword, value in SETTINGS.items()]
    lines.append("density 0.001")  # of the one-molecule template; the box sets its own
    return "\n".join(lines) + "\n"


def epsr_input_text(box: str) -> str:
    """A minimal ``.EPSR.inp``: EPSR supplies the defaults of every other keyword."""
    lines = [f"{box}.EPSR", f"fnameato {box}.ato", f"fnamepcof {box}.pcof"]
    lines.append(f"ndata {len(DATASETS)}")
    for index, (data_file, sample, _) in enumerate(DATASETS, start=1):
        lines += [
            f"data {index}",
            f"datafile {data_file}",
            f"wtsfile {sample}.NWTStot.wts",
            f"nrtype {DATA_TYPE}",
            f"rshmin {R_MIN_BACKGROUND}",
        ]
    lines.append("q")
    return "\n".join(lines) + "\n"


summary_lines: list[str] = []


def report(text: str = "") -> None:
    """Print a line and keep it for ``output/summary.txt``."""
    print(text)
    summary_lines.append(text)


# --- 1. An empty folder with only the measured data ----------------------------------
shutil.rmtree(OUTPUT, ignore_errors=True)
RUN_FOLDER.mkdir(parents=True)
for data_file, _, _ in DATASETS:
    shutil.copy2(MEASURED_DATA / data_file, RUN_FOLDER / data_file)

# --- 2. Molecule template and box -----------------------------------------------------
(RUN_FOLDER / "Water.mol").write_text(mol_file_text())
tools.run(RUN_FOLDER, "makemole", ["Water"])  # Water.mol -> Water.ato
# mixato asks: number of molecule types, the template, how many, atomic density, box name
tools.run(
    RUN_FOLDER,
    "mixato",
    answers=["1", "Water", str(N_MOLECULES), str(ATOMIC_DENSITY), "waterbox"],
)
tools.run(RUN_FOLDER, "randomise", ["waterbox"])  # random positions and orientations
(
    RUN_FOLDER / "Water.mol"
).unlink()  # the Rust fmole refuses to run while the .mol is present
tools.run(RUN_FOLDER, "fmole", ["waterbox", "100"])  # relax the molecule geometry

# --- 3. Neutron weights for each isotopic sample --------------------------------------
for _, sample, hydrogen_abundances in DATASETS:
    tools.write_keyword_file(
        RUN_FOLDER / f"{sample}.NWTS.dat",
        "nwts",
        {"fnameato": "waterbox.ato", "nsamples": 1, "ncomponent": len(ATOM_TYPES)},
        [
            {"atom": "OW", "atsymbol": "O", "iexchange": 0, "abundances": "0 1"},
            # iexchange 1: H exchanges between molecules (labile hydrogens)
            {
                "atom": "HW",
                "atsymbol": "H",
                "iexchange": 1,
                "abundances": hydrogen_abundances,
            },
        ],
    )
    tools.nwts(RUN_FOLDER, f"{sample}.NWTS.dat")  # -> <sample>.NWTStot.wts

# --- 4. EPSR control file --------------------------------------------------------------
(RUN_FOLDER / "waterbox.EPSR.inp").write_text(epsr_input_text("waterbox"))
project = epsr.Project(RUN_FOLDER, "waterbox.EPSR.inp")

# Record every EPSR iteration for the figures (step 7) when epsr.vis is available.
recorder = vis.TrajectoryRecorder(case="water") if vis is not None else None


def run(simulation: epsr.Simulation) -> epsr.Summary:
    """One EPSR run, recorded iteration by iteration when epsr.vis is available."""
    if recorder is None:
        return simulation.run()
    return vis.run_recorded(simulation, recorder)


# --- 5. EPSR: scratch start, equilibrate on the reference potential, refine -------------
first = project.simulation()  # no .pcof yet: EPSR creates it with default settings
report("Reference potential read back from the box (per atom type):")
for t in first.potential.type_params:
    report(
        f"  {t.label}: epsilon {t.epsilon:.4f} kJ/mol, sigma {t.sigma:.4f} Å,"
        f" charge {t.charge:+.3f} e, mass {t.mass:.4f} amu"
    )
report("Mixed pair parameters (Lorentz-Berthelot) and Coulomb prefactor:")
labels = list(ATOM_TYPES)
for i, a in enumerate(labels):
    for b in labels[i:]:
        pa, pb = ATOM_TYPES[a], ATOM_TYPES[b]
        report(
            f"  {a}-{b}: eps_ij {math.sqrt(pa['epsilon'] * pb['epsilon']):.4f} kJ/mol,"
            f" sigma_ij {(pa['sigma'] + pb['sigma']) / 2:.4f} Å,"
            f" 1390·q_i·q_j {1390 * pa['charge'] * pb['charge']:+.2f} kJ mol⁻¹ Å"
        )
run(first)
equilibrated = run(project.simulation(overrides={"niter": EQUILIBRATION_ITERATIONS}))
refinement = {
    "n_runs": 2,
    "overrides": {"niter": REFINEMENT_ITERATIONS},
    "initial_overrides": {"ereq": REQUESTED_EP_AMPLITUDE},
}
if recorder is None:
    history = epsr.refine(project, progress=False, **refinement)
else:
    history = vis.refine_recorded(project, recorder=recorder, **refinement)

# --- 6. Results -------------------------------------------------------------------------
simulation = project.simulation()  # reads the accumulated results back from the folder
report()
report(f"R-factor: {equilibrated.r_factor:.3e} on the reference potential")
report(
    f"          {history.r_factor[-1]:.3e} after refinement with the empirical potential"
)
for dataset_index, dataset in enumerate(simulation.datasets):
    fit = simulation.fit(dataset_index)
    report(
        f"  {dataset.name:24s} fitted over Q = {fit.q[0]:.2f} to {fit.q[-1]:.2f} 1/Å"
    )
r = simulation.r
for a, b, low, high, meaning in (
    ("OW", "OW", 2.0, 4.0, "nearest-neighbour waters"),
    ("OW", "HW", 1.5, 2.3, "hydrogen bond"),
):
    gr = simulation.gr_inter(a, b)
    window = (r > low) & (r < high)
    report(f"{a}-{b} first peak at {r[window][gr[window].argmax()]:.2f} Å ({meaning})")


def write_summary() -> None:
    (OUTPUT / "summary.txt").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )


if vis is None:
    report("(install the 'plot' extra for the figures)")
    write_summary()
    sys.exit()

# --- 7. Figures -------------------------------------------------------------------------
import matplotlib  # noqa: E402  (available: epsr.vis imported)

matplotlib.use("Agg")  # write files, no display needed
import matplotlib.pyplot as plt  # noqa: E402

trajectory = recorder.trajectory()
recorder.save(OUTPUT / "trajectory")  # reload later: vis.Trajectory.load(...)
report(f"\n{trajectory.nsteps} EPSR iterations recorded; figures in {FIGURE_FOLDER}")
vis.make_report(trajectory, FIGURE_FOLDER, verbose=False)  # 01-09 + summary tables
fit_quality = vis.metrics.fit_metrics(trajectory)
for dataset_index, name in enumerate(trajectory.datasets):
    report(
        f"  {name:24s} R² of the fit: {fit_quality['r2'][0, dataset_index]:7.3f}"
        f" (first iteration) -> {fit_quality['r2'][-1, dataset_index]:.3f} (last)"
    )

axes = epsr.plot.fits(simulation)  # final fit to each dataset
axes.figure.savefig(FIGURE_FOLDER / "fits.png", dpi=150)
plt.close(axes.figure)

# The reference potential defined in ATOM_TYPES, as EPSR uses it, and the final
# empirical potential on top of it.
figure, (left, right) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
epsr.plot.potentials(simulation, which="reference", ax=left)
left.set_xlim(0.0, 12.5)  # the reference potential is truncated at rmaxpt = 12 Å
left.set_ylim(-400.0, 600.0)  # the Coulomb terms are hundreds of kJ/mol at 2-3 Å
left.set_title("reference potential (LJ + Coulomb)")
epsr.plot.potentials(simulation, which="empirical", ax=right)
right.set_title("empirical potential after refinement")
figure.savefig(FIGURE_FOLDER / "12_reference_and_empirical_potential.png", dpi=150)
plt.close(figure)

# Partial g(r) against the experimental partials of Soper (2000), 298 K.
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
oxygen, hydrogen, _ = vis.metrics.neighbour_pairs(
    trajectory.positions[-1],
    trajectory.box(-1),
    trajectory.atom_labels,
    "OW",
    "HW",
    2.45,
)
molecule = trajectory.molecule_of_atom
hydrogen_bonds = (molecule[oxygen] != molecule[hydrogen]).sum()
report(
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
        report("  3D configuration + partial g(r): 11_configuration.png")
    except ImportError:
        report("  (install the 'viewer' extra for the 3D configuration image)")

write_summary()
