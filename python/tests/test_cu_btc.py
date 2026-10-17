"""Cu-BTC (HKUST-1) refined against real neutron data, through the Python wrapper.

The metal-organic framework Cu-BTC is the Dissolve example in
``experimental_data/cu-btc``: one cubic unit cell of 624 atoms (a = 26.3336 Å) and
the NIMROD data of the empty framework. Unlike the water and silica cases the box
is not built by ``mixato``: ``validation/cu-btc/make_pdb.py`` converts the unit cell
to a PDB file in which the framework is split into Cu2 paddlewheels and
benzene-tricarboxylate linkers, ``readpdb`` turns it into an EPSR configuration,
``changeato`` applies the Dissolve force field and ``nwts`` makes the neutron
weights. EPSR then runs with Bragg peaks up to 5 Å⁻¹, first on the reference
potential and then with the empirical potential (``validation/cu-btc/steps.txt``).

The original Fortran EPSR was run on the same steps; its folder after every EPSR
run is archived in ``validation/cu-btc/fortran_reference.tar.gz``
(``make_fortran_reference.py``). After each EPSR run through the wrapper, every
file of the folder is compared with that archive, so this module needs no Fortran
build. The full schedule (1 + 99 + 100 + 100 iterations) takes about 4 minutes.
"""

from __future__ import annotations

import shutil
import tarfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pytest

import epsr

from epsr_cases import Comparison, cli_binary, compare_dirs
from epsr_workflow import (
    CU_BTC,
    Step,
    WrapperStepResult,
    cli_command,
    coordination_number,
    first_peak_position,
    load_steps,
    prepare_case_folder,
    read_weights,
    run_step_externally,
    run_steps_with_wrapper,
)

pytestmark = [
    pytest.mark.experimental,
    pytest.mark.slow,
    pytest.mark.xdist_group("cu-btc"),
]

REFERENCE_ARCHIVE = CU_BTC.case_folder / "fortran_reference.tar.gz"
N_PADDLEWHEELS, ATOMS_PER_PADDLEWHEEL = 24, 2
N_LINKERS, ATOMS_PER_LINKER = 32, 18
CRYSTAL_COPPER_OXYGEN_BOND = 1.95  # Å, in the crystal structure
OXYGENS_PER_COPPER = 4  # square-planar coordination in the paddlewheel


@dataclass
class CuBtcRun:
    """The full Cu-BTC schedule through the wrapper, compared with the Fortran archive."""

    folder: Path  # the wrapper's folder after the last run
    epsr_results: list[WrapperStepResult] = field(default_factory=list)
    comparisons: dict[int, Comparison] = field(default_factory=dict)  # by step number

    @property
    def final(self) -> epsr.Simulation:
        simulation = self.epsr_results[-1].simulation
        assert simulation is not None
        return simulation

    def summary_after(self, run_index: int) -> epsr.Summary:
        summary = self.epsr_results[run_index].summary
        assert summary is not None
        return summary


def _extract_reference(destination: Path) -> Path:
    destination.mkdir(parents=True)
    with tarfile.open(REFERENCE_ARCHIVE) as archive:
        archive.extractall(destination, filter="data")
    return destination


@pytest.fixture(scope="module")
def cu_btc(tmp_path_factory: pytest.TempPathFactory) -> CuBtcRun:
    """Run every step of ``steps.txt`` through the wrapper; compare after each EPSR run."""
    if not CU_BTC.available or not REFERENCE_ARCHIVE.exists():
        pytest.skip(
            f"{CU_BTC.case_folder}, its data or its Fortran reference not found"
        )
    work = tmp_path_factory.mktemp("cu_btc")
    reference = _extract_reference(work / "fortran")
    pristine = prepare_case_folder(CU_BTC, work / "pristine")
    folder = prepare_case_folder(CU_BTC, work / "wrapper")

    run = CuBtcRun(folder)
    steps = load_steps(CU_BTC)
    step_results = run_steps_with_wrapper(steps, folder, CU_BTC.case_folder)
    for step_number, result in enumerate(step_results, start=1):
        if not result.step.is_epsr_run:
            continue
        run.epsr_results.append(result)
        fortran_snapshot = reference / f"after_step_{step_number:02d}"
        run.comparisons[step_number] = compare_dirs(fortran_snapshot, folder, pristine)
    return run


# --------------------------------------------------------------------- accuracy
def _wrapper_matches_cli(steps: list[Step], work: Path) -> Comparison | None:
    """Replay the schedule with the ``epsr-rs`` CLI and compare the final folders.

    Monte Carlo trajectories are chaotic, so identical final folders after 300
    iterations mean the wrapper and the CLI agreed at every step.
    """
    cli_executable = cli_binary()
    if cli_executable is None:
        return None
    wrapper_folder = prepare_case_folder(CU_BTC, work / "wrapper_again")
    cli_folder = prepare_case_folder(CU_BTC, work / "cli")
    for result in run_steps_with_wrapper(steps, wrapper_folder, CU_BTC.case_folder):
        run_step_externally(
            result.step, cli_folder, CU_BTC.case_folder, cli_command(cli_executable)
        )
    return compare_dirs(cli_folder, wrapper_folder)


def test_matches_fortran_after_every_run(cu_btc: CuBtcRun, tmp_path: Path) -> None:
    """Every file equals the original Fortran's after each of the four EPSR runs."""
    reports = {
        step_number: comparison.report()
        for step_number, comparison in cu_btc.comparisons.items()
        if not comparison.passed
    }
    close_files = [
        f"step {step_number}: {item.name}"
        for step_number, comparison in cu_btc.comparisons.items()
        for item in comparison.with_status("close")
    ]
    if close_files:
        # Not byte-identical to the Fortran: the wrapper must still equal the CLI, so
        # that the difference lies in the engine (libm rounding), not in the wrapper.
        versus_cli = _wrapper_matches_cli(load_steps(CU_BTC), tmp_path)
        if versus_cli is not None:
            assert (
                versus_cli.all_identical
            ), f"wrapper output differs from the epsr-rs CLI:\n{versus_cli.report()}"
        warnings.warn(
            "Cu-BTC: files that agree with the Fortran only in the last digits:\n  "
            + "\n  ".join(close_files),
            stacklevel=1,
        )
    assert not reports, "\n".join(
        f"after step {step_number}:\n{report}"
        for step_number, report in reports.items()
    )


# --------------------------------------------------------------------- setup
def test_framework_configuration(cu_btc: CuBtcRun) -> None:
    """readpdb splits the unit cell into 24 Cu2 paddlewheels and 32 BTC linkers."""
    simulation = cu_btc.final
    assert simulation.box_size == pytest.approx(CU_BTC.box_side, abs=1e-3)
    assert simulation.n_atoms == CU_BTC.n_atoms
    assert simulation.n_molecules == CU_BTC.n_molecules
    molecule_sizes = sorted(len(atom_range) for atom_range in simulation.molecule_atoms)
    expected_sizes = sorted(
        [ATOMS_PER_PADDLEWHEEL] * N_PADDLEWHEELS + [ATOMS_PER_LINKER] * N_LINKERS
    )
    assert molecule_sizes == expected_sizes
    assert simulation.atom_types == ["CU", "O", "C1", "C2", "C3", "H"]
    assert simulation.n_pairs == 21


def test_force_field_and_weights(cu_btc: CuBtcRun) -> None:
    """changeato applied the Dissolve Lennard-Jones parameters and charges."""
    parameters = {entry.label: entry for entry in cu_btc.final.potential.type_params}
    assert parameters["CU"].epsilon == pytest.approx(0.02092)
    assert parameters["CU"].sigma == pytest.approx(3.144)
    assert parameters["CU"].charge == pytest.approx(1.098)
    assert parameters["O"].charge == pytest.approx(-0.665)
    assert parameters["C1"].charge == pytest.approx(0.778)
    assert parameters["H"].mass == pytest.approx(1.008)
    total_charge = sum(entry.charge * entry.fraction for entry in parameters.values())
    assert total_charge == pytest.approx(0.0, abs=1e-5)  # the framework is neutral

    weights = read_weights(cu_btc.folder / "CuBTC.NWTStot.wts")
    assert len(weights) == 21
    # Natural hydrogen has a negative scattering length: pairs with one H atom get a
    # negative weight, H-H and every pair without H a positive one.
    for (first_type, second_type), weight in weights.items():
        n_hydrogens = (first_type == "H") + (second_type == "H")
        assert (weight < 0) == (n_hydrogens == 1), (first_type, second_type, weight)


def read_bragg_reflections(bragg_file: Path) -> tuple[int, list[float]]:
    """Number of reflections and the unique Bragg Q values listed in a ``.bragg`` file.

    Every line of the file starts with ``#``: a header, the line
    ``# Number of Bragg refections calculated: N``, then one line
    ``#  index  Q  multiplicity`` per unique Q value.
    """
    n_reflections, unique_q_values = 0, []
    listing_unique_values = False
    for line in bragg_file.read_text().splitlines():
        words = line.lstrip("# ").split()
        if line.startswith("# Number of Bragg refections calculated:"):
            n_reflections = int(words[-1])
        elif "Unique Bragg Q values" in line:
            listing_unique_values = True
        elif listing_unique_values and len(words) == 3:
            unique_q_values.append(float(words[1]))
    return n_reflections, unique_q_values


def test_bragg_peaks_are_calculated(cu_btc: CuBtcRun) -> None:
    """The crystal keywords make EPSR enumerate the reflections and write their intensities."""
    simulation = cu_btc.final
    n_reflections, unique_q_values = read_bragg_reflections(
        simulation.output_path("bragg")
    )
    assert (
        n_reflections > 1000
    )  # a 26 Å cell has many reflections below the Bragg Q limit
    assert len(unique_q_values) > 100
    # The list starts with the (0 0 0) reflection at Q = 0 and extends past qbraggmax,
    # because the Bragg g(r) truncation needs reflections beyond the fitted range.
    assert min(unique_q_values) == 0.0
    assert max(unique_q_values) > 5.0

    # .inten: the Bragg intensity of every atom pair on the Q grid up to qbraggmax = 5 Å⁻¹.
    intensities = np.loadtxt(simulation.output_path("inten"), comments="#", ndmin=2)
    q_step = simulation.control.qstep
    n_bragg_points = round(5.0 / q_step) + 1
    assert intensities.shape[0] == n_bragg_points
    assert intensities[0, 0] == 0.0 and intensities[-1, 0] == pytest.approx(
        5.0, abs=1e-5
    )
    assert np.max(np.abs(intensities[:, 1:])) > 0


# --------------------------------------------------------------------- refinement
def test_empirical_potential_improves_the_fit(cu_btc: CuBtcRun) -> None:
    scratch_start = cu_btc.summary_after(0)
    end_of_equilibration = cu_btc.summary_after(1)
    after_refinement = cu_btc.summary_after(-1)

    assert end_of_equilibration.ep_amplitude == 0.0
    assert end_of_equilibration.r_factor < scratch_start.r_factor
    assert after_refinement.r_factor < end_of_equilibration.r_factor
    requested = CU_BTC.requested_ep_amplitude
    assert after_refinement.ep_amplitude == pytest.approx(requested, rel=1e-5)
    assert cu_btc.final.potential.params.ereq == pytest.approx(requested)

    fit = cu_btc.final.fit(0)
    assert fit.q[0] < 0.5 and fit.q[-1] > 25.0  # the NIMROD data cover 0.25 to 50 Å⁻¹
    assert np.all(np.isfinite(fit.simulated))


# --------------------------------------------------------------------- structure
def test_framework_structure(cu_btc: CuBtcRun) -> None:
    """Molecules keep their geometry; the Cu-O coordination survives as intermolecular."""
    simulation = cu_btc.final
    r = simulation.r
    volume = simulation.box_size**3
    oxygen_density = N_LINKERS * 6 / volume

    # Bonded distances inside the molecules (readpdb bonds, held by EPSR):
    # the Cu-Cu of the paddlewheel and the linker bonds of the Dissolve species.
    assert first_peak_position(
        r, simulation.gr_intra("CU", "CU"), 2.0, 3.2
    ) == pytest.approx(2.63, abs=0.05)
    assert first_peak_position(
        r, simulation.gr_intra("C1", "O"), 1.0, 1.6
    ) == pytest.approx(1.255, abs=0.05)
    assert first_peak_position(
        r, simulation.gr_intra("C3", "H"), 0.5, 1.5
    ) == pytest.approx(0.93, abs=0.05)

    # The Cu-O coordination bond is intermolecular here (Lennard-Jones + Coulomb only),
    # so it is longer than in the crystal but each Cu still has its four oxygens.
    copper_oxygen_gr = simulation.gr_inter("CU", "O")
    bond_peak = first_peak_position(r, copper_oxygen_gr, 1.5, 2.8)
    assert (
        CRYSTAL_COPPER_OXYGEN_BOND - 0.2 < bond_peak < CRYSTAL_COPPER_OXYGEN_BOND + 0.5
    )
    oxygens_within_reach = coordination_number(r, copper_oxygen_gr, oxygen_density, 2.8)
    assert oxygens_within_reach == pytest.approx(OXYGENS_PER_COPPER, abs=0.8)

    # Neighbouring paddlewheels are 8.0 Å apart in the framework (a peak of the
    # intermolecular Cu-Cu g(r); further shells follow at 9.7 Å and beyond).
    copper_copper_gr = simulation.gr_inter("CU", "CU")
    paddlewheel_spacing = first_peak_position(r, copper_copper_gr, 7.0, 9.0)
    assert paddlewheel_spacing == pytest.approx(8.0, abs=0.3)
    assert float(np.max(copper_copper_gr[(r > 7.0) & (r < 9.0)])) > 1.0


def test_configuration_round_trip(cu_btc: CuBtcRun, tmp_path: Path) -> None:
    """The final .ato file re-read by a new Simulation gives the same configuration.

    Positions agree to the precision of the .ato file (8 significant digits).
    """
    copy = tmp_path / "copy"
    shutil.copytree(cu_btc.folder, copy)
    reread = epsr.Simulation(copy, CU_BTC.input_file)
    np.testing.assert_allclose(
        reread.positions, cu_btc.final.positions, rtol=1e-6, atol=1e-6
    )
    assert reread.molecule_atoms == cu_btc.final.molecule_atoms
    assert reread.atom_labels == cu_btc.final.atom_labels
