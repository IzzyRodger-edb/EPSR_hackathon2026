"""EPSR refinements against real neutron data, run through the Python wrapper.

Two systems from ``experimental_data/`` (the Dissolve examples), set up as in
``validation/`` and ``validation/PYTHON_WRAPPER_WORKFLOW.md``:

* **water**: 1000 H2O molecules at 0.1 atoms/Å³, fitted to the SANDALS H2O, D2O and
  HDO (50:50) data;
* **silica**: 1500 Si + 3000 O at 0.0664 atoms/Å³, fitted to the NIMROD SiO2 data.

Each case is built from scratch with the setup programs (molecule templates, box,
randomisation, geometry relaxation, neutron weights), then run on the reference
potential and refined with the empirical potential.

The default tests use a shortened schedule (10 instead of 99/100 iterations per
run) and check that

* the setup programs build the expected boxes and neutron weights,
* the experimental data are read and fitted,
* the structure is physically sensible and the empirical potential improves the fit,
* the wrapper writes exactly the same files as the ``epsr-rs`` CLI after every step.

``test_full_schedule_matches_fortran`` runs the complete schedule of
``validation/<case>/steps.txt`` through the wrapper and through the original Fortran
programs and compares every file after every EPSR run. It needs a Fortran build
(``rust/tools/build_fortran_reference.sh BUILD_DIR``) and is enabled with
``EPSR_FORTRAN_BUILD=BUILD_DIR``; it takes several minutes per case.
"""

from __future__ import annotations

import os
import warnings
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pytest

import epsr

from epsr_cases import cli_binary, compare_dirs
from epsr_workflow import (
    MIXATO_CASES,
    SILICA,
    WATER,
    ExperimentalCase,
    Step,
    WrapperRun,
    WrapperStepResult,
    cli_command,
    coordination_number,
    first_peak_position,
    fortran_command,
    load_steps,
    prepare_case_folder,
    read_weights,
    run_step_externally,
    run_steps_with_wrapper,
    shorten_schedule,
)

#: Iterations per EPSR run in the shortened schedule (the full one uses 99 and 100).
SHORT_EQUILIBRATION_ITERATIONS = 10
SHORT_REFINEMENT_ITERATIONS = 10

pytestmark = [pytest.mark.experimental, pytest.mark.slow]


def _case_parameter(case: ExperimentalCase) -> object:
    marks = [pytest.mark.xdist_group(f"experimental-{case.name}")]
    if not case.available:
        marks.append(
            pytest.mark.skip(reason=f"{case.case_folder} or its data not found")
        )
    return pytest.param(case, id=case.name, marks=marks)


CASE_PARAMETERS = [_case_parameter(case) for case in MIXATO_CASES]


@dataclass
class WorkflowRun(WrapperRun):
    """The shortened workflow of one case, run through the wrapper and the CLI side by side."""

    differences_from_cli: list[str] = field(default_factory=list)
    cli_available: bool = True


@dataclass
class SideBySideRun:
    """Outcome of running the steps through the wrapper and a reference program."""

    epsr_results: list[WrapperStepResult] = field(default_factory=list)
    differences: list[str] = field(default_factory=list)  # one entry per bad comparison
    close_files: list[str] = field(
        default_factory=list
    )  # agree in all but the last digits


def _run_workflow_against(
    case: ExperimentalCase,
    steps: list[Step],
    wrapper_folder: Path,
    reference_folder: Path | None,
    reference_command: Callable[[Step], list[str]] | None,
    require_identical: bool,
) -> SideBySideRun:
    """Run the steps through the wrapper and, optionally, a reference program side by side.

    After every EPSR run both folders are compared, requiring identical files or
    agreement within the last-digit tolerance (see ``docs/testing.md``).
    """
    prepare_case_folder(case, wrapper_folder)
    if reference_folder is not None:
        prepare_case_folder(case, reference_folder)

    outcome = SideBySideRun()
    for result in run_steps_with_wrapper(steps, wrapper_folder, case.case_folder):
        if reference_folder is not None and reference_command is not None:
            run_step_externally(
                result.step, reference_folder, case.case_folder, reference_command
            )
        if not result.step.is_epsr_run:
            continue  # set steps are applied to the files only by the reference side
        outcome.epsr_results.append(result)
        if reference_folder is None:
            continue

        comparison = compare_dirs(reference_folder, wrapper_folder)
        agrees = comparison.all_identical if require_identical else comparison.passed
        if not agrees:
            outcome.differences.append(f"after '{result.step}':\n{comparison.report()}")
        outcome.close_files += [
            f"{result.step}: {item.name}" for item in comparison.with_status("close")
        ]
    return outcome


@pytest.fixture(scope="module", params=CASE_PARAMETERS)
def workflow(
    request: pytest.FixtureRequest, tmp_path_factory: pytest.TempPathFactory
) -> WorkflowRun:
    """The shortened workflow of a case, through the wrapper and, side by side, the CLI."""
    case: ExperimentalCase = request.param
    steps = shorten_schedule(
        load_steps(case), SHORT_EQUILIBRATION_ITERATIONS, SHORT_REFINEMENT_ITERATIONS
    )
    work = tmp_path_factory.mktemp(case.name)
    cli_executable = cli_binary()
    outcome = _run_workflow_against(
        case,
        steps,
        wrapper_folder=work / "wrapper",
        reference_folder=work / "cli" if cli_executable else None,
        reference_command=cli_command(cli_executable) if cli_executable else None,
        require_identical=True,
    )
    return WorkflowRun(
        case=case,
        folder=work / "wrapper",
        epsr_results=outcome.epsr_results,
        differences_from_cli=outcome.differences,
        cli_available=cli_executable is not None,
    )


# --------------------------------------------------------------------- equivalence
def test_wrapper_matches_cli_after_every_step(workflow: WorkflowRun) -> None:
    """Setup programs, overrides and EPSR runs give the same files as the epsr-rs CLI.

    The CLI side applies the ``set`` steps by editing the files; the wrapper passes
    them as Simulation overrides. Identical files show both routes are equivalent.
    """
    if not workflow.cli_available:
        pytest.skip("epsr-rs CLI not available")
    assert not workflow.differences_from_cli, "\n".join(workflow.differences_from_cli)


# --------------------------------------------------------------------- setup
def test_box(workflow: WorkflowRun) -> None:
    """mixato builds a cubic box of the requested density (values from the workflow notes)."""
    simulation = workflow.final
    assert simulation.box_size == pytest.approx(workflow.case.box_side, abs=1e-3)
    assert simulation.n_atoms == workflow.case.n_atoms
    assert simulation.n_molecules == workflow.case.n_molecules
    density_from_box = simulation.n_atoms / simulation.box_size**3
    assert density_from_box == pytest.approx(simulation.density, rel=1e-5)


def test_neutron_weights(workflow: WorkflowRun) -> None:
    """nwts writes one weights file per isotopic sample, with the expected values."""
    if workflow.case == WATER:
        hdo_weights = read_weights(workflow.folder / "HDO.NWTStot.wts")
        # (1/2 (b_H + b_D))^2 x (2/3)^2, the sanity check of PYTHON_WRAPPER_WORKFLOW.md §2.5
        assert hdo_weights[("HW", "HW")] == pytest.approx(0.02147, abs=1e-5)
        h2o_weights = read_weights(workflow.folder / "H2O.NWTStot.wts")
        d2o_weights = read_weights(workflow.folder / "D2O.NWTStot.wts")
        # hydrogen has a negative scattering length, deuterium a positive one
        assert h2o_weights[("OW", "HW")] < 0 < d2o_weights[("OW", "HW")]
        assert h2o_weights[("OW", "OW")] == pytest.approx(
            d2o_weights[("OW", "OW")], rel=1e-6
        )
    else:
        silica_weights = read_weights(workflow.folder / "SiO2.NWTStot.wts")
        assert set(silica_weights) == {("Si", "Si"), ("Si", "O"), ("O", "O")}
        assert all(weight > 0 for weight in silica_weights.values())


# --------------------------------------------------------------------- data
def test_experimental_data_are_fitted(workflow: WorkflowRun) -> None:
    simulation = workflow.final
    data_file_names = [path.name for path in workflow.case.data_files]
    assert sorted(dataset.infili for dataset in simulation.datasets) == data_file_names

    qmin_value = simulation.get_parameter("qmin")
    assert qmin_value is not None
    q_minimum = float(qmin_value.split()[0])
    for dataset_index in range(simulation.n_datasets):
        fit = simulation.fit(dataset_index)
        assert fit.q[0] >= q_minimum - 1e-6
        assert fit.q[-1] > 25.0  # the data extend to high Q
        assert np.all(np.isfinite(fit.data)) and np.any(fit.data != 0)
        assert np.all(np.isfinite(fit.simulated))
        data_spread = float(np.std(fit.data))
        root_mean_square_residual = float(np.sqrt(np.mean(fit.residual**2)))
        assert root_mean_square_residual < data_spread  # the model captures the data


def test_empirical_potential_improves_the_fit(workflow: WorkflowRun) -> None:
    """R-factors drop once the empirical potential is switched on, within the ereq limit."""
    end_of_equilibration = workflow.summary_after(
        1
    )  # reference potential only (ereq 0)
    after_refinement = workflow.summary_after(-1)
    assert end_of_equilibration.ep_amplitude == 0.0
    assert after_refinement.r_factor < 0.75 * end_of_equilibration.r_factor
    for before, after in zip(
        end_of_equilibration.r_factor_per_dataset,
        after_refinement.r_factor_per_dataset,
        strict=True,
    ):
        assert after < before

    requested_amplitude = workflow.case.requested_ep_amplitude
    assert 0.0 < after_refinement.ep_amplitude <= requested_amplitude * (1 + 1e-5)
    assert workflow.final.potential.params.ereq == pytest.approx(requested_amplitude)
    pcof_file = workflow.folder / workflow.case.input_file.replace(".EPSR.inp", ".pcof")
    saved_ereq = epsr.read_keyword_file(pcof_file).values["ereq"]
    assert float(saved_ereq) == pytest.approx(
        requested_amplitude
    )  # the override was saved


def test_energy_relaxes(workflow: WorkflowRun) -> None:
    """The random start has huge overlap energies; by the end the energy is negative."""
    scratch_start = workflow.summary_after(0)
    after_refinement = workflow.summary_after(-1)
    assert scratch_start.energy > 1e6
    assert after_refinement.energy < 0


# --------------------------------------------------------------------- structure
def test_water_structure(workflow: WorkflowRun) -> None:
    if workflow.case != WATER:
        pytest.skip("water only")
    simulation = workflow.final
    r = simulation.r
    # Intramolecular geometry: O-H bond 0.976 Å, H-O-H angle 107.134°, so H...H 1.570 Å.
    assert first_peak_position(
        r, simulation.gr_intra("OW", "HW"), 0.5, 1.5
    ) == pytest.approx(0.976, abs=0.04)
    assert first_peak_position(
        r, simulation.gr_intra("HW", "HW"), 1.0, 2.0
    ) == pytest.approx(1.570, abs=0.04)
    # Liquid water: first O...O neighbour shell at about 2.8 Å.
    assert first_peak_position(
        r, simulation.gr_inter("OW", "OW"), 2.0, 4.0
    ) == pytest.approx(2.8, abs=0.1)
    # g(r) of every pair tends to 1 at large r.
    for pair in simulation.pair_labels:
        assert float(np.mean(simulation.gr_inter(pair)[-50:])) == pytest.approx(
            1.0, abs=0.05
        )


def test_silica_structure(workflow: WorkflowRun) -> None:
    if workflow.case != SILICA:
        pytest.skip("silica only")
    simulation = workflow.final
    r = simulation.r
    # Amorphous silica: Si-O about 1.6 Å, O-O about 2.6 Å, Si-Si about 3.1 Å.
    assert first_peak_position(
        r, simulation.gr_inter("Si", "O"), 1.0, 2.2
    ) == pytest.approx(1.6, abs=0.15)
    assert first_peak_position(
        r, simulation.gr_inter("O", "O"), 2.0, 3.2
    ) == pytest.approx(2.6, abs=0.15)
    assert first_peak_position(
        r, simulation.gr_inter("Si", "Si"), 2.5, 3.6
    ) == pytest.approx(3.1, abs=0.2)
    # A strong Si-O bond peak, and close contacts below 1 Å negligible next to it
    # (the small Si sigma of 1.03 Å allows a few).
    silicon_oxygen_gr = simulation.gr_inter("Si", "O")
    bond_peak_height = float(np.max(silicon_oxygen_gr[(r > 1.0) & (r < 2.2)]))
    assert bond_peak_height > 5.0
    assert float(np.max(silicon_oxygen_gr[r < 1.0])) < 0.01 * bond_peak_height
    for pair in simulation.pair_labels:
        assert float(np.mean(simulation.gr_inter(pair)[-50:])) == pytest.approx(
            1.0, abs=0.05
        )


# --------------------------------------------------------------------- full schedule
def _check_converged_water(simulation: epsr.Simulation) -> None:
    """Liquid-water structure after the full refinement against the three datasets."""
    r = simulation.r
    volume = simulation.box_size**3
    oxygen_density, hydrogen_density = 1000 / volume, 2000 / volume
    oxygen_hydrogen_gr = simulation.gr_inter("OW", "HW")
    # Hydrogen bonds: an intermolecular O...H peak at about 1.8 Å ...
    hydrogen_bond_peak = first_peak_position(r, oxygen_hydrogen_gr, 1.5, 2.3)
    assert hydrogen_bond_peak == pytest.approx(1.85, abs=0.15)
    # ... with between one and two hydrogen-bonded H around each O.
    bonded_hydrogens = coordination_number(
        r, oxygen_hydrogen_gr, hydrogen_density, 2.35
    )
    assert 1.2 < bonded_hydrogens < 2.2
    # About 4.5 water neighbours in the first shell (tetrahedral network).
    oxygen_oxygen_gr = simulation.gr_inter("OW", "OW")
    first_shell = coordination_number(r, oxygen_oxygen_gr, oxygen_density, 3.4)
    assert 4.0 < first_shell < 5.2


def _check_converged_silica(simulation: epsr.Simulation) -> None:
    """Amorphous-silica network after the full refinement."""
    r = simulation.r
    volume = simulation.box_size**3
    silicon_density, oxygen_density = 1500 / volume, 3000 / volume
    silicon_oxygen_gr = simulation.gr_inter("Si", "O")
    assert first_peak_position(r, silicon_oxygen_gr, 1.0, 2.2) == pytest.approx(
        1.61, abs=0.06
    )
    # SiO4 tetrahedra joined by bridging oxygens: about 4 O around each Si
    # and 2 Si around each O.
    oxygens_per_silicon = coordination_number(r, silicon_oxygen_gr, oxygen_density, 2.2)
    silicons_per_oxygen = coordination_number(
        r, silicon_oxygen_gr, silicon_density, 2.2
    )
    assert oxygens_per_silicon == pytest.approx(4.0, abs=0.5)
    assert silicons_per_oxygen == pytest.approx(2.0, abs=0.25)


CONVERGED_STRUCTURE_CHECKS = {
    WATER.name: _check_converged_water,
    SILICA.name: _check_converged_silica,
}

FORTRAN_BUILD = os.environ.get("EPSR_FORTRAN_BUILD", "")


@pytest.mark.skipif(
    not FORTRAN_BUILD,
    reason="set EPSR_FORTRAN_BUILD to a build from rust/tools/build_fortran_reference.sh",
)
@pytest.mark.parametrize(
    "case",
    [
        pytest.param(
            case, id=case.name, marks=pytest.mark.xdist_group(f"fortran-{case.name}")
        )
        for case in MIXATO_CASES
    ],
)
def test_full_schedule_matches_fortran(case: ExperimentalCase, tmp_path: Path) -> None:
    """The complete validation schedule: wrapper against the original Fortran, every EPSR run."""
    if not case.available:
        pytest.skip(f"{case.case_folder} or its data not found")
    fortran_build = Path(FORTRAN_BUILD)
    if not (fortran_build / "epsr").exists():
        pytest.fail(f"no Fortran epsr in {fortran_build}")

    outcome = _run_workflow_against(
        case,
        load_steps(case),
        wrapper_folder=tmp_path / "wrapper",
        reference_folder=tmp_path / "fortran",
        reference_command=fortran_command(fortran_build),
        require_identical=False,
    )
    assert not outcome.differences, "\n".join(outcome.differences)
    if outcome.close_files:
        warnings.warn(
            f"{case.name}: files that agree with the Fortran only in the last digits "
            f"(see docs/testing.md):\n  " + "\n  ".join(outcome.close_files),
            stacklevel=1,
        )

    full_run = WorkflowRun(case, tmp_path / "wrapper", outcome.epsr_results)
    end_of_equilibration = full_run.summary_after(1)
    after_refinement = full_run.summary_after(-1)
    assert after_refinement.r_factor < end_of_equilibration.r_factor
    CONVERGED_STRUCTURE_CHECKS[case.name](full_run.final)
