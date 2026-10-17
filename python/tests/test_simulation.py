"""Accessors, shapes, overrides, run-order rules and error mapping of ``epsr.Simulation``.

In the test case, atom types S, C, O and H belong to DMSO; O1 and H1 to water.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epsr
from epsr.simulation import render_overrides

from conftest import CompletedRun

DATA_FILES = ["D6DMSO_D2O.mint01", "D6DMSO_HDO.mint01", "D6DMSO_H2O.mint01"]
WEIGHTS_FILES = ["DMSODWaterD.NWTStot.wts", "DMSODWaterNull.NWTStot.wts", "DMSODWaterH.NWTStot.wts"]


# --------------------------------------------------------------- arrays and metadata
def test_grids_and_partials(completed_run: CompletedRun) -> None:
    simulation = completed_run.simulation
    nr, nq, n_pairs = simulation.nr, simulation.nq, simulation.n_pairs
    assert simulation.q.shape == (nq,) and simulation.r.shape == (nr,)
    assert simulation.q.dtype == np.float32 and simulation.r.dtype == np.float32
    np.testing.assert_allclose(np.diff(simulation.q), simulation.control.qstep, rtol=1e-4)

    expected_lengths = {"gr_inter": nr, "gr_intra": nr, "sq_inter": nq, "sq_intra": nq}
    for accessor_name, n_points in expected_lengths.items():
        accessor = getattr(simulation, accessor_name)
        all_pairs = accessor()
        assert all_pairs.shape == (n_pairs, n_points), accessor_name
        assert all_pairs.dtype == np.float32
        one_pair = accessor("O1", "H1")
        np.testing.assert_array_equal(one_pair, all_pairs[simulation.pair_index("O1", "H1")])

    assert simulation.current_gr().shape == (2, n_pairs, nr)
    assert simulation.current_sq().shape == (2, n_pairs, nq)
    # An intermolecular g(r) tends to 1 at large r.
    water_oxygen_gr = simulation.gr_inter("O1", "O1")
    assert float(np.mean(water_oxygen_gr[-50:])) == pytest.approx(1.0, abs=0.05)


def test_arrays_are_copies(completed_run: CompletedRun) -> None:
    simulation = completed_run.simulation
    returned_array = simulation.gr_inter()
    returned_array[:] = 0.0
    assert np.any(simulation.gr_inter() != 0.0)


def test_pair_labels_and_lookup(completed_run: CompletedRun) -> None:
    simulation = completed_run.simulation
    labels = simulation.pair_labels
    assert len(labels) == simulation.n_pairs == 21
    assert set(simulation.atom_types) == {atom_type for pair in labels for atom_type in pair}
    for index, (first_type, second_type) in enumerate(labels):
        assert simulation.pair_index(first_type, second_type) == index
        assert simulation.pair_index(second_type, first_type) == index
        assert simulation.pair_index((first_type, second_type)) == index
        assert simulation.pair_index(index) == index
    with pytest.raises(KeyError):
        simulation.pair_index("X", "Y")
    with pytest.raises(IndexError):
        simulation.pair_index(99)


def test_datasets_and_fits(completed_run: CompletedRun) -> None:
    simulation = completed_run.simulation
    datasets = simulation.datasets
    assert [dataset.infili for dataset in datasets] == DATA_FILES
    assert [Path(dataset.fnamedata).name for dataset in datasets] == DATA_FILES
    assert [dataset.name for dataset in datasets] == ["D6DMSO_D2O", "D6DMSO_HDO", "D6DMSO_H2O"]

    for dataset_index in range(simulation.n_datasets):
        first_q, last_q = simulation.data_q_range(dataset_index)
        assert 0 <= first_q < last_q < simulation.nq
        fit = simulation.fit(dataset_index)
        expected_shape = (last_q - first_q + 1,)
        assert fit.q.shape == fit.data.shape == fit.simulated.shape == fit.errors.shape
        assert fit.q.shape == expected_shape
        root_mean_square_residual = float(np.sqrt(np.mean(fit.residual**2)))
        assert root_mean_square_residual < 0.05  # the refined model fits the data closely

    assert len(simulation.summary.r_factor_per_dataset) == simulation.n_datasets
    with pytest.raises(IndexError):
        simulation.data_fq(3)
    np.testing.assert_array_equal(simulation.data_fq(-1), simulation.data_fq(2))


def test_configuration(completed_run: CompletedRun) -> None:
    simulation = completed_run.simulation
    atoms_per_dmso, atoms_per_water = 10, 3
    assert simulation.n_molecules == 600
    assert simulation.n_atoms == 200 * atoms_per_dmso + 400 * atoms_per_water

    positions = simulation.positions
    assert positions.shape == (simulation.n_atoms, 3) and positions.dtype == np.float32
    assert len(simulation.atom_labels) == simulation.n_atoms

    molecules = simulation.molecule_atoms
    assert len(molecules) == simulation.n_molecules
    assert molecules[0] == range(0, atoms_per_dmso)
    assert molecules[-1].stop == simulation.n_atoms

    assert simulation.box_size == pytest.approx(33.1213, abs=1e-3)
    assert simulation.density == pytest.approx(0.0880700, rel=1e-4)
    density_from_box = simulation.n_atoms / simulation.box_size**3
    assert density_from_box == pytest.approx(simulation.density, rel=1e-4)


def test_summary_and_control(completed_run: CompletedRun) -> None:
    simulation = completed_run.simulation
    summary = simulation.summary
    assert isinstance(summary, epsr.Summary)
    assert summary.iteration == 1
    assert summary.fortran["rfactot"] == pytest.approx(summary.r_factor)
    assert set(summary.as_dict()) >= {"energy", "pressure", "r_factor", "ep_amplitude"}

    control = simulation.control
    assert (control.nq, control.ntimes, control.niter) == (600, 5, 1)
    assert [Path(path).name for path in simulation.weights_files] == WEIGHTS_FILES
    assert simulation.output_path("g01") == completed_run.project.output("g01")
    assert simulation.output_path("g01").exists()
    assert simulation.completed and simulation.last_iteration == 1


# --------------------------------------------------------------- overrides
def test_render_overrides() -> None:
    assert render_overrides(None) == []
    assert render_overrides(["ntimes 3"]) == ["ntimes 3"]
    assert render_overrides({"ntimes": 3, "potfac": [1.0, 0.5], "flag": True}) == [
        "ntimes 3",
        "potfac 1.0 0.5",
        "flag T",
    ]
    lines = render_overrides({"data": {2: {"tweak": 1.05}}, "ereq": 12.5})
    assert lines == ["ereq 12.5", "data 2", "tweak 1.05"]
    with pytest.raises(ValueError):
        render_overrides({"data": {0: {"tweak": 1.0}}})
    with pytest.raises(TypeError):
        render_overrides({"ntimes": {"nested": 1}})


def test_dict_overrides_reach_engine(project: epsr.Project) -> None:
    overrides = {
        "ntimes": 2,
        "ereq": 15.0,
        "feedback": 0.8,
        "potfac": [0.5, 0.5],
        "data": {3: {"tweak": 1.1}},
    }
    simulation = project.simulation(overrides)
    assert simulation.control.ntimes == 2
    potential_parameters = simulation.potential.params
    assert potential_parameters.ereq == pytest.approx(15.0)
    assert potential_parameters.feedback == pytest.approx(0.8)
    assert simulation.get_parameter("ereq") == "15.0"
    assert simulation.datasets[2].tweak == pytest.approx(1.1)
    assert simulation.datasets[0].tweak == pytest.approx(1.0)


def test_overrides_are_saved_to_input(project: epsr.Project) -> None:
    project.simulation({"ntimes": 2}).run()
    assert project.input_file().values["ntimes"] == "2"
    assert project.simulation().control.ntimes == 2


# --------------------------------------------------------------- run-order rules
def test_iterate_rules(project: epsr.Project) -> None:
    simulation = project.simulation({"niter": 2})
    with pytest.raises(ValueError, match="out of range"):
        simulation.iterate(0)
    with pytest.raises(ValueError, match="in order"):
        simulation.iterate(2)

    first_summary = simulation.iterate(1)
    assert first_summary.iteration == 1 and not simulation.completed
    with pytest.raises(RuntimeError):
        simulation.run()

    simulation.iterate(2)
    assert simulation.completed
    with pytest.raises(ValueError, match="new Simulation"):
        simulation.iterate(3)
    assert list(simulation.iterations()) == []


def test_run_twice_is_refused(project: epsr.Project) -> None:
    simulation = project.simulation()
    simulation.run()
    with pytest.raises(RuntimeError):
        simulation.run()


# --------------------------------------------------------------- errors
def test_missing_folder(tmp_path: Path) -> None:
    with pytest.raises(epsr.EpsrIOError) as raised:
        epsr.Simulation(tmp_path / "nope", "x.EPSR.inp")
    assert isinstance(raised.value, OSError) and isinstance(raised.value, epsr.EpsrError)


def test_missing_input_file(project: epsr.Project) -> None:
    with pytest.raises(epsr.EpsrError):
        epsr.Simulation(project.folder, "missing.EPSR.inp")


def test_bad_keyword(project: epsr.Project) -> None:
    with pytest.raises(epsr.EpsrInputError) as raised:
        project.simulation(["no_such_keyword 1"])
    assert isinstance(raised.value, ValueError)


def test_inp_autodetect(project: epsr.Project, tmp_path: Path) -> None:
    assert epsr.Simulation(project.folder).inp_name == "200DMSO400Water.EPSR.inp"
    empty_folder = tmp_path / "empty"
    empty_folder.mkdir()
    with pytest.raises(FileNotFoundError):
        epsr.Simulation(empty_folder)


def test_truncated_configuration_is_an_epsr_error(project: epsr.Project) -> None:
    """An .ato file without the ecore/dcore line is rejected (the Fortran would hang)."""
    ato_file = project.folder / "200DMSO400Water.ato"
    first_line = ato_file.read_text(encoding="latin-1").splitlines()[0]
    ato_file.write_text(first_line + "\n", encoding="latin-1")
    with pytest.raises(epsr.EpsrError) as raised:
        project.simulation()
    # Truncating the file may surface as an input error or as "unsupported".
    assert isinstance(raised.value, (epsr.EpsrUnsupportedError, epsr.EpsrInputError))
    if isinstance(raised.value, epsr.EpsrUnsupportedError):
        assert isinstance(raised.value, NotImplementedError)
