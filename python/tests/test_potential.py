"""Reading and perturbing the potentials through ``Simulation.potential``.

The empirical potential (EP) of the test case is at its amplitude limit,
``ereq = 10``; several tests depend on that (see ``epsr/potential.py``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epsr

WATER_HYDROGEN_BOND = ("O1", "H1")


@pytest.fixture
def simulation(project: epsr.Project) -> epsr.Simulation:
    return project.simulation()


def test_read_potentials(simulation: epsr.Simulation) -> None:
    potential = simulation.potential
    table_shape = (simulation.n_pairs, simulation.nr)
    assert potential.empirical().shape == table_shape
    assert potential.reference().shape == table_shape
    assert potential.total().shape == table_shape
    np.testing.assert_array_equal(potential.r, simulation.r)

    # total = reference + empirical, to single precision
    np.testing.assert_allclose(
        potential.total(), potential.reference() + potential.empirical(), rtol=1e-6, atol=1e-4
    )

    pair_index = simulation.pair_index(*WATER_HYDROGEN_BOND)
    np.testing.assert_array_equal(
        potential.empirical(*WATER_HYDROGEN_BOND), potential.empirical()[pair_index]
    )
    np.testing.assert_array_equal(
        potential.reference(("H1", "O1")), potential.reference()[pair_index]
    )

    parameters = potential.params
    assert parameters.ncoeffp == 116
    assert parameters.ereq == pytest.approx(10.0)
    assert not parameters.gaussian
    assert potential.ep_coefficients().shape == (simulation.n_pairs, parameters.ncoeffp - 1)
    np.testing.assert_array_equal(potential.ep_flags(), np.ones(simulation.n_pairs, np.float32))
    assert potential.shared_pairs() == []


def test_type_params(simulation: epsr.Simulation) -> None:
    parameters_by_type = {
        parameters.label: parameters for parameters in simulation.potential.type_params
    }
    assert list(parameters_by_type) == simulation.atom_types
    assert parameters_by_type["S"].sigma == pytest.approx(3.56)
    assert parameters_by_type["S"].charge == pytest.approx(0.54)
    assert parameters_by_type["H1"].mass == pytest.approx(2.0)  # deuterium: D6-DMSO in D2O

    all_types = parameters_by_type.values()
    total_charge = sum(parameters.charge * parameters.fraction for parameters in all_types)
    assert total_charge == pytest.approx(0.0, abs=1e-5)  # neutral system
    assert sum(parameters.fraction for parameters in all_types) == pytest.approx(1.0, abs=1e-5)


def test_set_get_round_trip_is_exact(simulation: epsr.Simulation) -> None:
    potential = simulation.potential
    original_ep = potential.empirical()
    original_coefficients = potential.ep_coefficients()
    potential.set_ep_coefficients(original_coefficients)
    np.testing.assert_array_equal(potential.empirical(), original_ep)
    np.testing.assert_array_equal(potential.ep_coefficients(), original_coefficients)


def test_scale_down_is_linear(simulation: epsr.Simulation) -> None:
    potential = simulation.potential
    original_ep = potential.empirical()
    original_coefficients = potential.ep_coefficients()
    potential.scale_ep(0.5)
    np.testing.assert_allclose(potential.ep_coefficients(), 0.5 * original_coefficients, rtol=1e-6)
    np.testing.assert_allclose(potential.empirical(), 0.5 * original_ep, rtol=1e-4, atol=1e-5)


def test_ereq_limit_is_applied(simulation: epsr.Simulation) -> None:
    """The EP sits at ereq = 10, so the engine undoes scaling beyond it."""
    potential = simulation.potential
    original_coefficients = potential.ep_coefficients()
    potential.scale_ep(2.0)
    np.testing.assert_allclose(
        potential.ep_coefficients(), original_coefficients, rtol=1e-5, atol=1e-7
    )


def test_raising_ereq_allows_larger_ep(project: epsr.Project) -> None:
    original_coefficients = project.simulation().potential.ep_coefficients()
    potential = project.simulation({"ereq": 20.0}).potential
    potential.scale_ep(2.0)
    np.testing.assert_allclose(
        potential.ep_coefficients(), 2.0 * original_coefficients, rtol=1e-5, atol=1e-7
    )


def test_zero_ep(simulation: epsr.Simulation) -> None:
    potential = simulation.potential
    potential.zero_ep()
    assert not np.any(potential.empirical())
    assert not np.any(potential.ep_coefficients())
    np.testing.assert_array_equal(potential.total(), potential.reference())


def test_perturb_one_pair(simulation: epsr.Simulation) -> None:
    potential = simulation.potential
    potential.scale_ep(0.5)  # leave room below the ereq limit
    coefficients_before = potential.ep_coefficients()
    ep_before = potential.empirical()
    pair_index = simulation.pair_index(*WATER_HYDROGEN_BOND)

    delta = np.zeros(coefficients_before.shape[1], np.float32)
    delta[10:20] = 0.01
    potential.perturb_ep(WATER_HYDROGEN_BOND, delta)

    coefficients_after = potential.ep_coefficients()
    changed_pairs = np.nonzero(np.any(coefficients_after != coefficients_before, axis=1))[0]
    assert changed_pairs.tolist() == [pair_index]
    np.testing.assert_allclose(
        coefficients_after[pair_index] - coefficients_before[pair_index], delta, atol=1e-6
    )

    ep_after = potential.empirical()
    assert np.any(ep_after[pair_index] != ep_before[pair_index])
    np.testing.assert_array_equal(
        np.delete(ep_after, pair_index, axis=0), np.delete(ep_before, pair_index, axis=0)
    )


def test_set_flags_and_bad_shapes(simulation: epsr.Simulation) -> None:
    potential = simulation.potential
    potential.set_ep_flags(0.0)
    assert not np.any(potential.ep_flags())
    assert not np.any(potential.empirical())  # the flag multiplies the EP of each pair
    with pytest.raises(ValueError):
        potential.set_ep_coefficients(np.zeros((3, 3)))


def test_edits_persist_to_pcof(project: epsr.Project) -> None:
    # potfac 0 switches refinement off, so the coefficients are saved exactly as set.
    simulation = project.simulation({"potfac": [0.0, 0.0]})
    halved_coefficients = 0.5 * simulation.potential.ep_coefficients()
    simulation.potential.set_ep_coefficients(halved_coefficients)
    simulation.run()
    coefficients_read_back = project.simulation().potential.ep_coefficients()
    np.testing.assert_allclose(coefficients_read_back, halved_coefficients, rtol=1e-6, atol=1e-9)


def test_edited_ep_changes_the_simulation(project: epsr.Project, template: Path) -> None:
    unchanged_project = epsr.Project.from_template(template, project.folder.parent / "unchanged")
    unchanged_summary = unchanged_project.simulation().run()
    simulation = project.simulation()
    simulation.potential.zero_ep()
    edited_summary = simulation.run()
    assert edited_summary.energy != unchanged_summary.energy
    assert edited_summary.r_factor != unchanged_summary.r_factor


def test_refinement_updates_coefficients(simulation: epsr.Simulation) -> None:
    coefficients_before = simulation.potential.ep_coefficients()
    simulation.run()
    assert np.any(simulation.potential.ep_coefficients() != coefficients_before)
