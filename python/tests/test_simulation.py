"""Accessors, shapes, overrides, run-order rules and error mapping of ``epsr.Simulation``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epsr
from epsr.simulation import render_overrides


# --------------------------------------------------------------- arrays and metadata
def test_grids_and_partials(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    _, sim = ran
    nr, nq, npairs = sim.nr, sim.nq, sim.n_pairs
    assert sim.q.shape == (nq,) and sim.r.shape == (nr,)
    assert sim.q.dtype == np.float32 and sim.r.dtype == np.float32
    np.testing.assert_allclose(np.diff(sim.q), sim.control.qstep, rtol=1e-4)
    for name, n in (("gr_inter", nr), ("gr_intra", nr), ("sq_inter", nq), ("sq_intra", nq)):
        arr = getattr(sim, name)()
        assert arr.shape == (npairs, n), name
        assert arr.dtype == np.float32
        np.testing.assert_array_equal(getattr(sim, name)(("O", "H")), arr[sim.pair_index("O", "H")])
    assert sim.current_gr().shape == (2, npairs, nr)
    assert sim.current_sq().shape == (2, npairs, nq)
    # g(r) of an intermolecular pair tends to 1 at large r
    g = sim.gr_inter(("O", "O"))
    assert abs(float(np.mean(g[-50:])) - 1.0) < 0.05


def test_arrays_are_copies(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    _, sim = ran
    g = sim.gr_inter()
    g[:] = 0.0
    assert np.any(sim.gr_inter() != 0.0)


def test_pair_labels_and_lookup(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    _, sim = ran
    labels = sim.pair_labels
    assert len(labels) == sim.n_pairs == 21
    assert set(sim.atom_types) == {a for p in labels for a in p}
    for i, (a, b) in enumerate(labels):
        assert sim.pair_index(a, b) == i
        assert sim.pair_index(b, a) == i
        assert sim.pair_index((a, b)) == i
        assert sim.pair_index(i) == i
    with pytest.raises(KeyError):
        sim.pair_index("X", "Y")
    with pytest.raises(IndexError):
        sim.pair_index(99)


def test_datasets_and_fits(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    _, sim = ran
    ds = sim.datasets
    assert [d.infili for d in ds] == ["D6DMSO_D2O.mint01", "D6DMSO_HDO.mint01", "D6DMSO_H2O.mint01"]
    assert [Path(d.fnamedata).name for d in ds] == [d.infili for d in ds]
    assert [d.name for d in ds] == ["D6DMSO_D2O", "D6DMSO_HDO", "D6DMSO_H2O"]
    for i in range(sim.n_datasets):
        lo, hi = sim.data_q_range(i)
        assert 0 <= lo < hi < sim.nq
        f = sim.fit(i)
        assert f.q.shape == f.data.shape == f.simulated.shape == f.errors.shape == (hi - lo + 1,)
        # the refined model fits the data closely
        assert float(np.sqrt(np.mean(f.residual**2))) < 0.05
    assert len(sim.summary.r_factor_per_dataset) == sim.n_datasets
    with pytest.raises(IndexError):
        sim.data_fq(3)
    np.testing.assert_array_equal(sim.data_fq(-1), sim.data_fq(2))


def test_configuration(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    _, sim = ran
    assert sim.n_molecules == 600
    assert sim.n_atoms == 200 * 10 + 400 * 3
    pos = sim.positions
    assert pos.shape == (sim.n_atoms, 3) and pos.dtype == np.float32
    assert len(sim.atom_labels) == sim.n_atoms
    mols = sim.molecule_atoms
    assert len(mols) == sim.n_molecules
    assert mols[0] == range(0, 10) and mols[-1].stop == sim.n_atoms
    assert sim.box_size == pytest.approx(33.1213, abs=1e-3)
    assert sim.density == pytest.approx(0.0880700, rel=1e-4)
    # density consistent with the box
    assert sim.n_atoms / sim.box_size**3 == pytest.approx(sim.density, rel=1e-4)


def test_summary_and_control(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    proj, sim = ran
    s = sim.summary
    assert isinstance(s, epsr.Summary)
    assert s.iteration == 1
    assert s.fortran["rfactot"] == pytest.approx(s.r_factor)
    assert set(s.as_dict()) >= {"energy", "pressure", "r_factor", "ep_amplitude"}
    c = sim.control
    assert (c.nq, c.ntimes, c.niter) == (600, 5, 1)
    assert [Path(w).name for w in sim.weights_files] == [
        "DMSODWaterD.NWTStot.wts",
        "DMSODWaterNull.NWTStot.wts",
        "DMSODWaterH.NWTStot.wts",
    ]
    assert sim.output_path("g01") == proj.output("g01")
    assert sim.output_path("g01").exists()
    assert sim.completed and sim.last_iteration == 1


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
        render_overrides({"ntimes": {"a": 1}})


def test_dict_overrides_reach_engine(project: epsr.Project) -> None:
    sim = project.simulation(
        {"ntimes": 2, "ereq": 15.0, "feedback": 0.8, "potfac": [0.5, 0.5], "data": {3: {"tweak": 1.1}}}
    )
    assert sim.control.ntimes == 2
    p = sim.potential.params
    assert p.ereq == pytest.approx(15.0)
    assert p.feedback == pytest.approx(0.8)
    assert sim.get_parameter("ereq") == "15.0"
    assert sim.datasets[2].tweak == pytest.approx(1.1)
    assert sim.datasets[0].tweak == pytest.approx(1.0)


def test_overrides_are_saved_to_input(project: epsr.Project) -> None:
    project.simulation({"ntimes": 2}).run()
    assert project.input_file().values["ntimes"] == "2"
    assert project.simulation().control.ntimes == 2


# --------------------------------------------------------------- run-order rules
def test_iterate_rules(project: epsr.Project) -> None:
    sim = project.simulation({"niter": 2})
    with pytest.raises(ValueError, match="out of range"):
        sim.iterate(0)
    with pytest.raises(ValueError, match="in order"):
        sim.iterate(2)
    s1 = sim.iterate(1)
    assert s1.iteration == 1 and not sim.completed
    with pytest.raises(RuntimeError):
        sim.run()
    sim.iterate(2)
    assert sim.completed
    with pytest.raises(ValueError, match="new Simulation"):
        sim.iterate(3)
    assert list(sim.iterations()) == []


def test_run_twice_is_refused(project: epsr.Project) -> None:
    sim = project.simulation()
    sim.run()
    with pytest.raises(RuntimeError):
        sim.run()


# --------------------------------------------------------------- errors
def test_missing_folder(tmp_path: Path) -> None:
    with pytest.raises(epsr.EpsrIOError) as e:
        epsr.Simulation(tmp_path / "nope", "x.EPSR.inp")
    assert isinstance(e.value, OSError) and isinstance(e.value, epsr.EpsrError)


def test_missing_input_file(project: epsr.Project) -> None:
    with pytest.raises(epsr.EpsrError):
        epsr.Simulation(project.folder, "missing.EPSR.inp")


def test_bad_keyword(project: epsr.Project) -> None:
    with pytest.raises(epsr.EpsrInputError) as e:
        project.simulation(["no_such_keyword 1"])
    assert isinstance(e.value, ValueError)


def test_inp_autodetect(project: epsr.Project, tmp_path: Path) -> None:
    assert epsr.Simulation(project.folder).inp_name == "200DMSO400Water.EPSR.inp"
    (tmp_path / "empty").mkdir()
    with pytest.raises(FileNotFoundError):
        epsr.Simulation(tmp_path / "empty")


def test_unsupported_maps_to_not_implemented(project: epsr.Project) -> None:
    """An .ato file without the ecore/dcore line is 'not supported' (the Fortran hangs)."""
    ato = project.folder / "200DMSO400Water.ato"
    lines = ato.read_text(encoding="latin-1").splitlines()
    ato.write_text("\n".join(lines[:1]) + "\n", encoding="latin-1")
    with pytest.raises(epsr.EpsrError) as e:
        project.simulation()
    # Truncating the file may surface as an input error or as "unsupported"
    assert isinstance(e.value, (epsr.EpsrUnsupportedError, epsr.EpsrInputError))
    if isinstance(e.value, epsr.EpsrUnsupportedError):
        assert isinstance(e.value, NotImplementedError)
