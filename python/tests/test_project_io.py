"""``epsr.Project`` folder handling and the read-only parsers in ``epsr.io``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epsr
from epsr.io import ERG_COLUMNS, read_erg, read_keyword_file, read_table, read_table_series


def test_from_template_copies_and_protects(template: Path, tmp_path: Path) -> None:
    proj = epsr.Project.from_template(template, tmp_path / "w")
    names = {p.name for p in template.iterdir() if p.is_file()}
    assert {p.name for p in proj.folder.iterdir()} == names
    assert proj.inp_name == "200DMSO400Water.EPSR.inp"
    assert proj.run_name == "200DMSO400Water.EPSR"
    with pytest.raises(FileExistsError):
        epsr.Project.from_template(template, tmp_path / "w")
    (proj.folder / "marker").write_text("x")
    proj2 = epsr.Project.from_template(template, tmp_path / "w", overwrite=True)
    assert not (proj2.folder / "marker").exists()
    with pytest.raises(FileNotFoundError):
        epsr.Project.from_template(tmp_path / "missing", tmp_path / "w2")


def test_template_is_untouched_by_runs(template: Path, project: epsr.Project) -> None:
    before = (template / "200DMSO400Water.ato").read_bytes()
    project.simulation().run()
    assert (template / "200DMSO400Water.ato").read_bytes() == before
    assert (project.folder / "200DMSO400Water.ato").read_bytes() != before


def test_snapshot(project: epsr.Project) -> None:
    snap = project.snapshot("start")
    assert snap == project.folder.parent / "work.snapshots" / "start"
    assert (snap / project.inp_name).read_bytes() == (project.folder / project.inp_name).read_bytes()
    with pytest.raises(FileExistsError):
        project.snapshot("start")
    project.snapshot("start", overwrite=True)


def test_files(project: epsr.Project) -> None:
    f = project.files()
    assert f["ato"] == project.folder / "200DMSO400Water.ato"
    assert f["pcof"] == project.folder / "200DMSO400Water.pcof"
    assert [p.name for p in f["data"]] == ["D6DMSO_D2O.mint01", "D6DMSO_HDO.mint01", "D6DMSO_H2O.mint01"]
    assert len(f["weights"]) == 3
    assert all(p.exists() for p in f["data"] + f["weights"])


def test_read_keyword_file_inp(template: Path) -> None:
    kw = read_keyword_file(template / "200DMSO400Water.EPSR.inp")
    assert kw.title == "200DMSO400Water.EPSR"
    assert kw.values["ntimes"] == "5"
    assert kw.words("potfac") == ["1.000000", "1.000000"]
    assert kw.values["control-p"] == "none 0 0"
    assert kw.block_keyword == "data"
    assert len(kw.blocks) == 3
    assert kw.blocks[2]["datafile"] == "D6DMSO_H2O.mint01"
    assert kw.descriptions["ntimes"].startswith("Number of MC cycles")


def test_read_keyword_file_pcof_and_nwts(template: Path) -> None:
    pcof = read_keyword_file(template / "200DMSO400Water.pcof")
    assert float(pcof.values["ereq"]) == pytest.approx(10.0)
    assert pcof.values["Gaussian"] == "Poisson"
    assert pcof.trailer[0].strip() == "21"  # number of pairs follows the `q` line
    nw = read_keyword_file(template / "DMSODWaterD.NWTS.dat")
    assert nw.block_keyword == "component"
    assert len(nw.blocks) == int(nw.values["ncomponent"]) == 6
    assert nw.blocks[0]["atom"] == "S"


def test_read_table(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    proj, sim = ran
    g = read_table(proj.output("g01"))
    assert g.x_label == "r"
    assert g.labels == [f"{a}-{b}" for a, b in sim.pair_labels]
    assert set(g.errors) == set(g.columns)
    np.testing.assert_allclose(g["O-H"][: sim.nr], sim.gr_inter("O", "H"), atol=1e-6)
    u = proj.read_table("u01")
    assert u.labels == ["D6DMSO_D2O", "D6DMSO_HDO", "D6DMSO_H2O"]
    np.testing.assert_allclose(u["D6DMSO_D2O"], sim.simulated_fq(0), atol=1e-6)
    z = read_table(proj.output("z01"))  # values only, each pair listed twice
    assert not z.errors and "S-S#2" in z.columns
    assert len(read_table_series(proj.output("g01"))) == len(g)


def test_table_to_dataframe(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    pytest.importorskip("pandas")
    proj, _ = ran
    df = proj.read_table("u01").to_dataframe()
    assert list(df.columns)[:3] == ["Q", "D6DMSO_D2O", "D6DMSO_D2O err"]


def test_read_erg(template: Path) -> None:
    erg = read_erg(template / "200DMSO400Water.EPSR.erg")
    assert erg.shape[1] == len(ERG_COLUMNS) == 10
    assert erg.shape[0] > 100
    assert np.all(erg[:, 4] == pytest.approx(0.08807, rel=1e-4))  # density column


def test_read_table_rejects_non_table(template: Path) -> None:
    with pytest.raises(ValueError):
        read_table(template / "200DMSO400Water.EPSR.erg")
