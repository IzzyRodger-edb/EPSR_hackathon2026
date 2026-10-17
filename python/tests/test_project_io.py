"""``epsr.Project`` folder handling and the read-only parsers in ``epsr.io``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epsr
from epsr.io import ERG_COLUMNS, read_erg, read_keyword_file, read_table, read_table_series

from conftest import CompletedRun

DATA_FILES = ["D6DMSO_D2O.mint01", "D6DMSO_HDO.mint01", "D6DMSO_H2O.mint01"]
DATASET_NAMES = ["D6DMSO_D2O", "D6DMSO_HDO", "D6DMSO_H2O"]
DENSITY_COLUMN = ERG_COLUMNS.index("density")


def test_from_template_copies_and_protects(template: Path, tmp_path: Path) -> None:
    work_folder = tmp_path / "work"
    project = epsr.Project.from_template(template, work_folder)
    template_files = {path.name for path in template.iterdir() if path.is_file()}
    assert {path.name for path in project.folder.iterdir()} == template_files
    assert project.inp_name == "200DMSO400Water.EPSR.inp"
    assert project.run_name == "200DMSO400Water.EPSR"

    with pytest.raises(FileExistsError):
        epsr.Project.from_template(template, work_folder)

    marker = project.folder / "marker"
    marker.write_text("x")
    epsr.Project.from_template(template, work_folder, overwrite=True)
    assert not marker.exists()

    with pytest.raises(FileNotFoundError):
        epsr.Project.from_template(tmp_path / "missing", tmp_path / "other")


def test_template_is_untouched_by_runs(template: Path, project: epsr.Project) -> None:
    template_ato = template / "200DMSO400Water.ato"
    original_content = template_ato.read_bytes()
    project.simulation().run()
    assert template_ato.read_bytes() == original_content
    assert (project.folder / "200DMSO400Water.ato").read_bytes() != original_content


def test_snapshot(project: epsr.Project) -> None:
    snapshot_folder = project.snapshot("start")
    assert snapshot_folder == project.folder.parent / "work.snapshots" / "start"
    snapshot_input = (snapshot_folder / project.inp_name).read_bytes()
    assert snapshot_input == (project.folder / project.inp_name).read_bytes()
    with pytest.raises(FileExistsError):
        project.snapshot("start")
    project.snapshot("start", overwrite=True)


def test_files(project: epsr.Project) -> None:
    paths = project.files()
    assert paths["ato"] == project.folder / "200DMSO400Water.ato"
    assert paths["pcof"] == project.folder / "200DMSO400Water.pcof"
    assert [path.name for path in paths["data"]] == DATA_FILES
    assert len(paths["weights"]) == 3
    assert all(path.exists() for path in paths["data"] + paths["weights"])


def test_read_keyword_file_inp(template: Path) -> None:
    input_file = read_keyword_file(template / "200DMSO400Water.EPSR.inp")
    assert input_file.title == "200DMSO400Water.EPSR"
    assert input_file.values["ntimes"] == "5"
    assert input_file.words("potfac") == ["1.000000", "1.000000"]
    assert input_file.values["control-p"] == "none 0 0"
    assert input_file.block_keyword == "data"
    assert len(input_file.blocks) == 3
    assert input_file.blocks[2]["datafile"] == "D6DMSO_H2O.mint01"
    assert input_file.descriptions["ntimes"].startswith("Number of MC cycles")


def test_read_keyword_file_pcof_and_nwts(template: Path) -> None:
    pcof_file = read_keyword_file(template / "200DMSO400Water.pcof")
    assert float(pcof_file.values["ereq"]) == pytest.approx(10.0)
    assert pcof_file.values["Gaussian"] == "Poisson"
    assert pcof_file.trailer[0].strip() == "21"  # the number of pairs follows the `q` line

    nwts_file = read_keyword_file(template / "DMSODWaterD.NWTS.dat")
    assert nwts_file.block_keyword == "component"
    assert len(nwts_file.blocks) == int(nwts_file.values["ncomponent"]) == 6
    assert nwts_file.blocks[0]["atom"] == "S"


def test_read_table(completed_run: CompletedRun) -> None:
    project, simulation = completed_run
    gr_table = read_table(project.output("g01"))
    assert gr_table.x_label == "r"
    assert gr_table.labels == [f"{first}-{second}" for first, second in simulation.pair_labels]
    assert set(gr_table.errors) == set(gr_table.columns)
    np.testing.assert_allclose(
        gr_table["O1-H1"][: simulation.nr], simulation.gr_inter("O1", "H1"), atol=1e-6
    )

    total_fq_table = project.read_table("u01")
    assert total_fq_table.labels == DATASET_NAMES
    np.testing.assert_allclose(total_fq_table["D6DMSO_D2O"], simulation.simulated_fq(0), atol=1e-6)

    coordination_table = read_table(project.output("z01"))  # values only, each pair twice
    assert not coordination_table.errors and "S-S#2" in coordination_table.columns
    assert len(read_table_series(project.output("g01"))) == len(gr_table)


def test_table_to_dataframe(completed_run: CompletedRun) -> None:
    pytest.importorskip("pandas")
    dataframe = completed_run.project.read_table("u01").to_dataframe()
    assert list(dataframe.columns)[:3] == ["Q", "D6DMSO_D2O", "D6DMSO_D2O err"]


def test_read_erg(template: Path) -> None:
    erg_data = read_erg(template / "200DMSO400Water.EPSR.erg")
    assert erg_data.shape[1] == len(ERG_COLUMNS) == 10
    assert erg_data.shape[0] > 100
    assert np.all(erg_data[:, DENSITY_COLUMN] == pytest.approx(0.08807, rel=1e-4))


def test_read_table_rejects_non_table(template: Path) -> None:
    with pytest.raises(ValueError):
        read_table(template / "200DMSO400Water.EPSR.erg")
