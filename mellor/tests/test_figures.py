"""``mellor.figures`` on the example water run: every figure is a non-empty PNG."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest
from example_run import build_example_result, example_available

from mellor import figures
from mellor.figures import FIGURE_CAPTIONS
from mellor.runner import RunResult

pytestmark = pytest.mark.skipif(
    not example_available(), reason="toolkit example data or water data not found"
)

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
TRAJECTORY_FIGURES = {
    "01_convergence",
    "02_fq_fit",
    "03_fq_evolution",
    "04_residual_map",
    "05_partials",
    "06_potentials",
    "07_first_shell",
    "08_local_structure",
    "09_displacement",
}


@pytest.fixture(scope="module")
def example_figures(
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[RunResult, dict[str, Path]]:
    result = build_example_result(tmp_path_factory.mktemp("example_run"))
    return result, figures.make_figures(result)


def test_every_figure_is_a_png(
    example_figures: tuple[RunResult, dict[str, Path]],
) -> None:
    result, written = example_figures
    expected = TRAJECTORY_FIGURES | {
        "00_data",
        "10_final_outputs",
        "11_fit_final",
        "summary.csv",
        "summary.md",
    }
    assert set(written) == expected
    for name, path in written.items():
        assert path.parent == result.folder / "figures", name
        assert path.stat().st_size > 0, name
        if path.suffix == ".png":
            assert path.stem == name
            assert path.read_bytes()[:8] == PNG_MAGIC, name
            assert name in FIGURE_CAPTIONS, name
    summary = (result.folder / "figures" / "summary.md").read_text()
    assert "rfactor_total" in summary and "iter 300" in summary
    assert (
        (result.folder / "figures" / "summary.csv").read_text().startswith("iteration,")
    )


def test_captions_cover_every_figure() -> None:
    from epsrvis.report import FIGURES

    assert set(FIGURES) <= set(FIGURE_CAPTIONS)
    assert {"00_data", "10_final_outputs", "11_fit_final"} <= set(FIGURE_CAPTIONS)
    assert all(caption.endswith(".") for caption in FIGURE_CAPTIONS.values())


def test_out_dir_and_dpi(
    example_figures: tuple[RunResult, dict[str, Path]], tmp_path: Path
) -> None:
    result, default = example_figures
    written = figures.make_figures(result, tmp_path / "figs", dpi=50)
    assert written["11_fit_final"] == tmp_path / "figs" / "11_fit_final.png"
    assert set(written) == set(default)
    for name, path in written.items():
        assert path.stat().st_size > 0, name
    # a lower resolution gives a smaller file
    assert written["02_fq_fit"].stat().st_size < default["02_fq_fit"].stat().st_size


def test_without_trajectory(tmp_path: Path) -> None:
    result = build_example_result(tmp_path / "run", with_trajectory=False)
    written = figures.make_figures(result)
    assert set(written) == {"00_data", "10_final_outputs"}
    for path in written.values():
        assert path.read_bytes()[:8] == PNG_MAGIC
    assert not (tmp_path / "run" / "figures" / "summary.csv").exists()


def test_data_figure_falls_back_to_the_original_files(tmp_path: Path) -> None:
    result = build_example_result(
        tmp_path / "run", with_trajectory=False, with_epsr=False, with_data=False
    )
    written = figures.make_figures(result)
    assert set(written) == {"00_data"}
    assert written["00_data"].read_bytes()[:8] == PNG_MAGIC


def test_nothing_to_plot_is_logged_not_raised(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    result = build_example_result(
        tmp_path / "run", with_trajectory=False, with_epsr=False, with_data=False
    )
    for dataset in result.spec.datasets:
        dataset.file = str(tmp_path / "gone" / dataset.file_name)
    with caplog.at_level(logging.WARNING, logger="mellor.figures"):
        written = figures.make_figures(result)
    assert written == {}
    assert any("00_data skipped" in record.message for record in caplog.records)


def test_case_name_is_set_for_a_trajectory_without_one(tmp_path: Path) -> None:
    """A trajectory recorded without a case name gets the script's (water here)."""
    import json

    result = build_example_result(tmp_path / "run", with_epsr=False, with_data=False)
    # replace the symlinked trajectory by a copy whose meta.json has no case
    source = result.trajectory_path
    assert source is not None
    meta = json.loads((source / "meta.json").read_text())
    copy = tmp_path / "run" / "trajectory_copy"
    copy.mkdir()
    for path in source.iterdir():
        if path.name != "meta.json":
            (copy / path.name).symlink_to(path.resolve())
    meta["case"] = ""
    (copy / "meta.json").write_text(json.dumps(meta))
    result.manifest["trajectory"] = "trajectory_copy"
    traj = figures._trajectory(result)
    assert traj is not None and traj.case == "water"
    result.manifest["epsrvis_case"] = None
    result.spec.atom_types = [result.spec.atom_types[0]]  # no longer the water types
    traj = figures._trajectory(result)
    assert traj is not None and traj.case == "waterbox"
