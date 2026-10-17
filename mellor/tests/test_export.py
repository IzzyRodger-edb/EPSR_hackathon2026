"""Tidy-table export on the toolkit's example data (water, 300 iterations)."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import pytest

from mellor import export

EXAMPLE = (
    Path(__file__).resolve().parents[2]
    / "visualisation"
    / "python_visualisation_toolkit"
    / "example_H2O"
    / "data"
)


@pytest.fixture(scope="module")
def run_folder(tmp_path_factory: pytest.TempPathFactory) -> Path:
    if not (EXAMPLE / "trajectory" / "meta.json").is_file():
        pytest.skip("toolkit example data not found")
    folder = tmp_path_factory.mktemp("run")
    os.symlink(EXAMPLE / "trajectory", folder / "trajectory")
    os.symlink(EXAMPLE / "epsr_output", folder / "epsr")
    (folder / "manifest.json").write_text(
        json.dumps(
            {
                "epsr_folder": "epsr",
                "runs": [
                    {
                        "run": 1,
                        "stage": "scratch",
                        "niter": 1,
                        "ereq": 0.0,
                        "duration_s": 1.0,
                        "summary": {
                            "r_factor": 0.021,
                            "energy": -10.0,
                            "pressure": 1.0,
                            "ep_amplitude": 0.0,
                        },
                    },
                    {
                        "run": 2,
                        "stage": "reference potential",
                        "niter": 99,
                        "ereq": 0.0,
                        "duration_s": 99.0,
                        "summary": {
                            "r_factor": 6.1e-4,
                            "energy": -40.0,
                            "pressure": 5.0,
                            "ep_amplitude": 0.0,
                        },
                    },
                ],
            }
        )
    )
    return folder


def test_tables_from_trajectory(run_folder: Path) -> None:
    tables = export.tidy_tables(run_folder)
    assert set(tables) == set(export.TABLES)
    fits = tables["fits"]
    assert {row["dataset"] for row in fits} == {
        "SLS18502-D2O",
        "SLS18500-HDO5050",
        "SLS18498-H2O",
    }
    assert all(row["error"] is not None for row in fits)
    assert all(
        abs(row["residual"] - (row["data"] - row["simulated"])) < 1e-9 for row in fits
    )
    partials = tables["partials"]
    assert {row["pair"] for row in partials} == {"OW-OW", "OW-HW", "HW-HW"}
    history = tables["history"]
    assert len(history) == 300
    assert history[-1]["r_factor"] < history[0]["r_factor"]
    assert "r_factor[SLS18498-H2O]" in history[-1]
    assert len(tables["runs"]) == 2


def test_write_csv(run_folder: Path, tmp_path: Path) -> None:
    written = export.write_tables(run_folder, tmp_path / "tables")
    assert set(written) == set(export.TABLES)
    with (tmp_path / "tables" / "fits.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows and rows[0]["dataset"] and float(rows[0]["q"]) > 0


def test_fallback_to_epsr_files(tmp_path: Path) -> None:
    if not (EXAMPLE / "epsr_output").is_dir():
        pytest.skip("toolkit example data not found")
    os.symlink(EXAMPLE / "epsr_output", tmp_path / "epsr")
    tables = export.tidy_tables(tmp_path)
    assert tables["fits"] and tables["fits"][0]["error"] is None
    assert tables["partials"] and tables["potentials"]
    assert tables["history"] and tables["history"][0]["run"] is None
    assert tables["runs"] == []
