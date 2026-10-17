"""Tests of ``mellor.data`` on the measured ``.mint01`` files of the repository."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from mellor.data import (
    Mint01,
    describe,
    describe_all,
    guess_sample,
    plot_data,
    read_mint01,
)

ROOT = Path(__file__).resolve().parents[2]
WATER_DIR = ROOT / "experimental_data" / "water" / "data"
SILICA_FILE = ROOT / "experimental_data" / "silica_bulk" / "data" / "SiO2_NIMROD.mint01"
CU_BTC_FILE = ROOT / "experimental_data" / "cu-btc" / "data" / "Empty-CuBTC.mint01"
EPSR_TEST = ROOT / "EPSR-test"
WATER_FILES = sorted(WATER_DIR.glob("*.mint01"))
DMSO_FILES = sorted(EPSR_TEST.glob("D6DMSO_*.mint01"))
ALL_FILES = WATER_FILES + [SILICA_FILE, CU_BTC_FILE] + DMSO_FILES

HEADER = "# {source}\n# {title}\n#     1  {n}\n#    10\n# {user}\n# {date}\n# spec.bad\n# groups_def.dat\n"


def write_mint01(
    path: Path,
    title: str,
    rows: list[tuple[float, float, float]] | None = None,
    *,
    source: str = "SLS00001.mint01",
    user: str = "Tester",
    date: str = "01-JAN-2020 12:00:00",
) -> Path:
    """Write a small synthetic ``.mint01`` file in Gudrun's layout."""
    rows = (
        rows
        if rows is not None
        else [(0.1, -0.2, 0.01), (0.2, -0.1, 0.01), (0.3, 0.1, 0.02), (0.4, 0.0, 0.0)]
    )
    text = HEADER.format(source=source, title=title, n=len(rows), user=user, date=date)
    text += "".join(f"  {q:.7E}  {f:.7E}  {e:.7E}\n" for q, f, e in rows)
    path.write_text(text)
    return path


# --------------------------------------------------------------------------- real files
def test_water_h2o_file() -> None:
    data = read_mint01(WATER_DIR / "SLS18498-H2O.mint01")
    assert data.name == "SLS18498-H2O.mint01"
    assert data.source == "sls18498.mdcs01"
    assert data.title == "H2O in 1mm can 5 in position 5 : 25C"
    assert data.header[0] == data.source
    assert data.header[1] == data.title
    assert "SANDALS" in data.instrument
    assert data.run_number == "18498"
    assert data.user == "Bychkov"
    assert data.date == "28-FEB-2001 16:15:03"
    assert len(data.q) == 997
    assert data.n_points == 996  # the last row is 0 0 padding
    assert not data.valid[-1] and data.valid[:-1].all()
    assert data.q_range == pytest.approx((0.175, 49.925))
    assert data.q_step == pytest.approx(0.05)
    assert data.f[0] == pytest.approx(-0.3225306)
    assert data.err[0] == pytest.approx(6.120029e-3)


def test_water_runs_and_titles() -> None:
    runs = {read_mint01(path).run_number for path in WATER_FILES}
    assert runs == {"18498", "18500", "18502"}
    hdo = read_mint01(WATER_DIR / "SLS18500-HDO5050.mint01")
    assert hdo.title == "50/50 H2O/D2O in 2mm can 11 in position"
    d2o = read_mint01(WATER_DIR / "SLS18502-D2O.mint01")
    assert d2o.title.startswith("D2O in 4mm can")


def test_silica_file_padding() -> None:
    data = read_mint01(SILICA_FILE)
    assert data.source == "NIMROD00001122.mint01"
    assert "NIMROD" in data.instrument
    assert data.run_number == "1122"
    assert data.title == "3mm SiO2 flat plate in pos 9"
    assert data.user == "" and data.date == ""
    assert len(data.q) == 1002
    assert not data.valid[0]  # first row 0 0 0
    assert not data.valid[-1] and not data.valid[-2]  # trailing 0 0 padding
    assert data.n_points == 999
    assert data.q_range == pytest.approx((0.025, 49.925))
    assert data.q_step == pytest.approx(0.05)


def test_cu_btc_file() -> None:
    data = read_mint01(CU_BTC_FILE)
    assert data.run_number == "24401"
    assert "NIMROD" in data.instrument
    assert data.title == "CuBTC Pellets in Cell 959 on superstick"
    assert data.user == "Youngs TGA, Sheppard"
    assert data.date == "24-OCT-2013 10:59:00"
    assert data.n_points == len(data.q) - 1
    assert data.q_range[0] == pytest.approx(0.2534767)
    assert data.q_range[1] == pytest.approx(49.98219)
    assert 0.02 < data.q_step < 0.04  # non-uniform grid, median spacing


def test_dmso_files() -> None:
    datasets = [read_mint01(path) for path in DMSO_FILES]
    assert [data.run_number for data in datasets] == ["3863", "3850", "3871"]
    assert all("SANDALS" in data.instrument for data in datasets)
    assert datasets[0].title == "DMSO(2):2D2O IN 1MM CAN"
    assert datasets[2].title == 'DMSO(D) in 2"Null" waters in 1mm can'
    assert datasets[1].date == "9-SEP-1991 16:34:04"
    assert all(data.n_points == 594 for data in datasets)


@pytest.mark.parametrize("path", ALL_FILES, ids=[path.name for path in ALL_FILES])
def test_describe_every_file(path: Path) -> None:
    data = read_mint01(path)
    description = describe(data)
    json.dumps(description)  # JSON-serialisable
    expected_keys = {
        "file",
        "path",
        "source",
        "title",
        "instrument",
        "run_number",
        "user",
        "date",
        "n_points",
        "q_min",
        "q_max",
        "q_step",
        "f_min",
        "f_max",
        "f_high_q",
        "mean_error",
        "sample",
    }
    assert expected_keys <= set(description)
    assert description["file"] == path.name
    assert Path(description["path"]).is_absolute()
    assert description["n_points"] == int(data.valid.sum()) > 500
    assert description["q_min"] < 0.4 < 25.0 < description["q_max"]
    assert description["f_min"] <= description["f_high_q"] <= description["f_max"]
    assert abs(description["f_high_q"]) < 0.01  # the high-Q level is about zero
    assert description["mean_error"] > 0
    assert set(description["sample"]) >= {
        "keywords",
        "hydrogen",
        "hd_fraction",
        "candidates",
        "sample",
    }
    assert description["sample"][
        "candidates"
    ]  # every repository file names a known system


def test_describe_all_and_high_q_level() -> None:
    descriptions = describe_all([str(path) for path in WATER_FILES])
    assert [item["file"] for item in descriptions] == [
        path.name for path in WATER_FILES
    ]
    data = read_mint01(WATER_FILES[0])
    q, f = data.q[data.valid], data.f[data.valid]
    expected = f[q >= q.max() - 0.1 * (q.max() - q.min())].mean()
    assert descriptions[0]["f_high_q"] == pytest.approx(expected, abs=1e-6)


# --------------------------------------------------------------------------- sample hints
@pytest.mark.parametrize(
    ("title", "file_name", "hydrogen", "fraction", "candidates", "sample"),
    [
        ("H2O in 1mm can 5 in position 5 : 25C", "", "H", 0.0, ["water"], "H2O"),
        ("D2O in 4mm can 15 in position 9 : 25 C", "", "D", 1.0, ["water"], "D2O"),
        ("50/50 H2O/D2O in 2mm can 11 in position", "", "HD", 0.5, ["water"], "HDO"),
        ("HDO in a can", "", "HD", 0.5, ["water"], "HDO"),
        ("water 5050 mixture", "", "HD", 0.5, ["water"], "HDO"),
        ("50:50 H2O:D2O", "", "HD", 0.5, ["water"], "HDO"),
        ("75:25 D2O:H2O in can", "", "HD", 0.75, ["water"], "HDO"),
        ("25/75 H2O/D2O in can", "", "HD", 0.75, ["water"], "HDO"),
        ("liquid water at 298 K", "", None, None, ["water"], None),
        ("heavy water", "", "D", 1.0, ["water"], "D2O"),
        ("3mm SiO2 flat plate in pos 9", "", None, None, ["silica"], None),
        ("amorphous silica", "", None, None, ["silica"], None),
        ("CuBTC Pellets in Cell 959 on superstick", "", None, None, ["cu_btc"], None),
        ("Cu-BTC empty framework", "", None, None, ["cu_btc"], None),
        ("DMSO(2):2D2O IN 1MM CAN", "", "D", 1.0, ["dmso_water"], "D2O"),
        ("DMSO(2):2H2O IN 1MM CAN", "", "H", 0.0, ["dmso_water"], "H2O"),
        (
            'DMSO(D) in 2"Null" waters in 1mm can',
            "D6DMSO_HDO.mint01",
            "HD",
            0.36,
            ["dmso_water"],
            "Null",
        ),
        ("benzene in 1mm can", "", None, None, ["benzene"], None),
        ("C6H6 benzene", "", "H", 0.0, ["benzene"], "C6H6"),
        ("C6D6 benzene", "", "D", 1.0, ["benzene"], "C6D6"),
        ("5050 benzene", "", "HD", 0.5, ["benzene"], "HD50"),
        ("argon at 85 K", "", None, None, ["argon"], None),
        ("", "SLS18502-D2O.mint01", "D", 1.0, ["water"], "D2O"),  # file name fallback
        ("", "5050.mint01", "HD", 0.5, [], "HD50"),  # a ratio alone names no system
        ("unknown liquid", "", None, None, [], None),
        ("", "", None, None, [], None),
    ],
)
def test_guess_sample(
    title: str,
    file_name: str,
    hydrogen: str | None,
    fraction: float | None,
    candidates: list[str],
    sample: str | None,
) -> None:
    guess = guess_sample(title, file_name)
    assert guess["hydrogen"] == hydrogen
    assert guess["hd_fraction"] == (
        pytest.approx(fraction) if fraction is not None else None
    )
    assert guess["candidates"] == candidates
    assert guess["sample"] == sample
    assert isinstance(guess["keywords"], list)


def test_guess_sample_title_beats_file_name() -> None:
    guess = guess_sample("D2O in 4mm can", "SLS18498-H2O.mint01")
    assert guess["hydrogen"] == "D" and guess["sample"] == "D2O"


def test_guess_sample_keywords() -> None:
    assert guess_sample("50/50 H2O/D2O in 2mm can")["keywords"] == [
        "h2o",
        "d2o",
        "50:50",
        "water",
    ]
    assert guess_sample("3mm SiO2 flat plate")["keywords"] == ["silica"]
    assert "dmso" in guess_sample("DMSO(2):2D2O")["keywords"]


# --------------------------------------------------------------------------- synthetic files and errors
def test_synthetic_file_round_trip(tmp_path: Path) -> None:
    path = write_mint01(
        tmp_path / "GEM12345-D2O.mint01", "D2O test sample", source="GEM00012345.raw"
    )
    data = read_mint01(path)
    assert data.instrument == "GEM (ISIS)"
    assert data.run_number == "12345"
    assert data.user == "Tester"
    assert data.date == "01-JAN-2020 12:00:00"
    assert data.n_points == 3  # the 0.4 0 0 row is padding
    assert data.q_range == pytest.approx((0.1, 0.3))
    assert data.q_step == pytest.approx(0.1)
    assert describe(data)["sample"]["sample"] == "D2O"


def test_fortran_exponents_and_blank_lines(tmp_path: Path) -> None:
    path = tmp_path / "odd.mint01"
    path.write_text(
        "# unknown123.dat\n# title\n\n  1.0D-01  -2.0D-01  1.0D-02\n\n  2.0E-01  -1.0E-01  1.0E-02\n"
    )
    data = read_mint01(path)
    assert data.instrument == ""
    assert data.run_number == "123"
    assert data.user == "" and data.date == ""
    assert data.q.tolist() == pytest.approx([0.1, 0.2])
    assert data.f.tolist() == pytest.approx([-0.2, -0.1])


def test_date_found_anywhere_in_header(tmp_path: Path) -> None:
    path = tmp_path / "short.mint01"
    path.write_text("# src.dat\n# title\n# 02-MAR-2021 10:00:00\n  0.1  0.2  0.01\n")
    data = read_mint01(path)
    assert data.date == "02-MAR-2021 10:00:00"
    assert data.user == ""


def test_valid_mask_ignores_nan_rows() -> None:
    data = Mint01(
        path=Path("x.mint01"),
        source="",
        title="",
        header=[],
        q=np.array([0.1, 0.2, 0.3, 0.4]),
        f=np.array([1.0, np.nan, 0.0, 0.0]),
        err=np.array([0.1, 0.1, 0.1, 0.0]),
    )
    assert data.valid.tolist() == [True, False, True, False]
    assert data.n_points == 2
    assert data.instrument == "" and data.run_number == ""


def test_instrument_prefixes() -> None:
    def source(name: str) -> Mint01:
        return Mint01(Path("x"), name, "", [], np.zeros(0), np.zeros(0), np.zeros(0))

    assert source("POLARIS00099.raw").instrument == "POLARIS (ISIS)"
    assert source("POLARIS00099.raw").run_number == "99"
    assert source("nimrod00001122.mint01").instrument == "NIMROD (ISIS)"
    assert source("D4_000123.dat").instrument == "D4 (ILL)"
    assert source("D4_000123.dat").run_number == "123"
    assert source("mystery").run_number == ""


def test_read_errors(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="cannot read"):
        read_mint01(tmp_path / "missing.mint01")
    header_only = tmp_path / "header.mint01"
    header_only.write_text("# src\n# title\n")
    with pytest.raises(ValueError, match="no three-column data"):
        read_mint01(header_only)
    columns = tmp_path / "columns.mint01"
    columns.write_text("# src\n# title\n  0.1  0.2\n")
    with pytest.raises(ValueError, match="line 3.*expected 3 columns"):
        read_mint01(columns)
    text = tmp_path / "text.mint01"
    text.write_text("# src\n# title\n  0.1  abc  0.3\n")
    with pytest.raises(ValueError, match="line 3.*not a numeric"):
        read_mint01(text)
    zeros = tmp_path / "zeros.mint01"
    zeros.write_text("# src\n# title\n  0.1  0.0  0.0\n  0.2  0.0  0.0\n")
    with pytest.raises(ValueError, match="zero value and zero error"):
        read_mint01(zeros)


# --------------------------------------------------------------------------- figure
def test_plot_data_writes_png(tmp_path: Path) -> None:
    datasets = [read_mint01(path) for path in WATER_FILES]
    out = plot_data(datasets, tmp_path / "figures" / "00_data.png")
    assert out == tmp_path / "figures" / "00_data.png"
    assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert out.stat().st_size > 10_000
    single = plot_data(datasets[:1], tmp_path / "single.png", offset=1.0)
    assert single.read_bytes()[:4] == b"\x89PNG"


def test_plot_data_many_series_and_empty(tmp_path: Path) -> None:
    base = read_mint01(WATER_FILES[0])
    many = [base] * 10
    assert plot_data(many, tmp_path / "many.png").exists()
    with pytest.raises(ValueError, match="no data sets"):
        plot_data([], tmp_path / "none.png")
