"""``mellor.report``: the PDF report built from a hand-made run folder.

The fixture mirrors what ``runner.run_spec`` leaves behind (docs/architecture.md §3)
without running EPSR: the library's water specification for the three SANDALS files,
the generated script, the toolkit's example EPSR output folder as ``epsr/``, a
manifest as in §6, an analysis dict with the keys of ``analysis.analyse`` and a few
PNG figures named as in ``figures.FIGURE_CAPTIONS``.
"""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path
from typing import Any

import pytest
from pypdf import PdfReader
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import HRFlowable, Paragraph, Preformatted, Table
from specs import REPOSITORY, WATER_DATA

from mellor import data, inputs, library, report, script
from mellor.figures import FIGURE_CAPTIONS
from mellor.runner import RunResult

EXAMPLE_EPSR = (
    REPOSITORY
    / "visualisation"
    / "python_visualisation_toolkit"
    / "example_H2O"
    / "data"
    / "epsr_output"
)
WATER_FILES = [
    WATER_DATA / name
    for name in (
        "SLS18502-D2O.mint01",
        "SLS18500-HDO5050.mint01",
        "SLS18498-H2O.mint01",
    )
]
MODEL = "Qwen/Qwen3-4B-Instruct-2507"

DISCUSSION = """# Water at 25 °C

The refinement reproduces the three SANDALS data sets well: the final R-factor of
2.27e-4 is **two orders of magnitude** below the scratch start, and the *empirical
potential* reached the requested `ereq` of 3 kJ/mol.

## Structure

- The first OW–OW peak sits at 2.78 Å with 4.6 neighbours within 3.35 Å, the
  tetrahedral hydrogen-bonded network of liquid water.
- The OW–HW peak at 1.85 Å is the hydrogen bond.
  * Its coordination number of 1.9 means almost every hydrogen is bonded.
- Nothing above 10 Å⁻¹ in F(Q) is left unfitted.

1. The pressure of 5.4 kbar is typical for EPSR water.
2. The energy is flat over the last run.

> Every number above was computed by EPSR.

| Pair | r (Å) | n |
|:-----|------:|--:|
| OW-OW | 2.78 | 4.6 |
| OW-HW | 1.85 | 1.9 |

```python
print("verbatim & <code>")
```

---

See [Soper 2013](https://doi.org/10.1155/2013/279463) and ***bold italic*** text,
then a lone asterisk a * b * c that is not emphasis.
"""


# --------------------------------------------------------------------------- fixture
def _png(path: Path, title: str) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    figure = Figure(figsize=(6.0, 3.0), dpi=72)
    axes = figure.add_subplot(111)
    axes.plot([0.0, 1.0, 2.0], [0.0, 1.0, 0.0])
    axes.set_title(title)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path)
    return path


def _manifest(spec: Any) -> dict[str, Any]:
    """A manifest as the generated script writes it (docs/architecture.md §6)."""
    r_factors = [2.12e-2, 6.1e-4, 2.9e-4, 2.27e-4]
    per_dataset = [
        [2.4e-2, 1.9e-2, 2.0e-2],
        [7.0e-4, 5.5e-4, 5.9e-4],
        [3.1e-4, 2.7e-4, 2.9e-4],
        [2.4e-4, 2.1e-4, 2.3e-4],
    ]
    energies = [-12.3, -40.2, -41.1, -41.3]
    pressures = [1.2, 5.9, 5.5, 5.4]
    amplitudes = [0.0, 0.0, 2.8, 2.97]
    durations = [1.3, 99.0, 100.5, 101.2]
    runs = []
    for index, (stage, niter, ereq) in enumerate(spec.schedule.iterations()):
        runs.append(
            {
                "run": index + 1,
                "stage": stage,
                "niter": niter,
                "ereq": ereq,
                "overrides": [],
                "duration_s": durations[index],
                "summary": {
                    "r_factor": r_factors[index],
                    "r_factor_per_dataset": per_dataset[index],
                    "energy": energies[index],
                    "energy_inter": energies[index] - 5.0,
                    "energy_intra": 5.0,
                    "pressure": pressures[index],
                    "ep_amplitude": amplitudes[index],
                    "chisq": 1.4,
                },
            }
        )
    generated = inputs.generate(spec)
    return {
        "name": spec.name,
        "title": spec.title,
        "kind": spec.kind,
        "quick": False,
        "record": True,
        "complete": True,
        "created": "2026-10-17T14:02:11",
        "duration_s": sum(durations),
        "python": "3.12.3",
        "versions": {"epsr": "0.1.0", "epsrvis": "0.1.0", "mellor": "0.1.0"},
        "script": f"run_{spec.name}.py",
        "epsr_folder": "epsr",
        "trajectory": "trajectory",
        "inp_name": spec.inp_name,
        "epsrvis_case": "water",
        "temperature": 300.0,
        "box_side": 31.0723,
        "n_atoms": 3000,
        "n_molecules": 1000,
        "density": 0.1,
        "atom_types": ["OW", "HW"],
        "pairs": ["OW-OW", "OW-HW", "HW-HW"],
        "datasets": [
            {"label": d.label, "file": d.file_name, "wts": d.weights_file}
            for d in spec.datasets
        ],
        "setup_steps": generated.setup_steps,
        "schedule": [
            {"stage": stage, "niter": niter, "ereq": ereq}
            for stage, niter, ereq in spec.schedule.iterations()
        ],
        "runs": runs,
        "final": runs[-1]["summary"],
    }


def _analysis(result: RunResult) -> dict[str, Any]:
    """An analysis dict with the keys documented in ``analysis.analyse``."""
    spec, manifest = result.spec, result.manifest
    labels = [d.label for d in spec.datasets]
    runs = []
    for run in manifest["runs"]:
        summary = run["summary"]
        runs.append(
            {
                "run": run["run"],
                "stage": run["stage"],
                "niter": run["niter"],
                "ereq": run["ereq"],
                "r_factor": summary["r_factor"],
                "r_factor_per_dataset": dict(
                    zip(labels, summary["r_factor_per_dataset"], strict=True)
                ),
                "energy": summary["energy"],
                "energy_inter": summary["energy_inter"],
                "energy_intra": summary["energy_intra"],
                "pressure": summary["pressure"],
                "ep_amplitude": summary["ep_amplitude"],
                "chisq": summary["chisq"],
                "duration_s": run["duration_s"],
                "r2": {label: 0.99 for label in labels},
                "chi2": {label: 1.2 for label in labels},
            }
        )
    first, reference, last = runs[0], runs[1], runs[-1]
    fit = {
        label: {
            "file": dataset.file_name,
            "initial_r_factor": first["r_factor_per_dataset"][label],
            "reference_r_factor": reference["r_factor_per_dataset"][label],
            "final_r_factor": last["r_factor_per_dataset"][label],
            "improvement_factor": first["r_factor_per_dataset"][label]
            / last["r_factor_per_dataset"][label],
            "r2": 0.998,
            "chi2": 1.1,
            "q_min": 0.175,
            "q_max": 49.925,
            "n_points": 996,
        }
        for label, dataset in zip(labels, spec.datasets, strict=True)
    }
    shells = {
        "OW-OW": {
            "centre": "OW",
            "neighbour": "OW",
            "r_peak": 2.78,
            "g_peak": 3.1,
            "r_min": 3.35,
            "n_coord": 4.6,
            "n_coord_reverse": None,
            "partner_density": 0.0333,
            "description": "first peak at 2.78 Å (g = 3.10), minimum at 3.35 Å",
        },
        "OW-HW": {
            "centre": "OW",
            "neighbour": "HW",
            "r_peak": 1.85,
            "g_peak": 1.4,
            "r_min": 2.45,
            "n_coord": 1.9,
            "n_coord_reverse": 0.95,
            "partner_density": 0.0667,
            "description": "first peak at 1.85 Å (g = 1.40), minimum at 2.45 Å",
        },
    }
    return {
        "system": {
            "name": spec.name,
            "title": spec.title,
            "kind": spec.kind,
            "source": spec.source,
            "n_atoms": 3000,
            "n_molecules": 1000,
            "box_side": 31.0723,
            "density": 0.1,
            "mass_density_g_cm3": 1.107,
            "temperature": 300.0,
            "atom_types": [
                {"label": "OW", "count": 1000},
                {"label": "HW", "count": 2000},
            ],
            "pairs": ["OW-OW", "OW-HW", "HW-HW"],
            "datasets": [
                {"label": d.label, "file": d.file_name} for d in spec.datasets
            ],
        },
        "schedule": {
            "runs": manifest["schedule"],
            "total_iterations": 300,
            "quick": False,
            "ereq": 3.0,
            "ntimes": 5,
            "complete": True,
            "recorded": True,
            "setup_steps": manifest["setup_steps"],
        },
        "runs": runs,
        "fit": fit,
        "convergence": {
            "last_run": 4,
            "last_run_stage": "EP refinement",
            "ereq": 3.0,
            "ep_amplitude": 2.97,
            "ep_fraction_of_ereq": 0.99,
            "ep_at_limit": True,
            "energy": -41.3,
            "pressure": 5.4,
            "r_factor": 2.27e-4,
            "r_factor_initial": 2.12e-2,
            "r_factor_before_last_run": 2.9e-4,
            "r_factor_change_percent": -21.7,
            "trends": {
                "energy": {
                    "start": -41.2,
                    "end": -41.3,
                    "change": -0.1,
                    "change_percent": -0.2,
                },
                "pressure": {
                    "start": 5.6,
                    "end": 5.4,
                    "change": -0.2,
                    "change_percent": -3.6,
                },
                "ep_amplitude": {
                    "start": 2.9,
                    "end": 2.97,
                    "change": 0.07,
                    "change_percent": 2.4,
                },
                "r_factor": {
                    "start": 2.5e-4,
                    "end": 2.27e-4,
                    "change": -2.3e-5,
                    "change_percent": -9.2,
                },
            },
            "r_factor_rising": False,
            "energy_drift_percent": -0.2,
            "last_run_iterations": 100,
            "quarter_iterations": 25,
            "ep_limit_reached_after": 12,
            "verdict": "the EP amplitude reached the ereq limit of 3 kJ/mol after 12 iterations; "
            "the R-factor fell from 2.9e-04 to 2.3e-04 over the last run (-22 %).",
        },
        "shells": shells,
        "flags": [
            "the EP amplitude is at the ereq limit (2.97 of 3 kJ/mol): the fit may improve "
            "with a larger ereq",
        ],
        "figures": {},
    }


@pytest.fixture(scope="module")
def run_folder(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    if not EXAMPLE_EPSR.is_dir():
        pytest.skip("toolkit example data not found")
    folder = tmp_path_factory.mktemp("run")
    datasets = [data.read_mint01(path) for path in WATER_FILES]
    spec = library.make_spec("water", datasets)
    spec.validate()
    (folder / "data").mkdir()
    for path in WATER_FILES:
        shutil.copy(path, folder / "data" / path.name)
    spec.save(folder / "spec.json")
    generated = inputs.generate(spec)
    script.write_script(spec, generated, folder)
    shutil.copytree(EXAMPLE_EPSR, folder / "epsr")
    result = RunResult(
        folder=folder, spec=spec, manifest=_manifest(spec), duration_s=302.0
    )
    figures_dir = folder / "figures"
    figures = {
        "00_data": _png(figures_dir / "00_data.png", "data"),
        "01_convergence": _png(figures_dir / "01_convergence.png", "convergence"),
        "11_fit_final": _png(figures_dir / "11_fit_final.png", "fit"),
        "05_partials": figures_dir / "05_partials.png",  # never written: skipped
        "summary.csv": figures_dir / "summary.csv",  # not an image: skipped
        "09_displacement": figures_dir / "09_displacement.png",  # unreadable: skipped
    }
    (figures_dir / "summary.csv").write_text("pair,r_peak\nOW-OW,2.78\n")
    (figures_dir / "09_displacement.png").write_text("not a png")
    return {
        "result": result,
        "analysis": _analysis(result),
        "figures": figures,
        "descriptions": data.describe_all(WATER_FILES),
    }


def _text(path: Path) -> tuple[int, str]:
    reader = PdfReader(str(path))
    return len(reader.pages), "\n".join(page.extract_text() for page in reader.pages)


# --------------------------------------------------------------------------- build_pdf
def test_build_pdf_full(run_folder: dict[str, Any], tmp_path: Path) -> None:
    result: RunResult = run_folder["result"]
    out = tmp_path / "report.pdf"
    path = report.build_pdf(
        result,
        run_folder["analysis"],
        run_folder["figures"],
        DISCUSSION,
        out,
        model_name=MODEL,
        data_descriptions=run_folder["descriptions"],
    )
    assert path == out and path.is_file()
    n_pages, text = _text(path)
    assert n_pages >= 6
    # title page
    assert result.spec.title in text
    assert "Mellor" in text and MODEL in text
    # every section
    for heading in report.SECTION_HEADINGS:
        assert heading in text, heading
    # 1. the data (from the descriptions)
    assert "SANDALS" in text and "SLS18502-D2O.mint01" in text
    # 2. the model
    assert "OW" in text and "3.1655" in text  # σ of OW as in the specification
    # 3. the calculation
    assert "makemole Water" in text and "reference potential" in text
    # 4. results: run table value, fit table, shells, flags, figure captions
    assert "2.27e-4" in text and "2.12e-2" in text
    assert "Fit quality per data set" in text and "First coordination shells" in text
    assert "2.78" in text and "4.60" in text
    assert "at the ereq limit" in text
    assert FIGURE_CAPTIONS["11_fit_final"] in text
    assert FIGURE_CAPTIONS["01_convergence"].split(":")[0] in text
    assert FIGURE_CAPTIONS["05_partials"] not in text  # missing file: skipped
    assert FIGURE_CAPTIONS["09_displacement"] not in text  # unreadable: skipped
    # 5. discussion
    assert "Written by the language model" in text
    assert "hydrogen-bonded network" in text and "verbatim & <code>" in text
    # 6. reproducibility: command, layout, hashes
    assert f"python run_{result.spec.name}.py --out" in text
    assert "manifest.json" in text
    digest = hashlib.sha256(WATER_FILES[0].read_bytes()).hexdigest()
    assert digest in text
    # appendices
    assert "Water.mol" in text and result.spec.inp_name in text
    assert "class RecorderAdapter" in text and "def main(" in text
    # footer
    assert "page 1" in text and f"page {n_pages}" in text


def test_build_pdf_minimal(run_folder: dict[str, Any], tmp_path: Path) -> None:
    """An empty analysis, no figures, no discussion, no descriptions: still a report."""
    result: RunResult = run_folder["result"]
    path = report.build_pdf(result, {}, {}, "", tmp_path / "minimal.pdf")
    n_pages, text = _text(path)
    assert n_pages >= 3
    for heading in report.SECTION_HEADINGS:
        assert heading in text, heading
    assert "SANDALS" in text  # the data rows come from the files in <run>/data
    assert "No discussion was written" in text
    assert "model not named" in text
    assert "No figures are available" in text
    # the manifest's runs are still tabulated without an analysis dict
    assert "2.27e-4" in text and "EP refinement" in text


def test_build_pdf_default_path_and_escaping(
    run_folder: dict[str, Any], tmp_path: Path
) -> None:
    """``out_path`` None means <run>/report.pdf; ``&``, ``<`` and ``>`` survive."""
    original: RunResult = run_folder["result"]
    folder = tmp_path / "copy"
    shutil.copytree(original.folder, folder)
    spec = library.make_spec("water", [data.read_mint01(p) for p in WATER_FILES])
    spec.title = "Water: H2O & D2O <mixtures> at 25 °C"
    spec.notes = "Chosen because <title> says water & the files are SANDALS runs."
    result = RunResult(folder=folder, spec=spec, manifest=original.manifest)
    path = report.build_pdf(result, run_folder["analysis"], {}, "A & B < C > D")
    assert path == folder / "report.pdf"
    _n_pages, text = _text(path)
    assert "H2O & D2O <mixtures>" in text
    assert "<title> says water &" in text
    assert "A & B < C > D" in text


def test_build_pdf_without_script_or_data(tmp_path: Path) -> None:
    """A bare folder (no script, no data copies, no epsr folder): the report says so."""
    spec = library.make_spec("water", [data.read_mint01(p) for p in WATER_FILES])
    folder = tmp_path / "bare"
    folder.mkdir()
    result = RunResult(folder=folder, spec=spec, manifest={})
    path = report.build_pdf(result, {}, {}, "", tmp_path / "bare.pdf")
    _n_pages, text = _text(path)
    assert "No EPSR runs are recorded" in text
    assert "was not found in the run folder" in text  # Appendix B, regenerated script
    assert "class RecorderAdapter" in text
    assert "file not found" in text  # the hash table of spec.json / the script


# --------------------------------------------------------------------------- markdown
def _plain_text(flowables: list[Any]) -> str:
    parts = []
    for flowable in flowables:
        if isinstance(flowable, Paragraph):
            parts.append(flowable.getPlainText())
        elif isinstance(flowable, Preformatted):
            parts.append("\n".join(flowable.lines))
    return "\n".join(parts)


def test_markdown_to_flowables_every_construct() -> None:
    styles = report.make_styles()
    flowables = report.markdown_to_flowables(DISCUSSION, styles)
    assert len(flowables) > 0
    kinds = {type(flowable) for flowable in flowables}
    assert {Paragraph, Table, Preformatted, HRFlowable} <= kinds
    text = _plain_text(flowables)
    assert "Å" in text and "⁻¹" in text and "°C" in text
    paragraphs = [f for f in flowables if isinstance(f, Paragraph)]
    markup = "\n".join(p.text for p in paragraphs)
    assert "<b>two orders of magnitude</b>" in markup
    assert "<i>empirical potential</i>" in markup
    assert "<b><i>bold italic</i></b>" in markup
    assert 'face="DejaVuSansMono">ereq</font>' in markup
    assert '<a href="https://doi.org/10.1155/2013/279463"' in markup
    assert "a * b * c" in text  # a lone asterisk is not emphasis
    bullets = [p.bulletText for p in paragraphs if getattr(p, "bulletText", None)]
    assert "•" in bullets and "1." in bullets and "2." in bullets
    headings = [p for p in paragraphs if p.style.name in ("MdH1", "MdH2")]
    assert [h.getPlainText() for h in headings] == ["Water at 25 °C", "Structure"]
    code = [f for f in flowables if isinstance(f, Preformatted)]
    assert code and 'print("verbatim & <code>")' in code[0].lines[0]
    table = next(f for f in flowables if isinstance(f, Table))
    assert len(table._cellvalues) == 3  # header + two rows


def test_markdown_to_flowables_other_stylesheets() -> None:
    sample = (
        "Plain **bold** line with Å⁻¹.\n\n- item\n\n| a | b |\n|---|---|\n| 1 | 2 |\n"
    )
    assert report.markdown_to_flowables(sample, None)
    assert report.markdown_to_flowables(sample, getSampleStyleSheet())
    assert report.markdown_to_flowables("", report.make_styles()) == []


def test_markdown_unknown_syntax_falls_back_to_paragraphs() -> None:
    text = "<div>html block</div>\n\nterm\n: definition\n\n#hashtag not a heading"
    flowables = report.markdown_to_flowables(text, report.make_styles())
    assert all(isinstance(flowable, Paragraph) for flowable in flowables)
    assert "<div>html block</div>" in _plain_text(flowables)
    assert "#hashtag not a heading" in _plain_text(flowables)


# --------------------------------------------------------------------------- helpers
def test_number_formats() -> None:
    assert report.format_r_factor(2.27e-4) == "2.27e-4"
    assert report.format_r_factor(0.0212) == "2.12e-2"
    assert report.format_r_factor(1.5) == "1.50e0"
    assert report.format_r_factor(0) == "0"
    assert report.format_r_factor(None) == "–"
    assert report.format_r_factor("nan") == "–"
    assert report.format_fixed(-41.26, 1) == "-41.3"
    assert report.format_fixed(2.7849, 2) == "2.78"
    assert report.format_fixed(True, 2) == "–"
    assert report.format_general(3.165492) == "3.1655"
    assert report.format_general(-0.82) == "-0.82"
    assert report.format_int(996.0) == "996"
    assert report.escape_text("a & <b>") == "a &amp; &lt;b&gt;"
    assert report.escape_text(None) == ""


def test_register_fonts_uses_dejavu() -> None:
    fonts = report.register_fonts()
    assert fonts["sans"] == "DejaVuSans"
    assert fonts["bold"] == "DejaVuSans-Bold"
    assert fonts["mono"] == "DejaVuSansMono"
