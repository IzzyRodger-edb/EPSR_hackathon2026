"""Pins what the app assumes about Mellor's API.

The adapter is written against `mellor/docs/integration.md`. The doc and the code
can drift, and a mismatch would otherwise show up as a run that fails ten minutes
in. These tests read the real signatures instead, so a rename or a changed return
type fails here, in a second, with the name of what moved.

Skipped where Mellor is not importable (it needs the compiled EPSR engine), which
is also the case the app's demonstration mode exists for.
"""

from __future__ import annotations

import dataclasses
import inspect

import pytest

from epsrweb import mellor_adapter, plots, storage

mellor = pytest.importorskip("mellor", reason="mellor (and the EPSR engine) is not installed")


def _params(function) -> set[str]:
    return set(inspect.signature(function).parameters)


# ------------------------------------------------------------------ the agent


def test_the_app_agrees_that_mellor_is_available():
    status = mellor_adapter.mellor_status()
    assert status["available"] is True, status["reason"]


def test_session_takes_the_fields_the_worker_passes():
    from mellor.llm.tools import Session

    fields = {field.name for field in dataclasses.fields(Session)}
    assert {"work_dir", "data_files", "quick"} <= fields
    # Read back by the results page and by chat on a later request.
    assert {"analysis", "figures", "report_path", "spec_source"} <= fields
    assert callable(Session.load_existing)


def test_agent_takes_on_event_and_run_pipeline_takes_a_prompt():
    from mellor.llm.agent import Mellor

    assert "on_event" in _params(Mellor.__init__)
    assert "user_prompt" in _params(Mellor.run_pipeline)
    assert "user_message" in _params(Mellor.chat)


def test_pipeline_outcome_fields_the_app_stores():
    from mellor.llm.agent import PipelineOutcome

    fields = {field.name for field in dataclasses.fields(PipelineOutcome)}
    assert {"report_path", "final_message", "spec_source", "steps"} <= fields


def test_the_no_model_constants_exist():
    from mellor.llm import agent

    assert agent.DISCUSSION_NAME == "discussion.md"
    assert agent.TRANSCRIPT_NAME == "mellor_transcript.jsonl"
    assert agent.NO_MODEL and agent.NO_MODEL_NOTE
    # Both are offered as downloads by name.
    assert agent.DISCUSSION_NAME in __import__("epsrweb.runviews", fromlist=["x"]).SERVABLE_FILES
    assert agent.TRANSCRIPT_NAME in __import__("epsrweb.runviews", fromlist=["x"]).SERVABLE_FILES


def test_make_backend_accepts_the_strings_the_form_offers():
    from mellor.llm.backend import make_backend

    assert "spec" in _params(make_backend)
    offered = {value for value, _label, _token in mellor_adapter.BACKENDS}
    # `library` is the app's own name for "no model at all"; the rest go to Mellor.
    for value in offered - {"library"}:
        assert mellor_adapter.backend_string(value, "a/model").startswith(value)


# ----------------------------------------------------------------- the library


def test_known_systems_matches_the_fallback_list():
    from mellor import library

    assert set(mellor_adapter.FALLBACK_SYSTEMS) <= set(library.KNOWN_SYSTEMS)
    assert mellor_adapter.known_systems() == tuple(library.KNOWN_SYSTEMS)


def test_auto_spec_returns_a_spec_and_a_reason():
    """The doc implies a bare spec; the code returns `(spec, reason)`."""
    from mellor import library

    assert _params(library.auto_spec) >= {"datasets", "name"}
    annotation = inspect.signature(library.auto_spec).return_annotation
    assert "tuple" in str(annotation), annotation


def test_make_spec_takes_a_key_and_the_datasets():
    from mellor import library

    assert list(_params(library.make_spec))[:2] == ["key", "datasets"]


# ------------------------------------------------------------------ the runner


def test_run_spec_can_overwrite_a_folder_that_already_holds_app_files():
    """The run folder is never empty: it holds the uploads and the event log."""
    from mellor import runner

    assert {"quick", "record", "overwrite", "progress"} <= _params(runner.run_spec)
    # overwrite must remove only the products of an earlier run.
    for product in ("epsr", "trajectory", "figures", "tables", "analysis.json"):
        assert product in runner.RUN_PRODUCTS
    for kept in ("uploads", "events.jsonl", "result.json", plots.WEB_FIGURES_DIR):
        assert kept not in runner.RUN_PRODUCTS


# -------------------------------------------------------- analysis and output


def test_analysis_entry_points_and_keys_the_results_page_reads():
    from mellor import analysis

    assert _params(analysis.analyse) == {"result"}
    assert list(_params(analysis.save_analysis))[:2] == ["analysis", "path"]
    assert _params(analysis.analysis_markdown) == {"analysis"}
    # The keys `runviews.analysis_view` looks for are the ones `analyse` documents.
    documented = analysis.analyse.__doc__ or ""
    for key in ("runs", "fit", "convergence", "shells", "flags"):
        assert f"``{key}``" in documented, key


def test_figure_captions_cover_the_names_the_app_knows():
    from mellor import figures

    assert _params(figures.make_figures) >= {"result", "out_dir", "dpi"}
    for name in storage.FALLBACK_CAPTIONS:
        assert name in figures.FIGURE_CAPTIONS, name
    # The app merges Mellor's captions over its own, so every figure gets one.
    merged = mellor_adapter.figure_captions()
    assert set(figures.FIGURE_CAPTIONS) <= set(merged)


def test_build_pdf_takes_the_arguments_the_no_model_path_passes():
    from mellor import report

    names = list(_params(report.build_pdf))
    assert names[:4] == ["result", "analysis", "figures", "discussion_markdown"]
    assert "model_name" in names


def test_the_seaborn_figures_read_tables_that_export_writes():
    from mellor import export

    assert "folder" in _params(export.write_tables)
    needed = {"fits", "partials", "potentials", "history", "runs"}
    assert needed <= set(export.TABLES)


def test_read_mint01_mirrors_the_apps_validator(water_files):
    from mellor import data

    for path in water_files:
        theirs = data.read_mint01(path)
        ours = storage.validate_mint01(path.read_text(encoding="utf-8", errors="replace"),
                                       name=path.name)
        assert ours.title == theirs.title, (ours.title, theirs.title)
