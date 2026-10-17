"""The same contract, checked by reading Mellor's source instead of importing it.

`import mellor` pulls in the compiled EPSR engine, so
`test_mellor_contract.py` can only run where that has been built. These checks
parse the sibling package's source with `ast`, so the assumptions the adapter was
written against are verified on any machine -- including one that cannot build
the engine, which is where a silent drift would otherwise go unnoticed longest.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from epsrweb import mellor_adapter, plots, runviews, storage

MELLOR_SRC = Path(__file__).resolve().parents[2] / "mellor" / "src" / "mellor"

pytestmark = pytest.mark.skipif(
    not MELLOR_SRC.is_dir(), reason="the mellor source tree is not in this checkout"
)


def _module(relative: str) -> ast.Module:
    return ast.parse((MELLOR_SRC / relative).read_text(encoding="utf-8"))


def _function(relative: str, name: str) -> ast.FunctionDef:
    for node in ast.walk(_module(relative)):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
    raise AssertionError(f"{relative} has no function {name}")


def _parameters(relative: str, name: str) -> list[str]:
    node = _function(relative, name)
    args = node.args
    return [a.arg for a in (*args.posonlyargs, *args.args, *args.kwonlyargs)]


def _class_fields(relative: str, name: str) -> list[str]:
    for node in ast.walk(_module(relative)):
        if isinstance(node, ast.ClassDef) and node.name == name:
            return [
                item.target.id
                for item in node.body
                if isinstance(item, ast.AnnAssign) and isinstance(item.target, ast.Name)
            ]
    raise AssertionError(f"{relative} has no class {name}")


def _constant(relative: str, name: str) -> object:
    for node in _module(relative).body:
        targets = (
            [node.target] if isinstance(node, ast.AnnAssign)
            else getattr(node, "targets", [])
        )
        for target in targets:
            if isinstance(target, ast.Name) and target.id == name:
                return ast.literal_eval(node.value)
    raise AssertionError(f"{relative} has no constant {name}")


def _dict_keys(relative: str, name: str) -> list[str]:
    value = _constant(relative, name)
    assert isinstance(value, dict)
    return list(value)


# ------------------------------------------------------------------ the agent


def test_session_fields():
    fields = _class_fields("llm/tools.py", "Session")
    assert {"work_dir", "data_files", "quick"} <= set(fields)
    assert {"analysis", "figures", "report_path", "spec_source"} <= set(fields)


def test_pipeline_outcome_fields():
    assert {"report_path", "final_message", "spec_source", "steps"} <= set(
        _class_fields("llm/agent.py", "PipelineOutcome")
    )


def test_run_pipeline_and_chat_parameters():
    assert "user_prompt" in _parameters("llm/agent.py", "run_pipeline")
    assert "on_event" in _parameters("llm/agent.py", "__init__")
    assert "user_message" in _parameters("llm/agent.py", "chat")


def test_file_name_constants():
    assert _constant("llm/agent.py", "DISCUSSION_NAME") == "discussion.md"
    assert _constant("llm/agent.py", "TRANSCRIPT_NAME") == "mellor_transcript.jsonl"
    for name in ("discussion.md", "mellor_transcript.jsonl", "run.log"):
        assert name in runviews.SERVABLE_FILES


# --------------------------------------------- the three integration surprises


def test_auto_spec_returns_a_tuple_not_a_bare_spec():
    node = _function("library.py", "auto_spec")
    returns = ast.unparse(node.returns) if node.returns else ""
    assert returns.startswith("tuple"), returns


def test_run_spec_accepts_overwrite():
    """The run folder already holds the uploads and the event log when EPSR starts."""
    assert {"quick", "record", "overwrite", "progress"} <= set(
        _parameters("runner.py", "run_spec")
    )
    products = _constant("runner.py", "RUN_PRODUCTS")
    assert {"epsr", "trajectory", "figures", "tables", "analysis.json"} <= set(products)
    for kept in ("uploads", "events.jsonl", "result.json", plots.WEB_FIGURES_DIR):
        assert kept not in products


def test_build_pdf_argument_order():
    names = _parameters("report.py", "build_pdf")
    assert names[:4] == ["result", "analysis", "figures", "discussion_markdown"]
    assert "model_name" in names


# ------------------------------------------------------------- names and keys


def test_known_systems_cover_the_fallback_list():
    """`KNOWN_SYSTEMS` is built by comprehension, so collect the `key=` of each system."""
    keys = {
        keyword.value.value
        for node in ast.walk(_module("library.py"))
        if isinstance(node, ast.Call)
        for keyword in node.keywords
        if keyword.arg == "key"
        and isinstance(keyword.value, ast.Constant)
        and isinstance(keyword.value.value, str)
    }
    assert set(mellor_adapter.FALLBACK_SYSTEMS) <= keys, sorted(keys)


def test_every_figure_name_the_app_knows_still_exists():
    captions = _dict_keys("figures.py", "FIGURE_CAPTIONS")
    for name in storage.FALLBACK_CAPTIONS:
        assert name in captions, name


def test_the_tables_the_seaborn_figures_read():
    assert {"fits", "partials", "potentials", "history", "runs"} <= set(
        _dict_keys("export.py", "TABLES")
    )


def test_the_event_types_the_log_pane_renders():
    """Every `self._emit({"type": ...})` in the agent is a type `summarise` names."""
    from epsrweb import events

    emitted = {
        node.value.value
        for node in ast.walk(_module("llm/agent.py"))
        if isinstance(node, ast.Dict)
        for key, node.value in zip(node.keys, node.values)
        if isinstance(key, ast.Constant) and key.value == "type"
        and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str)
    }
    known = set(events.MELLOR_TYPES) | set(events.APP_TYPES)
    assert emitted <= known, f"unhandled event types: {sorted(emitted - known)}"


def test_analysis_documents_the_keys_the_results_page_reads():
    doc = ast.get_docstring(_function("analysis.py", "analyse")) or ""
    for key in ("runs", "fit", "convergence", "shells", "flags"):
        assert f"``{key}``" in doc, key
