"""The whole pipeline, end to end: a scripted backend plays the model and EPSR runs
the library's water system with the quick schedule (31 iterations, about a minute).

Nothing is faked: data, library, inputs, script, runner, analysis, figures and
report all do their real work on ``experimental_data/water``. The model's script
follows the setup procedure (inspect_data, library_spec, build_and_run), writes the
discussion and then answers one question through ``get_results``.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from example_run import WATER_FILES
from pypdf import PdfReader

from mellor import analysis as analysis_module
from mellor import report as report_module
from mellor.llm.agent import DISCUSSION_NAME, TRANSCRIPT_NAME, Mellor, PipelineOutcome
from mellor.llm.backend import AssistantMessage, ScriptedBackend, ToolCall
from mellor.llm.tools import Session
from mellor.spec import SystemSpec

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        not all(path.is_file() for path in WATER_FILES),
        reason="the water .mint01 files are not in the repository",
    ),
]

NARRATIVE = """## Discussion

The quick refinement of liquid water against the three SANDALS data sets ran to the
end of its schedule. The empirical potential reached the requested amplitude and the
R-factor fell during the refinement runs; the tetrahedral hydrogen-bonded network of
water is visible in the partial g(r). A quick schedule is a preview, not a converged
refinement.
"""
QUESTION = "What is the final R-factor?"


def scripted_model(messages: list[dict]) -> AssistantMessage:
    """The model's script, keyed on the last message of the conversation."""
    last = messages[-1]
    if last["role"] == "tool":
        content = json.loads(last["content"])
        name = last["name"]
        if name == "inspect_data":
            return AssistantMessage(
                "Three water samples; the library water system fits.",
                [ToolCall("library_spec", {"system": "water"})],
            )
        if name == "library_spec":
            # the spec of the previous tool result goes straight to build_and_run
            return AssistantMessage(
                "", [ToolCall("build_and_run", {"spec": content["spec"]})]
            )
        if name == "build_and_run":
            return AssistantMessage(
                "The run finished." if content.get("ok") else "The run failed."
            )
        if name == "get_results":
            r_factor = content["convergence"]["r_factor"]
            return AssistantMessage(f"The final R-factor is {r_factor:.3e}.")
    if last["role"] == "user":
        text = last["content"]
        if text.startswith("Set up and run an EPSR refinement"):
            return AssistantMessage(
                "Let me inspect the data.", [ToolCall("inspect_data", {})]
            )
        if 'Write the "Discussion" section' in text:
            return AssistantMessage(NARRATIVE)
        if f"User: {QUESTION}" in text:
            return AssistantMessage(
                "", [ToolCall("get_results", {"section": "convergence"})]
            )
    return AssistantMessage("I do not know.")


def _read_transcript(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _pdf_text(path: Path) -> str:
    return "\n".join(page.extract_text() for page in PdfReader(str(path)).pages)


# --------------------------------------------------------------------------- the run
@pytest.fixture(scope="module")
def pipeline(tmp_path_factory: pytest.TempPathFactory) -> SimpleNamespace:
    work_dir = tmp_path_factory.mktemp("pipeline") / "water"
    session = Session(work_dir=work_dir, data_files=list(WATER_FILES), quick=True)
    backend = ScriptedBackend(scripted_model)
    events: list[dict] = []
    agent = Mellor(backend, session, on_event=events.append)
    outcome = agent.run_pipeline("Liquid water at 25 C, three isotopic samples.")
    return SimpleNamespace(
        session=session, agent=agent, backend=backend, outcome=outcome, events=events
    )


def test_outcome_and_session(pipeline: SimpleNamespace) -> None:
    outcome, session = pipeline.outcome, pipeline.session
    assert isinstance(outcome, PipelineOutcome)
    assert (
        outcome.spec_source == "library:water"
        and session.spec_source == "library:water"
    )
    assert outcome.steps == 4  # inspect_data, library_spec, build_and_run, narrative
    assert pipeline.agent.spec_failures == 0 and pipeline.agent.notes == []
    assert outcome.report_path == session.work_dir / "report.pdf"
    assert (
        "Report:" in outcome.final_message and "quick schedule" in outcome.final_message
    )

    result = session.result
    assert result is not None and result.quick and result.complete
    assert result.folder == session.work_dir
    assert [run["niter"] for run in result.runs] == [1, 10, 10, 10]
    assert [run["ereq"] for run in result.runs] == [0.0, 0.0, 3.0, 3.0]
    assert 0 < result.final_summary["r_factor"] < result.runs[0]["summary"]["r_factor"]
    assert result.spec.source == "library:water" and result.spec.name == "waterbox"
    assert [dataset.label for dataset in result.spec.datasets] == ["D2O", "HDO", "H2O"]

    analysis = session.analysis
    assert analysis is not None
    assert set(analysis["fit"]) == {"D2O", "HDO", "H2O"}
    assert set(analysis["shells"]) == {"OW-OW", "OW-HW", "HW-HW"}
    assert analysis["schedule"]["total_iterations"] == 31
    assert analysis["convergence"]["r_factor"] == pytest.approx(
        result.final_summary["r_factor"], rel=1e-6
    )
    assert any("quick schedule" in flag for flag in analysis["flags"])
    assert "01_convergence" in session.figures and "00_data" in session.figures
    assert pipeline.agent.narrative_source == "model"
    assert session.data_descriptions is not None and len(session.data_descriptions) == 3
    assert session.data_descriptions[0]["path"] == str(WATER_FILES[0])


def test_run_folder(pipeline: SimpleNamespace) -> None:
    folder = pipeline.session.work_dir
    for name in (
        "spec.json",
        "run_waterbox.py",
        "data/SLS18502-D2O.mint01",
        "epsr/waterbox.EPSR.inp",
        "epsr/waterbox.EPSR.u01",
        "trajectory/meta.json",
        "manifest.json",
        "run.log",
        "analysis.json",
        "figures/00_data.png",
        "figures/01_convergence.png",
        "figures/11_fit_final.png",
        "figures/summary.md",
        "report.pdf",
        DISCUSSION_NAME,
        TRANSCRIPT_NAME,
    ):
        assert (folder / name).is_file(), name
    spec = SystemSpec.load(folder / "spec.json")
    assert spec.source == "library:water" and spec.validate(check_files=True) == []

    saved = analysis_module.load_analysis(folder / "analysis.json")
    assert saved["figures"]["01_convergence"] == str(
        folder / "figures" / "01_convergence.png"
    )
    assert (
        saved["runs"][-1]["r_factor"]
        == pipeline.session.analysis["runs"][-1]["r_factor"]
    )
    assert (folder / DISCUSSION_NAME).read_text() == NARRATIVE.rstrip("\n") + "\n"

    script = (folder / "run_waterbox.py").read_text()
    assert "import epsr" in script and "waterbox.EPSR.inp" in script
    assert "import mellor" not in script and "from mellor" not in script

    transcript = _read_transcript(folder / TRANSCRIPT_NAME)
    assert transcript[0]["role"] == "pipeline" and transcript[0]["content"] == "start"
    assert (
        transcript[0]["spec_given"] is False
        and transcript[0]["narrative_by_model"] is True
    )
    assert transcript[-1]["role"] == "pipeline" and transcript[-1]["content"] == "done"
    tool_entries = [entry["name"] for entry in transcript if entry["role"] == "tool"]
    assert tool_entries == ["inspect_data", "library_spec", "build_and_run"]
    narrative = [entry for entry in transcript if entry["role"] == "narrative"]
    assert len(narrative) == 1 and narrative[0]["source"] == "model"
    assert narrative[0]["model"] == "scripted" and narrative[0]["chars"] == len(
        NARRATIVE.strip()
    )


def test_report_pdf(pipeline: SimpleNamespace) -> None:
    text = _pdf_text(pipeline.outcome.report_path)
    assert "tetrahedral" in text and "Discussion" in text
    assert "Written by the language model" in text and "scripted" in text
    for heading in report_module.SECTION_HEADINGS:
        assert heading in text, heading
    assert "run_waterbox.py" in text and "waterbox.EPSR.inp" in text


def test_events(pipeline: SimpleNamespace) -> None:
    events = pipeline.events
    kinds = [event["type"] for event in events]
    phases = [event["phase"] for event in events if event["type"] == "phase"]
    assert phases == ["setup", "narrative", "report"]
    calls = [event["name"] for event in events if event["type"] == "tool_call"]
    assert calls == ["inspect_data", "library_spec", "build_and_run"]
    results = [event for event in events if event["type"] == "tool_result"]
    assert [event["ok"] for event in results] == [True, True, True]
    progress = [event["line"] for event in events if event["type"] == "progress"]
    assert progress[0].startswith("$ ") and any("run 4/4" in line for line in progress)
    assert any("makemole Water" in line for line in progress)
    assert kinds[-1] == "report" and events[-1]["path"] == str(
        pipeline.outcome.report_path
    )


def test_chat_answers_through_get_results(pipeline: SimpleNamespace) -> None:
    agent, session = pipeline.agent, pipeline.session
    answer = agent.chat(QUESTION)
    final = session.analysis["convergence"]["r_factor"]
    assert answer == f"The final R-factor is {final:.3e}."
    assert (
        agent.messages[-2]["role"] == "tool"
        and agent.messages[-2]["name"] == "get_results"
    )
    assert agent.steps == pipeline.outcome.steps + 2
    context = agent.messages[-4]["content"]
    assert "follow-up questions" in context and f"User: {QUESTION}" in context
    assert "figures/01_convergence.png" in context and "report.pdf" in context
    transcript = _read_transcript(session.work_dir / TRANSCRIPT_NAME)
    assert transcript[-1]["role"] == "chat" and transcript[-1]["content"] == "turn done"
