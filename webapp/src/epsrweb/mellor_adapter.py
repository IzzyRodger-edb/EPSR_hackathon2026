"""The only module that knows how to call Mellor.

Everything the app does with the agent goes through here, for two reasons. The
calls are fixed by `mellor/docs/integration.md`, so keeping them in one file
makes them easy to check against that contract. And Mellor is a sibling package
that needs a Rust build: when it is not importable this module falls back to
`demo`, so the app, its tests and the Overleaf path all work on a machine where
the engine has never been compiled.
"""

from __future__ import annotations

import importlib
import logging
from pathlib import Path
from typing import Callable, Sequence

from . import demo, plots, storage

LOG = logging.getLogger(__name__)

#: Backends the user may pick, as (value, label, needs a token) triples.
BACKENDS: tuple[tuple[str, str, bool], ...] = (
    ("auto", "Automatic (local model if present, else the hosted API)", False),
    ("hf-api", "Hugging Face Inference Providers (hosted, needs a token)", True),
    ("transformers", "Local model on this machine's GPU", False),
    ("library", "No language model: use a validated model system", False),
    ("scripted", "Scripted replies (smoke test, no model)", False),
)

#: Validated systems from `mellor.library.KNOWN_SYSTEMS`, used by the `library`
#: backend. Queried from the package when it is installed.
FALLBACK_SYSTEMS = ("water", "silica", "dmso_water", "benzene", "cu_btc")


def mellor_status() -> dict[str, object]:
    """Whether a real refinement can run here, and what is missing when it cannot.

    Three things have to be in place, and they fail differently often enough to be
    worth telling apart: the agent itself, the EPSR engine (`epsr`, a compiled Rust
    extension and the usual culprit) and the plotting toolkit (`epsrvis`).
    `mellor` imports both at module level, so one missing piece makes
    `import mellor` fail with a message that does not name the cause.
    """
    for module, hint in (
        ("epsr", "the EPSR engine is not installed: build it with "
                 "`uv sync --extra agent` in webapp/ (needs a Rust toolchain and a C compiler)"),
        ("epsrvis", "the plotting toolkit is not installed: "
                    "`uv sync --extra agent` in webapp/"),
        ("mellor", "the agent is not installed: `uv sync --extra agent` in webapp/"),
    ):
        if importlib.util.find_spec(module) is None:
            return {"available": False, "reason": hint, "missing": module}
    try:
        importlib.import_module("mellor.llm.agent")
    except Exception as exc:
        return {"available": False, "missing": "mellor",
                "reason": f"mellor is installed but could not be imported "
                          f"({type(exc).__name__}: {exc})"}
    import mellor

    return {"available": True, "reason": None, "missing": None,
            "version": getattr(mellor, "__version__", "unknown")}


def known_systems() -> tuple[str, ...]:
    try:
        from mellor import library

        systems = tuple(library.KNOWN_SYSTEMS)
        return systems or FALLBACK_SYSTEMS
    except Exception:
        return FALLBACK_SYSTEMS


def figure_captions() -> dict[str, str]:
    """Mellor's captions, falling back to our copy of them."""
    captions = dict(storage.FALLBACK_CAPTIONS)
    try:
        from mellor import figures

        captions.update(dict(figures.FIGURE_CAPTIONS))
    except Exception:
        pass
    captions.update(storage.SEABORN_CAPTIONS)
    return captions


def backend_string(backend: str, model: str | None) -> str:
    """Build the string `mellor.llm.backend.make_backend` expects."""
    if backend in {"hf-api", "transformers"} and model:
        return f"{backend}:{model}"
    return backend


def _run_with_agent(
    run_dir: Path, data_files: Sequence[Path], prompt: str | None, backend: str,
    model: str | None, quick: bool, emit: Callable[[dict], None],
) -> dict:
    """The documented pipeline: the model proposes the system, code runs EPSR."""
    from mellor.llm.agent import Mellor
    from mellor.llm.backend import make_backend
    from mellor.llm.tools import Session

    session = Session(work_dir=run_dir, data_files=list(data_files), quick=quick)
    agent = Mellor(make_backend(backend_string(backend, model)), session, on_event=emit)
    outcome = agent.run_pipeline(user_prompt=prompt)
    # `PipelineOutcome.report_path` is None when the report could not be built;
    # the run itself may still have succeeded, so this is not an error.
    report = outcome.report_path or (run_dir / "report.pdf")
    return {
        "report": str(report) if Path(report).is_file() else None,
        "spec_source": outcome.spec_source or "model",
        "message": outcome.final_message,
        "steps": outcome.steps,
        "demo": False,
    }


def _run_without_model(
    run_dir: Path, data_files: Sequence[Path], system: str | None, quick: bool,
    emit: Callable[[dict], None],
) -> dict:
    """The fixed-specification path, for the "no language model" option.

    integration.md section 4.4 in order: read the data, take a validated
    specification from the library (or let `auto_spec` pick one from the data
    titles), run it, analyse it, draw the figures, write the discussion and build
    the PDF. Deterministic, free, and the baseline an agent-proposed run can be
    compared against.
    """
    from mellor import analysis, data, figures, library, report, runner
    from mellor.llm.agent import DISCUSSION_NAME, NO_MODEL, NO_MODEL_NOTE

    def progress(line: object) -> None:
        emit({"type": "progress", "line": str(line).rstrip()})

    emit({"type": "phase", "phase": "setup"})
    datasets = [data.read_mint01(path) for path in data_files]
    if system:
        spec = library.make_spec(system, datasets)
        source = f"library:{system}"
        reason = f"the {system} template was chosen by the user"
    else:
        # auto_spec returns the specification AND the reason it picked it.
        spec, reason = library.auto_spec(datasets)
        source = spec.source or "library:auto"
    emit({"type": "note",
          "content": f"specification from {source}; no language model was used ({reason})"})
    for warning in spec.validate():
        emit({"type": "note", "content": f"specification warning: {warning}"})

    # overwrite=True because the folder is not empty: it already holds the
    # uploads, this app's event log and its worker log. run_spec only removes
    # the products of an earlier run (runner.RUN_PRODUCTS).
    result = runner.run_spec(spec, run_dir, quick=quick, overwrite=True, progress=progress)
    summary = analysis.analyse(result)
    analysis.save_analysis(summary, run_dir / "analysis.json")
    pngs = figures.make_figures(result)

    emit({"type": "phase", "phase": "report"})
    blank = "\n\n"
    discussion = f"_{NO_MODEL_NOTE}_{blank}{analysis.analysis_markdown(summary)}"
    (run_dir / DISCUSSION_NAME).write_text(discussion + "\n", encoding="utf-8")
    pdf = report.build_pdf(result, summary, pngs, discussion, model_name=NO_MODEL)
    return {
        "report": str(pdf),
        "spec_source": source,
        "message": f"Refined with the validated {system or spec.name} model system, "
                   "without a language model.",
        "demo": False,
    }


def run_pipeline(
    run_dir: Path,
    data_files: Sequence[Path],
    *,
    prompt: str | None = None,
    backend: str = "auto",
    model: str | None = None,
    quick: bool = False,
    system: str | None = None,
    emit: Callable[[dict], None] = lambda event: None,
) -> dict:
    """Run one refinement end to end and leave a complete run folder behind.

    Blocks for the whole run (minutes), so the caller must already be in a
    worker process. Returns the keys the app stores: `report`, `spec_source`,
    `message`, `demo`, plus `web_figures`.
    """
    status = mellor_status()
    if not status["available"]:
        emit({"type": "note",
              "content": f"Mellor is not installed ({status['reason']}); "
                         "writing a demonstration run instead."})
        outcome = demo.run_demo(run_dir, data_files, prompt=prompt, system=system,
                                quick=quick, emit=emit)
    elif backend == "library":
        outcome = _run_without_model(run_dir, data_files, system, quick, emit)
    else:
        outcome = _run_with_agent(run_dir, data_files, prompt, backend, model, quick, emit)

    if not outcome.get("demo"):
        # Tidy CSV tables for the app's own seaborn figures (integration.md 5.1).
        try:
            from mellor import export

            export.write_tables(run_dir)
            emit({"type": "progress", "line": "wrote tables/*.csv"})
        except Exception as exc:
            LOG.warning("export.write_tables failed: %s", exc)
            emit({"type": "note", "content": f"could not write the CSV tables ({exc}); "
                                             "the app's seaborn figures will be missing"})

    emit({"type": "phase", "phase": "figures"})
    try:
        outcome["web_figures"] = plots.make_seaborn_figures(run_dir)
        emit({"type": "progress",
              "line": f"drew {len(outcome['web_figures'])} seaborn figures"})
    except Exception as exc:
        LOG.warning("seaborn figures failed: %s", exc)
        outcome["web_figures"] = {}
        emit({"type": "note", "content": f"the seaborn figures failed ({exc})"})
    return outcome


def chat(run_dir: Path, message: str, *, backend: str = "auto", model: str | None = None) -> str:
    """Answer a follow-up question about a finished run.

    Rebuilds the session from its folder, as integration.md section 4.2 requires
    for a new process. The earlier conversation is not restored into the model's
    context; the results are, through its tools.
    """
    status = mellor_status()
    if not status["available"]:
        return _demo_answer(run_dir, message)
    from mellor.llm.agent import Mellor
    from mellor.llm.backend import make_backend
    from mellor.llm.tools import Session

    session = Session(work_dir=run_dir)
    session.load_existing()
    agent = Mellor(make_backend(backend_string(backend, model)), session)
    return str(agent.chat(message))


def _demo_answer(run_dir: Path, message: str) -> str:
    analysis = storage.read_json(run_dir / "analysis.json") or {}
    fit = analysis.get("fit") or []
    lines = [
        f"**{demo.BANNER}**",
        "",
        "There is no language model in this environment, so this is not an answer to "
        f"your question (*{message.strip()[:200]}*) but a summary of what the run "
        "folder holds:",
        "",
    ]
    for entry in fit:
        if isinstance(entry, dict):
            lines.append(
                f"- `{entry.get('dataset')}`: {entry.get('points')} points, "
                f"Q = {entry.get('q_min', 0):.2f} to {entry.get('q_max', 0):.1f} 1/A"
            )
    lines += ["", "Install `mellor` next to this app to ask real questions about a real "
              "refinement."]
    return "\n".join(lines)


def local_model_available() -> bool:
    """Whether a `transformers` backend could run here (torch is importable)."""
    return importlib.util.find_spec("torch") is not None


def backend_problem(backend: str, *, user_token: bool, server_token: bool) -> str | None:
    """Why this backend cannot run on this machine, or None if it can.

    `make_backend` raises at the point of use, which would cost the user a queued
    job and a failed run to find out. This says the same thing in the form.
    """
    if not mellor_status()["available"]:
        return None                      # demonstration mode serves every choice
    if backend in {"library", "scripted"}:
        return None
    has_token = user_token or server_token
    if backend == "transformers" and not local_model_available():
        return ("This server has no local model: torch is not installed. Choose the "
                "hosted backend, or 'No language model' to use a validated system.")
    if backend == "hf-api" and not has_token:
        return ("The hosted backend needs a Hugging Face token. Add one in Settings, "
                "or choose 'No language model' to use a validated system.")
    if backend == "auto" and not (local_model_available() or has_token):
        return ("No language model is available here: torch is not installed and no "
                "Hugging Face token is set. Add a token in Settings, or choose "
                "'No language model' to use a validated system.")
    return None
