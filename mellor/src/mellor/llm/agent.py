"""Mellor, the agent: a tool-call loop around a chat backend.

``Mellor.run_pipeline`` is the agentic workflow (docs/architecture.md §7):

1. setup: the model inspects the data and calls ``build_and_run`` with a valid
   specification (``max_spec_attempts`` validation failures, then the library
   fallback through ``library.auto_spec`` with the reason recorded in the session);
2. the run, the analysis and the figures happen inside the ``build_and_run`` tool;
3. narrative: the model writes the discussion from ``analysis_markdown`` (fallback:
   that markdown itself);
4. the report is built by code.

``run_pipeline(spec=...)`` skips step 1: a specification given by the user or built
from the library (``mellor run --spec`` / ``--system``) goes straight to
``build_and_run`` and the model only writes the discussion; ``narrative=False`` skips
the model there too (a run without any language model: the factual summary is used
and labelled as such).

``Mellor.chat`` continues the conversation with the same tools (questions about the
results, searching and downloading comparable data, running again).

Every model message, tool call and tool result is appended to
``<work_dir>/mellor_transcript.jsonl``; a ``narrative`` entry records who wrote the
discussion (``source`` model / fallback / summary and the model's name).
"""

from __future__ import annotations

import datetime as _dt
import json
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..spec import SpecError, SystemSpec
from .backend import ChatBackend
from .prompts import (
    SYSTEM_PROMPT,
    format_chat_context,
    format_narrative_task,
    format_setup_task,
)
from .tools import Session, ToolRegistry, compact_analysis, default_tools

__all__ = [
    "Mellor",
    "PipelineOutcome",
    "TRANSCRIPT_NAME",
    "DISCUSSION_NAME",
    "NO_MODEL",
    "NO_MODEL_NOTE",
]

TRANSCRIPT_NAME = "mellor_transcript.jsonl"
DISCUSSION_NAME = "discussion.md"
#: A narrative shorter than this is treated as "no narrative" (fallback text used).
MIN_NARRATIVE_CHARS = 80
#: The tools the model may call while writing the discussion.
NARRATIVE_TOOLS = ["get_results", "get_partial"]
#: Tool calls that count as a failed specification attempt.
SPEC_TOOLS = ("validate_spec", "build_and_run")
#: The ``model_name`` the report receives when no language model wrote the
#: discussion (``mellor.report`` renders names starting with ``none`` accordingly).
NO_MODEL = "none: the discussion is the factual summary computed by mellor.analysis"
#: The note that heads a discussion written without a language model.
NO_MODEL_NOTE = (
    "No language model was used for this run; this section is the factual summary "
    "computed from the run."
)


@dataclass
class PipelineOutcome:
    session: Session
    report_path: Path | None
    final_message: str
    spec_source: str
    steps: int  # model calls made


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


class Mellor:
    def __init__(
        self,
        backend: ChatBackend,
        session: Session,
        *,
        tools: ToolRegistry | None = None,
        max_steps: int = 30,
        max_spec_attempts: int = 4,
        on_event: Callable[[dict[str, Any]], None] | None = None,
    ) -> None:
        """``on_event`` receives dicts (``{"type": "assistant"|"tool_call"|"tool_result"|..., ...}``)
        for progress display; ``tools`` defaults to ``default_tools(session)``.

        ``max_steps`` bounds the model calls of one phase (the setup, the narrative,
        one chat turn); ``steps`` counts them over the whole conversation."""
        self.backend = backend
        self.session = session
        self.tools = tools if tools is not None else default_tools(session)
        self.max_steps = int(max_steps)
        self.max_spec_attempts = int(max_spec_attempts)
        self.on_event = on_event
        self.steps = 0
        self.spec_failures = 0
        self.notes: list[str] = []
        self.narrative: str = ""
        self.narrative_source: str = ""
        self.transcript_path = Path(session.work_dir) / TRANSCRIPT_NAME
        self._transcript: list[dict[str, Any]] = []
        self._messages: list[dict[str, Any]] = [
            {"role": "system", "content": SYSTEM_PROMPT}
        ]
        self._chat_context_added = False
        self._pipeline_ran = False
        if not session.model_name:
            session.model_name = getattr(backend, "name", "") or type(backend).__name__
        if session.progress is None:
            session.progress = self._progress

    # ------------------------------------------------------------------ bookkeeping
    @property
    def messages(self) -> list[dict[str, Any]]:
        """The conversation so far (system prompt first)."""
        return self._messages

    def _emit(self, event: dict[str, Any]) -> None:
        if self.on_event is None:
            return
        try:
            self.on_event(event)
        except Exception as error:  # a display problem must not stop the run
            self._record(
                {"role": "note", "content": f"on_event callback failed: {error}"}
            )

    def _progress(self, line: str) -> None:
        self._emit({"type": "progress", "line": line})

    def _record(self, entry: dict[str, Any]) -> None:
        """Append one line to the transcript (rewriting it if the file vanished, e.g.
        because the run folder was recreated)."""
        entry = {"time": _now(), **entry}
        self._transcript.append(entry)
        path = self.transcript_path
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.exists():
                with path.open("a", encoding="utf-8") as handle:
                    handle.write(
                        json.dumps(entry, ensure_ascii=False, default=str) + "\n"
                    )
            else:
                with path.open("w", encoding="utf-8") as handle:
                    for item in self._transcript:
                        handle.write(
                            json.dumps(item, ensure_ascii=False, default=str) + "\n"
                        )
        except OSError as error:
            self._emit(
                {"type": "note", "content": f"cannot write the transcript: {error}"}
            )

    def _append(self, message: dict[str, Any], **extra: Any) -> None:
        self._messages.append(message)
        self._record({**message, **extra})

    def _note(self, text: str) -> None:
        self.notes.append(text)
        self._record({"role": "note", "content": text})
        self._emit({"type": "note", "content": text})

    # ------------------------------------------------------------------ the pipeline
    def run_pipeline(
        self,
        user_prompt: str | None = None,
        *,
        spec: SystemSpec | Mapping[str, Any] | None = None,
        narrative: bool = True,
    ) -> PipelineOutcome:
        """The whole workflow (docs/architecture.md §7) on the session's data files.

        ``spec`` (a :class:`~mellor.spec.SystemSpec` or its dict) skips the setup
        phase: the specification is run through ``build_and_run`` as given (its
        ``source`` is kept, ``user`` when empty) and the model is told about it, so it
        can still write the discussion and answer questions. ``narrative=False``
        does not ask the model for the discussion either: the factual summary of the
        analysis is used, headed by :data:`NO_MODEL_NOTE` (a run without any
        language model).

        Raises ``ValueError``/``FileNotFoundError`` for missing data files,
        :class:`~mellor.spec.SpecError` when a given specification is invalid and
        ``RuntimeError`` when the EPSR run (or the library fallback) fails.
        """
        session = self.session
        if not session.data_files:
            raise ValueError("the session has no data files (.mint01) to refine")
        missing = [str(path) for path in session.data_files if not Path(path).is_file()]
        if missing:
            raise FileNotFoundError(f"data file(s) not found: {missing}")
        Path(session.work_dir).mkdir(parents=True, exist_ok=True)
        self._record(
            {
                "role": "pipeline",
                "content": "start",
                "data_files": [str(p) for p in session.data_files],
                "quick": session.quick,
                "model": self.session.model_name,
                "spec_given": spec is not None,
                "narrative_by_model": narrative,
            }
        )

        descriptions = self._describe_data()
        if spec is not None:
            self._run_given_spec(spec, user_prompt)
        else:
            listing = self._library_listing()
            self._emit({"type": "phase", "phase": "setup"})
            self.spec_failures = 0
            self._append(
                {
                    "role": "user",
                    "content": format_setup_task(
                        descriptions,
                        listing,
                        user_prompt,
                        str(session.work_dir),
                        session.quick,
                    ),
                }
            )
            self._loop(stop_after="build_and_run")
            if session.result is None:
                self._fallback()
        # the analysis summary enters the conversation with the narrative task;
        # without it the chat context repeats the summary for the model
        self._pipeline_ran = narrative

        self._emit({"type": "phase", "phase": "narrative"})
        narrative_text, source = self._narrative(ask_model=narrative)
        self._emit({"type": "phase", "phase": "report"})
        report_path = self._build_report(narrative_text, descriptions)
        final_message = self._final_message(source)
        self._record(
            {
                "role": "pipeline",
                "content": "done",
                "report": str(report_path),
                "steps": self.steps,
            }
        )
        return PipelineOutcome(
            session=session,
            report_path=report_path,
            final_message=final_message,
            spec_source=session.spec_source,
            steps=self.steps,
        )

    def _describe_data(self) -> list[dict[str, Any]]:
        from .. import data as data_module

        session = self.session
        try:
            files_arg: list[str | Path] = list(session.data_files)
            descriptions = [dict(item) for item in data_module.describe_all(files_arg)]
        except Exception as error:
            self._note(
                f"the data files could not be described ({type(error).__name__}: {error}); only their names are given to the model"
            )
            descriptions = [{"file": Path(path).name} for path in session.data_files]
        for description, path in zip(descriptions, session.data_files, strict=False):
            description.setdefault("file", Path(path).name)
            description.setdefault("path", str(path))
        session.data_descriptions = descriptions
        return descriptions

    def _library_listing(self) -> list[dict[str, Any]]:
        from .. import library as library_module

        try:
            return [dict(item) for item in library_module.describe_library()]
        except Exception as error:
            self._note(
                f"the library listing is unavailable ({type(error).__name__}: {error})"
            )
            return []

    def _run_given_spec(
        self, spec: SystemSpec | Mapping[str, Any], user_prompt: str | None
    ) -> None:
        """Run a specification that was given (``mellor run --spec`` / ``--system``)
        through ``build_and_run`` without asking the model, and tell the model."""
        data = spec.to_dict() if isinstance(spec, SystemSpec) else dict(spec)
        source = str(data.get("source") or "").strip() or "user"
        data["source"] = source
        self._emit({"type": "phase", "phase": "run"})
        self._record(
            {
                "role": "spec",
                "content": f"specification given ({source}); the model was not asked for the setup",
                "source": source,
                "name": data.get("name"),
            }
        )
        arguments = {"spec": data}
        self._emit(
            {"type": "tool_call", "name": "build_and_run", "arguments": arguments}
        )
        result = self.tools.call("build_and_run", arguments)
        text = self.tools.serialise(result)
        ok = isinstance(result, dict) and result.get("ok") is True
        self._emit(
            {
                "type": "tool_result",
                "name": "build_and_run",
                "ok": ok,
                "preview": text[:200],
            }
        )
        request = (user_prompt or "").strip()
        self._append(
            {
                "role": "user",
                "content": (
                    f"[Mellor] The specification was given ({source}) and run without "
                    "asking you for the setup."
                    + (f" The user's request: {request}" if request else "")
                    + f"\nResult of build_and_run:\n{text}"
                ),
            }
        )
        if not ok:
            if isinstance(result, dict) and result.get("errors"):
                raise SpecError([str(problem) for problem in result["errors"]])
            detail = result.get("error") if isinstance(result, dict) else result
            log_tail = (
                "\n".join(result.get("log_tail", []))
                if isinstance(result, dict)
                else ""
            )
            raise RuntimeError(
                f"the run of the given specification failed: {detail}\n{log_tail}".rstrip()
            )

    def _fallback(self) -> None:
        """Run the best library specification when the model produced none."""
        from .. import data as data_module
        from .. import library as library_module

        session = self.session
        if self.spec_failures:
            reason = f"the model produced {self.spec_failures} invalid specification(s)"
        else:
            reason = "the model stopped without a successful build_and_run"
        datasets = [data_module.read_mint01(path) for path in session.data_files]
        try:
            spec, why = library_module.auto_spec(datasets)
        except LookupError as error:
            raise RuntimeError(
                f"Mellor could not produce a valid specification ({reason}) and no "
                f"library system matches the data: {error}"
            ) from error
        key = spec.source.split(":")[-1] if spec.source else spec.name
        spec.source = f"fallback:library:{key}"
        note = f"Mellor fell back to the built-in library system '{key}' because {reason}. {why}".strip()
        spec.notes = f"{spec.notes}\n\n{note}" if spec.notes else note
        self._note(note)
        self._emit({"type": "fallback", "system": key, "reason": reason})
        session.events.append(
            {"type": "fallback", "system": key, "reason": reason, "time": _now()}
        )
        result = self.tools.call("build_and_run", {"spec": spec.to_dict()})
        text = self.tools.serialise(result)
        self._record(
            {
                "role": "fallback",
                "name": "build_and_run",
                "system": key,
                "content": text,
            }
        )
        self._append(
            {
                "role": "user",
                "content": f"[Mellor] The specification came from the library fallback ('{key}') because {reason}. Result of build_and_run:\n{text}",
            }
        )
        if not (isinstance(result, dict) and result.get("ok") is True):
            detail = (
                result.get("error") or result.get("errors")
                if isinstance(result, dict)
                else result
            )
            log_tail = (
                "\n".join(result.get("log_tail", []))
                if isinstance(result, dict)
                else ""
            )
            raise RuntimeError(
                f"the library fallback run failed: {detail}\n{log_tail}".rstrip()
            )

    def _narrative(self, ask_model: bool = True) -> tuple[str, str]:
        """``(discussion markdown, source)``; the source is ``model``, ``fallback``
        (the model was asked but wrote nothing usable) or ``summary`` (not asked)."""
        from .. import analysis as analysis_module

        session = self.session
        markdown = ""
        if session.analysis:
            try:
                markdown = analysis_module.analysis_markdown(session.analysis)
            except Exception as error:
                self._note(
                    f"analysis_markdown failed ({type(error).__name__}: {error}); the compact analysis is used"
                )
                markdown = (
                    "```json\n"
                    + json.dumps(compact_analysis(session.analysis), indent=1)
                    + "\n```"
                )
        else:
            markdown = "(the analysis of the run is not available; see the run table in the report)"
        if not ask_model:
            narrative, source = f"_{NO_MODEL_NOTE}_\n\n{markdown}", "summary"
        else:
            notes = session.result.spec.notes if session.result is not None else ""
            self._append(
                {"role": "user", "content": format_narrative_task(markdown, notes)}
            )
            try:
                text = self._loop(allowed_tools=NARRATIVE_TOOLS).strip()
            except (
                Exception
            ) as error:  # the run is done: the report must still be written
                self._note(
                    f"the model failed while writing the discussion ({type(error).__name__}: {error})"
                )
                text = ""
            if len(text) >= MIN_NARRATIVE_CHARS:
                narrative, source = text, "model"
            else:
                note = "The language model did not produce a discussion; this section is the automatic factual summary computed from the run."
                self._note(note)
                narrative, source = f"_{note}_\n\n{markdown}", "fallback"
        self.narrative, self.narrative_source = narrative, source
        self._record(
            {
                "role": "narrative",
                "source": source,
                "model": self._model_name(),
                "chars": len(narrative),
            }
        )
        session.events.append(
            {
                "type": "narrative",
                "source": source,
                "chars": len(narrative),
                "time": _now(),
            }
        )
        return narrative, source

    def _model_name(self) -> str:
        return str(getattr(self.backend, "name", "") or self.session.model_name)

    def _report_model_name(self) -> str:
        """What the report says wrote the discussion: the model, or :data:`NO_MODEL`."""
        return NO_MODEL if self.narrative_source == "summary" else self._model_name()

    def _build_report(self, narrative: str, descriptions: list[dict[str, Any]]) -> Path:
        from .. import report as report_module

        session = self.session
        if session.result is None:
            raise RuntimeError("no run to report on")
        folder = Path(session.result.folder)
        try:
            (folder / DISCUSSION_NAME).write_text(narrative + "\n", encoding="utf-8")
        except OSError as error:
            self._note(f"could not write {DISCUSSION_NAME}: {error}")
        path = report_module.build_pdf(
            session.result,
            session.analysis or {},
            session.figures,
            narrative,
            folder / "report.pdf",
            model_name=self._report_model_name(),
            data_descriptions=descriptions,
        )
        session.report_path = Path(path)
        session.events.append({"type": "report", "path": str(path), "time": _now()})
        self._emit({"type": "report", "path": str(path)})
        return session.report_path

    def _final_message(self, narrative_source: str) -> str:
        session = self.session
        result = session.result
        parts = []
        if result is not None:
            parts.append(
                f"EPSR refinement of {result.spec.title or result.spec.name} finished."
            )
        parts.append(f"Specification: {session.spec_source or 'unknown'}.")
        if session.analysis:
            flags = session.analysis.get("flags") or []
            if flags:
                parts.append("Flags: " + "; ".join(str(flag) for flag in flags) + ".")
        final = result.final_summary if result is not None else {}
        if final.get("r_factor") is not None:
            parts.append(f"Final R-factor {final['r_factor']:.4g}.")
        if narrative_source == "fallback":
            parts.append(
                "The discussion is the automatic summary (the model wrote none)."
            )
        elif narrative_source == "summary":
            parts.append(
                "The discussion is the factual summary (no language model was used)."
            )
        if session.report_path is not None:
            parts.append(f"Report: {session.report_path}")
        if self.notes:
            parts.append("Notes: " + " | ".join(self.notes))
        return " ".join(parts)

    # ------------------------------------------------------------------ chat
    def chat(self, user_message: str) -> str:
        """One user turn: tool loop until the model answers without tool calls (or ``finish``)."""
        content = user_message.strip()
        if not self._chat_context_added:
            content = self._chat_context() + "\n\nUser: " + content
            self._chat_context_added = True
        self._append({"role": "user", "content": content})
        self._emit({"type": "phase", "phase": "chat"})
        text = self._loop().strip()
        if not text:
            text = "Mellor gave no answer." + (
                " (" + "; ".join(self.notes[-1:]) + ")" if self.notes else ""
            )
        self._record({"role": "chat", "content": "turn done", "steps": self.steps})
        return text

    def _chat_context(self) -> str:
        from .. import analysis as analysis_module

        session = self.session
        markdown = None
        if session.analysis and not self._pipeline_ran:
            try:
                markdown = analysis_module.analysis_markdown(session.analysis)
            except Exception as error:
                self._note(
                    f"analysis_markdown failed ({type(error).__name__}: {error})"
                )
                markdown = (
                    "```json\n"
                    + json.dumps(compact_analysis(session.analysis), indent=1)
                    + "\n```"
                )
        files: list[str] = []
        folder = Path(session.work_dir)
        if folder.is_dir():
            for path in sorted(folder.iterdir()):
                if path.name.startswith("."):
                    continue
                files.append(path.name + ("/" if path.is_dir() else ""))
        files += [f"figures/{Path(path).name}" for path in session.figures.values()]
        report = str(session.report_path) if session.report_path else None
        return format_chat_context(markdown, files, report)

    # ------------------------------------------------------------------ the loop
    def _loop(
        self, allowed_tools: list[str] | None = None, *, stop_after: str | None = None
    ) -> str:
        """Call the backend, execute tool calls, repeat; return the model's final text.

        ``stop_after`` names a tool whose successful call ends the loop (used for
        ``build_and_run`` in the setup phase). Guards: ``max_steps``, repeated identical
        tool calls (the same name and arguments three times in a row ends the loop with
        a note), and an empty reply twice in a row."""
        schemas = self.tools.schemas(allowed_tools)
        phase_steps = 0
        last_key: tuple[str, str] | None = None
        repeats = 0
        empty_replies = 0
        final_text = ""
        while True:
            if phase_steps >= self.max_steps:
                self._note(
                    f"stopped: the limit of {self.max_steps} model calls was reached"
                )
                break
            reply = self.backend.chat(self._messages, schemas)
            self.steps += 1
            phase_steps += 1
            self._append(reply.to_message(), raw=reply.raw)
            self._emit(
                {
                    "type": "assistant",
                    "step": self.steps,
                    "content": reply.content,
                    "tool_calls": [
                        {"name": c.name, "arguments": c.arguments}
                        for c in reply.tool_calls
                    ],
                }
            )
            if not reply.tool_calls:
                malformed = "<tool_call>" in reply.content or "<tool_call>" in (
                    reply.raw or ""
                )
                if reply.content.strip() and not malformed:
                    final_text = reply.content
                    break
                empty_replies += 1
                if empty_replies >= 2:
                    self._note("stopped: two empty or malformed replies in a row")
                    break
                if malformed:
                    nudge = (
                        "Your tool call could not be parsed. Reply with exactly one "
                        '<tool_call> block containing valid JSON: {"name": "<tool>", "arguments": {...}}.'
                    )
                else:
                    nudge = "Your reply was empty. Call a tool, or answer in text."
                self._append({"role": "user", "content": nudge})
                continue
            empty_replies = 0
            stop = False
            for index, call in enumerate(reply.tool_calls, start=1):
                if not call.id:  # a backend may leave ids empty
                    call.id = f"call_{index}"
                key = (
                    call.name,
                    json.dumps(call.arguments, sort_keys=True, default=str),
                )
                repeats = repeats + 1 if key == last_key else 1
                last_key = key
                if repeats >= 3:
                    self._note(
                        f"stopped: the tool call {call.name} was repeated three times with the same arguments"
                    )
                    self._append(
                        {
                            "role": "tool",
                            "name": call.name,
                            "tool_call_id": call.id,
                            "content": json.dumps(
                                {
                                    "error": "this identical call was already made twice; Mellor stopped the loop"
                                }
                            ),
                        }
                    )
                    stop = True
                    break
                if allowed_tools is not None and call.name not in allowed_tools:
                    result: Any = {
                        "error": f"tool {call.name!r} is not available in this phase",
                        "hint": "available: " + ", ".join(allowed_tools),
                    }
                else:
                    self._emit(
                        {
                            "type": "tool_call",
                            "name": call.name,
                            "arguments": call.arguments,
                        }
                    )
                    result = self.tools.call(call.name, call.arguments)
                text = self.tools.serialise(result)
                self._append(
                    {
                        "role": "tool",
                        "name": call.name,
                        "tool_call_id": call.id,
                        "content": text,
                    }
                )
                failed = isinstance(result, dict) and "error" in result
                self._emit(
                    {
                        "type": "tool_result",
                        "name": call.name,
                        "ok": not failed,
                        "preview": text[:200],
                    }
                )
                if call.name == "finish":
                    final_text = str(call.arguments.get("message", "")).strip()
                    stop = True
                    break
                if (
                    stop_after
                    and call.name == stop_after
                    and isinstance(result, dict)
                    and result.get("ok") is True
                ):
                    final_text = reply.content
                    stop = True
                    break
                if call.name in SPEC_TOOLS and _is_invalid_spec_result(result):
                    self.spec_failures += 1
                    if (
                        stop_after == "build_and_run"
                        and self.spec_failures >= self.max_spec_attempts
                    ):
                        self._note(
                            f"stopped: {self.spec_failures} invalid specification(s) (limit {self.max_spec_attempts})"
                        )
                        stop = True
                        break
            if stop:
                break
        return final_text


def _is_invalid_spec_result(result: Any) -> bool:
    if not isinstance(result, dict):
        return False
    if result.get("valid") is False:
        return True
    return result.get("ok") is False and "errors" in result
