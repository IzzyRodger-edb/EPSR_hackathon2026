"""The progress stream of a run.

Mellor reports progress by calling `on_event` with dicts (integration.md
section 4.3). That callback runs in the worker *process*, while the HTTP
request that streams it runs in the web process, so the two are joined through
a file: the worker appends one JSON object per line to `events.jsonl` and the
SSE endpoint tails it. That also means a reload of the results page can replay
everything that happened, and a server restart loses nothing.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Iterator

EVENTS_FILE = "events.jsonl"
LOG_FILE = "run.log"

#: Event types Mellor emits, plus the ones the app adds for its own phases.
MELLOR_TYPES = ("phase", "assistant", "tool_call", "tool_result", "progress",
                "note", "fallback", "narrative", "report")
APP_TYPES = ("queued", "started", "finished", "failed", "cancelled")


class EventWriter:
    """Append-only writer, used from the worker process."""

    def __init__(self, run_dir: Path) -> None:
        self.path = run_dir / EVENTS_FILE
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._seq = _last_seq(self.path)

    def emit(self, event: dict[str, Any]) -> None:
        self._seq += 1
        record = {"seq": self._seq, "at": time.time(), **event}
        try:
            line = json.dumps(record, default=_stringify)
        except (TypeError, ValueError):
            line = json.dumps({"seq": self._seq, "at": time.time(), "type": "note",
                               "content": "an event could not be serialised"})
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")
            handle.flush()

    def __call__(self, event: dict[str, Any]) -> None:
        """So the writer can be passed straight to Mellor as `on_event`."""
        self.emit(event if isinstance(event, dict) else {"type": "note", "content": str(event)})


def _stringify(value: Any) -> str:
    return str(value)


def _last_seq(path: Path) -> int:
    last = 0
    for event in read_events(path.parent):
        last = max(last, int(event.get("seq", 0)))
    return last


def read_events(run_dir: Path, after_seq: int = 0) -> list[dict[str, Any]]:
    """Every event after `after_seq`. Partial trailing lines are skipped."""
    path = run_dir / EVENTS_FILE
    out: list[dict[str, Any]] = []
    try:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    event = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if isinstance(event, dict) and int(event.get("seq", 0)) > after_seq:
                    out.append(event)
    except OSError:
        return out
    return out


def follow_events(run_dir: Path, after_seq: int = 0) -> Iterator[dict[str, Any]]:
    """Yield what is already there, then whatever is appended next.

    The caller decides when to stop (the SSE endpoint stops when the job reaches
    a terminal state and the stream has caught up).
    """
    seq = after_seq
    while True:
        new = read_events(run_dir, seq)
        if not new:
            return
        for event in new:
            seq = int(event.get("seq", seq))
            yield event


def summarise(event: dict[str, Any]) -> str:
    """One line for the log pane. Keeps tool arguments short."""
    kind = str(event.get("type", "?"))
    if kind == "progress":
        return str(event.get("line", "")).rstrip()
    if kind == "phase":
        return f"-- phase: {event.get('phase')}"
    if kind == "tool_call":
        args = json.dumps(event.get("arguments", {}), default=_stringify)
        if len(args) > 300:
            args = args[:297] + "..."
        return f"-> {event.get('name')}({args})"
    if kind == "tool_result":
        result = str(event.get("result", ""))
        if len(result) > 300:
            result = result[:297] + "..."
        return f"<- {event.get('name')}: {result}"
    if kind == "assistant":
        content = str(event.get("content", "")).strip()
        return f"model: {content}" if content else "model: (tool calls)"
    if kind == "note":
        return f"note: {event.get('content')}"
    if kind == "fallback":
        return (f"note: the model's specification was unusable; fell back to the "
                f"library system {event.get('system')} ({event.get('reason')})")
    if kind == "report":
        return f"-- report written: {event.get('path')}"
    return f"{kind}: {event.get('content') or event.get('message') or ''}".rstrip()
