"""What runs in the chat process (see `chat.py`).

A loop over a pipe: receive a question, answer it with Mellor, send the answer
back. The backend is built once per backend string and kept, so the local model
is loaded on the first question and reused for every one after it, instead of
being loaded again for each message.
"""

from __future__ import annotations

import gc
import logging
import os
from pathlib import Path
from typing import Any


def serve(conn: Any) -> None:
    """Entry point of the chat process. Returns when the pipe closes."""
    from . import processes

    # Its own process group, so that stopping it also stops an EPSR rerun that
    # the `build_and_run` tool started.
    processes.become_group_leader()
    os.environ.setdefault("MPLBACKEND", "Agg")
    logging.basicConfig(level=logging.WARNING,
                        format="%(asctime)s %(levelname)s chat: %(message)s")

    from . import mellor_adapter

    key: str | None = None
    llm: Any = None
    while True:
        try:
            request = conn.recv()
        except (EOFError, OSError):
            return
        if request is None:
            return
        try:
            if request["backend"] != key:
                llm = None
                gc.collect()
                llm = mellor_adapter.make_chat_backend(request["backend"])
                key = request["backend"]
            answer = mellor_adapter.chat(Path(request["run_dir"]), request["message"], llm=llm)
            reply: dict[str, Any] = {"ok": True, "answer": answer}
        except Exception as exc:
            logging.exception("the question could not be answered")
            reply = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
        try:
            conn.send(reply)
        except (OSError, ValueError):
            return
