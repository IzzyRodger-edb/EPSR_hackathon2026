"""Follow-up questions, answered in one long-lived process that keeps its model.

A chat used to run Mellor inside the web server and build a new backend for
every message. With the local model that meant loading 8 GB of weights into the
server process for each question, and two questions at once meant two copies
on the GPU. Here instead:

- every question goes to a single chat process, started on the first question
  and kept running, so the model is loaded once and reused;
- one question is answered at a time; a second gets `ChatBusy` straight away;
- the local model is never on the GPU twice: a question that needs it is refused
  while a refinement holds it, and a refinement that needs it waits for the
  question in flight and then has the chat process stopped (`release`);
- after `chat_idle_seconds` without a question the process is stopped, which
  hands its RAM and GPU memory back.

The job manager owns this service and shares its lock with it, so "is the model
free?" and "take the model" are one step for both sides.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import threading
import time
from pathlib import Path
from typing import Any, Callable

from . import mellor_adapter, processes
from .config import Settings

LOG = logging.getLogger(__name__)


class ChatBusy(RuntimeError):
    """The question cannot be taken now; the message says why. Try again later."""


class ChatFailed(RuntimeError):
    """The chat process could not answer."""


class ChatService:
    def __init__(self, settings: Settings, *, lock: threading.RLock,
                 model_claimed: Callable[[], bool]) -> None:
        """`model_claimed()` is called under `lock` and is true while a refinement
        runs on the local model or waits for it."""
        self.settings = settings
        self._lock = lock
        self._model_claimed = model_claimed
        self._ctx = mp.get_context("spawn")
        self._process: Any = None
        self._conn: Any = None
        #: The backend string the process last answered with, and whether it is
        #: a model on this machine (and so holds GPU memory once loaded).
        self._key: str | None = None
        self._key_local = False
        self._busy = False
        self._last_used = 0.0

    # ------------------------------------------------------------------ asking

    def ask(self, run_dir: Path, message: str, *, backend: str, model: str | None) -> str:
        """Answer one question about a finished run. Blocks until it is answered."""
        if not mellor_adapter.mellor_status()["available"]:
            # Demonstration mode: a summary of the folder, no model involved.
            return mellor_adapter.chat(run_dir, message)

        key = mellor_adapter.backend_string(backend, model)
        local = mellor_adapter.uses_local_model(backend)
        with self._lock:
            if self._busy:
                raise ChatBusy("Another question is being answered. Ask again when it is done.")
            if local and self._model_claimed():
                raise ChatBusy("A refinement is using the local model. Ask again when it "
                               "has finished.")
            self._busy = True
        conn = None
        try:
            conn = self._connect(key, local)
            conn.send({"run_dir": str(run_dir), "message": message, "backend": key})
            reply = conn.recv()
        except (EOFError, OSError) as exc:
            # The process died, or was stopped for shutdown. This thread holds the
            # connection while a question is in flight, so it closes it.
            if conn is not None:
                conn.close()
            with self._lock:
                if self._conn is conn:
                    self._discard_locked()
            raise ChatFailed("the chat process stopped before it answered "
                             f"({type(exc).__name__})") from None
        finally:
            with self._lock:
                self._busy = False
                self._last_used = time.monotonic()
        if not reply.get("ok"):
            raise ChatFailed(reply.get("error") or "the chat process gave no answer")
        return str(reply.get("answer") or "")

    def _connect(self, key: str, local: bool) -> Any:
        """The connection to a chat process ready for `key`, starting one if needed.

        Switching away from a local model restarts the process, because stopping
        it is the only reliable way to hand back its GPU memory.
        """
        with self._lock:
            # The asking thread owns the connection here, so it may close it.
            if self._process is not None and not self._process.is_alive():
                self._discard_locked(close=True)
            if self._process is not None and self._key != key and self._key_local:
                LOG.info("chat: switching from %s to %s, restarting the chat process",
                         self._key, key)
                self._discard_locked(close=True)
            if self._process is None:
                from . import chat_worker

                parent, child = self._ctx.Pipe()
                process = self._ctx.Process(target=chat_worker.serve, args=(child,),
                                            name="epsr-chat", daemon=True)
                process.start()
                child.close()  # so that recv() sees EOF if the process dies
                self._process, self._conn = process, parent
                LOG.info("chat process started as pid %s", process.pid)
            self._key, self._key_local = key, local
            return self._conn

    # ---------------------------------------------------- called by the job manager

    @property
    def busy(self) -> bool:
        return self._busy

    def release(self) -> bool:
        """Free the local model for a refinement. False while a question needs it.

        Called by the job manager, under the shared lock, before it starts a run
        that will load the local model.
        """
        with self._lock:
            if self._busy and self._key_local:
                return False
            if not self._busy and self._key_local and self._process is not None:
                LOG.info("chat: stopping the chat process to free the model for a refinement")
                self._discard_locked()
            return True

    def reap_idle(self) -> None:
        """Stop the chat process once it has been idle for `chat_idle_seconds`."""
        with self._lock:
            if self._process is None or self._busy:
                return
            if time.monotonic() - self._last_used >= self.settings.chat_idle_seconds:
                LOG.info("chat: idle for %d s, stopping the chat process",
                         self.settings.chat_idle_seconds)
                self._discard_locked()

    def stop(self) -> None:
        """Stop the chat process, answering or not. Used at shutdown."""
        with self._lock:
            self._discard_locked()

    # ----------------------------------------------------------------- internals

    def _discard_locked(self, *, close: bool | None = None) -> None:
        """Stop the process and forget it. Call with the lock held.

        By default the connection is closed only when no question is in flight;
        while one is, the asking thread is blocked reading it and closes it itself
        once the read fails. `close=True` is for the asking thread.
        """
        if close is None:
            close = not self._busy
        process, conn = self._process, self._conn
        self._process = self._conn = None
        self._key, self._key_local = None, False
        if process is not None:
            processes.stop(process, grace=5.0)
        if conn is not None and close:
            conn.close()
