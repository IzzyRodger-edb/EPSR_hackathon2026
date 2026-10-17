"""The job queue: one OS process per refinement, with real cancellation.

`run_pipeline` blocks and the engine saturates the CPU, so jobs cannot run in
the web process. Mellor has no cancellation of its own; the documented way to
stop a run is to kill the worker and mark its folder (integration.md section 6).
A process pool would not allow that -- a submitted task cannot be interrupted --
so this manager starts a `multiprocessing.Process` per job and keeps a handle on
it.

The manager owns no results. It moves rows of the `runs` table between states
and lets the worker write the run folder.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import threading
from collections import deque
from pathlib import Path
from typing import Any

from . import db, events, storage, worker
from .config import Settings

LOG = logging.getLogger(__name__)

_POLL_SECONDS = 0.25


class JobManager:
    def __init__(self, settings: Settings) -> None:
        self.settings = settings
        # "spawn" everywhere, not only on Windows: a forked copy of the web
        # process would inherit its SQLite connections and its event loop.
        self._ctx = mp.get_context("spawn")
        self._queue: deque[dict[str, Any]] = deque()
        self._running: dict[str, Any] = {}
        self._job_meta: dict[str, dict[str, Any]] = {}
        self._cancelling: set[str] = set()
        self._lock = threading.RLock()
        self._wake = threading.Event()
        self._stopping = threading.Event()
        self._thread: threading.Thread | None = None

    # ------------------------------------------------------------------ lifecycle

    def start(self) -> None:
        self.recover()
        self._thread = threading.Thread(target=self._loop, name="epsrweb-jobs", daemon=True)
        self._thread.start()

    def stop(self, *, terminate_running: bool = False) -> None:
        self._stopping.set()
        self._wake.set()
        if terminate_running:
            with self._lock:
                for run_id, process in list(self._running.items()):
                    self._kill(run_id, process)
        if self._thread:
            self._thread.join(timeout=5.0)

    def recover(self) -> None:
        """Mark runs that were in flight when the server stopped.

        Their folders may hold a partial result, which is worth keeping and worth
        labelling, so they become `interrupted` rather than `failed`.
        """
        with db.session(self.settings.db_path) as conn:
            rows = conn.execute(
                "SELECT id FROM runs WHERE status IN ('queued', 'running')"
            ).fetchall()
            for row in rows:
                conn.execute(
                    "UPDATE runs SET status='interrupted', finished_at=?, "
                    "error='the server stopped while this run was in progress' WHERE id=?",
                    (storage.utcnow(), row["id"]),
                )
        if rows:
            LOG.warning("marked %d run(s) as interrupted at start-up", len(rows))

    # -------------------------------------------------------------------- queueing

    def submit(self, job: dict[str, Any]) -> None:
        """Queue a job whose `runs` row already exists with status `queued`."""
        events.EventWriter(Path(job["run_dir"])).emit(
            {"type": "queued", "content": "waiting for a free worker"}
        )
        with self._lock:
            self._queue.append(job)
        self._wake.set()

    def cancel(self, run_id: str) -> bool:
        """Stop a queued or running job. True if anything was stopped."""
        with self._lock:
            for job in list(self._queue):
                if job["run_id"] == run_id:
                    self._queue.remove(job)
                    self._finish(run_id, "cancelled", error="cancelled before it started")
                    return True
            process = self._running.get(run_id)
            if process is not None:
                self._cancelling.add(run_id)
                self._kill(run_id, process)
                return True
        return False

    def queue_position(self, run_id: str) -> int | None:
        with self._lock:
            for index, job in enumerate(self._queue, start=1):
                if job["run_id"] == run_id:
                    return index
        return None

    def active_for_user(self, user_id: int) -> int:
        """Queued plus running jobs of one user, for the per-user cap."""
        with self._lock:
            queued = sum(1 for job in self._queue if job.get("user_id") == user_id)
            running = sum(
                1 for run_id in self._running
                if self._job_meta.get(run_id, {}).get("user_id") == user_id
            )
        return queued + running

    # ------------------------------------------------------------------- scheduler

    def _loop(self) -> None:
        while not self._stopping.is_set():
            try:
                self._reap()
                self._dispatch()
            except Exception:
                LOG.exception("job scheduler error")
            self._wake.wait(_POLL_SECONDS)
            self._wake.clear()

    def _dispatch(self) -> None:
        with self._lock:
            while self._queue and len(self._running) < self.settings.max_workers:
                self._start(self._queue.popleft())

    def _start(self, job: dict[str, Any]) -> None:
        run_id = job["run_id"]
        process = self._ctx.Process(
            target=worker.run_job, args=(job,), name=f"epsr-run-{run_id}", daemon=False
        )
        try:
            process.start()
        except Exception as exc:
            LOG.exception("could not start a worker for %s", run_id)
            self._finish(run_id, "failed", error=f"could not start a worker: {exc}")
            return
        self._running[run_id] = process
        self._job_meta[run_id] = job
        with db.session(self.settings.db_path) as conn:
            conn.execute(
                "UPDATE runs SET status='running', started_at=? WHERE id=?",
                (storage.utcnow(), run_id),
            )
        LOG.info("run %s started as pid %s", run_id, process.pid)

    def _reap(self) -> None:
        finished = False
        with self._lock:
            for run_id, process in list(self._running.items()):
                if process.is_alive():
                    continue
                process.join(timeout=1.0)
                del self._running[run_id]
                job = self._job_meta.pop(run_id, {})
                self._settle(run_id, process.exitcode, Path(job.get("run_dir", ".")))
                finished = True
        if finished:
            self._wake.set()

    def _settle(self, run_id: str, exitcode: int | None, run_dir: Path) -> None:
        cancelled = run_id in self._cancelling
        self._cancelling.discard(run_id)
        result = worker.read_result(run_dir)

        if cancelled:
            self._finish(run_id, "cancelled", error="cancelled while running")
            events.EventWriter(run_dir).emit(
                {"type": "cancelled", "content": "stopped by the user"}
            )
            return
        if result and result.get("ok"):
            self._finish(run_id, "done", result=result)
            return
        error = (result or {}).get("error") or (
            f"the worker exited with code {exitcode} and wrote no result"
        )
        self._finish(run_id, "failed", error=error)

    def _kill(self, run_id: str, process: Any) -> None:
        LOG.info("terminating the worker of run %s (pid %s)", run_id, process.pid)
        process.terminate()
        process.join(timeout=10.0)
        if process.is_alive():
            process.kill()
            process.join(timeout=5.0)

    def _finish(self, run_id: str, status: str, *, result: dict | None = None,
                error: str | None = None) -> None:
        with db.session(self.settings.db_path) as conn:
            conn.execute(
                "UPDATE runs SET status=?, finished_at=?, error=?, "
                "spec_source=COALESCE(?, spec_source), "
                "final_message=COALESCE(?, final_message), "
                "demo=COALESCE(?, demo) WHERE id=?",
                (
                    status,
                    storage.utcnow(),
                    error,
                    (result or {}).get("spec_source"),
                    (result or {}).get("message"),
                    1 if (result or {}).get("demo") else None,
                    run_id,
                ),
            )
        LOG.info("run %s -> %s", run_id, status)
