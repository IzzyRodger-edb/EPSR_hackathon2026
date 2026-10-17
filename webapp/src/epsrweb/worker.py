"""What runs in the job process.

One process per job. `run_pipeline` is CPU-bound for minutes (the EPSR engine
uses every core) and Mellor offers no cancellation, so the only way to stop a
run is to kill its process -- which is exactly why each job gets one of its own
instead of a slot in a thread pool (integration.md section 6).

The process communicates two ways: events appended to `events.jsonl` while it
runs, and `result.json` plus the exit code when it is done. Nothing is returned
through a pipe, so the parent surviving a restart loses nothing.
"""

from __future__ import annotations

import json
import logging
import os
import sys
import traceback
from pathlib import Path
from typing import Any

RESULT_FILE = "result.json"
WORKER_LOG = "worker.log"


def _configure_logging(run_dir: Path) -> None:
    handler = logging.FileHandler(run_dir / WORKER_LOG, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    root = logging.getLogger()
    root.handlers = [handler]
    root.setLevel(logging.INFO)


def _apply_environment(job: dict[str, Any]) -> None:
    """Set only what this job needs, in its own process.

    `RAYON_NUM_THREADS` caps the engine's thread count so several jobs can share
    a machine; `HF_TOKEN` may be a per-user key. Neither is ever logged.
    """
    if job.get("rayon_threads"):
        os.environ["RAYON_NUM_THREADS"] = str(job["rayon_threads"])
    token = job.get("hf_token")
    if token:
        os.environ["HF_TOKEN"] = token
    os.environ.setdefault("MPLBACKEND", "Agg")


def run_job(job: dict[str, Any]) -> None:
    """Entry point of the child process. Exits non-zero on failure.

    `job` must be picklable: it carries strings and paths only, never a database
    connection or an open file.
    """
    run_dir = Path(job["run_dir"])
    run_dir.mkdir(parents=True, exist_ok=True)
    _configure_logging(run_dir)
    _apply_environment(job)

    # Imported here, not at module level: on Windows `multiprocessing` re-imports
    # this module in the child, and matplotlib/pandas should be paid for once, in
    # the child that uses them.
    from . import events, mellor_adapter

    writer = events.EventWriter(run_dir)
    writer.emit({"type": "started", "backend": job.get("backend"), "quick": job.get("quick")})
    logging.info("starting run %s in %s", job.get("run_id"), run_dir)

    try:
        outcome = mellor_adapter.run_pipeline(
            run_dir,
            [Path(p) for p in job["data_files"]],
            prompt=job.get("prompt") or None,
            backend=job.get("backend") or "auto",
            model=job.get("model") or None,
            quick=bool(job.get("quick")),
            system=job.get("system") or None,
            emit=writer,
        )
    except BaseException as exc:  # including KeyboardInterrupt from a terminate()
        detail = traceback.format_exc(limit=8)
        logging.error("run failed: %s", detail)
        _write_result(run_dir, {"ok": False, "error": f"{type(exc).__name__}: {exc}",
                                "traceback": detail})
        writer.emit({"type": "failed", "error": f"{type(exc).__name__}: {exc}"})
        sys.exit(1)

    _write_result(run_dir, {"ok": True, **outcome})
    writer.emit({"type": "finished", "spec_source": outcome.get("spec_source"),
                 "demo": bool(outcome.get("demo")), "message": outcome.get("message")})
    logging.info("run %s finished", job.get("run_id"))


def _write_result(run_dir: Path, payload: dict[str, Any]) -> None:
    try:
        (run_dir / RESULT_FILE).write_text(json.dumps(payload, indent=2, default=str),
                                          encoding="utf-8")
    except OSError:
        logging.exception("could not write %s", RESULT_FILE)


def read_result(run_dir: Path) -> dict[str, Any] | None:
    try:
        with (run_dir / RESULT_FILE).open("r", encoding="utf-8") as handle:
            loaded = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    return loaded if isinstance(loaded, dict) else None
