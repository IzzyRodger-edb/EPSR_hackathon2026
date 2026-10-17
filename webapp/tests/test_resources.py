"""The model, the chat process and child processes: who may use what, and when.

No model is loaded here. The chat service and the job manager are driven
directly with the model's availability patched, so the rules are checked in
milliseconds: one question at a time, the local model never on the GPU twice,
the right default model for each backend, and that stopping a worker stops what
it started.
"""

from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

import pytest

from epsrweb import mellor_adapter, processes
from epsrweb.chat import ChatBusy
from epsrweb.jobs import JobManager


def _pretend(monkeypatch, *, local: bool = True) -> None:
    monkeypatch.setattr(mellor_adapter, "mellor_status", lambda: {"available": True})
    monkeypatch.setattr(mellor_adapter, "local_model_available", lambda: local)


# ----------------------------------------------------------------- default models


def test_the_local_backend_defaults_to_the_small_model(settings, monkeypatch):
    _pretend(monkeypatch)
    assert mellor_adapter.choose_model(settings, "transformers", "") == settings.default_local_model
    assert mellor_adapter.choose_model(settings, "transformers", "a/b") == "a/b"
    assert mellor_adapter.choose_model(settings, "hf-api", "") == settings.default_model
    assert mellor_adapter.choose_model(settings, "hf-api", "a/b") == "a/b"


def test_auto_takes_the_default_for_where_it_runs(settings, monkeypatch):
    # Older rows carry the hosted default even when they ran locally; auto must
    # not try to load that onto the GPU.
    _pretend(monkeypatch, local=True)
    assert (mellor_adapter.choose_model(settings, "auto", settings.default_model)
            == settings.default_local_model)
    _pretend(monkeypatch, local=False)
    assert mellor_adapter.choose_model(settings, "auto", "a/b") == settings.default_model
    assert mellor_adapter.backend_string("auto", "a/b") == "auto:a/b"


def test_which_backends_use_the_local_model(monkeypatch):
    _pretend(monkeypatch, local=True)
    assert mellor_adapter.uses_local_model("auto")
    assert mellor_adapter.uses_local_model("transformers")
    assert not mellor_adapter.uses_local_model("hf-api")
    assert not mellor_adapter.uses_local_model("library")
    _pretend(monkeypatch, local=False)
    assert not mellor_adapter.uses_local_model("auto")


# ------------------------------------------------------------------- the chat rules


def test_a_second_question_is_refused_while_one_is_answered(settings, monkeypatch, tmp_path):
    _pretend(monkeypatch)
    jobs = JobManager(settings)
    jobs.chat._busy = True
    with pytest.raises(ChatBusy, match="Another question"):
        jobs.chat.ask(tmp_path, "hi", backend="hf-api", model="a/b")


def test_a_local_question_is_refused_while_a_run_holds_the_model(settings, monkeypatch, tmp_path):
    _pretend(monkeypatch)
    jobs = JobManager(settings)
    jobs._job_meta["running"] = {"run_id": "running", "backend": "auto"}
    with pytest.raises(ChatBusy, match="refinement"):
        jobs.chat.ask(tmp_path, "hi", backend="auto", model="a/b")
    assert not jobs.chat.busy


def test_a_local_run_waits_for_the_question_in_flight(settings, monkeypatch, tmp_path):
    _pretend(monkeypatch)
    jobs = JobManager(settings)
    started: list[str] = []
    monkeypatch.setattr(jobs, "_start", lambda job: started.append(job["run_id"]))
    jobs.chat._busy, jobs.chat._key_local = True, True
    jobs._queue.append({"run_id": "next", "backend": "auto"})

    jobs._dispatch()
    assert started == [] and jobs._waiting_for_model
    # Meanwhile new local questions are refused, so the run cannot be starved.
    jobs.chat._busy = False
    with pytest.raises(ChatBusy):
        jobs.chat.ask(tmp_path, "hi", backend="auto", model="a/b")

    jobs._dispatch()
    assert started == ["next"] and not jobs._waiting_for_model


def test_a_run_without_the_local_model_does_not_wait(settings, monkeypatch):
    _pretend(monkeypatch)
    jobs = JobManager(settings)
    started: list[str] = []
    monkeypatch.setattr(jobs, "_start", lambda job: started.append(job["run_id"]))
    jobs.chat._busy, jobs.chat._key_local = True, True
    jobs._queue.append({"run_id": "hosted", "backend": "hf-api"})
    jobs._dispatch()
    assert started == ["hosted"]


def test_an_idle_chat_process_is_stopped(settings):
    class Finished:
        pid = None

        def is_alive(self) -> bool:
            return False

    jobs = JobManager(settings)
    jobs.chat._process, jobs.chat._key, jobs.chat._key_local = Finished(), "auto:a/b", True
    jobs.chat._last_used = time.monotonic() - settings.chat_idle_seconds - 1
    jobs.chat.reap_idle()
    assert jobs.chat._process is None and not jobs.chat._key_local


def test_stop_can_be_called_twice(settings):
    jobs = JobManager(settings)
    jobs.stop(terminate_running=True)
    jobs.stop(terminate_running=True)


# ------------------------------------------------------------- stopping processes


def _leader_with_a_child(pid_file: str) -> None:
    """Stands in for a worker: leads its group and runs a long subprocess."""
    processes.become_group_leader()
    child = subprocess.Popen(["sleep", "60"])
    Path(pid_file).write_text(str(child.pid))
    time.sleep(60)


def _gone(pid: int) -> bool:
    try:
        state = Path(f"/proc/{pid}/stat").read_text().split(")")[-1].split()[0]
    except OSError:
        return True
    return state == "Z"


@pytest.mark.skipif(os.name != "posix" or not Path("/proc").is_dir(),
                    reason="process groups and /proc")
def test_stopping_a_worker_stops_the_subprocess_it_started(tmp_path):
    import multiprocessing as mp

    pid_file = tmp_path / "child.pid"
    worker = mp.get_context("spawn").Process(target=_leader_with_a_child, args=(str(pid_file),))
    worker.start()
    deadline = time.time() + 30
    while not pid_file.exists() and time.time() < deadline:
        time.sleep(0.05)
    child_pid = int(pid_file.read_text())

    processes.stop(worker, grace=5.0)

    deadline = time.time() + 5
    while not _gone(child_pid) and time.time() < deadline:
        time.sleep(0.05)
    assert not worker.is_alive()
    assert _gone(child_pid), "the worker's subprocess outlived it"
