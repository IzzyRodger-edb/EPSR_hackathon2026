"""Starting and stopping child processes together with everything they start.

A refinement worker runs the EPSR script as a subprocess of its own (Mellor's
`runner.run_spec`), and a chat can do the same through its `build_and_run` tool.
Terminating only the worker would leave that script running, using every core,
with nobody reading its output. So each child makes itself the leader of a new
process group as its first act, and is stopped by signalling the whole group.

On platforms without process groups (Windows) this falls back to signalling the
child alone, which is what `multiprocessing` does anyway.
"""

from __future__ import annotations

import os
import signal
from typing import Any

_POSIX = os.name == "posix"


def become_group_leader() -> None:
    """Put the calling process, and every process it starts, in a new group.

    Called first thing in a child. A side effect worth knowing: Ctrl+C in the
    terminal no longer reaches the child directly; the server stops it on the
    way out instead (`JobManager.stop`).
    """
    if _POSIX:
        try:
            os.setpgid(0, 0)
        except OSError:
            pass


def _signal_group(process: Any, sig: int) -> bool:
    """Signal the child's group, or the child alone if it has no group yet.

    Returns True when the group was reached. `killpg(pid)` can only reach a group
    led by this very child: a process group id is its leader's pid, and a pid is
    not reused while its group exists. A child that has not yet called
    `become_group_leader` leads no group, so the call fails and the child is
    signalled directly.
    """
    if _POSIX:
        try:
            os.killpg(process.pid, sig)
            return True
        except (ProcessLookupError, PermissionError):
            pass
    if sig == signal.SIGTERM:
        process.terminate()
    else:
        process.kill()
    return False


def stop(process: Any, *, grace: float = 10.0) -> None:
    """Terminate a `multiprocessing` child and everything it started.

    SIGTERM first, SIGKILL after `grace` seconds, then SIGKILL once more to the
    group for anything that outlived its parent.
    """
    if process.pid is None:
        return
    grouped = False
    if process.is_alive():
        grouped = _signal_group(process, signal.SIGTERM)
        process.join(grace)
    if process.is_alive():
        grouped = _signal_group(process, getattr(signal, "SIGKILL", signal.SIGTERM)) or grouped
        process.join(5.0)
    if grouped:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
