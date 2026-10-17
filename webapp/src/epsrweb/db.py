"""SQLite storage for accounts, runs and the Overleaf push log.

Mellor keeps no user state -- a run is entirely described by its folder
(integration.md §7.2) -- so this database only holds what the *app* adds:
who owns which folder, what state its job is in, and each user's Overleaf
credentials. Results are never copied into it.
"""

from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

#: Job lifecycle. `interrupted` is set at start-up for runs that were in flight
#: when the server stopped; their folder may hold a partial result.
STATUSES = ("queued", "running", "done", "failed", "cancelled", "interrupted")
TERMINAL_STATUSES = ("done", "failed", "cancelled", "interrupted")

_SCHEMA = """
CREATE TABLE IF NOT EXISTS users (
    id                    INTEGER PRIMARY KEY AUTOINCREMENT,
    email                 TEXT NOT NULL UNIQUE COLLATE NOCASE,
    password_hash         TEXT NOT NULL,
    created_at            TEXT NOT NULL,
    overleaf_token_enc    TEXT,
    overleaf_project_url  TEXT,
    overleaf_project_id   TEXT,
    hf_token_enc          TEXT
);

CREATE TABLE IF NOT EXISTS runs (
    id             TEXT PRIMARY KEY,
    user_id        INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    title          TEXT NOT NULL,
    status         TEXT NOT NULL,
    prompt         TEXT,
    backend        TEXT NOT NULL,
    model          TEXT,
    quick          INTEGER NOT NULL DEFAULT 0,
    system         TEXT,
    spec_source    TEXT,
    final_message  TEXT,
    error          TEXT,
    demo           INTEGER NOT NULL DEFAULT 0,
    created_at     TEXT NOT NULL,
    started_at     TEXT,
    finished_at    TEXT
);
CREATE INDEX IF NOT EXISTS runs_user_created ON runs(user_id, created_at DESC);

CREATE TABLE IF NOT EXISTS overleaf_pushes (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id      TEXT NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    user_id     INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    project_id  TEXT NOT NULL,
    figures     TEXT NOT NULL,
    status      TEXT NOT NULL,
    detail      TEXT,
    created_at  TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS pushes_run ON overleaf_pushes(run_id, created_at DESC);
"""


def connect(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(path, timeout=30.0)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA busy_timeout=30000")
    return conn


def init_db(path: Path) -> None:
    with connect(path) as conn:
        conn.executescript(_SCHEMA)


@contextmanager
def session(path: Path) -> Iterator[sqlite3.Connection]:
    """A connection per unit of work; commits on success, rolls back on error."""
    conn = connect(path)
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()
