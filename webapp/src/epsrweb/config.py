"""Configuration, read once from the environment.

Everything the app needs to know about the machine it runs on is here. The
defaults are chosen so that `uv run epsrweb` works on a laptop with no
environment set at all; everything that matters for a real deployment
(`EPSRWEB_SECRET_KEY`, the worker count, the data root) is overridable.
"""

from __future__ import annotations

import os
import secrets
from dataclasses import dataclass
from pathlib import Path

#: Where the app keeps its state when `EPSRWEB_DATA_ROOT` is unset: `webapp/var`.
_DEFAULT_ROOT = Path(__file__).resolve().parents[2] / "var"


def _flag(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def _int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if not raw:
        return default
    try:
        return int(raw)
    except ValueError:
        return default


@dataclass(frozen=True)
class Settings:
    data_root: Path
    secret_key: bytes
    #: Concurrently running refinements. One EPSR run uses every core
    #: (integration.md §6), so more than one worker per machine only pays off
    #: with `RAYON_NUM_THREADS` set as well.
    max_workers: int
    #: Per-user cap on simultaneous runs, enforced when a run is submitted.
    max_runs_per_user: int
    max_upload_bytes: int
    max_files_per_run: int
    rayon_threads: int | None
    default_backend: str
    #: The hosted model (`hf-api`), and the model run on this machine's GPU
    #: (`transformers`, and `auto` when torch is installed). The local default is
    #: the small model that fits a 12 GB card; see `mellor_adapter.choose_model`.
    default_model: str
    default_local_model: str
    #: Seconds the chat process keeps its model loaded after the last question.
    chat_idle_seconds: int
    allow_registration: bool
    cookie_secure: bool
    session_max_age: int
    #: Set `EPSRWEB_DISABLE_OVERLEAF_PUSH=1` to keep the Overleaf button on the
    #: zip fallback only (useful for a shared demo machine).
    overleaf_push_enabled: bool
    git_timeout: int

    @property
    def db_path(self) -> Path:
        return self.data_root / "epsrweb.db"

    @property
    def runs_root(self) -> Path:
        return self.data_root / "runs"

    def user_root(self, user_id: int) -> Path:
        """Per-user run root. This is the trust boundary (integration.md §7.2)."""
        return self.runs_root / str(user_id)


def _resolve_secret(data_root: Path) -> bytes:
    """Secret for session cookies and token encryption.

    From `EPSRWEB_SECRET_KEY` if set. Otherwise generated once and stored in
    `<data_root>/secret.key` with 0600 where the platform honours it, so that a
    development server keeps its logins and its stored Overleaf tokens across
    restarts. A deployment should set the variable.
    """
    env = os.environ.get("EPSRWEB_SECRET_KEY")
    if env:
        return env.encode("utf-8")
    key_file = data_root / "secret.key"
    if key_file.exists():
        return key_file.read_bytes().strip()
    data_root.mkdir(parents=True, exist_ok=True)
    key = secrets.token_hex(32).encode("ascii")
    key_file.write_bytes(key)
    try:
        key_file.chmod(0o600)
    except OSError:
        pass
    return key


def load_settings(data_root: Path | None = None) -> Settings:
    root = Path(os.environ.get("EPSRWEB_DATA_ROOT") or data_root or _DEFAULT_ROOT).resolve()
    root.mkdir(parents=True, exist_ok=True)
    return Settings(
        data_root=root,
        secret_key=_resolve_secret(root),
        max_workers=max(1, _int("EPSRWEB_WORKERS", 1)),
        max_runs_per_user=max(1, _int("EPSRWEB_MAX_RUNS_PER_USER", 1)),
        max_upload_bytes=_int("EPSRWEB_MAX_UPLOAD_BYTES", 8 * 1024 * 1024),
        max_files_per_run=_int("EPSRWEB_MAX_FILES", 6),
        rayon_threads=_int("RAYON_NUM_THREADS", 0) or None,
        default_backend=os.environ.get("MELLOR_BACKEND", "auto"),
        default_model=os.environ.get("MELLOR_MODEL", "Qwen/Qwen3-235B-A22B-Instruct-2507"),
        default_local_model=os.environ.get("EPSRWEB_LOCAL_MODEL", "Qwen/Qwen3-4B-Instruct-2507"),
        chat_idle_seconds=max(0, _int("EPSRWEB_CHAT_IDLE_SECONDS", 600)),
        allow_registration=_flag("EPSRWEB_ALLOW_REGISTRATION", True),
        cookie_secure=_flag("EPSRWEB_COOKIE_SECURE", False),
        session_max_age=_int("EPSRWEB_SESSION_MAX_AGE", 14 * 24 * 3600),
        overleaf_push_enabled=not _flag("EPSRWEB_DISABLE_OVERLEAF_PUSH", False),
        git_timeout=_int("EPSRWEB_GIT_TIMEOUT", 120),
    )
