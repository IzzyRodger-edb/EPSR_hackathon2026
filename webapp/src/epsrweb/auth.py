"""Accounts, sessions and the per-request current user.

Mellor stores nothing about users (integration.md section 7.2), so this is
entirely the app's. Sessions are signed cookies; there is no server-side session
table to keep clean. Each user owns one root under `runs/<user id>/`, and a
request may only ever address a run inside their own root.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from fastapi import HTTPException, Request, status
from fastapi.responses import RedirectResponse

from . import db, security, storage
from .config import Settings

SESSION_COOKIE = "epsrweb_session"


@dataclass(frozen=True)
class User:
    id: int
    email: str
    overleaf_project_url: str | None
    overleaf_project_id: str | None
    has_overleaf_token: bool
    has_hf_token: bool

    @classmethod
    def from_row(cls, row: sqlite3.Row) -> "User":
        return cls(
            id=int(row["id"]),
            email=str(row["email"]),
            overleaf_project_url=row["overleaf_project_url"],
            overleaf_project_id=row["overleaf_project_id"],
            has_overleaf_token=bool(row["overleaf_token_enc"]),
            has_hf_token=bool(row["hf_token_enc"]),
        )


class LoginRequired(Exception):
    """Raised for a page request without a session; handled as a redirect."""

    def __init__(self, next_url: str = "/") -> None:
        self.next_url = next_url


def create_user(settings: Settings, email: str, password: str) -> int:
    email = email.strip().lower()
    if "@" not in email or len(email) < 5:
        raise ValueError("Enter a valid email address.")
    problem = security.password_problem(password)
    if problem:
        raise ValueError(problem)
    with db.session(settings.db_path) as conn:
        try:
            cursor = conn.execute(
                "INSERT INTO users (email, password_hash, created_at) VALUES (?, ?, ?)",
                (email, security.hash_password(password), storage.utcnow()),
            )
        except sqlite3.IntegrityError as exc:
            raise ValueError("An account with that email address already exists.") from exc
        user_id = int(cursor.lastrowid or 0)
    settings.user_root(user_id).mkdir(parents=True, exist_ok=True)
    return user_id


def authenticate(settings: Settings, email: str, password: str) -> int | None:
    with db.session(settings.db_path) as conn:
        row = conn.execute(
            "SELECT id, password_hash FROM users WHERE email = ?", (email.strip().lower(),)
        ).fetchone()
    if row is None:
        # Spend the same time as a real verification, so the response time does
        # not reveal whether the address exists.
        security.verify_password(password, security.hash_password("placeholder"))
        return None
    if not security.verify_password(password, row["password_hash"]):
        return None
    return int(row["id"])


def load_user(settings: Settings, user_id: int) -> User | None:
    with db.session(settings.db_path) as conn:
        row = conn.execute("SELECT * FROM users WHERE id = ?", (user_id,)).fetchone()
    return User.from_row(row) if row else None


def current_user(request: Request) -> User | None:
    """The signed-in user, or None. Cached on the request."""
    if hasattr(request.state, "user"):
        return request.state.user
    settings: Settings = request.app.state.settings
    user_id = security.read_session_token(
        settings.secret_key, request.cookies.get(SESSION_COOKIE)
    )
    user = load_user(settings, user_id) if user_id else None
    request.state.user = user
    return user


def require_user(request: Request) -> User:
    """For pages: raises `LoginRequired`, which the app turns into a redirect."""
    user = current_user(request)
    if user is None:
        raise LoginRequired(request.url.path)
    return user


def require_user_api(request: Request) -> User:
    """For JSON and file endpoints: a plain 401."""
    user = current_user(request)
    if user is None:
        raise HTTPException(status.HTTP_401_UNAUTHORIZED, "Sign in first.")
    return user


def set_session(response: Any, settings: Settings, user_id: int) -> None:
    response.set_cookie(
        SESSION_COOKIE,
        security.make_session_token(settings.secret_key, user_id, settings.session_max_age),
        max_age=settings.session_max_age,
        httponly=True,
        samesite="lax",
        secure=settings.cookie_secure,
        path="/",
    )


def clear_session(response: Any) -> None:
    response.delete_cookie(SESSION_COOKIE, path="/")


def csrf_token(request: Request) -> str:
    """Bound to the session cookie, so it cannot be replayed into another one."""
    settings: Settings = request.app.state.settings
    cookie = request.cookies.get(SESSION_COOKIE)
    return security.make_csrf_token(settings.secret_key, cookie) if cookie else ""


def check_csrf(request: Request, submitted: str | None) -> None:
    settings: Settings = request.app.state.settings
    cookie = request.cookies.get(SESSION_COOKIE)
    if not security.check_csrf_token(settings.secret_key, cookie, submitted):
        raise HTTPException(status.HTTP_403_FORBIDDEN, "The form expired. Try again.")


def redirect(url: str, *, notice: str | None = None, error: str | None = None) -> RedirectResponse:
    """A redirect carrying a short message for the next page to display."""
    from urllib.parse import urlencode

    params = {k: v for k, v in (("notice", notice), ("error", error)) if v}
    target = f"{url}?{urlencode(params)}" if params else url
    return RedirectResponse(target, status_code=status.HTTP_303_SEE_OTHER)


def run_root(settings: Settings, user: User, run_id: str) -> Path:
    """The folder of one run, without checking that it exists."""
    return settings.user_root(user.id) / run_id
