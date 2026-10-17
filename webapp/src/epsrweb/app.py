"""The FastAPI application: pages, accounts and wiring.

Everything to do with runs lives in `runs_routes`; this module holds the app
factory, the pages the brief asks for (About and Instructions), the account
pages, and the start-up and shutdown of the job manager.
"""

from __future__ import annotations

import logging
import os
import threading
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, Form, Request, Response, status
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from . import auth, db, jobs as jobs_module, mellor_adapter, overleaf, runs_routes, security, storage
from .config import Settings, load_settings

LOG = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent


def create_app(settings: Settings | None = None) -> FastAPI:
    settings = settings or load_settings()
    db.init_db(settings.db_path)
    settings.runs_root.mkdir(parents=True, exist_ok=True)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        app.state.jobs.start()
        LOG.info("job manager started with %d worker slot(s)", settings.max_workers)
        yield
        # Stop the running refinements and the chat process with the server,
        # rather than leave them to finish (or to hold up the exit) on their own.
        app.state.jobs.stop(terminate_running=True)

    app = FastAPI(title="EPSR / Mellor", docs_url=None, redoc_url=None, lifespan=lifespan)
    app.state.settings = settings
    app.state.jobs = jobs_module.JobManager(settings)
    #: Set by `__main__` on Ctrl+C, so that long requests end instead of holding
    #: the shutdown up.
    app.state.shutting_down = threading.Event()

    templates = Jinja2Templates(directory=str(HERE / "templates"))
    templates.env.globals.update(
        csrf_token=auth.csrf_token,
        current_user=auth.current_user,
        mellor_status=mellor_adapter.mellor_status(),
        backends=mellor_adapter.BACKENDS,
        known_systems=mellor_adapter.known_systems(),
        default_model=settings.default_model,
        default_local_model=settings.default_local_model,
        max_files=settings.max_files_per_run,
        max_upload=storage.human_size(settings.max_upload_bytes),
        overleaf_push_enabled=settings.overleaf_push_enabled,
        allow_registration=settings.allow_registration,
    )
    app.state.templates = templates

    app.mount("/static", StaticFiles(directory=str(HERE / "static")), name="static")
    app.include_router(runs_routes.router)
    _register_pages(app)
    return app


def _register_pages(app: FastAPI) -> None:

    @app.exception_handler(auth.LoginRequired)
    async def _login_required(request: Request, exc: auth.LoginRequired) -> Response:
        from urllib.parse import quote

        return auth.redirect(f"/login?next={quote(exc.next_url)}")

    @app.get("/healthz", response_class=PlainTextResponse)
    async def healthz() -> str:
        status_info = mellor_adapter.mellor_status()
        return "ok mellor={}".format("yes" if status_info["available"] else "no")

    @app.get("/", response_class=HTMLResponse)
    async def home(request: Request) -> Response:
        user = auth.current_user(request)
        templates = request.app.state.templates
        if user is None:
            return templates.TemplateResponse(request, "landing.html", {})
        settings: Settings = request.app.state.settings
        with db.session(settings.db_path) as conn:
            rows = conn.execute(
                "SELECT * FROM runs WHERE user_id = ? ORDER BY created_at DESC LIMIT 100",
                (user.id,),
            ).fetchall()
        from . import runviews

        runs = [runviews.run_summary(row, settings.user_root(user.id) / row["id"])
                for row in rows]
        return templates.TemplateResponse(request, "dashboard.html", {
            "runs": runs,
            "active": request.app.state.jobs.active_for_user(user.id),
            "max_active": settings.max_runs_per_user,
        })

    @app.get("/new", response_class=HTMLResponse)
    async def new_run(request: Request) -> Response:
        auth.require_user(request)
        return request.app.state.templates.TemplateResponse(request, "new_run.html", {})

    @app.get("/about", response_class=HTMLResponse)
    async def about(request: Request) -> Response:
        return request.app.state.templates.TemplateResponse(request, "about.html", {})

    @app.get("/instructions", response_class=HTMLResponse)
    async def instructions(request: Request) -> Response:
        return request.app.state.templates.TemplateResponse(request, "instructions.html", {})

    # ------------------------------------------------------------------ accounts

    @app.get("/login", response_class=HTMLResponse)
    async def login_form(request: Request) -> Response:
        if auth.current_user(request):
            return auth.redirect("/")
        return request.app.state.templates.TemplateResponse(request, "login.html", {})

    @app.post("/login")
    async def login(request: Request, email: str = Form(...), password: str = Form(...),
                    next: str = Form(default="/")) -> Response:
        settings: Settings = request.app.state.settings
        user_id = auth.authenticate(settings, email, password)
        if user_id is None:
            return auth.redirect("/login", error="Wrong email address or password.")
        target = next if next.startswith("/") and not next.startswith("//") else "/"
        response = auth.redirect(target)
        auth.set_session(response, settings, user_id)
        settings.user_root(user_id).mkdir(parents=True, exist_ok=True)
        return response

    @app.get("/register", response_class=HTMLResponse)
    async def register_form(request: Request) -> Response:
        settings: Settings = request.app.state.settings
        if not settings.allow_registration:
            return auth.redirect("/login", error="Registration is closed on this server.")
        return request.app.state.templates.TemplateResponse(request, "register.html", {})

    @app.post("/register")
    async def register(request: Request, email: str = Form(...), password: str = Form(...),
                       confirm: str = Form(default="")) -> Response:
        settings: Settings = request.app.state.settings
        if not settings.allow_registration:
            return auth.redirect("/login", error="Registration is closed on this server.")
        if password != confirm:
            return auth.redirect("/register", error="The two passwords do not match.")
        try:
            user_id = auth.create_user(settings, email, password)
        except ValueError as exc:
            return auth.redirect("/register", error=str(exc))
        response = auth.redirect("/new", notice="Account created. Upload your data.")
        auth.set_session(response, settings, user_id)
        return response

    @app.post("/logout")
    async def logout(request: Request, csrf: str = Form(default="")) -> Response:
        auth.check_csrf(request, csrf)
        response = auth.redirect("/")
        auth.clear_session(response)
        return response

    # ------------------------------------------------------------------ settings

    @app.get("/settings", response_class=HTMLResponse)
    async def settings_page(request: Request) -> Response:
        user = auth.require_user(request)
        settings: Settings = request.app.state.settings
        with db.session(settings.db_path) as conn:
            pushes = conn.execute(
                "SELECT * FROM overleaf_pushes WHERE user_id = ? ORDER BY created_at DESC "
                "LIMIT 20", (user.id,),
            ).fetchall()
        return request.app.state.templates.TemplateResponse(request, "settings.html", {
            "user": user,
            "pushes": [dict(row) for row in pushes],
            "server_hf_token": bool(os.environ.get("HF_TOKEN")),
        })

    @app.post("/settings/overleaf")
    async def save_overleaf(request: Request, project: str = Form(default=""),
                            token: str = Form(default=""), csrf: str = Form(default="")) -> Response:
        user = auth.require_user(request)
        auth.check_csrf(request, csrf)
        settings: Settings = request.app.state.settings

        project_id = None
        if project.strip():
            try:
                project_id = overleaf.parse_project_id(project)
            except overleaf.OverleafError as exc:
                return auth.redirect("/settings", error=str(exc))

        with db.session(settings.db_path) as conn:
            if project_id:
                conn.execute(
                    "UPDATE users SET overleaf_project_url = ?, overleaf_project_id = ? "
                    "WHERE id = ?", (project.strip(), project_id, user.id),
                )
            if token.strip():
                conn.execute(
                    "UPDATE users SET overleaf_token_enc = ? WHERE id = ?",
                    (security.encrypt_secret(settings.secret_key, token.strip()), user.id),
                )
        return auth.redirect("/settings", notice="Overleaf settings saved.")

    @app.post("/settings/overleaf/forget")
    async def forget_overleaf(request: Request, csrf: str = Form(default="")) -> Response:
        user = auth.require_user(request)
        auth.check_csrf(request, csrf)
        settings: Settings = request.app.state.settings
        with db.session(settings.db_path) as conn:
            conn.execute(
                "UPDATE users SET overleaf_token_enc = NULL WHERE id = ?", (user.id,)
            )
        return auth.redirect("/settings", notice="The stored Overleaf token was deleted. "
                                                 "Revoke it in Overleaf as well.")

    @app.post("/settings/hf")
    async def save_hf(request: Request, token: str = Form(default=""),
                      csrf: str = Form(default="")) -> Response:
        user = auth.require_user(request)
        auth.check_csrf(request, csrf)
        settings: Settings = request.app.state.settings
        value = (
            security.encrypt_secret(settings.secret_key, token.strip())
            if token.strip() else None
        )
        with db.session(settings.db_path) as conn:
            conn.execute("UPDATE users SET hf_token_enc = ? WHERE id = ?", (value, user.id))
        return auth.redirect(
            "/settings",
            notice="Hugging Face token saved." if value else "Hugging Face token deleted.",
        )

    @app.exception_handler(status.HTTP_404_NOT_FOUND)
    async def not_found(request: Request, exc: Exception) -> Response:
        if request.url.path.startswith("/runs") and "application/json" in (
            request.headers.get("accept") or ""
        ):
            return JSONResponse({"detail": "not found"}, status_code=404)
        return request.app.state.templates.TemplateResponse(
            request, "error.html",
            {"code": 404, "message": "That page or run does not exist."}, status_code=404,
        )
