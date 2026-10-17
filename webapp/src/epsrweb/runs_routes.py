"""Everything under /runs: create, watch, read, publish, delete.

The HTTP shape follows the minimal mapping in integration.md section 7, with
the pages the brief asks for added on top.
"""

from __future__ import annotations

import asyncio
import io
import json
import logging
import os
import shutil
import sqlite3
import zipfile
from pathlib import Path

from fastapi import APIRouter, File, Form, HTTPException, Request, UploadFile, status
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response, StreamingResponse
from starlette.concurrency import run_in_threadpool

from . import auth, db, events, mellor_adapter, overleaf, runviews, security, storage
from .chat import ChatBusy, ChatFailed

LOG = logging.getLogger(__name__)
router = APIRouter()

_READ_CHUNK = 256 * 1024


# ----------------------------------------------------------------------- helpers


def _load_run(request: Request, user: auth.User, run_id: str) -> tuple[sqlite3.Row, Path]:
    """The run row and its folder, or 404 if it is not this user's."""
    settings = request.app.state.settings
    with db.session(settings.db_path) as conn:
        row = conn.execute(
            "SELECT * FROM runs WHERE id = ? AND user_id = ?", (run_id, user.id)
        ).fetchone()
    if row is None:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "No such run.")
    return row, auth.run_root(settings, user, run_id)


async def _save_uploads(
    files: list[UploadFile], target: Path, *, max_bytes: int, max_files: int
) -> tuple[list[Path], list[storage.DataFileInfo]]:
    """Write the uploads, rejecting anything that is not a usable data file.

    Validation happens here rather than in the worker so that a mistake costs a
    form redisplay instead of a queued job (integration.md section 6).
    """
    real = [f for f in files if f and f.filename]
    if not real:
        raise ValueError("Upload at least one .mint01 data file.")
    if len(real) > max_files:
        raise ValueError(f"At most {max_files} files per run; {len(real)} were given.")

    target.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []
    infos: list[storage.DataFileInfo] = []
    used: set[str] = set()

    for upload in real:
        name = storage.safe_data_filename(upload.filename or "data.mint01")
        stem, _, ext = name.rpartition(".")
        index = 2
        while name in used:
            name = f"{stem}-{index}.{ext}"
            index += 1
        used.add(name)

        path = target / name
        written = 0
        with path.open("wb") as handle:
            while True:
                chunk = await upload.read(_READ_CHUNK)
                if not chunk:
                    break
                written += len(chunk)
                if written > max_bytes:
                    handle.close()
                    path.unlink(missing_ok=True)
                    raise ValueError(
                        f"{upload.filename} is larger than "
                        f"{storage.human_size(max_bytes)}; a .mint01 file is normally "
                        "100 kB to 1 MB."
                    )
                handle.write(chunk)
        await upload.close()

        text = path.read_text(encoding="utf-8", errors="replace")
        try:
            infos.append(storage.validate_mint01(text, name=name))
        except ValueError:
            path.unlink(missing_ok=True)
            raise
        saved.append(path)
    return saved, infos


# ------------------------------------------------------------------ create a run


@router.post("/runs")
async def create_run(
    request: Request,
    files: list[UploadFile] = File(default=[]),
    title: str = Form(default=""),
    prompt: str = Form(default=""),
    backend: str = Form(default="auto"),
    model: str = Form(default=""),
    system: str = Form(default=""),
    quick: str = Form(default=""),
    csrf: str = Form(default=""),
) -> Response:
    user = auth.require_user(request)
    auth.check_csrf(request, csrf)
    settings = request.app.state.settings
    jobs = request.app.state.jobs

    if jobs.active_for_user(user.id) >= settings.max_runs_per_user:
        return auth.redirect(
            "/new",
            error=f"You already have {settings.max_runs_per_user} run(s) in progress. "
                  "One EPSR run uses every core of the machine, so wait for it to finish "
                  "or cancel it.",
        )
    # `make_backend` raises at the point of use, which would cost a queued job and a
    # failed run to discover. Say it in the form instead.
    problem = mellor_adapter.backend_problem(
        backend,
        user_token=user.has_hf_token,
        server_token=bool(os.environ.get("HF_TOKEN")),
    )
    if problem:
        return auth.redirect("/new", error=problem)

    run_id = storage.new_run_id()
    run_dir = settings.user_root(user.id) / run_id
    try:
        saved, infos = await _save_uploads(
            files, run_dir / "uploads",
            max_bytes=settings.max_upload_bytes, max_files=settings.max_files_per_run,
        )
    except ValueError as exc:
        shutil.rmtree(run_dir, ignore_errors=True)
        return auth.redirect("/new", error=str(exc))

    run_title = title.strip() or (infos[0].title if infos else run_id)
    model_id = mellor_adapter.choose_model(settings, backend, model)
    with db.session(settings.db_path) as conn:
        conn.execute(
            "INSERT INTO runs (id, user_id, title, status, prompt, backend, model, quick, "
            "system, created_at) VALUES (?, ?, ?, 'queued', ?, ?, ?, ?, ?, ?)",
            (run_id, user.id, run_title[:200], prompt.strip() or None, backend,
             model_id, 1 if quick else 0,
             system.strip() or None, storage.utcnow()),
        )

    writer = events.EventWriter(run_dir)
    for info in infos:
        writer.emit({"type": "progress",
                     "line": f"{info.name}: '{info.title}', {info.points} points, "
                             f"Q = {info.q_min:.3f} to {info.q_max:.2f} 1/A"})
        for warning in info.warnings:
            writer.emit({"type": "note", "content": f"{info.name}: {warning}"})

    hf_token = None
    if user.has_hf_token:
        with db.session(settings.db_path) as conn:
            row = conn.execute(
                "SELECT hf_token_enc FROM users WHERE id = ?", (user.id,)
            ).fetchone()
        if row and row["hf_token_enc"]:
            hf_token = security.decrypt_secret(settings.secret_key, row["hf_token_enc"])

    jobs.submit({
        "run_id": run_id,
        "user_id": user.id,
        "run_dir": str(run_dir),
        "data_files": [str(p) for p in saved],
        "prompt": prompt.strip(),
        "backend": backend,
        "model": model_id,
        "system": system.strip(),
        "quick": bool(quick),
        "rayon_threads": settings.rayon_threads,
        "hf_token": hf_token,
    })
    return auth.redirect(f"/runs/{run_id}")


# ------------------------------------------------------------------- the run page


@router.get("/runs/{run_id}", response_class=HTMLResponse)
async def run_page(request: Request, run_id: str) -> Response:
    user = auth.require_user(request)
    row, run_dir = _load_run(request, user, run_id)
    settings = request.app.state.settings
    jobs = request.app.state.jobs
    templates = request.app.state.templates

    return templates.TemplateResponse(
        request,
        "run.html",
        {
            "run": runviews.run_summary(row, run_dir),
            "figures": runviews.list_figures(run_dir),
            "analysis": runviews.analysis_view(run_dir),
            "downloads": runviews.downloads(run_dir),
            "log_lines": runviews.log_lines(run_dir),
            "discussion": (run_dir / "discussion.md").read_text(encoding="utf-8", errors="replace")
            if (run_dir / "discussion.md").is_file() else None,
            "queue_position": jobs.queue_position(run_id),
            "uploads": sorted(p.name for p in (run_dir / "uploads").glob("*"))
            if (run_dir / "uploads").is_dir() else [],
            "has_epsr_folder": (run_dir / "epsr").is_dir(),
            "overleaf_enabled": settings.overleaf_push_enabled,
        },
    )


@router.get("/runs/{run_id}/status")
async def run_status(request: Request, run_id: str) -> JSONResponse:
    user = auth.require_user_api(request)
    row, run_dir = _load_run(request, user, run_id)
    jobs = request.app.state.jobs
    return JSONResponse({
        "id": row["id"],
        "status": row["status"],
        "error": row["error"],
        "demo": bool(row["demo"]),
        "spec_source": row["spec_source"],
        "message": row["final_message"],
        "queue_position": jobs.queue_position(run_id),
        "figures": [f.key for f in runviews.list_figures(run_dir)],
    })


@router.get("/runs/{run_id}/events")
async def run_events(request: Request, run_id: str) -> StreamingResponse:
    """Mellor's `on_event` stream as Server-Sent Events.

    Stops once the job is in a terminal state and the file has been drained, so
    a finished run does not hold a connection open.
    """
    user = auth.require_user_api(request)
    row, run_dir = _load_run(request, user, run_id)
    settings = request.app.state.settings

    async def stream():
        seq = 0
        idle = 0
        while True:
            if await request.is_disconnected() or request.app.state.shutting_down.is_set():
                return
            new = await run_in_threadpool(events.read_events, run_dir, seq)
            for event in new:
                seq = int(event.get("seq", seq))
                payload = {"seq": seq, "type": event.get("type"),
                           "line": events.summarise(event)}
                yield f"data: {json.dumps(payload)}\n\n"
            with db.session(settings.db_path) as conn:
                current = conn.execute(
                    "SELECT status FROM runs WHERE id = ?", (run_id,)
                ).fetchone()
            state = current["status"] if current else "failed"
            if state in db.TERMINAL_STATUSES and not new:
                yield f"event: end\ndata: {json.dumps({'status': state})}\n\n"
                return
            idle += 1
            if idle % 30 == 0:
                yield ": keep-alive\n\n"
            await asyncio.sleep(0.5)

    return StreamingResponse(stream(), media_type="text/event-stream", headers={
        "Cache-Control": "no-cache", "X-Accel-Buffering": "no", "Connection": "keep-alive",
    })


# ------------------------------------------------------------------ serving files


def _servable(run_dir: Path, relative: str) -> Path:
    """Resolve a request path inside the run folder, or refuse it.

    Only the folders and files in `runviews` are reachable; `trajectory/` and
    `events.jsonl` are not, the first because of its size and the second because
    the page already streams it.
    """
    try:
        path = storage.resolve_within(run_dir, relative)
    except ValueError:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, "Bad path.")
    parts = Path(relative).as_posix().split("/")
    allowed = (
        parts[0] in runviews.SERVABLE_DIRS
        or relative in runviews.SERVABLE_FILES
        or (len(parts) == 1 and parts[0].startswith("run_") and parts[0].endswith(".py"))
    )
    if not allowed or not path.is_file():
        raise HTTPException(status.HTTP_404_NOT_FOUND, "No such file in this run.")
    return path


@router.get("/runs/{run_id}/files/{relative:path}")
async def run_file(request: Request, run_id: str, relative: str) -> FileResponse:
    user = auth.require_user_api(request)
    _row, run_dir = _load_run(request, user, run_id)
    path = _servable(run_dir, relative)
    # Figures are displayed inline; everything else downloads.
    inline = path.suffix.lower() in {".png", ".pdf", ".svg"}
    return FileResponse(
        path,
        filename=None if inline else path.name,
        headers={"Cache-Control": "private, max-age=60"},
    )


@router.get("/runs/{run_id}/epsr.zip")
async def epsr_archive(request: Request, run_id: str) -> Response:
    """The EPSR simulation folder: inputs, .ato, .pcof and every output file."""
    user = auth.require_user_api(request)
    _row, run_dir = _load_run(request, user, run_id)
    folder = run_dir / "epsr"
    if not folder.is_dir():
        raise HTTPException(status.HTTP_404_NOT_FOUND, "This run has no epsr/ folder.")
    size = storage.dir_size(folder)
    if size > 250 * 1024 * 1024:
        raise HTTPException(
            status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            f"The folder is {storage.human_size(size)}; copy it off the server instead.",
        )

    def build() -> bytes:
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
            for item in sorted(folder.rglob("*")):
                if item.is_file():
                    archive.write(item, f"epsr/{item.relative_to(folder).as_posix()}")
        return buffer.getvalue()

    payload = await run_in_threadpool(build)
    return Response(payload, media_type="application/zip", headers={
        "Content-Disposition": f'attachment; filename="{run_id}-epsr.zip"'
    })


# ------------------------------------------------------------- cancel and delete


@router.post("/runs/{run_id}/cancel")
async def cancel_run(request: Request, run_id: str, csrf: str = Form(default="")) -> Response:
    user = auth.require_user(request)
    auth.check_csrf(request, csrf)
    _load_run(request, user, run_id)
    stopped = request.app.state.jobs.cancel(run_id)
    return auth.redirect(
        f"/runs/{run_id}",
        notice="The run was stopped." if stopped else None,
        error=None if stopped else "That run is not queued or running.",
    )


@router.post("/runs/{run_id}/delete")
async def delete_run(request: Request, run_id: str, csrf: str = Form(default="")) -> Response:
    user = auth.require_user(request)
    auth.check_csrf(request, csrf)
    row, run_dir = _load_run(request, user, run_id)
    if row["status"] in {"queued", "running"}:
        request.app.state.jobs.cancel(run_id)
    settings = request.app.state.settings
    # Re-resolve under the user's root before removing anything.
    try:
        target = storage.resolve_within(settings.user_root(user.id), run_id)
    except ValueError:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, "Bad run id.")
    shutil.rmtree(target, ignore_errors=True)
    with db.session(settings.db_path) as conn:
        conn.execute("DELETE FROM runs WHERE id = ? AND user_id = ?", (run_id, user.id))
    return auth.redirect("/", notice=f"Run {run_id} and its folder were deleted.")


# -------------------------------------------------------------------------- chat


@router.post("/runs/{run_id}/chat")
async def chat(request: Request, run_id: str) -> JSONResponse:
    user = auth.require_user_api(request)
    body = await request.json()
    auth.check_csrf(request, request.headers.get("x-csrf-token"))
    message = str(body.get("message") or "").strip()
    if not message:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, "Ask something.")
    row, run_dir = _load_run(request, user, run_id)
    if row["status"] in {"queued", "running"}:
        raise HTTPException(status.HTTP_409_CONFLICT, "Wait for the run to finish.")

    settings = request.app.state.settings
    backend = row["backend"] if row["backend"] != "library" else "auto"
    try:
        # Answered in the chat process, which keeps its model loaded between
        # questions and takes one question at a time (chat.py).
        answer = await run_in_threadpool(
            request.app.state.jobs.chat.ask, run_dir, message,
            backend=backend, model=mellor_adapter.choose_model(settings, backend, row["model"]),
        )
    except ChatBusy as exc:
        raise HTTPException(status.HTTP_409_CONFLICT, str(exc))
    except ChatFailed as exc:
        if request.app.state.shutting_down.is_set():
            raise HTTPException(status.HTTP_503_SERVICE_UNAVAILABLE,
                                "The server is shutting down.")
        LOG.warning("chat failed for %s: %s", run_id, exc)
        raise HTTPException(status.HTTP_502_BAD_GATEWAY, f"{type(exc).__name__}: {exc}")
    except Exception as exc:
        LOG.exception("chat failed for %s", run_id)
        raise HTTPException(status.HTTP_502_BAD_GATEWAY, f"{type(exc).__name__}: {exc}")
    return JSONResponse({"answer": answer})


# -------------------------------------------------------------- publish figures


def _selections(
    run_dir: Path, run_id: str, keys: list[str], prefer: str
) -> tuple[list[tuple[Path, str, str]], list[str]]:
    """Resolve figure keys to (source, name in figures/, caption).

    Names are prefixed with the run id so that figures from two runs of the same
    system cannot overwrite each other in the Overleaf project
    (integration.md section 7.1). PDF is preferred for LaTeX, PNG is the
    fallback: Mellor's own figures are PNG only.
    """
    available = {figure.key: figure for figure in runviews.list_figures(run_dir)}
    selections: list[tuple[Path, str, str]] = []
    missing: list[str] = []
    for key in keys:
        figure = available.get(key)
        if figure is None:
            missing.append(key)
            continue
        relative = (figure.pdf or figure.png) if prefer == "pdf" else (figure.png or figure.pdf)
        if not relative:
            missing.append(key)
            continue
        source = run_dir / relative
        if not source.is_file():
            missing.append(key)
            continue
        selections.append((source, f"{run_id}_{key}{source.suffix}", figure.caption))
    return selections, missing


@router.post("/runs/{run_id}/overleaf")
async def push_to_overleaf(request: Request, run_id: str) -> JSONResponse:
    """Commit the selected figures to `figures/` in the user's Overleaf project."""
    user = auth.require_user_api(request)
    body = await request.json()
    auth.check_csrf(request, request.headers.get("x-csrf-token"))
    settings = request.app.state.settings
    if not settings.overleaf_push_enabled:
        raise HTTPException(status.HTTP_403_FORBIDDEN,
                            "Pushing to Overleaf is disabled on this server; "
                            "download the zip instead.")

    row, run_dir = _load_run(request, user, run_id)
    keys = [str(k) for k in (body.get("figures") or [])]
    if not keys:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, "Select at least one figure.")
    prefer = "pdf" if str(body.get("format", "pdf")) == "pdf" else "png"

    with db.session(settings.db_path) as conn:
        account = conn.execute(
            "SELECT overleaf_token_enc, overleaf_project_id FROM users WHERE id = ?",
            (user.id,),
        ).fetchone()
    token = security.decrypt_secret(settings.secret_key, account["overleaf_token_enc"]) \
        if account and account["overleaf_token_enc"] else None
    if not token:
        raise HTTPException(status.HTTP_400_BAD_REQUEST,
                            "No Overleaf token is stored. Add one in Settings.")
    try:
        project_id = overleaf.parse_project_id(
            str(body.get("project") or account["overleaf_project_id"] or "")
        )
    except overleaf.OverleafError as exc:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, str(exc))

    selections, missing = _selections(run_dir, run_id, keys, prefer)
    if not selections:
        raise HTTPException(status.HTTP_400_BAD_REQUEST,
                            "None of the selected figures exist on disk.")

    message = (f"Add EPSR figures from run {run_id}"
               + (f" ({row['title']})" if row["title"] else ""))
    try:
        result = await run_in_threadpool(
            overleaf.push_figures, project_id, token, selections,
            commit_message=message, author_name="EPSR web app",
            author_email=user.email, timeout=settings.git_timeout,
        )
        outcome, detail = "ok", result.detail
    except overleaf.OverleafError as exc:
        outcome, detail = "error", str(exc)
        result = None
    except Exception as exc:
        LOG.exception("overleaf push failed for %s", run_id)
        outcome, detail = "error", f"{type(exc).__name__}: {exc}"
        result = None

    with db.session(settings.db_path) as conn:
        conn.execute(
            "INSERT INTO overleaf_pushes (run_id, user_id, project_id, figures, status, "
            "detail, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)",
            (run_id, user.id, project_id, ", ".join(name for _s, name, _c in selections),
             outcome, detail[:2000], storage.utcnow()),
        )

    if outcome == "error":
        return JSONResponse({"ok": False, "detail": detail}, status_code=502)
    return JSONResponse({
        "ok": True,
        "detail": detail,
        "commit": result.commit if result else None,
        "pushed": result.pushed if result else [],
        "missing": missing,
        "snippets": [
            overleaf.latex_snippet(name, caption, Path(name).stem)
            for _source, name, caption in selections
        ],
    })


@router.post("/runs/{run_id}/figures.zip")
async def figures_zip(
    request: Request,
    run_id: str,
    figures: list[str] = Form(default=[]),
    figure_format: str = Form(default="pdf"),
    csrf: str = Form(default=""),
) -> Response:
    """The fallback for accounts without Overleaf's premium Git integration."""
    user = auth.require_user_api(request)
    auth.check_csrf(request, csrf)
    _row, run_dir = _load_run(request, user, run_id)
    if not figures:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, "Select at least one figure.")
    prefer = "pdf" if figure_format == "pdf" else "png"
    selections, _missing = _selections(run_dir, run_id, figures, prefer)
    if not selections:
        raise HTTPException(status.HTTP_400_BAD_REQUEST,
                            "None of the selected figures exist on disk.")
    payload = await run_in_threadpool(
        overleaf.build_zip, selections,
        project_note=f"Run {run_id}. Drag figures/ into the Overleaf file tree.",
    )
    return Response(payload, media_type="application/zip", headers={
        "Content-Disposition": f'attachment; filename="{run_id}-figures.zip"'
    })
