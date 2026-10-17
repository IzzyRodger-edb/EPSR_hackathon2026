"""Publishing selected figures to an Overleaf project.

Overleaf has no public upload API. The supported programmatic route is its Git
bridge, which the app drives server-side: clone the project, drop the chosen
figures into `figures/`, commit, push (integration.md section 7.1).

Two things follow from that and shape this module.

*The Git bridge is a premium Overleaf feature.* Users without it cannot be made
to work by better code, so `build_zip` produces the same figures plus a
`figures.tex` of ready-made snippets for drag-and-drop upload instead.

*One token grants write access to every project of an account.* So the token is
never passed on a command line (where it would show up in `ps`), never written
into `.git/config`, and never allowed into a log or an error page. It reaches git
through `GIT_ASKPASS` and the environment of the child process only, and
`redact` is applied to everything git prints before it is stored or shown.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import dataclass, field
from io import BytesIO
from pathlib import Path
from typing import Sequence

LOG = logging.getLogger(__name__)

GIT_HOST = "https://git.overleaf.com"
FIGURES_SUBDIR = "figures"
_TOKEN_ENV = "EPSRWEB_GIT_TOKEN"
#: Overleaf project ids are 24-character hex strings; allow a little either way
#: and nothing that could be read as a command-line option or a path.
_PROJECT_ID = re.compile(r"^[A-Za-z0-9]{8,40}$")


class OverleafError(RuntimeError):
    """Something went wrong that the user can be told about verbatim."""


@dataclass
class PushResult:
    pushed: list[str] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    commit: str | None = None
    detail: str = ""


def parse_project_id(value: str) -> str:
    """Accept a project URL, a git remote or a bare id; return the id.

    `https://www.overleaf.com/project/65f1c0ffee1234567890abcd` ->
    `65f1c0ffee1234567890abcd`.
    """
    text = (value or "").strip().rstrip("/")
    if not text:
        raise OverleafError("Give the Overleaf project URL or its id.")
    text = re.sub(r"[?#].*$", "", text)
    candidate = text.rsplit("/", 1)[-1]
    if not _PROJECT_ID.match(candidate):
        raise OverleafError(
            "That does not look like an Overleaf project. Paste the project URL, "
            "for example https://www.overleaf.com/project/65f1c0ffee1234567890abcd"
        )
    return candidate


def redact(text: str, token: str | None) -> str:
    """Remove a token, and anything shaped like one, from git's output."""
    cleaned = text or ""
    if token:
        cleaned = cleaned.replace(token, "***")
    # git sometimes echoes the remote with credentials in it
    cleaned = re.sub(r"https://[^@\s/]*:[^@\s/]*@", "https://***@", cleaned)
    return cleaned


def latex_snippet(filename: str, caption: str, label: str) -> str:
    """One ready-to-paste figure block for the user's document."""
    safe_caption = caption.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
    return (
        r"\begin{figure}[htbp]" "\n"
        r"  \centering" "\n"
        rf"  \includegraphics[width=\linewidth]{{{FIGURES_SUBDIR}/{filename}}}" "\n"
        rf"  \caption{{{safe_caption}}}" "\n"
        rf"  \label{{fig:{label}}}" "\n"
        r"\end{figure}" "\n"
    )


def build_zip(
    selections: Sequence[tuple[Path, str, str]],
    *,
    project_note: str = "",
) -> bytes:
    """A zip of the chosen figures plus `figures.tex`, for manual upload.

    `selections` is a sequence of (source file, name inside `figures/`, caption).
    """
    buffer = BytesIO()
    snippets: list[str] = []
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for source, name, caption in selections:
            archive.write(source, f"{FIGURES_SUBDIR}/{name}")
            snippets.append(latex_snippet(name, caption, Path(name).stem))
        header = (
            "% Figures exported from the EPSR / Mellor web app.\n"
            "% Upload the figures/ folder to your Overleaf project, then "
            r"\input" " this file\n"
            "% or copy the blocks you want into your document.\n"
        )
        if project_note:
            header += f"% {project_note}\n"
        archive.writestr("figures.tex", header + "\n" + "\n".join(snippets))
    return buffer.getvalue()


def _write_askpass(directory: Path) -> Path:
    """A helper git can run to obtain the password, without it touching argv.

    The remote URL carries the username (`git@`), so git asks only for the
    password and this one-liner answers it from the environment of the child.
    """
    if sys.platform.startswith("win"):
        script = directory / "askpass.bat"
        script.write_text(f"@echo %{_TOKEN_ENV}%\r\n", encoding="ascii")
        return script
    script = directory / "askpass.sh"
    script.write_text(f'#!/bin/sh\nprintf "%s" "${_TOKEN_ENV}"\n', encoding="ascii")
    script.chmod(script.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return script


def _git(args: Sequence[str], *, cwd: Path, env: dict[str, str], token: str,
         timeout: int) -> subprocess.CompletedProcess[str]:
    LOG.info("git %s", " ".join(args[:2]))
    try:
        completed = subprocess.run(
            ["git", *args], cwd=cwd, env=env, capture_output=True, text=True,
            timeout=timeout, check=False,
        )
    except FileNotFoundError as exc:
        raise OverleafError("git is not installed on the server.") from exc
    except subprocess.TimeoutExpired as exc:
        raise OverleafError(f"git {args[0]} timed out after {timeout} s.") from exc
    completed.stdout = redact(completed.stdout, token)
    completed.stderr = redact(completed.stderr, token)
    return completed


def push_figures(
    project_id: str,
    token: str,
    selections: Sequence[tuple[Path, str, str]],
    *,
    commit_message: str,
    author_name: str = "EPSR web app",
    author_email: str = "noreply@localhost",
    timeout: int = 120,
) -> PushResult:
    """Clone the project, add the chosen figures to `figures/`, commit and push.

    `selections` is a sequence of (source file, name inside `figures/`, caption).
    Names are expected to be prefixed with the run id by the caller, so two runs
    of the same system cannot overwrite each other's figures.
    """
    project_id = parse_project_id(project_id)
    if not token:
        raise OverleafError("No Overleaf token is stored for this account.")
    if not selections:
        raise OverleafError("No figures were selected.")

    remote = f"{GIT_HOST.replace('https://', 'https://git@')}/{project_id}"
    workspace = Path(tempfile.mkdtemp(prefix="epsrweb-overleaf-"))
    try:
        askpass = _write_askpass(workspace)
        env = {
            **os.environ,
            _TOKEN_ENV: token,
            "GIT_ASKPASS": str(askpass),
            "GIT_TERMINAL_PROMPT": "0",     # fail instead of hanging on a prompt
            "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_LFS_SKIP_SMUDGE": "1",
        }
        clone = workspace / "project"
        result = _git(["clone", "--depth", "1", remote, str(clone)], cwd=workspace, env=env,
                      token=token, timeout=timeout)
        if result.returncode != 0:
            raise OverleafError(_clone_hint(result.stderr))

        figures_dir = clone / FIGURES_SUBDIR
        figures_dir.mkdir(parents=True, exist_ok=True)
        pushed: list[str] = []
        for source, name, _caption in selections:
            shutil.copyfile(source, figures_dir / name)
            pushed.append(name)

        base = ["-c", f"user.name={author_name}", "-c", f"user.email={author_email}"]
        _git([*base, "add", "--", FIGURES_SUBDIR], cwd=clone, env=env, token=token,
             timeout=timeout)
        staged = _git(["diff", "--cached", "--name-only"], cwd=clone, env=env, token=token,
                      timeout=timeout)
        if not staged.stdout.strip():
            return PushResult(pushed=[], skipped=pushed, commit=None,
                              detail="Every selected figure is already in the project, "
                                     "byte for byte. Nothing was pushed.")

        commit = _git([*base, "commit", "-m", commit_message], cwd=clone, env=env, token=token,
                      timeout=timeout)
        if commit.returncode != 0:
            raise OverleafError(f"git commit failed: {commit.stderr.strip() or commit.stdout}")

        push = _git(["push", "origin", "HEAD"], cwd=clone, env=env, token=token, timeout=timeout)
        if push.returncode != 0:
            # Someone edited the project in the browser since the clone: rebase
            # once onto their work and try again (integration.md section 7.1).
            pull = _git([*base, "pull", "--rebase", "origin"], cwd=clone, env=env, token=token,
                        timeout=timeout)
            if pull.returncode != 0:
                raise OverleafError(
                    "The project moved on and could not be rebased automatically: "
                    f"{pull.stderr.strip() or pull.stdout.strip()}"
                )
            push = _git(["push", "origin", "HEAD"], cwd=clone, env=env, token=token,
                        timeout=timeout)
            if push.returncode != 0:
                raise OverleafError(f"git push failed: {push.stderr.strip() or push.stdout}")

        head = _git(["rev-parse", "--short", "HEAD"], cwd=clone, env=env, token=token,
                    timeout=timeout)
        return PushResult(pushed=pushed, skipped=[], commit=head.stdout.strip() or None,
                          detail=f"Pushed {len(pushed)} figure(s) to {FIGURES_SUBDIR}/ "
                                 f"in project {project_id}.")
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def _clone_hint(stderr: str) -> str:
    """Turn git's clone failure into something a user can act on."""
    text = (stderr or "").strip()
    lowered = text.lower()
    if "authentication failed" in lowered or "403" in lowered:
        return (
            "Overleaf refused the token. Check that it is a current Git authentication "
            "token and that your account has Git integration, which is an Overleaf "
            "premium feature. Use the zip download instead if it does not."
        )
    if "not found" in lowered or "404" in lowered:
        return ("Overleaf does not recognise that project id, or your account cannot "
                "reach it. Check the project URL.")
    if "could not resolve host" in lowered or "timed out" in lowered:
        return "The server could not reach git.overleaf.com."
    return f"git clone failed: {text}"
