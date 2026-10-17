"""The ``mellor`` command.

    mellor run DATA.mint01 [...] [--out DIR] [--name NAME] [--prompt TEXT]
               [--backend auto|transformers|hf-api|scripted] [--model ID]
               [--quick] [--spec spec.json | --system KEY] [--no-chat]
               [--overwrite] [--quiet]
        The whole workflow. With --spec or --system the model is not asked for the
        specification (it still writes the discussion unless --backend scripted);
        afterwards it drops into the chat unless --no-chat or stdin is not a terminal.
    mellor chat RUN_DIR [--backend ...] [--model ID] [--quiet]
        Questions about a finished run; `/quit` ends it. When stdin is not a terminal
        every line of it is a question and the answers are printed.
    mellor report RUN_DIR [--discussion FILE.md]
        Re-make the analysis, the figures and the PDF of a run without a model.
    mellor spec DATA.mint01 [...] --system KEY [-o spec.json] [--name NAME]
        Write a library specification for the data (to edit by hand; stdout by default).
    mellor search QUERY [--sources local,dissolve,zenodo] [--limit N] [--download URL --to DIR]
        Comparable data from the repository, the Dissolve examples and Zenodo.
    mellor info
        Backends available, GPU, default model, versions.

Exit status 0 on success, 1 on a handled error (message on stderr), 2 on usage errors.

Without a model (``--backend scripted``, or no backend available while ``--spec`` or
``--system`` is given) the run still happens: the discussion of the report is then
the factual summary computed by ``mellor.analysis``, labelled as such. The progress
of a run (the model's replies, tool calls and the EPSR log) is printed to stderr as
it happens unless ``--quiet``.
"""

from __future__ import annotations

import argparse
import contextlib
import datetime as _dt
import importlib.util
import json
import os
import platform
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import IO, Any

from ._version import __version__

__all__ = ["main", "build_parser", "CliError", "UsageError", "event_printer"]

BACKENDS = ("auto", "transformers", "hf-api", "scripted")
#: Where ``mellor run`` writes without ``--out``: ``work/<run name>``.
DEFAULT_WORK_ROOT = Path("work")
#: Files whose presence marks a folder as holding an earlier run.
RUN_MARKERS = ("manifest.json", "report.pdf", "epsr")
QUIT_COMMANDS = ("/quit", "/exit", "/q")
PROMPT = "you> "
ANSWER_PREFIX = "Mellor> "
#: Characters of a tool argument / reply shown in the progress stream.
PREVIEW_CHARS = 200


class CliError(Exception):
    """A handled error: its message goes to stderr and the command returns 1."""


class UsageError(CliError):
    """Arguments that argparse cannot check: the message goes to stderr, status 2."""


# --------------------------------------------------------------------------- parser
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mellor",
        description=(
            "Mellor: an LLM agent that sets up, runs, analyses and reports EPSR "
            "refinements of measured structure-factor data (.mint01 files)."
        ),
    )
    parser.add_argument("--version", action="version", version=f"mellor {__version__}")
    commands = parser.add_subparsers(dest="command", metavar="COMMAND")
    commands.required = True

    def backend_options(sub: argparse.ArgumentParser) -> None:
        sub.add_argument(
            "--backend",
            default=None,
            metavar="KIND",
            help="one of "
            + ", ".join(BACKENDS)
            + " (default auto, or MELLOR_BACKEND); transformers and hf-api "
            "accept :MODEL",
        )
        sub.add_argument(
            "--model", default=None, metavar="ID", help="Hugging Face model id"
        )
        sub.add_argument(
            "--quiet",
            action="store_true",
            help="do not print the progress (model replies, tool calls, EPSR log)",
        )

    run = commands.add_parser(
        "run",
        help="run the whole workflow on .mint01 files",
        description="Set up, run, analyse and report an EPSR refinement of the data.",
    )
    run.add_argument(
        "data", nargs="*", metavar="DATA.mint01", help="the measured data files"
    )
    run.add_argument("--out", metavar="DIR", help="run folder (default work/<name>)")
    run.add_argument("--name", metavar="NAME", help="run name (spec.name)")
    run.add_argument(
        "--prompt", metavar="TEXT", help="what the sample is and what you want"
    )
    backend_options(run)
    run.add_argument(
        "--quick",
        action="store_true",
        help="the short schedule (10 iterations per run) for a preview",
    )
    how = run.add_mutually_exclusive_group()
    how.add_argument(
        "--spec", metavar="FILE", help="run this specification (JSON) as given"
    )
    how.add_argument(
        "--system", metavar="KEY", help="run the library system KEY (see mellor info)"
    )
    run.add_argument(
        "--no-chat", action="store_true", help="do not start the chat afterwards"
    )
    run.add_argument(
        "--overwrite",
        action="store_true",
        help="replace the results of an earlier run in the run folder",
    )
    run.set_defaults(func=run_command)

    chat = commands.add_parser(
        "chat",
        help="questions about a finished run",
        description="Chat with Mellor about a finished run (/quit to leave).",
    )
    chat.add_argument("run_dir", metavar="RUN_DIR", help="a finished run folder")
    backend_options(chat)
    chat.set_defaults(func=chat_command)

    report = commands.add_parser(
        "report",
        help="re-make the analysis, figures and PDF of a run",
        description="Re-make the analysis, the figures and the PDF report of a run "
        "without a model.",
    )
    report.add_argument("run_dir", metavar="RUN_DIR", help="a finished run folder")
    report.add_argument(
        "--discussion",
        metavar="FILE.md",
        help="use this Markdown text as the discussion section",
    )
    report.set_defaults(func=report_command)

    spec = commands.add_parser(
        "spec",
        help="write a library specification for the data",
        description="Write the specification of a library system for the data files "
        "(JSON, to edit by hand and run with mellor run --spec).",
    )
    spec.add_argument(
        "data", nargs="+", metavar="DATA.mint01", help="the measured data files"
    )
    spec.add_argument("--system", required=True, metavar="KEY", help="library key")
    spec.add_argument(
        "-o", "--output", metavar="FILE", help="write here (default: stdout)"
    )
    spec.add_argument("--name", metavar="NAME", help="run name (spec.name)")
    spec.set_defaults(func=spec_command)

    search = commands.add_parser(
        "search",
        help="look for comparable experimental data",
        description="Comparable data from this repository, the Dissolve examples "
        "(GitHub) and Zenodo; --download fetches one file.",
    )
    search.add_argument("query", nargs="?", metavar="QUERY", help="sample keywords")
    search.add_argument(
        "--sources",
        metavar="LIST",
        default=None,
        help="comma-separated: local, dissolve, zenodo (default: all)",
    )
    search.add_argument(
        "--limit", type=int, default=8, metavar="N", help="results per source"
    )
    search.add_argument("--download", metavar="URL", help="download this file")
    search.add_argument("--to", metavar="DIR", help="folder for --download")
    search.set_defaults(func=search_command)

    info = commands.add_parser(
        "info",
        help="backends available, GPU, default model, versions",
        description="What this installation can do.",
    )
    info.set_defaults(func=info_command)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exit_:  # --help (0) or a usage error (2)
        return exit_.code if isinstance(exit_.code, int) else 2
    try:
        return int(args.func(args))
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130
    except UsageError as error:
        print(f"usage error: {error}", file=sys.stderr)
        return 2
    except (CliError, LookupError, ValueError, OSError, RuntimeError) as error:
        # SpecError is a ValueError, RunError a RuntimeError; make_backend raises
        # RuntimeError when no model is available and ValueError for a bad kind
        message = str(error).strip() or type(error).__name__
        print(f"error: {message}", file=sys.stderr)
        return 1


# --------------------------------------------------------------------------- progress
def _short(text: str, limit: int = PREVIEW_CHARS) -> str:
    text = " ".join(str(text).split())
    return text if len(text) <= limit else text[: limit - 3] + "..."


def _short_arguments(arguments: Any) -> str:
    """Tool arguments on one line; a specification is shown by its name only."""
    if not isinstance(arguments, dict):
        return _short(json.dumps(arguments, default=str))
    shown = dict(arguments)
    spec = shown.get("spec")
    if isinstance(spec, dict):
        shown["spec"] = f"<spec {spec.get('name') or '?'}>"
    return _short(json.dumps(shown, ensure_ascii=False, default=str))


def event_printer(stream: IO[str] | None = None) -> Callable[[dict[str, Any]], None]:
    """An ``on_event`` callback that prints Mellor's events (default: stderr).

    Phases, the model's replies, tool calls with their arguments, tool results
    (``ok`` or ``error`` with a preview), the EPSR log lines (``progress``), notes,
    fallbacks and the report path.
    """

    def show(event: dict[str, Any]) -> None:
        out = stream if stream is not None else sys.stderr
        kind = event.get("type")
        if kind == "phase":
            line = f"== {event.get('phase')} =="
        elif kind == "assistant":
            content = _short(event.get("content") or "")
            calls = event.get("tool_calls") or []
            if not content and not calls:
                content = "(empty reply)"
            line = f"[model] {content}" if content else ""
            if calls and not content:
                line = "[model] " + ", ".join(
                    f"calls {call.get('name')}" for call in calls
                )
        elif kind == "tool_call":
            line = f"[tool] {event.get('name')}({_short_arguments(event.get('arguments'))})"
        elif kind == "tool_result":
            status = "ok" if event.get("ok") else "error"
            line = f"[tool] {event.get('name')} -> {status}: {_short(event.get('preview') or '')}"
        elif kind == "progress":
            line = f"  {event.get('line', '')}"
        elif kind == "note":
            line = f"[note] {event.get('content')}"
        elif kind == "fallback":
            line = f"[fallback] library system '{event.get('system')}': {event.get('reason')}"
        elif kind == "report":
            line = f"[report] {event.get('path')}"
        else:
            line = f"[{kind}] {_short(json.dumps(event, default=str))}"
        if line:
            print(line, file=out)
            out.flush()

    return show


# --------------------------------------------------------------------------- helpers
def _check_data_files(paths: Sequence[Path]) -> list[Path]:
    files = [Path(path).expanduser() for path in paths]
    missing = [str(path) for path in files if not path.is_file()]
    if missing:
        raise CliError("data file(s) not found: " + ", ".join(missing))
    return [path.resolve() for path in files]


def _read_datasets(files: Sequence[Path]) -> list[Any]:
    from . import data as data_module

    return [data_module.read_mint01(path) for path in files]


def _library_spec(key: str, files: Sequence[Path], name: str | None) -> Any:
    """The library specification of ``key`` for ``files`` (``LookupError``/``ValueError``)."""
    from . import library as library_module

    key = key.strip().lower().replace("-", "_").replace(" ", "_")
    return library_module.make_spec(key, _read_datasets(files), name=name)


def _load_spec_file(path: str, name: str | None) -> dict[str, Any]:
    """A specification file as a dict with ``source`` ``user``; relative data-file
    paths that do not exist from the current folder are taken relative to the file."""
    spec_path = Path(path).expanduser()
    if not spec_path.is_file():
        raise CliError(f"specification file not found: {spec_path}")
    try:
        data = json.loads(spec_path.read_text(encoding="utf-8"))
    except ValueError as error:
        raise CliError(f"{spec_path} is not valid JSON: {error}") from error
    if not isinstance(data, dict):
        raise CliError(f"{spec_path} does not hold a JSON object")
    data["source"] = "user"
    if name:
        data["name"] = name
    for dataset in data.get("datasets") or []:
        if not isinstance(dataset, dict):
            continue
        file = str(dataset.get("file") or "")
        candidate = Path(file).expanduser()
        if file and not candidate.is_absolute() and not candidate.is_file():
            beside = spec_path.parent / candidate
            if beside.is_file():
                dataset["file"] = str(beside.resolve())
    return data


def _spec_name(spec: Any) -> str:
    if isinstance(spec, dict):
        return str(spec.get("name") or "")
    return str(getattr(spec, "name", "") or "")


def _work_dir(args: argparse.Namespace, spec: Any) -> Path:
    """The run folder: ``--out``, else ``work/<run name>`` (absolute)."""
    if args.out:
        folder = Path(args.out).expanduser()
    else:
        name = args.name or _spec_name(spec) or args.system
        if not name:
            name = "run-" + _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
        folder = DEFAULT_WORK_ROOT / name
    return Path(os.path.abspath(folder))


def _check_work_dir(folder: Path, overwrite: bool) -> None:
    if folder.exists() and not folder.is_dir():
        raise CliError(f"{folder} exists and is not a folder")
    if overwrite or not folder.is_dir():
        return
    present = [name for name in RUN_MARKERS if (folder / name).exists()]
    if present:
        raise CliError(
            f"{folder} already holds a run ({', '.join(present)}); choose another "
            "--out or pass --overwrite to replace its results"
        )


def _make_backend(args: argparse.Namespace) -> Any:
    from .llm import backend as backend_module

    options: dict[str, Any] = {}
    if args.model:
        options["model"] = args.model
    return backend_module.make_backend(args.backend, **options)


def _is_terminal(stream: Any) -> bool:
    try:
        return bool(stream.isatty())
    except (AttributeError, ValueError, OSError):
        return False


def _chat_loop(agent: Any, out: IO[str], err: IO[str]) -> None:
    """Read questions from the terminal until ``/quit`` (or EOF, Ctrl-C)."""
    with contextlib.suppress(ImportError):  # Windows without pyreadline
        import readline  # noqa: F401  (line editing and history for input())
    print("Ask Mellor about the run (/quit to leave).", file=out)
    while True:
        try:
            line = input(PROMPT)
        except (EOFError, KeyboardInterrupt):
            print(file=out)
            return
        text = line.strip()
        if not text:
            continue
        if text.lower() in QUIT_COMMANDS:
            return
        try:
            answer = agent.chat(text)
        except Exception as error:  # one failed turn must not end the session
            print(f"error: {type(error).__name__}: {error}", file=err)
            continue
        print(f"{ANSWER_PREFIX}{answer}\n", file=out)


def _chat_from_stream(agent: Any, questions: IO[str], out: IO[str]) -> None:
    """Non-interactive chat: every non-empty line of ``questions`` is a turn."""
    for line in questions:
        text = line.strip()
        if not text:
            continue
        if text.lower() in QUIT_COMMANDS:
            return
        print(f"{PROMPT}{text}", file=out)
        print(f"{ANSWER_PREFIX}{agent.chat(text)}\n", file=out)


# --------------------------------------------------------------------------- run
def run_command(args: argparse.Namespace) -> int:
    from .llm.agent import Mellor
    from .llm.backend import ScriptedBackend
    from .llm.tools import Session

    out, err = sys.stdout, sys.stderr
    spec: Any = None
    data_files = [Path(path) for path in args.data]
    if args.spec:
        spec = _load_spec_file(args.spec, args.name)
        if not data_files:
            data_files = [
                Path(str(item.get("file")))
                for item in spec.get("datasets") or []
                if isinstance(item, dict) and item.get("file")
            ]
    if not data_files:
        raise UsageError("no data files given (.mint01)")
    files = _check_data_files(data_files)
    if args.system:
        spec = _library_spec(args.system, files, args.name)

    # the backend: with a given specification a missing model is not fatal
    try:
        backend = _make_backend(args)
    except RuntimeError as error:
        if spec is None:
            raise
        print(
            f"note: no language model ({error}); the discussion will be the factual summary",
            file=err,
        )
        backend = ScriptedBackend([])
    use_model = not isinstance(backend, ScriptedBackend)
    prompt = args.prompt
    if spec is None and not use_model:
        from . import library as library_module

        try:
            spec, reason = library_module.auto_spec(
                _read_datasets(files), name=args.name
            )
        except LookupError as error:
            raise CliError(
                f"{error}; without a language model name the system with "
                "--system KEY or give a specification with --spec FILE"
            ) from error
        print(f"note: no language model; {reason}", file=err)
        if prompt:
            print("note: the prompt is not used without a language model", file=err)
    elif spec is None and args.name:
        prompt = (
            (prompt or "").strip() + f"\nUse '{args.name}' as the run name (spec.name)."
        ).strip()

    work_dir = _work_dir(args, spec)
    _check_work_dir(work_dir, args.overwrite)
    session = Session(work_dir=work_dir, data_files=files, quick=args.quick)
    agent = Mellor(
        backend, session, on_event=None if args.quiet else event_printer(err)
    )
    outcome = agent.run_pipeline(prompt, spec=spec, narrative=use_model)

    print(file=out)
    if session.result is not None:
        print(session.result.table(), file=out)
        print(file=out)
    message = outcome.final_message.replace(f"Report: {outcome.report_path}", "")
    print(" ".join(message.split()), file=out)
    print(f"Run folder: {session.work_dir}", file=out)
    print(f"Report: {outcome.report_path}", file=out)
    if use_model and not args.no_chat and _is_terminal(sys.stdin):
        print(file=out)
        _chat_loop(agent, out, err)
    return 0


# --------------------------------------------------------------------------- chat
def chat_command(args: argparse.Namespace) -> int:
    from .llm.agent import Mellor
    from .llm.tools import Session

    out, err = sys.stdout, sys.stderr
    session = Session(work_dir=Path(args.run_dir).expanduser())
    session.load_existing()
    backend = _make_backend(args)
    agent = Mellor(
        backend, session, on_event=None if args.quiet else event_printer(err)
    )
    if _is_terminal(sys.stdin):
        _chat_loop(agent, out, err)
    else:
        _chat_from_stream(agent, sys.stdin, out)
    return 0


# --------------------------------------------------------------------------- report
def _narrative_provenance(folder: Path) -> tuple[str, str]:
    """``(source, model name)`` of the discussion from the transcript (``narrative``
    entries; the pipeline start entry for older transcripts); empty when unknown."""
    from .llm.agent import TRANSCRIPT_NAME

    path = folder / TRANSCRIPT_NAME
    if not path.is_file():
        return "", ""
    source, model = "", ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            entry = json.loads(line)
        except ValueError:
            continue
        if not isinstance(entry, dict):
            continue
        if entry.get("role") == "narrative":
            source = str(entry.get("source") or "")
            model = str(entry.get("model") or "")
        elif entry.get("role") == "pipeline" and not model:
            model = str(entry.get("model") or "")
    return source, model


def _discussion_for_report(
    args: argparse.Namespace, folder: Path, summary: dict[str, Any]
) -> tuple[str, str]:
    """``(discussion markdown, model_name for the report)``."""
    from . import analysis as analysis_module
    from .llm.agent import DISCUSSION_NAME, NO_MODEL, NO_MODEL_NOTE

    if args.discussion:
        path = Path(args.discussion).expanduser()
        if not path.is_file():
            raise CliError(f"discussion file not found: {path}")
        text = path.read_text(encoding="utf-8")
        (folder / DISCUSSION_NAME).write_text(text, encoding="utf-8")
        return text, "none: the discussion was given with --discussion"
    existing = folder / DISCUSSION_NAME
    if existing.is_file():
        text = existing.read_text(encoding="utf-8")
        source, model = _narrative_provenance(folder)
        if source == "summary":
            return text, NO_MODEL
        return text, model
    text = f"_{NO_MODEL_NOTE}_\n\n{analysis_module.analysis_markdown(summary)}"
    (folder / DISCUSSION_NAME).write_text(text + "\n", encoding="utf-8")
    return text, NO_MODEL


def report_command(args: argparse.Namespace) -> int:
    from . import analysis as analysis_module
    from . import figures as figures_module
    from . import report as report_module
    from . import runner as runner_module

    out = sys.stdout
    folder = Path(args.run_dir).expanduser()
    result = runner_module.load_result(folder)
    folder = result.folder
    summary = analysis_module.analyse(result)
    figures = figures_module.make_figures(result)
    summary["figures"] = {name: str(path) for name, path in figures.items()}
    analysis_module.save_analysis(summary, folder / "analysis.json")
    discussion, model_name = _discussion_for_report(args, folder, summary)
    path = report_module.build_pdf(
        result,
        summary,
        figures,
        discussion,
        folder / "report.pdf",
        model_name=model_name,
    )
    print(result.table(), file=out)
    print(file=out)
    print(f"Analysis: {folder / 'analysis.json'}", file=out)
    print(f"Figures: {len(figures)} files in {folder / 'figures'}", file=out)
    print(f"Report: {path}", file=out)
    return 0


# --------------------------------------------------------------------------- spec
def spec_command(args: argparse.Namespace) -> int:
    files = _check_data_files([Path(path) for path in args.data])
    spec = _library_spec(args.system, files, args.name)
    warnings = spec.validate(check_files=True)
    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    text = spec.to_json() + "\n"
    if args.output:
        target = Path(args.output).expanduser()
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")
        print(
            f"wrote {target} ({spec.source}, {len(spec.datasets)} data set(s)); "
            f"run it with: mellor run --spec {target}",
            file=sys.stderr,
        )
    else:
        sys.stdout.write(text)
    return 0


# --------------------------------------------------------------------------- search
def _print_candidates(found: dict[str, Any], out: IO[str]) -> None:
    candidates = found.get("candidates") or []
    print(f"{len(candidates)} candidate(s) for {found.get('query')!r}", file=out)
    for number, candidate in enumerate(candidates, start=1):
        files = candidate.get("files") or []
        shown = ", ".join(str(name) for name in files[:6])
        if len(files) > 6:
            shown += f", ... ({len(files)} files)"
        print(
            f"{number:3d}. [{candidate.get('source')}] score {float(candidate.get('score') or 0):.1f}  "
            f"{_short(candidate.get('title') or '', 100)}",
            file=out,
        )
        description = _short(candidate.get("description") or "", 160)
        if description:
            print(f"      {description}", file=out)
        if shown:
            print(f"      files: {shown}", file=out)
        if candidate.get("url"):
            print(f"      url:   {candidate['url']}", file=out)
        extra = candidate.get("extra") or {}
        downloads = extra.get("downloads") if isinstance(extra, dict) else None
        if isinstance(downloads, dict) and downloads:
            first_name, first_url = next(iter(downloads.items()))
            print(f"      download: {first_url}  ({first_name})", file=out)
    for note in found.get("notes") or []:
        print(f"note: {note}", file=out)


def search_command(args: argparse.Namespace) -> int:
    from . import sources as sources_module

    out = sys.stdout
    if not args.query and not args.download:
        raise UsageError("give a QUERY, or --download URL --to DIR")
    if args.download and not args.to:
        raise UsageError("--download needs --to DIR (the folder to save into)")
    if args.query:
        found = sources_module.find_similar_data(
            args.query, sources=args.sources, limit=max(1, int(args.limit))
        )
        _print_candidates(found, out)
    if args.download:
        path = sources_module.download(args.download, Path(args.to).expanduser())
        print(f"downloaded {path} ({path.stat().st_size} bytes)", file=out)
    return 0


# --------------------------------------------------------------------------- info
def _distribution_version(name: str) -> str:
    from importlib import metadata

    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return "not installed"


def _installed(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def info_command(args: argparse.Namespace) -> int:
    from .llm.backend import DEFAULT_MODEL

    out = sys.stdout
    print(
        f"mellor {__version__}  (python {platform.python_version()}, {sys.executable})",
        file=out,
    )
    versions = ", ".join(
        f"{name} {_distribution_version(name)}"
        for name in (
            "epsr-py",
            "epsrvis",
            "numpy",
            "scipy",
            "matplotlib",
            "reportlab",
            "huggingface_hub",
        )
    )
    print(f"packages: {versions}", file=out)
    print("backends:", file=out)
    local = _installed("torch") and _installed("transformers")
    if local:
        line = (
            f"  transformers: available (torch {_distribution_version('torch')}, "
            f"transformers {_distribution_version('transformers')}; {_cuda_description()})"
        )
    else:
        line = "  transformers: not installed (uv sync --extra local)"
    print(line, file=out)
    token = os.environ.get("HF_TOKEN", "")
    print(
        (
            "  hf-api: HF_TOKEN set (Hugging Face Inference API)"
            if token
            else "  hf-api: HF_TOKEN not set"
        ),
        file=out,
    )
    print("  scripted: no model (tests; mellor run --system KEY)", file=out)
    configured = os.environ.get("MELLOR_BACKEND") or "auto"
    resolved = configured
    if configured == "auto":
        resolved = (
            "transformers"
            if local
            else (
                "hf-api" if token else "none (use --backend scripted with --system KEY)"
            )
        )
    print(
        f"default backend: {configured} -> {resolved}"
        + ("" if os.environ.get("MELLOR_BACKEND") else "  (MELLOR_BACKEND unset)"),
        file=out,
    )
    print(
        f"default model: {os.environ.get('MELLOR_MODEL') or DEFAULT_MODEL}"
        + ("" if os.environ.get("MELLOR_MODEL") else "  (MELLOR_MODEL unset)"),
        file=out,
    )
    try:
        from . import library as library_module

        print("library systems: " + ", ".join(library_module.KNOWN_SYSTEMS), file=out)
    except Exception as error:  # the info must come out whatever is broken
        print(f"library systems: unavailable ({error})", file=out)
    try:
        from . import sources as sources_module

        roots = sources_module.default_local_roots()
        print(
            "local data folders: " + (", ".join(str(root) for root in roots) or "none"),
            file=out,
        )
    except Exception as error:
        print(f"local data folders: unavailable ({error})", file=out)
    return 0


def _cuda_description() -> str:
    try:
        import torch
    except Exception as error:  # a broken torch install must not break `info`
        return f"torch cannot be imported: {type(error).__name__}: {error}"
    try:
        if not torch.cuda.is_available():
            return "CUDA not available, the model would run on the CPU"
        properties = torch.cuda.get_device_properties(0)
        return f"CUDA: {properties.name}, {properties.total_memory / 1e9:.1f} GB"
    except Exception as error:
        return f"CUDA state unknown: {type(error).__name__}: {error}"


if __name__ == "__main__":
    sys.exit(main())
