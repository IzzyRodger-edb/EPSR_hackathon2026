"""The tools Mellor can call, and the session state they act on.

Every tool is a plain function taking JSON-compatible keyword arguments and returning
a JSON-serialisable value (dicts, lists, numbers, strings). The registry turns them
into OpenAI-style function schemas for the model and executes calls, converting
exceptions into ``{"error": ..., "hint": ...}`` results so the model can recover.

Tools (name -> what it does):

* ``inspect_data(files=None)``: ``mellor.data.describe`` of the session's data files
  (or the given ones): title, instrument, Q range, sample hints.
* ``list_known_systems()``: ``mellor.library.describe_library``.
* ``library_spec(system, name=None, adjust=None)``: the specification of a known system
  for the session's data (``library.make_spec``); returns the spec JSON and the
  warnings, and stores it as the session's draft.
* ``validate_spec(spec)``: ``SystemSpec.from_dict`` + ``validate(check_files=True)``;
  returns ``{"valid": bool, "errors": [...], "warnings": [...], "summary": {...}}``.
* ``build_and_run(spec, quick=None)``: validate, then ``runner.run_spec`` into the
  session's work folder (``quick`` defaults to the session setting); on success runs
  ``analysis.analyse`` and ``figures.make_figures`` and returns the analysis summary
  (compact: runs table, fit, convergence, shells, flags); on failure the error and the
  tail of the log.
* ``get_results(section=None)``: the stored analysis dict, or one section of it.
* ``get_partial(pair, kind="gr_inter", r_max=8.0)``: the final partial g(r) or S(Q)
  of one pair from the trajectory, down-sampled to ~60 points, with the first-shell
  numbers.
* ``list_files(folder=None)``: the run folder listing (figures, report, script...).
* ``read_text_file(path, max_chars=6000)``: an input/output text file of the run.
* ``write_report(discussion_markdown)``: ``report.build_pdf`` with the given narrative.
* ``find_similar_data(query, sources=None, limit=8)``: ``mellor.sources.find_similar_data``.
* ``download_data(url)``: ``mellor.sources.download`` into the session's data folder.
* ``finish(message)``: ends the tool loop with a final message to the user.

The schemas must stay small: a 4B model handles about a dozen tools well.
"""

from __future__ import annotations

import contextlib
import datetime as _dt
import inspect
import json
import math
import time
import traceback
from collections.abc import Callable
from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any

from ..runner import RunError, RunResult
from ..spec import SPEC_JSON_SCHEMA, SpecError, SystemSpec

__all__ = [
    "Tool",
    "ToolRegistry",
    "Session",
    "default_tools",
    "ToolError",
    "SessionTools",
    "jsonable",
]

#: Name of the trajectory folder when the manifest does not say.
DEFAULT_TRAJECTORY_FOLDER = "trajectory"
#: Lines of the run log returned after a failed run.
LOG_TAIL_LINES = 40
#: Extensions of files ``read_text_file`` refuses (binary).
BINARY_SUFFIXES = {".png", ".pdf", ".gif", ".f32", ".npy", ".npz", ".pyc", ".so"}
#: The analysis sections ``get_results`` serves besides the analysis dict itself.
EXTRA_SECTIONS = ("spec", "manifest", "table", "figures")


class ToolError(Exception):
    """A tool could not do what was asked; the message goes back to the model."""

    def __init__(self, message: str, hint: str = ""):
        super().__init__(message)
        self.hint = hint


# --------------------------------------------------------------------------- JSON helpers
def jsonable(value: Any) -> Any:
    """``value`` with numpy scalars/arrays, paths, tuples, sets and dataclasses turned
    into plain JSON types; non-finite floats become ``None``."""
    if value is None or isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value) and not isinstance(value, type):
        return {f.name: jsonable(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [jsonable(item) for item in value]
    if hasattr(value, "tolist"):  # numpy arrays and scalars
        return jsonable(value.tolist())
    if hasattr(value, "item"):
        return jsonable(value.item())
    if isinstance(value, (_dt.datetime, _dt.date)):
        return value.isoformat()
    return str(value)


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def _short(value: Any, limit: int = 300) -> str:
    text = json.dumps(jsonable(value), ensure_ascii=False, default=str)
    return text if len(text) <= limit else text[: limit - 3] + "..."


# --------------------------------------------------------------------------- session
@dataclass
class Session:
    """What the agent knows and where it writes."""

    work_dir: Path
    data_files: list[Path] = field(default_factory=list)
    quick: bool = False
    draft_spec: SystemSpec | None = None
    result: RunResult | None = None
    analysis: dict[str, Any] | None = None
    figures: dict[str, Path] = field(default_factory=dict)
    report_path: Path | None = None
    model_name: str = ""
    spec_source: str = ""  # "model", "library:<key>", "user", "fallback:library:<key>"
    events: list[dict[str, Any]] = field(
        default_factory=list
    )  # tool calls and outcomes
    data_descriptions: list[dict[str, Any]] | None = None  # data.describe_all, cached
    progress: Callable[[str], None] | None = None  # receives the run log line by line

    def __post_init__(self) -> None:
        self.work_dir = Path(self.work_dir)
        self.data_files = [Path(path) for path in self.data_files]

    @property
    def data_folder(self) -> Path:
        return self.work_dir / "data"

    def load_existing(self) -> None:
        """Restore ``result``, ``analysis``, ``figures`` and ``report_path`` from ``work_dir``
        (for ``mellor chat`` on a finished run)."""
        from .. import analysis as analysis_module
        from .. import runner as runner_module

        folder = self.work_dir
        if not folder.is_dir():
            raise FileNotFoundError(f"run folder {folder} does not exist")
        manifest = folder / "manifest.json"
        analysis_path = folder / "analysis.json"
        report = folder / "report.pdf"
        if not (manifest.is_file() or analysis_path.is_file() or report.is_file()):
            raise FileNotFoundError(
                f"{folder} holds no finished run (no manifest.json, analysis.json or report.pdf)"
            )
        if manifest.is_file() and (folder / "spec.json").is_file():
            self.result = runner_module.load_result(folder)
            self.draft_spec = self.result.spec
            self.quick = bool(self.result.quick)
            if not self.spec_source:
                self.spec_source = self.result.spec.source or "user"
        if analysis_path.is_file():
            self.analysis = analysis_module.load_analysis(analysis_path)
        figures_dir = folder / "figures"
        if figures_dir.is_dir():
            for path in sorted(figures_dir.iterdir()):
                if path.suffix.lower() in (".png", ".csv", ".md", ".gif"):
                    key = (
                        path.stem
                        if path.suffix.lower() in (".png", ".gif")
                        else path.name
                    )
                    self.figures[key] = path
        if report.is_file():
            self.report_path = report
        if not self.data_files and self.data_folder.is_dir():
            self.data_files = sorted(self.data_folder.glob("*.mint01"))
        self.events.append(
            {"type": "load_existing", "folder": str(folder), "time": _now()}
        )


# --------------------------------------------------------------------------- registry
@dataclass
class Tool:
    name: str
    description: str
    parameters: dict[str, Any]  # JSON schema of the arguments
    func: Callable[..., Any]

    def schema(self) -> dict[str, Any]:
        return {
            "type": "function",
            "function": {
                "name": self.name,
                "description": self.description,
                "parameters": self.parameters,
            },
        }


class ToolRegistry:
    def __init__(self) -> None:
        self._tools: dict[str, Tool] = {}
        self.events: list[dict[str, Any]] | None = None  # where calls are recorded

    def register(self, tool: Tool) -> None:
        self._tools[tool.name] = tool

    def names(self) -> list[str]:
        return list(self._tools)

    def get(self, name: str) -> Tool | None:
        return self._tools.get(name)

    def __contains__(self, name: object) -> bool:
        return name in self._tools

    def __len__(self) -> int:
        return len(self._tools)

    def schemas(self, names: list[str] | None = None) -> list[dict[str, Any]]:
        if names is None:
            return [tool.schema() for tool in self._tools.values()]
        unknown = [name for name in names if name not in self._tools]
        if unknown:
            raise KeyError(f"unknown tool(s) {unknown}; registered: {self.names()}")
        wanted = set(names)
        return [tool.schema() for tool in self._tools.values() if tool.name in wanted]

    def call(self, name: str, arguments: dict[str, Any]) -> Any:
        """Run a tool; unknown names, bad arguments and exceptions become ``{"error": ...}``."""
        started = time.monotonic()
        result = self._call(name, arguments)
        if self.events is not None:
            error = result.get("error") if isinstance(result, dict) else None
            self.events.append(
                {
                    "type": "tool_call",
                    "tool": name,
                    "arguments": _short(arguments),
                    "ok": error is None,
                    "error": error,
                    "elapsed_s": round(time.monotonic() - started, 3),
                    "time": _now(),
                }
            )
        return result

    def _call(self, name: str, arguments: Any) -> Any:
        tool = self._tools.get(name)
        if tool is None:
            return {
                "error": f"unknown tool {name!r}",
                "hint": "available tools: " + ", ".join(self.names()),
            }
        if arguments is None:
            arguments = {}
        if isinstance(arguments, str):
            try:
                arguments = json.loads(arguments)
            except json.JSONDecodeError as error:
                return {
                    "error": f"arguments of {name} are not valid JSON: {error}",
                    "hint": "pass the arguments as a JSON object",
                }
        if not isinstance(arguments, dict):
            return {
                "error": f"arguments of {name} must be a JSON object, not {type(arguments).__name__}",
                "hint": "pass the arguments as a JSON object",
            }
        try:
            inspect.signature(tool.func).bind(**arguments)
        except TypeError as error:
            expected = list((tool.parameters.get("properties") or {}).keys())
            return {
                "error": f"bad arguments for {name}: {error}",
                "hint": f"{name} takes: {', '.join(expected) or 'no arguments'}",
            }
        try:
            result = tool.func(**arguments)
        except ToolError as error:
            return {
                "error": str(error),
                "hint": error.hint or "check the arguments and try again",
            }
        except SpecError as error:
            return {
                "error": "invalid specification",
                "errors": list(error.errors),
                "hint": "fix every listed problem and call validate_spec again",
            }
        except RunError as error:
            return {
                "error": str(error),
                "log_tail": _log_tail(error.log),
                "hint": "read the log tail, fix the specification and run again",
            }
        except Exception as error:  # any failure must reach the model
            return {
                "error": f"{type(error).__name__}: {error}{_where(error)}",
                "hint": "this tool failed; check the arguments or use another tool",
            }
        try:
            return jsonable(result)
        except Exception as error:
            return {
                "error": f"the result of {name} could not be serialised: {error}",
                "hint": "",
            }

    @staticmethod
    def serialise(result: Any, *, max_chars: int = 12_000) -> str:
        """``result`` as JSON text, truncated (with a note) to ``max_chars``."""
        text = json.dumps(jsonable(result), ensure_ascii=False, default=str)
        if len(text) <= max_chars:
            return text
        dropped = len(text) - max_chars
        return (
            text[:max_chars]
            + f"\n... [truncated: {dropped} more characters; ask for one section or file at a time]"
        )

    def call_json(
        self, name: str, arguments: dict[str, Any], *, max_chars: int = 12_000
    ) -> str:
        """``call`` serialised to JSON, truncated (with a note) to ``max_chars``."""
        return self.serialise(self.call(name, arguments), max_chars=max_chars)


def _log_tail(log: str, lines: int = LOG_TAIL_LINES) -> list[str]:
    return [line.rstrip() for line in (log or "").splitlines()[-lines:]]


def _where(error: BaseException) -> str:
    """' (in <function>, <file>:<line>)' of the innermost frame of an exception."""
    frames = traceback.extract_tb(error.__traceback__) if error.__traceback__ else []
    if not frames:
        return ""
    frame = frames[-1]
    return f" (in {frame.name}, {Path(frame.filename).name}:{frame.lineno})"


# --------------------------------------------------------------------------- the tools
class SessionTools:
    """The tool functions, bound to one :class:`Session`.

    Each public method is one tool (same names and arguments as the module
    docstring). The methods raise :class:`ToolError` for problems the model can fix;
    :func:`default_tools` wraps them in a :class:`ToolRegistry`.
    """

    def __init__(self, session: Session):
        self.session = session

    # ---- data
    def inspect_data(self, files: list[str] | None = None) -> dict[str, Any]:
        from .. import data as data_module

        paths = (
            [Path(item) for item in files] if files else list(self.session.data_files)
        )
        if not paths:
            raise ToolError(
                "the session has no data files", "pass 'files' with .mint01 paths"
            )
        missing = [str(path) for path in paths if not path.is_file()]
        if missing:
            raise ToolError(
                f"data file(s) not found: {missing}", "use the paths listed in the task"
            )
        files_arg: list[str | Path] = list(paths)
        descriptions = [dict(item) for item in data_module.describe_all(files_arg)]
        for description, path in zip(descriptions, paths, strict=False):
            description.setdefault("file", path.name)
            description.setdefault("path", str(path))
        if not files:
            self.session.data_descriptions = descriptions
        return {"n_files": len(descriptions), "files": descriptions}

    def list_known_systems(self) -> dict[str, Any]:
        from .. import library as library_module

        return {"systems": library_module.describe_library()}

    def library_spec(
        self, system: str, name: str | None = None, adjust: dict[str, Any] | None = None
    ) -> dict[str, Any]:
        from .. import library as library_module

        key = str(system).strip().lower().replace("-", "_").replace(" ", "_")
        if not self.session.data_files:
            raise ToolError(
                "the session has no data files to build a specification for"
            )
        if adjust is not None and not isinstance(adjust, dict):
            raise ToolError("'adjust' must be an object, e.g. {\"ereq\": 5}")
        datasets = self._read_datasets()
        try:
            spec = library_module.make_spec(key, datasets, name=name, **(adjust or {}))
        except (KeyError, LookupError) as error:
            known = list(getattr(library_module, "KNOWN_SYSTEMS", {}) or [])
            raise ToolError(
                f"unknown system {system!r}: {error}",
                f"known systems: {', '.join(known) or 'see list_known_systems'}",
            ) from error
        except TypeError as error:
            raise ToolError(
                f"library_spec({system!r}) rejected the adjustments: {error}"
            ) from error
        if not spec.source:
            spec.source = f"library:{key}"
        warnings = spec.validate(check_files=True)
        self.session.draft_spec = spec
        return {
            "system": key,
            "spec": spec.to_dict(),
            "warnings": warnings,
            "summary": spec_summary(spec),
            "next": "pass this spec (adapted if needed) to build_and_run",
        }

    def _read_datasets(self) -> list[Any]:
        from .. import data as data_module

        datasets = []
        for path in self.session.data_files:
            try:
                datasets.append(data_module.read_mint01(path))
            except (OSError, ValueError) as error:
                raise ToolError(f"cannot read {path}: {error}") from error
        return datasets

    # ---- specification
    def _spec_from_argument(self, spec: Any) -> SystemSpec:
        if isinstance(spec, str):
            try:
                spec = json.loads(spec)
            except json.JSONDecodeError as error:
                raise SpecError([f"the spec is not valid JSON: {error}"]) from error
        if (
            isinstance(spec, dict)
            and isinstance(spec.get("spec"), dict)
            and "schedule" not in spec
        ):
            spec = spec["spec"]  # the model echoed a library_spec result
        if not isinstance(spec, dict):
            raise SpecError(["the spec must be a JSON object"])
        parsed = SystemSpec.from_dict(spec)
        resolve_data_files(parsed, self.session)
        return parsed

    def validate_spec(self, spec: dict[str, Any]) -> dict[str, Any]:
        try:
            parsed = self._spec_from_argument(spec)
            warnings = parsed.validate(check_files=True)
        except SpecError as error:
            return {"valid": False, "errors": list(error.errors), "warnings": []}
        self.session.draft_spec = parsed
        return {
            "valid": True,
            "errors": [],
            "warnings": warnings,
            "summary": spec_summary(parsed),
        }

    def build_and_run(
        self, spec: dict[str, Any], quick: bool | None = None
    ) -> dict[str, Any]:
        from .. import analysis as analysis_module
        from .. import figures as figures_module
        from .. import runner as runner_module

        session = self.session
        try:
            parsed = self._spec_from_argument(spec)
            warnings = parsed.validate(check_files=True)
        except SpecError as error:
            return {
                "ok": False,
                "errors": list(error.errors),
                "hint": "fix every listed problem and call build_and_run again",
            }
        if not parsed.source:
            parsed.source = "model"
        session.draft_spec = parsed
        use_quick = session.quick if quick is None else bool(quick)
        folder = self._run_folder()
        try:
            result = runner_module.run_spec(
                parsed,
                folder,
                quick=use_quick,
                overwrite=True,
                progress=session.progress,
            )
        except RunError as error:
            session.events.append(
                {
                    "type": "run_failed",
                    "folder": str(folder),
                    "error": str(error),
                    "time": _now(),
                }
            )
            return {
                "ok": False,
                "error": str(error),
                "log_tail": _log_tail(error.log),
                "hint": "read the log tail; fix the specification (or the data files) and run again",
            }
        session.result = result
        session.spec_source = parsed.source
        session.analysis = None
        session.figures = {}
        session.report_path = None
        response: dict[str, Any] = {
            "ok": True,
            "folder": str(result.folder),
            "spec_source": parsed.source,
            "quick": use_quick,
            "warnings": warnings,
        }
        try:
            analysis = analysis_module.analyse(result)
            session.analysis = analysis
        except Exception as error:  # the run is done; report the problem
            response["analysis_error"] = f"{type(error).__name__}: {error}"
            analysis = None
        if analysis is not None:
            try:
                session.figures = dict(figures_module.make_figures(result))
            except Exception as error:
                response["figures_error"] = f"{type(error).__name__}: {error}"
            analysis["figures"] = {
                name: str(path) for name, path in session.figures.items()
            }
            try:
                analysis_module.save_analysis(analysis, result.folder / "analysis.json")
            except Exception as error:
                response["analysis_error"] = f"could not save analysis.json: {error}"
            response.update(compact_analysis(analysis))
            response["figures"] = sorted(session.figures)
        else:
            response.update(
                {
                    "runs": [jsonable(run) for run in result.runs],
                    "flags": [],
                    "figures": [],
                }
            )
        session.events.append(
            {
                "type": "run",
                "folder": str(result.folder),
                "spec_source": parsed.source,
                "quick": use_quick,
                "duration_s": result.duration_s,
                "time": _now(),
            }
        )
        return response

    def _run_folder(self) -> Path:
        """The session's work folder for the first run; ``rerun_<n>`` below it after a
        successful run, so the first results are kept."""
        session = self.session
        if session.result is None:
            return session.work_dir
        index = 1
        while (session.work_dir / f"rerun_{index}").exists():
            index += 1
        return session.work_dir / f"rerun_{index}"

    # ---- results
    def get_results(self, section: str | None = None) -> dict[str, Any]:
        session = self.session
        if session.analysis is None and session.result is None:
            raise ToolError("no results yet", "run build_and_run first")
        analysis = session.analysis or {}
        result = session.result
        if section is None:
            if analysis:
                return analysis
            if result is None:
                raise ToolError("no results yet", "run build_and_run first")
            return {"runs": [jsonable(run) for run in result.runs]}
        section = str(section).strip()
        if section in analysis:
            return {section: analysis[section]}
        if section == "spec" and session.result is not None:
            return {"spec": session.result.spec.to_dict()}
        if section == "manifest" and session.result is not None:
            return {"manifest": session.result.manifest}
        if section == "table" and session.result is not None:
            return {"table": session.result.table()}
        if section == "figures":
            return {
                "figures": {name: str(path) for name, path in session.figures.items()}
            }
        sections = list(analysis) + [s for s in EXTRA_SECTIONS if s not in analysis]
        raise ToolError(
            f"unknown section {section!r}", f"sections: {', '.join(sections)}"
        )

    def get_partial(
        self, pair: str, kind: str = "gr_inter", r_max: float = 8.0
    ) -> dict[str, Any]:
        from epsrvis import Trajectory, metrics

        session = self.session
        if session.result is None:
            raise ToolError("no run yet", "run build_and_run first")
        folder = trajectory_folder(session.result)
        if folder is None:
            raise ToolError(
                "the run has no recorded trajectory", "run again with record=True"
            )
        trajectory = Trajectory.load(folder)
        kinds = (
            "gr_inter",
            "gr_intra",
            "gr_total",
            "sq_inter",
            "sq_intra",
            "ep",
            "ref_pot",
        )
        kind = str(kind).strip().lower()
        if kind not in kinds:
            raise ToolError(f"unknown kind {kind!r}", f"kinds: {', '.join(kinds)}")
        index = _pair_index(trajectory.pairs, pair)
        step = trajectory.nsteps - 1
        if step < 0:
            raise ToolError("the trajectory is empty")
        if kind == "gr_total":
            values = trajectory.gr_total(step)[index]
        else:
            array = trajectory.arrays.get(kind)
            if array is None:
                raise ToolError(
                    f"{kind} was not recorded",
                    f"recorded: {', '.join(trajectory.arrays)}",
                )
            values = array[step][index]
        in_r = kind.startswith("gr") or kind in ("ep", "ref_pot")
        x = trajectory.r if in_r else trajectory.q
        values = [float(v) for v in values]
        if in_r:
            keep = [i for i, r in enumerate(x) if r <= float(r_max)]
        else:
            keep = list(range(len(x)))
        stride = max(1, len(keep) // 60)
        keep = keep[::stride]
        result: dict[str, Any] = {
            "pair": trajectory.pairs[index],
            "kind": kind,
            "iteration": int(trajectory.steps[step]),
            "x_label": "r (Å)" if in_r else "Q (Å^-1)",
            "y_label": {"ep": "U_EP (kJ/mol)", "ref_pot": "U_ref (kJ/mol)"}.get(
                kind, "g(r)" if in_r else "S(Q)"
            ),
            "x": [round(float(x[i]), 4) for i in keep],
            "y": [round(values[i], 5) for i in keep],
            "n_points_recorded": len(values),
        }
        if kind.startswith("gr"):
            coordination = metrics.running_coordination(
                trajectory, step, total=kind == "gr_total"
            )
            shell = metrics.first_shell(
                trajectory.r, _as_array(values), coordination[index]
            )
            result["first_shell"] = {
                "r_peak": shell["r_peak"],
                "g_peak": shell["g_peak"],
                "r_min": shell["r_min"],
                "coordination_number": shell["n_coord"],
                "note": "coordination number N(r_min): partners of the second type around the first"
                + (" (intermolecular only)" if kind == "gr_inter" else ""),
            }
        return jsonable(result)

    # ---- files
    def _allowed_roots(self) -> list[Path]:
        from .. import library as library_module

        roots = [self.session.work_dir.resolve()]
        roots += [Path(path).resolve().parent for path in self.session.data_files]
        with contextlib.suppress(Exception):  # only a convenience
            roots.append(library_module.repository_root().resolve())
        return roots

    def _resolve_path(self, path: str) -> Path:
        session = self.session
        candidate = Path(path).expanduser()
        candidates = (
            [candidate]
            if candidate.is_absolute()
            else [
                session.work_dir / candidate,
                session.work_dir / "epsr" / candidate,
                session.work_dir / "figures" / candidate,
                session.data_folder / candidate,
                Path.cwd() / candidate,
            ]
        )
        for item in candidates:
            if item.exists():
                resolved = item.resolve()
                if any(
                    resolved == root or root in resolved.parents
                    for root in self._allowed_roots()
                ):
                    return resolved
                raise ToolError(f"{path} is outside the run folder and the repository")
        raise ToolError(
            f"{path} not found", f"relative to {session.work_dir}; use list_files"
        )

    def list_files(self, folder: str | None = None) -> dict[str, Any]:
        root = self.session.work_dir if folder is None else self._resolve_path(folder)
        if not root.is_dir():
            raise ToolError(f"{root} is not a folder")
        entries: list[str] = []
        limit = 300
        for path in sorted(root.rglob("*")):
            if path.name.startswith(".") or "__pycache__" in path.parts:
                continue
            relative = path.relative_to(root)
            if len(relative.parts) > 2:
                continue
            if path.is_dir():
                entries.append(f"{relative}/")
            else:
                entries.append(f"{relative} ({_size(path.stat().st_size)})")
            if len(entries) >= limit:
                entries.append("... (listing truncated)")
                break
        return {"folder": str(root), "files": entries}

    def read_text_file(self, path: str, max_chars: int = 6000) -> dict[str, Any]:
        resolved = self._resolve_path(path)
        if resolved.is_dir():
            raise ToolError(f"{path} is a folder", "use list_files")
        if resolved.suffix.lower() in BINARY_SUFFIXES:
            raise ToolError(f"{path} is a binary file", "only text files can be read")
        max_chars = max(200, int(max_chars))
        text = resolved.read_text(errors="replace")
        truncated = len(text) > max_chars
        return {
            "path": str(resolved),
            "chars": len(text),
            "truncated": truncated,
            "text": text[:max_chars] + ("\n... [truncated]" if truncated else ""),
        }

    # ---- report
    def write_report(self, discussion_markdown: str) -> dict[str, Any]:
        from .. import report as report_module

        session = self.session
        if session.result is None:
            raise ToolError("no run to report on", "run build_and_run first")
        if not isinstance(discussion_markdown, str) or not discussion_markdown.strip():
            raise ToolError("discussion_markdown is empty")
        out_path = session.result.folder / "report.pdf"
        path = report_module.build_pdf(
            session.result,
            session.analysis or {},
            session.figures,
            discussion_markdown,
            out_path,
            model_name=session.model_name,
            data_descriptions=session.data_descriptions,
        )
        session.report_path = Path(path)
        session.events.append({"type": "report", "path": str(path), "time": _now()})
        return {"ok": True, "report": str(path)}

    # ---- comparable data
    def find_similar_data(
        self, query: str, sources: list[str] | str | None = None, limit: int = 8
    ) -> dict[str, Any]:
        from .. import sources as sources_module

        if not str(query).strip():
            raise ToolError(
                "query is empty",
                "describe the sample, e.g. 'D2O water neutron SANDALS'",
            )
        kwargs: dict[str, Any] = {"limit": max(1, int(limit))}
        if sources:
            names = sources.split(",") if isinstance(sources, str) else list(sources)
            kwargs["sources"] = tuple(
                str(name).strip().lower() for name in names if str(name).strip()
            )
        return sources_module.find_similar_data(str(query), **kwargs)

    def download_data(self, url: str) -> dict[str, Any]:
        from .. import sources as sources_module

        session = self.session
        session.data_folder.mkdir(parents=True, exist_ok=True)
        try:
            path = Path(sources_module.download(str(url), session.data_folder))
        except (OSError, ValueError, RuntimeError) as error:
            raise ToolError(f"download failed: {error}") from error
        added = path.suffix.lower() == ".mint01" and path not in session.data_files
        if added:
            session.data_files.append(path)
        return {
            "path": str(path),
            "bytes": path.stat().st_size if path.exists() else 0,
            "added_to_session": added,
            "data_files": [str(item) for item in session.data_files],
        }

    def finish(self, message: str = "") -> dict[str, Any]:
        return {"finished": True, "message": str(message)}


# --------------------------------------------------------------------------- helpers
def _size(n_bytes: int) -> str:
    if n_bytes < 1024:
        return f"{n_bytes} B"
    if n_bytes < 1024 * 1024:
        return f"{n_bytes / 1024:.0f} kB"
    return f"{n_bytes / 1024 / 1024:.1f} MB"


def _as_array(values: list[float]) -> Any:
    import numpy as np

    return np.asarray(values, dtype=float)


def _pair_index(pairs: list[str], pair: str) -> int:
    wanted = str(pair).strip().upper().replace("_", "-").replace(" ", "-")
    for index, label in enumerate(pairs):
        first, _, second = label.partition("-")
        if wanted in (label.upper(), f"{second}-{first}".upper()):
            return index
    raise ToolError(f"unknown pair {pair!r}", f"pairs: {', '.join(pairs)}")


def trajectory_folder(result: RunResult) -> Path | None:
    """The run's trajectory folder, or None when it was not recorded."""
    try:
        path = result.trajectory_path
    except NotImplementedError:
        path = None
    if path is None:
        path = result.folder / str(
            result.manifest.get("trajectory") or DEFAULT_TRAJECTORY_FOLDER
        )
    path = Path(path)
    return path if (path / "meta.json").is_file() else None


def resolve_data_files(spec: SystemSpec, session: Session) -> None:
    """Replace bare data-file names in ``spec`` by the session's paths (in place)."""
    known = {Path(path).name: Path(path) for path in session.data_files}
    for dataset in spec.datasets:
        if not dataset.file or Path(dataset.file).is_file():
            continue
        name = Path(dataset.file).name
        if name in known:
            dataset.file = str(known[name])
        elif (session.data_folder / name).is_file():
            dataset.file = str(session.data_folder / name)
    if spec.pdb is not None and spec.pdb.file and not Path(spec.pdb.file).is_file():
        for folder in (session.work_dir, session.data_folder):
            candidate = folder / Path(spec.pdb.file).name
            if candidate.is_file():
                spec.pdb.file = str(candidate)
                break


def spec_summary(spec: SystemSpec) -> dict[str, Any]:
    """A compact, JSON-serialisable summary of a specification for the model."""
    summary: dict[str, Any] = {
        "name": spec.name,
        "kind": spec.kind,
        "molecules": [
            {"name": m.name, "count": m.count, "atoms": m.n_atoms}
            for m in spec.molecules
        ],
        "atom_types": spec.labels,
        "n_atoms": spec.n_atoms,
        "n_molecules": spec.n_molecules,
        "density": spec.density,
        "temperature": spec.temperature,
        "datasets": [
            {
                "label": d.label,
                "file": d.file_name,
                "isotopes": {
                    k: [[i.mass_number, i.fraction] for i in v]
                    for k, v in d.isotopes.items()
                },
                "exchangeable": list(d.exchangeable),
            }
            for d in spec.datasets
        ],
        "schedule": {
            "runs": [list(run) for run in spec.schedule.iterations()],
            "total_iterations": spec.schedule.total_iterations(),
            "quick_total_iterations": spec.schedule.total_iterations(quick=True),
            "ereq": spec.schedule.ereq,
        },
        "epsr_keywords": dict(spec.epsr_keywords),
        "source": spec.source,
    }
    try:
        summary["box_side"] = round(spec.box_side, 3)
        summary["mass_density_g_cm3"] = round(spec.mass_density(), 4)
    except (ZeroDivisionError, KeyError, ValueError):
        pass
    return summary


def compact_analysis(analysis: dict[str, Any]) -> dict[str, Any]:
    """The parts of an analysis dict the model needs after a run."""
    runs = []
    for run in analysis.get("runs") or []:
        if not isinstance(run, dict):
            continue
        compact = {
            key: run[key]
            for key in (
                "run",
                "stage",
                "niter",
                "ereq",
                "r_factor",
                "r_factor_per_dataset",
                "energy",
                "pressure",
                "ep_amplitude",
                "duration_s",
            )
            if key in run
        }
        runs.append(compact)
    return jsonable(
        {
            "runs": runs,
            "fit": analysis.get("fit"),
            "convergence": analysis.get("convergence"),
            "shells": analysis.get("shells"),
            "flags": list(analysis.get("flags") or []),
        }
    )


# --------------------------------------------------------------------------- registry factory
_SPEC_PARAMETER = {
    "type": "object",
    "description": "The complete system specification (JSON object). Start from library_spec and change only what the data require.",
}


def default_tools(session: Session) -> ToolRegistry:
    """The registry of the tools listed in the module docstring, bound to ``session``."""
    tools = SessionTools(session)
    registry = ToolRegistry()
    registry.events = session.events

    def add(
        name: str,
        description: str,
        properties: dict[str, Any],
        required: list[str] | None = None,
    ) -> None:
        parameters: dict[str, Any] = {"type": "object", "properties": properties}
        if required:
            parameters["required"] = required
        registry.register(Tool(name, description, parameters, getattr(tools, name)))

    add(
        "inspect_data",
        "Describe the measured .mint01 data files: sample title, instrument, run number, Q range, number of points, high-Q level and sample hints (H/D labelling, candidate systems).",
        {
            "files": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Paths; default: the session's data files.",
            }
        },
    )
    add(
        "list_known_systems",
        "List the built-in library of validated EPSR systems (key, title, description, atom types, keywords).",
        {},
    )
    add(
        "library_spec",
        "Build the specification of a known system (by library key) for the session's data files, with the isotopic composition of every data set taken from its title. Returns the spec JSON, warnings and a summary.",
        {
            "system": {
                "type": "string",
                "description": "Library key, e.g. 'water', 'silica', 'dmso_water', 'benzene', 'argon', 'cu_btc'.",
            },
            "name": {"type": "string", "description": "Run name (optional)."},
            "adjust": {
                "type": "object",
                "description": 'Optional overrides, e.g. {"temperature": 298, "ereq": 3}.',
            },
        },
        ["system"],
    )
    add(
        "validate_spec",
        "Check a specification without running it: returns valid, errors (to fix), warnings and a summary (box, counts, schedule).",
        {"spec": _SPEC_PARAMETER},
        ["spec"],
    )
    add(
        "build_and_run",
        "Validate the specification, generate the EPSR inputs and the Python script, run the EPSR refinement, then analyse it and make the figures. Returns ok with the runs table, fit quality, convergence, first shells, flags and figure names; or ok=false with the errors.",
        {
            "spec": dict(SPEC_JSON_SCHEMA),
            "quick": {
                "type": "boolean",
                "description": "Short schedule for a preview (default: the session setting).",
            },
        },
        ["spec"],
    )
    add(
        "get_results",
        "The analysis of the finished run: all sections, or one of system, schedule, runs, fit, convergence, shells, flags, spec, manifest, table, figures.",
        {"section": {"type": "string", "description": "Section name (optional)."}},
    )
    add(
        "get_partial",
        "The final partial function of one atom pair from the recorded trajectory, down-sampled to about 60 points, with the first-shell numbers (peak position, height, first minimum, coordination number).",
        {
            "pair": {"type": "string", "description": "Pair label, e.g. 'OW-HW'."},
            "kind": {
                "type": "string",
                "enum": [
                    "gr_inter",
                    "gr_intra",
                    "gr_total",
                    "sq_inter",
                    "sq_intra",
                    "ep",
                    "ref_pot",
                ],
                "description": "Default gr_inter.",
            },
            "r_max": {
                "type": "number",
                "description": "Largest r (Å) returned for r-space functions (default 8).",
            },
        },
        ["pair"],
    )
    add(
        "list_files",
        "List the files of the run folder (inputs, outputs, figures, report, script).",
        {"folder": {"type": "string", "description": "Sub-folder or path (optional)."}},
    )
    add(
        "read_text_file",
        "Read a text file of the run (input files, .EPSR.inp, run.log, the Python script, figures/summary.md).",
        {
            "path": {
                "type": "string",
                "description": "File name or path, relative to the run folder.",
            },
            "max_chars": {"type": "integer", "description": "Default 6000."},
        },
        ["path"],
    )
    add(
        "write_report",
        "Build the PDF report of the finished run with the given discussion text (Markdown).",
        {
            "discussion_markdown": {
                "type": "string",
                "description": "The discussion section, Markdown.",
            }
        },
        ["discussion_markdown"],
    )
    add(
        "find_similar_data",
        "Search for comparable experimental data: .mint01 files in this repository, the Dissolve example data sets (GitHub) and Zenodo records.",
        {
            "query": {
                "type": "string",
                "description": "Sample keywords, e.g. 'D2O water neutron'.",
            },
            "sources": {
                "type": "array",
                "items": {"type": "string", "enum": ["local", "dissolve", "zenodo"]},
                "description": "Default: all three.",
            },
            "limit": {"type": "integer", "description": "Default 8."},
        },
        ["query"],
    )
    add(
        "download_data",
        "Download a data file (http/https URL) into the session's data folder; .mint01 files are added to the session.",
        {"url": {"type": "string"}},
        ["url"],
    )
    add(
        "finish",
        "End the conversation turn with a final message to the user.",
        {"message": {"type": "string"}},
        ["message"],
    )
    return registry
