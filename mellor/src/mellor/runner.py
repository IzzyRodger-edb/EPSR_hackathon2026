"""Running a specification: write the files, execute the generated script, collect the results.

The script runs in a subprocess (same interpreter) with the run folder as working
directory. This keeps the EPSR setup programs, which change the process's working
directory, out of the agent's process, and makes the run reproducible from the
command line: ``python run_<name>.py --out .``.

What :func:`run_spec` leaves in the folder (docs/architecture.md §3): ``data/`` with
copies of the measured files (and the PDB file), ``spec.json``, ``run_<name>.py``,
``run.log`` with the script's console output, and what the script writes:
``epsr/``, ``trajectory/`` and ``manifest.json``.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import threading
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any

from . import inputs as inputs_module
from . import script as script_module
from .spec import SpecError, SystemSpec

__all__ = [
    "RunResult",
    "run_spec",
    "load_result",
    "RunError",
    "RUN_PRODUCTS",
    "format_runs_table",
]

#: Everything a run (and the analysis of it) writes besides the data copies, the
#: specification and the script. ``run_spec(..., overwrite=True)`` removes these
#: before running, so that no result of an earlier run survives next to a new one.
RUN_PRODUCTS = (
    "epsr",
    "trajectory",
    "figures",
    "tables",
    "manifest.json",
    "run.log",
    "analysis.json",
    "report.pdf",
)

#: Where the script records the trajectory when the manifest does not say.
DEFAULT_TRAJECTORY_FOLDER = "trajectory"


class RunError(RuntimeError):
    """The script failed; ``log`` holds its output."""

    def __init__(self, message: str, log: str = ""):
        super().__init__(message)
        self.log = log


@dataclass
class RunResult:
    """A finished (or loaded) run folder; see docs/architecture.md §3."""

    folder: Path
    spec: SystemSpec
    manifest: dict[str, Any]
    duration_s: float = 0.0
    quick: bool = False

    def __post_init__(self) -> None:
        self.folder = Path(self.folder)

    @property
    def spec_path(self) -> Path:
        return self.folder / "spec.json"

    @property
    def script_path(self) -> Path:
        return self.folder / f"run_{self.spec.name}.py"

    @property
    def epsr_folder(self) -> Path:
        return self.folder / self.manifest.get("epsr_folder", "epsr")

    @property
    def trajectory_path(self) -> Path | None:
        """The trajectory folder, or None when the run was not recorded."""
        name = self.manifest.get("trajectory", DEFAULT_TRAJECTORY_FOLDER)
        if not name:
            return None
        return self.folder / str(name)

    @property
    def has_trajectory(self) -> bool:
        """True when a recorded trajectory (``meta.json``) is present in the folder."""
        path = self.trajectory_path
        return path is not None and (path / "meta.json").is_file()

    @property
    def log_path(self) -> Path:
        return self.folder / "run.log"

    @property
    def data_folder(self) -> Path:
        return self.folder / "data"

    @property
    def runs(self) -> list[dict[str, Any]]:
        return list(self.manifest.get("runs", []))

    @property
    def final_summary(self) -> dict[str, Any]:
        return dict(self.manifest.get("final", {}))

    @property
    def complete(self) -> bool:
        """True when the script finished every run of the schedule."""
        return bool(self.manifest.get("complete", bool(self.runs)))

    def table(self) -> str:
        """Plain-text table of the runs: stage, iterations, ereq, R-factor, energy, pressure, EP amplitude."""
        return format_runs_table(self.runs)


# --------------------------------------------------------------------------- the table
def _cell(value: Any, spec: str, width: int) -> str:
    """``value`` formatted with ``spec``; a dash when it is missing or not a number."""
    if value is None:
        return f"{'-':>{width}}"
    try:
        return f"{float(value):{spec}}"
    except (TypeError, ValueError):
        return f"{str(value):>{width}}"


def format_runs_table(runs: list[dict[str, Any]]) -> str:
    """The run summaries of a manifest as a plain-text table (one line per run)."""
    header = (
        f"{'run':>3}  {'stage':<20} {'niter':>5} {'ereq':>6} {'R-factor':>10}"
        f" {'energy':>10} {'pressure':>10} {'EP ampl':>8} {'time':>9}"
    )
    lines = [header]
    for run in runs:
        summary = run.get("summary") or {}
        lines.append(
            f"{_cell(run.get('run'), '>3.0f', 3)}  {str(run.get('stage', '')):<20}"
            f" {_cell(run.get('niter'), '>5.0f', 5)} {_cell(run.get('ereq'), '>6g', 6)}"
            f" {_cell(summary.get('r_factor'), '>10.3e', 10)}"
            f" {_cell(summary.get('energy'), '>10.4g', 10)}"
            f" {_cell(summary.get('pressure'), '>10.4g', 10)}"
            f" {_cell(summary.get('ep_amplitude'), '>8.3f', 8)}"
            f" {_cell(run.get('duration_s'), '>8.1f', 8)}s"
        )
    if not runs:
        lines.append("(no runs)")
    return "\n".join(lines)


# --------------------------------------------------------------------------- running
def run_spec(
    spec: SystemSpec,
    folder: str | Path,
    *,
    quick: bool = False,
    record: bool = True,
    overwrite: bool = False,
    timeout_s: float | None = None,
    progress: Callable[[str], None] | None = None,
) -> RunResult:
    """Run ``spec`` in ``folder`` and return the result.

    Steps: validate the spec (``check_files=True``); create the folder (refuse a
    non-empty one unless ``overwrite``); copy the data files to ``folder/data``
    (and the PDB file); write ``spec.json``; generate the inputs and write the script;
    run ``python run_<name>.py --out . [--quick] [--no-record]`` with the folder as
    working directory, streaming its output to ``run.log`` (and to ``progress`` line by
    line when given); on failure raise :class:`RunError` with the log; load
    ``manifest.json``.

    Raises :class:`~mellor.spec.SpecError` for an invalid specification or a missing
    data file, ``FileExistsError`` for a non-empty folder without ``overwrite``, and
    :class:`RunError` when the script fails, exceeds ``timeout_s`` or writes no
    manifest. With ``overwrite`` the products of an earlier run (:data:`RUN_PRODUCTS`)
    are removed first; the data copies and other files in the folder are kept.
    """
    warnings = spec.validate(check_files=True)
    folder = Path(os.path.abspath(Path(folder).expanduser()))
    _prepare_folder(folder, overwrite)
    generated = inputs_module.generate(spec)
    _copy_inputs(generated, folder / "data")
    spec.save(folder / "spec.json")
    script_path = script_module.write_script(spec, generated, folder)

    command = [sys.executable, "-I", script_path.name, "--out", "."]
    if quick:
        command.append("--quick")
    if not record:
        command.append("--no-record")
    preamble = [f"$ {' '.join(command)}  (in {folder})"]
    preamble += [f"warning: {warning}" for warning in warnings]

    started = time.perf_counter()
    log_text, returncode = _execute(
        command, folder, folder / "run.log", preamble, progress, timeout_s
    )
    duration = round(time.perf_counter() - started, 3)
    if returncode != 0:
        raise RunError(
            f"{script_path.name} failed with exit status {returncode} after "
            f"{duration:.1f} s; the output is in {folder / 'run.log'}",
            log=log_text,
        )
    manifest_path = folder / "manifest.json"
    if not manifest_path.is_file():
        raise RunError(
            f"{script_path.name} exited 0 but wrote no {manifest_path.name} in {folder}",
            log=log_text,
        )
    manifest = _read_json(manifest_path)
    return RunResult(
        folder=folder, spec=spec, manifest=manifest, duration_s=duration, quick=quick
    )


def load_result(folder: str | Path) -> RunResult:
    """A :class:`RunResult` for an existing run folder (``spec.json`` + ``manifest.json``)."""
    folder = Path(os.path.abspath(Path(folder).expanduser()))
    if not folder.is_dir():
        raise FileNotFoundError(f"run folder {folder} does not exist")
    spec_path = folder / "spec.json"
    manifest_path = folder / "manifest.json"
    for path in (spec_path, manifest_path):
        if not path.is_file():
            raise FileNotFoundError(
                f"{folder} is not a finished Mellor run folder: {path.name} is missing"
            )
    try:
        spec = SystemSpec.load(spec_path)
    except SpecError as error:
        raise SpecError(
            [f"{spec_path}: {problem}" for problem in error.errors]
        ) from error
    manifest = _read_json(manifest_path)
    duration = manifest.get("duration_s")
    return RunResult(
        folder=folder,
        spec=spec,
        manifest=manifest,
        duration_s=float(duration) if isinstance(duration, (int, float)) else 0.0,
        quick=bool(manifest.get("quick", False)),
    )


# --------------------------------------------------------------------------- helpers
def _read_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except ValueError as error:
        raise ValueError(f"{path} is not valid JSON: {error}") from error
    if not isinstance(data, dict):
        raise ValueError(f"{path} does not hold a JSON object")
    return data


def _prepare_folder(folder: Path, overwrite: bool) -> None:
    if folder.exists() and not folder.is_dir():
        raise NotADirectoryError(f"{folder} exists and is not a folder")
    folder.mkdir(parents=True, exist_ok=True)
    if not any(folder.iterdir()):
        return
    if not overwrite:
        raise FileExistsError(
            f"{folder} is not empty; choose another folder or pass overwrite=True "
            "(this removes the products of the earlier run: "
            + ", ".join(RUN_PRODUCTS)
            + ")"
        )
    for name in RUN_PRODUCTS:
        path = folder / name
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        elif path.exists() or path.is_symlink():
            path.unlink()


def _copy_inputs(generated: inputs_module.EpsrInputs, data_folder: Path) -> None:
    """Copy the measured data files (and the PDB file) into ``data_folder``."""
    sources = list(generated.data_files)
    if generated.pdb_file is not None:
        sources.append(generated.pdb_file)
    data_folder.mkdir(parents=True, exist_ok=True)
    for source in sources:
        if not source.is_file():
            raise SpecError([f"input file {source} not found"])
        target = data_folder / source.name
        if target.exists() and target.samefile(source):
            continue  # the specification already points at the copy in data/
        shutil.copy2(source, target)


def _pump(
    stream: IO[str],
    log: IO[str],
    lines: list[str],
    progress: Callable[[str], None] | None,
) -> None:
    """Copy the script's output, line by line, to the log file and ``progress``."""
    for raw in stream:
        line = raw.rstrip("\r\n")
        lines.append(line)
        log.write(line + "\n")
        log.flush()
        if progress is not None:
            try:
                progress(line)
            except Exception as error:  # a broken callback must not stop the run
                note = f"[runner] progress callback failed ({type(error).__name__}: {error}); not called again"
                lines.append(note)
                log.write(note + "\n")
                log.flush()
                progress = None


def _execute(
    command: list[str],
    cwd: Path,
    log_path: Path,
    preamble: list[str],
    progress: Callable[[str], None] | None,
    timeout_s: float | None,
) -> tuple[str, int]:
    """Run ``command`` in ``cwd``, streaming its output; ``(log text, exit status)``."""
    lines: list[str] = []
    with log_path.open("w", encoding="utf-8") as log:
        for line in preamble:
            lines.append(line)
            log.write(line + "\n")
            if progress is not None:
                progress(line)
        log.flush()
        process = subprocess.Popen(
            command,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        assert process.stdout is not None
        pump = threading.Thread(
            target=_pump, args=(process.stdout, log, lines, progress), daemon=True
        )
        pump.start()
        try:
            returncode = process.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
            pump.join()
            note = f"[runner] killed after {timeout_s} s (timeout)"
            lines.append(note)
            log.write(note + "\n")
            raise RunError(
                f"{command[2]} did not finish within {timeout_s} s and was killed; "
                f"the output is in {log_path}",
                log="\n".join(lines),
            ) from None
        pump.join()
    return "\n".join(lines), returncode
