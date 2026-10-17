"""Helpers shared by the wrapper tests: paths, output comparison, case replay.

The comparison reuses the masking rules of the Rust project's own tools
(``rust/tools/compare_outputs.py`` and ``rust/tools/make_reference_hashes.py``)
so that the wrapper is held to exactly the same standard as the ``epsr-rs`` CLI.
"""

from __future__ import annotations

import functools
import importlib.util
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType

import epsr

PYTHON_DIR = Path(__file__).resolve().parents[1]
REPOSITORY = PYTHON_DIR.parent
TEMPLATE = REPOSITORY / "EPSR-test"
TEST_PROBLEMS = REPOSITORY / "test_problems"
RUST_DIR = REPOSITORY / "rust"
TEMPLATE_INPUT_FILE = "200DMSO400Water.EPSR.inp"

#: Set EPSR_STRICT=1 to require byte-identical output everywhere (no numeric tolerance).
REQUIRE_IDENTICAL = os.environ.get("EPSR_STRICT", "") not in ("", "0")
#: Tolerances for files that are not byte-identical (last-digit libm rounding):
#: relative, and absolute as a fraction of the largest value in the same row or column.
RELATIVE_TOLERANCE = 1e-5
ABSOLUTE_TOLERANCE_FRACTION = 1e-5
IS_GLIBC = platform.system() == "Linux" and platform.libc_ver()[0] == "glibc"

#: Files of a case folder that are not inputs of the case.
CASE_FILES_TO_SKIP = {"commands.txt", "README.md", ".gz", ".pdf"}


def _load_rust_tool(module_name: str, path: Path) -> ModuleType:
    """Import a script from ``rust/tools`` without writing ``__pycache__`` into ``rust/``."""
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec and spec.loader, path
    module = importlib.util.module_from_spec(spec)
    previous_setting, sys.dont_write_bytecode = sys.dont_write_bytecode, True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous_setting
    return module


#: ``mask``, ``mask_more``, ``NUM`` and ``IGNORE`` of the Rust comparison script.
compare_outputs = _load_rust_tool(
    "compare_outputs", RUST_DIR / "tools" / "compare_outputs.py"
)
#: ``mask`` and ``fnv1a64`` used to make ``reference_hashes.txt``.
reference_hashes = _load_rust_tool(
    "make_reference_hashes", RUST_DIR / "tools" / "make_reference_hashes.py"
)


def copy_files(source: Path, destination: Path, skip: set[str] | None = None) -> Path:
    """Copy the regular files of ``source``, except names or suffixes listed in ``skip``."""
    destination.mkdir(parents=True, exist_ok=True)
    skip = skip or set()
    for source_file in source.iterdir():
        if (
            source_file.is_file()
            and source_file.name not in skip
            and source_file.suffix not in skip
        ):
            shutil.copy2(source_file, destination / source_file.name)
    return destination


# --------------------------------------------------------------------- comparison
@dataclass
class FileResult:
    name: str
    status: str  # "identical", "close", "different", "missing", "extra"
    detail: str = ""


def _as_number(token: str) -> float | None:
    """The value of a Fortran number token (``1.5E+00``, ``1.5D+00``...), or None."""
    if not compare_outputs.NUM.fullmatch(token):
        return None
    return float(token.replace("D", "E").replace("d", "e"))


def _largest_magnitude_per_column(lines: list[str]) -> dict[int, float]:
    """Largest magnitude of the numbers at each token position over a file."""
    largest: dict[int, float] = {}
    for line in lines:
        for column, token in enumerate(line.split()):
            value = _as_number(token)
            if value is not None:
                largest[column] = max(largest.get(column, 0.0), abs(value))
    return largest


def _last_digit_unit(token: str) -> float:
    """One unit in the last printed digit of a real token (``0.21800E-02`` gives 1e-7).

    Integers return 0: they must match exactly.
    """
    text = token.lstrip("+-").replace("D", "E").replace("d", "e")
    if "." not in text and "E" not in text and "e" not in text:
        return 0.0
    mantissa, _, exponent = text.partition("E") if "E" in text else text.partition("e")
    n_decimals = len(mantissa.split(".", 1)[1]) if "." in mantissa else 0
    return 10.0 ** (int(exponent or 0) - n_decimals)


def _line_difference(
    reference_line: str, candidate_line: str, column_magnitudes: dict[int, float]
) -> float | None:
    """How far two lines are from agreeing, in units of the tolerance.

    Returns a value <= 1 when the lines agree within tolerance, and None when
    their text (not just their numbers) differs. Lines are compared as
    whitespace-separated tokens because fixed-width Fortran fields shift when a
    sign appears (``-0.00000`` vs ``0.00000``). Two numbers agree when their
    difference is at most the largest of

    * ``RELATIVE_TOLERANCE`` times the larger magnitude (relative rounding),
    * one unit in the last printed digit (output printed with few digits),
    * ``ABSOLUTE_TOLERANCE_FRACTION`` times the largest magnitude in the same row
      or column (cancellation noise on values that are zero by symmetry).
    """
    reference_tokens, candidate_tokens = reference_line.split(), candidate_line.split()
    if len(reference_tokens) != len(candidate_tokens):
        return None
    row_values = [_as_number(token) for token in reference_tokens[1:]]
    row_magnitude = max(
        (abs(value) for value in row_values if value is not None), default=0.0
    )

    worst_ratio = 0.0
    token_pairs = zip(reference_tokens, candidate_tokens, strict=True)
    for column, (reference_token, candidate_token) in enumerate(token_pairs):
        if reference_token == candidate_token:
            continue
        reference_value = _as_number(reference_token)
        candidate_value = _as_number(candidate_token)
        if reference_value is None or candidate_value is None:
            return None
        if reference_value == candidate_value:  # e.g. -0.0 and 0.0
            continue
        tolerance = max(
            RELATIVE_TOLERANCE * max(abs(reference_value), abs(candidate_value)),
            _last_digit_unit(reference_token),
            _last_digit_unit(candidate_token),
            ABSOLUTE_TOLERANCE_FRACTION
            * max(column_magnitudes.get(column, 0.0), row_magnitude),
        )
        difference = abs(reference_value - candidate_value)
        ratio = difference / tolerance if tolerance > 0 else float("inf")
        worst_ratio = max(worst_ratio, ratio)
    return worst_ratio


def _masked_lines(name: str, content: bytes) -> list[str]:
    """File lines with the run-dependent parts masked (date, timings), as the Rust tools do."""
    lines = content.decode("latin-1").split("\n")
    return compare_outputs.mask_more(name, compare_outputs.mask(name, lines))


def compare_file(name: str, reference: Path, candidate: Path) -> FileResult:
    """Compare one output file: identical after masking, numerically close, or different."""
    reference_content, candidate_content = (
        reference.read_bytes(),
        candidate.read_bytes(),
    )
    if reference_content == candidate_content:
        return FileResult(name, "identical")
    reference_lines = _masked_lines(name, reference_content)
    candidate_lines = _masked_lines(name, candidate_content)
    if reference_lines == candidate_lines:
        return FileResult(name, "identical")
    if len(reference_lines) != len(candidate_lines):
        return FileResult(
            name,
            "different",
            f"line counts {len(reference_lines)} vs {len(candidate_lines)}",
        )

    column_magnitudes = _largest_magnitude_per_column(reference_lines)
    worst_ratio, n_different_lines = 0.0, 0
    line_pairs = zip(reference_lines, candidate_lines, strict=True)
    for line_index, (reference_line, candidate_line) in enumerate(line_pairs):
        if reference_line == candidate_line:
            continue
        n_different_lines += 1
        ratio = _line_difference(reference_line, candidate_line, column_magnitudes)
        if ratio is None:
            detail = f"line {line_index + 1}: {reference_line[:100]!r} vs {candidate_line[:100]!r}"
            return FileResult(name, "different", detail)
        worst_ratio = max(worst_ratio, ratio)

    status = "close" if worst_ratio <= 1.0 else "different"
    detail = f"{n_different_lines} line(s) differ in the last digits ({worst_ratio:.2f} x tolerance)"
    return FileResult(name, status, detail)


@dataclass
class Comparison:
    results: list[FileResult] = field(default_factory=list)

    def with_status(self, *statuses: str) -> list[FileResult]:
        return [result for result in self.results if result.status in statuses]

    @property
    def passed(self) -> bool:
        failing_statuses = ["different", "missing", "extra"]
        if REQUIRE_IDENTICAL:
            failing_statuses.append("close")
        return not self.with_status(*failing_statuses)

    @property
    def all_identical(self) -> bool:
        return all(result.status == "identical" for result in self.results)

    def report(self) -> str:
        return "\n".join(
            f"  {result.name:36s} {result.status} {result.detail}"
            for result in self.results
            if result.status != "identical"
        )


def changed_files(pristine: Path, work: Path) -> set[str]:
    """Files in ``work`` that are new or differ from ``pristine``."""
    changed = set()
    for work_file in work.iterdir():
        if not work_file.is_file():
            continue
        pristine_file = pristine / work_file.name
        if (
            not pristine_file.exists()
            or pristine_file.read_bytes() != work_file.read_bytes()
        ):
            changed.add(work_file.name)
    return changed


def compare_dirs(
    reference: Path, candidate: Path, pristine: Path | None = None
) -> Comparison:
    """Compare every file of ``reference`` with ``candidate``.

    With ``pristine`` (the inputs before the run), files the candidate changed
    that the reference did not are reported as ``extra``.
    """
    comparison = Comparison()
    ignored_names = set(getattr(compare_outputs, "IGNORE", set()))
    reference_names = set()
    for reference_file in sorted(reference.iterdir()):
        if not reference_file.is_file() or reference_file.name in ignored_names:
            continue
        reference_names.add(reference_file.name)
        candidate_file = candidate / reference_file.name
        if candidate_file.exists():
            result = compare_file(reference_file.name, reference_file, candidate_file)
        else:
            result = FileResult(reference_file.name, "missing")
        comparison.results.append(result)

    if pristine is not None:
        unexpected = (
            changed_files(pristine, candidate) - reference_names - ignored_names
        )
        for name in sorted(unexpected):
            comparison.results.append(
                FileResult(
                    name, "extra", "written by the wrapper but not by the Fortran"
                )
            )
    return comparison


def hash_file(name: str, path: Path) -> str:
    """The masked FNV-1a hash used by ``rust/epsr-core/tests/reference_hashes.txt``."""
    masked_content = reference_hashes.mask(name, path.read_bytes())
    return f"{reference_hashes.fnv1a64(masked_content):016x}"


# --------------------------------------------------------------------- case replay
@dataclass
class Command:
    """One line of a case's ``commands.txt``."""

    program: str
    args: list[str]
    stdin_file: str | None  # file name in the case folder


@dataclass
class Case:
    """One ``EPSR-test-ex*`` verification case."""

    path: Path
    commands: list[Command]
    repeat: int

    @property
    def name(self) -> str:
        return self.path.name

    @property
    def number(self) -> int:
        return int(re.sub(r"\D", "", self.name.rsplit("ex", 1)[1]))


def load_case(path: Path) -> Case:
    """Parse ``commands.txt`` the way ``rust/tools/run_ex_test.sh`` does."""
    commands, repeat = [], 1
    for raw_line in (path / "commands.txt").read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("repeat"):
            repeat = int(line.split()[1])
            continue
        stdin_file = None
        if "<" in line:
            line, stdin_file = (part.strip() for part in line.split("<", 1))
        words = shlex.split(line)
        commands.append(
            Command(program=words[0].lower(), args=words[1:], stdin_file=stdin_file)
        )
    return Case(path, commands, repeat)


def all_cases() -> list[Case]:
    command_files = TEST_PROBLEMS.glob("*/EPSR-test-ex*/commands.txt")
    cases = [load_case(command_file.parent) for command_file in command_files]
    return sorted(cases, key=lambda case: (case.path.parent.name, case.number))


def replay(case: Case, work: Path) -> list[str]:
    """Run every command of ``case`` through the Python wrapper in ``work``.

    ``epsr`` lines become ``epsr.Simulation(...).run()``; every other program goes
    through ``epsr.tools.run``. As in ``run_ex_test.sh``, a failing command is
    recorded and the replay continues (some cases deliberately exercise error
    exits of the Fortran programs). Returns the error messages.
    """
    errors = []
    for pass_number in range(1, case.repeat + 1):
        for command in case.commands:
            answers = ""
            if command.stdin_file:
                answers = (case.path / command.stdin_file).read_text(encoding="latin-1")
            try:
                if command.program == "epsr":
                    epsr.Simulation(work, command.args[0]).run()
                else:
                    epsr.tools.run(work, command.program, command.args, answers)
            except epsr.EpsrError as error:
                arguments = " ".join(command.args)
                errors.append(
                    f"pass {pass_number} {command.program} {arguments}: {error}"
                )
    return errors


def replay_cli(case: Case, work: Path, cli_executable: Path) -> list[str]:
    """Run ``case`` with the ``epsr-rs`` CLI exactly as ``rust/tools/run_ex_test.sh`` does
    (from inside ``work``, folder argument ``./``)."""
    errors = []
    for pass_number in range(1, case.repeat + 1):
        for command in case.commands:
            stdin_path = (
                case.path / command.stdin_file if command.stdin_file else os.devnull
            )
            with open(stdin_path, "rb") as stdin:
                completed = subprocess.run(
                    [
                        str(cli_executable),
                        "./",
                        command.program,
                        *command.args,
                        "--quiet",
                    ],
                    stdin=stdin,
                    capture_output=True,
                    cwd=work,
                )
            if completed.returncode:
                errors.append(
                    f"pass {pass_number} {command.program}: exit {completed.returncode}"
                )
    return errors


def extract_archive(archive_path: Path, destination: Path) -> Path:
    """Extract a .tar.gz archive into ``destination`` (created if needed)."""
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path) as archive:
        archive.extractall(destination, filter="data")
    return destination


def reference_dir(case: Case, destination: Path) -> Path:
    """Extract the case's ``fortran_reference.tar.gz`` into ``destination``."""
    return extract_archive(case.path / "fortran_reference.tar.gz", destination)


@functools.cache
def cli_binary() -> Path | None:
    """The ``epsr-rs`` CLI, rebuilt with cargo so that it matches the engine the wrapper uses.

    The wrapper's extension is rebuilt by ``uv`` whenever the engine sources change;
    the CLI is not, so it is built here (a no-op when it is up to date). Without
    cargo, an existing binary is used as it is. Cached per process.
    """
    executable = RUST_DIR / "target" / "release" / "epsr-rs"
    cargo = shutil.which("cargo") or str(Path.home() / ".cargo" / "bin" / "cargo")
    if not Path(cargo).exists():
        return executable if executable.exists() else None
    build_command = [
        cargo,
        "build",
        "--release",
        "--quiet",
        "--package",
        "epsr-cli",
        "--manifest-path",
        str(RUST_DIR / "Cargo.toml"),
    ]
    build = subprocess.run(build_command, capture_output=True)
    return executable if build.returncode == 0 and executable.exists() else None
