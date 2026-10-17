"""Helpers shared by the wrapper tests: paths, output comparison, case replay.

The comparison reuses the masking rules of the Rust project's own tools
(``rust/tools/compare_outputs.py`` and ``rust/tools/make_reference_hashes.py``)
so that the wrapper is held to exactly the same standard as the ``epsr-rs`` CLI.
"""

from __future__ import annotations

import importlib.util
import os
import platform
import re
import shlex
import shutil
import tarfile
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType

import epsr

PYTHON_DIR = Path(__file__).resolve().parents[1]
REPO = PYTHON_DIR.parent
TEMPLATE = REPO / "EPSR-test"
TEST_PROBLEMS = REPO / "test_problems"
RUST_DIR = REPO / "rust"
INP = "200DMSO400Water.EPSR.inp"

#: Set EPSR_STRICT=1 to require byte-identical output everywhere (no numeric tolerance).
STRICT = os.environ.get("EPSR_STRICT", "") not in ("", "0")
#: Relative tolerance for files that are not byte-identical (last-digit libm rounding).
RTOL = 1e-5
IS_GLIBC = platform.system() == "Linux" and platform.libc_ver()[0] == "glibc"


def _load(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader, path
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


#: ``mask``, ``mask_more``, ``NUM`` and ``IGNORE`` of the Rust comparison script.
compare_outputs = _load("compare_outputs", RUST_DIR / "tools" / "compare_outputs.py")
#: ``mask`` and ``fnv1a64`` used to make ``reference_hashes.txt``.
reference_hashes = _load("make_reference_hashes", RUST_DIR / "tools" / "make_reference_hashes.py")


def copy_files(src: Path, dst: Path, skip: set[str] | None = None) -> Path:
    dst.mkdir(parents=True, exist_ok=True)
    for p in src.iterdir():
        if p.is_file() and not (skip and (p.name in skip or p.suffix in skip)):
            shutil.copy2(p, dst / p.name)
    return dst


# --------------------------------------------------------------------- comparison
@dataclass
class FileResult:
    name: str
    status: str  # "identical", "close", "different", "missing", "extra"
    detail: str = ""


def _numbers(line: str) -> list[float]:
    return [float(x.replace("D", "E").replace("d", "e")) for x in compare_outputs.NUM.findall(line)]


def compare_file(name: str, ref: Path, cand: Path) -> FileResult:
    """Compare one output file: identical after masking, numerically close, or different."""
    a, b = ref.read_bytes(), cand.read_bytes()
    if a == b:
        return FileResult(name, "identical")
    la = compare_outputs.mask_more(name, compare_outputs.mask(name, a.decode("latin-1").split("\n")))
    lb = compare_outputs.mask_more(name, compare_outputs.mask(name, b.decode("latin-1").split("\n")))
    if la == lb:
        return FileResult(name, "identical")
    if len(la) != len(lb):
        return FileResult(name, "different", f"line counts {len(la)} vs {len(lb)}")
    worst = 0.0
    for i, (x, y) in enumerate(zip(la, lb)):
        if x == y:
            continue
        if compare_outputs.NUM.sub("#", x) != compare_outputs.NUM.sub("#", y):
            return FileResult(name, "different", f"line {i + 1}: {x[:100]!r} vs {y[:100]!r}")
        for u, v in zip(_numbers(x), _numbers(y)):
            s = max(abs(u), abs(v))
            if s > 0:
                worst = max(worst, abs(u - v) / s)
    if worst <= RTOL:
        return FileResult(name, "close", f"max relative difference {worst:.2e}")
    return FileResult(name, "different", f"max relative difference {worst:.2e}")


@dataclass
class Comparison:
    results: list[FileResult] = field(default_factory=list)

    def by_status(self, *status: str) -> list[FileResult]:
        return [r for r in self.results if r.status in status]

    @property
    def ok(self) -> bool:
        bad = ("different", "missing", "extra") + (("close",) if STRICT else ())
        return not self.by_status(*bad)

    def report(self) -> str:
        return "\n".join(f"  {r.name:36s} {r.status} {r.detail}" for r in self.results if r.status != "identical")


def changed_files(pristine: Path, work: Path) -> set[str]:
    """Files in ``work`` that are new or differ from ``pristine``."""
    out = set()
    for p in work.iterdir():
        if not p.is_file():
            continue
        q = pristine / p.name
        if not q.exists() or q.read_bytes() != p.read_bytes():
            out.add(p.name)
    return out


def compare_dirs(ref: Path, cand: Path, pristine: Path | None = None) -> Comparison:
    """Compare every file of ``ref`` with ``cand``; with ``pristine``, also flag
    files the candidate changed that the reference did not."""
    comp = Comparison()
    ignore = getattr(compare_outputs, "IGNORE", set())
    for p in sorted(ref.iterdir()):
        if not p.is_file() or p.name in ignore:
            continue
        q = cand / p.name
        if not q.exists():
            comp.results.append(FileResult(p.name, "missing"))
        else:
            comp.results.append(compare_file(p.name, p, q))
    if pristine is not None:
        for name in sorted(changed_files(pristine, cand) - {p.name for p in ref.iterdir()} - set(ignore)):
            comp.results.append(FileResult(name, "extra", "written by the wrapper but not by the Fortran"))
    return comp


def hash_file(name: str, path: Path) -> str:
    return f"{reference_hashes.fnv1a64(reference_hashes.mask(name, path.read_bytes())):016x}"


# --------------------------------------------------------------------- case replay
@dataclass
class Command:
    command: str
    args: list[str]
    stdin: str | None  # file name in the case folder


@dataclass
class Case:
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
    cmds, repeat = [], 1
    for raw in (path / "commands.txt").read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("repeat"):
            repeat = int(line.split()[1])
            continue
        stdin = None
        if "<" in line:
            line, stdin = (s.strip() for s in line.split("<", 1))
        words = shlex.split(line)
        cmds.append(Command(words[0].lower(), words[1:], stdin))
    return Case(path, cmds, repeat)


def all_cases() -> list[Case]:
    found = [load_case(p.parent) for p in TEST_PROBLEMS.glob("*/EPSR-test-ex*/commands.txt")]
    return sorted(found, key=lambda c: (c.path.parent.name, c.number))


CASE_SKIP = {"commands.txt", "README.md", ".gz", ".pdf"}


def replay(case: Case, work: Path) -> list[str]:
    """Run every command of ``case`` through the Python wrapper in ``work``.

    ``epsr`` lines become ``epsr.Simulation(...).run()``; every other program goes
    through ``epsr.tools.run``. As in ``run_ex_test.sh``, a failing command is
    recorded and the replay continues (some cases deliberately exercise error
    exits of the Fortran programs). Returns the error messages.
    """
    errors = []
    for k in range(1, case.repeat + 1):
        for c in case.commands:
            answers = (case.path / c.stdin).read_text(encoding="latin-1") if c.stdin else ""
            try:
                if c.command == "epsr":
                    epsr.Simulation(work, c.args[0]).run()
                else:
                    epsr.tools.run(work, c.command, c.args, answers)
            except epsr.EpsrError as e:
                errors.append(f"pass {k} {c.command} {' '.join(c.args)}: {e}")
    return errors


def reference_dir(case: Case, dst: Path) -> Path:
    """Extract the case's ``fortran_reference.tar.gz`` into ``dst``."""
    dst.mkdir(parents=True, exist_ok=True)
    with tarfile.open(case.path / "fortran_reference.tar.gz") as t:
        t.extractall(dst, filter="data")
    return dst


def cli_binary() -> Path | None:
    """The ``epsr-rs`` CLI, built on demand with cargo if necessary."""
    exe = RUST_DIR / "target" / "release" / "epsr-rs"
    if exe.exists():
        return exe
    cargo = shutil.which("cargo") or str(Path.home() / ".cargo" / "bin" / "cargo")
    if not Path(cargo).exists():
        return None
    import subprocess

    r = subprocess.run(
        [cargo, "build", "--release", "-q", "-p", "epsr-cli", "--manifest-path", str(RUST_DIR / "Cargo.toml")],
        capture_output=True,
    )
    return exe if r.returncode == 0 and exe.exists() else None
