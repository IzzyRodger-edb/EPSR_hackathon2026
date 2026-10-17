"""Replay every ``test_problems/*/EPSR-test-ex*`` verification case through the wrapper.

This is the Python counterpart of ``rust/tools/run_all_ex_tests.sh``. Each case
folder holds input files, ``commands.txt`` (one program call per line,
``repeat N`` runs the list N times) and ``fortran_reference.tar.gz``, the files
the original Fortran programs wrote. Here every ``epsr`` line runs through
``epsr.Simulation(...).run()`` and every other program through
``epsr.tools.run``; then every reference file is compared with the wrapper's.

A file passes when it is identical after the standard masking (date line of
``.out``, timing columns of ``fort.35``, see ``rust/tools/compare_outputs.py``).
Without ``EPSR_STRICT=1`` a file also passes when it differs only in the last
digits of numbers: within a relative 1e-5, or one unit in the last printed
digit, or 1e-5 of the largest value in the same row or column (for numbers that
are zero by symmetry and come out as float32 noise). Such differences come
from libm and float rounding that varies between platforms; to prove they are
not caused by the wrapper, every such case is also run through the ``epsr-rs``
CLI and the wrapper's files must be byte-identical to the CLI's. Files written
by the wrapper but not by the Fortran fail the test.

A few cases exercise Fortran bugs that the Rust port fixes on purpose, so their
engine output differs from the Fortran reference by design (``FORTRAN_BUGS.md`` at
the repository root; ``rust/docs/VALIDATION.md``, "Fortran bugs deliberately FIXED
in the port"). The engine team marks each such case with an "Expected divergence
from the Fortran" section in its ``README.md``. For those cases the wrapper's files
must be byte-identical to the ``epsr-rs`` CLI's instead, and the test checks that
the divergence from the Fortran is really present.

Select a subset with ``-k``, e.g. ``pytest tests/test_ex_cases.py -k "ex20 or ex50"``;
skip the whole module with ``-m "not slow"``.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import pytest

from epsr_cases import (
    CASE_FILES_TO_SKIP,
    Case,
    Comparison,
    all_cases,
    cli_binary,
    compare_dirs,
    copy_files,
    reference_dir,
    replay,
    replay_cli,
)

CASES = all_cases()
#: Heading the engine team adds to the README of a case that deliberately differs.
DIVERGENCE_HEADING = "Expected divergence from the Fortran"


def is_documented_divergence(case: Case) -> bool:
    readme = case.path / "README.md"
    return readme.exists() and DIVERGENCE_HEADING in readme.read_text()


def _compare_with_cli(
    case: Case, wrapper_output: Path, pristine: Path, tmp_path: Path
) -> Comparison | None:
    """Replay the case with the ``epsr-rs`` CLI and compare its files with the wrapper's."""
    cli_executable = cli_binary()
    if cli_executable is None:
        return None
    cli_output = copy_files(case.path, tmp_path / "cli", CASE_FILES_TO_SKIP)
    replay_cli(case, cli_output, cli_executable)
    return compare_dirs(cli_output, wrapper_output, pristine)


@pytest.mark.slow
@pytest.mark.parametrize("case", CASES, ids=[case.name for case in CASES])
def test_case_matches_fortran(case: Case, tmp_path: Path) -> None:
    wrapper_output = copy_files(case.path, tmp_path / "work", CASE_FILES_TO_SKIP)
    pristine = copy_files(case.path, tmp_path / "pristine", CASE_FILES_TO_SKIP)
    command_errors = replay(case, wrapper_output)
    fortran_reference = reference_dir(case, tmp_path / "reference")
    versus_fortran = compare_dirs(fortran_reference, wrapper_output, pristine)

    if is_documented_divergence(case):
        # The engine fixes a Fortran bug here on purpose: the wrapper must still write
        # exactly what the CLI writes, and the documented divergence must be present.
        versus_cli = _compare_with_cli(case, wrapper_output, pristine, tmp_path)
        if versus_cli is None:
            pytest.skip("the epsr-rs CLI is needed to check a documented divergence")
        assert (
            versus_cli.all_identical
        ), f"{case.name}: wrapper output differs from the epsr-rs CLI:\n{versus_cli.report()}"
        assert not versus_fortran.all_identical, (
            f"{case.name}: README documents a divergence from the Fortran, but the "
            "output is identical to the Fortran reference"
        )
        return

    close_files = versus_fortran.with_status("close")
    if versus_fortran.passed and not close_files:
        return  # byte-identical to the Fortran (after masking)

    # Not byte-identical to the Fortran: check whether the epsr-rs CLI writes exactly the
    # same files as the wrapper. If so, the difference lies in the engine, not the wrapper.
    versus_cli = _compare_with_cli(case, wrapper_output, pristine, tmp_path)
    if versus_cli is not None:
        assert (
            versus_cli.all_identical
        ), f"{case.name}: wrapper output differs from the epsr-rs CLI:\n{versus_cli.report()}"
        origin = "the epsr-rs CLI writes the same files"
    else:
        origin = "CLI not checked"

    detail = versus_fortran.report()
    if command_errors:
        detail += "\ncommand errors:\n  " + "\n  ".join(command_errors)
    assert (
        versus_fortran.passed
    ), f"{case.name} differs from the Fortran reference ({origin}):\n{detail}"
    warnings.warn(
        f"{case.name}: {len(close_files)} file(s) differ from the Fortran in the last digits "
        f"({origin}):\n{versus_fortran.report()}",
        stacklevel=1,
    )


def test_all_cases_found() -> None:
    assert len(CASES) == 115
    assert all((case.path / "fortran_reference.tar.gz").exists() for case in CASES)
    divergent = sorted(case.name for case in CASES if is_documented_divergence(case))
    assert divergent == [
        "EPSR-test-ex108",
        "EPSR-test-ex110",
        "EPSR-test-ex111",
        "EPSR-test-ex117",
        "EPSR-test-ex118",
    ], "the documented divergences changed; update docs/testing.md"
