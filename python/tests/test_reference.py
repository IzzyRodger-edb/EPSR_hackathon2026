"""Port of ``rust/epsr-core/tests/test_case.rs`` to the Python wrapper.

Runs the engine through ``epsr.Simulation`` on a copy of ``EPSR-test`` and checks
every output file against the hashes of the original Fortran program's output
(``rust/epsr-core/tests/reference_hashes.txt``), the same check the Rust test
makes. A second test proves the wrapper drives the engine exactly like the
``epsr-rs`` CLI, byte for byte, which holds on any platform.

Platform note: the reference potential table ``.o01`` holds values of order
1e13 at very small r. Their last printed digit depends on the libm rounding of
``powf``, which differs between glibc versions, so on some machines that one
file differs from the Fortran hash in the 7th significant digit (the Rust test
fails the same way there). That file is reported as a warning rather than a
failure unless ``EPSR_STRICT=1``; ``test_ex_cases.py`` checks its numbers
against the full Fortran reference within a relative tolerance of 1e-5.
"""

from __future__ import annotations

import subprocess
import warnings
from pathlib import Path

import numpy as np
import pytest

import epsr

from conftest import CompletedRun
from epsr_cases import (
    IS_GLIBC,
    REQUIRE_IDENTICAL,
    RUST_DIR,
    TEMPLATE_INPUT_FILE,
    cli_binary,
    compare_dirs,
    copy_files,
    hash_file,
)

#: Files whose last printed digit is sensitive to the platform's libm (see module docstring).
LIBM_SENSITIVE_FILES = {"200DMSO400Water.EPSR.o01"}
EXPECTED_R_FACTOR = 1.1478e-4  # from rust/epsr-core/tests/test_case.rs


def read_reference_hashes() -> list[tuple[str, str]]:
    """``(hash, file name)`` rows of the Rust test's reference file."""
    text = (RUST_DIR / "epsr-core" / "tests" / "reference_hashes.txt").read_text()
    rows = []
    for line in text.splitlines():
        if line.strip() and not line.startswith("#"):
            expected_hash, file_name = line.split(" ", 1)
            rows.append((expected_hash, file_name))
    return rows


def test_outputs_match_fortran_reference(completed_run: CompletedRun) -> None:
    reference_rows = read_reference_hashes()
    assert len(reference_rows) == 26

    mismatches, libm_sensitive_mismatches = [], []
    for expected_hash, file_name in reference_rows:
        output_file = completed_run.project.folder / file_name
        assert output_file.exists(), f"missing output {file_name}"
        if hash_file(file_name, output_file) == expected_hash:
            continue
        if file_name in LIBM_SENSITIVE_FILES and not REQUIRE_IDENTICAL:
            libm_sensitive_mismatches.append(file_name)
        else:
            mismatches.append(file_name)

    if not IS_GLIBC:
        pytest.skip(
            "bit-for-bit agreement with the Fortran is only expected on Linux/glibc"
        )
    assert not mismatches, (
        f"{len(mismatches)} of {len(reference_rows)} files differ "
        f"from the Fortran reference: {mismatches}"
    )
    if libm_sensitive_mismatches:
        warnings.warn(
            f"{libm_sensitive_mismatches} differ from the Fortran hash "
            "(libm rounding; see test_reference.py docstring)",
            stacklevel=1,
        )


def test_api_sanity(completed_run: CompletedRun) -> None:
    simulation = completed_run.simulation
    assert simulation.n_pairs == 21
    assert simulation.n_datasets == 3
    assert len(simulation.q) == 600
    assert simulation.summary.r_factor == pytest.approx(EXPECTED_R_FACTOR, abs=1e-7)


def test_overrides_are_applied(project: epsr.Project) -> None:
    """Port of the Rust test ``overrides_are_applied``."""
    simulation = project.simulation(["ntimes 2", "data 2", "tweak 1.05"])
    assert simulation.control.ntimes == 2
    assert simulation.get_parameter("tweak", 1) == "1.05"
    assert simulation.get_parameter("tweak", 0) == "1.0"
    assert simulation.datasets[1].tweak == pytest.approx(1.05, abs=1e-6)


def test_wrapper_matches_cli_byte_for_byte(template: Path, tmp_path: Path) -> None:
    """The wrapper and ``epsr-rs`` produce identical files from identical inputs."""
    cli_executable = cli_binary()
    if cli_executable is None:
        pytest.skip("epsr-rs CLI not available (cargo not found)")
    wrapper_folder = copy_files(template, tmp_path / "wrapper")
    cli_folder = copy_files(template, tmp_path / "cli")
    pristine = copy_files(template, tmp_path / "pristine")

    for _ in range(2):  # two chained runs
        epsr.Simulation(wrapper_folder, TEMPLATE_INPUT_FILE).run()
        subprocess.run(
            [
                str(cli_executable),
                str(cli_folder),
                "epsr",
                TEMPLATE_INPUT_FILE,
                "--quiet",
            ],
            check=True,
            capture_output=True,
        )

    comparison = compare_dirs(cli_folder, wrapper_folder, pristine)
    assert comparison.all_identical, comparison.report()


def test_iterations_equal_run(template: Path, tmp_path: Path) -> None:
    """Driving the iterations from Python writes the same files as ``run()``."""
    run_folder = copy_files(template, tmp_path / "run")
    iterate_folder = copy_files(template, tmp_path / "iterate")

    summary_from_run = epsr.Simulation(run_folder, overrides={"niter": 2}).run()
    stepped_simulation = epsr.Simulation(iterate_folder, overrides={"niter": 2})
    summaries = list(stepped_simulation.iterations())

    assert [summary.iteration for summary in summaries] == [1, 2]
    assert summaries[-1] == summary_from_run
    comparison = compare_dirs(run_folder, iterate_folder)
    assert comparison.all_identical, comparison.report()


def test_runs_are_deterministic(template: Path, tmp_path: Path) -> None:
    first_folder = copy_files(template, tmp_path / "first")
    second_folder = copy_files(template, tmp_path / "second")

    first_summary = epsr.Simulation(first_folder).run()
    second_summary = epsr.Simulation(second_folder).run()

    assert first_summary == second_summary
    for file_name in ("200DMSO400Water.pcof", "200DMSO400Water.ato"):
        assert (first_folder / file_name).read_bytes() == (
            second_folder / file_name
        ).read_bytes()
    np.testing.assert_array_equal(
        epsr.Simulation(first_folder).gr_inter(),
        epsr.Simulation(second_folder).gr_inter(),
    )
