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

import shutil
import subprocess
import warnings
from pathlib import Path

import numpy as np
import pytest
from epsr_cases import INP, IS_GLIBC, RUST_DIR, STRICT, cli_binary, compare_dirs, copy_files, hash_file

import epsr

#: Files whose last printed digit is sensitive to the platform's libm (see module docstring).
LIBM_SENSITIVE = {"200DMSO400Water.EPSR.o01"}


def reference_hashes() -> list[tuple[str, str]]:
    text = (RUST_DIR / "epsr-core" / "tests" / "reference_hashes.txt").read_text()
    rows = [line.split(" ", 1) for line in text.splitlines() if line.strip() and not line.startswith("#")]
    return [(h, name) for h, name in rows]


def test_outputs_match_fortran_reference(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    proj, _ = ran
    failures, sensitive = [], []
    rows = reference_hashes()
    assert len(rows) == 26
    for expected, name in rows:
        path = proj.folder / name
        assert path.exists(), f"missing output {name}"
        if hash_file(name, path) != expected:
            (sensitive if name in LIBM_SENSITIVE and not STRICT else failures).append(name)
    if not IS_GLIBC:
        pytest.skip("bit-for-bit agreement with the Fortran is only expected on Linux/glibc")
    assert not failures, f"{len(failures)} of {len(rows)} files differ from the Fortran reference: {failures}"
    if sensitive:
        warnings.warn(
            f"{sensitive} differ from the Fortran hash (libm rounding; see test_reference.py docstring)",
            stacklevel=1,
        )


def test_api_sanity(ran: tuple[epsr.Project, epsr.Simulation]) -> None:
    _, sim = ran
    assert sim.n_pairs == 21
    assert sim.n_datasets == 3
    assert len(sim.q) == 600
    assert abs(sim.summary.r_factor - 1.1478e-4) < 1e-7, sim.summary.r_factor


def test_overrides_are_applied(project: epsr.Project) -> None:
    """Port of the Rust test ``overrides_are_applied``."""
    sim = project.simulation(["ntimes 2", "data 2", "tweak 1.05"])
    assert sim.control.ntimes == 2
    assert sim.get_parameter("tweak", 1) == "1.05"
    assert sim.get_parameter("tweak", 0) == "1.0"
    assert abs(sim.datasets[1].tweak - 1.05) < 1e-6


def test_wrapper_matches_cli_byte_for_byte(template: Path, tmp_path: Path) -> None:
    """The wrapper and ``epsr-rs`` produce identical files from identical inputs."""
    exe = cli_binary()
    if exe is None:
        pytest.skip("epsr-rs CLI not available (cargo not found)")
    a = copy_files(template, tmp_path / "wrapper")
    b = copy_files(template, tmp_path / "cli")
    pristine = copy_files(template, tmp_path / "pristine")
    for _ in range(2):  # two chained runs
        epsr.Simulation(a, INP).run()
        subprocess.run([str(exe), str(b), "epsr", INP, "--quiet"], check=True, capture_output=True)
    comp = compare_dirs(b, a, pristine)
    not_identical = [r for r in comp.results if r.status != "identical"]
    assert not not_identical, comp.report()


def test_iterations_equal_run(template: Path, tmp_path: Path) -> None:
    """Driving the iterations from Python writes the same files as ``run()``."""
    a = copy_files(template, tmp_path / "run")
    b = copy_files(template, tmp_path / "iterate")
    s_run = epsr.Simulation(a, INP, overrides={"niter": 2}).run()
    sim = epsr.Simulation(b, INP, overrides={"niter": 2})
    summaries = list(sim.iterations())
    assert [s.iteration for s in summaries] == [1, 2]
    assert summaries[-1] == s_run
    comp = compare_dirs(a, b)
    assert all(r.status == "identical" for r in comp.results), comp.report()


def test_runs_are_deterministic(template: Path, tmp_path: Path) -> None:
    a = copy_files(template, tmp_path / "a")
    b = copy_files(template, tmp_path / "b")
    sa = epsr.Simulation(a).run()
    sb = epsr.Simulation(b).run()
    assert sa == sb
    assert (a / "200DMSO400Water.pcof").read_bytes() == (b / "200DMSO400Water.pcof").read_bytes()
    assert (a / "200DMSO400Water.ato").read_bytes() == (b / "200DMSO400Water.ato").read_bytes()
    np.testing.assert_array_equal(epsr.Simulation(a).gr_inter(), epsr.Simulation(b).gr_inter())
    shutil.rmtree(a)
