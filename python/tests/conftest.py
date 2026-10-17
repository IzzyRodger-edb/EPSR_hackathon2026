"""Shared fixtures for the wrapper tests."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import NamedTuple

import pytest

import epsr

from epsr_cases import TEMPLATE, TEMPLATE_INPUT_FILE


class CompletedRun(NamedTuple):
    """A copy of EPSR-test on which one EPSR run has finished."""

    project: epsr.Project
    simulation: epsr.Simulation


@pytest.fixture(scope="session")
def template() -> Path:
    if not (TEMPLATE / TEMPLATE_INPUT_FILE).exists():
        pytest.skip(f"test case folder {TEMPLATE} not found")
    return TEMPLATE


@pytest.fixture
def project(template: Path, tmp_path: Path) -> epsr.Project:
    """A fresh copy of EPSR-test (function scope: safe to modify)."""
    return epsr.Project.from_template(template, tmp_path / "work")


@pytest.fixture(scope="module")
def completed_run(template: Path, tmp_path_factory: pytest.TempPathFactory) -> CompletedRun:
    """One completed EPSR run on a copy of EPSR-test, shared by a test module (read-only use)."""
    run_project = epsr.Project.from_template(template, tmp_path_factory.mktemp("run") / "work")
    simulation = run_project.simulation()
    simulation.run()
    return CompletedRun(run_project, simulation)


def _file_contents(folder: Path) -> dict[str, bytes]:
    return {path.name: path.read_bytes() for path in sorted(folder.iterdir()) if path.is_file()}


@pytest.fixture(scope="session", autouse=True)
def template_is_never_modified(template: Path) -> Iterator[None]:
    """Fail the session if any test writes into the repository's EPSR-test folder."""
    contents_before = _file_contents(template)
    yield
    contents_after = _file_contents(template)
    all_names = contents_before.keys() | contents_after.keys()
    changed = sorted(
        name for name in all_names if contents_before.get(name) != contents_after.get(name)
    )
    assert not changed, f"tests modified the template folder {template}: {changed}"
