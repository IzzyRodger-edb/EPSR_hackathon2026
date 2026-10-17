"""Shared fixtures for the wrapper tests."""

from __future__ import annotations

from pathlib import Path

import pytest
from epsr_cases import INP, TEMPLATE

import epsr


@pytest.fixture(scope="session")
def template() -> Path:
    if not (TEMPLATE / INP).exists():
        pytest.skip(f"test case folder {TEMPLATE} not found")
    return TEMPLATE


@pytest.fixture
def project(template: Path, tmp_path: Path) -> epsr.Project:
    """A fresh copy of EPSR-test (function scope: safe to modify)."""
    return epsr.Project.from_template(template, tmp_path / "work")


@pytest.fixture(scope="module")
def ran(template: Path, tmp_path_factory: pytest.TempPathFactory) -> tuple[epsr.Project, epsr.Simulation]:
    """One completed EPSR run on a copy of EPSR-test, shared by a test module (read-only use)."""
    proj = epsr.Project.from_template(template, tmp_path_factory.mktemp("ran") / "work")
    sim = proj.simulation()
    sim.run()
    return proj, sim
