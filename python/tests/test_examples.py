"""The scripts in ``python/examples`` run end to end."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from epsr_cases import PYTHON_DIR

EXAMPLE_SCRIPTS = sorted((PYTHON_DIR / "examples").glob("0*.py"))
#: Command-line arguments that keep the longer examples short under test.
SHORT_RUN_ARGUMENTS = {
    "02_refinement_loop": ["3"],  # number of EPSR runs
    "05_experimental_water": ["--quick"],  # 10 iterations per EPSR run
}


@pytest.mark.parametrize(
    "script", EXAMPLE_SCRIPTS, ids=[script.stem for script in EXAMPLE_SCRIPTS]
)
def test_example_runs(script: Path, tmp_path: Path) -> None:
    environment = {
        **os.environ,
        "EPSR_EXAMPLES_WORK": str(tmp_path),
        "MPLBACKEND": "Agg",
    }
    command = [sys.executable, str(script), *SHORT_RUN_ARGUMENTS.get(script.stem, [])]
    completed = subprocess.run(
        command,
        cwd=PYTHON_DIR,
        env=environment,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert completed.stdout.strip()


def test_examples_found() -> None:
    assert [script.stem for script in EXAMPLE_SCRIPTS] == [
        "01_single_run",
        "02_refinement_loop",
        "03_perturb_potential",
        "04_tools",
        "05_experimental_water",
    ]
