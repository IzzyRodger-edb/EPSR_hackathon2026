"""The scripts in ``python/examples`` run end to end."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from epsr_cases import PYTHON_DIR

EXAMPLE_SCRIPTS = sorted((PYTHON_DIR / "examples").glob("0*.py"))
#: Number of runs passed to the refinement-loop example to keep the test short.
REFINEMENT_RUNS = "3"


@pytest.mark.parametrize(
    "script", EXAMPLE_SCRIPTS, ids=[script.stem for script in EXAMPLE_SCRIPTS]
)
def test_example_runs(script: Path, tmp_path: Path) -> None:
    environment = {
        **os.environ,
        "EPSR_EXAMPLES_WORK": str(tmp_path),
        "MPLBACKEND": "Agg",
    }
    command = [sys.executable, str(script)]
    if "refinement" in script.stem:
        command.append(REFINEMENT_RUNS)
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
    ]
