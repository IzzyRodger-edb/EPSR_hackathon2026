#!/usr/bin/env python3
"""Run the Cu-BTC case with the original Fortran EPSR and archive the reference results.

usage: make_fortran_reference.py FORTRAN_BUILD_DIR OUTPUT_ARCHIVE

FORTRAN_BUILD_DIR is made by rust/tools/build_fortran_reference.sh. The input files
of this folder are copied to a temporary folder and the steps of steps.txt are run
there exactly as validation/run_validation.sh runs them (each program from inside
the folder, `OMP_NUM_THREADS=1`, unlimited stack; `set` steps through
validation/set_keyword.py). After every EPSR run, every file of the folder is added
to OUTPUT_ARCHIVE (a .tar.gz) under `after_step_NN/`, NN being the step number.
python/tests/test_cu_btc.py compares the Python wrapper with these files.
"""

import os
import resource
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

CASE_FOLDER = Path(__file__).resolve().parent
SET_KEYWORD_SCRIPT = CASE_FOLDER.parent / "set_keyword.py"
NOT_INPUTS = {"steps.txt", "README.md", "make_pdb.py", "make_fortran_reference.py"}


def unlimited_stack():
    _, hard_limit = resource.getrlimit(resource.RLIMIT_STACK)
    resource.setrlimit(resource.RLIMIT_STACK, (hard_limit, hard_limit))


def load_steps():
    steps = []
    for raw_line in (CASE_FOLDER / "steps.txt").read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line:
            steps.append(line)
    return steps


def run_step(step, folder, fortran_build):
    kind, _, rest = step.partition(" ")
    if kind == "run":
        stdin_name = None
        if "<" in rest:
            rest, stdin_name = (part.strip() for part in rest.split("<", 1))
        program, *args = rest.split()
        stdin_path = CASE_FOLDER / stdin_name if stdin_name else Path(os.devnull)
        with open(stdin_path, "rb") as stdin:
            subprocess.run(
                [str(fortran_build / program), "./", program, *args],
                cwd=folder,
                stdin=stdin,
                capture_output=True,
                env={**os.environ, "OMP_NUM_THREADS": "1"},
                preexec_fn=unlimited_stack,
                check=False,
            )
        return program == "epsr"
    if kind == "set":
        file_name, keyword, *values = rest.split()
        subprocess.run(
            [sys.executable, str(SET_KEYWORD_SCRIPT), str(folder / file_name), keyword, *values],
            check=True,
        )
        return False
    if kind == "delete":
        for file_name in rest.split():
            (folder / file_name).unlink(missing_ok=True)
        return False
    raise SystemExit(f"unknown step: {step}")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    fortran_build, archive_path = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory() as temporary, tarfile.open(archive_path, "w:gz") as archive:
        folder = Path(temporary)
        for input_file in CASE_FOLDER.iterdir():
            is_input = (
                input_file.is_file()
                and input_file.name not in NOT_INPUTS
                and not input_file.name.endswith(".tar.gz")  # e.g. the archive being written
            )
            if is_input:
                shutil.copy2(input_file, folder / input_file.name)
        for step_number, step in enumerate(load_steps(), start=1):
            if run_step(step, folder, fortran_build):
                snapshot = f"after_step_{step_number:02d}"
                for output_file in sorted(folder.iterdir()):
                    archive.add(output_file, arcname=f"{snapshot}/{output_file.name}")
                print(f"step {step_number:2d}  {step}  -> {snapshot}/")
    print(f"wrote {archive_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
