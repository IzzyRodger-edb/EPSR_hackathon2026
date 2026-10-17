"""Helpers for the experimental-data tests: build a case from real data and run its steps.

The cases are those of ``validation/``: water and amorphous silica, built from
scratch with ``mixato`` (``test_experimental_data.py``), and the Cu-BTC crystal,
whose unit cell is read from a PDB file with ``readpdb`` (``test_cu_btc.py``). Each
case folder holds the setup inputs and ``steps.txt``, the calculation step by step
(see ``validation/run_validation.sh``). The neutron data themselves are taken from
``experimental_data/<case>/data/``.

A step list can be run in two ways:

* through the Python wrapper (:func:`run_steps_with_wrapper`): setup programs with
  ``epsr.tools.run``, EPSR runs with ``epsr.Simulation``. A ``set`` step that changes
  a keyword of the ``.EPSR.inp`` or ``.pcof`` file becomes an override of the next
  Simulation, which is how a Python user changes parameters;
* with an external program (:func:`run_step_externally`): the ``epsr-rs`` CLI or the
  original Fortran build, with ``set`` steps applied by ``validation/set_keyword.py``,
  exactly as ``validation/run_validation.sh`` does.
"""

from __future__ import annotations

import os
import resource
import shutil
import subprocess
import sys
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

import epsr

from epsr_cases import REPOSITORY, copy_files

EXPERIMENTAL_DATA = REPOSITORY / "experimental_data"
VALIDATION = REPOSITORY / "validation"
SET_KEYWORD_SCRIPT = VALIDATION / "set_keyword.py"

#: Files of a validation case folder that are not inputs of the calculation.
CASE_FILES_TO_SKIP = {"steps.txt", "README.md"}
#: Suffixes of case files that are not inputs: scripts, reference archives, figures,
#: and the data files, which are taken from experimental_data/ instead.
CASE_SUFFIXES_TO_SKIP = {".py", ".gz", ".png", ".mint01"}
#: Files whose keywords are EPSR parameters, which the wrapper sets as overrides.
EPSR_KEYWORD_FILE_SUFFIXES = (".EPSR.inp", ".pcof")


@dataclass(frozen=True)
class ExperimentalCase:
    """One validation system fitted to real neutron data."""

    name: str  # folder name in validation/ and experimental_data/
    input_file: str  # the .EPSR.inp file
    requested_ep_amplitude: float  # ereq set before the empirical-potential refinement
    box_side: float  # Å, side of the cubic box
    n_atoms: int
    n_molecules: int

    @property
    def case_folder(self) -> Path:
        return VALIDATION / self.name

    @property
    def data_folder(self) -> Path:
        return EXPERIMENTAL_DATA / self.name / "data"

    @property
    def data_files(self) -> list[Path]:
        return sorted(self.data_folder.glob("*.mint01"))

    @property
    def available(self) -> bool:
        return (self.case_folder / "steps.txt").exists() and bool(self.data_files)


WATER = ExperimentalCase(
    name="water",
    input_file="waterbox.EPSR.inp",
    requested_ep_amplitude=3.0,
    box_side=31.0723,
    n_atoms=3000,
    n_molecules=1000,
)
SILICA = ExperimentalCase(
    name="silica_bulk",
    input_file="silica.EPSR.inp",
    requested_ep_amplitude=50.0,
    box_side=40.7707,
    n_atoms=4500,
    n_molecules=4500,  # single-atom "molecules"
)
CU_BTC = ExperimentalCase(
    name="cu-btc",
    input_file="CuBTC.EPSR.inp",
    requested_ep_amplitude=10.0,
    box_side=26.3336,  # one unit cell of the crystal, from the PDB (not built by mixato)
    n_atoms=624,
    n_molecules=56,  # 24 Cu2 paddlewheels + 32 BTC linkers
)
#: The cases whose box is built with mixato (test_experimental_data.py); Cu-BTC is not.
MIXATO_CASES = [WATER, SILICA]


@dataclass(frozen=True)
class Step:
    """One line of ``steps.txt``."""

    kind: str  # "run", "set" or "delete"
    #: run: program and arguments; set: file, keyword, values; delete: files
    words: tuple[str, ...]
    stdin_file: str | None = None  # run: file with the answers, in the case folder

    @property
    def program(self) -> str:
        return self.words[0]

    @property
    def args(self) -> list[str]:
        return list(self.words[1:])

    @property
    def is_epsr_run(self) -> bool:
        return self.kind == "run" and self.program == "epsr"

    def __str__(self) -> str:
        text = f"{self.kind} {' '.join(self.words)}"
        return text + (f" < {self.stdin_file}" if self.stdin_file else "")


def load_steps(case: ExperimentalCase) -> list[Step]:
    """Parse ``steps.txt`` the way ``validation/run_validation.sh`` does."""
    steps = []
    for raw_line in (case.case_folder / "steps.txt").read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        kind, _, rest = line.partition(" ")
        stdin_file = None
        if kind == "run" and "<" in rest:
            rest, stdin_file = (part.strip() for part in rest.split("<", 1))
        steps.append(Step(kind, tuple(rest.split()), stdin_file))
    return steps


def shorten_schedule(
    steps: list[Step], equilibration_iterations: int, refinement_iterations: int
) -> list[Step]:
    """The same steps with fewer EPSR iterations per run.

    ``steps.txt`` sets ``niter 99`` before the equilibration run and ``niter 100``
    before the refinement runs; those values are replaced.
    """
    replacements = {
        "99": str(equilibration_iterations),
        "100": str(refinement_iterations),
    }
    shortened = []
    for step in steps:
        if step.kind == "set" and step.words[1] == "niter":
            new_value = replacements.get(step.words[2], step.words[2])
            step = replace(step, words=(step.words[0], "niter", new_value))
        shortened.append(step)
    return shortened


def prepare_case_folder(case: ExperimentalCase, destination: Path) -> Path:
    """Copy the case's setup inputs and the experimental data into a new folder."""
    destination.mkdir(parents=True)  # fails if the folder exists: inputs must be fresh
    copy_files(
        case.case_folder, destination, CASE_FILES_TO_SKIP | CASE_SUFFIXES_TO_SKIP
    )
    for data_file in case.data_files:
        shutil.copy2(data_file, destination / data_file.name)
    return destination


# --------------------------------------------------------------------- the wrapper
@dataclass
class WrapperStepResult:
    """What the wrapper produced for one step."""

    step: Step
    simulation: epsr.Simulation | None = None  # set for EPSR runs
    summary: epsr.Summary | None = None


@dataclass
class WrapperRun:
    """A case run through the wrapper: its folder and the result of every EPSR run."""

    case: ExperimentalCase
    folder: Path
    epsr_results: list[WrapperStepResult] = field(default_factory=list)

    @property
    def final(self) -> epsr.Simulation:
        """The Simulation of the last EPSR run."""
        simulation = self.epsr_results[-1].simulation
        assert simulation is not None
        return simulation

    def summary_after(self, run_index: int) -> epsr.Summary:
        """The summary of EPSR run ``run_index`` (0-based; negative counts from the end)."""
        summary = self.epsr_results[run_index].summary
        assert summary is not None
        return summary


def run_steps_with_wrapper(
    steps: list[Step], folder: Path, case_folder: Path
) -> Iterator[WrapperStepResult]:
    """Run the steps through the Python wrapper, yielding a result after each step.

    ``set`` steps on EPSR keyword files are collected and passed as overrides to the
    next EPSR run; EPSR saves them back to its files at the end of that run.
    """
    pending_overrides: list[str] = []
    for step in steps:
        if step.kind == "run" and step.is_epsr_run:
            simulation = epsr.Simulation(
                folder, step.args[0], overrides=pending_overrides
            )
            summary = simulation.run()
            pending_overrides = []
            yield WrapperStepResult(step, simulation, summary)
            continue

        if step.kind == "run":
            answers = (
                (case_folder / step.stdin_file).read_text() if step.stdin_file else ""
            )
            epsr.tools.run(folder, step.program, step.args, answers)
        elif step.kind == "set":
            file_name, keyword, *values = step.words
            if not file_name.endswith(EPSR_KEYWORD_FILE_SUFFIXES):
                raise ValueError(f"cannot turn {step} into an override")
            pending_overrides.append(f"{keyword} {' '.join(values)}")
        elif step.kind == "delete":
            for file_name in step.words:
                (folder / file_name).unlink(missing_ok=True)
        else:
            raise ValueError(f"unknown step {step}")
        yield WrapperStepResult(step)


# --------------------------------------------------------------------- external programs
def _unlimited_stack() -> None:
    """The Fortran programs need a large stack (``ulimit -s unlimited``)."""
    soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_STACK)
    resource.setrlimit(resource.RLIMIT_STACK, (hard_limit, hard_limit))


def cli_command(cli_executable: Path) -> Callable[[Step], list[str]]:
    """Command line of a step for the ``epsr-rs`` CLI."""

    def command(step: Step) -> list[str]:
        return [str(cli_executable), "./", step.program, *step.args, "--quiet"]

    return command


def fortran_command(fortran_build: Path) -> Callable[[Step], list[str]]:
    """Command line of a step for the Fortran programs built by ``build_fortran_reference.sh``."""

    def command(step: Step) -> list[str]:
        return [str(fortran_build / step.program), "./", step.program, *step.args]

    return command


def run_step_externally(
    step: Step,
    folder: Path,
    case_folder: Path,
    command_for_step: Callable[[Step], list[str]],
) -> None:
    """Run one step with an external program, as ``validation/run_validation.sh`` does."""
    if step.kind == "run":
        stdin_path = (
            case_folder / step.stdin_file if step.stdin_file else Path(os.devnull)
        )
        environment = {**os.environ, "OMP_NUM_THREADS": "1"}
        with open(stdin_path, "rb") as stdin:
            subprocess.run(
                command_for_step(step),
                cwd=folder,
                stdin=stdin,
                capture_output=True,
                env=environment,
                preexec_fn=_unlimited_stack,
                check=False,  # like run_validation.sh: a failing program shows up in the comparison
            )
    elif step.kind == "set":
        file_name, keyword, *values = step.words
        subprocess.run(
            [
                sys.executable,
                str(SET_KEYWORD_SCRIPT),
                str(folder / file_name),
                keyword,
                *values,
            ],
            check=True,
            capture_output=True,
        )
    elif step.kind == "delete":
        for file_name in step.words:
            (folder / file_name).unlink(missing_ok=True)
    else:
        raise ValueError(f"unknown step {step}")


# --------------------------------------------------------------------- reading results
def read_weights(path: Path) -> dict[tuple[str, str], float]:
    """Intermolecular weight of each atom-type pair from a ``.NWTStot.wts`` file.

    The file has a header line, then for each pair a line with the two labels and a
    line with the weights (the first is the intermolecular one).
    """
    lines = path.read_text().splitlines()[1:]
    weights = {}
    for label_line, value_line in zip(lines[0::2], lines[1::2], strict=True):
        first_type, second_type = label_line.split()
        weights[(first_type, second_type)] = float(value_line.split()[0])
    return weights


#: Marker lines of a .bragg file; "refections" is the engine's own spelling (sic).
BRAGG_COUNT_MARKER = "# Number of Bragg refections calculated:"
BRAGG_UNIQUE_VALUES_MARKER = "Unique Bragg Q values"


def read_bragg_reflections(bragg_file: Path) -> tuple[int, list[float]]:
    """Number of reflections and the unique Bragg Q values listed in a ``.bragg`` file.

    Every line of the file starts with ``#``. After the lattice header come the
    line ``# Number of Bragg refections calculated: N``, the heading
    ``# Unique Bragg Q values and multiplicities`` and then one line
    ``#  index  Q  multiplicity`` per unique Q value, in the order the engine
    enumerated the reflections (not sorted; the (0 0 0) reflection sits mid-list).
    """
    n_reflections, unique_q_values = 0, []
    listing_unique_values = False
    for line in bragg_file.read_text().splitlines():
        words = line.lstrip("# ").split()
        if line.startswith(BRAGG_COUNT_MARKER):
            n_reflections = int(words[-1])
        elif BRAGG_UNIQUE_VALUES_MARKER in line:
            listing_unique_values = True
        elif listing_unique_values and len(words) == 3:
            unique_q_values.append(float(words[1]))
    return n_reflections, unique_q_values


def first_peak_position(r: object, gr: object, r_min: float, r_max: float) -> float:
    """Position of the highest point of ``gr`` between ``r_min`` and ``r_max``."""
    r_values, gr_values = np.asarray(r), np.asarray(gr)
    window = (r_values > r_min) & (r_values < r_max)
    return float(r_values[window][np.argmax(gr_values[window])])


def coordination_number(
    r: object, gr: object, partner_density: float, cutoff: float
) -> float:
    """Average number of partner atoms within ``cutoff``: 4 pi rho integral r^2 g(r) dr.

    The integral uses the trapezoidal rule because the EPSR r grid is not uniform
    (its first interval is a quarter of the others).
    """
    r_values = np.asarray(r, dtype=float)
    gr_values = np.asarray(gr, dtype=float)
    inside = r_values <= cutoff
    r_inside = r_values[inside]
    integrand = r_inside**2 * gr_values[inside]
    integral = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r_inside))
    return float(4 * np.pi * partner_density * integral)
