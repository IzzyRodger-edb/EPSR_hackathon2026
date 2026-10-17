"""The standalone Python script: the "Python input file" a researcher can check and run.

``render_script`` writes ``run_<name>.py``. The script is self-contained: it embeds
the generated input texts verbatim, imports only ``epsr`` (the wrapper) and
``epsrvis`` (the toolkit) plus the standard library, and reproduces the calculation
step by step with the wrapper API. Mellor itself runs this script
(:mod:`mellor.runner`), so the researcher re-runs exactly what produced the results.

Script contract (``python run_<name>.py --out DIR [--quick] [--no-record]``):

1. creates ``DIR/epsr`` (fails if it exists and is not empty), copies the data files
   from ``DIR/data`` (or the paths embedded in the script) and writes the input files;
2. setup: ``makemole`` per template (molecular/atomic), ``mixato`` with the answers,
   ``randomise`` (if scheduled), deletes the ``.mol`` files, ``fmole <box> <n>``
   (molecular only), ``nwts`` per data set — or ``readpdb`` + ``changeato`` + ``nwts``
   for a pdb system; every step via ``epsr.tools.run``;
3. the EPSR runs of the schedule: a fresh ``epsr.Simulation`` per run with the
   overrides (``niter``, ``ereq``), driven iteration by iteration through
   ``epsrvis.TrajectoryRecorder``/``record_run`` with an adapter that maps the wrapper
   names (``control.niter``, ``potential.empirical()``...) to the recorder's duck-typed
   interface; ``--no-record`` uses ``run()`` instead;
4. writes ``DIR/trajectory`` and ``DIR/manifest.json`` (see docs/architecture.md §6),
   prints a table of the run summaries and exits 0.

The script has a configuration block at the top (name, schedule, overrides) so the
researcher can see and change the numbers without reading the rest.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from . import __version__
from .inputs import EpsrInputs
from .spec import SystemSpec

__all__ = [
    "render_script",
    "write_script",
    "script_name",
    "epsrvis_case",
    "EPSRVIS_CASES",
    "SCRIPT_BODY",
]

#: Written into every script (``MELLOR_VERSION``) and the manifest.
MELLOR_VERSION = __version__

#: The systems ``epsrvis.systems`` knows (case name -> its atom-type labels). Other
#: systems get the toolkit's generic defaults under the run name.
EPSRVIS_CASES: dict[str, frozenset[str]] = {
    "water": frozenset({"OW", "HW"}),
    "silica_bulk": frozenset({"Si", "O"}),
}


def script_name(spec: SystemSpec) -> str:
    return f"run_{spec.name}.py"


def epsrvis_case(spec: SystemSpec) -> str:
    """The ``epsrvis`` case name of the system: ``water`` or ``silica_bulk`` when the
    atom types are exactly those systems' (their colours, bonds and structural
    analyses apply), otherwise the run name (generic defaults)."""
    labels = frozenset(spec.labels)
    for case, case_labels in EPSRVIS_CASES.items():
        if labels == case_labels:
            return case
    return spec.name


def render_script(spec: SystemSpec, inputs: EpsrInputs) -> str:
    """The script text. ``inputs`` must come from ``inputs.generate(spec)``."""
    if inputs.spec is not spec and inputs.spec.name != spec.name:
        raise ValueError(
            f"inputs were generated for {inputs.spec.name!r}, not for {spec.name!r}"
        )
    if not inputs.inp_text or not inputs.inp_name:
        raise ValueError("inputs have no .EPSR.inp text: use inputs.generate(spec)")
    return _header(spec, inputs) + _configuration(spec, inputs) + SCRIPT_BODY


def write_script(spec: SystemSpec, inputs: EpsrInputs, folder: str | Path) -> Path:
    """Write ``run_<name>.py`` into ``folder`` and return its path."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / script_name(spec)
    path.write_text(render_script(spec, inputs), encoding="utf-8")
    return path


# --------------------------------------------------------------------------- rendering
def _docstring_safe(text: str) -> str:
    """Text that can sit inside a triple-double-quoted string."""
    return text.replace("\\", "/").replace('"""', "'''")


def _text_literal(text: str) -> str:
    """A Python literal of a file text: a verbatim triple-quoted block when the text
    allows it (no backslashes, no ``'''``, newline-terminated), ``repr`` otherwise."""
    if "\\" in text or "'''" in text or not text.endswith("\n"):
        return repr(text)
    return "'''\\\n" + text + "'''"


def _rows_literal(rows: list[tuple[Any, ...]]) -> str:
    """A list of tuples, one tuple per line."""
    if not rows:
        return "[]"
    body = "".join(
        "    (" + ", ".join(repr(item) for item in row) + "),\n" for row in rows
    )
    return "[\n" + body + "]"


def _dict_literal(items: dict[str, str]) -> str:
    """A dict of texts (keys sorted as given), one entry per line or block."""
    if not items:
        return "{}"
    body = "".join(
        f"    {key!r}: {_text_literal(value)},\n" for key, value in items.items()
    )
    return "{\n" + body + "}"


def _plain_dict_literal(items: dict[str, Any]) -> str:
    if not items:
        return "{}"
    body = "".join(f"    {key!r}: {value!r},\n" for key, value in items.items())
    return "{\n" + body + "}"


def _header(spec: SystemSpec, inputs: EpsrInputs) -> str:
    schedule = spec.schedule.iterations()
    runs = ", ".join(
        f"{stage} ({niter} it., ereq {ereq:g})" for stage, niter, ereq in schedule
    )
    steps = "; ".join(inputs.setup_steps)
    lines = [
        "#!/usr/bin/env python3",
        f'"""{spec.name}: EPSR refinement of {_docstring_safe(spec.title) or spec.name}.',
        "",
        f"Written by Mellor {MELLOR_VERSION} on {datetime.now():%Y-%m-%d %H:%M} from the system",
        "specification (spec.json). This is the Python input file of the calculation: it",
        "contains every number and every EPSR input file, and runs the refinement with the",
        "``epsr`` wrapper of the Rust EPSR engine and the ``epsrvis`` toolkit, nothing else.",
        "",
        "    python "
        + script_name(spec)
        + " --out DIR [--quick] [--no-record] [--verbose]",
        "",
        "1. DIR/epsr: the input files embedded below are written and the measured data are",
        "   copied from DIR/data (or from the original locations listed in DATA_SOURCES);",
        f"2. setup with the EPSR programs: {_docstring_safe(steps)};",
        f"3. the EPSR runs of the schedule, each a fresh epsr.Simulation: {runs};",
        "   every iteration is recorded into DIR/trajectory (epsrvis folder format);",
        "4. DIR/manifest.json: the summary of every run (R-factor, energy, pressure, EP",
        "   amplitude, timings, versions) and a table on the console.",
        "",
        "Change the numbers in the configuration block to run a variant; the EPSR input",
        "files are the same text EPSR reads (list-directed, so spacing is free).",
    ]
    if spec.notes:
        lines += [
            "",
            "Notes from the specification:",
            _docstring_safe(spec.notes).rstrip(),
        ]
    lines += ['"""', ""]
    return "\n".join(lines) + "\n"


def _configuration(spec: SystemSpec, inputs: EpsrInputs) -> str:
    schedule = spec.schedule
    datasets = [
        (dataset.label, dataset.file_name, dataset.weights_file)
        for dataset in spec.datasets
    ]
    data_sources = {path.name: str(path) for path in inputs.data_files}
    molecules = [(molecule.name, molecule.count) for molecule in spec.molecules]
    pdb_file = inputs.pdb_name if spec.kind == "pdb" else None
    pdb_source = str(inputs.pdb_file) if inputs.pdb_file is not None else None
    lines = [
        "from __future__ import annotations",
        "",
        "import argparse",
        "import json",
        "import platform",
        "import shutil",
        "import sys",
        "import time",
        "from datetime import datetime",
        "from pathlib import Path",
        "",
        "import epsr",
        "import epsrvis",
        "",
        "# " + "=" * 76 + " #",
        "# Configuration: everything a researcher may want to check or change.          #",
        "# " + "=" * 76 + " #",
        f"MELLOR_VERSION = {MELLOR_VERSION!r}",
        f"NAME = {spec.name!r}  # run name: <NAME>.ato, <NAME>.pcof and <NAME>.EPSR.* files",
        f"TITLE = {spec.title!r}",
        f"KIND = {spec.kind!r}  # molecular | atomic | pdb",
        f"INP_NAME = {inputs.inp_name!r}",
        f"EPSRVIS_CASE = {epsrvis_case(spec)!r}  # system settings of epsrvis (figures)",
        f"TEMPERATURE_K = {float(spec.temperature)!r}",
        f"DENSITY = {float(spec.density)!r}  # atoms/A^3 (mixato); the box side follows",
        f"ATOM_TYPES = {spec.labels!r}",
        "",
        "# Molecule templates in mixato order: (template name, number of molecules).",
        f"MOLECULES = {_rows_literal(molecules)}",
        f"MIXATO_ANSWERS = {inputs.mixato_answers!r}",
        f"RANDOMISE = {bool(schedule.randomise)!r}  # randomise molecule positions and orientations",
        f"FMOLE_ITERATIONS = {int(schedule.fmole_iterations)!r}  # geometry relaxation shakes (0 = none)",
        "",
        "# pdb systems: the unit cell is read with readpdb (crystal positions kept) and the",
        "# force field applied with changeato; the PDB file is copied as PDB_FILE.",
        f"PDB_FILE = {pdb_file!r}",
        f"PDB_SOURCE = {pdb_source!r}",
        f"READPDB_ANSWERS = {inputs.readpdb_answers!r}",
        f"CHANGEATO_COMMANDS = {inputs.changeato_commands!r}",
        "",
        "# Data sets in .EPSR.inp order: (label, data file, weights file written by nwts).",
        f"DATASETS = {_rows_literal(datasets)}",
        "# Where the data files came from (used when DIR/data holds no copy).",
        f"DATA_SOURCES = {_plain_dict_literal(data_sources)}",
        "",
        "# EPSR runs: (stage, niter, ereq [kJ/mol]). Each run is one epsr invocation on the",
        "# files the previous run left. The scratch run creates the .pcof (ereq 0); ereq is",
        "# set as an override before the first EP refinement run and niter whenever it",
        "# changes (EPSR saves overrides back into its files).",
        f"SCHEDULE = {_rows_literal(schedule.iterations())}",
        "# The schedule of --quick (previews, tests).",
        f"QUICK_SCHEDULE = {_rows_literal(schedule.iterations(quick=True))}",
        "",
        "# The setup programs in order (the steps.txt syntax of validation/).",
        f"SETUP_STEPS = {inputs.setup_steps!r}",
        "",
        "# The EPSR input files, verbatim (.mol templates, .NWTS.dat weights, .EPSR.inp).",
        f"INPUT_FILES = {_dict_literal(inputs.files())}",
        "# The standard input of the interactive programs, also written as *.stdin files.",
        f"ANSWER_FILES = {_dict_literal(inputs.answer_files())}",
        "",
        "# Folder layout under --out.",
        'EPSR_FOLDER = "epsr"',
        'DATA_FOLDER = "data"',
        'TRAJECTORY_FOLDER = "trajectory"',
        'MANIFEST_FILE = "manifest.json"',
    ]
    return "\n".join(lines) + "\n"


#: The fixed part of every script: the recorder adapter, setup, runs, manifest.
SCRIPT_BODY = '''
# ============================================================================ #
# The calculation. Nothing below needs editing.                                #
# ============================================================================ #


def log(message: str = "") -> None:
    print(message, flush=True)


class RecorderAdapter:
    """``epsr.Simulation`` seen through the duck-typed interface of ``epsrvis.record_run``.

    The recorder looks every quantity up by name, as a property or a method (see
    ``epsrvis.recorder.SimulationLike``); this class maps the wrapper's names onto
    them. The wrapper writes the output files itself on the last ``iterate`` call,
    so no ``write_fort35`` is exposed and the folder ends up as after ``run()``.
    """

    def __init__(self, simulation: epsr.Simulation) -> None:
        self.simulation = simulation

    # ---- the run
    def iterate(self, iteration: int) -> epsr.Summary:
        return self.simulation.iterate(iteration)

    @property
    def niter(self) -> int:
        return self.simulation.control.niter

    # ---- grids, sizes and labels
    @property
    def q(self):
        return self.simulation.q

    @property
    def r(self):
        return self.simulation.r

    @property
    def n_datasets(self) -> int:
        return self.simulation.n_datasets

    @property
    def pair_labels(self) -> list[tuple[str, str]]:
        return self.simulation.pair_labels

    @property
    def atom_types(self) -> list[str]:
        return self.simulation.atom_types

    @property
    def atom_labels(self) -> list[str]:
        return self.simulation.atom_labels

    @property
    def molecule_atoms(self) -> list[range]:
        return self.simulation.molecule_atoms

    def get_parameter(self, name: str, dataset: int | None = None) -> str | None:
        return self.simulation.get_parameter(name, dataset)

    # ---- data and fits, per data set (0-based)
    def data_fq(self, dataset: int):
        return self.simulation.data_fq(dataset)

    def data_errors(self, dataset: int):
        return self.simulation.data_errors(dataset)

    def simulated_fq(self, dataset: int):
        return self.simulation.simulated_fq(dataset)

    def data_q_range(self, dataset: int) -> tuple[int, int]:
        return self.simulation.data_q_range(dataset)

    # ---- partials and potentials, all pairs at once: (n_pairs, n)
    def gr_inter(self):
        return self.simulation.gr_inter()

    def gr_intra(self):
        return self.simulation.gr_intra()

    def sq_inter(self):
        return self.simulation.sq_inter()

    def sq_intra(self):
        return self.simulation.sq_intra()

    def empirical_potential(self):
        return self.simulation.potential.empirical()

    def reference_potential(self):
        return self.simulation.potential.reference()

    # ---- configuration
    @property
    def positions(self):
        return self.simulation.positions

    @property
    def box_size(self) -> float:
        return self.simulation.box_size

    # ---- convergence numbers
    @property
    def summary(self) -> dict:
        """The ``Summary`` fields under the wrapper's names and under the Fortran names
        (``foldsv``, ``sentb``, ``epcup``, ``epccoeff``...), so that every alias the
        recorder knows resolves."""
        summary = self.simulation.summary
        values = dict(summary.fortran)
        values.update(summary.as_dict())
        return values


def summary_dict(summary: epsr.Summary) -> dict:
    """The descriptive fields of a ``Summary`` as plain JSON values."""
    values = summary.as_dict()
    values["r_factor_per_dataset"] = [float(v) for v in summary.r_factor_per_dataset]
    return values


def prepare_folder(out: Path) -> Path:
    """DIR/epsr with the input files, the answer files and copies of the data."""
    folder = out / EPSR_FOLDER
    if folder.exists() and any(folder.iterdir()):
        raise SystemExit(
            f"error: {folder} exists and is not empty (choose another --out or remove it)"
        )
    folder.mkdir(parents=True, exist_ok=True)
    for name, text in INPUT_FILES.items():
        (folder / name).write_text(text, encoding="latin-1")
    for name, text in ANSWER_FILES.items():
        (folder / name).write_text(text, encoding="latin-1")
    copies = [(data_file, DATA_SOURCES.get(data_file)) for _, data_file, _ in DATASETS]
    if PDB_FILE:
        copies.append((PDB_FILE, PDB_SOURCE))
    for name, source in copies:
        candidates = [out / DATA_FOLDER / name]
        if source:
            candidates += [out / DATA_FOLDER / Path(source).name, Path(source)]
        found = next((path for path in candidates if path.is_file()), None)
        if found is None:
            raise SystemExit(
                f"error: {name} not found in {out / DATA_FOLDER}"
                + (f" nor at {source}" if source else "")
            )
        shutil.copy2(found, folder / name)
    return folder


def run_setup(folder: Path, verbose: bool) -> None:
    """The EPSR setup programs (SETUP_STEPS) through ``epsr.tools.run``."""

    def tool(command: str, args: list[str], answers: list[str] | None = None) -> None:
        log(f"  {command} {' '.join(args)}".rstrip())
        try:
            epsr.tools.run(folder, command, args, answers or [], verbose=verbose)
        except epsr.EpsrError as error:
            raise SystemExit(f"error: {command} {' '.join(args)} failed: {error}") from error

    def expect(name: str, program: str) -> None:
        if not (folder / name).is_file():
            raise SystemExit(f"error: {program} did not write {name} in {folder}")

    if KIND == "pdb":
        tool("readpdb", [PDB_FILE], READPDB_ANSWERS)
        expect(f"{NAME}.ato", f"readpdb {PDB_FILE}")
        tool("changeato", [f"{NAME}.ato"], CHANGEATO_COMMANDS)
    else:
        for template, _count in MOLECULES:
            tool("makemole", [template])
            expect(f"{template}.ato", f"makemole {template}")
        tool("mixato", [], MIXATO_ANSWERS)
        expect(f"{NAME}.ato", "mixato")
        if RANDOMISE:
            tool("randomise", [NAME])
        for template, _count in MOLECULES:
            # fmole regenerates a template whose .mol is present, which the Rust port
            # does not support (validation/README.md); the template .ato is what counts.
            (folder / f"{template}.mol").unlink(missing_ok=True)
        if KIND == "molecular" and FMOLE_ITERATIONS > 0:
            tool("fmole", [NAME, str(FMOLE_ITERATIONS)])
    for label, _data_file, weights_file in DATASETS:
        tool("nwts", [f"{label}.NWTS.dat"])
        expect(weights_file, f"nwts {label}.NWTS.dat")


def system_facts(simulation: epsr.Simulation) -> dict:
    """Sizes and labels as the engine read them."""
    return {
        "box_side": float(simulation.box_size),
        "n_atoms": int(simulation.n_atoms),
        "n_molecules": int(simulation.n_molecules),
        "density": float(simulation.density),
        "atom_types": list(simulation.atom_types),
        "pairs": [f"{first}-{second}" for first, second in simulation.pair_labels],
    }


def write_manifest(out: Path, run_info: dict, system: dict, runs: list[dict]) -> Path:
    manifest = {
        "name": NAME,
        "title": TITLE,
        "kind": KIND,
        **run_info,
        "python": platform.python_version(),
        "versions": {
            "epsr": epsr.__version__,
            "epsrvis": epsrvis.__version__,
            "mellor": MELLOR_VERSION,
        },
        "script": Path(__file__).name,
        "epsr_folder": EPSR_FOLDER,
        "trajectory": TRAJECTORY_FOLDER if run_info["record"] else None,
        "inp_name": INP_NAME,
        "epsrvis_case": EPSRVIS_CASE,
        "temperature": TEMPERATURE_K,
        **system,
        "datasets": [
            {"label": label, "file": data_file, "wts": weights_file}
            for label, data_file, weights_file in DATASETS
        ],
        "setup_steps": SETUP_STEPS,
        "schedule": [
            {"stage": stage, "niter": niter, "ereq": ereq}
            for stage, niter, ereq in run_info["schedule_runs"]
        ],
        "runs": runs,
        "final": runs[-1]["summary"] if runs else {},
    }
    manifest.pop("schedule_runs")
    path = out / MANIFEST_FILE
    path.write_text(json.dumps(manifest, indent=2) + "\\n")
    return path


def run_schedule(
    out: Path, folder: Path, schedule: list, record: bool, verbose: bool, run_info: dict
) -> list[dict]:
    """The EPSR runs, one fresh Simulation each; manifest and trajectory after every run."""
    recorder = epsrvis.TrajectoryRecorder(case=EPSRVIS_CASE) if record else None
    runs: list[dict] = []
    system: dict = {}
    # the scratch-start .EPSR.inp has niter 1 and there is no .pcof yet (ereq 0)
    previous_niter, previous_ereq = 1, 0.0
    for number, (stage, niter, ereq) in enumerate(schedule, start=1):
        overrides = []
        if ereq != previous_ereq:
            overrides.append(f"ereq {ereq!r}")
        if niter != previous_niter:
            overrides.append(f"niter {niter}")
        log(
            f"run {number}/{len(schedule)}: {stage}, {niter} iteration(s), ereq {ereq:g}"
            + (f", overrides {overrides}" if overrides else "")
        )
        started = time.perf_counter()
        simulation = epsr.Simulation(folder, INP_NAME, overrides=overrides, verbose=verbose)
        if not system:
            system = system_facts(simulation)
        if recorder is not None:
            epsrvis.record_run(RecorderAdapter(simulation), recorder, ereq=ereq)
        else:
            simulation.run()
        duration = time.perf_counter() - started
        summary = simulation.summary
        runs.append(
            {
                "run": number,
                "stage": stage,
                "niter": niter,
                "ereq": ereq,
                "overrides": overrides,
                "duration_s": round(duration, 3),
                "summary": summary_dict(summary),
            }
        )
        previous_niter, previous_ereq = niter, ereq
        log(
            f"  R-factor {summary.r_factor:.4e}, energy {summary.energy:.3f} kJ/mol,"
            f" pressure {summary.pressure:.3f} kbar, EP {summary.ep_amplitude:.3f} kJ/mol,"
            f" {duration:.1f} s"
        )
        if recorder is not None:
            recorder.save(out / TRAJECTORY_FOLDER)
        run_info["complete"] = number == len(schedule)
        run_info["duration_s"] = round(time.perf_counter() - run_info["t0"], 3)
        write_manifest(out, {k: v for k, v in run_info.items() if k != "t0"}, system, runs)
    return runs


def table(runs: list[dict]) -> str:
    """A plain-text table of the run summaries."""
    lines = [
        f"{'run':>3}  {'stage':<20} {'niter':>5} {'ereq':>6} {'R-factor':>10}"
        f" {'energy':>9} {'pressure':>9} {'EP ampl':>8} {'time':>8}"
    ]
    for run in runs:
        summary = run["summary"]
        lines.append(
            f"{run['run']:>3}  {run['stage']:<20} {run['niter']:>5} {run['ereq']:>6g}"
            f" {summary['r_factor']:>10.3e} {summary['energy']:>9.4g}"
            f" {summary['pressure']:>9.4g} {summary['ep_amplitude']:>8.3f}"
            f" {run['duration_s']:>7.1f}s"
        )
    return "\\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=f"EPSR refinement {NAME}: {TITLE}")
    parser.add_argument(
        "--out",
        default=".",
        help="run folder: DIR/epsr, DIR/trajectory and DIR/manifest.json are written "
        "(default: the current folder)",
    )
    parser.add_argument(
        "--quick", action="store_true", help="run QUICK_SCHEDULE instead of SCHEDULE"
    )
    parser.add_argument(
        "--no-record",
        action="store_true",
        help="run each EPSR run with Simulation.run() and write no trajectory",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="print the engine's own progress messages"
    )
    args = parser.parse_args(argv)
    out = Path(args.out).resolve()
    schedule = QUICK_SCHEDULE if args.quick else SCHEDULE
    record = not args.no_record
    run_info = {
        "quick": bool(args.quick),
        "record": record,
        "created": datetime.now().isoformat(timespec="seconds"),
        "schedule_runs": schedule,
        "complete": False,
        "duration_s": 0.0,
        "t0": time.perf_counter(),
    }

    log(f"{NAME}: {TITLE}")
    log(
        f"epsr {epsr.__version__}, epsrvis {epsrvis.__version__}, "
        f"mellor {MELLOR_VERSION}, python {platform.python_version()}"
    )
    log(f"writing {out / EPSR_FOLDER}")
    folder = prepare_folder(out)
    log("setup:")
    run_setup(folder, args.verbose)
    total = sum(niter for _, niter, _ in schedule)
    log(
        f"EPSR: {len(schedule)} runs, {total} iterations"
        + (" (quick schedule)" if args.quick else "")
        + ("" if record else " (not recorded)")
    )
    runs = run_schedule(out, folder, schedule, record, args.verbose, run_info)
    log()
    log(table(runs))
    log()
    log(f"done in {time.perf_counter() - run_info['t0']:.1f} s: {out / MANIFEST_FILE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
'''
