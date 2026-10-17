# Testing

The suite in `python/tests/` checks the wrapper against the same references the
Rust engine is validated with, plus wrapper-specific behaviour.

```bash
cd python
uv run pytest -n auto                 # everything: 217 tests (2 opt-in), about 5.5 minutes on 4 cores
uv run pytest -m "not slow"           # skip case replays and experimental-data runs: about 1 minute
uv run pytest -m experimental          # only the refinements of real neutron data
uv run pytest tests/test_ex_cases.py -k "ex20 or ex50"    # selected verification cases
EPSR_STRICT=1 uv run pytest           # require byte-identical output everywhere
```

`-n auto` uses pytest-xdist (in the dev dependency group). The `epsr-rs` CLI is
built with cargo on first use for the cross-checks; if cargo is unavailable those
checks are skipped.

## Equivalence with the Rust tests

| Rust test | wrapper test |
|---|---|
| `rust/epsr-core/tests/test_case.rs::dmso_water_test_case_matches_fortran` | `test_reference.py::test_outputs_match_fortran_reference` (same 26 hashes, same masking) and `test_api_sanity` |
| `test_case.rs::overrides_are_applied` | `test_reference.py::test_overrides_are_applied` |
| `test_case.rs::epsr_is_send` | `test_threads.py` (runs on threads with the GIL released) |
| `rust/tools/run_all_ex_tests.sh` (115 `EPSR-test-ex*` cases against the Fortran) | `test_ex_cases.py::test_case_matches_fortran[...]` |

`test_ex_cases.py` replays each case's `commands.txt` through the wrapper:
`epsr` lines through `epsr.Simulation(...).run()`, every other program through
`epsr.tools.run`, with `repeat N` honoured. It then compares every file in the
case's `fortran_reference.tar.gz` with the wrapper's output, using the masking
rules of `rust/tools/compare_outputs.py`. Files the wrapper writes that the
Fortran did not also fail the test.

## Wrapper-specific tests

| file | what it checks |
|---|---|
| `test_reference.py` | wrapper output byte-identical to the `epsr-rs` CLI over two chained runs; `iterations()` writes the same files as `run()`; determinism |
| `test_simulation.py` | shapes, dtypes, copies, pair lookup, datasets, fits, configuration, overrides (dict and list), run-order rules, error mapping |
| `test_potential.py` | EP get/set round trip is exact; scaling is linear below `ereq` and capped above it; raising `ereq` lifts the cap; zero, perturb, flags; edits persist to `.pcof`; edits change the simulation |
| `test_project_io.py` | template copying and protection, snapshots, `files()`, keyword-file and table parsers agree with the arrays |
| `test_refine.py` | the run loop, `.erg` growth, callbacks, early stopping, overrides, snapshots, CSV and pandas export |
| `test_tools.py` | `nwts` rebuilds the shipped `.wts` files byte for byte, keyword tables, `write_keyword_file` round trip, interactive answers, the deliberate `makeato` crash maps to `EpsrUnsupportedError` |
| `test_threads.py` | GIL release, parallel runs equal sequential runs, `refine_many`, tools from several threads; byte-identical results with 1 and 4 engine (rayon) threads for a cubic and a triclinic box; concurrent cubic and triclinic runs sharing the rayon pool equal their sequential runs |
| `test_examples.py` | the five example scripts run (the long ones in their short mode) |
| `test_plot.py` | the plotting helpers draw |
| `test_experimental_data.py` | water and silica built from scratch and refined against real neutron data; see the next section |
| `test_cu_btc.py` | the Cu-BTC crystal read from a PDB file, refined with Bragg peaks against real neutron data, compared with the archived Fortran results after every EPSR run; see the next section |

A session fixture in `conftest.py` fails the run if any test modifies the
repository's `EPSR-test` folder.

## Tests on real experimental data

`test_experimental_data.py` refines two systems against real neutron data from
`experimental_data/` (the Dissolve examples), set up as described in
`validation/README.md` and `validation/PYTHON_WRAPPER_WORKFLOW.md`:

| case | system | data |
|---|---|---|
| water | 1000 H2O, 0.1 atoms/Å³, 300 K | SANDALS H2O, D2O and HDO (50:50), `experimental_data/water/data/*.mint01` |
| silica_bulk | 1500 Si + 3000 O, 0.0664 atoms/Å³, 298 K | NIMROD SiO2, `experimental_data/silica_bulk/data/SiO2_NIMROD.mint01` |

Each case is built from scratch by following `validation/<case>/steps.txt`
through the wrapper. The setup programs (`makemole`, `mixato`, `randomise`,
`fmole`, `nwts`) run through `epsr.tools.run`. The EPSR runs (scratch start,
equilibration on the reference potential, then refinement with the empirical
potential) run through `epsr.Simulation`. Where `steps.txt` edits a keyword
(`niter` in the `.EPSR.inp`, `ereq` in the `.pcof`), the wrapper passes it as
an override to the next `Simulation`, as a Python user would. The helpers are
in `tests/epsr_workflow.py`.

**Default tests** use a shortened schedule: 10 instead of 99 or 100 iterations
per run, about 30 s per case. They check that:

* the boxes have the expected side, atom count and density, and the neutron
  weights are right (HDO HW–HW weight 0.02147, the sanity value of the workflow
  notes; H2O and D2O have O–H weights of opposite sign);
* every data file is read and fitted over its Q range;
* switching on the empirical potential lowers the overall and per-dataset
  R-factors, the potential stays within `ereq`, and the `ereq` override is saved
  to the `.pcof` file;
* the structure is physical: water O–H bond 0.976 Å and H–H 1.570 Å inside the
  molecule, O–O first neighbours at 2.8 Å; silica Si–O 1.6 Å, O–O 2.6 Å,
  Si–Si 3.1 Å; every g(r) tends to 1;
* the wrapper writes **byte-identical** files to the `epsr-rs` CLI after every
  EPSR run. The CLI side applies the keyword edits to the files with
  `validation/set_keyword.py`, so this also shows that overrides and file edits
  are equivalent.

**Full schedule against the original Fortran** (`test_full_schedule_matches_fortran`)
runs the complete `steps.txt` (1 + 99 + 100 + 100 iterations) through the wrapper
and through the Fortran programs, and compares every file after every EPSR run
with the tolerance policy below. It is skipped unless a Fortran build is given:

```bash
rust/tools/build_fortran_reference.sh ~/epsr-work/fortran-build
cd python
EPSR_FORTRAN_BUILD=~/epsr-work/fortran-build uv run pytest -n 2 \
    tests/test_experimental_data.py -k full_schedule
```

It then checks the converged structure: the water hydrogen-bond peak and
coordination numbers, and the silica SiO4 network.

Results on the development machine (17 October 2026, Fortran built with
`rust/tools/build_fortran_reference.sh`): every file of both cases was
**byte-identical** to the Fortran after every EPSR run (39 files for water, 36 for
silica). The two runs together take about 10 minutes on two cores. After the full
schedule:

| | water | silica |
|---|---|---|
| raw R-factor after 100 / 200 / 300 iterations | 6.0e-4 / 2.3e-4 / 2.3e-4 | 1.25e-3 / 3.7e-4 / 1.25e-4 |
| empirical-potential amplitude (ereq) | 3.0 (3) | 23.7 (50) |
| first peaks | O–O 2.79 Å, O···H 1.80 Å | Si–O 1.59 Å, O–O 2.55 Å, Si–Si 3.09 Å |
| coordination | 4.6 O around O (3.4 Å), 1.6 H-bonded H per O | 3.8 O per Si, 1.9 Si per O (2.2 Å) |

Coordination numbers integrate g(r) with the trapezoidal rule because the EPSR r
grid is not uniform: its first interval is a quarter of the others.

`validation/<case>/gr_agreement.png` shows the partial g(r) of each case from the
three codes side by side (original Fortran, Rust CLI, Python wrapper) with the
differences underneath, which are zero; `validation/plot_gr_agreement.py` makes them.

### Cu-BTC: a crystal from a PDB file, against archived Fortran results

`test_cu_btc.py` runs the Cu-BTC (HKUST-1) case of `validation/cu-btc/` through the
wrapper: the unit cell of the Dissolve example in `experimental_data/cu-btc`,
converted to a PDB file, read with `readpdb` (24 Cu₂ paddlewheels and 32
benzene-tricarboxylate linkers), given the Dissolve force field with `changeato`,
weighted with `nwts`, and refined against the NIMROD data of the empty framework
with Bragg peaks calculated up to 5 Å⁻¹. The schedule is the full one of
`steps.txt` (1 + 99 + 100 + 100 iterations, about 4 minutes).

The original Fortran's folder after each of the four EPSR runs is archived in
`validation/cu-btc/fortran_reference.tar.gz` (made by `make_fortran_reference.py`
from a Fortran build), so the comparison needs no Fortran build and runs in the
default suite. After every EPSR run, every file of the wrapper's folder is compared
with the archive under the tolerance policy below; files that agree only in the
last digits trigger a replay with the `epsr-rs` CLI, whose final folder must be
byte-identical to the wrapper's.

Results on the development machine: **every file byte-identical to the Fortran
after every run** (34 to 36 files per snapshot, including the Bragg files
`.bragg` and `.inten`). Further checks: the molecules and atom types `readpdb`
creates, the force field `changeato` applies, the 21 neutron weights, the Bragg
files, the R-factor falling from 3.1×10⁻² (scratch start) to 2.3×10⁻² (reference
potential) to 1.6×10⁻² (empirical potential at `ereq` 10), and the structure:
intramolecular Cu–Cu 2.64 Å, C1–O 1.26 Å and C3–H 0.93 Å, about four O around each
Cu, paddlewheels 8 Å apart. `validation/cu-btc/README.md` explains the conversion
and its approximations.

## Tolerance policy

A file passes when it is identical after masking (the date line of `.out` and the
timing columns of `fort.35`). Otherwise, unless `EPSR_STRICT=1`, it passes when
every line has the same text and each differing number satisfies one of:

* relative difference at most 1e-5;
* difference at most one unit in the last printed digit (for example
  `0.21800E-02` against `0.21801E-02`);
* difference at most 1e-5 of the largest value in the same row or column
  (values that are zero by symmetry and come out as float32 noise, such as
  `0.0000000E+00` against `-0.2523153E-08`).

Integers must match exactly. Lines are compared token by token, because a sign
shifts fixed-width Fortran fields (`-0.00000` against `0.00000`).

Every case that passes only within tolerance is also replayed through the
`epsr-rs` CLI, and the wrapper's files must be **byte-identical** to the CLI's.
That proves that any remaining difference from the Fortran comes from the engine
on this platform, not from the wrapper.

## Results on the development machine

Linux (WSL2), glibc 2.39, Rust 1.99, Python 3.12, 17 October 2026:

| | cases |
|---|---|
| byte-identical to the Fortran reference | 56 |
| within tolerance of the Fortran, and byte-identical to the CLI | 59 |
| different | 0 |

The 59 cases within tolerance are all 49 cases that run the `epsr` engine (the
`.o01` reference-potential table always, plus last-digit differences in `.erg`,
`.terg`, `.out`, `delfit.dat` and a few result tables in some cases), nine of the
`sharm`/`sdf`/`sdfcube`/`mapgr` cases ex70–ex79 (ex76 is identical), and
`torangles` (ex54). The other 56 tool cases are byte-identical.

The Rust engine's own regression test shows the same behaviour on this machine:
`cargo test --release` in `rust/` reports `200DMSO400Water.EPSR.o01` as differing
from the Fortran hash. That table holds values of order 1e13 at very small r,
and their 7th significant digit depends on the libm rounding of `powf`, which
differs between glibc versions. `test_reference.py` reports this one file as a
warning instead of a failure unless `EPSR_STRICT=1`; `test_ex_cases.py` checks its
numbers within the tolerance above.

## Adding a test

* For a new accessor, compare it with the output file that holds the same
  quantity, as `test_project_io.py::test_read_table` does for `.g01` and `.u01`.
* Use the `project` fixture (a fresh copy of `EPSR-test`) for anything that runs or
  writes; the module-scoped `ran` fixture holds one completed run for read-only checks.
* `tests/epsr_cases.py` has the comparison helpers (`compare_dirs`,
  `compare_file`), the case replay (`replay`, `replay_cli`) and `cli_binary()`.
