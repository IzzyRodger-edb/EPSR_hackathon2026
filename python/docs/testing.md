# Testing

The suite in `python/tests/` checks the wrapper against the same references the
Rust engine is validated with, plus wrapper-specific behaviour.

```bash
cd python
uv run pytest -n auto                 # everything: 188 tests, about 1.5 minutes on 8 cores
uv run pytest -m "not slow"           # skip the 115 case replays: about 1 minute
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
| `test_threads.py` | GIL release, parallel runs equal sequential runs, `refine_many`, tools from several threads |
| `test_examples.py` | the four example scripts run |
| `test_plot.py` | the plotting helpers draw |

A session fixture in `conftest.py` fails the run if any test modifies the
repository's `EPSR-test` folder.

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
