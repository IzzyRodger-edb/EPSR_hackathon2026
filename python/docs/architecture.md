# Architecture of the wrapper

For maintainers. The wrapper has two layers on top of the unchanged Rust engine.

```
 Python user code
        │
 epsr (pure Python)        epsr/*.py     descriptive names, dataclasses, 2-D arrays, Project,
        │                                refine loop, Potential, tools helpers, io parsers, plots
 epsr._core (PyO3)         src/*.rs      thin bindings with Fortran names; float32 copies;
        │                                error mapping; GIL released around engine calls
 epsr-core (Rust)          ../rust/epsr-core   the engine, used through its public API only
```

## Ground rules

* **The engine is not modified.** Other people work on `rust/` concurrently. The
  bindings use only items that are already `pub`: `Epsr` and its `api.rs`
  accessors, the public fields `st`, `ctl`, `data`, `errs`, `wav`, `qwind`,
  `summary`, and the methods `generate_ep`, `multiply_pot_coeffs`,
  `tools::run_tool`, `symtab_tools::ALL`. If an engine change ever becomes
  unavoidable, add a new accessor at the end of `rust/epsr-core/src/api.rs` and
  record it in `CHANGES_TO_RUST.md`.
* **No physics in Python.** Python reads results, sets inputs and drives runs.
  Parsing EPSR's text outputs with numpy is the only computation it does.
* **EPSR semantics are kept.** One `Simulation` is one invocation of the Fortran
  `epsr` program; files are updated in place; iterations run in order up to `niter`.

## Build and packaging

* `python/` is its own Cargo workspace root (`[workspace]` in `Cargo.toml`) with a
  path dependency on `../rust/epsr-core`, so `cargo test` in `rust/` never needs
  PyO3 or Python headers. The engine compiles into `python/target/`.
* PyO3 0.30 with the `abi3-py310` feature: one wheel serves Python ≥ 3.10. The
  `extension-module` feature no longer exists in PyO3 0.30; maturin configures it.
* maturin is the build backend (`[tool.maturin] module-name = "epsr._core"`, mixed
  Rust/Python layout with `epsr/` next to `pyproject.toml`). uv manages the
  environment and `uv.lock`; `[tool.uv] cache-keys` include the engine sources so
  `uv sync` rebuilds after engine changes.
* `[profile.dev.package."*"] opt-level = 3` keeps debug builds usable: the engine
  is very slow unoptimised. maturin builds with the release profile by default.

## The binding layer (`src/`)

| file | contents |
|---|---|
| `lib.rs` | module definition: `Simulation`, `run_tool`, `tool_commands`, `tool_tables`, `tool_keywords`, `VERSION_DATE` |
| `sim.rs` | `#[pyclass] Simulation` wrapping `epsr_core::Epsr` plus the number of the last completed iteration |
| `tools.rs` | `run_tool` (`ToolContext` + `tools::run_tool`), keyword tables from `symtab_tools::ALL` |
| `convert.rs` | error mapping, slice to numpy helpers |

Design points:

* **Arrays are copies.** `array_1d`/`array_2d` copy the used part of the 1-based Fortran
  columns (`F2::col(j)[..n]`) into new numpy arrays. Sizes are small (21 pairs ×
  600 points for the test case), and copies avoid lifetime problems.
* **GIL.** `Simulation.__new__`, `run`, `iterate` and `run_tool` run inside
  `Python::detach` (PyO3's former `allow_threads`). `Epsr` is `Send`, which the
  engine's own test `epsr_is_send` checks.
* **Engine threads.** The engine parallelises its inner loops with rayon (see
  `PARALLELISM_PLAN.md` at the repository root). The wrapper needs no code for
  this: the rayon work happens inside the detached call. All simulations of a
  process share rayon's global pool, sized by `RAYON_NUM_THREADS`. The engine keeps
  the box geometry in a thread-local variable that every parallel task re-syncs;
  `test_threads.py` checks thread-count invariance and concurrent runs with a
  cubic and a triclinic box through the wrapper.
* **Run order.** `sim.rs` tracks `last_iteration` so that `run` happens once and
  `iterate(iteration)` only accepts `last_iteration + 1 ≤ niter`. After the last iteration
  `iterate` also writes `fort.35`, so stepping through iterations produces the
  same files as `run()`.
* **Errors.** `EpsrError::{Io, Input, Unsupported}` map to the Python classes in
  `epsr/errors.py`. They are defined in Python rather than with
  `create_exception!` because they use multiple inheritance with `OSError`,
  `ValueError` and `NotImplementedError`. `convert::to_py` imports
  `epsr.errors` when raising.
* **EP writes.** `set_ep_coefficients` writes `pot.pcoefp[(ic, irefpcof[ip])]` for
  every pair and calls `generate_ep(false)`, the same call the engine makes after
  each refinement; `scale_ep` calls `multiply_pot_coeffs` and then
  `generate_ep(false)`. Every write therefore goes through the engine's `ereq`
  amplitude limit.

### Name mapping

| Python (`epsr`) | `_core` | engine |
|---|---|---|
| `Simulation(...)` | `Simulation(folder, inp, extra_lines, verbose)` | `Epsr::setup_with_lines` |
| `run`, `iterate` | `run`, `iterate` | `Epsr::run`, `Epsr::iterate` (+ `write_fort35`) |
| `summary` | `summary` (dict) | `Epsr.summary` |
| `control`, `datasets` | dicts | `Epsr.ctl`, `Epsr.data` |
| `q`, `r` | `q_grid`, `r_grid` | `api::q_grid`, `api::r_grid` |
| `gr_inter` ... `sq_intra` | same, 2-D | `api::gr_inter(ip)` ... per pair |
| `current_gr`, `current_sq` | same | `Epsr.wav`, `Epsr.qwind` |
| `data_fq`, `simulated_fq`, `data_q_range` | same | `api::...` |
| `data_errors` | same | `Epsr.errs` |
| `positions`, `box_size`, `density` | `atom_positions`, `box_size`, `rho` | `api::...`, `st.dist.rho` |
| `Potential.empirical/reference/total` | `empirical_potential`, `reference_potential`, `total_potential` | `xavs`, `ecoul`, `xavtot` |
| `Potential.ep_coefficients` | `ep_coefficients` (with flag column) | `api::ep_coefficients` (`pcoefp`) |
| `Potential.set_ep_coefficients`, `scale_ep` | same | `pcoefp`, `generate_ep`, `multiply_pot_coeffs` |
| `Potential.params`, `type_params` | `pot_params`, `type_params` | `st.pot.*`, `st.mol.{epscor,sigcor,charge,amassa,cfrac}` |
| `tools.run` | `run_tool` | `tools::run_tool(ToolContext)` |
| `tools.keywords` | `tool_keywords` | `symtab_tools::ALL` |

## The Python layer (`epsr/`)

| module | role |
|---|---|
| `simulation.py` | `Simulation`, the `Summary`/`Control`/`Dataset`/`Fit` dataclasses, `render_overrides` |
| `potential.py` | `Potential`, `PotParams`, `TypeParams`; checks shared `.pcof` columns before writes |
| `project.py` | `Project`, `copy_folder`: folder copies, snapshots, output paths |
| `refine.py` | `refine`, `refine_many`, `History` |
| `tools.py` | `run` (switches to the folder under a lock), keyword tables, `write_keyword_file`, typed helpers |
| `io.py` | read-only parsers for result tables, `.erg` and keyword files |
| `plot.py` | matplotlib helpers, imported lazily |
| `errors.py` | exception classes |
| `_core.pyi` | type stubs for the extension |

### Why `tools.run` changes the working directory

The box-building programs read `system_commands.txt` and write `dirlist.txt` in
the *current working directory*, exactly like the Fortran (see the header of
`rust/epsr-core/src/tools/setup_common.rs`): the EPSR shell runs every program
from inside the simulation folder. Without a change of directory, `mixato` and
friends cannot find their files and `dirlist.txt` lands wherever Python was
started. `tools.run` therefore switches to the folder, passes `./` (as the Rust
verification script `rust/tools/run_ex_test.sh` does, which also makes the
`writexyz` headers match the Fortran reference), and restores the directory
afterwards. Because the working directory is process-wide, a re-entrant lock
serialises tool runs. `Simulation` does not depend on the working directory.

### Overrides

`Epsr::setup_with_lines` processes extra lines after reading `.EPSR.inp` and
before the run parameters are parsed from the keyword table, so overrides reach
`Control` and the potential parameters (`ereq` → `absolute_energy`). Applying
lines to an existing object later would leave `Control` out of sync with the
table, so overrides are only accepted at construction.

## Extending the wrapper

To expose another engine quantity:

1. Check whether it is reachable through public fields (`self.engine.st.<struct>.<field>` in `sim.rs`;
   see the state tables in `rust/docs/ARCHITECTURE.md`). Remember that the arrays
   are 1-based and column-major with fixed leading dimensions: slice to the used length.
2. Add a method to `src/sim.rs` returning a copy, add its signature to
   `epsr/_core.pyi`, and wrap it with a descriptive name in the Python layer.
3. Add a test. If the quantity is also written to an output file, compare with
   that file (see `test_project_io.py::test_read_table`).
4. Rebuild with `uv run maturin develop --release --uv`.

## Code style

* **Names.** Use descriptive names for variables, parameters and helpers
  (`pair_index`, `dataset`, `coefficients`), never single letters or
  abbreviations. Fortran names appear only where they mirror the engine: the
  `_core` method names, the keys of the dicts it returns, and dataclass fields
  that copy engine structs (`Control`, `Dataset`), each documented with a comment.
* **Formatting.** isort orders imports and black lays out the code, both with a
  line length of 88 (black's default) (`[tool.isort]`, `[tool.black]` in `pyproject.toml`). Run
  `uv run isort . && uv run black .` before committing; ruff lints (it does not
  format) and mypy checks types. The Rust bindings use `cargo fmt` and
  `cargo clippy`.
* **Tests.** Every change is followed by `uv run pytest -n auto`.

## Deviations from PLAN.md

The implementation follows `PLAN.md` except where the engine's actual behaviour
required otherwise:

| plan | implementation | reason |
|---|---|---|
| exceptions created in Rust with `create_exception!` | defined in `epsr/errors.py`, raised from Rust | they need multiple inheritance (`EpsrIOError` is also an `OSError`, ...) |
| `tools.run` calls the engine directly | switches to the folder under a lock and passes `./` | the box-building programs use the current working directory, as in the Fortran |
| test: `scale_ep(2.0)` doubles the EP | doubling is undone when the EP is already at `ereq`; tests check linear scaling below the limit and the cap above it | `generate_ep` applies the `ereq` amplitude limit after every rebuild |
| examples use pair `("O", "H")` for water | `("O1", "H1")` | in the test case `O`/`H` are DMSO atoms, `O1`/`H1` water |
| example 04 runs `makeato` | runs `nwts`, `partials` and `writexyz` | `makeato` deliberately reproduces a Fortran crash (empty `.ato`) |
| PyO3 0.22 sketch, `extension-module` feature | PyO3 0.30 (`Python::detach`), no such feature | current PyO3; maturin sets the extension-module configuration |
| — | added `_core.ep_columns`, `_core.tool_tables`, `Potential.shared_pairs`, `refine(initial_overrides=...)` | safe EP writes for pairs sharing `.pcof` coefficients; listing keyword tables; one-off overrides such as `ireset 1` |
| — | typed tool helpers raise `FileNotFoundError` for a missing keyword file | the engine, like the Fortran, silently does nothing |
| — | the last `iterate()` also writes `fort.35` | so that stepping through iterations writes the same files as `run()` |
| test suite: the reference test and `nwts` | also replays all 115 `EPSR-test-ex*` cases, with a tolerance policy and CLI cross-check | the user asked for wrapper equivalents of the Rust tests |
