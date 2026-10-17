# Task: build the `epsr` Python wrapper for the Rust EPSR engine

You are implementing a Python package that wraps the Rust translation of EPSR
(Empirical Potential Structure Refinement) in this repository. The full specification is
`python/PLAN.md`. Read it completely before you write any code, then read the documents it
lists in its first section, in that order. The plan takes precedence over
`rust/docs/PYTHON_WRAPPER_NOTES.md` wherever the two differ.

## What to build

Everything in `python/PLAN.md`: the `python/` maturin project with the PyO3 extension
`epsr._core`, the pure-Python `epsr` package (`Project`, `Simulation`, `Potential`,
`refine`, `History`, `io`, `tools`, `plot`), the four example scripts, the pytest suite,
`README.md`, `.gitignore`, `uv.lock`, and an empty `CHANGES_TO_RUST.md`. Work through the
phases in the plan's section 7 in order; each phase has a "done when" criterion.

## Hard constraints

1. **Do not change the structure, file names or function names of `rust/epsr-core` or
   `rust/epsr-cli`.** Other agents are editing that code concurrently. Do not add members
   to `rust/Cargo.toml`. Expect to make zero edits there; everything you need is already
   `pub`. If an engine change is truly unavoidable, add one new method at the end of
   `rust/epsr-core/src/api.rs` and record it in `python/CHANGES_TO_RUST.md`.
2. **Package management is `uv` only**, with maturin as the build backend. Commit
   `uv.lock`. No pip, venv, poetry or conda anywhere.
3. **Use existing libraries** (numpy required; pandas, matplotlib, tqdm optional; pytest,
   ruff, mypy dev). Do not reimplement physics, CSV writing, progress bars or plotting.
4. **Preserve EPSR semantics.** One `Simulation` = one Fortran `epsr` invocation. Never
   call `iterate` past `niter`. Files are modified in place; the folder copy helper is the
   documented path.
5. Do not reopen the decisions listed in the plan's section 8.

## Environment facts

* Rust 1.99 is installed at `~/.cargo/bin` but is not on `PATH` in non-login shells.
* `uv` is **not** installed. Install it with the official standalone installer
  (`curl -LsSf https://astral.sh/uv/install.sh | sh`); it lands in `~/.local/bin`.
* System Python is 3.12 with numpy 1.26; headers are in `/usr/include/python3.12`.
* The reference test case is `EPSR-test/` at the repository root (200 DMSO + 400 water,
  three neutron datasets). Never run on it directly; copy it first. The engine's own
  regression test, `rust/epsr-core/tests/test_case.rs`, shows the expected values.
* Additional verified cases for the auxiliary tools are under
  `test_problems/200DMSO400Water/EPSR-test-ex*/` (each has `README.md`, `commands.txt`
  and any `*.stdin` answer files).

## How to verify

* `uv run pytest -q` in `python/` is green, including the byte-for-byte reference test
  (`test_reference.py`) and the `nwts` tool test that must reproduce the shipped
  `*.NWTStot.wts` files exactly.
* All four scripts in `python/examples/` run end to end against a copy of `EPSR-test`.
* `cargo test --release --manifest-path rust/Cargo.toml` still passes and
  `git status rust/` is clean.
* `uv run ruff check` and `uv run mypy epsr` are clean.

## How to work

* Commit after each phase with a message naming the phase. Do not push.
* Keep the Rust binding layer thin; put ergonomics in Python.
* When something in the plan turns out to be wrong about the engine (a field that is not
  `pub`, a layout that differs), fix the plan's assumption in the simplest way that keeps
  constraint 1, and note it in your final report.
* Final report: what was built, test results with output, anything left out and why, and
  the contents of `CHANGES_TO_RUST.md` if it is not empty.
