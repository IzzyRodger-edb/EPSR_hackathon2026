# epsr-py: Python interface to the Rust EPSR engine

`epsr` lets you run and inspect **EPSR** (Empirical Potential Structure Refinement,
A. K. Soper) from Python. The heavy lifting is done by the Rust translation of the
original Fortran in [`../rust`](../rust), which reproduces the Fortran output
byte for byte; this package adds a Pythonic API, numpy arrays, a refinement loop,
potential perturbation, the 35 setup and analysis tools, and plotting helpers.

```python
import epsr

proj = epsr.Project.from_template("../EPSR-test", "work/dmso")  # copy the folder first
sim = proj.simulation(overrides={"ntimes": 5})                  # read all inputs
sim.run()                                                       # one EPSR run (Rust, GIL released)

sim.summary.r_factor                 # 1.1478e-04
sim.gr_inter("O1", "H1")             # partial g(r) as a float32 numpy array
sim.potential.empirical()            # (n_pairs, nr) empirical potential, kJ/mol
hist = epsr.refine(proj, n_runs=20)  # repeated runs with convergence history
```

## Installation

The project is managed with [uv](https://docs.astral.sh/uv/) and built with
[maturin](https://www.maturin.rs/). You need a Rust toolchain (stable, ≥ 1.85) and
Python ≥ 3.10.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh     # if uv is not installed
export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
cd python
uv sync --all-extras          # creates .venv, compiles the Rust engine, installs epsr
uv run python examples/01_single_run.py
```

| command | purpose |
|---|---|
| `uv sync` | install `epsr` with numpy only |
| `uv sync --all-extras` | also matplotlib, pandas and tqdm |
| `uv run maturin develop --release --uv` | fast rebuild after editing `src/*.rs` |
| `uv run pytest -n auto` | the full test suite (about 1.5 minutes on 8 cores) |
| `uv run pytest -m "not slow"` | skip the 115 verification-case replays (about 1 minute) |
| `uv run isort . && uv run black .` | format the Python code (import order, then layout) |
| `uv run ruff check . && uv run mypy` | lint and type check |
| `uv build` | wheel and sdist in `dist/` |

`uv sync` rebuilds the extension automatically when the Rust engine or the
bindings change.

## Things to know before you start

* **EPSR rewrites its input files.** Every run overwrites `.ato`, `.pcof`,
  `.EPSR.inp` and `.EPSR.inpa`, appends to `.erg` and writes all result files in the
  simulation folder. Even creating a `Simulation` rewrites the `.qdr` file. This is
  how EPSR carries its state from run to run. Always work on a copy:
  `epsr.Project.from_template(original, work)`.
* **One `Simulation` is one EPSR run.** It reads the folder, runs `niter`
  iterations (Monte Carlo cycles plus an empirical-potential update) and saves the
  state. For the next run create a new `Simulation`, or let `epsr.refine` do the loop.
* **Overrides are sticky.** Parameters passed as `overrides` are written back to the
  input files, so they persist into later runs.
* **Threads.** `run()` releases the GIL. Simulations on different folders can run in
  parallel threads (`epsr.refine_many`). Never run two on the same folder.
* **Precision.** Arrays are float32 copies, the precision EPSR computes in.
* **Platform.** Bit-for-bit agreement with the Fortran holds on Linux with glibc for
  most files; a few files can differ in the last printed digit on other libm
  versions. See [docs/testing.md](docs/testing.md).

## What is in the package

| module | contents |
|---|---|
| `epsr.Simulation` | one run: `run`, `iterate`, `iterations`, `summary`, grids `q`/`r`, partials `gr_inter`/`gr_intra`/`sq_inter`/`sq_intra`, data fits `data_fq`/`simulated_fq`/`fit`, configuration `positions`, `potential` |
| `epsr.Potential` | `sim.potential`: empirical, reference and total potentials; EP coefficients get/set/scale/zero/perturb; potential and atom-type parameters |
| `epsr.Project` | a simulation folder: `from_template`, `simulation`, `snapshot`, `read_table`, `history`, `files` |
| `epsr.refine`, `epsr.refine_many`, `epsr.History` | repeated runs, convergence history, CSV/pandas export |
| `epsr.tools` | the 35 auxiliary programs (`nwts`, `partials`, `makeato`, `mixato`, `fmole`, `sharm`, ...), keyword tables, `write_keyword_file` |
| `epsr.io` | read-only parsers for result tables (`.g01`, `.u01`, ...), `.erg`, keyword files |
| `epsr.plot` | matplotlib helpers: `fits`, `partials`, `potentials`, `history` |
| `epsr.errors` | `EpsrError` and its subclasses `EpsrIOError` (an `OSError`), `EpsrInputError` (a `ValueError`) and `EpsrUnsupportedError` (a `NotImplementedError`) |

## Documentation

* [docs/user_guide.md](docs/user_guide.md): workflow, refinement, perturbing the potential, tools, parallel runs.
* [docs/api.md](docs/api.md): reference for every public class and function.
* [docs/architecture.md](docs/architecture.md): how the wrapper is built and how to extend it.
* [docs/testing.md](docs/testing.md): the test suite, tolerances and verification results.
* [examples/](examples): four runnable scripts: a single run, a refinement loop, potential perturbation, tools.

## Layout

```
python/
├── pyproject.toml, uv.lock     project definition (maturin backend, uv lock file)
├── Cargo.toml, src/*.rs        PyO3 bindings: the extension module epsr._core
├── epsr/                       the Python package
├── examples/                   runnable scripts
├── tests/                      pytest suite (wrapper equivalents of the Rust tests)
├── docs/                       documentation
├── PLAN.md, PROMPT.md          the specification this package was built from
└── CHANGES_TO_RUST.md          changes made to rust/ for the wrapper (none)
```

The Rust engine under `../rust` is used unchanged through its public API.
