# epsr-py: Python interface to the Rust EPSR engine

`epsr` lets you run and inspect **EPSR** (Empirical Potential Structure Refinement,
A. K. Soper) from Python. The heavy lifting is done by the Rust translation of the
original Fortran in [`../rust`](../rust), which reproduces the Fortran output
byte for byte; this package adds a Pythonic API, numpy arrays, a refinement loop,
potential perturbation, the 35 setup and analysis tools, and plotting helpers.

## Example: refining liquid water against real neutron data

This builds 1000 water molecules from scratch and refines them against the SANDALS
neutron data for H2O, D2O and HDO (50:50) in
[`../experimental_data/water/data`](../experimental_data/water/data). The setup
inputs (molecule template, neutron-weight files, EPSR input file) are those of
[`../validation/water`](../validation/water). Run it from the `python/` folder
after [installing](#installation):

```bash
uv run python examples/05_experimental_water.py            # full schedule, about 4 minutes
uv run python examples/05_experimental_water.py --quick    # 10 iterations per run, about 30 s
```

The core of [examples/05_experimental_water.py](examples/05_experimental_water.py):

```python
import shutil
from pathlib import Path

import epsr
from epsr import tools

work = Path("work/water")

# 1. Setup inputs plus the measured data, in a fresh folder.
epsr.copy_folder("../validation/water", work, overwrite=True)
for data_file in Path("../experimental_data/water/data").glob("*.mint01"):
    shutil.copy2(data_file, work / data_file.name)

# 2. Build the box with the EPSR setup programs.
tools.run(work, "makemole", ["Water"])
tools.run(work, "mixato", answers=["1", "Water", "1000", "0.1", "waterbox"])
tools.run(work, "randomise", ["waterbox"])
(work / "Water.mol").unlink()  # fmole refuses to run while the .mol is present
tools.run(work, "fmole", ["waterbox", "100"])
for sample in ("H2O", "D2O", "HDO"):
    tools.nwts(work, f"{sample}.NWTS.dat")

# 3. Scratch start, then equilibrate on the reference potential.
project = epsr.Project(work, "waterbox.EPSR.inp")
project.simulation().run()
project.simulation(overrides={"niter": 99}).run()

# 4. Switch on the empirical potential (ereq 3 kJ/mol) and refine twice.
history = epsr.refine(
    project,
    n_runs=2,
    overrides={"niter": 100},
    initial_overrides={"ereq": 3.0},
)

# 5. Inspect the result.
simulation = project.simulation()
print(history.r_factor[-1])           # overall R-factor of the fit
simulation.fit(0)                     # data and simulated F(Q) of the D2O dataset
simulation.gr_inter("OW", "HW")       # O...H partial g(r): hydrogen-bond peak at 1.8 Å
simulation.potential.empirical()      # the refined empirical potential, kJ/mol
epsr.plot.fits(simulation)            # data and fits, with the `plot` extra
```

Output of the full schedule:

```
R-factor: 6.136e-04 on the reference potential
          2.268e-04 after refinement with the empirical potential
  SLS18502-D2O       fitted over Q = 0.20 to 29.90 1/Å
  SLS18500-HDO5050   fitted over Q = 0.20 to 29.90 1/Å
  SLS18498-H2O       fitted over Q = 0.20 to 29.90 1/Å
OW-OW first peak at 2.79 Å (nearest-neighbour waters)
OW-HW first peak at 1.80 Å (hydrogen bond)
```

The same workflow for amorphous silica, and a check that every output file is
identical to the original Fortran EPSR, are described in
[docs/testing.md](docs/testing.md#tests-on-real-experimental-data).

## Installation

The project is managed with [uv](https://docs.astral.sh/uv/) and built with
[maturin](https://www.maturin.rs/). You need a Rust toolchain (stable, ≥ 1.85) and
Python ≥ 3.10.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh     # if uv is not installed
export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
cd python
uv sync --all-extras          # creates .venv, compiles the Rust engine, installs epsr
uv run python examples/05_experimental_water.py --quick   # the example above
```

| command | purpose |
|---|---|
| `uv sync` | install `epsr` with numpy only |
| `uv sync --all-extras` | also matplotlib, pandas and tqdm |
| `uv run maturin develop --release --uv` | fast rebuild after editing `src/*.rs` |
| `uv run pytest -n auto` | the full test suite (about 4.5 minutes on 4 cores) |
| `uv run pytest -m "not slow"` | skip the 115 verification-case replays and the experimental-data refinements (about 1 minute) |
| `uv run pytest -m experimental` | refine real water and silica neutron data from `experimental_data/` (about 2 minutes) |
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
* **Threads.** The engine uses all CPU cores inside each run; set `RAYON_NUM_THREADS`
  to limit it. Results are identical at any thread count. `run()` releases the GIL,
  so simulations on different folders can also run in parallel threads
  (`epsr.refine_many`). Never run two on the same folder.
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
* [examples/](examples): five runnable scripts: a single run on the test case, a refinement loop, potential perturbation, tools, and the water refinement against experimental data shown above.

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
