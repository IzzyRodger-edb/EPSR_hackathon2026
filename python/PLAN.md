# Plan: Python wrapper for the Rust EPSR engine

This document is the specification for building the `epsr` Python package on top of
the Rust engine in `rust/epsr-core`. It is written for the agent that will implement
it. Read it fully before writing code, then read, in this order:

1. `rust/docs/PYTHON_WRAPPER_NOTES.md` — the engine author's notes for the wrapper (API table, side effects, determinism).
2. `rust/docs/ARCHITECTURE.md` — the state model and the data flow of one iteration.
3. `rust/epsr-core/src/api.rs` — the 0-based accessors you will bind.
4. `rust/epsr-core/src/epsr.rs` lines 30–170 and 580–600 — `Dataset`, `Control`, `Summary`, `Epsr`, `setup_with_lines`, `iterate`, `run`.
5. `rust/epsr-core/src/tools/mod.rs` — `ToolContext`, `COMMANDS`, `run_tool`.
6. `rust/epsr-core/tests/test_case.rs` and `rust/tools/make_reference_hashes.py` — the regression test you will port to pytest.

## Relationship to `rust/docs/PYTHON_WRAPPER_NOTES.md`

The engine author's notes are the source for the Rust API table, the in-place file
semantics, determinism, the GIL/`Send` argument and the testing approach, and this plan
follows them. Where the two differ, **this plan takes precedence**; the differences are
deliberate:

| topic | notes say | this plan | why |
|---|---|---|---|
| crate location | `rust/epsr-py` | `python/` with its own Cargo root | keeps `rust/` untouched for the other agents; the notes themselves ask for it to be outside the workspace's default members |
| user entry point | `epsr.Simulation(folder, inp)` | same, plus `Project.from_template(...).simulation()` | makes the mandatory folder copy the default path |
| `summary` | a dict | `Summary` dataclass with descriptive names and a `.fortran` raw dict | IDE discoverability; Fortran names stay reachable |
| `run_many` | loop in Rust, optional | loop in Python (`refine`) | callbacks, early stop and snapshots are easier in Python; the GIL is released during each run anyway |
| potential perturbation | not covered | §4.3 `Potential` (in-memory EP coefficient edits via the public `pcoefp` + `generate_ep`) | required by the project brief |
| `read_table`, `prepare` helpers | suggested | `io.read_table`, `Project.from_template` | same intent, placed in modules |

Everything else in the notes (array copies as float32, error mapping, `tools::run_tool`
wrapping, pytest port of the hash test, the `.dat` keyword-file helper) is adopted as is.

## 0. Ground rules

* **Do not modify the structure or function names of `rust/epsr-core` or `rust/epsr-cli`.**
  Other agents work on that code concurrently. Do not rename, move, split or change the
  signature of anything there, and do not add `rust/` workspace members.
* **Every Rust binding lives in the new crate `python/` (package `epsr-py`).** All state
  needed for the features below is already `pub` (`Epsr.st`, `st.pot`, `st.mol`, `ctl`,
  `data`, `datai`, `errs`, `wav`, `qwind`, `summary`) and the methods you need
  (`generate_ep`, `generate_rp`, `multiply_pot_coeffs`, `moloutc`, `apply_input_lines`)
  are `pub`. The expectation is **zero edits to `epsr-core`**. If something is genuinely
  impossible without an engine change, add a *new* accessor method at the bottom of
  `rust/epsr-core/src/api.rs` (additive only), and list it in `python/CHANGES_TO_RUST.md`.
* **Do not reimplement physics in Python.** The engine is validated byte-for-byte against
  Fortran. Python reads results, sets inputs, and drives runs. The one exception is parsing
  text output files (`.erg`, `.?01`) with numpy, which is harmless.
* **Preserve EPSR's semantics.** One `Simulation` object = one Fortran `epsr` invocation
  (`setup` + `niter` iterations). Files in the folder are updated in place at the end;
  that *is* how EPSR carries state between runs. Never call `iterate` past `niter` on the
  same object; make a new `Simulation` for the next run.
* **Keep the Rust layer thin and the Python layer ergonomic.** The extension module
  `epsr._core` mirrors the Rust API almost 1:1 (Fortran names allowed). The pure-Python
  package `epsr` provides descriptive names, dataclasses, numpy 2-D arrays, the
  refinement loop, folder management and plotting.
* **Use existing Python libraries rather than hand-rolling.** Required: `numpy`
  (arrays, `loadtxt`/`genfromtxt` for the text outputs). Optional extras: `pandas`
  (`History.to_dataframe()`, `Table.to_dataframe()`), `matplotlib` (`epsr.plot`),
  `tqdm` (progress bar in `refine` when installed). Dev: `pytest`, `maturin`, `ruff`,
  `mypy`. Do not write a CSV writer, an argument parser, a progress bar or a plotting
  layer of your own when one of these does the job. Do not add heavy dependencies
  (`scipy`, `xarray`, `ase`, …) unless a feature in this plan needs them, and then
  only as an optional extra.
* **Package management is `uv`.** The project is defined in `pyproject.toml` with the
  maturin build backend; `uv` manages the virtual environment, the lock file and the
  dev/optional dependency groups. Commit `uv.lock`. Never use `pip`, `venv`, `poetry`
  or `conda` in scripts, docs or CI.

## 1. Target user experience

These scripts define "done". Write them as `python/examples/*.py` and make them run.

```python
# examples/01_single_run.py
import epsr

proj = epsr.Project.from_template("EPSR-test", "work/dmso")   # copies the folder
sim = proj.simulation(overrides={"ntimes": 5, "niter": 1})     # Epsr::setup_with_lines
sim.run()                                                      # releases the GIL
print(sim.summary.r_factor, sim.summary.energy, sim.summary.pressure)

q, r = sim.q, sim.r                                            # float32 numpy arrays
for i, ds in enumerate(sim.datasets):
    data, fit = sim.data_fq(i), sim.simulated_fq(i)            # (nq,) each
g = sim.gr_inter()                                             # (n_pairs, nr)
u = sim.potential.empirical()                                  # (n_pairs, nr)  kJ/mol
sim.pair_labels                                                # [("S","S"), ("S","C"), ...]
sim.positions                                                  # (n_atoms, 3) Å
```

```python
# examples/02_refinement_loop.py
import epsr

proj = epsr.Project.from_template("EPSR-test", "work/refine")
hist = epsr.refine(
    proj, n_runs=50,
    overrides={"ntimes": 5},
    callback=lambda k, sim: print(k, sim.summary.r_factor),
    until=lambda h: h.r_factor[-1] < 1e-4,
)
hist.r_factor, hist.energy, hist.ep_amplitude    # numpy arrays, one entry per run
hist.to_csv("work/refine/history.csv")
epsr.plot.history(hist)                          # matplotlib (optional dependency)
```

```python
# examples/03_perturb_potential.py
import epsr, numpy as np

proj = epsr.Project.from_template("EPSR-test", "work/perturb")

# (a) the standard EPSR knobs, applied as .EPSR.inp / .pcof keyword overrides
sim = proj.simulation(overrides={"ereq": 15.0, "feedback": 0.8, "potfac": [1.0, 1.0]})

# (b) direct manipulation of the empirical potential, in memory, before running
pot = sim.potential
coef = pot.ep_coefficients()            # (n_pairs, ncoeffp-1) float32, flag column removed
coef[pot.pair_index("O", "H")] *= 0.5   # damp one pair
pot.set_ep_coefficients(coef)           # writes pcoefp, calls generate_ep(false)
pot.scale_ep(0.9)                       # multiply_pot_coeffs + generate_ep
pot.zero_ep()                           # same as ireset 2
pot.empirical("O", "H")                 # (nr,) current EP table for that pair
pot.reference("O", "H")                 # LJ + Coulomb table
pot.total("O", "H")                     # reference + empirical (xavtot)

sim.run()                               # the modified EP is refined further and saved to .pcof

# (c) fix the EP and only run Monte Carlo (potfac 0 disables refinement)
sim2 = proj.simulation(overrides={"potfac": [0.0, 0.0]})
sim2.run()
```

```python
# examples/04_tools.py
import epsr

proj = epsr.Project.from_template("EPSR-test", "work/tools")
epsr.tools.run(proj.folder, "nwts", ["DMSODWaterD.NWTS.dat"])            # keyword-file tool
epsr.tools.run(proj.folder, "makeato", ["Ow"], answers=["0.65", "3.16", "15.9994", "-0.8476", "O", "300", "0.1"])
epsr.tools.write_keyword_file(proj.folder / "x.PARTIALS.dat", "partials", {"fnameato": "200DMSO400Water.ato", ...})
epsr.tools.commands()                                                     # list of 35 names
```

## 2. Repository layout (new files only)

```
python/                          the wrapper: one maturin project, its own Cargo root
├── PLAN.md                      this file
├── CHANGES_TO_RUST.md           empty unless an api.rs accessor had to be added
├── README.md                    install, quick start, side-effects warning, API summary
├── pyproject.toml               [build-system] maturin; [project] deps numpy; extras plot/pandas;
│                                [dependency-groups] dev (pytest, maturin, ruff, mypy); [tool.maturin]; [tool.uv]
├── uv.lock                      committed; produced by `uv sync`
├── .python-version              "3.12" (what is installed here; abi3 wheels still serve ≥ 3.10)
├── Cargo.toml                   package "epsr-py", [lib] name "_core", crate-type cdylib
│                                depends on epsr-core = { path = "../rust/epsr-core" }
│                                contains an empty [workspace] table (own root)
├── src/
│   └── lib.rs                   PyO3 bindings (thin). May be split into
│       ├── sim.rs               PySimulation
│       ├── convert.rs           error mapping, F1/F2 -> numpy helpers
│       └── tools.rs             run_tool binding, keyword tables
├── epsr/                        pure-Python package  (import epsr)
│   ├── __init__.py              re-exports: Project, Simulation, Summary, refine, History, tools, plot, errors
│   ├── _core.pyi                type stubs for the extension module
│   ├── errors.py                EpsrError, EpsrInputError(ValueError), EpsrUnsupportedError(NotImplementedError)
│   ├── project.py               Project: folder, run name, template copy, snapshots, file readers
│   ├── simulation.py            Simulation, Summary, Dataset, Control dataclasses, Potential handle
│   ├── potential.py             Potential class (EP / reference manipulation)
│   ├── refine.py                refine(), History
│   ├── io.py                    read_table(.?01), read_erg, read_inp/pcof keyword files (read-only)
│   ├── tools.py                 run(), commands(), keywords(), write_keyword_file()
│   └── plot.py                  matplotlib helpers (lazy import)
├── examples/                    the four scripts of §1
└── tests/
    ├── conftest.py              fixtures: template folder path, tmp copy, "slow" marker
    ├── test_reference.py        pytest port of rust/epsr-core/tests/test_case.rs
    ├── test_simulation.py       accessors, shapes, dtypes, overrides, errors
    ├── test_potential.py        EP get/set/scale/zero round trips and persistence to .pcof
    ├── test_refine.py           loop, history, callback early stop
    ├── test_tools.py            nwts reproduces the shipped .wts files byte-identically
    └── test_threads.py          GIL release; two folders in parallel threads give identical results to sequential
```

Why a separate Cargo root under `python/` rather than a member of `rust/Cargo.toml`:
it needs Python headers and PyO3, which must not become a requirement of `cargo test` in
`rust/`, and it keeps the engine's workspace untouched. Cargo accepts a path dependency on a
package that belongs to another workspace; the engine is simply compiled into
`python/target/`.

## 3. The Rust binding layer (`python/src/lib.rs`)

Toolchain: PyO3 (latest stable; the notes were written against 0.22, use the current
`Bound` API) + the `numpy` crate of the matching version + maturin ≥ 1.5. Enable
`abi3-py310` so one wheel serves every Python ≥ 3.10. Rust 1.99 is installed at
`~/.cargo/bin` (not on `PATH` in non-login shells). numpy 1.26 is installed; maturin and
pytest are not (`pip install maturin pytest`).

### 3.1 Module `epsr._core`

Expose exactly these. Names on the Rust side keep the engine's names so the mapping to
`api.rs` is obvious; the Python package renames.

**`class Simulation`** (`#[pyclass]`, holds `epsr_core::Epsr`; `Epsr` is `Send`, so no
`unsendable`).

| method / property | Rust | notes |
|---|---|---|
| `__new__(folder: str, inp_name: str, extra_lines: list[str] = [], verbose: bool = False)` | `Epsr::setup_with_lines` | map errors (§3.3) |
| `run()` | `Epsr::run` | inside `py.allow_threads` |
| `iterate(iter: int) -> dict` | `Epsr::iterate` | inside `allow_threads`; raise `ValueError` if `iter < 1` or `iter > niter` |
| `write_fort35()` | `Epsr::write_fort35` | |
| `summary -> dict` | `sim.summary` | keys = Fortran names (`iter, foldsv, foldsb, energd, pressd, energabs, rfactot, rfacraw, rfacd (list), epcup, epccoeff, sentb, chisq`) |
| `control -> dict` | `sim.ctl` | all `Control` fields except `pressure`; add `pressure` as a nested dict only if trivial |
| `datasets -> list[dict]` | `sim.data` | all `Dataset` fields |
| `n_datasets, n_pairs, nq, nr, n_molecules, n_atoms` | same | ints |
| `atom_types() -> list[str]`, `pair_labels() -> list[tuple[str,str]]`, `atom_labels()`, `molecule_atoms() -> list[tuple[int,int]]` | same | 0-based half-open ranges |
| `q_grid(), r_grid()` | same | 1-D float32 copies |
| `gr_inter(), gr_intra(), sq_inter(), sq_intra()` | per-pair accessors | return **2-D** `(n_pairs, nr|nq)` float32; build with `PyArray2::from_vec2` |
| `data_fq(id), simulated_fq(id)` | same | 1-D float32 |
| `data_errors(id)` | `sim.errs.col(id+1)[..nq]` | 1-D float32 |
| `data_q_range(id) -> (int, int)` | same | |
| `current_gr(), current_sq()` | `sim.wav`, `sim.qwind` | the *current-configuration* histograms, index 1 = inter, 2 = intra; return `(2, n_pairs, n)` float32. Check the F3 layout in `epsr.rs` before indexing. Optional, implement last. |
| `empirical_potential(), reference_potential(), total_potential()` | `api::empirical_potential`, `api::reference_potential`, `st.pot.xavtot.col(ip+1)[..nr]` | 2-D `(n_pairs, nr)` |
| `ep_coefficients() -> ndarray (n_pairs, ncoeffp)` | `api::ep_coefficients` | column 0 is the refinement flag `pcoefp(1)`; keep it here, strip it in Python |
| `set_ep_coefficients(coef: ndarray (n_pairs, ncoeffp))` | write `st.pot.pcoefp[(ic, irefpcof[ip+1])]` for `ic in 1..=ncoeffp`, then `st.generate_ep(false)` | validate shape; **note that several pairs can share one `irefpcof` column** (pairs matched to the same `.pcof` entry) — write in pair order and document that the last write wins |
| `scale_ep(factor: float)` | `st.multiply_pot_coeffs(f)` then `st.generate_ep(false)` | |
| `regenerate_ep()` | `st.generate_ep(false)` | |
| `pot_params() -> dict` | `st.pot.{absolute_energy, potfac, potfacamp, refpotfac, feedback, ncoeffp, pdmin, pdstep, pdmax, psigma2, gaussian}` | read-only snapshot |
| `type_params() -> list[dict]` | `st.mol.{type_, epscor, sigcor, charge, amassa}` indexed `1..=ntype` | read-only; LJ ε, σ, charge, mass per atom type |
| `box_size()`, `atom_positions() -> (n_atoms,3)` | same | |
| `get_parameter(name, dataset=None) -> str | None` | same | |
| `weights_files() -> list[str]` | same | |
| `rho -> float` | `st.dist.rho` | |

**`function run_tool(folder: str, command: str, args: list[str], answers: str, verbose: bool = False)`**
→ `ToolContext::new(folder, args, TextFile::from_string(answers), verbose)`,
`tools::run_tool`. Inside `allow_threads`.

**`function tool_commands() -> list[str]`** → `tools::COMMANDS`.

**`function tool_keywords(name: str) -> dict`** → iterate
`epsr_core::symtab_tools::ALL`, match on `ToolTable.name` (case-insensitive), return
`{"name", "extinp", "prompt", "symbols": [(symbol, type, iflag, default, description), ...]}`.
Read the `ToolTable` definition in `symbols.rs` ~line 580 for the exact field types.

**`VERSION_DATE: str`** → `epsr_core::VERSION_DATE`.

### 3.2 Array conversions

All engine arrays are `f32`. Return **copies** as `float32` numpy arrays (the notes
confirm this is cheap: 21 pairs × 600 points). Never hand out views into `Epsr` memory.
For `F2` columns use the existing `col(j)` and slice to the used length as `api.rs` does.
Write one helper `rows_to_array2(py, Vec<&[f32]>) -> Bound<PyArray2<f32>>`.

### 3.3 Error mapping

Create three Python exception types in the extension (`create_exception!`) and export them:

| Rust | Python |
|---|---|
| `EpsrError::Io(s)` | `EpsrIOError(OSError)` |
| `EpsrError::Input(s)` | `EpsrInputError(ValueError)` |
| `EpsrError::Unsupported(s)` | `EpsrUnsupportedError(NotImplementedError)` |

The message is the Rust `Display` string, which already prefixes `"not supported by the
Rust port: ..."`. The pure-Python `errors.py` re-exports them under a common base
`EpsrError` (a mixin added via `__bases__`-free approach: define the three in Rust with a
common Rust-side base `EpsrError(Exception)` so `except epsr.EpsrError` works).

### 3.4 GIL and threads

* `run`, `iterate`, `run_tool` release the GIL (`py.allow_threads`). `Epsr: Send` is
  asserted by the engine's own test `epsr_is_send`.
* The engine keeps one piece of thread-local state: the active box geometry in
  `neighbour.rs` (`GEOM`), which `Epsr::iterate` re-syncs at its start via
  `sync_lattice`. Consequence: a `Simulation` may be used from any thread between calls,
  and two `Simulation`s on **different folders** may run on two threads at once. Two
  objects on the **same folder** must never run concurrently (they overwrite each other's
  files). Document this; `test_threads.py` checks the first claim.
* Progress messages (`verbose=True`) go to Rust's stdout through `println!`, not through
  Python's `sys.stdout`; they cannot be captured by Python logging without an engine
  change. Default `verbose=False`, say so in the docstring, and do not add a log hook to
  `state.rs`.

## 4. The pure-Python package (`python/epsr/`)

### 4.1 `project.py` — `Project`

Owns a folder and knows the EPSR naming conventions.

```python
class Project:
    folder: Path
    inp_name: str            # e.g. "200DMSO400Water.EPSR.inp" (auto-detected: the unique *.EPSR.inp)
    run_name: str            # "200DMSO400Water.EPSR" (first word of line 1 of the .inp file)

    @classmethod
    def from_template(cls, template: PathLike, work: PathLike, *, overwrite=False) -> "Project"
        # copies regular files only (no subfolders), refuses to overwrite unless asked
    def simulation(self, overrides=None, verbose=False) -> Simulation
    def snapshot(self, label: str) -> Path         # copy the folder to <folder>.snapshots/<label>/
    def output(self, ext: str) -> Path             # ext "g01" -> folder / f"{run_name}.g01"
    def read_table(self, ext: str) -> Table        # see io.py
    def history(self) -> History                   # parsed .erg
    def files(self) -> dict[str, Path]             # inp, ato, pcof, erg, out, data files, weights files
```

The copy helper is mandatory: the notes say loudly that every run overwrites `.ato`,
`.pcof`, `.EPSR.inp`, `.EPSR.inpa` and appends to `.erg`. Put the same warning in the
`Simulation` docstring and the README.

### 4.2 `simulation.py` — `Simulation`, `Summary`, `Dataset`, `Control`

`Simulation` wraps `_core.Simulation`. Constructor signature:

```python
Simulation(folder, inp_name, *, overrides: dict | list[str] | None = None, verbose=False)
```

`overrides` accepts the raw `.EPSR.inp` line form (`["ntimes 10", "data 2", "tweak 1.05"]`)
or a dict that is rendered to lines: scalar → `"key value"`, list/tuple → values joined
by spaces (`{"potfac": [1.0, 1.0]}` → `"potfac 1.0 1.0"`), and a nested
`{"data": {2: {"tweak": 1.05}}}` → `["data 2", "tweak 1.05"]`. Keywords are those of
`rust/epsr-core/src/symtab.rs` (engine) and the `.pcof` header (`ereq`, `feedback`,
`potfac`, `refpotfac`, … — these also go through the same symbol table; verify with
`get_parameter` in a test). Remember overrides are **sticky**: the engine writes them back
to `.EPSR.inp` / `.pcof` at the end of the run.

Methods and properties (descriptive names; keep the Fortran name in each docstring):

* `run()`, `iterate(i)`, `iterations()` — a generator yielding `Summary` for
  `i = 1..niter`; raises if the object has already completed its run.
* `summary: Summary` (frozen dataclass):
  `iteration (iter)`, `energy_inter (foldsv)`, `energy_intra (foldsb)`, `energy (energd)`,
  `pressure (pressd)`, `ep_amplitude (energabs)`, `r_factor (rfactot)`,
  `r_factor_raw (rfacraw)`, `r_factor_per_dataset (rfacd)`, `ep_up (epcup)`,
  `ep_coefficient (epccoeff)`, `gr_roughness (sentb)`, `chisq`. Provide `.as_dict()` and
  `.fortran` (the raw dict).
* `control: Control`, `datasets: list[Dataset]` — dataclasses mirroring the dicts.
* `q`, `r`, `nq`, `nr`, `n_pairs`, `n_datasets`, `n_atoms`, `n_molecules`, `box_size`,
  `density (rho)`, `pair_labels`, `atom_types`, `atom_labels`, `molecule_atoms`,
  `positions`.
* `pair_index(a, b)` — order-insensitive lookup into `pair_labels`.
* `gr_inter(pair=None)`, `gr_intra`, `sq_inter`, `sq_intra` — 2-D array, or the 1-D row
  when `pair` (index or label tuple) is given.
* `data_fq(i)`, `simulated_fq(i)`, `data_errors(i)`, `data_q_range(i)`, `fit(i)` → a
  small `Fit` dataclass `(q, data, errors, simulated, q_range)` restricted to the fitted
  range, for plotting.
* `potential: Potential` (§4.3).
* `folder`, `inp_name`, `run_name`, `completed: bool`.

### 4.3 `potential.py` — `Potential`

A handle onto a live `Simulation`; no state of its own. This is the "perturbation" API.

```python
class Potential:
    # read
    def empirical(self, pair=None) -> ndarray       # xavs      (n_pairs, nr) or (nr,)
    def reference(self, pair=None) -> ndarray       # ecoul
    def total(self, pair=None) -> ndarray           # xavtot
    def ep_coefficients(self) -> ndarray            # (n_pairs, ncoeffp-1): flag column stripped
    def ep_flags(self) -> ndarray                   # (n_pairs,) the pcoefp(1) refinement flags
    params: PotParams                               # ereq (absolute_energy), potfac, refpotfac, feedback, ncoeffp, pdstep, pdmax, psigma2, gaussian
    type_params: list[TypeParams]                   # label, epsilon, sigma, charge, mass
    r: ndarray
    # write (in memory; persisted to .pcof when the run finishes)
    def set_ep_coefficients(self, coef)             # (n_pairs, ncoeffp-1); reattaches flags
    def scale_ep(self, factor)                      # multiply_pot_coeffs
    def zero_ep(self)                               # scale_ep(0.0); equivalent to `ireset 2`
    def perturb_ep(self, pair, delta_coef)          # add to one pair's coefficients
    def pair_index(self, a, b)
```

Design notes for the implementer:

* The EP in r-space is **derived** from the coefficients through the basis in
  `generate_ep` (Poisson or Gaussian functions with spacing `pdstep`). Users therefore
  perturb coefficients, not the r-space table. Do not offer a `set_empirical(r_table)`;
  if asked for, it would require a least-squares projection onto the basis and belongs to a
  later phase.
* The reference potential's parameters (LJ ε/σ, charge) live in the `.ato` file. Expose
  them read-only in phase 3. Writing them (edit `st.mol.epscor/sigcor/charge`, call
  `st.generate_rp()`) is feasible and additive but is **phase 6 / optional**; it also
  changes what `moloutc` writes back to `.ato`.
* The standard EPSR way to perturb the potential globally is the `ereq`/`potfac`/
  `feedback`/`refpotfac` keywords and per-dataset `efilereq`/`tweak`; these go through
  `overrides` at construction. `Potential.params` only *reads* them.
* After any in-memory change the engine must regenerate the tables before
  `iterate` (`fmolec1` uses `xavtot`). The `_core` setters do that; never expose a raw
  `pcoefp` write without the regenerate.

### 4.4 `refine.py` — `refine()`, `History`

```python
def refine(project, n_runs, *, overrides=None, callback=None, until=None, verbose=False,
           snapshot_every=None) -> History
```

Faithful loop from the notes §3: each run constructs a fresh `Simulation` (so
`.EPSR.inp` and the accumulators are re-read), calls `run()`, appends `summary` to the
history, then calls `callback(k, sim)`; stop early if `callback` returns `True` or
`until(history)` is true. `overrides` are applied to every run (they are sticky anyway;
applying them each time is harmless and explicit). `snapshot_every=n` calls
`project.snapshot(f"run{k:04d}")`. Show a `tqdm` progress bar when `tqdm` is
importable and `progress=True` (default: `True` when stdout is a TTY).

`History`: columns as numpy arrays (`run`, `iteration`, `energy`, `pressure`,
`r_factor`, `r_factor_per_dataset` (2-D), `ep_amplitude`, `chisq`, …), `append(summary)`,
`to_csv(path)` (via `numpy.savetxt` with a header), `to_dataframe()` (pandas, optional
extra; raise a clear `ImportError` otherwise), `from_erg(path)` (parse the `.erg` file: 10 columns written by
`epsr.rs` ~line 845 as `energd, pressd, rfacraw, energabs, rho, epcup, epccoeff, aparab,
bparab, cparab`, one line per iteration), `__len__`, `__getitem__`.

Optional (phase 5): `refine_many(projects, n_runs, workers=N)` using
`ThreadPoolExecutor`, valid because the GIL is released and the folders differ.

### 4.5 `io.py`

Read-only parsers for EPSR text outputs; numpy only.

* `read_table(path) -> Table(x, columns: dict[str, ndarray], errors: dict[str, ndarray], labels)`
  for every `.?01` column file; parse the numeric block with `numpy.loadtxt(comments="#")`
  and offer `Table.to_dataframe()` (pandas optional). Line 1 is a `#` header with the x name and the pair /
  dataset labels (see `EPSR-test/200DMSO400Water.EPSR.g01`); numeric rows are x followed
  by value/error pairs. Check `rust/epsr-core/src/output.rs` for which extensions have
  error columns and the table of extensions.
* `read_erg(path) -> ndarray (n, 10)`.
* `read_keyword_file(path) -> dict` for `.EPSR.inp`, `.pcof` headers, `*.NWTS.dat`
  etc.: keyword, values (split on whitespace, description ignored). Needed only for
  inspection and for `write_keyword_file` round trips; the engine remains the writer.

### 4.6 `tools.py`

```python
def commands() -> list[str]
def run(folder, command, args=(), answers: str | Sequence[str] = "", verbose=False) -> None
def keywords(command) -> ToolKeywords           # from _core.tool_keywords
def write_keyword_file(path, command, values: dict, title=None) -> Path
```

`answers` is what the Fortran program reads from standard input; a sequence is joined
with newlines and a trailing newline is added. `write_keyword_file` lays the file out as
the Fortran does: title line, then `keyword<pad>value(s)<pad>description` lines in table
order with defaults for anything not given, then `q`. Copy the column layout from an
existing file (e.g. `EPSR-test/DMSODWaterD.NWTS.dat`) and from `symbols.rs::save_inp_file`;
the tool rewrites the file itself afterwards, so exact spacing is not critical, but the
`component n` / `data n` block structure is — look at how `dataprompt` is used in the
tables. Tools with a `component`/`data` block take a nested dict like the engine overrides.

Typed convenience wrappers for the most common tools are welcome but not required:
`nwts(folder, dat_file)`, `partials(folder, dat_file)`, `coord(...)`, `writexyz(folder, ato, out)`.

### 4.7 `plot.py`

Lazy-import matplotlib. Three functions, each taking an optional `ax`:
`fits(sim, datasets=None)` (data vs simulated F(Q), offset per dataset),
`partials(sim, pairs=None, kind="gr_inter")`, `potentials(sim, pairs=None)`, and
`history(hist)`. Nothing fancy; this replaces the unported PGPLOT front ends.

## 5. Build, install, test commands

`uv` is not installed on this machine yet; install it first with the official standalone
installer (`curl -LsSf https://astral.sh/uv/install.sh | sh`) and put `~/.local/bin` on
`PATH` alongside `~/.cargo/bin`.

```bash
export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
cd python
uv sync --all-extras                        # creates .venv, resolves uv.lock, builds the extension
                                            # (maturin is the build backend, so `uv sync` compiles ../rust/epsr-core)
uv run maturin develop --release --uv       # fast rebuild after editing src/lib.rs (optimised)
uv run pytest -q                            # reference test ≈ 3 s per run
uv run python examples/01_single_run.py
uv build                                    # wheel + sdist via maturin, into python/dist/
cargo test --release --manifest-path ../rust/Cargo.toml   # must still pass, untouched
```

Notes for the `pyproject.toml`:

* `[build-system] requires = ["maturin>=1.5,<2"]`, `build-backend = "maturin"`.
* `[tool.maturin] module-name = "epsr._core"`, `features = ["pyo3/extension-module"]`,
  and `python-source` left unset (the `epsr/` package sits next to `pyproject.toml`).
* `[tool.uv] cache-keys = [{file = "pyproject.toml"}, {file = "Cargo.toml"}, {file = "src/**/*.rs"}, {file = "../rust/epsr-core/src/**/*.rs"}]`
  so `uv sync` rebuilds the extension when the Rust engine changes.
* Set `[tool.uv] package = true` (the default for a project with a build backend) so
  `uv sync` installs `epsr` itself into the environment; a debug build is fine for
  development, but run the reference test with a release build
  (`uv run maturin develop --release --uv`) because the engine is slow unoptimised.
* Add `uv run` task equivalents to `README.md`; no Makefile needed.

Add `python/target/`, `python/.venv/`, `python/dist/`, `python/**/__pycache__/`,
`*.so`, `work/` to a `python/.gitignore`. Do **not** ignore `uv.lock`.

## 6. Tests (pytest) and acceptance criteria

1. **`test_reference.py`** — copy `EPSR-test` to a tmp folder, `Simulation(...).run()`,
   mask and FNV-1a-hash every file listed in `rust/epsr-core/tests/reference_hashes.txt`
   (reuse `mask`/`fnv1a64` from `rust/tools/make_reference_hashes.py` by importing the
   file, or copy the two functions), compare. Also assert `n_pairs == 21`,
   `n_datasets == 3`, `len(q) == 600`, `abs(summary.r_factor - 1.1478e-4) < 1e-7`.
   Mark `skipif` not Linux/glibc for the hash part (libm caveat in the notes).
2. **`test_simulation.py`** — shapes `(21, nr)` / `(21, 600)`, dtype float32,
   `pair_index` symmetry, `overrides` dict → `get_parameter("tweak", 1) == "1.05"`,
   `iterate(niter+1)` raises, missing folder raises `EpsrIOError`, bad keyword raises
   `EpsrInputError`, `iterations()` yields `niter` summaries and the last equals `summary`.
3. **`test_potential.py`** — `set_ep_coefficients(ep_coefficients())` leaves
   `empirical()` bit-identical; `scale_ep(2.0)` doubles `ep_coefficients()` and
   `empirical()` (basis is linear in the coefficients); `zero_ep()` gives an all-zero
   `empirical()` and `total() == reference()`; after `run()` the `.pcof` file reflects the
   edit (reopen in a new `Simulation`, compare coefficients).
4. **`test_refine.py`** — 3 runs on a template copy: history length 3, `.erg` gained 3
   lines, `callback` early stop works, `until` works, snapshots exist.
5. **`test_tools.py`** — run `nwts` on each `DMSODWater*.NWTS.dat` in a copy of
   `EPSR-test`; the produced `*.NWTStot.wts` must be byte-identical to the shipped ones
   (this is what `test_problems/200DMSO400Water/EPSR-test-ex20/README.md` states).
   `tool_commands()` has 35 entries. `write_keyword_file` + `read_keyword_file` round trip.
6. **`test_threads.py`** — two template copies refined once each on two threads give the
   same summaries as sequentially; while one `run()` executes, the main thread can do work
   (GIL released) — assert elapsed-time overlap or use a counter thread.
7. **Determinism** — two fresh copies of the template give identical `summary` and
   identical `.pcof` bytes after one run.

The four example scripts must run end to end from `python/` against `../EPSR-test`.

## 7. Implementation order

| phase | deliverable | done when |
|---|---|---|
| 0 | `Cargo.toml`, `pyproject.toml`, `src/lib.rs` with `Simulation.__new__/run/summary`, `maturin develop` succeeds | `import epsr._core` works; one run reproduces `r_factor` |
| 1 | all `_core.Simulation` accessors of §3.1 except `set_*`/`current_*`; `_core.pyi`; `errors.py`; `simulation.py` dataclasses and array properties | `test_reference.py`, `test_simulation.py` pass |
| 2 | `project.py`, `refine.py`, `io.py`, examples 01–02 | `test_refine.py` passes; `History.from_erg` matches `History` built from summaries for `rfacraw`/`energd` |
| 3 | `_core` setters, `potential.py`, example 03 | `test_potential.py` passes |
| 4 | `tools.py`, `_core.run_tool`/`tool_keywords`, example 04 | `test_tools.py` passes |
| 5 | `plot.py`, `README.md`, `test_threads.py`, `refine_many` | all tests green; `cargo test` in `rust/` still green and `git status rust/` clean |
| 6 (optional) | `current_gr/current_sq`, writable reference-potential parameters, `set_empirical` via basis projection | only if time permits; keep additive |

## 8. Decisions already made (do not reopen)

* Package/import name `epsr`; extension module `epsr._core`; distribution name `epsr-py`.
* PyO3 + maturin + the `numpy` crate; arrays returned as float32 copies; 2-D
  `(n_pairs, n)` for pair-resolved quantities.
* One `Simulation` per EPSR invocation; no `iterate` beyond `niter`; new object per run.
* Overrides are passed at construction (`setup_with_lines`), not applied to a live object:
  `Control` is parsed from the symbol table only during `setup`, so a later
  `apply_input_lines` would desynchronise `ctl` from the table.
* Potential perturbation acts on EP coefficients in memory and is persisted by the engine's
  own `.pcof` writer; Python never writes `.pcof` or `.ato`.
* Rust engine is not modified; any unavoidable addition goes at the end of `api.rs` and is
  recorded in `python/CHANGES_TO_RUST.md`.
* `uv` is the only package/environment manager; `uv.lock` is committed; maturin is the
  build backend invoked through `uv`.
* Existing libraries over home-grown code: numpy required; pandas, matplotlib, tqdm as
  optional extras; pytest, ruff, mypy in the dev group. No other runtime dependencies
  without a feature in this plan that needs them.

## 9. Known caveats to document in `README.md`

* Files are modified in place; use `Project.from_template`.
* Bit-for-bit reproducibility of the Fortran holds on Linux/glibc; elsewhere results agree
  closely but not necessarily to the last bit.
* Inputs the Fortran itself cannot handle raise `EpsrUnsupportedError`; the list is in
  `rust/docs/FORTRAN_MAP.md` ("Still unsupported").
* `verbose=True` prints from Rust directly to the process stdout.
* The plotting programs of the EPSR suite are not ported; use `epsr.plot`.
* `makeato` reproduces a Fortran crash (empty `.ato`, error exit) by design, see
  `rust/docs/TOOLS.md`.
