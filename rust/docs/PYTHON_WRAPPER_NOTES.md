# Notes for building the Python wrapper

These notes are for the developer (or agent) who will wrap `epsr-core` for
Python.  Read `ARCHITECTURE.md` first for the engine's structure.

## 1. What a Python user needs to do with EPSR

The essential EPSR workflow, which the wrapper should make easy:

1. **Run refinement iterations** on a prepared simulation folder
   (`.EPSR.inp`, `.ato`, `.pcof`, weights `.wts`, data files), usually many
   times in a row, with optional parameter changes between runs
   (`potfac`/feedback, `ereq`, `ntimes`, ...).
2. **Monitor convergence**: energy per molecule, pressure, R-factor per
   dataset, the EP amplitude (`energabs`), the R-space coefficient.
3. **Inspect results** as arrays: data vs. simulated F(Q) per dataset,
   partial g(r) and S(Q) per atom pair, the empirical and reference
   potentials, atom coordinates.
4. **Plot** them (the Fortran PGPLOT/gnuplot tools are not ported; matplotlib
   on the arrays replaces them).

## 2. The Rust API to wrap

Everything is on `epsr_core::Epsr` (see `src/epsr.rs` and `src/api.rs`):

| Rust | meaning |
|------|---------|
| `Epsr::setup(dir, inp_name, verbose) -> Result<Epsr>` | read all inputs (one "invocation" of the Fortran `epsr`) |
| `Epsr::setup_with_lines(dir, inp_name, &[String], verbose)` | same, then apply extra `.EPSR.inp` lines, e.g. `"ntimes 10"`, `"data 2"`, `"tweak 1.05"` |
| `sim.iterate(iter) -> Result<&Summary>` | one EPSR iteration; output files are written when `iter == sim.ctl.niter` |
| `sim.run() -> Result<()>` | `iterate(1..=niter)` + `fort.35` — exactly one Fortran `epsr` invocation |
| `sim.summary` (`Summary`) | `foldsv` (inter energy/molecule, running), `foldsb` (intra), `energd`, `pressd`, `energabs`, `rfactot`, `rfacraw`, `rfacd[dataset]`, `epcup`, `epccoeff`, `sentb`, `chisq`, `iter` |
| `sim.ctl` (`Control`) | run parameters read from the input (`ntimes`, `niter`, `nq`, `qstep`, ...) |
| `q_grid()`, `r_grid()` | Q (Å⁻¹) and r (Å) grids |
| `data_fq(id)`, `simulated_fq(id)`, `data_q_range(id)` | measured and simulated total F(Q) of dataset `id`, fitted range |
| `gr_inter(ip)`, `gr_intra(ip)`, `sq_inter(ip)`, `sq_intra(ip)` | running-average partials of pair `ip` |
| `empirical_potential(ip)`, `reference_potential(ip)`, `ep_coefficients(ip)` | potentials (kJ/mol) on the r grid; EP coefficients |
| `pair_labels()`, `atom_types()`, `n_pairs()`, `n_datasets()` | labelling |
| `atom_positions()` (Å), `atom_labels()`, `molecule_atoms()`, `box_size()` | configuration |
| `get_parameter(name, dataset)` | current value of an input keyword |
| `epsr_core::EpsrError` | `Io`, `Input`, `Unsupported` |

All accessor slices are 0-based and of the used length.  Pair index `ip`
follows `pair_labels()`; dataset index `id` follows the `data n` blocks.

Lower-level state is public (`sim.st.mol`, `sim.st.dist`, `sim.st.pot`, ...)
but uses Fortran layout (1-based `F1`/`F2`, fixed leading dimensions) —
prefer adding accessors to `api.rs` over exposing it directly.

## 3. Important behaviour to preserve / expose

* **Files are modified in place.**  Every `run()` overwrites `.ato`, `.pcof`,
  `.EPSR.inp`, `.EPSR.inpa`, appends to `.erg`, writes `delfit.dat` and
  `fort.35` in the simulation folder.  This *is* how EPSR carries state
  between runs (accumulated averages, EP, RNG state).  Offer a helper that
  copies a template folder first, and document the side effects loudly.
* **A "run" is setup + run.**  The Fortran program is normally invoked
  repeatedly by a script (EPSR's `runepsr`/shell loop).  The faithful Python
  equivalent is

  ```python
  for k in range(n_runs):
      sim = epsr.Simulation(folder, "200DMSO400Water.EPSR.inp")  # Epsr::setup
      sim.run()                                                  # Epsr::run
      history.append(sim.summary)
  ```

  Do not call `iterate` past `ctl.niter` on the same object: the Fortran
  never does that, and output files are only written on the last iteration.
  (A `run_many(n)` convenience doing the loop above in Rust is cheap to add.)
* **Determinism.**  Results depend only on the folder contents (the RNG state
  is stored in the `.ato` file).  The same folder always gives the same
  numbers, bit for bit, on any platform with IEEE floats and glibc-compatible
  libm (`sinf`, `cosf`, `expf`, `logf`, `powf`, `atanf`, `sinhf` and double
  `sin`, `cos`, `log`).  On other
  platforms (macOS, Windows) results should agree closely but not
  necessarily bit-for-bit, because those libm implementations may round
  differently.
* **Unsupported features** raise `EpsrError::Unsupported` — map it to
  `NotImplementedError` with the message, so users know which EPSR feature
  their input needs.
* **Paths.**  `dir` is the folder; `inp_name` a file name inside it.  The run
  name is taken from the first line of the `.EPSR.inp` file.

## 4. Suggested implementation

A separate crate `rust/epsr-py` using **PyO3** + **maturin** and the
**numpy** crate:

```toml
# rust/epsr-py/Cargo.toml
[package]
name = "epsr-py"
version = "0.1.0"
edition = "2021"
[lib]
name = "epsr"            # Python module name
crate-type = ["cdylib"]
[dependencies]
epsr-core = { path = "../epsr-core" }
pyo3 = { version = "0.22", features = ["extension-module"] }
numpy = "0.22"
```

Keep it **out of the workspace's default members** (or in its own
workspace), so that `cargo build`/`cargo test` of the engine does not
require Python headers.

Sketch (PyO3 0.22 "Bound" API; adjust to the version used):

```rust
use numpy::{PyArray1, PyArray2, ToPyArray};
use pyo3::{exceptions::*, prelude::*};

fn to_py(e: epsr_core::EpsrError) -> PyErr {
    match e {
        epsr_core::EpsrError::Io(s) => PyIOError::new_err(s),
        epsr_core::EpsrError::Input(s) => PyValueError::new_err(s),
        epsr_core::EpsrError::Unsupported(s) => PyNotImplementedError::new_err(s),
    }
}

#[pyclass]
struct Simulation { inner: epsr_core::Epsr }

#[pymethods]
impl Simulation {
    #[new]
    #[pyo3(signature = (folder, inp_name, overrides=None, verbose=false))]
    fn new(folder: std::path::PathBuf, inp_name: &str, overrides: Option<Vec<String>>, verbose: bool) -> PyResult<Self> {
        let lines = overrides.unwrap_or_default();
        epsr_core::Epsr::setup_with_lines(&folder, inp_name, &lines, verbose)
            .map(|inner| Simulation { inner }).map_err(to_py)
    }
    /// One full EPSR invocation (niter iterations); releases the GIL.
    fn run(&mut self, py: Python<'_>) -> PyResult<()> {
        py.allow_threads(|| self.inner.run()).map_err(to_py)
    }
    #[getter]
    fn q<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> { self.inner.q_grid().to_pyarray_bound(py) }
    #[getter]
    fn r<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> { self.inner.r_grid().to_pyarray_bound(py) }
    /// (n_pairs, n_r) array of intermolecular partial g(r)
    fn gr_inter<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f32>> {
        let rows: Vec<Vec<f32>> = (0..self.inner.n_pairs()).map(|ip| self.inner.gr_inter(ip).to_vec()).collect();
        PyArray2::from_vec2_bound(py, &rows).unwrap()
    }
    // ... summary as a dict, pair_labels as list of tuples, etc.
}

#[pymodule]
fn epsr(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Simulation>()?;
    Ok(())
}
```

Build/test: `pip install maturin numpy && maturin develop --release -m rust/epsr-py/Cargo.toml`.

`Epsr` contains only owned data (no `Rc`/raw pointers), so it is `Send`
(checked by the test `epsr_is_send` in `tests/test_case.rs`) and can be
stored in a `#[pyclass]` and used inside `allow_threads`.

Returning copies (`to_pyarray`) is simplest and cheap (e.g. 21 pairs × 600
points for the test case).  Arrays are `float32`, which is the precision EPSR computes in.

## 5. Suggested Python API

```python
import epsr
sim = epsr.Simulation("work/", "200DMSO400Water.EPSR.inp", overrides=["ntimes 10"])
sim.run()
sim.summary            # dict: energy, pressure, r_factor, r_factor_per_dataset, ep_amplitude, ...
sim.q, sim.r           # grids
sim.data_fq(0), sim.simulated_fq(0)
sim.pair_labels        # [("S","S"), ("S","C"), ...]
sim.gr_inter()         # (n_pairs, n_r)
sim.empirical_potential()  # (n_pairs, n_r)
sim.positions          # (n_atoms, 3) Å
```

plus helpers in pure Python:
* `epsr.prepare(template_dir, work_dir)` — copy a folder before running;
* `epsr.read_table(path)` — read any `.?01` output file
  (`numpy.loadtxt(path, comments="#")`; columns are x, then value/error
  pairs per pair or dataset; the header line gives the labels);
* `epsr.refine(folder, inp, n_runs, callback=None)` — the loop of §3.

## 6. Testing the wrapper

* Port `epsr-core/tests/test_case.rs` to pytest: copy `EPSR-test` to a temp
  folder, run once, hash the masked outputs (`tools/make_reference_hashes.py`
  contains the masking + FNV-1a hash in Python), compare with
  `epsr-core/tests/reference_hashes.txt`.  This proves the wrapper drives
  the engine exactly like the CLI.
* Check `sim.summary["r_factor"] ≈ 1.1478e-4` after one run of the test case.
* `tools/chain_test.sh` remains the way to compare against the Fortran.

## 7. The auxiliary programs (setup and analysis tools)

All 35 setup and analysis programs are ported. [TOOLS.md](TOOLS.md) lists them with their test cases.
They share a single entry point:

```rust
use epsr_core::fortran::listread::TextFile;
use epsr_core::tools::{run_tool, ToolContext, COMMANDS};

let stdin = TextFile::from_string("answers\nto\nprompts\n");   // what the Fortran reads from unit 5
let mut ctx = ToolContext::new(Path::new(folder), vec!["DMSOWater.NWTS.dat".into()], stdin, false);
run_tool("nwts", &mut ctx)?;                                      // writes its files into `folder`
```

* Wrap this as a single function, `epsr.run_tool(folder, command, args=[], answers="")`.
  `COMMANDS` gives the valid names.
* Keyword-file programs (`nwts`, `xwts`, `partials`, `coord`, `sharm`, ...) take the `.dat` file name as their
  argument. They rewrite it with the defaults filled in, as the Fortran does. A Python helper can create the `.dat` file
  from a dict, using the keyword tables in `symtab_tools.rs`.
* Interactive programs (`makeato`, `mixato`, `fmole`, `addato`, `changeato`, `makelattice`, `readxyz`, ...) read
  their answers from `stdin`. The `stdin` files of the corresponding EPSR-test-ex cases show the expected order of answers.
* Outputs are text files in the folder. Parse them as described in §5.
* The tools are not thread-safe with respect to each other in one folder, because they write files. They hold no
  global state apart from the thread-local lattice of `neighbour.rs`, so different folders can run in
  parallel threads.

## 8. Known gaps relevant to Python users

* The plotting programs (PGPLOT/gnuplot drivers) are not ported. Plot the arrays with matplotlib.
* Inputs on which the Fortran itself fails raise `Unsupported` (list in `FORTRAN_MAP.md`).
