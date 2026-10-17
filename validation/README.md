# Validation of the Rust port on experimental data: water and amorphous silica

This folder checks the Rust port (`rust/`) against the original Fortran EPSR on two
systems fitted to real neutron data. Both are taken from the Dissolve examples in
`experimental_data/`:

| case | system | data |
|------|--------|------|
| [water/](water/) | 1000 H₂O, 0.1 atoms/Å³, 300 K | SANDALS: H₂O, D₂O, HDO (50:50) — `SLS18498/18500/18502.mint01` |
| [silica_bulk/](silica_bulk/) | 1500 Si + 3000 O, 0.0664 atoms/Å³, 298 K | NIMROD: `SiO2_NIMROD.mint01` |

Each case is built from scratch with the EPSR setup programs (molecule → box →
weights → `.EPSR.inp`). It is then equilibrated and refined with the empirical
potential. Every step is run by the Fortran and by the Rust program on separate
copies of the folder, and **all files of the two folders are compared after every
step**.

For how these calculations map onto the Rust API (for the Python wrapper), see
[PYTHON_WRAPPER_WORKFLOW.md](PYTHON_WRAPPER_WORKFLOW.md).

## Results

RESULTS_PLACEHOLDER

## How to reproduce

```bash
# 1. build the Fortran reference (gfortran) and the Rust port
rust/tools/build_fortran_reference.sh ~/epsr-work/fortran-build
(cd rust && cargo build --release)

# 2. run a case: Fortran vs Rust CLI, every file compared after every step
validation/run_validation.sh validation/water       ~/epsr-work/fortran-build ~/epsr-work/val-water
validation/run_validation.sh validation/silica_bulk ~/epsr-work/fortran-build ~/epsr-work/val-silica

# 3. optional: run the same steps through the epsr-core library API, then compare with the Fortran folder
(cd validation/api_driver && cargo build --release)
validation/api_driver/target/release/epsr-api-driver validation/water ~/epsr-work/api-water
python3 rust/tools/compare_outputs.py ~/epsr-work/val-water/fortran ~/epsr-work/api-water --quiet
```

Each case takes about 25 minutes on one core (Fortran and Rust run one after the
other). Put the work folders on a real disk: `/tmp` on this machine is a small tmpfs.
Don't edit `run_validation.sh` while it is running; bash reads scripts as it goes.

`WORK_DIR/summary.txt` has one line per step. `WORK_DIR/logs/` holds the console
output of every program and the full comparison after every step.

## Files

| file | purpose |
|------|---------|
| `run_validation.sh` | driver: runs `steps.txt` in a Fortran and a Rust copy of the case, compares after each step |
| `set_keyword.py` | sets a keyword in an EPSR keyword file (the edit a user makes between runs) |
| `<case>/steps.txt` | the calculation, step by step (`run`, `set`, `delete`) |
| `<case>/*.mol`, `*.NWTS.dat`, `*.EPSR.inp`, `mixato.stdin`, `*.mint01` | inputs |
| `api_driver/` | standalone crate: the same steps through `epsr-core` (`run_tool`, `Epsr::setup`/`run`) |
| `PYTHON_WRAPPER_WORKFLOW.md` | the calculation documented call by call, with pitfalls |

## From the Dissolve inputs to EPSR inputs

The files in `experimental_data/` are Dissolve inputs, not EPSR inputs. The
following was carried over:

* **Water.** Geometry (O–H 0.976 Å, H–O–H 107.134°), LJ (ε 0.6503 kJ/mol,
  σ 3.165492 Å on O) and charges (−0.82/+0.41 e). Also the population (1000),
  density (0.1 atoms/Å³), temperature (300 K), the three isotopic samples (H₂O;
  D₂O; HDO 50:50 with exchangeable H) and `EReq 3`.
* **Silica.** LJ parameters and charges (Si: ε 0.175, σ 1.03, +2; O: ε 0.165,
  σ 3.5, −1), populations 1500/3000, density 0.0664 atoms/Å³, 298 K, Q_min 0.1 Å⁻¹
  (Dissolve `XMin 0.1`) and `EReq 50`.

What differs from the Dissolve setups:

* **The silica box is cubic.** Dissolve uses a 90/90/120° cell, but `mixato`
  builds cubic boxes. The density and atom count are the same.
* **EPSR's own algorithm and defaults are used**, not Dissolve's (MC with
  `ntimes 5`, `feedback 0.9`, Poisson EP functions, `nq 600`, `qstep 0.05`,
  `fwhmq 0.02`, ...). These values come from the EPSR test case. Dissolve settings
  without an EPSR counterpart (MD moves, `Weighting 0.1`, broadening, window
  functions) are not reproduced.
* **The schedule is short.** It runs 100 iterations on the reference potential,
  then 2 × 100 with the empirical potential. That is enough to compare the two codes
  over 300 coupled iterations, but not a converged production refinement.
* `rshmin` (background-subtraction radius) is 0.7 Å for water (as in the EPSR test
  case) and 1.2 Å for silica (below the 1.6 Å Si–O distance).
