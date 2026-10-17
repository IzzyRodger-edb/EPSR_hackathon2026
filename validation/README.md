# Validation of the Rust port on experimental data: water, amorphous silica and Cu-BTC

This folder checks the Rust port (`rust/`) against the original Fortran EPSR on three
systems fitted to real neutron data. All are taken from the Dissolve examples in
`experimental_data/`:

| case | system | data |
|------|--------|------|
| [water/](water/) | 1000 H₂O, 0.1 atoms/Å³, 300 K | SANDALS: H₂O, D₂O, HDO (50:50) — `SLS18498/18500/18502.mint01` |
| [silica_bulk/](silica_bulk/) | 1500 Si + 3000 O, 0.0664 atoms/Å³, 298 K | NIMROD: `SiO2_NIMROD.mint01` |
| [cu-btc/](cu-btc/) | Cu-BTC (HKUST-1) metal-organic framework, one 624-atom unit cell, 300 K | NIMROD: `Empty-CuBTC.mint01` |

The liquid cases are built from scratch with the EPSR setup programs (molecule →
box → weights → `.EPSR.inp`). Cu-BTC is a crystal: its unit cell is converted to a
PDB file and read with `readpdb`, the force field is applied with `changeato`, and
EPSR calculates its Bragg peaks (see [cu-btc/README.md](cu-btc/README.md)). Each
case is then equilibrated and refined with the empirical potential. Every step is run by the Fortran and by the Rust program on separate
copies of the folder, and **all files of the two folders are compared after every
step**.

For how these calculations map onto the Rust API (for the Python wrapper), see
[PYTHON_WRAPPER_WORKFLOW.md](PYTHON_WRAPPER_WORKFLOW.md).

## Results

**All three cases: every file is byte-identical between Fortran and Rust after every
step.** The only masking is the date on line 2 of `.EPSR.out` and the timing
columns of `fort.35`.

| case | steps | EPSR iterations | files compared after the last step | result |
|------|-------|-----------------|------------------------------------|--------|
| water | 15 (7 setup, 4 EPSR runs, 4 edits) | 300 | 39 | ALL STEPS IDENTICAL |
| silica_bulk | 12 (5 setup, 4 EPSR runs, 3 edits) | 300 | 36 | ALL STEPS IDENTICAL |
| cu-btc | 10 (3 setup, 4 EPSR runs, 3 edits) | 300 | 36 | ALL STEPS IDENTICAL |

That covers the setup programs (`makemole`, `mixato`, `randomise`, `fmole`, `nwts`,
and for Cu-BTC `readpdb` and `changeato`), a scratch start without a `.pcof`, 99
iterations on the reference potential, and two 100-iteration runs with the
empirical potential, the second restarting from the first's files. Cu-BTC also
exercises the Bragg-peak calculation (`.bragg`, `.inten`) on every iteration. MC trajectories are chaotic, so a one-ulp difference anywhere would
have grown into visible differences within a few iterations.

The same steps run through the **library API** (`api_driver/`, i.e. `run_tool` and
`Epsr::setup` + `run`, as a Python wrapper would call them) also give a folder
identical to the Fortran one. The only extra files are the driver's own
`api_fit_<n>.txt`.

State after the last run (from `sim.summary`; identical to the Fortran `.erg`):

| case | energy / molecule [kJ/mol] | EPSR pressure | EP amplitude `energabs` | R-factor total | R-factor per data set |
|------|------------------|----------|-----------|---------|-------------|
| water | −44.89 | 3.87 | 3.00 (= `ereq`) | 2.22×10⁻⁴ | D₂O 2.13×10⁻⁴, HDO 2.30×10⁻⁴, H₂O 2.24×10⁻⁴ |
| silica | −845.3 | 108.7 | 23.7 (rising towards `ereq` 50) | 2.30×10⁻⁴ | 2.30×10⁻⁴ |
| cu-btc | −558.0 | 2.6 | 10.0 (= `ereq`) | 1.64×10⁻² | 1.64×10⁻² |

R-factor history: water fell from 2.1×10⁻² (random start) to 6.1×10⁻⁴ after the
reference-potential stage, then to 2.2×10⁻⁴ with the EP. Silica went from
3.5×10⁻³ to 1.0×10⁻³, then to 2.3×10⁻⁴. The water energy of −45 kJ/mol is typical
for this kind of water model. The silica schedule is too short for a converged
refinement: the EP amplitude was still rising and the pressure is high. Cu-BTC went
from 3.1×10⁻² to 2.3×10⁻², then to 1.6×10⁻²; its fit stays far poorer than the
liquids' because the Cu–O coordination bonds are only Lennard-Jones and Coulomb
interactions between molecules (see its README).

Speed: Rust took 1.2–1.6× as long as the Fortran. Each code ran single-threaded and the two codes alternated, but the two cases sometimes ran at the same time on this multi-core machine.
Water: 1.32 vs 1.62 s per iteration on the reference potential, 1.55 vs 2.4 s with
the EP. Silica: 2.05 vs 2.7 s, then 2.1 vs 3.1–3.2 s. Cu-BTC (624 atoms with Bragg peaks):
0.35 vs 0.42 s on the reference potential, 0.40 vs 0.78 s with the EP.

### Agreement plots

`plot_gr_agreement.py` overlays the intermolecular partial g(r) after the full
schedule as written by the Fortran, by the Rust CLI and by the Python wrapper,
with the differences from the Fortran underneath. They are zero for every pair
of every case: the three `.g01` files are the same file.

| case | plot |
|---|---|
| water | [water/gr_agreement.png](water/gr_agreement.png) |
| silica_bulk | [silica_bulk/gr_agreement.png](silica_bulk/gr_agreement.png) |
| cu-btc | [cu-btc/gr_agreement.png](cu-btc/gr_agreement.png) |

To regenerate one, point the script at the three result folders (the Fortran and
Rust folders of `run_validation.sh`, and the wrapper's folder from
`python/tests/test_experimental_data.py` or `test_cu_btc.py`), from the `python/`
folder so that the wrapper is importable:

```bash
cd python
uv run python ../validation/plot_gr_agreement.py water WORK/fortran WORK/rust WORK/wrapper ../validation/water/gr_agreement.png
```

### Discrepancies and issues found

1. **`fmole` with the molecule's `.mol` file present: not ported (Rust stops).**
   The Fortran `fmole` first runs `makemole` on every template whose `.mol` is in
   the folder. That is the normal situation, because the user just ran `makemole`.
   The Rust port stops with `not supported by the Rust port: fmole: template
   ./Water.mol exists ...`. This is a known, documented gap
   (`rust/epsr-core/src/tools/fmole.rs`, header comment). Porting it needs care:
   that `makemole` runs inside `fmole` and writes `fmole`'s random-number state into
   the template `.ato`. On the first attempt the Rust `fmole` exited with that error,
   so `Water.ato` and `waterbox.ato` differed from the Fortran from then on.
   *Workaround used:* delete `Water.mol` before `fmole` (step 4 of `water/steps.txt`).
   The Rust code was not changed.
2. **The setup programs depend on the current working directory (faithful, but
   matters for the wrapper).** `mixato` (and the other interactive programs) read
   `system_commands.txt`, write `dirlist.txt` and open typed file names in the
   *current* directory, not in the folder passed to them. Run from another
   directory, both the Fortran and the Rust `mixato` report "No .ato files found"
   and write no box. A wrapper must `chdir` into the folder. Setting up several
   folders in parallel therefore needs separate processes, not threads. This
   contradicts a remark in `rust/docs/PYTHON_WRAPPER_NOTES.md` §7; see
   PYTHON_WRAPPER_WORKFLOW.md §6.1.
3. **Unexplained one-off failure of the API driver (not reproduced).** One early
   run of both cases through `api_driver`, launched in parallel, stopped at the
   second `Epsr::setup` with `Io("No such file or directory (os error 2)")`. The
   message names no file. Both cases then ran through the same driver binary
   without error, one of them concurrently with other jobs. The error could have
   come from the driver or from the engine; I couldn't tell which. It is listed
   here because an `EpsrError::Io` without a file name is hard to diagnose. Adding
   the path to those messages would help wrapper users.

Faithful Fortran behaviours worth knowing (reproduced by the port, not bugs in it):
a `.pcof` with zero pairs switches Coulomb off; the default `.pcof` has `ereq = 0`;
`ireset`/`iinit` are reset to 0 after a run; `nrtype 5` treats Q as bin edges.
See PYTHON_WRAPPER_WORKFLOW.md §6.

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
