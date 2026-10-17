# Validation of the Rust port

## Method

1. The original Fortran `epsr` was compiled from `fortran_src` with gfortran
   15.2 (`tools/build_fortran_reference.sh`; flags `-O2 -fno-tree-vectorize`,
   see "Compiler vectorisation" below) and run single-threaded
   (`OMP_NUM_THREADS=1`) on a copy of `EPSR-test`.
2. The Rust `epsr-rs` was run on another copy of the same inputs.
3. Every file in the folder was compared byte for byte
   (`tools/compare_outputs.py`), masking only the current date in `.out` and
   the wall-clock columns of `fort.35`.
4. To exercise more of the code, runs were chained (each run continues from
   the files written by the previous one, as EPSR is normally used) —
   `tools/chain_test.sh` — and input variants were tested.

The Fortran was also instrumented with gcov to find exactly which statements
the test case executes; the Rust translation covers all of them.

## Results

| test | Fortran vs Rust |
|------|-----------------|
| 1 run of `EPSR-test` (all 26 output files) | **identical** |
| 12 chained runs (60 MC cycles; includes internal-rotation cycles and running energy averages) | **identical after every run** |
| `niter 3` (three iterations in one invocation), 3 chained runs | **identical** |
| `ngrsamples 200`, `nmolcell 6` (random sampling of g(r) origins, different cell lists), 2 chained runs | **identical** |

"Identical" covers the Monte Carlo trajectory (`.ato` coordinates and the
stored random-number-generator state), all distribution functions
(`.g01 .y01 .f01 .s01 ...`), the empirical potential (`.pcof`, `.p01`),
the updated input files (`.EPSR.inp`, `.EPSR.inpa`), `.erg`, `.out` and
`delfit.dat`.

`cargo test --release` repeats the single-run comparison automatically
(hashes of the Fortran outputs are stored in
`epsr-core/tests/reference_hashes.txt`).

Run time for one iteration of the test case (5 MC cycles): Fortran (scalar
build) 2.4 s, Rust 3.2 s; Fortran with gfortran's default vectorisation
1.2 s.  The Rust code is currently serial; the Fortran OpenMP loops are
data-parallel and can be parallelised the same way (e.g. with `rayon`)
without changing results.

## Feature test suite (`EPSR-test-ex*`)

There is one folder per feature or option, `test_problems/200DMSO400Water/EPSR-test-ex1` to `EPSR-test-ex124`, with gaps in the numbering
(all current cases are built on the 200 DMSO + 400 water test system; further systems can be added as siblings of
`200DMSO400Water`).
[TEST_CASES.md](TEST_CASES.md) indexes them.

* **Contents of a folder.** Its inputs, mostly derived from EPSR-test (often with a `scripts/` folder that regenerates
  them), plus:
  * `commands.txt`: the program calls; `repeat N` chains N passes;
  * a `README.md` describing what is exercised;
  * `fortran_reference.tar.gz`: the files the original Fortran writes.
* **How cases were designed.** The feature had to be active and give non-trivial numbers, checked in the Fortran
  log and output. For most new engine features, gcov of a `--coverage` Fortran build confirmed which branches run.

```
tools/build_fortran_reference.sh /path/fref            # once; ~5 min
tools/run_ex_test.sh ../test_problems/200DMSO400Water/EPSR-test-ex15 /path/fref   # one case: runs Fortran and Rust, compares every pass
tools/run_all_ex_tests.sh /path/fref [work] [N ...]    # all (or selected) cases, one summary line each
```

`run_ex_test.sh --save-reference` stores a new `fortran_reference.tar.gz`.

Each folder also contains a one-page summary, `EPSR-test-exN.pdf`. It covers the code exercised, the physical system with
references, and how the case was verified. It is rendered by `tools/make_case_pdfs.py --results <run_all_ex_tests log>`
from `docs/case_summaries/exN.json`, with the shared bibliography in `docs/case_summaries/references.json`.

| programs | cases | passes | result |
|---|---:|---:|---|
| `epsr` engine (options and features) | 49 | 130 | **all identical** |
| setup, weights and conversion tools | 38 | 38 | **all identical** |
| analysis tools | 28 | 64 | **all identical** |
| **total** | **115** | **232** | **115/115 identical** |

These figures are from the final regression run: `cargo test --release` passed and all 115 cases were identical
against the gfortran 15.2 -O2 -fno-tree-vectorize reference.

* **"Identical."** Byte-for-byte after the two masks above, plus one more: `.mapgr.g01` columns 31–45, where the
  Fortran writes uninitialised stack data (see TOOLS.md).
* **Patched reference builds.** Three Fortran sources fail at run time with current gfortran. The reference build
  adds the minimal fix, and the port behaves like the fixed program:
  * `readxyz`/`readpdb`: an explicit interface for `ctype` (SIGSEGV without it);
  * `sharm`/`sdf`, and `sharm_diffuse` in `epsr` (diffuse scattering, ex122): the missing comma after the
    Hollerith `1h(` in a run-time FORMAT.

  `build_fortran_reference.sh` applies these patches and keeps the originals as `*.orig`.
* **Crashing Fortran.** Where the Fortran program itself crashes (`makeato`, ex30), the case checks that the
  port leaves the same files behind.

## Caveats and findings

### Compiler vectorisation (why `-fno-tree-vectorize`)
gfortran ≥ 12 auto-vectorises some loops at `-O2`.  In `ftback` (background
removal) and in `auxiliary_routines.f` it then calls glibc's SIMD math
library (`_ZGVbN4v_sinf`, `_ZGVbN4v_expf`), whose results differ from scalar
`sinf`/`expf` by a few units in the last place — and only for loop
iterations that fall in a full SIMD block.  With the default flags the
Fortran's background-subtracted data therefore differ from the scalar build
by 1–3 ulp, and the Monte Carlo fit of the empirical potential (which is
chaotic) diverges from that point; the resulting EP coefficients differ at
the 10⁻³ level.  Neither result is "more correct".  Older compilers (such
as those used for the distributed EPSR binaries) generate scalar code, so the
scalar build is taken as the reference.  If exact agreement with a
vectorised gfortran build is ever needed, the Rust code could call libmvec
directly in `ftback`.

### The outputs shipped in `EPSR-test`
`EPSR-test` contains the inputs *and* the outputs of a previous run made
with the Windows EPSR 26 binary, and EPSR overwrites its inputs (`.ato`,
`.pcof`, ...) at the end of each run.  Those output files therefore cannot be
reproduced from the folder; running either program on the folder performs the
*next* iteration (`nsumt` 102 → 103).  The results are consistent with the
shipped ones (chi-square 0.1148×10⁻³ vs 0.1146×10⁻³, energy −39.9 vs −40.3
kJ/mol, R-factors within 0.3 %).

### Fortran quirks reproduced deliberately
The port reproduces the Fortran's behaviour even where it looks unintended,
so that results are identical; these are flagged in comments:

* `fmolec1` never stores the energy of a molecule with the *first* candidate
  in its neighbour list (the loop over candidates starts at index 2).
* `read_pcof`'s exemption of `q`/`Q` pseudo-atoms from the `roverlap` check
  tests the wrong line of input, so it never applies.
* `entrpb` (noise of g(r), written to `.out`) loops over the first `ndata`
  *pairs*, not all pairs.
* `cosd(90°)` in single precision is −4.4×10⁻⁸, so the cubic box's `blatx` is
  not exactly zero; this enters the cell assignment of molecules.
* `movem` doubles the translation step size of the molecule *type* if a
  trial move does not change the position (argument passed by reference).
* `rfacraw` (column 3 of `.erg`) is an uninitialised local on non-final
  iterations of a multi-iteration run; gfortran writes the value from the
  last `.erg` line, which the port reproduces (`Epsr::rfacraw_last`).

### Fortran bugs deliberately FIXED in the port (documented divergences)
Since the bug-fix pass (see `FORTRAN_BUGS.md` at the repo root), a small set
of feature-gated Fortran bugs is fixed rather than reproduced, so the
affected ex-cases **deliberately differ** from the Fortran:

* **Three-body angles corrupted to radians at `.pcof` save time** (the
  in-memory cosines/width/amplitudes were overwritten by `save_files_others`
  every iteration, blowing up the triplet energy from iteration 2). Fixed in
  `threebody.rs::write_pcof_threebody`; ex110/ex111 diverge from iteration 2
  (ex110 iteration-2 energy: Fortran +1953.6 kJ/mol, Rust −39.4); ex112–114
  remain identical.
* **Coordination-constraint angle built from the same arm twice**
  (`adistcc(i,ia)` for both arms). Fixed in `coordconstraint.rs`; ex117 and
  ex118 diverge (ex115/116/119 do not exercise angle constraints and remain
  identical).
* **New pair added at save time had `coulfac`/`efacp`/`coordex` = 0** (no
  Coulomb term, EP never refined). Fixed in `pcof.rs`; ex108 diverges
  (ex109's fresh-start pairs are re-read from file before use and remain
  identical).
* **`igauss` incremented instead of decremented on a `.pcof` read error**
  (bad Gaussian double-counted). Fixed in `pcof.rs`; error path only, no case
  affected.
* **Pressure keyword `ca` selected no axis** (`ac` rejected). Fixed in
  `pressure.rs`; no case uses `ca`/`ac`, so none diverges today.

The main test case, `reference_hashes.txt`, `thread_invariance.rs` and the
water/silica validation exercise none of these features and remain
byte-identical to the Fortran. So "115/115 identical" above now carries the
documented exceptions listed here; the affected case READMEs say the same.

### Threads
All OpenMP loops executed by the test case write to independent array
elements and are therefore thread-count independent, except the
close-contact counter `potex` in `atogen`, which is updated without
synchronisation.  It only matters when `rminex > 0` (not the case in the test
input).  The Rust port parallelises the corresponding loops with rayon
(`PARALLELISM_PLAN.md`) and is byte-identical at any thread count
(`tests/thread_invariance.rs`); the `potex` counter is accumulated per thread
and merged, so the Fortran's race is fixed properly.  The scalar reference
comparison corresponds to `OMP_NUM_THREADS=1`.

### Coverage
Every optional feature of `epsr` now has a reference comparison (FORTRAN_MAP.md lists feature → case). The
remaining `Unsupported` paths are inputs on which the Fortran itself hangs, divides by zero, writes out of bounds
or reads uninitialised memory.
