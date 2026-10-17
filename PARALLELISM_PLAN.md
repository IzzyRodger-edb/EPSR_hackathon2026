# Parallelism Plan for the EPSR Rust Port

This document identifies where parallelism will actually improve the performance of the
Rust EPSR engine (`rust/epsr-core`, `rust/epsr-cli`), ranks the candidates by likely
speedup versus risk, and records the constraints every change must respect. It is an
analysis only — no code has been changed.

## Summary

- The workload is **entirely CPU-bound**: tight numeric loops over atom pairs, Q-points
  and r-bins. There is no meaningful blocking I/O, so **async/tokio has no role here**;
  the right tool is data parallelism (**rayon**, approved as the project's first
  external dependency — `rust/docs/VALIDATION.md` already anticipates exactly this).
- The governing constraint is **byte-identical output**: the test suite
  (`rust/epsr-core/tests/test_case.rs`) FNV-hashes every output file against
  `tests/reference_hashes.txt` with zero tolerance, and 232 ex-test passes are
  byte-for-byte against the Fortran reference. Every candidate below has a recipe that
  preserves bit-exactness; results must be identical at any thread count.
- A second hard constraint: **the code structure is frozen**. Other agents are working
  on this code concurrently. Function names, signatures and observable behaviour must
  not change. All parallelism goes *inside* existing function bodies (internal parallel
  loops, per-thread local buffers) — no renames, no signature changes, no splitting or
  moving functions.
- The highest-value targets are the `atogen` O(N²) pair histogram, the `fit_expansion`
  chi² evaluation, and the `gensq` Fourier transform. The Monte Carlo move loop itself
  is inherently sequential (single RNG stream) and must be left alone.

## Constraints that shape every change

### 1. One RNG stream, and its draw order is output

`Ran1` (`rust/epsr-core/src/rng.rs:27-65`) is a single Park–Miller generator living in
`State.rng`; its full state is saved into the `.ato` file, so the global draw order is
part of the validated output. Consumers, in draw order per iteration: MC moves
(`montecarlo.rs`), `atogen` origin sampling (`structure.rs:521`), optional lattice
refinement (`bragg.rs:966-984`), `fit_expansion` (`potential.rs:404-416`),
`automatic_ereq` (`epsr.rs:928-942`). Any phase that changes how many values are drawn,
or in what order, changes everything downstream. **No RNG call may move onto a worker
thread.**

### 2. f32 reduction order must be preserved

All accumulations are `f32` (some mixed with `f64`) in a fixed sequential order; a 1–3
ulp change cascades into visible differences through the chaotic EP fit
(`VALIDATION.md`). Three patterns are safe and cover every candidate:

- **(a) Disjoint writes** — iterations write different array elements; no reduction at all.
- **(b) Exact integer counts** — histogram increments of +1.0/+2.0 are exact in `f32`
  below 2²⁴, so per-thread histograms merged in any order are bit-identical.
- **(c) Parallel partials, serial reduce** — compute independent per-element values in
  parallel into an indexed buffer, then sum serially in the original index order.

### 3. The `GEOM` thread-local (the one real trap)

`thread_local! static GEOM: Cell<BoxGeom>` (`neighbour.rs:43-47`) is read by
`minimagec` (`neighbour.rs:52-65`) on every minimum-image calculation, and
`sync_lattice` (`neighbour.rs:69-73`) sets it **only on the calling thread**. A rayon
worker sees the default geometry — which happens to be correct for cubic boxes
(box-fraction units) but is **silently wrong for non-cubic boxes**.

Prerequisite fix, before any parallel region that calls `minimagec`: initialise the
thread-local on every worker — rayon's pool `start_handler` / `broadcast` calling
`sync_lattice`, or re-sync at the top of each parallel closure (it is a cheap `Cell`
write). Do **not** change `minimagec`'s signature (frozen-structure rule). Add a
non-cubic ex-case run to the verification of the first parallel region.

### 4. Shared scratch inside `&mut self` kernels

Several kernels use fields of `State` as scratch: `en_change` writes `erg.ergnew`;
`press` uses `erg.ergold/ergnew/nbref` per molecule; `atogen` uses
`iz1values`/`icrefvalues` and bumps `pot.potex`, `csampled`, `cg.wavinter`. Parallel
versions need per-thread local buffers for these. (Note: `VALIDATION.md:136-141` records
that the Fortran OpenMP build actually **races** on `potex`; the Rust version can fix
that correctly with per-thread counts, since the counts are exact — pattern (b).)

Good news otherwise: there is no `Rc`, `RefCell`, `static mut` or `Mutex` anywhere in
`src`; `State` is owned `Vec`s and scalars and is `Send` (tested in
`test_case.rs:84-88`).

## Step 0: measure before implementing

The iteration loop already logs a `timing>` line after each phase
(`epsr.rs:621,646,700,796,888`). Capture that profile on the test case
(`EPSR-test/200DMSO400Water`) first and confirm the ranking below — the estimates come
from loop bounds (600 molecules / ~3200 atoms, nq=600, ~550 r-bins, 3 datasets, 116 EP
terms, ntimes=5), not measurements. Reported serial Rust time: ~3.2 s per iteration.

**Measured** (2026-10-17, Windows 11, 8 cores, serial release build, 1 iteration,
5.35 s total):

| Phase | s | Parallelisable content |
|---|---|---|
| Monte Carlo cycles | 1.59 | none (sequential RNG chain) |
| output files | 1.32 | `w_intr` formatting; last iteration only, amortises over `niter` |
| potential refinement | 1.16 | candidate #2 (`fit_expansion` chi²) |
| structure accumulation + pressure | 0.69 | candidates #1 (`atogen`), #3 (`gensq`), #5 (`press`) |
| neighbour lists + energies | 0.44 | candidate #4 (`fmolec1`) |

So ~2.3 s of the 5.35 s is addressable by candidates #1–#5 on this case; the MC phase
is the serial floor. `gensq` (#3, done) is a small slice of its 0.69 s phase — it was
implemented first as the low-risk proving ground, not for its absolute gain.

### Validation caveat on Windows

The byte-for-byte reference hashes **do not pass on Windows even for unmodified
serial HEAD** (15 of 26 files differ; the `trig_intrinsics_match_gfortran` unit test
fails by 1 ulp) — exactly the glibc-libm caveat in VALIDATION.md. On Windows, validate
parallel changes by **self-consistency** instead: run unmodified serial HEAD, the
parallel build at `RAYON_NUM_THREADS=1`, and at full thread count on copies of the test
case; all output files must be byte-identical across the three runs (mask the `.out`
date line and `fort.35` timing columns). Final reference-hash validation still happens
on a glibc Linux machine.

## Prerequisites

1. Add `rayon` to `epsr-core` (**approved** — first external dependency).
2. Fix the `GEOM` thread-local as described above.
3. Thread-count control: use rayon's default pool (respects `RAYON_NUM_THREADS`), so
   single-threaded runs remain available for debugging and validation.

## Ranked candidates

All are CPU-bound and suit rayon. Suggested implementation order: **#3 first** (lowest
risk, proves the harness and the GEOM fix is not even needed there), then #1, #2, #4, #5.

### 1. `atogen` intermolecular histogram — high gain, medium risk

`structure.rs:500-596`. Brute-force upper-triangle loop over all atom pairs
(~5M distance evaluations with `sqrt`+`nint` on the test case, **no cell list**), once
per iteration (more if lattice refinement is on). Likely a dominant phase.

- **Parallel axis:** origin molecule `im`.
- **Bit-identical recipe:** the one RNG draw per `im` (`structure.rs:521`) does not
  depend on positions, so draw all `nmol` values serially up front; the
  `molsampled`-dependent accept decisions also stay serial using those pre-drawn values.
  Then run the per-`im` pair work in parallel with per-thread `wav`/`potex`/`csampled`/
  `wavinter` histograms (all exact counts — pattern (b), merge in any order) and
  per-thread `iz1values`/`icrefvalues` scratch. Requires the GEOM fix.
- **Risk:** medium — several scratch arrays to privatise and the RNG pre-draw to get
  exactly right, but the merge is order-insensitive by construction.

### 2. `fit_expansion` chi² evaluation — high gain, medium risk

`potential.rs:409-413` and `428-431`, inside the EP-coefficient Monte Carlo fit.
Possibly the single largest phase: ≥1000 steps × 3 datasets per iteration, each step
evaluating chi² over ~600 Q-rows × 116 terms (0.2–3 G multiply-adds per iteration).

- **Parallel axis:** Q-row `i` inside one chi² evaluation. The step chain itself is a
  sequential Markov chain (RNG + `wtold` updates) and **stays serial**.
- **Bit-identical recipe:** `expand(i)` rows are independent (the loop-carried clamp is
  within each `i`); compute residual² per `i` in parallel into a buffer, then sum
  serially in index order — pattern (c).
- **Risk/caveat:** granularity. ~70k FLOPs per step is close to rayon's fork-join
  overhead, so use the persistent pool with chunked `par_iter` and **measure**; expect
  ~2–3× on this phase, not linear scaling. If overhead wins, fall back to leaving it
  serial — the measurement in Step 0 decides.

### 3. `gensq` Fourier transform (g(r) → S(Q)) — medium gain, lowest risk — **DONE**

`structure.rs:638-665`. npairs × 600 Q × ~550 r × 2 serial sums, once per iteration.
`&self`, writes are disjoint per `(i, ic)` into `sofq`, inner sums stay serial —
pattern (a)+(c) for free.

- **Parallel axis:** pair `ic` (or flattened `(ic, i)`), plain `par_iter`.
- **Implemented 2026-10-17**: `into_par_iter` over `ic`, per-pair column collected in
  order, serial write-back; rayon 1.12 added to `epsr-core`. Verified byte-identical
  across serial HEAD / 1 thread / 8 threads (self-consistency protocol above).

### 4. `fmolec1` all-pairs stored energies — medium gain, medium risk

`energy.rs:439-554`. Full molecule-pair energy pass over the cell neighbourhood, 1–2×
per iteration (plus after shake cycles).

- **Parallel axis:** cell / origin molecule `im`; `en_mole_pair` (`energy.rs:98`) is
  `&self` and independent per pair.
- **Bit-identical recipe:** compute pair energies in parallel into per-`im` local
  buffers; then a serial pass calls `set_ergpair` (`neighbour.rs:322-345`) in the
  original order — insertion order defines the neighbour order the MC loop later sees,
  so it must be preserved — and does the serial `ftest`/`ftestb` sums. Preserve the
  reproduced Fortran quirk that the store loop starts at `inn=2` (`energy.rs:485-495`).
  Requires the GEOM fix.
- **Risk:** medium — the serial ordering pass must be exactly faithful.

### 5. `press` energy/virial recomputation — medium gain, medium risk

`energy.rs:596-699`. Recomputes all stored-neighbour pair energies and the virial, once
per iteration; similar cost to `fmolec1`.

- **Parallel axis:** molecule `im`.
- **Bit-identical recipe:** per-thread scratch replacing the shared
  `self.erg.ergold/ergnew/nbref`; per-`im` partials for `sumc` (`f32`) and `presst`
  (`f64`) reduced serially in `im` order — pattern (c). Requires the GEOM fix.

### 6. `stog` (S(Q) → g(r)) — small gain, low risk

`structure.rs:125-152`. Called 4 times in `write_results`, last iteration only.
Parallel over r-bin `j` with a per-thread `sum` buffer; `gofr` writes are disjoint.
Worth doing only because it is nearly free once rayon is in.

### 7. Conditional targets (not exercised by the main test case)

- **`bragginten` structure-factor amplitudes** (`bragg.rs:440-452`), only when
  `nbragg > 0`: swap the loops so hkl index `ib` is outer (disjoint writes per `ib`,
  atom order within each `ib` preserved → bit-identical). The subsequent broadening
  scatter (`bragg.rs:491-535`) is order-dependent — leave serial.
- **Three-body** (`find_triplets`, `threebody.rs:250-405`), only when `nthreeb > 0`:
  brute force O(n1·n2·n3), but it appends to `tripletlist` in order — parallelise only
  with an order-preserving gather, and only if a three-body case shows it matters.
- **`generate_ep`/`generate_rp`** (`potential.rs`): independent per pair but cheap;
  touch only if Step 0 measurements say otherwise.

## Not worth it / unsafe — leave serial

| Region | Why |
|---|---|
| `molfit` MC molecule loop (`montecarlo.rs:232-437`) | RNG Markov chain; accepted moves mutate positions and `erg.ergpair` read by the next molecule; shared step-size adaptation and accept counters. Parallelising (e.g. checkerboard) changes the trajectory and breaks Fortran validation. |
| `fit_expansion` step chain (`potential.rs:396-446`) | Sequential Markov chain (only the chi² *inside* a step is parallel — candidate #2). |
| `invert_cwt` (`weights.rs:404-522`) | RNG-driven inversion; setup-only. |
| Lattice refinement outer loop (`bragg.rs:861-1012`) | Sequential MC over lattice parameters. |
| `en_change` per-move neighbour loop (`energy.rs:381-414`) | Independent in principle, but ~1.5k table lookups (~µs) per call — below rayon's fork-join overhead. Pointless. |
| Any I/O (async/tokio) | Inputs are whole-file reads at setup (`fortran/listread.rs:229-237`); outputs are buffered and written once (`fortran/output.rs:31-52`). ~6 MB of text per run against ~3.2 s of compute — I/O is a rounding error. The per-row `fmt` formatting cost in `w_intr` (`output.rs:52-86`) is CPU, runs once at the end, and is not worth threading. |

## Aside (not parallelism): per-move timer overhead

`en_change`/`en_change_intra` call `Instant::now()`/`elapsed()` on **every MC move**
(`energy.rs:360,367,374,432`) to accumulate a timer that almost always adds 0. This is
measurable serial overhead on the hottest path and could be cheapened independently of
this plan — flagged here, out of scope.

## Verification

For every parallel change, in order:

1. `cargo test` — the reference hashes in `tests/reference_hashes.txt` must pass
   **unchanged** (no regenerated references).
2. Run the test case at `RAYON_NUM_THREADS=1` and at the machine's core count; output
   hashes must be identical at any thread count.
3. Compare before/after `timing>` phase lines to confirm the change actually pays — if
   a region shows no gain (likely risk for #2), revert it to serial.
4. After the first region that calls `minimagec` in parallel (#1, #4 or #5): run a
   non-cubic-box case from `test_problems/.../EPSR-test-ex*` to exercise the GEOM
   worker-sync fix.
5. The full ex-suite (`tools/run_all_ex_tests.sh`, 232 passes) before merging.
