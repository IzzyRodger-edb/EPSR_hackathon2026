# Known Bugs in the EPSR Fortran and Their Status in the Rust Port

The Rust port was validated byte-for-byte against the original Fortran, which
means it initially reproduced the Fortran's bugs on purpose. This document
catalogues every documented bug (marked in the Fortran source or found during
the port), says which have now been **fixed in the Rust** (a deliberate,
documented divergence from the Fortran), and gives the possible fix and blast
radius for the ones still reproduced.

Ground rules for fixing more of these:
- A fix must not change function names, signatures, or module structure
  (concurrent work: Python wrapper, validation harness).
- The blast radius below says which outputs change. Fixing an "always active"
  bug changes **every** engine output: all 26 reference hashes,
  the 232 ex-test passes, and the published water/silica validation numbers
  (`validation/README.md`) would all need regeneration in a coordinated way.
- Live Fortran-vs-Rust comparisons (`run_ex_test.sh`, `chain_test.sh`,
  `validation/run_validation.sh`) will report differences for any fixed case;
  those differences are the fix working. The affected case READMEs say so.

## 1. Fixed in the Rust port (deliberate divergence from the Fortran)

### 1.1 Three-body angles corrupted to radians after every save (the big one)
- **Fortran:** `read_and_save_files.f` — `read_pcof` stores each requested
  triplet angle as a **cosine** (`cos(piconv*deg)`, lines 1369-70) and forces
  `threebawidth = 1.0` (1348-49). `save_files_others`, which runs at the end
  of **every iteration**, re-parses the same `tb` line and overwrites the
  in-memory arrays with the angle in **radians** (816-17), the width in
  radians (797-98), and the amplitude without `abs()`. The `.pcof` file on
  disk stays correct; only memory is corrupted.
- **Symptom:** from iteration 2 of any multi-iteration run with `tb` lines,
  `en_mole_triplets`/`pr_mole_triplets` compare cos θ against radians. The
  energy minimum moves to an unreachable value and the scale blows up
  (ex110: ~0.2 → ~2000 kJ/mol per molecule between iterations 1 and 2),
  corrupting MC acceptance, energy/molecule, and pressure. A fresh `epsr`
  invocation re-reads the `.pcof` and is correct again for its first iteration.
- **Fix (Rust):** `threebody.rs::write_pcof_threebody` keeps the re-parsed
  values local (same arithmetic, so the written line is byte-identical) and
  leaves memory in `read_pcof`'s conventions. This also restores the `abs()`
  on amplitudes that the Fortran's save path dropped.
- **Divergence (verified):** ex110, ex111 differ from the Fortran from
  iteration 2. Measured on ex110: the Fortran's iteration-2 energy/molecule is
  +1953.6 kJ/mol (corrupted) while the fixed Rust stays at −39.4, consistent
  with iteration 1; ex111: +688.9 vs −38.5. ex112–114 (niter 1) remain
  byte-identical.

### 1.2 Coordination-constraint angle built from the same arm twice
- **Fortran:** `coordconstraint_mod.f90:623` uses `adistcc(i,ia)` where the
  comment above it says "separation of ia and ja2", i.e. index `j`. Every
  `tb`-angle coordination constraint used one arm's length for both arms.
- **Fix (Rust):** `coordconstraint.rs` (`find_all_constraints` angle pass) now
  uses `adistcc[(j, ia)]`. Divergence (verified): ex117 and ex118; ex115,
  ex116 and ex119 do not exercise angle constraints numerically and remain
  byte-identical.

### 1.3 New pair added at save time gets no Coulomb term and no EP refinement
- **Fortran:** `save_pcof` leaves `coulfac`, `coordex`, `efacp` of a pair that
  first appears at save time as unset COMMON data (0.0).
- **Fix (Rust):** `pcof.rs` gives the new pair the same defaults a pair read
  from file gets (`coulfac = 1.0`, `efacp = 1.0`, `coordex = 0.0`).
  Divergence (verified): ex108 (the Fortran's next-pass refinement starts at
  +412.7 kJ/mol without the pair's Coulomb term; the fixed Rust at −51.4).
  ex109 (fresh start, no `.pcof`) is unchanged: its created pairs are re-read
  from file — acquiring the same defaults — before anything uses them.

### 1.4 Malformed Gaussian in `.pcof` double-counted on save
- **Fortran:** `read_and_save_files.f:748-750` does `igauss = igauss + 1` on a
  read error in `save_files_others`, where `read_pcof` (1166-67) decrements.
- **Fix (Rust):** `pcof.rs` decrements, matching the read path. Error-path
  only; no test case exercises it.

### 1.5 Pressure-control keyword `ca` selects no axis
- **Fortran:** the accepted-keyword list contains `ca` but the axis branch
  tests `ac`, so `ca` is accepted yet controls nothing (`npcontrols = 0`) and
  `ac` is rejected outright.
- **Fix (Rust):** `pressure.rs` accepts both `ca` and `ac`, selecting the c+a
  axes. Divergence: none today — no test case uses `ca`/`ac` (ex105 uses
  `abc` and `a` and stays byte-identical).

## 2. Still reproduced by design — always active (fixing changes every output)

These are genuine defects, but each one perturbs the Monte Carlo trajectory of
**every** run, so fixing any of them invalidates all byte-identity validation
at once. Fix them, if at all, together and with a coordinated regeneration of
`reference_hashes.txt`, the ex-test expectations, and the validation numbers.

| Bug | Location (Fortran / Rust) | Possible fix | Effect of the bug |
|---|---|---|---|
| First neighbour's pair energy never stored: the store loop in `fmolec1` starts at `inn=2` | `ato_manipulation_rev4.f:2451-59` / `energy.rs` (store loop after the parallel pair pass) | start the loop at `inn=1` | One candidate pair per molecule missing from the stored energy list every iteration |
| Whole-molecule translation step permanently doubled when a move rounds to no change | `molfit.f:727-28` / `montecarlo.rs:54-57,328-330` | double a local copy, not the per-type step | Step-size adaptation drifts upward spuriously |
| `entrpb` noise loops over the first `nfile` pairs instead of all pairs | `epsr_standalone_rev3.f:3357-59` / `structure.rs` (`entrpb`) | loop to `npairs` | Noise figure in `.out` computed from the wrong columns |
| Single-precision `cosd(90°)` = −4.4e-8 skews a nominally cubic box | gfortran intrinsic / `ato.rs:44-51` | treat 90° exactly (return 0.0) | `blatx` not exactly 0; minimum-image uses the general-box path |
| `rfacraw` uninitialised on non-final iterations (`.erg` column 3) | `epsr_standalone_rev3.f:2317-26` / `epsr.rs` (`rfacraw_last`) | compute it every iteration | Garbage/stale value in `.erg` until the last iteration |
| q/Q atoms' exemption from the `roverlap` minimum never applies (tests the number line, not the labels) | `read_and_save_files.f:1183-89` / `pcof.rs:258-264` | test `tlab1`/`tlab2` | q/Q pairs get a minimum radius they shouldn't |
| `generate_mp` reads uninitialised `it`/`jt`, so container-force lists are always empty there | / `potential.rs:561-563` | derive `it`/`jt` from the pair index | Container forces never contribute to the molecular potential |

## 3. Still reproduced — feature-gated or inert

| Bug | Location | Notes / possible fix |
|---|---|---|
| `ireset 1` misaligns the `tb` symbol slot, so the save-time re-parse is skipped and in-memory cosines are written as angles (−57.30, width 57.30) | `threebody.rs:180` (`nwords > 12` guard); ex113 | Fix = address symbols by name, not slot index. Left as is: the written file is wrong either way in the Fortran, and ex113 documents it |
| `nthreebangsmax` set from the wrong index at save time | `read_and_save_files.f:811`; Rust reports `unsupported` (`threebody.rs`) | Unreachable in practice; Rust already fails loudly instead of corrupting |
| `get_triplet_pot_cons` subtracts cosines from radians (iteration 1) | `threebody_mod.f90:283` / `threebody.rs:423-424` | Print-only value (`consthreeb`); fix = convert consistently |
| Coordination-constraint stale-slot family: locals not reset, search stops at slot 2, stale slot 1 read, `save_coord_energies` only after rotations | `coordconstraint.rs:466, 585, 611, 968, 1039-41` (ex115–119) | Each is a small index/reset fix; fix together with a dedicated regeneration of the cc ex-cases |
| `awidth` field of a `tb` line ignored on read (forced to 1.0) | `read_and_save_files.f:1348-49` / `threebody.rs:110` | Arguably intended (width folded into amplitudes); honouring the field would change the tb energy scale |
| `get_direct_bonds` assumes the angle constraint is the triangle's longest side | `ato_manipulation_rev4.f:2530` (marked `!Hence this method would not work...`) | Wrong for bond angles under ~60°; needs a real topological test |
| `changepotfac` branch compares against the wrong absolute energy | `epsr_standalone_rev3.f:6043-44` (marked `!CAUTION!`) / `potential.rs:814` | Unreachable (never called with `.true.`) |
| `phix/phiy/phiz` transformed for untethered molecules; `rho`/`xstar` not recomputed after a box change | `pressure.rs:145` and ex105 README | Needs physics review before fixing |
| `moloutc` crash for some inputs (SIGSEGV in Fortran) | `ato.rs:1011-19` reproduces as empty file + error | Rust already degrades gracefully |

The analysis tools (`tools/` modules) reproduce a further set of
uninitialised-data quirks, e.g. `epsrwts.rs:97`, `xwts.rs:142`,
`nwts.rs:84,261`, `dockato.rs:332-340`, `sdfcube.rs:12-15`, `readpdb.rs:19`,
`bonds.rs:101`, `fmole.rs:21` and others — each marked with an `NB:` comment
at the site. They affect standalone tools, not the engine; fix on demand with
a per-tool reference regeneration.

## 4. Historical notes (nothing to do)

- `epsr_standalone_rev3.f:122-153` describes bugs fixed in 2006–2007
  (atom-type forwarding, pressure sign, negative X-ray form factors).
- `molfit.f:85` and `:119-120` are stale warnings about code that was since
  corrected (`centreofmassi` after accept; neighbour interactions included).
- `epsr_standalone_rev3.f:367`: `efilereq` must not be zero for every dataset
  when `feedback > 0` — a documented input restriction, not checked anywhere.

## 5. How a fix here gets validated

1. The main test case, `tests/reference_hashes.txt`, `thread_invariance.rs`,
   and the water/silica validation use no `tb` lines, coordination
   constraints, pressure control, or save-time pairs — they must stay green
   unchanged (checked on glibc Linux; the hashes do not hold on Windows libm).
2. Run `rust/tools/run_ex_test.sh` for the affected ex-cases against the
   unpatched Fortran reference: the diff set must match exactly the cases the
   fix predicts, and the affected READMEs must say why.
3. A full `rust/tools/run_all_ex_tests.sh` sweep must show divergence **only**
   in the documented cases.
