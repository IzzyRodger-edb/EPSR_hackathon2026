# EPSR-test-ex121 — Bragg peaks with Lorentzian broadening, hklqmin cut, triclinic crystal

**Program:** `epsr` (Bragg peak calculation; Rust `epsr-core/src/bragg.rs`).

**System:** triclinic crystal of tethered one-atom molecules built with `makelattice`
from `tric.cell` (a,b,c = 4.3, 4.8, 5.4 A; alpha, beta, gamma = 80, 100, 110 deg;
S at (0,0,0), O1 at (0.5,0.45,0.55), no charges) replicated 7 x 7 x 6 -> 588 atoms,
box 30.1 x 33.6 x 32.4 A. The atom labels S and O1 exist in the EPSR-test `.pcof` and
weights files, which are used unchanged (only the S-S, S-O1, O1-O1 pairs are matched).
`scripts/make_case.sh FORTRAN_BUILD_DIR OUT` regenerates all inputs (verified byte-identical).

**Keywords changed from EPSR-test:** `iinit 1`, `ntimes 2`, `qwidthqmax -0.02 4.0`
(negative width = Lorentzian peak shape), `mplicities 7 7 6 1 1 2` (6 words: usea/useb/usec
read, the value 2 is reset to 1), `hklqmin 1.0 8.0 0.0` (reflections below 1 A^-1 are not
broadened, rminbragg 8 A, no Debye-Waller factor).

**What is exercised:** triclinic reciprocal lattice (non-zero off-diagonal components, the
k/l ranges depending on h and k), Lorentzian branch (nfirst = 1), `hklqmin` cut (nlast = 0),
debyewaller = 0 branch, tethered molecules in `crystfrac`, non-cubic box.

**Verification:** ALL IDENTICAL (both passes).
