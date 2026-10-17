# EPSR-test-ex120 — Bragg peaks, Gaussian broadening, Debye-Waller factor (orthorhombic box)

**Program:** `epsr` (`epsr_standalone_rev3.f` + `bragginten` in `ato_manipulation_rev4.f`;
Rust: `epsr-core/src/bragg.rs`).

**System:** molecular crystal of DMSO and water: a 6.84 A cubic unit cell with 2 DMSO + 4 water
(`mixato` -> `randomise` -> 200 EPSR Monte Carlo cycles with `potfac 0` to remove overlaps),
replicated 5 x 4 x 4 with `makelatticeato` -> orthorhombic 34.2 x 27.4 x 27.4 A box, 480
molecules / 2560 atoms (written as a general, non-cubic box). Atom types are those of
EPSR-test, so the EPSR-test `.pcof`, weights and data files are used unchanged (they are
only numerical input here). `scripts/make_case.sh FORTRAN_BUILD_DIR OUT` regenerates all
inputs (verified byte-identical).

**Keywords changed from EPSR-test:** `iinit 1`, `ntimes 2`,
`qwidthqmax 0.01 3.0` (Gaussian width 0.01 + fwhmq*Q, Bragg calculation up to 3 A^-1),
`mplicities 5 4 4` (unit cell = box/(5,4,4); usea=useb=usec=1 by default),
`hklqmin 0.0 10.0 0.1` (all reflections, rminbragg 10 A, Debye-Waller 0.1).

**What is exercised:** keyword parsing (3-word `mplicities`, `nbragg`/`nbraggup`);
`bragginten`: reciprocal lattice, hkl enumeration (153 reflections, 8 unique Q values),
non-cubic `crystfrac`, amplitudes/intensities per pair, Gaussian broadening (both the
`Q,Qb > 0` branch and the Q = 0 branch), single-atom subtraction, `.bragg` file; in
`runepsr17`: `.inten` file, division by the number of unit cells, subtraction of the
intramolecular S(Q), Debye-Waller factor, truncation convolution (`rminbragg`), and the
replacement of the intermolecular partial S(Q) up to `qbraggmax` (which then feeds the
accumulated S(Q), the EP refinement and every output file). `repeat 2` re-reads all
written files.

**Verification:** `run_ex_test.sh` against the Fortran reference: ALL IDENTICAL (both passes).
