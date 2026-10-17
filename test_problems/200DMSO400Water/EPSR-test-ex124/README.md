# EPSR-test-ex124 — lattice refinement of a cubic box, Bragg peaks with nbraggup = nq

**Program:** `epsr` (lattice refinement, cubic branch; Bragg peaks; Rust `bragg.rs`).

**System:** DMSO/water crystal, the relaxed 6.84 A unit cell of EPSR-test-ex120 replicated
4 x 4 x 4 with `makelatticeato` -> cubic 27.4 A box (384 molecules, 2048 atoms).
`scripts/make_case.sh` regenerates the inputs (verified byte-identical).

**Keywords changed from EPSR-test:** `iinit 1`, `ntimes 2`, `nq 100` (Q up to 5 A^-1),
`qwidthqmax 0.01 4.5`, `mplicities 4 4 4`, `hklqmin 0.5 5.0 0.05`, `nlatrefine 2`,
`latrefine 0.3 0 0 0 0 0`.

**What is exercised:** cubic `crystfrac` branch; `nbraggup` limited to `nq`
((4.5 + 3/5)/0.05 + 1 = 103 > 100); `hklqmin` cut with Gaussian peaks; cubic branch of the
lattice refinement (one `ran1` call, b = c = a, angles 90). For a cubic box the trial change
of `a` is undone exactly by the constant-density rescaling, so the fit does not change and
every trial is accepted (`newfit <= oldfit`); on the last step `ilatrefine` is advanced to
skip the recomputation — this is the Fortran behaviour.

**Verification:** ALL IDENTICAL (both passes).
