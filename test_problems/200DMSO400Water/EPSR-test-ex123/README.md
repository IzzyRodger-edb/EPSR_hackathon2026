# EPSR-test-ex123 — lattice refinement of a triclinic crystal (with Bragg peaks)

**Program:** `epsr`: the `do while (ilatrefine.lt.nlatrefine+1)` loop of `runepsr17`,
`change_boxsize`, `convertphibthetacphic`, `convertalphabetagamma` (Rust `bragg.rs`).

**System:** the triclinic S/O1 crystal of EPSR-test-ex121 (588 tethered atoms, see
`scripts/make_triclinic_crystal.sh`); EPSR-test `.pcof`/weights/data used unchanged.
`scripts/make_case.sh` regenerates the inputs (verified byte-identical).

**Keywords changed from EPSR-test:** `iinit 1`, `ntimes 2`, `qwidthqmax 0.02 3.0`,
`mplicities 7 7 6`, `hklqmin 0.0 8.0 0.05`, `nlatrefine 3`,
`latrefine 0.2 0.2 0.2 1.0 1.0 1.0`.

**What is exercised:** per refinement step: `atogen` + `gensq` + Bragg correction, the raw
R-factor against the three datasets, acceptance/rejection of the trial lattice, random trial
a, b, c, alpha, beta, gamma (non-cubic branch, 6 `ran1` calls), `convertalphabetagamma`,
`change_boxsize` (gfortran intrinsic `cosd`/`sind`, new `dmax`/`cell`, fractional
coordinates of molecule centres and tether points), the constant-density rescaling
(`vol**(1/3)`), and `(un)normalise_coordinates` with the changed `dmax`. In pass 1 the
trials of steps 1 and 2 are accepted and the last one is rejected (box restored without
recomputing the structure); in pass 2 step 1 is rejected, 2 accepted, 3 rejected. The refined
lattice is written to the `.ato` file and read back in pass 2.

**Fortran quirks reproduced:** `clatnew` is rescaled with `blatfac`; when the last trial is
rejected the box is restored but S(Q) is not recomputed (the accumulated structure belongs
to the rejected box).

**Verification:** ALL IDENTICAL (both passes).
