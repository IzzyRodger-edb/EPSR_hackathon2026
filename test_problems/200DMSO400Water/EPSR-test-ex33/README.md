# EPSR-test-ex33 — `mixato` + `fmole` (standard route)

**Programs:** `mixato`, `fmole` (`fortran_src/fmole.f`: `fmole`, `initpt_fmole`, `update_ato`
plus `molfit`, `fmolec1`, `adjust_intra_wts`, `generate_mp`, cell lists).

**Inputs:** `DMSO.ato`, `WaterMolecule.ato` (from `EPSR-test`), `mixato.stdin`.
No `.mol` files are present, so the `makemole` call that `fmole` makes for each
template does nothing (as in the Fortran).

**Commands:**
1. `mixato < mixato.stdin` — 10 DMSO + 20 water, density 0.1, output `mixbox`.
2. `fmole mixbox.ato 50` — 50 intramolecular shakes of a 30-molecule box (rotational
   groups switched off because nmol > 1); templates `DMSO.ato` and `WaterMolecule.ato`
   are read by `update_ato`.
3. `fmole DMSO 200` — single molecule: atom shakes **and** internal rotations of the two
   CH3 rotational groups (Metropolis test on intra + inter energy); the whole-molecule and rotation step sizes are zeroed during the run and restored, so the `.ato` step-size line is unchanged.

**Expected outputs:** `mixbox.ato`, updated `DMSO.ato`, `dirlist.txt`.

**Verification:** ALL IDENTICAL (coordinates, step sizes and random-number state).
