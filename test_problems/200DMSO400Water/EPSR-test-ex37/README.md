# EPSR-test-ex37 — `makelattice`

**Program:** `makelattice` (`fortran_src/makelattice.f`): builds a crystal of tethered
one-atom molecules from a unit-cell file.

**Inputs (unit-cell files written for this test):**
* `NaCl.cell` — rock salt, a = 5.64 Å, `fractional` coordinates, 8 atoms, 2 atom types,
  `ecore dcore` line and a random-number-state line.
* `tric.cell` — triclinic cell (80°, 100°, 110°: exercises `convertalphabetagamma`),
  `absolute` coordinates, atom type line without `iexchange`, no ecore line, no
  random-number line (generator reset).
* `ortho.cell` — no coordinate-type line (count read from line 3), gamma = 120°.
* `ml1.stdin`, `ml2.stdin`, `ml3.stdin`.

**Commands:**
1. `makelattice NaCl.cell 2 2 2 < ml1.stdin` (output name with extension).
2. `makelattice < ml2.stdin` — file and `na` on one stdin line, `nb nc` prompted.
3. `makelattice < ml3.stdin` — file typed, `na nb nc` prompted.

**Fortran details reproduced:** the box is always written as non-cubic (`cubic` is
never set); `cosd/sind` are the gfortran *intrinsics* (not EPSR's own functions, because
they are not declared `external`); `natomsph` is implicitly INTEGER, so `ctype` stores
the bit pattern of integer 1 (1.4E-45) as `natomsphere` (visible in the `.ato` atom-type
lines and the coarse-grained density 0).

**Expected outputs:** `NaCl_222.ato`, `tric.ato`, `ortho.ato`.

**Verification:** ALL IDENTICAL.
