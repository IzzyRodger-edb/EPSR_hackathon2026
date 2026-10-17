# EPSR-test-ex43 — writexyz with pre-existing output files

Program: `writexyz` (`fortran_src/writexyz.f`).

Inputs: `DMSO.ato`, `WaterMolecule.ato`, `200DMSO400Water.ato` (EPSR-test) and pre-existing output files:
`DMSOappend.xyz` (one complete block followed by an incomplete one), `WaterMoleculeappend.xyz` (two complete
blocks plus extra blank lines), `200DMSO400Water.nblock` containing 7.

* `writexyz DMSO.ato append 1 5 5 -5 5 0 90 0 0` — the last block of the existing file is incomplete, so the
  Fortran rewinds: the file is rewritten with block 00001 only.
* `writexyz WaterMolecule.ato append 0 50 50 -50 50 10 0 0 0` — appended as block 00003 after the existing
  blocks (the trailing blank lines are kept).
* `writexyz 200DMSO400Water.ato noappend 999 6 7 -8 9 -20 30 -40 4 abc 3` — counter 7 -> file
  `200DMSO400Water00008.xyz`; centre molecule beyond nmol (centre = origin); the exclusion list stops at the
  non-integer word `abc` (type 4 excluded, 3 not).

Verification: `ALL IDENTICAL` against the original Fortran reference binary.
