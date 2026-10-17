# EPSR-test-ex42 — writexyz, append and numbered-file modes

Program: `writexyz` (`fortran_src/writexyz.f`): writes (part of) a configuration as `.xyz` blocks, optionally
centred on a molecule, restricted to a slab, rotated by Euler angles, with atom types excluded.

Inputs: `200DMSO400Water.ato`, `DMSO.ato` from EPSR-test; `writexyz_dmso.stdin`.

* `writexyz 200DMSO400Water.ato append 5 8 8 -6 6 30 45 60 4` (twice) — centre molecule 5, |x|,|y| <= 8,
  -6 <= z <= 6 after rotation (30,45,60), atom type 4 (DMSO H) excluded; minimum image, `rotatecoords`;
  the second run appends block 00002 to `200DMSO400Waterappend.xyz` (reading of existing blocks).
* `writexyz 200DMSO400Water noappend 0 100 100 -100 100 0 0 0 0` (twice) — whole box, numbered files
  `200DMSO400Water00001.xyz`, `...00002.xyz`, counter file `200DMSO400Water.nblock`.
* `writexyz DMSO.ato < writexyz_dmso.stdin` — every parameter from stdin (append y, centre 1, distances on two
  lines, angles, exclusion list "3 4"; atoms of the central molecule are always written).

Expected outputs: `200DMSO400Waterappend.xyz`, `200DMSO400Water0000[12].xyz`, `200DMSO400Water.nblock`,
`DMSOappend.xyz`.

Verification: `ALL IDENTICAL` against the original Fortran reference binary.
