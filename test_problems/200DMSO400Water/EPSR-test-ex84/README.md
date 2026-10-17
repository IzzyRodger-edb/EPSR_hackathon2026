# EPSR-test-ex84 — dockato: water docked into the 200 DMSO + 400 water box

Inputs: `WaterMolecule.ato`, `200DMSO400Water.ato` (EPSR-test). One molecule per run (ndock = 1, see note).
* run 1: overlap factor 1.0 — hundreds of candidate water sites are rejected by `test_overlap` (minimum image) before one is accepted;
* run 2: factor 0.6; run 3: overlaps allowed, fixed dihedral 60° (DIH group).
Note: with a multi-molecule docked file the Fortran extends `nal(imv)` of the docked molecule although the new atoms are at the end of the atom list, so the molecule atom ranges in the output overlap (Fortran behaviour, reproduced). With ndock > 1 on this box a later failed search sets `imv` to a molecule without molecule type and the Fortran writes `moltemplate(0)` (out of bounds); the port refuses that case, hence ndock = 1.
Outputs: `boxwater.ato`, `boxwater2.ato`, `boxwater3.ato`.

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-exNN <fortran build dir> [work dir] [--save-reference]
```
Each line of `commands.txt` runs `dockato ./ dockato <args> [< answers]`. Command-line words are, in order:
docking .ato file, atomd, atomdr, atomdrr, ndock, radmax, docked .ato file, atomv, atomvr, atomvrr,
dminrr, dmaxrr, overlaps (factor / `yes` = allow all / `nooverlaps` = 1.0), output file. Missing trailing
words are read from standard input. The only output is the new `.ato` file (stored in
`fortran_reference.tar.gz`, made with the gfortran -O2 -fno-tree-vectorize reference build).

## Verification
Rust port (`epsr-rs ... dockato`) vs Fortran reference: **ALL IDENTICAL** (byte-for-byte).
