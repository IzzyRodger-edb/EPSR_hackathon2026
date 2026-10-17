# EPSR-test-ex85 — dockato: water docked onto a DMSO of the box

Inputs: `WaterOH110.ato` (see ex83), `200DMSO400Water.ato`.
* run 1: C–H site of the box DMSO, random 0–360°, factor 1.0; run 2: overlaps allowed, fixed dihedral 90° referenced to S.
Exercises insertion of the new ROT/DIH group among the 400 existing groups of the box.
Outputs: `boxdmso.ato`, `boxdmso2.ato`.

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
