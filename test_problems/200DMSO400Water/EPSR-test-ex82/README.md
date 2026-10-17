# EPSR-test-ex82 — dockato: no / one-sided second rotation

Inputs: `WaterMolecule.ato`.
* run 1: atomdrr = atomvrr = `0`, angles 0 0, `nooverlaps` (factor 1.0): no second rotation; since dminrr = dmaxrr the Fortran still creates DIH groups whose atom pair is (0,0) (it reads `nbond(0)`, the preceding COMMON element, = 0) — reproduced.
* run 2: atomvrr = 0 only: the docking molecule is simply rotated by the random angle (10–170°).
Outputs: `waternorot.ato`, `waterhalf.ato`.

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
