# EPSR-test-ex81 — dockato: fixed dihedral angle (DIH groups)

Inputs: `WaterMolecule.ato`. All parameters on the command line.
dminrr = dmaxrr → a DIH rotational group with the requested dihedral is created at each bond (exercises the DIH writing added to `moloutc`).
* run 1: 120°, overlaps allowed, 4 molecules; run 2: 75°, overlap factor 0.8 (2 of 4 placements accepted).
Outputs: `waterdih.ato`, `waterdih2.ato`.

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
