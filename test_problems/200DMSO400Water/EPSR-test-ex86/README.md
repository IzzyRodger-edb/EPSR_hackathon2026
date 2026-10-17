# EPSR-test-ex86 — dockato: Qv bonding sites

Inputs: `WaterQv.ato` (created): water with an extra dummy atom `Qv` bonded to O1 at 1.45 Å (Qv atoms are invisible, qtype 1), plus its atom type.
The Qv of the docking molecule is placed on the docked O1, the docking O1 is rotated onto the docked Qv, both Qv atoms are removed (O–O bond of 1.45 Å). The Qv type disappears (`clean_atom_types`); in run 1 the second docking finds no site (restore path), but its `ftype` call has re-created the Qv type, which stays in the output (Fortran behaviour).
* run 1: random 0–360°, factor 1.0, ndock 2; run 2: fixed dihedral 110° (DIH), `nooverlaps`, answers from `fixed.stdin`.
Outputs: `peroxide.ato`, `peroxide_dih.ato`.

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
