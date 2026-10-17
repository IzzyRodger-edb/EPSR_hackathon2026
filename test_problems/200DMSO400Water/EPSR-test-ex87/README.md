# EPSR-test-ex87 — dockato: Qv site as first atom / first atom type

Inputs: `WaterQvA.ato` (created): like ex86 but Qv is atom 1 and the first atom type, so removing it shifts the remaining types (shift branch of `clean_atom_types`; the re-created Qv type gets the stale parameters of the last slot — Fortran behaviour, reproduced).
* run 1: random 20–340°, factor 1.0, ndock 3 (only 1 possible); run 2: rotation only on the docking side (atomvrr = 0), overlaps allowed.
Outputs: `peroxideA.ato`, `peroxideB.ato`.

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
