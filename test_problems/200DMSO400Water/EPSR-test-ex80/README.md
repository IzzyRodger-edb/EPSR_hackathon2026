# EPSR-test-ex80 — dockato: water onto water, answers from standard input

Inputs: `WaterMolecule.ato` (copied from EPSR-test).
* run 1: docking parameters on the command line, docked-site parameters (O1 H1 H1, dihedral range 30–90°, overlap factor 1.0, output name) from `dock_random.stdin`. The docking H1 is placed on the docked O1, the docking O1 is rotated onto the docked H1, both H's removed and an O–O bond formed; random second (dihedral) rotation using ran1; overlap test rejects some placements (2 of 3 added). New ROT groups, box size rescaled.
* run 2: overlaps allowed (`yes`), dihedral range −60–60°, blank line before the output name (re-prompt), name given with `.ato` extension; 6 molecules chained.
Outputs: `waterchain.ato`, `waterchain_yes.ato`.

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
