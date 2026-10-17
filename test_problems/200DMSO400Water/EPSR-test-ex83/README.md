# EPSR-test-ex83 — dockato: water onto DMSO (docked molecule with ROT groups)

Inputs: `DMSO.ato` (EPSR-test) and `WaterOH110.ato` = `WaterMolecule.ato` with the O–H bond entries changed to 1.10 Å (dockato requires equal bond lengths for the two site pairs, DMSO C–H = 1.10 Å).
The water H1 is docked onto DMSO C, its O1 onto the C–H hydrogen: C–O bond, existing DMSO ROT groups are shifted and a new ROT group inserted.
* run 1: site answers from `dmso_random.stdin` (C H H, 30–90°, factor 1.0), 3 molecules;
* run 2: fixed dihedral 180° referenced to S, overlaps allowed → DIH groups;
* run 3: no second rotation, factor 0.9, 6 molecules.
Outputs: `dmsowater.ato`, `dmsodih.ato`, `dmsonorot.ato`.

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
