# EPSR-test-ex27 — `changeato`, navigation and molecule removal

Program: `changeato`. `MixD.ato` is a copy of `EPSR-test/200DMSO400Water.ato`,
`WaterE.ato` of `EPSR-test/WaterMolecule.ato`.

| stdin | what it tests |
|---|---|
| `caseD` (MixD) | `remove`: a blank answer, 3 DMSO molecules from the *start* of the file (the Fortran does not renumber the bond partners / rotation-group atoms of the molecules that move down — reproduced exactly), 395 waters; then keyword navigation with `d`/`u` and bare values, a non-numeric value, `help`, `atom 4`, `atom 99` (clamped), `default`, an unknown command, `save`/`y`, `q` |
| `caseE` (WaterE) | `atlabel` of the H type, epsilon of O, `label` of all H atoms (fraction 1.0), temp, `bye`, save at the final prompt with overwrite |

Verification: **ALL IDENTICAL** (reference in `fortran_reference.tar.gz`).
