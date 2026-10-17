# EPSR-test-ex28 — `makemole` on the EPSR-test molecules

Program: `makemole` (subroutine `makemole` in `ato_manipulation_rev4.f`).
Inputs: `DMSO.mol` and `WaterMolecule.mol` from `EPSR-test` (`.gmol` type: atoms
with coordinates and neighbours, bonds, angles, `rot` groups, `changelabel`,
`potential`, temperatures, `ecoredcore`, `density`, `pseudofac`). The second
command gives the file name without extension.

Outputs: `DMSO.ato`, `WaterMolecule.ato` (single-molecule configurations with the
bond/angle distance constraints and rotational groups); the `.mol` files are
rewritten (unchanged).

Verification: **ALL IDENTICAL** (reference in `fortran_reference.tar.gz`).
