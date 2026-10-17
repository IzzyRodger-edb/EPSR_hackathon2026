# EPSR-test-ex13 — tethered molecules and dihedral groups

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.ato`: the first 20 DMSO molecules are tethered to atom 1 (`T00001`, tether position = initial centre), tolerance `tol 0.5`; the methyl `ROT` groups of all 200 DMSO molecules are turned into `DIH` groups (dihedral 180°; no `ROT` groups remain). Exercises `minitc`/`moloutc` for tethered molecules and DIH groups, `en_tethered`, `check_tethered` (move rejection), `en_dihedral`/`en_mole_dihedral`, `ftethered` in the `.out` file and wrapping of tethered molecules.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex13 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
