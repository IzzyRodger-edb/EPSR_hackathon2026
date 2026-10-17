# EPSR-test-ex114 — three-body potentials with large triplet lists

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.pcof`: appended after the pair potentials
```
tb O1  O1  O1  0.5 20.0 0.1 6.0 6.0 0 0 0 109.5 1.0 70.0 0.3
tb H1  O1  O1  1.0 15.0 0.1 5.0 3.0 0 0 0 100.0 1.0
tb O1  S   O1  1.0 15.0 0.1 4.0 4.0 0 0 0 90.0 1.0
```
Long distance ranges, all atoms on different molecules: more than 10000 triplets in a group (the Fortran
reallocates `tripletlist`), more than 4 triplet neighbours per molecule (reallocation of `threebneighbmol`),
DMSO-centred triplets with water neighbours.  `commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2`.

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex114 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
