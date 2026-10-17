# EPSR-test-ex1 — exponential repulsive terms

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.pcof`: three pairs get `rminex pcoefex coordex` values (S-S 4.0/0/0.1, two water pairs 2.6/0.5/0 and 1.5/0.2/0.05), so the exponential repulsive core of `generate_rp`, `repulsivept` (coordination-number feedback on `pcoefex`) and the corresponding `potex` sums of `atogen` are active.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex1 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
