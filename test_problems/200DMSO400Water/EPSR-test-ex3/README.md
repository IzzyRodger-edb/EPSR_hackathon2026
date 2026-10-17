# EPSR-test-ex3 — unbiased `ereq` stepping

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.pcof`: `ereqstep 2.0 1.0`, `bias 0` — automatic EP amplitude control with unbiased steps, 6 chained runs.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 6` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex3 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
