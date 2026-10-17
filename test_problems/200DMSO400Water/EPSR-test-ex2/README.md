# EPSR-test-ex2 — automatic `ereq` control with a ceiling

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.pcof`: `ereqmax 14.0`, `ereqstep 1.0 1.0` — the EP amplitude is stepped automatically (`automatic_ereq`) over 6 chained runs (10 → 12); the `ereqmax` ceiling is checked every iteration but not reached.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 6` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex2 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
