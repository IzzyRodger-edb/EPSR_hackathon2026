# EPSR-test-ex5 — complete reset (`ireset 1`)

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.EPSR.inp`: `ireset 1` — EP coefficients reset, `pressreq/ereqstep/erequnit` zeroed, accumulators restarted, default pair potentials; after the first run the Fortran writes `ireset 0` back so passes 2-3 continue normally.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex5 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
