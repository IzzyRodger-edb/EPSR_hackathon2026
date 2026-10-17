# EPSR-test-ex8 — molecule size factor

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.EPSR.inp`: `sizefactor 2.0 0.9 0.0` — molecules start shrunk by a factor 2 and are expanded by 0.9-steps over the runs (WCA reference potential `generate_mp`, `scalefacmol`, `molfit` rescaling); 4 chained runs.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 4` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex8 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
