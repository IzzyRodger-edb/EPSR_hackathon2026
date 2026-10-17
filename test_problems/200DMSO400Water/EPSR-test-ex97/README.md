# EPSR-test-ex97 — Porod scattering (porodpwr, porodadd)

Program: `epsr`. Based on EPSR-test; only `200DMSO400Water.EPSR.inp` differs (`scripts/make_inp.sh`): `porodpwr 0.5 3.5` (Porod Q and power), `porodadd` = 0.25, 0.10, -0.5.

Exercises: search of `nporodq` (Q index nearest to porodq, here 11), `porodamp(in)=datai(nporodq,in)+porodadd(in)` (log: 0.05630, 0.08910, -0.36131 — dataset 3 negative, so no Porod term for it), the Porod term `porodamp*(porodq/q)**porodpwr` added to the background-subtraction difference `dts` and to the simulated fit `soqism`/`ufile` (.u01, R-factors).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex97 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in all 3 passes (byte-for-byte after masking the `.out` date line and the fort.35 timings).
