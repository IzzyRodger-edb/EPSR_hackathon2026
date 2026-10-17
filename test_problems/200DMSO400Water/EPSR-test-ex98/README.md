# EPSR-test-ex98 — revlorch (Q-space Lorch broadening)

Program: `epsr`. Based on EPSR-test; only `200DMSO400Water.EPSR.inp` differs (`scripts/make_inp.sh`): `revlorch 0.05 0.01` (`lorchqwidth`, `lorchqwidthq`).

Exercises the broadening of the r->Q transform matrix `foqols` in `inita`: `bwidth=rr*lorchlim*(lorchqwidth+lorchqwidthq*qq)`, `lorchlim=4.4934/rmaxgr`, `sums=sums*pofr(bwidth)` (function `pofr`, single precision). All simulated S(Q) output (.q01/.u01/.x01/...) changes relative to the base case.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex98 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in all 3 passes (byte-for-byte after masking the `.out` date line and the fort.35 timings).
