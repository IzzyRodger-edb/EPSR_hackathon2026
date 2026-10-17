# EPSR-test-ex96 — szeros /= 0 (extrapolation of the data to Q=0)

Program: `epsr`. Based on EPSR-test; only `200DMSO400Water.EPSR.inp` differs (`scripts/make_inp.sh`): `szeros` = -0.30 (data 1), -0.12 (data 2), 0.0 (data 3, mixed case).

Exercises (`epsr_standalone_rev3.f`): `szerostemp = szeros - sum(ciwt*soqintra(1,:))`; the linear extrapolation `conss0+grads0*(q-q(1))-sum` below `nqmin` in the background-subtraction difference `dts` (before `ftback`); in `ptgnsq` the extrapolated `qtemp` below max(qminp/qminm, q(nqmin)) (visible in the first columns of `delfit.dat` for datasets 1 and 2, zero for dataset 3) and the szeros/=0 R-factor branch (1..nq, log only). The `fit_expansion` call is the same in both Fortran branches. `nqmin` = 8 > 1, so the extrapolation branch is taken.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex96 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in all 3 passes (byte-for-byte after masking the `.out` date line and the fort.35 timings).
