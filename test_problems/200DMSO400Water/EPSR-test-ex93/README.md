# EPSR-test-ex93 — Gaussian representation of the empirical potential (`Gaussian T`)

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.
Inputs made with `scripts/make_case.py`:

* `.pcof`: `Gaussian T` (instead of `Poisson`). The EP coefficients of EPSR-test (fitted with the Poisson basis) are kept and are
  reinterpreted in the Gaussian basis, so the EP changes strongly (EP energy ≈ 99 vs 3.7 in the Poisson run).

Exercised (gcov-checked): `inita` (`T` → `gaussian=.true.`, symbol value rewritten as `Gaussian`, so passes 2-3 test the `Gaussian`
spelling); `generate_ep` with `calc_gaussr_grad` (r-space gradients of Gaussians at `d=m*pdstep`, width `psigma2`); `ptgnsq` with
`calc_gaussq` (Q-space basis `sin(Qd)/(Qd)*exp(-(Q*sigma)**2/2)`).

Note: `rbroad`, `dmax` and `rho` (in `calc_gaussq`) are unused by the Gaussian basis routines; `rmin=0` so `ncfirst` is never touched.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex93 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in all 3 passes (byte-for-byte after masking the `.out` date line and the fort.35 timings).
Branch coverage of pass 1 checked with a gcov build of the Fortran.
