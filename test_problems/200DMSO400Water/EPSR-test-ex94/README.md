# EPSR-test-ex94 — `nsmoop 0`, `refpotfac 0.8`, `rbroad 0.05`

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.
Inputs made with `scripts/make_case.py`:

* `.EPSR.inp`: `nsmoop 0` (no smooth background subtraction from the data: `ftback` is not called);
* `.pcof`: `refpotfac 0.80000` (reference LJ + Coulomb potential scaled by 0.8 in `generate_rp`), `rbroad 0.05000000`
  (decay `exp(-rbroad*r)` of the Poisson EP basis functions in `calc_pwrexpr_grad`, both in `expon` and the gradient factor
  `n/r-1/sigma-rbroad`). The values are written back to the `.pcof`/`.EPSR.inp` (`rbroad` list-directed: `5.00000007E-02`).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex94 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in all 3 passes (byte-for-byte after masking the `.out` date line and the fort.35 timings).
Branch coverage of pass 1 checked with a gcov build of the Fortran.
