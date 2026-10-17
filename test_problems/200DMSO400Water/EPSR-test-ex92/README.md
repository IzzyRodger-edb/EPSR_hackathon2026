# EPSR-test-ex92 — exponential repulsive core and container force with gradient limit (`expecf` 2nd value > 0)

Program: `epsr`. Based on EPSR-test; inputs made by `scripts/make_ex15.py ../../../EPSR-test . 4.0 0.5` (as EPSR-test-ex15, container force on Q-O1,
`iinit 1`) followed by `scripts/make_case.py`:

* `.pcof`: `expecf 0.300 1.500`; `.pcof` repulsive-core lines as in EPSR-test-ex1: S-S `4.0 0 0.1`, O1-O1 `2.6 0.5 0`, O1-H1 `1.5 0.2 0.05` (`rminex pcoefex coordex`).

Exercised (gcov-checked): `generate_rp` exponential core with `gradlimit=1.5`: linear term `pcoefex*gradlimit*(1+ddiff-log(gradlimit))`
for `ddiff >= log(gradlimit)` and the `exp(ddiff)` term otherwise; `en_containerforce` exponential branch with `gradlimit=1.5`
(345 linear / 2055 exponential evaluations). (EPSR-test-ex15 only has the `gradlimit=1, loggradlimit=0` default of the container force
and no exponential core.)

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex92 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in all 3 passes (byte-for-byte after masking the `.out` date line and the fort.35 timings).
Branch coverage of pass 1 checked with a gcov build of the Fortran.
