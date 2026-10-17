# EPSR-test-ex91 — harmonic repulsive core and harmonic container force, with gradient limit

Program: `epsr`. Based on EPSR-test; inputs made by `scripts/make_ex15.py ../../../EPSR-test . 4.0 0.5` (as EPSR-test-ex15: Q site at every DMSO S,
28 pairs, Q-O1 container force `rmaxex=4.0 pcoefexmax=0.5`, `iinit 1`) followed by `scripts/make_case.py`:

* `.pcof`: `reppottype harmonic`, `expecf 0.300 1.500` (gradient limit 1.5); `.pcof` repulsive-core lines as in EPSR-test-ex1: S-S `4.0 0 0.1`, O1-O1 `2.6 0.5 0`, O1-H1 `1.5 0.2 0.05` (`rminex pcoefex coordex`).

Exercised (gcov-checked): `generate_rp` harmonic core, both the linear (`ddiff >= gradlimit`: `0.5*pcoefex*gradlimit*(2*ddiff-gradlimit)`) and
the quadratic branch; `en_containerforce` harmonic branch (`rtest=rmaxex-expecfuse`), linear (`ddiff > gradlimit`, 281 calls) and quadratic
(1580 calls) parts, in `fmolec1` / `en_change`.

Fortran quirk (reproduced): the harmonic core uses `ddiff >= gradlimit`, the harmonic container force `ddiff > gradlimit`; for
`expecfpwr <= 0` the container force uses `gradlimit=1` but the core uses `gradlimit=0` (no limit).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex91 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in all 3 passes (byte-for-byte after masking the `.out` date line and the fort.35 timings).
Branch coverage of pass 1 checked with a gcov build of the Fortran.
