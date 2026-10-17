# EPSR-test-ex90 — harmonic repulsive core (`reppottype harmonic`)

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.
Inputs made with `scripts/make_case.py` (keyword/line edits of the EPSR-test files):

* `.pcof`: `reppottype harmonic`; `.pcof` repulsive-core lines as in EPSR-test-ex1: S-S `4.0 0 0.1`, O1-O1 `2.6 0.5 0`, O1-H1 `1.5 0.2 0.05` (`rminex pcoefex coordex`).

Exercised: `inita` sets `harmonic=.true.` (the symbol value stays `harmonic` and is written back); `generate_rp` builds the
Gaussian/harmonic repulsive core `0.5*pcoefex*ddiff**2`, `ddiff=(rminex+expecfuse-r)/expecfuse`, `expecfuse=expecf*sqrt(rminex)`
(quadratic branch only: gradient limit `expecf` second value = 0); `repulsivept` drives `pcoefex` of S-S from 0 (coordination feedback).
Results differ from EPSR-test-ex1 (same pcof with the exponential core).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex90 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in all 3 passes (byte-for-byte after masking the `.out` date line and the fort.35 timings).
Branch coverage of pass 1 checked with a gcov build of the Fortran.
