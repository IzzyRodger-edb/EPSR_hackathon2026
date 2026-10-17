# EPSR-test-ex11 — additional Gaussian potentials

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.pcof`: Gaussian `(r, w, h)` triples on three pairs — two Gaussians on one pair, and a *negative* r on another, which holds the EP of that pair fixed (`efacp = 0`). Exercises the Gaussian branch of `read_pcof`/`generate_rp` and the `(3(1x,f10.5),N(1x,f9.4,1x,f9.4,1x,e11.4))` output format.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex11 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
