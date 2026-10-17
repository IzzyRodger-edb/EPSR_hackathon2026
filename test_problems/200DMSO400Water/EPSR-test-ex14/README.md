# EPSR-test-ex14 — non-cubic (triclinic) box

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.ato`: the cubic header is replaced by the 9-number cell header (a, b, c = 33.1213 Å and the cell angles in EPSR's `phib thetac phic` convention, 80°, 10°, 20°), so `minimagec` uses the general-cell branch and `Lattice` / `sync_lattice` the non-orthogonal metric; the atomic number density becomes 0.0908 Å⁻³. All other inputs are the EPSR-test files.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex14 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
