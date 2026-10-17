# EPSR-test-ex111 — three-body potentials: several groups, several angles

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.pcof`: appended after the pair potentials
```
tb O1  O1  O1  20.0 15.0 0.1 3.3 3.3 0 0 0 109.5 1.0 60.0 0.5
tb H1  O1  H1  3.0 5.0 0.1 1.2 1.2 1 1 1 104.5 1.0
some comment line that is ignored
tb O1 H1
tb O1  H1  O   10.0 10.0 0.2 1.2 2.6 1 0 0 180.0 1.0 160.0 -0.5 140.0 0.25 120.0 0.1 100.0 0.05
```
* water O···O···O network angle (all three atoms on different molecules, two requested angles: nearest-angle
  selection, `im1 /= im2` neighbour bookkeeping);
* purely intramolecular H1–O1–H1 angle (same molecule flags 1 1 1: no triplet neighbours);
* lines that are not `tb` with > 13 words are ignored (as in the Fortran);
* water O1–H1···O(DMSO) with five angles (first and third types differ → `ia3 /= ia1` rule); a negative
  amplitude is made positive by `read_pcof` (`abs`), but `save_files_others` re-reads it without `abs`, so it
  is negative in the second iteration and is written back as -0.50.  With 5 angles the Fortran builds the
  format `(...,3(1x,i1),*(1x,f7.2))` (`i1.1` of 10 gives `*`, an unlimited repeat) — same output as `10(...)`.

`.EPSR.inp`: `niter 2`.  `commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3`.
Three-body energy per molecule in the last `.out`: ~712 kJ/mol (second iteration, radians quirk — see ex110).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex111 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).

## Expected divergence from the Fortran (deliberate bug fix)

Since the bug-fix pass (see FORTRAN_BUGS.md at the repo root), this case
deliberately differs from the Fortran reference because of the three-body fix: the in-memory triplet angles/width/amplitudes are no longer overwritten with radians at .pcof save time, so iterations 2+ use the correct cosine convention (Fortran iteration-2 energy +688.9 kJ/mol vs fixed Rust -38.5).
The difference is the fix working; run_ex_test.sh reporting SOME FILES DIFFER
here is expected.
