# EPSR-test-ex75: sdfcube water–water (identical), chained 6 times

Input configuration: `200DMSO400Water.ato` copied from `EPSR-test` (200 DMSO + 400 water).
Input files were written by hand from the keyword tables (`setup_symbols_*`), one `keyword value` per line, terminated by `q`.

Program: `sdfcube`; `repeat 6`.

* `water.SDFCUBE.dat`: SDF of water around water on a 17×17×17 grid
  (`mpixel 8`, `rmax 6`).  Pass 1 creates `water.SDFCUBE.cube`; pass 6
  (nsumt 5) reads the step size and the old cube (`initcube`, `outcube`) and
  accumulates (`fc1*old + fc*new`; the Fortran's uninitialised `fc1` is folded
  to 1.0 by gfortran — reproduced).
* For identical molecules the Fortran reads `averagemol` from the `axiss1`
  value (missing `iref=iref+2`), i.e. here a single molecule is written to the
  xyz files — reproduced.

## Verification

`rust/tools/run_ex_test.sh <this folder> <FORTRAN_BUILD_DIR> <work> ` →
**ALL IDENTICAL** for every pass (Fortran reference = `fortran_reference.tar.gz`).
