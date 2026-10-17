# EPSR-test-ex105 — pressure control by box size (`control-p abc` and `control-p a`)

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets).

Two runs in the same folder (`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `epsr boxa.EPSR.inp`, `repeat 2`):

* `200DMSO400Water`: `control-p abc 0.0 0.01`, `ntimes 2`, `niter 2`. After every iteration
  (`runepsr17`, "Adjust boxsize to control pressure") `rat = 0.01*clip(pressd - 0.0, ±1)`,
  `dmaxfac = (1+rat)**(1/3)` and `change_boxsize` scales a, b and c: the molecule centres (and
  `phix/phiy/phiz` of *every* molecule) are converted to fractional coordinates of the old cell
  (`generate_fractional_coordinates`, with minimum image), the lattice is rebuilt with the gfortran
  COSD/SIND intrinsics, the centres converted back, `dmax` = min. perpendicular face distance,
  `cell = dmax`; then `uniform` recomputes `ycoh`. The box stays cubic (33.1213 -> 33.4521 Å after 4 iterations).
* `boxa` (made by `scripts/clone_run.py`): same, with `control-p a 0.0 0.01` (`npcontrols = 1`):
  only a is scaled, the box becomes non-cubic (`cubic = .false.`), the `.ato` is written with the
  9-number cell header and pass 2 reads it back through the general-cell branch of `minitc`.

Fortran quirks reproduced (`epsr-core/src/pressure.rs`): `rho` is not recomputed after the box change
(only on the next read of the `.ato`), so `uniform` recomputes the same `ycoh`; `xstar = x/dmax` (r grid
in box units) is not recomputed either; phix/phiy/phiz of non-tethered molecules are transformed too; the
keyword list accepts `ca` but the branch tests `ac`, so `ca` selects no axis and `ac` is rejected.

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex105 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).
Inputs are generated from the base EPSR-test files by `scripts/` (see `scripts/make_inputs.*`).

## Verification
gcov of the Fortran: the `dmaxfac`, `change_boxsize` and `uniform` lines of the pressure-control block are
executed in every iteration. Rust port vs Fortran reference: **ALL IDENTICAL** in both passes.
