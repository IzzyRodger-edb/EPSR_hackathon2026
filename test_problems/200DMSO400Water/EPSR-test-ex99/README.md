# EPSR-test-ex99 — miscellaneous dataset / run options

Program: `epsr`. Based on EPSR-test; only `200DMSO400Water.EPSR.inp` differs (`scripts/make_inp.sh`):

* `qmin 0.60 0.80`: qminm raises `nqmin` from 8 to 13 for all datasets, qminp=0.8 sets the minimum Q of the EP fit in `ptgnsq`;
* `rshmin` 0.5 / 0.9 / 1.1 (background subtraction `ftback` and the minimum radius of the EP perturbation, log `Minimum radius integer of potential: 19 / 33 / 39`);
* `tweak` 1.03 / 1.0 / 0.97 (data scaling);
* `efilereq` 1.0 / 0.5 / 0.0 (`efacm = efilereq/mean(nonzero efilereq)` scales the dataset rows in the inversion matrix (`get_wts`; the inversion is recalculated whenever efacm changes); dataset 3 gets zero weight in the EP);
* `fwhm 0.03` (Q-independent resolution term in `inita`);
* `rspcrmin 1.50`, `rspcfrac 0.3` (R-space coefficient `epcup` / `epccoeff` and the `.erg` trend line);
* `qinscalefc 0.98 0.5` (qscaleinputfactor < 1, so no coarse graining: only read and written back; selfbroadpwr is only used with coarse graining).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex99 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in all 3 passes (byte-for-byte after masking the `.out` date line and the fort.35 timings).
