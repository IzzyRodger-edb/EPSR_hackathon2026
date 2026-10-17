# EPSR-test-ex9 — one X-ray data set

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.EPSR.inp`: `iinit 1`; dataset 3 replaced by X-ray data `DMSOWater_xray.mint01` with weights `DMSOWater.XWTS.wts` made by `xwts` from `DMSOWater.XWTS.dat` (and `f0_WaasKirf.dat` form factors). The X-ray data files are copies of the neutron data files (`DMSOWater_xray.mint01` = `D6DMSO_H2O.mint01`): the case tests the X-ray code path, not a physical X-ray fit. Exercises `read_xray_weights`, `wxtype`, `wxcreate` (writes `X-raywts.dat`, `X-rayf0x.dat`), the X-ray branches of `iterate`/`ptgnsq` and the Q interpolation of the data.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex9 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
