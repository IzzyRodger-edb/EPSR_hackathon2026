# EPSR-test-ex10 — two X-ray data sets

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

As ex9 but datasets 2 and 3 are X-ray (`D6DMSO_HDO_xray.mint01` + `DMSOWaterA.XWTS.wts`, `DMSOWater_xray.mint01` + `DMSOWaterB.XWTS.wts`), `nxfile = 2` (the `MXFILE=2` limit). The X-ray data files are copies of the neutron data files (`D6DMSO_HDO_xray.mint01` = `D6DMSO_HDO.mint01`, `DMSOWater_xray.mint01` = `D6DMSO_H2O.mint01`): the case tests the X-ray code path, not a physical X-ray fit. `wxcreate` rewrites `X-raywts.dat`/`X-rayf0x.dat` for each X-ray set, so the files of the last set remain.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex10 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
