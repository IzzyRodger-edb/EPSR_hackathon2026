# EPSR-test-ex76: sdfcube DMSO–water (distinct), nsumt = 5 without existing cube file

Input configuration: `200DMSO400Water.ato` copied from `EPSR-test` (200 DMSO + 400 water).
Input files were written by hand from the keyword tables (`setup_symbols_*`), one `keyword value` per line, terminated by `q`.

Program: `sdfcube`.

* `dmsowat.SDFCUBE.dat`: water (O1 centre) around DMSO (centre S/O/C,
  axes `z 4`, `x 2 3`), 21×21×21 grid, `averagemol n`, `nsumt 5` with no cube
  file (initcube resets `nsumt` to 0; the saved value is 1).

## Verification

`rust/tools/run_ex_test.sh <this folder> <FORTRAN_BUILD_DIR> <work> ` →
**ALL IDENTICAL** for every pass (Fortran reference = `fortran_reference.tar.gz`).
