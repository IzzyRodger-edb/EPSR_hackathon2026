# EPSR-test-ex50 — `partials`, all atom-type pairs, chained accumulation

Program: `partials` (`fortran_src/partials.f`), input `all.PARTIALS.dat`
(keyword file written by hand from the defaults printed by the Fortran program).

* Configuration: `200DMSO400Water.ato` (copy of EPSR-test), 6 atom types -> 21 pairs.
* No centre atom (`atom-c/z/r = *`): intermolecular g(r) of every pair, nz = 300, `rnormalise 1.0`.
* `nsumt 0` and `repeat 2`: the second run reads `all.PARTIALS.g01` back (r_int) and
  forms the running average and spread (outgofz).

Outputs compared: `all.PARTIALS.g01`, the rewritten `all.PARTIALS.dat` (nz, nsumt, ndist
updated), `fort.14` (pair list — written to unit 14 after `save_files_others` has closed it).

Verification: `run_ex_test.sh` -> ALL IDENTICAL (passes 1 and 2).
