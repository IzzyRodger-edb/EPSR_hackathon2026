# EPSR-test-ex52 — `coord` (coordination-number distributions)

Program: `coord` (`fortran_src/coord.f`). Four inputs, list repeated twice (chained runs):

* `hbond.COORD.dat`: 3 distributions (O1-H1 1.4-2.4, S-O1 2.5-4.5, O-H1 1.2-2.5 A), no centre
  atom, nsumt 0 -> accumulation over the two passes.
* `slabz.COORD.dat`: `atom-z O1`, O1 atoms within |z| in [0,1) A of an O1 atom, nsumt -1.
* `cylr.COORD.dat`: `atom-r S`, xy-distance window around S.
* `centreS.COORD.dat` (last, so its `fort.14` is the one kept): `atom-c S`, two distributions
  with `rminrmaxc` windows; the `fort.14` line uses format reversion (8 values, 4(1x,f10.5)).

Outputs compared: `*.COORD.n01` (histogram + mean bond length columns), rewritten `.dat`
files, `fort.14` (mean coordination number and rms). The `fort.14` of `hbond` (ltypc = 0
branch) was also checked identical in a separate single-command run.

Verification: ALL IDENTICAL (passes 1 and 2).
