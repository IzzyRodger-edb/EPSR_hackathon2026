# EPSR-test-ex108 — atom pairs missing from the `.pcof` file; `ncoeffp` given in the `.EPSR.inp`

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets). Inputs made by
`scripts/make_inputs.py`: the pairs S-H1 and O1-H1 and the `ncoeffp` header line are removed from the
`.pcof` (19 pairs left); `ncoeffp 60` is added to the `.EPSR.inp` before `fnameato`/`fnamepcof`;
`ntimes 2`, `niter 2`. `commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2`.

* `ncoeffp` in the input file: `check_symbol` reads it into `ncoeffp` (capping at mcoeff = 1000 is tested in
  ex109). Since the `.pcof` has no `ncoeffp` line, `read_pcof` uses this value: only 60 EP coefficients per
  pair are read (the rest of each record is skipped, coefficients 61..116 stay 0); `inita` then resets
  `ncoeffp` to 116 from rminpt/rmaxpt and the `.pcof` is written with `ncoeffp 116`. The keyword is not
  written back to the `.EPSR.inp` (it belongs to the `.pcof` part of the symbol table).
* Missing pairs: `read_pcof` creates a new set of coefficients for S-H1 and O1-H1 (`npairspcof` 20, 21):
  `pcoefp = (1, 0, ...)`, `rminex = max(rminfac*dcor, roverlap)` (0 here), new pair symbol with a
  list-directed `rminex`.

Fortran quirks reproduced: for a new pair `coulfac`, `coordex` and `efacp` are not set (COMMON, i.e. 0):
the pair has **no Coulomb term** (`coulfac 0.00000` is then written to the `.pcof` and read back in the next
run) and its EP is **not refined** in the run that creates it (efacp = 0); it is refined from pass 2 on.

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex108 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).
Inputs are generated from the base EPSR-test files by `scripts/` (see `scripts/make_inputs.*`).

## Verification
gcov: `npairspcof=npairspcof+1` executed twice, `check_symbol` ncoeffp block once.
Rust vs Fortran: **ALL IDENTICAL** in both passes.

## Expected divergence from the Fortran (deliberate bug fix)

Since the bug-fix pass (see FORTRAN_BUGS.md at the repo root), this case
deliberately differs from the Fortran reference because of the new-pair fix: a pair added at save time now gets coulfac/efacp/coordex defaults (1.0/1.0/0.0) instead of unset zeros, so later passes include its Coulomb term and refine its EP.
The difference is the fix working; run_ex_test.sh reporting SOME FILES DIFFER
here is expected.
