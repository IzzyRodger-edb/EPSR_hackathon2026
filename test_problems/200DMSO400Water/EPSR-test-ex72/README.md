# EPSR-test-ex72: sharm with nsumt = 5 (no previous file), the orientation test, mapgr and mapgrentropy

Input configuration: `200DMSO400Water.ato` copied from `EPSR-test` (200 DMSO + 400 water).
Input files were written by hand from the keyword tables (`setup_symbols_*`), one `keyword value` per line, terminated by `q`.

Programs: `sharm` (2 inputs), `mapgr`, `mapgrentropy`; `repeat 2`.

* `watdmso.SHARM.dat`: central water, second DMSO (molecule 2 = DMSO occurs
  first in the .ato, so the "molecule 2 at origin" pass is used), axes
  `z,x` and `-z,x`, l2 = 0 2 3, l = 0..5, n1step 2, n2step 1, `rmax 30`
  (> half box: capped to 0.5·box), `averagemol n`, and `nsumt 5` with no
  existing `.h01` (outgofzs_sharm: read fails → `nsum = -1`, no averaging).
  Pass 2 (nsumt 6) skips the calculation.
* `orient.SHARM.dat`: `test 10 20 30 40 50 60 70 80` — the orientation test
  (fixed phi/theta/chi for both molecules and the lab vector; nsumt is not
  incremented, so pass 2 recomputes), axes `y,z` and `y,-z`, extra blanks in
  `atom-s` (reformatted by `check_symbol`).
* `watdmso.mapgr.txt`: views with grtype 3 and 2, no smoothing.
* `orient.mapgrentropy.dat`: on the test coefficients, `ident -1`,
  `mirrorsym y`, smoothing `1.0 0 1 2`; `ran1.dat` carried to pass 2.

## Fortran reference build (important for `sharm` / `sdf`)

With gfortran >= 14 (the reference uses 15.2) the original `sharm` and `sdf`
stop with `Fortran runtime error: Missing comma between descriptors` when they
write the coefficient column titles: `run_sharm`/`run_sdf` build the run-time
format `'(1h('//il1fmt//',1h,,...'`, i.e. `(1h(i1,1h,,i1,...)`, with no comma
after the Hollerith `1h(`; older gfortran accepted it.  The original programs
therefore only write the `-first.xyz`/`-second.xyz` files and never the `.h01`
coefficients or the updated input file, so nothing downstream (mapgr,
mapgrentropy) can run.

`fortran_reference.tar.gz` was therefore generated with `sharm`/`sdf` rebuilt
from copies of `sharm_rev2.f`/`sdf_rev2.f` in which only that format is changed to
`'(1h(,'//il1fmt//...` (one inserted comma; same output as older compilers),
with the same flags as `build_fortran_reference.sh`; all other programs are the
unmodified reference binaries:

```
sed "s/formatout='(1h('\/\/il1fmt/formatout='(1h(,'\/\/il1fmt/" sharm_rev2.f > patched/sharm_rev2.f   # same for sdf_rev2.f
gfortran -c -O2 -fno-backslash -ffixed-line-length-none -Wno-unused-label -fopenmp -w -fno-tree-vectorize sharm_rev2.f
gfortran -fopenmp -o sharm <objects of the sharm_rev2 line of fortran_link_lists.txt>
```

The Rust port writes the titles `(l1,l2,l;n1,n2)` the intended way.  Run the
comparison with a FORTRAN_BUILD_DIR that holds these patched `sharm`/`sdf`
(with the unpatched binaries the case reports differences for every file after
the xyz files).

## Masking rule (`mapgr`)

`mapgrquad.f` never sets `grwerr(i,1)` (the "error" column of the first,
`g_centres(r)` column of `<run>.mapgr.g01`), so the Fortran writes
uninitialised stack contents there (denormals, NaN, ...); the Rust port writes
0.  `compare_outputs.py` masks characters 31-45 of the data lines of
`*.mapgr.g01` files (that column only).

## Verification

`rust/tools/run_ex_test.sh <this folder> <FORTRAN_BUILD_DIR> <work> ` →
**ALL IDENTICAL** for every pass (Fortran reference = `fortran_reference.tar.gz`).
