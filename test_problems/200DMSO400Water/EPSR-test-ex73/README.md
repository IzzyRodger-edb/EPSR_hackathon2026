# EPSR-test-ex73: sdf water–water (identical), chained 6 times, with mapgr

Input configuration: `200DMSO400Water.ato` copied from `EPSR-test` (200 DMSO + 400 water).
Input files were written by hand from the keyword tables (`setup_symbols_*`), one `keyword value` per line, terminated by `q`.

Programs: `sdf`, `mapgr`; `repeat 6` (calculation in pass 1 and pass 6
with accumulation).

* `water.SDF.dat`: "molecules" found by atom types and distances: O1 H1 H1
  with pair distances 0.976 0.976 1.54 Å ±20 %; the second-molecule
  distances/axes in the file are overwritten by `check_symbol` (atom-s =
  atom-c).  Note: the Fortran checks the remaining (H–H) distance with
  `zdiff = ytemp(i2) - ztemp(i1)`, so only 235 of the 400 waters are accepted
  — reproduced by the port.
* `water.mapgr.txt`: two full orientational views.

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
