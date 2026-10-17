# EPSR-test-ex77: sharm DMSO–DMSO with high l (≤ 6, L ≤ 12) → mapgr → mapgrentropy

Input configuration: `200DMSO400Water.ato` copied from `EPSR-test` (200 DMSO + 400 water).
Input files were written by hand from the keyword tables (`setup_symbols_*`), one `keyword value` per line, terminated by `q`.

Programs: `sharm`, `mapgr`, `mapgrentropy`.

* `dmso.SHARM.dat`: identical DMSO molecules (centre = S, C, C, O mass
  centre, axes `-z 4`, `y 2`), l1,l2 = 0 2 4 6, l = 0..12 (even),
  n1step = n2step = 3 → 181 coefficients (`.h01`, `.h02`): large
  log-factorial arguments in the Clebsch–Gordan coefficients and rotation
  matrices (precision test of `cgcoef`/`dlm`).
* `dmso.mapgr.txt`: 3 views (radial, two orientational), smoothing
  `3.0 2 1 1`; `dmso.mapgrentropy.dat`: 40 trials.

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
