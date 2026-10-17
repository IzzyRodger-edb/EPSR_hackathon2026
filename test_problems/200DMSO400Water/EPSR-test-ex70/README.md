# EPSR-test-ex70: sharm → mapgr → mapgrentropy, water–water (identical molecules), chained 6 times

Input configuration: `200DMSO400Water.ato` copied from `EPSR-test` (200 DMSO + 400 water).
Input files were written by hand from the keyword tables (`setup_symbols_*`), one `keyword value` per line, terminated by `q`.

Programs: `sharm`, `mapgr`, `mapgrentropy` (`commands.txt`, `repeat 6`).

* `water.SHARM.dat`: water–water spherical-harmonic coefficients from
  `200DMSO400Water.ato` (400 water molecules): centre = mass centre of O1/H1,
  z axis along the H–H midpoint (`axisc1 z 2 3`), second axis from H 2
  (`axisc2 y 2`); `atom-s` identical to `atom-c` (→ `ldistinct = 0`, axis copy
  in `check_symbol`, only l2 ≤ l1 terms); l1,l2 = 0 1 2, l = 0..4,
  `n1step = n2step = 2` (C2v water: n = 0, ±2) → 38 coefficients.
* The calculation is done when `nsumt` is a multiple of 5: pass 1 (nsumt 0,
  new accumulation) and pass 6 (nsumt 5: the previous `.h01` is read back and
  averaged, `fc1 = fc2 = 0.5`); passes 2–5 only update `nsumt` in the input file.
* `water.mapgr.txt`: g(r,ω) for 3 views — radial (grtype 1), spatial (2),
  full orientational (3); radial/spatial/orientational smoothing 1/1/1;
  `ncoeffs`/`ident` taken from `water.SHARM.dat` (`find_plot_par`).
* `water.mapgrentropy.dat`: entropy estimate from 200 random orientations;
  `ran1.dat` (random-number state) is created on pass 1 and carried through
  the chained runs.

Outputs: `water.SHARM.h01`, `water.SHARM-first.xyz`, `water.SHARM-second.xyz`,
`water.SHARM.dat` (updated), `water.mapgr.g01`, `water.mapgrentropy.out`,
`water.mapgrentropy.s01`, `negativegr.txt`, `ran1.dat`.

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
