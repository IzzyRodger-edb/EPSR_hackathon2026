# EPSR-test-ex122 — diffuse scattering (`diffuse`, `sharm_diffuse`)

**Program:** `epsr`: `sharm_diffuse` / `initialise_diffuse_values` (`ato_manipulation_rev4.f`),
called after the Bragg-peak calculation (Rust `epsr-core/src/bragg.rs`).

**System:** the DMSO/water crystal of EPSR-test-ex120 (orthorhombic box, 5 x 4 x 4 cubic
6.84 A unit cells, 2560 atoms, 21 atom-type pairs). Inputs = outputs of a pre-run of the
Fortran `epsr` with the keywords below and `iinit 1` (accumulators, `.ato`, `.erg` and the
`<run>dNN.DIFFUSE.h01` files), then `nsumt 1`, keeping only the `d01` and `d02` coefficient
files. `scripts/make_case.sh FORTRAN_BUILD_DIR OUT` regenerates all inputs (verified
byte-identical with the reference build described below).

**Keywords (from EPSR-test):** `ntimes 2`, `qwidthqmax 0.01 2.2`, `mplicities 5 4 4 1 1 0`
(only (h,k,0) reflections), `hklqmin 0.0 10.0 0.1`, `diffuse 2 2 1 4 2 1`
(non-Bragg = h or k odd; l = 0, 2, 4, m = 0..l), `nsumt 1`, `iinit 0`.

**What is exercised:** `diffuse` keyword parsing; `bragginten` keeping the hkl arrays; the
spherical-harmonic coefficients (`dlm`, intrinsic `cosd`/`sind`, `atan2d`, inversion image
for h > 0), averaging per Q bin, `smood`, the titles, and `outgofzs_sharm` accumulation:
with `nsumt 1` the `d01` and `d02` files are read and averaged (nsums 1 -> 2 -> 3), the
missing `d03` file sets `nsums = -1`, which switches the running averages of the whole run to
the `nesum` mode (Fortran quirk: `nsums` is the COMMON variable passed to `outgofzs_sharm`).
Pass 2 starts with `nsumt -1`.

**Limits (Fortran behaviour):**
* `sharm_diffuse` indexes `intenhkl(nhkl,npairs)` as `intenhkl(ip,ib)`; this stays inside the
  array only while nhkl <= ~npairs. qbraggmax is chosen so that nhkl = 15 <= 21 pairs (with
  more reflections the Fortran reads/writes past the array: abort/segfault). The Rust port
  reads the same (transposed) storage element and reports the out-of-range case as unsupported.
* The original source builds the title FORMAT `(1h(i1,1h,,...` without a comma after `1h(`;
  with gfortran >= 14 the original `epsr` stops at that statement ("Fortran runtime error:
  Missing comma between descriptors") after writing `.bragg`/`.inten`. As for `sharm` and
  `sdf`, the reference is built with the comma inserted (`tools/build_fortran_reference.sh`
  now patches `ato_manipulation_rev4.f`; the patch only affects `sharm_diffuse`, called by
  `epsr` only). `fortran_reference.tar.gz` was made with that build.

**Verification:** ALL IDENTICAL (both passes) against the patched reference build. Against
an unpatched build (e.g. the old fref3) the Fortran stops in `sharm_diffuse` in pass 1.
