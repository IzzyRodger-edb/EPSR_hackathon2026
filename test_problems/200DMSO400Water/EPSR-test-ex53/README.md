# EPSR-test-ex53 — `triangles` (bond-angle distributions)

Program: `triangles` (`fortran_src/triangles.f`). Four inputs, list repeated twice:

* `slabz.TRI.dat`: `atom-z O`, first atom = centre type (z-only distance for 1-2), ltype 1 and 2.
* `cylr.TRI.dat`: `atom-r S`, ltype 1 and 2.
* `centreS.TRI.dat`: `atom-c S` with S as first atom (damin/damax replaced by rminrmaxc), plus a
  distribution with O1 as first atom restricted to middle atoms near S.
* `water.TRI.dat` (last; nsumt 0 -> accumulated): O1-O1-O1 (ltype 1), O1-O1-<H1 bonded to
  middle O1> (ltype 2, bond-averaged third atom), S-O1-H1 (ltype 1).

Outputs compared: `*.TRI.c01`, rewritten `.dat` files, `fort.14` (order parameter per
distribution, list-directed). `fort.14` of the other three inputs was also checked identical
in separate runs.

Note: `triangles.f` calls `sind`/`cosd` without declaring them EXTERNAL, so gfortran 15 uses its
*intrinsic* SIND/COSD (`_gfortran_sind_r4`), not the routines of `auxiliary_routines.f`; the Rust
port reproduces the libgfortran algorithm (`tools/anal_common.rs`).

Verification: ALL IDENTICAL (passes 1 and 2).
