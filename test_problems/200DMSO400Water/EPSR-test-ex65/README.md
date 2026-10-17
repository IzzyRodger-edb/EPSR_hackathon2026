# EPSR-test-ex65 — `rings`, default nsize

Program: `rings`, input `ow.RINGS.dat`, `repeat 2`, `nsumt 0`.

Water O1–H1···O* (water or DMSO acceptor) rings; `nsize` and `ndist` not given, so the computed default `nsize` = `mring` = 25 ("Maximum ring length (max 25)") of `setup_symbols_rings` is used (59 rings found).

Inputs: `200DMSO400Water.ato` copied unchanged from `EPSR-test`; the keyword file was written by hand (the Fortran program creates no file when the named input file does not exist; keywords/defaults from `setup_symbols_generic.f`).

Run: `rust/tools/run_ex_test.sh <this folder> <fortran build dir> [work dir]` (Fortran references in `fortran_reference.tar.gz`, made with `--save-reference`).

Verification: Rust port (`epsr-rs ./ rings ...`) vs gfortran reference: **ALL IDENTICAL** on every pass (all written files byte-identical).
