# EPSR-test-ex67 — `lsi`, two distributions, wildcards, rtrimax, failed restart

Program: `lsi`, input `mix.LSI.dat`, `repeat 2`.

Distribution 1: two criteria `O1 H1 O1` (water donor) and `O* S H1` (DMSO/water O with S in the same molecule, to water H1), angle limits 30° and 60°, `rtrimax 3.5`, initial `nbonded 1.5 0.5` (two values: spread squared on input) and `nnonbonded 0.3` (one value). Distribution 2: S–C···O1 (water molecules have no S: `iwmin = 0` path), `rtrimax 4.5` (the last value read is the one used, as in the Fortran). `nsumt 2` but no `mix.LSI.n01` exists, so `outgofz` fails to read it and sets nsumt = −1 (no accumulation thereafter); fc = 0.5 in run 1. Because `find_lsi` zeroes the whole `gofz` array at the start of every distribution, the columns of distribution 1 in `mix.LSI.n01` are all zero — a quirk of the Fortran reproduced by the port.

Inputs: `200DMSO400Water.ato` copied unchanged from `EPSR-test`; the keyword file was written by hand (the Fortran program creates no file when the named input file does not exist; keywords/defaults from `setup_symbols_generic.f`).

Run: `rust/tools/run_ex_test.sh <this folder> <fortran build dir> [work dir]` (Fortran references in `fortran_reference.tar.gz`, made with `--save-reference`).

Verification: Rust port (`epsr-rs ./ lsi ...`) vs gfortran reference: **ALL IDENTICAL** on every pass (all written files byte-identical).
