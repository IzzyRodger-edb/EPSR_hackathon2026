# EPSR-test-ex64 — `rings`, water H-bond rings, two distributions

Program: `rings` (`rings_sp_criterion.f` + `chains_routines.f`), input `hbring.RINGS.dat`, `repeat 2` with `nsumt 0` (accumulation).

Distribution 1: water O1–H1···O1 H-bonds with angle ≤ 40°; distribution 2: two criteria `O1 H1 O*` and `O1 H1 S`. `nsize 12` (so `nspmax = 7`). Exercises the shortest-path ring search with the "no short cut" test (`nspvaluet` route counts), the shared `spmin` bookkeeping, and the histograms of neighbour numbers and ring sizes (41 and 90 rings found).

Outputs: `hbring.RINGS.n01`, rewritten `hbring.RINGS.dat`, `fort.14`.

Inputs: `200DMSO400Water.ato` copied unchanged from `EPSR-test`; the keyword file was written by hand (the Fortran program creates no file when the named input file does not exist; keywords/defaults from `setup_symbols_generic.f`).

Run: `rust/tools/run_ex_test.sh <this folder> <fortran build dir> [work dir]` (Fortran references in `fortran_reference.tar.gz`, made with `--save-reference`).

Verification: Rust port (`epsr-rs ./ rings ...`) vs gfortran reference: **ALL IDENTICAL** on every pass (all written files byte-identical).
