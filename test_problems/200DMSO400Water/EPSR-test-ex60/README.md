# EPSR-test-ex60 — `clusters`, three distributions, accumulated over 2 runs

Program: `clusters` (`clusters.f`), input `hb.CLUSTERS.dat`, `commands.txt`: `repeat 2` × `clusters hb.CLUSTERS.dat`.

Distributions: (1) water O1–H1 hydrogen bonds 1.4–2.4 Å (19 clusters, 88 % of molecules in percolating clusters), (2) DMSO O – water H1 1.4–2.5 Å (molecules of both types take part), (3) S–O1 3.0–4.0 Å contacts. `nsize 60`, `nsumt 0`, so run 1 starts an accumulation and run 2 reads `hb.CLUSTERS.n01` back and averages (`outgofz` accumulation branch, error column non-zero).

Outputs: `hb.CLUSTERS.n01` (row 1 = percolating fraction, rows 2.. = normalised cluster-size histogram), rewritten `hb.CLUSTERS.dat` (nsumt 1, 2) and `fort.14` — the Fortran writes the "fraction containing 2 or more molecules" lines to unit 14 *after* `save_files_others` has closed it, so they go to `fort.14` (in the working directory; the Rust port writes it to the case folder).

Inputs: `200DMSO400Water.ato` copied unchanged from `EPSR-test`; the keyword file was written by hand (the Fortran program creates no file when the named input file does not exist; keywords/defaults from `setup_symbols_generic.f`).

Run: `rust/tools/run_ex_test.sh <this folder> <fortran build dir> [work dir]` (Fortran references in `fortran_reference.tar.gz`, made with `--save-reference`).

Verification: Rust port (`epsr-rs ./ clusters ...`) vs gfortran reference: **ALL IDENTICAL** on every pass (all written files byte-identical).
