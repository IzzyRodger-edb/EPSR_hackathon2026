# EPSR-test-ex61 — `clusters`, overflow bin and no accumulation

Program: `clusters`, input `oo.CLUSTERS.dat`, `repeat 2`.

Distributions: (1) water O1–O1 2.5–3.3 Å, (2) DMSO C – water O1 3.0–3.6 Å. `nsize 10` is smaller than the largest clusters, so the overflow branch (`gofz(nz,id)` collects all larger clusters) is used; `nsumt` is not given (default −1) so no accumulation takes place (fc1 = 1, nsumt stays −1 in the rewritten file). Outputs as in ex60 (`oo.CLUSTERS.n01`, `oo.CLUSTERS.dat`, `fort.14`).

Inputs: `200DMSO400Water.ato` copied unchanged from `EPSR-test`; the keyword file was written by hand (the Fortran program creates no file when the named input file does not exist; keywords/defaults from `setup_symbols_generic.f`).

Run: `rust/tools/run_ex_test.sh <this folder> <fortran build dir> [work dir]` (Fortran references in `fortran_reference.tar.gz`, made with `--save-reference`).

Verification: Rust port (`epsr-rs ./ clusters ...`) vs gfortran reference: **ALL IDENTICAL** on every pass (all written files byte-identical).
