# EPSR-test-ex66 — `lsi`, water local structure index, 3 chained runs

Program: `lsi` (`lsi.f` + `chains_routines.f`), input `water.LSI.dat`, `repeat 3`, `nsumt 0`.

Water O1–H1···O1 bond criterion; LSI range −1..2 Å with 100 bins. Run 1 and 2 use fc = 1 (the Fortran uses the *input* nsumt: 1/nsumt only for nsumt > 0), run 3 fc = 0.5 for the running averages `nbonded`, `nnonbonded`, `avbond` written back to the `.dat` file (format `2(1x,e12.5)`; second value = spread). DMSO molecules have no O1 atom, so their 20 "neighbours" have no matching atom pair (`iwmin = ia2min = 0`); the Fortran then reads `dfs1(id,0)` and `xa(0)` out of bounds — emulated in the port (see lsi.rs).

Outputs: `water.LSI.n01` (10 columns: bonded / non-bonded neighbour counts, LSI, three neighbour-triangle angle distributions, r3..r6), rewritten `water.LSI.dat`.

Inputs: `200DMSO400Water.ato` copied unchanged from `EPSR-test`; the keyword file was written by hand (the Fortran program creates no file when the named input file does not exist; keywords/defaults from `setup_symbols_generic.f`).

Run: `rust/tools/run_ex_test.sh <this folder> <fortran build dir> [work dir]` (Fortran references in `fortran_reference.tar.gz`, made with `--save-reference`).

Verification: Rust port (`epsr-rs ./ lsi ...`) vs gfortran reference: **ALL IDENTICAL** on every pass (all written files byte-identical).
