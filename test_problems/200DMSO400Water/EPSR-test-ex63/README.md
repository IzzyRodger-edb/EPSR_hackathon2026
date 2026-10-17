# EPSR-test-ex63 — `chains`, default values and nsize cap

Program: `chains`, input `ss.CHAINS.dat`, `repeat 2`.

Single S(–O)···S criterion between DMSO molecules (water molecules have no neighbours and count as chains of length 1). `nsize 300` exceeds `mchain` = 200 and is capped; `nsumt`, `ndist`, angle keywords are not given, so the defaults are used (including the computed `nsize` text "Maximum chain length (max 200)" of `setup_symbols_chains`, visible in the rewritten `.dat`); `nsumt -1` means no accumulation.

Inputs: `200DMSO400Water.ato` copied unchanged from `EPSR-test`; the keyword file was written by hand (the Fortran program creates no file when the named input file does not exist; keywords/defaults from `setup_symbols_generic.f`).

Run: `rust/tools/run_ex_test.sh <this folder> <fortran build dir> [work dir]` (Fortran references in `fortran_reference.tar.gz`, made with `--save-reference`).

Verification: Rust port (`epsr-rs ./ chains ...`) vs gfortran reference: **ALL IDENTICAL** on every pass (all written files byte-identical).
