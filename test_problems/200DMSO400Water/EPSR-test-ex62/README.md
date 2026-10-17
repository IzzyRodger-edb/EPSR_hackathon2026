# EPSR-test-ex62 — `chains`, H-bond chains, multi-criterion distribution

Program: `chains` (`chains_sp_criterion.f` + `chains_routines.f`), input `hbchain.CHAINS.dat`, `repeat 2` with `nsumt 0` (accumulation over the two runs).

Distribution 1: water–water H-bond O1(–H1)···O1 with O–H 0.8–1.2, O···O 2.5–3.3, H···O 1.5–2.4 Å and H–O···O angle 0–30° (angle test active). Distribution 2: two criteria (`iw` loop) — `O1 H1 O*` (`*` wildcard: matches DMSO O and water O1) and `O1 H1 S`, with angle ranges 0–35° (tested) and 0–180° (angle test skipped because the cosine range is ≥ 1.9). Exercises the bond search, `shortest_path`, the unimodal shortest-path chain search, histograms of neighbour numbers and chain lengths.

Outputs: `hbchain.CHAINS.n01`, rewritten `hbchain.CHAINS.dat`, `fort.14` (mean and spread of each histogram; written after unit 14 was closed, see ex60).

Inputs: `200DMSO400Water.ato` copied unchanged from `EPSR-test`; the keyword file was written by hand (the Fortran program creates no file when the named input file does not exist; keywords/defaults from `setup_symbols_generic.f`).

Run: `rust/tools/run_ex_test.sh <this folder> <fortran build dir> [work dir]` (Fortran references in `fortran_reference.tar.gz`, made with `--save-reference`).

Verification: Rust port (`epsr-rs ./ chains ...`) vs gfortran reference: **ALL IDENTICAL** on every pass (all written files byte-identical).
