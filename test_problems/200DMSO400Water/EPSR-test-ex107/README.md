# EPSR-test-ex107 — `ref_intra` (EP and reference potential on the intramolecular structure) + `control-p bc`

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets).

`.EPSR.inp`: `ref_intra 0.5 0.3`, `control-p bc 0.0 0.01`, `ntimes 2`, `niter 2`.
`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2`.

* `refine_intra_fac = 0.5`: `en_atom_intra` adds `0.5*(xavs + xavs1*dtdif)` (the empirical potential) to the
  harmonic energy of bonded pairs (when `iz <= npt2` and both atoms are not Q-type) and to the
  intramolecular potential of non-bonded pairs.
* `coul_intra_fac = 0.3`: `generate_rp` replaces the intramolecular potential `xav` of every pair by
  `0.3*ecoul` (full reference potential: LJ + Coulomb + Gaussians + repulsive core) — note `xav(npt)` is
  still zeroed afterwards and the gradient `xav1` is computed before that, as in the Fortran.
* `control-p bc` (`npcontrols = 2`, `dmaxfac = (1+rat)**0.5`) scales b and c: the box becomes non-cubic
  (a = 33.1213, b = c = 33.6477 Å after 4 iterations).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex107 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).
Inputs are generated from the base EPSR-test files by `scripts/` (see `scripts/make_inputs.*`).

## Verification
gcov of molfit.f: the bonded `refine_intra` line is executed 51600 times, the non-bonded one 37800 times;
`xav(j,ic)=coul_intra_fac*ecoul(j,ic)` 69678 times. Rust vs Fortran: **ALL IDENTICAL** in both passes.
