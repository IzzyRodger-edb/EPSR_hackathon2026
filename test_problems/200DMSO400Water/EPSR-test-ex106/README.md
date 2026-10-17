# EPSR-test-ex106 — pressure control through the repulsive cores (`control-p p`)

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets).

`.EPSR.inp`: `control-p p 0.0 0.01`, `ntimes 2`, `niter 2`. `.pcof`: repulsive cores for two pairs
(`rminex pcoefex coordex`): O-O1 `2.6 0.5 0` and O1-O1 `2.7 0 0` (all other pairs 0, as in the base case).
`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2`.

With `p`, after each iteration `pfac = 1/(1+0.01*clip(pressd))` and (`potfac > 0`) `multiply_pcoefex(pfac)`
multiplies `rminex` of every pair with `pcoefex /= 0`, `rminex > 0` (and `|pcoefex| > 1e-4` or `pfac > 1`),
rewrites the pair symbol value with a list-directed write (`rminex` followed by the remaining words of the old
value) and regenerates the reference potential (`generate_rp`, exponential repulsive core). This also
exercises the `rminex > 0` branches of `repulsivept` (the O1-O1 core starts with `pcoefex = 0` and is grown
by `repulsivept`; O-O1 rminex 2.6 -> 2.4985, O1-O1 2.7 -> 2.5946 after 4 iterations).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex106 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).
Inputs are generated from the base EPSR-test files by `scripts/` (see `scripts/make_inputs.*`).

## Verification
gcov: `multiply_pcoefex` called in every iteration, `rminex(irefp)=rminex(irefp)*epfact` executed 8 times
(the one-word branch of the symbol rewrite is not reached). Rust vs Fortran: **ALL IDENTICAL** in both passes.
