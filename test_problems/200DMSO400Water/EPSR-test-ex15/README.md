# EPSR-test-ex15 — container forces (`rmaxex`/`pcoefexmax`)

Program: `epsr`. Based on EPSR-test; inputs made by `scripts/make_ex15.py ../../../EPSR-test . 4.0 0.5` (run inside this folder):

* every DMSO molecule gets an unbonded, massless `Q` site at its S position (`qtype=1`), plus the `Q` atom type;
* the `.pcof` file gets the seven `Q-*` pairs (28 pairs); the `Q-O1` line carries 5 numbers,
  `rminex pcoefex coordex rmaxex=4.0 pcoefexmax=0.5`, which switches on a container force
  (`containerforces_mod.f90`): each water O1 further than 4 Å from its nearest Q site is pulled back with an
  exponential energy `pcoefexmax*exp((r-rmaxex)/(expecf*sqrt(rmaxex)))`;
* `iinit 1` (the accumulators of EPSR-test have 21 pairs).

Exercised: `read_pcof` 5-number lines, `get_containerforce_list` (from `generate_rp`), `en_containerforce`
in `fmolec1` and `en_change`, `dist_containerforce` in `atogen`, `repulsivept` (`pcoefexmax := potexmax`, the number of
contained atoms beyond `rmaxex` — 98 after pass 1, so pass 2 runs with a much stronger force), and writing
`rmaxex`/`pcoefexmax` back to the `.pcof` file and symbol table.

Fortran quirk reproduced: `ergold_container` is never updated after an accepted move (only by `fmolec1`).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2`.

Verification: Rust vs Fortran reference (gfortran -O2 -fno-tree-vectorize) — **ALL IDENTICAL**, both passes.
