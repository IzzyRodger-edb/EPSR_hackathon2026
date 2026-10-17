# EPSR-test-ex48 — alignmole on a box (zero displacement), stdin input

Program: `alignmole` (`fortran_src/alignmole.f`).

* `box.ato` (built with a one-off helper script that was not kept; molecule order DWDWWD): 6 DMSO/water molecules whose first atom sits exactly
  at the molecule centre. With axis atoms `1 1` the displacement is identically 0, so every rotation is
  accepted while the molecules still rotate (ran1 sequence, `rotatm`, step-size adaptation every niter/10
  iterations, several molecules):
  * `alignmole box.ato y 30 1 1`;
  * `alignmole box.ato x 20 7 < box_atoms.stdin` — a single axis atom on the command line is not enough, the
    list is read again from stdin (`1 1`).
* `atoms.ato` — three one-atom molecules (`rotatm` does nothing); axis and number of iterations from stdin
  (`atoms.stdin`), all atoms used.

Expected outputs: `box.ato` (after both runs), `atoms.ato`.

**Limitation of the Fortran.** When a trial rotation is rejected, `alignmole.f` restores the molecule from the work arrays `xr, yr, zr`, which are allocated but never assigned: it copies uninitialised heap memory into the coordinates (in practice zeros and leftovers of earlier allocations, `**************` in the `.ato`). Such runs cannot be reproduced; the Rust translation uses zero-filled work arrays and prints a warning. The cases here are therefore designed so that no rotation is ever rejected. Also `niter < 10` makes the Fortran divide by zero (`naccum = niter/10`); the Rust port reports an input error.

Verification: `ALL IDENTICAL` against the original Fortran reference binary.
