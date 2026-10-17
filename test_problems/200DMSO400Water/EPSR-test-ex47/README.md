# EPSR-test-ex47 — alignmole, single water molecules (all moves accepted)

Program: `alignmole` (`fortran_src/alignmole.f`, standalone program; command name `alignmole`): Monte Carlo
rotation (`rotatm`, `ran1`) of every molecule to minimise the mean square distance of its atoms from an axis.

Inputs (built from `WaterMolecule.ato` with a one-off helper script that was not kept, step size `stepr` = 0.05, random
generator re-seeded with `idum = -seed`):

* `waterz.ato` — molecule in the xy plane (seed 61): `alignmole waterz.ato z 10 all`;
* `waterx.ato` — molecule in the yz plane (seed 329): `alignmole waterx x 10 1 2 3` (axis-atom list branch).

Starting at the maximum of the displacement, and with these seeds (found by search), all 10 trial rotations
decrease the displacement; the step size is adapted at iteration 10 and restored before writing.

Expected outputs: `waterz.ato`, `waterx.ato` (rotated coordinates, advanced random-number state).

**Limitation of the Fortran.** When a trial rotation is rejected, `alignmole.f` restores the molecule from the work arrays `xr, yr, zr`, which are allocated but never assigned: it copies uninitialised heap memory into the coordinates (in practice zeros and leftovers of earlier allocations, `**************` in the `.ato`). Such runs cannot be reproduced; the Rust translation uses zero-filled work arrays and prints a warning. The cases here are therefore designed so that no rotation is ever rejected. Also `niter < 10` makes the Fortran divide by zero (`naccum = niter/10`); the Rust port reports an input error.

Verification: `ALL IDENTICAL` against the original Fortran reference binary.
