# EPSR-test-ex22 — `nwts`, second-order differences

Program: `nwts`, three samples (`nsamples 3`): sample 3 is a mixture of samples 1
and 2 (`xfrac`), and the XX, XY (=XH) and YY (=HH) partial weights are written.
Inputs made from `EPSR-test/DMSODWaterD.NWTS.dat`.

| input | options |
|---|---|
| `SecondDiffA` | per molecule, `xfrac 0.4`, `normtot 0`, H and H1 substituted (H1 exchangeable) |
| `SecondDiffB` | per atom, `xfrac 0.65`, `normtot 2`, C (13C), H and H1 substituted, no exchange (intramolecular YY/XY/XX branches for non-exchanging substituted atoms) |
| `SecondDiffC` | `nsamples 4` (clamped to 3 by `check_symbol`), `normtot 1` |

Outputs: `<run>tot1/2/3.wts`, `<run>xx.wts`, `<run>xy.wts`, `<run>yy.wts`, the
rewritten `.NWTS.dat` files.

Verification: **ALL IDENTICAL** (reference in `fortran_reference.tar.gz`).
