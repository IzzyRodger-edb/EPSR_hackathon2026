# EPSR-test-ex23 — `xwts` (X-ray form factors)

Program: `xwts` (`fortran_src/xwts.f`, keyword file `.XWTS.dat`). The form-factor
table `f0_WaasKirf.dat` (copied from `EPSR-test`) is read from the current folder.

| command | what it tests |
|---|---|
| `Xray1` | input with only `fnameato`: the `fnameato` branch of `check_symbol` fills the components (labels and default form-factor symbols from the `.ato` file) |
| `Xray2` | per molecule, `normtot 1`, `iafformaff 2`, ionic/valence form factors (`Cval`, `O1-`, `O2-`, `Hh+`) and non-zero MAFF q/delta values |
| `Xray3 < xwts3.stdin` | `output 0` and `normtot 5` (clamped by `check_symbol`); unknown symbols `Zz` and `Xx` re-typed on standard input (one wrong answer `Qq` first) |

Outputs: `<run>.wts` (form-factor coefficients, format `11(1x,f11.6)` with
reversion) and the rewritten `.XWTS.dat` files.

Verification: **ALL IDENTICAL** (reference in `fortran_reference.tar.gz`).
