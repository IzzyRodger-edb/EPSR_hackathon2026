# EPSR-test-ex21 — `nwts`, first-order differences and interactive corrections

Program: `nwts`. Inputs made from `EPSR-test/DMSODWaterD.NWTS.dat` by changing
keyword values; `WaterX.ato` is `EPSR-test/WaterMolecule.ato` with the chemical
symbol of O1 changed to the unknown `Ox`.

| command | what it tests |
|---|---|
| `FirstDiffA` | 2 samples, per molecule (`output 2`), `normtot 1`, `normdif 1`; H with two isotopes, H and H1 substituted, `iexchange 1` requested for an O atom (forced to 0 by `check_symbol`) |
| `FirstDiffB` | 2 samples, `output 5` (clamped to 2), `normtot 2` (<b^2> normalisation), `normdif 0`; substitution of S (34S) and H (two isotopes in sample 2), O1 with a zero-abundance isotope, no exchange |
| `FirstDiffC < diffc.stdin` | unknown mass numbers 99 and 77 in `abundances`/`abund2`: the program asks for new mass numbers on standard input |
| `WaterX < waterx.stdin` | input file with only `fnameato`: all component entries are filled from the `.ato` file (defaults); the unknown chemical symbol `Ox` is re-typed on standard input and the corrected symbol written back to `WaterX.ato` |

Outputs: `<run>tot1.wts`, `<run>tot2.wts`, `<run>dif.wts`, `WaterX.NWTStot.wts`,
rewritten `.NWTS.dat` files and `WaterX.ato`.

Verification: **ALL IDENTICAL** (reference in `fortran_reference.tar.gz`).
