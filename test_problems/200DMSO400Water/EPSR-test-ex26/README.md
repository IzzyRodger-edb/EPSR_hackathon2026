# EPSR-test-ex26 — `changeato` (editing .ato files), parameters and bonds

Program: `changeato` (`fortran_src/changeato.f` + the interactive keyword loop and
the changeato routines of `process_symbols.f`). Commands come from standard input.
`MixA.ato`/`MixC.ato` are copies of `EPSR-test/200DMSO400Water.ato`, `DMSOB.ato` of
`EPSR-test/DMSO.ato`.

| stdin | what it tests |
|---|---|
| `caseA` (MixA) | all global parameters (temp, tol, step sizes, vibtemp, angtemp, dihtemp, ecore, dcore), `rho 0.09` (box rescaled through `change_boxsize`, molecule centres moved), `box cubic`, an attempt to change the fixed `natomtypes`, a non-numeric value (rejected), `n` (revert and re-apply the previous value), per-type values of atom type 2 (epsilon, sigma, atmass, charge, qradius, atsymbol with an unknown symbol rejected, iexchange, natomsph, radminsph, radmaxsph, cgmolfac < 1 → 1), `atlabel` of type 6, `d` + bare value, `save` + overwrite `y`, `e` |
| `caseB` (DMSOB) | `bond` (change a bond length), `addbond`, `removebond`, `label` (relabel half of the H atoms at random — uses the random-number state of the `.ato` file), `rho`; exit with `e`, save at the prompt, refuse overwriting and give the new name `DMSOnew` |
| `caseC` (MixC) | `remove` all water molecules, change a charge, `forcesave`, `q` |

Outputs: `MixA.ato`, `DMSOnew.ato`, `MixC.ato`.

Verification: **ALL IDENTICAL** (reference in `fortran_reference.tar.gz`).
