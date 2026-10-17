# EPSR-test-ex24 — `epsrwts` (interactive neutron weights)

Program: `epsrwts` (`fortran_src/epsrwts.f`); all answers come from standard input
(`case*.stdin`), the `.ato` file is `EPSR-test/200DMSO400Water.ato`.

| stdin | what it tests |
|---|---|
| `case1` | 1 sample, per atom, D on H and exchangeable H1, not normalised → `D1tot.wts` |
| `case2` | per molecule; invalid sample count (5) re-asked; two-isotope C and H (the Fortran normalises only the last abundance — reproduced); invalid exchange flag (7) re-asked; unknown mass number (99) and a non-numeric answer re-asked; 2 samples with H and H1 substituted; base name with leading blanks and extra words; normalised/not normalised totals and difference → `HD2tot1/2.wts`, `HD2dif.wts` |
| `case3` | 3 samples (second-order difference), per molecule, `fracmix 0.5` → `X3tot1/2/3.wts`, `X3xx/xy/yy.wts` |
| `case4` | 3 samples, per atom, non-exchanging H and H1 with the same two isotopes (isotope-resolved intramolecular H-H weights), S and H substituted, `fracmix 0.3` → `X4*.wts` |

The `.ato` file is rewritten after each run.

Verification: **ALL IDENTICAL** (reference in `fortran_reference.tar.gz`).
