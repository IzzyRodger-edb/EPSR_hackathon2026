# EPSR-test-ex25 — `epsrwtsx` (interactive X-ray weights)

Program: `epsrwtsx` (`fortran_src/epsrwtsx.f`); answers from standard input,
`f0_WaasKirf.dat` read from the current folder.

| stdin | what it tests |
|---|---|
| `case1` | per atom, neutral-atom form factors, IAFF → `W1X-ray.wts` |
| `case2` | `nrho 3` (invalid → per molecule, written as -1), valence/ionic form factors, an unknown symbol (`Xx`) re-asked, MAFF q/delta (one pair comma-separated), `normx 2`, `iafformaff 2`, blank line before the base name → `W2X-ray.wts` |
| `case3` | per molecule, `.ato` file named without extension → `W3X-ray.wts` |

Verification: **ALL IDENTICAL** (reference in `fortran_reference.tar.gz`).
