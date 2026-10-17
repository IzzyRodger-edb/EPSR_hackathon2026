# Energy logs used by `make_readme_figures.py`

Column 1 of each file is EPSR's energy per molecule (kJ/mol) after each iteration (`.erg` files).

| file | run |
|---|---|
| `water_fortran.erg`, `water_rust.erg` | the 300-iteration water validation (`validation/water`, `run_validation.sh`), Fortran and Rust; byte-identical |
| `chaos_fortran.erg`, `chaos_rust.erg` | 30 further iterations from the final validated state (`niter 30`), Fortran and Rust; every output file identical |
| `chaos_nudged.erg` | the same 30 iterations in Rust after moving molecule 1 by +1×10⁻⁶ Å in x in `waterbox.ato` (`-0.11350319E+01` → `-0.11350309E+01`) |

In the nudged run the energy matches for 12 iterations, then differs from iteration 13 on. A
nudge of 1×10⁻⁷ Å (one single-precision step) did not change any Monte Carlo decision within
30 iterations; nudges of 10⁻⁵ to 10⁻² Å diverged within 1–2 iterations.
