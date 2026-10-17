# EPSR-test-ex110 — three-body potential: water O–H···O hydrogen bond

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.pcof`: one three-body line appended after the pair potentials
```
tb O1  H1  O1  20.0 10.0 0.2 1.2 2.6 1 0 0 180.0 1.0
```
i.e. triplets O1(donor)–H1(centre)···O1(acceptor) with O1–H1 on the same molecule (< 1.2 Å),
H1···O1 on different molecules (< 2.6 Å), requested angle 180° (amplitude 20).
`.EPSR.inp`: `niter 2`.  `commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2`.

Exercises: `read_pcof` three-body branch (symbol `tb`), `find_triplets` (called by `fmolec1` and `press`,
same-molecule and different-molecule flags, the `ltyp1 = ltyp3 → ia3 > ia1` rule, triplet-neighbour lists),
`en_mole_triplets` in `fmolec1`, in `en_change` (moved molecule + its triplet neighbours) and in `press`,
`save_threeb_energies` after accepted whole-molecule moves and internal rotations, `pr_mole_triplets`
(three-body pressure), `get_triplet_pot_cons`, normalisation of `threebrmax`, the `.out` lines
"Three-body energy/pressure" and "Energy for each triplet", and writing the `tb` line in `save_files_others`.

Fortran quirk reproduced: `read_pcof` stores the cosines of the angles and sets the angular width to 1.0,
but `save_files_others` (end of every iteration) re-reads the `tb` symbol and stores the angles and the width
in radians.  In the second iteration of each pass the triplet energy therefore jumps from ~0.2 to ~2000 kJ/mol
per molecule (cos θ is compared with π rad, width 0.1745), which visibly changes the Monte Carlo.

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex110 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).

## Expected divergence from the Fortran (deliberate bug fix)

Since the bug-fix pass (see FORTRAN_BUGS.md at the repo root), this case
deliberately differs from the Fortran reference because of the three-body fix: the in-memory triplet angles/width/amplitudes are no longer overwritten with radians at .pcof save time, so iterations 2+ use the correct cosine convention (Fortran iteration-2 energy +1953.6 kJ/mol vs fixed Rust -39.4).
The difference is the fix working; run_ex_test.sh reporting SOME FILES DIFFER
here is expected.
