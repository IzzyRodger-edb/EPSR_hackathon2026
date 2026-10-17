# EPSR-test-ex118 — coordination-number constraints: angles + intramolecular angles

Program: `epsr`. All other inputs are the EPSR-test files unchanged (`ntimes 5`, `niter 1`, `rotfreq 5`), made by copying the 16 base input files of EPSR-test and appending the `coordcon` line(s) below to the end of `200DMSO400Water.pcof` (after the pair potentials). A `coordcon` line is `coordcon atc nc cctemp ccdistort ccgaussweight ccgausswidth ats1 ds1 dsmax1 [ats2 ds2 dsmax2 ...] [angles wt distort gw gwidth ang1 ...] [anglesintra wt distort gw gwidth ats ang ...]`; it is read by `read_pcof` (`init_constraint`), its lengths are scaled by `normalise_coordinates`, and it is rewritten by `save_files_others` as a 256-character record (so pass 2 and 3 re-read the program's own format).

`.pcof` line:
```
coordcon O1 2 5.0 1.0 0.5 1.0 H1 1.8 2.5 angles 1.5 0.9 0.0 0.0 109.5 anglesintra 1.0 1.0 0.2 0.5 H1 104.5 H1 120.0
```

Water O1 with 2 H1 partners, `angles` (one angle, `tbgaussweight=0`) and `anglesintra 1.0 1.0 0.2 0.5 H1 104.5 H1 120.0`: for each of the two intra angles, every own H1 of an O1 molecule (donor limit = number of H1 per molecule) is bonded (`wtcc=3`) to the near neighbour of O1 whose separation best matches the angle, subject to the donor/acceptor bookkeeping (no mutual donor-acceptor pairs, both passes of the donor search; 600/596 acceptors/donors). Exercised: the full `init_constraint` parser (`angles` followed by `anglesintra`), `hastbintra` part of `find_constraints`, `wtcc=3` parameters in `en_coordcon`, and the `.pcof` record with both angle parts. Constraint energy per molecule: -11.64 kJ/mol.

Common to all coordination-constraint cases (ex115–ex119): `find_constraints` is run by `fmolec1` (which also writes `coordcon.txt`) and by `press`; `en_mole_coord` is added to the molecule energy in `fmolec1`, `en_change` (MC moves) and `press`; atom pairs that are constraint-bonded (`check_coord_bonded`) use the reference potential `epot` without the extra Gaussians / repulsive core (`en_atom_inter_scale_coordextra`, and `epot/epot1` in `press`); `.out` reports `fcoordconstraints`.

Fortran quirks reproduced: `save_coord_energies(im)` is only called after accepted *internal-rotation* moves, so after an accepted whole-molecule move `ergold_coord(im)` keeps its old value (from `fmolec1`); only bonds of the moved molecule's own atoms enter `en_mole_coord(im)`; in pass 2 of `find_constraints` the search for the atom to remove stops at slot 2; the bond arrays keep their contents between calls (stale slot 1 is read in pass 3 for an atom whose bonds were all removed).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3`.

Verification: Rust vs Fortran reference (gfortran -O2 -fno-tree-vectorize, `fortran_reference.tar.gz`) — **ALL IDENTICAL**, all 3 passes (including `coordcon.txt`).

## Expected divergence from the Fortran (deliberate bug fix)

Since the bug-fix pass (see FORTRAN_BUGS.md at the repo root), this case
deliberately differs from the Fortran reference because of the coordination-constraint fix: the tb-angle distance now uses adistcc(j,ia) for the second arm instead of the same arm twice.
The difference is the fix working; run_ex_test.sh reporting SOME FILES DIFFER
here is expected.
