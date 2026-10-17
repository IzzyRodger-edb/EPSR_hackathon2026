# EPSR-test-ex117 — coordination-number constraints: three-body angles

Program: `epsr`. All other inputs are the EPSR-test files unchanged (`ntimes 5`, `niter 1`, `rotfreq 5`), made by copying the 16 base input files of EPSR-test and appending the `coordcon` line(s) below to the end of `200DMSO400Water.pcof` (after the pair potentials). A `coordcon` line is `coordcon atc nc cctemp ccdistort ccgaussweight ccgausswidth ats1 ds1 dsmax1 [ats2 ds2 dsmax2 ...] [angles wt distort gw gwidth ang1 ...] [anglesintra wt distort gw gwidth ats ang ...]`; it is read by `read_pcof` (`init_constraint`), its lengths are scaled by `normalise_coordinates`, and it is rewritten by `save_files_others` as a 256-character record (so pass 2 and 3 re-read the program's own format).

`.pcof` line:
```
coordcon O1 2 5.0 1.0 0.5 1.0 H1 1.8 2.5 angles 2.0 1.0 0.3 0.5 104.5 109.5
```

Water O1 with 2 H1 partners plus `angles 2.0 1.0 0.3 0.5 104.5 109.5`: after the direct bonds, every pair of neighbours (ja1, ja2) of an O1 atom gets a derived bond (`wtcc=2`, negative `tbondcc`) whose length follows from the law of cosines with the angle closest to the actual separation (two angles -> the `ntbangles>1` search), unless the molecules are already neighbours. These bonds use `en_gauss` (Gaussian term `tbgaussweight=0.3`). NB reproduced: the Fortran uses `adistcc(i,ia)` for both distances `tb1` and `tb2`; the derived bonds are appended (array reallocation of the 4 bond arrays). The `.pcof` record carries the `angles` part (Fortran FORMAT with an extra `)` that gfortran ignores). Constraint energy per molecule: -8.770 kJ/mol.

Common to all coordination-constraint cases (ex115–ex119): `find_constraints` is run by `fmolec1` (which also writes `coordcon.txt`) and by `press`; `en_mole_coord` is added to the molecule energy in `fmolec1`, `en_change` (MC moves) and `press`; atom pairs that are constraint-bonded (`check_coord_bonded`) use the reference potential `epot` without the extra Gaussians / repulsive core (`en_atom_inter_scale_coordextra`, and `epot/epot1` in `press`); `.out` reports `fcoordconstraints`.

Fortran quirks reproduced: `save_coord_energies(im)` is only called after accepted *internal-rotation* moves, so after an accepted whole-molecule move `ergold_coord(im)` keeps its old value (from `fmolec1`); only bonds of the moved molecule's own atoms enter `en_mole_coord(im)`; in pass 2 of `find_constraints` the search for the atom to remove stops at slot 2; the bond arrays keep their contents between calls (stale slot 1 is read in pass 3 for an atom whose bonds were all removed).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3`.

Verification: Rust vs Fortran reference (gfortran -O2 -fno-tree-vectorize, `fortran_reference.tar.gz`) — **ALL IDENTICAL**, all 3 passes (including `coordcon.txt`).

## Expected divergence from the Fortran (deliberate bug fix)

Since the bug-fix pass (see FORTRAN_BUGS.md at the repo root), this case
deliberately differs from the Fortran reference because of the coordination-constraint fix: the tb-angle distance now uses adistcc(j,ia) for the second arm instead of the same arm twice.
The difference is the fix working; run_ex_test.sh reporting SOME FILES DIFFER
here is expected.
