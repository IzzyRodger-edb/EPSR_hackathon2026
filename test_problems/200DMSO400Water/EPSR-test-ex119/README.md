# EPSR-test-ex119 — coordination-number constraints: over-bonded secondary atoms

Program: `epsr`. All other inputs are the EPSR-test files unchanged (`ntimes 5`, `niter 1`, `rotfreq 5`), made by copying the 16 base input files of EPSR-test and appending the `coordcon` line(s) below to the end of `200DMSO400Water.pcof` (after the pair potentials). A `coordcon` line is `coordcon atc nc cctemp ccdistort ccgaussweight ccgausswidth ats1 ds1 dsmax1 [ats2 ds2 dsmax2 ...] [angles wt distort gw gwidth ang1 ...] [anglesintra wt distort gw gwidth ats ang ...]`; it is read by `read_pcof` (`init_constraint`), its lengths are scaled by `normalise_coordinates`, and it is rewritten by `save_files_others` as a 256-character record (so pass 2 and 3 re-read the program's own format).

`.pcof` line:
```
coordcon H1 1 3.0 1.0 0.2 0.5 O 1.9 2.8
```

Every water H1 (800 atoms) must bond to `nc=1` DMSO O (200 atoms) at 1.9/2.8 Å, so `alternatenc = nint(800/200) = 4`. Many H1 pick the same nearest O, so the DMSO O lists overflow `nbondccsize` (reallocation in pass 1), pass 2 removes the farthest H1 from over-bonded O (search from the end, shift of the remaining entries), and pass 3 re-bonds the H1 left with no bond to the nearest O that still has room (reading the stale slot 1 of the emptied H1 list, as the Fortran does). Finally 100 % of the atoms have the expected number of bonds. `coordcon.txt` has 4-bond records. The constraint involves both molecule types (DMSO rotations -> `save_coord_energies`). Constraint energy per molecule: -1.954 kJ/mol.

Common to all coordination-constraint cases (ex115–ex119): `find_constraints` is run by `fmolec1` (which also writes `coordcon.txt`) and by `press`; `en_mole_coord` is added to the molecule energy in `fmolec1`, `en_change` (MC moves) and `press`; atom pairs that are constraint-bonded (`check_coord_bonded`) use the reference potential `epot` without the extra Gaussians / repulsive core (`en_atom_inter_scale_coordextra`, and `epot/epot1` in `press`); `.out` reports `fcoordconstraints`.

Fortran quirks reproduced: `save_coord_energies(im)` is only called after accepted *internal-rotation* moves, so after an accepted whole-molecule move `ergold_coord(im)` keeps its old value (from `fmolec1`); only bonds of the moved molecule's own atoms enter `en_mole_coord(im)`; in pass 2 of `find_constraints` the search for the atom to remove stops at slot 2; the bond arrays keep their contents between calls (stale slot 1 is read in pass 3 for an atom whose bonds were all removed).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3`.

Verification: Rust vs Fortran reference (gfortran -O2 -fno-tree-vectorize, `fortran_reference.tar.gz`) — **ALL IDENTICAL**, all 3 passes (including `coordcon.txt`).
