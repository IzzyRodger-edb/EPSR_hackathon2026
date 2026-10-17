# EPSR-test-ex116 — coordination-number constraints: two constraints, two secondary types

Program: `epsr`. All other inputs are the EPSR-test files unchanged (`ntimes 5`, `niter 1`, `rotfreq 5`), made by copying the 16 base input files of EPSR-test and appending the `coordcon` line(s) below to the end of `200DMSO400Water.pcof` (after the pair potentials). A `coordcon` line is `coordcon atc nc cctemp ccdistort ccgaussweight ccgausswidth ats1 ds1 dsmax1 [ats2 ds2 dsmax2 ...] [angles wt distort gw gwidth ang1 ...] [anglesintra wt distort gw gwidth ats ang ...]`; it is read by `read_pcof` (`init_constraint`), its lengths are scaled by `normalise_coordinates`, and it is rewritten by `save_files_others` as a 256-character record (so pass 2 and 3 re-read the program's own format).

`.pcof` lines:
```
coordcon O1 3 4.0 0.8 0.3 0.5 H1 1.85 2.6 O1 2.8 3.4
coordcon O 2 6.0 1.2 0.0 0.0 H1 1.9 2.7
```

Two constraints: water O1 with `nc=3` neighbours of two types (H1 at 1.85/2.6 Å, O1 at 2.8/3.4 Å; O1-O1 bonds are same-type, so no reverse bond is stored), `alternatenc = nint(3*400/1200) = 1`; and DMSO O with `nc=2` H1 neighbours (`alternatenc = nint(2*200/800) = nint(0.5) = 1`) with `ccgaussweight=0` (no extra term). DMSO molecules have rotational groups, so accepted internal rotations call `save_coord_energies` (passes 2 and 3 diverge if it is omitted). Exercised in addition to ex115: `ncoordconstraints=2` in `check_coord_bonded`/`en_mole_coord`/`writecoordcons`, `ns=2` (type search loops, the second `ats` entry in the `.pcof` record), the reverse-bond search through the second constraint. Constraint energy per molecule: -11.80 kJ/mol.

Common to all coordination-constraint cases (ex115–ex119): `find_constraints` is run by `fmolec1` (which also writes `coordcon.txt`) and by `press`; `en_mole_coord` is added to the molecule energy in `fmolec1`, `en_change` (MC moves) and `press`; atom pairs that are constraint-bonded (`check_coord_bonded`) use the reference potential `epot` without the extra Gaussians / repulsive core (`en_atom_inter_scale_coordextra`, and `epot/epot1` in `press`); `.out` reports `fcoordconstraints`.

Fortran quirks reproduced: `save_coord_energies(im)` is only called after accepted *internal-rotation* moves, so after an accepted whole-molecule move `ergold_coord(im)` keeps its old value (from `fmolec1`); only bonds of the moved molecule's own atoms enter `en_mole_coord(im)`; in pass 2 of `find_constraints` the search for the atom to remove stops at slot 2; the bond arrays keep their contents between calls (stale slot 1 is read in pass 3 for an atom whose bonds were all removed).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3`.

Verification: Rust vs Fortran reference (gfortran -O2 -fno-tree-vectorize, `fortran_reference.tar.gz`) — **ALL IDENTICAL**, all 3 passes (including `coordcon.txt`).
