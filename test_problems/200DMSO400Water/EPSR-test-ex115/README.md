# EPSR-test-ex115 — coordination-number constraints: water O1 accepting 2 H1

Program: `epsr`. All other inputs are the EPSR-test files unchanged (`ntimes 5`, `niter 1`, `rotfreq 5`), made by copying the 16 base input files of EPSR-test and appending the `coordcon` line(s) below to the end of `200DMSO400Water.pcof` (after the pair potentials). A `coordcon` line is `coordcon atc nc cctemp ccdistort ccgaussweight ccgausswidth ats1 ds1 dsmax1 [ats2 ds2 dsmax2 ...] [angles wt distort gw gwidth ang1 ...] [anglesintra wt distort gw gwidth ats ang ...]`; it is read by `read_pcof` (`init_constraint`), its lengths are scaled by `normalise_coordinates`, and it is rewritten by `save_files_others` as a 256-character record (so pass 2 and 3 re-read the program's own format).

`.pcof` line:
```
coordcon O1 2 5.0 1.0 0.5 1.0 H1 1.8 2.5
```

Each water oxygen (O1, 400 atoms) is constrained to 2 hydrogen-bond partners H1 (on other molecules) at 1.8 Å (max. 2.5 Å); `alternatenc = nint(2*cfrac(O1)/cfrac(H1)) = 1`, so each H1 takes one O1. Morse-like energy `en_harmonic` with depth `cctemp=5`, `ccdistort=1`, plus the harmonic/linear term (`ccgaussweight=0.5`, `ccgausswidth=1`, both branches |args|<1 and >=1 taken). Exercised: passes 1-3 of `find_constraints` (nearest-neighbour lists on distinct molecules, ordered insertion of the reverse bonds, removal of over-bonded atoms, completion of under-bonded O1 atoms: 97.7 % of the atoms end with the expected number of bonds, 66 % of the bonds are beyond `dsmax`). Coordination-constraint energy per molecule (`.out`, pass 3): -7.546 kJ/mol.

Common to all coordination-constraint cases (ex115–ex119): `find_constraints` is run by `fmolec1` (which also writes `coordcon.txt`) and by `press`; `en_mole_coord` is added to the molecule energy in `fmolec1`, `en_change` (MC moves) and `press`; atom pairs that are constraint-bonded (`check_coord_bonded`) use the reference potential `epot` without the extra Gaussians / repulsive core (`en_atom_inter_scale_coordextra`, and `epot/epot1` in `press`); `.out` reports `fcoordconstraints`.

Fortran quirks reproduced: `save_coord_energies(im)` is only called after accepted *internal-rotation* moves, so after an accepted whole-molecule move `ergold_coord(im)` keeps its old value (from `fmolec1`); only bonds of the moved molecule's own atoms enter `en_mole_coord(im)`; in pass 2 of `find_constraints` the search for the atom to remove stops at slot 2; the bond arrays keep their contents between calls (stale slot 1 is read in pass 3 for an atom whose bonds were all removed).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3`.

Verification: Rust vs Fortran reference (gfortran -O2 -fno-tree-vectorize, `fortran_reference.tar.gz`) — **ALL IDENTICAL**, all 3 passes (including `coordcon.txt`).

The `coordex`/`pcoefex` coordination-number control of the exponential repulsive term (`repulsivept`) is already exercised by EPSR-test-ex1 (`coordex` 0.1 and 0.05 in its `.pcof`), so it is not repeated here.
