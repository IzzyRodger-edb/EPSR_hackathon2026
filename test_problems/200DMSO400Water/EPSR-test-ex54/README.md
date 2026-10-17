# EPSR-test-ex54 — `torangles` (torsion-angle distributions)

Program: `torangles` (`fortran_src/torangles.f`), input `dmso.TOR.dat`, repeat 2 (nsumt 0 ->
accumulation):

* DMSO (molecules starting with S): O4-S1-C2-H5 torsion (`ltyp3 4 4`, `ltyp4 5 5`) and
  C2-S1-C3-<mid(H8,H9)> (`ltyp4 8 9`, mid-point of two atoms).
* Water (O1): degenerate O-H-H "torsion" (`ltyp3 3 3`, `ltyp4 3 3`).

Outputs compared: `dmso.TOR.c01`, `dmso.TOR.dat`, `toranglesvalues.txt` (angle of each
molecule, run-time format `(1x,i6,n(1x,e14.7))`).

Verification: ALL IDENTICAL (passes 1 and 2).
