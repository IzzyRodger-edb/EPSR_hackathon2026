# EPSR-test-ex58 — `local_density` (multi-bond cluster analysis)

Program: `local_density` (`fortran_src/local_density.f`, keyword table of `clusters`,
input `hbnet.CLUSTERS.dat`), repeat 2 (nsumt 0 -> accumulation):

* cluster 1: single bond definition O1-H1 1.2-2.4 A (water H-bond network);
* cluster 2: two definitions (O1-H1, O-H1) — water and DMSO-water H-bonds;
* cluster 3: three definitions (S-S 0-5, C-O 2-3.2, O1-H1 1.5-2) — mixed clusters.

Checks parsing of word lists for atom1/atom2/rmin/rmax, the composed column titles
("O1-H1, O-H1"), cluster labelling, percolation fraction (row 1) and the
fraction-in-clusters lines in `fort.14`.

Outputs compared: `hbnet.CLUSTERS.n01`, `hbnet.CLUSTERS.dat`, `fort.14`.

Verification: ALL IDENTICAL (passes 1 and 2).
