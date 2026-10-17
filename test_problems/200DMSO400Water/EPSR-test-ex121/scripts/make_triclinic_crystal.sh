#!/usr/bin/env bash
# Build the triclinic crystal of tethered S / O1 atoms used by EPSR-test-ex121/123.
# usage: make_triclinic_crystal.sh FORTRAN_BUILD_DIR WORK_DIR
# tric.cell: a,b,c = 4.3 4.8 5.4 A, alpha,beta,gamma = 80 100 110 deg, 2 atoms
# (S at 0,0,0 and O1 at 0.5,0.45,0.55), atom types S and O1 (labels that exist
# in the EPSR-test .pcof and weights files), no charges.
# makelattice tric.cell 7 7 6 -> tric.ato (588 tethered one-atom molecules).
set -euo pipefail
FB=$(realpath "$1"); WD=$2
mkdir -p "$WD"; cd "$WD"
ulimit -s unlimited || true
export OMP_NUM_THREADS=1
cat > tric.cell <<'CELL'
4.3 4.8 5.4
80.0 100.0 110.0
fractional
2
S 0.0 0.0 0.0
O1 0.5 0.45 0.55
S S 0
0.8 3.5 32.06 0.0 0.0
O1 O 0
0.65 3.0 15.999 0.0 0.0
0.5 2.0
CELL
echo tric.ato > ml.stdin
"$FB/makelattice" ./ makelattice tric.cell 7 7 6 < ml.stdin > makelattice.log 2>&1
echo "built $WD/tric.ato"
