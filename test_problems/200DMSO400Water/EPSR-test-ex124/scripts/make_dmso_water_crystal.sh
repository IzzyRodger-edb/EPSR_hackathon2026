#!/usr/bin/env bash
# Build the DMSO/water molecular crystal used by EPSR-test-ex120/122/124.
# usage: make_dmso_water_crystal.sh FORTRAN_BUILD_DIR WORK_DIR
#   FORTRAN_BUILD_DIR: reference programs (mixato, randomise, epsr, makelatticeato)
#   EPSR_TEST: folder with the base EPSR-test files (default ../../../../EPSR-test, i.e. EPSR-test in the repository root)
# Steps:
#  1. mixato: 2 DMSO + 4 water at 0.1 atoms/A^3 -> cell.ato (cubic, 6.84 A; all
#     molecules at the same centre)
#  2. randomise cell -> random positions/orientations
#  3. epsr (EPSR-test inputs, potfac 0, iinit 1, ntimes 200) on the single cell
#     to remove overlaps -> cellr.ato
#  4. makelatticeato cellr 5 4 4 lattr n  -> orthorhombic 34.2 x 27.4 x 27.4 A box
#     (480 molecules, 2560 atoms)  [ex120, ex122]
#     makelatticeato cellr 4 4 4 lattc n  -> cubic 27.4 A box (384 molecules)  [ex124]
set -euo pipefail
FB=$(realpath "$1"); WD=$2
HERE=$(cd "$(dirname "$0")" && pwd)
B=${EPSR_TEST:-$HERE/../../../../EPSR-test}
mkdir -p "$WD"; cd "$WD"
ulimit -s unlimited || true
export OMP_NUM_THREADS=1
cp "$B/DMSO.ato" "$B/WaterMolecule.ato" "$B/system_commands.txt" .
printf '2\ne\nDMSO\n2\ne\nWaterMolecule\n4\n0.1\ncell\n' > mixato.stdin
"$FB/mixato" ./ mixato < mixato.stdin > mixato.log 2>&1
"$FB/randomise" ./ randomise cell > randomise.log 2>&1
for f in 200DMSO400Water.EPSR.inp 200DMSO400Water.EPSR.inpa 200DMSO400Water.pcof 200DMSO400Water.EPSR.erg \
  D6DMSO_D2O.mint01 D6DMSO_HDO.mint01 D6DMSO_H2O.mint01 DMSODWaterD.NWTStot.wts DMSODWaterNull.NWTStot.wts DMSODWaterH.NWTStot.wts; do
  cp "$B/$f" .
done
cp cell.ato 200DMSO400Water.ato
python3 "$HERE/setkw.py" 200DMSO400Water.EPSR.inp iinit=1 ntimes=200 potfac="0.0 0.0"
"$FB/epsr" ./ epsr 200DMSO400Water.EPSR.inp > relax.log 2>&1
cp 200DMSO400Water.ato cellr.ato
"$FB/makelatticeato" ./ makelatticeato cellr 5 4 4 lattr n > mla1.log 2>&1
"$FB/makelatticeato" ./ makelatticeato cellr 4 4 4 lattc n > mla2.log 2>&1
echo "built $WD/lattr.ato and $WD/lattc.ato"
