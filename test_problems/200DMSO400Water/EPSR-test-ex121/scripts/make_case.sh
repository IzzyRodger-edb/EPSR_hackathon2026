#!/usr/bin/env bash
# Regenerate the inputs of EPSR-test-ex121 into OUT_DIR.
# usage: make_case.sh FORTRAN_BUILD_DIR OUT_DIR
set -euo pipefail
FB=$(realpath "$1"); OUT=$2; HERE=$(cd "$(dirname "$0")" && pwd); B=$HERE/../../../../EPSR-test
mkdir -p "$OUT"
for f in 200DMSO400Water.EPSR.inp 200DMSO400Water.EPSR.inpa 200DMSO400Water.pcof 200DMSO400Water.EPSR.erg 200DMSO400Water.EPSR.g01 200DMSO400Water.EPSR.y01 200DMSO400Water.EPSR.f01 200DMSO400Water.EPSR.s01 D6DMSO_D2O.mint01 D6DMSO_HDO.mint01 D6DMSO_H2O.mint01 DMSODWaterD.NWTStot.wts DMSODWaterNull.NWTStot.wts DMSODWaterH.NWTStot.wts system_commands.txt; do cp "$B/$f" "$OUT/"; done
EPSR_TEST=$B "$HERE/make_triclinic_crystal.sh" "$FB" "$OUT/crystal_build"
cp "$OUT/crystal_build/tric.ato" "$OUT/200DMSO400Water.ato"
python3 "$HERE/setkw.py" "$OUT/200DMSO400Water.EPSR.inp" iinit=1 ntimes=2 "qwidthqmax=-0.02 4.0" "mplicities=7 7 6 1 1 2" "hklqmin=1.0 8.0 0.0"
