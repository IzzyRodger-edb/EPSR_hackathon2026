#!/usr/bin/env bash
# Regenerate the inputs of EPSR-test-ex122 into OUT_DIR.
# usage: make_case.sh FORTRAN_BUILD_DIR OUT_DIR
# FORTRAN_BUILD_DIR must contain an epsr built with the corrected sharm_diffuse
# title FORMAT (tools/build_fortran_reference.sh does this); the original source
# stops with "Missing comma between descriptors" when `diffuse` is active.
# 1. DMSO/water crystal lattr.ato (see make_dmso_water_crystal.sh)
# 2. pre-run: EPSR-test inputs + lattr.ato, iinit 1, ntimes 2, Bragg + diffuse
#    keywords (below) -> accumulators (.inpa, .g01 ...), .ato, .erg and the
#    three <run>dNN.DIFFUSE.h01 files
# 3. case inputs = pre-run outputs with nsumt 1 (iinit is 0 after the pre-run),
#    keeping only the d01 and d02 .DIFFUSE.h01 files (d03 missing on purpose).
set -euo pipefail
FB=$(realpath "$1"); OUT=$2; HERE=$(cd "$(dirname "$0")" && pwd); B=$HERE/../../../../EPSR-test
P="$OUT/prerun"; mkdir -p "$P"
for f in 200DMSO400Water.EPSR.inp 200DMSO400Water.EPSR.inpa 200DMSO400Water.pcof 200DMSO400Water.EPSR.erg 200DMSO400Water.EPSR.g01 200DMSO400Water.EPSR.y01 200DMSO400Water.EPSR.f01 200DMSO400Water.EPSR.s01 D6DMSO_D2O.mint01 D6DMSO_HDO.mint01 D6DMSO_H2O.mint01 DMSODWaterD.NWTStot.wts DMSODWaterNull.NWTStot.wts DMSODWaterH.NWTStot.wts system_commands.txt; do cp "$B/$f" "$P/"; done
EPSR_TEST=$B "$HERE/make_dmso_water_crystal.sh" "$FB" "$OUT/crystal_build"
cp "$OUT/crystal_build/lattr.ato" "$P/200DMSO400Water.ato"
python3 "$HERE/setkw.py" "$P/200DMSO400Water.EPSR.inp" iinit=1 ntimes=2 "qwidthqmax=0.01 2.2" "mplicities=5 4 4 1 1 0" "hklqmin=0.0 10.0 0.1" "diffuse=2 2 1 4 2 1"
( cd "$P" && ulimit -s unlimited; OMP_NUM_THREADS=1 "$FB/epsr" ./ epsr 200DMSO400Water.EPSR.inp > prerun.log 2>&1 )
cd "$P"
cp 200DMSO400Water.EPSR.inp 200DMSO400Water.EPSR.inpa 200DMSO400Water.pcof 200DMSO400Water.EPSR.erg 200DMSO400Water.EPSR.g01 200DMSO400Water.EPSR.y01 200DMSO400Water.EPSR.f01 200DMSO400Water.EPSR.s01 200DMSO400Water.ato D6DMSO_D2O.mint01 D6DMSO_HDO.mint01 D6DMSO_H2O.mint01 DMSODWaterD.NWTStot.wts DMSODWaterNull.NWTStot.wts DMSODWaterH.NWTStot.wts system_commands.txt 200DMSO400Waterd01.DIFFUSE.h01 200DMSO400Waterd02.DIFFUSE.h01 "$OUT/"
python3 "$HERE/setkw.py" "$OUT/200DMSO400Water.EPSR.inp" nsumt=1
