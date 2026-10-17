#!/usr/bin/env bash
# Run the Fortran reference and the Rust port side by side for N chained EPSR
# runs (each run continues from the files written by the previous one) and
# compare all output files after every run.
#
# usage: chain_test.sh FORTRAN_EPSR_BINARY PRISTINE_TEST_DIR WORK_DIR [N]
#
#   FORTRAN_EPSR_BINARY  epsr executable built by tools/build_fortran_reference.sh
#   PRISTINE_TEST_DIR    folder with the test inputs (copied, never modified)
#   WORK_DIR             scratch folder (created/overwritten)
#   N                    number of chained runs (default 10)
set -euo pipefail
FBIN=$(realpath "$1"); SRC=$(realpath "$2"); WORK=$3; N=${4:-10}
HERE=$(cd "$(dirname "$0")" && pwd)
RBIN="$HERE/../target/release/epsr-rs"
INP=$(cd "$SRC" && ls *.EPSR.inp | head -1)
rm -rf "$WORK"; mkdir -p "$WORK"
cp -r "$SRC" "$WORK/fortran"; cp -r "$SRC" "$WORK/rust"
ulimit -s unlimited || true
status=0
for k in $(seq 1 "$N"); do
  ( cd "$WORK/fortran" && OMP_NUM_THREADS=1 "$FBIN" ./ epsr "$INP" < /dev/null > run.log 2>&1 )
  ( cd "$WORK/rust" && "$RBIN" ./ epsr "$INP" --quiet > run.log 2>&1 )
  if python3 "$HERE/compare_outputs.py" "$WORK/fortran" "$WORK/rust" --quiet > "$WORK/compare_$k.txt"; then
    echo "run $k: ALL IDENTICAL   $(grep -m1 'Chi-square' "$WORK/rust/"*.EPSR.out | awk '{print "chi2="$2}')"
  else
    echo "run $k: DIFFERENCES (see $WORK/compare_$k.txt)"; cat "$WORK/compare_$k.txt"; status=1; break
  fi
done
exit $status
