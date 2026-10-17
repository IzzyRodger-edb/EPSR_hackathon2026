#!/usr/bin/env bash
# Run one EPSR-test-ex* verification case with the Fortran reference programs
# and with the Rust port, and compare every file they write.
#
# usage: run_ex_test.sh EX_DIR FORTRAN_BUILD_DIR [WORK_DIR] [--save-reference]
#
# EX_DIR contains the input files plus `commands.txt`.  Each non-comment line of
# commands.txt is one program invocation, run in the case folder:
#
#     <command> [args...] [< stdin_file]
#
# e.g.  `nwts DMSODWaterD.NWTS.dat`   or   `makeato < makeato.stdin`.
# A line `repeat N` runs the whole command list N times (chained runs).
# The Fortran program FORTRAN_BUILD_DIR/<command> (built by
# build_fortran_reference.sh) and `epsr-rs <folder> <command> ...` are each run
# on their own copy of EX_DIR.  With --save-reference the files created or
# changed by the Fortran run are stored in EX_DIR/fortran_reference.tar.gz.
set -euo pipefail
EX=$(realpath "$1"); FB=$(realpath "$2")
WORK=${3:-$(mktemp -d)}
SAVE=0; for a in "$@"; do [ "$a" = "--save-reference" ] && SAVE=1; done
[ "${3:-}" = "--save-reference" ] && WORK=$(mktemp -d)
HERE=$(cd "$(dirname "$0")" && pwd)
RBIN="$HERE/../target/release/epsr-rs"
rm -rf "$WORK"; mkdir -p "$WORK/fortran" "$WORK/rust" "$WORK/pristine" "$WORK/logs"
for d in fortran rust pristine; do
  (cd "$EX" && find . -maxdepth 1 -type f ! -name commands.txt ! -name README.md ! -name 'fortran_reference.tar.gz' ! -name 'EPSR-test-ex*.pdf' -exec cp {} "$WORK/$d/" \;)
done
REPEAT=$(grep -E '^repeat ' "$EX/commands.txt" | awk '{print $2}' | head -1 || true); REPEAT=${REPEAT:-1}
ulimit -s unlimited || true
n=0
for k in $(seq 1 "$REPEAT"); do
  while IFS= read -r line; do
    line=${line%%#*}; [ -z "${line// }" ] && continue
    case $line in repeat*) continue;; esac
    n=$((n+1))
    stdin=/dev/null
    if [[ $line == *"<"* ]]; then stdin="$WORK/pristine/$(echo "${line#*<}" | xargs)"; line=${line%%<*}; fi
    read -r cmd args <<< "$line"
    ( cd "$WORK/fortran" && OMP_NUM_THREADS=1 timeout 3600 "$FB/$cmd" ./ "$cmd" $args < "$stdin" > "$WORK/logs/f$n.log" 2>&1 ) || echo "  (fortran $cmd exited with $?)"
    ( cd "$WORK/rust" && timeout 3600 "$RBIN" ./ "$cmd" $args --quiet < "$stdin" > "$WORK/logs/r$n.log" 2>&1 ) || { echo "  (rust $cmd exited with $?)"; tail -3 "$WORK/logs/r$n.log"; }
  done < "$EX/commands.txt"
  echo "== $(basename "$EX") pass $k:"
  python3 "$HERE/compare_outputs.py" "$WORK/fortran" "$WORK/rust" --quiet || true
done
if [ $SAVE = 1 ]; then
  ( cd "$WORK/fortran" && files=$(for f in *; do [ -f "$f" ] || continue; if ! cmp -s "$f" "$WORK/pristine/$f" 2>/dev/null; then echo "$f"; fi; done); tar czf "$EX/fortran_reference.tar.gz" $files )
  echo "saved $(tar tzf "$EX/fortran_reference.tar.gz" | wc -l) reference files to $EX/fortran_reference.tar.gz"
fi
echo "work folder: $WORK"
