#!/usr/bin/env bash
# Run one validation case with the original Fortran EPSR programs and with the
# Rust port (epsr-rs), on separate copies of the case folder, and compare every
# file in the two folders after each step.
#
# usage: run_validation.sh CASE_DIR FORTRAN_BUILD_DIR WORK_DIR
#
#   CASE_DIR           folder with the input files and steps.txt (e.g. validation/water)
#   FORTRAN_BUILD_DIR  built by rust/tools/build_fortran_reference.sh
#   WORK_DIR           created/overwritten; gets fortran/, rust/, logs/ and summary.txt
#
# steps.txt, one step per line (# starts a comment):
#   run <command> [args...] [< stdin_file]   run an EPSR program in both folders:
#                                            Fortran: FORTRAN_BUILD_DIR/<command> ./ <command> args
#                                            Rust:    epsr-rs ./ <command> args --quiet
#   set <file> <keyword> <value...>          edit a keyword in both folders (set_keyword.py)
#   delete <file...>                         remove files from both folders
# After every step all files of the two folders are compared (rust/tools/compare_outputs.py).
# Exit status: 0 if the folders were identical after every step.
set -uo pipefail
CASE=$(realpath "$1"); FB=$(realpath "$2"); WORK=$3
HERE=$(cd "$(dirname "$0")" && pwd)
RBIN="$HERE/../rust/target/release/epsr-rs"
CMP="$HERE/../rust/tools/compare_outputs.py"
[ -x "$RBIN" ] || { echo "build the Rust port first: (cd rust && cargo build --release)"; exit 2; }
rm -rf "$WORK"; mkdir -p "$WORK/fortran" "$WORK/rust" "$WORK/logs"
WORK=$(realpath "$WORK")
for d in fortran rust; do
  (cd "$CASE" && find . -maxdepth 1 -type f ! -name steps.txt ! -name README.md -exec cp {} "$WORK/$d/" \;)
done
ulimit -s unlimited || true
SUMMARY="$WORK/summary.txt"; : > "$SUMMARY"
n=0; fail=0
while IFS= read -r line || [ -n "$line" ]; do
  line=${line%%#*}; [ -z "${line// }" ] && continue
  n=$((n+1))
  read -r kind rest <<< "$line"
  case $kind in
    run)
      stdin=/dev/null; cmdline=$rest
      if [[ $rest == *"<"* ]]; then stdin="$CASE/$(echo "${rest#*<}" | xargs)"; cmdline=${rest%%<*}; fi
      read -r cmd args <<< "$cmdline"
      t0=$(date +%s.%N)
      ( cd "$WORK/fortran" && OMP_NUM_THREADS=1 "$FB/$cmd" ./ "$cmd" $args < "$stdin" > "$WORK/logs/step$n.fortran.log" 2>&1 ) \
        || echo "  (fortran $cmd exited with $?)" | tee -a "$SUMMARY"
      t1=$(date +%s.%N)
      ( cd "$WORK/rust" && "$RBIN" ./ "$cmd" $args --quiet < "$stdin" > "$WORK/logs/step$n.rust.log" 2>&1 ) \
        || echo "  (rust $cmd exited with $?)" | tee -a "$SUMMARY"
      t2=$(date +%s.%N)
      timing=$(printf "fortran %.1fs, rust %.1fs" "$(echo "$t1 - $t0" | bc)" "$(echo "$t2 - $t1" | bc)")
      ;;
    set)
      for d in fortran rust; do python3 "$HERE/set_keyword.py" "$WORK/$d/"$rest || fail=1; done
      timing="edit"
      ;;
    delete)
      for d in fortran rust; do (cd "$WORK/$d" && rm -f $rest); done
      timing="edit"
      ;;
    *) echo "unknown step: $line"; exit 2;;
  esac
  python3 "$CMP" "$WORK/fortran" "$WORK/rust" --quiet > "$WORK/logs/step$n.compare.txt" 2>&1
  status=$(tail -1 "$WORK/logs/step$n.compare.txt")
  [ "$status" = "ALL IDENTICAL" ] || fail=1
  printf "step %2d  %-55s %-28s %s\n" "$n" "$kind $rest" "($timing)" "$status" | tee -a "$SUMMARY"
  [ "$status" = "ALL IDENTICAL" ] || sed '$d' "$WORK/logs/step$n.compare.txt" | sed 's/^/    /' | tee -a "$SUMMARY"
done < "$CASE/steps.txt"
echo "== $(basename "$CASE"): $([ $fail = 0 ] && echo 'ALL STEPS IDENTICAL' || echo 'DIFFERENCES FOUND')" | tee -a "$SUMMARY"
exit $fail
