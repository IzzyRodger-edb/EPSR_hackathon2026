#!/usr/bin/env bash
# Run every EPSR-test-ex* verification case (or those given) against the
# Fortran reference programs and print one summary line per case.
#
# usage: run_all_ex_tests.sh FORTRAN_BUILD_DIR [WORK_DIR] [case numbers...]
#   (cases are the folders test_problems/<system>/EPSR-test-ex<N> of the repository;
#    $EPSR_ROOT overrides the repository root)
#   e.g. run_all_ex_tests.sh /tmp/epsr-fortran              (all cases)
#        run_all_ex_tests.sh /tmp/epsr-fortran /tmp/w 1 2 50 (selected cases)
# Exit status 0 if every pass of every case is ALL IDENTICAL.
set -uo pipefail
FB=$(realpath "$1"); WORK=${2:-$(mktemp -d)}; shift; [ $# -gt 0 ] && shift
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=${EPSR_ROOT:-$(cd "$HERE/../.." && pwd)}  # repository root
TP="$ROOT/test_problems"
if [ $# -gt 0 ]; then
  cases=$(for n in "$@"; do d=$(ls -d "$TP"/*/EPSR-test-ex"$n" 2>/dev/null | head -1)
          if [ -n "$d" ]; then echo "$d"; else echo "no test case ex$n under $TP" >&2; fi; done)
else cases=$(ls -d "$TP"/*/EPSR-test-ex* | sort -V); fi
fail=0; total=0
for d in $cases; do
  [ -f "$d/commands.txt" ] || continue
  name=$(basename "$d"); total=$((total+1))
  out=$("$HERE/run_ex_test.sh" "$d" "$FB" "$WORK/$name" 2>&1)
  passes=$(echo "$out" | grep -c '^== ')
  ok=$(echo "$out" | grep -c '^ALL IDENTICAL')
  if [ "$passes" -gt 0 ] && [ "$ok" = "$passes" ]; then
    echo "$name: ALL IDENTICAL ($passes pass(es))"
  else
    echo "$name: DIFFERENCES ($ok/$passes identical) — see $WORK/$name"; fail=$((fail+1))
    echo "$out" | grep -v "^ALL IDENTICAL\|^== \|^work folder" | head -8 | sed 's/^/    /'
  fi
done
echo "$((total-fail))/$total cases identical"
[ $fail = 0 ]
