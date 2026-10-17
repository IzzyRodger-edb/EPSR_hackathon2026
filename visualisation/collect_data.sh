#!/usr/bin/env bash
# Record a validation case with the per-iteration recorder (or reuse a finished
# recording) and copy what the visualisation tools need into visualisation/data/<case>/.
#
# usage: visualisation/collect_data.sh CASE [WORK_DIR] [VALIDATED_RUST_DIR]
#   CASE                water | silica_bulk   (a folder in validation/)
#   WORK_DIR            where the recorder runs (default ~/epsr-work/rec-<CASE>);
#                       if WORK_DIR/trajectory already exists the recording is reused
#   VALIDATED_RUST_DIR  optional: the Rust folder of validation/run_validation.sh
#                       (e.g. ~/epsr-work/val-water/rust); the recorder's EPSR output
#                       files are compared with it (they must be identical)
#
# Takes ~10 min (water) / ~30 min (silica) when it has to record.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(dirname "$here")"
case_name="$1"
work="${2:-$HOME/epsr-work/rec-$case_name}"
validated="${3:-}"

if [[ ! -f "$work/trajectory/history.tsv" ]]; then
    (cd "$here/recorder" && cargo build --release)
    ulimit -s unlimited
    OMP_NUM_THREADS=1 "$here/recorder/target/release/epsr-recorder" "$root/validation/$case_name" "$work"
fi

if [[ -n "$validated" ]]; then
    echo "== comparing the recorder's EPSR files with $validated"
    python3 "$root/rust/tools/compare_outputs.py" "$validated" "$work" --quiet
fi

dest="$here/data/$case_name"
rm -rf "$dest"
mkdir -p "$dest/epsr_output"
cp -r "$work/trajectory" "$dest/"
# final EPSR output of the last run: result tables, history, configuration, potentials, input
for f in "$work"/*.EPSR.* "$work"/*.ato "$work"/*.pcof "$work"/*.wts "$work"/*.NWTS.dat; do
    [[ -f "$f" ]] && cp "$f" "$dest/epsr_output/"
done
echo "copied to $dest ($(du -sh "$dest" | cut -f1))"
