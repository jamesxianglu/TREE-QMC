#!/usr/bin/env bash
#
# Sweep the corner-row resolution-margin threshold and write runs that
# results/analysis/edge_metrics.py can score directly.
#
# Usage:  ./sweep_resolution.sh <group> <margin>[,<margin>...] [nproc]
# Example: ./sweep_resolution.sh n25 0,0.005,0.01,0.02,0.04 6
#
# A margin of 0 disables the test and reproduces the shipped corner-row run,
# which is the paired baseline every other column is compared against.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"
ROOT="$(dirname "$SOURCE_DIR")"
BINARY="$SOURCE_DIR/build/tree-qmc"

GROUP="${1:?group, e.g. n25}"
MARGINS="${2:?comma-separated margins}"
NPROC="${3:-6}"

ORACLE="${ORACLE:-t1+cf|maj}"
TAU="${TAU:-0.15,0.60}"
ALPHA="${ALPHA:-0.001}"
SAMPLES="${SAMPLES:-20}"

export R_HOME="${R_HOME:-$(R RHOME)}"
export PATH="$HOME/.juliaup/bin:$PATH"

DATA="$ROOT/data/camus-dataset/$GROUP"
REF="$ROOT/data/refinements/$GROUP"
# Per-dataset, matching run.meta's "per-dataset (iqtree_500.nwk, else g_500.nwk)":
# n25 carries iqtree_500 for only some replicates.

run_one() {
    local rep="$1" margin="$2" label="$3"
    local out="$ROOT/results/runs/$label/$GROUP/$rep"
    mkdir -p "$out"
    [[ -s "$out/estimate.nwk" ]] && return 0
    local gt="iqtree_500.nwk"
    [[ -f "$DATA/$rep/$gt" ]] || gt="g_500.nwk"
    local flag=()
    [[ "$margin" != "0" ]] && flag=(--cornerrow-resolution-margin "$margin")
    local t0=$SECONDS
    "$BINARY" -i "$DATA/$rep/$gt" \
        --blobsearchonly "$REF/$rep/astral4-rooted.tre" \
        --blob --cornerrow-blob --oracle "$ORACLE" --corner-tau "$TAU" \
        "${flag[@]}" --resolution-samples "$SAMPLES" \
        --query-alpha "$ALPHA" --override -o "$out/estimate.nwk" \
        > "$out/run.log" 2>&1
    local rc=$?
    echo "wall_clock_seconds" > "$out/row.csv"
    echo "$((SECONDS - t0))" >> "$out/row.csv"
    [[ $rc -ne 0 ]] && echo "FAILED $GROUP/$rep margin=$margin rc=$rc" >&2
    return 0
}
export -f run_one
export ROOT BINARY GROUP DATA REF ORACLE TAU ALPHA SAMPLES

for margin in ${MARGINS//,/ }; do
    label="cr_res_${margin//./p}"
    mkdir -p "$ROOT/results/runs/$label"
    cat > "$ROOT/results/runs/$label/run.meta" <<EOF
method=cornerrow
args=--oracle $ORACLE --corner-tau $TAU --cornerrow-resolution-margin $margin --resolution-samples $SAMPLES --query-alpha $ALPHA
gene_trees=per-dataset (iqtree_500.nwk, else g_500.nwk)
refinement=data/refinements/<group>/<rep>/astral4-rooted.tre
runner=local
label=$label
EOF
    echo "=== $label on $GROUP (${NPROC} way) ==="
    ls "$REF" | sort | xargs -P "$NPROC" -I{} bash -c 'run_one "$@"' _ {} "$margin" "$label"
done
echo "DONE $GROUP $MARGINS"
