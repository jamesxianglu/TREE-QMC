#!/usr/bin/env bash
#
# Assess row_sweep_test_idx on one dataset and report where it is weak.
#
# Usage:
#   ./run_split_test.sh <dataset-path> <delta[,delta2,...]> [gene_trees] [extra flags]
#
# Examples:
#   ./run_split_test.sh n25/10 0.15
#   ./run_split_test.sh n15/00 0.05,0.15,0.3 g_true.nwk
#   ./run_split_test.sh n25/10 0.15 iqtree_500.nwk --neg-mult 4 --seed 7
#
# <dataset-path> is relative to data/camus-dataset (e.g. n25/10).
# gene_trees defaults to iqtree_500.nwk.
#
set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Make julia (PhyloNetworks) and R (embedded in tree-qmc) reachable.
export PATH="$HOME/.juliaup/bin:/opt/homebrew/bin:$PATH"
if command -v R >/dev/null 2>&1; then
    export R_HOME="$(R RHOME)"
fi

if [ "$#" -lt 2 ]; then
    sed -n '3,18p' "$0"
    exit 1
fi

# Ensure the tree-qmc binary exists; build it into ../build if missing.
SRC="$(dirname "$HERE")"
BIN="$SRC/build/tree-qmc"
if [ ! -x "$BIN" ]; then
    echo "tree-qmc binary not found at $BIN -- building..."
    cmake -S "$SRC" -B "$SRC/build" -DCMAKE_BUILD_TYPE=Release
    cmake --build "$SRC/build" -j4
    if [ ! -x "$BIN" ]; then
        echo "ERROR: build did not produce $BIN" >&2
        exit 1
    fi
    echo "Built $BIN"
fi

DATASET="$1"; shift
DELTAS="$1"; shift
GT="iqtree_500.nwk"
if [ "$#" -gt 0 ] && [ "${1#-}" = "$1" ]; then GT="$1"; shift; fi   # optional gene-tree file

IFS=',' read -ra DARR <<< "$DELTAS"
for d in "${DARR[@]}"; do
    python3 "$HERE/analyze_split.py" "$DATASET" "$d" "$GT" "$@"
    echo
done
