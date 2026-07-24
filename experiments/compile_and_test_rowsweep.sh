#!/usr/bin/env bash
#
# Build TREE-QMC, evaluate row_sweep_test_idx as a split classifier, construct
# a tree of blobs from a binary refinement, and compare it with the true ToB.
#
# Usage:
#   ./compile_and_test_rowsweep.sh DATASET DELTA [QUERY_ALPHA] [GENE_TREES] [REFINEMENT]
#
# Example:
#   ./compile_and_test_rowsweep.sh n15/00 0.15 0.0005
#   ./compile_and_test_rowsweep.sh n15/00 0.15 0.0005 g_true.nwk my_refinement.nwk
#
# DATASET is relative to data/camus-dataset. If REFINEMENT is omitted, the
# script first estimates a binary species tree with TREE-QMC and uses it as the
# refinement tested by the row-sweep constructor.

set -euo pipefail

if [[ $# -lt 2 || $# -gt 5 ]]; then
    sed -n '3,16p' "$0"
    exit 1
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"
DATA_ROOT="$(dirname "$SOURCE_DIR")/data/camus-dataset"
BUILD_DIR="$SOURCE_DIR/build"
BINARY="$BUILD_DIR/tree-qmc"

DATASET="$1"
DELTA="$2"
QUERY_ALPHA="${3:-0.0005}"
GENE_TREES="${4:-iqtree_500.nwk}"
PROVIDED_REFINEMENT="${5:-}"
DATASET_DIR="$DATA_ROOT/$DATASET"
GENE_TREE_PATH="$DATASET_DIR/$GENE_TREES"
TRUE_NETWORK="$DATASET_DIR/true_net.nwk"

if [[ ! -f "$GENE_TREE_PATH" || ! -f "$TRUE_NETWORK" ]]; then
    echo "ERROR: missing dataset files under $DATASET_DIR" >&2
    exit 1
fi
if [[ -n "$PROVIDED_REFINEMENT" && ! -f "$PROVIDED_REFINEMENT" ]]; then
    echo "ERROR: refinement does not exist: $PROVIDED_REFINEMENT" >&2
    exit 1
fi

export PATH="$HOME/.juliaup/bin:/opt/homebrew/bin:$PATH"
if command -v R >/dev/null 2>&1; then
    export R_HOME="$(R RHOME)"
fi

TAG="${DATASET//\//_}_d${DELTA}_a${QUERY_ALPHA}"
OUTDIR="$HERE/results/compile_test/$TAG"
CLASSIFICATION_DIR="$OUTDIR/classification"
mkdir -p "$CLASSIFICATION_DIR"

echo "== Configuring and compiling TREE-QMC =="
cmake -S "$SOURCE_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" -j4

if [[ -n "$PROVIDED_REFINEMENT" ]]; then
    REFINEMENT="$PROVIDED_REFINEMENT"
    echo "Using supplied refinement: $REFINEMENT"
else
    REFINEMENT="$OUTDIR/refinement.nwk"
    echo "Estimating binary TREE-QMC refinement: $REFINEMENT"
    "$BINARY" -i "$GENE_TREE_PATH" --override -o "$REFINEMENT" \
        > "$OUTDIR/refinement.log" 2>&1
fi

INFERRED_TOB="$OUTDIR/rowsweep_tob.nwk"
echo "Constructing row-sweep ToB: $INFERRED_TOB"
"$BINARY" -i "$GENE_TREE_PATH" \
    --blobsearchonly "$REFINEMENT" \
    --blob --rowsweep-blob \
    --delta "$DELTA" --query-alpha "$QUERY_ALPHA" \
    --override -o "$INFERRED_TOB" \
    > "$OUTDIR/constructor.log" 2>&1

julia --startup-file=no "$HERE/compare_inferred_tob.jl" \
    "$TRUE_NETWORK" "$INFERRED_TOB" \
    "$OUTDIR/true_tob.nwk"

echo
echo "-- inferred tree --"
sed -n '1p' "$INFERRED_TOB"
echo
echo "Results: $OUTDIR"
