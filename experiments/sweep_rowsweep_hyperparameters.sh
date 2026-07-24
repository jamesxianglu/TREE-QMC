#!/usr/bin/env bash
#
# Empirically tune row-sweep's delta and query-alpha on a reproducible random
# sample of 5 n15 and 3 n25 CAMUS datasets. The only persistent output is
# comparisons.tsv; refinements, trees, and logs are temporary.
#
# Usage:
#   ./sweep_rowsweep_hyperparameters.sh [OUTPUT_DIR]
#
# Optional environment overrides:
#   SEED=1
#   GENE_TREES=iqtree_500.nwk
#   DELTAS="0.01 0.05 0.10 0.15 0.20 0.25 0.30"
#   QUERY_ALPHAS="0.000001 0.000005 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05"

set -euo pipefail

if [[ $# -gt 1 ]]; then
    sed -n '3,15p' "$0"
    exit 1
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"
DATA_ROOT="$(dirname "$SOURCE_DIR")/data/camus-dataset"
BUILD_DIR="$SOURCE_DIR/build"
BINARY="$BUILD_DIR/tree-qmc"

SEED="${SEED:-1}"
GENE_TREES="${GENE_TREES:-iqtree_500.nwk}"
DELTAS="${DELTAS:-0.01 0.05 0.10 0.15 0.20 0.25 0.30}"
QUERY_ALPHAS="${QUERY_ALPHAS:-0.000001 0.000005 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05}"
OUTDIR="${1:-$HERE/results/rowsweep_hyperparameter_sweep/seed${SEED}}"
COMPARISONS="$OUTDIR/comparisons.tsv"

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/tree-qmc-sweep.XXXXXX")"
cleanup() {
    rm -rf -- "$WORKDIR"
}
trap cleanup EXIT

SELECTED="$WORKDIR/selected_datasets.txt"
MANIFEST="$WORKDIR/manifest.tsv"

export PATH="$HOME/.juliaup/bin:/opt/homebrew/bin:$PATH"
if command -v R >/dev/null 2>&1; then
    export R_HOME="$(R RHOME)"
fi

sample_group() {
    local group="$1"
    local count="$2"
    local seed="$3"

    find "$DATA_ROOT/$group" -mindepth 1 -maxdepth 1 -type d | sort | \
        awk -v group="$group" -v gene_trees="$GENE_TREES" \
            -v count="$count" -v seed="$seed" '
            BEGIN { srand(seed) }
            {
                true_network = $0 "/true_net.nwk"
                gene_tree_file = $0 "/" gene_trees
                check_true = "test -f \"" true_network "\""
                check_genes = "test -f \"" gene_tree_file "\""
                if (system(check_true) != 0 || system(check_genes) != 0) next
                seen++
                component_count = split($0, path_components, "/")
                candidate = group "/" path_components[component_count]
                if (seen <= count) {
                    chosen[seen] = candidate
                } else {
                    replacement = int(rand() * seen) + 1
                    if (replacement <= count) chosen[replacement] = candidate
                }
            }
            END {
                if (seen < count) exit 2
                for (i = 1; i <= count; i++) print chosen[i]
            }
        '
}

: > "$SELECTED"
sample_group n15 5 "$SEED" >> "$SELECTED"
sample_group n25 3 "$((SEED + 1))" >> "$SELECTED"

if [[ "$(wc -l < "$SELECTED" | tr -d ' ')" != "8" ]]; then
    echo "ERROR: expected 8 selected datasets" >&2
    exit 1
fi

echo "Selected datasets (seed=$SEED):"
sed 's/^/  /' "$SELECTED"

echo "Configuring TREE-QMC..."
if ! cmake -S "$SOURCE_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release \
    > "$WORKDIR/configure.log" 2>&1; then
    cat "$WORKDIR/configure.log" >&2
    exit 1
fi
echo "Building TREE-QMC..."
if ! cmake --build "$BUILD_DIR" -j4 > "$WORKDIR/build.log" 2>&1; then
    cat "$WORKDIR/build.log" >&2
    exit 1
fi

printf 'dataset\tdelta\tquery_alpha\ttrue_network\tinferred_tree\ttrue_tob_output\n' > "$MANIFEST"
delta_count="$(awk '{print NF}' <<< "$DELTAS")"
alpha_count="$(awk '{print NF}' <<< "$QUERY_ALPHAS")"
total_runs="$((8 * delta_count * alpha_count))"
run_index=0

while IFS= read -r dataset; do
    dataset_dir="$DATA_ROOT/$dataset"
    dataset_tag="${dataset//\//_}"
    dataset_work="$WORKDIR/$dataset_tag"
    gene_tree_path="$dataset_dir/$GENE_TREES"
    true_network="$dataset_dir/true_net.nwk"
    refinement="$dataset_work/refinement.nwk"
    true_tob="$dataset_work/true_tob.nwk"
    mkdir -p "$dataset_work"

    echo "[$dataset] estimating refinement..."
    if ! "$BINARY" -i "$gene_tree_path" --override -o "$refinement" \
        > "$dataset_work/refinement.log" 2>&1; then
        cat "$dataset_work/refinement.log" >&2
        exit 1
    fi

    for delta in $DELTAS; do
        for query_alpha in $QUERY_ALPHAS; do
            tag="d${delta}_a${query_alpha}"
            inferred_tob="$dataset_work/tob_${tag}.nwk"
            constructor_log="$dataset_work/${tag}.log"

            run_index="$((run_index + 1))"
            echo "[$run_index/$total_runs] dataset=$dataset delta=$delta query_alpha=$query_alpha"
            if ! "$BINARY" -i "$gene_tree_path" \
                --blobsearchonly "$refinement" \
                --blob --rowsweep-blob \
                --delta "$delta" --query-alpha "$query_alpha" \
                --override -o "$inferred_tob" \
                > "$constructor_log" 2>&1; then
                cat "$constructor_log" >&2
                exit 1
            fi

            printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$dataset" "$delta" "$query_alpha" "$true_network" \
                "$inferred_tob" "$true_tob" >> "$MANIFEST"
        done
    done
done < "$SELECTED"

mkdir -p "$OUTDIR"
echo "Computing FN, FP, and normalized RF..."
julia --startup-file=no "$HERE/compare_inferred_tob.jl" \
    --batch "$MANIFEST" "$COMPARISONS"

echo "Saved: $COMPARISONS"
