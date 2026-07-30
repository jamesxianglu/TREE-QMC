#!/usr/bin/env bash
#
# Merge the per-dataset CSV shards written by run_cornerrow.sbatch into a single
# results/cornerrow.csv, and report any datasets that failed or never finished.
#
# Usage:
#   ./collect_cornerrow.sh [n15 ...]
#
# With no arguments, every group under results/cornerrow_trees is collected.
# This is safe to run repeatedly and while jobs are still in flight; it always
# rewrites the merged CSV from whatever shards currently exist.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"
PROJECT_ROOT="$(dirname "$SOURCE_DIR")"
RESULTS_DIR="$PROJECT_ROOT/results"
TREE_OUTPUT_ROOT="$RESULTS_DIR/cornerrow_trees"
RESULTS_CSV="$RESULTS_DIR/cornerrow.csv"

CSV_HEADER="taxa,network_id,true_tob,estimated_tob,method,k,heavy,tau,query_alpha,seed,fn,fp,rf,true_internal,estimated_internal,fnr,fpr,mean_rate,wall_clock_seconds,max_rss_kb,node"

if [[ ! -d "$TREE_OUTPUT_ROOT" ]]; then
    echo "ERROR: no output tree found at $TREE_OUTPUT_ROOT" >&2
    exit 1
fi

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/cornerrow-collect.XXXXXX")"
trap 'rm -rf -- "$WORKDIR"' EXIT
SHARDS="$WORKDIR/shards.txt"

if [[ $# -gt 0 ]]; then
    : > "$SHARDS"
    for group in "$@"; do
        group_dir="$TREE_OUTPUT_ROOT/$group"
        if [[ ! -d "$group_dir" ]]; then
            echo "ERROR: no results for group '$group'" >&2
            exit 1
        fi
        find "$group_dir" -mindepth 2 -maxdepth 2 -name cornerrow.csv | sort >> "$SHARDS"
    done
    search_roots=()
    for group in "$@"; do search_roots+=( "$TREE_OUTPUT_ROOT/$group" ); done
    # Dataset directories sit one level below a group directory, two below the
    # tree root. Getting this depth wrong silently drops PENDING datasets from
    # the completeness report instead of erroring.
    dataset_depth=1
else
    find "$TREE_OUTPUT_ROOT" -mindepth 3 -maxdepth 3 -name cornerrow.csv | sort > "$SHARDS"
    search_roots=( "$TREE_OUTPUT_ROOT" )
    dataset_depth=2
fi

if [[ ! -s "$SHARDS" ]]; then
    echo "ERROR: no result shards found; are the jobs still running?" >&2
    exit 1
fi

# Never silently overwrite a merged CSV whose columns are not the ones we are
# about to write: results/cornerrow.csv predates the FNR/FPR columns, and the
# local evaluate_corner_row.sh writes a same-named file.
if [[ -s "$RESULTS_CSV" ]]; then
    existing_header="$(sed -n '1p' "$RESULTS_CSV")"
    if [[ "$existing_header" != "$CSV_HEADER" ]]; then
        echo "ERROR: $RESULTS_CSV has a different header and would be overwritten" >&2
        echo "  expected: $CSV_HEADER" >&2
        echo "  found:    $existing_header" >&2
        echo "Move the old file aside (e.g. to cornerrow-legacy.csv) and re-run." >&2
        exit 1
    fi
fi

mkdir -p "$RESULTS_DIR"
printf '%s\n' "$CSV_HEADER" > "$RESULTS_CSV"

merged=0
while IFS= read -r shard; do
    shard_header="$(sed -n '1p' "$shard")"
    if [[ "$shard_header" != "$CSV_HEADER" ]]; then
        echo "ERROR: $shard has an unexpected header" >&2
        echo "  expected: $CSV_HEADER" >&2
        echo "  found:    $shard_header" >&2
        exit 1
    fi
    row_count="$(( $(wc -l < "$shard") - 1 ))"
    if [[ "$row_count" -ne 1 ]]; then
        echo "ERROR: $shard has $row_count data rows, expected 1" >&2
        exit 1
    fi
    sed -n '2p' "$shard" >> "$RESULTS_CSV"
    merged="$((merged + 1))"
done < "$SHARDS"

echo "Merged $merged result(s) into $RESULTS_CSV"

# Report datasets that started but did not produce a result.
failed=0
while IFS= read -r marker; do
    dataset_dir="$(dirname "$marker")"
    dataset="$(basename "$(dirname "$dataset_dir")")/$(basename "$dataset_dir")"
    echo "  FAILED  $dataset: $(cat "$marker")"
    failed="$((failed + 1))"
done < <(find "${search_roots[@]}" -name cornerrow.failed | sort)

missing=0
while IFS= read -r dataset_dir; do
    if [[ ! -s "$dataset_dir/cornerrow.csv" && ! -f "$dataset_dir/cornerrow.failed" ]]; then
        dataset="$(basename "$(dirname "$dataset_dir")")/$(basename "$dataset_dir")"
        echo "  PENDING $dataset (queued, running, or killed by SLURM)"
        missing="$((missing + 1))"
    fi
done < <(find "${search_roots[@]}" -mindepth "$dataset_depth" -maxdepth "$dataset_depth" -type d | sort)

if [[ "$failed" -gt 0 || "$missing" -gt 0 ]]; then
    echo
    echo "$failed failed, $missing incomplete. Re-run submit_cornerrow.sh to fill the gaps;"
    echo "it skips datasets that already have a result."
fi
