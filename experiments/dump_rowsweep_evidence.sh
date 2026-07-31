#!/usr/bin/env bash
#
# Record the quartet evidence the row sweep consumes, for offline study of
# alternative decision rules. This runs the same queries as evaluate_rowsweep.sh
# but stores qCF counts and T1 p-values instead of contracting anything, so it
# costs about one row-sweep run per replicate and changes no result.
#
# Usage:
#   ./dump_rowsweep_evidence.sh n15 [n25 ...]
#   ./dump_rowsweep_evidence.sh n50 0-9        # IDs/ranges follow their group
#
# Environment:
#   ALL_TARGETS=1   also fill in the two T1 p-values the sweep never asks for
#                   (roughly 2x slower; only affordable for small groups)
#   ALL_ANCHORS=1   query every pair from the smaller side, not just the pairs
#                   containing the arbitrary anchor R[0]; lets an offline rule
#                   re-anchor or average over anchors, at a factor |R|/2 cost
#   QUERY_ALPHA     recorded in the dump; p-values themselves are alpha-free
#   GENE_TREES      basename of the gene tree file to use. The default mirrors
#                   evaluate_rowsweep.sh (iqtree_500.nwk, else g_500.nwk), but
#                   results/rowsweep.csv was produced with iqtree_500.nwk only
#                   for n15 and g_500.nwk from n25 up, so set this explicitly
#                   when a dump must match a published row.
#
set -euo pipefail

if [[ $# -eq 0 ]]; then
    sed -n '3,14p' "$0"
    exit 1
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"
PROJECT_ROOT="$(dirname "$SOURCE_DIR")"
DATA_ROOT="$PROJECT_ROOT/data/camus-dataset"
REFINEMENT_ROOT="$PROJECT_ROOT/data/refinements"
BINARY="$SOURCE_DIR/build/tree-qmc"
EVIDENCE_ROOT="${EVIDENCE_ROOT:-$PROJECT_ROOT/results/rowsweep_evidence}"

QUERY_ALPHA="${QUERY_ALPHA:-0.001}"
ALL_TARGETS_FLAG=""
[[ "${ALL_TARGETS:-0}" == "1" ]] && ALL_TARGETS_FLAG="--rowsweep-dump-all-targets"
ALL_ANCHORS_FLAG=""
[[ "${ALL_ANCHORS:-0}" == "1" ]] && ALL_ANCHORS_FLAG="--rowsweep-dump-all-anchors"

export PATH="$HOME/.juliaup/bin:/opt/homebrew/bin:$PATH"
if command -v R >/dev/null 2>&1; then
    export R_HOME="$(R RHOME)"
fi
[[ -x "$BINARY" ]] || { echo "ERROR: build tree-qmc first ($BINARY)" >&2; exit 1; }

# Expand "group [id|range ...]" arguments into a work list.
WORK="$(mktemp)"
trap 'rm -f "$WORK"' EXIT
group=""
for argument in "$@"; do
    for token in ${argument//,/ }; do
        if [[ "$token" =~ ^n[0-9]+$ ]]; then
            group="$token"
            # A bare group means every replicate that has a refinement.
            find "$REFINEMENT_ROOT/$group" -mindepth 1 -maxdepth 1 -type d \
                -exec basename {} \; | sort | sed "s|^|$group	|" >> "$WORK"
            continue
        fi
        [[ -n "$group" ]] || { echo "ERROR: '$token' must follow a group" >&2; exit 1; }
        # An explicit selector replaces the group's default expansion.
        grep -v "^$group	" "$WORK" > "$WORK.tmp" || true
        mv "$WORK.tmp" "$WORK"
        start="${token%%-*}"; end="${token#*-}"; [[ "$token" == *-* ]] || end="$start"
        for ((i = 10#$start; i <= 10#$end; i++)); do
            printf '%s\t%0*d\n' "$group" "${#start}" "$i" >> "$WORK"
        done
    done
done

count=0
while IFS=$'\t' read -r group rep; do
    dataset_dir="$DATA_ROOT/$group/$rep"
    refinement="$REFINEMENT_ROOT/$group/$rep/astral4-rooted.tre"
    if [[ -n "${GENE_TREES:-}" ]]; then
        gene_trees="$dataset_dir/$GENE_TREES"
        [[ -f "$gene_trees" ]] || { echo "ERROR: no $GENE_TREES for $group/$rep" >&2; exit 1; }
    elif [[ -f "$dataset_dir/iqtree_500.nwk" ]]; then
        gene_trees="$dataset_dir/iqtree_500.nwk"
    elif [[ -f "$dataset_dir/g_500.nwk" ]]; then
        gene_trees="$dataset_dir/g_500.nwk"
    else
        echo "ERROR: no gene trees for $group/$rep" >&2
        exit 1
    fi
    [[ -f "$refinement" ]] || { echo "ERROR: no refinement for $group/$rep" >&2; exit 1; }

    out_dir="$EVIDENCE_ROOT/$group/$rep"
    mkdir -p "$out_dir"
    if [[ -s "$out_dir/dump.quartets.bin" && -s "$out_dir/dump.meta.json" ]]; then
        echo "skip $group/$rep (already dumped)"
        continue
    fi

    count=$((count + 1))
    echo "[$count] dumping $group/$rep"
    start_time="$(python3 -c 'import time; print(time.monotonic())')"
    if ! "$BINARY" -i "$gene_trees" \
        --blobsearchonly "$refinement" \
        --blob --rowsweep-blob \
        --query-alpha "$QUERY_ALPHA" \
        --rowsweep-dump "$out_dir/dump" $ALL_TARGETS_FLAG $ALL_ANCHORS_FLAG \
        --override -o /dev/null > "$out_dir/dump.log" 2>&1; then
        tail -20 "$out_dir/dump.log" >&2
        exit 1
    fi
    end_time="$(python3 -c 'import time; print(time.monotonic())')"
    awk -v s="$start_time" -v e="$end_time" 'BEGIN { printf "  %.1fs\n", e - s }'
done < "$WORK"

echo "Dumped $count replicates to $EVIDENCE_ROOT"
