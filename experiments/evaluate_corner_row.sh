#!/usr/bin/env bash
#
# Evaluate corner-row contraction on every dataset with a precomputed refinement
# in one or more leaf-count groups, then append the results to
# results/cornerrow.csv.
#
# Usage:
#   ./evaluate_corner_row.sh n15 [n25 ...]
#   ./evaluate_corner_row.sh n15,n25
#   ./evaluate_corner_row.sh n15 10-19
#
# Numeric dataset IDs and inclusive ranges apply to the preceding leaf-count
# group. For example, "n15 10-19" evaluates n15/10 through n15/19.
#
# The default parameters are k=0 (meaning the largest legal row sample size,
# n-3), heavy=1, tau=0.3125 (which is (1+delta)/4 at the tuned delta=0.25),
# query-alpha=0.001, and seed=20250729. Each may be overridden through the
# environment, e.g. CORNER_HEAVY=2 CORNER_TAU=0.25 ./evaluate_corner_row.sh n15.
# CORNER_HEAVY=m makes every row draw min(m*k, row population) tuples.
#
# For each dataset, iqtree_500.nwk is preferred; g_500.nwk is a fallback.
#
# Each row records FN and FP counts, the normalized RF, and the two rates
#   fnr = fn / (internal branches of the true tree of blobs)
#   fpr = fp / (internal branches of the estimated tree of blobs)
# together with their unweighted mean and both denominators.
#
# To run this on a SLURM cluster instead, one job per dataset, use
# submit_cornerrow.sh followed by collect_cornerrow.sh.

set -euo pipefail

if [[ $# -eq 0 ]]; then
    sed -n '3,29p' "$0"
    exit 1
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"
PROJECT_ROOT="$(dirname "$SOURCE_DIR")"
DATA_ROOT="$PROJECT_ROOT/data/camus-dataset"
REFINEMENT_ROOT="$PROJECT_ROOT/data/refinements"
BUILD_DIR="$SOURCE_DIR/build"
BINARY="$BUILD_DIR/tree-qmc"
RESULTS_DIR="$PROJECT_ROOT/results"
RESULTS_CSV="$RESULTS_DIR/cornerrow.csv"
TREE_OUTPUT_ROOT="$RESULTS_DIR/cornerrow_trees"

K="${CORNER_K:-0}"
HEAVY="${CORNER_HEAVY:-1}"
TAU="${CORNER_TAU:-0.3125}"
QUERY_ALPHA="${CORNER_QUERY_ALPHA:-0.001}"
SEED="${CORNER_SEED:-20250729}"
METHOD="CORNERROW"
CSV_HEADER="taxa,network_id,true_tob,estimated_tob,method,k,heavy,tau,query_alpha,seed,fn,fp,rf,true_internal,estimated_internal,fnr,fpr,mean_rate,wall_clock_seconds"

if [[ ! "$K" =~ ^[0-9]+$ ]]; then
    echo "ERROR: CORNER_K must be a non-negative integer: $K" >&2
    exit 1
fi
if [[ ! "$HEAVY" =~ ^[0-9]+$ ]] || [[ "$HEAVY" -lt 1 ]]; then
    echo "ERROR: CORNER_HEAVY must be a positive integer: $HEAVY" >&2
    exit 1
fi
if [[ ! "$SEED" =~ ^[0-9]+$ ]]; then
    echo "ERROR: CORNER_SEED must be a non-negative integer: $SEED" >&2
    exit 1
fi

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/cornerrow-evaluation.XXXXXX")"
cleanup() {
    rm -rf -- "$WORKDIR"
}
trap cleanup EXIT

MANIFEST="$WORKDIR/manifest.tsv"
COMPARISONS="$WORKDIR/comparisons.tsv"
GROUPS="$WORKDIR/groups.txt"
FILTERS="$WORKDIR/filters.tsv"
DATASET_DIRS="$WORKDIR/dataset_dirs.txt"

# Accept leaf-count groups plus optional dataset IDs/ranges associated with the
# most recently named group. Commas may be used as separators as well.
: > "$GROUPS"
: > "$FILTERS"
current_group=""
for argument in "$@"; do
    argument="${argument//,/ }"
    for token in $argument; do
        if [[ "$token" =~ ^n[0-9]+$ ]]; then
            current_group="$token"
            if ! grep -qxF "$current_group" "$GROUPS"; then
                printf '%s\n' "$current_group" >> "$GROUPS"
            fi
            continue
        fi

        if [[ ! "$token" =~ ^[0-9]+(-[0-9]+)?$ ]]; then
            echo "ERROR: invalid argument '$token' (expected e.g. n15, 10, or 10-19)" >&2
            exit 1
        fi
        if [[ -z "$current_group" ]]; then
            echo "ERROR: dataset selector '$token' must follow a leaf-count group" >&2
            exit 1
        fi

        range_start="${token%%-*}"
        if [[ "$token" == *-* ]]; then
            range_end="${token#*-}"
        else
            range_end="$range_start"
        fi
        start_number="$((10#$range_start))"
        end_number="$((10#$range_end))"
        if [[ "$start_number" -gt "$end_number" ]]; then
            echo "ERROR: dataset range must be ascending: $token" >&2
            exit 1
        fi

        width="${#range_start}"
        if [[ "${#range_end}" -gt "$width" ]]; then
            width="${#range_end}"
        fi
        for ((dataset_number = start_number; dataset_number <= end_number; dataset_number++)); do
            printf -v selected_id "%0${width}d" "$dataset_number"
            filter_line="${current_group}"$'\t'"${selected_id}"
            if ! grep -qxF "$filter_line" "$FILTERS"; then
                printf '%s\n' "$filter_line" >> "$FILTERS"
            fi
        done
    done
done

if [[ ! -s "$GROUPS" ]]; then
    echo "ERROR: no leaf-count groups were provided" >&2
    exit 1
fi

export PATH="$HOME/.juliaup/bin:/opt/homebrew/bin:$PATH"
if command -v R >/dev/null 2>&1; then
    export R_HOME="$(R RHOME)"
fi

for command_name in cmake julia python3; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $command_name" >&2
        exit 1
    fi
done

if [[ -s "$RESULTS_CSV" ]]; then
    existing_header="$(sed -n '1p' "$RESULTS_CSV")"
    if [[ "$existing_header" != "$CSV_HEADER" ]]; then
        echo "ERROR: $RESULTS_CSV has an unexpected header" >&2
        echo "  expected: $CSV_HEADER" >&2
        echo "  found:    $existing_header" >&2
        echo "Rows written before the FNR/FPR columns were added have the shorter" >&2
        echo "header; move that file aside and re-run to regenerate it." >&2
        exit 1
    fi
fi

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
if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: build did not produce $BINARY" >&2
    exit 1
fi

mkdir -p "$RESULTS_DIR"
if [[ ! -s "$RESULTS_CSV" ]]; then
    printf '%s\n' "$CSV_HEADER" > "$RESULTS_CSV"
fi

echo "Parameters: k=$K heavy=$HEAVY tau=$TAU query-alpha=$QUERY_ALPHA seed=$SEED"

dataset_count=0
appended=0
while IFS= read -r group; do
    dataset_group="$DATA_ROOT/$group"
    refinement_group="$REFINEMENT_ROOT/$group"
    if [[ ! -d "$dataset_group" ]]; then
        echo "ERROR: dataset group does not exist: $dataset_group" >&2
        exit 1
    fi
    if [[ ! -d "$refinement_group" ]]; then
        echo "ERROR: refinement group does not exist: $refinement_group" >&2
        exit 1
    fi

    : > "$DATASET_DIRS"
    if grep -qF "${group}"$'\t' "$FILTERS"; then
        while IFS=$'\t' read -r filter_group filter_id; do
            if [[ "$filter_group" == "$group" ]]; then
                refinement_dataset_dir="$refinement_group/$filter_id"
                if [[ ! -d "$refinement_dataset_dir" ]]; then
                    echo "ERROR: requested refinement dataset does not exist: $group/$filter_id" >&2
                    exit 1
                fi
                printf '%s\n' "$refinement_dataset_dir" >> "$DATASET_DIRS"
            fi
        done < "$FILTERS"
    else
        find "$refinement_group" -mindepth 1 -maxdepth 1 -type d | sort \
            > "$DATASET_DIRS"
    fi

    found_in_group=0
    while IFS= read -r refinement_dataset_dir; do
        found_in_group=1
        network_id="$(basename "$refinement_dataset_dir")"
        dataset="$group/$network_id"
        dataset_dir="$dataset_group/$network_id"
        refinement="$refinement_dataset_dir/astral4-rooted.tre"
        true_network="$dataset_dir/true_net.nwk"

        if [[ -f "$dataset_dir/iqtree_500.nwk" ]]; then
            gene_tree_path="$dataset_dir/iqtree_500.nwk"
        elif [[ -f "$dataset_dir/g_500.nwk" ]]; then
            gene_tree_path="$dataset_dir/g_500.nwk"
        else
            echo "ERROR: neither g_500.nwk nor iqtree_500.nwk exists for $dataset" >&2
            exit 1
        fi
        if [[ ! -f "$refinement" ]]; then
            echo "ERROR: missing refinement for $dataset: $refinement" >&2
            exit 1
        fi
        if [[ ! -f "$true_network" ]]; then
            echo "ERROR: missing true network for $dataset: $true_network" >&2
            exit 1
        fi

        output_dir="$TREE_OUTPUT_ROOT/$group/$network_id"
        inferred_tob="$output_dir/cornerrow.nwk"
        true_tob="$output_dir/true_tob.nwk"
        constructor_log="$output_dir/cornerrow.log"
        mkdir -p "$output_dir"

        dataset_count="$((dataset_count + 1))"
        echo "[$dataset_count] $dataset ($(basename "$gene_tree_path"))"
        start_time="$(python3 -c 'import time; print(time.monotonic())')"
        if ! "$BINARY" -i "$gene_tree_path" \
            --blobsearchonly "$refinement" \
            --blob --cornerrow-blob \
            --corner-k "$K" --heavy-sampling "$HEAVY" \
            --corner-tau "$TAU" --corner-seed "$SEED" \
            --query-alpha "$QUERY_ALPHA" \
            --override -o "$inferred_tob" \
            > "$constructor_log" 2>&1; then
            cat "$constructor_log" >&2
            exit 1
        fi
        end_time="$(python3 -c 'import time; print(time.monotonic())')"
        elapsed_seconds="$(awk -v start="$start_time" -v end="$end_time" \
            'BEGIN { printf "%.6f", end - start }')"

        # Analyze and persist this dataset before starting the next one. Using a
        # one-row manifest retains the comparator's machine-readable interface;
        # its two label columns carry tau and query alpha here.
        printf 'dataset\tdelta\tquery_alpha\ttrue_network\tinferred_tree\ttrue_tob_output\n' \
            > "$MANIFEST"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$dataset" "$TAU" "$QUERY_ALPHA" "$true_network" \
            "$inferred_tob" "$true_tob" >> "$MANIFEST"

        echo "  Computing FN, FP, normalized RF, FNR, and FPR..."
        julia --startup-file=no "$HERE/compare_inferred_tob.jl" \
            --batch --rates "$MANIFEST" "$COMPARISONS"

        exec 3< "$COMPARISONS"
        if ! IFS= read -r comparison_header <&3 || \
            ! IFS=$'\t' read -r comparison_dataset _ _ fn fp rf true_internal \
                estimated_internal fnr fpr mean_rate <&3; then
            echo "ERROR: missing comparison row for $dataset" >&2
            exit 1
        fi
        if [[ "$comparison_header" != $'dataset\tdelta\tquery_alpha\tfn\tfp\tnormalized_rf\ttrue_internal\testimated_internal\tfnr\tfpr\tmean_rate' ]]; then
            echo "ERROR: comparison output has an unexpected header" >&2
            exit 1
        fi
        if IFS= read -r extra_comparison <&3; then
            echo "ERROR: expected one comparison row for $dataset" >&2
            exit 1
        fi
        exec 3<&-
        if [[ "$comparison_dataset" != "$dataset" ]]; then
            echo "ERROR: comparison dataset mismatch: $comparison_dataset != $dataset" >&2
            exit 1
        fi

        taxa="${dataset%%/*}"
        printf '%s,%s,"%s","%s",%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "$taxa" "$network_id" "$true_tob" "$inferred_tob" "$METHOD" \
            "$K" "$HEAVY" "$TAU" "$QUERY_ALPHA" "$SEED" "$fn" "$fp" "$rf" \
            "$true_internal" "$estimated_internal" "$fnr" "$fpr" "$mean_rate" \
            "$elapsed_seconds" \
            >> "$RESULTS_CSV"
        appended="$((appended + 1))"
        echo "  Appended result to $RESULTS_CSV"
    done < "$DATASET_DIRS"

    if [[ "$found_in_group" -eq 0 ]]; then
        echo "ERROR: no refinement datasets found under $refinement_group" >&2
        exit 1
    fi
done < "$GROUPS"
if [[ "$appended" -ne "$dataset_count" ]]; then
    echo "ERROR: expected to append $dataset_count rows, got $appended" >&2
    exit 1
fi

echo "Appended $appended rows to $RESULTS_CSV"
