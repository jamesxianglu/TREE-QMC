#!/usr/bin/env bash
#
# Submit one SLURM job per dataset to evaluate corner-row contraction on every
# dataset with a precomputed refinement in one or more leaf-count groups. Each
# job writes its own CSV shard; run collect_cornerrow.sh afterwards to merge them
# into results/cornerrow.csv.
#
# Usage:
#   ./submit_cornerrow.sh [--dry-run] [--force] [--no-build] n15 [n25 ...]
#   ./submit_cornerrow.sh n15,n25
#   ./submit_cornerrow.sh n15 10-19
#
# Numeric dataset IDs and inclusive ranges apply to the preceding leaf-count
# group. For example, "n15 10-19" evaluates n15/10 through n15/19.
#
# Options:
#   --dry-run   print the sbatch commands without submitting
#   --force     re-submit datasets that already have a result CSV
#   --no-build  skip cmake configure/build and reuse the existing binary
#
# The default parameters are k=0 (meaning the largest legal row sample size,
# n-3), heavy=1, tau=0.3125 (which is (1+delta)/4 at the tuned delta=0.25),
# query-alpha=0.001, and seed=20250729. Each may be overridden through the
# environment, e.g. CORNER_HEAVY=2 ./submit_cornerrow.sh n50. Note that heavy
# sampling has been measured to change nothing: on n15/00-04, heavy=2 and
# heavy=1000 give bit-identical FN/FP/RF to heavy=1.
#
# For each dataset, iqtree_500.nwk is preferred; g_500.nwk is used as a fallback.

set -euo pipefail

if [[ $# -eq 0 ]]; then
    sed -n '3,28p' "$0"
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
TREE_OUTPUT_ROOT="$RESULTS_DIR/cornerrow_trees"
LOG_ROOT="$RESULTS_DIR/cornerrow_logs"
BIN_SNAPSHOT_ROOT="$RESULTS_DIR/cornerrow_bin"
SBATCH_SCRIPT="$HERE/run_cornerrow.sbatch"
COMPARE_JL="$HERE/compare_inferred_tob.jl"

K="${CORNER_K:-0}"
HEAVY="${CORNER_HEAVY:-1}"
TAU="${CORNER_TAU:-0.3125}"
QUERY_ALPHA="${CORNER_QUERY_ALPHA:-0.001}"
SEED="${CORNER_SEED:-20250729}"

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

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/cornerrow-submit.XXXXXX")"
cleanup() {
    rm -rf -- "$WORKDIR"
}
trap cleanup EXIT

GROUPS="$WORKDIR/groups.txt"
FILTERS="$WORKDIR/filters.tsv"
DATASET_DIRS="$WORKDIR/dataset_dirs.txt"

DRY_RUN=0
FORCE=0
DO_BUILD=1

# ---------------------------------------------------------------------------
# Argument parsing. Flags first, then leaf-count groups plus optional dataset
# IDs/ranges associated with the most recently named group. Commas may be used
# as separators as well.
# ---------------------------------------------------------------------------
: > "$GROUPS"
: > "$FILTERS"
current_group=""
for argument in "$@"; do
    case "$argument" in
        --dry-run)  DRY_RUN=1;   continue ;;
        --force)    FORCE=1;     continue ;;
        --no-build) DO_BUILD=0;  continue ;;
        --*)
            echo "ERROR: unknown option '$argument'" >&2
            exit 1
            ;;
    esac

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

# ---------------------------------------------------------------------------
# Environment. Adjust these to whatever the cluster actually needs; the old
# /opt/homebrew path was macOS-only and will silently do nothing here.
#
# sbatch exports the submitting environment to jobs by default, so PATH,
# JULIA_BIN, and R_HOME set here are visible inside each job. If your site sets
# --export=NONE in a system-wide plugstack, set them in the sbatch script too.
# ---------------------------------------------------------------------------
# module load julia cmake gcc R    # <- uncomment/edit for your cluster
export PATH="$HOME/.juliaup/bin:$PATH"
if command -v R >/dev/null 2>&1; then
    R_HOME="$(R RHOME)"
    export R_HOME
fi
export JULIA_BIN="${JULIA_BIN:-julia}"

for command_name in sbatch julia; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $command_name" >&2
        exit 1
    fi
done
if [[ "$DO_BUILD" -eq 1 ]] && ! command -v cmake >/dev/null 2>&1; then
    echo "ERROR: required command not found: cmake" >&2
    exit 1
fi
if [[ ! -f "$SBATCH_SCRIPT" ]]; then
    echo "ERROR: missing job script: $SBATCH_SCRIPT" >&2
    exit 1
fi
if [[ ! -f "$COMPARE_JL" ]]; then
    echo "ERROR: missing comparator: $COMPARE_JL" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Build exactly once, here, before anything is queued. Jobs must never build:
# concurrent jobs would race on the same object files.
# ---------------------------------------------------------------------------
if [[ "$DO_BUILD" -eq 1 ]]; then
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
fi
if [[ ! -x "$BINARY" ]]; then
    echo "ERROR: no binary at $BINARY (drop --no-build, or build first)" >&2
    exit 1
fi

# Pin the binary for this run. Jobs may sit in the queue for hours; if you
# rebuild in the meantime, queued jobs would otherwise silently execute
# different code than the ones that already ran.
RUN_STAMP="$(date +%Y%m%dT%H%M%S)"
RUN_BINARY="$BIN_SNAPSHOT_ROOT/tree-qmc-$RUN_STAMP"
mkdir -p "$BIN_SNAPSHOT_ROOT" "$LOG_ROOT"
cp "$BINARY" "$RUN_BINARY"
chmod +x "$RUN_BINARY"

# Provenance: what code, what parameters, when. Cheap now, invaluable later.
{
    echo "run_stamp:    $RUN_STAMP"
    echo "submitted_by: $(whoami)@$(hostname)"
    echo "binary:       $RUN_BINARY"
    echo "corner_k:     $K"
    echo "heavy:        $HEAVY"
    echo "corner_tau:   $TAU"
    echo "query_alpha:  $QUERY_ALPHA"
    echo "corner_seed:  $SEED"
    echo "git_commit:   $(git -C "$SOURCE_DIR" rev-parse HEAD 2>/dev/null || echo unavailable)"
    echo "git_dirty:    $(git -C "$SOURCE_DIR" status --porcelain 2>/dev/null | wc -l) modified files"
    echo "groups:       $(tr '\n' ' ' < "$GROUPS")"
} > "$BIN_SNAPSHOT_ROOT/run-$RUN_STAMP.info"

echo "Parameters: k=$K heavy=$HEAVY tau=$TAU query-alpha=$QUERY_ALPHA seed=$SEED"

# ---------------------------------------------------------------------------
# Enumerate datasets and submit.
# ---------------------------------------------------------------------------
submitted=0
skipped=0
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

        # Resumability: a completed dataset already has its shard. Re-running
        # the submitter after a partial failure fills only the gaps.
        if [[ "$FORCE" -eq 0 && -s "$output_dir/cornerrow.csv" ]]; then
            skipped="$((skipped + 1))"
            continue
        fi

        mkdir -p "$output_dir" "$LOG_ROOT/$group"
        job_log="$LOG_ROOT/$group/cornerrow_${network_id}_%j.log"

        if [[ "$DRY_RUN" -eq 1 ]]; then
            echo "[dry-run] $dataset ($(basename "$gene_tree_path"))"
        else
            sbatch --job-name="cornerrow_${group}_${network_id}" \
                --output="$job_log" \
                --error="$job_log" \
                "$SBATCH_SCRIPT" \
                "$group" "$network_id" "$gene_tree_path" "$refinement" \
                "$true_network" "$output_dir" "$RUN_BINARY" "$COMPARE_JL" \
                "$K" "$HEAVY" "$TAU" "$QUERY_ALPHA" "$SEED"
        fi
        submitted="$((submitted + 1))"
    done < "$DATASET_DIRS"

    if [[ "$found_in_group" -eq 0 ]]; then
        echo "ERROR: no refinement datasets found under $refinement_group" >&2
        exit 1
    fi
done < "$GROUPS"

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "Dry run: $submitted job(s) would be submitted, $skipped already complete."
else
    echo "Submitted $submitted job(s); skipped $skipped already complete."
    echo "Watch:    squeue -u \"$(whoami)\""
    echo "History:  sacct -X --starttime=today --format=JobID,JobName%30,State,Elapsed,MaxRSS,ExitCode"
    echo "Collect:  $HERE/collect_cornerrow.sh"
fi
