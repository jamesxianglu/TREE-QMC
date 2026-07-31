#!/usr/bin/env bash
#
# Run one tree-of-blobs experiment: one method, one configuration, over a chosen
# set of datasets, either on SLURM or locally. Every experiment lives in its own
# run directory named after the method and its parameters, and each dataset
# leaves a one-row CSV shard there; collect_tob.sh merges them.
#
# The method, its flags, the datasets and the runner have to be typed out: a run
# is only reproducible if its configuration was written down. The gene tree file
# is the one exception, since it is a property of the dataset rather than of the
# experiment; whichever file is used is recorded in every row of the output.
#
# Usage:
#   ./submit_tob.sh --method METHOD --args 'FLAGS' --datasets 'SPEC' \
#                   --runner (local|slurm) [options]
#
# Required:
#   --method METHOD      rowsweep | cornerrow
#   --args 'FLAGS'       tree-qmc flags for this configuration, quoted as one
#                        string, e.g.
#                          '--oracle t1 --delta 0.25 --query-alpha 0.001'
#                          '--oracle t1+cf|maj --rowsweep-tau 0.45,0.60 \
#                           --rowsweep-heavy 1,2 --query-alpha 0.001'
#
# Options:
#   --gene-trees FILE    basename inside each dataset directory. Defaults to
#                        iqtree_500.nwk, falling back to g_500.nwk per dataset
#                        when the first is absent. Which file was actually used
#                        changes the result, so the resolved name is recorded in
#                        every CSV row rather than left to be inferred.
#   --datasets 'SPEC'    space-separated groups, each optionally narrowed by
#                        ':ids', e.g. 'n15 n25 n50:00-09 n100:07'
#   --runner RUNNER      local | slurm
#
#   --label NAME         run directory name; defaults to a slug built from the
#                        method and flags
#   --out-root DIR       where run directories live (default results/runs)
#   --dry-run            print what would happen and stop
#   --force              redo datasets that already have a shard
#   --no-build           reuse the existing binary instead of rebuilding
#
# Example:
#   ./submit_tob.sh --method rowsweep --runner slurm --gene-trees g_500.nwk \
#       --datasets 'n15 n25 n50 n100 n150 n200' \
#       --args '--oracle t1+cf|maj --rowsweep-tau 0.45,0.60 --rowsweep-heavy 1,2 --query-alpha 0.001'
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"
PROJECT_ROOT="$(dirname "$SOURCE_DIR")"
DATA_ROOT="$PROJECT_ROOT/data/camus-dataset"
REFINEMENT_ROOT="$PROJECT_ROOT/data/refinements"
BUILD_DIR="$SOURCE_DIR/build"
TRUE_TOB_DIR="$PROJECT_ROOT/results/true_tob"
RUNNER_SCRIPT="$HERE/run_tob.sbatch"
TOB_JL="$HERE/compute_tree_of_blob.jl"
COMPARE_PY="$PROJECT_ROOT/results/analysis/compare_two_trees.py"

usage() { sed -n '3,45p' "$0"; exit "${1:-1}"; }
[[ $# -eq 0 ]] && usage

METHOD=""; METHOD_ARGS=""; GENE_TREES=""; DATASETS=""; RUNNER=""
LABEL=""; OUT_ROOT="$PROJECT_ROOT/results/runs"
DRY_RUN=0; FORCE=0; BUILD=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --method)      METHOD="${2:-}";      shift 2 ;;
        --args)        METHOD_ARGS="${2:-}"; shift 2 ;;
        --gene-trees)  GENE_TREES="${2:-}";  shift 2 ;;
        --datasets)    DATASETS="${2:-}";    shift 2 ;;
        --runner)      RUNNER="${2:-}";      shift 2 ;;
        --label)       LABEL="${2:-}";       shift 2 ;;
        --out-root)    OUT_ROOT="${2:-}";    shift 2 ;;
        --dry-run)     DRY_RUN=1;            shift ;;
        --force)       FORCE=1;              shift ;;
        --no-build)    BUILD=0;              shift ;;
        -h|--help)     usage 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage ;;
    esac
done

missing=""
[[ -n "$METHOD" ]]      || missing="$missing --method"
[[ -n "$METHOD_ARGS" ]] || missing="$missing --args"
[[ -n "$DATASETS" ]]    || missing="$missing --datasets"
[[ -n "$RUNNER" ]]      || missing="$missing --runner"
if [[ -n "$missing" ]]; then
    echo "ERROR: required and not supplied:$missing" >&2
    exit 1
fi
case "$METHOD" in rowsweep|cornerrow) ;; *)
    echo "ERROR: --method must be rowsweep or cornerrow, got '$METHOD'" >&2; exit 1 ;;
esac
case "$RUNNER" in local|slurm) ;; *)
    echo "ERROR: --runner must be local or slurm, got '$RUNNER'" >&2; exit 1 ;;
esac
if [[ "$RUNNER" == slurm ]] && ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: --runner slurm but sbatch is not on PATH" >&2
    exit 1
fi

# --- run identity ----------------------------------------------------------
# The slug has to survive being a directory name and still say what was run, so
# flag punctuation collapses to dashes: "--oracle t1+cf|maj --rowsweep-tau
# 0.45,0.60" becomes "oracle-t1+cf-maj-rowsweep-tau-0.45-0.60".
slugify() {
    printf '%s' "$1" |
        sed -e 's/--//g' -e 's/[^A-Za-z0-9._+-]\{1,\}/-/g' \
            -e 's/-\{2,\}/-/g' -e 's/^-//' -e 's/-$//'
}
if [[ -z "$LABEL" ]]; then
    LABEL="${METHOD}__$(slugify "$METHOD_ARGS")"
fi
RUN_DIR="$OUT_ROOT/$LABEL"

# --- dataset expansion -----------------------------------------------------
WORK="$(mktemp)"
trap 'rm -f "$WORK"' EXIT
for token in $DATASETS; do
    group="${token%%:*}"
    [[ "$group" =~ ^n[0-9]+$ ]] || { echo "ERROR: bad group in --datasets: $token" >&2; exit 1; }
    [[ -d "$REFINEMENT_ROOT/$group" ]] || { echo "ERROR: no refinements for $group" >&2; exit 1; }
    if [[ "$token" == *:* ]]; then
        selector="${token#*:}"
        for part in ${selector//,/ }; do
            start="${part%%-*}"; end="${part#*-}"; [[ "$part" == *-* ]] || end="$start"
            [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ ]] ||
                { echo "ERROR: bad id selector: $part" >&2; exit 1; }
            for ((i = 10#$start; i <= 10#$end; i++)); do
                printf '%s\t%0*d\n' "$group" "${#start}" "$i" >> "$WORK"
            done
        done
    else
        find "$REFINEMENT_ROOT/$group" -mindepth 1 -maxdepth 1 -type d \
            -exec basename {} \; | sort | sed "s|^|$group\t|" >> "$WORK"
    fi
done
[[ -s "$WORK" ]] || { echo "ERROR: --datasets selected nothing" >&2; exit 1; }

echo "method      $METHOD"
echo "args        $METHOD_ARGS"
echo "gene trees  ${GENE_TREES:-iqtree_500.nwk, else g_500.nwk}"
echo "runner      $RUNNER"
echo "run dir     $RUN_DIR"
echo "datasets    $(wc -l < "$WORK" | tr -d ' ')"

# --- build -----------------------------------------------------------------
BINARY="$RUN_DIR/tree-qmc"
if [[ "$DRY_RUN" -eq 0 ]]; then
    mkdir -p "$RUN_DIR"
    if [[ "$BUILD" -eq 1 ]]; then
        echo "Building TREE-QMC..."
        cmake -S "$SOURCE_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release \
            > "$RUN_DIR/configure.log" 2>&1
        cmake --build "$BUILD_DIR" -j4 > "$RUN_DIR/build.log" 2>&1
    fi
    [[ -x "$BUILD_DIR/tree-qmc" ]] || { echo "ERROR: no binary at $BUILD_DIR/tree-qmc" >&2; exit 1; }
    # Pin the binary to the run so a later rebuild cannot change what these
    # numbers mean.
    cp -f "$BUILD_DIR/tree-qmc" "$BINARY"
    {
        echo "method=$METHOD"
        echo "args=$METHOD_ARGS"
        echo "gene_trees=${GENE_TREES:-auto}"
        echo "runner=$RUNNER"
        echo "label=$LABEL"
        echo "submitted=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } > "$RUN_DIR/run.meta"
fi

PYTHON="$PROJECT_ROOT/.venv/bin/python"
[[ -x "$PYTHON" ]] || PYTHON="python3"

submitted=0
skipped=0
while IFS=$'\t' read -r group rep; do
    dataset_dir="$DATA_ROOT/$group/$rep"
    refinement="$REFINEMENT_ROOT/$group/$rep/astral4-rooted.tre"
    true_network="$dataset_dir/true_net.nwk"
    out_dir="$RUN_DIR/$group/$rep"

    # Resolve the gene trees per dataset: an explicit --gene-trees must exist,
    # otherwise prefer the estimated trees and fall back to g_500.nwk.
    if [[ -n "$GENE_TREES" ]]; then
        gene_trees="$dataset_dir/$GENE_TREES"
        [[ -f "$gene_trees" ]] || { echo "ERROR: no $GENE_TREES for $group/$rep" >&2; exit 1; }
    elif [[ -f "$dataset_dir/iqtree_500.nwk" ]]; then
        gene_trees="$dataset_dir/iqtree_500.nwk"
    elif [[ -f "$dataset_dir/g_500.nwk" ]]; then
        gene_trees="$dataset_dir/g_500.nwk"
    else
        echo "ERROR: neither iqtree_500.nwk nor g_500.nwk for $group/$rep" >&2
        exit 1
    fi

    for f in "$refinement" "$true_network"; do
        [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
    done
    if [[ -s "$out_dir/row.csv" && "$FORCE" -eq 0 ]]; then
        skipped=$((skipped + 1))
        continue
    fi

    args=("$METHOD" "$METHOD_ARGS" "$group" "$rep" "$gene_trees" "$refinement"
          "$true_network" "$out_dir" "$TRUE_TOB_DIR" "$BINARY" "$TOB_JL"
          "$COMPARE_PY" "$PYTHON")
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "would run $group/$rep"
        submitted=$((submitted + 1))
        continue
    fi
    mkdir -p "$out_dir"
    if [[ "$RUNNER" == slurm ]]; then
        sbatch --job-name="${METHOD}-${group}-${rep}" \
               --output="$out_dir/slurm-%j.out" \
               "$RUNNER_SCRIPT" "${args[@]}" > /dev/null
    else
        bash "$RUNNER_SCRIPT" "${args[@]}"
    fi
    submitted=$((submitted + 1))
done < "$WORK"

echo "$submitted dataset(s) $( [[ "$RUNNER" == slurm ]] && echo submitted || echo run ), $skipped already done"
[[ "$DRY_RUN" -eq 1 ]] && exit 0
echo "Collect with: $HERE/collect_tob.sh --run-dir $RUN_DIR"
