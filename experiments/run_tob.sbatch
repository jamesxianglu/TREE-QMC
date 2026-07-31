#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --qos=highmem
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --constraint=EPYC-7313
#SBATCH --exclusive
#SBATCH --mem-bind=local
#
# Reconstruct a tree of blobs for ONE dataset with ONE method configuration and
# write a one-row CSV shard. Method-agnostic: everything that distinguishes one
# experiment from another arrives in METHOD_ARGS.
#
# Not meant to be run by hand -- submit_tob.sh invokes it, via sbatch or
# directly for a local run. Shards are merged afterwards by collect_tob.sh.
#
# Positional arguments, all required:
#   $1  METHOD        rowsweep | cornerrow
#   $2  METHOD_ARGS   tree-qmc flags for this configuration, as one string
#   $3  GROUP         leaf-count group, e.g. n15
#   $4  NETWORK_ID    dataset id within the group, e.g. 07
#   $5  GENE_TREES    gene tree file for this dataset
#   $6  REFINEMENT    astral4-rooted.tre for this dataset
#   $7  TRUE_NETWORK  true_net.nwk for this dataset
#   $8  OUT_DIR       per-dataset output directory
#   $9  TRUE_TOB_DIR  shared cache of ground-truth trees of blobs
#   $10 BINARY        tree-qmc binary (a run-pinned snapshot, not the build dir)
#   $11 TOB_JL        compute_tree_of_blob.jl
#   $12 COMPARE_PY    compare_two_trees.py
#   $13 PYTHON        interpreter that can import treeswift
#
# Peak memory is captured with GNU time when it is available (set GNU_TIME to
# point at it, e.g. gtime on macOS); it is reported as NA otherwise.
#
# Resource notes (verified against the CBCB cluster, July 2026):
#   --constraint=EPYC-7313 pins the CPU model so runtimes are comparable across
#   datasets; --exclusive keeps co-tenants from contending for cache and memory
#   bandwidth; --mem-bind=local keeps allocations NUMA-local.
#
set -euo pipefail

if [[ $# -ne 13 ]]; then
    sed -n '13,33p' "$0" >&2
    exit 1
fi

METHOD="$1"; METHOD_ARGS="$2"; GROUP="$3"; NETWORK_ID="$4"
GENE_TREES="$5"; REFINEMENT="$6"; TRUE_NETWORK="$7"; OUT_DIR="$8"
TRUE_TOB_DIR="$9"; BINARY="${10}"; TOB_JL="${11}"; COMPARE_PY="${12}"; PYTHON="${13}"

case "$METHOD" in
    rowsweep)  METHOD_FLAG="--rowsweep-blob"  ;;
    cornerrow) METHOD_FLAG="--cornerrow-blob" ;;
    *) echo "ERROR: unknown method: $METHOD (expected rowsweep or cornerrow)" >&2; exit 1 ;;
esac

export PATH="$HOME/.juliaup/bin:/opt/homebrew/bin:$PATH"
command -v R >/dev/null 2>&1 && export R_HOME="$(R RHOME)"

for f in "$GENE_TREES" "$REFINEMENT" "$TRUE_NETWORK" "$BINARY" "$TOB_JL" "$COMPARE_PY"; do
    [[ -e "$f" ]] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done

mkdir -p "$OUT_DIR" "$TRUE_TOB_DIR/$GROUP/$NETWORK_ID"
ESTIMATE="$OUT_DIR/estimate.nwk"
LOG="$OUT_DIR/run.log"
SHARD="$OUT_DIR/row.csv"
TRUE_TOB="$TRUE_TOB_DIR/$GROUP/$NETWORK_ID/true_tob.nwk"

# The ground truth depends only on the dataset, so it is computed once and
# shared by every configuration.
if [[ ! -s "$TRUE_TOB" ]]; then
    echo "Computing the true tree of blobs for $GROUP/$NETWORK_ID"
    julia --startup-file=no "$TOB_JL" "$TRUE_NETWORK" "$TRUE_TOB" >> "$LOG" 2>&1
fi
[[ -s "$TRUE_TOB" ]] || { echo "ERROR: no true tree of blobs for $GROUP/$NETWORK_ID" >&2; exit 1; }

echo "$METHOD $GROUP/$NETWORK_ID: $METHOD_ARGS"
USAGE_LOG="$OUT_DIR/usage.log"
TIMER="${GNU_TIME:-/usr/bin/time}"
TIMER_ARGS=()
if [[ -x "$TIMER" ]] && "$TIMER" -v true >/dev/null 2>&1; then
    TIMER_ARGS=("$TIMER" -v -o "$USAGE_LOG")
fi

START="$(python3 -c 'import time; print(time.monotonic())')"
# METHOD_ARGS is deliberately word-split: it carries several flags.
# shellcheck disable=SC2086
if ! ${TIMER_ARGS[@]+"${TIMER_ARGS[@]}"} "$BINARY" -i "$GENE_TREES" \
        --blobsearchonly "$REFINEMENT" \
        --blob $METHOD_FLAG $METHOD_ARGS \
        --override -o "$ESTIMATE" >> "$LOG" 2>&1; then
    tail -20 "$LOG" >&2
    exit 1
fi
END="$(python3 -c 'import time; print(time.monotonic())')"
ELAPSED="$(awk -v s="$START" -v e="$END" 'BEGIN { printf "%.3f", e - s }')"
MAX_RSS_KB="NA"
if [[ -s "$USAGE_LOG" ]]; then
    MAX_RSS_KB="$(awk -F': ' '/Maximum resident set size/{print $NF}' "$USAGE_LOG")"
    MAX_RSS_KB="${MAX_RSS_KB:-NA}"
fi

# compare_two_trees.py prints n_true,n_esti,nl,i1,i2,fn,fp
STATS="$("$PYTHON" "$COMPARE_PY" -t1 "$TRUE_TOB" -t2 "$ESTIMATE")"

ORACLE="$(printf '%s\n' "$METHOD_ARGS" | sed -n 's/.*--oracle[= ]\{1,\}\([^ ]*\).*/\1/p')"
[[ -n "$ORACLE" ]] || ORACLE="t1"

"$PYTHON" - "$METHOD" "$ORACLE" "$METHOD_ARGS" "$GROUP" "$NETWORK_ID" \
    "$(basename "$GENE_TREES")" "$STATS" "$ELAPSED" "$MAX_RSS_KB" "$(hostname)" \
    > "$SHARD" <<'PYEOF'
import csv
import sys

(method, oracle, args, taxa, network_id, gene_trees, stats, seconds, max_rss_kb,
 host) = sys.argv[1:11]
n_true, n_esti, nl, i1, i2, fn, fp = (int(v) for v in stats.split(","))


def rate(num, den):
    # A star tree holds no internal edge, so it can neither miss nor invent one.
    if den == 0:
        return f"{0:.10f}" if num == 0 else "NA"
    return f"{num / den:.10f}"


fnr, fpr = rate(fn, i1), rate(fp, i2)
avg = "NA" if "NA" in (fnr, fpr) else f"{(float(fnr) + float(fpr)) / 2:.10f}"
csv.writer(sys.stdout, lineterminator="\n").writerow([
    method, oracle, args, taxa, network_id, gene_trees, nl, i1, i2, fn, fp,
    fnr, fpr, avg, fn + fp, rate(fn + fp, 2 * nl - 6), seconds, max_rss_kb, host])
PYEOF

echo "Wrote $SHARD"
