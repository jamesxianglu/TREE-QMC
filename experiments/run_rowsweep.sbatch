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
# Evaluate row-sweep on a SINGLE dataset and write a per-dataset CSV shard.
#
# This script is not meant to be run by hand; submit_rowsweep.sh invokes it via
# sbatch, once per dataset. Results are collected afterwards by
# collect_rowsweep.sh.
#
# Positional arguments (all required, all absolute paths where applicable):
#   $1  GROUP         leaf-count group, e.g. n15
#   $2  NETWORK_ID    dataset id within the group, e.g. 07
#   $3  GENE_TREES    gene tree file (g_500.nwk or iqtree_500.nwk)
#   $4  REFINEMENT    astral4-rooted.tre for this dataset
#   $5  TRUE_NETWORK  true_net.nwk for this dataset
#   $6  OUT_DIR       per-dataset output directory
#   $7  BINARY        tree-qmc binary (a run-pinned snapshot, not the build dir)
#   $8  COMPARE_JL    path to compare_inferred_tob.jl
#   $9  DELTA
#   $10 QUERY_ALPHA
#
# Resource notes (verified against the CBCB cluster, July 2026):
#   --constraint=EPYC-7313 pins the CPU model so runtimes are comparable across
#   datasets. It resolves to cbcb[00-20]: 21 nodes, 32 CPUs and ~2 TB each.
#   --exclusive prevents co-tenant jobs from contending for cache and memory
#   bandwidth. --mem-bind=local keeps allocations NUMA-local. Together these are
#   what make the reported wall clock a property of the algorithm rather than of
#   whatever else happened to be running.
#
#   --qos=highmem is REQUIRED, not redundant. --exclusive allocates all 32 CPUs
#   of the node, but the default QOS caps a job at cpu=4,mem=32G, so without
#   this line the job pends forever on QOSMaxCpuPerJobLimit instead of ever
#   being scheduled. highmem allows cpu=128,mem=2T with a 21-day walltime.
#
#   --mem=256G is a cgroup ceiling, not a reservation; --exclusive already
#   reserves the whole node. Raise it (up to ~2T) if row-sweep is ever OOM-
#   killed on the larger leaf counts. Re-verify any of this with:
#     sinfo -p cbcb -o "%20N %10c %10m %50f" | grep 7313
#     sacctmgr show qos format=Name,MaxTRESPerJob%30,MaxWall
#     sacctmgr show assoc user=$USER format=Account,User,QOS%40

# Deliberately NOT -e: we want to record failures rather than die silently.
set -uo pipefail

if [[ $# -ne 10 ]]; then
    echo "ERROR: expected 10 arguments, got $#" >&2
    exit 2
fi

GROUP="$1"
NETWORK_ID="$2"
GENE_TREES="$3"
REFINEMENT="$4"
TRUE_NETWORK="$5"
OUT_DIR="$6"
BINARY="$7"
COMPARE_JL="$8"
DELTA="$9"
QUERY_ALPHA="${10}"

METHOD="ROWSWEEP"
DATASET="$GROUP/$NETWORK_ID"
TIMER="/usr/bin/time"
JULIA_BIN="${JULIA_BIN:-julia}"

CSV_HEADER="taxa,network_id,true_tob,estimated_tob,method,delta,query_alpha,fn,fp,rf,wall_clock_seconds,max_rss_kb,node"

INFERRED_TOB="$OUT_DIR/rowsweep.nwk"
TRUE_TOB="$OUT_DIR/true_tob.nwk"
CONSTRUCTOR_LOG="$OUT_DIR/rowsweep.log"
USAGE_LOG="$OUT_DIR/rowsweep.usage.log"
MANIFEST="$OUT_DIR/manifest.tsv"
COMPARISONS="$OUT_DIR/comparisons.tsv"
CSV="$OUT_DIR/rowsweep.csv"
FAILED_MARKER="$OUT_DIR/rowsweep.failed"

mkdir -p "$OUT_DIR"
# Clear stale state so a re-run never leaves a half-old result behind.
rm -f "$CSV" "$FAILED_MARKER"

fail() {
    printf '%s\n' "$*" > "$FAILED_MARKER"
    echo "ERROR [$DATASET]: $*" >&2
    exit 1
}

# ---------------------------------------------------------------------------
# Preflight. Cheap, and turns a confusing 3-hour-later failure into a fast one.
# ---------------------------------------------------------------------------
[[ -x "$BINARY" ]]        || fail "binary is not executable: $BINARY"
[[ -f "$GENE_TREES" ]]    || fail "missing gene trees: $GENE_TREES"
[[ -f "$REFINEMENT" ]]    || fail "missing refinement: $REFINEMENT"
[[ -f "$TRUE_NETWORK" ]]  || fail "missing true network: $TRUE_NETWORK"
[[ -f "$COMPARE_JL" ]]    || fail "missing comparator: $COMPARE_JL"
[[ -x "$TIMER" ]]         || fail "GNU time not found at $TIMER (the bash builtin does not support -v)"
command -v "$JULIA_BIN" >/dev/null 2>&1 || fail "julia not on PATH inside the job (set JULIA_BIN or fix the module environment)"

echo "=== $DATASET on $(hostname) ==="
echo "cpu:        $(awk -F': ' '/model name/{print $2; exit}' /proc/cpuinfo 2>/dev/null)"
echo "gene trees: $(basename "$GENE_TREES")"
echo "delta:      $DELTA"
echo "alpha:      $QUERY_ALPHA"

# ---------------------------------------------------------------------------
# Run the algorithm under GNU time. This is the only timed region.
# ---------------------------------------------------------------------------
"$TIMER" -v -o "$USAGE_LOG" \
    "$BINARY" -i "$GENE_TREES" \
        --blobsearchonly "$REFINEMENT" \
        --blob --rowsweep-blob \
        --delta "$DELTA" --query-alpha "$QUERY_ALPHA" \
        --override -o "$INFERRED_TOB" \
    > "$CONSTRUCTOR_LOG" 2>&1
run_status=$?

if [[ "$run_status" -ne 0 ]]; then
    cat "$CONSTRUCTOR_LOG" >&2
    fail "tree-qmc exited with status $run_status"
fi
[[ -s "$INFERRED_TOB" ]] || fail "tree-qmc produced no output at $INFERRED_TOB"

# ---------------------------------------------------------------------------
# Parse GNU time's report.
#
# The elapsed line is either h:mm:ss.ss or m:ss.ss, so it has to be split on
# colons and reassembled rather than read as a float.
# ---------------------------------------------------------------------------
elapsed_field="$(awk -F': ' '/Elapsed \(wall clock\) time/{print $NF}' "$USAGE_LOG")"
[[ -n "$elapsed_field" ]] || fail "no elapsed time line in $USAGE_LOG"

IFS=: read -r part1 part2 part3 <<< "$elapsed_field"
if [[ -n "${part3:-}" ]]; then
    WALL_CLOCK="$(awk -v h="$part1" -v m="$part2" -v s="$part3" \
        'BEGIN { printf "%.3f", h*3600 + m*60 + s }')"
else
    WALL_CLOCK="$(awk -v m="$part1" -v s="$part2" \
        'BEGIN { printf "%.3f", m*60 + s }')"
fi

MAX_RSS_KB="$(awk -F': ' '/Maximum resident set size/{print $NF}' "$USAGE_LOG")"
MAX_RSS_KB="${MAX_RSS_KB:-NA}"

echo "  wall clock: ${WALL_CLOCK}s"
echo "  max RSS:    ${MAX_RSS_KB} KB"

# ---------------------------------------------------------------------------
# Compare against the true network. The one-row manifest preserves the
# comparator's existing machine-readable batch interface.
# ---------------------------------------------------------------------------
printf 'dataset\tdelta\tquery_alpha\ttrue_network\tinferred_tree\ttrue_tob_output\n' \
    > "$MANIFEST"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$DATASET" "$DELTA" "$QUERY_ALPHA" "$TRUE_NETWORK" \
    "$INFERRED_TOB" "$TRUE_TOB" >> "$MANIFEST"

echo "  computing FN, FP, and normalized RF..."
if ! "$JULIA_BIN" --startup-file=no "$COMPARE_JL" --batch "$MANIFEST" "$COMPARISONS"; then
    fail "comparator failed"
fi

exec 3< "$COMPARISONS" || fail "cannot read $COMPARISONS"
if ! IFS= read -r comparison_header <&3 || \
   ! IFS=$'\t' read -r comparison_dataset _ _ fn fp rf <&3; then
    fail "missing comparison row"
fi
if [[ "$comparison_header" != $'dataset\tdelta\tquery_alpha\tfn\tfp\tnormalized_rf' ]]; then
    fail "comparison output has an unexpected header: $comparison_header"
fi
if IFS= read -r _extra <&3; then
    fail "expected exactly one comparison row"
fi
exec 3<&-
[[ "$comparison_dataset" == "$DATASET" ]] || \
    fail "comparison dataset mismatch: $comparison_dataset != $DATASET"

# ---------------------------------------------------------------------------
# Write this dataset's shard. One file per job: no concurrent appends, so
# nothing can interleave on the shared filesystem.
# ---------------------------------------------------------------------------
{
    printf '%s\n' "$CSV_HEADER"
    printf '%s,%s,"%s","%s",%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$GROUP" "$NETWORK_ID" "$TRUE_TOB" "$INFERRED_TOB" "$METHOD" \
        "$DELTA" "$QUERY_ALPHA" "$fn" "$fp" "$rf" \
        "$WALL_CLOCK" "$MAX_RSS_KB" "$(hostname)"
} > "$CSV"

echo "  wrote $CSV"
