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
#   $1  METHOD        rowsweep | cornerrow | branchcut
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
#   $13 PYTHON        interpreter that can import treeswift, as an absolute path
#   $14 MODULES       space-separated modulefiles to load first, "" for none
#
# Peak memory is captured with GNU time when it is available (set GNU_TIME to
# point at it, e.g. gtime on macOS); it is reported as NA otherwise.
#
# Resource notes (verified against the CBCB cluster, July 2026):
#   --constraint=EPYC-7313 pins the CPU model so runtimes are comparable across
#   datasets; --exclusive keeps co-tenants from contending for cache and memory
#   bandwidth; --mem-bind=local keeps allocations NUMA-local.
#
#   This is a whole reserved node per dataset. It buys comparable wall clocks at
#   the price of throughput, and a sweep of a few hundred datasets will sit in
#   the queue for a long time -- and be rejected outright if the account is not
#   allowed the qos or the constraint matches no node. Override per submission
#   with submit_tob.sh --sbatch-args, e.g.
#       --sbatch-args '--mem=16G --qos=normal --constraint= --oversubscribe'
#   Options given to sbatch on the command line beat these directives.
#
set -euo pipefail

if [[ $# -ne 14 ]]; then
    sed -n '13,34p' "$0" >&2
    exit 1
fi

METHOD="$1"; METHOD_ARGS="$2"; GROUP="$3"; NETWORK_ID="$4"
GENE_TREES="$5"; REFINEMENT="$6"; TRUE_NETWORK="$7"; OUT_DIR="$8"
TRUE_TOB_DIR="$9"; BINARY="${10}"; TOB_JL="${11}"; COMPARE_PY="${12}"; PYTHON="${13}"
MODULES="${14}"

case "$METHOD" in
    rowsweep)  METHOD_FLAG="--rowsweep-blob"  ;;
    cornerrow) METHOD_FLAG="--cornerrow-blob" ;;
    branchcut) METHOD_FLAG="--branchcut-blob" ;;
    *) echo "ERROR: unknown method: $METHOD (expected rowsweep, cornerrow or branchcut)" >&2; exit 1 ;;
esac

export PATH="$HOME/.juliaup/bin:/opt/homebrew/bin:$PATH"

# --- modules ----------------------------------------------------------------
# submit_tob.sh --modules, forwarded so the common case needs no site_env.sh.
# `module` is a shell function and this shell is not interactive, so its init
# script has to be sourced first; those scripts are not written to survive
# `set -u`, hence the fence. Kept in step with the same block in submit_tob.sh,
# which loads these before checking the interpreter -- the point of loading them
# in both places is that what was checked there is what runs here.
if [[ -n "$MODULES" ]]; then
    set +u
    if ! type module >/dev/null 2>&1; then
        for init in /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash \
                    /etc/profile.d/modules.sh; do
            # shellcheck disable=SC1090
            [[ -f "$init" ]] && { source "$init"; break; }
        done
    fi
    if ! type module >/dev/null 2>&1; then
        set -u
        echo "ERROR: modules requested ($MODULES) but there is no module command" >&2
        echo "       on $(hostname)." >&2
        exit 1
    fi
    for m in $MODULES; do
        # A module that did not load must stop the job. Carrying on runs the
        # dataset under whatever was on PATH anyway and writes a shard that
        # looks exactly like the others.
        if ! module load "$m"; then
            set -u
            echo "ERROR: module load $m failed on $(hostname)" >&2
            exit 1
        fi
    done
    set -u
fi

# Whatever a particular machine needs to make julia and python findable --
# module loads, a different PATH -- belongs to that machine and not to this
# script, which is why it lives in an optional file beside it rather than in the
# lines above. Create experiments/site_env.sh on the cluster (see
# site_env.sh.example); on a laptop it is simply absent.
#
# The path is derived from TOB_JL rather than from BASH_SOURCE: sbatch copies
# this script into the node's spool directory before running it, so
# BASH_SOURCE[0] here is /var/spool/.../slurm_script, not the checkout.
SITE_ENV="$(dirname "$TOB_JL")/site_env.sh"
if [[ -f "$SITE_ENV" ]]; then
    echo "sourcing $SITE_ENV"
    # shellcheck disable=SC1090
    source "$SITE_ENV"
fi

command -v R >/dev/null 2>&1 && export R_HOME="$(R RHOME)"

for f in "$GENE_TREES" "$REFINEMENT" "$TRUE_NETWORK" "$BINARY" "$TOB_JL" "$COMPARE_PY"; do
    [[ -e "$f" ]] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done

# The interpreter is only needed at the end, but it is checked here, because the
# end is hours of tree-qmc away. This used to surface as a traceback after the
# estimate had already been computed and thrown away with the job.
#
# submit_tob.sh ran this same check before submitting. Reaching it here means
# the two machines disagree -- an interpreter under a prefix the compute nodes
# do not mount, or a module that was loaded on the login node and not passed to
# --modules -- so the message says that rather than repeating the install line.
if ! python_error="$("$PYTHON" -c 'import treeswift' 2>&1)"; then
    echo "ERROR: $PYTHON cannot import treeswift on $(hostname):" >&2
    printf '%s\n' "$python_error" | sed 's/^/         /' >&2
    echo "       This interpreter passed the same check on the submitting node," >&2
    echo "       so the two machines differ. Usual causes: the path is on storage" >&2
    echo "       local to the login node, or it needs a module this job did not" >&2
    echo "       load -- pass it with submit_tob.sh --modules." >&2
    exit 1
fi

mkdir -p "$OUT_DIR" "$TRUE_TOB_DIR/$GROUP/$NETWORK_ID"
ESTIMATE="$OUT_DIR/estimate.nwk"
LOG="$OUT_DIR/run.log"
SHARD="$OUT_DIR/row.csv"
TRUE_TOB="$TRUE_TOB_DIR/$GROUP/$NETWORK_ID/true_tob.nwk"
# A rerun is a fresh attempt, so the log describes that attempt and not the
# concatenation of every attempt so far.
: > "$LOG"

# The ground truth depends only on the dataset, so it is computed once and
# shared by every configuration.
if [[ ! -s "$TRUE_TOB" ]]; then
    command -v julia >/dev/null 2>&1 || {
        echo "ERROR: julia is not on PATH on $(hostname); it is needed to build" >&2
        echo "       the true tree of blobs. Load the module in this script or" >&2
        echo "       precompute $TRUE_TOB." >&2
        exit 1
    }
    echo "Computing the true tree of blobs for $GROUP/$NETWORK_ID"
    if ! julia --startup-file=no "$TOB_JL" "$TRUE_NETWORK" "$TRUE_TOB" >> "$LOG" 2>&1; then
        tail -20 "$LOG" >&2
        exit 1
    fi
fi
[[ -s "$TRUE_TOB" ]] || { echo "ERROR: no true tree of blobs for $GROUP/$NETWORK_ID" >&2; exit 1; }

echo "$METHOD $GROUP/$NETWORK_ID: $METHOD_ARGS"
USAGE_LOG="$OUT_DIR/usage.log"
TIMER="${GNU_TIME:-/usr/bin/time}"
TIMER_ARGS=()
if [[ -x "$TIMER" ]] && "$TIMER" -v true >/dev/null 2>&1; then
    TIMER_ARGS=("$TIMER" -v -o "$USAGE_LOG")
fi

START="$("$PYTHON" -c 'import time; print(time.monotonic())')"
# METHOD_ARGS is deliberately word-split: it carries several flags.
# shellcheck disable=SC2086
if ! ${TIMER_ARGS[@]+"${TIMER_ARGS[@]}"} "$BINARY" -i "$GENE_TREES" \
        --blobsearchonly "$REFINEMENT" \
        --blob $METHOD_FLAG $METHOD_ARGS \
        --override -o "$ESTIMATE" >> "$LOG" 2>&1; then
    tail -20 "$LOG" >&2
    exit 1
fi
END="$("$PYTHON" -c 'import time; print(time.monotonic())')"
[[ -s "$ESTIMATE" ]] || { echo "ERROR: $BINARY exited 0 but wrote no tree to $ESTIMATE" >&2; exit 1; }
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
