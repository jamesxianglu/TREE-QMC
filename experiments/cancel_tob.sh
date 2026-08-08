#!/usr/bin/env bash
#
# Cancel the jobs of one run, and optionally clear its outputs so the run can be
# resubmitted from scratch.
#
# Only the job ids that this run recorded in jobids.txt are cancelled. A run
# directory is the only thing that knows which jobs belong to it, and
# "scancel -u $USER" would take out unrelated work alongside them.
#
# Usage:
#   ./cancel_tob.sh --run-dir DIR [options]
#
# Required:
#   --run-dir DIR   a run directory written by submit_tob.sh, or just its label,
#                   which is resolved against results/runs
#
# Options:
#   --reset         after cancelling, delete the run directory so a resubmission
#                   redoes every dataset. The cached true trees of blobs live
#                   outside the run directory, under results/true_tob, and are
#                   deliberately kept: they depend on the dataset alone, not on
#                   the configuration, and recomputing them costs a julia run
#                   per dataset for no change in the answer.
#   --yes           do not ask before deleting
#   --wait SECONDS  how long to wait for the queue to drain (default 60)
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$HERE")")"
OUT_ROOT="$PROJECT_ROOT/results/runs"

usage() { sed -n '3,27p' "$0"; exit "${1:-1}"; }
[[ $# -eq 0 ]] && usage

RUN_DIR=""; RESET=0; ASSUME_YES=0; WAIT_FOR=60
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-dir) RUN_DIR="${2:-}";  shift 2 ;;
        --reset)   RESET=1;           shift ;;
        --yes)     ASSUME_YES=1;      shift ;;
        --wait)    WAIT_FOR="${2:-}"; shift 2 ;;
        -h|--help) usage 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage ;;
    esac
done
[[ -n "$RUN_DIR" ]] || { echo "ERROR: --run-dir is required" >&2; exit 1; }

# Run labels are long enough that pasting the full path is where this goes
# wrong, so a bare label is resolved against the default out-root, and a path
# that does not resolve is answered with the list of labels that do.
if [[ ! -d "$RUN_DIR" && -d "$OUT_ROOT/$RUN_DIR" ]]; then
    RUN_DIR="$OUT_ROOT/$RUN_DIR"
fi
if [[ ! -d "$RUN_DIR" ]]; then
    echo "ERROR: no such run directory: $RUN_DIR" >&2
    parent="$(dirname "$RUN_DIR")"
    [[ -d "$parent" ]] || parent="$OUT_ROOT"
    if [[ -d "$parent" ]]; then
        echo >&2
        echo "Run directories under $parent:" >&2
        find "$parent" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; |
            sort | sed 's/^/  /' >&2
        echo >&2
        echo "--run-dir also accepts a bare label from $OUT_ROOT." >&2
    else
        echo "  ($parent does not exist either)" >&2
    fi
    exit 1
fi
# A run directory always has a run.meta. Requiring it here is what keeps --reset
# from being pointed at an arbitrary path.
[[ -f "$RUN_DIR/run.meta" ]] ||
    { echo "ERROR: $RUN_DIR has no run.meta; refusing to touch it" >&2; exit 1; }

# --- cancel ----------------------------------------------------------------
if [[ -s "$RUN_DIR/jobids.txt" ]]; then
    # submit_tob.sh records "jobid<TAB>group/rep"; --parsable emits "jobid" or,
    # under multi-cluster, "jobid;cluster".
    # Not mapfile: macOS still ships bash 3.2, and this script is worth being
    # able to run on the laptop that edits it.
    ids=()
    while IFS= read -r id; do ids+=("$id"); done < <(
        cut -f1 "$RUN_DIR/jobids.txt" | cut -d';' -f1 | grep -E '^[0-9]+$')
    echo "${#ids[@]} job id(s) recorded for this run"
    if command -v squeue >/dev/null 2>&1; then
        live="$(squeue -h -j "$(IFS=,; echo "${ids[*]}")" -o '%i' 2>/dev/null | wc -l | tr -d ' ')"
        echo "$live still in the queue"
    fi
    if command -v scancel >/dev/null 2>&1 && [[ "${#ids[@]}" -gt 0 ]]; then
        scancel "${ids[@]}" 2>/dev/null || true
        echo "scancel sent"
        waited=0
        while [[ "$waited" -lt "$WAIT_FOR" ]]; do
            left="$(squeue -h -j "$(IFS=,; echo "${ids[*]}")" -o '%i' 2>/dev/null | wc -l | tr -d ' ')"
            [[ "$left" -eq 0 ]] && break
            sleep 5
            waited=$((waited + 5))
        done
        left="$(squeue -h -j "$(IFS=,; echo "${ids[*]}")" -o '%i' 2>/dev/null | wc -l | tr -d ' ')"
        if [[ "$left" -eq 0 ]]; then
            echo "queue is clear"
        else
            echo "WARNING: $left job(s) still listed after ${WAIT_FOR}s." >&2
            echo "         Cancelling jobs is not instant; check squeue before resubmitting," >&2
            echo "         or a survivor will write into the new run directory." >&2
        fi
    fi
else
    echo "no jobids.txt in $RUN_DIR; nothing to cancel"
fi

# --- reset -----------------------------------------------------------------
[[ "$RESET" -eq 1 ]] || { echo "Left $RUN_DIR in place (pass --reset to clear it)."; exit 0; }

shards="$(find "$RUN_DIR" -mindepth 3 -maxdepth 3 -name row.csv | wc -l | tr -d ' ')"
echo
echo "About to delete $RUN_DIR"
echo "  $shards completed shard(s) will be lost"
echo "  results/true_tob is outside this directory and is kept"
if [[ "$ASSUME_YES" -eq 0 ]]; then
    read -r -p "Type the run label to confirm: " typed
    [[ "$typed" == "$(basename "$RUN_DIR")" ]] || { echo "Not confirmed; nothing deleted." >&2; exit 1; }
fi
rm -rf "$RUN_DIR"
echo "deleted $RUN_DIR"
echo "Resubmit with the same submit_tob.sh command."
