#!/usr/bin/env bash
#
# Run one branch-cut configuration over one group locally, into
# results/runs/<label>/<group>/<rep>/estimate.nwk so edge_metrics.py and
# final_compare.py can score it directly.
#
#   ./bc_local.sh <label> <group> [nproc] -- <tree-qmc flags...>
#
# Existing estimate.nwk files are kept, so a re-run only fills gaps.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"; ROOT="$(dirname "$SOURCE_DIR")"
BINARY="${BC_BINARY:-$SOURCE_DIR/build/tree-qmc}"

LABEL="${1:?label}"; GROUP="${2:?group}"; NPROC="${3:-7}"
shift 3
[[ "${1:-}" == "--" ]] && shift
FLAGS=("$@")

export R_HOME="${R_HOME:-$(R RHOME)}"
export PATH="$HOME/.juliaup/bin:$PATH"
DATA="$ROOT/data/camus-dataset/$GROUP"
REF="${BC_REF:-$ROOT/data/refinements/$GROUP}"
REF_FILE="${BC_REF_FILE:-astral4-rooted.tre}"
SCORE_ROOT="${BC_SCORE_OUT:-}"

run_one() {
    local rep="$1"
    local out="$ROOT/results/runs/$LABEL/$GROUP/$rep"; mkdir -p "$out"
    [[ -s "$out/estimate.nwk" ]] && return 0
    local gt="${BC_GENES:-}"
    if [[ -z "$gt" ]]; then
        gt="iqtree_500.nwk"; [[ -f "$DATA/$rep/$gt" ]] || gt="g_500.nwk"
    fi
    local extra=()
    if [[ -n "$SCORE_ROOT" ]]; then
        mkdir -p "$SCORE_ROOT/$GROUP"
        extra=(--branchcut-score-out "$SCORE_ROOT/$GROUP/${rep}.tsv")
        [[ -n "${BC_QUAD_OUT:-}" ]] && extra+=(--branchcut-quad-out "$SCORE_ROOT/$GROUP/${rep}.quad.tsv")
    fi
    local t0=$SECONDS
    "$BINARY" -i "$DATA/$rep/$gt" --blobsearchonly "$REF/$rep/$REF_FILE" \
        --blob --branchcut-blob "${FLAGS[@]}" "${extra[@]}" \
        --override -o "$out/estimate.nwk" > "$out/run.log" 2>&1
    local rc=$?
    printf 'wall_clock_seconds\n%d\n' "$((SECONDS-t0))" > "$out/row.csv"
    [[ $rc -ne 0 ]] && echo "FAILED $GROUP/$rep rc=$rc" >&2
    return 0
}
export -f run_one
export ROOT BINARY GROUP DATA REF REF_FILE LABEL SCORE_ROOT
export BC_QUAD_OUT="${BC_QUAD_OUT:-}"
export BC_GENES="${BC_GENES:-}"
export FLAGS_STR="${FLAGS[*]}"
# xargs cannot carry an array; re-expand from the string inside the subshell.
mkdir -p "$ROOT/results/runs/$LABEL"
printf 'method=branchcut\nargs=%s\nrunner=local\nlabel=%s\nrefinement=%s/<rep>/%s\n' \
    "${FLAGS[*]}" "$LABEL" "$REF" "$REF_FILE" > "$ROOT/results/runs/$LABEL/run.meta"
echo "=== $LABEL on $GROUP: ${FLAGS[*]} ==="
# REPS="00 01 ..." restricts the run to a subset of replicates.
if [[ -n "${REPS:-}" ]]; then printf '%s\n' $REPS; else ls "$REF" | sort; fi | xargs -P "$NPROC" -I{} bash -c 'FLAGS=($FLAGS_STR); run_one "$@"' _ {}
echo "DONE $LABEL $GROUP"
