#!/usr/bin/env bash
# Branch cut with the FAITHFUL cycle-cover sampling of Lem. cycle-cover
# (--branchcut-cycles h) against the older uniform random sampling
# (--branchcut-samples s). Budget-matched at h == s.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"; ROOT="$(dirname "$SOURCE_DIR")"
BINARY="$SOURCE_DIR/build/tree-qmc"
GROUP="${1:?group}"; SPECS="${2:?comma-separated, e.g. cyc4,cyc8,cyc16,rnd8}"; NPROC="${3:-7}"
ORACLE="${ORACLE:-sym+cf|maj}"; TAU="${TAU:-0.15,0.60}"; ALPHA="${ALPHA:-0.001}"
MINSUP="${MINSUP:-4}"
export R_HOME="${R_HOME:-$(R RHOME)}"; export PATH="$HOME/.juliaup/bin:$PATH"
DATA="$ROOT/data/camus-dataset/$GROUP"; REF="$ROOT/data/refinements/$GROUP"
run_one() {
    local rep="$1" spec="$2" label="$3"
    local out="$ROOT/results/runs/$label/$GROUP/$rep"; mkdir -p "$out"
    [[ -s "$out/estimate.nwk" ]] && return 0
    local gt="iqtree_500.nwk"; [[ -f "$DATA/$rep/$gt" ]] || gt="g_500.nwk"
    local flag=()
    case "$spec" in
        cyc*) flag=(--branchcut-cycles "${spec#cyc}") ;;
        rnd*) flag=(--branchcut-samples "${spec#rnd}") ;;
    esac
    local t0=$SECONDS
    "$BINARY" -i "$DATA/$rep/$gt" --blobsearchonly "$REF/$rep/astral4-rooted.tre" \
        --blob --branchcut-blob --oracle "$ORACLE" --branchcut-tau "$TAU" \
        --branchcut-min-support "$MINSUP" "${flag[@]}" \
        --query-alpha "$ALPHA" --override -o "$out/estimate.nwk" > "$out/run.log" 2>&1
    local rc=$?
    printf 'wall_clock_seconds\n%d\n' "$((SECONDS-t0))" > "$out/row.csv"
    [[ $rc -ne 0 ]] && echo "FAILED $GROUP/$rep $spec rc=$rc" >&2
    return 0
}
export -f run_one; export ROOT BINARY GROUP DATA REF ORACLE TAU ALPHA MINSUP
for spec in ${SPECS//,/ }; do
    label="bc_${spec}"; mkdir -p "$ROOT/results/runs/$label"
    printf 'method=branchcut\nargs=--oracle %s --branchcut-tau %s --branchcut-min-support %s (%s)\ngene_trees=per-dataset\nrunner=local\nlabel=%s\n' \
        "$ORACLE" "$TAU" "$MINSUP" "$spec" "$label" > "$ROOT/results/runs/$label/run.meta"
    echo "=== $label on $GROUP ==="
    ls "$REF" | sort | xargs -P "$NPROC" -I{} bash -c 'run_one "$@"' _ {} "$spec" "$label"
done
echo "DONE $GROUP $SPECS"
