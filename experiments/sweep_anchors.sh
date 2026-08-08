#!/usr/bin/env bash
# Sweep --rowsweep-anchors, which repeats the row sweep from several taxa of R
# and rejects if any anchor fires. Already implemented; never swept.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"; ROOT="$(dirname "$SOURCE_DIR")"
BINARY="$SOURCE_DIR/build/tree-qmc"
GROUP="${1:?group}"; ANCHORS_LIST="${2:?comma-separated anchor counts}"; NPROC="${3:-7}"
ORACLE="${ORACLE:-t1+cf|maj}"; TAU="${TAU:-0.45,0.60}"; HEAVY="${HEAVY:-1,2}"; ALPHA="${ALPHA:-0.001}"
export R_HOME="${R_HOME:-$(R RHOME)}"; export PATH="$HOME/.juliaup/bin:$PATH"
DATA="$ROOT/data/camus-dataset/$GROUP"; REF="$ROOT/data/refinements/$GROUP"
run_one() {
    local rep="$1" na="$2" label="$3"
    local out="$ROOT/results/runs/$label/$GROUP/$rep"; mkdir -p "$out"
    [[ -s "$out/estimate.nwk" ]] && return 0
    local gt="iqtree_500.nwk"; [[ -f "$DATA/$rep/$gt" ]] || gt="g_500.nwk"
    local t0=$SECONDS
    "$BINARY" -i "$DATA/$rep/$gt" --blobsearchonly "$REF/$rep/astral4-rooted.tre" \
        --blob --rowsweep-blob --oracle "$ORACLE" --rowsweep-tau "$TAU" \
        --rowsweep-heavy "$HEAVY" --rowsweep-anchors "$na" \
        --query-alpha "$ALPHA" --override -o "$out/estimate.nwk" > "$out/run.log" 2>&1
    local rc=$?
    printf 'wall_clock_seconds\n%d\n' "$((SECONDS-t0))" > "$out/row.csv"
    [[ $rc -ne 0 ]] && echo "FAILED $GROUP/$rep anchors=$na rc=$rc" >&2
    return 0
}
export -f run_one; export ROOT BINARY GROUP DATA REF ORACLE TAU HEAVY ALPHA
for na in ${ANCHORS_LIST//,/ }; do
    label="rs_anch_${na}"; mkdir -p "$ROOT/results/runs/$label"
    cat > "$ROOT/results/runs/$label/run.meta" <<META
method=rowsweep
args=--oracle $ORACLE --rowsweep-tau $TAU --rowsweep-heavy $HEAVY --rowsweep-anchors $na --query-alpha $ALPHA
gene_trees=per-dataset (iqtree_500.nwk, else g_500.nwk)
refinement=data/refinements/<group>/<rep>/astral4-rooted.tre
runner=local
label=$label
META
    echo "=== $label on $GROUP ==="
    ls "$REF" | sort | xargs -P "$NPROC" -I{} bash -c 'run_one "$@"' _ {} "$na" "$label"
done
echo "DONE $GROUP $ANCHORS_LIST"
