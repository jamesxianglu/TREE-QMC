#!/usr/bin/env bash
#
# Merge the per-dataset CSV shards of one experiment into a single results CSV.
# The output is named after the run, which is named after the method and its
# parameters, so a directory listing of results/ says which configurations have
# been measured.
#
# Usage:
#   ./collect_tob.sh --run-dir DIR [--out FILE]
#
# Required:
#   --run-dir DIR   a run directory written by submit_tob.sh
#
# Options:
#   --out FILE      output path (default results/<run label>.csv)
#   --allow-partial write the CSV even if some datasets have no shard
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$HERE")")"

usage() { sed -n '3,18p' "$0"; exit "${1:-1}"; }
[[ $# -eq 0 ]] && usage

RUN_DIR=""; OUT=""; ALLOW_PARTIAL=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-dir)        RUN_DIR="${2:-}"; shift 2 ;;
        --out)            OUT="${2:-}";     shift 2 ;;
        --allow-partial)  ALLOW_PARTIAL=1;  shift ;;
        -h|--help)        usage 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage ;;
    esac
done
[[ -n "$RUN_DIR" ]] || { echo "ERROR: --run-dir is required" >&2; exit 1; }
[[ -d "$RUN_DIR" ]] || { echo "ERROR: no such run directory: $RUN_DIR" >&2; exit 1; }
[[ -f "$RUN_DIR/run.meta" ]] || { echo "ERROR: $RUN_DIR has no run.meta" >&2; exit 1; }

LABEL="$(sed -n 's/^label=//p' "$RUN_DIR/run.meta")"
[[ -n "$LABEL" ]] || LABEL="$(basename "$RUN_DIR")"
[[ -n "$OUT" ]] || OUT="$PROJECT_ROOT/results/$LABEL.csv"

# A dataset directory without a shard is a job that failed or has not finished;
# silently dropping it would understate the error of whichever configuration
# happened to crash.
expected=0
missing=0
while IFS= read -r d; do
    expected=$((expected + 1))
    [[ -s "$d/row.csv" ]] || { missing=$((missing + 1)); echo "missing: ${d#$RUN_DIR/}" >&2; }
done < <(find "$RUN_DIR" -mindepth 2 -maxdepth 2 -type d | sort)

if [[ "$missing" -gt 0 && "$ALLOW_PARTIAL" -eq 0 ]]; then
    echo "ERROR: $missing of $expected datasets have no shard; rerun them or pass --allow-partial" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUT")"
{
    echo "method,oracle,args,taxa,network_id,gene_trees,nl,i1,i2,fn,fp,fnr,fpr,fnr_fpr_avg,err,nrf,wall_clock_seconds,max_rss_kb,host"
    # Group order by taxon count, then dataset id, so the file reads naturally.
    find "$RUN_DIR" -mindepth 3 -maxdepth 3 -name row.csv |
        awk -F/ '{ g=$(NF-2); sub(/^n/, "", g); print g "\t" $(NF-1) "\t" $0 }' |
        sort -k1,1n -k2,2 | cut -f3 | xargs -I{} cat {}
} > "$OUT"

rows=$(($(wc -l < "$OUT") - 1))
echo "wrote $OUT ($rows rows, $missing missing)"
