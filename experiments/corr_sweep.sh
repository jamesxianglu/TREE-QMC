#!/usr/bin/env bash
# Split-sample corroboration sweep: --branchcut-corroborate g for g in 1..4,
# on top of the shipped branch-cut configuration, with fixed streams so the
# comparison is paired on the query set.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NPROC="${NPROC:-6}"
BASE=(--oracle 'sym+cf|maj' --branchcut-tau 0.15,0.60 --branchcut-cycles 32
      --branchcut-resolution-margin 0.01 --query-alpha 0.001
      --branchcut-fixed-streams)

for g in "${GROUPS:-n15 n25 n50}"; do :; done
GROUPS="${GROUPS:-n15 n25 n50}"
CORRS="${CORRS:-0 2 3 4}"

for group in $GROUPS; do
  for c in $CORRS; do
    if [[ "$c" == "0" ]]; then
      "$HERE/bc_local.sh" "c_base" "$group" "$NPROC" -- "${BASE[@]}"
    else
      "$HERE/bc_local.sh" "c_corr$c" "$group" "$NPROC" -- "${BASE[@]}" \
          --branchcut-corroborate "$c"
    fi
  done
done
echo "SWEEP DONE"
