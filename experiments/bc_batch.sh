#!/usr/bin/env bash
# Run a list of branch-cut configurations over a list of groups.
#   ./bc_batch.sh "n15 n25" <<'EOF'
#   label1 : --flag a --flag b
#   label2 : --flag c
#   EOF
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BC_GROUPS="${1:?groups}"
NPROC="${2:-8}"
while IFS= read -r line; do
    [[ -z "${line// }" || "$line" == \#* ]] && continue
    label="${line%%:*}"; flags="${line#*:}"
    label="${label// /}"
    for g in $BC_GROUPS; do
        eval "\"$HERE/bc_local.sh\" \"$label\" \"$g\" \"$NPROC\" -- $flags" \
            > /dev/null 2>&1
        echo "done $label $g"
    done
done
