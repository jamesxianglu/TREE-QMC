#!/usr/bin/env bash
#
# Run one tree-of-blobs experiment: one method, one configuration, over a chosen
# set of datasets, either on SLURM or locally. Every experiment lives in its own
# run directory named after the method and its parameters, and each dataset
# leaves a one-row CSV shard there; collect_tob.sh merges them.
#
# The method, its flags, the datasets and the runner have to be typed out: a run
# is only reproducible if its configuration was written down. The gene tree file
# is the one exception, and is never typed: it is a property of the dataset
# rather than of the experiment, so it is resolved per dataset -- iqtree_500.nwk
# when the dataset has it, g_500.nwk otherwise -- and the file that was actually
# used is recorded in every row of the output. See the note under --datasets
# about what that costs.
#
# Usage:
#   ./submit_tob.sh --method METHOD --args 'FLAGS' --datasets 'SPEC' \
#                   --runner (local|slurm) [options]
#
# Required:
#   --method METHOD      rowsweep | cornerrow | branchcut
#   --args 'FLAGS'       tree-qmc flags for this configuration, quoted as one
#                        string, e.g.
#                          '--oracle t1 --delta 0.25 --query-alpha 0.001'
#                          '--oracle t1+cf|maj --rowsweep-tau 0.45,0.60 \
#                           --rowsweep-heavy 1,2 --query-alpha 0.001'
#
# Options:
#   --datasets 'SPEC'    space-separated groups, each optionally narrowed by
#                        ':ids', e.g. 'n15 n25 n50:00-09 n100:07'
#
#                        The gene tree file is not selectable: each dataset uses
#                        its own iqtree_500.nwk, or g_500.nwk when it has no
#                        iqtree_500.nwk. In this data that is uniform per group
#                        except for n25, where 20 of the 50 replicates have both
#                        files and so take iqtree_500.nwk while the other 30 take
#                        g_500.nwk. A selection that lands on more than one file
#                        is not a like-for-like comparison, and since there is no
#                        flag to correct it, the only correction is to narrow
#                        --datasets to replicates that share one file. The
#                        breakdown is printed before anything is submitted.
#   --runner RUNNER      local | slurm
#
#   --label NAME         run directory name; defaults to a slug built from the
#                        method and flags
#   --out-root DIR       where run directories live (default results/runs)
#   --refinement DIR/FILE
#                        the binary refinement to contract, as a directory laid
#                        out <group>/<rep>/<file>. Default
#                        data/refinements/astral4-rooted.tre -- the ASTRAL-IV
#                        species tree, which is *not* a refinement of the true
#                        tree of blobs. Pass
#                        data/refinements_true/true-refinement.tre (written by
#                        results/analysis/make_true_refinement.py) to run the
#                        theory's own setting, where every spurious edge is
#                        blob-internal by construction. Recorded in run.meta.
#   --sbatch-args 'STR'  extra sbatch options, prepended to the ones this script
#                        sets, so the resource policy baked into run_tob.sbatch
#                        can be overridden without editing it, e.g.
#                          --sbatch-args '--mem=16G --qos=normal'
#   --skip-missing       drop datasets whose gene tree file is absent instead of
#                        refusing to run. Off by default: a silently shrinking
#                        dataset set is how a configuration ends up compared
#                        against a different set of networks than its rival.
#   --python PATH        interpreter for compare_two_trees.py, which needs
#                        treeswift. Defaults to <project root>/.venv/bin/python
#                        if that exists and to python3 otherwise -- note that
#                        this is the project root, not this directory. Whichever
#                        is chosen is printed with its version before anything
#                        runs.
#
#                        A bare name is resolved to an absolute path here and
#                        the path is what the jobs are given. A name is resolved
#                        against the PATH of whichever machine looks it up, and
#                        that is not the same machine that runs the job: on a
#                        cluster whose login node has a python3 the compute
#                        nodes do not, "python3" passes every check here and
#                        then runs as some other interpreter there.
#   --modules 'LIST'     space-separated modulefiles to load inside each job
#                        before anything else, e.g. --modules 'Python3/3.11.11'.
#                        A job starts from a bare environment, so a module
#                        loaded in your login shell is not loaded there. The
#                        same modules are loaded for the checks below, so what
#                        is verified here is what the jobs will run. Versions
#                        are pinned deliberately: a site upgrade must not change
#                        interpreters halfway through a sweep.
#
# Neither --python nor --modules normally has to be given. Both describe the
# machine rather than the experiment, so they belong in experiments/site.conf,
# which is sourced here when it exists and is overridden by the flags:
#
#     # experiments/site.conf -- not shared between machines
#     PYTHON_OPT=/fs/cbcb-lab/ekmolloy/jameslu/.venv/bin/python
#     MODULES='Python3/3.11.11'
#
# The values it supplies are printed before the run and written to run.meta, so
# a run that leaned on it records the same facts as one that spelled them out.
# (experiments/site_env.sh is a different file, sourced inside the job for
# anything site.conf cannot express; it is not needed for the above.)
#   --status-wait SECS   after a slurm submission, wait this long and report what
#                        sacct says happened (default 20, 0 to skip). sbatch
#                        accepting a job says nothing about whether it ran.
#   --dry-run            print what would happen and stop
#   --force              redo datasets that already have a shard
#   --no-build           reuse the existing binary instead of rebuilding
#
# Example:
#   ./submit_tob.sh --method rowsweep --runner slurm \
#       --datasets 'n50 n100 n150 n200' \
#       --args '--oracle t1+cf|maj --rowsweep-tau 0.45,0.60 --rowsweep-heavy 1,2 --query-alpha 0.001'
#
# n15 and n25 are left out of that example on purpose. Those four groups are
# g_500.nwk throughout, so they compare against each other; n15 is iqtree_500.nwk
# and n25 is split between the two, so each belongs to its own row set.
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(dirname "$HERE")"
PROJECT_ROOT="$(dirname "$SOURCE_DIR")"
DATA_ROOT="$PROJECT_ROOT/data/camus-dataset"
REFINEMENT_ROOT="$PROJECT_ROOT/data/refinements"
REFINEMENT_NAME="astral4-rooted.tre"
REFINEMENT_SPEC=""
BUILD_DIR="$SOURCE_DIR/build"
TRUE_TOB_DIR="$PROJECT_ROOT/results/true_tob"
RUNNER_SCRIPT="$HERE/run_tob.sbatch"
TOB_JL="$HERE/compute_tree_of_blob.jl"
COMPARE_PY="$PROJECT_ROOT/results/analysis/compare_two_trees.py"

# The header above is the usage message, so it is printed rather than restated:
# a second copy is a copy that goes stale.
usage() { sed -n '3,/^set -euo/p' "$0" | sed '$d'; exit "${1:-1}"; }
[[ $# -eq 0 ]] && usage

METHOD=""; METHOD_ARGS=""; DATASETS=""; RUNNER=""
LABEL=""; OUT_ROOT="$PROJECT_ROOT/results/runs"; SBATCH_ARGS=""; PYTHON_OPT=""
MODULES=""
DRY_RUN=0; FORCE=0; BUILD=1; SKIP_MISSING=0; STATUS_WAIT=20

# --- per-machine defaults ---------------------------------------------------
# Which interpreter to use and which modules to load are facts about the machine,
# not about the experiment, and they do not change between runs -- so they are
# read from a file on that machine rather than retyped on every submission. An
# option that has to be repeated correctly every time is an option that will
# eventually be left out, and leaving out --modules does not fail loudly; it runs
# under a different interpreter and writes shards that look like all the others.
#
# Sourced before the options are parsed, so a flag on the command line still
# wins. Whatever it sets is printed below and recorded in run.meta, which is what
# keeps a defaulted run as reproducible as a fully typed one.
SITE_CONF="$HERE/site.conf"
if [[ -f "$SITE_CONF" ]]; then
    # shellcheck disable=SC1090
    source "$SITE_CONF"
fi
SITE_PYTHON="$PYTHON_OPT"; SITE_MODULES="$MODULES"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --method)       METHOD="${2:-}";      shift 2 ;;
        --args)         METHOD_ARGS="${2:-}"; shift 2 ;;
        --datasets)     DATASETS="${2:-}";    shift 2 ;;
        --runner)       RUNNER="${2:-}";      shift 2 ;;
        --label)        LABEL="${2:-}";       shift 2 ;;
        --out-root)     OUT_ROOT="${2:-}";    shift 2 ;;
        --refinement)   REFINEMENT_SPEC="${2:-}"; shift 2 ;;
        --sbatch-args)  SBATCH_ARGS="${2:-}"; shift 2 ;;
        --python)       PYTHON_OPT="${2:-}";  shift 2 ;;
        --modules)      MODULES="${2:-}";     shift 2 ;;
        --skip-missing) SKIP_MISSING=1;       shift ;;
        --status-wait)  STATUS_WAIT="${2:-}"; shift 2 ;;
        --dry-run)      DRY_RUN=1;            shift ;;
        --force)        FORCE=1;              shift ;;
        --no-build)     BUILD=0;              shift ;;
        -h|--help)      usage 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage ;;
    esac
done

if [[ -n "$REFINEMENT_SPEC" ]]; then
    # "DIR/FILE": the last component is the per-replicate file name, the rest is
    # the root holding <group>/<rep>. A bare directory keeps the default name.
    if [[ "$REFINEMENT_SPEC" == */*.tre || "$REFINEMENT_SPEC" == */*.nwk ]]; then
        REFINEMENT_NAME="$(basename "$REFINEMENT_SPEC")"
        REFINEMENT_SPEC="$(dirname "$REFINEMENT_SPEC")"
    fi
    [[ "$REFINEMENT_SPEC" = /* ]] || REFINEMENT_SPEC="$PROJECT_ROOT/$REFINEMENT_SPEC"
    REFINEMENT_ROOT="$REFINEMENT_SPEC"
fi

missing=""
[[ -n "$METHOD" ]]      || missing="$missing --method"
[[ -n "$METHOD_ARGS" ]] || missing="$missing --args"
[[ -n "$DATASETS" ]]    || missing="$missing --datasets"
[[ -n "$RUNNER" ]]      || missing="$missing --runner"
if [[ -n "$missing" ]]; then
    echo "ERROR: required and not supplied:$missing" >&2
    exit 1
fi
case "$METHOD" in rowsweep|cornerrow|branchcut) ;; *)
    echo "ERROR: --method must be rowsweep, cornerrow or branchcut, got '$METHOD'" >&2; exit 1 ;;
esac
case "$RUNNER" in local|slurm) ;; *)
    echo "ERROR: --runner must be local or slurm, got '$RUNNER'" >&2; exit 1 ;;
esac
if [[ "$RUNNER" == slurm ]] && ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: --runner slurm but sbatch is not on PATH" >&2
    exit 1
fi

# --- run identity ----------------------------------------------------------
# The slug has to survive being a directory name and still say what was run, so
# flag punctuation collapses to dashes: "--oracle t1+cf|maj --rowsweep-tau
# 0.45,0.60" becomes "oracle-t1+cf-maj-rowsweep-tau-0.45-0.60".
slugify() {
    printf '%s' "$1" |
        sed -e 's/--//g' -e 's/[^A-Za-z0-9._+-]\{1,\}/-/g' \
            -e 's/-\{2,\}/-/g' -e 's/^-//' -e 's/-$//'
}
if [[ -z "$LABEL" ]]; then
    # The gene tree file is not part of the label. It used to be, back when it
    # was chosen on the command line and two sweeps of one configuration could
    # differ by nothing else -- the second would then write its shards into the
    # first one's directory. Now the file follows from the dataset, so a label
    # and a dataset list determine it, and there is nothing left to disambiguate.
    LABEL="${METHOD}__$(slugify "$METHOD_ARGS")"
fi
RUN_DIR="$OUT_ROOT/$LABEL"

# --- dataset expansion -----------------------------------------------------
WORK="$(mktemp)"
JOBS="$(mktemp)"
trap 'rm -f "$WORK" "$JOBS"' EXIT
for token in $DATASETS; do
    group="${token%%:*}"
    [[ "$group" =~ ^n[0-9]+$ ]] || { echo "ERROR: bad group in --datasets: $token" >&2; exit 1; }
    [[ -d "$REFINEMENT_ROOT/$group" ]] || { echo "ERROR: no refinements for $group" >&2; exit 1; }
    if [[ "$token" == *:* ]]; then
        selector="${token#*:}"
        for part in ${selector//,/ }; do
            start="${part%%-*}"; end="${part#*-}"; [[ "$part" == *-* ]] || end="$start"
            [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ ]] ||
                { echo "ERROR: bad id selector: $part" >&2; exit 1; }
            for ((i = 10#$start; i <= 10#$end; i++)); do
                printf '%s\t%0*d\n' "$group" "${#start}" "$i" >> "$WORK"
            done
        done
    else
        find "$REFINEMENT_ROOT/$group" -mindepth 1 -maxdepth 1 -type d \
            -exec basename {} \; | sort | sed "s|^|$group\t|" >> "$WORK"
    fi
done
[[ -s "$WORK" ]] || { echo "ERROR: --datasets selected nothing" >&2; exit 1; }

# --- the module environment -------------------------------------------------
# Loaded here as well as in the job, so the interpreter checked below is the one
# the jobs will run. Kept in step with the block in run_tob.sbatch; the two
# scripts each stand alone, which is worth ten duplicated lines.
#
# `module` is a shell function, absent from non-interactive shells until its
# init script is sourced, and those init scripts are not written to survive
# `set -u`, hence the fence.
load_modules() {
    [[ -n "$MODULES" ]] || return 0
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
        echo "ERROR: --modules given but there is no module command on $(hostname)" >&2
        return 1
    fi
    local m
    for m in $MODULES; do
        if ! module load "$m"; then
            set -u
            echo "ERROR: module load $m failed on $(hostname)" >&2
            return 1
        fi
    done
    set -u
}
# A failed load must stop the run. Loading nothing and carrying on is how a
# sweep ends up running under whatever interpreter happened to be on PATH --
# which is the failure this option exists to prevent.
load_modules || exit 1

# Picking an interpreter by guessing at a path is how a forgotten virtualenv
# gets used for months without anyone noticing, so the choice is always stated
# out loud below, and --python overrides it.
if [[ -n "$PYTHON_OPT" ]]; then
    PYTHON="$PYTHON_OPT"
    # Named apart so the summary says where the choice came from. "--python" for
    # something never typed sends the next reader looking through their shell
    # history for a flag that is not there.
    if [[ "$PYTHON_OPT" == "$SITE_PYTHON" ]]; then
        PYTHON_SOURCE="site.conf"
    else
        PYTHON_SOURCE="--python"
    fi
elif [[ -x "$PROJECT_ROOT/.venv/bin/python" ]]; then
    PYTHON="$PROJECT_ROOT/.venv/bin/python"
    PYTHON_SOURCE="\$PROJECT_ROOT/.venv"
else
    PYTHON="python3"
    PYTHON_SOURCE="PATH"
fi

# Resolve a bare name to a path now, and send the path to the jobs. A name is
# resolved by whoever looks it up: "python3" checked here is the login node's
# python3, and the same word on a compute node with a different PATH is a
# different interpreter that was never checked. That difference does not
# announce itself -- it surfaces as an error from deep inside a library, on a
# machine you are not logged into.
if [[ "$PYTHON" != */* ]]; then
    if resolved="$(command -v "$PYTHON" 2>/dev/null)"; then
        PYTHON_SOURCE="$PYTHON_SOURCE, resolved to a path here"
        PYTHON="$resolved"
    else
        echo "ERROR: no '$PYTHON' on PATH (via $PYTHON_SOURCE)" >&2
        echo "       Pass --python /absolute/path/to/python." >&2
        exit 1
    fi
fi

# --- preflight: the tools the jobs will need --------------------------------
# These live on the submitting side but are used on the compute side, and a job
# that cannot find them dies in milliseconds. Checking here turns one error
# message into one error message instead of into a few hundred failed jobs.
toolerr=0
for f in "$TOB_JL" "$COMPARE_PY"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; toolerr=1; }
done
if ! PYTHON_VERSION="$("$PYTHON" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])' 2>/dev/null)"; then
    echo "ERROR: cannot run $PYTHON (found via $PYTHON_SOURCE)" >&2
    PYTHON_VERSION="?"
    toolerr=1
elif ! "$PYTHON" -c 'import sys; sys.exit(sys.version_info < (3, 6))' 2>/dev/null; then
    # Checked separately from the import because the import failure it causes
    # does not look like a version problem: current treeswift uses f-strings, so
    # an old interpreter reports a SyntaxError from a line inside the library,
    # which reads as a corrupt install rather than as the wrong python.
    echo "ERROR: $PYTHON is python $PYTHON_VERSION (found via $PYTHON_SOURCE)," >&2
    echo "       and treeswift needs 3.6 or newer. Point --python at a newer" >&2
    echo "       interpreter; installing treeswift into this one cannot help." >&2
    toolerr=1
elif ! import_error="$("$PYTHON" -c 'import treeswift' 2>&1)"; then
    # Naming the interpreter and its version matters more than naming the
    # module: the usual cause is that this is not the interpreter that pip was
    # pointed at, and a version that looks nothing like the one just installed
    # into says so at a glance. The import's own error follows, since a broken
    # install and an absent one need different fixes.
    echo "ERROR: $PYTHON (python $PYTHON_VERSION, found via $PYTHON_SOURCE)" >&2
    echo "       cannot import treeswift, which compare_two_trees.py needs." >&2
    printf '%s\n' "$import_error" | sed 's/^/         /' >&2
    echo "       Install it into that interpreter:" >&2
    echo "         $PYTHON -m pip install treeswift" >&2
    echo "       or point elsewhere with --python /path/to/python." >&2
    toolerr=1
fi

# An interpreter under a machine-local prefix is the one failure this preflight
# cannot see: it runs fine here and is simply absent on the compute node, where
# the job falls back to whatever else answers to that name. Checked by prefix
# rather than by asking a compute node, since one probe job per submission costs
# more than it is worth -- but say so, because the symptom is unrecognisable.
if [[ "$RUNNER" == slurm ]]; then
    case "$PYTHON" in
        /opt/*|/usr/local/*|/tmp/*|/var/*|/scratch/*)
            echo "WARNING: $PYTHON is under a prefix that is often local to the machine" >&2
            echo "         it is on, and the jobs run elsewhere. Confirm the compute" >&2
            echo "         nodes have it:" >&2
            echo "           srun --partition=cbcb --account=cbcb $PYTHON -c 'import treeswift'" >&2
            echo "         An interpreter on shared storage avoids the question." >&2
            ;;
    esac
fi
if ! PATH="$HOME/.juliaup/bin:/opt/homebrew/bin:$PATH" command -v julia >/dev/null 2>&1; then
    echo "WARNING: julia is not on PATH here. It is only needed for datasets whose" >&2
    echo "         true tree of blobs is not yet cached under $TRUE_TOB_DIR," >&2
    echo "         and what matters is the compute node, not this one -- but if it" >&2
    echo "         is missing there too, every such job will fail immediately." >&2
fi
[[ "$toolerr" -eq 0 ]] || exit 1

# --- preflight: the datasets ------------------------------------------------
# Every input is resolved and checked before anything is built or submitted. A
# missing file used to abort partway through the submission loop, which left a
# run half in the queue and reported nothing about the datasets it never
# reached; the whole selection is validated first so the answer is all-or-
# nothing and every problem is listed at once.
BAD="$(mktemp)"
trap 'rm -f "$WORK" "$WORK.inputs" "$JOBS" "$BAD"' EXIT

while IFS=$'\t' read -r group rep; do
    dataset_dir="$DATA_ROOT/$group/$rep"
    refinement="$REFINEMENT_ROOT/$group/$rep/$REFINEMENT_NAME"
    true_network="$dataset_dir/true_net.nwk"

    # The gene trees follow from the dataset and are never passed in: prefer the
    # estimated trees, fall back to the true ones. Preferring iqtree_500.nwk is
    # what makes n25 inhomogeneous, since only some of its replicates have it.
    gene_trees=""
    if [[ -f "$dataset_dir/iqtree_500.nwk" ]]; then
        gene_trees="$dataset_dir/iqtree_500.nwk"
    elif [[ -f "$dataset_dir/g_500.nwk" ]]; then
        gene_trees="$dataset_dir/g_500.nwk"
    fi

    trouble=""
    [[ -n "$gene_trees" ]] ||
        trouble="no iqtree_500.nwk or g_500.nwk in $dataset_dir"
    [[ -f "$refinement" ]]   || trouble="${trouble:+$trouble; }missing $refinement"
    [[ -f "$true_network" ]] || trouble="${trouble:+$trouble; }missing $true_network"
    if [[ -n "$trouble" ]]; then
        printf '%s/%s\t%s\n' "$group" "$rep" "$trouble" >> "$BAD"
        continue
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$group" "$rep" "$gene_trees" "$refinement" "$true_network" >> "$JOBS"
done < "$WORK"

nbad="$(wc -l < "$BAD" | tr -d ' ')"
njobs="$(wc -l < "$JOBS" | tr -d ' ')"

echo "method      $METHOD"
echo "args        $METHOD_ARGS"
echo "gene trees  per dataset: iqtree_500.nwk, else g_500.nwk"
echo "python      $PYTHON (python $PYTHON_VERSION, via $PYTHON_SOURCE)"
if [[ -n "$MODULES" && "$MODULES" == "$SITE_MODULES" ]]; then
    echo "modules     $MODULES (site.conf)"
else
    echo "modules     ${MODULES:-none}"
fi
echo "runner      $RUNNER"
echo "refinement  $REFINEMENT_ROOT/<group>/<rep>/$REFINEMENT_NAME"
echo "run dir     $RUN_DIR"
echo "datasets    $njobs usable of $(wc -l < "$WORK" | tr -d ' ') selected"

if [[ "$nbad" -gt 0 ]]; then
    echo
    echo "$nbad selected dataset(s) cannot be run:" >&2
    head -20 "$BAD" >&2
    [[ "$nbad" -gt 20 ]] && echo "  ... and $((nbad - 20)) more" >&2
    if [[ "$SKIP_MISSING" -eq 0 ]]; then
        echo >&2
        echo "Nothing was submitted. Narrow --datasets to what exists, or pass" >&2
        echo "--skip-missing to run the rest anyway." >&2
        exit 1
    fi
    echo "--skip-missing given: running the remaining $njobs." >&2
fi
[[ "$njobs" -gt 0 ]] || { echo "ERROR: no runnable datasets" >&2; exit 1; }

# Mixing estimated and true gene trees across a comparison silently changes what
# the error numbers mean, and resolving the file per dataset does exactly that
# whenever a selection spans datasets that do not have the same files. Since the
# file is no longer selectable, this warning is the only place it is caught, so
# it reports the split per group rather than in total: a run that is one file per
# group is a set of comparable groups, and one that splits inside a group is not
# comparable even against itself. That distinction is invisible in the totals.
awk -F'\t' '{ n = split($3, p, "/"); print $1 "\t" p[n] }' "$JOBS" |
    sort | uniq -c | awk '{ print $2, $3, $1 }' > "$WORK.inputs"
if [[ "$(awk '{ print $2 }' "$WORK.inputs" | sort -u | wc -l | tr -d ' ')" -gt 1 ]]; then
    echo
    echo "WARNING: this selection does not use one gene tree file throughout." >&2
    echo "         The error numbers mean different things across the split, so" >&2
    echo "         the rows are not directly comparable:" >&2
    awk '{ printf "           %-6s %-16s %3d replicate(s)\n", $1, $2, $3 }' \
        "$WORK.inputs" >&2
    # Named because it is the only group in this data that splits, and a split
    # inside one group is the case that looks like a normal run and is not one.
    if [[ "$(awk '$1 == "n25"' "$WORK.inputs" | wc -l | tr -d ' ')" -gt 1 ]]; then
        echo "         n25 is split within the group: only 20 of its 50 replicates" >&2
        echo "         have iqtree_500.nwk, and those 20 take it. Narrow --datasets" >&2
        echo "         to one side of that split before comparing within n25." >&2
    fi
fi

# --- build -----------------------------------------------------------------
BINARY="$RUN_DIR/tree-qmc"
if [[ "$DRY_RUN" -eq 0 ]]; then
    mkdir -p "$RUN_DIR"
    if [[ "$BUILD" -eq 1 ]]; then
        echo "Building TREE-QMC..."
        # cmake's output goes to a log, so on failure the log is what has to be
        # shown; without this the script died here with an empty terminal.
        if ! cmake -S "$SOURCE_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release \
                > "$RUN_DIR/configure.log" 2>&1; then
            echo "ERROR: cmake configure failed; tail of $RUN_DIR/configure.log:" >&2
            tail -20 "$RUN_DIR/configure.log" >&2
            exit 1
        fi
        if ! cmake --build "$BUILD_DIR" -j4 > "$RUN_DIR/build.log" 2>&1; then
            echo "ERROR: build failed; tail of $RUN_DIR/build.log:" >&2
            tail -20 "$RUN_DIR/build.log" >&2
            exit 1
        fi
    fi
    [[ -x "$BUILD_DIR/tree-qmc" ]] || { echo "ERROR: no binary at $BUILD_DIR/tree-qmc" >&2; exit 1; }
    # Pin the binary to the run so a later rebuild cannot change what these
    # numbers mean.
    cp -f "$BUILD_DIR/tree-qmc" "$BINARY"
    {
        echo "method=$METHOD"
        echo "args=$METHOD_ARGS"
        echo "gene_trees=per-dataset (iqtree_500.nwk, else g_500.nwk)"
        # Which refinement was contracted decides what the error means: the
        # ASTRAL tree is not a refinement of the true tree of blobs, and a
        # correct one makes every spurious edge blob-internal by construction.
        echo "refinement=$REFINEMENT_ROOT/<group>/<rep>/$REFINEMENT_NAME"
        echo "runner=$RUNNER"
        # Which interpreter and which modules produced these numbers is part of
        # the run: the same flags under a different python are a different run.
        echo "python=$PYTHON ($PYTHON_VERSION)"
        echo "modules=${MODULES:-none}"
        echo "label=$LABEL"
        echo "submitted=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } > "$RUN_DIR/run.meta"
    # The dataset list is part of the run: collect_tob.sh needs to know what was
    # meant to exist, not just what happens to be on disk.
    cut -f1,2 "$JOBS" > "$RUN_DIR/datasets.tsv"
    : > "$RUN_DIR/jobids.txt"
fi

submitted=0
skipped=0
failed=0
# The dataset table is read on fd 3: the runner and everything it starts inherit
# stdin, and a local run would otherwise let julia or tree-qmc swallow the
# remaining lines.
while IFS=$'\t' read -r group rep gene_trees refinement true_network <&3; do
    out_dir="$RUN_DIR/$group/$rep"
    if [[ -s "$out_dir/row.csv" && "$FORCE" -eq 0 ]]; then
        skipped=$((skipped + 1))
        continue
    fi

    args=("$METHOD" "$METHOD_ARGS" "$group" "$rep" "$gene_trees" "$refinement"
          "$true_network" "$out_dir" "$TRUE_TOB_DIR" "$BINARY" "$TOB_JL"
          "$COMPARE_PY" "$PYTHON" "$MODULES")
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "would run $group/$rep with $(basename "$gene_trees")"
        submitted=$((submitted + 1))
        continue
    fi
    mkdir -p "$out_dir"
    if [[ "$RUNNER" == slurm ]]; then
        # One rejected job must not abort the other 269, and the job id is the
        # only handle on a submitted job, so it is kept rather than discarded.
        if jobid="$(sbatch --parsable ${SBATCH_ARGS} \
                       --job-name="${METHOD}-${group}-${rep}" \
                       --output="$out_dir/slurm-%j.out" \
                       "$RUNNER_SCRIPT" "${args[@]}" 2>&1)"; then
            printf '%s\t%s/%s\n' "$jobid" "$group" "$rep" >> "$RUN_DIR/jobids.txt"
            submitted=$((submitted + 1))
        else
            echo "sbatch refused $group/$rep: $jobid" >&2
            failed=$((failed + 1))
            [[ "$failed" -ge 5 ]] && { echo "ERROR: 5 rejections, stopping" >&2; break; }
        fi
    else
        bash "$RUNNER_SCRIPT" "${args[@]}" < /dev/null || {
            echo "run failed for $group/$rep" >&2
            failed=$((failed + 1))
        }
        submitted=$((submitted + 1))
    fi
done 3< "$JOBS"

echo "$submitted dataset(s) $( [[ "$RUNNER" == slurm ]] && echo submitted || echo run ), $skipped already done, $failed failed"
[[ "$DRY_RUN" -eq 1 ]] && exit 0

# --- did they actually land? ------------------------------------------------
# sbatch returning a job id says only that the scheduler accepted the script. A
# job whose script dies on its first line is accepted just the same, and a few
# hundred of them clear the queue faster than squeue can be typed -- which reads
# exactly like nothing was ever submitted. So look, rather than report success
# and leave.
if [[ "$RUNNER" == slurm && "$submitted" -gt 0 ]]; then
    echo "Job ids in $RUN_DIR/jobids.txt"
    if [[ "$STATUS_WAIT" -gt 0 ]] && command -v sacct >/dev/null 2>&1; then
        echo "Waiting ${STATUS_WAIT}s, then asking sacct what became of them..."
        sleep "$STATUS_WAIT"
        # --parsable prints "jobid" or, under multi-cluster, "jobid;cluster".
        ids="$(cut -f1 "$RUN_DIR/jobids.txt" | cut -d';' -f1 | paste -sd, -)"
        echo "states after ${STATUS_WAIT}s:"
        sacct -j "$ids" -X -n -o State%-20 | awk 'NF {print $1}' |
            sort | uniq -c | sed 's/^/  /'
        done_now="$(find "$RUN_DIR" -mindepth 3 -maxdepth 3 -name row.csv | wc -l | tr -d ' ')"
        echo "  $done_now of $submitted shard(s) written so far"
        if grep -qE 'FAILED|CANCELLED|TIMEOUT|NODE_FAIL|OUT_OF_ME' \
               <(sacct -j "$ids" -X -n -o State 2>/dev/null); then
            first_out="$(find "$RUN_DIR" -name 'slurm-*.out' -size +0 | sort | head -1)"
            echo
            echo "Some jobs have already ended badly. First job output:" >&2
            [[ -n "$first_out" ]] && { echo "  $first_out" >&2; sed 's/^/  /' "$first_out" | head -25 >&2; }
        fi
    fi
fi
echo "Collect with: $HERE/collect_tob.sh --run-dir $RUN_DIR"
[[ "$failed" -eq 0 ]]
