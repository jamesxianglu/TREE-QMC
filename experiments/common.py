"""Shared helpers for the row-sweep experiments: path/env resolution, the
PhyloNetworks tree-of-blobs call, and running/scoring tree-qmc."""

import os
import shutil
import subprocess

from tree_of_blobs import splits_from_tob_tree, tree_of_blobs_splits, split_blob_adjacency

HERE = os.path.dirname(os.path.abspath(__file__))


def camus_root():
    """Locate the data/camus-dataset directory by searching upward from here."""
    d = HERE
    for _ in range(6):
        cand = os.path.join(d, "data", "camus-dataset")
        if os.path.isdir(cand):
            return cand
        d = os.path.dirname(d)
    raise FileNotFoundError("could not find data/camus-dataset above " + HERE)


def default_binary():
    return os.path.join(os.path.dirname(HERE), "build", "tree-qmc")


def r_home():
    if os.environ.get("R_HOME"):
        return os.environ["R_HOME"]
    for cand in ("/opt/homebrew/bin/R", "/usr/local/bin/R", shutil.which("R") or "R"):
        try:
            out = subprocess.run([cand, "RHOME"], capture_output=True, text=True)
            if out.returncode == 0 and out.stdout.strip():
                return out.stdout.strip()
        except FileNotFoundError:
            continue
    return None


def julia_bin():
    for cand in (shutil.which("julia"),
                 os.path.expanduser("~/.juliaup/bin/julia"),
                 "/opt/homebrew/bin/julia"):
        if cand and os.path.exists(cand):
            return cand
    return "julia"


def child_env():
    env = dict(os.environ)
    extra = "/opt/homebrew/bin:" + os.path.expanduser("~/.juliaup/bin")
    env["PATH"] = extra + ":" + env.get("PATH", "")
    rh = r_home()
    if rh:
        env["R_HOME"] = rh
    return env


def compute_tob_splits(net_path, tob_out_path):
    """Compute tree-of-blobs splits for a network.

    Primary path: PhyloNetworks `treeofblobs` (trusted reference) via
    compute_tree_of_blob.jl, then parse the resulting tree's splits.
    Fallback: the built-in bridge finder (only if Julia/PhyloNetworks fails).

    Returns (all_leaves, splits, source_str).
    """
    jl = os.path.join(HERE, "compute_tree_of_blob.jl")
    cmd = [julia_bin(), "--startup-file=no", jl, net_path, tob_out_path]
    proc = subprocess.run(cmd, env=child_env(), capture_output=True, text=True)
    if proc.returncode == 0 and os.path.exists(tob_out_path) and os.path.getsize(tob_out_path):
        all_leaves, splits = splits_from_tob_tree(open(tob_out_path).read())
        return all_leaves, splits, "PhyloNetworks.treeofblobs"
    err = (proc.stderr or proc.stdout or "").strip().replace("\n", " ")[:300]
    all_leaves, splits = tree_of_blobs_splits(open(net_path).read())
    return all_leaves, splits, f"bridge-fallback (julia failed: {err})"


def run_rowsweep(binary, gene_trees, bip_file, out_file, delta, log_file):
    cmd = [binary, "-i", gene_trees, "--rowsweep", bip_file,
           "--rowsweep_out", out_file, "--delta", str(delta)]
    with open(log_file, "w") as log:
        proc = subprocess.run(cmd, env=child_env(), stdout=log, stderr=subprocess.STDOUT)
    return proc.returncode, cmd


def read_labels(bip_file):
    labels, sides = {}, {}
    with open(bip_file) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            labels[p[0]] = int(p[3])
            sides[p[0]] = (p[1], p[2])
    return labels, sides


def read_predictions(out_file):
    preds = {}
    with open(out_file) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or p[0] == "id":
                continue
            preds[p[0]] = int(p[1])
    return preds


def score(labels, preds):
    n_pos = n_neg = fn = fp = missing = 0
    fn_ids, fp_ids = [], []
    for rid, lab in labels.items():
        if rid not in preds:
            missing += 1
            continue
        p = preds[rid]
        if lab == 1:
            n_pos += 1
            if p == 0:
                fn += 1; fn_ids.append(rid)
        else:
            n_neg += 1
            if p == 1:
                fp += 1; fp_ids.append(rid)
    return dict(n_pos=n_pos, n_neg=n_neg, fn=fn, fp=fp, missing=missing,
               fn_ids=fn_ids, fp_ids=fp_ids)
