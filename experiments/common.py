"""Shared path, execution, TSV, and scoring helpers for row-sweep experiments."""

import csv
import os
import shutil
import subprocess
from dataclasses import dataclass

from tree_of_blobs import splits_from_tob_tree, tree_of_blobs_splits

HERE = os.path.dirname(os.path.abspath(__file__))


@dataclass(frozen=True)
class SplitTestScore:
    """Classification counts and the IDs of misclassified bipartitions."""

    positive_count: int
    negative_count: int
    false_negative_ids: tuple
    false_positive_ids: tuple
    missing_count: int

    @property
    def false_negative_count(self):
        return len(self.false_negative_ids)

    @property
    def false_positive_count(self):
        return len(self.false_positive_ids)

    @property
    def false_negative_rate(self):
        return _safe_rate(self.false_negative_count, self.positive_count)

    @property
    def false_positive_rate(self):
        return _safe_rate(self.false_positive_count, self.negative_count)

    @property
    def accuracy(self):
        total = self.positive_count + self.negative_count
        errors = self.false_negative_count + self.false_positive_count
        return (total - errors) / total if total else float("nan")


def _safe_rate(numerator, denominator):
    return numerator / denominator if denominator else float("nan")


def camus_root():
    """Locate the data/camus-dataset directory by searching upward from here."""
    search_directory = HERE
    for _ in range(6):
        candidate = os.path.join(search_directory, "data", "camus-dataset")
        if os.path.isdir(candidate):
            return candidate
        search_directory = os.path.dirname(search_directory)
    raise FileNotFoundError("could not find data/camus-dataset above " + HERE)


def default_binary():
    return os.path.join(os.path.dirname(HERE), "build", "tree-qmc")


def r_home():
    if os.environ.get("R_HOME"):
        return os.environ["R_HOME"]
    candidates = ("/opt/homebrew/bin/R", "/usr/local/bin/R", shutil.which("R") or "R")
    for candidate in candidates:
        try:
            result = subprocess.run([candidate, "RHOME"], capture_output=True, text=True)
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout.strip()
        except FileNotFoundError:
            continue
    return None


def julia_bin():
    candidates = (
        shutil.which("julia"),
        os.path.expanduser("~/.juliaup/bin/julia"),
        "/opt/homebrew/bin/julia",
    )
    for candidate in candidates:
        if candidate and os.path.exists(candidate):
            return candidate
    return "julia"


def child_env():
    env = dict(os.environ)
    extra = "/opt/homebrew/bin:" + os.path.expanduser("~/.juliaup/bin")
    env["PATH"] = extra + ":" + env.get("PATH", "")
    r_home_path = r_home()
    if r_home_path:
        env["R_HOME"] = r_home_path
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
        with open(tob_out_path) as tob_file:
            all_leaves, splits = splits_from_tob_tree(tob_file.read())
        return all_leaves, splits, "PhyloNetworks.treeofblobs"

    julia_error = (proc.stderr or proc.stdout or "").strip().replace("\n", " ")[:300]
    with open(net_path) as network_file:
        all_leaves, splits = tree_of_blobs_splits(network_file.read())
    return all_leaves, splits, f"bridge-fallback (julia failed: {julia_error})"


def run_rowsweep(binary, gene_trees, bip_file, out_file, delta, log_file,
                 query_alpha=0.05):
    """Run one row-sweep command and return its process exit code."""
    cmd = [binary, "-i", gene_trees, "--rowsweep", bip_file,
           "--rowsweep_out", out_file, "--delta", str(delta),
           "--query-alpha", str(query_alpha)]
    with open(log_file, "w") as log:
        proc = subprocess.run(cmd, env=child_env(), stdout=log, stderr=subprocess.STDOUT)
    return proc.returncode


def log_tail(log_file, line_count=20):
    """Return the last few lines of a log, or an empty string if it is absent."""
    if not os.path.exists(log_file):
        return ""
    with open(log_file) as log:
        return "".join(log.readlines()[-line_count:])


def read_labels(bip_file):
    """Read row-sweep ground-truth labels and the two sides of each split."""
    labels, sides = {}, {}
    with open(bip_file) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row_id = row["id"]
            labels[row_id] = int(row["label"])
            sides[row_id] = (row["A"], row["B"])
    return labels, sides


def read_predictions(out_file):
    """Read the binary ACCEPT=1 / REJECT=0 predictions emitted by tree-qmc."""
    preds = {}
    with open(out_file) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            preds[row["id"]] = int(row["prediction"])
    return preds


def score(labels, preds):
    """Compare predictions with labels and return named classification metrics."""
    positive_count = 0
    negative_count = 0
    missing_count = 0
    false_negative_ids = []
    false_positive_ids = []

    for row_id, label in labels.items():
        if row_id not in preds:
            missing_count += 1
            continue

        prediction = preds[row_id]
        if label == 1:
            positive_count += 1
            if prediction == 0:
                false_negative_ids.append(row_id)
        else:
            negative_count += 1
            if prediction == 1:
                false_positive_ids.append(row_id)

    return SplitTestScore(
        positive_count=positive_count,
        negative_count=negative_count,
        false_negative_ids=tuple(false_negative_ids),
        false_positive_ids=tuple(false_positive_ids),
        missing_count=missing_count,
    )
