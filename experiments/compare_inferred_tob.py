#!/usr/bin/env python3
"""Compare an inferred tree's non-trivial splits with the true tree of blobs."""

import argparse
import os

from common import compute_tob_splits
from tree_of_blobs import splits_from_tob_tree


def safe_ratio(numerator, denominator):
    return numerator / denominator if denominator else float("nan")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("true_network")
    parser.add_argument("inferred_tree")
    parser.add_argument("--truth-output", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.truth_output)), exist_ok=True)
    true_leaves, true_splits, source = compute_tob_splits(
        args.true_network, args.truth_output
    )
    with open(args.inferred_tree) as inferred_file:
        inferred_leaves, inferred_splits = splits_from_tob_tree(
            inferred_file.read()
        )

    if true_leaves != inferred_leaves:
        missing = sorted(true_leaves - inferred_leaves)
        extra = sorted(inferred_leaves - true_leaves)
        parser.error(f"taxon mismatch: missing={missing}, extra={extra}")

    tp = len(inferred_splits & true_splits)
    fp = len(inferred_splits - true_splits)
    fn = len(true_splits - inferred_splits)
    precision = safe_ratio(tp, tp + fp)
    recall = safe_ratio(tp, tp + fn)
    f1 = safe_ratio(2 * precision * recall, precision + recall)

    print("\n-- constructed tree split comparison --")
    print(f"  ground truth source : {source}")
    print(f"  true splits         : {len(true_splits)}")
    print(f"  inferred splits     : {len(inferred_splits)}")
    print(f"  true positives      : {tp}")
    print(f"  false positives     : {fp}")
    print(f"  false negatives     : {fn}")
    print(f"  precision           : {precision:.1%}")
    print(f"  recall              : {recall:.1%}")
    print(f"  F1                  : {f1:.1%}")
    print("  accuracy            : n/a (no negative split set is defined for a tree)")


if __name__ == "__main__":
    main()
