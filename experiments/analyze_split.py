#!/usr/bin/env python3
"""Analyse row_sweep_test_idx on ONE dataset at ONE delta, with a weakness report.

Ground truth = PhyloNetworks' tree of blobs. For every true split (positive) and
a mixed set of non-splits (negatives) we run tree-qmc --rowsweep and report
false-positive / false-negative rates, then break the errors down structurally
(blob adjacency, split size, hard-vs-random negatives) to expose where and why
the algorithm fails.
"""

import argparse
import math
import os
from dataclasses import dataclass

from common import (camus_root, compute_tob_splits, default_binary, log_tail,
                    read_labels, read_predictions, run_rowsweep, score)
from gen_bipartitions import build_bipartition_file
from tree_of_blobs import leaf_blob_adjacency, split_blob_adjacency


@dataclass(frozen=True)
class ExperimentFiles:
    network: str
    gene_trees: str
    tob: str
    bipartitions: str
    predictions: str
    log: str
    report: str


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="dataset path relative to camus-dataset, e.g. n25/10")
    parser.add_argument("delta", type=float)
    parser.add_argument("gene_trees", nargs="?", default="iqtree_500.nwk")
    parser.add_argument("--binary", default=default_binary())
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--neg-mult", type=int, default=3)
    parser.add_argument("--min-neg", type=int, default=40)
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--eps", type=float, default=0.05,
                        help="epsilon for the theoretical row-sweep bound")
    parser.add_argument("--query-alpha", type=float, default=0.05,
                        help="per-query significance level for the T1 test")
    args = parser.parse_args()

    if not 0 <= args.query_alpha <= 1:
        parser.error("--query-alpha must be between 0 and 1")
    if not 0 < args.eps < 1:
        parser.error("--eps must be between 0 and 1")
    return args


def experiment_files(args):
    dataset_dir = os.path.join(camus_root(), args.path)
    output_dir = args.outdir or os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "results",
        args.path.replace("/", "_"),
    )
    os.makedirs(output_dir, exist_ok=True)

    tag = f"d{args.delta}_a{args.query_alpha}"
    return ExperimentFiles(
        network=os.path.join(dataset_dir, "true_net.nwk"),
        gene_trees=os.path.join(dataset_dir, args.gene_trees),
        tob=os.path.join(output_dir, "tob.nwk"),
        bipartitions=os.path.join(output_dir, "bips.tsv"),
        predictions=os.path.join(output_dir, f"pred_{tag}.tsv"),
        log=os.path.join(output_dir, f"run_{tag}.log"),
        report=os.path.join(output_dir, f"report_{tag}.txt"),
    )


def split_key(side_text):
    """Parse the canonical (non-reference) A side stored in the TSV."""
    return frozenset(taxon for taxon in side_text.split(",") if taxon)


def format_taxa(taxa):
    return "{" + ",".join(sorted(taxa, key=lambda taxon: (len(taxon), taxon))) + "}"


def smaller_side(split, all_leaves):
    complement = all_leaves - split
    return split if len(split) <= len(complement) else complement


def larger_side_size(split, leaf_count):
    return max(len(split), leaf_count - len(split))


def yes_no_unknown(value):
    if value is None:
        return "?"
    return "YES" if value else "no"


def numeric_row_id(row_id):
    return int(row_id[1:])


def nearest_true_split(candidate, true_splits):
    """Return the true split with minimum symmetric difference from candidate."""
    if not true_splits:
        return frozenset(), len(candidate)
    nearest = min(true_splits, key=lambda split: len(candidate ^ split))
    return nearest, len(candidate ^ nearest)


def theory_requirement(leaf_count, epsilon, delta):
    """Minimum s-1 in the proof; infinite outside its delta < 1/3 regime."""
    if delta >= 1 / 3:
        return float("inf")
    return 8 * math.log(leaf_count * leaf_count / epsilon) / (1 - 3 * delta) ** 2


def positive_side_sizes(labels, sides, leaf_count):
    sizes = []
    for row_id, label in labels.items():
        if label == 1:
            sizes.append(larger_side_size(split_key(sides[row_id][0]), leaf_count))
    return sizes


def append_false_negative_details(lines, result, sides, all_leaves, blob_adjacency):
    blob_adjacent_count = 0
    if not result.false_negative_ids:
        return blob_adjacent_count

    lines.append("-- FALSE NEGATIVES (rejected true splits) --")
    lines.append(f"  {'id':<5} {'|small|':>7} {'|large|=s':>9} {'blob-adj':>9}   smaller side")
    for row_id in sorted(result.false_negative_ids):
        split = split_key(sides[row_id][0])
        small_side = smaller_side(split, all_leaves)
        is_blob_adjacent = blob_adjacency.get(split)
        blob_adjacent_count += int(bool(is_blob_adjacent))
        lines.append(
            f"  {row_id:<5} {len(small_side):>7} "
            f"{larger_side_size(split, len(all_leaves)):>9} "
            f"{yes_no_unknown(is_blob_adjacent):>9}   {format_taxa(small_side)}"
        )
    lines.append("")
    return blob_adjacent_count


def append_false_positive_details(lines, result, sides, all_leaves, true_splits,
                                  leaf_adjacency, hard_negative_count):
    hard_count = 0
    moved_at_blob_count = 0
    if not result.false_positive_ids:
        return hard_count, moved_at_blob_count

    lines.append("-- FALSE POSITIVES (accepted non-splits) --")
    lines.append(
        f"  {'id':<5} {'size':>4} {'type':>6} {'symdiff':>7} "
        f"{'moved@blob':>10}   nearest true split"
    )
    for row_id in sorted(result.false_positive_ids, key=numeric_row_id):
        candidate = split_key(sides[row_id][0])
        candidate_small_side = smaller_side(candidate, all_leaves)
        is_hard = numeric_row_id(row_id) < hard_negative_count
        hard_count += int(is_hard)

        nearest, symmetric_difference = nearest_true_split(candidate, true_splits)
        moved_taxa = candidate ^ nearest
        moved_at_blob = any(leaf_adjacency.get(taxon, False) for taxon in moved_taxa)
        moved_at_blob_count += int(moved_at_blob)
        nearest_small_side = smaller_side(nearest, all_leaves)

        lines.append(
            f"  {row_id:<5} {len(candidate_small_side):>4} "
            f"{('hard' if is_hard else 'random'):>6} {symmetric_difference:>7} "
            f"{('YES' if moved_at_blob else 'no'):>10}   {format_taxa(nearest_small_side)}"
        )
    lines.append("")
    return hard_count, moved_at_blob_count


def build_report(args, all_leaves, true_splits, source, bipartition_stats,
                 labels, sides, result, blob_adjacency, leaf_adjacency):
    leaf_count = len(all_leaves)
    theta = (1 + args.delta) / 4
    required_s_minus_one = theory_requirement(leaf_count, args.eps, args.delta)
    larger_side_sizes = positive_side_sizes(labels, sides, leaf_count)
    minimum_s = min(larger_side_sizes) if larger_side_sizes else 0
    bound_satisfied = minimum_s - 1 >= required_s_minus_one
    bound_status = (
        "SATISFIED" if bound_satisfied
        else "NOT satisfied -> worst-case bound vacuous; relies on realized noise << delta"
    )

    lines = [
        "=" * 70,
        f"row_sweep_test_idx  |  {args.path}  |  delta = {args.delta}  "
        f"(theta = {theta:.4f})  |  T1 alpha = {args.query_alpha}",
        "=" * 70,
        f"leaves                 : {leaf_count}",
        f"gene trees             : {args.gene_trees}",
        f"tree-of-blobs source   : {source}",
        f"positives (true splits): {bipartition_stats.positive_count}",
        f"negatives (non-splits) : {bipartition_stats.negative_count}  "
        f"({bipartition_stats.hard_negative_count} hard + "
        f"{bipartition_stats.random_negative_count} random)",
        "",
        f"-- theory reference (Thm, eps={args.eps:.2f}) --",
        f"  larger-side size s ranges over positives; min s = {minimum_s}",
        f"  guarantee needs  s-1 >= {required_s_minus_one:.0f}   ({bound_status})",
        "",
        "-- results --",
        f"  FALSE NEGATIVES : {result.false_negative_count:>2}/{result.positive_count:<2}  "
        f"({result.false_negative_rate:6.1%})   (true split wrongly REJECTED)",
        f"  FALSE POSITIVES : {result.false_positive_count:>2}/{result.negative_count:<2}  "
        f"({result.false_positive_rate:6.1%})   (non-split wrongly ACCEPTED)",
        f"  accuracy        : {result.accuracy:6.1%}"
        + (f"   [{result.missing_count} predictions missing]" if result.missing_count else ""),
        "",
    ]

    false_negative_blob_count = append_false_negative_details(
        lines, result, sides, all_leaves, blob_adjacency
    )
    false_positive_hard_count, false_positive_moved_blob_count = append_false_positive_details(
        lines, result, sides, all_leaves, true_splits, leaf_adjacency,
        bipartition_stats.hard_negative_count,
    )

    lines.append("-- weakness summary --")
    if result.false_negative_count:
        small_false_negative_count = sum(
            len(smaller_side(split_key(sides[row_id][0]), all_leaves)) <= 3
            for row_id in result.false_negative_ids
        )
        lines.append(
            f"  * {false_negative_blob_count}/{result.false_negative_count} false negatives "
            "are blob-adjacent splits (true clades hanging off a reticulation cycle)."
        )
        lines.append(
            f"  * {small_false_negative_count}/{result.false_negative_count} false negatives "
            "have a small side <= 3 taxa."
        )
    else:
        lines.append("  * no false negatives at this delta.")

    if result.false_positive_count:
        lines.append(
            f"  * {false_positive_hard_count}/{result.false_positive_count} false positives "
            "are HARD negatives (small perturbations of a real split), not random ones."
        )
        lines.append(
            f"  * {false_positive_moved_blob_count}/{result.false_positive_count} false positives "
            "differ from the nearest true split by a blob-adjacent taxon."
        )
    else:
        lines.append("  * no false positives at this delta.")

    return "\n".join(lines)


def main():
    args = parse_args()
    files = experiment_files(args)
    for required_file in (files.network, files.gene_trees, args.binary):
        if not os.path.exists(required_file):
            raise SystemExit(f"ERROR: missing {required_file}")

    all_leaves, true_splits, source = compute_tob_splits(files.network, files.tob)
    with open(files.network) as network_file:
        network_text = network_file.read()
    blob_adjacency = split_blob_adjacency(network_text)[1]
    leaf_adjacency = leaf_blob_adjacency(network_text)

    bipartition_stats = build_bipartition_file(
        all_leaves, true_splits, files.bipartitions,
        seed=args.seed, neg_mult=args.neg_mult, min_neg=args.min_neg,
    )
    return_code = run_rowsweep(
        args.binary, files.gene_trees, files.bipartitions, files.predictions,
        args.delta, files.log, args.query_alpha,
    )
    if return_code != 0 or not os.path.exists(files.predictions):
        raise SystemExit(
            f"ERROR: tree-qmc failed (rc={return_code}). log tail:\n{log_tail(files.log)}"
        )

    labels, sides = read_labels(files.bipartitions)
    result = score(labels, read_predictions(files.predictions))
    report = build_report(
        args, all_leaves, true_splits, source, bipartition_stats,
        labels, sides, result, blob_adjacency, leaf_adjacency,
    )
    print(report)
    with open(files.report, "w") as report_file:
        report_file.write(report + "\n")


if __name__ == "__main__":
    main()
