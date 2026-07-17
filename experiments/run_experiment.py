#!/usr/bin/env python3
"""
Sweep row_sweep_test_idx (via tree-qmc --rowsweep) over many datasets and deltas.

Ground truth = PhyloNetworks' tree of blobs (bridge fallback if Julia fails).
For each dataset: positives = true non-trivial ToB splits; negatives = mixed
hard + random non-splits. Scoring: a positive predicted REJECT is a false
negative; a negative predicted ACCEPT is a false positive. Reports pooled
FNR/FPR per delta plus a per-(dataset, delta) CSV.

Usage:
    python3 run_experiment.py [--group n15] [--datasets N] [--gene-trees FILE]
                              [--deltas 0.05,...] [--query-alpha A]
                              [--seed S] [--outdir DIR]
"""

import argparse
import csv
import os
from dataclasses import dataclass

from common import (camus_root, compute_tob_splits, default_binary, log_tail,
                    read_labels, read_predictions, run_rowsweep, score)
from gen_bipartitions import build_bipartition_file

DEFAULT_DELTAS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
HERE = os.path.dirname(os.path.abspath(__file__))


@dataclass(frozen=True)
class PreparedDataset:
    name: str
    gene_trees: str
    bipartitions: str
    labels: dict


def select_datasets(group_dir, limit):
    datasets = []
    for name in sorted(os.listdir(group_dir)):
        dataset_dir = os.path.join(group_dir, name)
        network = os.path.join(dataset_dir, "true_net.nwk")
        if os.path.isdir(dataset_dir) and os.path.exists(network):
            datasets.append(name)
        if len(datasets) >= limit:
            break
    return datasets


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default=default_binary())
    parser.add_argument("--group", default="n15")
    parser.add_argument("--datasets", type=int, default=10)
    parser.add_argument("--gene-trees", default="iqtree_500.nwk")
    parser.add_argument("--deltas", default=",".join(str(d) for d in DEFAULT_DELTAS))
    parser.add_argument("--query-alpha", type=float, default=0.05,
                        help="per-query significance level for the T1 test")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--neg-mult", type=int, default=3)
    parser.add_argument("--min-neg", type=int, default=30)
    parser.add_argument("--outdir", default=os.path.join(HERE, "results", "sweep"))
    args = parser.parse_args()
    if not 0 <= args.query_alpha <= 1:
        parser.error("--query-alpha must be between 0 and 1")
    return args


def ensure_output_directories(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for subdirectory in ("bips", "preds", "logs", "tob"):
        os.makedirs(os.path.join(output_dir, subdirectory), exist_ok=True)


def prepare_datasets(args, group_dir, dataset_names):
    prepared = []
    for name in dataset_names:
        network = os.path.join(group_dir, name, "true_net.nwk")
        gene_trees = os.path.join(group_dir, name, args.gene_trees)
        tob = os.path.join(args.outdir, "tob", f"{name}.nwk")
        bipartitions = os.path.join(args.outdir, "bips", f"{name}.tsv")

        all_leaves, true_splits, source = compute_tob_splits(network, tob)
        stats = build_bipartition_file(
            all_leaves, true_splits, bipartitions,
            seed=args.seed + int(name),
            neg_mult=args.neg_mult,
            min_neg=args.min_neg,
        )
        labels = read_labels(bipartitions)[0]
        prepared.append(PreparedDataset(name, gene_trees, bipartitions, labels))
        print(
            f"  {args.group}/{name}: {stats.positive_count} pos, "
            f"{stats.negative_count} neg "
            f"({stats.hard_negative_count}h+{stats.random_negative_count}r)  "
            f"[{source.split(' ')[0]}]"
        )
    return prepared


def evaluate_datasets(args, datasets, deltas):
    results = []
    for dataset in datasets:
        if not os.path.exists(dataset.gene_trees):
            print(f"  WARN: {dataset.gene_trees} missing; skipping {dataset.name}")
            continue

        for delta in deltas:
            run_tag = f"d{delta}_a{args.query_alpha}"
            prediction_file = os.path.join(
                args.outdir, "preds", f"{dataset.name}_{run_tag}.tsv"
            )
            log_file = os.path.join(args.outdir, "logs", f"{dataset.name}_{run_tag}.log")
            return_code = run_rowsweep(
                args.binary, dataset.gene_trees, dataset.bipartitions,
                prediction_file, delta, log_file, args.query_alpha,
            )
            if return_code != 0 or not os.path.exists(prediction_file):
                raise SystemExit(
                    f"ERROR: tree-qmc failed (rc={return_code}) "
                    f"{dataset.name} d={delta}\n{log_tail(log_file)}"
                )

            result = score(dataset.labels, read_predictions(prediction_file))
            results.append({
                "dataset": dataset.name,
                "delta": delta,
                "query_alpha": args.query_alpha,
                "n_pos": result.positive_count,
                "n_neg": result.negative_count,
                "fn": result.false_negative_count,
                "fp": result.false_positive_count,
                "missing": result.missing_count,
                "fn_rate": result.false_negative_rate,
                "fp_rate": result.false_positive_rate,
            })
            print(
                f"  {args.group}/{dataset.name} d={delta:<4}: "
                f"FN {result.false_negative_count}/{result.positive_count} "
                f"({result.false_negative_rate:.1%})  "
                f"FP {result.false_positive_count}/{result.negative_count} "
                f"({result.false_positive_rate:.1%})"
            )
    return results


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize_by_delta(results, deltas):
    print("\n=== POOLED across datasets (per delta) ===")
    print(f"{'delta':>6} | {'FN/pos':>10} {'FNR':>7} | {'FP/neg':>10} {'FPR':>7}")
    print("-" * 50)

    summary = []
    for delta in deltas:
        matching_results = [row for row in results if row["delta"] == delta]
        if not matching_results:
            continue

        false_negatives = sum(row["fn"] for row in matching_results)
        positives = sum(row["n_pos"] for row in matching_results)
        false_positives = sum(row["fp"] for row in matching_results)
        negatives = sum(row["n_neg"] for row in matching_results)
        false_negative_rate = false_negatives / positives if positives else float("nan")
        false_positive_rate = false_positives / negatives if negatives else float("nan")
        summary.append({
            "delta": delta,
            "query_alpha": matching_results[0]["query_alpha"],
            "FN": false_negatives,
            "n_pos": positives,
            "FP": false_positives,
            "n_neg": negatives,
            "fn_rate": false_negative_rate,
            "fp_rate": false_positive_rate,
        })
        print(
            f"{delta:>6} | {f'{false_negatives}/{positives}':>10} "
            f"{false_negative_rate:>6.1%} | {f'{false_positives}/{negatives}':>10} "
            f"{false_positive_rate:>6.1%}"
        )
    return summary


def main():
    args = parse_args()

    if not os.path.exists(args.binary):
        raise SystemExit(f"ERROR: tree-qmc binary not found at {args.binary}")

    group_dir = os.path.join(camus_root(), args.group)
    deltas = [float(x) for x in args.deltas.split(",")]
    dataset_names = select_datasets(group_dir, args.datasets)
    ensure_output_directories(args.outdir)

    print(f"Group: {args.group}   Datasets: {dataset_names}")
    print(f"Deltas: {deltas}   T1 alpha: {args.query_alpha}   "
          f"Gene trees: {args.gene_trees}\n")

    datasets = prepare_datasets(args, group_dir, dataset_names)
    results = evaluate_datasets(args, datasets, deltas)
    detail_fields = [
        "dataset", "delta", "query_alpha", "n_pos", "n_neg", "fn", "fp", "missing",
        "fn_rate", "fp_rate",
    ]
    alpha_tag = f"a{args.query_alpha}"
    write_csv(os.path.join(args.outdir, f"detail_{alpha_tag}.csv"), results, detail_fields)

    summary = summarize_by_delta(results, deltas)
    summary_fields = [
        "delta", "query_alpha", "FN", "n_pos", "FP", "n_neg", "fn_rate", "fp_rate",
    ]
    write_csv(
        os.path.join(args.outdir, f"summary_by_delta_{alpha_tag}.csv"),
        summary,
        summary_fields,
    )
    print(f"\nWrote results under {args.outdir}")


if __name__ == "__main__":
    main()
