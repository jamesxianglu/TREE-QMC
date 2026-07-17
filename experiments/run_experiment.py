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
                              [--deltas 0.05,...] [--seed S] [--outdir DIR]
"""

import argparse
import csv
import os

from common import (camus_root, default_binary, compute_tob_splits, run_rowsweep,
                    read_labels, read_predictions, score)
from gen_bipartitions import build_bipartition_file

DEFAULT_DELTAS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
HERE = os.path.dirname(os.path.abspath(__file__))


def select_datasets(group_dir, n):
    ds = []
    for d in sorted(os.listdir(group_dir)):
        p = os.path.join(group_dir, d)
        if os.path.isdir(p) and os.path.exists(os.path.join(p, "true_net.nwk")):
            ds.append(d)
        if len(ds) >= n:
            break
    return ds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=default_binary())
    ap.add_argument("--group", default="n15")
    ap.add_argument("--datasets", type=int, default=10)
    ap.add_argument("--gene-trees", default="iqtree_500.nwk")
    ap.add_argument("--deltas", default=",".join(str(d) for d in DEFAULT_DELTAS))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--neg-mult", type=int, default=3)
    ap.add_argument("--min-neg", type=int, default=30)
    ap.add_argument("--outdir", default=os.path.join(HERE, "results", "sweep"))
    args = ap.parse_args()

    if not os.path.exists(args.binary):
        raise SystemExit(f"ERROR: tree-qmc binary not found at {args.binary}")

    group_dir = os.path.join(camus_root(), args.group)
    deltas = [float(x) for x in args.deltas.split(",")]
    datasets = select_datasets(group_dir, args.datasets)
    os.makedirs(args.outdir, exist_ok=True)
    for sub in ("bips", "preds", "logs", "tob"):
        os.makedirs(os.path.join(args.outdir, sub), exist_ok=True)

    print(f"Group: {args.group}   Datasets: {datasets}")
    print(f"Deltas: {deltas}   Gene trees: {args.gene_trees}\n")

    labels_by_ds = {}
    for d in datasets:
        net = os.path.join(group_dir, d, "true_net.nwk")
        tob = os.path.join(args.outdir, "tob", f"{d}.nwk")
        bip = os.path.join(args.outdir, "bips", f"{d}.tsv")
        all_leaves, splits, source = compute_tob_splits(net, tob)
        stats = build_bipartition_file(all_leaves, splits, bip,
                                       seed=args.seed + int(d),
                                       neg_mult=args.neg_mult, min_neg=args.min_neg)
        labels_by_ds[d] = read_labels(bip)[0]
        print(f"  {args.group}/{d}: {stats['n_positive']} pos, {stats['n_negative']} neg "
              f"({stats['n_hard_neg']}h+{stats['n_rand_neg']}r)  [{source.split(' ')[0]}]")

    rows = []
    for d in datasets:
        gt = os.path.join(group_dir, d, args.gene_trees)
        if not os.path.exists(gt):
            print(f"  WARN: {gt} missing; skipping {d}")
            continue
        bip = os.path.join(args.outdir, "bips", f"{d}.tsv")
        for delta in deltas:
            out = os.path.join(args.outdir, "preds", f"{d}_d{delta}.tsv")
            log = os.path.join(args.outdir, "logs", f"{d}_d{delta}.log")
            rc, cmd = run_rowsweep(args.binary, gt, bip, out, delta, log)
            if rc != 0 or not os.path.exists(out):
                tail = "".join(open(log).readlines()[-20:]) if os.path.exists(log) else ""
                raise SystemExit(f"ERROR: tree-qmc failed (rc={rc}) {d} d={delta}\n{tail}")
            s = score(labels_by_ds[d], read_predictions(out))
            fnr = s["fn"] / s["n_pos"] if s["n_pos"] else float("nan")
            fpr = s["fp"] / s["n_neg"] if s["n_neg"] else float("nan")
            rows.append(dict(dataset=d, delta=delta, n_pos=s["n_pos"], n_neg=s["n_neg"],
                            fn=s["fn"], fp=s["fp"], missing=s["missing"],
                            fn_rate=fnr, fp_rate=fpr))
            print(f"  {args.group}/{d} d={delta:<4}: FN {s['fn']}/{s['n_pos']} ({fnr:.1%})  "
                  f"FP {s['fp']}/{s['n_neg']} ({fpr:.1%})")

    with open(os.path.join(args.outdir, "detail.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dataset", "delta", "n_pos", "n_neg",
                                          "fn", "fp", "missing", "fn_rate", "fp_rate"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("\n=== POOLED across datasets (per delta) ===")
    print(f"{'delta':>6} | {'FN/pos':>10} {'FNR':>7} | {'FP/neg':>10} {'FPR':>7}")
    print("-" * 50)
    summary = []
    for delta in deltas:
        rs = [r for r in rows if r["delta"] == delta]
        if not rs:
            continue
        FN = sum(r["fn"] for r in rs); P = sum(r["n_pos"] for r in rs)
        FP = sum(r["fp"] for r in rs); Nn = sum(r["n_neg"] for r in rs)
        fnr = FN / P if P else float("nan")
        fpr = FP / Nn if Nn else float("nan")
        summary.append(dict(delta=delta, FN=FN, n_pos=P, FP=FP, n_neg=Nn, fn_rate=fnr, fp_rate=fpr))
        print(f"{delta:>6} | {f'{FN}/{P}':>10} {fnr:>6.1%} | {f'{FP}/{Nn}':>10} {fpr:>6.1%}")

    with open(os.path.join(args.outdir, "summary_by_delta.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["delta", "FN", "n_pos", "FP", "n_neg", "fn_rate", "fp_rate"])
        w.writeheader()
        for r in summary:
            w.writerow(r)
    print(f"\nWrote results under {args.outdir}")


if __name__ == "__main__":
    main()
