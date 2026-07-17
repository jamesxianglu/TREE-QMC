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

from common import (camus_root, default_binary, compute_tob_splits, run_rowsweep,
                    read_labels, read_predictions, score)
from tree_of_blobs import split_blob_adjacency, leaf_blob_adjacency
from gen_bipartitions import build_bipartition_file


def key_of(sideA):
    """Canonical split key = the (non-OUT) A side, as stored in the bip file."""
    return frozenset(t for t in sideA.split(",") if t)


def taxa(s):
    return "{" + ",".join(sorted(s, key=lambda x: (len(x), x))) + "}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="dataset path relative to camus-dataset, e.g. n25/10")
    ap.add_argument("delta", type=float)
    ap.add_argument("gene_trees", nargs="?", default="iqtree_500.nwk")
    ap.add_argument("--binary", default=default_binary())
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--neg-mult", type=int, default=3)
    ap.add_argument("--min-neg", type=int, default=40)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--eps", type=float, default=0.05, help="epsilon for the theory bound")
    args = ap.parse_args()

    root = camus_root()
    ds_dir = os.path.join(root, args.path)
    net = os.path.join(ds_dir, "true_net.nwk")
    gt = os.path.join(ds_dir, args.gene_trees)
    for f in (net, gt):
        if not os.path.exists(f):
            raise SystemExit(f"ERROR: missing {f}")
    if not os.path.exists(args.binary):
        raise SystemExit(f"ERROR: tree-qmc binary not found at {args.binary}")

    outdir = args.outdir or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "results", args.path.replace("/", "_"))
    os.makedirs(outdir, exist_ok=True)
    tob_file = os.path.join(outdir, "tob.nwk")
    bip_file = os.path.join(outdir, "bips.tsv")
    out_file = os.path.join(outdir, f"pred_d{args.delta}.tsv")
    log_file = os.path.join(outdir, f"run_d{args.delta}.log")

    # 1) ground truth
    all_leaves, splits, source = compute_tob_splits(net, tob_file)
    blob_adj = split_blob_adjacency(open(net).read())[1]      # {split_key: bool}
    leaf_adj = leaf_blob_adjacency(open(net).read())          # {taxon: bool}
    n = len(all_leaves)

    # 2) bipartitions + run
    stats = build_bipartition_file(all_leaves, splits, bip_file, seed=args.seed,
                                   neg_mult=args.neg_mult, min_neg=args.min_neg)
    rc, cmd = run_rowsweep(args.binary, gt, bip_file, out_file, args.delta, log_file)
    if rc != 0 or not os.path.exists(out_file):
        tail = "".join(open(log_file).readlines()[-20:]) if os.path.exists(log_file) else ""
        raise SystemExit(f"ERROR: tree-qmc failed (rc={rc}). log tail:\n{tail}")

    labels, sides = read_labels(bip_file)
    preds = read_predictions(out_file)
    s = score(labels, preds)

    # theory reference: s-1 >= 8 ln(n^2/eps) / (1-3 delta)^2
    theta = (1 + args.delta) / 4
    if args.delta < 1 / 3:
        req = 8 * math.log(n * n / args.eps) / (1 - 3 * args.delta) ** 2
    else:
        req = float("inf")
    pos_larger = [max(len(key_of(sides[r][0])), n - len(key_of(sides[r][0])))
                  for r, l in labels.items() if l == 1]
    min_s = min(pos_larger) if pos_larger else 0

    # ---- report ----
    L = []
    L.append("=" * 70)
    L.append(f"row_sweep_test_idx  |  {args.path}  |  delta = {args.delta}  (theta = {theta:.4f})")
    L.append("=" * 70)
    L.append(f"leaves                 : {n}")
    L.append(f"gene trees             : {args.gene_trees}")
    L.append(f"tree-of-blobs source   : {source}")
    L.append(f"positives (true splits): {stats['n_positive']}")
    L.append(f"negatives (non-splits) : {stats['n_negative']}  "
             f"({stats['n_hard_neg']} hard + {stats['n_rand_neg']} random)")
    L.append("")
    L.append("-- theory reference (Thm, eps=%.2f) --" % args.eps)
    L.append(f"  larger-side size s ranges over positives; min s = {min_s}")
    L.append(f"  guarantee needs  s-1 >= {req:.0f}   "
             + ("(SATISFIED)" if min_s - 1 >= req else "(NOT satisfied -> worst-case bound vacuous; "
                                                       "relies on realized noise << delta)"))
    L.append("")
    acc = (s['n_pos'] + s['n_neg'] - s['fn'] - s['fp']) / max(1, s['n_pos'] + s['n_neg'])
    fnr = s['fn'] / s['n_pos'] if s['n_pos'] else float('nan')
    fpr = s['fp'] / s['n_neg'] if s['n_neg'] else float('nan')
    L.append("-- results --")
    L.append(f"  FALSE NEGATIVES : {s['fn']:>2}/{s['n_pos']:<2}  ({fnr:6.1%})   "
             "(true split wrongly REJECTED)")
    L.append(f"  FALSE POSITIVES : {s['fp']:>2}/{s['n_neg']:<2}  ({fpr:6.1%})   "
             "(non-split wrongly ACCEPTED)")
    L.append(f"  accuracy        : {acc:6.1%}"
             + (f"   [{s['missing']} predictions missing]" if s['missing'] else ""))
    L.append("")

    # FN detail
    fn_blob = 0
    if s['fn_ids']:
        L.append("-- FALSE NEGATIVES (rejected true splits) --")
        L.append(f"  {'id':<5} {'|small|':>7} {'|large|=s':>9} {'blob-adj':>9}   smaller side")
        for rid in sorted(s['fn_ids']):
            k = key_of(sides[rid][0])
            small = k if len(k) <= n - len(k) else all_leaves - k
            ba = blob_adj.get(k, None)
            fn_blob += 1 if ba else 0
            L.append(f"  {rid:<5} {len(small):>7} {max(len(k), n-len(k)):>9} "
                     f"{('YES' if ba else ('no' if ba is not None else '?')):>9}   {taxa(small)}")
        L.append("")

    # FP detail
    fp_hard = fp_movedblob = 0
    if s['fp_ids']:
        L.append("-- FALSE POSITIVES (accepted non-splits) --")
        L.append(f"  {'id':<5} {'size':>4} {'type':>6} {'symdiff':>7} {'moved@blob':>10}   nearest true split")
        for rid in sorted(s['fp_ids'], key=lambda r: int(r[1:])):
            k = key_of(sides[rid][0])
            small = k if len(k) <= n - len(k) else all_leaves - k
            is_hard = int(rid[1:]) < stats['n_hard_neg']
            fp_hard += 1 if is_hard else 0
            # nearest true split by symmetric-difference size
            best, bestd = None, 10 ** 9
            for sp in splits:
                d = len(k ^ sp)
                if d < bestd:
                    best, bestd = sp, d
            moved = (k ^ best) if best is not None else set()
            moved_at_blob = any(leaf_adj.get(t, False) for t in moved)
            fp_movedblob += 1 if moved_at_blob else 0
            nb = best if (best is None or len(best) <= n - len(best)) else all_leaves - best
            L.append(f"  {rid:<5} {len(small):>4} {('hard' if is_hard else 'random'):>6} "
                     f"{bestd:>7} {('YES' if moved_at_blob else 'no'):>10}   {taxa(nb)}")
        L.append("")

    # weakness summary
    L.append("-- weakness summary --")
    if s['fn']:
        small_fn = sum(1 for rid in s['fn_ids']
                       if min(len(key_of(sides[rid][0])), n - len(key_of(sides[rid][0]))) <= 3)
        L.append(f"  * {fn_blob}/{s['fn']} false negatives are blob-adjacent splits "
                 "(true clades hanging off a reticulation cycle).")
        L.append(f"  * {small_fn}/{s['fn']} false negatives have a small side <= 3 taxa.")
    else:
        L.append("  * no false negatives at this delta.")
    if s['fp']:
        L.append(f"  * {fp_hard}/{s['fp']} false positives are HARD negatives "
                 "(small perturbations of a real split), not random ones.")
        L.append(f"  * {fp_movedblob}/{s['fp']} false positives differ from the nearest true "
                 "split by a blob-adjacent taxon.")
    else:
        L.append("  * no false positives at this delta.")

    report = "\n".join(L)
    print(report)
    with open(os.path.join(outdir, f"report_d{args.delta}.txt"), "w") as f:
        f.write(report + "\n")


if __name__ == "__main__":
    main()
