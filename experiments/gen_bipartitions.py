"""
Generate a labelled bipartition test set for one dataset's true network.

Positives : the non-trivial splits of the tree of blobs (should be ACCEPTed).
Negatives : a MIXED set of bipartitions that are NOT tree-of-blobs splits
            (should be REJECTed), combining
              - hard negatives  : true splits perturbed by moving 1-2 taxa across
                                  the split (structurally "close" to a real split)
              - random negatives : random bipartitions with both sides >= 2 taxa

Output: a TSV with header `id  A  B  label` where
   id    : unique id (p* = positive, n* = negative)
   A, B  : comma-separated taxon names for the two sides
   label : 1 if the bipartition is a true ToB split, 0 otherwise
tree-qmc's --rowsweep reader consumes only the first three columns (id, A, B)
and ignores the trailing `label`, so the same file drives the tool and scores it.
"""

import random

from tree_of_blobs import tree_of_blobs_splits


def _canon(side, all_leaves, ref):
    """Canonical key for a bipartition = the side that does not contain `ref`."""
    side = frozenset(side)
    return frozenset(all_leaves - side) if ref in side else side


def build_bipartition_file(all_leaves, splits, out_path, seed=0, neg_mult=3, min_neg=30):
    """Write the labelled bipartition TSV for one network. Returns a stats dict.

    `all_leaves` / `splits` are the taxon set and the set of true non-trivial
    tree-of-blobs splits (each a frozenset of the side without the reference
    taxon), supplied by the caller -- typically from PhyloNetworks' treeofblobs.
    """
    rng = random.Random(seed)
    all_leaves = frozenset(all_leaves)
    splits = set(splits)
    ref = "OUT" if "OUT" in all_leaves else min(all_leaves)
    n = len(all_leaves)
    ingroup = sorted(all_leaves - {ref}, key=lambda x: (len(x), x))

    # ----- positives: the true tree-of-blobs splits -----
    positives = sorted(splits, key=lambda s: (len(s), sorted(s)))

    # ----- negatives -----
    target_neg = max(neg_mult * len(positives), min_neg)
    used = set(splits)                       # canonical keys we must not reproduce
    hard, rand = [], []

    # hard negatives: perturb each positive by moving 1 taxon in either direction
    pos_pool = list(positives)
    rng.shuffle(pos_pool)
    for S in pos_pool:
        if len(hard) >= target_neg // 2:
            break
        S = set(S)
        C = set(all_leaves) - S              # complement (contains ref)
        moves = []
        # move one taxon OUT of S  (needs |S| >= 3 so the new side keeps >= 2)
        if len(S) >= 3:
            x = rng.choice(sorted(S))
            moves.append(S - {x})
        # move one non-ref taxon INTO S from the complement
        cand = sorted(C - {ref})
        if cand and len(C) >= 3:             # complement keeps >= 2 after removal
            y = rng.choice(cand)
            moves.append(S | {y})
        for newS in moves:
            if not (2 <= len(newS) <= n - 2):
                continue
            key = _canon(newS, all_leaves, ref)
            if key in used:
                continue
            used.add(key)
            hard.append(key)

    # random negatives: random subsets of the ingroup, both sides >= 2
    attempts = 0
    while len(hard) + len(rand) < target_neg and attempts < 20000:
        attempts += 1
        k = rng.randint(2, n - 2)            # size of the non-ref side
        subset = frozenset(rng.sample(ingroup, min(k, len(ingroup))))
        if not (2 <= len(subset) <= n - 2):
            continue
        key = _canon(subset, all_leaves, ref)
        if key in used:
            continue
        used.add(key)
        rand.append(key)

    # ----- write -----
    rows = []
    for i, key in enumerate(positives):
        side = set(key)
        rows.append((f"p{i}", side, all_leaves - side, 1))
    negatives = hard + rand
    for i, key in enumerate(negatives):
        side = set(key)
        rows.append((f"n{i}", side, all_leaves - side, 0))

    def fmt(taxa):
        return ",".join(sorted(taxa, key=lambda x: (len(x), x)))

    with open(out_path, "w") as f:
        f.write("id\tA\tB\tlabel\n")
        for rid, A, B, lab in rows:
            f.write(f"{rid}\t{fmt(A)}\t{fmt(B)}\t{lab}\n")

    return {
        "n_leaves": n,
        "n_positive": len(positives),
        "n_negative": len(negatives),
        "n_hard_neg": len(hard),
        "n_rand_neg": len(rand),
    }


if __name__ == "__main__":
    import sys, os
    d = sys.argv[1] if len(sys.argv) > 1 else "00"
    net = os.path.join(os.path.dirname(__file__), "..", "..",
                       "data", "camus-dataset", "n15", d, "true_net.nwk")
    all_leaves, splits = tree_of_blobs_splits(open(net).read())
    out = f"/tmp/bips_{d}.tsv"
    stats = build_bipartition_file(all_leaves, splits, out, seed=1)
    print(stats)
    print(open(out).read())
