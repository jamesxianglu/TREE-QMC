"""
Generate a labelled bipartition test set for one dataset's true network.

Positives : the non-trivial splits of the tree of blobs (should be ACCEPTed).
Negatives : a MIXED set of bipartitions that are NOT tree-of-blobs splits
            (should be REJECTed), combining
              - hard negatives  : true splits perturbed by moving one taxon across
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
from dataclasses import dataclass

from tree_of_blobs import tree_of_blobs_splits


@dataclass(frozen=True)
class BipartitionSetStats:
    leaf_count: int
    positive_count: int
    negative_count: int
    hard_negative_count: int
    random_negative_count: int


def taxon_sort_key(taxon):
    return len(taxon), taxon


def _canonical_side(side, all_leaves, reference_taxon):
    """Canonical key = the side that does not contain the reference taxon."""
    side = frozenset(side)
    return frozenset(all_leaves - side) if reference_taxon in side else side


def _format_taxa(taxa):
    return ",".join(sorted(taxa, key=taxon_sort_key))


def build_bipartition_file(all_leaves, splits, out_path, seed=0, neg_mult=3, min_neg=30):
    """Write the labelled bipartition TSV and return named row counts.

    `all_leaves` / `splits` are the taxon set and the set of true non-trivial
    tree-of-blobs splits (each a frozenset of the side without the reference
    taxon), supplied by the caller -- typically from PhyloNetworks' treeofblobs.
    """
    rng = random.Random(seed)
    all_leaves = frozenset(all_leaves)
    splits = set(splits)
    reference_taxon = "OUT" if "OUT" in all_leaves else min(all_leaves)
    leaf_count = len(all_leaves)
    ingroup = sorted(all_leaves - {reference_taxon}, key=taxon_sort_key)

    # ----- positives: the true tree-of-blobs splits -----
    positives = sorted(splits, key=lambda side: (len(side), sorted(side)))

    # ----- negatives -----
    target_negative_count = max(neg_mult * len(positives), min_neg)
    used_splits = set(splits)                # never emit a duplicate or true split
    hard_negatives = []
    random_negatives = []

    # hard negatives: perturb each positive by moving 1 taxon in either direction
    positive_pool = list(positives)
    rng.shuffle(positive_pool)
    for positive_side in positive_pool:
        if len(hard_negatives) >= target_negative_count // 2:
            break
        positive_side = set(positive_side)
        complement = set(all_leaves) - positive_side
        perturbed_sides = []

        # Move one taxon out, while leaving at least two taxa on this side.
        if len(positive_side) >= 3:
            moved_taxon = rng.choice(sorted(positive_side))
            perturbed_sides.append(positive_side - {moved_taxon})

        # Move one non-reference taxon in, while leaving two in the complement.
        movable_taxa = sorted(complement - {reference_taxon})
        if movable_taxa and len(complement) >= 3:
            moved_taxon = rng.choice(movable_taxa)
            perturbed_sides.append(positive_side | {moved_taxon})

        for perturbed_side in perturbed_sides:
            if not (2 <= len(perturbed_side) <= leaf_count - 2):
                continue
            key = _canonical_side(perturbed_side, all_leaves, reference_taxon)
            if key in used_splits:
                continue
            used_splits.add(key)
            hard_negatives.append(key)

    # random negatives: random subsets of the ingroup, both sides >= 2
    attempts = 0
    while (len(hard_negatives) + len(random_negatives) < target_negative_count
           and attempts < 20000):
        attempts += 1
        side_size = rng.randint(2, leaf_count - 2)
        subset = frozenset(rng.sample(ingroup, side_size))
        if not (2 <= len(subset) <= leaf_count - 2):
            continue
        key = _canonical_side(subset, all_leaves, reference_taxon)
        if key in used_splits:
            continue
        used_splits.add(key)
        random_negatives.append(key)

    # ----- write -----
    rows = []
    for index, key in enumerate(positives):
        side = set(key)
        rows.append((f"p{index}", side, all_leaves - side, 1))
    negatives = hard_negatives + random_negatives
    for index, key in enumerate(negatives):
        side = set(key)
        rows.append((f"n{index}", side, all_leaves - side, 0))

    with open(out_path, "w") as f:
        f.write("id\tA\tB\tlabel\n")
        for row_id, side_a, side_b, label in rows:
            f.write(
                f"{row_id}\t{_format_taxa(side_a)}\t{_format_taxa(side_b)}\t{label}\n"
            )

    return BipartitionSetStats(
        leaf_count=leaf_count,
        positive_count=len(positives),
        negative_count=len(negatives),
        hard_negative_count=len(hard_negatives),
        random_negative_count=len(random_negatives),
    )


if __name__ == "__main__":
    import os
    import sys

    dataset = sys.argv[1] if len(sys.argv) > 1 else "00"
    net = os.path.join(os.path.dirname(__file__), "..", "..",
                       "data", "camus-dataset", "n15", dataset, "true_net.nwk")
    with open(net) as network_file:
        all_leaves, splits = tree_of_blobs_splits(network_file.read())
    out = f"/tmp/bips_{dataset}.tsv"
    stats = build_bipartition_file(all_leaves, splits, out, seed=1)
    print(stats)
    with open(out) as output_file:
        print(output_file.read())
