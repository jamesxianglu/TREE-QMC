"""
Compute the non-trivial splits of the TREE OF BLOBS of a phylogenetic network
given in extended Newick format (with #H reticulation nodes).

Background
----------
The tree of blobs T(N) of a network N is obtained by contracting every blob
(maximal biconnected subgraph, i.e. every 2-edge-connected component that
contains a cycle) of N to a single node.  The non-trivial splits of T(N) are
in one-to-one correspondence with the CUT EDGES (bridges) of the underlying
undirected graph of N: removing a bridge splits the leaf set into two groups.
A split is "non-trivial" when both sides contain >= 2 leaves.

These bridge-induced bipartitions are exactly the bipartitions that
`row_sweep_test_idx` should ACCEPT; every other bipartition (both sides >= 2)
should be REJECTed.

This module only parses the network and finds bridges -- no external deps.
"""

import re
import sys
from collections import defaultdict


class NetworkParseError(Exception):
    pass


def parse_extended_newick(nwk):
    """Parse an extended-Newick string into an undirected multigraph.

    Returns (edges, leaf_names) where
      edges      : list of (u, v) integer node-id pairs (one per graph edge)
      leaf_names : dict node_id -> taxon name, for leaf nodes only
    Reticulation (#H) nodes that appear twice are merged into a single node so
    that both parent edges connect to the same graph vertex (creating a cycle).
    """
    s = nwk.strip()
    if s.endswith(";"):
        s = s[:-1]

    nodes = {}                 # id -> {"name": str|None, "hybrid": str|None, "children": [ids]}
    hybrid_map = {}            # hybrid tag -> node id
    counter = [0]

    def new_node():
        nid = counter[0]
        counter[0] += 1
        nodes[nid] = {"name": None, "hybrid": None, "children": []}
        return nid

    i = [0]
    n = len(s)

    def parse_subtree():
        children = []
        if s[i[0]] == "(":
            i[0] += 1
            while True:
                children.append(parse_subtree())
                c = s[i[0]]
                if c == ",":
                    i[0] += 1
                    continue
                if c == ")":
                    i[0] += 1
                    break
                raise NetworkParseError(f"expected ',' or ')' at pos {i[0]}: ...{s[i[0]:i[0]+20]}")

        # Read the label token: everything up to ',', ')', ':', or end.
        start = i[0]
        while i[0] < n and s[i[0]] not in ",():":
            i[0] += 1
        label = s[start:i[0]]

        # Skip branch-length / support / gamma annotations ":len:sup:gamma".
        if i[0] < n and s[i[0]] == ":":
            while i[0] < n and s[i[0]] not in ",()":
                i[0] += 1

        # Split label into taxon-name part and optional hybrid tag.
        name = None
        hybrid = None
        if "#" in label:
            base, tag = label.split("#", 1)
            name = base if base else None
            hybrid = tag
        else:
            name = label if label else None

        if hybrid is not None:
            if hybrid in hybrid_map:
                nid = hybrid_map[hybrid]          # merge with the earlier occurrence
            else:
                nid = new_node()
                hybrid_map[hybrid] = nid
                nodes[nid]["hybrid"] = hybrid
            if name and nodes[nid]["name"] is None:
                nodes[nid]["name"] = name
        else:
            nid = new_node()
            nodes[nid]["name"] = name

        # Attach the children discovered in THIS occurrence (the bare reference
        # occurrence of a hybrid has no children, giving it a 2nd parent edge).
        for ch in children:
            nodes[nid]["children"].append(ch)
        return nid

    parse_subtree()
    if i[0] != n:
        raise NetworkParseError(f"trailing characters after root at pos {i[0]}: {s[i[0]:]}")

    # Build undirected edge list from parent->child relations.
    edges = []
    for pid, node in nodes.items():
        for ch in node["children"]:
            edges.append((pid, ch))

    # Leaves = named nodes with no children.
    leaf_names = {
        nid: node["name"]
        for nid, node in nodes.items()
        if node["name"] is not None and len(node["children"]) == 0
    }
    return edges, leaf_names


def find_bridges(num_nodes, edges):
    """Return the set of edge indices that are bridges (cut edges).

    Uses the standard low-link DFS.  Parallel edges are handled by skipping the
    specific edge index we arrived on, not merely the parent node.
    """
    adj = defaultdict(list)                 # node -> list of (neighbor, edge_index)
    for idx, (u, v) in enumerate(edges):
        adj[u].append((v, idx))
        adj[v].append((u, idx))

    disc = [-1] * num_nodes
    low = [0] * num_nodes
    bridges = set()
    timer = [0]

    sys.setrecursionlimit(10000)

    def dfs(u, in_edge):
        disc[u] = low[u] = timer[0]
        timer[0] += 1
        for v, eidx in adj[u]:
            if eidx == in_edge:
                continue
            if disc[v] == -1:
                dfs(v, eidx)
                low[u] = min(low[u], low[v])
                if low[v] > disc[u]:
                    bridges.add(eidx)
            else:
                low[u] = min(low[u], disc[v])

    for start in range(num_nodes):
        if disc[start] == -1:
            dfs(start, -1)
    return bridges


def _reachable_leaves(edges, exclude_edge_idx, source, leaf_ids):
    """Leaves reachable from `source` in the graph with one edge removed."""
    adj = defaultdict(list)
    for idx, (u, v) in enumerate(edges):
        if idx == exclude_edge_idx:
            continue
        adj[u].append(v)
        adj[v].append(u)
    seen = {source}
    stack = [source]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return {leaf_ids[nid] for nid in seen if nid in leaf_ids}


def tree_of_blobs_splits(nwk):
    """Compute the network's leaf set and its non-trivial tree-of-blobs splits.

    Returns (all_leaves, splits) where
      all_leaves : frozenset of every taxon name
      splits     : set of frozensets.  Each frozenset is one SIDE of a
                   non-trivial split, canonicalised as the side that does NOT
                   contain the outgroup 'OUT' (or, if 'OUT' is absent, the
                   lexicographically smaller side).  Both sides have >= 2 taxa.
    """
    edges, leaf_names = parse_extended_newick(nwk)
    num_nodes = 1 + max((max(u, v) for u, v in edges), default=-1)
    all_leaves = frozenset(leaf_names.values())
    ref = "OUT" if "OUT" in all_leaves else min(all_leaves)

    bridges = find_bridges(num_nodes, edges)
    splits = set()
    for eidx in bridges:
        u, v = edges[eidx]
        side = _reachable_leaves(edges, eidx, v, leaf_names)
        other = all_leaves - side
        if len(side) < 2 or len(other) < 2:
            continue                          # trivial split
        key = frozenset(other if ref in side else side)
        splits.add(key)
    return all_leaves, splits


def splits_from_tob_tree(tob_nwk):
    """Non-trivial splits of an already-computed tree of blobs (a plain tree).

    Used to consume the output of PhyloNetworks' `treeofblobs` (the trusted
    reference).  Since a tree of blobs is a tree, every edge is a bridge, so
    this is `tree_of_blobs_splits` applied to that tree.
    """
    return tree_of_blobs_splits(tob_nwk)


def split_blob_adjacency(nwk):
    """Map each non-trivial network split -> whether its cut edge touches a blob.

    An edge is inside a blob iff it is *not* a bridge, so the "blob nodes" are
    exactly the endpoints of non-bridge edges.  A split's inducing bridge is
    "blob-adjacent" when one of its endpoints is a blob node, i.e. the clade
    hangs directly off a reticulation cycle -- the structurally hard case.

    Returns (all_leaves, {split_key: blob_adjacent_bool}).
    """
    edges, leaf_names = parse_extended_newick(nwk)
    num_nodes = 1 + max((max(u, v) for u, v in edges), default=-1)
    all_leaves = frozenset(leaf_names.values())
    ref = "OUT" if "OUT" in all_leaves else min(all_leaves)

    bridges = find_bridges(num_nodes, edges)
    blob_nodes = set()
    for idx, (u, v) in enumerate(edges):
        if idx not in bridges:
            blob_nodes.add(u)
            blob_nodes.add(v)

    info = {}
    for eidx in bridges:
        u, v = edges[eidx]
        side = _reachable_leaves(edges, eidx, v, leaf_names)
        other = all_leaves - side
        if len(side) < 2 or len(other) < 2:
            continue
        key = frozenset(other if ref in side else side)
        info[key] = (u in blob_nodes) or (v in blob_nodes)
    return all_leaves, info


def leaf_blob_adjacency(nwk):
    """Return {taxon: bool} where True means the leaf attaches directly to a
    blob node (its pendant edge lands on a reticulation cycle)."""
    edges, leaf_names = parse_extended_newick(nwk)
    num_nodes = 1 + max((max(u, v) for u, v in edges), default=-1)
    bridges = find_bridges(num_nodes, edges)
    blob_nodes = set()
    for idx, (u, v) in enumerate(edges):
        if idx not in bridges:
            blob_nodes.add(u); blob_nodes.add(v)
    # neighbour of each leaf
    nbr = {}
    for u, v in edges:
        if u in leaf_names:
            nbr[u] = v
        if v in leaf_names:
            nbr[v] = u
    return {leaf_names[nid]: (nbr.get(nid) in blob_nodes) for nid in leaf_names}


if __name__ == "__main__":
    import os
    base = os.path.join(os.path.dirname(__file__), "..", "..", "data", "camus-dataset", "n15")
    for d in sorted(os.listdir(base))[:5]:
        net = os.path.join(base, d, "true_net.nwk")
        if not os.path.exists(net):
            continue
        leaves, splits = tree_of_blobs_splits(open(net).read())
        print(f"n15/{d}: {len(leaves)} leaves, {len(splits)} non-trivial ToB splits")
        for sp in sorted(splits, key=lambda s: (len(s), sorted(s))):
            print("   ", "{" + ",".join(sorted(sp, key=lambda x: (len(x), x))) + "}")
