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
    next_node_id = 0

    def new_node():
        nonlocal next_node_id
        nid = next_node_id
        next_node_id += 1
        nodes[nid] = {"name": None, "hybrid": None, "children": []}
        return nid

    position = 0
    n = len(s)

    def parse_subtree():
        nonlocal position
        if position >= n:
            raise NetworkParseError("unexpected end of extended Newick string")

        children = []
        if s[position] == "(":
            position += 1
            while True:
                children.append(parse_subtree())
                if position >= n:
                    raise NetworkParseError("unexpected end while reading child list")
                delimiter = s[position]
                if delimiter == ",":
                    position += 1
                    continue
                if delimiter == ")":
                    position += 1
                    break
                context = s[position:position + 20]
                raise NetworkParseError(
                    f"expected ',' or ')' at position {position}: ...{context}"
                )

        # Read the label token: everything up to ',', ')', ':', or end.
        label_start = position
        while position < n and s[position] not in ",():":
            position += 1
        label = s[label_start:position]

        # Skip branch-length / support / gamma annotations ":len:sup:gamma".
        if position < n and s[position] == ":":
            while position < n and s[position] not in ",()":
                position += 1

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
    if position != n:
        raise NetworkParseError(
            f"trailing characters after root at position {position}: {s[position:]}"
        )

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
    adjacency = defaultdict(list)           # node -> list of (neighbor, edge_index)
    for edge_index, (u, v) in enumerate(edges):
        adjacency[u].append((v, edge_index))
        adjacency[v].append((u, edge_index))

    discovery_time = [-1] * num_nodes
    low = [0] * num_nodes
    bridges = set()
    next_discovery_time = 0

    sys.setrecursionlimit(10000)

    def dfs(node, incoming_edge):
        nonlocal next_discovery_time
        discovery_time[node] = low[node] = next_discovery_time
        next_discovery_time += 1
        for neighbor, edge_index in adjacency[node]:
            if edge_index == incoming_edge:
                continue
            if discovery_time[neighbor] == -1:
                dfs(neighbor, edge_index)
                low[node] = min(low[node], low[neighbor])
                if low[neighbor] > discovery_time[node]:
                    bridges.add(edge_index)
            else:
                low[node] = min(low[node], discovery_time[neighbor])

    for start in range(num_nodes):
        if discovery_time[start] == -1:
            dfs(start, -1)
    return bridges


def _reachable_leaves(edges, exclude_edge_idx, source, leaf_ids):
    """Leaves reachable from `source` in the graph with one edge removed."""
    adjacency = defaultdict(list)
    for edge_index, (u, v) in enumerate(edges):
        if edge_index == exclude_edge_idx:
            continue
        adjacency[u].append(v)
        adjacency[v].append(u)
    seen = {source}
    stack = [source]
    while stack:
        x = stack.pop()
        for neighbor in adjacency[x]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return {leaf_ids[node_id] for node_id in seen if node_id in leaf_ids}


def _reference_taxon(all_leaves):
    """Use OUT when present; otherwise choose a deterministic reference taxon."""
    return "OUT" if "OUT" in all_leaves else min(all_leaves)


def _canonical_split(side, all_leaves, reference_taxon):
    """Represent a split by the side that excludes the reference taxon."""
    other = all_leaves - side
    return frozenset(other if reference_taxon in side else side)


def _blob_nodes(edges, bridges):
    """Nodes incident to a non-bridge edge are precisely the nodes in blobs."""
    nodes = set()
    for edge_index, (u, v) in enumerate(edges):
        if edge_index not in bridges:
            nodes.update((u, v))
    return nodes


def tree_of_blobs_splits(nwk):
    """Compute the network's leaf set and its non-trivial tree-of-blobs splits.

    Returns (all_leaves, splits) where
      all_leaves : frozenset of every taxon name
      splits     : set of frozensets.  Each frozenset is one SIDE of a
                   non-trivial split, canonicalised as the side that does NOT
                   contain the outgroup 'OUT' (or, if 'OUT' is absent, the
                   lexicographically first taxon). Both sides have >= 2 taxa.
    """
    edges, leaf_names = parse_extended_newick(nwk)
    num_nodes = 1 + max((max(u, v) for u, v in edges), default=-1)
    all_leaves = frozenset(leaf_names.values())
    reference_taxon = _reference_taxon(all_leaves)

    bridges = find_bridges(num_nodes, edges)
    splits = set()
    for edge_index in bridges:
        _, v = edges[edge_index]
        side = _reachable_leaves(edges, edge_index, v, leaf_names)
        other = all_leaves - side
        if len(side) < 2 or len(other) < 2:
            continue                          # trivial split
        key = _canonical_split(side, all_leaves, reference_taxon)
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
    reference_taxon = _reference_taxon(all_leaves)

    bridges = find_bridges(num_nodes, edges)
    blob_nodes = _blob_nodes(edges, bridges)

    info = {}
    for edge_index in bridges:
        u, v = edges[edge_index]
        side = _reachable_leaves(edges, edge_index, v, leaf_names)
        other = all_leaves - side
        if len(side) < 2 or len(other) < 2:
            continue
        key = _canonical_split(side, all_leaves, reference_taxon)
        info[key] = (u in blob_nodes) or (v in blob_nodes)
    return all_leaves, info


def leaf_blob_adjacency(nwk):
    """Return whether each leaf attaches directly to a blob node.

    Equivalently, the leaf's pendant edge lands on a reticulation cycle.
    """
    edges, leaf_names = parse_extended_newick(nwk)
    num_nodes = 1 + max((max(u, v) for u, v in edges), default=-1)
    bridges = find_bridges(num_nodes, edges)
    blob_nodes = _blob_nodes(edges, bridges)

    # Record the graph neighbor attached to each leaf.
    leaf_neighbors = {}
    for u, v in edges:
        if u in leaf_names:
            leaf_neighbors[u] = v
        if v in leaf_names:
            leaf_neighbors[v] = u
    return {
        leaf_names[node_id]: leaf_neighbors.get(node_id) in blob_nodes
        for node_id in leaf_names
    }


if __name__ == "__main__":
    import os
    base = os.path.join(os.path.dirname(__file__), "..", "..", "data", "camus-dataset", "n15")
    for dataset in sorted(os.listdir(base))[:5]:
        net = os.path.join(base, dataset, "true_net.nwk")
        if not os.path.exists(net):
            continue
        with open(net) as network_file:
            leaves, splits = tree_of_blobs_splits(network_file.read())
        print(f"n15/{dataset}: {len(leaves)} leaves, {len(splits)} non-trivial ToB splits")
        for split in sorted(splits, key=lambda side: (len(side), sorted(side))):
            ordered_taxa = sorted(split, key=lambda taxon: (len(taxon), taxon))
            print("   ", "{" + ",".join(ordered_taxa) + "}")
