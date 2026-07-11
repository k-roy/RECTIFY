"""UMI / position clustering primitives for the cdna pipeline.

Two graph-clustering algorithms used at Stage 1:
  * connected-components (legacy, pre-v1.13) — chains independent molecules at
    polyA hotspots; retained for ``--umi-clustering=components`` comparison.
  * umi_tools-style directional clustering (v1.13+) — breaks transitivity via
    the 2× rule, resolving polyA-hotspot mega-buckets into individual molecule
    clusters.

Parallel implementations on integer positions are used at Stage 2 to group
clusters into isoforms (5'/3' position grouping with a stricter 2× rule).

Also provides :func:`canonical_umi` / :func:`umi_canon_placeholder` helpers
that pick the canonical UMI of a cluster for tagging.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Tuple

from rapidfuzz.distance import Levenshtein

from .read_info import ReadInfo


def umi_components(umis: List[str], max_edit: int) -> List[List[int]]:
    """Connected-components clustering of UMIs by Levenshtein ≤ max_edit.

    LEGACY (pre-v1.13). Chains independent molecules together when many UMIs
    in a bucket are within Lev≤max_edit (the "polyA-hotspot chain-merge"
    pathology). Retained for comparison / debugging via --umi-clustering=components.
    """
    n = len(umis)
    if n == 0:
        return []
    parent = list(range(n))
    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    for i in range(n):
        for j in range(i + 1, n):
            if Levenshtein.distance(umis[i], umis[j], score_cutoff=max_edit) <= max_edit:
                union(i, j)
    comps: defaultdict[int, list] = defaultdict(list)
    for i in range(n):
        comps[find(i)].append(i)
    return list(comps.values())


def umi_components_directional(umis: List[str], max_edit: int) -> List[List[int]]:
    """umi_tools-style directional clustering (Smith, Heger & Sudbery 2017, Genome Res).

    Build a directed graph over UNIQUE UMIs in the bucket:
      node a → b iff Lev(a, b) ≤ max_edit AND count(a) ≥ 2*count(b) − 1.

    Then walk from the highest-count unvisited node via outgoing edges (BFS);
    everything reached forms one molecule cluster. Repeat with next highest-count
    unvisited node.

    The 2× rule constrains transitive merging: a small-count UMI (likely an error
    variant of a deeper UMI nearby) gets absorbed by its parent, but cannot itself
    absorb other UMIs of comparable depth. This breaks the chain-merge that
    connected-components suffers at polyA hotspots where thousands of independent
    molecules' UMIs land within Lev ≤ 3 of each other.

    For singleton-dominated buckets (every UMI has count 1, the 2× rule reduces
    to "1 ≥ 1", allowing edges between any two Lev-≤-max_edit singletons), the
    behavior collapses back to connected-components — that's fine; directional
    helps where it's needed (multi-read buckets with mixed coverage levels).
    """
    n = len(umis)
    if n == 0: return []
    if n == 1: return [[0]]

    # Group read indices by exact UMI → unique UMI list + counts.
    umi_to_read_indices: Dict[str, List[int]] = defaultdict(list)
    for i, u in enumerate(umis):
        umi_to_read_indices[u].append(i)
    unique = list(umi_to_read_indices.keys())
    counts = [len(umi_to_read_indices[u]) for u in unique]
    n_uniq = len(unique)

    # Order unique UMIs by count desc (ties broken by UMI lex order for determinism).
    order_desc = sorted(range(n_uniq), key=lambda i: (-counts[i], unique[i]))
    pos_in_order = [0] * n_uniq
    for rank, idx in enumerate(order_desc):
        pos_in_order[idx] = rank

    # --- Candidate-neighbour finding ---------------------------------------
    # The directional edge predicate (below) is:  a→b iff dist(a,b) ≤ max_edit
    # AND count(a) ≥ 2·count(b) − 1. The naive way to find the dist-≤-max_edit
    # candidates is all-pairs Levenshtein — O(n_uniq²), which HANGS on ultra-deep
    # single-position buckets (e.g. PGK1's ~21k-read CPA pileup, ~4.6e8 pairs).
    #
    # FAST PATH (max_edit == 1 AND all UNIQUE umis the same length): for equal-length
    # strings Levenshtein-1 ⟺ Hamming-1 (a single insertion/deletion changes the
    # length, so the only length-preserving edit is a substitution). Hamming-≤1
    # pairs are found in O(n_uniq · L) via masked-key hashing: two umis that agree
    # everywhere except (at most) one position share a "position-p-masked" key.
    # This is PROVABLY IDENTICAL to the all-pairs Levenshtein(score_cutoff=1) result,
    # just without the quadratic blow-up. cdna Type-1 umis are always exactly
    # UMI_LEN chars (see extract_read_info), so this path is the norm.
    uniq_lens = {len(u) for u in unique}
    fast_hamming = (max_edit == 1 and len(uniq_lens) == 1 and next(iter(uniq_lens)) >= 1)

    candidate_neighbours: Dict[int, set] = defaultdict(set)
    if fast_hamming:
        L = next(iter(uniq_lens))
        for p in range(L):
            key_buckets: Dict[str, List[int]] = defaultdict(list)
            for idx, u in enumerate(unique):
                key_buckets[u[:p] + u[p + 1:]].append(idx)
            for group in key_buckets.values():
                if len(group) > 1:
                    for a in group:
                        for b in group:
                            if a != b:
                                candidate_neighbours[a].add(b)

    # Build adjacency: directed edge a→b for b later in order_desc (count(b) ≤
    # count(a)) whose count satisfies the 2× rule (count(b) ≤ (count(a)+1)//2) and
    # whose distance to a is ≤ max_edit.
    adj: Dict[int, List[int]] = defaultdict(list)
    for ai_pos, a in enumerate(order_desc):
        ca = counts[a]
        max_cb = (ca + 1) // 2  # need cb ≤ this for 2× rule
        ua = unique[a]
        if fast_hamming:
            # Only Hamming-1 candidates can be children; filter to later-in-order +
            # 2× rule. Sort by order rank for deterministic adjacency.
            cand = [b for b in candidate_neighbours[a]
                    if pos_in_order[b] > ai_pos and counts[b] <= max_cb]
            for b in sorted(cand, key=lambda b: pos_in_order[b]):
                adj[a].append(b)
        else:
            # Fallback (max_edit > 1 or mixed-length umis): all-pairs Levenshtein.
            # counts in order_desc are non-increasing; skip until first b with cb ≤ max_cb.
            bi_start = ai_pos + 1
            while bi_start < n_uniq and counts[order_desc[bi_start]] > max_cb:
                bi_start += 1
            for bi_pos in range(bi_start, n_uniq):
                b = order_desc[bi_pos]
                d = Levenshtein.distance(ua, unique[b], score_cutoff=max_edit)
                if d <= max_edit:
                    adj[a].append(b)

    # BFS from each unvisited highest-count root.
    visited: set = set()
    clusters: List[List[int]] = []
    for root in order_desc:
        if root in visited: continue
        component_uniqs = []
        queue = [root]
        while queue:
            u = queue.pop()
            if u in visited: continue
            visited.add(u)
            component_uniqs.append(u)
            for v in adj[u]:
                if v not in visited:
                    queue.append(v)
        # Expand back to read indices.
        read_indices: List[int] = []
        for uid in component_uniqs:
            read_indices.extend(umi_to_read_indices[unique[uid]])
        clusters.append(read_indices)
    return clusters


def position_components_directional(positions: List[int], weights: List[int],
                                     max_dist: int) -> List[List[int]]:
    """v1.17: directional clustering on integer positions for isoform grouping.

    Parallel to umi_components_directional but with two differences:
      1. Distance = abs(pos[a] - pos[b]) (integer), not Levenshtein.
      2. STRICT 2× rule: edge a→b iff w(a) ≥ 2·w(b) AND |pos[a]-pos[b]| ≤ max_dist.
         (Not the umi_tools "n ≥ 2m − 1" — strict prevents bridge-chains where
         count=1 valley positions link two real peaks through a chain of edges.)

    `positions` and `weights` are parallel arrays — one entry per input cluster.
    Multiple input clusters at the same position get their weights summed and
    are treated as one "node" for the directional graph; this is exactly the
    behavior we want (a peak's strength = total read support, not count of
    sub-clusters).

    Returns list of components, each as a list of indices into the input arrays.
    """
    n = len(positions)
    if n == 0: return []
    if n == 1: return [[0]]

    # Group input indices by exact position; weight per node = SUM of input weights.
    pos_to_indices: Dict[int, List[int]] = defaultdict(list)
    for i, p in enumerate(positions):
        pos_to_indices[p].append(i)
    unique = sorted(pos_to_indices.keys())  # sorted by position (helps locality)
    node_weight = [sum(weights[i] for i in pos_to_indices[p]) for p in unique]
    n_nodes = len(unique)

    # Order nodes by weight desc; ties broken by position for determinism.
    order_desc = sorted(range(n_nodes), key=lambda i: (-node_weight[i], unique[i]))

    # Build adjacency. STRICT 2× rule + max_dist.
    adj: Dict[int, List[int]] = defaultdict(list)
    for ai_pos, a in enumerate(order_desc):
        wa = node_weight[a]
        pa = unique[a]
        # Iterate later positions in count-desc order. Need w(b) ≤ wa/2 (strict)
        # AND |pa - pb| ≤ max_dist.
        for bi_pos in range(ai_pos + 1, n_nodes):
            b = order_desc[bi_pos]
            if node_weight[b] * 2 > wa:
                continue  # strict 2× rule fails
            if abs(unique[b] - pa) > max_dist:
                continue
            adj[a].append(b)

    # BFS from each unvisited highest-weight root.
    visited: set = set()
    clusters: List[List[int]] = []
    for root in order_desc:
        if root in visited: continue
        component_nodes = []
        queue = [root]
        while queue:
            u = queue.pop()
            if u in visited: continue
            visited.add(u)
            component_nodes.append(u)
            for v in adj[u]:
                if v not in visited:
                    queue.append(v)
        # Expand back to input-cluster indices.
        cluster_indices: List[int] = []
        for nid in component_nodes:
            cluster_indices.extend(pos_to_indices[unique[nid]])
        clusters.append(cluster_indices)
    return clusters


def position_components(reads: List[ReadInfo]) -> List[List[int]]:
    """v1.15: exact (aln_start, aln_end) grouping for Type-2 (SSP-less) reads.

    No UMI is available for Type-2 reads, so we collapse exact-position duplicates
    only. This is conservative — molecules sharing the SAME aln_start AND aln_end
    are treated as PCR/sequencing duplicates of one decay/truncation event.
    Wobble of even ±1 bp produces a new cluster; that's fine because real distinct
    decay-intermediate species typically have well-defined termini.
    """
    pos_to_indices: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i, r in enumerate(reads):
        pos_to_indices[(r.aln_start, r.aln_end)].append(i)
    return list(pos_to_indices.values())


def canonical_umi(cluster_umis: List[str]) -> str:
    """Pick the canonical UMI for a cluster (most frequent; ties → first lexicographically)."""
    counts: defaultdict[str, int] = defaultdict(int)
    for u in cluster_umis:
        counts[u] += 1
    return min(counts, key=lambda u: (-counts[u], u))


def umi_canon_placeholder(cluster: List[ReadInfo]) -> str:
    """Return the canonical UMI of a Type-1 cluster (most-frequent UMI string).
    For Type-2 clusters this returns empty string."""
    counts: defaultdict[str, int] = defaultdict(int)
    for r in cluster:
        counts[r.umi] += 1
    if not counts: return ""
    return min(counts, key=lambda u: (-counts[u], u))
