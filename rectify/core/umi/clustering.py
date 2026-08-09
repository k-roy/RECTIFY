"""Protocol-agnostic UMI / position clustering primitives.

Moved here from ``rectify/core/cdna/umi.py`` (v1.13+) because these algorithms are
not cDNA-specific: they operate on plain ``List[str]`` + Levenshtein and know
nothing about ONT, pysam, or polyA. ``core/cdna/umi.py`` re-exports them, so the
cdna pipeline's public API is unchanged.

The ReadInfo-coupled helpers (``position_components``, ``umi_canon_placeholder``)
deliberately stay in ``core/cdna/umi.py`` — they depend on the cdna ReadInfo
dataclass and have no meaning for short reads.

Algorithms:
  * connected-components (legacy, pre-v1.13) — chains independent molecules at
    hotspots; retained for ``--umi-clustering=components`` comparison.
  * umi_tools-style directional clustering (v1.13+) — breaks transitivity via
    the 2x rule (Smith, Heger & Sudbery 2017, Genome Res 27:491-499).
  * directional clustering on integer positions, for isoform grouping.
"""
from __future__ import annotations

import logging
import os
from collections import defaultdict
from typing import Dict, List, Optional

from rapidfuzz.distance import Levenshtein

logger = logging.getLogger(__name__)

# 651: a pip sdist fallback silently installs rapidfuzz as pure Python (Root-Is-Purelib), which
# runs this module's distance computations interpreted (~50x slower; cost 200 CPU-h on H2 before
# it was caught — planning/619). Warn loudly at import so it can never be silent again.
_HAS_CPP_BACKEND = "cpp" in getattr(Levenshtein.distance, "__module__", "")
if not _HAS_CPP_BACKEND:
    logger.warning(
        "rapidfuzz C++ backend unavailable (Levenshtein.distance resolves to %s) — UMI "
        "clustering will run ~50x slower. Fix: pip install --force-reinstall "
        "--only-binary=:all: rapidfuzz",
        getattr(Levenshtein.distance, "__module__", "?"))


def _levenshtein_neighbour_sets(unique: List[str], max_edit: int
                                ) -> Optional[Dict[int, set]]:
    """Exact undirected neighbour sets {i -> {j}} for Lev(u_i, u_j) <= max_edit.

    651: computes the pair set ONCE per bucket with rapidfuzz ``process.cdist`` — the same
    predicate and cutoff as the per-pair loop, but batched on the C++ side (row-chunked so the
    dense block stays ~<=32 MB; entries above the cutoff come back as cutoff+1 and are masked
    out). The result feeds BOTH the directional adjacency build and ``_split_into_stars``, which
    previously paid a second quadratic distance pass on percolated components.

    Returns None when the batched path is unavailable (pure-Python rapidfuzz backend, or numpy /
    process missing) — callers fall back to the original per-pair loop.

    ``RECTIFY_UMI_CDIST_WORKERS`` (default 1) sets cdist's thread count. It stays 1 by default
    because the production caller (`correct-cdna --workers N`) already runs one process per
    region; raise it only for serial single-region runs.
    """
    if not _HAS_CPP_BACKEND:
        return None
    try:
        import numpy as np
        from rapidfuzz import process as _rf_process
    except ImportError:
        return None
    neighbours: Dict[int, set] = {}
    n = len(unique)
    if n < 2:
        return neighbours
    workers = 1
    env = os.environ.get("RECTIFY_UMI_CDIST_WORKERS")
    if env:
        try:
            workers = max(1, int(env))
        except ValueError:
            pass
    block_rows = max(1, min(n, (32 << 20) // n))
    for r0 in range(0, n, block_rows):
        rows = unique[r0:r0 + block_rows]
        cols = unique[r0:]  # upper triangle only (plus a negligible in-block lower wedge)
        block = _rf_process.cdist(rows, cols, scorer=Levenshtein.distance,
                                  score_cutoff=max_edit, dtype=np.uint8,
                                  workers=workers)
        ii, jj = np.nonzero(block <= max_edit)
        for il, jl in zip(ii.tolist(), jj.tolist()):
            i = r0 + il
            j = r0 + jl
            if j > i:
                neighbours.setdefault(i, set()).add(j)
                neighbours.setdefault(j, set()).add(i)
    return neighbours



def umi_components(umis: List[str], max_edit: int) -> List[List[int]]:
    """Connected-components clustering of UMIs by Levenshtein <= max_edit.

    LEGACY (pre-v1.13). Chains independent molecules together when many UMIs
    in a bucket are within Lev<=max_edit (the "polyA-hotspot chain-merge"
    pathology). Retained for comparison / debugging via --umi-clustering=components.

    NOTE (short-read callers): this is single-linkage. It is the documented
    transitivity trap (fgbio calls the equivalent ``--strategy=edit``): a~b and
    b~c groups a with c even when dist(a,c) > max_edit. Prefer
    :func:`umi_components_directional` unless you are deliberately comparing.
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
      node a -> b iff Lev(a, b) <= max_edit AND count(a) >= 2*count(b) - 1.

    Then walk from the highest-count unvisited node via outgoing edges (BFS);
    everything reached forms one molecule cluster. Repeat with next highest-count
    unvisited node.

    The 2x rule constrains transitive merging: a small-count UMI (likely an error
    variant of a deeper UMI nearby) gets absorbed by its parent, but cannot itself
    absorb other UMIs of comparable depth. This breaks the chain-merge that
    connected-components suffers at polyA hotspots where thousands of independent
    molecules' UMIs land within Lev <= 3 of each other.

    WHY THE COUNT GRADIENT IS THE REAL DISCRIMINATOR (matters when choosing
    max_edit): a UMI produced by a *sequencing error* off a true UMI is necessarily
    RARER than its parent, because it arises from a subset of that parent's reads.
    Two genuinely distinct molecules have NO expected count relationship. So the
    2x rule is a likelihood-ratio-like test separating "error off a parent" from
    "two real molecules that happen to be neighbours" -- it, not the distance
    threshold, is what makes max_edit=1 safe. It does NOT scale to protect
    max_edit=2 equally (the Hamming ball of a 12-mer grows 37 -> 631).

    For singleton-dominated buckets (every UMI has count 1, the 2x rule reduces
    to "1 >= 1", allowing edges between any two Lev-<=-max_edit singletons), the
    behavior collapses back to connected-components -- that's fine; directional
    helps where it's needed (multi-read buckets with mixed coverage levels).
    """
    n = len(umis)
    if n == 0: return []
    if n == 1: return [[0]]

    # Group read indices by exact UMI -> unique UMI list + counts.
    umi_to_read_indices: Dict[str, List[int]] = defaultdict(list)
    for i, u in enumerate(umis):
        umi_to_read_indices[u].append(i)
    unique = list(umi_to_read_indices.keys())
    counts = [len(umi_to_read_indices[u]) for u in unique]
    n_uniq = len(unique)

    # NOTE (planning/633, /634): an earlier guard here refused merging when a BUCKET was
    # almost all singletons ("no amplification evidence"). It was REMOVED because it is
    # refuted: healthy buckets reach 0.9997 distinct/read — HIGHER than broken ones — so no
    # cut on that ratio separates them, and the per-CLUSTER version fails for the same reason
    # (median distinct/reads ~1.0 at EVERY cluster size, because at ONT error rates on a 27-nt
    # UMI two genuine duplicates usually do NOT share an exact UMI). Counting distinct UMIs
    # cannot distinguish the two ways a cluster acquires them. The GEOMETRY can — see
    # `_split_into_stars` below, which is the replacement.

    # Order unique UMIs by count desc (ties broken by UMI lex order for determinism).
    order_desc = sorted(range(n_uniq), key=lambda i: (-counts[i], unique[i]))
    pos_in_order = [0] * n_uniq
    for rank, idx in enumerate(order_desc):
        pos_in_order[idx] = rank

    # --- Candidate-neighbour finding ---------------------------------------
    # The directional edge predicate (below) is:  a->b iff dist(a,b) <= max_edit
    # AND count(a) >= 2*count(b) - 1. The naive way to find the dist-<=-max_edit
    # candidates is all-pairs Levenshtein -- O(n_uniq^2), which HANGS on ultra-deep
    # single-position buckets (e.g. PGK1's ~21k-read CPA pileup, ~4.6e8 pairs).
    #
    # FAST PATH (max_edit == 1 AND all UNIQUE umis the same length): for equal-length
    # strings Levenshtein-1 <=> Hamming-1 (a single insertion/deletion changes the
    # length, so the only length-preserving edit is a substitution). Hamming-<=1
    # pairs are found in O(n_uniq * L) via masked-key hashing: two umis that agree
    # everywhere except (at most) one position share a "position-p-masked" key.
    # This is PROVABLY IDENTICAL to the all-pairs Levenshtein(score_cutoff=1) result,
    # just without the quadratic blow-up.
    #
    # Fixed-length random N-mer UMIs (Illumina: Lexogen CORALL N12, and most
    # standard chemistries) satisfy the length precondition by construction, so
    # short-read callers using the recommended max_edit=1 always take this path.
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

    # 651: for the ed>1 path, precompute the exact undirected pair set once (batched C++
    # cdist). Reused twice: adjacency filter below, and _split_into_stars membership (which
    # previously re-computed distances hub-by-hub — a second quadratic pass on percolated
    # components). None -> per-pair fallback loop. On the fast_hamming path the masked-key
    # candidate sets ARE the exact Lev<=1 neighbour sets (equal length: Lev 1 <=> Hamming 1).
    neighbour_sets: Optional[Dict[int, set]] = None
    if not fast_hamming:
        neighbour_sets = _levenshtein_neighbour_sets(unique, max_edit)
    star_neighbours: Optional[Dict[int, set]] = (
        candidate_neighbours if fast_hamming else neighbour_sets)

    # Build adjacency: directed edge a->b for b later in order_desc (count(b) <=
    # count(a)) whose count satisfies the 2x rule (count(b) <= (count(a)+1)//2) and
    # whose distance to a is <= max_edit.
    adj: Dict[int, List[int]] = defaultdict(list)
    for ai_pos, a in enumerate(order_desc):
        ca = counts[a]
        max_cb = (ca + 1) // 2  # need cb <= this for 2x rule
        ua = unique[a]
        if fast_hamming:
            # Only Hamming-1 candidates can be children; filter to later-in-order +
            # 2x rule. Sort by order rank for deterministic adjacency.
            cand = [b for b in candidate_neighbours[a]
                    if pos_in_order[b] > ai_pos and counts[b] <= max_cb]
            for b in sorted(cand, key=lambda b: pos_in_order[b]):
                adj[a].append(b)
        elif neighbour_sets is not None:
            # 651: exact pair set precomputed by batched cdist — filter to later-in-order +
            # 2x rule; no per-pair Python-loop distance calls. Same edge set as the loop below.
            cand = [b for b in neighbour_sets.get(a, ())
                    if pos_in_order[b] > ai_pos and counts[b] <= max_cb]
            for b in sorted(cand, key=lambda b: pos_in_order[b]):
                adj[a].append(b)
        else:
            # Fallback (batched path unavailable): all-pairs Levenshtein.
            # counts in order_desc are non-increasing; skip until first b with cb <= max_cb.
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
        # ---- 635: ENFORCE THE STAR INVARIANT (radius, not chain depth) ----------------
        # A cluster asserts that all members are error copies of ONE true UMI, so every
        # member must lie within max_edit of an OBSERVED centre. BFS does not enforce that:
        # it admits transitive chains a~b~c whose ends are >= 2*max_edit apart.
        for star in _split_into_stars(component_uniqs, unique, counts, max_edit,
                                      neighbours=star_neighbours):
            read_indices: List[int] = []
            for uid in star:
                read_indices.extend(umi_to_read_indices[unique[uid]])
            clusters.append(read_indices)
    return clusters


def _split_into_stars(component: List[int], unique: List[str], counts: List[int],
                      max_edit: int,
                      neighbours: Optional[Dict[int, set]] = None) -> List[List[int]]:
    """Partition a connected component into STARS, each centred on an OBSERVED member.

    WHY. Connected-components/BFS answers "is there a path?", but a UMI cluster claims
    something stronger: that every member is an error copy of ONE true sequence. Under that
    claim every member must be within `max_edit` of a single centre — and that centre must be
    something we actually observed, not an inferred consensus. BFS instead admits chains
    a~b~c where a and c are 2*max_edit apart or more, which merges distinct molecules.

    MEASURED (planning/635b — matched control, same BAM, same reads, ed the ONLY variable;
    radius = edit distance from each member to the cluster consensus):

        ed=1 : median radius 0, max 1, 100.0% of members within ed<=3   -> STARS
        ed=3 : median radius 6, max 16,  22.3% within ed<=3             -> CHAINS
               and it degrades with size: 36% at 10-30 reads -> 0.0% above 1,000 reads

    So the chaining is caused by the edit radius allowing transitive closure, and this
    restores the invariant directly rather than trying to predict which buckets will break
    (three count-based predictors were tried and refuted — planning/626, /633, /634).

    GREEDY, and a NO-OP on components that are already valid stars: take the highest-count
    remaining member as the hub, absorb everything within `max_edit` of it, repeat. If the
    component was a genuine star the first pass absorbs all of it and output is unchanged —
    which is why this is safe to apply unconditionally, including at max_edit=1.

    Deterministic: hubs are chosen by (-count, umi) so ties break on UMI lex order.

    651: when the caller already holds the exact neighbour sets (batched cdist, or the
    fast_hamming masked-key candidates at ed=1), membership is a set lookup and NO distances are
    re-computed here. With ``neighbours=None`` the original per-pair verification runs unchanged.
    """
    if len(component) <= 1:
        return [component]
    remaining = sorted(component, key=lambda i: (-counts[i], unique[i]))
    stars: List[List[int]] = []
    while remaining:
        hub = remaining[0]
        hub_umi = unique[hub]
        hub_neigh = neighbours.get(hub, ()) if neighbours is not None else None
        star = [hub]
        rest: List[int] = []
        for i in remaining[1:]:
            if (i in hub_neigh) if hub_neigh is not None else (
                    Levenshtein.distance(hub_umi, unique[i],
                                         score_cutoff=max_edit) <= max_edit):
                star.append(i)
            else:
                rest.append(i)
        stars.append(star)
        remaining = rest
    return stars


def position_components_directional(positions: List[int], weights: List[int],
                                     max_dist: int) -> List[List[int]]:
    """v1.17: directional clustering on integer positions for isoform grouping.

    Parallel to umi_components_directional but with two differences:
      1. Distance = abs(pos[a] - pos[b]) (integer), not Levenshtein.
      2. STRICT 2x rule: edge a->b iff w(a) >= 2*w(b) AND |pos[a]-pos[b]| <= max_dist.
         (Not the umi_tools "n >= 2m - 1" -- strict prevents bridge-chains where
         count=1 valley positions link two real peaks through a chain of edges.)

    `positions` and `weights` are parallel arrays -- one entry per input cluster.
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

    # Build adjacency. STRICT 2x rule + max_dist.
    adj: Dict[int, List[int]] = defaultdict(list)
    for ai_pos, a in enumerate(order_desc):
        wa = node_weight[a]
        pa = unique[a]
        # Iterate later positions in count-desc order. Need w(b) <= wa/2 (strict)
        # AND |pa - pb| <= max_dist.
        for bi_pos in range(ai_pos + 1, n_nodes):
            b = order_desc[bi_pos]
            if node_weight[b] * 2 > wa:
                continue  # strict 2x rule fails
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


def canonical_umi(cluster_umis: List[str]) -> str:
    """Pick the canonical UMI for a cluster (most frequent; ties -> first lexicographically)."""
    counts: defaultdict[str, int] = defaultdict(int)
    for u in cluster_umis:
        counts[u] += 1
    return min(counts, key=lambda u: (-counts[u], u))
