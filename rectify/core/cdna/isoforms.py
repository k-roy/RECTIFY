"""v1.17 / v2: Stage-2 isoform grouping + Type-1↔Type-2 same-molecule reconciliation.

  * :func:`assign_isoforms` — within each (gene, orient, read_type) group, run
    directional clustering on 5' and 3' positions (Type-1) or 3' only (Type-2).
    Emits per-isoform records with modal termini, tail length, and the cluster
    membership list.
  * :func:`reconcile_t1_t2_pairs` — link Type-1 ↔ Type-2 clusters of the same
    molecule by termini-distance one-to-one greedy match (replaces the legacy
    cross-orient pair check, which had no biological substrate).
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .read_info import ReadInfo
from .umi import position_components_directional, umi_canon_placeholder


def assign_isoforms(clusters: List[List[ReadInfo]],
                     cluster_xg: Dict[int, Optional[str]],
                     tol5: int, tol3: int
                     ) -> Tuple[Dict[int, str], List[dict]]:
    """v1.17: group Stage-1 clusters into isoform clusters.

    Within each (gene, orient, read_type) group:
      Type-1 (SSP+UMI captured, 5' = real mRNA 5' end):
        - directional clustering on 5' positions (tol5, strict 2× rule)
        - directional clustering on 3' positions (tol3, strict 2× rule)
        - isoform = (gene, orient, type=1, 5'-group, 3'-group)
      Type-2 (SSP-less, 5' is random truncation point):
        - skip 5' clustering (positions are noisy)
        - directional clustering on 3' positions (tol3, strict 2× rule)
        - isoform = (gene, orient, type=2, "*", 3'-group)

    Weights = cluster.n_reads (PCR-collapsed support).

    Returns:
      cluster_to_isoform_id: dict cluster_id → isoform_id string
      isoform_records: list of dicts (one per isoform) with aggregated stats
    """
    cluster_to_iso: Dict[int, str] = {}
    iso_records: List[dict] = []

    # Group cluster indices by (gene, orient, read_type)
    groups: Dict[Tuple[str, str, int], List[int]] = defaultdict(list)
    for cid, c in enumerate(clusters):
        gene = cluster_xg.get(cid)
        if gene is None: continue  # skip unannotated
        groups[(gene, c[0].orient, c[0].read_type)].append(cid)

    for (gene, orient, rt), cids in groups.items():
        # 5'/3' position per cluster (mRNA-relative coords).
        # orient=fwd: 5' = aln_start, 3' = aln_end-1.
        # orient=rev: 5' = aln_end-1, 3' = aln_start.
        pos5: List[int] = []
        pos3: List[int] = []
        weights: List[int] = []
        for cid in cids:
            r0 = clusters[cid][0]
            # v1.19: use pos5_corrected (bridge walk-forward applied for Type-1;
            # raw aln-side position for Type-2).
            pos5.append(r0.pos5_corrected)
            if orient == "fwd":
                pos3.append(r0.aln_end - 1)
            else:
                pos3.append(r0.aln_start)
            weights.append(len(clusters[cid]))

        # 3' clustering (used by both Type 1 and Type 2)
        comps3 = position_components_directional(pos3, weights, tol3)
        cid_to_3g = {}
        for g3_idx, comp in enumerate(comps3):
            for local_idx in comp:
                cid_to_3g[cids[local_idx]] = g3_idx

        if rt == 1:
            comps5 = position_components_directional(pos5, weights, tol5)
            cid_to_5g = {}
            for g5_idx, comp in enumerate(comps5):
                for local_idx in comp:
                    cid_to_5g[cids[local_idx]] = g5_idx
            # Build isoform IDs and group clusters by (5'-group, 3'-group)
            iso_buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
            for cid in cids:
                iso_buckets[(cid_to_5g[cid], cid_to_3g[cid])].append(cid)
            for (g5, g3), iso_cids in iso_buckets.items():
                iso_id = f"{gene}_t1_5g{g5}_3g{g3}"
                _emit_isoform(iso_id, gene, orient, 1, iso_cids, clusters,
                              cluster_to_iso, iso_records, has_5g=True)
        else:
            # Type-2: only 3' grouping
            iso_buckets3: Dict[int, List[int]] = defaultdict(list)
            for cid in cids:
                iso_buckets3[cid_to_3g[cid]].append(cid)
            for g3, iso_cids in iso_buckets3.items():
                iso_id = f"{gene}_t2_3g{g3}"
                _emit_isoform(iso_id, gene, orient, 2, iso_cids, clusters,
                              cluster_to_iso, iso_records, has_5g=False)

    return cluster_to_iso, iso_records


def reconcile_t1_t2_pairs(clusters: List[List[ReadInfo]],
                            cluster_xg: Dict[int, Optional[str]],
                            tol5: int, tol3: int
                            ) -> Tuple[Dict[int, int], List[dict]]:
    """v2: link Type-1 ↔ Type-2 clusters of the SAME molecule (same gene + same orient
    + both 5' and 3' termini within tol bp).

    Biology: a single dsDNA molecule produces some reads from Strand A (SSP captured
    at reliable basecalled-3' → Type-1, UMI extractable) and some from Strand B
    where SSP at basecalled-5' gets truncated → Type-2 (no UMI). The two read
    populations end up in different Stage-1 clusters (Type-1 in UMI cluster,
    Type-2 in position cluster). Reconciliation rescues the link.

    Cross-orient pairs (the v1.13 Stage-2 idea) don't exist in PCB114 chemistry —
    both physical PCR strands of the same molecule give the same orient label
    after BAM normalization. The Type-1↔Type-2 same-orient pair IS the real
    biological substrate.

    Matching is one-to-one greedy: each Type-2 cluster pairs with its nearest
    Type-1 candidate (closest sum |Δ5'|+|Δ3'|). Each Type-1 can be claimed by
    at most one Type-2 (first-match-by-distance wins).

    Returns:
      cluster_xl: dict cluster_id → partner_cluster_id (bidirectional links)
      pair_records: list of dicts (one per pair) for the TSV manifest
    """
    cluster_xl: Dict[int, int] = {}
    pair_records: List[dict] = []

    # Group cluster indices by (gene, orient)
    groups: Dict[Tuple[str, str], Dict[int, List[int]]] = defaultdict(
        lambda: {1: [], 2: []})
    for cid, c in enumerate(clusters):
        gene = cluster_xg.get(cid)
        if gene is None: continue
        groups[(gene, c[0].orient)][c[0].read_type].append(cid)

    for (gene, orient), by_type in groups.items():
        t1_cids = by_type[1]
        t2_cids = by_type[2]
        if not t1_cids or not t2_cids:
            continue
        # Build (5', 3') for each
        def termini(cid: int) -> Tuple[int, int]:
            r0 = clusters[cid][0]
            # v1.19 pos5_corrected (bridge walk-forward); 3' is raw polyA-side
            if orient == "fwd":
                return (r0.pos5_corrected, r0.aln_end - 1)
            return (r0.pos5_corrected, r0.aln_start)

        t1_termini = {cid: termini(cid) for cid in t1_cids}
        t2_termini = {cid: termini(cid) for cid in t2_cids}

        # Greedy match: for each Type-2, find closest available Type-1 within tol
        claimed_t1: set = set()
        # Sort Type-2 by n_reads desc (claim higher-coverage Type-2 first)
        t2_sorted = sorted(t2_cids, key=lambda c: -len(clusters[c]))
        for t2 in t2_sorted:
            p5_2, p3_2 = t2_termini[t2]
            best_t1 = None
            best_dist = None
            for t1 in t1_cids:
                if t1 in claimed_t1: continue
                p5_1, p3_1 = t1_termini[t1]
                d5 = abs(p5_1 - p5_2)
                d3 = abs(p3_1 - p3_2)
                if d5 > tol5 or d3 > tol3: continue
                dist = d5 + d3
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    best_t1 = t1
            if best_t1 is None: continue
            claimed_t1.add(best_t1)
            cluster_xl[best_t1] = t2   # Type-1 → Type-2 partner cid
            cluster_xl[t2] = best_t1   # Type-2 → Type-1 partner cid
            p5_1, p3_1 = t1_termini[best_t1]
            pair_records.append({
                "t1_cid": best_t1, "t2_cid": t2,
                "gene": gene, "orient": orient,
                "t1_pos5": p5_1, "t1_pos3": p3_1,
                "t2_pos5": p5_2, "t2_pos3": p3_2,
                "d5": abs(p5_1 - p5_2), "d3": abs(p3_1 - p3_2),
                "t1_n_reads": len(clusters[best_t1]),
                "t2_n_reads": len(clusters[t2]),
                "t1_umi": umi_canon_placeholder(clusters[best_t1]),
            })
    return cluster_xl, pair_records


def _emit_isoform(iso_id: str, gene: str, orient: str, read_type: int,
                   iso_cids: List[int], clusters: List[List[ReadInfo]],
                   cluster_to_iso: Dict[int, str], iso_records: List[dict],
                   has_5g: bool) -> None:
    pos5_vals = []; pos3_vals = []; n_reads_total = 0
    tail_lens: List[int] = []
    chrom = clusters[iso_cids[0]][0].chrom
    read_subtype = clusters[iso_cids[0]][0].read_subtype
    for cid in iso_cids:
        cluster_to_iso[cid] = iso_id
        r0 = clusters[cid][0]
        pos5_vals.append(r0.pos5_corrected)
        pos3_vals.append(r0.aln_end - 1 if orient == "fwd" else r0.aln_start)
        n_reads_total += len(clusters[cid])
        # Each read's tail_len weighted by appearing in the cluster
        for r in clusters[cid]:
            tail_lens.append(r.tail_len)
    def _modal(xs):
        if not xs: return -1
        c: defaultdict[int, int] = defaultdict(int)
        for x in xs: c[x] += 1
        return max(c, key=lambda k: (c[k], -k))
    def _median(xs):
        if not xs: return 0
        s = sorted(xs); return s[len(s) // 2]
    iso_records.append({
        "isoform_id": iso_id,
        "gene": gene,
        "chrom": chrom,
        "orient": orient,
        "read_type": read_type,
        "read_subtype": read_subtype,
        "n_clusters": len(iso_cids),
        "n_reads_total": n_reads_total,
        "pos5_modal": _modal(pos5_vals) if has_5g else -1,
        "pos3_modal": _modal(pos3_vals),
        "tail_len_median": _median(tail_lens),
        "cluster_ids": ",".join(str(c) for c in iso_cids),
    })
