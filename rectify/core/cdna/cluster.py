"""Stage-1 cluster_reads: anchor-bucket reads, then UMI/position-cluster within
each bucket. Also exposes pick_representative for singleton / consensus-fallback.
"""
from __future__ import annotations

from collections import defaultdict
from typing import List, Tuple

import pysam

from .read_info import ReadInfo
from .umi import (
    position_components,
    umi_components,
    umi_components_directional,
)


def cluster_reads(reads: List[ReadInfo], anchor_window: int, max_edit: int,
                  per_cluster_cap: int,
                  clustering_method: str = "directional"
                  ) -> Tuple[List[List[ReadInfo]], dict]:
    """Stage-1 clustering: anchor-bucket reads, then UMI-cluster within each bucket.

    clustering_method:
      - "directional" (v1.13 default): umi_tools-style; breaks polyA-hotspot
        chain-merge. Per-bucket cap is advisory only (large buckets logged but
        still processed).
      - "components" (legacy): connected-components; needs `per_cluster_cap` to
        drop oversized buckets to avoid catastrophic chain-merge.
    """
    # Bucket key includes read_type (v1.15) — Type-1 and Type-2 reads never
    # mix within a bucket, so they get different clustering algorithms.
    anchor_buckets: dict = defaultdict(list)
    for r in reads:
        bucket = (r.chrom, r.anchor // anchor_window, r.orient, r.read_type)
        anchor_buckets[bucket].append(r)

    stats = dict(
        anchor_buckets=len(anchor_buckets),
        buckets_dropped_polyA_pileup=0,
        reads_in_dropped_buckets=0,
        biggest_bucket_size=0,
        molecule_clusters=0,
        molecule_clusters_size_1=0,
        molecule_clusters_size_2=0,
        molecule_clusters_size_ge_3=0,
        type1_clusters=0,
        type2_clusters=0,
        type1_reads=sum(1 for r in reads if r.read_type == 1),
        type2_reads=sum(1 for r in reads if r.read_type == 2),
    )

    umi_cluster_fn = umi_components_directional if clustering_method == "directional" else umi_components

    out_clusters: List[List[ReadInfo]] = []
    for bucket_key, bucket_reads in anchor_buckets.items():
        read_type = bucket_key[3]
        bsize = len(bucket_reads)
        if bsize > stats["biggest_bucket_size"]:
            stats["biggest_bucket_size"] = bsize
        # Legacy components path still respects the cap to avoid chain-merge (Type-1 only).
        if read_type == 1 and clustering_method == "components" and bsize > per_cluster_cap:
            stats["buckets_dropped_polyA_pileup"] += 1
            stats["reads_in_dropped_buckets"] += bsize
            continue
        if read_type == 1:
            umis = [r.umi for r in bucket_reads]
            comps = umi_cluster_fn(umis, max_edit)
        else:
            comps = position_components(bucket_reads)
        for comp in comps:
            cluster = [bucket_reads[i] for i in comp]
            out_clusters.append(cluster)
            sz = len(cluster)
            if sz == 1: stats["molecule_clusters_size_1"] += 1
            elif sz == 2: stats["molecule_clusters_size_2"] += 1
            else: stats["molecule_clusters_size_ge_3"] += 1
            if read_type == 1: stats["type1_clusters"] += 1
            else: stats["type2_clusters"] += 1
    stats["molecule_clusters"] = len(out_clusters)
    return out_clusters, stats


def pick_representative(reads_in_cluster: List[pysam.AlignedSegment]) -> pysam.AlignedSegment:
    """Choose the cluster representative (umi_tools-style longest-and-best).

    Used for singletons (size 1) and as fallback when pileup consensus fails.
    """
    def key(r: pysam.AlignedSegment) -> tuple:
        ref_len = (r.reference_end or 0) - r.reference_start
        return (-ref_len, -r.mapping_quality, r.query_name)
    return min(reads_in_cluster, key=key)
