"""Type-2 (SSP-less, no-UMI) reads must NOT be deduplicated.

Regression guard for a real defect. Stage 1 used to route Type-2 buckets through
``position_components``, collapsing exact ``(aln_start, aln_end)`` matches and calling them
PCR duplicates. Type-2 reads carry **no UMI**, so there is no evidence by which to call any
two of them the same molecule — the collapse measured **positional concentration, not
amplification**, which is the exact substitution the Chanfreau DRS/cDNA invariants forbid.

Measured on the 827 splicing_aa cohort (18 libraries): the old path removed **51.5 % of all
Type-2 reads** (13,292,754 -> 6,450,950), and the rate scaled with DEPTH — 4-6 % on the ~50 k
libraries vs 44-57 % on the multi-million-read ones — while true UMI-measured PCR duplication
on the SAME libraries was 24-41 %. Duplication does not scale with depth like that; positional
crowding does. The excess was distinct molecules merged away, understating Type-2 abundance
roughly two-fold in a depth-dependent manner.
"""
from __future__ import annotations

from dataclasses import dataclass

from rectify.core.cdna.cluster import cluster_reads


@dataclass
class _R:
    chrom: str
    anchor: int
    orient: str
    read_type: int
    umi: str
    aln_start: int
    aln_end: int
    read_id: str = "r"
    xf_tier: int = 0
    tail_len: int = 0
    is_reverse: bool = False
    read_subtype: str = "umi_not_captured"


def _t2(n, start=100, end=500):
    """n Type-2 reads sharing IDENTICAL coordinates — the worst case for the old path."""
    return [_R("I", start, "fwd", 2, "", start, end, read_id=f"t2_{i}") for i in range(n)]


def test_type2_reads_are_never_collapsed_by_default():
    reads = _t2(10)
    clusters, stats = cluster_reads(reads, anchor_window=25, max_edit=3, per_cluster_cap=200)
    assert len(clusters) == 10, "co-located Type-2 reads must stay 10 separate observations"
    assert all(len(c) == 1 for c in clusters)
    assert stats["type2_clusters"] == 10
    assert stats["type2_reads"] == 10
    assert stats.get("type2_reads_collapsed", 0) == 0


def test_legacy_position_collapse_is_opt_in_and_still_works():
    """The old behaviour stays reachable ONLY to reproduce pre-fix outputs."""
    reads = _t2(10)
    clusters, stats = cluster_reads(reads, anchor_window=25, max_edit=3, per_cluster_cap=200,
                                    type2_collapse="position")
    assert len(clusters) == 1, "legacy path collapses identical coordinates into one"
    assert stats["type2_reads_collapsed"] == 9


def test_type1_umi_dedup_is_unaffected_by_the_type2_change():
    """The fix must not touch Type-1: UMI-anchored dedup there is correct and stays."""
    u = "TTACGTTTCCGTTTAAGTTTCCCTTT"
    reads = [_R("I", 100, "fwd", 1, u, 100, 500, read_id=f"t1_{i}",
                read_subtype="umi_captured_fwd") for i in range(4)]
    clusters, stats = cluster_reads(reads, anchor_window=25, max_edit=3, per_cluster_cap=200)
    assert len(clusters) == 1, "identical UMIs are real PCR siblings and SHOULD merge"
    assert stats["type1_clusters"] == 1
    assert stats["type1_reads"] == 4


def test_mixed_bucket_keeps_type1_dedup_and_type2_observations():
    u = "TTACGTTTCCGTTTAAGTTTCCCTTT"
    reads = ([_R("I", 100, "fwd", 1, u, 100, 500, read_id=f"t1_{i}",
                 read_subtype="umi_captured_fwd") for i in range(3)]
             + _t2(5))
    clusters, stats = cluster_reads(reads, anchor_window=25, max_edit=3, per_cluster_cap=200)
    # 1 Type-1 molecule (3 siblings merged) + 5 uncollapsed Type-2 observations
    assert stats["type1_clusters"] == 1
    assert stats["type2_clusters"] == 5
    assert len(clusters) == 6
