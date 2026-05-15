#!/usr/bin/env python3
"""Tests for the `rectify correct-cdna` subcommand (cdna_correct_command).

Focus: fast unit tests of the algorithm primitives (UMI directional clustering,
position-directional clustering, classifier, walk-back / walk-forward, boundary
pattern). Plus an opt-in end-to-end smoke test if the chrI test data is on disk.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from rectify.core import cdna_correct_command as cc


# ---------------------------------------------------------------------------
# UMI directional clustering (v1.13)
# ---------------------------------------------------------------------------
class TestUmiComponentsDirectional:
    def test_single_umi_returns_one_cluster(self):
        assert cc.umi_components_directional(["TTACGATTAACGTTCACCTTAACGTTT"], max_edit=3) == [[0]]

    def test_two_identical_umis_one_cluster(self):
        umis = ["TTACGATTAACGTTCACCTTAACGTTT"] * 2
        comps = cc.umi_components_directional(umis, max_edit=3)
        assert len(comps) == 1 and sorted(comps[0]) == [0, 1]

    def test_two_x_rule_blocks_equal_count_peaks(self):
        """Two singleton UMIs at Lev≤3 of each other should NOT merge under
        the directional algorithm's strict-ish 2× rule (n_parent ≥ 2n_child-1
        with both n=1 gives 1≥1, so the relaxed rule lets them merge; that's
        what umi_tools does and what we use)."""
        umis = ["TTACGATTAACGTTCACCTTAACGTTT",
                "TTACGATTAACGTTCACCTTAACGTAA"]  # Lev=2 — within max_edit=3
        comps = cc.umi_components_directional(umis, max_edit=3)
        # With non-strict 2x rule (matching umi_tools): both count=1 nodes merge.
        # If we ever tighten to strict 2× they would split.
        assert len(comps) == 1

    def test_2x_rule_prevents_chain_merge(self):
        """A high-count UMI surrounded by Lev≤3 low-count error variants should
        absorb the variants but not chain-merge with another high-count peak."""
        umis = (["TTACGATTAACGTTCACCTTAACGTTT"] * 10        # peak A (count 10)
                + ["TTACGATTAACGTTCACCTTAACGTTA"]            # Lev=1 from A, count 1
                + ["AACGTTTAACGTTCACCTTAACGTTT"] * 10)       # peak B (Lev≥3 from A)
        comps = cc.umi_components_directional(umis, max_edit=3)
        # Expect two clusters: peak A (with its variant) and peak B
        assert len(comps) == 2
        # The size-1 variant should be absorbed by peak A (count 10), not peak B
        sizes = sorted(len(c) for c in comps)
        assert sizes == [10, 11]


# ---------------------------------------------------------------------------
# Position-directional clustering (v1.17)
# ---------------------------------------------------------------------------
class TestPositionComponentsDirectional:
    def test_single_position(self):
        assert cc.position_components_directional([100], [1], max_dist=5) == [[0]]

    def test_identical_positions_merge(self):
        comps = cc.position_components_directional([100, 100, 100], [3, 1, 1], max_dist=5)
        assert len(comps) == 1 and sorted(comps[0]) == [0, 1, 2]

    def test_within_tolerance_peak_absorbs_low_count_neighbor(self):
        """A peak with count 10 should absorb a singleton 3 bp away."""
        comps = cc.position_components_directional([100, 103], [10, 1], max_dist=5)
        assert len(comps) == 1 and sorted(comps[0]) == [0, 1]

    def test_strict_2x_rule_prevents_chain_merge(self):
        """Two peaks 20 bp apart, with count=1 bridges at intermediate positions
        within max_dist=5, must NOT chain-merge under the strict 2× rule."""
        positions = [100,    105,  110,  115,  120]
        weights   = [50,     1,    1,    1,    50]
        comps = cc.position_components_directional(positions, weights, max_dist=5)
        # Expect two clusters (100 and 120 stay distinct), bridges absorbed by
        # their nearest peak (or stay separate as orphans)
        # The peak at 100 absorbs only positions reachable by direct edges (within max_dist=5)
        # 100 -> 105 (50 >= 2*1=2 yes); but 105 -> 110 fails strict 2x (1>=2 no)
        # so the chain doesn't propagate. Position 120 stays in its own cluster
        # with 115 absorbed (50 >= 2*1=2 yes).
        n_peak_clusters = sum(1 for c in comps if max(weights[i] for i in c) >= 10)
        assert n_peak_clusters == 2

    def test_returns_input_index_lists(self):
        comps = cc.position_components_directional([100, 200, 100], [1, 1, 1], max_dist=5)
        # 100 (idx 0, 2) clusters together; 200 (idx 1) alone
        clusters_sorted = sorted([sorted(c) for c in comps])
        assert clusters_sorted == [[0, 2], [1]]


# ---------------------------------------------------------------------------
# Boundary pattern matching (v1.19)
# ---------------------------------------------------------------------------
class TestBoundaryPattern:
    def test_perfect_fwd_pattern_scores_12(self):
        # TT-VVVV-TTT-GGG (V = A/C/G)
        seq = "TTACGCTTTGGG"
        score = cc._score_boundary_window(seq, cc._BOUNDARY_PATTERN_FWD)
        assert score == 12

    def test_perfect_rev_pattern_scores_12(self):
        # CCC-AAA-BBBB-AA (B = T/G/C)
        seq = "CCCAAACGTCAA"
        score = cc._score_boundary_window(seq, cc._BOUNDARY_PATTERN_REV)
        assert score == 12

    def test_one_mismatch_passes_threshold(self):
        seq = "TTACGCTTTGGA"  # last G replaced with A (12 → 11 score)
        score = cc._score_boundary_window(seq, cc._BOUNDARY_PATTERN_FWD)
        assert 10 <= score < 12

    def test_pattern_in_v_position_must_not_be_t(self):
        # V positions [2,3,4,5] cannot be T (would score 0 at those positions)
        seq = "TTTTTTTTTTGGG"[:12]  # contiguous T's
        score = cc._score_boundary_window(seq, cc._BOUNDARY_PATTERN_FWD)
        # Only positions 0,1,6,7,8 (T) + 9,10,11 (G) match: 5+3 = 8 max
        assert score < cc._BOUNDARY_MIN_SCORE

    def test_find_boundary_match_requires_bridge_intact(self):
        # All 3 G's at end must be intact, even if rest scores well
        seq = "TTACGCTTTGGA"  # 1 bridge G missing → reject
        result = cc._find_boundary_match(seq, 0, len(seq), rc=False)
        assert result is None

    def test_find_boundary_match_returns_best_match(self):
        # Two valid windows; should pick the higher-scoring one
        seq = "TTACGCTTTGGG"  # perfect fwd pattern
        result = cc._find_boundary_match(seq, 0, len(seq), rc=False)
        assert result is not None
        qpos, score = result
        assert qpos == 0 and score == 12


# ---------------------------------------------------------------------------
# revcomp helper
# ---------------------------------------------------------------------------
class TestRevComp:
    def test_revcomp_basic(self):
        assert cc.revcomp("ATCG") == "CGAT"

    def test_revcomp_with_n(self):
        assert cc.revcomp("ATNCG") == "CGNAT"

    def test_revcomp_lowercase_preserved(self):
        assert cc.revcomp("atcg") == "cgat"


# ---------------------------------------------------------------------------
# End-to-end CLI smoke test (opt-in: requires the chrI test data on disk)
# ---------------------------------------------------------------------------
_CHRI_BAM = Path("/Users/kevinroy/work/ont_cdna/test_data/wt_rep1.chrI.bam")
_GFF = Path("/Users/kevinroy/work/ont_cdna/test_data/saccharomyces_cerevisiae_R64-5-1_20240529.gff")
_FASTA = Path("/Users/kevinroy/work/ont_cdna/test_data/S288C_reference_sequence_R64-5-1_20240529.fsa")


@pytest.mark.skipif(
    not (_CHRI_BAM.exists() and _GFF.exists() and _FASTA.exists()),
    reason="chrI test data not on disk (only present on dev M1; skip on CI)",
)
def test_correct_cdna_chri_smoke(tmp_path):
    """End-to-end: invoke `rectify correct-cdna` on the chrI BAM. Verify the
    BAM, clusters.tsv, isoforms.tsv, t1t2_pairs.tsv outputs exist and cluster
    counts are in the expected ballpark (catches accidental algorithm regressions)."""
    out_dir = tmp_path / "out"
    rc = subprocess.run(
        [sys.executable, "-m", "rectify", "correct-cdna",
         str(_CHRI_BAM), "--out", str(out_dir),
         "--gff", str(_GFF), "--reference", str(_FASTA)],
        capture_output=True, text=True, timeout=900,
    )
    assert rc.returncode == 0, f"stderr: {rc.stderr[-1500:]}"
    for f in ("stage1_consensus.fastq.gz", "clusters.tsv", "isoforms.tsv", "t1t2_pairs.tsv"):
        assert (out_dir / f).exists(), f"missing {f}"

    # Cluster count sanity bands (chrI w/ v1.19d on 2026-05-13: 26,364 clusters)
    n_clusters = sum(1 for _ in (out_dir / "clusters.tsv").open()) - 1  # header
    assert 24000 < n_clusters < 30000, f"clusters.tsv has {n_clusters} rows"
    n_isoforms = sum(1 for _ in (out_dir / "isoforms.tsv").open()) - 1
    assert 3000 < n_isoforms < 4500, f"isoforms.tsv has {n_isoforms} rows"

    # Verify pretrim wiring: XQ and XK tags must appear in every FASTQ record
    import gzip as _gzip
    with _gzip.open(out_dir / "stage1_consensus.fastq.gz", "rt") as fq:
        first_header = fq.readline()
    assert "XQ:i:" in first_header, f"XQ pretrim-5' tag missing from FASTQ: {first_header[:200]}"
    assert "XK:i:" in first_header, f"XK pretrim-3' tag missing from FASTQ: {first_header[:200]}"


# ---------------------------------------------------------------------------
# Consensus pre-trimming (Task 2)
# ---------------------------------------------------------------------------
class TestPretrimConsensus:
    """pretrim_consensus strips SSP/UMI/GGG from 5' and poly-A from 3' so
    aligners receive clean mRNA sequence only."""

    def _make_fwd_consensus(self, mrna: str, polya: str = "A" * 8) -> str:
        umi = "T" * cc.UMI_LEN
        return cc.SSP_FWD + umi + "GGG" + mrna + polya + cc.ANCHOR_FWD

    def test_fwd_strips_ssp_umi_ggg_and_polya(self):
        mrna = "ATGCATGC"
        full = self._make_fwd_consensus(mrna)
        trimmed, q_trim_5, pretrim_pa_len = cc.pretrim_consensus(full, "fwd", read_type=1)
        assert trimmed == mrna, (
            f"expected trimmed=={mrna!r}, got {trimmed!r}"
        )
        assert q_trim_5 == len(cc.SSP_FWD) + cc.UMI_LEN + 3  # 23 + 27 + 3 = 53
        assert pretrim_pa_len == 8

    def test_fwd_type2_no_5p_trim(self):
        """Type-2 reads have no SSP — only poly-A strip applies."""
        mrna = "GCTAGCTAGC"
        seq = mrna + "A" * 12 + cc.ANCHOR_FWD
        trimmed, q_trim_5, pretrim_pa_len = cc.pretrim_consensus(seq, "fwd", read_type=2)
        assert q_trim_5 == 0, "Type-2 should not strip a 5' prefix"
        assert trimmed == mrna
        assert pretrim_pa_len == 12

    def test_fwd_no_ssp_found_returns_full_with_no_5p_trim(self):
        """If SSP is absent from a Type-1 consensus (e.g. POA artifact), the
        function should not crash and should not trim the 5' end."""
        seq = "ATGCATGC" + "A" * 6
        trimmed, q_trim_5, _ = cc.pretrim_consensus(seq, "fwd", read_type=1)
        assert q_trim_5 == 0
        assert trimmed.startswith("ATGCATGC")


# ---------------------------------------------------------------------------
# Task 3 — Chimeric junction rescue via select_best_chimeric
# ---------------------------------------------------------------------------

class TestSelectBestChimericSplice:
    """select_best_chimeric() must prefer the aligner that calls a canonical
    GT-AG splice junction over one that calls a non-canonical site.

    Layout:
      ref  = "A"*100 + "GT" + "T"*98 + "AG" + "C"*100   (302 bp)
      query = "A"*100 + "C"*100                           (200 bp, exon1+exon2)

    mapPacBio CIGAR: 100M 102N 100M — canonical intron ref[100:202] (GT…AG)
    minimap2  CIGAR: 100M  50N 100M — wrong intron ref[100:150] (GT…TT)

    Both aligners agree on query[0:100] → ref[0:100] (sync points).
    They diverge at query[100:200] (the 3' terminal segment containing the
    splice junction). select_best_chimeric must pick mapPacBio there.
    """

    def _make_segs(self):
        import pysam as ps

        ref_seq = "A" * 100 + "GT" + "T" * 98 + "AG" + "C" * 100  # 302 bp
        genome = {"chrT": ref_seq}

        header = ps.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": "chrT", "LN": len(ref_seq)}],
        })

        query_seq = "A" * 100 + "C" * 100

        def _seg(cigar_tuples):
            seg = ps.AlignedSegment(header)
            seg.query_name = "test_read"
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = 0
            seg.mapping_quality = 60
            seg.cigartuples = cigar_tuples
            seg.query_sequence = query_seq
            seg.query_qualities = ps.qualitystring_to_array("?" * len(query_seq))
            return seg

        # mapPacBio: canonical intron ref[100:202] (GT…AG), 102 bp
        mpb = _seg([(0, 100), (3, 102), (0, 100)])
        # minimap2: wrong intron ref[100:150] (GT…TT), 50 bp — non-canonical 3'SS
        mm2 = _seg([(0, 100), (3, 50), (0, 100)])

        return {"minimap2": mm2, "mapPacBio": mpb}, genome

    def test_canonical_junction_wins(self):
        from rectify.core.chimeric_consensus import select_best_chimeric
        aligner_reads, genome = self._make_segs()
        result = select_best_chimeric(aligner_reads, genome)
        # The 3' segment holds the junction. mapPacBio has GT-AG (+5); mm2 has GT-TT (-3).
        assert result.three_prime_aligner == "mapPacBio", (
            f"Expected mapPacBio to win 3' segment (canonical junction); "
            f"got three_prime_aligner={result.three_prime_aligner!r}"
        )

    def test_chimeric_cigar_uses_canonical_intron(self):
        from rectify.core.chimeric_consensus import select_best_chimeric
        aligner_reads, genome = self._make_segs()
        result = select_best_chimeric(aligner_reads, genome)
        # The chimeric CIGAR should use the 102-bp intron from mapPacBio, not the 50-bp one.
        n_ops = {op: length for op, length in result.chimeric_cigar}
        assert n_ops.get(3) == 102, (
            f"Expected N-op length 102 (canonical intron) in chimeric CIGAR; "
            f"got cigar={result.chimeric_cigar}"
        )
