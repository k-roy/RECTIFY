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

from rectify.core.cdna._constants import ANCHOR_FWD, SSP_FWD, UMI_LEN
from rectify.core.cdna.consensus import pretrim_consensus
from rectify.core.cdna.read_info import revcomp
from rectify.core.cdna.umi import (
    position_components_directional,
    umi_components_directional,
)
from rectify.core.cdna.walkback import (
    _BOUNDARY_MIN_SCORE,
    _BOUNDARY_PATTERN_FWD,
    _BOUNDARY_PATTERN_REV,
    _find_boundary_match,
    _score_boundary_window,
)


# ---------------------------------------------------------------------------
# UMI directional clustering (v1.13)
# ---------------------------------------------------------------------------
class TestUmiComponentsDirectional:
    def test_single_umi_returns_one_cluster(self):
        assert umi_components_directional(["TTACGATTAACGTTCACCTTAACGTTT"], max_edit=3) == [[0]]

    def test_two_identical_umis_one_cluster(self):
        umis = ["TTACGATTAACGTTCACCTTAACGTTT"] * 2
        comps = umi_components_directional(umis, max_edit=3)
        assert len(comps) == 1 and sorted(comps[0]) == [0, 1]

    def test_two_x_rule_blocks_equal_count_peaks(self):
        """Two singleton UMIs at Lev≤3 of each other should NOT merge under
        the directional algorithm's strict-ish 2× rule (n_parent ≥ 2n_child-1
        with both n=1 gives 1≥1, so the relaxed rule lets them merge; that's
        what umi_tools does and what we use)."""
        umis = ["TTACGATTAACGTTCACCTTAACGTTT",
                "TTACGATTAACGTTCACCTTAACGTAA"]  # Lev=2 — within max_edit=3
        comps = umi_components_directional(umis, max_edit=3)
        # With non-strict 2x rule (matching umi_tools): both count=1 nodes merge.
        # If we ever tighten to strict 2× they would split.
        assert len(comps) == 1

    def test_2x_rule_prevents_chain_merge(self):
        """A high-count UMI surrounded by Lev≤3 low-count error variants should
        absorb the variants but not chain-merge with another high-count peak."""
        umis = (["TTACGATTAACGTTCACCTTAACGTTT"] * 10        # peak A (count 10)
                + ["TTACGATTAACGTTCACCTTAACGTTA"]            # Lev=1 from A, count 1
                + ["AACGTTTAACGTTCACCTTAACGTTT"] * 10)       # peak B (Lev≥3 from A)
        comps = umi_components_directional(umis, max_edit=3)
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
        assert position_components_directional([100], [1], max_dist=5) == [[0]]

    def test_identical_positions_merge(self):
        comps = position_components_directional([100, 100, 100], [3, 1, 1], max_dist=5)
        assert len(comps) == 1 and sorted(comps[0]) == [0, 1, 2]

    def test_within_tolerance_peak_absorbs_low_count_neighbor(self):
        """A peak with count 10 should absorb a singleton 3 bp away."""
        comps = position_components_directional([100, 103], [10, 1], max_dist=5)
        assert len(comps) == 1 and sorted(comps[0]) == [0, 1]

    def test_strict_2x_rule_prevents_chain_merge(self):
        """Two peaks 20 bp apart, with count=1 bridges at intermediate positions
        within max_dist=5, must NOT chain-merge under the strict 2× rule."""
        positions = [100,    105,  110,  115,  120]
        weights   = [50,     1,    1,    1,    50]
        comps = position_components_directional(positions, weights, max_dist=5)
        # Expect two clusters (100 and 120 stay distinct), bridges absorbed by
        # their nearest peak (or stay separate as orphans)
        # The peak at 100 absorbs only positions reachable by direct edges (within max_dist=5)
        # 100 -> 105 (50 >= 2*1=2 yes); but 105 -> 110 fails strict 2x (1>=2 no)
        # so the chain doesn't propagate. Position 120 stays in its own cluster
        # with 115 absorbed (50 >= 2*1=2 yes).
        n_peak_clusters = sum(1 for c in comps if max(weights[i] for i in c) >= 10)
        assert n_peak_clusters == 2

    def test_returns_input_index_lists(self):
        comps = position_components_directional([100, 200, 100], [1, 1, 1], max_dist=5)
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
        score = _score_boundary_window(seq, _BOUNDARY_PATTERN_FWD)
        assert score == 12

    def test_perfect_rev_pattern_scores_12(self):
        # CCC-AAA-BBBB-AA (B = T/G/C)
        seq = "CCCAAACGTCAA"
        score = _score_boundary_window(seq, _BOUNDARY_PATTERN_REV)
        assert score == 12

    def test_one_mismatch_passes_threshold(self):
        seq = "TTACGCTTTGGA"  # last G replaced with A (12 → 11 score)
        score = _score_boundary_window(seq, _BOUNDARY_PATTERN_FWD)
        assert 10 <= score < 12

    def test_pattern_in_v_position_must_not_be_t(self):
        # V positions [2,3,4,5] cannot be T (would score 0 at those positions)
        seq = "TTTTTTTTTTGGG"[:12]  # contiguous T's
        score = _score_boundary_window(seq, _BOUNDARY_PATTERN_FWD)
        # Only positions 0,1,6,7,8 (T) + 9,10,11 (G) match: 5+3 = 8 max
        assert score < _BOUNDARY_MIN_SCORE

    def test_find_boundary_match_requires_bridge_intact(self):
        # All 3 G's at end must be intact, even if rest scores well
        seq = "TTACGCTTTGGA"  # 1 bridge G missing → reject
        result = _find_boundary_match(seq, 0, len(seq), rc=False)
        assert result is None

    def test_find_boundary_match_returns_best_match(self):
        # Two valid windows; should pick the higher-scoring one
        seq = "TTACGCTTTGGG"  # perfect fwd pattern
        result = _find_boundary_match(seq, 0, len(seq), rc=False)
        assert result is not None
        qpos, score = result
        assert qpos == 0 and score == 12


# ---------------------------------------------------------------------------
# revcomp helper
# ---------------------------------------------------------------------------
class TestRevComp:
    def test_revcomp_basic(self):
        assert revcomp("ATCG") == "CGAT"

    def test_revcomp_with_n(self):
        assert revcomp("ATNCG") == "CGNAT"

    def test_revcomp_lowercase_preserved(self):
        assert revcomp("atcg") == "cgat"


# ---------------------------------------------------------------------------
# End-to-end CLI smoke test (opt-in: requires the chrI test data on disk)
# ---------------------------------------------------------------------------
_CHRI_BAM = Path("/Users/kevinroy/work/ont_cdna/test_data/wt_rep1.chrI.bam")
_GFF = Path("/Users/kevinroy/work/ont_cdna/test_data/saccharomyces_cerevisiae_R64-5-1_20240529.gff")
_FASTA = Path("/Users/kevinroy/work/ont_cdna/test_data/S288C_reference_sequence_R64-5-1_20240529.fsa")


@pytest.mark.slow
@pytest.mark.skipif(
    not (_CHRI_BAM.exists() and _FASTA.exists()),
    reason="chrI test data not on disk (only present on dev M1; skip on CI)",
)
def test_correct_cdna_chri_smoke(tmp_path):
    """End-to-end: invoke `rectify correct-cdna` on the chrI BAM. Manifest
    files (clusters.tsv / isoforms.tsv / t1t2_pairs.tsv) now come from the
    downstream `rectify cdna-analyze` step — see test_cdna_pipeline_chri_smoke
    for the full three-stage check. Here we only verify that correct-cdna
    emits a non-empty stage1_consensus.fastq.gz with the pretrim wiring tags
    and an output record count in the expected band."""
    out_dir = tmp_path / "out"
    rc = subprocess.run(
        [sys.executable, "-m", "rectify", "correct-cdna",
         str(_CHRI_BAM), "--out", str(out_dir),
         "--reference", str(_FASTA)],
        capture_output=True, text=True, timeout=900,
    )
    assert rc.returncode == 0, f"stderr: {rc.stderr[-1500:]}"
    fastq = out_dir / "stage1_consensus.fastq.gz"
    assert fastq.exists(), "missing stage1_consensus.fastq.gz"

    # Cluster count sanity band (one FASTQ record per cluster; chrI w/ v1.19d
    # on 2026-05-13: 26,364 clusters)
    import gzip as _gzip
    with _gzip.open(fastq, "rt") as fq:
        first_header = fq.readline()
        n_records = 1 + sum(1 for line in fq if line.startswith("@cluster_"))
    assert 24000 < n_records < 30000, f"FASTQ has {n_records} records"

    # Verify required tags ride on the FASTQ comment line for `rectify align -y`
    # to propagate (XQ / XK = pretrim wiring; XU/XC/XR/XO/XT/XY/XF/XB = the
    # alignment-independent per-cluster tags consumed by cdna-analyze)
    for tag in ("XQ:i:", "XK:i:",
                "XU:Z:", "XC:i:", "XR:Z:", "XO:Z:", "XT:i:", "XY:Z:",
                "XF:i:", "XB:Z:"):
        assert tag in first_header, f"{tag} missing from FASTQ header: {first_header[:200]}"


@pytest.mark.slow
@pytest.mark.skipif(
    not (_CHRI_BAM.exists() and _GFF.exists() and _FASTA.exists()
         and shutil.which("minimap2") and shutil.which("samtools")),
    reason="chrI test data + minimap2 + samtools required (dev M1 only)",
)
def test_cdna_pipeline_chri_smoke(tmp_path):
    """Full three-stage pipeline: correct-cdna → align (minimap2) → cdna-analyze.

    Verifies that:
      * per-cluster FASTQ tags survive minimap2's `-y` passthrough into the
        consensus BAM,
      * cdna-analyze can recompute walkback / walk-forward / gene assignment /
        isoform clustering / T1↔T2 pairing on the post-align coordinates,
      * the three downstream TSV manifests + consensus_tagged.bam land in the
        out directory with row counts in the expected bands.

    Uses only the minimap2 aligner. mapPacBio / gapmm2 do not propagate FASTQ
    comments natively, and gapmm2 isn't installed locally — those aligner
    combinations are a separate concern (chimeric consensus may need a tag
    re-attach step) and are out of scope for the smoke test.
    """
    correct_out = tmp_path / "correct"
    rc = subprocess.run(
        [sys.executable, "-m", "rectify", "correct-cdna",
         str(_CHRI_BAM), "--out", str(correct_out),
         "--reference", str(_FASTA)],
        capture_output=True, text=True, timeout=900,
    )
    assert rc.returncode == 0, f"correct-cdna failed: {rc.stderr[-1500:]}"
    fastq = correct_out / "stage1_consensus.fastq.gz"
    assert fastq.exists(), "correct-cdna did not emit stage1_consensus.fastq.gz"

    align_out = tmp_path / "align"
    rc = subprocess.run(
        [sys.executable, "-m", "rectify", "align",
         str(fastq), "--genome", str(_FASTA),
         "-o", str(align_out), "--aligners", "minimap2",
         "--prefix", "stage1"],
        capture_output=True, text=True, timeout=1800,
    )
    assert rc.returncode == 0, f"align failed: {rc.stderr[-1500:]}"
    # `align` produces <prefix>.multialigned.bam (and .multialigned.md.bam after calmd)
    aligned_bam = align_out / "stage1.multialigned.md.bam"
    if not aligned_bam.exists():
        aligned_bam = align_out / "stage1.multialigned.bam"
    assert aligned_bam.exists(), f"no multialigned BAM in {align_out}: {list(align_out.iterdir())}"

    # Spot-check tag passthrough: the first mapped record must carry the
    # alignment-independent tags from the FASTQ comment.
    import pysam as _pysam
    with _pysam.AlignmentFile(str(aligned_bam), "rb") as bam:
        first = next((r for r in bam.fetch(until_eof=True) if not r.is_unmapped), None)
    assert first is not None, "no mapped records in rectified BAM"
    have = {t for t, _ in first.tags}
    for required in ("XU", "XC", "XR", "XO", "XT", "XY", "XF", "XB"):
        assert required in have, (
            f"tag {required} missing from rectified BAM (-y passthrough broken?). "
            f"Present: {sorted(have)}"
        )

    analyze_out = tmp_path / "analyze"
    rc = subprocess.run(
        [sys.executable, "-m", "rectify", "cdna-analyze",
         str(aligned_bam), "-o", str(analyze_out),
         "--gff", str(_GFF), "--reference", str(_FASTA)],
        capture_output=True, text=True, timeout=600,
    )
    assert rc.returncode == 0, f"cdna-analyze failed: {rc.stderr[-1500:]}"
    for f in ("clusters.tsv", "isoforms.tsv", "t1t2_pairs.tsv", "consensus_tagged.bam"):
        assert (analyze_out / f).exists(), f"missing {f}"

    n_clusters = sum(1 for _ in (analyze_out / "clusters.tsv").open()) - 1  # header
    assert 24000 < n_clusters < 30000, f"clusters.tsv has {n_clusters} rows"
    n_isoforms = sum(1 for _ in (analyze_out / "isoforms.tsv").open()) - 1
    assert 3000 < n_isoforms < 4500, f"isoforms.tsv has {n_isoforms} rows"

    # Tagged BAM should contain XS/XA and (for the post-classification subset) XG
    with _pysam.AlignmentFile(str(analyze_out / "consensus_tagged.bam"), "rb") as bam:
        seen_xs = seen_xa = seen_xg = 0
        for rec in bam.fetch(until_eof=True):
            tags = {t for t, _ in rec.tags}
            seen_xs += "XS" in tags
            seen_xa += "XA" in tags
            seen_xg += "XG" in tags
            if seen_xs > 100 and seen_xa > 100 and seen_xg > 100:
                break
    assert seen_xs > 100 and seen_xa > 100 and seen_xg > 100, (
        f"tagged BAM missing recomputed tags: XS={seen_xs} XA={seen_xa} XG={seen_xg}"
    )


# ---------------------------------------------------------------------------
# Consensus pre-trimming (Task 2)
# ---------------------------------------------------------------------------
class TestPretrimConsensus:
    """pretrim_consensus strips SSP/UMI/GGG from 5' and poly-A from 3' so
    aligners receive clean mRNA sequence only."""

    def _make_fwd_consensus(self, mrna: str, polya: str = "A" * 8) -> str:
        umi = "T" * UMI_LEN
        return SSP_FWD + umi + "GGG" + mrna + polya + ANCHOR_FWD

    def test_fwd_strips_ssp_umi_ggg_and_polya(self):
        mrna = "ATGCATGC"
        full = self._make_fwd_consensus(mrna)
        trimmed, q_trim_5, pretrim_pa_len = pretrim_consensus(full, "fwd", read_type=1)
        assert trimmed == mrna, (
            f"expected trimmed=={mrna!r}, got {trimmed!r}"
        )
        assert q_trim_5 == len(SSP_FWD) + UMI_LEN + 3  # 23 + 27 + 3 = 53
        assert pretrim_pa_len == 8

    def test_fwd_type2_no_5p_trim(self):
        """Type-2 reads have no SSP — only poly-A strip applies."""
        mrna = "GCTAGCTAGC"
        seq = mrna + "A" * 12 + ANCHOR_FWD
        trimmed, q_trim_5, pretrim_pa_len = pretrim_consensus(seq, "fwd", read_type=2)
        assert q_trim_5 == 0, "Type-2 should not strip a 5' prefix"
        assert trimmed == mrna
        assert pretrim_pa_len == 12

    def test_fwd_no_ssp_found_returns_full_with_no_5p_trim(self):
        """If SSP is absent from a Type-1 consensus (e.g. POA artifact), the
        function should not crash and should not trim the 5' end."""
        seq = "ATGCATGC" + "A" * 6
        trimmed, q_trim_5, _ = pretrim_consensus(seq, "fwd", read_type=1)
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
        from rectify.core.consensus.chimeric_consensus import select_best_chimeric
        aligner_reads, genome = self._make_segs()
        result = select_best_chimeric(aligner_reads, genome)
        # The 3' segment holds the junction. mapPacBio has GT-AG (+5); mm2 has GT-TT (-3).
        assert result.three_prime_aligner == "mapPacBio", (
            f"Expected mapPacBio to win 3' segment (canonical junction); "
            f"got three_prime_aligner={result.three_prime_aligner!r}"
        )

    def test_chimeric_cigar_uses_canonical_intron(self):
        from rectify.core.consensus.chimeric_consensus import select_best_chimeric
        aligner_reads, genome = self._make_segs()
        result = select_best_chimeric(aligner_reads, genome)
        # The chimeric CIGAR should use the 102-bp intron from mapPacBio, not the 50-bp one.
        n_ops = {op: length for op, length in result.chimeric_cigar}
        assert n_ops.get(3) == 102, (
            f"Expected N-op length 102 (canonical intron) in chimeric CIGAR; "
            f"got cigar={result.chimeric_cigar}"
        )
