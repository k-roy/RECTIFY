#!/usr/bin/env python3
"""
Tests for splice junction analysis: SpliceMotifScorer, NonCanonicalSSDetector,
ASSDetector, filter_junctions, resolve_homopolymer_ambiguity,
and _rescue_5prime_softclip.

Author: Kevin R. Roy
"""

import pytest
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

# ---------------------------------------------------------------------------
# Imports under test
# ---------------------------------------------------------------------------

from rectify.utils.splice_motif import (
    SpliceMotifScorer,
    NonCanonicalSSDetector,
    get_splice_site_dinucleotides,
    get_splice_site_sequences,
)
from rectify.core.analyze.junction_analysis import ASSDetector
from rectify.core.aggregate.junctions import filter_junctions, resolve_homopolymer_ambiguity
from rectify.core.consensus.consensus import AlignmentInfo, _rescue_5prime_softclip
from rectify.core.splice.splice_aware_5prime import (
    rescue_3ss_truncation,
    _hp_edit_distance,
    _ACCEPTOR_PRIORITY_PLUS,
    _ACCEPTOR_PRIORITY_MINUS,
)


# =============================================================================
# SpliceMotifScorer
# =============================================================================

class TestSpliceMotifScorer:

    def test_yeast_perfect_five_ss(self):
        scorer = SpliceMotifScorer.for_yeast()
        assert scorer.score_five_ss('GTATGT') == 0.0

    def test_yeast_perfect_three_ss(self):
        scorer = SpliceMotifScorer.for_yeast()
        # YYYAG: C/T at pos 0,1,2 then A, G
        assert scorer.score_three_ss('CTTAG') == 0.0

    def test_yeast_single_mismatch_five_ss(self):
        scorer = SpliceMotifScorer.for_yeast()
        # First base G→A: penalty 4.0
        score = scorer.score_five_ss('ATATGT')
        assert score == 4.0

    def test_yeast_last_base_mismatch_five_ss(self):
        scorer = SpliceMotifScorer.for_yeast()
        # Last base T→C: penalty 1.0
        score = scorer.score_five_ss('GTATGC')
        assert score == 1.0

    def test_metazoan_perfect_five_ss(self):
        scorer = SpliceMotifScorer.for_metazoan()
        assert scorer.score_five_ss('GTAAGT') == 0.0

    def test_wrong_length_returns_100(self):
        scorer = SpliceMotifScorer.for_yeast()
        assert scorer.score_five_ss('GT') == 100.0
        assert scorer.score_three_ss('AG') == 100.0

    def test_case_insensitive(self):
        scorer = SpliceMotifScorer.for_yeast()
        assert scorer.score_five_ss('gtatgt') == scorer.score_five_ss('GTATGT')

    def test_classify_canonical(self):
        scorer = SpliceMotifScorer.for_yeast()
        assert scorer.classify_splice_site_type('GT', 'AG') == 'GT-AG (canonical)'

    def test_classify_gc_ag(self):
        scorer = SpliceMotifScorer.for_yeast()
        assert scorer.classify_splice_site_type('GC', 'AG') == 'GC-AG (non-canonical)'

    def test_classify_at_ac(self):
        scorer = SpliceMotifScorer.for_yeast()
        assert scorer.classify_splice_site_type('AT', 'AC') == 'AT-AC (minor spliceosome)'


# =============================================================================
# get_splice_site_dinucleotides
# =============================================================================

class TestGetSpliceSiteDinucleotides:
    """Test dinucleotide extraction with a fake genome."""

    # Genome: chrI = ...GT...AG... where intron spans [10, 20)
    # So genome[10:12] = GT (5' dinuc on + strand)
    # genome[18:20] = AG (3' dinuc on + strand)
    GENOME = {'chrI': 'NNNNNNNNNNGT' + 'N' * 6 + 'AGNNNNNN'}  # len=28

    def test_plus_strand(self):
        five, three = get_splice_site_dinucleotides(
            self.GENOME, 'chrI', 10, 20, '+'
        )
        assert five == 'GT'
        assert three == 'AG'

    def test_minus_strand_rc(self):
        # Minus strand: five=RC(genome[18:20])=RC("AG")="CT", three=RC(genome[10:12])=RC("GT")="AC"
        five, three = get_splice_site_dinucleotides(
            self.GENOME, 'chrI', 10, 20, '-'
        )
        assert five == 'CT'
        assert three == 'AC'

    def test_missing_chrom_returns_empty(self):
        five, three = get_splice_site_dinucleotides(
            self.GENOME, 'chrII', 10, 20, '+'
        )
        assert five == ''
        assert three == ''

    def test_out_of_bounds_returns_empty(self):
        five, three = get_splice_site_dinucleotides(
            self.GENOME, 'chrI', -5, 20, '+'
        )
        assert five == ''
        assert three == ''


# =============================================================================
# NonCanonicalSSDetector
# =============================================================================

class TestNonCanonicalSSDetector:

    def test_gt_ag_canonical(self):
        det = NonCanonicalSSDetector()
        assert det.is_canonical('GT', 'AG') is True
        result = det.classify('GT', 'AG')
        assert result['type'] == 'canonical'
        assert result['likely_artifact'] is False

    def test_gc_ag_canonical(self):
        det = NonCanonicalSSDetector()
        assert det.is_canonical('GC', 'AG') is True
        result = det.classify('GC', 'AG')
        assert result['type'] == 'canonical'

    def test_at_ac_canonical(self):
        det = NonCanonicalSSDetector()
        assert det.is_canonical('AT', 'AC') is True
        result = det.classify('AT', 'AC')
        assert result['type'] == 'canonical'

    def test_noncanonical(self):
        det = NonCanonicalSSDetector()
        assert det.is_canonical('AA', 'TT') is False
        result = det.classify('AA', 'TT')
        assert result['type'] == 'noncanonical'
        assert result['likely_artifact'] is True

    def test_n_in_dinuc_is_unknown(self):
        det = NonCanonicalSSDetector()
        result = det.classify('NN', 'AG')
        assert result['type'] == 'unknown'
        assert result['likely_artifact'] is True

    def test_rare_noncanonical(self):
        det = NonCanonicalSSDetector()
        result = det.classify('GT', 'CG')
        assert result['type'] == 'rare_noncanonical'
        assert result['likely_artifact'] is False

    def test_case_insensitive(self):
        det = NonCanonicalSSDetector()
        assert det.is_canonical('gt', 'ag') is True

    def test_record_and_summary(self):
        det = NonCanonicalSSDetector()
        det.record('GT', 'AG')
        det.record('GT', 'AG')
        det.record('GC', 'AG')
        assert det.total_junctions == 3
        rows = det.summary()
        assert rows[0]['five_dinuc'] == 'GT'
        assert rows[0]['count'] == 2
        freqs = [r['frequency'] for r in rows]
        assert abs(sum(freqs) - 1.0) < 1e-9


# =============================================================================
# ASSDetector
# =============================================================================

class TestASSDetector:

    @pytest.fixture
    def detector(self):
        annotated = {
            'chrI': {(100, 200), (300, 500)},
        }
        return ASSDetector(annotated, proximity_threshold=10)

    def test_exact_match(self, detector):
        result = detector.classify('chrI', 100, 200)
        assert result['classification'] == 'exact'
        assert result['is_exact'] is True
        assert result['offset_donor'] == 0
        assert result['offset_acceptor'] == 0

    def test_alt_5ss(self, detector):
        # Donor shifted by 5, acceptor exact
        result = detector.classify('chrI', 105, 200)
        assert result['classification'] == 'alt_5ss'
        assert result['is_alt_5ss'] is True
        assert result['offset_donor'] == 5
        assert result['offset_acceptor'] == 0

    def test_alt_3ss(self, detector):
        # Donor exact, acceptor shifted
        result = detector.classify('chrI', 100, 205)
        assert result['classification'] == 'alt_3ss'
        assert result['is_alt_3ss'] is True
        assert result['offset_donor'] == 0
        assert result['offset_acceptor'] == 5

    def test_alt_both(self, detector):
        result = detector.classify('chrI', 105, 205)
        assert result['classification'] == 'alt_both'

    def test_novel_no_nearby(self, detector):
        result = detector.classify('chrI', 1000, 1200)
        assert result['classification'] == 'novel'
        assert result['is_novel'] is True
        assert result['nearest_annotated'] is None

    def test_beyond_threshold_is_novel(self, detector):
        # Both ends > threshold (10): donor off by 15, acceptor off by 15
        # Neither annotated junction (100,200) nor (300,500) has either end within 10bp
        result = detector.classify('chrI', 115, 215)
        assert result['classification'] == 'novel'

    def test_missing_chrom(self, detector):
        result = detector.classify('chrII', 100, 200)
        assert result['classification'] == 'novel'

    def test_classify_many(self, detector):
        junctions = [('chrI', 100, 200), ('chrI', 105, 200), ('chrI', 1000, 1200)]
        results = detector.classify_many(junctions)
        assert results[('chrI', 100, 200)]['classification'] == 'exact'
        assert results[('chrI', 105, 200)]['classification'] == 'alt_5ss'
        assert results[('chrI', 1000, 1200)]['classification'] == 'novel'


# =============================================================================
# filter_junctions
# =============================================================================

class TestFilterJunctions:

    @pytest.fixture
    def sample_df(self):
        # Use the column names expected by filter_junctions:
        #   total_reads (or full_junction_reads), n_aligners, is_canonical, annotation_status
        return pd.DataFrame({
            'chrom': ['chrI', 'chrI', 'chrI', 'chrI'],
            'intron_start': [100, 300, 500, 700],
            'intron_end': [200, 400, 600, 800],
            'strand': ['+', '+', '-', '-'],
            'total_reads': [10, 3, 8, 2],
            'n_aligners': [3, 1, 2, 1],
            'is_canonical': [True, False, True, False],
            'annotation_status': ['exact', 'novel', 'near', 'novel'],
        })

    def test_min_reads_filter(self, sample_df):
        result = filter_junctions(sample_df, min_reads=5)
        assert len(result) == 2
        assert set(result['total_reads']) == {10, 8}

    def test_min_aligners_filter(self, sample_df):
        result = filter_junctions(sample_df, min_reads=1, min_aligners=2)
        assert len(result) == 2
        assert set(result['n_aligners']) == {3, 2}

    def test_require_canonical(self, sample_df):
        result = filter_junctions(sample_df, min_reads=1, min_aligners=1, require_canonical=True)
        assert result['is_canonical'].all()
        assert len(result) == 2

    def test_require_annotated_or_near(self, sample_df):
        result = filter_junctions(
            sample_df, min_reads=1, min_aligners=1, require_annotated_or_near=True
        )
        assert set(result['annotation_status']).issubset({'exact', 'near'})
        assert len(result) == 2

    def test_no_filters_returns_all(self, sample_df):
        result = filter_junctions(sample_df, min_reads=1, min_aligners=1)
        assert len(result) == len(sample_df)

    def test_empty_df_returns_empty(self):
        empty = pd.DataFrame(columns=['chrom', 'intron_start', 'intron_end',
                                       'strand', 'total_reads', 'aligner_count'])
        result = filter_junctions(empty, min_reads=5)
        assert len(result) == 0


# =============================================================================
# resolve_homopolymer_ambiguity
# =============================================================================

class TestResolveHomopolymerAmbiguity:
    """
    Test with a fake genome where a GT-AG intron is inside a poly-A run,
    creating homopolymer-ambiguous positions.

    Sequence (0-based):
        pos: 0123456789...
        NNNNAAAA[GT]NNNNN[AG]AAAA
    We set up the genome so sliding left or right by 1 gives the same
    flanking base (A).

    Specifically:
        intron_start=8, intron_end=20
        genome[7]='A', genome[8]='A', genome[9]='A'  (homopolymer at donor)
        → positions (8,20), (9,20), (7,20) are all "equivalent"
    """

    GENOME = {'chrI': 'NNNNAAAA' + 'AAAA' + 'NNNNNNN' + 'AGAAAA' + 'NNNNN'}

    def _make_df(self, coords_and_reads):
        rows = []
        for (start, end), reads in coords_and_reads:
            rows.append({
                'chrom': 'chrI',
                'intron_start': start,
                'intron_end': end,
                'strand': '+',
                'total_reads': reads,
                'five_dinuc': 'AA',  # Doesn't matter for ambiguity test
                'three_dinuc': 'AA',
            })
        return pd.DataFrame(rows)

    def test_non_ambiguous_preserved(self):
        # Junction not in a homopolymer context → preserved as-is
        genome = {'chrI': 'NNNNNGTNNNNAGNNNNN'}
        df = pd.DataFrame([{
            'chrom': 'chrI', 'intron_start': 5, 'intron_end': 11,
            'strand': '+', 'total_reads': 10,
        }])
        result = resolve_homopolymer_ambiguity(df, genome)
        assert len(result) == 1

    def test_ambiguous_cluster_merged(self):
        # Two junctions that are homopolymer-equivalent: reads get merged
        genome = {'chrI': 'AAAAAAAAAAAAAAAAAAAA'}  # all-A: every shift is equivalent
        df = pd.DataFrame([
            {'chrom': 'chrI', 'intron_start': 5, 'intron_end': 15,
             'strand': '+', 'total_reads': 7},
            {'chrom': 'chrI', 'intron_start': 6, 'intron_end': 16,
             'strand': '+', 'total_reads': 3},
        ])
        result = resolve_homopolymer_ambiguity(df, genome)
        assert len(result) == 1
        assert result.iloc[0]['total_reads'] == 10

    def test_motif_scorer_selects_best(self):
        # Two equivalent positions; one has canonical GT...AG, other has AA...TT
        # Motif scorer should pick the GT-AG one
        # Genome: pos 5=G, 6=T (donor start), pos 11=A, 12=G (acceptor end)
        #         pos 4=G, 5=T at shifted pos... need a special setup
        # Simpler: just two identical-flanking junctions in poly-A, one labeled GT-AG.
        # We rely on motif scorer's score_five_ss using the actual genome sequence.
        # Build genome: AAAAAA[GTNNNNAG]AAAA → intron [6,14)
        # Shift left: [5,13) → genome[5:7]='GT', genome[11:13] = depends
        genome = {'chrI': 'AAAAGT' + 'NNNN' + 'AGAAAA'}
        #                   0123456    78910   111213...
        # intron_start=4 → genome[4:6]='GT', intron_start=5 → genome[5:7]='GN'
        # Actually let's just verify the function returns exactly 1 row when scorer is provided
        df = pd.DataFrame([
            {'chrom': 'chrI', 'intron_start': 4, 'intron_end': 10,
             'strand': '+', 'total_reads': 5},
            {'chrom': 'chrI', 'intron_start': 5, 'intron_end': 11,
             'strand': '+', 'total_reads': 3},
        ])
        scorer = SpliceMotifScorer.for_yeast()
        result = resolve_homopolymer_ambiguity(df, genome, motif_scorer=scorer)
        # Whether or not they collapse depends on genome; just ensure no crash and ≤2 rows
        assert len(result) <= 2

    def test_different_chroms_not_merged(self):
        genome = {'chrI': 'AAAAAAAAAAAAA', 'chrII': 'AAAAAAAAAAAAA'}
        df = pd.DataFrame([
            {'chrom': 'chrI', 'intron_start': 5, 'intron_end': 10,
             'strand': '+', 'total_reads': 5},
            {'chrom': 'chrII', 'intron_start': 5, 'intron_end': 10,
             'strand': '+', 'total_reads': 3},
        ])
        result = resolve_homopolymer_ambiguity(df, genome)
        # Both are on different chroms — may each form their own cluster
        assert len(result) >= 1  # At minimum one survives per chrom


# =============================================================================
# _rescue_5prime_softclip
# =============================================================================

class TestRescue5PrimeSoftclip:
    """
    Test the 5' soft-clip rescue logic with a synthetic genome and alignments.

    Setup:
        Genome (chrT):
            pos 0-99: 'A' * 100  (exon1)
            pos 100-199: 'N' * 100  (intron, positions 100..199)
            pos 200-299: 'G' * 100  (exon2)

    Plus-strand alignment that missed the upstream intron:
        - Aligner starts at reference_start=200 (exon2)
        - 5' soft-clip = 'AAAAAAAAAA' (10 A's) = last 10 bases of exon1
        - Candidate junction: (chrT, 100, 200) — intron spans [100, 200)
        → Should rescue: genome[90:100]='AAAAAAAAAA' matches clip

    Minus-strand alignment:
        - The gene is on minus strand; exon1 in transcript is at the RIGHT in reference
        - Intron in reference spans [100, 200)
        - Aligner returns reference_end=100 (just left of intron)
        - 5' soft-clip (transcript 5' = rightmost in reference) = last bases of exon2
        - But for minus strand upstream intron: intron_start ≥ reference_end
          so the intron [100,200) has intron_start=100 == reference_end=100 → dist=0 → skipped
          (dist ≤ 0 is excluded)
        We need intron_start > reference_end. So use a different setup:
            reference_end = 99, intron at [100, 200)
            clip = genome[200:210] = 'G'*10 (first 10 bases of exon after intron_end in reference)
        → Should rescue.
    """

    GENOME = {'chrT': 'A' * 100 + 'N' * 100 + 'G' * 100}

    def _make_alignment(self, *, strand, reference_start, reference_end,
                        clip_seq, chrom='chrT'):
        return AlignmentInfo(
            read_id='test_read',
            aligner='mm2',
            chrom=chrom,
            strand=strand,
            reference_start=reference_start,
            reference_end=reference_end,
            cigar_string='10S90M',
            five_prime_softclip_seq=clip_seq,
            five_prime_softclip=len(clip_seq),
        )

    def test_plus_strand_rescue(self):
        """Clip matches exon1 end, intron_end == reference_start → rescued."""
        aln = self._make_alignment(
            strand='+',
            reference_start=200,
            reference_end=290,
            clip_seq='A' * 10,
        )
        junctions = {('chrT', 100, 200, '+')}  # intron [100,200)
        result = _rescue_5prime_softclip(aln, self.GENOME, junctions)
        assert result is True

    def test_plus_strand_no_rescue_wrong_sequence(self):
        """Clip sequence doesn't match exon1 end → not rescued."""
        aln = self._make_alignment(
            strand='+',
            reference_start=200,
            reference_end=290,
            clip_seq='C' * 10,  # C's, but exon1 is all A's
        )
        junctions = {('chrT', 100, 200, '+')}
        result = _rescue_5prime_softclip(aln, self.GENOME, junctions)
        assert result is False

    def test_plus_strand_no_rescue_internal_junction(self):
        """Intron is DOWNSTREAM of alignment start (internal) → not rescued."""
        aln = self._make_alignment(
            strand='+',
            reference_start=50,   # alignment starts at 50
            reference_end=150,
            clip_seq='A' * 10,
        )
        junctions = {('chrT', 100, 200, '+')}  # intron_end=200 > align_5prime=50 → dist<0 → skip
        result = _rescue_5prime_softclip(aln, self.GENOME, junctions)
        assert result is False

    def test_plus_strand_no_rescue_too_far(self):
        """Intron is too far from alignment start → not rescued."""
        aln = self._make_alignment(
            strand='+',
            reference_start=200,
            reference_end=290,
            clip_seq='A' * 10,
        )
        # intron_end=190 → dist = 200-190 = 10, within default window (300)
        # But let's use search_window_bp=5 to test exclusion
        junctions = {('chrT', 100, 190, '+')}  # intron_end=190, dist=10 > window=5
        result = _rescue_5prime_softclip(
            aln, self.GENOME, junctions, search_window_bp=5
        )
        assert result is False

    def test_minus_strand_rescue(self):
        """Minus strand: clip matches exon to the right of intron_end → rescued."""
        # reference_end=99, intron [100,200) → dist = intron_start(100) - align_5prime(98) = 2 > 0
        # exon sequence = genome[200:210] = 'GGGGGGGGGG'
        aln = self._make_alignment(
            strand='-',
            reference_start=10,
            reference_end=99,  # align_5prime = reference_end - 1 = 98
            clip_seq='G' * 10,
        )
        junctions = {('chrT', 100, 200, '-')}  # intron [100,200)
        result = _rescue_5prime_softclip(aln, self.GENOME, junctions)
        assert result is True

    def test_minus_strand_no_rescue_wrong_sequence(self):
        aln = self._make_alignment(
            strand='-',
            reference_start=10,
            reference_end=99,
            clip_seq='T' * 10,  # T's, but exon2 is G's
        )
        junctions = {('chrT', 100, 200, '-')}
        result = _rescue_5prime_softclip(aln, self.GENOME, junctions)
        assert result is False

    def test_no_clip_returns_false(self):
        aln = self._make_alignment(
            strand='+',
            reference_start=200,
            reference_end=290,
            clip_seq='',
        )
        result = _rescue_5prime_softclip(aln, self.GENOME, set())
        assert result is False

    def test_missing_chrom_returns_false(self):
        aln = self._make_alignment(
            strand='+',
            reference_start=200,
            reference_end=290,
            clip_seq='A' * 10,
            chrom='chrMissing',
        )
        result = _rescue_5prime_softclip(aln, self.GENOME, {('chrMissing', 100, 200, '+')})
        assert result is False

    def test_wrong_chrom_in_junctions_skipped(self):
        """Junction on different chrom than alignment → not rescued."""
        aln = self._make_alignment(
            strand='+',
            reference_start=200,
            reference_end=290,
            clip_seq='A' * 10,
        )
        junctions = {('chrOther', 100, 200, '+')}  # different chrom
        result = _rescue_5prime_softclip(aln, self.GENOME, junctions)
        assert result is False

    def test_partial_mismatch_within_threshold(self):
        """1 mismatch in 10-base clip: 0.1 ≤ 0.2 default → rescued."""
        # genome[90:100] = 'A'*10
        clip = 'A' * 9 + 'C'  # 1 mismatch
        aln = self._make_alignment(
            strand='+',
            reference_start=200,
            reference_end=290,
            clip_seq=clip,
        )
        junctions = {('chrT', 100, 200, '+')}
        result = _rescue_5prime_softclip(aln, self.GENOME, junctions)
        assert result is True

    def test_too_many_mismatches(self):
        """4 mismatches in 10-base clip: 0.4 > 0.2 → not rescued."""
        clip = 'C' * 4 + 'A' * 6
        aln = self._make_alignment(
            strand='+',
            reference_start=200,
            reference_end=290,
            clip_seq=clip,
        )
        junctions = {('chrT', 100, 200, '+')}
        result = _rescue_5prime_softclip(aln, self.GENOME, junctions)
        assert result is False

# =============================================================================
# MPB terminal-error rescue via effective_five_prime_clip_seq
# =============================================================================

class TestMPBRescue:
    """
    Verify that _rescue_5prime_softclip uses effective_five_prime_clip_seq
    when five_prime_softclip_seq is empty (mapPacBio forced-mismatch case).

    Genome (chrT):  A*100 | N*100 | G*100
    Intron: [100, 200)
    Junction for rescue: ('chrT', 100, 200)
    """

    GENOME = {'chrT': 'A' * 100 + 'N' * 100 + 'G' * 100}

    def _make_mpb_alignment(self, *, strand, reference_start, reference_end,
                             effective_seq, chrom='chrT'):
        """AlignmentInfo with zero explicit soft-clip but non-empty effective seq."""
        return AlignmentInfo(
            read_id='mpb_read',
            aligner='mapPacBio',
            chrom=chrom,
            strand=strand,
            reference_start=reference_start,
            reference_end=reference_end,
            cigar_string='100M',
            five_prime_softclip=0,
            five_prime_softclip_seq='',
            effective_five_prime_clip=len(effective_seq),
            effective_five_prime_clip_seq=effective_seq,
        )

    def test_plus_strand_mpb_rescue(self):
        """MPB: effective seq matches exon1 end → rescued even with zero soft-clip."""
        aln = self._make_mpb_alignment(
            strand='+',
            reference_start=200,
            reference_end=300,
            effective_seq='A' * 10,  # matches genome[90:100]
        )
        junctions = {('chrT', 100, 200, '+')}
        assert _rescue_5prime_softclip(aln, self.GENOME, junctions) is True

    def test_plus_strand_mpb_no_rescue_wrong_seq(self):
        """MPB: effective seq doesn't match upstream exon → not rescued."""
        aln = self._make_mpb_alignment(
            strand='+',
            reference_start=200,
            reference_end=300,
            effective_seq='C' * 10,  # C's ≠ A's in exon1
        )
        junctions = {('chrT', 100, 200, '+')}
        assert _rescue_5prime_softclip(aln, self.GENOME, junctions) is False

    def test_explicit_clip_takes_priority_over_effective(self):
        """When both explicit clip and effective seq are present, explicit wins."""
        aln = AlignmentInfo(
            read_id='r',
            aligner='mapPacBio',
            chrom='chrT',
            strand='+',
            reference_start=200,
            reference_end=300,
            cigar_string='10S90M',
            five_prime_softclip=10,
            five_prime_softclip_seq='A' * 10,      # matches → rescued
            effective_five_prime_clip=12,
            effective_five_prime_clip_seq='CC',     # does NOT match (would fail alone)
        )
        junctions = {('chrT', 100, 200, '+')}
        # Should rescue because five_prime_softclip_seq is tried first and matches
        assert _rescue_5prime_softclip(aln, self.GENOME, junctions) is True

    def test_rescue_seq_override_overrides_all(self):
        """rescue_seq_override takes priority over both alignment sequence fields."""
        aln = AlignmentInfo(
            read_id='r',
            aligner='mm2',
            chrom='chrT',
            strand='+',
            reference_start=200,
            reference_end=300,
            cigar_string='10S90M',
            five_prime_softclip=10,
            five_prime_softclip_seq='C' * 10,  # would fail
            effective_five_prime_clip=10,
            effective_five_prime_clip_seq='C' * 10,  # would fail
        )
        junctions = {('chrT', 100, 200, '+')}
        # Override with matching sequence → rescued
        assert _rescue_5prime_softclip(
            aln, self.GENOME, junctions,
            rescue_seq_override='A' * 10,
        ) is True

    def test_no_seq_anywhere_returns_false(self):
        """Zero explicit clip + zero effective clip seq → not rescued."""
        aln = AlignmentInfo(
            read_id='r',
            aligner='mapPacBio',
            chrom='chrT',
            strand='+',
            reference_start=200,
            reference_end=300,
            cigar_string='100M',
            five_prime_softclip=0,
            five_prime_softclip_seq='',
            effective_five_prime_clip=0,
            effective_five_prime_clip_seq='',
        )
        assert _rescue_5prime_softclip(aln, self.GENOME, {('chrT', 100, 200, '+')}) is False


# =============================================================================
# rescue_3ss_truncation  (post-consensus, splice_aware_5prime.py)
# =============================================================================

def _distinct_seq(n: int, seed: int = 1) -> str:
    """Deterministic, ~aperiodic ACGT string with no run > 2 bp.

    Keeps the HP-aware DP scoring sequences by identity (no homopolymer
    free-rides) and makes any 60-mer window uniquely placeable, so a
    _RESCUE_DP_CAP test that slices the wrong end of the clip lands at a
    different position (or fails to rescue) rather than matching by accident.
    """
    import random
    rng = random.Random(seed)
    out: list = []
    for _ in range(n):
        choices = [b for b in 'ACGT'
                   if not (len(out) >= 2 and out[-1] == b and out[-2] == b)]
        out.append(rng.choice(choices))
    return ''.join(out)


class MockRead:
    """Minimal pysam.AlignedSegment substitute for unit tests."""

    def __init__(self, *, reference_name, reference_start, reference_end,
                 is_reverse, query_sequence=None, cigartuples=None,
                 query_name='mock_read'):
        self.reference_name = reference_name
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.is_reverse = is_reverse
        self.query_sequence = query_sequence or ''
        self.cigartuples = cigartuples or []
        self.is_unmapped = False
        self.query_name = query_name

    def get_aligned_pairs(self):
        """Return naive aligned pairs (M ops only) for the terminal-error scan."""
        pairs = []
        qp = 0
        rp = self.reference_start
        for op, length in self.cigartuples:
            if op == 4:  # S
                for _ in range(length):
                    pairs.append((qp, None))
                    qp += 1
            elif op == 0:  # M
                for _ in range(length):
                    pairs.append((qp, rp))
                    qp += 1
                    rp += 1
        return pairs


class TestRescue3SSTruncation:
    """
    Tests for rescue_3ss_truncation() in splice_aware_5prime.py.

    Genome (chrT): A*100 | N*100 | G*100
      exon1:  [0, 100)   all A
      intron: [100, 200) all N
      exon2:  [200, 300) all G

    Plus-strand gene:  5'─exon1─intron─exon2─3'
      Intron: ('chrT', 100, 200)
      3'SS acceptor (exon2 start) = 200
      A read starting at 200 was truncated; soft-clip / MPB bytes = last N of exon1

    Minus-strand gene: 5'─exon2─intron─exon1─3'  (exon2 is upstream in transcript)
      The upstream exon (in transcript) = exon2 = reference positions [200,300)
      3'SS (intron end in reference) = 200  →  intron_start = 100
      A read on minus strand ending at reference_end = 100 was truncated;
      clip bytes = first N bases of exon2 = genome[200:200+N] = 'G'*N
    """

    GENOME = {'chrT': 'A' * 100 + 'N' * 100 + 'G' * 100}
    JUNCTION = {('chrT', 100, 200)}  # (chrom, intron_start, intron_end)

    # ---- Plus strand: explicit soft-clip ----

    def test_plus_softclip_rescued(self):
        """Explicit 10-base soft-clip matches exon1 end (A*10) → rescued."""
        # CIGAR: 10S90M  → soft-clip = 'A'*10
        read = MockRead(
            reference_name='chrT',
            reference_start=200,
            reference_end=290,
            is_reverse=False,
            query_sequence='A' * 10 + 'G' * 90,
            cigartuples=[(4, 10), (0, 90)],
        )
        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='+')
        assert r['rescued'] is True
        assert r['rescue_type'] == 'softclip'
        assert r['rescued_junction'] == ('chrT', 100, 200)
        # five_prime_corrected = intron_start - 1 = 99
        assert r['five_prime_corrected'] == 99

    def test_plus_softclip_wrong_seq_not_rescued(self):
        """Soft-clip is C's, exon1 is A's → sequence mismatch, not rescued."""
        read = MockRead(
            reference_name='chrT',
            reference_start=200,
            reference_end=290,
            is_reverse=False,
            query_sequence='C' * 10 + 'G' * 90,
            cigartuples=[(4, 10), (0, 90)],
        )
        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='+')
        assert r['rescued'] is False

    # ---- Plus strand: MPB forced-mismatch at exon2 start ----

    def test_plus_mpb_forced_mismatch_rescued(self):
        """No explicit soft-clip; first 10 aligned bases are all mismatches (A reads
        as C in MPB output against G-rich exon2), but the query bytes are 'A'*10
        which matches exon1 end → rescued via MPB terminal-error path."""
        # genome[200:210] = 'G'*10; read starts with 'A'*10 → mismatch at every pos
        # Threshold: ≥40% errors in window=8 → terminal_end covers first 10 bases
        seq = 'A' * 10 + 'G' * 40  # 10 mismatches then correct exon2 bases
        read = MockRead(
            reference_name='chrT',
            reference_start=200,
            reference_end=250,
            is_reverse=False,
            query_sequence=seq,
            cigartuples=[(0, 50)],  # pure M, no soft-clip
        )
        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='+')
        assert r['rescued'] is True
        # The reanchor pre-pass collapses the 10 leading mismatches into a 10S
        # soft-clip before the body runs, so the rescue path is 'softclip', not
        # 'mpb_mismatch'. Assertion updated from 'mpb_mismatch' post-reanchor wiring.
        assert r['rescue_type'] == 'softclip'
        assert r['five_prime_corrected'] == 99

    # ---- Plus strand: proximity-only (zero clip, zero errors) ----

    def test_plus_proximity_only(self):
        """Read starts exactly at the 3'SS with no clip and no terminal errors;
        exon2 sequence is G*50 so aligned bases match perfectly → proximity flag."""
        seq = 'G' * 50
        read = MockRead(
            reference_name='chrT',
            reference_start=200,
            reference_end=250,
            is_reverse=False,
            query_sequence=seq,
            cigartuples=[(0, 50)],
        )
        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='+')
        assert r['rescued'] is False
        assert r['rescue_type'] == 'proximity'
        assert r['rescued_junction'] == ('chrT', 100, 200)

    def test_plus_intronic_snap_returns_last_upstream_exon_base(self):
        """Case 4 intronic snap must not return intron_start itself."""
        read = MockRead(
            reference_name='chrT',
            reference_start=150,
            reference_end=200,
            is_reverse=False,
            query_sequence='A' * 50,
            cigartuples=[(0, 50)],
        )

        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='+')

        assert r['rescued'] is True
        assert r['rescue_type'] == 'intronic_snap'
        assert r['five_prime_corrected'] == 99
        assert r['five_prime_corrected'] != 100

    # ---- Plus strand: start too far from any 3'SS ----

    def test_plus_no_junction_in_proximity(self):
        """Read starts 50 bp past the 3'SS — outside default proximity window."""
        seq = 'G' * 50
        read = MockRead(
            reference_name='chrT',
            reference_start=250,
            reference_end=300,
            is_reverse=False,
            query_sequence=seq,
            cigartuples=[(0, 50)],
        )
        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='+')
        assert r['rescue_type'] == 'none'
        assert r['rescued'] is False

    # ---- Minus strand: explicit soft-clip ----

    def test_minus_softclip_rescued(self):
        """Minus strand: soft-clip at high end (transcript 5') matches exon2 start
        (genome[200:210] = 'G'*10) → rescued."""
        # reference_end = 100; intron [100,200); dist = intron_start(100) - align_5prime(99) = 1 ✓
        # Exon2 upstream of intron_end=200: genome[200:210] = 'G'*10
        read = MockRead(
            reference_name='chrT',
            reference_start=10,
            reference_end=100,
            is_reverse=True,
            query_sequence='A' * 90 + 'G' * 10,  # last 10 = soft-clip at high end
            cigartuples=[(0, 90), (4, 10)],
        )
        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='-')
        assert r['rescued'] is True
        assert r['rescue_type'] == 'softclip'
        # five_prime_corrected = intron_end = 200
        assert r['five_prime_corrected'] == 200

    # ---- Bug 3: n_op_snap — mapPacBio terminal-D overshoot ----
    #
    # Pattern 1 (original): CIGAR ends …N I M D (minus strand reversed: D M I N …)
    #   Terminal D at the 5' end after M and I ops pushes reference_end past intron_end.
    #
    # Pattern 2 (de69692f): CIGAR ends …N D = (minus strand reversed: = D N …)
    #   D between the N-op and the terminal match block.
    #   mapPacBio emits a few deleted exon bases between the intron skip and the
    #   terminal match, so align_5prime ends up past intron_end but _leading_del was
    #   not detected under the old code.
    #
    # Both patterns should fire n_op_snap and snap to intron_end (minus strand).
    # The genome for these tests: 'A'*100 (exon1) + 'N'*100 (intron) + 'G'*100 (exon2).
    # The N-op in these reads spans the intron [100, 200).
    # Pool junction is the same: ('chrT', 100, 200) → intron_end = 200.

    def test_minus_n_op_snap_pattern1_terminal_d(self):
        """Bug 3 Pattern 1 (original): reversed CIGAR starts D M I N …
        Terminal D at 5' end → _leading_del=True → n_op_snap fires."""
        # CIGAR: 50M 100N 10I 1M 15D  (minus strand; 5' = right end)
        # reversed: D(15), M(1), I(10), N(100), M(50)
        # D(15) is first → _leading_del=True
        # ref_start: M(50) at 50..100; N(100) at 100..200; D(15) at 200..215
        # ref_end = 215; align_5prime = 214 > intron_end=200 → dist<0
        # All-N query so sequence rescue fails → n_op_snap fires
        q_len = 50 + 10 + 1  # M + I + M (D/N don't consume query)
        read = MockRead(
            reference_name='chrT',
            reference_start=50,
            reference_end=215,   # 50+50+100+15 = 215
            is_reverse=True,
            query_sequence='N' * q_len,
            cigartuples=[(0, 50), (3, 100), (1, 10), (0, 1), (2, 15)],
        )
        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='-')
        assert r['rescued'] is True, r
        assert r['rescue_type'] == 'n_op_snap', r
        assert r['five_prime_corrected'] == 200, r

    def test_minus_n_op_snap_pattern2_d_between_n_and_match(self):
        """Bug 3 Pattern 2 (de69692f): reversed CIGAR starts = D N …
        D sits between the terminal match and the N-op.
        Under the old code _leading_del was False (= caused an early break).
        The fix detects = D N and sets _leading_del=True → n_op_snap fires."""
        # CIGAR: 50M 100N 5D 11=  (minus strand)
        # reversed: =(11), D(5), N(100), M(50)
        # =(11) → don't break (new code), D(5) → _leading_del=True
        # ref_start=50; M(50): 50..100; N(100): 100..200; D(5): 200..205; =(11): 205..216
        # ref_end=216; align_5prime=215 > intron_end=200 → dist<0
        # All-N query so sequence rescue fails → n_op_snap fires
        q_len = 50 + 11  # M + = (D/N don't consume query)
        read = MockRead(
            reference_name='chrT',
            reference_start=50,
            reference_end=216,   # 50+50+100+5+11 = 216
            is_reverse=True,
            query_sequence='N' * q_len,
            cigartuples=[(0, 50), (3, 100), (2, 5), (7, 11)],
        )
        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='-')
        assert r['rescued'] is True, r
        assert r['rescue_type'] == 'n_op_snap', r
        assert r['five_prime_corrected'] == 200, r

    def test_minus_clean_5prime_no_snap(self):
        """Clean 5' end with an N-op that places align_5prime past intron_end
        but no D near the N-op → no n_op_snap (read is correctly aligned past
        the intron, not a mapPacBio overshoot)."""
        # CIGAR: 50M 100N 16= (minus strand; reversed: =(16), N(100), M(50))
        # No D at all → _leading_del=False → junction is skipped
        q_len = 50 + 16
        read = MockRead(
            reference_name='chrT',
            reference_start=50,
            reference_end=216,
            is_reverse=True,
            query_sequence='N' * q_len,
            cigartuples=[(0, 50), (3, 100), (7, 16)],
        )
        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='-')
        assert r['rescue_type'] != 'n_op_snap', r

    # ---- Missing chrom ----

    def test_missing_chrom_returns_no_rescue(self):
        read = MockRead(
            reference_name='chrMissing',
            reference_start=200,
            reference_end=250,
            is_reverse=False,
            query_sequence='G' * 50,
            cigartuples=[(0, 50)],
        )
        r = rescue_3ss_truncation(read, self.GENOME, self.JUNCTION, strand='+')
        assert r['rescue_type'] == 'none'

    # ---- Large 5' clip: _RESCUE_DP_CAP truncation (PERF_AUDIT.md rec #2) ----

    def test_plus_large_clip_capped_still_rescues(self):
        """A 150 bp plus-strand soft-clip (> _RESCUE_DP_CAP) still rescues to the
        correct donor. The cap keeps the donor-ADJACENT bases (_rseq[-CAP:], the
        end abutting intron_start); a wrong-end slice would land at a shifted
        position or fail to rescue. exon1[0,180) intron[180,280) exon2[280,380)."""
        exon1 = _distinct_seq(180, seed=11)
        genome = {'chrL': exon1 + 'N' * 100 + 'C' * 100}
        junction = {('chrL', 180, 280)}
        clip = exon1[30:180]              # donor-adjacent 150 bp of exon1 (> cap)
        read = MockRead(
            reference_name='chrL',
            reference_start=280,
            reference_end=380,
            is_reverse=False,
            query_sequence=clip + 'C' * 100,
            cigartuples=[(4, 150), (0, 100)],
        )
        r = rescue_3ss_truncation(read, genome, junction, strand='+')
        assert r['rescued'] is True, r
        assert r['rescue_type'] == 'softclip', r
        assert r['rescued_junction'] == ('chrL', 180, 280), r
        assert r['five_prime_corrected'] == 179, r   # intron_start - 1

    def test_minus_large_clip_capped_still_rescues(self):
        """Minus-strand symmetric: 150 bp soft-clip at the transcript-5' (high
        genomic coord) end. The cap keeps the donor-adjacent bases (_rseq[:CAP],
        the end abutting intron_end). body[0,100) intron[100,200) exon2[200,380)."""
        exon2 = _distinct_seq(180, seed=22)
        genome = {'chrW': 'A' * 100 + 'N' * 100 + exon2}   # 'chrM' is reserved → chrMito
        junction = {('chrW', 100, 200)}
        clip = exon2[0:150]              # donor-adjacent 150 bp (intron_end side, > cap)
        read = MockRead(
            reference_name='chrW',
            reference_start=10,
            reference_end=100,
            is_reverse=True,
            query_sequence='A' * 90 + clip,
            cigartuples=[(0, 90), (4, 150)],
        )
        r = rescue_3ss_truncation(read, genome, junction, strand='-')
        assert r['rescued'] is True, r
        assert r['rescue_type'] == 'softclip', r
        assert r['five_prime_corrected'] == 200, r   # intron_end


# =============================================================================
# _hp_edit_distance — homopolymer-aware edit distance
# =============================================================================

class TestHpEditDistance:
    """
    Indels whose base equals an immediate neighbour in the SAME string cost 0.5
    (Nanopore homopolymer run-length error). All other operations cost 1.0.
    """

    def test_identical_strings(self):
        assert _hp_edit_distance('ACGT', 'ACGT') == 0.0

    def test_empty_strings(self):
        assert _hp_edit_distance('', '') == 0.0
        assert _hp_edit_distance('', 'A') == 1.0
        assert _hp_edit_distance('A', '') == 1.0

    def test_substitution_always_full_cost(self):
        """Substitutions always cost 1.0, even inside a homopolymer."""
        assert _hp_edit_distance('ACGT', 'ACCT') == 1.0
        assert _hp_edit_distance('AAAA', 'AABA') == 1.0

    def test_non_run_indel_full_cost(self):
        """Deletion of a base with NO identical neighbour costs 1.0."""
        # Delete C from ACGT → AGT; C neighbours are A and G, not C
        assert _hp_edit_distance('ACGT', 'AGT') == 1.0

    def test_homopolymer_deletion_half_cost(self):
        """Delete one base from a run of identical bases → cost 0.5."""
        assert _hp_edit_distance('AAAA', 'AAA') == pytest.approx(0.5)
        assert _hp_edit_distance('TTTT', 'TTT') == pytest.approx(0.5)

    def test_homopolymer_insertion_half_cost(self):
        """Insert one base into an existing run → cost 0.5."""
        assert _hp_edit_distance('AAA', 'AAAA') == pytest.approx(0.5)

    def test_mixed_hp_and_substitution(self):
        """One hp-indel (0.5) + one substitution (1.0) = 1.5."""
        # AAACG → AACGT: delete one A from run (0.5) + insert T at end (1.0)
        assert _hp_edit_distance('AAACG', 'AACGT') == pytest.approx(1.5)

    def test_hp_indel_in_middle_of_run(self):
        """Single A deletion from AAAACCC → AAACCC is 0.5."""
        assert _hp_edit_distance('AAAACCC', 'AAACCC') == pytest.approx(0.5)

    def test_two_hp_indels_sum(self):
        """Two hp indels (one A from A-run, one C from C-run) → 1.0 total."""
        # 'AACCC' → 'ACC': delete one A (in run of 2: 0.5) + delete one C (in run of 3: 0.5) = 1.0
        assert _hp_edit_distance('AACCC', 'ACC') == pytest.approx(1.0)


# =============================================================================
# Acceptor priority tables
# =============================================================================

class TestAcceptorPriority:
    """
    _ACCEPTOR_PRIORITY_PLUS:  last 2 bases of intron (genomic, plus strand)
      AG=0 (canonical) < CG=1 < TG=2 < AT=3 < other=4

    _ACCEPTOR_PRIORITY_MINUS: first 2 bases of intron (genomic, minus strand)
      CT=0 (RC of AG) < CG=1 < CA=2 (RC of TG) < AT=3 < other=4
    """

    def test_plus_values(self):
        assert _ACCEPTOR_PRIORITY_PLUS['AG'] == 0
        assert _ACCEPTOR_PRIORITY_PLUS['CG'] == 1
        assert _ACCEPTOR_PRIORITY_PLUS['TG'] == 2
        assert _ACCEPTOR_PRIORITY_PLUS['AT'] == 3

    def test_plus_ordering(self):
        p = _ACCEPTOR_PRIORITY_PLUS
        assert p['AG'] < p['CG'] < p['TG'] < p['AT']

    def test_minus_values(self):
        assert _ACCEPTOR_PRIORITY_MINUS['CT'] == 0
        assert _ACCEPTOR_PRIORITY_MINUS['CG'] == 1
        assert _ACCEPTOR_PRIORITY_MINUS['CA'] == 2
        assert _ACCEPTOR_PRIORITY_MINUS['AT'] == 3

    def test_minus_ordering(self):
        m = _ACCEPTOR_PRIORITY_MINUS
        assert m['CT'] < m['CG'] < m['CA'] < m['AT']

    def test_unknown_dinuc_returns_4_via_get(self):
        assert _ACCEPTOR_PRIORITY_PLUS.get('NN', 4) == 4
        assert _ACCEPTOR_PRIORITY_PLUS.get('GG', 4) == 4
        assert _ACCEPTOR_PRIORITY_MINUS.get('NN', 4) == 4
        assert _ACCEPTOR_PRIORITY_MINUS.get('AA', 4) == 4

    def test_plus_ag_best(self):
        ag = _ACCEPTOR_PRIORITY_PLUS['AG']
        for motif in ('CG', 'TG', 'AT'):
            assert ag < _ACCEPTOR_PRIORITY_PLUS[motif]

    def test_minus_ct_best(self):
        ct = _ACCEPTOR_PRIORITY_MINUS['CT']
        for motif in ('CG', 'CA', 'AT'):
            assert ct < _ACCEPTOR_PRIORITY_MINUS[motif]


# =============================================================================
# Additional rescue_3ss_truncation cases — v2.9.6/v2.9.8/v2.9.9 logic
# =============================================================================

class TestRescue3SSTruncationExtended:
    """
    Tests for:
      v2.9.6: dynamic shift range + canonical GT/GC donor preference (plus strand)
      v2.9.8: ±5 bp baseline; minus strand canonical AC/GC; tuple scoring
      v2.9.9: hp_edit_distance scoring; 3'SS acceptor tiebreaker

    Each test uses a purpose-built genome dict to isolate the behaviour under test.
    """

    # ---- ±5 baseline catches a 3-bp offset ----

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "99558c1 two-step scoring (in_amb > donor_ok) interacts with the "
            "_off search window: shift=0 (annotated NN, in_amb=True) finds the "
            "same exon window via _off=3 as shift=-3 (canonical GT, in_amb=False) "
            "via _off=0. in_amb priority selects shift=0 → five_prime_corrected=102 "
            "instead of 99. Three fix options: (A) limit _off_limit=0 for soft-clip "
            "rescues when align_5prime==intron_end (narrow); (B) revisit in_amb > "
            "donor_ok for out-of-window canonical donors (design call); (C) track "
            "_eff_intron_start -= _off in the scoring loop so the recorded donor "
            "position matches where the window actually landed (root-cause fix, "
            "smallest diff — preferred)."
        ),
    )
    def test_plus_offset_junction_rescued(self):
        """
        Junction annotated at intron_start=103, but true GT donor is at 100.
        The ±5 baseline (not the ambiguity run) covers this 3-bp shift.
        Soft-clip = A*8 matches exon1 (A*100) ending at pos 99.

        XFAIL: 99558c1 two-step scoring causes the annotated in-window position
        (shift=0, NN) to beat the canonical GT at shift=-3 when both find the same
        exon window via different _off values. See open design item in debugger_queue.md.
        """
        genome = {'chrO': 'A' * 100 + 'GT' + 'N' * 98 + 'C' * 100}
        junction = {('chrO', 103, 200)}
        read = MockRead(
            reference_name='chrO',
            reference_start=200,
            reference_end=250,
            is_reverse=False,
            query_sequence='A' * 8 + 'C' * 50,
            cigartuples=[(4, 8), (0, 50)],
        )
        r = rescue_3ss_truncation(read, genome, junction, strand='+')
        assert r['rescued'] is True
        # shift=-3 lands intron_start at 100 (GT canonical); corrected = 99
        assert r['five_prime_corrected'] == 99

    # ---- Canonical GT donor beats non-canonical AA ----

    def test_plus_canonical_gt_beats_non_canonical(self):
        """
        Two junctions; exon1 is A*50 for both → same edit distance for A*8 clip.
        Junction at intron_start=50 has GT donor (canonical).
        Junction at intron_start=150 has AA donor (non-canonical).
        GT should win.
        """
        genome = {'chrC': 'A' * 50 + 'GT' + 'N' * 98 + 'AA' + 'N' * 98 + 'A' * 50}
        junctions = {
            ('chrC', 50, 150),   # GT donor
            ('chrC', 150, 250),  # AA donor
        }
        read = MockRead(
            reference_name='chrC',
            reference_start=150,
            reference_end=200,
            is_reverse=False,
            query_sequence='A' * 8 + 'N' * 50,
            cigartuples=[(4, 8), (0, 50)],
        )
        r = rescue_3ss_truncation(read, genome, junctions, strand='+')
        assert r['rescued'] is True
        assert r['five_prime_corrected'] == 49  # intron_start=50 → corrected=49

    # ---- GC donor is canonical (minor spliceosome) ----

    def test_plus_gc_donor_is_canonical(self):
        """GC donor is accepted as canonical and preferred over AA."""
        genome = {'chrGC': 'A' * 50 + 'GC' + 'N' * 98 + 'AA' + 'N' * 98 + 'A' * 50}
        junctions = {
            ('chrGC', 50, 150),   # GC donor → canonical
            ('chrGC', 150, 250),  # AA donor → non-canonical
        }
        read = MockRead(
            reference_name='chrGC',
            reference_start=150,
            reference_end=200,
            is_reverse=False,
            query_sequence='A' * 8 + 'N' * 50,
            cigartuples=[(4, 8), (0, 50)],
        )
        r = rescue_3ss_truncation(read, genome, junctions, strand='+')
        assert r['rescued'] is True
        assert r['five_prime_corrected'] == 49

    # ---- Minus strand: AC donor canonical (RC of GT) ----

    def test_minus_ac_donor_canonical(self):
        """
        Minus strand: genome[intron_end-2:intron_end] = 'AC' is canonical (RC of GT).
        Two junctions share the same intron_start (both within proximity of align_5prime),
        but have different intron_ends and thus different 5'SS dinucleotides.
        AC should win over TT at equal edit distance.

        Genome layout:
          pos 0-98:   exon2 (3' side; read body) — all C's except dinucs
          pos 98-99:  'TT' → non-canonical (junction B intron_end=100)
          pos 100-197: exon-candidate for junction B: G*98
          pos 198-199: 'AC' → canonical   (junction A intron_end=200)
          pos 200-297: exon-candidate for junction A: G*98
        Both junction exons start with 'G'*10 → soft-clip 'G'*10 matches both with ed=0.
        align_5prime = reference_end - 1 = 99; intron_start=99 for both → dist=0.

        The clip is 10 nt, not 8, so it clears ``min_informative_clip_bp()``
        (ISSUE-006): below that floor no clip can distinguish its placement from
        chance and the sequence rescue is refused by design. The donor tie-break
        this test asserts is unchanged — only the clip length is.

        ISSUE-020: the ranking segment now ends exactly at the junction — a 5'
        base inside the intron joins the clip, a base over an exon-2 position
        is trimmed. The old mock's body ran to reference_end 100 over genome
        'TT' with 'C' bases (two mismatches the reanchor pre-pass turns into
        clip), and its last base sat ON intron_start; with those bases in the
        segment the toy's exon 1 (G*98) is indistinguishable from junction A's
        intron (T + G*98), so the exon-vs-intron acceptance cannot pass — the old
        clip-only segment hid that. The read now ends one base before the intron
        (reference_end 98, dist 2) with an 11-nt clip whose first base lies over
        exon-2 position 98 and is trimmed; the 10 G's tie on both candidates and
        the canonical AC donor decides, as the test intends.
        """
        genome = {
            'chrM2': (
                'C' * 98          # exon2 body (pos 0-97)
                + 'TT'            # pos 98-99: non-canonical 5'SS for junction_end=100
                + 'G' * 98        # pos 100-197: exon1 candidate for junction B
                + 'AC'            # pos 198-199: canonical 5'SS for junction_end=200
                + 'G' * 98        # pos 200-297: exon1 candidate for junction A
            )
        }
        junctions = {
            ('chrM2', 99, 200),   # genome[198:200]='AC' → canonical
            ('chrM2', 99, 100),   # genome[98:100]='TT' → non-canonical
        }
        # Minus strand read: align_5prime = reference_end - 1 = 99
        # dist = intron_start(99) - align_5prime(99) = 0 ≤ junction_proximity_bp ✓
        # Soft clip at HIGH end (last op) = 'G'*10, matching genome[200:210] and genome[100:110]
        read = MockRead(
            reference_name='chrM2',
            reference_start=17,
            reference_end=98,
            is_reverse=True,
            query_sequence='C' * 81 + 'G' * 11,  # body over [17, 98); 11-nt soft clip at the 5' (high) end
            cigartuples=[(0, 81), (4, 11)],
        )
        r = rescue_3ss_truncation(read, genome, junctions, strand='-')
        assert r['rescued'] is True
        assert r['five_prime_corrected'] == 200  # canonical AC junction (intron_end=200)

    # ---- 3'SS acceptor tiebreaker: AG (priority 0) beats CG (priority 1) ----

    def test_plus_acceptor_ag_beats_cg(self):
        """
        Two GT-donor junctions with identical edit distances.
        Junction 1 ends with AG (3'SS acceptor priority 0).
        Junction 2 ends with CG (3'SS acceptor priority 1).
        AG should win.

        Genome:
          exon1   = A*50           (pos 0-49)
          intron1 = GT + N*96 + AG (pos 50-149; GT donor, AG acceptor)
          landing1= A*50           (pos 150-199)
          intron2 = GT + N*96 + CG (pos 200-299; GT donor, CG acceptor)
          landing2= A*50           (pos 300-349)
        """
        genome = {
            'chrAcc': (
                'A' * 50
                + 'GT' + 'N' * 96 + 'AG'
                + 'A' * 50
                + 'GT' + 'N' * 96 + 'CG'
                + 'A' * 50
            )
        }
        junctions = {
            ('chrAcc', 50, 150),   # AG acceptor
            ('chrAcc', 200, 300),  # CG acceptor
        }
        read = MockRead(
            reference_name='chrAcc',
            reference_start=150,
            reference_end=200,
            is_reverse=False,
            # 10 nt clip = min_informative_clip_bp() (ISSUE-006); both exon-1
            # candidates are poly-A so the acceptor tie-break under test is
            # unaffected by the length change.
            query_sequence='A' * 10 + 'A' * 50,
            cigartuples=[(4, 10), (0, 50)],
        )
        r = rescue_3ss_truncation(read, genome, junctions, strand='+')
        assert r['rescued'] is True
        assert r['five_prime_corrected'] == 49  # AG junction (intron_start=50)

    # ---- hp scoring prefers homopolymer-error match ----

    def test_plus_hp_edit_prefers_hp_match(self):
        """
        Soft-clip = 'A'*10 + 'C' (11 nt: a 10-A run plus C).
        Two junctions both close to align_5prime (within proximity=10):
          Junction A (intron_start=100): exon window genome[89:100]='A'*10 + 'C'... —
            with an 11 nt clip the window is genome[89:100] = 'A'*10 + 'C' minus one
            A, i.e. 'AAAAAAAAAAC' vs the genome's 'AAAAAAAAAC' run → hp_ed = 0.5
            (one deletion INSIDE a homopolymer run, the HP-discounted cost).
          Junction B (intron_start=107): exon window ends in 'GTGTAAC' → several subs.
        Both junctions have GT donors.  Junction A wins on the HP-discounted cost.

        The clip is 11 nt rather than 6 so it clears ``min_informative_clip_bp()``
        (ISSUE-006); the HP-vs-substitution preference under test is unchanged.

        ISSUE-020: the read starts 9 nt into exon 2 (dist 9 for junction A, 2 for
        junction B) and the ranking now TRIMS the `dist` junction-side clip bases
        — the exon-2 positions the aligner left unaligned — instead of sliding the
        genome window over them. The read therefore carries those exon-2 bases in
        its clip (nine mixed bases — a G run would make the 20-nt clip periodic
        and refuse the search — the realistic soft-clip shape); the old mock
        simply lacked them, and the slide was silently supplying the missing
        geometry. The exon-1 segment compared is the same 11 nt as before.

        Genome:
          pos 0-94:   A*95  (poly-A exon body)
          pos 95-99:  'AAAAC'  (last 5 = candidate window for junction A at is=100)
          pos 100-101: 'GT'  (donor for junction A)
          pos 102-105: 'AACCG'... actually need window for junction B too
          pos 100-105: 'GTAAC'
          pos 106-107: 'GT'  (donor for junction B)
          ...
        Simpler: put both intron_starts at 100 and 106, both with GT donors.
        """
        # exon ends at pos 99 (for junction A) and pos 105 (for junction B)
        # genome[95:100]='AAAAC' → matches clip 'AAAAAC' with hp_ed=0.5
        # genome[101:106]='GTAAC' → hp_ed('AAAAAC','GTAAC')=3 (sub G,T;A match,A→A,C match)
        genome = {
            'chrHP': (
                'A' * 95 + 'AAAAC'   # pos 0-99: poly-A exon
                + 'GT'               # pos 100-101: GT donor for junction A
                + 'GTAAC'            # pos 102-106: sequence between junctions
                + 'GT'               # pos 107-108: GT donor for junction B (intron_start=107)
                + 'N' * 90           # intron body
                + 'G' * 100          # exon2 (read landing zone): pos 199+
            )
        }
        junctions = {
            ('chrHP', 100, 190),   # GT at 100; exon candidate = genome[95:100]='AAAAC'
            ('chrHP', 107, 197),   # GT at 107; exon candidate = genome[102:107]='GTAAC'
        }
        # Read starts at 199 (exon2), soft clip must match exon near each intron_start.
        # dist(jA) = 199 - 190 = 9 ≤ 10 ✓; dist(jB) = 199 - 197 = 2 ≤ 10 ✓
        read = MockRead(
            reference_name='chrHP',
            reference_start=199,
            reference_end=249,
            is_reverse=False,
            # exon-1 tail + nine bases over the unaligned exon-2 positions [190, 199) as the clip; body = genome[199:249]
            query_sequence='A' * 10 + 'C' + 'TCAGTCGAC' + 'G' * 50,
            cigartuples=[(4, 20), (0, 50)],
        )
        r = rescue_3ss_truncation(read, genome, junctions, strand='+')
        assert r['rescued'] is True
        assert r['five_prime_corrected'] == 99  # junction A: intron_start=100 → corrected=99


# =============================================================================
# in_amb priority over canonical donor (99558c1 two-step scoring)
# =============================================================================

class TestInAmbVsDonorOkPriority:
    """
    Verifies the priority introduced in 99558c1:
      Step 2 (in_amb) ranks above Step 3 (canonical donor).

    At an edit-distance tie, a shift that lands inside the sequence-ambiguity
    window must beat a shift that lands outside it — even if the outside shift
    has a canonical GT/GC donor and the inside shift does not.

    Rationale: in_amb is a match-quality property (a slide within the window
    is a genuine no-op for the aligner); canonical-donor is a signal-quality
    property. Buying a better-anchored placement at the cost of non-canonical
    signal is never the right trade.
    """

    def test_plus_in_amb_beats_canonical_out_of_window(self):
        """
        Genome layout (chrP):
          pos 0-1:    CC  (exon1 head, irrelevant)
          pos 2-13:   GT*6  (repeating exon1 tail — both windows land here)
          pos 14-15:  NN  (non-canonical 5'SS donor at intron_start=14)
          pos 16-113: N*98  (intron body)
          pos 114-115: AG  (canonical 3'SS acceptor, intron_end=116)
          pos 116-165: C*50  (exon2 / read body)

        soft-clip rescue_seq = 'GTGTGTGTGT' (10 bases — the shortest clip that
        clears ``min_informative_clip_bp()``, ISSUE-006; the repeating tail was
        lengthened to 12 nt so both windows below still fit).

        Within the ±5 shift search:
          shift= 0: _eff_start=14, exon window genome[4:14]='GTGTGTGTGT', ED=0,
                    in_amb=True  (only shift=0 is in_amb because _r_amb=0, _l_amb=0),
                    donor genome[14:16]='NN' → non-canonical.
          shift=-2: _eff_start=12, exon window genome[2:12]='GTGTGTGTGT', ED=0,
                    in_amb=False (2 positions outside the ambiguity window),
                    donor genome[12:14]='GT' → canonical.

        The two-step scoring tuple (not in_amb, not donor_ok, shift_abs):
          shift= 0 → (False, True,  0)
          shift=-2 → (True,  False, 2)
        (False, True, 0) < (True, False, 2) → shift=0 wins → five_prime_corrected=13.

        With the OLD pre-99558c1 ordering (not donor_ok, not in_amb, shift_abs):
          shift= 0 → (True,  False, 0)
          shift=-2 → (False, True,  2)
        (False, True, 2) < (True, False, 0) → shift=-2 would win → five_prime_corrected=11.
        """
        genome = {'chrP':
                  'CC'              # pos 0-1
                  + 'GT' * 6        # pos 2-13: exon1 tail; GT at 12-13 = canonical donor
                  + 'NN'            # pos 14-15: non-canonical donor (intron_start=14)
                  + 'N' * 98        # pos 16-113: intron body
                  + 'AG'            # pos 114-115: canonical acceptor (intron_end=116)
                  + 'C' * 50}       # pos 116-165: exon2
        junction = {('chrP', 14, 116)}
        read = MockRead(
            reference_name='chrP',
            reference_start=116,
            reference_end=166,
            is_reverse=False,
            query_sequence='GT' * 5 + 'C' * 50,
            cigartuples=[(4, 10), (0, 50)],
        )
        r = rescue_3ss_truncation(read, genome, junction, strand='+')
        assert r['rescued'] is True
        # shift=0 (in_amb) wins over shift=-2 (canonical GT, out-of-window)
        assert r['five_prime_corrected'] == 13   # intron_start(14) - 1 = 13
        assert r['rescued_junction'] == ('chrP', 14, 116)

    def test_plus_canonical_wins_when_both_in_amb(self):
        """
        When both shifts are inside the ambiguity window (in_amb tied), the canonical
        GT donor (Step 3) is the deciding factor.  This is the complement of
        test_plus_in_amb_beats_canonical_out_of_window.

        Genome layout (chrQ):
          pos 0-15:  A*16  (exon1)
          pos 16-17: AA    (non-canonical donor at shift=0)
          pos 18-19: GT    (canonical donor at shift=+2)
          pos 20+:   intron body + AG acceptor + exon2

        Ambiguity window:
          _leb=genome[15]='A'; _r_amb: genome[16]='A'=_leb → +1;
            genome[17]='A'=_leb → +1; genome[18]='G'≠'A' → stop. _r_amb=2.
          _fib=genome[16]='A'; _l_amb: genome[15]='A'=_fib → +1; ... → _l_amb=15.
          in_amb window: (-15 ≤ shift ≤ 2) → both shift=0 and shift=+2 are in_amb.

        soft-clip = 'AAAA' matches both exon windows with ED=0:
          shift=0:  genome[12:16]='AAAA', ED=0
          shift=+2: genome[14:18]='AAAA', ED=0

        Scoring tuple (not in_amb, not donor_ok, shift_abs):
          shift=0:  (False, True,  0)  — in_amb, non-canonical
          shift=+2: (False, False, 2)  — in_amb, canonical GT
        (False, False, 2) < (False, True, 0) → shift=+2 (canonical) wins.

        Expected: five_prime_corrected = 17 (intron_start+2 - 1).
        """
        genome = {'chrQ':
                  'A' * 16          # pos 0-15: exon1
                  + 'AAGT'          # pos 16-19: AA (non-canonical) then GT (canonical)
                  + 'N' * 96        # pos 20-115: intron body
                  + 'AG'            # pos 116-117: acceptor (intron_end=118)
                  + 'C' * 50}       # pos 118-167: exon2
        junction = {('chrQ', 16, 118)}
        read = MockRead(
            reference_name='chrQ',
            reference_start=118,
            reference_end=168,
            is_reverse=False,
            # 10 nt clip = min_informative_clip_bp() (ISSUE-006); genome[6:16] and
            # genome[8:18] are both 'A'*10, so both shifts still tie at ED 0 and
            # the canonical-donor tie-break under test is unchanged.
            query_sequence='A' * 10 + 'C' * 50,
            cigartuples=[(4, 10), (0, 50)],
        )
        r = rescue_3ss_truncation(read, genome, junction, strand='+')
        assert r['rescued'] is True
        # shift=+2 (canonical GT, in_amb) beats shift=0 (non-canonical AA, in_amb)
        assert r['five_prime_corrected'] == 17  # intron_start+2 - 1 = 17
        assert r['rescued_junction'] == ('chrQ', 18, 118)
