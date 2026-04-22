#!/usr/bin/env python3
"""
Unit tests for consensus.py — select_best_alignment and related functions.

Covers:
- extract_junctions_from_cigar: CIGAR N-op parsing for plus and minus strand reads
- check_canonical_splice_sites: GT/AG motif detection in genome sequence
- score_alignment: 5' clip and A-tract depth penalty arithmetic
- select_best_alignment: winner selection, tiebreaker logic (3' agreement,
  annotated junction count, canonical motif count)
- Real-data tests using wt_by4742_rep1 BAM (skipped when not available)
"""

import gzip
import os
from dataclasses import replace
from pathlib import Path
from typing import Dict, List, Tuple

import pytest
import pysam

from rectify.core.consensus import (
    AlignmentInfo,
    ConsensusResult,
    check_canonical_splice_sites,
    extract_junctions_from_cigar,
    score_alignment,
    select_best_alignment,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_read(
    chrom: str,
    start: int,
    cigar: List[Tuple[int, int]],
    strand: str,
    seq: str,
) -> pysam.AlignedSegment:
    hdr = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6'},
        'SQ': [{'SN': chrom, 'LN': 3_000_000}],
    })
    read = pysam.AlignedSegment(hdr)
    read.query_name = 'test_read'
    read.reference_name = chrom
    read.reference_start = start
    read.cigartuples = cigar
    read.is_reverse = (strand == '-')
    read.is_unmapped = False
    read.is_secondary = False
    read.is_supplementary = False
    read.mapping_quality = 60
    read.query_sequence = seq
    return read


def _make_alignment_info(
    read_id: str = 'r1',
    aligner: str = 'minimap2',
    chrom: str = 'chrI',
    strand: str = '+',
    start: int = 1000,
    end: int = 2000,
    junctions: List[Tuple[int, int]] = None,
    five_prime_softclip: int = 0,
    three_prime_softclip: int = 0,
    effective_five_prime_clip: int = 0,
    effective_three_prime_clip: int = 0,
    five_prime_softclip_seq: str = '',
    three_prime_atract_depth: int = 0,
    corrected_3prime: int = None,
    canonical_count: int = 0,
) -> AlignmentInfo:
    if junctions is None:
        junctions = []
    return AlignmentInfo(
        read_id=read_id,
        aligner=aligner,
        chrom=chrom,
        strand=strand,
        reference_start=start,
        reference_end=end,
        cigar_string='100M',
        mapq=60,
        junctions=junctions,
        five_prime_softclip=five_prime_softclip,
        three_prime_softclip=three_prime_softclip,
        effective_five_prime_clip=effective_five_prime_clip,
        effective_three_prime_clip=effective_three_prime_clip,
        five_prime_softclip_seq=five_prime_softclip_seq,
        three_prime_atract_depth=three_prime_atract_depth,
        corrected_3prime=corrected_3prime if corrected_3prime is not None else end - 1,
        canonical_count=canonical_count,
    )


# ---------------------------------------------------------------------------
# extract_junctions_from_cigar
# ---------------------------------------------------------------------------

class TestExtractJunctionsFromCigar:

    def test_no_junctions_pure_match(self):
        read = _make_read('chrI', 1000, [(0, 100)], '+', 'A' * 100)
        assert extract_junctions_from_cigar(read) == []

    def test_single_intron_plus_strand(self):
        # 50M 200N 50M → intron at (1050, 1250)
        read = _make_read('chrI', 1000, [(0, 50), (3, 200), (0, 50)], '+', 'A' * 100)
        junctions = extract_junctions_from_cigar(read)
        assert len(junctions) == 1
        assert junctions[0] == (1050, 1250)

    def test_two_introns(self):
        # 30M 100N 30M 200N 30M
        read = _make_read('chrI', 500, [(0, 30), (3, 100), (0, 30), (3, 200), (0, 30)], '+', 'A' * 90)
        junctions = extract_junctions_from_cigar(read)
        assert junctions[0] == (530, 630)
        assert junctions[1] == (660, 860)

    def test_minus_strand_intron(self):
        # Minus strand read — junction coordinates are still genomic
        read = _make_read('chrI', 2000, [(0, 50), (3, 586), (0, 50)], '-', 'T' * 100)
        junctions = extract_junctions_from_cigar(read)
        assert len(junctions) == 1
        assert junctions[0] == (2050, 2636)

    def test_softclip_not_counted_as_junction(self):
        # 5S 50M 100N 50M 3S
        read = _make_read('chrI', 1000, [(4, 5), (0, 50), (3, 100), (0, 50), (4, 3)], '+', 'A' * 108)
        junctions = extract_junctions_from_cigar(read)
        assert len(junctions) == 1
        assert junctions[0] == (1050, 1150)

    def test_deletion_not_counted_as_junction(self):
        # Deletions (op=2) consume reference but are not junctions
        read = _make_read('chrI', 1000, [(0, 50), (2, 5), (0, 50)], '+', 'A' * 100)
        assert extract_junctions_from_cigar(read) == []

    def test_no_cigar_returns_empty(self):
        read = _make_read('chrI', 1000, [(0, 50)], '+', 'A' * 50)
        read.cigartuples = None
        assert extract_junctions_from_cigar(read) == []


# ---------------------------------------------------------------------------
# check_canonical_splice_sites
# ---------------------------------------------------------------------------

class TestCheckCanonicalSpliceSites:

    def test_gt_ag_is_canonical(self):
        # Build genome with GT..AG around a junction at (10, 20)
        seq = 'A' * 10 + 'GT' + 'A' * 6 + 'AG' + 'A' * 100
        genome = {'chrI': seq}
        canon, non_canon = check_canonical_splice_sites([(10, 20)], 'chrI', genome)
        assert canon == 1
        assert non_canon == 0

    def test_gc_ag_is_canonical(self):
        # GC at 5'SS is also in CANONICAL_5SS
        seq = 'A' * 10 + 'GC' + 'A' * 6 + 'AG' + 'A' * 100
        genome = {'chrI': seq}
        canon, non_canon = check_canonical_splice_sites([(10, 20)], 'chrI', genome)
        assert canon == 1
        assert non_canon == 0

    def test_non_canonical_detected(self):
        seq = 'A' * 10 + 'AT' + 'A' * 6 + 'TT' + 'A' * 100
        genome = {'chrI': seq}
        canon, non_canon = check_canonical_splice_sites([(10, 20)], 'chrI', genome)
        assert canon == 0
        assert non_canon == 1

    def test_mixed_canonical_and_non_canonical(self):
        # canonical junction at (10, 20):
        #   5'SS = seq[10:12] = 'GT'
        #   3'SS = seq[end-2:end] = seq[18:20] = 'AG'
        # non-canonical junction at (50, 62):
        #   5'SS = seq[50:52] = 'AT'
        #   3'SS = seq[60:62] = 'TT'
        # Layout: 10xA | GT | 6xA | AG | 30xA | AT | 8xA | TT | 38xA = 100 bp
        seq = 'A' * 10 + 'GT' + 'A' * 6 + 'AG' + 'A' * 30 + 'AT' + 'A' * 8 + 'TT' + 'A' * 38
        assert len(seq) == 100
        genome = {'chrI': seq}
        canon, non_canon = check_canonical_splice_sites([(10, 20), (50, 62)], 'chrI', genome)
        assert canon == 1
        assert non_canon == 1

    def test_empty_junctions(self):
        genome = {'chrI': 'A' * 100}
        canon, non_canon = check_canonical_splice_sites([], 'chrI', genome)
        assert canon == 0
        assert non_canon == 0

    def test_chrom_not_in_genome(self):
        genome = {'chrII': 'A' * 100}
        canon, non_canon = check_canonical_splice_sites([(10, 20)], 'chrI', genome)
        assert canon == 0
        assert non_canon == 0

    def test_yal030w_real_gtag(self):
        """YAL030W intron at chrI:87387-87500 (0-based) has canonical GT-AG."""
        from rectify.data import get_bundled_genome_path
        import gzip as _gz
        genome_path = get_bundled_genome_path('saccharomyces_cerevisiae')
        seq = ''
        with _gz.open(genome_path, 'rt') as fh:
            for line in fh:
                if line.startswith('>') and 'chrI' in line:
                    continue
                if line.startswith('>'):
                    if seq:
                        break
                    continue
                seq += line.strip()
        genome = {'chrI': seq}
        canon, non_canon = check_canonical_splice_sites([(87387, 87500)], 'chrI', genome)
        assert canon == 1, "YAL030W intron should be GT-AG canonical"
        assert non_canon == 0


# ---------------------------------------------------------------------------
# score_alignment
# ---------------------------------------------------------------------------

class TestScoreAlignment:

    def _genome(self, chrom='chrI', size=1_000_000, content='A'):
        return {chrom: content * size}

    def test_clean_alignment_scores_zero(self):
        al = _make_alignment_info(effective_five_prime_clip=0,
                                   three_prime_atract_depth=0,
                                   effective_three_prime_clip=0)
        score = score_alignment(al, self._genome())
        assert score == 0.0

    def test_five_prime_clip_penalty_is_minus_2_per_base(self):
        al = _make_alignment_info(effective_five_prime_clip=5,
                                   three_prime_atract_depth=0,
                                   effective_three_prime_clip=0)
        score = score_alignment(al, self._genome())
        assert score == -10.0  # 5 bases × -2

    def test_atract_depth_penalty_capped_at_10(self):
        al = _make_alignment_info(effective_five_prime_clip=0,
                                   three_prime_atract_depth=15,
                                   effective_three_prime_clip=0)
        score = score_alignment(al, self._genome())
        assert score == -10.0  # capped at 10

    def test_atract_depth_penalty_uncapped_below_10(self):
        al = _make_alignment_info(effective_five_prime_clip=0,
                                   three_prime_atract_depth=7,
                                   effective_three_prime_clip=0)
        score = score_alignment(al, self._genome())
        assert score == -7.0

    def test_three_prime_clip_penalty_minus_2_per_base_capped_at_10(self):
        al = _make_alignment_info(effective_five_prime_clip=0,
                                   three_prime_atract_depth=0,
                                   effective_three_prime_clip=3)
        score = score_alignment(al, self._genome())
        assert score == -6.0  # 3 × -2

    def test_three_prime_clip_penalty_capped(self):
        al = _make_alignment_info(effective_five_prime_clip=0,
                                   three_prime_atract_depth=0,
                                   effective_three_prime_clip=20)
        score = score_alignment(al, self._genome())
        assert score == -10.0  # capped

    def test_combined_penalties(self):
        al = _make_alignment_info(effective_five_prime_clip=3,
                                   three_prime_atract_depth=4,
                                   effective_three_prime_clip=2)
        score = score_alignment(al, self._genome())
        # 5' penalty: 3 × -2 = -6; atract: -4; 3' clip: 2 × -2 = -4 → total -14
        assert score == -14.0

    def test_score_written_back_to_alignment(self):
        al = _make_alignment_info(effective_five_prime_clip=2,
                                   three_prime_atract_depth=0,
                                   effective_three_prime_clip=0)
        score_alignment(al, self._genome())
        assert al.junction_score == -4.0


# ---------------------------------------------------------------------------
# select_best_alignment
# ---------------------------------------------------------------------------

class TestSelectBestAlignment:

    def _genome(self, chrom='chrI', size=1_000_000):
        return {chrom: 'A' * size}

    def test_single_aligner_always_wins(self):
        al = _make_alignment_info(aligner='minimap2', effective_five_prime_clip=5)
        result = select_best_alignment({'minimap2': al}, self._genome())
        assert result.best_aligner == 'minimap2'
        assert result.confidence == 'high'

    def test_clean_alignment_beats_clipped(self):
        clean = _make_alignment_info(aligner='minimap2',
                                      effective_five_prime_clip=0,
                                      three_prime_atract_depth=0,
                                      effective_three_prime_clip=0)
        clipped = _make_alignment_info(aligner='mapPacBio',
                                        effective_five_prime_clip=10,
                                        three_prime_atract_depth=0,
                                        effective_three_prime_clip=0)
        result = select_best_alignment({'minimap2': clean, 'mapPacBio': clipped},
                                        self._genome())
        assert result.best_aligner == 'minimap2'

    def test_was_5prime_rescued_flag(self):
        clean = _make_alignment_info(aligner='minimap2',
                                      effective_five_prime_clip=0)
        clipped = _make_alignment_info(aligner='mapPacBio',
                                        effective_five_prime_clip=8)
        result = select_best_alignment({'minimap2': clean, 'mapPacBio': clipped},
                                        self._genome())
        assert result.was_5prime_rescued is True

    def test_no_5prime_rescue_when_scores_equal_and_no_clip_difference(self):
        al1 = _make_alignment_info(aligner='a1', effective_five_prime_clip=0)
        al2 = _make_alignment_info(aligner='a2', effective_five_prime_clip=0)
        result = select_best_alignment({'a1': al1, 'a2': al2}, self._genome())
        assert result.was_5prime_rescued is False

    def test_tiebreaker_uses_annotated_junction_count(self):
        # Both aligners have equal clip/atract score but different annotated junction counts.
        # al_annotated has 1 annotated junction; al_novel has none.
        intron_start, intron_end = 500, 700
        annotated = {('chrI', intron_start, intron_end, '+')}

        al_annotated = _make_alignment_info(
            aligner='minimap2',
            junctions=[(intron_start, intron_end)],
            effective_five_prime_clip=0,
            canonical_count=1,
        )
        al_novel = _make_alignment_info(
            aligner='mapPacBio',
            junctions=[(600, 800)],   # different, unannotated junction
            effective_five_prime_clip=0,
            canonical_count=0,
        )
        result = select_best_alignment(
            {'minimap2': al_annotated, 'mapPacBio': al_novel},
            self._genome(),
            annotated_junctions=annotated,
        )
        assert result.best_aligner == 'minimap2'

    def test_tiebreaker_uses_3prime_agreement(self):
        # Three aligners: two agree on corrected 3', one disagrees.
        # The two that agree have equal junction scores.
        pos_majority = 1099
        al1 = _make_alignment_info(aligner='a1', corrected_3prime=pos_majority, effective_five_prime_clip=0)
        al2 = _make_alignment_info(aligner='a2', corrected_3prime=pos_majority, effective_five_prime_clip=0)
        al3 = _make_alignment_info(aligner='a3', corrected_3prime=1050, effective_five_prime_clip=0)
        result = select_best_alignment(
            {'a1': al1, 'a2': al2, 'a3': al3},
            self._genome(),
        )
        # Winner must be a1 or a2 (both agree with majority 3' end)
        assert result.best_aligner in ('a1', 'a2')

    def test_confidence_high_when_all_agree(self):
        junctions = [(500, 700)]
        al1 = _make_alignment_info(aligner='a1', junctions=junctions, effective_five_prime_clip=0)
        al2 = _make_alignment_info(aligner='a2', junctions=junctions, effective_five_prime_clip=0)
        result = select_best_alignment({'a1': al1, 'a2': al2}, self._genome())
        assert result.confidence == 'high'

    def test_confidence_medium_when_two_of_three_agree(self):
        junctions_shared = [(500, 700)]
        al1 = _make_alignment_info(aligner='a1', junctions=junctions_shared, effective_five_prime_clip=0)
        al2 = _make_alignment_info(aligner='a2', junctions=junctions_shared, effective_five_prime_clip=0)
        al3 = _make_alignment_info(aligner='a3', junctions=[(600, 800)], effective_five_prime_clip=0)
        result = select_best_alignment({'a1': al1, 'a2': al2, 'a3': al3}, self._genome())
        # The winner is one of a1/a2; two aligners agree
        assert result.confidence in ('high', 'medium')
        assert result.n_aligners_agree >= 2

    def test_empty_alignments_returns_none_best(self):
        result = select_best_alignment({}, self._genome())
        assert result.best_aligner == 'none'
        assert result.best_alignment is None

    def test_aligners_compared_list_matches_input(self):
        al1 = _make_alignment_info(aligner='minimap2')
        al2 = _make_alignment_info(aligner='mapPacBio')
        result = select_best_alignment({'minimap2': al1, 'mapPacBio': al2}, self._genome())
        assert set(result.aligners_compared) == {'minimap2', 'mapPacBio'}

    def test_is_best_flag_set_on_winner(self):
        al = _make_alignment_info(aligner='minimap2')
        result = select_best_alignment({'minimap2': al}, self._genome())
        assert result.best_alignment.is_best is True


# ---------------------------------------------------------------------------
# Real-data tests — wt_by4742_rep1 BAM (skipped when not available)
# ---------------------------------------------------------------------------

_REAL_BAM = Path(
    # Set RECTIFY_TEST_BAM env var or place file at this path for real-data tests
    os.environ.get('RECTIFY_TEST_BAM', '/path/to/wt_by4742_rep1.sorted.bam')
)
_real_bam = pytest.mark.skipif(not _REAL_BAM.exists(), reason='Real BAM not available')


@pytest.fixture(scope='module')
def yal030w_reads():
    """Return a list of real primary reads spanning the YAL030W intron region."""
    if not _REAL_BAM.exists():
        pytest.skip('Real BAM not available')
    reads = []
    with pysam.AlignmentFile(str(_REAL_BAM), 'rb') as bam:
        for read in bam.fetch('chrI', 86500, 88500):
            if not (read.is_secondary or read.is_supplementary or read.is_unmapped):
                reads.append(read)
    return reads


class TestWithRealBam:
    """Tests using real wt_by4742_rep1 reads at the YAL030W locus.

    YAL030W has a single canonical intron at chrI:87387-87500 (0-based, GT-AG).
    """

    @_real_bam
    def test_region_contains_reads(self, yal030w_reads):
        assert len(yal030w_reads) > 0, "Expected reads near YAL030W"

    @_real_bam
    def test_extract_junctions_finds_yal030w_intron(self, yal030w_reads):
        """At least some reads should show the expected intron junction."""
        found_intron = False
        for read in yal030w_reads:
            for junc_start, junc_end in extract_junctions_from_cigar(read):
                # Allow 5 bp tolerance for soft-clip / trimming variation
                if abs(junc_start - 87387) <= 5 and abs(junc_end - 87500) <= 5:
                    found_intron = True
                    break
            if found_intron:
                break
        assert found_intron, (
            "Expected at least one read with junction near YAL030W intron (87387-87500)"
        )

    @_real_bam
    def test_junctions_are_non_overlapping(self, yal030w_reads):
        """Junctions within any single read must be non-overlapping and ordered."""
        for read in yal030w_reads[:100]:
            junctions = extract_junctions_from_cigar(read)
            for i in range(1, len(junctions)):
                assert junctions[i][0] >= junctions[i - 1][1], (
                    f"Overlapping junctions in read {read.query_name}: "
                    f"{junctions[i-1]} vs {junctions[i]}"
                )

    @_real_bam
    def test_select_best_alignment_with_real_junction(self, yal030w_reads):
        """select_best_alignment must prefer the aligner that captures the
        YAL030W junction over a hypothetical aligner that soft-clips it.

        We construct two AlignmentInfo objects from the first spliced real read:
        - 'minimap2': uses the real junction (as observed)
        - 'mapPacBio': simulated as having a 20 bp 5' soft-clip instead
        """
        # Find first read with a junction
        spliced = None
        for read in yal030w_reads:
            if any(op == 3 for op, _ in (read.cigartuples or [])):
                spliced = read
                break
        if spliced is None:
            pytest.skip("No spliced read found in test region")

        junctions = extract_junctions_from_cigar(spliced)
        chrom = spliced.reference_name
        strand = '-' if spliced.is_reverse else '+'

        al_spliced = _make_alignment_info(
            read_id=spliced.query_name,
            aligner='minimap2',
            chrom=chrom,
            strand=strand,
            start=spliced.reference_start,
            end=spliced.reference_end,
            junctions=junctions,
            effective_five_prime_clip=0,
            three_prime_atract_depth=0,
        )
        al_clipped = _make_alignment_info(
            read_id=spliced.query_name,
            aligner='mapPacBio',
            chrom=chrom,
            strand=strand,
            start=spliced.reference_start,
            end=spliced.reference_end,
            junctions=[],
            effective_five_prime_clip=20,  # simulated: soft-clips the exon instead of splicing
            three_prime_atract_depth=0,
        )

        result = select_best_alignment(
            {'minimap2': al_spliced, 'mapPacBio': al_clipped},
            {chrom: 'A' * 3_000_000},
        )
        assert result.best_aligner == 'minimap2', (
            "Aligner with the splice junction should beat the one that soft-clips"
        )
        assert result.was_5prime_rescued is True
