#!/usr/bin/env python3
"""
Unit tests for terminal_exon_refiner.py.

Covers:
- SpliceSiteIndex add/lookup methods
- load_splice_sites_from_gff (plus and minus strand, 1-based→0-based conversion)
- merge_splice_indices (annotation-takes-precedence rule)
- detect_junction_truncated_reads (plus and minus strand truncation detection)
- detect_partial_junction_crossings (soft-clip boundary detection, both strands)
- get_soft_clip_info (left/right clip extraction)
- simple_align (best-position search, mismatch tolerance)
- detect_mismatch_clusters (cluster detection logic)
"""

import io
import gzip
import os
import struct
import tempfile
import textwrap
from pathlib import Path
from typing import List, Tuple
from unittest.mock import Mock

import pytest
import pysam

from rectify.core.terminal_exon_refiner import (
    SpliceSite,
    SpliceSiteIndex,
    load_splice_sites_from_gff,
    merge_splice_indices,
    detect_junction_truncated_reads,
    detect_partial_junction_crossings,
    get_soft_clip_info,
    simple_align,
    detect_mismatch_clusters,
)


# ---------------------------------------------------------------------------
# Helpers for building synthetic BAM files
# ---------------------------------------------------------------------------

def _make_bam_read(
    read_id: str,
    chrom: str,
    start: int,
    cigar: List[Tuple[int, int]],
    strand: str,
    seq: str,
) -> pysam.AlignedSegment:
    """Build a synthetic pysam.AlignedSegment without writing to disk."""
    hdr = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6'},
        'SQ': [{'SN': chrom, 'LN': 1_000_000}],
    })
    read = pysam.AlignedSegment(hdr)
    read.query_name = read_id
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


def _write_bam(tmp_path: Path, reads: List[pysam.AlignedSegment], chrom: str = 'chrI') -> str:
    """Write reads to a temp BAM file and return its path as a string."""
    bam_path = str(tmp_path / 'test.bam')
    hdr = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6'},
        'SQ': [{'SN': chrom, 'LN': 1_000_000}],
    })
    with pysam.AlignmentFile(bam_path, 'wb', header=hdr) as bam:
        for r in reads:
            bam.write(r)
    pysam.sort('-o', bam_path + '.sorted.bam', bam_path)
    pysam.index(bam_path + '.sorted.bam')
    return bam_path + '.sorted.bam'


def _write_gff(tmp_path: Path, lines: str) -> str:
    """Write GFF lines to a temp .gff file, return path."""
    p = tmp_path / 'ann.gff'
    p.write_text(lines)
    return str(p)


# ---------------------------------------------------------------------------
# SpliceSiteIndex
# ---------------------------------------------------------------------------

class TestSpliceSiteIndex:

    def test_add_and_retrieve_5ss(self):
        idx = SpliceSiteIndex()
        site = SpliceSite(chrom='chrI', position=100, strand='+', site_type='5ss')
        idx.add_site(site)
        assert 'chrI' in idx.five_ss
        assert 100 in idx.five_ss['chrI']

    def test_add_and_retrieve_3ss(self):
        idx = SpliceSiteIndex()
        site = SpliceSite(chrom='chrI', position=200, strand='-', site_type='3ss')
        idx.add_site(site)
        assert 200 in idx.three_ss['chrI']

    def test_get_nearby_5ss_within_distance(self):
        idx = SpliceSiteIndex()
        idx.add_site(SpliceSite('chrI', 100, '+', '5ss'))
        idx.add_site(SpliceSite('chrI', 160, '+', '5ss'))
        nearby = idx.get_nearby_5ss('chrI', 105, max_distance=10)
        assert len(nearby) == 1
        assert nearby[0].position == 100

    def test_get_nearby_5ss_sorted_by_distance(self):
        idx = SpliceSiteIndex()
        idx.add_site(SpliceSite('chrI', 100, '+', '5ss'))
        idx.add_site(SpliceSite('chrI', 108, '+', '5ss'))
        nearby = idx.get_nearby_5ss('chrI', 105, max_distance=20)
        assert nearby[0].position == 108  # 3 bp away, closer than 100 (5 bp)

    def test_get_nearby_empty_chrom(self):
        idx = SpliceSiteIndex()
        result = idx.get_nearby_5ss('chrVII', 500, max_distance=50)
        assert result == []

    def test_5ss_and_3ss_stored_separately(self):
        idx = SpliceSiteIndex()
        idx.add_site(SpliceSite('chrI', 100, '+', '5ss'))
        idx.add_site(SpliceSite('chrI', 100, '+', '3ss'))
        assert 100 in idx.five_ss.get('chrI', {})
        assert 100 in idx.three_ss.get('chrI', {})


# ---------------------------------------------------------------------------
# load_splice_sites_from_gff
# ---------------------------------------------------------------------------

class TestLoadSpliceSitesFromGff:

    def _intron_gff(self, chrom, start_1based, end_1based, strand, name='intron1'):
        """Build a minimal 1-based GFF intron record."""
        return (
            f'{chrom}\t.\tintron\t{start_1based}\t{end_1based}\t.\t{strand}\t.\t'
            f'Name={name}\n'
        )

    def test_plus_strand_5ss_position(self, tmp_path):
        # Intron at 1-based [101, 200] on + strand
        # 0-based: [100, 200)  → 5'SS at 100, 3'SS at 198
        gff = _write_gff(tmp_path, self._intron_gff('chrI', 101, 200, '+'))
        idx = load_splice_sites_from_gff(gff)
        assert 100 in idx.five_ss.get('chrI', {}), "5'SS should be at intron start (0-based)"

    def test_plus_strand_3ss_position(self, tmp_path):
        gff = _write_gff(tmp_path, self._intron_gff('chrI', 101, 200, '+'))
        idx = load_splice_sites_from_gff(gff)
        assert 198 in idx.three_ss.get('chrI', {}), "3'SS should be 2 bp before intron end (0-based)"

    def test_minus_strand_5ss_at_intron_end(self, tmp_path):
        # Minus strand: 5'SS is at high coordinate (intron end on genomic coords)
        # Intron 1-based [101, 200] → 0-based [100, 200)
        # 5'SS = end_0 - 2 = 198;  3'SS = start_0 = 100
        gff = _write_gff(tmp_path, self._intron_gff('chrI', 101, 200, '-'))
        idx = load_splice_sites_from_gff(gff)
        assert 198 in idx.five_ss.get('chrI', {}), "Minus 5'SS should be at intron genomic end"

    def test_minus_strand_3ss_at_intron_start(self, tmp_path):
        gff = _write_gff(tmp_path, self._intron_gff('chrI', 101, 200, '-'))
        idx = load_splice_sites_from_gff(gff)
        assert 100 in idx.three_ss.get('chrI', {}), "Minus 3'SS should be at intron genomic start"

    def test_comment_lines_skipped(self, tmp_path):
        content = '# comment\n' + self._intron_gff('chrI', 101, 200, '+')
        gff = _write_gff(tmp_path, content)
        idx = load_splice_sites_from_gff(gff)
        assert idx.five_ss  # Should still load the intron

    def test_non_intron_features_skipped(self, tmp_path):
        content = 'chrI\t.\texon\t101\t200\t.\t+\t.\tName=exon1\n'
        gff = _write_gff(tmp_path, content)
        idx = load_splice_sites_from_gff(gff)
        assert not idx.five_ss
        assert not idx.three_ss

    def test_gz_gff(self, tmp_path):
        line = self._intron_gff('chrI', 101, 200, '+')
        p = tmp_path / 'ann.gff.gz'
        with gzip.open(str(p), 'wt') as f:
            f.write(line)
        idx = load_splice_sites_from_gff(str(p))
        assert idx.five_ss


# ---------------------------------------------------------------------------
# merge_splice_indices
# ---------------------------------------------------------------------------

class TestMergeSpliceIndices:

    def test_merge_combines_sites(self):
        a = SpliceSiteIndex()
        a.add_site(SpliceSite('chrI', 100, '+', '5ss', intron_id='A'))
        b = SpliceSiteIndex()
        b.add_site(SpliceSite('chrI', 200, '+', '5ss', intron_id='B'))
        merged = merge_splice_indices(a, b)
        assert 100 in merged.five_ss['chrI']
        assert 200 in merged.five_ss['chrI']

    def test_first_index_takes_precedence(self):
        """Annotation sites should not be overwritten by novel sites."""
        a = SpliceSiteIndex()
        a.add_site(SpliceSite('chrI', 100, '+', '5ss', intron_id='annotated'))
        b = SpliceSiteIndex()
        b.add_site(SpliceSite('chrI', 100, '+', '5ss', intron_id='novel'))
        merged = merge_splice_indices(a, b)
        assert merged.five_ss['chrI'][100].intron_id == 'annotated'

    def test_merge_empty_indices(self):
        a = SpliceSiteIndex()
        b = SpliceSiteIndex()
        merged = merge_splice_indices(a, b)
        assert not merged.five_ss
        assert not merged.three_ss

    def test_merge_three_indices(self):
        a, b, c = SpliceSiteIndex(), SpliceSiteIndex(), SpliceSiteIndex()
        a.add_site(SpliceSite('chrI', 100, '+', '5ss', intron_id='a'))
        b.add_site(SpliceSite('chrI', 200, '+', '5ss', intron_id='b'))
        c.add_site(SpliceSite('chrI', 300, '+', '5ss', intron_id='c'))
        merged = merge_splice_indices(a, b, c)
        assert {100, 200, 300} == set(merged.five_ss['chrI'].keys())


# ---------------------------------------------------------------------------
# get_soft_clip_info
# ---------------------------------------------------------------------------

class TestGetSoftClipInfo:

    def _make_read(self, cigar, seq, strand='+'):
        return _make_bam_read('r1', 'chrI', 1000, cigar, strand, seq)

    def test_no_clips(self):
        read = self._make_read([(0, 10)], 'A' * 10)
        info = get_soft_clip_info(read)
        assert info['left']['length'] == 0
        assert info['right']['length'] == 0

    def test_left_clip(self):
        # 5S10M
        read = self._make_read([(4, 5), (0, 10)], 'C' * 5 + 'A' * 10)
        info = get_soft_clip_info(read)
        assert info['left']['length'] == 5
        assert info['left']['sequence'] == 'CCCCC'
        assert info['right']['length'] == 0

    def test_right_clip(self):
        # 10M5S
        read = self._make_read([(0, 10), (4, 5)], 'A' * 10 + 'T' * 5)
        info = get_soft_clip_info(read)
        assert info['right']['length'] == 5
        assert info['right']['sequence'] == 'TTTTT'
        assert info['left']['length'] == 0

    def test_both_clips(self):
        # 3S8M4S
        read = self._make_read([(4, 3), (0, 8), (4, 4)], 'G' * 3 + 'A' * 8 + 'C' * 4)
        info = get_soft_clip_info(read)
        assert info['left']['length'] == 3
        assert info['right']['length'] == 4

    def test_no_cigar(self):
        read = self._make_read([], '')
        read.cigartuples = None
        info = get_soft_clip_info(read)
        assert info['left']['length'] == 0


# ---------------------------------------------------------------------------
# simple_align
# ---------------------------------------------------------------------------

class TestSimpleAlign:

    def test_perfect_match(self):
        pos, score, mismatches = simple_align('ACGT', 'XXXXXACGTYYYYY')
        assert pos == 5
        assert mismatches == 0

    def test_match_at_start(self):
        pos, score, mismatches = simple_align('ACG', 'ACGTTTT')
        assert pos == 0

    def test_no_match_returns_minus_one(self):
        pos, score, mismatches = simple_align('AAAA', 'TTTTTTTT', max_mismatches=0)
        assert pos == -1

    def test_single_mismatch_within_tolerance(self):
        pos, score, mismatches = simple_align('ACGT', 'XXXTCGT', max_mismatches=1)
        # ACGT vs TCGT at pos 3: 1 mismatch
        assert pos != -1
        assert mismatches <= 1

    def test_case_insensitive(self):
        pos, score, mismatches = simple_align('acgt', 'XXACGTXX')
        assert pos != -1 and mismatches == 0

    def test_empty_query(self):
        pos, score, mismatches = simple_align('', 'ACGT')
        assert pos == -1

    def test_query_longer_than_target(self):
        pos, score, mismatches = simple_align('ACGTACGT', 'ACG')
        assert pos == -1


# ---------------------------------------------------------------------------
# detect_mismatch_clusters
# ---------------------------------------------------------------------------

class TestDetectMismatchClusters:

    def test_no_cigar_returns_empty(self):
        read = _make_bam_read('r1', 'chrI', 1000, [], '+', '')
        read.cigartuples = None
        clusters = detect_mismatch_clusters(read)
        assert clusters == []

    def test_pure_match_no_cluster(self):
        read = _make_bam_read('r1', 'chrI', 1000, [(0, 100)], '+', 'A' * 100)
        clusters = detect_mismatch_clusters(read)
        assert clusters == []

    def test_clustered_deletions_detected(self):
        # 20M 3D 20M 3D 20M in a tight window — should trigger a cluster
        cigar = [(0, 5), (2, 2), (0, 5), (2, 2), (0, 5), (2, 2), (0, 5)]
        seq_len = 5 + 5 + 5 + 5  # M ops only
        read = _make_bam_read('r1', 'chrI', 1000, cigar, '+', 'A' * seq_len)
        clusters = detect_mismatch_clusters(read, window_size=20, min_events=3, min_density=0.1)
        assert len(clusters) >= 1

    def test_sparse_events_no_cluster(self):
        # Two deletions far apart
        cigar = [(0, 50), (2, 2), (0, 200), (2, 2), (0, 50)]
        seq = 'A' * 100
        read = _make_bam_read('r1', 'chrI', 1000, cigar, '+', seq)
        clusters = detect_mismatch_clusters(read, window_size=15, min_events=3)
        assert clusters == []

    def test_cluster_has_expected_fields(self):
        cigar = [(0, 3), (2, 1), (0, 3), (2, 1), (0, 3), (2, 1), (0, 3)]
        seq = 'A' * 12
        read = _make_bam_read('r1', 'chrI', 500, cigar, '-', seq)
        clusters = detect_mismatch_clusters(read, window_size=15, min_events=3, min_density=0.1)
        if clusters:
            c = clusters[0]
            assert c.chrom == 'chrI'
            assert c.strand == '-'
            assert c.start >= 500
            assert c.n_deletions >= 1


# ---------------------------------------------------------------------------
# detect_junction_truncated_reads
# ---------------------------------------------------------------------------

class TestDetectJunctionTruncatedReads:

    def _make_splice_index(self, chrom, intron_start, intron_end, strand):
        """Build an index with one intron on the given strand."""
        idx = SpliceSiteIndex()
        if strand == '+':
            idx.add_site(SpliceSite(chrom, intron_start, strand, '5ss', intron_id='i1'))
            idx.add_site(SpliceSite(chrom, intron_end - 2, strand, '3ss', intron_id='i1'))
        else:
            idx.add_site(SpliceSite(chrom, intron_end - 2, strand, '5ss', intron_id='i1'))
            idx.add_site(SpliceSite(chrom, intron_start, strand, '3ss', intron_id='i1'))
        return idx

    def test_plus_strand_truncated_at_3ss(self, tmp_path):
        # Intron: [500, 700).  3'SS (acceptor) is at 698.
        # A read that ends at 698 with no junction should be detected.
        intron_end = 700
        three_ss_pos = intron_end - 2  # 698

        idx = self._make_splice_index('chrI', 500, intron_end, '+')

        # Plus-strand read that ends AT the 3'SS (reference_end = three_ss_pos + 1 = 699)
        # reference_end is exclusive, so to have read_end == 698+1 we set start=400, 299M
        read = _make_bam_read('r1', 'chrI', 400, [(0, 299)], '+', 'A' * 299)
        # read.reference_end == 400 + 299 = 699, so read_end = 699
        # abs(699 - 698) = 1 ≤ tolerance=2 → should detect
        bam_path = _write_bam(tmp_path, [read])
        results = detect_junction_truncated_reads(bam_path, idx, tolerance=2)
        assert results['stats']['truncated_at_3ss'] >= 1

    def test_minus_strand_truncated_at_3ss(self, tmp_path):
        # Minus strand 3'SS is at intron_start (low coordinate).
        # Intron [500, 700): 3'SS at 500.
        # A minus-strand read that starts at 500 (read_start = 500) with no junction.
        intron_start = 500
        idx = self._make_splice_index('chrI', intron_start, 700, '-')

        # read_start == intron_start, length 100
        read = _make_bam_read('r1', 'chrI', intron_start, [(0, 100)], '-', 'T' * 100)
        bam_path = _write_bam(tmp_path, [read])
        results = detect_junction_truncated_reads(bam_path, idx, tolerance=2)
        assert results['stats']['truncated_at_3ss'] >= 1

    def test_read_with_junction_not_truncated(self, tmp_path):
        # A read that crosses the junction (has N op) should NOT appear in results.
        idx = self._make_splice_index('chrI', 500, 700, '+')
        # 100M 200N 100M → properly spliced read
        read = _make_bam_read('r1', 'chrI', 400, [(0, 100), (3, 200), (0, 100)], '+', 'A' * 200)
        bam_path = _write_bam(tmp_path, [read])
        results = detect_junction_truncated_reads(bam_path, idx, tolerance=10)
        assert results['stats']['truncated_at_3ss'] == 0

    def test_read_far_from_ss_not_detected(self, tmp_path):
        idx = self._make_splice_index('chrI', 500, 700, '+')
        # Read ends far from any splice site
        read = _make_bam_read('r1', 'chrI', 100, [(0, 50)], '+', 'A' * 50)
        bam_path = _write_bam(tmp_path, [read])
        results = detect_junction_truncated_reads(bam_path, idx, tolerance=5)
        assert results['stats']['truncated_at_3ss'] == 0

    def test_result_has_expected_keys(self, tmp_path):
        idx = self._make_splice_index('chrI', 500, 700, '+')
        bam_path = _write_bam(tmp_path, [])
        results = detect_junction_truncated_reads(bam_path, idx)
        assert 'truncated_at_3ss' in results
        assert 'truncated_at_5ss' in results
        assert 'stats' in results


# ---------------------------------------------------------------------------
# detect_partial_junction_crossings — boundary and strand checks
# ---------------------------------------------------------------------------

class TestDetectPartialJunctionCrossings:

    def _genome(self, chrom='chrI', size=1_000_000):
        return {chrom: 'A' * size}

    def _make_index_plus_with_3ss_at(self, pos, chrom='chrI'):
        idx = SpliceSiteIndex()
        idx.add_site(SpliceSite(chrom, pos, '+', '3ss', intron_id='i1'))
        idx.add_site(SpliceSite(chrom, pos - 100, '+', '5ss', intron_id='i1'))
        return idx

    def test_plus_strand_read_at_boundary_detected(self, tmp_path):
        # 3'SS at position 500.  Read starts at 500 with a left soft-clip.
        # 5S50M  → reference_start = 500, left clip length = 5 ≥ 1
        three_ss_pos = 500
        idx = self._make_index_plus_with_3ss_at(three_ss_pos)
        read = _make_bam_read('r1', 'chrI', three_ss_pos, [(4, 5), (0, 50)], '+', 'G' * 5 + 'A' * 50)
        bam_path = _write_bam(tmp_path, [read])
        genome = self._genome()
        results = detect_partial_junction_crossings(
            bam_path, genome, idx, boundary_tolerance=3, min_clip_length=1
        )
        assert results['stats']['reads_at_splice_boundary'] >= 1

    def test_read_without_softclip_not_counted_as_crossing(self, tmp_path):
        three_ss_pos = 500
        idx = self._make_index_plus_with_3ss_at(three_ss_pos)
        # No soft-clip: pure 50M
        read = _make_bam_read('r1', 'chrI', three_ss_pos, [(0, 50)], '+', 'A' * 50)
        bam_path = _write_bam(tmp_path, [read])
        genome = self._genome()
        results = detect_partial_junction_crossings(
            bam_path, genome, idx, boundary_tolerance=3, min_clip_length=1
        )
        # reads_at_splice_boundary may be 1 but reads_with_softclip_at_boundary should be 0
        assert results['stats']['reads_with_softclip_at_boundary'] == 0

    def test_read_with_junction_skipped(self, tmp_path):
        three_ss_pos = 500
        idx = self._make_index_plus_with_3ss_at(three_ss_pos)
        # Read has an N op → already spliced, should be skipped
        read = _make_bam_read('r1', 'chrI', 400, [(0, 50), (3, 100), (0, 50)], '+', 'A' * 100)
        bam_path = _write_bam(tmp_path, [read])
        genome = self._genome()
        results = detect_partial_junction_crossings(
            bam_path, genome, idx, boundary_tolerance=10, min_clip_length=1
        )
        assert results['stats']['total_reads'] == 1
        # Properly spliced reads skip the boundary check entirely
        reads_checked = results['stats'].get('reads_at_splice_boundary', 0)
        assert reads_checked == 0

    def test_result_contains_expected_stat_keys(self, tmp_path):
        idx = SpliceSiteIndex()
        bam_path = _write_bam(tmp_path, [])
        genome = self._genome()
        results = detect_partial_junction_crossings(bam_path, genome, idx)
        expected = {'total_reads', 'reads_at_splice_boundary', 'reads_with_softclip_at_boundary'}
        assert expected.issubset(set(results['stats'].keys()))


# ---------------------------------------------------------------------------
# Real-data tests — bundled GFF (always available) + real BAM (skip if absent)
# ---------------------------------------------------------------------------

from pathlib import Path as _Path

# Bundled GFF is part of the package and always present
def _bundled_gff() -> str:
    from rectify.data import get_bundled_annotation_path
    return get_bundled_annotation_path('saccharomyces_cerevisiae')


# Real wt_by4742_rep1 BAM — present on Oak but not in CI
_REAL_BAM = _Path(
    # Set RECTIFY_TEST_BAM env var or place file at this path for real-data tests
    os.environ.get('RECTIFY_TEST_BAM', '/path/to/wt_by4742_rep1.sorted.bam')
)
_real_bam = pytest.mark.skipif(not _REAL_BAM.exists(), reason='Real BAM not available')


class TestLoadSpliceSitesFromBundledGff:
    """Verify load_splice_sites_from_gff produces correct positions from the
    bundled S. cerevisiae annotation (R64-5-1).  No real BAM required."""

    @pytest.fixture(scope='class')
    def splice_index(self):
        return load_splice_sites_from_gff(_bundled_gff())

    def test_loads_expected_intron_count(self, splice_index):
        # R64-5-1 has 402 annotated introns but 17 share a splice-site position
        # with another intron (alternative 5'SS / 3'SS).  The dict stores one
        # site per position, giving 385 unique entries.
        n_5ss = sum(len(v) for v in splice_index.five_ss.values())
        n_3ss = sum(len(v) for v in splice_index.three_ss.values())
        assert n_5ss == 385
        assert n_3ss == 385

    def test_chrom_names_match_bam_convention(self, splice_index):
        # Bundled GFF uses 'chrI', 'chrII', etc. — must match aligned BAM
        assert 'chrI' in splice_index.five_ss, "Expected 'chrI' not 'I' in splice index"

    def test_yAL030W_plus_strand_5ss_position(self, splice_index):
        # YAL030W_intron: 1-based [87388, 87500] on + strand
        # 0-based 5'SS = 87387 (intron start, 0-based)
        five_ss_positions = splice_index.five_ss.get('chrI', {})
        assert 87387 in five_ss_positions, (
            f"Expected 5'SS at 87387 for YAL030W; found positions: "
            f"{sorted(five_ss_positions.keys())[:10]}"
        )
        assert five_ss_positions[87387].strand == '+'

    def test_yAL030W_plus_strand_3ss_position(self, splice_index):
        # 3'SS = end_0 - 2 = 87500 - 2 = 87498
        three_ss_positions = splice_index.three_ss.get('chrI', {})
        assert 87498 in three_ss_positions, (
            f"Expected 3'SS at 87498 for YAL030W; found: {sorted(three_ss_positions.keys())[:10]}"
        )

    def test_yAL001C_minus_strand_5ss_position(self, splice_index):
        # YAL001C_intron: 1-based [151007, 151096] on - strand
        # Minus 5'SS = end_0 - 2 = 151096 - 2 = 151094
        five_ss_positions = splice_index.five_ss.get('chrI', {})
        assert 151094 in five_ss_positions, (
            f"Expected minus-strand 5'SS at 151094 for YAL001C"
        )
        assert five_ss_positions[151094].strand == '-'

    def test_yAL001C_minus_strand_3ss_position(self, splice_index):
        # Minus 3'SS = start_0 = 151007 - 1 = 151006
        three_ss_positions = splice_index.three_ss.get('chrI', {})
        assert 151006 in three_ss_positions, (
            f"Expected minus-strand 3'SS at 151006 for YAL001C"
        )
        assert three_ss_positions[151006].strand == '-'

    def test_intron_id_field_populated(self, splice_index):
        # intron_id should be set from the Name= attribute
        site = splice_index.five_ss['chrI'][87387]
        assert site.intron_id, "intron_id should not be empty"
        assert 'YAL030W' in site.intron_id


@pytest.fixture(scope='module')
def yal030w_region_bam(tmp_path_factory):
    """Extract chrI:86500-88500 from the real BAM into a temp BAM.

    Skipped automatically when the real BAM is not present.
    """
    if not _REAL_BAM.exists():
        pytest.skip('Real BAM not available')
    tmp_path = tmp_path_factory.mktemp('real_bam')
    out_bam = str(tmp_path / 'yal030w_region.bam')
    with pysam.AlignmentFile(str(_REAL_BAM), 'rb') as src:
        with pysam.AlignmentFile(out_bam, 'wb', header=src.header) as dst:
            for read in src.fetch('chrI', 86500, 88500):
                if not (read.is_secondary or read.is_supplementary or read.is_unmapped):
                    dst.write(read)
    pysam.sort('-o', out_bam + '.sorted.bam', out_bam)
    pysam.index(out_bam + '.sorted.bam')
    return out_bam + '.sorted.bam'


class TestWithRealBam:
    """Tests using real wt_by4742_rep1 reads at the YAL030W locus (chrI:87387-87500)."""

    @_real_bam
    def test_real_reads_have_junctions_at_yal030w(self, yal030w_region_bam):
        """Reads spanning the YAL030W intron should carry N operations."""
        n_spliced = 0
        with pysam.AlignmentFile(yal030w_region_bam, 'rb') as bam:
            for read in bam.fetch('chrI', 87387, 87500):
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                has_n = any(op == 3 for op, _ in (read.cigartuples or []))
                if has_n:
                    n_spliced += 1
        assert n_spliced > 0, (
            "Expected at least one spliced read spanning the YAL030W intron"
        )

    @_real_bam
    def test_spliced_reads_not_detected_as_truncated(self, yal030w_region_bam):
        """Properly spliced reads crossing the YAL030W intron must NOT appear
        in detect_junction_truncated_reads results."""
        splice_index = load_splice_sites_from_gff(_bundled_gff())
        results = detect_junction_truncated_reads(
            yal030w_region_bam, splice_index, tolerance=5
        )
        # Collect read IDs detected as truncated
        truncated_ids = {r['read_id'] for r in results['truncated_at_3ss']}
        truncated_ids |= {r['read_id'] for r in results['truncated_at_5ss']}

        # Any read in truncated_ids must not have a junction (N op) — verify
        with pysam.AlignmentFile(yal030w_region_bam, 'rb') as bam:
            for read in bam.fetch():
                if read.query_name in truncated_ids:
                    has_n = any(op == 3 for op, _ in (read.cigartuples or []))
                    assert not has_n, (
                        f"Read {read.query_name} has N op but was classified as truncated"
                    )

    @_real_bam
    def test_gff_splice_sites_cover_real_intron(self):
        """The bundled GFF splice index must contain the YAL030W intron position."""
        splice_index = load_splice_sites_from_gff(_bundled_gff())
        assert 87387 in splice_index.five_ss.get('chrI', {}), \
            "YAL030W 5'SS (chrI:87387) must be in the splice index"
        assert 87498 in splice_index.three_ss.get('chrI', {}), \
            "YAL030W 3'SS (chrI:87498) must be in the splice index"

    @_real_bam
    def test_get_soft_clip_info_on_real_reads(self, yal030w_region_bam):
        """get_soft_clip_info must return valid dicts for real reads (no crashes)."""
        with pysam.AlignmentFile(yal030w_region_bam, 'rb') as bam:
            reads_checked = 0
            for read in bam.fetch('chrI', 87000, 88000):
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                info = get_soft_clip_info(read)
                assert 'left' in info and 'right' in info
                assert 'length' in info['left'] and 'length' in info['right']
                reads_checked += 1
                if reads_checked >= 50:
                    break
        assert reads_checked > 0, "No reads found in test region"
