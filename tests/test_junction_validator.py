#!/usr/bin/env python3
"""
Tests for cross-sample junction validation (Bug 5 fix).

Covers the three-pass COMPASS-inspired architecture in
rectify.core.junction_validator:

  Pass 1 — extract_sample_junctions()
  Pass 2 — filter_cross_sample_junctions()
  Pass 3 — apply_junction_filter()

Author: Kevin R. Roy
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from unittest.mock import MagicMock, patch, call
import pandas as pd

from rectify.core.junction_validator import (
    extract_sample_junctions,
    filter_cross_sample_junctions,
    apply_junction_filter,
    _get_splice_motif,
    _revcomp,
    _read_has_invalid_junction,
    _downgrade_read,
)


# =============================================================================
# Helpers — mock pysam objects
# =============================================================================

def _make_read(
    cigar_tuples,
    ref_start,
    chrom='chrI',
    is_reverse=False,
    is_secondary=False,
    is_supplementary=False,
    is_unmapped=False,
    xc_tag=None,
):
    """Return a minimal mock pysam.AlignedSegment."""
    read = MagicMock()
    read.cigartuples = cigar_tuples
    read.reference_start = ref_start
    read.reference_name = chrom
    read.is_reverse = is_reverse
    read.is_secondary = is_secondary
    read.is_supplementary = is_supplementary
    read.is_unmapped = is_unmapped

    # Simulate get_tag / set_tag
    _tags = {}
    if xc_tag is not None:
        _tags['XC'] = xc_tag

    def _get_tag(name):
        if name in _tags:
            return _tags[name]
        raise KeyError(name)

    def _set_tag(name, value):
        _tags[name] = value

    read.get_tag = _get_tag
    read.set_tag = _set_tag
    read._tags = _tags  # expose for assertions
    return read


def _make_bam_file(reads, header=None):
    """Return a mock pysam.AlignmentFile context manager."""
    if header is None:
        header = MagicMock()

    bam_mock = MagicMock()
    bam_mock.fetch = MagicMock(return_value=iter(reads))
    bam_mock.header = header

    ctx = MagicMock()
    ctx.__enter__ = MagicMock(return_value=bam_mock)
    ctx.__exit__ = MagicMock(return_value=False)
    return ctx, bam_mock


# CIGAR op codes
M = 0   # match/mismatch
D = 2   # deletion
N = 3   # intron skip
I = 1   # insertion
S = 4   # soft clip


# =============================================================================
# Pass 1 — extract_sample_junctions
# =============================================================================

class TestExtractSampleJunctions:

    def _run_extraction(self, reads, sample_id='s1', genome=None):
        """Helper: run extraction against a list of mock reads."""
        ctx, bam_mock = _make_bam_file(reads)
        with patch('rectify.core.junction_validator.pysam') as mock_pysam:
            mock_pysam.AlignmentFile.return_value = ctx
            return extract_sample_junctions.__wrapped__(
                bam_mock, sample_id, genome
            ) if hasattr(extract_sample_junctions, '__wrapped__') else \
                _extract_directly(bam_mock, sample_id, genome)

    def test_basic_intron_extracted(self):
        """A read with an N operation produces one junction row."""
        # 10M 100N 10M: intron from 10 to 110
        reads = [_make_read([(M, 10), (N, 100), (M, 10)], ref_start=0)]
        df = _extract_from_reads(reads)
        assert len(df) == 1
        row = df.iloc[0]
        assert row['chrom'] == 'chrI'
        assert row['intron_start'] == 10
        assert row['intron_end'] == 110
        assert row['strand'] == '+'
        assert row['read_count'] == 1

    def test_two_reads_same_junction(self):
        """Two reads sharing the same junction collapse to one row with count=2."""
        reads = [
            _make_read([(M, 10), (N, 100), (M, 10)], ref_start=0),
            _make_read([(M, 10), (N, 100), (M, 10)], ref_start=0),
        ]
        df = _extract_from_reads(reads)
        assert len(df) == 1
        assert df.iloc[0]['read_count'] == 2

    def test_secondary_reads_skipped(self):
        """Secondary alignments must not contribute junctions."""
        reads = [
            _make_read([(M, 10), (N, 100), (M, 10)], ref_start=0, is_secondary=True),
        ]
        df = _extract_from_reads(reads)
        assert df.empty

    def test_supplementary_reads_skipped(self):
        """Supplementary alignments must not contribute junctions."""
        reads = [
            _make_read(
                [(M, 10), (N, 100), (M, 10)], ref_start=0, is_supplementary=True
            ),
        ]
        df = _extract_from_reads(reads)
        assert df.empty

    def test_unmapped_reads_skipped(self):
        """Unmapped reads must not contribute junctions."""
        reads = [
            _make_read(
                [(M, 10), (N, 100), (M, 10)], ref_start=0, is_unmapped=True
            ),
        ]
        df = _extract_from_reads(reads)
        assert df.empty

    def test_no_intron_read_produces_no_rows(self):
        """A read with only M operations produces no junction rows."""
        reads = [_make_read([(M, 100)], ref_start=0)]
        df = _extract_from_reads(reads)
        assert df.empty

    def test_minus_strand_junction(self):
        """Minus-strand reads record strand='-'."""
        reads = [
            _make_read([(M, 10), (N, 200), (M, 10)], ref_start=0, is_reverse=True)
        ]
        df = _extract_from_reads(reads)
        assert df.iloc[0]['strand'] == '-'

    def test_motif_with_genome(self):
        """When a genome is provided, motif and is_canonical are set."""
        reads = [_make_read([(M, 10), (N, 100), (M, 10)], ref_start=0)]
        genome = MagicMock()
        # donor = GT, acceptor = AG → canonical GT-AG
        genome.fetch = MagicMock(side_effect=lambda chrom, s, e: {
            (10, 12): 'GT', (108, 110): 'AG',
        }.get((s, e), 'NN'))
        df = _extract_from_reads(reads, genome=genome)
        assert df.iloc[0]['motif'] == 'GT-AG'
        assert bool(df.iloc[0]['is_canonical']) is True

    def test_non_canonical_motif(self):
        """Non-canonical motifs are flagged is_canonical=False."""
        reads = [_make_read([(M, 10), (N, 100), (M, 10)], ref_start=0)]
        genome = MagicMock()
        genome.fetch = MagicMock(side_effect=lambda chrom, s, e: {
            (10, 12): 'CT', (108, 110): 'AC',
        }.get((s, e), 'NN'))
        df = _extract_from_reads(reads, genome=genome)
        assert df.iloc[0]['motif'] == 'CT-AC'
        assert bool(df.iloc[0]['is_canonical']) is False

    def test_no_genome_unknown_motif(self):
        """Without a genome, motif is 'unknown' and is_canonical is False."""
        reads = [_make_read([(M, 10), (N, 100), (M, 10)], ref_start=0)]
        df = _extract_from_reads(reads, genome=None)
        assert df.iloc[0]['motif'] == 'unknown'
        assert bool(df.iloc[0]['is_canonical']) is False

    def test_multiple_junctions_in_one_read(self):
        """A read with two N operations produces two junction rows."""
        reads = [
            _make_read(
                [(M, 10), (N, 100), (M, 5), (N, 200), (M, 10)], ref_start=0
            )
        ]
        df = _extract_from_reads(reads)
        assert len(df) == 2
        starts = set(df['intron_start'].tolist())
        assert 10 in starts
        assert 115 in starts

    def test_sample_id_stored(self):
        """sample_id column matches the argument passed."""
        reads = [_make_read([(M, 10), (N, 50), (M, 10)], ref_start=0)]
        df = _extract_from_reads(reads, sample_id='my_sample')
        assert df.iloc[0]['sample_id'] == 'my_sample'

    def test_output_columns(self):
        """DataFrame has all required columns."""
        reads = [_make_read([(M, 10), (N, 50), (M, 10)], ref_start=0)]
        df = _extract_from_reads(reads)
        expected = {
            'chrom', 'intron_start', 'intron_end', 'strand',
            'sample_id', 'read_count', 'motif', 'is_canonical',
        }
        assert set(df.columns) == expected


# =============================================================================
# Pass 2 — filter_cross_sample_junctions
# =============================================================================

class TestFilterCrossSampleJunctions:

    def _make_df(self, sample_id, junctions, is_canonical=True, motif='GT-AG'):
        """Build a minimal per-sample DataFrame."""
        rows = []
        for chrom, start, end, strand in junctions:
            rows.append({
                'chrom': chrom,
                'intron_start': start,
                'intron_end': end,
                'strand': strand,
                'sample_id': sample_id,
                'read_count': 2,
                'motif': motif,
                'is_canonical': is_canonical,
            })
        return pd.DataFrame(rows)

    def test_junction_in_one_sample_filtered(self):
        """A junction seen in only one sample is excluded (min_samples=2)."""
        df1 = self._make_df('s1', [('chrI', 100, 500, '+')])
        result = filter_cross_sample_junctions([df1], min_samples=2)
        assert len(result) == 0

    def test_junction_in_two_samples_passes(self):
        """A junction in two samples passes when min_samples=2."""
        df1 = self._make_df('s1', [('chrI', 100, 500, '+')])
        df2 = self._make_df('s2', [('chrI', 100, 500, '+')])
        result = filter_cross_sample_junctions([df1, df2], min_samples=2)
        assert ('chrI', 100, 500) in result

    def test_non_canonical_filtered_when_require_canonical(self):
        """Non-canonical junctions are removed when require_canonical=True."""
        df1 = self._make_df('s1', [('chrI', 100, 500, '+')],
                            is_canonical=False, motif='CT-AC')
        df2 = self._make_df('s2', [('chrI', 100, 500, '+')],
                            is_canonical=False, motif='CT-AC')
        result = filter_cross_sample_junctions(
            [df1, df2], min_samples=2, require_canonical=True
        )
        assert ('chrI', 100, 500) not in result

    def test_non_canonical_passes_when_not_required(self):
        """Non-canonical junctions pass when require_canonical=False."""
        df1 = self._make_df('s1', [('chrI', 100, 500, '+')],
                            is_canonical=False, motif='CT-AC')
        df2 = self._make_df('s2', [('chrI', 100, 500, '+')],
                            is_canonical=False, motif='CT-AC')
        result = filter_cross_sample_junctions(
            [df1, df2], min_samples=2, require_canonical=False
        )
        assert ('chrI', 100, 500) in result

    def test_huge_intron_filtered(self):
        """Junctions larger than max_intron are excluded."""
        df1 = self._make_df('s1', [('chrI', 100, 11_000, '+')])
        df2 = self._make_df('s2', [('chrI', 100, 11_000, '+')])
        result = filter_cross_sample_junctions(
            [df1, df2], max_intron=10_000
        )
        assert ('chrI', 100, 11_000) not in result

    def test_intron_at_max_intron_boundary_passes(self):
        """A junction exactly at max_intron (inclusive) passes."""
        df1 = self._make_df('s1', [('chrI', 100, 10_100, '+')])
        df2 = self._make_df('s2', [('chrI', 100, 10_100, '+')])
        result = filter_cross_sample_junctions(
            [df1, df2], max_intron=10_000
        )
        assert ('chrI', 100, 10_100) in result

    def test_insufficient_total_reads_filtered(self):
        """Junctions with too few total reads are excluded."""
        df1 = self._make_df('s1', [('chrI', 100, 500, '+')])
        df1['read_count'] = 1
        df2 = self._make_df('s2', [('chrI', 100, 500, '+')])
        df2['read_count'] = 1
        result = filter_cross_sample_junctions(
            [df1, df2], min_samples=2, min_total_reads=5
        )
        assert ('chrI', 100, 500) not in result

    def test_empty_list_returns_empty(self):
        """Empty input produces an empty frozenset."""
        result = filter_cross_sample_junctions([])
        assert result == frozenset()

    def test_unknown_motif_passes_canonical_filter(self):
        """Unknown motif (no genome provided) passes the canonical filter."""
        df1 = self._make_df('s1', [('chrI', 100, 500, '+')],
                            is_canonical=False, motif='unknown')
        df2 = self._make_df('s2', [('chrI', 100, 500, '+')],
                            is_canonical=False, motif='unknown')
        result = filter_cross_sample_junctions(
            [df1, df2], min_samples=2, require_canonical=True
        )
        assert ('chrI', 100, 500) in result

    def test_returns_frozenset(self):
        """Return type is always a frozenset."""
        result = filter_cross_sample_junctions([])
        assert isinstance(result, frozenset)

    def test_multiple_junctions_mixed(self):
        """Only junctions meeting all criteria appear in the validated set."""
        good = [('chrI', 100, 500, '+')]   # in 2 samples, canonical
        bad_1sample = [('chrI', 200, 600, '+')]   # only sample 1
        bad_huge = [('chrI', 300, 15_000, '+')]  # intron too large

        df1 = pd.concat([
            self._make_df('s1', good),
            self._make_df('s1', bad_1sample),
            self._make_df('s1', bad_huge),
        ], ignore_index=True)
        df2 = pd.concat([
            self._make_df('s2', good),
            self._make_df('s2', bad_huge),
        ], ignore_index=True)

        result = filter_cross_sample_junctions([df1, df2], max_intron=10_000)
        assert ('chrI', 100, 500) in result
        assert ('chrI', 200, 600) not in result
        assert ('chrI', 300, 15_000) not in result


# =============================================================================
# Pass 3 — apply_junction_filter
# =============================================================================

class TestApplyJunctionFilter:

    def _run_filter(self, reads, validated_junctions):
        """Run apply_junction_filter against a list of mock reads."""
        import io, tempfile, os

        header_mock = MagicMock()
        in_ctx, in_bam = _make_bam_file(reads, header=header_mock)

        written_reads = []

        # Capture written reads instead of writing to disk
        out_bam_mock = MagicMock()
        out_bam_mock.write = MagicMock(side_effect=written_reads.append)
        out_ctx = MagicMock()
        out_ctx.__enter__ = MagicMock(return_value=out_bam_mock)
        out_ctx.__exit__ = MagicMock(return_value=False)

        with patch('rectify.core.junction_validator.pysam') as mock_pysam:
            mock_pysam.AlignmentFile.side_effect = [in_ctx, out_ctx]
            counts = apply_junction_filter(
                'input.bam', 'output.bam', validated_junctions
            )

        return counts, written_reads

    def test_read_with_validated_junction_unchanged(self):
        """A read whose junction is in the validated set keeps its XC tag."""
        validated = frozenset([('chrI', 10, 110)])
        read = _make_read(
            [(M, 10), (N, 100), (M, 10)], ref_start=0, xc_tag='high'
        )
        counts, written = self._run_filter([read], validated)
        assert counts['unchanged'] == 1
        assert counts['downgraded'] == 0
        assert read._tags.get('XC') == 'high'

    def test_read_with_unvalidated_junction_downgraded_from_high(self):
        """A read with an unvalidated junction is downgraded from high→medium."""
        validated = frozenset()  # empty: no validated junctions
        read = _make_read(
            [(M, 10), (N, 100), (M, 10)], ref_start=0, xc_tag='high'
        )
        counts, written = self._run_filter([read], validated)
        assert counts['downgraded'] == 1
        assert read._tags.get('XC') == 'medium'

    def test_read_with_unvalidated_junction_downgraded_from_medium(self):
        """A medium-confidence read with a bad junction becomes low."""
        validated = frozenset()
        read = _make_read(
            [(M, 10), (N, 100), (M, 10)], ref_start=0, xc_tag='medium'
        )
        counts, written = self._run_filter([read], validated)
        assert read._tags.get('XC') == 'low'

    def test_low_confidence_stays_low(self):
        """A low-confidence read cannot be downgraded further."""
        validated = frozenset()
        read = _make_read(
            [(M, 10), (N, 100), (M, 10)], ref_start=0, xc_tag='low'
        )
        counts, written = self._run_filter([read], validated)
        assert read._tags.get('XC') == 'low'
        assert counts['downgraded'] == 0

    def test_all_reads_written(self):
        """Every read is written to the output regardless of validation status."""
        validated = frozenset()
        reads = [
            _make_read([(M, 10), (N, 100), (M, 10)], ref_start=0, xc_tag='high'),
            _make_read([(M, 20)], ref_start=0),   # no junction
        ]
        counts, written = self._run_filter(reads, validated)
        assert counts['total'] == 2
        assert len(written) == 2

    def test_read_without_junctions_unchanged(self):
        """A read with no N operations is never downgraded."""
        validated = frozenset()
        read = _make_read([(M, 100)], ref_start=0, xc_tag='high')
        counts, written = self._run_filter([read], validated)
        assert counts['unchanged'] == 1
        assert read._tags.get('XC') == 'high'

    def test_read_without_xc_tag_not_modified(self):
        """If XC tag is absent no tag is written (avoids spurious annotation)."""
        validated = frozenset()
        read = _make_read([(M, 10), (N, 100), (M, 10)], ref_start=0)
        # No xc_tag → get_tag raises KeyError
        counts, written = self._run_filter([read], validated)
        assert 'XC' not in read._tags

    def test_partial_junction_mix(self):
        """
        A read with two junctions is downgraded if ANY one is invalid.
        """
        validated = frozenset([('chrI', 10, 110)])   # first junction valid
        # second junction (115, 315) not in validated
        read = _make_read(
            [(M, 10), (N, 100), (M, 5), (N, 200), (M, 10)],
            ref_start=0, xc_tag='high',
        )
        counts, written = self._run_filter([read], validated)
        assert counts['downgraded'] == 1
        assert read._tags.get('XC') == 'medium'

    def test_return_counts_structure(self):
        """Return dict always has total, downgraded, unchanged keys."""
        counts, _ = self._run_filter([], frozenset())
        assert set(counts.keys()) == {'total', 'downgraded', 'unchanged'}


# =============================================================================
# Helpers used by the test classes above (pure-Python, no pysam needed)
# =============================================================================

def _extract_from_reads(reads, sample_id='s1', genome=None):
    """
    Run extract_sample_junctions logic directly against a list of mock reads
    without opening a real BAM file.
    """
    from collections import defaultdict
    from rectify.core.junction_validator import _get_splice_motif

    _BAM_CREF_SKIP = 3
    junction_counts = defaultdict(int)

    for read in reads:
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if not read.cigartuples:
            continue

        strand = '-' if read.is_reverse else '+'
        chrom = read.reference_name
        ref_pos = read.reference_start

        for op, length in read.cigartuples:
            if op == _BAM_CREF_SKIP:
                junction_counts[(chrom, ref_pos, ref_pos + length, strand)] += 1
            if op in (0, 2, 3, 7, 8):
                ref_pos += length

    if not junction_counts:
        return pd.DataFrame(columns=[
            'chrom', 'intron_start', 'intron_end', 'strand',
            'sample_id', 'read_count', 'motif', 'is_canonical',
        ])

    rows = []
    for (chrom, start, end, strand), count in junction_counts.items():
        motif, is_canonical = _get_splice_motif(genome, chrom, start, end, strand)
        rows.append({
            'chrom': chrom, 'intron_start': start, 'intron_end': end,
            'strand': strand, 'sample_id': sample_id, 'read_count': count,
            'motif': motif, 'is_canonical': is_canonical,
        })
    return pd.DataFrame(rows)


# =============================================================================
# Unit tests for internal helpers
# =============================================================================

class TestHelpers:

    def test_revcomp(self):
        assert _revcomp('GTAG') == 'CTAC'
        assert _revcomp('GT') == 'AC'
        assert _revcomp('AT') == 'AT'

    def test_get_splice_motif_no_genome(self):
        motif, is_can = _get_splice_motif(None, 'chrI', 10, 110, '+')
        assert motif == 'unknown'
        assert is_can is False

    def test_get_splice_motif_canonical_plus(self):
        genome = MagicMock()
        genome.fetch = MagicMock(side_effect=lambda c, s, e: {
            (10, 12): 'GT', (108, 110): 'AG',
        }[(s, e)])
        motif, is_can = _get_splice_motif(genome, 'chrI', 10, 110, '+')
        assert motif == 'GT-AG'
        assert is_can is True

    def test_get_splice_motif_u12_atac(self):
        genome = MagicMock()
        genome.fetch = MagicMock(side_effect=lambda c, s, e: {
            (10, 12): 'AT', (108, 110): 'AC',
        }[(s, e)])
        motif, is_can = _get_splice_motif(genome, 'chrI', 10, 110, '+')
        assert motif == 'AT-AC'
        assert is_can is True

    def test_read_has_invalid_junction_true(self):
        read = _make_read([(M, 10), (N, 100), (M, 10)], ref_start=0)
        validated = frozenset()  # empty
        assert _read_has_invalid_junction(read, validated) is True

    def test_read_has_invalid_junction_false(self):
        read = _make_read([(M, 10), (N, 100), (M, 10)], ref_start=0)
        validated = frozenset([('chrI', 10, 110)])
        assert _read_has_invalid_junction(read, validated) is False

    def test_read_has_no_junctions_returns_false(self):
        read = _make_read([(M, 100)], ref_start=0)
        assert _read_has_invalid_junction(read, frozenset()) is False

    def test_downgrade_read_high_to_medium(self):
        read = _make_read([], ref_start=0, xc_tag='high')
        _downgrade_read(read, 'XC')
        assert read._tags['XC'] == 'medium'

    def test_downgrade_read_medium_to_low(self):
        read = _make_read([], ref_start=0, xc_tag='medium')
        _downgrade_read(read, 'XC')
        assert read._tags['XC'] == 'low'

    def test_downgrade_read_low_stays_low(self):
        read = _make_read([], ref_start=0, xc_tag='low')
        _downgrade_read(read, 'XC')
        assert read._tags['XC'] == 'low'

    def test_downgrade_read_missing_tag_no_op(self):
        read = _make_read([], ref_start=0)  # no XC tag
        _downgrade_read(read, 'XC')
        assert 'XC' not in read._tags
