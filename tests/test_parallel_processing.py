#!/usr/bin/env python3
"""
Tests for parallel BAM processing.

Author: Kevin R. Roy
Date: 2026-03
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
import tempfile
import os
import pysam

# Test imports
from rectify.core.bam.regions import (
    get_processing_regions,
    find_coverage_gaps,
)
from rectify.core.bam.parallel import (
    _is_region_checkpoint_complete,
    _prescan_bam_once,
    _process_region_worker,
    _rebuild_output_from_region_files,
    _region_done_path,
    _region_tsv_path,
    _regions_from_intervals,
    _write_region_results_atomic,
)
from rectify.core.bam.processing_stats import ProcessingStats


class TestGetProcessingRegions:
    """Tests for region splitting logic."""

    def test_small_chromosomes_single_region(self):
        """Small chromosomes should be single regions."""
        # Mock BAM with small chromosomes
        with patch('pysam.AlignmentFile') as mock_bam:
            mock_bam.return_value.__enter__ = Mock(return_value=mock_bam.return_value)
            mock_bam.return_value.__exit__ = Mock(return_value=False)
            mock_bam.return_value.references = ['chrI', 'chrII']
            mock_bam.return_value.lengths = [100000, 200000]  # Both under max_region_size

            with patch('rectify.core.bam.regions.find_coverage_gaps', return_value=[]):
                regions = get_processing_regions(
                    'test.bam',
                    min_gap_size=10000,
                    max_region_size=500000
                )

                # Both chromosomes should be single regions
                assert len(regions) == 2
                assert regions[0] == ('chrI', 0, 100000)
                assert regions[1] == ('chrII', 0, 200000)


class TestFindCoverageGaps:
    """Tests for coverage gap detection."""

    def test_empty_bam_returns_full_gap(self):
        """Empty BAM should return one gap covering whole chromosome."""
        with patch('pysam.AlignmentFile') as mock_bam:
            mock_bam.return_value.__enter__ = Mock(return_value=mock_bam.return_value)
            mock_bam.return_value.__exit__ = Mock(return_value=False)
            mock_bam.return_value.references = ['chrI']
            mock_bam.return_value.lengths = [100000]
            mock_bam.return_value.fetch.return_value = []  # No reads

            gaps = find_coverage_gaps('test.bam', 'chrI', min_gap_size=10000)

            # Should return one gap covering the whole chromosome
            assert len(gaps) == 1
            assert gaps[0] == (0, 100000)


class TestProcessRegionWorker:
    """Tests for region worker function."""

    def test_skips_unmapped_reads(self):
        """Worker should skip unmapped reads."""
        # Create mock read
        mock_read = Mock()
        mock_read.is_unmapped = True
        mock_read.is_secondary = False
        mock_read.is_supplementary = False

        with patch('pysam.AlignmentFile') as mock_bam:
            mock_bam.return_value.__enter__ = Mock(return_value=mock_bam.return_value)
            mock_bam.return_value.__exit__ = Mock(return_value=False)
            mock_bam.return_value.fetch.return_value = [mock_read]

            results = _process_region_worker(
                region=('chrI', 0, 100000),
                bam_path='test.bam',
                genome={'chrI': 'ATCG' * 25000},  # Mock genome
                apply_atract=False,
                apply_ag_mispriming=False,
                apply_polya_trim=False,
                apply_indel_correction=False,
                netseq_dir=None
            )

            # Should return empty results (unmapped read skipped)
            assert len(results) == 0

    def test_skips_secondary_alignments(self):
        """Worker should skip secondary alignments."""
        mock_read = Mock()
        mock_read.is_unmapped = False
        mock_read.is_secondary = True
        mock_read.is_supplementary = False

        with patch('pysam.AlignmentFile') as mock_bam:
            mock_bam.return_value.__enter__ = Mock(return_value=mock_bam.return_value)
            mock_bam.return_value.__exit__ = Mock(return_value=False)
            mock_bam.return_value.fetch.return_value = [mock_read]

            results = _process_region_worker(
                region=('chrI', 0, 100000),
                bam_path='test.bam',
                genome={'chrI': 'ATCG' * 25000},
                apply_atract=False,
                apply_ag_mispriming=False,
                apply_polya_trim=False,
                apply_indel_correction=False,
                netseq_dir=None
            )

            assert len(results) == 0


class TestStreamingOutput:
    """Tests for streaming output mode."""

    def test_generates_summary_from_stats(self):
        """Summary generation from stats dict should work."""
        from rectify.core.bam.output import generate_summary_from_stats

        stats = {
            'total_reads': 1000,
            'with_ambiguity': 250,
            'position_corrected': 100,
            'by_confidence': {'high': 800, 'medium': 150, 'low': 50},
        }

        report = generate_summary_from_stats(stats)

        assert 'Total reads:' in report
        assert '1,000' in report
        assert '25.0%' in report  # 250/1000 with ambiguity
        assert 'high' in report
        assert 'Streaming Mode' in report

    def test_empty_stats_handled(self):
        """Empty stats should return appropriate message."""
        from rectify.core.bam.output import generate_summary_from_stats

        stats = {
            'total_reads': 0,
            'with_ambiguity': 0,
            'position_corrected': 0,
            'by_confidence': {'high': 0, 'medium': 0, 'low': 0},
        }

        report = generate_summary_from_stats(stats)
        assert 'No reads processed' in report


class TestWriteResultsChunk:
    """Tests for chunk writing."""

    def test_writes_correct_format(self):
        """Chunk writer should produce correct TSV format."""
        from rectify.core.bam.parallel import _write_results_chunk
        from io import StringIO

        results = [{
            'read_id': 'read1',
            'chrom': 'chrI',
            'strand': '+',
            'original_3prime': 100,
            'corrected_3prime': 95,
            'ambiguity_min': 90,
            'ambiguity_max': 100,
            'ambiguity_range': 10,
            'correction_applied': ['atract_ambiguity'],
            'confidence': 'high',
            'qc_flags': ['PASS'],
        }]

        output = StringIO()
        _write_results_chunk(output, results)

        written = output.getvalue()
        assert 'read1' in written
        assert 'chrI' in written
        assert '+' in written
        assert 'atract_ambiguity' in written
        assert 'PASS' in written

    def test_chunk_writer_matches_canonical_schema_width(self):
        """Streaming chunks must carry every column needed by BAM CIGAR surgery."""
        from rectify.core.bam.output import CORRECTION_TSV_HEADER
        from rectify.core.bam.parallel import _write_results_chunk
        from io import StringIO

        result = {
            'read_id': 'read1',
            'chrom': 'chrI',
            'strand': '+',
            'original_3prime': 100,
            'corrected_3prime': 95,
            'ambiguity_min': 95,
            'ambiguity_max': 100,
            'ambiguity_range': 5,
            'correction_applied': ['overcall_rescue'],
            'confidence': 'high',
            'qc_flags': ['PASS'],
            'oc_homopolymer_extension': 3,
            'oc_overcall_count': 2,
            'oc_terminal_base': 'G',
            'five_prime_upstream_trim': 4,
            'reanchor_clip_len': 5,
        }

        output = StringIO()
        _write_results_chunk(output, [result])
        fields = output.getvalue().rstrip('\n').split('\t')

        assert len(fields) == len(CORRECTION_TSV_HEADER)
        by_col = dict(zip(CORRECTION_TSV_HEADER, fields))
        assert by_col['oc_homopolymer_extension'] == '3'
        assert by_col['oc_overcall_count'] == '2'
        assert by_col['oc_terminal_base'] == 'G'
        assert by_col['five_prime_upstream_trim'] == '4'
        assert by_col['reanchor_clip_len'] == '5'


class TestCheckpointRegionFiles:
    """Tests for checkpointed parallel output durability."""

    def _result(self, read_id, chrom='chrI', corrected=100, fraction=1.0):
        return {
            'read_id': read_id,
            'chrom': chrom,
            'strand': '+',
            'original_3prime': corrected,
            'corrected_3prime': corrected,
            'ambiguity_min': corrected,
            'ambiguity_max': corrected,
            'ambiguity_range': 0,
            'correction_applied': [],
            'confidence': 'high',
            'qc_flags': ['PASS'],
            'fraction': fraction,
        }

    def test_region_done_requires_committed_tsv(self, tmp_path):
        """A .done sentinel alone must not be trusted as a completed region."""
        done_path = _region_done_path(tmp_path, 0)
        done_path.write_text('done\n')
        tmp_region = tmp_path / 'region_0000.tsv.tmp'
        tmp_region.write_text('partial row\n')

        assert not _is_region_checkpoint_complete(tmp_path, 0)

        _region_tsv_path(tmp_path, 0).write_text('')
        assert _is_region_checkpoint_complete(tmp_path, 0)

    def test_region_results_commit_tsv_stats_then_done(self, tmp_path):
        """Region commits should leave a trusted TSV, stats sidecar, and sentinel."""
        _write_region_results_atomic(tmp_path, 2, [self._result('read2')])

        assert _is_region_checkpoint_complete(tmp_path, 2)
        assert _region_tsv_path(tmp_path, 2).read_text().startswith('read2\tchrI\t+')
        assert not (tmp_path / 'region_0002.tsv.tmp').exists()
        assert _region_done_path(tmp_path, 2).read_text() == 'done\n'

    def test_rebuild_output_orders_regions_and_position_counts(self, tmp_path):
        """Final output is rebuilt in region order, independent of commit order."""
        _write_region_results_atomic(tmp_path, 1, [self._result('readB', corrected=200)])
        _write_region_results_atomic(tmp_path, 0, [self._result('readA', corrected=100)])

        prescan_stats = ProcessingStats(total_reads_in_bam=2)
        output_path = tmp_path / 'corrected_reads.tsv'

        stats, pos_counts = _rebuild_output_from_region_files(
            tmp_path,
            str(output_path),
            2,
            prescan_stats,
        )

        lines = output_path.read_text().splitlines()
        assert lines[0].startswith('read_id\tchrom\tstrand')
        assert lines[1].startswith('readA\tchrI\t+')
        assert lines[2].startswith('readB\tchrI\t+')
        assert stats.total_reads_in_bam == 2
        assert stats.reads_processed == 2
        assert pos_counts == {
            ('chrI', 100, '+'): 1.0,
            ('chrI', 200, '+'): 1.0,
        }


class TestParallelPrescan:
    def test_regions_from_intervals_splits_at_large_gaps(self):
        regions = _regions_from_intervals(
            {'chrI': 250000},
            {'chrI': [(100, 200), (50000, 50100), (200000, 200100)]},
            min_gap_size=10000,
        )

        assert regions == [
            ('chrI', 0, 200),
            ('chrI', 50000, 50100),
            ('chrI', 200000, 200100),
        ]

    def test_prescan_bam_once_counts_stats_and_plans_regions(self, tmp_path):
        bam_path = tmp_path / 'reads.bam'
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [{'SN': 'chrI', 'LN': 200000}],
        }
        with pysam.AlignmentFile(str(bam_path), 'wb', header=header) as out:
            read = pysam.AlignedSegment()
            read.query_name = 'read1'
            read.query_sequence = 'A' * 20
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 10
            read.mapping_quality = 60
            read.cigartuples = [(0, 20)]
            read.query_qualities = pysam.qualitystring_to_array('I' * 20)
            out.write(read)

            unmapped = pysam.AlignedSegment()
            unmapped.query_name = 'unmapped'
            unmapped.query_sequence = 'A' * 10
            unmapped.flag = 4
            unmapped.reference_id = -1
            unmapped.reference_start = -1
            unmapped.mapping_quality = 0
            unmapped.cigartuples = None
            unmapped.query_qualities = pysam.qualitystring_to_array('I' * 10)
            out.write(unmapped)

        stats, regions, rescue = _prescan_bam_once(
            str(bam_path),
            {'chrI': 'A' * 200000},
            min_gap_size=10000,
            variant_aware=False,
        )

        assert stats.total_reads_in_bam == 2
        assert stats.reads_unmapped == 1
        assert regions == [('chrI', 0, 30)]
        assert rescue is None
