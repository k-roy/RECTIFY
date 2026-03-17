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

# Test imports
from rectify.core.bam_processor import (
    get_processing_regions,
    find_coverage_gaps,
    _process_region_worker,
)


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

            with patch('rectify.core.bam_processor.find_coverage_gaps', return_value=[]):
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
        from rectify.core.bam_processor import generate_summary_from_stats

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
        from rectify.core.bam_processor import generate_summary_from_stats

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
        from rectify.core.bam_processor import _write_results_chunk
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
