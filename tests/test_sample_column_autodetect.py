"""
Tests for sample column auto-detection in load_corrected_positions.

Bug 2: Pipeline fails with ValueError when input data uses "replicate" (or other
common variants) instead of "sample" as the column name.

The fix adds auto-detection in both the small-file and large-file (chunked) paths
of load_corrected_positions / _load_large_file_chunked.
"""

import io
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from rectify.core.analyze_command import load_corrected_positions


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_CHROMS = ['chrI', 'chrI', 'chrII', 'chrII', 'chrI']
_STRANDS = ['+', '+', '-', '-', '+']
_POSITIONS = [100, 200, 300, 400, 500]


def _make_tsv(sample_col_name: str, tmp_path: Path) -> Path:
    """Write a minimal corrected-positions TSV with the given sample column name."""
    df = pd.DataFrame({
        'chrom': _CHROMS,
        'strand': _STRANDS,
        'corrected_position': _POSITIONS,
        sample_col_name: ['s1', 's1', 's2', 's2', 's1'],
    })
    tsv = tmp_path / f"positions_{sample_col_name}.tsv"
    df.to_csv(tsv, sep='\t', index=False)
    return tsv


# ---------------------------------------------------------------------------
# Tests: small-file path (< 0.5 GB — always true for these tiny fixtures)
# ---------------------------------------------------------------------------

class TestSampleColumnAutodetect:

    def test_explicit_sample_column_works(self, tmp_path):
        """Baseline: explicit --sample-column sample still works."""
        tsv = _make_tsv('sample', tmp_path)
        df = load_corrected_positions(str(tsv), sample_column='sample',
                                      normalize_chroms=False)
        assert 'sample' in df.columns
        assert set(df['sample']) == {'s1', 's2'}

    def test_replicate_column_autodetected(self, tmp_path):
        """Bug 2 scenario: data uses 'replicate'; default --sample-column sample."""
        tsv = _make_tsv('replicate', tmp_path)
        # Should NOT raise — auto-detects 'replicate'
        df = load_corrected_positions(str(tsv), sample_column='sample',
                                      normalize_chroms=False)
        # Column name in the returned DataFrame is the detected name
        assert 'replicate' in df.columns
        assert set(df['replicate']) == {'s1', 's2'}

    def test_sample_id_column_autodetected(self, tmp_path):
        """Data uses 'sample_id'; default --sample-column sample."""
        tsv = _make_tsv('sample_id', tmp_path)
        df = load_corrected_positions(str(tsv), sample_column='sample',
                                      normalize_chroms=False)
        assert 'sample_id' in df.columns
        assert set(df['sample_id']) == {'s1', 's2'}

    def test_sample_name_column_autodetected(self, tmp_path):
        """Data uses 'sample_name'; default --sample-column sample."""
        tsv = _make_tsv('sample_name', tmp_path)
        df = load_corrected_positions(str(tsv), sample_column='sample',
                                      normalize_chroms=False)
        assert 'sample_name' in df.columns

    def test_condition_column_autodetected(self, tmp_path):
        """Data uses 'condition' (last fallback); default --sample-column sample."""
        tsv = _make_tsv('condition', tmp_path)
        df = load_corrected_positions(str(tsv), sample_column='sample',
                                      normalize_chroms=False)
        assert 'condition' in df.columns

    def test_no_recognisable_sample_column_raises(self, tmp_path):
        """When no known sample column is present, a clear ValueError is raised."""
        df = pd.DataFrame({
            'chrom': _CHROMS,
            'strand': _STRANDS,
            'corrected_position': _POSITIONS,
            'experiment': ['e1'] * 5,   # not in the fallback list
        })
        tsv = tmp_path / "no_sample_col.tsv"
        df.to_csv(tsv, sep='\t', index=False)

        with pytest.raises(ValueError, match="not found"):
            load_corrected_positions(str(tsv), sample_column='sample',
                                     normalize_chroms=False)

    def test_explicit_replicate_arg_works_directly(self, tmp_path):
        """Passing --sample-column replicate explicitly still works as before."""
        tsv = _make_tsv('replicate', tmp_path)
        df = load_corrected_positions(str(tsv), sample_column='replicate',
                                      normalize_chroms=False)
        assert 'replicate' in df.columns
