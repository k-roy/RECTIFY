"""
Unit tests for rectify.visualize.metagene strand handling.

Tests cover:
- Correct plus strand signal extraction (peak at 0)
- Correct minus strand signal extraction (peak matches plus strand)
- Detection of missing signal reversal (the canonical strand bug)
- Detection of off-by-one center coordinate errors
- verify_strand_balance() pass/fail/warn behavior
- compute_center_profile() end-to-end with strand verification

Synthetic test data design:
- Create a PositionIndex with a signal spike at a known genomic position
- For plus strand loci: spike at (center + offset) → profile should peak at offset
- For minus strand loci: spike at (center - offset) genomically → after correct
  handling (flip + reverse), profile should also peak at +offset from center

Author: Kevin R. Roy
"""

import numpy as np
import pandas as pd
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from rectify.visualize.metagene import (
    PositionIndex,
    MetagenePipeline,
    StrandOrientationError,
    verify_strand_balance,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_index_with_spike(chrom, strand, spike_position, count=100):
    """Build a PositionIndex with a single signal spike at spike_position."""
    df = pd.DataFrame([{
        'chrom': chrom,
        'strand': strand,
        'position': spike_position,
    }] * count)
    return PositionIndex(df, position_col='position')


def make_profiles_with_peak_at(peak_offset, profile_length, n_loci=50, noise=0.01):
    """
    Build synthetic per-locus profiles (2D array) with peak at peak_offset.
    peak_offset is an index into the profile array (not a bp position).
    """
    rng = np.random.default_rng(42)
    profiles = rng.uniform(0, noise, size=(n_loci, profile_length))
    profiles[:, peak_offset] += 1.0  # spike
    return profiles


# ---------------------------------------------------------------------------
# verify_strand_balance() tests
# ---------------------------------------------------------------------------

def test_verify_strand_balance_passes():
    """Both strands peak at same position — no error."""
    x = np.arange(-50, 51)
    # Both peak at index 50 → x[50] = 0
    profiles_plus = make_profiles_with_peak_at(50, len(x))
    profiles_minus = make_profiles_with_peak_at(50, len(x))

    result = verify_strand_balance(
        profiles_plus, profiles_minus, x,
        tolerance=2, raise_on_fail=True,
    )

    assert result['passed'] is True
    assert result['plus_peak'] == 0
    assert result['minus_peak'] == 0
    assert result['difference'] == 0
    assert result['n_plus'] == 50
    assert result['n_minus'] == 50


def test_verify_strand_balance_raises_on_fail():
    """Peaks diverge by 10 bp — StrandOrientationError raised."""
    x = np.arange(-50, 51)
    # Plus peaks at x[50] = 0, minus peaks at x[60] = +10
    profiles_plus = make_profiles_with_peak_at(50, len(x))
    profiles_minus = make_profiles_with_peak_at(60, len(x))

    with pytest.raises(StrandOrientationError, match="Strand orientation bug"):
        verify_strand_balance(
            profiles_plus, profiles_minus, x,
            tolerance=2, raise_on_fail=True,
        )


def test_verify_strand_balance_warn_only():
    """Peaks diverge but raise_on_fail=False — returns dict with passed=False."""
    x = np.arange(-50, 51)
    profiles_plus = make_profiles_with_peak_at(50, len(x))
    profiles_minus = make_profiles_with_peak_at(65, len(x))  # +15 bp off

    result = verify_strand_balance(
        profiles_plus, profiles_minus, x,
        tolerance=2, raise_on_fail=False,
    )

    assert result['passed'] is False
    assert result['difference'] > 2


def test_verify_strand_balance_within_tolerance():
    """Peaks differ by exactly 2 bp — should pass (tolerance=2)."""
    x = np.arange(-50, 51)
    profiles_plus = make_profiles_with_peak_at(50, len(x))
    profiles_minus = make_profiles_with_peak_at(52, len(x))  # 2 bp off

    result = verify_strand_balance(
        profiles_plus, profiles_minus, x,
        tolerance=2, raise_on_fail=True,
    )

    assert result['passed'] is True
    assert result['difference'] == 2


# ---------------------------------------------------------------------------
# compute_center_profile() extraction correctness tests
# ---------------------------------------------------------------------------

def test_plus_strand_extraction_correct():
    """
    Plus strand: signal spike at (center + 0), i.e. exactly at center.
    Expected: profile peak at x=0.
    """
    chrom = 'chrI'
    center = 1000
    spike_pos = center  # spike at center → should appear at x=0

    index = make_index_with_spike(chrom, '+', spike_pos)
    loci = [{'chrom': chrom, 'strand': '+', 'center': center}]

    pipeline = MetagenePipeline()
    result = pipeline.compute_center_profile(
        loci, index,
        window=(-50, 50),
        normalize=False,
        verify_strands=False,  # single strand, skip verification
    )

    peak_idx = np.argmax(result['profile'])
    peak_x = result['x'][peak_idx]
    assert peak_x == 0, f"Expected peak at x=0, got x={peak_x}"


def test_minus_strand_extraction_correct():
    """
    Minus strand: signal spike at genomic position (center - 0) = center.
    After correct strand handling (flip + reverse), should appear at x=0.

    For minus strand with window (-50, 50):
    - Genomic extraction: [center - 50, center + 50]
    - After [::-1]: position at genomic center maps to profile center (index 50)
    - x[50] = 0 for window (-50, 50)
    """
    chrom = 'chrI'
    center = 1000
    spike_pos = center  # genomic position at center

    index = make_index_with_spike(chrom, '-', spike_pos)
    loci = [{'chrom': chrom, 'strand': '-', 'center': center}]

    pipeline = MetagenePipeline()
    result = pipeline.compute_center_profile(
        loci, index,
        window=(-50, 50),
        normalize=False,
        verify_strands=False,
    )

    peak_idx = np.argmax(result['profile'])
    peak_x = result['x'][peak_idx]
    assert peak_x == 0, f"Expected peak at x=0 for minus strand, got x={peak_x}"


def test_plus_and_minus_strand_peaks_match():
    """
    Both strands have their signal spike at center.
    compute_center_profile() should produce matching peaks AND pass
    strand verification automatically.
    """
    chrom = 'chrIV'
    n_loci = 30
    window = (-100, 100)

    # Build index with spikes at centers for both strands
    rows = []
    loci = []
    for i in range(n_loci):
        center = 5000 + i * 500
        rows.append({'chrom': chrom, 'strand': '+', 'position': center})
        loci.append({'chrom': chrom, 'strand': '+', 'center': center})
        rows.append({'chrom': chrom, 'strand': '-', 'position': center})
        loci.append({'chrom': chrom, 'strand': '-', 'center': center})

    df = pd.DataFrame(rows * 10)  # 10 reads per spike
    index = PositionIndex(df, position_col='position')

    pipeline = MetagenePipeline()
    # verify_strands=True is default — should NOT raise
    result = pipeline.compute_center_profile(loci, index, window=window, normalize=False)

    sv = result['strand_verification']
    assert sv is not None
    assert sv['passed'] is True, (
        f"Strand verification failed: plus={sv['plus_peak']}, minus={sv['minus_peak']}"
    )
    assert sv['plus_peak'] == 0
    assert sv['minus_peak'] == 0


def test_minus_strand_missing_reversal_detected():
    """
    Simulate the canonical strand bug: extract minus strand WITHOUT [::-1].
    verify_strand_balance() must detect and raise StrandOrientationError.

    We construct profiles manually: plus peaks at x=0, minus "peaks" at x=100
    (what you'd see if signals weren't reversed — upstream becomes downstream).
    """
    x = np.arange(-100, 101)

    # Plus: correctly peaks at x=0 (index 100)
    profiles_plus = make_profiles_with_peak_at(100, len(x), n_loci=30)

    # Minus: buggy — peaks at x=+100 (index 200), as if no reversal applied
    profiles_minus = make_profiles_with_peak_at(200, len(x), n_loci=30)

    with pytest.raises(StrandOrientationError):
        verify_strand_balance(
            profiles_plus, profiles_minus, x,
            tolerance=2, raise_on_fail=True,
        )


def test_off_by_one_center_detected():
    """
    Off-by-one in center coordinate: minus strand loci have center shifted by 5 bp.
    This mimics the 1-based vs 0-based coordinate confusion.
    verify_strand_balance() must detect the >2 bp peak mismatch.
    """
    chrom = 'chrII'
    n_loci = 20
    window = (-50, 50)
    offset = 5  # deliberate off-by-one (5 bp instead of 0 or 1)

    rows = []
    plus_loci = []
    minus_loci_correct = []
    minus_loci_buggy = []

    for i in range(n_loci):
        center = 3000 + i * 300
        # Plus: spike at center, center is correct
        rows.extend([{'chrom': chrom, 'strand': '+', 'position': center}] * 10)
        plus_loci.append({'chrom': chrom, 'strand': '+', 'center': center})
        # Minus: spike also at center (genomic)
        rows.extend([{'chrom': chrom, 'strand': '-', 'position': center}] * 10)
        minus_loci_correct.append({'chrom': chrom, 'strand': '-', 'center': center})
        # Buggy minus loci: center is shifted by offset
        minus_loci_buggy.append({'chrom': chrom, 'strand': '-', 'center': center + offset})

    df = pd.DataFrame(rows)
    index = PositionIndex(df, position_col='position')
    pipeline = MetagenePipeline()
    x = np.arange(window[0], window[1] + 1)

    # Extract plus profiles (correct)
    plus_result = pipeline.compute_center_profile(
        plus_loci, index, window=window, normalize=False, verify_strands=False
    )

    # Extract minus profiles with buggy centers
    minus_result = pipeline.compute_center_profile(
        minus_loci_buggy, index, window=window, normalize=False, verify_strands=False
    )

    # Manual verification — should fail due to the 5 bp center offset
    with pytest.raises(StrandOrientationError):
        verify_strand_balance(
            plus_result['profiles_plus'],
            minus_result['profiles_minus'],
            x,
            tolerance=2,
            raise_on_fail=True,
        )
