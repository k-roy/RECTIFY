#!/usr/bin/env python3
"""
Tests for apa_detection module.

Basic smoke tests verifying detect_apa_isoforms, build_gene_apa_profiles,
quantify_apa_usage, and summarize_apa_detection on a small synthetic dataset.
"""

import pytest
from rectify.core.analyze.apa_detection import (
    detect_apa_isoforms,
    build_gene_apa_profiles,
    quantify_apa_usage,
    summarize_apa_detection,
    HAS_PANDAS,
)
from rectify.core.unified_record import UnifiedReadRecord


@pytest.fixture
def tfc3_records():
    """Six reads for YAL001C (TFC3), two distinct TES clusters."""
    return [
        # TES cluster at ~3000
        UnifiedReadRecord(
            read_id="read_001", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=3000, three_prime_corrected=2995,
            alignment_start=1000, alignment_end=3000,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
        UnifiedReadRecord(
            read_id="read_002", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=3000, three_prime_corrected=2998,
            alignment_start=1000, alignment_end=3000,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
        UnifiedReadRecord(
            read_id="read_003", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=3000, three_prime_corrected=3002,
            alignment_start=1000, alignment_end=3000,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
        # TES cluster at ~2500 (APA site)
        UnifiedReadRecord(
            read_id="read_004", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=2500, three_prime_corrected=2495,
            alignment_start=1000, alignment_end=2500,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
        UnifiedReadRecord(
            read_id="read_005", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=2500, three_prime_corrected=2498,
            alignment_start=1000, alignment_end=2500,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
        UnifiedReadRecord(
            read_id="read_006", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=2500, three_prime_corrected=2502,
            alignment_start=1000, alignment_end=2500,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
    ]


@pytest.fixture
def tfc3_attributions():
    return {f"read_00{i}": ("YAL001C", "TFC3") for i in range(1, 7)}


def test_detect_apa_isoforms_basic(tfc3_records, tfc3_attributions):
    isoforms = detect_apa_isoforms(
        tfc3_records, tfc3_attributions, min_reads_per_isoform=3
    )
    assert len(isoforms) == 2
    read_counts = sorted(iso.n_reads for iso in isoforms)
    assert read_counts == [3, 3]


def test_detect_apa_isoforms_min_reads_filter(tfc3_records, tfc3_attributions):
    # Raise min_reads so only one cluster qualifies (neither does since both have exactly 3)
    isoforms = detect_apa_isoforms(
        tfc3_records, tfc3_attributions, min_reads_per_isoform=4
    )
    assert len(isoforms) == 0


def test_build_gene_apa_profiles(tfc3_records, tfc3_attributions):
    isoforms = detect_apa_isoforms(
        tfc3_records, tfc3_attributions, min_reads_per_isoform=3
    )
    profiles = build_gene_apa_profiles(isoforms)
    assert "YAL001C" in profiles
    profile = profiles["YAL001C"]
    assert profile.n_isoforms == 2
    assert profile.has_apa is True
    usage = profile.get_tes_usage()
    total_fraction = sum(usage.values())
    assert abs(total_fraction - 1.0) < 0.01


@pytest.mark.skipif(not HAS_PANDAS, reason="pandas not available")
def test_quantify_apa_usage(tfc3_records, tfc3_attributions):
    isoforms = detect_apa_isoforms(
        tfc3_records, tfc3_attributions, min_reads_per_isoform=3
    )
    df = quantify_apa_usage(isoforms)
    assert not df.empty
    assert "isoform_id" in df.columns


def test_summarize_apa_detection(tfc3_records, tfc3_attributions):
    isoforms = detect_apa_isoforms(
        tfc3_records, tfc3_attributions, min_reads_per_isoform=3
    )
    profiles = build_gene_apa_profiles(isoforms)
    summary = summarize_apa_detection(isoforms, profiles)
    assert "n_isoforms" in summary
    assert summary["n_isoforms"] == 2
