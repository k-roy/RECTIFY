#!/usr/bin/env python3
"""
Tests for RECTIFY configuration module.
"""

import pytest
from rectify import config


class TestConfig:
    """Test configuration parameters."""

    def test_version(self):
        """Test version string."""
        assert config.__version__ == "2.5.0"

    def test_chromosome_mappings(self):
        """Test chromosome name mappings."""
        assert config.CHROM_TO_GENOME['chrI'] == 'ref|NC_001133|'
        assert config.GENOME_TO_CHROM['ref|NC_001133|'] == 'chrI'

        # Check all chromosomes
        assert len(config.CHROM_TO_GENOME) == 17  # 16 + mito
        assert len(config.GENOME_TO_CHROM) == 17
        assert len(config.CHROM_SIZES) == 17

    def test_shift_corrections(self):
        """Test shift correction lookup."""
        assert config.get_shift_from_acount(0) == 0.0
        assert config.get_shift_from_acount(4) == 1.0
        assert config.get_shift_from_acount(10) == 3.8
        assert config.get_shift_from_acount(15) == 3.8  # Saturated

    def test_thresholds(self):
        """Test threshold values are valid."""
        assert 0.0 <= config.POLYA_RICHNESS_THRESHOLD <= 1.0
        assert 0.0 <= config.AG_RICHNESS_THRESHOLD <= 1.0
        assert 0.0 <= config.INDEL_FLANK_A_THRESHOLD <= 1.0
        assert 0.0 <= config.NETSEQ_PEAK_THRESHOLD <= 1.0


class TestHelperFunctions:
    """Test configuration helper functions."""

    def test_get_shift_from_acount(self):
        """Test shift correction retrieval."""
        # Known values
        assert config.get_shift_from_acount(0) == 0.0
        assert config.get_shift_from_acount(7) == 2.6

        # Saturated value
        assert config.get_shift_from_acount(100) == config.DEFAULT_MAX_SHIFT

    def test_validate_config(self):
        """Test configuration validation."""
        # Should not raise exception
        config.validate_config()
