#!/usr/bin/env python3
"""
Configuration constants for RECTIFY.

This module centralizes:
- Chromosome mappings (standard ↔ NCBI format)
- Shift correction parameters (from NET-seq analysis)
- Poly(A) tail model parameters
- Indel detection parameters
- NET-seq refinement parameters

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict

# =============================================================================
# Version
# =============================================================================

__version__ = "2.5.0"

# =============================================================================
# Chromosome Mappings
# =============================================================================

# Standard chromosome names (UCSC/Roman numeral format - canonical for RECTIFY)
# All BAM files aligned with RECTIFY bundled genome will use these names
CANONICAL_CHROMS: list[str] = [
    'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII',
    'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI',
    'chrMito',
]

# NCBI to canonical (chrI) format mapping (for legacy BAM files)
NCBI_TO_CHROM: Dict[str, str] = {
    'ref|NC_001133|': 'chrI',
    'ref|NC_001134|': 'chrII',
    'ref|NC_001135|': 'chrIII',
    'ref|NC_001136|': 'chrIV',
    'ref|NC_001137|': 'chrV',
    'ref|NC_001138|': 'chrVI',
    'ref|NC_001139|': 'chrVII',
    'ref|NC_001140|': 'chrVIII',
    'ref|NC_001141|': 'chrIX',
    'ref|NC_001142|': 'chrX',
    'ref|NC_001143|': 'chrXI',
    'ref|NC_001144|': 'chrXII',
    'ref|NC_001145|': 'chrXIII',
    'ref|NC_001146|': 'chrXIV',
    'ref|NC_001147|': 'chrXV',
    'ref|NC_001148|': 'chrXVI',
    'ref|NC_001224|': 'chrMito',
}

# Canonical to NCBI format mapping (for compatibility with legacy workflows)
CHROM_TO_NCBI: Dict[str, str] = {v: k for k, v in NCBI_TO_CHROM.items()}

# Legacy aliases (deprecated - use NCBI_TO_CHROM and CHROM_TO_NCBI instead)
CHROM_TO_GENOME: Dict[str, str] = CHROM_TO_NCBI
GENOME_TO_CHROM: Dict[str, str] = NCBI_TO_CHROM

# Chromosome sizes (sacCer3/R64)
CHROM_SIZES: Dict[str, int] = {
    'chrI': 230218,
    'chrII': 813184,
    'chrIII': 316620,
    'chrIV': 1531933,
    'chrV': 576874,
    'chrVI': 270161,
    'chrVII': 1090940,
    'chrVIII': 562643,
    'chrIX': 439888,
    'chrX': 745751,
    'chrXI': 666816,
    'chrXII': 1078177,
    'chrXIII': 924431,
    'chrXIV': 784333,
    'chrXV': 1091291,
    'chrXVI': 948066,
    'chrMito': 85779,
}

# =============================================================================
# Shift Correction Parameters (from NET-seq Analysis)
# =============================================================================

# These represent how much the apparent position is shifted RIGHTWARD (downstream)
# due to long poly(A) tails (>15 A's for nanopore) aligning to genomic A's/T's
# To correct: shift LEFTWARD (upstream) by this amount
SHIFT_CORRECTIONS_BY_ACOUNT: Dict[int, float] = {
    0: 0.0,   # No downstream A's - no shift
    1: 0.2,   # Minimal shift
    2: 0.3,
    3: 0.4,
    4: 1.0,   # Moderate shift
    5: 1.3,
    6: 1.7,
    7: 2.6,   # Strong shift
    8: 2.8,
    9: 2.9,
    10: 3.8,  # Maximum observed shift (saturated)
}

# For A-counts > 10, use the 10A value (saturated)
DEFAULT_MAX_SHIFT: float = 3.8

# =============================================================================
# Poly(A) Tail Model Parameters
# =============================================================================

# Minimum tail length for nanopore oligo-dT priming
MIN_POLYA_LENGTH: int = 15

# A-richness threshold for identifying poly(A) tails
POLYA_RICHNESS_THRESHOLD: float = 0.8  # 80% A content

# Window size for A-richness calculation
POLYA_WINDOW_SIZE: int = 10

# RTA adapter patterns to detect
ADAPTER_POLY_T_MIN: int = 6          # Minimum poly(T) length
ADAPTER_TC_MOTIFS = ['TC', 'TCTC', 'TCT']

# Scoring thresholds for tail classification
TAIL_SCORE_HIGH: float = 0.7    # Definite poly(A) tail
TAIL_SCORE_LOW: float = 0.4     # Uncertain

# =============================================================================
# Indel Detection Parameters
# =============================================================================

# Maximum distance from 3' end to consider deletions as artifacts
INDEL_SEARCH_WINDOW: int = 20  # bp

# Maximum deletion size to consider as artifact
INDEL_MAX_SIZE: int = 3  # bp

# Minimum A-richness in flanking regions to classify as artifact
INDEL_FLANK_A_THRESHOLD: float = 0.7  # 70% A content

# Minimum flank length to check for A-richness
INDEL_MIN_FLANK_LENGTH: int = 5  # bp

# =============================================================================
# Ambiguity Range Parameters
# =============================================================================

# Maximum shift to consider when calculating ambiguity ranges
MAX_AMBIGUITY_SHIFT: int = 5  # bp

# Window size for downstream A-count calculation
DOWNSTREAM_WINDOW_SIZE: int = 10  # bp

# =============================================================================
# NET-seq Refinement Parameters
# =============================================================================

# Agreement threshold for combining WT and dst1 NET-seq data
NETSEQ_AGREEMENT_THRESHOLD: float = 0.90  # 90% agreement at ±1bp

# Window around CPA for peak detection
NETSEQ_PEAK_WINDOW: int = 5  # ±5bp

# Minimum signal threshold for peak calling (relative to max)
NETSEQ_PEAK_THRESHOLD: float = 0.5  # 50% of max signal

# Minimum peak signal for high confidence
NETSEQ_SIGNAL_HIGH: float = 1.0
NETSEQ_SIGNAL_MEDIUM: float = 0.5

# Maximum distance between peaks to consider them "close"
NETSEQ_PEAK_CLOSE_DISTANCE: int = 2  # bp

# =============================================================================
# NET-seq Oligo-A Parameters (Empirically Derived from 0A Sites)
# =============================================================================
#
# These constants were derived from analysis of oligo-adenylated soft-clips at
# 0A sites (CPA positions with no downstream genomic adenosines) in:
#   - GSE25107 (Churchman 2011): 17.8% reads with detectable oligo-A
#   - GSE159603 (Couvillion 2022, hexamer-trimmed): 6.2% reads with oligo-A
#
# Both datasets show identical oligo-A length distributions at 0A sites,
# validating the robustness of this empirical model.

# Mean oligo-A tail length at 0A sites (no downstream A-tract ambiguity)
# Derived from 5000+ 0A CPA sites, both 2011 and 2022 NET-seq datasets
NETSEQ_OLIGO_A_MEAN_0A: float = 5.52  # bp

# Median oligo-A tail length at 0A sites
NETSEQ_OLIGO_A_MEDIAN_0A: float = 5.0  # bp

# Fraction of reads with detectable oligo-A at 0A sites (capture rate)
# Lower rate in 2022 due to hexamer trimming removing some short tails
NETSEQ_OLIGO_A_CAPTURE_RATE_2011: float = 0.178  # 17.8%
NETSEQ_OLIGO_A_CAPTURE_RATE_2022: float = 0.062  # 6.2%

# 0A Point-Spread-Function parameters
# At 0A sites: ~54% of signal at true CPA, ~46% spreads downstream
NETSEQ_PSF_POSITION_0_FRACTION: float = 0.54  # Signal at true CPA position
NETSEQ_PSF_DOWNSTREAM_FRACTION: float = 0.46  # Signal spreading downstream

# Peak shape statistics by A-tract length (from peak_shape_statistics.tsv)
# Format: (right_signal_fraction, asymmetry)
# Asymmetry > 0 indicates rightward (downstream) spreading
NETSEQ_PEAK_SHAPE_0A = (0.61, 0.26)      # Moderate rightward spread
NETSEQ_PEAK_SHAPE_1_3A = (0.408, -0.16)  # Slight leftward bias
NETSEQ_PEAK_SHAPE_4_6A = (0.472, -0.02)  # Balanced
NETSEQ_PEAK_SHAPE_7PLUS_A = (0.92, 0.85) # Severe rightward spread (deconvolution target)

# NNLS deconvolution regularization strength
# Higher values smooth the solution; lower values preserve peaks
NETSEQ_DECONV_REGULARIZATION: float = 0.01

# =============================================================================
# AG Mispriming Parameters (from original RECTIFY)
# =============================================================================

# Window size for AG-richness calculation
AG_RICHNESS_WINDOW: int = 50  # bp downstream

# AG content threshold for flagging likely mispriming
AG_RICHNESS_THRESHOLD: float = 0.65  # 65% A+G content

# Minimum window size if near chromosome end
AG_RICHNESS_MIN_WINDOW: int = 20  # bp

# =============================================================================
# Helper Functions
# =============================================================================

def get_shift_from_acount(a_count: int) -> float:
    """
    Get shift correction for given A-count.

    Args:
        a_count: Number of A's in downstream window

    Returns:
        Shift correction in bp (how much to shift leftward/upstream)
    """
    if a_count <= 10:
        return SHIFT_CORRECTIONS_BY_ACOUNT.get(a_count, 0.0)
    else:
        return DEFAULT_MAX_SHIFT


def validate_config():
    """Validate configuration parameters."""
    assert len(CHROM_TO_GENOME) == len(GENOME_TO_CHROM), "Chromosome mapping mismatch"
    assert len(CHROM_TO_GENOME) == len(CHROM_SIZES), "Chromosome size mapping incomplete"
    assert 0.0 <= POLYA_RICHNESS_THRESHOLD <= 1.0, "Invalid poly(A) richness threshold"
    assert 0.0 <= AG_RICHNESS_THRESHOLD <= 1.0, "Invalid AG richness threshold"
    assert MIN_POLYA_LENGTH > 0, "Invalid minimum poly(A) length"
    assert INDEL_SEARCH_WINDOW > 0, "Invalid indel search window"


# Run validation on import
validate_config()
