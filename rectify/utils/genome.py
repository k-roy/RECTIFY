#!/usr/bin/env python3
"""
Genome utilities for RECTIFY.

This module provides functions for:
- Loading and querying genome sequences
- Detecting A-tracts and T-tracts
- Calculating downstream A-counts
- Strand-aware sequence operations
- Chromosome name conversion

Consolidated from roadblocks project:
- scripts/02_nanopore_polishing/utils/genome_utils.py
- scripts/utils/sequence_utils.py

Author: Kevin R. Roy
Date: 2026-03-09
"""

from pathlib import Path
from typing import Dict, Optional, Tuple
from Bio import SeqIO

from ..config import (
    CHROM_TO_GENOME, GENOME_TO_CHROM, CHROM_SIZES,
    POLYA_RICHNESS_THRESHOLD
)

# =============================================================================
# Sequence Complement
# =============================================================================

# Use str.maketrans for efficient string translation
# This is much faster than per-character dict lookup
COMPLEMENT_TABLE = str.maketrans('ATGCatgc', 'TACGtacg')


def complement(seq: str) -> str:
    """
    Return complement of DNA sequence.

    Uses optimized str.translate for O(n) performance.

    Args:
        seq: DNA sequence (case-insensitive)

    Returns:
        Complemented sequence (preserves case, unknown bases unchanged)
    """
    return seq.translate(COMPLEMENT_TABLE)


def reverse_complement(seq: str) -> str:
    """
    Return reverse complement of DNA sequence.

    Args:
        seq: DNA sequence (case-insensitive)

    Returns:
        Reverse complemented sequence
    """
    return complement(seq)[::-1]


# =============================================================================
# Genome Loading
# =============================================================================

def load_genome(genome_path: Path) -> Dict[str, str]:
    """
    Load reference genome into memory.

    Handles both plain FASTA and gzip-compressed FASTA (.gz).

    A pickle sidecar cache (<genome>.pkl) is written alongside the FASTA on
    first load and reused on subsequent calls as long as it is newer than the
    FASTA file.  This saves 10-120 s per sample on large genomes.

    Args:
        genome_path: Path to genome FASTA file (.fa, .fasta, .fsa, or .gz)

    Returns:
        Dict mapping chromosome name (NCBI format) to sequence string
    """
    import gzip as _gzip
    import pickle as _pickle
    import logging as _logging

    _log = _logging.getLogger(__name__)
    genome_path = Path(genome_path)
    pickle_path = genome_path.with_suffix('.pkl')

    # Try cache first
    try:
        if (pickle_path.exists()
                and pickle_path.stat().st_mtime >= genome_path.stat().st_mtime):
            _log.debug("Loading genome from pickle cache: %s", pickle_path)
            with open(pickle_path, 'rb') as _fh:
                genome = _pickle.load(_fh)
            _log.debug("  Loaded %d chromosomes from cache", len(genome))
            return genome
    except Exception as _exc:
        _log.debug("Genome pickle cache unusable (%s); loading from FASTA", _exc)

    print(f"Loading genome from {genome_path}...")
    genome = {}

    if genome_path.suffix == '.gz':
        with _gzip.open(str(genome_path), 'rt') as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                genome[record.id] = str(record.seq).upper()
    else:
        for record in SeqIO.parse(str(genome_path), 'fasta'):
            genome[record.id] = str(record.seq).upper()

    print(f"  Loaded {len(genome)} chromosomes")

    # Write pickle cache for next call
    try:
        with open(pickle_path, 'wb') as _fh:
            _pickle.dump(genome, _fh, protocol=_pickle.HIGHEST_PROTOCOL)
        _log.debug("Wrote genome pickle cache: %s", pickle_path)
    except Exception as _exc:
        _log.debug("Could not write genome pickle cache (%s); continuing without cache", _exc)

    return genome


# =============================================================================
# Sequence Fetching
# =============================================================================

def fetch_genomic_sequence(genome: Dict[str, str],
                           chrom: str,
                           start: int,
                           end: int,
                           strand: str = '+') -> str:
    """
    Fetch genomic sequence in specified orientation.

    Args:
        genome: Genome dict from load_genome()
        chrom: Chromosome name (standard format: chrI, chrII, etc.)
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        strand: '+' or '-'

    Returns:
        Sequence string, or empty if invalid coordinates
    """
    # Get sequence (try canonical name first, fall back to NCBI format)
    seq = genome.get(chrom) or genome.get(CHROM_TO_GENOME.get(chrom, ''))
    if seq is None:
        return ''

    # Clamp coordinates
    start = max(0, start)
    end = min(len(seq), end)

    if start >= end:
        return ''

    # Extract sequence
    genomic_seq = seq[start:end]

    # Reverse complement if minus strand
    if strand == '-':
        genomic_seq = reverse_complement(genomic_seq)

    return genomic_seq


def get_downstream_sequence(genome: Dict[str, str],
                            chrom: str,
                            position: int,
                            strand: str,
                            downstream_bp: int = 10) -> str:
    """
    Get downstream sequence in gene orientation.

    + strand: Look RIGHT of position (higher genomic coords)
    - strand: Look LEFT of position (lower genomic coords), return in RNA orientation

    Args:
        genome: Genome dict from load_genome()
        chrom: Chromosome name (standard format)
        position: Position (0-based)
        strand: '+' or '-'
        downstream_bp: Number of bases downstream

    Returns:
        Downstream sequence in gene orientation (RNA 5'→3')
    """
    seq = genome.get(chrom) or genome.get(CHROM_TO_GENOME.get(chrom, ''))
    if seq is None:
        return ''

    if strand == '+':
        # Downstream = RIGHT = higher coords
        down_start = position + 1
        down_end = min(len(seq), position + 1 + downstream_bp)
        downstream_seq = seq[down_start:down_end]
    else:
        # Downstream in gene coords = LEFT = lower genomic coords
        down_start = max(0, position - downstream_bp)
        down_end = position
        genomic_seq = seq[down_start:down_end]
        # Reverse complement to get RNA orientation
        downstream_seq = reverse_complement(genomic_seq)

    return downstream_seq


def get_downstream_acount(genome: Dict[str, str],
                         chrom: str,
                         pos: int,
                         strand: str,
                         downstream_bp: int = 10) -> Optional[int]:
    """
    Get A-count in DOWNSTREAM region where poly(A) tail would align.

    + strand: Look RIGHT of CPA (higher genomic coords) for genomic A's
    - strand: Look LEFT of CPA (lower genomic coords) for genomic T's (RNA A's)

    Args:
        genome: Genome dict from load_genome()
        chrom: Chromosome name (standard format)
        pos: CPA position (0-based)
        strand: '+' or '-'
        downstream_bp: Number of bases to check downstream

    Returns:
        Number of A's in downstream window, or None if invalid
    """
    downstream_seq = get_downstream_sequence(genome, chrom, pos, strand, downstream_bp)

    if len(downstream_seq) < downstream_bp:
        return None

    return downstream_seq.count('A')


# =============================================================================
# A-tract Detection
# =============================================================================

def is_atract(seq: str, strand: str = '+', threshold: float = POLYA_RICHNESS_THRESHOLD) -> bool:
    """
    Check if sequence is A-tract (+ strand) or T-tract (- strand).

    Args:
        seq: Genomic sequence (in genomic orientation)
        strand: Gene strand
        threshold: Minimum fraction of A's (or T's) to classify as tract

    Returns:
        True if sequence is A-tract (+ strand) or T-tract (- strand)
    """
    if len(seq) == 0:
        return False

    if strand == '+':
        # Check for A-tract in genomic sequence
        a_frac = seq.count('A') / len(seq)
        return a_frac >= threshold
    else:
        # Check for T-tract in genomic sequence (complement of RNA A-tract)
        t_frac = seq.count('T') / len(seq)
        return t_frac >= threshold


def count_contiguous_a_tract(seq: str, strand: str = '+') -> int:
    """
    Count length of contiguous A-tract at beginning of sequence.

    Args:
        seq: Genomic sequence (in genomic orientation)
        strand: Gene strand

    Returns:
        Length of contiguous A's (+ strand) or T's (- strand) from start
    """
    target_base = 'A' if strand == '+' else 'T'
    count = 0

    for base in seq:
        if base == target_base:
            count += 1
        else:
            break

    return count


# NOTE: find_atract_boundaries has been consolidated to rectify.core.atract_detector
# Use: from rectify.core.atract_detector import find_atract_boundaries


# =============================================================================
# Chromosome Name Conversion
# =============================================================================

def standardize_chrom_name(chrom: str) -> str:
    """
    Standardize chromosome name to chrI, chrII, etc. format.

    Handles various input formats:
    - chrI, chrII, etc. (already standard)
    - chr1, chr2, etc. (numeric)
    - I, II, etc. (Roman numerals only)
    - ref|NC_001133|, etc. (NCBI format)

    Args:
        chrom: Chromosome name in any format

    Returns:
        Standardized chromosome name (chrI, chrII, etc.)
    """
    # Already standard format
    if chrom in CHROM_SIZES:
        return chrom

    # NCBI format
    if chrom in GENOME_TO_CHROM:
        return GENOME_TO_CHROM[chrom]

    # Numeric format (chr1 -> chrI)
    roman_map = {
        '1': 'I', '2': 'II', '3': 'III', '4': 'IV', '5': 'V',
        '6': 'VI', '7': 'VII', '8': 'VIII', '9': 'IX', '10': 'X',
        '11': 'XI', '12': 'XII', '13': 'XIII', '14': 'XIV', '15': 'XV', '16': 'XVI',
        'M': 'Mito', 'mt': 'Mito', 'MT': 'Mito', 'mito': 'Mito'
    }

    if chrom.startswith('chr'):
        num = chrom[3:]
        if num in roman_map:
            return 'chr' + roman_map[num]

    # Roman numeral only
    if chrom in ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII',
                 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'Mito']:
        return 'chr' + chrom

    # Return as-is if unrecognized
    return chrom


def validate_coordinates(chrom: str, start: int, end: int) -> bool:
    """
    Validate that coordinates are within chromosome bounds.

    Args:
        chrom: Chromosome name (standard format)
        start: Start position (0-based or 1-based, depending on context)
        end: End position

    Returns:
        True if coordinates are valid, False otherwise
    """
    if chrom not in CHROM_SIZES:
        return False
    if start < 0 or end < start:
        return False
    if end > CHROM_SIZES[chrom]:
        return False
    return True


def clamp_position(chrom: str, position: int) -> int:
    """
    Clamp position to valid chromosome range.

    Args:
        chrom: Chromosome name (standard format)
        position: Position to clamp (0-based)

    Returns:
        Clamped position within [0, chrom_size)
    """
    chrom_size = CHROM_SIZES.get(chrom, float('inf'))
    return max(0, min(position, chrom_size - 1))


# =============================================================================
# Position Utilities
# =============================================================================

def get_read_3prime_position(read, strand: str) -> int:
    """
    Get 3' end genomic position (strand-aware).

    Returns the position where the read alignment ends on the 3' side,
    which corresponds to the potential CPA site. Does NOT include soft-clipped
    bases (poly-A tails) since those extend beyond the true 3' end.

    Args:
        read: pysam.AlignedSegment
        strand: Gene strand ('+' or '-')

    Returns:
        Genomic coordinate of RNA 3' end (0-based, inclusive)
    """
    if strand == '+':
        # 3' end is at right (higher coordinate)
        # reference_end is 0-based exclusive, so subtract 1 for inclusive
        return read.reference_end - 1
    else:
        # 3' end is at left (lower coordinate)
        # reference_start is 0-based inclusive
        return read.reference_start
