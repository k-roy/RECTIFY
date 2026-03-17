#!/usr/bin/env python3
"""
Alignment utilities for RECTIFY.

This module provides functions for:
- CIGAR parsing and operations
- Coordinate conversion (read ↔ genomic)
- Soft-clip extraction
- Deletion detection and analysis

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import List, Tuple, Dict, Optional
import pysam

# CIGAR operation codes
CIGAR_OPS = {
    0: 'M',  # Match/Mismatch
    1: 'I',  # Insertion
    2: 'D',  # Deletion
    3: 'N',  # Skipped region (splicing)
    4: 'S',  # Soft clip
    5: 'H',  # Hard clip
    6: 'P',  # Padding
    7: '=',  # Sequence match
    8: 'X',  # Sequence mismatch
}

# Operations that consume reference
CONSUMES_REF = {0, 2, 3, 7, 8}  # M, D, N, =, X

# Operations that consume query
CONSUMES_QUERY = {0, 1, 4, 7, 8}  # M, I, S, =, X


# =============================================================================
# CIGAR Parsing
# =============================================================================

def parse_cigar(cigar_tuples: List[Tuple[int, int]]) -> str:
    """
    Convert pysam CIGAR tuples to CIGAR string.

    Args:
        cigar_tuples: List of (operation, length) tuples from pysam

    Returns:
        CIGAR string (e.g., "100M", "90M10S")
    """
    if not cigar_tuples:
        return ""

    cigar_str = ""
    for op, length in cigar_tuples:
        cigar_str += f"{length}{CIGAR_OPS[op]}"

    return cigar_str


def get_cigar_stats(cigar_tuples: List[Tuple[int, int]]) -> Dict[str, int]:
    """
    Get statistics from CIGAR string.

    Args:
        cigar_tuples: List of (operation, length) tuples

    Returns:
        Dict with counts for each operation type
    """
    stats = {op: 0 for op in CIGAR_OPS.values()}

    for op_code, length in cigar_tuples:
        op_char = CIGAR_OPS[op_code]
        stats[op_char] += length

    return stats


# =============================================================================
# Soft-clip Extraction
# =============================================================================

def extract_soft_clips(read: pysam.AlignedSegment) -> List[Dict]:
    """
    Extract soft-clipped sequences from read.

    Args:
        read: pysam AlignedSegment

    Returns:
        List of soft-clip dicts with keys:
            - side: 'left' or 'right'
            - seq: soft-clipped sequence
            - qual: soft-clipped quality scores
            - start: genomic start position of clip
            - length: length of clip
    """
    soft_clips = []

    if not read.cigartuples:
        return soft_clips

    read_seq = read.query_sequence
    read_qual = read.query_qualities

    # Check left (5') soft-clip
    if read.cigartuples[0][0] == 4:  # Soft-clip
        clip_length = read.cigartuples[0][1]
        soft_clips.append({
            'side': 'left',
            'seq': read_seq[:clip_length],
            'qual': read_qual[:clip_length] if read_qual else None,
            'start': read.reference_start - clip_length,
            'length': clip_length,
        })

    # Check right (3') soft-clip
    if read.cigartuples[-1][0] == 4:  # Soft-clip
        clip_length = read.cigartuples[-1][1]
        soft_clips.append({
            'side': 'right',
            'seq': read_seq[-clip_length:],
            'qual': read_qual[-clip_length:] if read_qual else None,
            'start': read.reference_end,
            'length': clip_length,
        })

    return soft_clips


def has_soft_clip(read: pysam.AlignedSegment, side: str = 'any') -> bool:
    """
    Check if read has soft-clipped bases.

    Args:
        read: pysam AlignedSegment
        side: 'left', 'right', or 'any'

    Returns:
        True if read has soft-clip on specified side
    """
    if not read.cigartuples:
        return False

    if side in ('left', 'any'):
        if read.cigartuples[0][0] == 4:
            return True

    if side in ('right', 'any'):
        if read.cigartuples[-1][0] == 4:
            return True

    return False


# =============================================================================
# Deletion Detection
# =============================================================================

def extract_deletions(read: pysam.AlignedSegment) -> List[Dict]:
    """
    Extract all deletions from CIGAR string.

    Args:
        read: pysam AlignedSegment

    Returns:
        List of deletion dicts with keys:
            - ref_pos: genomic position of deletion start
            - read_pos: position in read (0-based)
            - length: deletion length (bp)
            - ref_seq: deleted sequence (if available)
            - distance_from_3prime: distance from 3' end of alignment
    """
    if not read.cigartuples:
        return []

    deletions = []
    ref_pos = read.reference_start
    query_pos = 0

    # Try to get reference sequence (if MD tag available)
    read_seq = read.query_sequence or ''
    ref_seq = read.get_reference_sequence() if hasattr(read, 'get_reference_sequence') else None

    for i, (op, length) in enumerate(read.cigartuples):
        if op == 2:  # Deletion
            # Calculate distance from 3' end
            remaining_ops = read.cigartuples[i+1:]
            distance_from_3prime = sum(
                length for op, length in remaining_ops
                if op in CONSUMES_REF
            )

            # Get deleted sequence if possible
            deleted_seq = ''
            if ref_seq:
                # Calculate position in reference sequence
                ref_offset = ref_pos - read.reference_start
                deleted_seq = ref_seq[ref_offset:ref_offset + length]

            deletions.append({
                'ref_pos': ref_pos,
                'read_pos': query_pos,
                'length': length,
                'ref_seq': deleted_seq,
                'distance_from_3prime': distance_from_3prime,
                'cigar_index': i,
            })

        # Update positions
        if op in CONSUMES_REF:
            ref_pos += length
        if op in CONSUMES_QUERY:
            query_pos += length

    return deletions


def extract_insertions(read: pysam.AlignedSegment) -> List[Dict]:
    """
    Extract all insertions from CIGAR string.

    Args:
        read: pysam AlignedSegment

    Returns:
        List of insertion dicts with keys:
            - ref_pos: genomic position where insertion occurs
            - read_pos: position in read (0-based)
            - length: insertion length (bp)
            - seq: inserted sequence
            - distance_from_3prime: distance from 3' end of alignment
    """
    if not read.cigartuples:
        return []

    insertions = []
    ref_pos = read.reference_start
    query_pos = 0
    read_seq = read.query_sequence or ''

    for i, (op, length) in enumerate(read.cigartuples):
        if op == 1:  # Insertion
            # Calculate distance from 3' end
            remaining_ops = read.cigartuples[i+1:]
            distance_from_3prime = sum(
                length for op, length in remaining_ops
                if op in CONSUMES_REF
            )

            # Get inserted sequence
            inserted_seq = read_seq[query_pos:query_pos + length]

            insertions.append({
                'ref_pos': ref_pos,
                'read_pos': query_pos,
                'length': length,
                'seq': inserted_seq,
                'distance_from_3prime': distance_from_3prime,
                'cigar_index': i,
            })

        # Update positions
        if op in CONSUMES_REF:
            ref_pos += length
        if op in CONSUMES_QUERY:
            query_pos += length

    return insertions


def get_deletions_near_3prime(read: pysam.AlignedSegment,
                              window: int = 20) -> List[Dict]:
    """
    Get deletions near 3' end of read alignment.

    Args:
        read: pysam AlignedSegment
        window: Maximum distance from 3' end (bp)

    Returns:
        List of deletion dicts within window of 3' end
    """
    all_deletions = extract_deletions(read)

    near_3prime = [
        d for d in all_deletions
        if d['distance_from_3prime'] <= window
    ]

    return near_3prime


# =============================================================================
# Coordinate Conversion
# =============================================================================

def read_to_genomic_coord(read: pysam.AlignedSegment,
                         read_pos: int) -> Optional[int]:
    """
    Convert read coordinate to genomic coordinate.

    Args:
        read: pysam AlignedSegment
        read_pos: Position in read (0-based)

    Returns:
        Genomic coordinate (0-based), or None if position is in soft-clip/insertion
    """
    if not read.cigartuples:
        return None

    ref_pos = read.reference_start
    query_pos = 0

    for op, length in read.cigartuples:
        # Check if read_pos is in this operation
        if query_pos <= read_pos < query_pos + length:
            if op in CONSUMES_REF and op in CONSUMES_QUERY:  # M, =, X
                offset = read_pos - query_pos
                return ref_pos + offset
            else:
                # Position is in soft-clip or insertion
                return None

        # Update positions
        if op in CONSUMES_REF:
            ref_pos += length
        if op in CONSUMES_QUERY:
            query_pos += length

    return None


def genomic_to_read_coord(read: pysam.AlignedSegment,
                          genomic_pos: int) -> Optional[int]:
    """
    Convert genomic coordinate to read coordinate.

    Args:
        read: pysam AlignedSegment
        genomic_pos: Genomic position (0-based)

    Returns:
        Read coordinate (0-based), or None if position is in deletion/intron
    """
    if not read.cigartuples:
        return None

    ref_pos = read.reference_start
    query_pos = 0

    for op, length in read.cigartuples:
        # Check if genomic_pos is in this operation
        if ref_pos <= genomic_pos < ref_pos + length:
            if op in CONSUMES_REF and op in CONSUMES_QUERY:  # M, =, X
                offset = genomic_pos - ref_pos
                return query_pos + offset
            else:
                # Position is in deletion or intron
                return None

        # Update positions
        if op in CONSUMES_REF:
            ref_pos += length
        if op in CONSUMES_QUERY:
            query_pos += length

    return None


# =============================================================================
# Alignment Metrics
# =============================================================================

def get_alignment_length(cigar_tuples: List[Tuple[int, int]]) -> int:
    """
    Get total alignment length (reference span).

    Args:
        cigar_tuples: List of (operation, length) tuples

    Returns:
        Total length on reference
    """
    return sum(length for op, length in cigar_tuples if op in CONSUMES_REF)


def get_query_length(cigar_tuples: List[Tuple[int, int]]) -> int:
    """
    Get total query length (read length including soft-clips).

    Args:
        cigar_tuples: List of (operation, length) tuples

    Returns:
        Total length on query
    """
    return sum(length for op, length in cigar_tuples if op in CONSUMES_QUERY)


def get_match_length(cigar_tuples: List[Tuple[int, int]]) -> int:
    """
    Get total match length (M, =, X operations).

    Args:
        cigar_tuples: List of (operation, length) tuples

    Returns:
        Total match length
    """
    return sum(
        length for op, length in cigar_tuples
        if op in (0, 7, 8)  # M, =, X
    )


def get_alignment_identity(read: pysam.AlignedSegment) -> float:
    """
    Calculate alignment identity (fraction of matches).

    Args:
        read: pysam AlignedSegment

    Returns:
        Alignment identity (0.0-1.0)
    """
    if not read.cigartuples:
        return 0.0

    match_length = get_match_length(read.cigartuples)
    query_length = get_query_length(read.cigartuples)

    if query_length == 0:
        return 0.0

    return match_length / query_length


# =============================================================================
# Strand Utilities
# =============================================================================

def is_reverse_strand(read: pysam.AlignedSegment) -> bool:
    """
    Check if read is on reverse strand.

    Args:
        read: pysam AlignedSegment

    Returns:
        True if read is reverse strand
    """
    return read.is_reverse


def get_strand_symbol(read: pysam.AlignedSegment) -> str:
    """
    Get strand symbol ('+' or '-').

    Args:
        read: pysam AlignedSegment

    Returns:
        '+' for forward, '-' for reverse
    """
    return '-' if read.is_reverse else '+'
