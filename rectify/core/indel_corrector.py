#!/usr/bin/env python3
"""
Indel artifact correction for RECTIFY.

This module implements detection and correction of insertion/deletion artifacts
at A-tracts. These artifacts occur when aligners try to force-align poly(A) tails
to genomic A-tracts, creating spurious indels.

Problem: Aligners (minimap2, BWA, STAR) create deletion artifacts when aligning
poly(A) tails to genomic sequences. Example:
    Read:    ATTAAAAAA
    Genomic: ATTTAAAAA
    Aligned: ATT--A--AAAA (aligner removes genomic TT)

This shifts the apparent 3' end position and must be corrected.

Two correction approaches are provided:
1. `detect_indel_artifacts()` + `correct_position_for_indels()`:
   Detects small indels near 3' end and corrects for them

2. `find_polya_boundary()`:
   Genome-aware algorithm that walks backwards from mapped 3' end,
   comparing read to genome, to find the TRUE CPA (cleavage and
   polyadenylation) site where they agree on a non-A base.

The `find_polya_boundary()` approach is recommended for nanopore data as it
handles complex poly(A) alignment artifacts including multiple deletions,
mismatches (base-calling errors), and genomic A-tracts.

Author: Kevin R. Roy
Date: 2026-03-09
Updated: 2026-03-16 - Added genome-aware poly-A boundary detection
"""

from typing import Dict, List, Optional
import pysam

from ..config import (
    INDEL_MAX_SIZE,
    INDEL_SEARCH_WINDOW,
    INDEL_FLANK_A_THRESHOLD,
    INDEL_MIN_FLANK_LENGTH,
)
from ..utils.alignment import extract_deletions, extract_insertions
from ..utils.genome import fetch_genomic_sequence
from ..config import CHROM_TO_GENOME


# =============================================================================
# Genome-aware Poly-A Boundary Detection
# =============================================================================

def find_polya_boundary(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Dict[str, str],
    min_polya_len: int = 5
) -> Optional[Dict]:
    """
    Find the true 3' end by comparing mRNA and genome sequences.

    This is the recommended approach for nanopore data. It handles complex
    poly(A) alignment artifacts by walking backwards from the mapped 3' end
    and finding the first position where the read and genome AGREE on a non-A base.

    Algorithm:
    1. Build aligned positions from CIGAR (read_base, genome_base pairs)
    2. Start at the mapped 3' end
    3. Walk BACKWARDS toward the mRNA body
    4. Find the FIRST non-A base where genome and mRNA AGREE
    5. That's the true CPA

    Example (+ strand):
        mRNA:    CTGACGATGAAAAAGAAAATAAACAA-AAAAAAA
        genome:  CTGACGATGAAAAAGAAAATAAA-AAGAAAAAAA
                                   ^   ^
                                   T   C (mismatch: genome=A, read=C -> seq error)

        Walking backwards from 3' end:
        - Positions in poly-A region: A=A (match), A=A (match), etc.
        - Position 23: C vs A (MISMATCH) -> sequencing error, keep scanning
        - Position 19: T vs T (MATCH on non-A) -> TRUE CPA!

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'
        genome: Dict of chromosome sequences (NCBI format keys)
        min_polya_len: Minimum poly-A length to consider for correction (default 5)

    Returns:
        Dict with:
            - corrected_pos: True CPA position (0-based)
            - original_pos: Original mapped 3' position
            - polya_aligned_bp: Number of bp in the aligned poly-A region
            - correction_bp: How much the position was corrected
        Or None if no poly-A alignment artifacts found
    """
    cigar = read.cigartuples
    seq = read.query_sequence
    if not cigar or not seq:
        return None

    # Get genomic sequence
    chrom = read.reference_name
    genome_chrom = CHROM_TO_GENOME.get(chrom, chrom)
    genome_seq = genome.get(genome_chrom, '')
    if not genome_seq:
        return None

    # CIGAR ops: M=0, I=1, D=2, N=3, S=4, H=5

    if strand == '+':
        # For + strand, 3' end is at right side
        # Walk BACKWARDS from 3' end to find first non-A agreement

        ref_pos = read.reference_start
        read_pos = 0

        # Skip 5' soft-clip
        if cigar[0][0] == 4:
            read_pos = cigar[0][1]

        # Build list of aligned positions with both read and ref bases
        # (read_pos, ref_pos, read_base, genome_base)
        aligned_positions = []

        for op, length in cigar:
            if op == 4:  # Soft-clip
                continue
            elif op == 0 or op == 7 or op == 8:  # M, =, X - match/mismatch
                for i in range(length):
                    if ref_pos < len(genome_seq):
                        read_base = seq[read_pos].upper() if read_pos < len(seq) else 'N'
                        genome_base = genome_seq[ref_pos].upper()
                        aligned_positions.append((read_pos, ref_pos, read_base, genome_base))
                    read_pos += 1
                    ref_pos += 1
            elif op == 1:  # Insertion - consumes query only
                read_pos += length
            elif op == 2:  # Deletion - consumes reference only
                # Mark deletions (read has no base, genome has base)
                for i in range(length):
                    if ref_pos < len(genome_seq):
                        genome_base = genome_seq[ref_pos].upper()
                        aligned_positions.append((None, ref_pos, None, genome_base))
                    ref_pos += 1

        if not aligned_positions:
            return None

        # Walk BACKWARDS from the 3' end (last aligned position)
        # Find the FIRST non-A base where genome and read AGREE
        raw_3prime = read.reference_end - 1
        true_cpa_ref_pos = None

        for i in range(len(aligned_positions) - 1, -1, -1):
            rp, refp, rb, gb = aligned_positions[i]
            if rp is None:  # Skip deletions
                continue
            # Check: do genome and read agree on a non-A base?
            if rb == gb and gb != 'A':
                true_cpa_ref_pos = refp
                break

        if true_cpa_ref_pos is not None:
            polya_len = raw_3prime - true_cpa_ref_pos
            # Only return correction if we found significant poly-A
            if polya_len >= min_polya_len:
                return {
                    'corrected_pos': true_cpa_ref_pos,
                    'original_pos': raw_3prime,
                    'polya_aligned_bp': polya_len,
                    'correction_bp': raw_3prime - true_cpa_ref_pos,
                }

    else:  # strand == '-'
        # For - strand, 3' end is at left side (reference_start)
        # Poly-A appears as T's in genomic coordinates
        # Walk FORWARD from 3' end to find first non-T agreement

        ref_pos = read.reference_start
        read_pos = 0

        # Skip 3' soft-clip (at left for - strand)
        if cigar[0][0] == 4:
            read_pos = cigar[0][1]

        # Build aligned positions
        aligned_positions = []

        for op, length in cigar:
            if op == 4:  # Soft-clip
                continue
            elif op == 0 or op == 7 or op == 8:  # M, =, X
                for i in range(length):
                    if ref_pos < len(genome_seq):
                        read_base = seq[read_pos].upper() if read_pos < len(seq) else 'N'
                        genome_base = genome_seq[ref_pos].upper()
                        aligned_positions.append((read_pos, ref_pos, read_base, genome_base))
                    read_pos += 1
                    ref_pos += 1
            elif op == 1:  # Insertion
                read_pos += length
            elif op == 2:  # Deletion
                for i in range(length):
                    if ref_pos < len(genome_seq):
                        genome_base = genome_seq[ref_pos].upper()
                        aligned_positions.append((None, ref_pos, None, genome_base))
                    ref_pos += 1

        if not aligned_positions:
            return None

        # Walk FORWARD from the 3' end (first aligned position)
        # Find the FIRST non-T base where genome and read AGREE
        raw_3prime = read.reference_start
        true_cpa_ref_pos = None

        for i in range(len(aligned_positions)):
            rp, refp, rb, gb = aligned_positions[i]
            if rp is None:  # Skip deletions
                continue
            # Check: do genome and read agree on a non-T base?
            if rb == gb and gb != 'T':
                true_cpa_ref_pos = refp
                break

        if true_cpa_ref_pos is not None:
            polya_len = true_cpa_ref_pos - raw_3prime
            # Only return correction if we found significant poly-A
            if polya_len >= min_polya_len:
                return {
                    'corrected_pos': true_cpa_ref_pos,
                    'original_pos': raw_3prime,
                    'polya_aligned_bp': polya_len,
                    'correction_bp': true_cpa_ref_pos - raw_3prime,
                }

    return None


# =============================================================================
# Small Indel Detection (Original Approach)
# =============================================================================

def is_atract_deletion(
    deletion: Dict,
    genome_seq: str,
    read_seq: str,
    strand: str
) -> bool:
    """
    Determine if a deletion is an A-tract alignment artifact.

    Artifact criteria (all aligners):
    1. Deletion size ≤ INDEL_MAX_SIZE (default 3bp)
    2. Deleted bases are NOT A's (+ strand) or T's (- strand)
    3. Flanking regions in read are A-rich (>INDEL_FLANK_A_THRESHOLD)

    Args:
        deletion: Deletion dict with 'pos', 'length', 'ref_seq' keys
        genome_seq: Genomic sequence at deletion site
        read_seq: Read sequence flanking deletion
        strand: Gene strand ('+' or '-')

    Returns:
        True if deletion is likely an A-tract artifact
    """
    # Check deletion size
    if deletion['length'] > INDEL_MAX_SIZE:
        return False

    # Get deleted bases
    deleted_bases = deletion.get('ref_seq', '').upper()
    if not deleted_bases:
        return False

    # Check if deleted bases are A's/T's (these are expected in A-tracts)
    target_base = 'A' if strand == '+' else 'T'
    if all(b == target_base for b in deleted_bases):
        # Deletion is all A's/T's - this is expected, not an artifact
        return False

    # Check flanking regions are A/T-rich (depending on strand)
    # This indicates the read has poly(A) sequence at this position
    # For minus strand reads, poly(A) appears as poly(T) in genomic orientation
    if len(read_seq) < INDEL_MIN_FLANK_LENGTH * 2:
        return False

    # Get flanks around deletion position in read
    del_pos = deletion['read_pos']
    left_start = max(0, del_pos - INDEL_MIN_FLANK_LENGTH)
    left_flank = read_seq[left_start:del_pos]
    right_flank = read_seq[del_pos:del_pos + INDEL_MIN_FLANK_LENGTH]

    # Check richness of flanks (A for + strand, T for - strand)
    flank_base = 'A' if strand == '+' else 'T'

    if len(left_flank) >= INDEL_MIN_FLANK_LENGTH:
        left_richness = left_flank.count(flank_base) / len(left_flank)
        if left_richness < INDEL_FLANK_A_THRESHOLD:
            return False

    if len(right_flank) >= INDEL_MIN_FLANK_LENGTH:
        right_richness = right_flank.count(flank_base) / len(right_flank)
        if right_richness < INDEL_FLANK_A_THRESHOLD:
            return False

    return True


def detect_indel_artifacts(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Optional[Dict[str, str]] = None
) -> List[Dict]:
    """
    Detect indel artifacts near the 3' end of read.

    Scans CIGAR for deletions/insertions within INDEL_SEARCH_WINDOW of 3' end.
    Classifies each as artifact or legitimate variation.

    Args:
        read: pysam AlignedSegment
        strand: Gene strand ('+' or '-')
        genome: Optional genome dict for sequence validation

    Returns:
        List of artifact dicts with:
            - type: 'deletion' or 'insertion'
            - pos: Position in read (0-based)
            - ref_pos: Position in reference
            - length: Indel size
            - is_artifact: True if classified as artifact
    """
    artifacts = []

    # Determine 3' end position in read
    read_len = read.query_length
    if strand == '+':
        # 3' end is at right (end of read)
        three_prime_pos = read_len
        search_start = max(0, three_prime_pos - INDEL_SEARCH_WINDOW)
        search_end = three_prime_pos
    else:
        # 3' end is at left (start of read)
        three_prime_pos = 0
        search_start = 0
        search_end = min(read_len, INDEL_SEARCH_WINDOW)

    # Extract deletions
    deletions = extract_deletions(read)
    for deletion in deletions:
        # Check if deletion is near 3' end
        del_pos = deletion['read_pos']
        if search_start <= del_pos <= search_end:
            # Check if it's an artifact
            read_seq = read.query_sequence
            genome_seq = ''
            if genome:
                # Get genomic sequence for validation
                chrom = read.reference_name
                ref_pos = deletion['ref_pos']
                genome_seq = fetch_genomic_sequence(
                    genome, chrom, ref_pos - 10, ref_pos + 10
                )

            is_artifact = is_atract_deletion(
                deletion, genome_seq, read_seq, strand
            )

            artifacts.append({
                'type': 'deletion',
                'pos': del_pos,
                'ref_pos': deletion['ref_pos'],
                'length': deletion['length'],
                'ref_seq': deletion.get('ref_seq', ''),
                'is_artifact': is_artifact,
            })

    # Extract insertions (less common but can occur)
    insertions = extract_insertions(read)
    for insertion in insertions:
        ins_pos = insertion['read_pos']
        if search_start <= ins_pos <= search_end:
            # Insertions near 3' end are suspicious if small
            is_artifact = insertion['length'] <= INDEL_MAX_SIZE

            artifacts.append({
                'type': 'insertion',
                'pos': ins_pos,
                'ref_pos': insertion['ref_pos'],
                'length': insertion['length'],
                'seq': insertion.get('seq', ''),
                'is_artifact': is_artifact,
            })

    return artifacts


def correct_position_for_indels(
    original_pos: int,
    indel_artifacts: List[Dict],
    strand: str
) -> Dict:
    """
    Correct 3' end position by accounting for DELETION artifacts only.

    IMPORTANT: Only deletions affect reference coordinates. Insertions add bases
    to the read but do NOT shift the reference position. Therefore, only deletion
    artifacts require position correction.

    For deletions in homopolymer regions (nanopore error):
    - The aligner skips genomic bases that should be in the alignment
    - This shifts reference_end to a higher coordinate than the true 3' end
    - Correction: move 3' position UPSTREAM by deletion length

    Args:
        original_pos: Original 3' end position (0-based)
        indel_artifacts: List of artifact dicts from detect_indel_artifacts()
        strand: Gene strand ('+' or '-')

    Returns:
        Dict with:
            - corrected_position: Position after removing artifacts
            - correction_bp: Number of bp corrected
            - n_deletions: Number of deletion artifacts corrected
            - n_insertions: Number of insertion artifacts (flagged only, no correction)
    """
    correction_bp = 0
    n_deletions = 0
    n_insertions = 0

    for artifact in indel_artifacts:
        if not artifact['is_artifact']:
            continue

        if artifact['type'] == 'deletion':
            # Deletion artifact: aligner skipped genomic bases
            # The reference coordinates are shifted downstream by deletion length
            # Correct by moving position UPSTREAM:
            #   + strand: upstream = leftward = subtract
            #   - strand: upstream = rightward = add
            if strand == '+':
                correction_bp -= artifact['length']
            else:
                correction_bp += artifact['length']
            n_deletions += 1

        elif artifact['type'] == 'insertion':
            # Insertion artifact: extra bases in read that aren't in reference
            # IMPORTANT: Insertions do NOT affect reference coordinates!
            # The reference_start and reference_end are unchanged by insertions.
            # We count them for QC but do NOT apply position correction.
            n_insertions += 1
            # No correction_bp change for insertions

    return {
        'corrected_position': original_pos + correction_bp,
        'correction_bp': correction_bp,
        'n_deletions': n_deletions,
        'n_insertions': n_insertions,
    }


def correct_indels_from_read(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Optional[Dict[str, str]] = None
) -> Dict:
    """
    Full indel artifact correction for a single read.

    This function:
    1. Detects indel artifacts near 3' end
    2. Calculates corrected 3' end position
    3. Returns summary statistics

    Args:
        read: pysam AlignedSegment
        strand: Gene strand ('+' or '-')
        genome: Optional genome dict for validation

    Returns:
        Dict with:
            - original_3prime: Original 3' end position (0-based)
            - corrected_3prime: Corrected position
            - correction_bp: Position shift (bp)
            - has_artifacts: True if artifacts detected
            - artifacts: List of artifact dicts
    """
    # Get original 3' end position
    if strand == '+':
        original_3prime = read.reference_end - 1  # 0-based inclusive
    else:
        original_3prime = read.reference_start  # 0-based

    # Detect artifacts
    artifacts = detect_indel_artifacts(read, strand, genome)

    # Count artifacts
    artifact_list = [a for a in artifacts if a['is_artifact']]
    has_artifacts = len(artifact_list) > 0

    # Correct position
    if has_artifacts:
        correction = correct_position_for_indels(
            original_3prime, artifact_list, strand
        )
        corrected_3prime = correction['corrected_position']
        correction_bp = correction['correction_bp']
    else:
        corrected_3prime = original_3prime
        correction_bp = 0

    return {
        'original_3prime': original_3prime,
        'corrected_3prime': corrected_3prime,
        'correction_bp': correction_bp,
        'has_artifacts': has_artifacts,
        'artifacts': artifact_list,
    }


def calculate_indel_statistics(indel_results: list) -> Dict:
    """
    Calculate summary statistics for indel correction.

    Args:
        indel_results: List of result dicts from correct_indels_from_read()

    Returns:
        Dict with summary statistics
    """
    if not indel_results:
        return {
            'total': 0,
            'with_artifacts': 0,
            'artifact_rate': 0.0,
            'mean_correction': 0.0,
        }

    import numpy as np

    with_artifacts = [r for r in indel_results if r['has_artifacts']]
    corrections = [abs(r['correction_bp']) for r in with_artifacts]

    return {
        'total': len(indel_results),
        'with_artifacts': len(with_artifacts),
        'artifact_rate': len(with_artifacts) / len(indel_results),
        'mean_correction': np.mean(corrections) if corrections else 0.0,
        'median_correction': np.median(corrections) if corrections else 0.0,
        'max_correction': max(corrections) if corrections else 0,
    }


def format_indel_report(stats: Dict) -> str:
    """
    Format indel correction statistics as human-readable report.

    Args:
        stats: Statistics dict from calculate_indel_statistics()

    Returns:
        Formatted report string
    """
    report = []
    report.append("=" * 60)
    report.append("Indel Artifact Correction Summary")
    report.append("=" * 60)
    report.append("")

    report.append("Overall:")
    report.append(f"  Total reads:            {stats['total']:,}")
    report.append(f"  With artifacts:         {stats['with_artifacts']:,} ({stats['artifact_rate']:.1%})")
    report.append("")

    if stats['with_artifacts'] > 0:
        report.append("Position Corrections:")
        report.append(f"  Mean:                   {stats['mean_correction']:.1f} bp")
        report.append(f"  Median:                 {stats['median_correction']:.1f} bp")
        report.append(f"  Maximum:                {stats['max_correction']} bp")

    report.append("")
    report.append("=" * 60)

    return "\n".join(report)
