#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Splice-Aware 5' End Correction for RECTIFY

This module corrects the 5' end position of reads based on splice junction information.
Many introns are close to the TSS, and aligners may place the 5' end within an intron
rather than at the true transcription start site.

Problem:
    When the first exon is very short, the read's leftmost aligned position may be:
    - Within the first exon (correct 5' end)
    - Within an intron (needs correction to exon boundary)

Solution:
    1. Check if the read's 5' end falls within an annotated intron
    2. If so, shift the 5' end to the nearest upstream exon boundary
    3. Use consensus junctions to validate/refine the correction

Coordinate System:
    - 0-based, half-open coordinates (consistent with pysam/BED)
    - Plus strand: 5' end = leftmost position
    - Minus strand: 5' end = rightmost position

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-24
"""

from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
import pysam

try:
    from intervaltree import IntervalTree
    HAS_INTERVALTREE = True
except ImportError:
    HAS_INTERVALTREE = False


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class FivePrimeCorrection:
    """Result of 5' end correction."""
    five_prime_raw: int             # Original alignment position
    five_prime_corrected: int       # After splice-aware correction
    first_exon_start: Optional[int] # Start of first exon (may be None)
    starts_in_intron: bool          # True if correction was applied
    correction_bp: int              # Bases shifted (0 if no correction)
    correction_reason: str          # Explanation of correction


# =============================================================================
# Helper Functions
# =============================================================================

def get_read_5prime_position(read: pysam.AlignedSegment, strand: Optional[str] = None) -> int:
    """
    Get the 5' end genomic position from a read.

    The 5' end is the transcription start site (TSS) end of the RNA molecule.

    Args:
        read: pysam AlignedSegment
        strand: Optional strand override ('+' or '-')

    Returns:
        5' end position (0-based)
    """
    if strand is None:
        strand = '-' if read.is_reverse else '+'

    if strand == '+':
        return read.reference_start
    else:
        return read.reference_end - 1


def build_intron_interval_tree(
    annotation_df: 'pd.DataFrame',
    feature_types: Optional[List[str]] = None,
) -> Dict[Tuple[str, str], 'IntervalTree']:
    """
    Build IntervalTree index for introns from annotation.

    This identifies intron regions (gaps between CDS/exon features).

    Args:
        annotation_df: DataFrame with gene annotations
            Required columns: chrom, start, end, strand
            Optional columns: gene_id, feature_type
        feature_types: Feature types to use for exon boundaries (default: ['CDS', 'exon'])

    Returns:
        Dict mapping (chrom, strand) to IntervalTree of introns
    """
    if not HAS_INTERVALTREE:
        raise ImportError("intervaltree is required. Install with: pip install intervaltree")

    import pandas as pd
    from collections import defaultdict

    if feature_types is None:
        feature_types = ['CDS', 'exon', 'mRNA']

    # Filter to relevant feature types
    if 'feature_type' in annotation_df.columns:
        df = annotation_df[annotation_df['feature_type'].isin(feature_types)].copy()
    else:
        df = annotation_df.copy()

    # Build exon intervals per gene
    intron_trees = defaultdict(IntervalTree)

    # Group by gene and strand
    if 'gene_id' in df.columns:
        group_cols = ['chrom', 'strand', 'gene_id']
    else:
        group_cols = ['chrom', 'strand']

    for group_key, group_df in df.groupby(group_cols):
        if len(group_cols) == 3:
            chrom, strand, gene_id = group_key
        else:
            chrom, strand = group_key
            gene_id = "unknown"

        # Sort exons by position
        exons = group_df.sort_values('start')

        # Find introns (gaps between consecutive exons)
        exon_list = list(zip(exons['start'], exons['end']))

        for i in range(len(exon_list) - 1):
            _, exon1_end = exon_list[i]
            exon2_start, _ = exon_list[i + 1]

            # Intron is the gap between exons
            if exon2_start > exon1_end:
                intron_start = exon1_end
                intron_end = exon2_start

                key = (chrom, strand)
                intron_trees[key][intron_start:intron_end] = {
                    'gene_id': gene_id,
                    'upstream_exon_end': exon1_end,
                    'downstream_exon_start': exon2_start,
                    'intron_index': i,
                }

    return dict(intron_trees)


def find_overlapping_introns(
    chrom: str,
    position: int,
    strand: str,
    intron_trees: Dict[Tuple[str, str], 'IntervalTree'],
) -> List[Dict]:
    """
    Find introns that overlap a position.

    Args:
        chrom: Chromosome name
        position: Genomic position (0-based)
        strand: Strand ('+' or '-')
        intron_trees: Dict from build_intron_interval_tree()

    Returns:
        List of intron info dicts
    """
    key = (chrom, strand)

    if key not in intron_trees:
        return []

    tree = intron_trees[key]
    overlaps = tree[position]

    results = []
    for interval in overlaps:
        intron_info = interval.data.copy()
        intron_info['intron_start'] = interval.begin
        intron_info['intron_end'] = interval.end
        results.append(intron_info)

    return results


# =============================================================================
# Core Correction Functions
# =============================================================================

def correct_5prime_for_splicing(
    read: pysam.AlignedSegment,
    intron_trees: Optional[Dict[Tuple[str, str], 'IntervalTree']] = None,
    read_junctions: Optional[List[Tuple[int, int]]] = None,
    strand: Optional[str] = None,
) -> FivePrimeCorrection:
    """
    Correct 5' end position based on splice junction information.

    If the read starts within an intron (near TSS), the true 5' end
    should be at the upstream exon boundary, not the raw alignment start.

    Args:
        read: pysam AlignedSegment
        intron_trees: Dict from build_intron_interval_tree() (optional)
        read_junctions: List of (start, end) tuples for junctions in the read
        strand: Strand override, or None to infer from read

    Returns:
        FivePrimeCorrection with correction details
    """
    if strand is None:
        strand = '-' if read.is_reverse else '+'

    chrom = read.reference_name
    five_prime_raw = get_read_5prime_position(read, strand)

    # Default: no correction
    result = FivePrimeCorrection(
        five_prime_raw=five_prime_raw,
        five_prime_corrected=five_prime_raw,
        first_exon_start=None,
        starts_in_intron=False,
        correction_bp=0,
        correction_reason="",
    )

    # Strategy 1: Use read junctions to infer 5' correction
    # If the first junction in the read is very close to the 5' end,
    # the 5' end might be in an intron
    if read_junctions:
        if strand == '+':
            # For plus strand, 5' is leftmost
            # Check if first junction starts very close to 5' end
            first_junction_start = min(j[0] for j in read_junctions)
            first_junction_end = min(j[1] for j in read_junctions if j[0] == first_junction_start)

            distance_to_first_junction = first_junction_start - five_prime_raw

            # If 5' end is within the first junction, correct to upstream exon
            if five_prime_raw >= first_junction_start and five_prime_raw < first_junction_end:
                result.starts_in_intron = True
                result.five_prime_corrected = first_junction_start
                result.first_exon_start = first_junction_start
                result.correction_bp = five_prime_raw - first_junction_start
                result.correction_reason = "5' end within read junction - shifted to exon boundary"
                return result

        else:
            # For minus strand, 5' is rightmost
            # Check if last junction ends very close to 5' end
            last_junction_end = max(j[1] for j in read_junctions)
            last_junction_start = max(j[0] for j in read_junctions if j[1] == last_junction_end)

            # If 5' end is within the last junction, correct to downstream exon
            if five_prime_raw >= last_junction_start and five_prime_raw < last_junction_end:
                result.starts_in_intron = True
                result.five_prime_corrected = last_junction_end - 1
                result.first_exon_start = last_junction_end
                result.correction_bp = last_junction_end - 1 - five_prime_raw
                result.correction_reason = "5' end within read junction - shifted to exon boundary"
                return result

    # Strategy 2: Use annotation introns
    if intron_trees:
        overlapping_introns = find_overlapping_introns(chrom, five_prime_raw, strand, intron_trees)

        if overlapping_introns:
            # Find the most relevant intron (closest to 5' end)
            if strand == '+':
                # For plus strand, we want the intron that starts closest to the 5' end
                # and shift to the downstream exon (higher coordinate)
                best_intron = min(overlapping_introns, key=lambda x: x['intron_start'])
                corrected_pos = best_intron['downstream_exon_start']
            else:
                # For minus strand, we want the intron that ends closest to the 5' end
                # and shift to the upstream exon (lower coordinate)
                best_intron = max(overlapping_introns, key=lambda x: x['intron_end'])
                corrected_pos = best_intron['upstream_exon_end'] - 1

            result.starts_in_intron = True
            result.five_prime_corrected = corrected_pos
            result.first_exon_start = corrected_pos if strand == '+' else corrected_pos + 1
            result.correction_bp = abs(corrected_pos - five_prime_raw)
            result.correction_reason = f"5' end in annotated intron ({best_intron.get('gene_id', 'unknown')}) - shifted to exon boundary"

    return result


def correct_5prime_batch(
    reads_with_junctions: List[Tuple[pysam.AlignedSegment, List[Tuple[int, int]]]],
    intron_trees: Optional[Dict[Tuple[str, str], 'IntervalTree']] = None,
) -> List[FivePrimeCorrection]:
    """
    Correct 5' end positions for a batch of reads.

    Args:
        reads_with_junctions: List of (read, junctions) tuples
        intron_trees: Dict from build_intron_interval_tree()

    Returns:
        List of FivePrimeCorrection objects
    """
    corrections = []

    for read, junctions in reads_with_junctions:
        correction = correct_5prime_for_splicing(
            read,
            intron_trees=intron_trees,
            read_junctions=junctions,
        )
        corrections.append(correction)

    return corrections


# =============================================================================
# Validation
# =============================================================================

def validate_5prime_correction():
    """
    Print validation examples for 5' end correction.
    """
    print("=" * 60)
    print("SPLICE-AWARE 5' END CORRECTION VALIDATION")
    print("=" * 60)

    # Create mock read for testing
    class MockRead:
        def __init__(self, start, end, is_reverse, reference_name='chrI'):
            self.reference_start = start
            self.reference_end = end
            self.is_reverse = is_reverse
            self.reference_name = reference_name

    # Test case 1: Plus strand, 5' end in junction
    print("\nTest 1: Plus strand read with 5' end in junction")
    read1 = MockRead(start=1050, end=2000, is_reverse=False)
    # Junction from 1000-1100, so 5' end (1050) is in the junction
    junctions1 = [(1000, 1100)]
    result1 = correct_5prime_for_splicing(read1, read_junctions=junctions1)
    print(f"  Raw 5' end: {result1.five_prime_raw}")
    print(f"  Corrected 5' end: {result1.five_prime_corrected}")
    print(f"  Starts in intron: {result1.starts_in_intron}")
    print(f"  Correction: {result1.correction_bp} bp")
    print(f"  Reason: {result1.correction_reason}")
    if result1.starts_in_intron and result1.five_prime_corrected == 1000:
        print("  PASSED")
    else:
        print("  FAILED")

    # Test case 2: Plus strand, 5' end NOT in junction
    print("\nTest 2: Plus strand read with 5' end NOT in junction")
    read2 = MockRead(start=900, end=2000, is_reverse=False)
    # Junction from 1000-1100, so 5' end (900) is before the junction
    junctions2 = [(1000, 1100)]
    result2 = correct_5prime_for_splicing(read2, read_junctions=junctions2)
    print(f"  Raw 5' end: {result2.five_prime_raw}")
    print(f"  Corrected 5' end: {result2.five_prime_corrected}")
    print(f"  Starts in intron: {result2.starts_in_intron}")
    if not result2.starts_in_intron and result2.five_prime_corrected == result2.five_prime_raw:
        print("  PASSED")
    else:
        print("  FAILED")

    # Test case 3: Minus strand, 5' end in junction
    print("\nTest 3: Minus strand read with 5' end in junction")
    read3 = MockRead(start=1000, end=2050, is_reverse=True)
    # Junction from 2000-2100, so 5' end (2050-1=2049) is in the junction
    junctions3 = [(2000, 2100)]
    result3 = correct_5prime_for_splicing(read3, read_junctions=junctions3)
    print(f"  Raw 5' end: {result3.five_prime_raw}")
    print(f"  Corrected 5' end: {result3.five_prime_corrected}")
    print(f"  Starts in intron: {result3.starts_in_intron}")
    print(f"  Correction: {result3.correction_bp} bp")
    if result3.starts_in_intron and result3.five_prime_corrected == 2099:
        print("  PASSED")
    else:
        print("  FAILED")

    print("\n" + "=" * 60)
    print("VALIDATION COMPLETE")
    print("=" * 60)


if __name__ == '__main__':
    validate_5prime_correction()
