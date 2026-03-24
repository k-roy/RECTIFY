#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gene Attribution Module for RECTIFY

This module attributes 3' end positions to genes based on the coverage of read bodies.
For each 3' end (cluster), we track what fraction of reads overlap which CDS/ncRNA features.

## Coordinate System Documentation

RECTIFY uses **0-based, half-open coordinates** consistent with pysam and BED format:
- `start` is inclusive (first base is at position `start`)
- `end` is exclusive (last base is at position `end - 1`)

### 5' vs 3' End Positions (RNA perspective)

For an RNA molecule transcribed from plus strand DNA:
```
    5' -----> 3'  (RNA, transcript direction)

    Genomic:  100       110       120       130       140
              |---------|---------|---------|---------|
    Read:         [=================]
              start=105          end=135 (exclusive)

    5' end = reference_start = 105 (leftmost aligned base)
    3' end = reference_end - 1 = 134 (rightmost aligned base)
```

For an RNA molecule transcribed from minus strand DNA:
```
    3' <----- 5'  (RNA, transcript direction is RIGHT to LEFT)

    Genomic:  100       110       120       130       140
              |---------|---------|---------|---------|
    Read:         [=================]
              start=105          end=135 (exclusive)

    5' end = reference_end - 1 = 134 (rightmost = 5' for minus strand)
    3' end = reference_start = 105 (leftmost = 3' for minus strand)
```

### Summary Table

| Strand | 5' end position | 3' end position |
|--------|-----------------|-----------------|
| +      | reference_start | reference_end-1 |
| -      | reference_end-1 | reference_start |

### Validation Examples (S. cerevisiae)

Example 1: YAL038W (CDC19, plus strand)
- Gene: chrI:71,787-73,289 (GFF 1-based: 71788-73289)
- Expected 3' ends: Near 73,289
- Expected 5' ends: Near 71,787

Example 2: YAL001C (TFC3, minus strand)
- Gene: chrI:151,097-151,166 (GFF 1-based)
- Expected 3' ends: Near 151,097 (lower coord = 3' for minus)
- Expected 5' ends: Near 151,166 (higher coord = 5' for minus)

Author: Kevin R. Roy
Date: 2026-03-24
"""

from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import pysam


def get_read_5prime_position(read: pysam.AlignedSegment, strand: Optional[str] = None) -> int:
    """
    Get the 5' end genomic position from a read.

    The 5' end is the transcription start site (TSS) end of the RNA molecule.
    This is the OPPOSITE end from the 3' end (CPA site).

    Args:
        read: pysam AlignedSegment object
        strand: Optional strand override ('+' or '-'). If None, uses read.is_reverse.

    Returns:
        5' end position (0-based genomic coordinate)

    Coordinate Details:
        - pysam reference_start: 0-based leftmost aligned position (inclusive)
        - pysam reference_end: 0-based rightmost position + 1 (exclusive)

        Plus strand (+): RNA 5'→3' matches genomic left→right
            5' end = reference_start (leftmost = 5' end)

        Minus strand (-): RNA 5'→3' matches genomic right→left
            5' end = reference_end - 1 (rightmost = 5' end)

    Example (plus strand):
        Genomic:  [100]-----[109]  (reference_start=100, reference_end=110)
        RNA:      5'=========3'
        5' end = 100

    Example (minus strand):
        Genomic:  [100]-----[109]  (reference_start=100, reference_end=110)
        RNA:      3'=========5'
        5' end = 109
    """
    if strand is None:
        strand = '-' if read.is_reverse else '+'

    if strand == '+':
        return read.reference_start  # Leftmost position (0-based)
    else:
        return read.reference_end - 1  # Rightmost position (0-based inclusive)


def get_read_3prime_position(read: pysam.AlignedSegment, strand: Optional[str] = None) -> int:
    """
    Get the 3' end genomic position from a read.

    The 3' end is the cleavage and polyadenylation (CPA) site of the RNA molecule.

    Args:
        read: pysam AlignedSegment object
        strand: Optional strand override ('+' or '-'). If None, uses read.is_reverse.

    Returns:
        3' end position (0-based genomic coordinate)

    Coordinate Details:
        Plus strand (+): 3' end = reference_end - 1 (rightmost)
        Minus strand (-): 3' end = reference_start (leftmost)
    """
    if strand is None:
        strand = '-' if read.is_reverse else '+'

    if strand == '+':
        return read.reference_end - 1  # Rightmost position
    else:
        return read.reference_start  # Leftmost position


def get_read_body_interval(read: pysam.AlignedSegment) -> Tuple[int, int]:
    """
    Get the aligned interval (read body) from a read.

    Args:
        read: pysam AlignedSegment object

    Returns:
        Tuple of (start, end) in 0-based half-open coordinates
        The interval [start, end) covers all aligned bases.
    """
    return (read.reference_start, read.reference_end)


def build_cds_interval_tree(
    annotation_df: pd.DataFrame,
    feature_types: Optional[List[str]] = None,
) -> Dict[Tuple[str, str], 'IntervalTree']:
    """
    Build IntervalTree index for CDS and ncRNA features.

    Args:
        annotation_df: DataFrame with gene annotations
            Required columns: chrom, start, end, strand
            Optional columns: gene_id, gene_name, feature_type
        feature_types: List of feature types to include (default: ['CDS', 'gene'])

    Returns:
        Dict mapping (chrom, strand) to IntervalTree
        Each interval stores: {'gene_id': str, 'gene_name': str, 'feature_type': str}

    Note:
        The IntervalTree uses half-open intervals [start, end).
        For 1-based GFF annotations, we convert: start_0based = start_gff - 1
    """
    from intervaltree import IntervalTree

    if feature_types is None:
        feature_types = ['CDS', 'gene', 'mRNA', 'ncRNA', 'tRNA', 'snoRNA', 'snRNA']

    trees = defaultdict(IntervalTree)

    for _, row in annotation_df.iterrows():
        chrom = row['chrom']
        strand = row.get('strand', '+')

        # Handle both 0-based and 1-based annotations
        # If 'start' looks like GFF (1-based), convert to 0-based
        start = int(row['start'])
        end = int(row['end'])

        # GFF is 1-based inclusive, we need 0-based half-open
        # If start is 1-based, convert: start_0 = start - 1, end stays same (becomes exclusive)
        # Most annotation_df from RECTIFY's load_annotation are already 0-based

        gene_id = row.get('gene_id', '')
        gene_name = row.get('gene_name', gene_id)
        feature_type = row.get('feature_type', 'gene')

        if feature_type not in feature_types:
            continue

        key = (chrom, strand)

        # IntervalTree uses half-open [start, end)
        # Add 1 to end for inclusive end coordinate if needed
        if end <= start:
            continue  # Skip invalid intervals

        trees[key][start:end] = {
            'gene_id': gene_id,
            'gene_name': gene_name,
            'feature_type': feature_type,
            'start': start,
            'end': end,
        }

    return dict(trees)


def find_overlapping_genes(
    chrom: str,
    start: int,
    end: int,
    strand: str,
    interval_trees: Dict[Tuple[str, str], 'IntervalTree'],
) -> List[Dict]:
    """
    Find all genes that overlap with a genomic interval.

    Args:
        chrom: Chromosome name
        start: Start position (0-based)
        end: End position (0-based exclusive)
        strand: Strand ('+' or '-')
        interval_trees: Dict from build_cds_interval_tree()

    Returns:
        List of gene info dicts with keys: gene_id, gene_name, feature_type, overlap_bp
    """
    key = (chrom, strand)

    if key not in interval_trees:
        return []

    tree = interval_trees[key]
    overlaps = tree[start:end]

    results = []
    for interval in overlaps:
        gene_info = interval.data.copy()

        # Calculate overlap size
        overlap_start = max(start, interval.begin)
        overlap_end = min(end, interval.end)
        gene_info['overlap_bp'] = max(0, overlap_end - overlap_start)

        results.append(gene_info)

    return results


def compute_read_gene_attribution(
    read: pysam.AlignedSegment,
    interval_trees: Dict[Tuple[str, str], 'IntervalTree'],
    chrom: Optional[str] = None,
) -> List[str]:
    """
    Determine which genes a read body overlaps.

    Args:
        read: pysam AlignedSegment
        interval_trees: Dict from build_cds_interval_tree()
        chrom: Optional chromosome name override

    Returns:
        List of gene_id strings that the read body overlaps.
        Empty list if no overlap (intergenic).

    Example:
        If read spans genes A and B: ['YBR296C', 'YBR297W']
        If intergenic: []
    """
    if chrom is None:
        chrom = read.reference_name

    strand = '-' if read.is_reverse else '+'
    start, end = get_read_body_interval(read)

    overlaps = find_overlapping_genes(chrom, start, end, strand, interval_trees)

    # Return gene IDs sorted by overlap size (largest first)
    overlaps.sort(key=lambda x: x['overlap_bp'], reverse=True)

    return [g['gene_id'] for g in overlaps if g['gene_id']]


def aggregate_attributions_for_3prime_end(
    attributions: List[Tuple[str, ...]],
    min_fraction: float = 0.01,
) -> str:
    """
    Aggregate gene attributions from multiple reads sharing a 3' end.

    This computes the fraction of reads overlapping each unique combination
    of genes.

    Args:
        attributions: List of tuples, each containing gene IDs for one read.
            Example: [('geneA',), ('geneA', 'geneB'), ('geneA',), ()]
        min_fraction: Minimum fraction to report (default 1%)

    Returns:
        Attribution string in format: "geneA,50%;geneA,geneB,25%;none,25%"

        The percentages always sum to 100% (within rounding).
        "none" is used for intergenic reads (no gene overlap).

    Example:
        Input: [('YBR296C',), ('YBR296C',), ('YBR296C', 'YBR297W'), ()]
        Output: "YBR296C,50%;YBR296C,YBR297W,25%;none,25%"
    """
    if not attributions:
        return "none,100%"

    # Count occurrences of each unique gene combination
    counts = defaultdict(int)
    total = len(attributions)

    for genes in attributions:
        if not genes:
            key = ('none',)
        else:
            # Sort gene IDs for consistent ordering
            key = tuple(sorted(genes))
        counts[key] += 1

    # Convert to fractions and format
    parts = []
    for genes, count in sorted(counts.items(), key=lambda x: -x[1]):
        fraction = count / total
        if fraction >= min_fraction:
            pct = round(fraction * 100)
            gene_str = ','.join(genes)
            parts.append(f"{gene_str},{pct}%")

    # Ensure we have at least one entry
    if not parts:
        parts = ["none,100%"]

    return ';'.join(parts)


def compute_cluster_attribution(
    cluster_reads: List[pysam.AlignedSegment],
    interval_trees: Dict[Tuple[str, str], 'IntervalTree'],
    chrom: Optional[str] = None,
) -> Dict[str, any]:
    """
    Compute gene attribution for a cluster of reads sharing a 3' end region.

    Args:
        cluster_reads: List of reads in the cluster
        interval_trees: Dict from build_cds_interval_tree()
        chrom: Optional chromosome name override

    Returns:
        Dict with:
            'attribution_string': Formatted attribution (e.g., "geneA,75%;geneB,25%")
            'n_reads': Total number of reads
            'n_intergenic': Number of intergenic reads
            'primary_gene': Most common gene (or 'none')
            'primary_fraction': Fraction of reads with primary gene
    """
    attributions = []

    for read in cluster_reads:
        genes = compute_read_gene_attribution(read, interval_trees, chrom)
        attributions.append(tuple(genes))

    attribution_string = aggregate_attributions_for_3prime_end(attributions)

    # Parse to find primary gene
    n_intergenic = sum(1 for a in attributions if not a)

    # Find most common single-gene attribution
    gene_counts = defaultdict(int)
    for genes in attributions:
        if genes:
            for g in genes:
                gene_counts[g] += 1

    if gene_counts:
        primary_gene = max(gene_counts.keys(), key=lambda g: gene_counts[g])
        primary_fraction = gene_counts[primary_gene] / len(attributions)
    else:
        primary_gene = 'none'
        primary_fraction = 1.0

    return {
        'attribution_string': attribution_string,
        'n_reads': len(cluster_reads),
        'n_intergenic': n_intergenic,
        'primary_gene': primary_gene,
        'primary_fraction': primary_fraction,
    }


# ============================================================================
# Validation Functions
# ============================================================================

def validate_strand_handling():
    """
    Print validation examples for strand handling.

    Run this to verify coordinate logic is correct.
    """
    print("=" * 60)
    print("STRAND HANDLING VALIDATION")
    print("=" * 60)

    # Create mock read objects for testing
    class MockRead:
        def __init__(self, start, end, is_reverse):
            self.reference_start = start
            self.reference_end = end
            self.is_reverse = is_reverse
            self.reference_name = 'chrI'

    # Test case 1: Plus strand read
    read_plus = MockRead(start=100, end=200, is_reverse=False)
    print("\nPlus strand read (forward):")
    print(f"  reference_start = {read_plus.reference_start}")
    print(f"  reference_end = {read_plus.reference_end}")
    print(f"  5' end = {get_read_5prime_position(read_plus)} (should be 100)")
    print(f"  3' end = {get_read_3prime_position(read_plus)} (should be 199)")

    # Test case 2: Minus strand read
    read_minus = MockRead(start=100, end=200, is_reverse=True)
    print("\nMinus strand read (reverse):")
    print(f"  reference_start = {read_minus.reference_start}")
    print(f"  reference_end = {read_minus.reference_end}")
    print(f"  5' end = {get_read_5prime_position(read_minus)} (should be 199)")
    print(f"  3' end = {get_read_3prime_position(read_minus)} (should be 100)")

    # Verify correctness
    assert get_read_5prime_position(read_plus) == 100, "Plus strand 5' end wrong"
    assert get_read_3prime_position(read_plus) == 199, "Plus strand 3' end wrong"
    assert get_read_5prime_position(read_minus) == 199, "Minus strand 5' end wrong"
    assert get_read_3prime_position(read_minus) == 100, "Minus strand 3' end wrong"

    print("\n✓ All strand handling tests passed!")


def validate_attribution_format():
    """
    Print validation examples for attribution string format.
    """
    print("\n" + "=" * 60)
    print("ATTRIBUTION FORMAT VALIDATION")
    print("=" * 60)

    # Test case 1: Single gene
    attrs1 = [('YBR296C',), ('YBR296C',), ('YBR296C',), ('YBR296C',)]
    result1 = aggregate_attributions_for_3prime_end(attrs1)
    print(f"\n4 reads all in YBR296C:")
    print(f"  Input: {attrs1}")
    print(f"  Output: {result1}")
    assert result1 == "YBR296C,100%", f"Expected 'YBR296C,100%', got '{result1}'"

    # Test case 2: Mixed genes
    attrs2 = [('YBR296C',), ('YBR296C',), ('YBR297W',), ()]
    result2 = aggregate_attributions_for_3prime_end(attrs2)
    print(f"\n2 in YBR296C, 1 in YBR297W, 1 intergenic:")
    print(f"  Input: {attrs2}")
    print(f"  Output: {result2}")
    # Should be "YBR296C,50%;YBR297W,25%;none,25%"
    assert "YBR296C,50%" in result2, f"Expected YBR296C,50% in '{result2}'"
    assert "none,25%" in result2, f"Expected none,25% in '{result2}'"

    # Test case 3: Spanning multiple genes
    attrs3 = [('YBR296C', 'YBR297W'), ('YBR296C', 'YBR297W'), ('YBR296C',), ('YBR296C',)]
    result3 = aggregate_attributions_for_3prime_end(attrs3)
    print(f"\n2 spanning both genes, 2 in YBR296C only:")
    print(f"  Input: {attrs3}")
    print(f"  Output: {result3}")

    # Test case 4: All intergenic
    attrs4 = [(), (), (), ()]
    result4 = aggregate_attributions_for_3prime_end(attrs4)
    print(f"\n4 intergenic reads:")
    print(f"  Input: {attrs4}")
    print(f"  Output: {result4}")
    assert result4 == "none,100%", f"Expected 'none,100%', got '{result4}'"

    print("\n✓ All attribution format tests passed!")


if __name__ == '__main__':
    print("Running gene_attribution.py validation...\n")
    validate_strand_handling()
    validate_attribution_format()
    print("\n" + "=" * 60)
    print("ALL VALIDATION TESTS PASSED")
    print("=" * 60)
