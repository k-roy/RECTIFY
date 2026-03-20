#!/usr/bin/env python3
"""
Chromosome name utilities for RECTIFY.

Provides standardized chromosome name handling across different reference
genome naming conventions (UCSC, NCBI RefSeq, Ensembl).

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
"""

from typing import Dict, Optional, List
import re

# Standard chromosome mappings for S. cerevisiae (yeast)
# Maps from various formats to a canonical format

# UCSC (chrI) to NCBI RefSeq (ref|NC_001133|)
YEAST_CHR_TO_NCBI = {
    'chrI': 'ref|NC_001133|',
    'chrII': 'ref|NC_001134|',
    'chrIII': 'ref|NC_001135|',
    'chrIV': 'ref|NC_001136|',
    'chrV': 'ref|NC_001137|',
    'chrVI': 'ref|NC_001138|',
    'chrVII': 'ref|NC_001139|',
    'chrVIII': 'ref|NC_001140|',
    'chrIX': 'ref|NC_001141|',
    'chrX': 'ref|NC_001142|',
    'chrXI': 'ref|NC_001143|',
    'chrXII': 'ref|NC_001144|',
    'chrXIII': 'ref|NC_001145|',
    'chrXIV': 'ref|NC_001146|',
    'chrXV': 'ref|NC_001147|',
    'chrXVI': 'ref|NC_001148|',
    'chrMito': 'ref|NC_001224|',
    'chrM': 'ref|NC_001224|',
    'chrMT': 'ref|NC_001224|',
    # Simple Roman numerals (Ensembl format without 'chr' prefix)
    'I': 'ref|NC_001133|',
    'II': 'ref|NC_001134|',
    'III': 'ref|NC_001135|',
    'IV': 'ref|NC_001136|',
    'V': 'ref|NC_001137|',
    'VI': 'ref|NC_001138|',
    'VII': 'ref|NC_001139|',
    'VIII': 'ref|NC_001140|',
    'IX': 'ref|NC_001141|',
    'X': 'ref|NC_001142|',
    'XI': 'ref|NC_001143|',
    'XII': 'ref|NC_001144|',
    'XIII': 'ref|NC_001145|',
    'XIV': 'ref|NC_001146|',
    'XV': 'ref|NC_001147|',
    'XVI': 'ref|NC_001148|',
    'Mito': 'ref|NC_001224|',
    'MT': 'ref|NC_001224|',
}

# NCBI RefSeq (ref|NC_001133|) to UCSC (chrI)
YEAST_NCBI_TO_CHR = {v: k for k, v in YEAST_CHR_TO_NCBI.items() if not k.startswith('chrM')}
YEAST_NCBI_TO_CHR['ref|NC_001224|'] = 'chrMito'

# SGD format (BK006935.2) to NCBI
YEAST_SGD_TO_NCBI = {
    'BK006935.2': 'ref|NC_001133|',  # chrI
    'BK006936.2': 'ref|NC_001134|',  # chrII
    'BK006937.2': 'ref|NC_001135|',  # chrIII
    'BK006938.2': 'ref|NC_001136|',  # chrIV
    'BK006939.2': 'ref|NC_001137|',  # chrV
    'BK006940.2': 'ref|NC_001138|',  # chrVI
    'BK006941.2': 'ref|NC_001139|',  # chrVII
    'BK006942.2': 'ref|NC_001140|',  # chrVIII
    'BK006943.2': 'ref|NC_001141|',  # chrIX
    'BK006944.2': 'ref|NC_001142|',  # chrX
    'BK006945.2': 'ref|NC_001143|',  # chrXI
    'BK006946.2': 'ref|NC_001144|',  # chrXII
    'BK006947.2': 'ref|NC_001145|',  # chrXIII
    'BK006948.2': 'ref|NC_001146|',  # chrXIV
    'BK006949.2': 'ref|NC_001147|',  # chrXV
    'BK006934.2': 'ref|NC_001148|',  # chrXVI
    'AJ011856.1': 'ref|NC_001224|',  # chrMito
}


def detect_chromosome_format(chrom_names: List[str]) -> str:
    """
    Detect the chromosome naming format from a list of chromosome names.

    Args:
        chrom_names: List of chromosome names from a BAM/annotation file

    Returns:
        Format string: 'ucsc' (chrI), 'ncbi' (ref|NC_001133|),
                      'sgd' (BK006935.2), or 'unknown'
    """
    sample = set(chrom_names[:20])  # Check first 20

    if any(c.startswith('ref|NC_') for c in sample):
        return 'ncbi'
    elif any(c.startswith('chr') for c in sample):
        return 'ucsc'
    elif any(c.startswith('BK00') for c in sample):
        return 'sgd'
    else:
        return 'unknown'


def normalize_chromosome(
    chrom: str,
    target_format: str = 'ucsc',
    organism: str = 'yeast',
) -> str:
    """
    Normalize a chromosome name to a target format.

    Args:
        chrom: Input chromosome name
        target_format: Target format ('ncbi', 'ucsc', or 'passthrough')
        organism: Organism ('yeast' supported)

    Returns:
        Normalized chromosome name
    """
    if target_format == 'passthrough':
        return chrom

    if organism != 'yeast':
        # For other organisms, just return as-is
        return chrom

    if target_format == 'ncbi':
        # Convert to NCBI RefSeq format
        if chrom.startswith('ref|NC_'):
            return chrom
        elif chrom in YEAST_CHR_TO_NCBI:
            return YEAST_CHR_TO_NCBI[chrom]
        elif chrom in YEAST_SGD_TO_NCBI:
            return YEAST_SGD_TO_NCBI[chrom]
        else:
            return chrom

    elif target_format == 'ucsc':
        # Convert to UCSC format
        if chrom.startswith('chr'):
            return chrom
        elif chrom in YEAST_NCBI_TO_CHR:
            return YEAST_NCBI_TO_CHR[chrom]
        else:
            # Try reverse SGD lookup
            for sgd, ncbi in YEAST_SGD_TO_NCBI.items():
                if chrom == sgd:
                    return YEAST_NCBI_TO_CHR.get(ncbi, chrom)
            return chrom

    return chrom


def build_chromosome_map(
    source_names: List[str],
    target_names: List[str],
) -> Dict[str, str]:
    """
    Build a chromosome name mapping between two reference genomes.

    Args:
        source_names: Chromosome names from source (e.g., BAM file)
        target_names: Chromosome names from target (e.g., annotation)

    Returns:
        Dict mapping source names to target names
    """
    chrom_map = {}

    source_format = detect_chromosome_format(source_names)
    target_format = detect_chromosome_format(target_names)

    if source_format == target_format:
        # Same format - identity mapping
        for name in source_names:
            if name in target_names:
                chrom_map[name] = name
    else:
        # Different formats - need to convert
        for source in source_names:
            # Normalize to NCBI as intermediate
            ncbi = normalize_chromosome(source, 'ncbi')

            # Then convert to target format
            if target_format == 'ucsc':
                target = normalize_chromosome(ncbi, 'ucsc')
            elif target_format == 'ncbi':
                target = ncbi
            else:
                target = source

            if target in target_names:
                chrom_map[source] = target

    return chrom_map


def normalize_dataframe_chromosomes(
    df,
    chrom_col: str = 'chrom',
    target_format: str = 'ucsc',
    organism: str = 'yeast',
    inplace: bool = False,
):
    """
    Normalize chromosome names in a pandas DataFrame.

    Args:
        df: DataFrame with chromosome column
        chrom_col: Column name for chromosome
        target_format: Target format ('ncbi', 'ucsc')
        organism: Organism ('yeast')
        inplace: Modify in place or return copy

    Returns:
        DataFrame with normalized chromosome names (if not inplace)
    """
    if not inplace:
        df = df.copy()

    df[chrom_col] = df[chrom_col].apply(
        lambda x: normalize_chromosome(x, target_format, organism)
    )

    if not inplace:
        return df
