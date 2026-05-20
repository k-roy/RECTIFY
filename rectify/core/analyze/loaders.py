"""
Data loaders for the RECTIFY analyze pipeline.

Functions for loading corrected-positions TSVs (with optional chunked
aggregation for files >0.5 GB), GTF/GFF/TSV annotations, and the
per-sample compact position index (`*_index.bed.gz`).
"""

from collections import defaultdict
from pathlib import Path
from typing import Optional

import pandas as pd


def load_corrected_positions(
    filepath: str,
    sample_column: str,
    normalize_chroms: bool = True,
    chrom_format: str = 'passthrough',
    chunk_size: int = 1_000_000,
    max_rows: int = None,
) -> pd.DataFrame:
    """Load corrected positions from TSV file with optional chunked loading.

    Args:
        filepath: Path to TSV file with corrected positions
        sample_column: Column name for sample identifier
        normalize_chroms: Whether to normalize chromosome names
        chrom_format: Target chromosome format ('passthrough', 'ucsc', 'ncbi', 'sgd'); default passthrough preserves input names
        chunk_size: Rows per chunk for large files (default: 1M)
        max_rows: Maximum rows to load (for testing)

    Returns:
        DataFrame with corrected positions
    """
    from ...utils.chromosome import normalize_dataframe_chromosomes
    import os

    # Check file size to determine loading strategy
    file_size = os.path.getsize(filepath)
    file_size_gb = file_size / (1024**3)

    if file_size_gb > 0.5:
        # Large file - use chunked loading with aggregation
        print(f"  Large file ({file_size_gb:.1f} GB) - using chunked loading...")
        return _load_large_file_chunked(
            filepath, sample_column, normalize_chroms, chrom_format, chunk_size, max_rows
        )

    # Standard loading for smaller files
    df = pd.read_csv(filepath, sep='\t', nrows=max_rows)

    # Handle different position column names
    position_col = None
    for col in ['corrected_position', 'corrected_3prime', 'position']:
        if col in df.columns:
            position_col = col
            break

    if position_col is None:
        raise ValueError("No position column found (tried: corrected_position, corrected_3prime, position)")

    # Standardize to 'corrected_position' if needed
    if position_col != 'corrected_position':
        df['corrected_position'] = df[position_col]
        print(f"  Using '{position_col}' as position column")

    # Auto-detect sample column if the specified name is not found
    if sample_column not in df.columns:
        _alt_sample_cols = ['sample', 'replicate', 'sample_id', 'sample_name', 'condition']
        _detected = next((c for c in _alt_sample_cols if c in df.columns), None)
        if _detected:
            print(f"  Auto-detected sample column '{_detected}' (requested '{sample_column}' not found)")
            sample_column = _detected
        else:
            raise ValueError(
                f"Sample column '{sample_column}' not found. "
                f"Tried fallbacks: {_alt_sample_cols}. "
                f"Available columns: {list(df.columns)}"
            )

    required_cols = ['chrom', 'strand', 'corrected_position', sample_column]
    missing = [c for c in required_cols if c not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Normalize chromosome names for consistent analysis
    if normalize_chroms:
        df = normalize_dataframe_chromosomes(df, 'chrom', chrom_format)
        print(f"  Normalized chromosome names to {chrom_format} format")

    # Drop columns not needed downstream to reduce memory usage
    _keep_cols = {'chrom', 'strand', 'corrected_position', sample_column}
    for _opt in ('fraction', 'alignment_start', 'alignment_end', 'count'):
        if _opt in df.columns:
            _keep_cols.add(_opt)
    _drop_cols = [c for c in df.columns if c not in _keep_cols]
    if _drop_cols:
        df = df.drop(columns=_drop_cols)

    return df


def _load_large_file_chunked(
    filepath: str,
    sample_column: str,
    normalize_chroms: bool,
    chrom_format: str,
    chunk_size: int,
    max_rows: int,
) -> pd.DataFrame:
    """Load large file in chunks, aggregating counts by position.

    Instead of loading all 200M+ rows, aggregate counts by:
    (chrom, strand, position, sample) -> count

    This reduces memory by ~100x for most datasets.
    Uses vectorized pandas operations instead of iterrows for speed.
    """
    from ...utils.chromosome import normalize_chromosome

    print(f"  Aggregating counts by position (chunk_size={chunk_size:,})...")

    # First pass: determine position column name
    header = pd.read_csv(filepath, sep='\t', nrows=0)
    position_col = None
    for col in ['corrected_position', 'corrected_3prime', 'position']:
        if col in header.columns:
            position_col = col
            break

    if position_col is None:
        raise ValueError("No position column found")

    # Check for fraction column (proportional assignment)
    has_fraction = 'fraction' in header.columns

    # Build usecols list — only load columns actually needed
    _usecols_set = {position_col, sample_column, 'chrom', 'strand'}
    if has_fraction:
        _usecols_set.add('fraction')
    for _optional_col in ('alignment_start', 'alignment_end'):
        if _optional_col in header.columns:
            _usecols_set.add(_optional_col)
    # Preserve a stable order: required cols first, then optional
    _usecols = [c for c in [position_col, sample_column, 'chrom', 'strand', 'fraction', 'alignment_start', 'alignment_end'] if c in _usecols_set]

    # Aggregated DataFrames from each chunk
    aggregated_chunks = []

    total_rows = 0
    for chunk_num, chunk in enumerate(pd.read_csv(filepath, sep='\t', chunksize=chunk_size, usecols=_usecols)):
        if max_rows and total_rows >= max_rows:
            break

        # Normalize chromosomes in chunk (vectorized)
        if normalize_chroms:
            chunk['chrom'] = chunk['chrom'].map(
                lambda x: normalize_chromosome(x, chrom_format)
            )

        # VECTORIZED AGGREGATION: Use groupby instead of iterrows
        group_cols = ['chrom', 'strand', position_col, sample_column]

        if has_fraction:
            # Sum fractions for each position/sample combination
            chunk_agg = chunk.groupby(group_cols)['fraction'].sum().reset_index()
            chunk_agg.columns = ['chrom', 'strand', 'corrected_position', sample_column, 'count']
        else:
            # Count occurrences
            chunk_agg = chunk.groupby(group_cols).size().reset_index(name='count')
            chunk_agg.columns = ['chrom', 'strand', 'corrected_position', sample_column, 'count']

        aggregated_chunks.append(chunk_agg)

        total_rows += len(chunk)
        if chunk_num % 10 == 0:
            print(f"    Processed {total_rows:,} rows...", flush=True)

    # Combine all chunk aggregations and re-aggregate
    print(f"  Combining {len(aggregated_chunks)} chunks...")
    combined = pd.concat(aggregated_chunks, ignore_index=True)

    # Final aggregation to merge duplicate keys across chunks
    group_cols = ['chrom', 'strand', 'corrected_position', sample_column]
    final = combined.groupby(group_cols)['count'].sum().reset_index()

    print(f"  Aggregated {total_rows:,} rows into {len(final):,} position/sample combinations")

    return final


def load_annotation(filepath: str, normalize_chroms: bool = True, chrom_format: str = 'passthrough') -> pd.DataFrame:
    """Load gene annotation file (GTF or TSV).

    Args:
        filepath: Path to annotation file
        normalize_chroms: Whether to normalize chromosome names
        chrom_format: Target chromosome format ('passthrough', 'ucsc', 'ncbi', 'sgd'); default passthrough preserves input names

    Returns:
        DataFrame with gene annotations
    """
    from ...utils.chromosome import normalize_dataframe_chromosomes

    filepath_str = str(filepath)
    if any(filepath_str.endswith(ext) for ext in ('.gtf', '.gff', '.gtf.gz', '.gff.gz', '.gff3', '.gff3.gz')):
        # Parse GTF/GFF
        df = _parse_gtf(filepath)
    else:
        # Assume TSV
        df = pd.read_csv(filepath, sep='\t')

    # Normalize chromosome names to match position data
    if normalize_chroms and 'chrom' in df.columns:
        df = normalize_dataframe_chromosomes(df, 'chrom', chrom_format)

    return df


def load_position_index(
    tsv_path: str,
    sample_id: str,
    normalize_chroms: bool = True,
    chrom_format: str = 'passthrough',
) -> Optional[pd.DataFrame]:
    """Load compact position index (_index.bed.gz) if it exists.

    Returns None if the index file is missing.
    """
    from ...utils.chromosome import normalize_dataframe_chromosomes

    base = Path(tsv_path)
    stem = base.name
    for _suffix in ('.tsv.gz', '.tsv'):
        if stem.endswith(_suffix):
            stem = stem[:-len(_suffix)]
            break
    index_path = base.parent / f"{stem}_index.bed.gz"
    if not index_path.exists():
        return None

    df = pd.read_csv(str(index_path), sep='\t', compression='gzip')
    df = df.rename(columns={'corrected_3prime': 'corrected_position'})
    df['sample'] = sample_id
    # Older indices (pre-2026-05-19) lack count_ag_rich; default to 0.0
    # so downstream code can rely on the column always being present.
    if 'count_ag_rich' not in df.columns:
        df['count_ag_rich'] = 0.0

    if normalize_chroms and 'chrom' in df.columns:
        df = normalize_dataframe_chromosomes(df, 'chrom', chrom_format)

    return df


def _parse_gtf(filepath: str) -> pd.DataFrame:
    """Parse GTF/GFF3 file to extract gene coordinates.

    Supports both GTF format (key "value") and GFF3 format (key=value).
    For GFF3, extracts:
        - ID -> gene_id
        - Name -> systematic gene name
        - gene -> common gene name (preferred for display)
    """
    genes = []

    import gzip as _gzip
    _open = _gzip.open if str(filepath).endswith('.gz') else open
    with _open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            if fields[2] != 'gene':
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1  # GFF is 1-based inclusive → 0-based half-open
            end = int(fields[4])        # GFF end is inclusive = 0-based exclusive (no change)
            strand = fields[6]

            # Parse attributes - handle both GTF and GFF3 formats
            attrs = {}
            attr_string = fields[8]

            # Detect format: GFF3 uses key=value, GTF uses key "value"
            if '=' in attr_string:
                # GFF3 format: ID=YAL069W;Name=YAL069W;gene=PAU8
                for attr in attr_string.split(';'):
                    attr = attr.strip()
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        # URL decode common escapes
                        value = value.replace('%20', ' ').replace('%3B', ';').replace('%2C', ',')
                        attrs[key] = value
            else:
                # GTF format: gene_id "YAL069W"; gene_name "PAU8"
                for attr in attr_string.split(';'):
                    attr = attr.strip()
                    if ' ' in attr:
                        key, value = attr.split(' ', 1)
                        attrs[key] = value.strip('"')

            # Extract gene_id - try GFF3 keys first, then GTF
            gene_id = attrs.get('ID') or attrs.get('gene_id') or f'{chrom}_{start}'

            # Extract gene_name - prefer common name, fall back to systematic name
            # GFF3: 'gene' has common name, 'Name' has systematic name
            # GTF: 'gene_name' has the name
            gene_name = attrs.get('gene') or attrs.get('gene_name') or attrs.get('Name') or gene_id

            genes.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': gene_id,
                'gene_name': gene_name,
            })

    return pd.DataFrame(genes)
