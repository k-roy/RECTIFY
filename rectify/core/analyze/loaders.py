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


_MANIFEST_HEADER_COLS = ['region_id', 'chrom', 'start', 'end', 'tsv_path', 'n_rows', 'sha256']


def _is_manifest(filepath: str) -> bool:
    """Return True if filepath is a corrected_reads.manifest.tsv (Commit B format)."""
    try:
        with open(filepath) as fh:
            first_line = fh.readline().rstrip('\n')
        return first_line.split('\t') == _MANIFEST_HEADER_COLS
    except OSError:
        return False


def _load_manifest_as_dataframe(
    manifest_path: str,
    sample_column: str,
    normalize_chroms: bool = True,
    chrom_format: str = 'passthrough',
    chunk_size: int = 1_000_000,
    max_rows: int = None,
) -> pd.DataFrame:
    """Load a corrected_reads.manifest.tsv by concatenating its region TSVs.

    Iterates region TSVs in manifest order and concatenates them into one
    DataFrame. Commit D will restructure analyze to stream partials; this
    in-memory concat is intentionally simple for now.
    """
    import os
    from ...core.bam.tsv_partition import load_manifest
    entries = load_manifest(Path(manifest_path))
    frames = []
    rows_loaded = 0
    for entry in entries:
        region_tsv = entry['tsv_path']  # absolute Path
        if not region_tsv.exists():
            raise FileNotFoundError(
                f"Manifest references missing region TSV: {region_tsv}"
            )
        df_part = pd.read_csv(str(region_tsv), sep='\t')
        if max_rows is not None:
            remaining = max_rows - rows_loaded
            if remaining <= 0:
                break
            df_part = df_part.iloc[:remaining]
        frames.append(df_part)
        rows_loaded += len(df_part)
        if max_rows is not None and rows_loaded >= max_rows:
            break
    if not frames:
        # Return an empty DataFrame with sensible columns
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    # Apply the same post-processing that load_corrected_positions does below.
    return _postprocess_positions_df(combined, sample_column, normalize_chroms, chrom_format)


def _postprocess_positions_df(
    df: pd.DataFrame,
    sample_column: str,
    normalize_chroms: bool,
    chrom_format: str,
) -> pd.DataFrame:
    """Apply position-column normalisation, sample-column auto-detect, and column
    trimming to df.  Mirrors the post-processing in the single-TSV code path."""
    from ...utils.chromosome import normalize_dataframe_chromosomes

    position_col = _select_position_column(df.columns)
    _validate_corrected_position_columns(df)

    if position_col is None:
        raise ValueError(
            "No position column found (tried: corrected_3prime, corrected_position, position)"
        )

    if position_col != 'corrected_position':
        df['corrected_position'] = df[position_col]
        print(f"  Using '{position_col}' as position column")

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


def _select_position_column(columns) -> Optional[str]:
    for col in ['corrected_3prime', 'corrected_position', 'position']:
        if col in columns:
            return col
    return None


def _validate_corrected_position_columns(df: pd.DataFrame) -> None:
    if 'corrected_3prime' not in df.columns or 'corrected_position' not in df.columns:
        return
    lhs = pd.to_numeric(df['corrected_3prime'], errors='coerce')
    rhs = pd.to_numeric(df['corrected_position'], errors='coerce')
    mismatch = lhs.notna() & rhs.notna() & (lhs != rhs)
    if mismatch.any():
        raise ValueError(
            "TSV contains both corrected_3prime and corrected_position with "
            f"{int(mismatch.sum())} differing rows; refusing to choose silently"
        )


def load_corrected_positions(
    filepath: str,
    sample_column: str,
    normalize_chroms: bool = True,
    chrom_format: str = 'passthrough',
    chunk_size: int = 1_000_000,
    max_rows: int = None,
) -> pd.DataFrame:
    """Load corrected positions from TSV file with optional chunked loading.

    Auto-detects whether ``filepath`` is a manifest
    (``corrected_reads.manifest.tsv``) or a single concatenated TSV.  Both
    forms produce the same DataFrame so downstream code is unchanged.

    Args:
        filepath: Path to TSV file with corrected positions OR a manifest.
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

    # Manifest auto-detect (Commit B, 2026-05-20)
    if _is_manifest(filepath):
        print(f"  Detected manifest format — loading per-region TSVs...")
        return _load_manifest_as_dataframe(
            filepath, sample_column, normalize_chroms, chrom_format, chunk_size, max_rows
        )

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
    position_col = _select_position_column(df.columns)
    _validate_corrected_position_columns(df)

    if position_col is None:
        raise ValueError("No position column found (tried: corrected_3prime, corrected_position, position)")

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
    position_col = _select_position_column(header.columns)

    if position_col is None:
        raise ValueError("No position column found")

    # Check for fraction column (proportional assignment)
    has_fraction = 'fraction' in header.columns

    # Build usecols list — only load columns actually needed
    _usecols_set = {position_col, sample_column, 'chrom', 'strand'}
    if 'corrected_3prime' in header.columns and 'corrected_position' in header.columns:
        _usecols_set.update(['corrected_3prime', 'corrected_position'])
    if has_fraction:
        _usecols_set.add('fraction')
    for _optional_col in ('alignment_start', 'alignment_end'):
        if _optional_col in header.columns:
            _usecols_set.add(_optional_col)
    # Preserve a stable order: required cols first, then optional
    _usecols = [c for c in [position_col, sample_column, 'chrom', 'strand', 'corrected_3prime', 'corrected_position', 'fraction', 'alignment_start', 'alignment_end'] if c in _usecols_set]

    # Aggregated DataFrames from each chunk
    aggregated_chunks = []

    total_rows = 0
    for chunk_num, chunk in enumerate(pd.read_csv(filepath, sep='\t', chunksize=chunk_size, usecols=_usecols)):
        if max_rows and total_rows >= max_rows:
            break

        # Normalize chromosomes in chunk (vectorized)
        _validate_corrected_position_columns(chunk)
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


def _parse_annotation_attributes(attr_string: str) -> dict:
    """Parse GTF/GFF3 attribute text into a dict."""
    from urllib.parse import unquote

    attrs = {}
    if '=' in attr_string:
        for attr in attr_string.split(';'):
            attr = attr.strip()
            if '=' not in attr:
                continue
            key, value = attr.split('=', 1)
            attrs[key] = unquote(value.strip().strip('"'))
    else:
        for attr in attr_string.split(';'):
            attr = attr.strip()
            if ' ' not in attr:
                continue
            key, value = attr.split(' ', 1)
            attrs[key] = value.strip().strip('"')
    return attrs


def _normalize_feature_type(feature: str) -> str:
    """Normalize common GTF/GFF feature spellings to downstream categories."""
    raw = str(feature or '').strip()
    key = raw.lower().replace('-', '_')
    feature_map = {
        'three_prime_utr': 'UTR3',
        '3utr': 'UTR3',
        '3_utr': 'UTR3',
        'utr3': 'UTR3',
        "3'utr": 'UTR3',
        'five_prime_utr': 'UTR5',
        '5utr': 'UTR5',
        '5_utr': 'UTR5',
        'utr5': 'UTR5',
        "5'utr": 'UTR5',
        'cds': 'CDS',
        'snorna_gene': 'snoRNA',
        'snorna': 'snoRNA',
        'snrna_gene': 'snRNA',
        'snrna': 'snRNA',
        'trna_gene': 'tRNA',
        'trna': 'tRNA',
        'rrna_gene': 'rRNA',
        'rrna': 'rRNA',
        'ncrna_gene': 'ncRNA',
        'ncrna': 'ncRNA',
        'cut': 'CUT',
        'sut': 'SUT',
        'xut': 'XUT',
    }
    return feature_map.get(key, raw)


def _parse_gtf(filepath: str) -> pd.DataFrame:
    """Parse GTF/GFF3 annotation, preserving gene and feature-level rows."""
    records = []
    import gzip as _gzip
    _open = _gzip.open if str(filepath).endswith('.gz') else open
    with _open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            try:
                start = int(fields[3]) - 1
                end = int(fields[4])
            except ValueError:
                continue
            if start < 0 or end <= start:
                continue
            strand = fields[6]
            raw_feature = fields[2]
            attrs = _parse_annotation_attributes(fields[8])
            parent_id = (attrs.get('Parent') or '').split(',')[0] or None
            record_id = attrs.get('ID') or attrs.get('gene_id') or attrs.get('transcript_id')
            feature_type = _normalize_feature_type(raw_feature)

            records.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': attrs.get('gene_id') or (record_id if feature_type == 'gene' else None),
                'gene_name': attrs.get('gene_name') or attrs.get('gene') or attrs.get('Name'),
                'feature_type': feature_type,
                'raw_feature': raw_feature,
                'source': fields[1],
                'transcript_id': attrs.get('transcript_id'),
                'parent_id': parent_id,
                '_record_id': record_id,
            })

    if not records:
        return pd.DataFrame(columns=[
            'chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name',
            'feature_type', 'raw_feature', 'source', 'transcript_id',
            'parent_id',
        ])

    id_to_gene = {}
    id_to_name = {}
    for rec in records:
        rec_id = rec.get('_record_id')
        if not rec_id:
            continue
        gene_id = rec.get('gene_id') or rec_id
        gene_name = rec.get('gene_name') or gene_id
        id_to_gene[rec_id] = gene_id
        id_to_name[rec_id] = gene_name

    for rec in records:
        parent_id = rec.get('parent_id')
        rec_id = rec.get('_record_id')
        if not rec.get('gene_id'):
            rec['gene_id'] = (
                id_to_gene.get(parent_id)
                or id_to_gene.get(rec_id)
                or parent_id
                or rec_id
                or f"{rec['chrom']}_{rec['start']}"
            )
        if not rec.get('gene_name'):
            rec['gene_name'] = (
                id_to_name.get(parent_id)
                or id_to_name.get(rec_id)
                or rec['gene_id']
            )

    df = pd.DataFrame(records)
    return df.drop(columns=['_record_id'])
