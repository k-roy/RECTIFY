"""
Generic locus loading and signal index construction for metagene analysis.

Provides project-independent functions to build locus lists and PositionIndex
objects from common genomic data sources:

    Locus sources (→ List[Dict]):
        loci_from_tsv          — any TSV/CSV with explicit center column
        loci_from_bed          — BED-format annotation files
        loci_from_gff          — GFF3 feature annotations (genes, TSSs, CPAs, …)
        loci_from_motif_scan   — regex/IUPAC motif scan over genome sequences
        loci_from_pickle       — backward-compat for dict-based PKL caches

    Signal sources (→ PositionIndex):
        position_index_from_tsv     — position TSVs (RECTIFY, NET-seq, …)
        position_index_from_bigwig  — BigWig files (NET-seq, ChIP-seq, PAR-CLIP, …)

All locus dicts produced have at minimum:
    {'chrom': str, 'strand': str, 'center': int}   # center is always 0-based

Center semantics (in 5'→3' transcription direction):
    'start'     TSS / 5' end (gene start for + strand, gene end for - strand)
    'end'       TES / 3' end (gene end for + strand, gene start for - strand)
    'tss'       alias for 'start'
    'tes'       alias for 'end'
    'midpoint'  geometric midpoint

IUPAC nucleotide codes supported by loci_from_motif_scan:
    N=[ACGT] R=[AG] Y=[CT] S=[GC] W=[AT] K=[GT] M=[AC]
    B=[CGT] D=[AGT] H=[ACT] V=[ACG]

Author: Kevin R. Roy
"""

import re
import pickle
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .metagene import PositionIndex


# ---------------------------------------------------------------------------
# IUPAC helpers
# ---------------------------------------------------------------------------

_IUPAC_MAP = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]',
    'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]',
    'N': '[ACGT]',
}

_COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVNacgtryswkmbdhvn',
                             'TGCAYRSWMKVHDBNtgcayrswmkvhdbn')


def iupac_to_regex(motif: str) -> str:
    """Convert an IUPAC nucleotide string to a regex pattern.

    Supports: A C G T U R Y S W K M B D H V N (case-insensitive).
    Non-IUPAC characters (e.g. quantifiers like {5,}) are passed through.

    Example:
        >>> iupac_to_regex('TATAAA')
        'TATAAA'
        >>> iupac_to_regex('WTATAWW')
        '[AT]TATA[AT][AT]'
    """
    result = []
    i = 0
    motif_upper = motif.upper()
    while i < len(motif_upper):
        ch = motif_upper[i]
        if ch in _IUPAC_MAP:
            result.append(_IUPAC_MAP[ch])
        else:
            # Pass through regex meta-characters and quantifiers unchanged
            result.append(ch)
        i += 1
    return ''.join(result)


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _center_from_match(match, chr_len: int, strand: str, center: str) -> int:
    """Compute 0-based center position from a regex match object.

    For plus strand matches the position is the genomic position directly.
    For minus strand matches (found by scanning the reverse complement), the
    genomic position must be flipped back: pos = chr_len - 1 - rc_pos.

    Center semantics (in 5'→3' transcription direction):
        'start'/'tss'  → 5' end of match
        'end'/'tes'    → 3' end of match (last base, inclusive, 0-based)
        'midpoint'     → geometric midpoint (integer division)
    """
    center = center.lower()
    if center in ('start', 'tss'):
        if strand == '+':
            return match.start()
        else:
            return chr_len - 1 - match.start()
    elif center in ('end', 'tes'):
        if strand == '+':
            return match.end() - 1
        else:
            return chr_len - match.end()
    elif center == 'midpoint':
        mid_rc = (match.start() + match.end() - 1) // 2
        if strand == '+':
            return mid_rc
        else:
            return chr_len - 1 - mid_rc
    else:
        raise ValueError(
            f"Unknown center='{center}'. Use 'start', 'end', 'tss', 'tes', or 'midpoint'."
        )


# ---------------------------------------------------------------------------
# GFF3 attribute parsing
# ---------------------------------------------------------------------------

def _parse_gff3_attributes(attr_str: str) -> Dict[str, str]:
    """Parse GFF3 attributes string into a dict.

    Handles: key=value; pairs with URI-encoded characters.
    """
    attrs = {}
    for part in attr_str.rstrip(';').split(';'):
        part = part.strip()
        if '=' in part:
            key, _, value = part.partition('=')
            attrs[key.strip()] = value.strip()
    return attrs


# ---------------------------------------------------------------------------
# Locus loading functions
# ---------------------------------------------------------------------------

def loci_from_tsv(
    tsv_path: Union[str, Path],
    chrom_col: str = 'chrom',
    strand_col: str = 'strand',
    center_col: str = 'center',
    sep: str = '\t',
    additional_cols: Optional[List[str]] = None,
    comment: str = '#',
) -> List[Dict]:
    """
    Load loci from a TSV/CSV file with an explicit center column.

    The simplest loader: wraps any tabular file that already has
    chrom, strand, and center (0-based) columns.

    Args:
        tsv_path: Path to TSV/CSV file
        chrom_col: Column name for chromosome
        strand_col: Column name for strand ('+' or '-')
        center_col: Column name for 0-based center position
        sep: Column separator (default: tab)
        additional_cols: Extra columns to include in locus dicts
        comment: Lines starting with this character are ignored

    Returns:
        List of dicts with 'chrom', 'strand', 'center' plus any additional_cols
    """
    usecols = [chrom_col, strand_col, center_col]
    if additional_cols:
        usecols.extend(c for c in additional_cols if c not in usecols)

    df = pd.read_csv(tsv_path, sep=sep, usecols=usecols, comment=comment)

    loci = []
    for _, row in df.iterrows():
        locus = {
            'chrom': row[chrom_col],
            'strand': row[strand_col],
            'center': int(row[center_col]),
        }
        if additional_cols:
            for col in additional_cols:
                locus[col] = row[col]
        loci.append(locus)

    return loci


def loci_from_bed(
    bed_path: Union[str, Path],
    center: str = 'start',
    default_strand: str = '+',
) -> List[Dict]:
    """
    Load loci from a BED file.

    BED coordinates are 0-based half-open: [start, end).
    Strand is taken from column 5 if present; otherwise uses default_strand.

    Center semantics (in 5'→3' direction):
        'start'/'tss'  → 5' end of feature
        'end'/'tes'    → 3' end (last base, 0-based inclusive)
        'midpoint'     → (start + end - 1) // 2

    Args:
        bed_path: Path to BED file
        center: How to define center position ('start', 'end', 'midpoint',
                'tss', 'tes')
        default_strand: Strand to use if BED file has < 6 columns

    Returns:
        List of dicts with 'chrom', 'strand', 'center', 'name' (if present),
        'bed_start', 'bed_end'
    """
    df = pd.read_csv(
        bed_path, sep='\t', header=None, comment='#',
        names=['chrom', 'start', 'end', 'name', 'score', 'strand'],
        dtype=str,
    )

    has_strand = len(df.columns) >= 6 and df['strand'].dropna().isin(['+', '-']).any()
    center_key = center.lower()

    loci = []
    for _, row in df.iterrows():
        try:
            bed_start = int(row['start'])
            bed_end = int(row['end'])
        except (ValueError, KeyError):
            continue

        strand = row.get('strand', default_strand) if has_strand else default_strand
        if strand not in ('+', '-'):
            strand = default_strand

        if center_key in ('start', 'tss'):
            c = bed_start if strand == '+' else bed_end - 1
        elif center_key in ('end', 'tes'):
            c = bed_end - 1 if strand == '+' else bed_start
        elif center_key == 'midpoint':
            c = (bed_start + bed_end - 1) // 2
        else:
            raise ValueError(
                f"Unknown center='{center}'. Use 'start', 'end', 'tss', 'tes', or 'midpoint'."
            )

        locus: Dict = {
            'chrom': row['chrom'],
            'strand': strand,
            'center': c,
            'feature_start': bed_start,
            'feature_end': bed_end,
        }
        if 'name' in row and pd.notna(row['name']):
            locus['name'] = row['name']
        loci.append(locus)

    return loci


def loci_from_gff(
    gff_path: Union[str, Path],
    feature_type: Optional[str] = None,
    center: str = 'start',
    chrom_col: int = 0,
    feature_col: int = 2,
    start_col: int = 3,
    end_col: int = 4,
    strand_col: int = 6,
    attr_col: int = 8,
    name_attr: Optional[str] = 'Name',
    id_attr: Optional[str] = 'ID',
) -> List[Dict]:
    """
    Load loci from a GFF3 annotation file.

    GFF3 coordinates are 1-based, both endpoints inclusive.
    Converted to 0-based internally (center is always 0-based).

    Center semantics (in 5'→3' direction):
        'start'/'tss'  → TSS / 5' end (feature start for +, feature end for -)
        'end'/'tes'    → TES / 3' end (feature end for +, feature start for -)
        'midpoint'     → geometric midpoint

    Args:
        gff_path: Path to GFF3 file (plain or gzip-compressed)
        feature_type: GFF3 feature type to include (e.g. 'gene', 'mRNA',
                      'CDS', 'three_prime_UTR'). None = include all features.
        center: Center definition ('start', 'end', 'tss', 'tes', 'midpoint')
        chrom_col/feature_col/etc.: Column indices (default: standard GFF3)
        name_attr: GFF3 attribute key for feature name (None = skip)
        id_attr: GFF3 attribute key for feature ID (None = skip)

    Returns:
        List of dicts with 'chrom', 'strand', 'center', 'feature_start',
        'feature_end' (both 0-based), plus 'name' and/or 'id' if found.
    """
    import gzip

    gff_path = Path(gff_path)
    opener = gzip.open if gff_path.suffix == '.gz' else open

    center_key = center.lower()
    loci = []

    with opener(gff_path, 'rt') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue

            if feature_type is not None and cols[feature_col] != feature_type:
                continue

            strand = cols[strand_col]
            if strand not in ('+', '-'):
                continue

            try:
                # GFF3 is 1-based inclusive; convert to 0-based inclusive
                gff_start_1 = int(cols[start_col])
                gff_end_1 = int(cols[end_col])
                feat_start_0 = gff_start_1 - 1   # 0-based start (inclusive)
                feat_end_0 = gff_end_1 - 1        # 0-based end (inclusive)
            except ValueError:
                continue

            if center_key in ('start', 'tss'):
                c = feat_start_0 if strand == '+' else feat_end_0
            elif center_key in ('end', 'tes'):
                c = feat_end_0 if strand == '+' else feat_start_0
            elif center_key == 'midpoint':
                c = (feat_start_0 + feat_end_0) // 2
            else:
                raise ValueError(
                    f"Unknown center='{center}'. Use 'start', 'end', 'tss', 'tes', or 'midpoint'."
                )

            locus: Dict = {
                'chrom': cols[chrom_col],
                'strand': strand,
                'center': c,
                'feature_start': feat_start_0,
                'feature_end': feat_end_0,
            }

            # Parse attributes if requested
            if (name_attr or id_attr) and len(cols) > attr_col:
                attrs = _parse_gff3_attributes(cols[attr_col])
                if name_attr and name_attr in attrs:
                    locus['name'] = attrs[name_attr]
                if id_attr and id_attr in attrs:
                    locus['id'] = attrs[id_attr]

            loci.append(locus)

    return loci


def loci_from_motif_scan(
    sequences: Union[Dict[str, str], Union[str, Path]],
    motif: str,
    iupac: bool = False,
    center: str = 'start',
    strand: str = 'both',
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
    chroms: Optional[List[str]] = None,
    ignore_case: bool = True,
) -> List[Dict]:
    """
    Scan genome sequences for a motif and return match positions as loci.

    This is the generic version of TRT/A-tract scanning: any motif pattern
    can be used (T-tracts, G-quadruplexes, TATA boxes, poly-A signals, etc.).

    Minus strand scanning works by reverse-complementing the chromosome
    sequence and searching for the same pattern, then converting positions
    back to genomic coordinates. This ensures the motif matches the minus
    strand transcript sequence.

    Args:
        sequences: Either a dict of {chrom: sequence_str} or a path to a
                   FASTA file (plain or gzip-compressed)
        motif: Pattern to search for. Interpreted as:
                   - regex if iupac=False (e.g. r'T{5,}', r'TATAAA')
                   - IUPAC string if iupac=True (e.g. 'WTATAWW', 'TTTTTT')
        iupac: If True, convert IUPAC codes in motif to regex before scanning
        center: Center position relative to match ('start', 'end', 'midpoint',
                'tss', 'tes'). Defined in 5'→3' transcription direction.
        strand: Which strand(s) to scan:
                    'both' — scan + strand and RC for - strand matches
                    '+' — plus strand only
                    '-' — minus strand only (scans RC)
        min_length: Minimum match length to include (useful for variable-length
                    patterns like r'T{5,}')
        max_length: Maximum match length to include
        chroms: If provided, scan only these chromosomes
        ignore_case: If True (default), pattern matching is case-insensitive

    Returns:
        List of dicts with keys:
            'chrom', 'strand', 'center',      # always present
            'match_start', 'match_end',        # 0-based genomic, inclusive
            'match_length', 'match_seq'        # match details

    Example:
        # T-tract metagene (equivalent to TRT analysis)
        genome = load_genome("S288C_genome.fasta")
        loci = loci_from_motif_scan(genome, r'T{5,}', center='start')

        # TATA box metagene on promoter windows
        loci = loci_from_motif_scan(genome, 'WTATAWW', iupac=True, center='midpoint')

        # G-quadruplex sites (G4 motif)
        loci = loci_from_motif_scan(genome, r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}')
    """
    # Load genome if path provided
    if not isinstance(sequences, dict):
        from rectify.utils.genome import load_genome
        sequences = load_genome(Path(sequences))

    # Build regex pattern
    pattern_str = iupac_to_regex(motif) if iupac else motif
    flags = re.IGNORECASE if ignore_case else 0
    pattern = re.compile(pattern_str, flags)

    if chroms is not None:
        sequences = {k: v for k, v in sequences.items() if k in chroms}

    loci = []

    for chrom, seq in sequences.items():
        chr_len = len(seq)
        seq_upper = seq.upper()

        # Plus strand scan
        if strand in ('+', 'both'):
            for match in pattern.finditer(seq_upper):
                mlen = match.end() - match.start()
                if min_length is not None and mlen < min_length:
                    continue
                if max_length is not None and mlen > max_length:
                    continue
                c = _center_from_match(match, chr_len, '+', center)
                loci.append({
                    'chrom': chrom,
                    'strand': '+',
                    'center': c,
                    'match_start': match.start(),
                    'match_end': match.end() - 1,
                    'match_length': mlen,
                    'match_seq': match.group(),
                })

        # Minus strand scan: search the reverse complement
        if strand in ('-', 'both'):
            rc_seq = _reverse_complement(seq_upper)
            for match in pattern.finditer(rc_seq):
                mlen = match.end() - match.start()
                if min_length is not None and mlen < min_length:
                    continue
                if max_length is not None and mlen > max_length:
                    continue
                # Convert RC coordinates → genomic coordinates
                # RC position p → genomic position (chr_len - 1 - p)
                genomic_end_inclusive = chr_len - 1 - match.start()    # 5' (highest coord)
                genomic_start_inclusive = chr_len - match.end()         # 3' (lowest coord)
                c = _center_from_match(match, chr_len, '-', center)
                loci.append({
                    'chrom': chrom,
                    'strand': '-',
                    'center': c,
                    'match_start': genomic_start_inclusive,
                    'match_end': genomic_end_inclusive,
                    'match_length': mlen,
                    'match_seq': _reverse_complement(match.group()),  # in genomic orientation
                })

    return loci


def loci_from_pickle(
    cache_path: Union[str, Path],
    chrom_key: str = 'chrom',
    strand_key: str = 'strand',
    center_key: str = 'center',
    additional_keys: Optional[List[str]] = None,
) -> List[Dict]:
    """
    Load loci from a pickle file containing a list of dicts or a dict-with-key.

    Backward-compatible loader for existing PKL locus caches. Normalizes the
    loaded data to the standard locus format {'chrom', 'strand', 'center'}.

    Args:
        cache_path: Path to .pkl file
        chrom_key: Key for chromosome in source dicts
        strand_key: Key for strand in source dicts
        center_key: Key for center position in source dicts
        additional_keys: Extra keys to include from source dicts

    Returns:
        List of dicts with 'chrom', 'strand', 'center' plus any additional_keys
    """
    cache_path = Path(cache_path)
    with open(cache_path, 'rb') as f:
        raw = pickle.load(f)

    # Handle dict-with-loci-key and direct list formats
    if isinstance(raw, dict):
        if 'loci' in raw:
            source = raw['loci']
        else:
            raise ValueError(
                f"Pickle {cache_path.name} is a dict but has no 'loci' key. "
                f"Available keys: {list(raw.keys())}"
            )
    elif isinstance(raw, list):
        source = raw
    else:
        raise ValueError(
            f"Unexpected pickle format in {cache_path.name}: "
            f"expected list or dict, got {type(raw).__name__}"
        )

    loci = []
    for item in source:
        locus: Dict = {
            'chrom': item[chrom_key],
            'strand': item[strand_key],
            'center': int(item[center_key]),
        }
        if additional_keys:
            for key in additional_keys:
                if key in item:
                    locus[key] = item[key]
        loci.append(locus)

    return loci


# ---------------------------------------------------------------------------
# Signal index construction
# ---------------------------------------------------------------------------

def position_index_from_tsv(
    tsv_paths: Union[Union[str, Path], List[Union[str, Path]]],
    position_col: str = 'position',
    chrom_col: str = 'chrom',
    strand_col: str = 'strand',
    count_col: Optional[str] = None,
    sep: str = '\t',
) -> Tuple[PositionIndex, int]:
    """
    Build a PositionIndex from one or more TSV files of genomic positions.

    Each row in the TSV represents one read (or one count if count_col provided).
    Missing files are skipped with a warning.

    Args:
        tsv_paths: Single path or list of paths to TSV files
        position_col: Column name for genomic position (0-based)
        chrom_col: Column name for chromosome
        strand_col: Column name for strand
        count_col: Optional column with pre-computed counts. If provided,
                   each row contributes count_col reads instead of 1.
        sep: Column separator

    Returns:
        Tuple of (PositionIndex, total_reads)

    Example:
        # RECTIFY corrected 3' ends
        index, n = position_index_from_tsv(
            ["wt_rep1.tsv", "wt_rep2.tsv"],
            position_col='corrected_3prime'
        )

        # NET-seq 3' end positions
        index, n = position_index_from_tsv("netseq_wt.tsv", position_col='position')
    """
    if isinstance(tsv_paths, (str, Path)):
        tsv_paths = [tsv_paths]

    usecols = [chrom_col, strand_col, position_col]
    if count_col:
        usecols.append(count_col)

    dfs = []
    for path in tsv_paths:
        path = Path(path)
        if path.exists():
            df = pd.read_csv(path, sep=sep, usecols=usecols)
            dfs.append(df)
        else:
            print(f"WARNING: TSV file not found: {path}")

    if not dfs:
        raise ValueError(f"No TSV files could be loaded from {len(list(tsv_paths))} paths.")

    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.rename(columns={
        position_col: 'position',
        chrom_col: 'chrom',
        strand_col: 'strand',
    })
    if count_col:
        combined = combined.rename(columns={count_col: 'count'})

    total_reads = int(combined['count'].sum()) if count_col else len(combined)
    index = PositionIndex(
        combined,
        position_col='position',
        count_col='count' if count_col else None,
    )

    print(f"Built PositionIndex: {total_reads:,} reads from {len(dfs)} file(s)")
    return index, total_reads


def position_index_from_bigwig(
    bw_path: Union[str, Path],
    strand: str = '+',
    chroms: Optional[List[str]] = None,
    round_values: bool = True,
) -> Tuple[PositionIndex, int]:
    """
    Build a PositionIndex from a BigWig file.

    BigWig files typically store continuous signal values. Each non-zero
    position is treated as a read count (rounded to int by default). This
    works well for BigWigs generated from discrete read endpoints (NET-seq,
    PAR-CLIP, RECTIFY) but less so for smoothed/continuous ChIP-seq signal.

    Requires pyBigWig: ``pip install pyBigWig``

    Args:
        bw_path: Path to BigWig file
        strand: Strand label to assign to all positions ('+' or '-')
                BigWig files are strand-agnostic — you must specify the strand
                based on which file you're loading (e.g., the plus-strand BW)
        chroms: If provided, load only these chromosomes
        round_values: If True (default), round float values to int counts

    Returns:
        Tuple of (PositionIndex, total_reads)
    """
    try:
        import pyBigWig
    except ImportError:
        raise ImportError(
            "pyBigWig is required for position_index_from_bigwig. "
            "Install with: conda install -c bioconda pybigwig"
        )

    bw_path = Path(bw_path)
    bw = pyBigWig.open(str(bw_path))

    rows = []
    chrom_sizes = bw.chroms()
    if chroms:
        chrom_sizes = {c: s for c, s in chrom_sizes.items() if c in chroms}

    for chrom, size in chrom_sizes.items():
        try:
            values = bw.values(chrom, 0, size, numpy=True)
        except RuntimeError:
            continue

        if values is None:
            continue

        nz_positions = np.where(~np.isnan(values) & (values > 0))[0]
        nz_counts = values[nz_positions]

        if round_values:
            nz_counts = np.round(nz_counts).astype(int)

        for pos, count in zip(nz_positions, nz_counts):
            rows.append({
                'chrom': chrom,
                'strand': strand,
                'position': int(pos),
                'count': int(count),
            })

    bw.close()

    if not rows:
        raise ValueError(f"No signal found in BigWig: {bw_path}")

    df = pd.DataFrame(rows)
    total_reads = int(df['count'].sum())
    index = PositionIndex(df, position_col='position', count_col='count')

    print(f"Built PositionIndex from BigWig: {total_reads:,} counts, "
          f"strand='{strand}', chroms={len(chrom_sizes)}")
    return index, total_reads
