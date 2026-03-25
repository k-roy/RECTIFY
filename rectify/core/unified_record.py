#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified Read Record Module for RECTIFY

This module defines the unified per-read data structure that captures transcript ends,
splice junctions, and soft clips in a single pass. This enables efficient downstream
derivation of multiple analysis types without repeated BAM reprocessing.

Data Structure:
    Each read is represented as a UnifiedReadRecord with:
    - Transcript end positions (5' and 3', raw and corrected)
    - Splice junction information (from multi-aligner consensus)
    - Soft clip details (5' and 3' ends, poly(A) information)
    - Alignment quality metrics

Output Formats:
    - Parquet: Primary format for analysis (columnar, efficient, ~3-5x smaller)
    - TSV: Secondary format for debugging and inspection (JSON-encoded junctions)

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-24
"""

from dataclasses import dataclass, asdict, field
from typing import List, Tuple, Dict, Optional, Any
from pathlib import Path
import json
import pandas as pd
import numpy as np

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAS_PYARROW = True
except ImportError:
    HAS_PYARROW = False


@dataclass
class UnifiedReadRecord:
    """
    Per-read record with all extracted features.

    This is the fundamental unit of the unified data structure, capturing:
    - Both transcript ends (5' and 3')
    - Splice junction information
    - Soft clip details
    - Alignment quality metrics

    Coordinate System:
        All positions use 0-based, half-open coordinates consistent with pysam/BED.

        For strand-aware interpretation:
        - Plus strand: 5' = leftmost, 3' = rightmost
        - Minus strand: 5' = rightmost, 3' = leftmost
    """
    # Read identification
    read_id: str
    chrom: str
    strand: str

    # 5' end positions (TSS)
    five_prime_raw: int              # Raw leftmost/rightmost aligned position
    five_prime_corrected: int        # After splice-aware correction
    first_exon_start: Optional[int]  # Start of first exon (may differ from raw)
    starts_in_intron: bool           # True if read starts within intron

    # 3' end positions (CPA)
    three_prime_raw: int             # Raw mapped 3' end
    three_prime_corrected: int       # After poly(A) correction

    # Alignment span
    alignment_start: int             # Leftmost aligned position
    alignment_end: int               # Rightmost position + 1 (exclusive)

    # Splice junctions
    junctions: List[Tuple[int, int]] = field(default_factory=list)
    n_junctions: int = 0
    junctions_filtered: int = 0      # Number filtered as poly(A) artifacts

    # Soft clips
    five_prime_soft_clip_length: int = 0
    three_prime_soft_clip_length: int = 0
    five_prime_soft_clip_seq: str = ""
    three_prime_soft_clip_seq: str = ""

    # Poly(A) information
    polya_length: int = 0            # Total observed poly(A) length
    aligned_a_length: int = 0        # A's in aligned region
    soft_clip_a_length: int = 0      # A's in soft clip

    # Alignment quality
    mapq: int = 0
    cigar_summary: str = ""          # e.g., "100M2I50M1D30M"
    alignment_identity: float = 0.0

    # Consensus information (from multi-aligner)
    best_aligner: str = ""           # Which aligner was selected
    consensus_confidence: str = ""   # high/medium/low
    n_aligners_support: int = 0      # Number of aligners supporting junctions

    # Count (for collapsed reads)
    count: int = 1

    def to_dict(self) -> Dict[str, Any]:
        """Convert record to dictionary for DataFrame construction."""
        d = asdict(self)
        # Convert junctions list to string for TSV compatibility
        d['junctions_str'] = self.junctions_to_string()
        return d

    def junctions_to_string(self) -> str:
        """Convert junctions list to compact string format."""
        if not self.junctions:
            return ""
        return ";".join(f"{s}-{e}" for s, e in self.junctions)

    @staticmethod
    def junctions_from_string(s: str) -> List[Tuple[int, int]]:
        """Parse junctions from string format."""
        if not s:
            return []
        junctions = []
        for part in s.split(";"):
            if "-" in part:
                start, end = part.split("-")
                junctions.append((int(start), int(end)))
        return junctions


def records_to_dataframe(records: List[UnifiedReadRecord]) -> pd.DataFrame:
    """
    Convert list of UnifiedReadRecord to pandas DataFrame.

    Args:
        records: List of UnifiedReadRecord objects

    Returns:
        DataFrame with all record fields as columns
    """
    if not records:
        return pd.DataFrame()

    data = [r.to_dict() for r in records]
    df = pd.DataFrame(data)

    return df


def dataframe_to_records(df: pd.DataFrame) -> List[UnifiedReadRecord]:
    """
    Convert DataFrame back to list of UnifiedReadRecord.

    Args:
        df: DataFrame with record columns

    Returns:
        List of UnifiedReadRecord objects
    """
    records = []
    for _, row in df.iterrows():
        # Parse junctions from string if present
        if 'junctions_str' in row:
            junctions = UnifiedReadRecord.junctions_from_string(row['junctions_str'])
        elif 'junctions' in row and isinstance(row['junctions'], str):
            junctions = UnifiedReadRecord.junctions_from_string(row['junctions'])
        elif 'junctions' in row:
            junctions = row['junctions']
        else:
            junctions = []

        record = UnifiedReadRecord(
            read_id=row['read_id'],
            chrom=row['chrom'],
            strand=row['strand'],
            five_prime_raw=int(row['five_prime_raw']),
            five_prime_corrected=int(row['five_prime_corrected']),
            first_exon_start=int(row['first_exon_start']) if pd.notna(row.get('first_exon_start')) else None,
            starts_in_intron=bool(row['starts_in_intron']),
            three_prime_raw=int(row['three_prime_raw']),
            three_prime_corrected=int(row['three_prime_corrected']),
            alignment_start=int(row['alignment_start']),
            alignment_end=int(row['alignment_end']),
            junctions=junctions,
            n_junctions=int(row.get('n_junctions', len(junctions))),
            junctions_filtered=int(row.get('junctions_filtered', 0)),
            five_prime_soft_clip_length=int(row.get('five_prime_soft_clip_length', 0)),
            three_prime_soft_clip_length=int(row.get('three_prime_soft_clip_length', 0)),
            five_prime_soft_clip_seq=str(row.get('five_prime_soft_clip_seq', '')),
            three_prime_soft_clip_seq=str(row.get('three_prime_soft_clip_seq', '')),
            polya_length=int(row.get('polya_length', 0)),
            aligned_a_length=int(row.get('aligned_a_length', 0)),
            soft_clip_a_length=int(row.get('soft_clip_a_length', 0)),
            mapq=int(row.get('mapq', 0)),
            cigar_summary=str(row.get('cigar_summary', '')),
            alignment_identity=float(row.get('alignment_identity', 0.0)),
            best_aligner=str(row.get('best_aligner', '')),
            consensus_confidence=str(row.get('consensus_confidence', '')),
            n_aligners_support=int(row.get('n_aligners_support', 0)),
            count=int(row.get('count', 1)),
        )
        records.append(record)

    return records


# =============================================================================
# Export Functions
# =============================================================================

def export_to_parquet(
    records: List[UnifiedReadRecord],
    output_path: Path,
    compression: str = 'snappy',
) -> Path:
    """
    Export records to Parquet format.

    This is the primary format for analysis - columnar, efficient, ~3-5x smaller.

    Args:
        records: List of UnifiedReadRecord objects
        output_path: Path to output file
        compression: Compression codec ('snappy', 'gzip', 'zstd', or None)

    Returns:
        Path to written file
    """
    if not HAS_PYARROW:
        raise ImportError("pyarrow is required for Parquet export. Install with: pip install pyarrow")

    df = records_to_dataframe(records)

    # Convert junctions list to JSON string for Parquet storage
    df['junctions_json'] = df['junctions'].apply(
        lambda x: json.dumps(x) if isinstance(x, list) else json.dumps([])
    )
    df = df.drop(columns=['junctions'])

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_parquet(output_path, compression=compression, index=False)

    return output_path


def export_to_tsv(
    records: List[UnifiedReadRecord],
    output_path: Path,
    include_sequences: bool = False,
) -> Path:
    """
    Export records to TSV format.

    This is the secondary format for debugging and inspection.

    Args:
        records: List of UnifiedReadRecord objects
        output_path: Path to output file
        include_sequences: If True, include soft clip sequences (larger files)

    Returns:
        Path to written file
    """
    df = records_to_dataframe(records)

    # Use junction string format instead of list
    if 'junctions' in df.columns:
        df = df.drop(columns=['junctions'])

    if not include_sequences:
        df = df.drop(columns=['five_prime_soft_clip_seq', 'three_prime_soft_clip_seq'], errors='ignore')

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_path, sep='\t', index=False)

    return output_path


def export_both(
    records: List[UnifiedReadRecord],
    output_dir: Path,
    sample_name: str,
    include_sequences: bool = False,
) -> Dict[str, Path]:
    """
    Export records to both Parquet and TSV formats.

    Args:
        records: List of UnifiedReadRecord objects
        output_dir: Directory for output files
        sample_name: Sample name for file naming
        include_sequences: If True, include soft clip sequences in TSV

    Returns:
        Dict with 'parquet' and 'tsv' paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {}

    # Parquet (primary)
    parquet_path = output_dir / f"{sample_name}.unified_reads.parquet"
    paths['parquet'] = export_to_parquet(records, parquet_path)

    # TSV (secondary)
    tsv_path = output_dir / f"{sample_name}.unified_reads.tsv"
    paths['tsv'] = export_to_tsv(records, tsv_path, include_sequences=include_sequences)

    return paths


# =============================================================================
# Import Functions
# =============================================================================

def import_from_parquet(input_path: Path) -> List[UnifiedReadRecord]:
    """
    Import records from Parquet format.

    Args:
        input_path: Path to Parquet file

    Returns:
        List of UnifiedReadRecord objects
    """
    if not HAS_PYARROW:
        raise ImportError("pyarrow is required for Parquet import. Install with: pip install pyarrow")

    df = pd.read_parquet(input_path)

    # Parse junctions from JSON
    if 'junctions_json' in df.columns:
        df['junctions'] = df['junctions_json'].apply(
            lambda x: json.loads(x) if x else []
        )
        df['junctions_str'] = df['junctions'].apply(
            lambda x: ";".join(f"{s}-{e}" for s, e in x) if x else ""
        )
        df = df.drop(columns=['junctions_json'])

    return dataframe_to_records(df)


def import_from_tsv(input_path: Path) -> List[UnifiedReadRecord]:
    """
    Import records from TSV format.

    Args:
        input_path: Path to TSV file

    Returns:
        List of UnifiedReadRecord objects
    """
    df = pd.read_csv(input_path, sep='\t')

    return dataframe_to_records(df)


def import_unified_records(input_path: Path) -> List[UnifiedReadRecord]:
    """
    Import records from either Parquet or TSV format (auto-detected).

    Args:
        input_path: Path to input file (.parquet or .tsv)

    Returns:
        List of UnifiedReadRecord objects
    """
    input_path = Path(input_path)

    if input_path.suffix == '.parquet':
        return import_from_parquet(input_path)
    elif input_path.suffix in ('.tsv', '.txt'):
        return import_from_tsv(input_path)
    else:
        raise ValueError(f"Unknown file format: {input_path.suffix}")


# =============================================================================
# Validation
# =============================================================================

def validate_record(record: UnifiedReadRecord) -> List[str]:
    """
    Validate a UnifiedReadRecord for consistency.

    Args:
        record: Record to validate

    Returns:
        List of validation error messages (empty if valid)
    """
    errors = []

    # Check strand
    if record.strand not in ('+', '-'):
        errors.append(f"Invalid strand: {record.strand}")

    # Check coordinate consistency
    if record.alignment_start > record.alignment_end:
        errors.append(f"alignment_start ({record.alignment_start}) > alignment_end ({record.alignment_end})")

    # Check 5' position is within alignment
    if record.strand == '+':
        if record.five_prime_raw < record.alignment_start:
            errors.append(f"five_prime_raw ({record.five_prime_raw}) < alignment_start")
    else:
        if record.five_prime_raw > record.alignment_end - 1:
            errors.append(f"five_prime_raw ({record.five_prime_raw}) > alignment_end - 1")

    # Check junction count consistency
    if record.n_junctions != len(record.junctions):
        errors.append(f"n_junctions ({record.n_junctions}) != len(junctions) ({len(record.junctions)})")

    # Check junction coordinates
    for i, (start, end) in enumerate(record.junctions):
        if start >= end:
            errors.append(f"Junction {i}: start ({start}) >= end ({end})")
        if start < record.alignment_start or end > record.alignment_end:
            errors.append(f"Junction {i} outside alignment span")

    return errors


if __name__ == '__main__':
    # Quick test
    print("Testing UnifiedReadRecord...")

    record = UnifiedReadRecord(
        read_id="test_read_001",
        chrom="chrI",
        strand="+",
        five_prime_raw=1000,
        five_prime_corrected=1000,
        first_exon_start=1000,
        starts_in_intron=False,
        three_prime_raw=2000,
        three_prime_corrected=1995,
        alignment_start=1000,
        alignment_end=2000,
        junctions=[(1200, 1400), (1600, 1800)],
        n_junctions=2,
        junctions_filtered=0,
        polya_length=25,
        mapq=60,
        best_aligner="mapPacBio",
        consensus_confidence="high",
    )

    # Validate
    errors = validate_record(record)
    if errors:
        print(f"Validation errors: {errors}")
    else:
        print("Record validated successfully!")

    # Test round-trip
    df = records_to_dataframe([record])
    print(f"\nDataFrame columns: {list(df.columns)}")

    records_back = dataframe_to_records(df)
    print(f"Round-trip successful: {records_back[0].read_id == record.read_id}")

    print("\nUnifiedReadRecord module ready!")
