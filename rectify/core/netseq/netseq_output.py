"""
NET-seq output generation for RECTIFY.

Generates:
- Parquet: Read-level output (reuses UnifiedReadRecord export)
- BedGraph: Position-level RPM-normalized signal
- BigWig: Position-level RPM-normalized signal (if pyBigWig available)

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np

from ...config import CHROM_SIZES
from ..unified_record import UnifiedReadRecord, export_to_parquet, records_to_dataframe

# Try to import pyBigWig
try:
    import pyBigWig
    HAS_PYBIGWIG = True
except ImportError:
    HAS_PYBIGWIG = False


def calculate_rpm_factor(total_reads: int) -> float:
    """
    Calculate reads-per-million normalization factor.

    Args:
        total_reads: Total number of reads in sample

    Returns:
        RPM factor (1e6 / total_reads)
    """
    return 1e6 / total_reads if total_reads > 0 else 1.0


def position_counts_to_dataframe(
    position_counts: Dict[Tuple[str, str, int], float],
) -> pd.DataFrame:
    """
    Convert position counts dict to DataFrame.

    Args:
        position_counts: Dict mapping (chrom, strand, position) -> count

    Returns:
        DataFrame with chrom, strand, position, count columns
    """
    if not position_counts:
        return pd.DataFrame(columns=['chrom', 'strand', 'position', 'count'])

    data = [
        {'chrom': chrom, 'strand': strand, 'position': pos, 'count': count}
        for (chrom, strand, pos), count in position_counts.items()
    ]

    df = pd.DataFrame(data)
    return df.sort_values(['chrom', 'strand', 'position']).reset_index(drop=True)


def write_bedgraph(
    position_counts: Dict[Tuple[str, str, int], float],
    output_path: Path,
    strand: str,
    total_reads: int,
    normalize_rpm: bool = True,
    track_name: Optional[str] = None,
    track_description: Optional[str] = None,
) -> Path:
    """
    Write strand-specific bedGraph from position counts.

    BedGraph format: chrom<tab>start(0-based)<tab>end(exclusive)<tab>value

    Args:
        position_counts: (chrom, strand, pos) -> count dict
        output_path: Output file path
        strand: '+' or '-'
        total_reads: Total reads for RPM normalization
        normalize_rpm: Apply RPM normalization
        track_name: Track name for UCSC header
        track_description: Track description for UCSC header

    Returns:
        Path to written file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rpm_factor = calculate_rpm_factor(total_reads) if normalize_rpm else 1.0

    # Filter to strand and sort
    strand_data = [
        (chrom, pos, count * rpm_factor)
        for (chrom, s, pos), count in position_counts.items()
        if s == strand
    ]

    if not strand_data:
        # Write empty file with header only
        with open(output_path, 'w') as f:
            if track_name:
                f.write(f'track type=bedGraph name="{track_name}" description="{track_description or track_name}"\n')
        return output_path

    # Sort by chrom then position
    strand_data.sort(key=lambda x: (x[0], x[1]))

    with open(output_path, 'w') as f:
        # Write track header
        if track_name:
            f.write(f'track type=bedGraph name="{track_name}" description="{track_description or track_name}"\n')

        # Write data
        for chrom, pos, value in strand_data:
            # BedGraph: 0-based, half-open coordinates
            f.write(f"{chrom}\t{pos}\t{pos + 1}\t{value:.4f}\n")

    return output_path


def write_bigwig(
    position_counts: Dict[Tuple[str, str, int], float],
    output_path: Path,
    strand: str,
    total_reads: int,
    chrom_sizes: Dict[str, int] = None,
    normalize_rpm: bool = True,
) -> Optional[Path]:
    """
    Write strand-specific bigWig from position counts.

    Args:
        position_counts: (chrom, strand, pos) -> count dict
        output_path: Output file path
        strand: '+' or '-'
        total_reads: Total reads for RPM normalization
        chrom_sizes: Chromosome sizes dict (default: from config)
        normalize_rpm: Apply RPM normalization

    Returns:
        Path to written file, or None if pyBigWig not available
    """
    if not HAS_PYBIGWIG:
        return None

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if chrom_sizes is None:
        chrom_sizes = CHROM_SIZES

    rpm_factor = calculate_rpm_factor(total_reads) if normalize_rpm else 1.0

    # Filter to strand and organize by chromosome
    chrom_data: Dict[str, List[Tuple[int, float]]] = {}
    for (chrom, s, pos), count in position_counts.items():
        if s != strand:
            continue
        if chrom not in chrom_data:
            chrom_data[chrom] = []
        chrom_data[chrom].append((pos, count * rpm_factor))

    if not chrom_data:
        return None

    # Create bigWig file
    # Need to provide chromosome sizes in order
    ordered_chroms = sorted(chrom_sizes.keys(), key=lambda x: (len(x), x))
    header = [(chrom, chrom_sizes[chrom]) for chrom in ordered_chroms if chrom in chrom_data]

    if not header:
        return None

    bw = pyBigWig.open(str(output_path), "w")
    bw.addHeader(header)

    # Add data for each chromosome
    for chrom in ordered_chroms:
        if chrom not in chrom_data:
            continue

        # Sort positions and convert to arrays
        positions = sorted(chrom_data[chrom], key=lambda x: x[0])
        starts = [p[0] for p in positions]
        ends = [p[0] + 1 for p in positions]
        values = [p[1] for p in positions]

        bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)

    bw.close()

    return output_path


def export_netseq_results(
    records: List[UnifiedReadRecord],
    raw_counts: Dict[Tuple[str, str, int], float],
    deconv_counts: Optional[Dict[Tuple[str, str, int], float]],
    output_dir: Path,
    sample_name: str,
    export_parquet: bool = True,
    export_bedgraph: bool = True,
    export_bigwig: bool = True,
    normalize_rpm: bool = True,
    total_reads: Optional[int] = None,
    corrected_counts: Optional[Dict[Tuple[str, str, int], float]] = None,
) -> Dict[str, Path]:
    """
    Export all NET-seq output files.

    Creates:
    - {sample_name}.unified_reads.parquet (read-level)
    - {sample_name}.raw.plus.bedgraph (RAW alignment terminus, plus strand)
    - {sample_name}.raw.minus.bedgraph (RAW alignment terminus, minus strand)
    - {sample_name}.corrected.plus.bedgraph (after tail walkback + junction rescue)
    - {sample_name}.corrected.minus.bedgraph
    - {sample_name}.deconv.plus.bedgraph (deconvolved, plus strand)
    - {sample_name}.deconv.minus.bedgraph (deconvolved, minus strand)
    - (Same with .bw extension if bigwig enabled)

    ``raw`` and ``corrected`` are BOTH written whenever ``corrected_counts`` is supplied, so a run
    can be diffed against the pre-correction track without a second pass. ``raw`` keeps exactly the
    contents it had before the correction pass existed.

    Args:
        records: List of UnifiedReadRecord
        raw_counts: Raw position counts (three_prime_raw)
        deconv_counts: Deconvolved position counts (None to skip)
        output_dir: Output directory
        sample_name: Sample name for output files
        export_parquet: Generate parquet file
        export_bedgraph: Generate bedgraph files
        export_bigwig: Generate bigwig files
        normalize_rpm: Apply RPM normalization
        total_reads: Read total for RPM normalization. Defaults to len(records); pass it
            explicitly when records were streamed (not materialised) into raw_counts.
        corrected_counts: Corrected-3'-end position counts (three_prime_corrected). None to skip.

    Returns:
        Dict mapping output type to path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs: Dict[str, Path] = {}

    # Calculate total reads for normalization
    if total_reads is None:
        total_reads = len(records) if records else int(sum(raw_counts.values()))

    # Export parquet (read-level)
    if export_parquet and records:
        parquet_path = output_dir / f"{sample_name}.unified_reads.parquet"
        export_to_parquet(records, parquet_path)
        outputs['parquet'] = parquet_path
        print(f"  Exported parquet: {parquet_path.name} ({len(records):,} records)")

    # Export raw bedgraph
    if export_bedgraph and raw_counts:
        for strand, strand_name in [('+', 'plus'), ('-', 'minus')]:
            bg_path = output_dir / f"{sample_name}.raw.{strand_name}.bedgraph"
            track_name = f"{sample_name}_raw_{strand_name}"
            write_bedgraph(
                raw_counts, bg_path, strand, total_reads,
                normalize_rpm=normalize_rpm,
                track_name=track_name,
            )
            outputs[f'raw_{strand_name}_bedgraph'] = bg_path
            n_positions = sum(1 for (_, s, _) in raw_counts if s == strand)
            print(f"  Exported bedgraph: {bg_path.name} ({n_positions:,} positions)")

    # Export corrected bedgraph (tail walkback + donor-side junction rescue)
    if export_bedgraph and corrected_counts:
        for strand, strand_name in [('+', 'plus'), ('-', 'minus')]:
            bg_path = output_dir / f"{sample_name}.corrected.{strand_name}.bedgraph"
            track_name = f"{sample_name}_corrected_{strand_name}"
            write_bedgraph(
                corrected_counts, bg_path, strand, total_reads,
                normalize_rpm=normalize_rpm,
                track_name=track_name,
            )
            outputs[f'corrected_{strand_name}_bedgraph'] = bg_path
            n_positions = sum(1 for (_, s, _) in corrected_counts if s == strand)
            print(f"  Exported bedgraph: {bg_path.name} ({n_positions:,} positions)")

    # Export deconvolved bedgraph
    if export_bedgraph and deconv_counts:
        for strand, strand_name in [('+', 'plus'), ('-', 'minus')]:
            bg_path = output_dir / f"{sample_name}.deconv.{strand_name}.bedgraph"
            track_name = f"{sample_name}_deconv_{strand_name}"
            write_bedgraph(
                deconv_counts, bg_path, strand, total_reads,
                normalize_rpm=normalize_rpm,
                track_name=track_name,
            )
            outputs[f'deconv_{strand_name}_bedgraph'] = bg_path
            n_positions = sum(1 for (_, s, _) in deconv_counts if s == strand)
            print(f"  Exported bedgraph: {bg_path.name} ({n_positions:,} positions)")

    # Export bigwig files
    if export_bigwig and not HAS_PYBIGWIG:
        print("  WARNING: --output-format bigwig requested but pyBigWig is not importable -- "
              "NO .bw files were written (install with: pip install pyBigWig)")
    if export_bigwig and HAS_PYBIGWIG:
        if raw_counts:
            for strand, strand_name in [('+', 'plus'), ('-', 'minus')]:
                bw_path = output_dir / f"{sample_name}.raw.{strand_name}.bw"
                result = write_bigwig(
                    raw_counts, bw_path, strand, total_reads,
                    normalize_rpm=normalize_rpm,
                )
                if result:
                    outputs[f'raw_{strand_name}_bigwig'] = result
                    print(f"  Exported bigwig: {bw_path.name}")

        if corrected_counts:
            for strand, strand_name in [('+', 'plus'), ('-', 'minus')]:
                bw_path = output_dir / f"{sample_name}.corrected.{strand_name}.bw"
                result = write_bigwig(
                    corrected_counts, bw_path, strand, total_reads,
                    normalize_rpm=normalize_rpm,
                )
                if result:
                    outputs[f'corrected_{strand_name}_bigwig'] = result
                    print(f"  Exported bigwig: {bw_path.name}")

        if deconv_counts:
            for strand, strand_name in [('+', 'plus'), ('-', 'minus')]:
                bw_path = output_dir / f"{sample_name}.deconv.{strand_name}.bw"
                result = write_bigwig(
                    deconv_counts, bw_path, strand, total_reads,
                    normalize_rpm=normalize_rpm,
                )
                if result:
                    outputs[f'deconv_{strand_name}_bigwig'] = result
                    print(f"  Exported bigwig: {bw_path.name}")

    return outputs


def write_exclusion_stats(
    total_reads: int,
    n_excluded_rdna: int,
    n_excluded_tRNA: int,
    n_excluded_other: int,
    output_path: Path,
) -> Path:
    """
    Write exclusion statistics summary.

    Args:
        total_reads: Total reads processed
        n_excluded_rdna: Reads excluded for rDNA
        n_excluded_tRNA: Reads excluded for tRNA genes
        n_excluded_other: Reads excluded for other reasons
        output_path: Output file path

    Returns:
        Path to written file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_excluded = n_excluded_rdna + n_excluded_tRNA + n_excluded_other
    n_kept = total_reads - n_excluded

    with open(output_path, 'w') as f:
        f.write("category\tcount\tpercent\n")
        f.write(f"total_reads\t{total_reads}\t100.0\n")
        f.write(f"kept\t{n_kept}\t{100*n_kept/total_reads:.2f}\n")
        f.write(f"excluded_total\t{n_excluded}\t{100*n_excluded/total_reads:.2f}\n")
        f.write(f"excluded_rdna\t{n_excluded_rdna}\t{100*n_excluded_rdna/total_reads:.2f}\n")
        f.write(f"excluded_tRNA\t{n_excluded_tRNA}\t{100*n_excluded_tRNA/total_reads:.2f}\n")
        f.write(f"excluded_other\t{n_excluded_other}\t{100*n_excluded_other/total_reads:.2f}\n")

    return output_path
