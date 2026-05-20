#!/usr/bin/env python3
"""
Region planning for parallel BAM processing.

Splits chromosomes at coverage gaps (or uniformly when no gaps exist)
to produce balanced work units for the parallel correction pipeline.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import pysam


@dataclass(frozen=True)
class RegionPlan:
    """A single region's work plan for parallel BAM processing.

    Carries an id (stable across runs given the same input + min_gap_size)
    and a per-region temp dir so workers can write outputs without colliding.
    """
    region_id: str          # e.g. "r000", "r001"; left-padded so lex-sort == coord-sort
    chrom: str
    start: int              # 0-based, inclusive
    end: int                # 0-based, exclusive (pysam fetch convention)
    tmp_dir: Path           # per-region scratch dir (workers write region_NNN.bam etc. here)

    @property
    def region_bam(self) -> Path:
        return self.tmp_dir / f"{self.region_id}.bam"

    @property
    def ok_sentinel(self) -> Path:
        """Marker file written ONLY after region_bam is fully sorted + fsync'd.
        Existence -> this region's work is durably complete."""
        return self.tmp_dir / f"{self.region_id}.bam.ok"


def plan_regions(
    bam_path: str,
    tmp_dir: Path,
    min_gap_size: int = 10000,
    max_region_size: int = 100000,
) -> List[RegionPlan]:
    """Wrap get_processing_regions() into RegionPlan objects.

    Region ids are assigned in coordinate order using fixed-width zero-padding
    so that lex-sort matches coord-sort (matters for samtools merge of
    'region_*.bam' globbed in shell order).
    """
    raw = get_processing_regions(bam_path, min_gap_size=min_gap_size, max_region_size=max_region_size)
    width = max(3, len(str(max(len(raw) - 1, 0))))  # at least 3 digits, more if needed
    plans = []
    for i, (chrom, start, end) in enumerate(raw):
        plans.append(RegionPlan(
            region_id=f"r{i:0{width}d}",
            chrom=chrom,
            start=start,
            end=end,
            tmp_dir=tmp_dir,
        ))
    return plans


def find_coverage_gaps(
    bam_path: str,
    chrom: str,
    min_gap_size: int = 10000,
) -> List[Tuple[int, int]]:
    """
    Find gaps in read coverage that can be used to split chromosome.

    Large chromosomes can be split at coverage deserts where no reads exist,
    allowing independent parallel processing of regions.

    Args:
        bam_path: Path to BAM file
        chrom: Chromosome name
        min_gap_size: Minimum gap size in bp to consider a split point

    Returns:
        List of (start, end) tuples representing gap regions
    """
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        # Get chromosome length
        chrom_length = None
        for idx, ref in enumerate(bam.references):
            if ref == chrom:
                chrom_length = bam.lengths[idx]
                break

        if chrom_length is None:
            return []

        # Collect read intervals and merge to find covered regions
        # Uses O(reads) memory instead of O(genome_positions)
        intervals = []
        for read in bam.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            intervals.append((read.reference_start, read.reference_end))

    # Use index to efficiently find gaps
    # Strategy: sample positions and find where no reads map
    gaps = []
    last_covered = 0

    if not intervals:
        return [(0, chrom_length)]

    # Sort and merge overlapping intervals
    intervals.sort()
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])

    # Find gaps between merged intervals
    prev_end = 0
    for seg_start, seg_end in merged:
        if seg_start - prev_end >= min_gap_size:
            gaps.append((prev_end, seg_start))
        prev_end = seg_end

    # Check final gap
    if chrom_length - prev_end >= min_gap_size:
        gaps.append((prev_end, chrom_length))

    return gaps


def get_processing_regions(
    bam_path: str,
    min_gap_size: int = 10000,
    max_region_size: int = 100000,
) -> List[Tuple[str, int, int]]:
    """
    Get list of regions for parallel processing.

    Splits large chromosomes at coverage gaps to create balanced work units.
    Each region is capped at max_region_size reference bases to limit per-worker
    result list size and prevent OOM when many workers run in parallel.
    The 100k default limits peak worker memory to ~150 MB on a 1552x-coverage
    nanopore DRS yeast dataset (vs ~1 GB at 500k). Workers with very high
    local coverage (e.g. RPL genes at 50x) stay under ~300 MB.

    Args:
        bam_path: Path to BAM file
        min_gap_size: Minimum gap to use as split point
        max_region_size: Target maximum region size

    Returns:
        List of (chrom, start, end) tuples representing processing regions
    """
    bam = pysam.AlignmentFile(bam_path, 'rb')
    chromosomes = list(bam.references)
    chrom_lengths = dict(zip(bam.references, bam.lengths))
    bam.close()

    regions = []

    for chrom in chromosomes:
        chrom_len = chrom_lengths[chrom]

        # For small chromosomes, process as single region
        if chrom_len <= max_region_size:
            regions.append((chrom, 0, chrom_len))
            continue

        # For large chromosomes, split at gaps
        gaps = find_coverage_gaps(bam_path, chrom, min_gap_size)

        if not gaps:
            # No gaps found — split uniformly by coordinate so workers don't
            # accumulate millions of reads in one list.
            n_sub = max(1, (chrom_len + max_region_size - 1) // max_region_size)
            sub_size = chrom_len // n_sub
            for i in range(n_sub):
                sub_start = i * sub_size
                sub_end = chrom_len if i == n_sub - 1 else (i + 1) * sub_size
                regions.append((chrom, sub_start, sub_end))
            continue

        # Create regions between gaps
        current_start = 0
        for gap_start, gap_end in gaps:
            if gap_start > current_start:
                regions.append((chrom, current_start, gap_start))
            current_start = gap_end

        # Final region after last gap
        if current_start < chrom_len:
            regions.append((chrom, current_start, chrom_len))

    return regions
