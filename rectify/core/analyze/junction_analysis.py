#!/usr/bin/env python3
"""
Junction Analysis Module for RECTIFY

Implements IsoQuant-style splice site tolerance for junction matching.
Long-read sequencing often has small errors at splice boundaries due to
basecalling artifacts and alignment uncertainty. This module allows
flexible matching of read junctions to known/annotated junctions within
a configurable tolerance.

Key Features:
- Splice site tolerance (default: ±10bp)
- Junction coordinate correction (snap to known)
- Novel junction detection with tolerance-aware deduplication

Reference:
    IsoQuant: Accurate isoform discovery with long reads
    (Nature Biotechnology, 2022)

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-25
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple, Union
from collections import defaultdict
import numpy as np

# Import UnifiedReadRecord for type hints
from rectify.core.unified_record import UnifiedReadRecord


# Default tolerance for splice site matching
DEFAULT_SPLICE_SITE_TOLERANCE = 10  # bp


# =============================================================================
# Junction Matching with Tolerance
# =============================================================================

def match_junction_with_tolerance(
    read_junction: Tuple[int, int],
    known_junctions: Union[List[Tuple[int, int]], Set[Tuple[int, int]]],
    tolerance: int = DEFAULT_SPLICE_SITE_TOLERANCE,
) -> Optional[Tuple[int, int]]:
    """
    Match a read junction to a known junction within tolerance.

    Long reads have systematic splice site shifts due to:
    - Basecalling errors near splice boundaries
    - Alignment artifacts at GT-AG motifs
    - Homopolymer runs near splice sites

    This function finds the closest known junction within ±tolerance bp
    of both donor and acceptor sites.

    Args:
        read_junction: (donor, acceptor) from read alignment
        known_junctions: Collection of known (donor, acceptor) junctions
        tolerance: Maximum bp distance for matching (default: 10)

    Returns:
        Matched known junction, or None if no match within tolerance
    """
    donor, acceptor = read_junction

    best_match = None
    best_distance = float('inf')

    for known_donor, known_acceptor in known_junctions:
        donor_dist = abs(donor - known_donor)
        acceptor_dist = abs(acceptor - known_acceptor)

        # Both ends must be within tolerance
        if donor_dist <= tolerance and acceptor_dist <= tolerance:
            # Use sum of distances as match quality metric
            total_dist = donor_dist + acceptor_dist

            if total_dist < best_distance:
                best_distance = total_dist
                best_match = (known_donor, known_acceptor)

    return best_match


def match_junction_to_index(
    read_junction: Tuple[int, int],
    junction_index: Dict[Tuple[int, int], int],
    tolerance: int = DEFAULT_SPLICE_SITE_TOLERANCE,
) -> Optional[int]:
    """
    Match a read junction to a junction index with tolerance.

    More efficient when you have many reads to match against a fixed
    set of junctions.

    Args:
        read_junction: (donor, acceptor) from read alignment
        junction_index: Dictionary mapping (donor, acceptor) to index
        tolerance: Maximum bp distance for matching

    Returns:
        Index of matched junction, or None if no match
    """
    matched = match_junction_with_tolerance(
        read_junction,
        junction_index.keys(),
        tolerance
    )

    if matched:
        return junction_index[matched]
    return None


# =============================================================================
# Junction Coordinate Correction
# =============================================================================

def correct_junction_coordinates(
    record: UnifiedReadRecord,
    known_junctions: Set[Tuple[int, int]],
    tolerance: int = DEFAULT_SPLICE_SITE_TOLERANCE,
) -> UnifiedReadRecord:
    """
    Correct read junction coordinates to match nearby known junctions.

    If a read junction is within ±tolerance of a known junction,
    snap it to the known coordinates. This reduces false novel
    junction calls caused by alignment noise.

    Args:
        record: UnifiedReadRecord to correct
        known_junctions: Set of known (donor, acceptor) junctions
        tolerance: Maximum bp distance for snapping

    Returns:
        New UnifiedReadRecord with corrected junction coordinates
    """
    if not record.junctions:
        return record

    corrected_junctions = []
    n_corrected = 0

    for junction in record.junctions:
        matched = match_junction_with_tolerance(junction, known_junctions, tolerance)

        if matched:
            corrected_junctions.append(matched)
            if matched != junction:
                n_corrected += 1
        else:
            # Keep as novel (no matching known junction)
            corrected_junctions.append(junction)

    # Return new record if any changes were made
    if n_corrected > 0:
        return UnifiedReadRecord(
            read_id=record.read_id,
            chrom=record.chrom,
            strand=record.strand,
            five_prime_raw=record.five_prime_raw,
            five_prime_corrected=record.five_prime_corrected,
            first_exon_start=record.first_exon_start,
            starts_in_intron=record.starts_in_intron,
            three_prime_raw=record.three_prime_raw,
            three_prime_corrected=record.three_prime_corrected,
            alignment_start=record.alignment_start,
            alignment_end=record.alignment_end,
            junctions=corrected_junctions,
            n_junctions=len(corrected_junctions),
            junctions_filtered=record.junctions_filtered,
            five_prime_soft_clip_length=record.five_prime_soft_clip_length,
            three_prime_soft_clip_length=record.three_prime_soft_clip_length,
            five_prime_soft_clip_seq=record.five_prime_soft_clip_seq,
            three_prime_soft_clip_seq=record.three_prime_soft_clip_seq,
            polya_length=record.polya_length,
            aligned_a_length=record.aligned_a_length,
            soft_clip_a_length=record.soft_clip_a_length,
            mapq=record.mapq,
            cigar_summary=record.cigar_summary,
            alignment_identity=record.alignment_identity,
            best_aligner=record.best_aligner,
            consensus_confidence=record.consensus_confidence,
            n_aligners_support=record.n_aligners_support,
            count=record.count,
        )

    return record


def correct_all_junction_coordinates(
    records: List[UnifiedReadRecord],
    known_junctions: Set[Tuple[int, int]],
    tolerance: int = DEFAULT_SPLICE_SITE_TOLERANCE,
) -> Tuple[List[UnifiedReadRecord], Dict[str, int]]:
    """
    Correct junction coordinates for all records.

    Args:
        records: List of UnifiedReadRecord objects
        known_junctions: Set of known (donor, acceptor) junctions
        tolerance: Maximum bp distance for snapping

    Returns:
        Tuple of:
        - List of corrected records
        - Statistics dictionary
    """
    corrected_records = []
    stats = {
        'total_reads': len(records),
        'reads_with_junctions': 0,
        'reads_corrected': 0,
        'junctions_total': 0,
        'junctions_matched': 0,
        'junctions_novel': 0,
    }

    for record in records:
        if record.junctions:
            stats['reads_with_junctions'] += 1
            stats['junctions_total'] += len(record.junctions)

        corrected = correct_junction_coordinates(record, known_junctions, tolerance)

        # Check if any junctions were corrected
        if corrected.junctions != record.junctions:
            stats['reads_corrected'] += 1

        # Count matched vs novel
        for junction in corrected.junctions:
            if junction in known_junctions:
                stats['junctions_matched'] += 1
            else:
                stats['junctions_novel'] += 1

        corrected_records.append(corrected)

    return corrected_records, stats


# =============================================================================
# Novel Junction Deduplication
# =============================================================================

def deduplicate_novel_junctions(
    novel_junctions: List[Tuple[int, int]],
    tolerance: int = DEFAULT_SPLICE_SITE_TOLERANCE,
) -> List[Tuple[int, int]]:
    """
    Deduplicate novel junctions that are within tolerance of each other.

    Multiple reads may call the same novel junction at slightly different
    coordinates. This function merges such junctions to a consensus position.

    Args:
        novel_junctions: List of (donor, acceptor) novel junctions
        tolerance: Maximum bp distance for merging

    Returns:
        List of deduplicated (donor, acceptor) junctions
    """
    if not novel_junctions:
        return []

    # Sort by donor position
    sorted_junctions = sorted(novel_junctions, key=lambda x: (x[0], x[1]))

    # Cluster nearby junctions
    clusters = []
    current_cluster = [sorted_junctions[0]]

    for junction in sorted_junctions[1:]:
        # Check if junction is close to current cluster
        cluster_donor = np.mean([j[0] for j in current_cluster])
        cluster_acceptor = np.mean([j[1] for j in current_cluster])

        if (abs(junction[0] - cluster_donor) <= tolerance and
            abs(junction[1] - cluster_acceptor) <= tolerance):
            current_cluster.append(junction)
        else:
            clusters.append(current_cluster)
            current_cluster = [junction]

    clusters.append(current_cluster)

    # Get consensus position for each cluster
    deduplicated = []
    for cluster in clusters:
        # Use median for robustness
        consensus_donor = int(np.median([j[0] for j in cluster]))
        consensus_acceptor = int(np.median([j[1] for j in cluster]))
        deduplicated.append((consensus_donor, consensus_acceptor))

    return deduplicated


# =============================================================================
# Junction Signature Analysis
# =============================================================================

@dataclass
class JunctionSignature:
    """
    A unique junction pattern (set of junctions) representing an isoform.

    Attributes:
        junctions: Sorted tuple of (donor, acceptor) junctions
        chrom: Chromosome
        strand: Strand
        n_reads: Number of reads with this exact pattern
        read_ids: List of supporting read IDs
    """
    junctions: Tuple[Tuple[int, int], ...]
    chrom: str
    strand: str
    n_reads: int = 0
    read_ids: List[str] = None

    def __post_init__(self):
        if self.read_ids is None:
            self.read_ids = []

    @property
    def n_junctions(self) -> int:
        return len(self.junctions)

    @property
    def signature_str(self) -> str:
        """String representation of junction pattern."""
        if not self.junctions:
            return "intronless"
        return ";".join(f"{d}-{a}" for d, a in self.junctions)


def extract_junction_signatures(
    records: List[UnifiedReadRecord],
    correct_to_known: bool = True,
    known_junctions: Optional[Set[Tuple[int, int]]] = None,
    tolerance: int = DEFAULT_SPLICE_SITE_TOLERANCE,
) -> Dict[Tuple, JunctionSignature]:
    """
    Extract unique junction signatures (isoform patterns) from reads.

    Groups reads by their complete junction pattern after optional
    coordinate correction.

    Args:
        records: List of UnifiedReadRecord objects
        correct_to_known: If True, correct junctions to known before grouping
        known_junctions: Set of known junctions for correction
        tolerance: Tolerance for junction matching

    Returns:
        Dictionary mapping (chrom, strand, junction_tuple) to JunctionSignature
    """
    # Optionally correct junctions first
    if correct_to_known and known_junctions:
        records, _ = correct_all_junction_coordinates(
            records, known_junctions, tolerance
        )

    # Group by junction signature
    signatures = defaultdict(lambda: {'reads': [], 'n': 0})

    for record in records:
        # Create sorted tuple of junctions
        junction_tuple = tuple(sorted(record.junctions))
        key = (record.chrom, record.strand, junction_tuple)

        signatures[key]['reads'].append(record.read_id)
        signatures[key]['n'] += record.count

    # Convert to JunctionSignature objects
    result = {}
    for (chrom, strand, junction_tuple), data in signatures.items():
        result[(chrom, strand, junction_tuple)] = JunctionSignature(
            junctions=junction_tuple,
            chrom=chrom,
            strand=strand,
            n_reads=data['n'],
            read_ids=data['reads'],
        )

    return result


# =============================================================================
# Junction Index Building
# =============================================================================

def build_junction_index_from_records(
    records: List[UnifiedReadRecord],
    min_reads: int = 1,
) -> Dict[str, Set[Tuple[int, int]]]:
    """
    Build a junction index from read records.

    Collects all unique junctions per chromosome, optionally filtering
    by minimum read support.

    Args:
        records: List of UnifiedReadRecord objects
        min_reads: Minimum reads required to include junction

    Returns:
        Dictionary mapping chrom to set of (donor, acceptor) junctions
    """
    # Count junctions per chrom
    junction_counts = defaultdict(lambda: defaultdict(int))

    for record in records:
        for junction in record.junctions:
            junction_counts[record.chrom][junction] += record.count

    # Build index with filtering
    index = {}
    for chrom, junctions in junction_counts.items():
        index[chrom] = {
            junction for junction, count in junctions.items()
            if count >= min_reads
        }

    return index


def build_known_junction_index(
    annotation_path: str,
    file_format: str = 'auto',
) -> Dict[str, Set[Tuple[int, int]]]:
    """
    Build junction index from annotation file.

    Args:
        annotation_path: Path to GFF/GTF or BED file
        file_format: 'gff', 'bed', or 'auto' (detect from extension)

    Returns:
        Dictionary mapping chrom to set of (donor, acceptor) junctions
    """
    # Auto-detect format
    if file_format == 'auto':
        if annotation_path.endswith(('.gff', '.gff3', '.gtf')):
            file_format = 'gff'
        elif annotation_path.endswith('.bed'):
            file_format = 'bed'
        else:
            raise ValueError(f"Cannot auto-detect format for: {annotation_path}")

    index = defaultdict(set)

    if file_format in ('gff', 'gtf'):
        # Parse GFF/GTF for exon features, derive junctions
        current_transcript = None
        exons = []

        with open(annotation_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                chrom = parts[0]
                feature = parts[2]
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])
                strand = parts[6]
                attrs = parts[8]

                if feature == 'exon':
                    # Extract transcript ID
                    transcript_id = None
                    for attr in attrs.split(';'):
                        attr = attr.strip()
                        if attr.startswith('transcript_id'):
                            transcript_id = attr.split('"')[1] if '"' in attr else attr.split('=')[1]
                            break
                        elif attr.startswith('Parent='):
                            transcript_id = attr.split('=')[1]
                            break

                    if transcript_id != current_transcript:
                        # Process previous transcript
                        if current_transcript and len(exons) > 1:
                            exons.sort()
                            for i in range(len(exons) - 1):
                                # Junction = end of exon i to start of exon i+1
                                donor = exons[i][1]
                                acceptor = exons[i + 1][0]
                                index[exons[i][2]].add((donor, acceptor))

                        current_transcript = transcript_id
                        exons = []

                    exons.append((start, end, chrom))

            # Process last transcript
            if current_transcript and len(exons) > 1:
                exons.sort()
                for i in range(len(exons) - 1):
                    donor = exons[i][1]
                    acceptor = exons[i + 1][0]
                    index[exons[i][2]].add((donor, acceptor))

    elif file_format == 'bed':
        # Parse BED file (direct junction coordinates)
        with open(annotation_path, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('track'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue

                chrom = parts[0]
                start = int(parts[1])  # Already 0-based
                end = int(parts[2])

                index[chrom].add((start, end))

    return dict(index)


# =============================================================================
# Summary and Export
# =============================================================================

def summarize_junction_analysis(
    records: List[UnifiedReadRecord],
    known_junctions: Set[Tuple[int, int]],
    tolerance: int = DEFAULT_SPLICE_SITE_TOLERANCE,
) -> Dict:
    """
    Generate summary statistics for junction analysis.

    Args:
        records: List of UnifiedReadRecord objects
        known_junctions: Set of known junctions
        tolerance: Tolerance used for matching

    Returns:
        Dictionary with summary statistics
    """
    stats = {
        'tolerance_bp': tolerance,
        'total_reads': len(records),
        'reads_with_junctions': 0,
        'unique_junctions_raw': set(),
        'unique_junctions_corrected': set(),
        'junctions_matched_exact': 0,
        'junctions_matched_tolerance': 0,
        'junctions_novel': 0,
    }

    for record in records:
        if record.junctions:
            stats['reads_with_junctions'] += 1

            for junction in record.junctions:
                stats['unique_junctions_raw'].add(junction)

                # Check matching
                if junction in known_junctions:
                    stats['junctions_matched_exact'] += 1
                    stats['unique_junctions_corrected'].add(junction)
                else:
                    matched = match_junction_with_tolerance(
                        junction, known_junctions, tolerance
                    )
                    if matched:
                        stats['junctions_matched_tolerance'] += 1
                        stats['unique_junctions_corrected'].add(matched)
                    else:
                        stats['junctions_novel'] += 1
                        stats['unique_junctions_corrected'].add(junction)

    # Convert sets to counts
    stats['unique_junctions_raw'] = len(stats['unique_junctions_raw'])
    stats['unique_junctions_corrected'] = len(stats['unique_junctions_corrected'])

    return stats


if __name__ == '__main__':
    # Quick test
    print("Testing junction_analysis module...")

    # Create known junctions
    known = {(1200, 1400), (2000, 2200), (3000, 3300)}

    # Test junction with slight offset
    read_junction = (1198, 1402)  # 2bp off on each end

    matched = match_junction_with_tolerance(read_junction, known, tolerance=5)
    print(f"\nMatching {read_junction} with tolerance=5:")
    print(f"  Matched to: {matched}")

    # Test with tolerance=1 (should not match)
    matched_strict = match_junction_with_tolerance(read_junction, known, tolerance=1)
    print(f"\nMatching {read_junction} with tolerance=1:")
    print(f"  Matched to: {matched_strict}")

    # Test deduplication
    novel_junctions = [
        (5000, 5200),
        (5002, 5198),  # Close to first
        (5001, 5201),  # Close to first
        (6000, 6300),  # Different junction
    ]

    deduplicated = deduplicate_novel_junctions(novel_junctions, tolerance=10)
    print(f"\nDeduplicating {len(novel_junctions)} junctions:")
    print(f"  Result: {deduplicated}")

    print("\njunction_analysis module ready!")
