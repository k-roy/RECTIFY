"""
NET-seq BAM processing for RECTIFY.

Extracts 3' ends from NET-seq BAM files with:
- Strand-aware soft-clip oligo(A) detection
- CIGAR-based alignment parsing
- Chromosome name standardization

Unlike nanopore processing, NET-seq reads are short (~40-76bp) and represent
the 3' end of nascent RNA. The aligned 3' end of the read IS the 3' end
of the nascent transcript.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Generator
import re

import pysam

from ..config import NCBI_TO_CHROM
from ..utils.alignment import extract_soft_clips, parse_cigar
from .unified_record import UnifiedReadRecord
from .exclusion_regions import ExclusionRegionDetector


# Extended NCBI mapping for full contig names with metadata
def extract_ncbi_accession(contig_name: str) -> str:
    """
    Extract NCBI accession from full contig name.

    Handles formats like:
    - 'ref|NC_001133| [org=Saccharomyces cerevisiae] ...'
    - 'ref|NC_001133|'
    - 'NC_001133'

    Returns:
        NCBI accession in 'ref|NC_XXXXXX|' format, or original name if not recognized
    """
    # Match ref|NC_XXXXXX| pattern
    match = re.match(r'(ref\|NC_\d+\|)', contig_name)
    if match:
        return match.group(1)

    # Match bare NC_XXXXXX
    match = re.match(r'(NC_\d+)', contig_name)
    if match:
        return f"ref|{match.group(1)}|"

    return contig_name


def standardize_netseq_chrom(contig_name: str) -> str:
    """
    Standardize NET-seq BAM chromosome name to chrI, chrII, etc. format.

    Args:
        contig_name: Full contig name from BAM header

    Returns:
        Standardized chromosome name (chrI, chrII, etc.)
    """
    # Extract NCBI accession
    ncbi_acc = extract_ncbi_accession(contig_name)

    # Look up in NCBI_TO_CHROM mapping
    if ncbi_acc in NCBI_TO_CHROM:
        return NCBI_TO_CHROM[ncbi_acc]

    # Try without the pipes
    ncbi_bare = ncbi_acc.replace('ref|', '').replace('|', '')
    for key, value in NCBI_TO_CHROM.items():
        if ncbi_bare in key:
            return value

    # Return as-is if unrecognized
    return contig_name


def get_netseq_3prime_position(read: pysam.AlignedSegment) -> Tuple[int, str]:
    """
    Get 3' end position and strand from NET-seq read.

    NET-seq captures nascent RNA 3' ends. The read represents the
    terminal fragment, so:
    - Forward strand (FLAG 0): 3' end = right side = reference_end - 1
    - Reverse strand (FLAG 16): 3' end = left side = reference_start

    Args:
        read: pysam AlignedSegment

    Returns:
        Tuple of (position, strand) where position is 0-based
    """
    strand = '-' if read.is_reverse else '+'

    if strand == '+':
        # Plus strand: 3' end is rightmost aligned base
        position = read.reference_end - 1
    else:
        # Minus strand: 3' end is leftmost aligned base
        position = read.reference_start

    return position, strand


def detect_oligo_a_in_softclip(
    read: pysam.AlignedSegment,
    strand: str,
    min_a_fraction: float = 0.8,
) -> Dict:
    """
    Detect oligo(A) content in 3' soft-clipped sequence (strand-aware).

    NET-seq reads can have oligo(A) tails in soft-clips due to nascent
    transcript adenylation. This is strand-aware:
    - Plus strand: 3' end is RIGHT side, oligo(A) in right soft-clip
    - Minus strand: 3' end is LEFT side, oligo(A) in left soft-clip
      (soft-clip sequence contains T's which = A's in RNA)

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'
        min_a_fraction: Minimum A (or T for minus) content to classify as oligo(A)

    Returns:
        Dict with:
            - soft_clip_a_length: Number of A's (RNA orientation) in 3' soft-clip
            - three_prime_soft_clip_seq: Full soft-clip sequence (genomic orientation)
            - three_prime_soft_clip_length: Total soft-clip length
            - has_oligo_a: Boolean indicating if clip is A-rich
    """
    clips = extract_soft_clips(read)

    # Find the soft-clip at the biological 3' end
    if strand == '+':
        clip_side = 'right'  # 3' is right for plus strand
    else:
        clip_side = 'left'   # 3' is left for minus strand

    three_prime_clip = next((c for c in clips if c['side'] == clip_side), None)

    if three_prime_clip is None or not three_prime_clip['seq']:
        return {
            'soft_clip_a_length': 0,
            'three_prime_soft_clip_seq': '',
            'three_prime_soft_clip_length': 0,
            'has_oligo_a': False,
        }

    seq = three_prime_clip['seq']
    clip_length = len(seq)

    # Count A's in RNA orientation
    # For minus strand, T in genomic = A in RNA
    if strand == '-':
        a_count = seq.upper().count('T')
    else:
        a_count = seq.upper().count('A')

    a_fraction = a_count / clip_length if clip_length > 0 else 0
    is_a_rich = a_fraction >= min_a_fraction

    return {
        'soft_clip_a_length': a_count if is_a_rich else 0,
        'three_prime_soft_clip_seq': seq,
        'three_prime_soft_clip_length': clip_length,
        'has_oligo_a': is_a_rich,
    }


def get_5prime_softclip_info(
    read: pysam.AlignedSegment,
    strand: str,
) -> Dict:
    """
    Get 5' soft-clip information (for completeness).

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'

    Returns:
        Dict with 5' soft-clip length and sequence
    """
    clips = extract_soft_clips(read)

    # 5' end is opposite of 3' end
    if strand == '+':
        clip_side = 'left'   # 5' is left for plus strand
    else:
        clip_side = 'right'  # 5' is right for minus strand

    five_prime_clip = next((c for c in clips if c['side'] == clip_side), None)

    if five_prime_clip is None:
        return {
            'five_prime_soft_clip_length': 0,
            'five_prime_soft_clip_seq': '',
        }

    return {
        'five_prime_soft_clip_length': five_prime_clip['length'],
        'five_prime_soft_clip_seq': five_prime_clip['seq'] or '',
    }


def process_netseq_read(
    read: pysam.AlignedSegment,
    chrom_std: str,
    min_a_fraction: float = 0.8,
) -> UnifiedReadRecord:
    """
    Process single NET-seq read into UnifiedReadRecord.

    Creates a record with:
    - three_prime_raw: Raw 3' position from alignment
    - three_prime_corrected: Same as raw (deconvolution is position-level)
    - soft_clip_a_length: Oligo(A) in 3' soft-clip
    - All other standard fields

    Args:
        read: pysam AlignedSegment
        chrom_std: Standardized chromosome name
        min_a_fraction: Minimum A fraction for oligo(A) detection

    Returns:
        UnifiedReadRecord with all fields populated
    """
    # Get 3' position and strand
    three_prime_pos, strand = get_netseq_3prime_position(read)

    # Get 5' position
    if strand == '+':
        five_prime_pos = read.reference_start
    else:
        five_prime_pos = read.reference_end - 1

    # Detect oligo(A) in 3' soft-clip
    oligo_a_info = detect_oligo_a_in_softclip(read, strand, min_a_fraction)

    # Get 5' soft-clip info
    five_prime_info = get_5prime_softclip_info(read, strand)

    # Parse CIGAR
    cigar_str = parse_cigar(read.cigartuples) if read.cigartuples else ''

    # Create UnifiedReadRecord
    record = UnifiedReadRecord(
        read_id=read.query_name,
        chrom=chrom_std,
        strand=strand,

        # 5' end positions
        five_prime_raw=five_prime_pos,
        five_prime_corrected=five_prime_pos,  # No correction for NET-seq
        first_exon_start=None,
        starts_in_intron=False,

        # 3' end positions
        three_prime_raw=three_prime_pos,
        three_prime_corrected=three_prime_pos,  # Updated by deconvolution later

        # Alignment span
        alignment_start=read.reference_start,
        alignment_end=read.reference_end,

        # Splice junctions (NET-seq reads are short, rarely spliced)
        junctions=[],
        n_junctions=0,
        junctions_filtered=0,

        # Soft clips
        five_prime_soft_clip_length=five_prime_info['five_prime_soft_clip_length'],
        three_prime_soft_clip_length=oligo_a_info['three_prime_soft_clip_length'],
        five_prime_soft_clip_seq=five_prime_info['five_prime_soft_clip_seq'],
        three_prime_soft_clip_seq=oligo_a_info['three_prime_soft_clip_seq'],

        # Poly(A) / Oligo(A) information
        polya_length=oligo_a_info['soft_clip_a_length'],
        aligned_a_length=0,  # Not computed for NET-seq
        soft_clip_a_length=oligo_a_info['soft_clip_a_length'],

        # Alignment quality
        mapq=read.mapping_quality,
        cigar_summary=cigar_str,
        alignment_identity=0.0,  # Not computed for NET-seq

        # Consensus info (not applicable for NET-seq)
        best_aligner='bbmap',
        consensus_confidence='',
        n_aligners_support=0,

        # Count
        count=1,
    )

    return record


def process_netseq_bam(
    bam_path: Path,
    exclusion_detector: Optional[ExclusionRegionDetector] = None,
    min_mapq: int = 0,
    min_a_fraction: float = 0.8,
    max_reads: Optional[int] = None,
    show_progress: bool = True,
    progress_interval: int = 1_000_000,
) -> Generator[UnifiedReadRecord, None, None]:
    """
    Stream process NET-seq BAM file.

    Yields UnifiedReadRecord for each valid read, excluding:
    - Unmapped reads
    - Secondary/supplementary alignments
    - Reads in exclusion regions (if detector provided)
    - Reads below min_mapq

    Args:
        bam_path: Path to BAM file
        exclusion_detector: Optional ExclusionRegionDetector
        min_mapq: Minimum mapping quality (default 0)
        min_a_fraction: Minimum A fraction for oligo(A) detection
        max_reads: Limit for testing (None = process all)
        show_progress: Print progress messages
        progress_interval: Print progress every N reads

    Yields:
        UnifiedReadRecord for each valid read
    """
    bam_path = Path(bam_path)

    if show_progress:
        print(f"Processing NET-seq BAM: {bam_path.name}")

    # Build chromosome name mapping from BAM header
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    chrom_map = {}
    for sq in bam.header['SQ']:
        contig_name = sq['SN']
        chrom_map[contig_name] = standardize_netseq_chrom(contig_name)

    if show_progress:
        print(f"  Chromosome mapping: {len(chrom_map)} contigs")
        # Show a few examples
        for i, (orig, std) in enumerate(chrom_map.items()):
            if i < 3:
                print(f"    {orig[:50]}... -> {std}")

    # Process reads
    n_total = 0
    n_skipped_unmapped = 0
    n_skipped_secondary = 0
    n_skipped_mapq = 0
    n_skipped_excluded = 0
    n_yielded = 0

    for read in bam:
        n_total += 1

        # Progress reporting
        if show_progress and n_total % progress_interval == 0:
            print(f"  Processed {n_total:,} reads, yielded {n_yielded:,}...")

        # Skip invalid reads
        if read.is_unmapped:
            n_skipped_unmapped += 1
            continue
        if read.is_secondary or read.is_supplementary:
            n_skipped_secondary += 1
            continue
        if read.mapping_quality < min_mapq:
            n_skipped_mapq += 1
            continue

        # Get standardized chromosome name
        contig_name = read.reference_name
        chrom_std = chrom_map.get(contig_name, contig_name)

        # Get 3' position for exclusion check
        three_prime_pos, strand = get_netseq_3prime_position(read)

        # Check exclusion regions
        if exclusion_detector and exclusion_detector.is_excluded(chrom_std, three_prime_pos):
            n_skipped_excluded += 1
            continue

        # Process read
        record = process_netseq_read(read, chrom_std, min_a_fraction)
        n_yielded += 1

        yield record

        # Check max_reads limit
        if max_reads and n_yielded >= max_reads:
            if show_progress:
                print(f"  Reached max_reads limit ({max_reads})")
            break

    bam.close()

    if show_progress:
        print(f"  Completed: {n_total:,} total reads")
        print(f"    Yielded: {n_yielded:,}")
        print(f"    Skipped unmapped: {n_skipped_unmapped:,}")
        print(f"    Skipped secondary/supp: {n_skipped_secondary:,}")
        print(f"    Skipped low MAPQ: {n_skipped_mapq:,}")
        print(f"    Skipped excluded regions: {n_skipped_excluded:,}")


def aggregate_positions(
    records: Generator[UnifiedReadRecord, None, None],
    position_field: str = 'three_prime_raw',
) -> Dict[Tuple[str, str, int], int]:
    """
    Aggregate records into position counts.

    Args:
        records: Generator of UnifiedReadRecord
        position_field: Field to use for position ('three_prime_raw' or 'three_prime_corrected')

    Returns:
        Dict mapping (chrom, strand, position) -> count
    """
    counts: Dict[Tuple[str, str, int], int] = {}

    for record in records:
        position = getattr(record, position_field)
        key = (record.chrom, record.strand, position)
        counts[key] = counts.get(key, 0) + record.count

    return counts
