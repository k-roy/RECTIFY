"""
Multi-Aligner Consensus Module for RECTIFY.

Compares alignments from multiple aligners and selects the best one per read
based on junction quality, soft-clip rescue, and false junction detection.

Scoring priorities:
1. Prefer alignments that splice through junctions vs soft-clipping (5' rescue)
2. Prefer canonical splice sites (GT-AG)
3. Prefer junctions supported by multiple aligners
4. Remove 3' false junctions from poly(A) artifacts

Author: Kevin R. Roy
"""

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from pathlib import Path

import pysam
import numpy as np

logger = logging.getLogger(__name__)


# Canonical splice site dinucleotides
CANONICAL_5SS = {'GT', 'GC'}  # 5' splice site (donor)
CANONICAL_3SS = {'AG'}  # 3' splice site (acceptor)


@dataclass
class AlignmentInfo:
    """Stores key info about a read's alignment from one aligner."""
    read_id: str
    aligner: str
    chrom: str
    strand: str
    reference_start: int
    reference_end: int
    cigar_string: str
    mapq: int

    # Junction information
    junctions: List[Tuple[int, int]] = field(default_factory=list)  # (start, end) pairs

    # Soft-clip info (5' and 3')
    five_prime_softclip: int = 0
    three_prime_softclip: int = 0

    # Quality scores
    junction_score: float = 0.0
    canonical_count: int = 0
    non_canonical_count: int = 0

    # Flags
    has_false_3prime_junction: bool = False
    is_best: bool = False


@dataclass
class ConsensusResult:
    """Result of consensus selection for a read."""
    read_id: str
    best_aligner: str
    best_alignment: Optional[AlignmentInfo]
    aligners_compared: List[str]

    # Consensus metrics
    n_aligners_agree: int = 0  # Number of aligners with same junctions
    confidence: str = ""  # 'high', 'medium', 'low'

    # Rescue info
    was_5prime_rescued: bool = False  # Did we pick an alignment that spliced vs clipped
    false_junction_removed: bool = False

    # All alignments for this read (for debugging)
    all_alignments: Dict[str, AlignmentInfo] = field(default_factory=dict)


def extract_junctions_from_cigar(read: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    """
    Extract splice junctions from CIGAR N operations.

    Returns list of (intron_start, intron_end) in 0-based coords.
    intron_end is exclusive.
    """
    if not read.cigartuples:
        return []

    junctions = []
    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op == 3:  # N = skipped region (intron)
            junctions.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op in (0, 2, 7, 8):  # M, D, =, X consume reference
            ref_pos += length

    return junctions


def get_softclip_lengths(read: pysam.AlignedSegment) -> Tuple[int, int]:
    """
    Get 5' and 3' soft-clip lengths for a read.

    For plus strand: 5' is left, 3' is right
    For minus strand: 5' is right, 3' is left (in aligned orientation)

    Returns: (five_prime_clip, three_prime_clip)
    """
    if not read.cigartuples:
        return (0, 0)

    left_clip = 0
    right_clip = 0

    # Left soft-clip
    if read.cigartuples[0][0] == 4:  # S = soft-clip
        left_clip = read.cigartuples[0][1]

    # Right soft-clip
    if read.cigartuples[-1][0] == 4:
        right_clip = read.cigartuples[-1][1]

    # Adjust for strand
    if read.is_reverse:
        return (right_clip, left_clip)  # 5' is on right for minus strand
    else:
        return (left_clip, right_clip)  # 5' is on left for plus strand


def check_canonical_splice_sites(
    junctions: List[Tuple[int, int]],
    chrom: str,
    genome: Dict[str, str],
) -> Tuple[int, int]:
    """
    Count canonical vs non-canonical splice sites for junctions.

    Returns: (canonical_count, non_canonical_count)
    """
    if chrom not in genome:
        return (0, 0)

    seq = genome[chrom]
    canonical = 0
    non_canonical = 0

    for start, end in junctions:
        if start < 0 or end > len(seq):
            continue

        # 5'SS: first 2 bases of intron
        five_ss = seq[start:start + 2].upper()
        # 3'SS: last 2 bases of intron
        three_ss = seq[end - 2:end].upper()

        if five_ss in CANONICAL_5SS and three_ss in CANONICAL_3SS:
            canonical += 1
        else:
            non_canonical += 1

    return (canonical, non_canonical)


def detect_false_3prime_junction(
    read: pysam.AlignedSegment,
    junctions: List[Tuple[int, int]],
    genome: Dict[str, str],
    min_polya_for_false: int = 3,
) -> bool:
    """
    Detect if the 3'-most junction is a false junction from poly(A) artifacts.

    False junction pattern:
    - Junction is at the 3' end of the read
    - Region before junction ends with A's
    - Region after junction starts with A's
    - Read has no substantial sequence after the junction (just poly(A))

    Returns True if the 3'-most junction appears to be a poly(A) artifact.
    """
    if not junctions:
        return False

    chrom = read.reference_name
    if chrom not in genome:
        return False

    seq = genome[chrom]
    strand = '-' if read.is_reverse else '+'

    # Get 3'-most junction (last for + strand, first for - strand)
    if strand == '+':
        junc_start, junc_end = junctions[-1]
    else:
        junc_start, junc_end = junctions[0]

    # Check if this junction is near the 3' end of the read
    read_3prime = read.reference_end if strand == '+' else read.reference_start

    if strand == '+':
        # For plus strand, junction should be near reference_end
        dist_to_3prime = read.reference_end - junc_end
    else:
        # For minus strand, junction should be near reference_start
        dist_to_3prime = junc_start - read.reference_start

    # If junction is not near 3' end, it's likely real
    if dist_to_3prime > 50:
        return False

    # Check for A-tract pattern
    # Before junction (exon side)
    if strand == '+':
        before_seq = seq[max(0, junc_start - 10):junc_start]
    else:
        before_seq = seq[junc_end:min(len(seq), junc_end + 10)]

    # Count trailing A's before junction
    trailing_a = 0
    for base in reversed(before_seq):
        if base.upper() == 'A':
            trailing_a += 1
        else:
            break

    # After junction
    if strand == '+':
        after_seq = seq[junc_end:min(len(seq), junc_end + 10)]
    else:
        after_seq = seq[max(0, junc_start - 10):junc_start]

    # Count leading A's after junction
    leading_a = 0
    for base in after_seq:
        if base.upper() == 'A':
            leading_a += 1
        else:
            break

    # Pattern for false junction: A's on both sides
    if trailing_a >= min_polya_for_false and leading_a >= min_polya_for_false:
        logger.debug(f"Detected false 3' junction at {chrom}:{junc_start}-{junc_end}")
        return True

    return False


def score_alignment(
    alignment: AlignmentInfo,
    genome: Dict[str, str],
    annotated_junctions: Optional[Set[Tuple[str, int, int]]] = None,
) -> float:
    """
    Score an alignment based on junction quality.

    Scoring factors:
    - +10 per canonical junction (GT-AG)
    - +5 per annotated junction
    - -5 per non-canonical junction
    - -2 per soft-clipped base at 5' end (prefer spliced alignments)

    Note: False 3' junctions from poly(A) artifacts are handled by the
    walk back correction step, which eats through aligned A's and discards
    spurious N operations to find the true CPA site.
    """
    score = 0.0

    # Canonical junctions bonus
    score += alignment.canonical_count * 10

    # Non-canonical junction penalty
    score -= alignment.non_canonical_count * 5

    # 5' soft-clip penalty (prefer alignments that splice through)
    score -= alignment.five_prime_softclip * 2

    # Annotated junction bonus
    if annotated_junctions and alignment.junctions:
        for junc in alignment.junctions:
            if (alignment.chrom, junc[0], junc[1]) in annotated_junctions:
                score += 5

    alignment.junction_score = score
    return score


def extract_alignment_info(
    read: pysam.AlignedSegment,
    aligner: str,
    genome: Dict[str, str],
) -> AlignmentInfo:
    """Extract alignment info from a pysam read."""

    junctions = extract_junctions_from_cigar(read)
    five_clip, three_clip = get_softclip_lengths(read)

    chrom = read.reference_name
    canonical, non_canonical = check_canonical_splice_sites(junctions, chrom, genome)

    # Note: False 3' junction detection is available via detect_false_3prime_junction()
    # but not used in scoring since walk back correction handles this case.

    return AlignmentInfo(
        read_id=read.query_name,
        aligner=aligner,
        chrom=chrom,
        strand='-' if read.is_reverse else '+',
        reference_start=read.reference_start,
        reference_end=read.reference_end,
        cigar_string=read.cigarstring or "",
        mapq=read.mapping_quality,
        junctions=junctions,
        five_prime_softclip=five_clip,
        three_prime_softclip=three_clip,
        canonical_count=canonical,
        non_canonical_count=non_canonical,
        has_false_3prime_junction=False,  # Not used; walk back handles this
    )


def select_best_alignment(
    alignments: Dict[str, AlignmentInfo],
    genome: Dict[str, str],
    annotated_junctions: Optional[Set[Tuple[str, int, int]]] = None,
) -> ConsensusResult:
    """
    Select the best alignment from multiple aligners for a single read.

    Args:
        alignments: Dict mapping aligner name to AlignmentInfo
        genome: Dict mapping chrom to sequence
        annotated_junctions: Optional set of (chrom, start, end) for annotated junctions

    Returns:
        ConsensusResult with best alignment selected
    """
    if not alignments:
        return ConsensusResult(
            read_id="",
            best_aligner="none",
            best_alignment=None,
            aligners_compared=[],
        )

    read_id = list(alignments.values())[0].read_id

    # Score all alignments
    for aligner, alignment in alignments.items():
        score_alignment(alignment, genome, annotated_junctions)

    # Select best by score
    best_aligner = max(alignments.keys(), key=lambda a: alignments[a].junction_score)
    best_alignment = alignments[best_aligner]
    best_alignment.is_best = True

    # Check for 5' rescue
    # Did we pick an alignment that spliced through vs one that soft-clipped?
    was_rescued = False
    if len(alignments) > 1:
        min_5clip = min(a.five_prime_softclip for a in alignments.values())
        max_5clip = max(a.five_prime_softclip for a in alignments.values())
        if max_5clip > min_5clip and best_alignment.five_prime_softclip == min_5clip:
            was_rescued = True

    # Count junction agreement across aligners
    junction_sets = {}
    for aligner, alignment in alignments.items():
        junc_key = tuple(sorted(alignment.junctions))
        junction_sets[aligner] = junc_key

    # Count how many aligners agree with the best alignment's junctions
    best_juncs = junction_sets[best_aligner]
    n_agree = sum(1 for jset in junction_sets.values() if jset == best_juncs)

    # Confidence based on agreement
    if n_agree == len(alignments) and len(alignments) >= 2:
        confidence = 'high'
    elif n_agree >= 2:
        confidence = 'medium'
    else:
        confidence = 'low'

    return ConsensusResult(
        read_id=read_id,
        best_aligner=best_aligner,
        best_alignment=best_alignment,
        aligners_compared=list(alignments.keys()),
        n_aligners_agree=n_agree,
        confidence=confidence,
        was_5prime_rescued=was_rescued,
        false_junction_removed=False,  # Not tracked; walk back correction handles this
        all_alignments=alignments,
    )


def run_consensus_selection(
    bam_paths: Dict[str, str],
    genome: Dict[str, str],
    output_bam: str,
    annotated_junctions: Optional[Set[Tuple[str, int, int]]] = None,
    write_all_to_tag: bool = True,
) -> Dict[str, int]:
    """
    Run consensus selection across multiple BAM files.

    Args:
        bam_paths: Dict mapping aligner name to BAM path
        genome: Dict mapping chrom to sequence
        output_bam: Output path for consensus BAM
        annotated_junctions: Optional set of annotated junctions
        write_all_to_tag: If True, write all aligner info to BAM tags

    Returns:
        Summary statistics dict
    """
    # Open all BAM files
    bams = {}
    for aligner, path in bam_paths.items():
        bams[aligner] = pysam.AlignmentFile(path, 'rb')

    # Get header from first BAM (should be identical across aligners)
    first_bam = list(bams.values())[0]
    header = first_bam.header.to_dict()

    # Add program group for RECTIFY consensus
    if 'PG' not in header:
        header['PG'] = []
    header['PG'].append({
        'ID': 'RECTIFY',
        'PN': 'RECTIFY',
        'VN': '2.0',
        'CL': f'consensus selection from {",".join(bam_paths.keys())}',
    })

    # Open output BAM
    out_bam = pysam.AlignmentFile(output_bam, 'wb', header=pysam.AlignmentHeader.from_dict(header))

    # Build read index for all BAMs
    # Note: This assumes BAMs are sorted and indexed
    logger.info("Building read indices...")
    read_alignments: Dict[str, Dict[str, pysam.AlignedSegment]] = defaultdict(dict)

    for aligner, bam in bams.items():
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            read_alignments[read.query_name][aligner] = read

    logger.info(f"Found {len(read_alignments)} unique reads across {len(bams)} aligners")

    # Process each read
    stats = {
        'total_reads': 0,
        'consensus_high': 0,
        'consensus_medium': 0,
        'consensus_low': 0,
        '5prime_rescued': 0,
        'by_aligner': defaultdict(int),
    }

    for read_id, aligner_reads in read_alignments.items():
        stats['total_reads'] += 1

        # Extract alignment info
        alignments = {}
        for aligner, read in aligner_reads.items():
            alignments[aligner] = extract_alignment_info(read, aligner, genome)

        # Select best
        result = select_best_alignment(alignments, genome, annotated_junctions)

        # Update stats
        if result.confidence == 'high':
            stats['consensus_high'] += 1
        elif result.confidence == 'medium':
            stats['consensus_medium'] += 1
        else:
            stats['consensus_low'] += 1

        if result.was_5prime_rescued:
            stats['5prime_rescued'] += 1

        # Note: False 3' junctions are handled by walk back correction,
        # not during consensus selection.

        stats['by_aligner'][result.best_aligner] += 1

        # Write best alignment to output
        if result.best_aligner in aligner_reads:
            best_read = aligner_reads[result.best_aligner]

            # Add consensus tags
            best_read.set_tag('XA', result.best_aligner)  # Best aligner
            best_read.set_tag('XC', result.confidence)  # Confidence
            best_read.set_tag('XN', result.n_aligners_agree)  # N aligners agree

            if result.was_5prime_rescued:
                best_read.set_tag('XR', 1)  # Was 5' rescued

            if result.false_junction_removed:
                best_read.set_tag('XF', 1)  # Had false junction

            out_bam.write(best_read)

    # Close files
    for bam in bams.values():
        bam.close()
    out_bam.close()

    # Index output
    pysam.index(output_bam)

    # Log summary
    logger.info(f"\nConsensus selection complete:")
    logger.info(f"  Total reads: {stats['total_reads']}")
    logger.info(f"  High confidence: {stats['consensus_high']}")
    logger.info(f"  Medium confidence: {stats['consensus_medium']}")
    logger.info(f"  Low confidence: {stats['consensus_low']}")
    logger.info(f"  5' rescued: {stats['5prime_rescued']}")
    logger.info(f"  By aligner: {dict(stats['by_aligner'])}")

    return stats


def load_annotated_junctions(annotation_path: str) -> Set[Tuple[str, int, int]]:
    """
    Load annotated junctions from GFF/GTF file.

    Returns set of (chrom, intron_start, intron_end) tuples.
    """
    junctions = set()

    with open(annotation_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2].lower()

            # Look for intron features
            if feature_type == 'intron':
                chrom = parts[0]
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])  # Already exclusive in GFF end
                junctions.add((chrom, start, end))

    logger.info(f"Loaded {len(junctions)} annotated junctions from {annotation_path}")
    return junctions
