"""
Multi-Aligner Consensus Module for RECTIFY.

Compares alignments from multiple aligners and selects the best one per read
based on junction quality, soft-clip rescue, and false junction detection.

Scoring priorities:
1. Prefer alignments that splice through junctions vs soft-clipping (5' rescue)
2. Prefer alignments whose 3' end lands outside a downstream A-tract (3' end quality)
3. Prefer junctions supported by multiple aligners
4. Tiebreaker: prefer aligner whose corrected 3' position agrees with majority
5. Tiebreaker: canonical splice site motifs (GT/AG) and annotated junctions

Note: A-tract 3' correction is applied to each aligner pre-scoring using genome
sequence only. Full indel correction (MD-tag dependent) is applied post-consensus
as a refinement step.

Author: Kevin R. Roy
"""

import logging
import os
import hashlib
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

    # 3' end A-tract correction (computed pre-consensus using genome sequence only;
    # full indel correction requiring MD tags is applied post-consensus)
    corrected_3prime: Optional[int] = None   # estimated true CPA position
    three_prime_atract_depth: int = 0        # A's downstream of raw 3' end (0 = clean landing)

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
    - -2 per soft-clipped base at 5' end (prefer spliced alignments)

    Neither canonical splice site motifs (GT/AG) nor annotated junction
    matches are scored here, to avoid biasing against novel junctions.
    Both are used only as tiebreakers in select_best_alignment().

    Note: False 3' junctions from poly(A) artifacts are handled by the
    walk back correction step, which eats through aligned A's and discards
    spurious N operations to find the true CPA site.
    """
    score = 0.0

    # NOTE: Canonical junction motifs (GT-AG) are deliberately not scored here.
    # They are used only as a tiebreaker in select_best_alignment() to avoid
    # biasing against novel non-canonical junctions.

    # 5' soft-clip penalty (prefer alignments that splice through)
    score -= alignment.five_prime_softclip * 2

    # 3' A-tract depth penalty (prefer alignments landing closer to the true CPA).
    # Each downstream A the aligner runs into costs 1 point — same scale as 5' penalty
    # per base, but capped at 10 to avoid overwhelming junction scoring.
    score -= min(alignment.three_prime_atract_depth, 10)

    # NOTE: Annotated junction matches are deliberately not scored here.
    # They are used only as a tiebreaker in select_best_alignment() to avoid
    # biasing against novel unannotated junctions.

    alignment.junction_score = score
    return score


def extract_alignment_info(
    read: pysam.AlignedSegment,
    aligner: str,
    genome: Dict[str, str],
) -> AlignmentInfo:
    """Extract alignment info from a pysam read.

    Computes corrected_3prime pre-consensus using A-tract detection (genome-only,
    no MD tags required). Full indel correction (MD-dependent) is applied
    post-consensus as a refinement step.
    """
    from .atract_detector import calculate_atract_ambiguity

    junctions = extract_junctions_from_cigar(read)
    five_clip, three_clip = get_softclip_lengths(read)

    chrom = read.reference_name
    strand = '-' if read.is_reverse else '+'
    canonical, non_canonical = check_canonical_splice_sites(junctions, chrom, genome)

    # Estimate corrected 3' end using A-tract ambiguity detection.
    # Raw 3' end: reference_end - 1 for + strand, reference_start for - strand.
    raw_3prime = (read.reference_end - 1) if strand == '+' else read.reference_start
    corrected_3prime = raw_3prime
    atract_depth = 0

    chrom_std = chrom
    if genome.get(chrom_std) is None:
        # Try standardized name
        from ..utils.genome import standardize_chrom_name
        chrom_std = standardize_chrom_name(chrom, genome) or chrom

    if genome.get(chrom_std) is not None:
        try:
            atract = calculate_atract_ambiguity(
                genome, chrom_std, raw_3prime, strand, downstream_bp=10
            )
            atract_depth = atract.get('downstream_a_count') or 0
            # Best-guess corrected position: ambiguity_min for +, ambiguity_max for -
            if strand == '+':
                corrected_3prime = atract.get('ambiguity_min', raw_3prime)
            else:
                corrected_3prime = atract.get('ambiguity_max', raw_3prime)
        except Exception:
            pass  # Non-fatal; raw position used

    return AlignmentInfo(
        read_id=read.query_name,
        aligner=aligner,
        chrom=chrom,
        strand=strand,
        reference_start=read.reference_start,
        reference_end=read.reference_end,
        cigar_string=read.cigarstring or "",
        mapq=read.mapping_quality,
        junctions=junctions,
        five_prime_softclip=five_clip,
        three_prime_softclip=three_clip,
        corrected_3prime=corrected_3prime,
        three_prime_atract_depth=atract_depth,
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

    # Select best by score, using annotation and canonical motifs as tiebreakers
    max_score = max(a.junction_score for a in alignments.values())
    tied_aligners = [name for name, a in alignments.items()
                     if a.junction_score == max_score]

    if len(tied_aligners) == 1:
        best_aligner = tied_aligners[0]
    else:
        # Tiebreaker 1: prefer alignment whose corrected 3' end agrees with majority
        all_corrected = [a.corrected_3prime for a in alignments.values()
                         if a.corrected_3prime is not None]
        def _count_3prime_agreement(aligner_name):
            pos = alignments[aligner_name].corrected_3prime
            return sum(1 for p in all_corrected if p == pos) if pos is not None else 0

        # Tiebreaker 2: prefer alignment with more annotated junctions
        # Tiebreaker 3: prefer alignment with more canonical splice sites (GT/AG)
        def _tiebreak_key(aligner_name):
            a = alignments[aligner_name]
            n_annotated = 0
            if annotated_junctions and a.junctions:
                n_annotated = sum(
                    1 for junc in a.junctions
                    if (a.chrom, junc[0], junc[1]) in annotated_junctions
                )
            return (_count_3prime_agreement(aligner_name), n_annotated, a.canonical_count)

        best_aligner = max(tied_aligners, key=_tiebreak_key)

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
    if n_agree == len(alignments):
        confidence = 'high'  # All aligners agree (including single-aligner case)
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



def _ensure_name_sorted(bam_path: str) -> str:
    """
    Ensure a BAM file is name-sorted. If not, create a name-sorted copy.

    Returns path to name-sorted BAM (may be same as input if already sorted).
    """
    bam = pysam.AlignmentFile(bam_path, 'rb')
    header = bam.header.to_dict()
    bam.close()

    sort_order = header.get('HD', {}).get('SO', 'unknown')
    if sort_order == 'queryname':
        logger.debug(f"BAM already name-sorted: {bam_path}")
        return bam_path

    sorted_path = bam_path.replace('.bam', '.namesorted.bam')
    if os.path.exists(sorted_path):
        if os.path.getmtime(sorted_path) > os.path.getmtime(bam_path):
            logger.info(f"Using existing name-sorted BAM: {sorted_path}")
            return sorted_path

    logger.info(f"Name-sorting BAM: {bam_path} -> {sorted_path}")
    pysam.sort('-n', '-o', sorted_path, bam_path)
    return sorted_path


def _read_id_hash(read_id: str, n_buckets: int) -> int:
    """Deterministic hash of read_id for SLURM array splitting."""
    h = hashlib.md5(read_id.encode()).hexdigest()
    return int(h, 16) % n_buckets


def _filtered_read_iterator(bam: pysam.AlignmentFile):
    """Yield only primary, mapped reads from a BAM file."""
    for read in bam:
        if not (read.is_unmapped or read.is_secondary or read.is_supplementary):
            yield read


def _iter_name_grouped_bams(bam_paths: Dict[str, str]):
    """
    K-way merge across name-sorted BAMs, yielding all alignments per read.

    Memory: O(n_aligners) per read instead of O(total_reads * n_aligners).
    """
    bams = {}
    iterators = {}
    for aligner, path in bam_paths.items():
        bam = pysam.AlignmentFile(path, 'rb')
        bams[aligner] = bam
        iterators[aligner] = _filtered_read_iterator(bam)

    current_reads = {}
    for aligner, it in iterators.items():
        try:
            current_reads[aligner] = next(it)
        except StopIteration:
            current_reads[aligner] = None

    try:
        while any(r is not None for r in current_reads.values()):
            min_read_id = min(
                r.query_name for r in current_reads.values() if r is not None
            )
            group = {}
            for aligner in list(current_reads.keys()):
                read = current_reads[aligner]
                if read is not None and read.query_name == min_read_id:
                    group[aligner] = read
                    try:
                        current_reads[aligner] = next(iterators[aligner])
                    except StopIteration:
                        current_reads[aligner] = None
            yield min_read_id, group
    finally:
        for bam in bams.values():
            bam.close()


def _process_and_write_batch(read_batch, raw_read_batch, genome, annotated_junctions, out_bam, stats):
    """Process a batch of reads and write best alignments to output BAM."""
    for i, (read_id, alignments) in enumerate(read_batch):
        result = select_best_alignment(alignments, genome, annotated_junctions)
        if result.confidence == 'high':
            stats['consensus_high'] += 1
        elif result.confidence == 'medium':
            stats['consensus_medium'] += 1
        else:
            stats['consensus_low'] += 1
        if result.was_5prime_rescued:
            stats['5prime_rescued'] += 1
        stats['by_aligner'][result.best_aligner] += 1

        _, aligner_reads = raw_read_batch[i]
        if result.best_aligner in aligner_reads:
            best_read = aligner_reads[result.best_aligner]
            best_read.set_tag('XA', result.best_aligner)
            best_read.set_tag('XC', result.confidence)
            best_read.set_tag('XN', result.n_aligners_agree)
            if result.was_5prime_rescued:
                best_read.set_tag('XR', 1)
            if result.false_junction_removed:
                best_read.set_tag('XF', 1)
            out_bam.write(best_read)


def run_consensus_selection(
    bam_paths: Dict[str, str],
    genome: Dict[str, str],
    output_bam: str,
    annotated_junctions: Optional[Set[Tuple[str, int, int]]] = None,
    write_all_to_tag: bool = True,
    n_workers: int = 0,
    batch_size: int = 10000,
    slurm_array_task: Optional[int] = None,
    slurm_array_total: Optional[int] = None,
) -> Dict[str, int]:
    """
    Run consensus selection across multiple BAM files.

    Streams through name-sorted BAMs to avoid loading all reads into memory.
    Supports SLURM array job splitting for cluster-scale parallelism.

    Memory usage: O(batch_size * n_aligners) instead of O(total_reads * n_aligners).

    Args:
        bam_paths: Dict mapping aligner name to BAM path
        genome: Dict mapping chrom to sequence
        output_bam: Output path for consensus BAM
        annotated_junctions: Optional set of annotated junctions
        write_all_to_tag: If True, write all aligner info to BAM tags
        n_workers: Number of worker processes (0 = auto-detect, 1 = single-threaded)
        batch_size: Number of read groups to accumulate before processing
        slurm_array_task: Current SLURM array task ID (0-indexed).
                          When set, only reads where
                          hash(read_id) % slurm_array_total == slurm_array_task
                          are processed.
        slurm_array_total: Total number of SLURM array tasks.

    Returns:
        Summary statistics dict
    """
    from ..slurm import get_available_cpus, get_slurm_info

    # Auto-detect SLURM array settings from environment
    if slurm_array_task is None and slurm_array_total is None:
        slurm_info = get_slurm_info()
        if slurm_info.get('array_task_id') is not None:
            try:
                slurm_array_task = int(slurm_info['array_task_id'])
                slurm_array_total = int(os.environ.get(
                    'SLURM_ARRAY_TASK_COUNT',
                    os.environ.get('SLURM_ARRAY_TASK_MAX', '0')
                ))
                if slurm_array_total > 0:
                    task_min = int(os.environ.get('SLURM_ARRAY_TASK_MIN', '0'))
                    task_step = int(os.environ.get('SLURM_ARRAY_TASK_STEP', '1'))
                    if 'SLURM_ARRAY_TASK_COUNT' not in os.environ:
                        slurm_array_total = (slurm_array_total - task_min) // task_step + 1
                    logger.info(
                        f"SLURM array detected: task {slurm_array_task} of {slurm_array_total}"
                    )
                else:
                    slurm_array_task = None
                    slurm_array_total = None
            except (ValueError, TypeError):
                slurm_array_task = None
                slurm_array_total = None

    use_slurm_filter = (
        slurm_array_task is not None and
        slurm_array_total is not None and
        slurm_array_total > 1
    )

    # Auto-detect workers
    if n_workers <= 0:
        n_workers = get_available_cpus()

    # Ensure BAMs are name-sorted
    logger.info("Ensuring BAMs are name-sorted...")
    sorted_bam_paths = {}
    for aligner, path in bam_paths.items():
        sorted_bam_paths[aligner] = _ensure_name_sorted(path)

    # Get header from first BAM
    first_bam_path = list(sorted_bam_paths.values())[0]
    first_bam = pysam.AlignmentFile(first_bam_path, 'rb')
    header = first_bam.header.to_dict()
    first_bam.close()

    # Add program group for RECTIFY consensus
    if 'PG' not in header:
        header['PG'] = []
    header['PG'].append({
        'ID': 'RECTIFY',
        'PN': 'RECTIFY',
        'VN': '2.0',
        'CL': f'consensus selection from {",".join(bam_paths.keys())}',
    })

    # Modify output path for SLURM array tasks
    if use_slurm_filter:
        base, ext = os.path.splitext(output_bam)
        output_bam = f"{base}.task{slurm_array_task}{ext}"
        logger.info(f"SLURM array task {slurm_array_task}: writing to {output_bam}")

    # Open output BAM
    out_bam = pysam.AlignmentFile(
        output_bam, 'wb',
        header=pysam.AlignmentHeader.from_dict(header)
    )

    # Initialize stats
    stats = {
        'total_reads': 0,
        'reads_skipped_slurm_filter': 0,
        'consensus_high': 0,
        'consensus_medium': 0,
        'consensus_low': 0,
        '5prime_rescued': 0,
        'by_aligner': defaultdict(int),
    }

    # Stream through name-sorted BAMs
    logger.info(f"Streaming consensus selection (batch_size={batch_size})...")
    if use_slurm_filter:
        logger.info(
            f"  SLURM array filter: task {slurm_array_task}/{slurm_array_total}"
        )

    # Accumulate batches for processing
    read_batch = []
    raw_read_batch = []
    n_batches = 0

    for read_id, aligner_reads in _iter_name_grouped_bams(sorted_bam_paths):
        # SLURM array filtering
        if use_slurm_filter:
            if _read_id_hash(read_id, slurm_array_total) != slurm_array_task:
                stats['reads_skipped_slurm_filter'] += 1
                continue

        stats['total_reads'] += 1

        # Extract alignment info for scoring
        alignments = {}
        for aligner, read in aligner_reads.items():
            alignments[aligner] = extract_alignment_info(read, aligner, genome)

        read_batch.append((read_id, alignments))
        raw_read_batch.append((read_id, aligner_reads))

        # Process batch when full
        if len(read_batch) >= batch_size:
            _process_and_write_batch(
                read_batch, raw_read_batch, genome,
                annotated_junctions, out_bam, stats
            )
            read_batch = []
            raw_read_batch = []
            n_batches += 1

            if stats['total_reads'] % 100000 == 0:
                logger.info(f"  Processed {stats['total_reads']:,} reads...")

    # Process remaining reads
    if read_batch:
        _process_and_write_batch(
            read_batch, raw_read_batch, genome,
            annotated_junctions, out_bam, stats
        )
        n_batches += 1

    # Close output
    out_bam.close()

    # Sort output by coordinate and index
    logger.info("Coordinate-sorting output BAM...")
    sorted_output = output_bam.replace('.bam', '.sorted.bam')
    pysam.sort('-o', sorted_output, output_bam)
    os.replace(sorted_output, output_bam)
    pysam.index(output_bam)

    # Log summary
    logger.info(f"\nConsensus selection complete:")
    logger.info(f"  Total reads processed: {stats['total_reads']}")
    if use_slurm_filter:
        logger.info(f"  Reads skipped (other SLURM tasks): {stats['reads_skipped_slurm_filter']}")
    logger.info(f"  High confidence: {stats['consensus_high']}")
    logger.info(f"  Medium confidence: {stats['consensus_medium']}")
    logger.info(f"  Low confidence: {stats['consensus_low']}")
    logger.info(f"  5' rescued: {stats['5prime_rescued']}")
    logger.info(f"  By aligner: {dict(stats['by_aligner'])}")
    logger.info(f"  Batches processed: {n_batches}")

    return stats


def merge_slurm_array_bams(
    output_bam_pattern: str,
    n_tasks: int,
    merged_output: str,
):
    """
    Merge BAM files from SLURM array tasks into a single output.

    Call this after all array tasks have completed.

    Args:
        output_bam_pattern: Pattern with {task} placeholder
        n_tasks: Number of array tasks
        merged_output: Path for merged output BAM
    """
    task_bams = []
    for task_id in range(n_tasks):
        bam_path = output_bam_pattern.format(task=task_id)
        if os.path.exists(bam_path):
            task_bams.append(bam_path)
        else:
            logger.warning(f"Missing SLURM array task BAM: {bam_path}")

    if not task_bams:
        raise FileNotFoundError("No SLURM array task BAMs found")

    logger.info(f"Merging {len(task_bams)} SLURM array task BAMs...")
    pysam.merge('-f', merged_output, *task_bams)

    sorted_output = merged_output.replace('.bam', '.sorted.bam')
    pysam.sort('-o', sorted_output, merged_output)
    os.replace(sorted_output, merged_output)
    pysam.index(merged_output)

    logger.info(f"Merged output: {merged_output}")

    for bam_path in task_bams:
        idx_path = bam_path + '.bai'
        if os.path.exists(idx_path):
            os.remove(idx_path)
        os.remove(bam_path)

    logger.info("SLURM array merge complete")


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
