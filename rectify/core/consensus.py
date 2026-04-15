"""
Multi-Aligner Consensus Module for RECTIFY.

Compares alignments from multiple aligners and selects the best one per read
based on junction quality, soft-clip rescue, and false junction detection.

Scoring priorities:
1. Prefer alignments that splice through junctions vs soft-clipping (5' rescue):
   - If a 5' soft-clip can be explained by a missed intron (the clipped bases
     match the upstream exon end, checked via edit distance), the soft-clip is
     "rescued" and carries no penalty. This avoids penalizing aligners that
     correctly identify a junction but soft-clip a few upstream exon bases,
     relative to aligners that simply start mapping AFTER the junction.
   - Per-read junction pool = annotated junctions UNION all aligners' observed
     junctions for this read.
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
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from pathlib import Path

import pysam
import numpy as np

logger = logging.getLogger(__name__)

# Module-level A-tract result cache: (chrom, pos, strand) → atract dict
# Many reads share the same 3' endpoint; this avoids repeated sequence lookups.
_atract_cache: Dict[Tuple[str, int, str], dict] = {}

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

    # Soft-clipped sequence at the 5' end (for sequence-based rescue scoring)
    five_prime_softclip_seq: str = ""

    # Query bases in the terminal mismatch/indel region at the 5' end.
    # mapPacBio forces mismatches/indels at splice junction boundaries instead of
    # soft-clipping; those forced-error bases are structurally equivalent to a
    # soft-clip and should be used for sequence-based rescue if five_prime_softclip_seq
    # is empty.  Length = effective_five_prime_clip - five_prime_softclip.
    effective_five_prime_clip_seq: str = ""

    # Quality scores
    junction_score: float = 0.0
    canonical_count: int = 0
    non_canonical_count: int = 0

    # Effective 5' clip for scoring: max(explicit_soft_clip, terminal_error_region).
    # Some aligners (mapPacBio) substitute forced mismatches/indels for soft-clips at
    # splice junction boundaries — identical structural situation, different CIGAR encoding.
    # This field ensures fair scoring regardless of the aligner's clipping policy.
    effective_five_prime_clip: int = 0

    # Effective 3' clip for scoring: non-poly(A) terminal errors at the 3' end.
    # Aligners should only soft-clip the poly(A) tail; clipping or force-mismatching
    # real exon sequence means the aligner stopped before the true 3' end.
    # Complementary to the A-tract depth penalty (which catches going too far INTO poly(A)).
    effective_three_prime_clip: int = 0

    # Mismatches + indels within junction_window_bp of any splice junction.
    # Cleaner aligners (mapPacBio) score 0; aligners with forced errors near
    # junctions (minimap2, gapmm2) accumulate a penalty here.
    # Float because homopolymer-context errors are weighted 0.5 (nanopore DRS
    # commonly undercalls homopolymer lengths, making such errors less diagnostic
    # of misalignment than non-homopolymer mismatches).
    junction_proximity_errors: float = 0.0

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
    n_tied_score: int = 1      # Number of aligners with equal top junction score
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


def _rescue_5prime_softclip(
    alignment: AlignmentInfo,
    genome: Dict[str, str],
    candidate_junctions: Set[Tuple[str, int, int, str]],
    max_edit_frac: float = 0.2,
    search_window_bp: int = 300,
    rescue_seq_override: str = "",
) -> bool:
    """Check whether a 5' soft-clip (or MPB terminal error region) is explained
    by a missed upstream intron.

    Algorithm (2-pass rescue):
      Pass 1 — candidate junction pool is built by the caller from annotated
               junctions UNION junctions detected by any aligner for this read.
      Pass 2 — For each candidate junction that is UPSTREAM of the alignment's
               5' end (i.e., intron_end ≤ align_5prime for + strand, or
               intron_start ≥ reference_end for − strand), within search_window_bp:
                 1. Fetch the upstream exon-end sequence (last N bases of the exon
                    before the intron donor, where N = rescue sequence length).
                 2. Compute edit distance between rescue query bases and that
                    reference sequence.
                 3. If distance / clip_len ≤ max_edit_frac → rescue.

    Rescue sequence priority:
      rescue_seq_override (explicit caller-supplied)
        > five_prime_softclip_seq (explicit soft-clip, e.g. minimap2/gapmm2)
        > effective_five_prime_clip_seq (MPB forced-mismatch/indel terminal region)

    Coordinate conventions:
      + strand: alignment 5' end = reference_start (leftmost mapped base).
                Upstream intron has intron_end ≤ reference_start.
                Exon upstream of intron: genome[intron_start - N : intron_start].
      − strand: alignment 5' end = reference_end − 1 (rightmost mapped base).
                Upstream intron (in transcript) has intron_start ≥ reference_end.
                Exon upstream of intron (in transcript) = genome[intron_end : intron_end + N].
                BAM stores reverse-strand query_sequence in reference orientation,
                so no reverse-complement is needed for the comparison.

    Args:
        alignment: AlignmentInfo with five_prime_softclip_seq (and optionally
            effective_five_prime_clip_seq) populated
        genome: Dict mapping chrom to sequence
        candidate_junctions: Set of (chrom, intron_start, intron_end, strand)
        max_edit_frac: Max edit distance / clip_len to declare a rescue
        search_window_bp: Max distance from 5' alignment boundary to intron edge
        rescue_seq_override: Explicit sequence to use instead of alignment fields.
            When supplied, skips the field-selection logic entirely.

    Returns:
        True if the clip/error region is explained by a missed intron (no penalty),
        False if unexplained (apply normal penalty).
    """
    # Priority: explicit override > soft-clip seq > MPB terminal error seq
    clip_seq = (
        rescue_seq_override
        or alignment.five_prime_softclip_seq
        or alignment.effective_five_prime_clip_seq
    )
    if not clip_seq:
        return False

    clip_len = len(clip_seq)
    chrom = alignment.chrom

    if chrom not in genome:
        return False

    from .spikein_filter import edit_distance
    genome_seq = genome[chrom]
    clip_seq_upper = clip_seq.upper()  # Hoist out of per-junction loop

    if alignment.strand == '+':
        # 5' alignment start = leftmost mapped base
        align_5prime = alignment.reference_start
        for (j_chrom, intron_start, intron_end, j_strand) in candidate_junctions:
            if j_chrom != chrom:
                continue
            # Upstream intron: its 3'SS (intron_end) must be at or before align_5prime.
            # Directional: intron_end ≤ align_5prime (not abs, to exclude internal junctions).
            dist = align_5prime - intron_end
            if dist < 0 or dist > search_window_bp:
                continue
            # Exon1-end sequence: last clip_len bases of exon before the intron donor
            exon_end_start = intron_start - clip_len
            if exon_end_start < 0:
                continue
            exon_end_seq = genome_seq[exon_end_start:intron_start].upper()
            ed = edit_distance(clip_seq_upper, exon_end_seq)
            if ed / clip_len <= max_edit_frac:
                logger.debug(
                    f"5' rescue (+): {chrom}:{align_5prime} clip={clip_len}bp "
                    f"matches exon1 end before intron {intron_start}-{intron_end} "
                    f"(edit={ed}/{clip_len})"
                )
                return True
    else:
        # Minus strand: 5' end is at reference_end - 1 (rightmost mapped base).
        # In transcript orientation, "upstream" means higher reference coordinates.
        align_5prime = alignment.reference_end - 1
        reference_end = alignment.reference_end
        for (j_chrom, intron_start, intron_end, j_strand) in candidate_junctions:
            if j_chrom != chrom:
                continue
            # Upstream intron (in transcript): intron_start must be ≥ reference_end.
            # Directional: intron_start > align_5prime (excludes internal junctions).
            dist = intron_start - align_5prime
            if dist == 0 or dist > search_window_bp:
                continue
            # Exon upstream of intron (in transcript) = exon to the right of intron_end
            # in reference. BAM reverse-strand query_sequence is in reference orientation,
            # so compare directly without reverse-complementing.
            exon_end_seq = genome_seq[intron_end:intron_end + clip_len].upper()
            if len(exon_end_seq) < clip_len:
                continue
            ed = edit_distance(clip_seq_upper, exon_end_seq)
            if ed / clip_len <= max_edit_frac:
                logger.debug(
                    f"5' rescue (−): {chrom}:{align_5prime} clip={clip_len}bp "
                    f"matches exon upstream of intron {intron_start}-{intron_end} "
                    f"(edit={ed}/{clip_len})"
                )
                return True

    return False


def _get_effective_5prime_clip(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    scan_bp: int = 25,
    window: int = 8,
    error_threshold: float = 0.40,
    min_errors: int = 2,
) -> int:
    """
    Compute effective 5' clip length including both explicit soft-clips and
    terminal mismatch/indel regions.

    Some aligners (notably mapPacBio) substitute forced mismatches/indels for
    soft-clips at splice junction boundaries — the read cannot be aligned there
    due to a missed or partially-resolved intron, but instead of clipping, the
    aligner records mismatches. This is structurally the same situation as a
    soft-clip and should receive the same scoring penalty.

    Algorithm:
      1. Collect per-position error flags (mismatch or insertion) for the first
         `scan_bp` aligned bases from the 5' end (after explicit soft-clips).
      2. Scan with a sliding window of size `window`. The terminal error region
         extends as long as the leading window(s) have density >= error_threshold.
         Stop at the first clean window (greedy from 5' end).
      3. Return max(explicit_soft_clip, explicit_soft_clip + terminal_error_length),
         provided at least `min_errors` errors were found in the terminal region.

    The explicit soft-clip sequence in `five_prime_softclip_seq` is NOT updated —
    the sequence-based rescue uses that field and operates on the clip sequence only.

    Args:
        read: pysam aligned read
        genome: chromosome → sequence dict (for mismatch detection without MD tags)
        scan_bp: aligned positions to scan from 5' end (after explicit soft-clips)
        window: sliding window size for density estimation
        error_threshold: fraction of errors per window to qualify as terminal region
        min_errors: minimum total errors to apply terminal clip extension

    Returns:
        Effective 5' clip length (>= explicit soft-clip).
    """
    five_clip, _ = get_softclip_lengths(read)

    if not read.query_sequence or not read.cigartuples:
        return five_clip

    chrom = read.reference_name
    if not chrom or chrom not in genome:
        return five_clip

    ref_seq = genome[chrom]
    query_seq = read.query_sequence
    query_len = len(query_seq)

    try:
        pairs = read.get_aligned_pairs()
    except Exception as e:
        logger.warning(
            "_get_effective_5prime_clip: get_aligned_pairs failed for read %s: %s",
            getattr(read, "query_name", "<unknown>"),
            e,
        )
        return five_clip

    # Build error array for the first `scan_bp` aligned bases from the 5' end.
    # Skip explicit soft-clip positions; only examine the true aligned region.
    errors: List[int] = []
    if read.is_reverse:
        # 5' end is at high query coordinates; scan backward from cutoff
        cutoff_qp = query_len - 1 - five_clip  # last pos inside aligned region
        for qp, rp in reversed(pairs):
            if qp is None:
                continue
            if qp > cutoff_qp:
                continue  # still in soft-clip region
            if len(errors) >= scan_bp:
                break
            if rp is None:
                errors.append(1)  # insertion into read = mismatch-equivalent
            elif rp < len(ref_seq):
                ref_b = ref_seq[rp].upper()
                read_b = query_seq[qp].upper()
                errors.append(1 if (ref_b != 'N' and read_b != ref_b) else 0)
            else:
                errors.append(0)
    else:
        # 5' end is at low query coordinates; scan forward from cutoff
        cutoff_qp = five_clip  # first pos inside aligned region
        for qp, rp in pairs:
            if qp is None:
                continue
            if qp < cutoff_qp:
                continue  # still in soft-clip region
            if len(errors) >= scan_bp:
                break
            if rp is None:
                errors.append(1)
            elif rp < len(ref_seq):
                ref_b = ref_seq[rp].upper()
                read_b = query_seq[qp].upper()
                errors.append(1 if (ref_b != 'N' and read_b != ref_b) else 0)
            else:
                errors.append(0)

    if len(errors) < window:
        return five_clip

    # Greedy scan from 5' end: extend terminal boundary while leading windows are
    # high-error. Stop at the first clean window — terminal errors are contiguous.
    terminal_end = 0
    for i in range(len(errors) - window + 1):
        w = errors[i:i + window]
        n_err = sum(w)
        if n_err / window >= error_threshold:
            terminal_end = i + window
        else:
            break  # first clean window; stop extending

    if terminal_end == 0:
        return five_clip

    total_errors = sum(errors[:terminal_end])
    if total_errors < min_errors:
        return five_clip

    return five_clip + terminal_end


def _get_effective_3prime_clip(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    scan_bp: int = 25,
    window: int = 8,
    error_threshold: float = 0.40,
    min_errors: int = 2,
) -> int:
    """
    Compute effective 3' clip length for non-poly(A) terminal errors.

    Aligners should only soft-clip the poly(A) tail at the 3' end. Clipping or
    force-mismatching real exon sequence (non-A on + strand, non-T on - strand in
    read coordinates) means the aligner stopped before the true 3' end and should
    be penalized.

    Key distinction from 5' end: poly(A) base errors/clips are EXPECTED at the 3'
    end and are already handled by the A-tract depth penalty. Only non-A terminal
    errors are counted here.

    Algorithm:
      1. Collect per-position error flags for the first `scan_bp` aligned bases
         from the 3' end, scanning inward. Mark a position as an error only if it
         is a mismatch/insertion AND the read base is non-A (+ strand) or non-T
         (- strand in BAM coords, where poly(A) appears as poly(T)).
      2. Apply the same greedy sliding-window scan as _get_effective_5prime_clip.
      3. Return max(explicit_3prime_clip, explicit_3prime_clip + terminal_end).

    The explicit 3' soft-clip is also filtered: only non-A/T bases within the
    clip sequence count toward the penalty.

    Args:
        read: pysam aligned read
        genome: chromosome → sequence dict
        scan_bp: aligned positions to scan inward from 3' end
        window: sliding window size
        error_threshold: fraction of errors per window (non-A only)
        min_errors: minimum non-A errors to trigger penalty

    Returns:
        Effective 3' clip length (>= explicit soft-clip).
    """
    _, three_clip = get_softclip_lengths(read)

    if not read.query_sequence or not read.cigartuples:
        return three_clip

    chrom = read.reference_name
    if not chrom or chrom not in genome:
        return three_clip

    ref_seq = genome[chrom]
    query_seq = read.query_sequence
    query_len = len(query_seq)

    # Poly(A) base in read coordinates:
    # - Plus strand: A in query = expected at 3' end (poly(A) tail)
    # - Minus strand: poly(A) tail is stored as poly(T) in the reverse-complemented BAM query
    polya_base = 'T' if read.is_reverse else 'A'

    try:
        pairs = read.get_aligned_pairs()
    except Exception:
        return three_clip

    # Build error array for the first `scan_bp` aligned bases from the 3' end.
    # Scan inward (3' → 5') and flag only non-poly(A) errors.
    errors: List[int] = []
    if read.is_reverse:
        # 3' RNA end is at the LEFT (low query positions) for minus strand.
        # Explicit 3' clip occupies query_seq[:three_clip]; scan just inside.
        cutoff_qp = three_clip  # first pos inside aligned region
        for qp, rp in pairs:
            if qp is None:
                continue
            if qp < cutoff_qp:
                continue
            if len(errors) >= scan_bp:
                break
            read_b = query_seq[qp].upper()
            if rp is None:
                # Insertion — count as error only if not a poly(A)-base insertion
                errors.append(0 if read_b == polya_base else 1)
            elif rp < len(ref_seq):
                ref_b = ref_seq[rp].upper()
                is_mismatch = ref_b != 'N' and read_b != ref_b
                # Only count non-polyA mismatches
                errors.append(1 if (is_mismatch and read_b != polya_base) else 0)
            else:
                errors.append(0)
    else:
        # 3' RNA end is at the RIGHT (high query positions) for plus strand.
        # Explicit 3' clip occupies query_seq[-three_clip:]; scan just inside.
        cutoff_qp = query_len - 1 - three_clip  # last pos inside aligned region
        for qp, rp in reversed(pairs):
            if qp is None:
                continue
            if qp > cutoff_qp:
                continue
            if len(errors) >= scan_bp:
                break
            read_b = query_seq[qp].upper()
            if rp is None:
                errors.append(0 if read_b == polya_base else 1)
            elif rp < len(ref_seq):
                ref_b = ref_seq[rp].upper()
                is_mismatch = ref_b != 'N' and read_b != ref_b
                errors.append(1 if (is_mismatch and read_b != polya_base) else 0)
            else:
                errors.append(0)

    if len(errors) < window:
        return three_clip

    # Same greedy scan as _get_effective_5prime_clip
    terminal_end = 0
    for i in range(len(errors) - window + 1):
        w = errors[i:i + window]
        n_err = sum(w)
        if n_err / window >= error_threshold:
            terminal_end = i + window
        else:
            break

    if terminal_end == 0:
        return three_clip

    total_errors = sum(errors[:terminal_end])
    if total_errors < min_errors:
        return three_clip

    return max(three_clip, three_clip + terminal_end)


def _is_homopolymer_position(ref_seq: str, rp: int, min_run: int = 3) -> bool:
    """Return True if reference position rp is within a homopolymer run >= min_run."""
    if rp < 0 or rp >= len(ref_seq):
        return False
    base = ref_seq[rp].upper()
    if base == 'N':
        return False
    left = rp
    while left > 0 and ref_seq[left - 1].upper() == base:
        left -= 1
    right = rp + 1
    while right < len(ref_seq) and ref_seq[right].upper() == base:
        right += 1
    return (right - left) >= min_run


def _count_junction_proximity_errors(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    junction_window_bp: int = 5,
) -> float:
    """
    Count mismatches and indels within junction_window_bp of any splice junction.

    For each N operation (intron), inspect the `junction_window_bp` aligned
    bases on the exon side of both junction boundaries.  Errors (mismatches,
    insertions, deletions) in those windows indicate an aligner that placed
    indels or mismatches right at the junction rather than resolving the
    splice cleanly.

    Insertions are attributed to the preceding aligned reference position;
    deletions are attributed to their reference position directly.

    Errors at homopolymer reference positions are weighted 0.5 instead of 1.0.
    Nanopore DRS commonly undercalls homopolymer lengths; such errors are
    sequencing artifacts rather than evidence of misalignment.

    Returns:
        Total weighted error count summed across all junction-proximal windows.
    """
    if not read.cigartuples or not read.query_sequence:
        return 0

    chrom = read.reference_name
    if not chrom or chrom not in genome:
        return 0

    junctions = extract_junctions_from_cigar(read)
    if not junctions:
        return 0

    ref_seq = genome[chrom]
    query_seq = read.query_sequence

    # Build set of reference positions within junction_window_bp of any boundary.
    prox: set = set()
    for junc_start, junc_end in junctions:
        for rp in range(max(0, junc_start - junction_window_bp), junc_start):
            prox.add(rp)
        for rp in range(junc_end, min(len(ref_seq), junc_end + junction_window_bp)):
            prox.add(rp)

    try:
        pairs = read.get_aligned_pairs()
    except Exception:
        return 0

    errors: float = 0.0
    prev_rp: Optional[int] = None

    for qp, rp in pairs:
        if rp is None and qp is not None:
            # Insertion into read — attribute to preceding ref position
            if prev_rp is not None and prev_rp in prox:
                weight = 0.5 if _is_homopolymer_position(ref_seq, prev_rp) else 1.0
                errors += weight
        elif qp is None and rp is not None:
            # Deletion from read
            if rp in prox:
                weight = 0.5 if _is_homopolymer_position(ref_seq, rp) else 1.0
                errors += weight
            prev_rp = rp
        elif qp is not None and rp is not None:
            # Aligned pair — check for mismatch
            if rp in prox and rp < len(ref_seq):
                ref_b = ref_seq[rp].upper()
                read_b = query_seq[qp].upper()
                if ref_b != 'N' and read_b != ref_b:
                    weight = 0.5 if _is_homopolymer_position(ref_seq, rp) else 1.0
                    errors += weight
            prev_rp = rp

    return errors


def score_alignment(
    alignment: AlignmentInfo,
    genome: Dict[str, str],
    candidate_junctions: Optional[Set[Tuple[str, int, int, str]]] = None,
) -> float:
    """
    Score an alignment based on junction quality.

    Scoring factors:
    - 5' soft-clip penalty: -2 per unexplained clipped base; 0 if sequence-rescue
      confirms the aligner found an intron but could not align the upstream exon end
    - 3' A-tract depth penalty: -1 per downstream A (capped at 10)

    Neither canonical splice site motifs (GT/AG) nor annotated junction
    matches are scored here, to avoid biasing against novel junctions.
    Both are used only as tiebreakers in select_best_alignment().

    Note: False 3' junctions from poly(A) artifacts are handled by the
    walk back correction step, which eats through aligned A's and discards
    spurious N operations to find the true CPA site.

    Args:
        alignment: AlignmentInfo for this aligner
        genome: Dict mapping chrom to sequence
        candidate_junctions: Per-read junction pool (annotated + all aligners'
            observed junctions). Used for 5' soft-clip rescue. When provided,
            a soft-clip that matches an upstream exon end is NOT penalized.
    """
    score = 0.0

    # 5' terminal penalty: penalize both explicit soft-clips and terminal mismatch/
    # indel regions equivalently. `effective_five_prime_clip` = max(explicit_clip,
    # terminal_error_length) so aligners that force mismatches instead of clipping
    # (mapPacBio) are scored on equal footing with aligners that soft-clip (minimap2).
    #
    # Rescue uses whichever sequence is available: explicit soft-clip bytes
    # (minimap2/gapmm2/uLTRA) or MPB's forced-mismatch/indel terminal region
    # (effective_five_prime_clip_seq). Both represent the same structural
    # situation — bases that should align to the upstream exon across the intron.
    effective_clip = alignment.effective_five_prime_clip
    if effective_clip > 0:
        rescue_seq = (
            alignment.five_prime_softclip_seq
            or alignment.effective_five_prime_clip_seq
        )
        clip_rescued = (
            bool(rescue_seq)
            and candidate_junctions is not None
            and genome
            and _rescue_5prime_softclip(alignment, genome, candidate_junctions)
        )
        if not clip_rescued:
            score -= effective_clip * 2

    # 3' A-tract depth penalty (prefer alignments landing closer to the true CPA).
    # Each downstream A the aligner runs into costs 1 point — same scale as 5' penalty
    # per base, but capped at 10 to avoid overwhelming junction scoring.
    score -= min(alignment.three_prime_atract_depth, 10)

    # 3' non-poly(A) terminal error penalty: penalizes aligners that stop before the
    # true 3' end. Aligners should only soft-clip the poly(A) tail; non-A/T clipping or
    # forced mismatches at the 3' end indicate missed exon coverage.
    # Same -2/base scale as the 5' penalty; capped at 10 to avoid overwhelming junction scores.
    if alignment.effective_three_prime_clip > 0:
        score -= min(alignment.effective_three_prime_clip * 2, 10)

    # Junction-proximity mismatch/indel penalty: penalize aligners that place
    # errors right at splice junction boundaries.  -1 per proximal error, capped
    # at 10 to avoid overwhelming the junction/clip scores.  This favors aligners
    # (e.g. mapPacBio) that produce clean junctions with no forced errors in the
    # flanking exon sequence, relative to those (e.g. minimap2, gapmm2) that
    # accumulate mismatches/indels within a few bp of each junction.
    if alignment.junction_proximity_errors > 0:
        score -= min(alignment.junction_proximity_errors, 10)

    # NOTE: Canonical junction motifs (GT-AG) and annotated junction matches are
    # deliberately not scored here. They are used only as tiebreakers in
    # select_best_alignment() to avoid biasing against novel non-canonical junctions.

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

    # Effective 5' clip for scoring: includes terminal mismatch/indel region.
    # Must be computed before the AlignmentInfo is constructed.
    effective_five_clip = _get_effective_5prime_clip(read, genome)

    # Effective 3' clip for scoring: non-poly(A) terminal errors at 3' end.
    effective_three_clip = _get_effective_3prime_clip(read, genome)

    # Extract soft-clipped bases at the 5' end for sequence-based rescue
    five_prime_seq = ""
    if five_clip > 0 and read.query_sequence:
        if not read.is_reverse:
            five_prime_seq = read.query_sequence[:five_clip]
        else:
            five_prime_seq = read.query_sequence[-five_clip:]

    # Extract the terminal mismatch/indel region for mapPacBio rescue.
    # MPB forces mismatches/indels instead of soft-clipping at splice boundaries;
    # effective_five_clip > five_clip means there are terminal alignment errors.
    # We extract those query bases so _rescue_5prime_softclip can compare them
    # against the upstream exon end — same test as for an explicit soft-clip.
    effective_five_prime_seq = ""
    terminal_error_len = effective_five_clip - five_clip
    if terminal_error_len > 0 and read.query_sequence:
        qlen = len(read.query_sequence)
        if not read.is_reverse:
            # Plus strand: terminal errors are the first `terminal_error_len` aligned bases
            # (immediately after the explicit soft-clip region, if any)
            region_start = five_clip
            region_end = five_clip + terminal_error_len
            if region_end <= qlen:
                effective_five_prime_seq = read.query_sequence[region_start:region_end]
        else:
            # Minus strand: terminal errors are at the high end of query_sequence
            # (immediately before the explicit soft-clip region, if any)
            region_end = qlen - five_clip
            region_start = region_end - terminal_error_len
            if region_start >= 0:
                effective_five_prime_seq = read.query_sequence[region_start:region_end]

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
        chrom_std = standardize_chrom_name(chrom) or chrom

    if genome.get(chrom_std) is not None:
        try:
            # Cache A-tract results by (chrom, pos, strand) — many reads share 3' ends
            _cache_key = (chrom_std, raw_3prime, strand)
            if _cache_key not in _atract_cache:
                _atract_cache[_cache_key] = calculate_atract_ambiguity(
                    genome, chrom_std, raw_3prime, strand, downstream_bp=10
                )
            atract = _atract_cache[_cache_key]
            atract_depth = atract.get('downstream_a_count') or 0
            # Best-guess corrected position: ambiguity_min for +, ambiguity_max for -
            if strand == '+':
                corrected_3prime = atract.get('ambiguity_min', raw_3prime)
            else:
                corrected_3prime = atract.get('ambiguity_max', raw_3prime)
        except Exception:
            pass  # Non-fatal; raw position used

    # Count mismatches/indels within 5 bp of each junction boundary.
    # Used by score_alignment() to prefer aligners with clean junction handling.
    junction_prox_errors = _count_junction_proximity_errors(read, genome)

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
        five_prime_softclip_seq=five_prime_seq,
        effective_five_prime_clip=effective_five_clip,
        effective_five_prime_clip_seq=effective_five_prime_seq,
        effective_three_prime_clip=effective_three_clip,
        corrected_3prime=corrected_3prime,
        three_prime_atract_depth=atract_depth,
        canonical_count=canonical,
        non_canonical_count=non_canonical,
        junction_proximity_errors=junction_prox_errors,
        has_false_3prime_junction=False,  # Not used; walk back handles this
    )


def select_best_alignment(
    alignments: Dict[str, AlignmentInfo],
    genome: Dict[str, str],
    annotated_junctions: Optional[Set[Tuple[str, int, int, str]]] = None,
) -> ConsensusResult:
    """
    Select the best alignment from multiple aligners for a single read.

    Args:
        alignments: Dict mapping aligner name to AlignmentInfo
        genome: Dict mapping chrom to sequence
        annotated_junctions: Optional set of (chrom, start, end, strand) for annotated junctions

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

    # Build per-read junction pool for 5' soft-clip rescue.
    # Pool = annotated junctions UNION all junctions observed by any aligner for
    # this read. Using all aligners' junctions ensures that if any aligner correctly
    # identifies an intron, soft-clips at that intron in other aligners are rescued.
    chrom_for_read = list(alignments.values())[0].chrom
    candidate_junctions: Set[Tuple[str, int, int, str]] = set()
    if annotated_junctions:
        # Only keep annotated junctions on the same chrom to avoid scanning everything
        for j in annotated_junctions:
            if j[0] == chrom_for_read:
                candidate_junctions.add(j)
    for alignment in alignments.values():
        for junc_start, junc_end in alignment.junctions:
            candidate_junctions.add((alignment.chrom, junc_start, junc_end, alignment.strand))

    # Score all alignments
    for aligner, alignment in alignments.items():
        alignment.junction_score = score_alignment(alignment, genome, candidate_junctions)

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
                    if (a.chrom, junc[0], junc[1], a.strand) in annotated_junctions
                )
            return (_count_3prime_agreement(aligner_name), n_annotated, a.canonical_count)

        best_aligner = max(tied_aligners, key=_tiebreak_key)

    best_alignment = alignments[best_aligner]
    best_alignment.is_best = True

    # Check for 5' rescue using effective clip (includes terminal mismatch regions)
    # Did we pick an alignment that spliced through vs one that soft-clipped/had terminal errors?
    was_rescued = False
    if len(alignments) > 1:
        min_5clip = min(a.effective_five_prime_clip for a in alignments.values())
        max_5clip = max(a.effective_five_prime_clip for a in alignments.values())
        if max_5clip > min_5clip and best_alignment.effective_five_prime_clip == min_5clip:
            was_rescued = True

    # Count junction agreement across aligners
    junction_sets = {}
    for aligner, alignment in alignments.items():
        junc_key = (alignment.strand, tuple(sorted(alignment.junctions)))
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
        n_tied_score=len(tied_aligners),
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
    # Cap samtools sort memory: without -m, samtools uses 768MB × all threads.
    # With 5 BAMs sorted sequentially, Python's allocator retains each peak
    # in RSS, compounding to ~60GB on a 16-core node. 1G per sort is ample
    # for typical per-sample BAM sizes.
    pysam.sort('-n', '-m', '1G', '-o', sorted_path, bam_path)
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


def _natural_sort_key(s: str) -> list:
    """Key for natural (version) sort matching samtools queryname:natural order.

    Samtools natural sort compares runs of digits as integers rather than
    lexicographically.  Example: ``98297e97`` sorts before ``0633141e``
    because 98297 < 633141, even though ``'9' > '0'`` lexicographically.
    The K-way merge must use the same ordering as the BAM iterators so that
    reads present in only a subset of aligners do not desynchronise the merge.
    """
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', s)]


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
                (r.query_name for r in current_reads.values() if r is not None),
                key=_natural_sort_key,
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


def _cigar_query_length(read: pysam.AlignedSegment) -> int:
    """Return the total number of query-consuming bases implied by the CIGAR."""
    if not read.cigartuples:
        return 0
    # ops that consume query: M=0, I=1, S=4, =7, X=8
    query_ops = {0, 1, 4, 7, 8}
    return sum(length for op, length in read.cigartuples if op in query_ops)


def _restore_sequence_from_aligner_reads(
    best_read: pysam.AlignedSegment,
    aligner_reads: Dict[str, pysam.AlignedSegment],
) -> None:
    """Copy query sequence and quality scores to best_read from another aligner.

    gapmm2 outputs PAF which carries no read sequence, so _paf_to_bam() leaves
    query_sequence=None on every gapmm2 BAM record.  When gapmm2 wins consensus
    selection the output BAM would contain SEQ=* records that break all
    downstream steps (indel correction, poly-A trimming, etc.).

    This function looks through the other aligners' reads for the same read_id
    and copies the first non-None sequence whose length matches the CIGAR's
    expected query length.  Donors with mismatched lengths (e.g. hard-clipped
    records from deSALT) are skipped to prevent samtools from rejecting the
    BAM with "CIGAR and query sequence lengths differ".

    If no donor with the correct length is found, best_read is left as SEQ=*
    and a warning is logged.

    Args:
        best_read: The winning pysam.AlignedSegment (modified in place).
        aligner_reads: Dict mapping aligner name to pysam.AlignedSegment for
                       the same read_id.
    """
    expected_len = _cigar_query_length(best_read)
    for donor_read in aligner_reads.values():
        seq = donor_read.query_sequence
        if seq is None:
            continue
        if expected_len > 0 and len(seq) != expected_len:
            logger.debug(
                f"Skipping donor for '{best_read.query_name}': "
                f"sequence length {len(seq)} != CIGAR query length {expected_len}"
            )
            continue
        best_read.query_sequence = seq
        best_read.query_qualities = donor_read.query_qualities
        return
    logger.warning(
        f"No aligner has query_sequence for read '{best_read.query_name}'; "
        "writing SEQ=* record"
    )


def _process_and_write_batch(read_batch, raw_read_batch, genome, annotated_junctions, out_bam, stats, use_chimeric=False):
    """Process a batch of reads and write best alignments to output BAM."""
    if use_chimeric:
        from .chimeric_consensus import select_best_chimeric, build_chimeric_read

    for i, (read_id, alignments) in enumerate(read_batch):
        _, aligner_reads = raw_read_batch[i]

        if use_chimeric:
            chimeric_result = select_best_chimeric(aligner_reads, genome, annotated_junctions)

            # Pick a template read with a valid sequence (gapmm2 yields None).
            # Verify sequence length matches the chimeric CIGAR's query length
            # to prevent the "CIGAR and query sequence lengths differ" crash.
            query_ops = {0, 1, 4, 7, 8}  # M, I, S, =, X
            expected_len = sum(
                length for op, length in chimeric_result.chimeric_cigar if op in query_ops
            ) if chimeric_result.chimeric_cigar else 0
            template = None
            for r in aligner_reads.values():
                seq = r.query_sequence
                if seq is not None and (expected_len == 0 or len(seq) == expected_len):
                    template = r
                    break
            if template is None:
                logger.warning(
                    f"No valid template read for chimeric assembly of '{read_id}'; skipping"
                )
                continue

            out_read = build_chimeric_read(
                template_read=template,
                ref_start=chimeric_result.chimeric_ref_start,
                cigar_tuples=chimeric_result.chimeric_cigar,
                chimeric_result=chimeric_result,
                header=out_bam.header,
            )
            out_read.flag &= ~0x900  # enforce primary

            if chimeric_result.confidence == 'high':
                stats['consensus_high'] += 1
            elif chimeric_result.confidence == 'medium':
                stats['consensus_medium'] += 1
            else:
                stats['consensus_low'] += 1
            if chimeric_result.is_chimeric:
                stats['chimeric_reads'] += 1
            for _pos, winner, _qs, _qe in (chimeric_result.segment_winners or []):
                stats['by_aligner'][winner] += 1
            unique_winners = frozenset(
                w[1] for w in (chimeric_result.segment_winners or [])
            )
            stats['by_aligner_combo'][unique_winners] += 1

            out_bam.write(out_read)

        else:
            result = select_best_alignment(alignments, genome, annotated_junctions)
            if result.confidence == 'high':
                stats['consensus_high'] += 1
            elif result.confidence == 'medium':
                stats['consensus_medium'] += 1
            else:
                stats['consensus_low'] += 1
            if result.was_5prime_rescued:
                stats['5prime_rescued'] += 1
            if result.n_tied_score > 1:
                stats['tied_score'] += 1
            stats['by_aligner'][result.best_aligner] += 1
            stats['by_aligner_combo'][frozenset(result.aligners_compared)] += 1

            if result.best_aligner in aligner_reads:
                best_read = aligner_reads[result.best_aligner]

                # Enforce exactly one primary per read: clear secondary (0x100) and
                # supplementary (0x800) bits so the winning record is always primary.
                best_read.flag &= ~0x900

                # gapmm2 PAF→BAM conversion does not preserve read sequences;
                # restore SEQ from another aligner's record for the same read.
                if best_read.query_sequence is None:
                    _restore_sequence_from_aligner_reads(best_read, aligner_reads)

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
    annotated_junctions: Optional[Set[Tuple[str, int, int, str]]] = None,
    write_all_to_tag: bool = True,
    n_workers: int = 0,
    batch_size: int = 10000,
    slurm_array_task: Optional[int] = None,
    slurm_array_total: Optional[int] = None,
    use_chimeric: bool = False,
) -> Dict[str, int]:
    """
    Run consensus selection across multiple BAM files.

    Streams through name-sorted BAMs to avoid loading all reads into memory.
    Supports SLURM array job splitting for cluster-scale parallelism.

    Memory usage: O(batch_size * n_aligners) instead of O(total_reads * n_aligners).

    Args:
        bam_paths: Dict mapping aligner name to BAM path
        genome: Dict mapping chrom to sequence
        output_bam: Output path for rectified BAM
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

    # Auto-detect SLURM array settings from environment.
    # Only activates when RECTIFY_CONSENSUS_ARRAY_MODE=1 is explicitly set.
    # In run-all mode, SLURM array indices are for sample parallelism, not
    # read-level partitioning — do not auto-activate there.
    if slurm_array_task is None and slurm_array_total is None:
        slurm_info = get_slurm_info()
        if (slurm_info.get('array_task_id') is not None
                and os.environ.get('RECTIFY_CONSENSUS_ARRAY_MODE') == '1'):
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
    import time as _time
    _t_total = _time.perf_counter()
    logger.info("Ensuring BAMs are name-sorted...")
    _t_ns = _time.perf_counter()
    sorted_bam_paths = {}
    for aligner, path in bam_paths.items():
        sorted_bam_paths[aligner] = _ensure_name_sorted(path)
    logger.info(f"[TIMING] Name-sort: {_time.perf_counter() - _t_ns:.1f}s")

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
        'tied_score': 0,
        'chimeric_reads': 0,
        'by_aligner': defaultdict(int),
        'by_aligner_combo': defaultdict(int),  # frozenset of available aligners → count
    }

    try:
        # Stream through name-sorted BAMs
        _t_stream = _time.perf_counter()
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
                    annotated_junctions, out_bam, stats,
                    use_chimeric=use_chimeric,
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
                annotated_junctions, out_bam, stats,
                use_chimeric=use_chimeric,
            )
            n_batches += 1

    except Exception:
        out_bam.close()
        # Remove partial output BAM so callers don't see an incomplete file
        try:
            os.unlink(output_bam)
        except OSError:
            pass
        raise

    # Close output
    out_bam.close()
    logger.info(f"[TIMING] Streaming ({stats['total_reads']:,} reads, {n_batches} batches): {_time.perf_counter() - _t_stream:.1f}s")
    if stats['total_reads'] > 0:
        _reads_per_sec = stats['total_reads'] / max(_time.perf_counter() - _t_stream, 0.001)
        logger.info(f"[TIMING] Throughput: {_reads_per_sec:,.0f} reads/sec")

    # Sort output by coordinate and index
    _t_sort = _time.perf_counter()
    logger.info("Coordinate-sorting output BAM...")
    sorted_output = output_bam.replace('.bam', '.sorted.bam')
    pysam.sort('-m', '1G', '-o', sorted_output, output_bam)
    os.replace(sorted_output, output_bam)
    pysam.index(output_bam)
    logger.info(f"[TIMING] Coordinate-sort + index: {_time.perf_counter() - _t_sort:.1f}s")

    # Log summary
    logger.info(f"\nConsensus selection complete:")
    logger.info(f"  Total reads processed: {stats['total_reads']}")
    if use_slurm_filter:
        logger.info(f"  Reads skipped (other SLURM tasks): {stats['reads_skipped_slurm_filter']}")
    logger.info(f"  High confidence: {stats['consensus_high']}")
    logger.info(f"  Medium confidence: {stats['consensus_medium']}")
    logger.info(f"  Low confidence: {stats['consensus_low']}")
    logger.info(f"  5' rescued: {stats['5prime_rescued']}")
    logger.info(f"  Tied score (tiebreaker used): {stats['tied_score']}")
    if use_chimeric:
        logger.info(f"  Chimeric reads (multi-aligner segments): {stats['chimeric_reads']}")
    logger.info(f"  By aligner: {dict(stats['by_aligner'])}")
    logger.info(f"  By aligner combo: { {'+'.join(sorted(k)): v for k, v in stats['by_aligner_combo'].items()} }")
    logger.info(f"  Batches processed: {n_batches}")
    logger.info(f"[TIMING] run_consensus_selection total: {_time.perf_counter() - _t_total:.1f}s")

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
    pysam.sort('-m', '1G', '-o', sorted_output, merged_output)
    os.replace(sorted_output, merged_output)
    pysam.index(merged_output)

    logger.info(f"Merged output: {merged_output}")

    for bam_path in task_bams:
        idx_path = bam_path + '.bai'
        if os.path.exists(idx_path):
            os.remove(idx_path)
        os.remove(bam_path)

    logger.info("SLURM array merge complete")


def load_annotated_junctions(annotation_path: str) -> Set[Tuple[str, int, int, str]]:
    """
    Load annotated junctions from GFF/GTF file.

    Returns set of (chrom, intron_start, intron_end, strand) tuples where chrom is in
    standardized canonical format (chrI, chrII, etc.) so that junction lookups
    match the standardized chrom names used during correction.
    """
    from ..utils.genome import standardize_chrom_name

    junctions = set()

    import gzip as _gzip
    _open = _gzip.open if str(annotation_path).endswith('.gz') else open
    with _open(annotation_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2].lower()

            # Look for intron features
            if feature_type == 'intron':
                chrom = standardize_chrom_name(parts[0])
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])  # Already exclusive in GFF end
                strand = parts[6] if parts[6] in ('+', '-') else '+'
                junctions.add((chrom, start, end, strand))

    if len(junctions) == 0:
        logger.warning(
            "load_annotated_junctions: 0 junctions loaded from %s. "
            "Check that the file exists, is readable, and contains 'intron' "
            "feature records (column 3). Junction-guided scoring will be disabled.",
            annotation_path,
        )
    else:
        logger.info(f"Loaded {len(junctions)} annotated junctions from {annotation_path}")
    return junctions
