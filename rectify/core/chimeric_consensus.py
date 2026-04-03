"""
Chimeric Consensus Module for RECTIFY.

Instead of picking one aligner's entire alignment, this module finds "sync points"
where multiple aligners agree on query→reference mapping, scores segments between
those sync points independently, and constructs a chimeric alignment taking the
best segment from each aligner.

Key insight: DRS reads are long (1-10+ kb). One aligner may handle the 5' end
better (splice-through vs soft-clip), another may find better junctions in the
middle, and a third may handle the 3' poly(A) boundary more cleanly. There's no
reason to take all-or-nothing from a single aligner.

Segment scoring:
  - 5' terminal: penalize soft-clipping (prefer splice-through to real exon)
  - Interior: prefer canonical splice junctions, annotated junction matches,
              fewer mismatches
  - 3' terminal: prefer soft-clipping (clean poly(A) boundary) over alignment
                  into genomic A-tracts that create false junctions

Author: Kevin R. Roy
"""

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set

import pysam
import numpy as np

logger = logging.getLogger(__name__)

# Canonical splice site dinucleotides (reuse from consensus.py)
CANONICAL_5SS = {'GT', 'GC'}
CANONICAL_3SS = {'AG'}


# ============================================================================
# Data structures
# ============================================================================

@dataclass
class CigarEvent:
    """A single CIGAR operation with position tracking."""
    op: int          # CIGAR operation code
    length: int      # Operation length
    q_start: int     # Query position at start (for query-consuming ops)
    q_end: int       # Query position at end (exclusive)
    r_start: int     # Reference position at start (for ref-consuming ops)
    r_end: int       # Reference position at end (exclusive)


@dataclass
class SegmentScore:
    """Score for one aligner's contribution to a segment."""
    aligner: str
    score: float
    n_matches: int = 0
    n_mismatches: int = 0
    n_insertions: int = 0
    n_deletions: int = 0
    n_softclip: int = 0
    n_junctions: int = 0
    n_canonical_junctions: int = 0
    n_annotated_junctions: int = 0
    has_false_3prime_junction: bool = False


@dataclass
class ChimericSegment:
    """A segment of the read between sync points."""
    q_start: int             # Start query position (inclusive)
    q_end: int               # End query position (exclusive)
    position: str            # 'five_prime', 'interior', 'three_prime'
    winning_aligner: str = ""
    scores: Dict[str, SegmentScore] = field(default_factory=dict)
    cigar_events: Dict[str, List[CigarEvent]] = field(default_factory=dict)


@dataclass
class ChimericResult:
    """Result of chimeric consensus selection for a read."""
    read_id: str
    is_chimeric: bool        # True if segments came from different aligners
    segment_winners: List[Tuple[str, str, int, int]]  # [(position, aligner, q_start, q_end)]
    chimeric_cigar: List[Tuple[int, int]]  # Final CIGAR tuples
    chimeric_ref_start: int
    confidence: str          # 'high', 'medium', 'low'
    n_segments: int = 0
    n_aligners_used: int = 0  # How many distinct aligners contributed
    # Per-region statistics
    five_prime_aligner: str = ""
    interior_aligners: List[str] = field(default_factory=list)
    three_prime_aligner: str = ""
    # All scores for stats
    all_segment_scores: List[ChimericSegment] = field(default_factory=list)


# ============================================================================
# Alignment path analysis
# ============================================================================

def cigar_to_events(cigar_tuples: List[Tuple[int, int]], ref_start: int) -> List[CigarEvent]:
    """
    Convert CIGAR tuples to positional events.

    Each event tracks the query and reference position ranges consumed
    by that CIGAR operation.

    Args:
        cigar_tuples: List of (op, length) from pysam.cigartuples
        ref_start: Reference start position (read.reference_start)

    Returns:
        List of CigarEvent objects
    """
    events = []
    qpos = 0
    rpos = ref_start

    for op, length in cigar_tuples:
        if op in (0, 7, 8):  # M, =, X: consume both query and ref
            events.append(CigarEvent(op, length, qpos, qpos + length, rpos, rpos + length))
            qpos += length
            rpos += length
        elif op == 1:  # I: consume query only
            events.append(CigarEvent(1, length, qpos, qpos + length, rpos, rpos))
            qpos += length
        elif op == 2:  # D: consume reference only
            events.append(CigarEvent(2, length, qpos, qpos, rpos, rpos + length))
            rpos += length
        elif op == 3:  # N: splice/skip reference
            events.append(CigarEvent(3, length, qpos, qpos, rpos, rpos + length))
            rpos += length
        elif op == 4:  # S: soft clip (query only)
            events.append(CigarEvent(4, length, qpos, qpos + length, rpos, rpos))
            qpos += length
        elif op == 5:  # H: hard clip (neither)
            events.append(CigarEvent(5, length, qpos, qpos, rpos, rpos))

    return events


def build_query_ref_map(read: pysam.AlignedSegment) -> Dict[int, int]:
    """
    Build a mapping from query position to reference position for aligned bases.

    Only includes positions where the query base is actually aligned to a
    reference base (M, =, X operations). Soft-clipped and inserted bases
    are excluded.

    Args:
        read: pysam AlignedSegment

    Returns:
        Dict mapping query_pos -> ref_pos for all aligned positions
    """
    qr_map = {}
    qpos = 0
    rpos = read.reference_start

    for op, length in read.cigartuples:
        if op in (0, 7, 8):  # M, =, X
            for i in range(length):
                qr_map[qpos + i] = rpos + i
            qpos += length
            rpos += length
        elif op == 1:  # I
            qpos += length
        elif op in (2, 3):  # D, N
            rpos += length
        elif op == 4:  # S
            qpos += length
        elif op == 5:  # H
            pass

    return qr_map


# ============================================================================
# Sync point detection
# ============================================================================

def find_sync_points(
    qr_maps: Dict[str, Dict[int, int]],
) -> List[Tuple[int, int]]:
    """
    Find query positions where ALL aligners map to the same reference position.

    These "sync points" are natural boundaries for chimeric selection: at these
    positions, we can switch from one aligner's CIGAR to another's without any
    discontinuity in the reference coordinate.

    Args:
        qr_maps: Dict mapping aligner_name -> {query_pos: ref_pos}

    Returns:
        Sorted list of (query_pos, ref_pos) sync points
    """
    if not qr_maps or len(qr_maps) < 2:
        return []

    aligner_names = list(qr_maps.keys())

    # Find query positions present in ALL aligners
    common_qpos = set(qr_maps[aligner_names[0]].keys())
    for name in aligner_names[1:]:
        common_qpos &= set(qr_maps[name].keys())

    if not common_qpos:
        return []

    # Filter to positions where all aligners agree on ref_pos
    sync_points = []
    for qpos in sorted(common_qpos):
        ref_positions = {qr_maps[name][qpos] for name in aligner_names}
        if len(ref_positions) == 1:
            sync_points.append((qpos, ref_positions.pop()))

    return sync_points


def compress_sync_runs(
    sync_points: List[Tuple[int, int]],
) -> List[Tuple[int, int, int, int]]:
    """
    Compress consecutive sync points into runs of agreement.

    Returns list of (q_start, q_end, r_start, r_end) for contiguous
    agreement regions. The gaps between these runs are the "disagreement
    segments" where aligners diverge and chimeric selection matters.

    Args:
        sync_points: Sorted list of (query_pos, ref_pos) from find_sync_points

    Returns:
        List of (q_start, q_end_exclusive, r_start, r_end_exclusive) runs
    """
    if not sync_points:
        return []

    runs = []
    run_q_start, run_r_start = sync_points[0]
    prev_q, prev_r = sync_points[0]

    for q, r in sync_points[1:]:
        if q == prev_q + 1 and r == prev_r + 1:
            # Contiguous with previous
            prev_q, prev_r = q, r
        else:
            # Gap: end current run, start new one
            runs.append((run_q_start, prev_q + 1, run_r_start, prev_r + 1))
            run_q_start, run_r_start = q, r
            prev_q, prev_r = q, r

    # Close last run
    runs.append((run_q_start, prev_q + 1, run_r_start, prev_r + 1))
    return runs


# ============================================================================
# Segment identification
# ============================================================================

def identify_segments(
    agreement_runs: List[Tuple[int, int, int, int]],
    read_length: int,
) -> List[Tuple[int, int, str]]:
    """
    Identify disagreement segments between agreement runs.

    Between each pair of agreement runs, there's a disagreement segment
    where aligners diverge. The first and last disagreement segments are
    labeled as terminal (5' or 3' depending on strand), and interior ones
    are labeled as 'interior'.

    Also includes the agreement runs themselves as 'agreement' segments
    (these don't need scoring — any aligner's CIGAR is fine).

    Args:
        agreement_runs: From compress_sync_runs()
        read_length: Total query length

    Returns:
        List of (q_start, q_end, segment_type) where segment_type is
        'five_prime', 'three_prime', 'interior', or 'agreement'
    """
    if not agreement_runs:
        # No agreement at all — treat entire read as one segment
        return [(0, read_length, 'interior')]

    segments = []

    # 5' disagreement (before first agreement)
    first_agree_start = agreement_runs[0][0]
    if first_agree_start > 0:
        segments.append((0, first_agree_start, 'five_prime'))

    # Interleave agreement runs and interior disagreements
    for i, (aq_start, aq_end, _, _) in enumerate(agreement_runs):
        segments.append((aq_start, aq_end, 'agreement'))

        if i < len(agreement_runs) - 1:
            next_start = agreement_runs[i + 1][0]
            if aq_end < next_start:
                segments.append((aq_end, next_start, 'interior'))

    # 3' disagreement (after last agreement)
    last_agree_end = agreement_runs[-1][1]
    if last_agree_end < read_length:
        segments.append((last_agree_end, read_length, 'three_prime'))

    return segments


# ============================================================================
# CIGAR extraction for query sub-ranges
# ============================================================================

def extract_events_for_query_range(
    events: List[CigarEvent],
    q_start: int,
    q_end: int,
) -> List[CigarEvent]:
    """
    Extract and trim CIGAR events for query range [q_start, q_end).

    Handles partial overlap: if an M block spans [40, 80) and we want
    [50, 70), we extract an M of length 20 at the correct ref position.

    For reference-only operations (D, N), includes them if the adjacent
    query position falls within our range.

    Args:
        events: Full list of CigarEvent for an alignment
        q_start: Start query position (inclusive)
        q_end: End query position (exclusive)

    Returns:
        List of CigarEvent objects for the specified query range
    """
    result = []

    for ev in events:
        if ev.op in (0, 7, 8, 1, 4):  # Query-consuming operations
            # Calculate overlap with [q_start, q_end)
            ov_start = max(ev.q_start, q_start)
            ov_end = min(ev.q_end, q_end)

            if ov_start < ov_end:
                trim_left = ov_start - ev.q_start
                new_length = ov_end - ov_start

                if ev.op in (0, 7, 8):  # Also consume reference
                    new_r_start = ev.r_start + trim_left
                    new_r_end = new_r_start + new_length
                else:
                    # I or S: no reference movement
                    new_r_start = ev.r_start
                    new_r_end = ev.r_end

                result.append(CigarEvent(
                    ev.op, new_length,
                    ov_start, ov_end,
                    new_r_start, new_r_end,
                ))

        elif ev.op in (2, 3):  # D, N: reference-only operations
            # Include if the query position (same before and after) is
            # strictly inside our range. At boundaries, the next segment
            # will pick it up.
            if q_start <= ev.q_start < q_end:
                result.append(CigarEvent(
                    ev.op, ev.length,
                    ev.q_start, ev.q_end,
                    ev.r_start, ev.r_end,
                ))

    return result


# ============================================================================
# Segment scoring
# ============================================================================

def score_segment(
    events: List[CigarEvent],
    position: str,
    chrom: str,
    genome: Dict[str, str],
    annotated_junctions: Optional[Set[Tuple[str, int, int]]] = None,
) -> SegmentScore:
    """
    Score a segment from one aligner based on its position in the read.

    Scoring varies by position:
      - five_prime: penalize soft-clipping heavily (-3/base), reward aligned bases
      - interior: reward canonical junctions (+5), annotated junctions (+8),
                  penalize non-canonical junctions (-3), penalize mismatches
      - three_prime: REWARD soft-clipping (+2/base) — clean poly(A) boundary
                     is better than extending into genomic A-tracts. Penalize
                     false 3' junctions from A-tract alignment.

    Args:
        events: CigarEvent list for this segment from one aligner
        position: 'five_prime', 'interior', 'three_prime', or 'agreement'
        chrom: Chromosome name
        genome: Dict mapping chrom -> sequence
        annotated_junctions: Optional set of (chrom, start, end) tuples

    Returns:
        SegmentScore with detailed breakdown
    """
    score_obj = SegmentScore(aligner="", score=0.0)
    score = 0.0

    seq = genome.get(chrom, "")

    for ev in events:
        if ev.op in (0, 7, 8):  # M/=/X: aligned bases
            score_obj.n_matches += ev.length
            # Small reward per aligned base in terminal segments
            if position in ('five_prime', 'three_prime'):
                score += ev.length * 0.5

        elif ev.op == 1:  # Insertion
            score_obj.n_insertions += ev.length
            score -= ev.length * 1.0

        elif ev.op == 2:  # Deletion
            score_obj.n_deletions += ev.length
            score -= ev.length * 0.5

        elif ev.op == 4:  # Soft clip
            score_obj.n_softclip += ev.length

            if position == 'five_prime':
                # 5' soft-clipping is BAD: aligner gave up on mapping these bases
                # which likely correspond to a real exon (especially for DRS 5' ends
                # with short overhangs at splice junctions)
                score -= ev.length * 3.0

            elif position == 'three_prime':
                # 3' soft-clipping is GOOD: clean poly(A) boundary.
                # Aligners that extend through the poly(A) into a genomic A-tract
                # create false junctions and incorrect 3' end positions.
                score += ev.length * 2.0

            else:  # interior
                # Interior soft-clip shouldn't happen, penalize mildly
                score -= ev.length * 1.0

        elif ev.op == 3:  # N = splice junction
            score_obj.n_junctions += 1
            junc_start = ev.r_start
            junc_end = ev.r_end

            # Check canonical splice site motifs
            if seq and junc_start >= 0 and junc_end <= len(seq):
                five_ss = seq[junc_start:junc_start + 2].upper()
                three_ss = seq[junc_end - 2:junc_end].upper()

                if five_ss in CANONICAL_5SS and three_ss in CANONICAL_3SS:
                    score_obj.n_canonical_junctions += 1
                    score += 5.0
                else:
                    score -= 3.0

            # Check annotated junction match
            if annotated_junctions and (chrom, junc_start, junc_end) in annotated_junctions:
                score_obj.n_annotated_junctions += 1
                score += 8.0

            # Detect false 3' junctions (A-tract artifacts)
            if position == 'three_prime' and seq:
                is_false = _is_false_3prime_junction(
                    junc_start, junc_end, seq, min_a_tract=3
                )
                if is_false:
                    score_obj.has_false_3prime_junction = True
                    score -= 15.0  # Heavy penalty for false junction

    score_obj.score = score
    return score_obj


def _is_false_3prime_junction(
    junc_start: int,
    junc_end: int,
    seq: str,
    min_a_tract: int = 3,
) -> bool:
    """
    Check if a junction looks like a poly(A) artifact.

    Pattern: A-tract before the junction AND A-tract after = likely false
    junction where the aligner "jumped" over a genomic region to continue
    aligning to an A-rich stretch.
    """
    if junc_start < 0 or junc_end > len(seq):
        return False

    # Count trailing A's before junction
    before = seq[max(0, junc_start - 10):junc_start]
    trailing_a = 0
    for base in reversed(before):
        if base.upper() == 'A':
            trailing_a += 1
        else:
            break

    # Count leading A's after junction
    after = seq[junc_end:min(len(seq), junc_end + 10)]
    leading_a = 0
    for base in after:
        if base.upper() == 'A':
            leading_a += 1
        else:
            break

    return trailing_a >= min_a_tract and leading_a >= min_a_tract


# ============================================================================
# Chimeric CIGAR construction
# ============================================================================

def build_chimeric_cigar(
    segments: List[ChimericSegment],
    aligner_events: Dict[str, List[CigarEvent]],
    aligner_reads: Dict[str, pysam.AlignedSegment],
) -> Tuple[Optional[int], List[Tuple[int, int]]]:
    """
    Construct a chimeric CIGAR from segments won by different aligners.

    At sync-point boundaries, all aligners have the same (query_pos, ref_pos),
    so we can seamlessly concatenate CIGAR segments from different aligners.

    Reference continuity is enforced as the CIGAR is built:
    - If a reference-consuming event starts ahead of the current reference
      position (gap), an N (splice/skip) bridge is inserted automatically.
    - If a reference regression is detected (new event's r_start < cur_ref),
      the chimeric assembly is invalid. (None, []) is returned as a sentinel
      so the caller can fall back to single-aligner selection.

    Args:
        segments: List of ChimericSegment with winning_aligner set
        aligner_events: Dict of aligner -> full CigarEvent list
        aligner_reads: Dict of aligner -> pysam read

    Returns:
        (reference_start, cigar_tuples) for the chimeric alignment, or
        (None, []) if reference continuity cannot be maintained.
    """
    chimeric_ops: List[Tuple[int, int]] = []
    ref_start: Optional[int] = None
    cur_ref: Optional[int] = None  # tracks current reference position

    for seg in segments:
        winner = seg.winning_aligner
        if not winner:
            continue

        # Get events for this segment from the winning aligner
        seg_events = extract_events_for_query_range(
            aligner_events[winner], seg.q_start, seg.q_end
        )

        for ev in seg_events:
            if ev.op in (0, 7, 8, 2, 3):  # reference-consuming: M/=/X, D, N
                if cur_ref is None:
                    # First reference-consuming event anchors ref_start
                    ref_start = ev.r_start
                    cur_ref = ev.r_start
                elif ev.r_start > cur_ref:
                    # Reference gap between segments: insert N bridge.
                    # This replaces the phantom insertion that would otherwise
                    # appear when a new winner's events start ahead of cur_ref.
                    gap = ev.r_start - cur_ref
                    chimeric_ops.append((3, gap))  # op 3 = N (skip/intron)
                    logger.debug(
                        "Inserted %d-bp N bridge between chimeric segments "
                        "(prev ref=%d, next ref=%d)",
                        gap, cur_ref, ev.r_start,
                    )
                elif ev.r_start < cur_ref:
                    # Reference regression: chimeric assembly is geometrically
                    # invalid (two segments map to overlapping reference regions).
                    # Return sentinel so caller can fall back.
                    logger.debug(
                        "Reference regression in chimeric CIGAR: "
                        "cur_ref=%d, ev.r_start=%d — falling back",
                        cur_ref, ev.r_start,
                    )
                    return None, []

                chimeric_ops.append((ev.op, ev.length))
                cur_ref = ev.r_end  # advance cur_ref by this event's ref span

            else:
                # Query-only events (I=1, S=4): do not advance reference.
                chimeric_ops.append((ev.op, ev.length))

    if ref_start is None:
        return None, []

    # Merge adjacent operations of the same type
    merged = _merge_cigar_ops(chimeric_ops)

    return ref_start, merged


def _validate_chimeric_cigar(
    cigar_tuples: List[Tuple[int, int]],
    read_length: int,
    max_insertion: Optional[int] = None,
    max_intron: int = 10_000,
) -> bool:
    """
    Validate a chimeric CIGAR for biological plausibility.

    Rejects CIGARs with:
    - Any single I operation longer than max_insertion (default: read_length // 4,
      minimum 100 bp). A 947-bp insertion in a 1,869-bp read is a stitching artifact,
      not a real insertion.
    - Any single N operation longer than max_intron (default: 10,000 bp). Introns
      spanning tens of kilobases indicate cross-segment reference jumps that were
      bridged with a phantom N rather than a real splice junction.

    Args:
        cigar_tuples: List of (op, length) CIGAR tuples
        read_length: Full query length (used to compute max_insertion threshold)
        max_insertion: Override for maximum allowed I length
        max_intron: Maximum allowed N length (default 10,000 bp)

    Returns:
        True if the CIGAR passes all checks, False otherwise
    """
    if read_length is None:
        return False  # Can't validate without read length; reject to trigger fallback
    if max_insertion is None:
        max_insertion = max(read_length // 4, 100)

    for op, length in cigar_tuples:
        if op == 1 and length > max_insertion:  # I
            logger.debug(
                "Chimeric CIGAR rejected: %dI exceeds max insertion %d",
                length, max_insertion,
            )
            return False
        if op == 3 and length > max_intron:  # N
            logger.debug(
                "Chimeric CIGAR rejected: %dN exceeds max intron %d",
                length, max_intron,
            )
            return False

    return True


def _merge_cigar_ops(ops: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge adjacent CIGAR operations of the same type."""
    if not ops:
        return []

    merged = [ops[0]]
    for op, length in ops[1:]:
        if op == merged[-1][0]:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))

    # Remove zero-length operations
    return [(op, length) for op, length in merged if length > 0]


# ============================================================================
# Chimeric read construction
# ============================================================================

def build_chimeric_read(
    template_read: pysam.AlignedSegment,
    ref_start: int,
    cigar_tuples: List[Tuple[int, int]],
    chimeric_result: 'ChimericResult',
    header: pysam.AlignmentHeader,
) -> pysam.AlignedSegment:
    """
    Construct a new pysam.AlignedSegment from chimeric selection results.

    Uses the template read for sequence, quality, and basic flags.
    Replaces CIGAR and reference_start with the chimeric version.
    Adds custom tags documenting the chimeric selection.

    Args:
        template_read: Any aligner's read for this query (for seq/qual)
        ref_start: Chimeric reference start position
        cigar_tuples: Chimeric CIGAR tuples
        chimeric_result: ChimericResult with metadata
        header: BAM header for the output file

    Returns:
        New pysam.AlignedSegment ready to write
    """
    out = pysam.AlignedSegment(header)
    out.query_name = template_read.query_name
    out.query_sequence = template_read.query_sequence
    out.query_qualities = template_read.query_qualities
    out.flag = template_read.flag
    out.reference_id = template_read.reference_id
    out.reference_start = ref_start
    out.cigar = cigar_tuples
    out.mapping_quality = template_read.mapping_quality

    # Custom tags
    # XA: aligner(s) used — comma-separated if chimeric
    aligner_set = list(dict.fromkeys(
        w[1] for w in chimeric_result.segment_winners
    ))
    out.set_tag('XA', ','.join(aligner_set))

    # XC: confidence level
    out.set_tag('XC', chimeric_result.confidence)

    # XK: chimeric flag (1 = chimeric, 0 = single aligner)
    out.set_tag('XK', 1 if chimeric_result.is_chimeric else 0)

    # XS: number of segments
    out.set_tag('XS', chimeric_result.n_segments)

    # XU: number of unique aligners used
    out.set_tag('XU', chimeric_result.n_aligners_used)

    # X5: 5' segment aligner
    if chimeric_result.five_prime_aligner:
        out.set_tag('X5', chimeric_result.five_prime_aligner)

    # X3: 3' segment aligner
    if chimeric_result.three_prime_aligner:
        out.set_tag('X3', chimeric_result.three_prime_aligner)

    # XI: interior segment aligners (comma-separated)
    if chimeric_result.interior_aligners:
        out.set_tag('XI', ','.join(chimeric_result.interior_aligners))

    return out


# ============================================================================
# Main chimeric selection
# ============================================================================

def select_best_chimeric(
    aligner_reads: Dict[str, pysam.AlignedSegment],
    genome: Dict[str, str],
    annotated_junctions: Optional[Set[Tuple[str, int, int]]] = None,
    min_sync_fraction: float = 0.05,
    max_intron: int = 10_000,
) -> ChimericResult:
    """
    Select the best chimeric alignment from multiple aligners for a single read.

    Algorithm:
    1. Build query→ref maps for each aligner
    2. Find sync points where all aligners agree on mapping
    3. Compress sync points into agreement runs
    4. Identify disagreement segments between runs
    5. Score each aligner's version of each disagreement segment
    6. Pick the winning aligner per segment
    7. Return chimeric result (may be single-aligner if one wins everything)

    Falls back to simple best-alignment selection when:
    - Only 1 aligner produced an alignment
    - Aligners map to different chromosomes or strands
    - Too few sync points to meaningfully segment

    Args:
        aligner_reads: Dict mapping aligner_name -> pysam.AlignedSegment
        genome: Dict mapping chrom -> sequence
        annotated_junctions: Optional set of (chrom, start, end)
        min_sync_fraction: Minimum fraction of read that must be in sync
                           points for chimeric selection (default 5%)

    Returns:
        ChimericResult with segment winners and chimeric CIGAR
    """
    if not aligner_reads:
        return _empty_result()

    aligner_names = list(aligner_reads.keys())
    read_id = aligner_reads[aligner_names[0]].query_name
    read_length = aligner_reads[aligner_names[0]].query_length

    # ---- Single aligner: no comparison needed ----
    if len(aligner_reads) == 1:
        return _single_aligner_result(aligner_names[0], aligner_reads, read_id)

    # ---- Check preconditions for chimeric selection ----
    # All aligners must map to same chrom and strand
    chroms = {r.reference_name for r in aligner_reads.values()}
    strands = {r.is_reverse for r in aligner_reads.values()}

    if len(chroms) > 1 or len(strands) > 1:
        # Different chrom/strand: fall back to simple scoring
        return _fallback_simple_selection(aligner_reads, genome, annotated_junctions)

    chrom = chroms.pop()
    is_reverse = strands.pop()

    # ---- Build query→ref maps ----
    qr_maps = {}
    all_events = {}
    for name, read in aligner_reads.items():
        qr_maps[name] = build_query_ref_map(read)
        all_events[name] = cigar_to_events(read.cigartuples, read.reference_start)

    # ---- Find sync points ----
    sync_points = find_sync_points(qr_maps)

    # Check minimum sync fraction
    if len(sync_points) < read_length * min_sync_fraction:
        return _fallback_simple_selection(aligner_reads, genome, annotated_junctions)

    # ---- Compress into agreement runs and identify segments ----
    agreement_runs = compress_sync_runs(sync_points)
    raw_segments = identify_segments(agreement_runs, read_length)

    # Adjust terminal labels for strand
    # For minus strand: the 5' end in biological terms is at the RIGHT (high coords)
    # but in query coordinates, query[0] is always the 5' end of the *read* as sequenced.
    # DRS reads are sequenced 3'→5', so query[0] is actually the 3' end biologically.
    # For plus strand: query[0] ≈ 5' end (TSS side)
    # For minus strand: query[0] ≈ 3' end (CPA side)
    if is_reverse:
        # Swap 5' and 3' labels
        adjusted = []
        for q_start, q_end, seg_type in raw_segments:
            if seg_type == 'five_prime':
                adjusted.append((q_start, q_end, 'three_prime'))
            elif seg_type == 'three_prime':
                adjusted.append((q_start, q_end, 'five_prime'))
            else:
                adjusted.append((q_start, q_end, seg_type))
        raw_segments = adjusted

    # ---- Score each aligner per segment ----
    chimeric_segments = []
    for q_start, q_end, seg_type in raw_segments:
        seg = ChimericSegment(
            q_start=q_start, q_end=q_end, position=seg_type
        )

        if seg_type == 'agreement':
            # Any aligner is fine; pick first
            seg.winning_aligner = aligner_names[0]
            seg.cigar_events = {
                name: extract_events_for_query_range(all_events[name], q_start, q_end)
                for name in aligner_names
            }
        else:
            # Score each aligner for this segment
            best_score = float('-inf')
            best_aligner = aligner_names[0]

            for name in aligner_names:
                seg_events = extract_events_for_query_range(
                    all_events[name], q_start, q_end
                )
                seg.cigar_events[name] = seg_events

                score_result = score_segment(
                    seg_events, seg_type, chrom, genome, annotated_junctions
                )
                score_result.aligner = name
                seg.scores[name] = score_result

                if score_result.score > best_score:
                    best_score = score_result.score
                    best_aligner = name

            seg.winning_aligner = best_aligner

        chimeric_segments.append(seg)

    # ---- Build chimeric CIGAR ----
    ref_start, cigar_tuples = build_chimeric_cigar(
        chimeric_segments, all_events, aligner_reads
    )

    # Validate: fall back to simple selection if the chimeric CIGAR is
    # geometrically invalid (reference regression → sentinel None) or
    # biologically implausible (giant phantom insertions or huge N bridges).
    if ref_start is None or not _validate_chimeric_cigar(
        cigar_tuples, read_length, max_intron=max_intron
    ):
        logger.debug(
            "Read %s: chimeric CIGAR failed validation — falling back to "
            "single-aligner selection",
            read_id,
        )
        return _fallback_simple_selection(aligner_reads, genome, annotated_junctions)

    # ---- Collect results ----
    segment_winners = [
        (seg.position, seg.winning_aligner, seg.q_start, seg.q_end)
        for seg in chimeric_segments
    ]

    unique_winners = set(seg.winning_aligner for seg in chimeric_segments)
    is_chimeric = len(unique_winners) > 1

    # Per-region aligner tracking
    five_prime_aligner = ""
    three_prime_aligner = ""
    interior_aligners = []

    for seg in chimeric_segments:
        if seg.position == 'five_prime':
            five_prime_aligner = seg.winning_aligner
        elif seg.position == 'three_prime':
            three_prime_aligner = seg.winning_aligner
        elif seg.position == 'interior':
            interior_aligners.append(seg.winning_aligner)

    # Confidence based on how much of the read is in agreement
    total_agree_bases = sum(
        seg.q_end - seg.q_start
        for seg in chimeric_segments
        if seg.position == 'agreement'
    )
    agree_fraction = total_agree_bases / read_length if read_length > 0 else 0

    if agree_fraction > 0.8:
        confidence = 'high'
    elif agree_fraction > 0.5:
        confidence = 'medium'
    else:
        confidence = 'low'

    return ChimericResult(
        read_id=read_id,
        is_chimeric=is_chimeric,
        segment_winners=segment_winners,
        chimeric_cigar=cigar_tuples,
        chimeric_ref_start=ref_start,
        confidence=confidence,
        n_segments=len(chimeric_segments),
        n_aligners_used=len(unique_winners),
        five_prime_aligner=five_prime_aligner,
        interior_aligners=interior_aligners,
        three_prime_aligner=three_prime_aligner,
        all_segment_scores=chimeric_segments,
    )


# ============================================================================
# Fallback functions
# ============================================================================

def _empty_result() -> ChimericResult:
    """Return an empty result for no-alignment case."""
    return ChimericResult(
        read_id="", is_chimeric=False, segment_winners=[],
        chimeric_cigar=[], chimeric_ref_start=0, confidence='low',
    )


def _single_aligner_result(
    aligner_name: str,
    aligner_reads: Dict[str, pysam.AlignedSegment],
    read_id: str,
) -> ChimericResult:
    """Result when only one aligner produced an alignment."""
    read = aligner_reads[aligner_name]
    return ChimericResult(
        read_id=read_id,
        is_chimeric=False,
        segment_winners=[('whole', aligner_name, 0, read.query_length)],
        chimeric_cigar=list(read.cigartuples),
        chimeric_ref_start=read.reference_start,
        confidence='high',
        n_segments=1,
        n_aligners_used=1,
        five_prime_aligner=aligner_name,
        interior_aligners=[aligner_name],
        three_prime_aligner=aligner_name,
    )


def _fallback_simple_selection(
    aligner_reads: Dict[str, pysam.AlignedSegment],
    genome: Dict[str, str],
    annotated_junctions: Optional[Set[Tuple[str, int, int]]] = None,
) -> ChimericResult:
    """
    Fallback to simple whole-alignment scoring when chimeric isn't possible.

    Uses the existing score_alignment logic: penalize 5' soft-clipping,
    tiebreak on annotated junctions and canonical motifs.
    """
    from .consensus import extract_alignment_info, score_alignment

    if not aligner_reads:
        raise ValueError("_fallback_simple_selection called with empty aligner_reads")

    best_score = float('-inf')
    best_aligner = list(aligner_reads.keys())[0]

    for name, read in aligner_reads.items():
        info = extract_alignment_info(read, name, genome)
        score = score_alignment(info, genome, annotated_junctions)
        if score > best_score:
            best_score = score
            best_aligner = name

    best_read = aligner_reads[best_aligner]
    return ChimericResult(
        read_id=best_read.query_name,
        is_chimeric=False,
        segment_winners=[('whole', best_aligner, 0, best_read.query_length)],
        chimeric_cigar=list(best_read.cigartuples),
        chimeric_ref_start=best_read.reference_start,
        confidence='medium',
        n_segments=1,
        n_aligners_used=1,
        five_prime_aligner=best_aligner,
        interior_aligners=[best_aligner],
        three_prime_aligner=best_aligner,
    )


# ============================================================================
# Statistics aggregation
# ============================================================================

@dataclass
class ChimericStats:
    """Aggregated statistics across all reads for chimeric consensus."""
    total_reads: int = 0
    chimeric_reads: int = 0           # Reads where >1 aligner contributed
    single_aligner_reads: int = 0     # Reads won entirely by one aligner
    fallback_reads: int = 0           # Reads that couldn't use chimeric

    # Per-region aligner contribution counts
    five_prime_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    three_prime_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    interior_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    overall_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))

    # Per-region base counts (how many bases each aligner contributed)
    five_prime_bases_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    three_prime_bases_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    interior_bases_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))

    # Confidence distribution
    confidence_high: int = 0
    confidence_medium: int = 0
    confidence_low: int = 0

    # Segment statistics
    total_segments: int = 0
    avg_segments_per_read: float = 0.0
    avg_agreement_fraction: float = 0.0

    def update(self, result: ChimericResult):
        """Update stats with one read's chimeric result."""
        self.total_reads += 1

        if result.is_chimeric:
            self.chimeric_reads += 1
        elif result.n_segments == 1 and result.segment_winners:
            if result.segment_winners[0][0] == 'whole':
                self.single_aligner_reads += 1
            else:
                self.single_aligner_reads += 1
        else:
            self.single_aligner_reads += 1

        # Confidence
        if result.confidence == 'high':
            self.confidence_high += 1
        elif result.confidence == 'medium':
            self.confidence_medium += 1
        else:
            self.confidence_low += 1

        # Segments
        self.total_segments += result.n_segments

        # Per-region tracking
        for position, aligner, q_start, q_end in result.segment_winners:
            n_bases = q_end - q_start
            self.overall_by_aligner[aligner] += 1

            if position == 'five_prime':
                self.five_prime_by_aligner[aligner] += 1
                self.five_prime_bases_by_aligner[aligner] += n_bases
            elif position == 'three_prime':
                self.three_prime_by_aligner[aligner] += 1
                self.three_prime_bases_by_aligner[aligner] += n_bases
            elif position == 'interior':
                self.interior_by_aligner[aligner] += 1
                self.interior_bases_by_aligner[aligner] += n_bases

    def finalize(self):
        """Compute derived statistics."""
        if self.total_reads > 0:
            self.avg_segments_per_read = self.total_segments / self.total_reads

    def summary_dict(self) -> Dict:
        """Return a summary dict suitable for JSON serialization."""
        self.finalize()
        return {
            'total_reads': self.total_reads,
            'chimeric_reads': self.chimeric_reads,
            'single_aligner_reads': self.single_aligner_reads,
            'fallback_reads': self.fallback_reads,
            'confidence': {
                'high': self.confidence_high,
                'medium': self.confidence_medium,
                'low': self.confidence_low,
            },
            'avg_segments_per_read': round(self.avg_segments_per_read, 2),
            'five_prime_segment_wins': dict(self.five_prime_by_aligner),
            'three_prime_segment_wins': dict(self.three_prime_by_aligner),
            'interior_segment_wins': dict(self.interior_by_aligner),
            'five_prime_bases': dict(self.five_prime_bases_by_aligner),
            'three_prime_bases': dict(self.three_prime_bases_by_aligner),
            'interior_bases': dict(self.interior_bases_by_aligner),
            'overall_segment_wins': dict(self.overall_by_aligner),
        }

    def log_summary(self):
        """Log a human-readable summary."""
        self.finalize()
        logger.info("\n" + "=" * 60)
        logger.info("CHIMERIC CONSENSUS STATISTICS")
        logger.info("=" * 60)
        logger.info(f"  Total reads:         {self.total_reads:,}")
        logger.info(f"  Chimeric reads:      {self.chimeric_reads:,} "
                     f"({100*self.chimeric_reads/max(1,self.total_reads):.1f}%)")
        logger.info(f"  Single-aligner:      {self.single_aligner_reads:,}")
        logger.info(f"  Fallback (simple):   {self.fallback_reads:,}")
        logger.info(f"  Avg segments/read:   {self.avg_segments_per_read:.1f}")
        logger.info("")
        logger.info("  Confidence: high={}, medium={}, low={}".format(
            self.confidence_high, self.confidence_medium, self.confidence_low))
        logger.info("")

        # Per-region tables
        all_aligners = sorted(set(
            list(self.five_prime_by_aligner.keys()) +
            list(self.three_prime_by_aligner.keys()) +
            list(self.interior_by_aligner.keys())
        ))

        if all_aligners:
            logger.info("  Segment wins by region:")
            header = f"    {'Aligner':<15} {'5-prime':>10} {'Interior':>10} {'3-prime':>10}"
            logger.info(header)
            logger.info("    " + "-" * 47)
            for a in all_aligners:
                logger.info(f"    {a:<15} {self.five_prime_by_aligner.get(a,0):>10,} "
                             f"{self.interior_by_aligner.get(a,0):>10,} "
                             f"{self.three_prime_by_aligner.get(a,0):>10,}")

            logger.info("")
            logger.info("  Bases contributed by region:")
            header = f"    {'Aligner':<15} {'5-prime':>10} {'Interior':>10} {'3-prime':>10}"
            logger.info(header)
            logger.info("    " + "-" * 47)
            for a in all_aligners:
                logger.info(f"    {a:<15} {self.five_prime_bases_by_aligner.get(a,0):>10,} "
                             f"{self.interior_bases_by_aligner.get(a,0):>10,} "
                             f"{self.three_prime_bases_by_aligner.get(a,0):>10,}")

        logger.info("=" * 60)
