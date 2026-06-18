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

# Splice junctions are coordinate-ambiguous: when the bases flanking the donor
# and acceptor repeat, the SAME spliced product can be written with the intron
# shifted by a few bp. A junction (intron [start,end)) slides left by 1 with an
# identical spliced product iff seq[start-1] == seq[end-1] (the base leaving the
# 3' end of exon1 equals the base leaving the 3' end of the intron); it slides
# right by 1 iff seq[start] == seq[end]. This is the same base-equality slide
# used by junction_refiner._apply_junction_replacement's pure-slide fast path
# and the _l_amb/_r_amb windows in splice_aware_5prime. Comparing junctions by
# exact coordinate (the old behaviour) charges an aligner that placed the SAME
# junction one bp over as a mismatch — so all junction comparisons below
# normalize to the LEFTMOST equivalent coordinate first.
_JUNCTION_AMBIGUITY_MAX_SHIFT = 15  # matches splice_aware_5prime.MAX_SS_SHIFT
# The annotated bonus is GATED on read support: only trust an annotated match
# when the read's own alignment anchors the junction with at least this many
# contiguous matched bases on the SHORTER flank. Otherwise the catalog would
# override read evidence (forcing a short-overhang junction onto an annotated
# site and suppressing genuine nearby novel junctions). Matches the human
# anchor-gate K (corrected_consensus _MIN_JUNCTION_ANCHOR / K=10).
_ANNOTATED_SUPPORT_MIN_ANCHOR = 10


def normalize_junction(start: int, end: int, seq: str,
                       max_shift: int = _JUNCTION_AMBIGUITY_MAX_SHIFT) -> Tuple[int, int]:
    """Left-normalize a junction to the leftmost ambiguity-equivalent coordinate.

    Slides the intron left while ``seq[start-1] == seq[end-1]`` (the base-equality
    slide that preserves the spliced product). Two junctions that are
    ambiguity-equivalent normalize to the same (start, end), so a normalized set
    membership test is an ambiguity-aware match. Bounded by ``max_shift``.
    """
    shifts = 0
    while (shifts < max_shift and start > 0 and end > 0
           and start - 1 < len(seq) and end - 1 < len(seq)
           and seq[start - 1].upper() == seq[end - 1].upper()):
        start -= 1
        end -= 1
        shifts += 1
    return start, end


def junction_ambiguity_window(start: int, end: int, seq: str,
                              max_shift: int = _JUNCTION_AMBIGUITY_MAX_SHIFT) -> Tuple[int, int]:
    """Return (l_amb, r_amb): how many bp the junction may slide left / right
    while preserving the spliced product (the sequence-ambiguity window)."""
    l_amb = 0
    while (l_amb < max_shift and start - 1 - l_amb >= 0 and end - 1 - l_amb >= 0
           and seq[start - 1 - l_amb].upper() == seq[end - 1 - l_amb].upper()):
        l_amb += 1
    r_amb = 0
    while (r_amb < max_shift and start + r_amb < len(seq) and end + r_amb < len(seq)
           and seq[start + r_amb].upper() == seq[end + r_amb].upper()):
        r_amb += 1
    return l_amb, r_amb


def _full_junction_anchor(full_events, junc_start: int, junc_end: int):
    """True flanking anchor (min of the two contiguous matched runs) for a
    junction, measured from the aligner's FULL alignment — NOT the per-segment
    clipped events. Segmentation puts the agreed flanking exon into agreement
    runs and bounds the interior disagreement segment tightly around the
    junction, so a segment-local anchor under-counts the real read support and
    the support gate would fail spuriously. Returns None if the junction isn't
    found in full_events (caller falls back to segment-local)."""
    if not full_events:
        return None
    for i, ev in enumerate(full_events):
        if ev.op == 3 and ev.r_start == junc_start and ev.r_end == junc_end:
            left = (full_events[i - 1].length
                    if i > 0 and full_events[i - 1].op in (0, 7, 8) else 0)
            right = (full_events[i + 1].length
                     if i + 1 < len(full_events) and full_events[i + 1].op in (0, 7, 8) else 0)
            return min(left, right)
    return None


def _canonical_within_window(start: int, end: int, seq: str,
                             l_amb: int, r_amb: int) -> bool:
    """True if ANY ambiguity-equivalent placement of the junction yields a
    canonical (GT/GC..AG) motif. The biologically real junction is the canonical
    placement, so an aligner that landed one bp off a canonical site within the
    window still gets canonical credit."""
    for s in range(-l_amb, r_amb + 1):
        js, je = start + s, end + s
        if js < 0 or je > len(seq) or je - 2 < 0:
            continue
        five_ss = seq[js:js + 2].upper()
        three_ss = seq[je - 2:je].upper()
        if five_ss in CANONICAL_5SS and three_ss in CANONICAL_3SS:
            return True
    return False


# Cache the normalized annotation set by identity: load_annotated_junctions
# returns one set per run reused across every read, so normalizing (and 4→3
# tuple projection) once instead of per-read keeps the hot path cheap.
_ANNOT_NORM_CACHE: Dict[int, Tuple[object, Set[Tuple[str, int, int]]]] = {}


def _normalized_annotation_set(annotated_junctions, genome):
    """Left-normalize the annotation set to leftmost ambiguity-equivalent
    coordinates (and project any 4-tuples down to (chrom,start,end)), cached by
    object identity so it runs ONCE per run, not once per read."""
    if not annotated_junctions:
        return annotated_junctions
    key = id(annotated_junctions)
    hit = _ANNOT_NORM_CACHE.get(key)
    if hit is not None and hit[0] is annotated_junctions:
        return hit[1]
    norm = set()
    for j in annotated_junctions:
        chrom, s, e = j[0], j[1], j[2]
        seq = genome.get(chrom, "")
        if seq:
            s, e = normalize_junction(s, e, seq)
        norm.add((chrom, s, e))
    _ANNOT_NORM_CACHE.clear()  # single annotation set per run; bound the cache
    _ANNOT_NORM_CACHE[key] = (annotated_junctions, norm)
    return norm


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
    # Canonical (GT-AG within the ambiguity window) but NOT in the annotation set:
    # a de-novo novel-junction candidate. Tracked so the annotated-bonus tuning
    # can't silently suppress discovery (a novel canonical junction must not be
    # scored as "wrong" merely for being absent from the catalog).
    n_novel_canonical_junctions: int = 0
    # Annotated match that FAILED the read-support gate (short anchor) — the
    # catalog wanted it but the read didn't earn it; bonus withheld.
    n_annotated_unsupported: int = 0
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
    # True when this result came from a fallback path (no usable sync points /
    # CIGAR build failed) rather than genuine multi-aligner stitching.  Used by
    # ChimericStats to separate real segment wins from degenerate single-aligner
    # pass-throughs so a measurement can't be confounded by fallback rate.
    is_fallback: bool = False


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
    min_aligners: Optional[int] = None,
) -> List[Tuple[int, int]]:
    """
    Find query positions where aligners agree on a reference position.

    These "sync points" are natural boundaries for chimeric selection: at these
    positions, we can switch from one aligner's CIGAR to another's without any
    discontinuity in the reference coordinate.

    Args:
        qr_maps: Dict mapping aligner_name -> {query_pos: ref_pos}
        min_aligners: Seam stringency.
            * None (default) → require ALL aligners present at the position to
              agree (the original behaviour; byte-identical).
            * int k → a sync point forms where ≥ k aligners agree on the SAME
              ref position. MUST be a STRICT majority (k > n/2) so the agreed
              ref is unique — otherwise two disjoint groups could each hit k
              and the seam ref would be ambiguous (phantom-bridge risk). Values
              that are not a strict majority are raised to one.

    Returns:
        Sorted list of (query_pos, ref_pos) sync points. Under k-of-n, ref_pos
        is the majority-agreed reference position.
    """
    if not qr_maps or len(qr_maps) < 2:
        return []

    aligner_names = list(qr_maps.keys())
    n = len(aligner_names)

    # ---- ALL-agree mode (default): unchanged, byte-identical ----
    if min_aligners is None:
        common_qpos = set(qr_maps[aligner_names[0]].keys())
        for name in aligner_names[1:]:
            common_qpos &= set(qr_maps[name].keys())
        if not common_qpos:
            return []
        sync_points = []
        for qpos in sorted(common_qpos):
            ref_positions = {qr_maps[name][qpos] for name in aligner_names}
            if len(ref_positions) == 1:
                sync_points.append((qpos, ref_positions.pop()))
        return sync_points

    # ---- k-of-n mode: strict-majority agreement on a unique ref ----
    # Enforce strict majority so the agreed ref position is unique.
    k = max(int(min_aligners), n // 2 + 1)
    if k > n:
        return []

    # Positions present in at least k aligners are the only candidates.
    qpos_counts: Dict[int, int] = defaultdict(int)
    for name in aligner_names:
        for qpos in qr_maps[name]:
            qpos_counts[qpos] += 1

    sync_points = []
    for qpos in sorted(q for q, c in qpos_counts.items() if c >= k):
        ref_votes: Dict[int, int] = defaultdict(int)
        for name in aligner_names:
            r = qr_maps[name].get(qpos)
            if r is not None:
                ref_votes[r] += 1
        # Strict majority (k > n/2) guarantees at most one ref can reach k.
        best_ref, best_votes = max(ref_votes.items(), key=lambda kv: kv[1])
        if best_votes >= k:
            sync_points.append((qpos, best_ref))

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
    annotated_min_anchor: int = _ANNOTATED_SUPPORT_MIN_ANCHOR,
    aligner_full_events: Optional[List[CigarEvent]] = None,
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

    Junction scoring is AMBIGUITY-AWARE: the canonical-motif check slides within
    the junction's sequence-ambiguity window, and the annotated match normalizes
    both the read junction and the catalog to the leftmost equivalent coordinate
    (so the same junction written one bp over is not charged as a mismatch). The
    annotated +8 bonus is GATED on read support (a minimum flanking anchor) so
    the catalog breaks ties rather than overriding read evidence — a junction the
    aligner FORCED onto a catalogued site with a short overhang earns no bonus,
    and a genuine novel canonical junction nearby is not suppressed.

    Args:
        events: CigarEvent list for this segment from one aligner
        position: 'five_prime', 'interior', 'three_prime', or 'agreement'
        chrom: Chromosome name
        genome: Dict mapping chrom -> sequence
        annotated_junctions: Optional set of (chrom, start, end) tuples. MUST be
            left-normalized (see normalize_junction); select_best_chimeric does
            this once per run.
        annotated_min_anchor: minimum contiguous matched bases on the SHORTER
            flank of a junction for its annotated bonus to apply.

    Returns:
        SegmentScore with detailed breakdown
    """
    score_obj = SegmentScore(aligner="", score=0.0)
    score = 0.0

    seq = genome.get(chrom, "")

    for idx, ev in enumerate(events):
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

            # Flanking anchor = contiguous matched bases immediately on each side
            # (read support for the junction). Prefer the FULL alignment's anchor
            # — the per-segment events are clipped tightly around the disagreement
            # so a segment-local anchor under-counts real support and the gate
            # would fail spuriously. Fall back to segment-local if the junction
            # isn't located in the full events.
            min_anchor = _full_junction_anchor(aligner_full_events, junc_start, junc_end)
            if min_anchor is None:
                left_anchor = (events[idx - 1].length
                               if idx > 0 and events[idx - 1].op in (0, 7, 8) else 0)
                right_anchor = (events[idx + 1].length
                                if idx + 1 < len(events) and events[idx + 1].op in (0, 7, 8) else 0)
                min_anchor = min(left_anchor, right_anchor)

            # Canonical motif, ambiguity-aware (any equivalent placement GT-AG).
            is_canonical = False
            if seq and junc_start >= 0 and junc_end <= len(seq):
                l_amb, r_amb = junction_ambiguity_window(junc_start, junc_end, seq)
                is_canonical = _canonical_within_window(
                    junc_start, junc_end, seq, l_amb, r_amb
                )
                if is_canonical:
                    score_obj.n_canonical_junctions += 1
                    score += 5.0
                else:
                    score -= 3.0

            # Annotated match, ambiguity-aware: normalize the read junction to the
            # leftmost equivalent coordinate and test against the (pre-normalized)
            # catalog. Bonus is GATED on read support.
            is_annotated = False
            if annotated_junctions and seq:
                njs, nje = normalize_junction(junc_start, junc_end, seq)
                is_annotated = (chrom, njs, nje) in annotated_junctions
            if is_annotated:
                if min_anchor >= annotated_min_anchor:
                    score_obj.n_annotated_junctions += 1
                    score += 8.0
                else:
                    # Catalog wanted it but the read didn't anchor it — withhold
                    # the bonus so annotation can't override weak read evidence.
                    score_obj.n_annotated_unsupported += 1
            elif is_canonical:
                # Canonical but not catalogued = de-novo novel-junction candidate.
                # Must NOT be treated as wrong; tracked so the bonus tuning can be
                # audited for discovery suppression.
                score_obj.n_novel_canonical_junctions += 1

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
                    # Reference gap between segments.
                    gap = ev.r_start - cur_ref
                    if ev.op in (2, 3):
                        # D or N event at a segment boundary: the gap is
                        # contiguous with this reference-only operation.
                        # Extend the event to cover the gap rather than
                        # inserting a separate N bridge, which would either
                        # duplicate the op (N bridge + N event) or produce
                        # a phantom N before a D event (malformed CIGAR).
                        ev = CigarEvent(
                            ev.op, ev.length + gap,
                            ev.q_start, ev.q_end,
                            cur_ref, ev.r_end,
                        )
                        logger.debug(
                            "Extended %d-bp gap into %s event at chimeric "
                            "boundary (prev ref=%d, next ref=%d)",
                            gap, "D" if ev.op == 2 else "N", cur_ref, ev.r_end,
                        )
                    else:
                        # M/=/X event: gap requires a separate bridge op.
                        # Small gaps (≤ MAX_STITCH_DELETION) at chimeric stitch
                        # points are almost certainly artefactual deletions, not
                        # genuine intron skips.  Use D (op 2) so they render as
                        # deletions in IGV rather than as phantom introns.
                        _MAX_STITCH_DELETION = 10
                        if gap <= _MAX_STITCH_DELETION:
                            chimeric_ops.append((2, gap))  # D (deletion)
                            logger.debug(
                                "Inserted %d-bp D bridge at chimeric stitch "
                                "(prev ref=%d, next ref=%d)",
                                gap, cur_ref, ev.r_start,
                            )
                        else:
                            chimeric_ops.append((3, gap))  # N (genuine intron skip)
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

    # Custom tags — lowercase second-letter to avoid colliding with the
    # cDNA pipeline's X[upper] tags (XU=UMI, XC=cluster_size, XA=tail_len,
    # XS=sense/antisense, XK=3' pre-trim length, XI=isoform_id, ...).
    # Xa: aligner(s) used — comma-separated if chimeric
    aligner_set = list(dict.fromkeys(
        w[1] for w in chimeric_result.segment_winners
    ))
    out.set_tag('Xa', ','.join(aligner_set))

    # Xc: confidence level
    out.set_tag('Xc', chimeric_result.confidence)

    # Xz: chimeric flag (1 = chimeric, 0 = single aligner)
    out.set_tag('Xz', 1 if chimeric_result.is_chimeric else 0)

    # Xg: number of segments
    out.set_tag('Xg', chimeric_result.n_segments)

    # Xm: number of unique aligners used (multi-aligner / "merged")
    out.set_tag('Xm', chimeric_result.n_aligners_used)

    # Xq: 5' segment aligner
    if chimeric_result.five_prime_aligner:
        out.set_tag('Xq', chimeric_result.five_prime_aligner)

    # Xw: 3' segment aligner
    if chimeric_result.three_prime_aligner:
        out.set_tag('Xw', chimeric_result.three_prime_aligner)

    # Xy: interior segment aligners (comma-separated)
    if chimeric_result.interior_aligners:
        out.set_tag('Xy', ','.join(chimeric_result.interior_aligners))

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

    # ---- Normalize annotated junctions for the segment scorer ----
    # score_segment tests 3-tuple (chrom, start, end) membership, but
    # load_annotated_junctions (and the fallback / winner-take-all scorer in
    # scoring.py) use 4-tuples (chrom, start, end, strand). Passing the 4-tuple
    # set straight through meant the membership test NEVER matched, so the
    # annotated +8 reward was silently dead in the chimeric path for every
    # organism. Project to 3-tuples here (strand-agnostic — identical intron
    # coordinates on opposite strands are vanishingly rare). The original set is
    # still handed to _fallback_simple_selection, whose scorer unpacks strand.
    annotated_for_segments = _normalized_annotation_set(annotated_junctions, genome)

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
                    seg_events, seg_type, chrom, genome, annotated_for_segments,
                    aligner_full_events=all_events[name],
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
        is_fallback=True,
    )


# ============================================================================
# Statistics aggregation
# ============================================================================

def _classify_term(winner: 'SegmentScore', loser: 'SegmentScore') -> str:
    """Classify the term that separated `winner` from `loser` on one contested
    segment.  Checked in priority order so the dominant motif/structure term is
    attributed before the residue-level fallback:

      annotated    — winner placed more annotated junctions than the loser
      canonical    — winner placed more canonical (GT-AG) junctions
      false_3prime — loser carried a spurious 3' junction the winner did not
      softclip     — loser soft-clipped more (the "went for it" overhang term)
      indel        — loser carried more insertions+deletions
      residue      — decided on matches/mismatches only (no structural term)

    This lets a downstream measurement see whether GMAP's segment wins/losses
    ride on the annotated (+8) term — which is DEAD when annotated_junctions is
    empty — versus the canonical (+5/-3) term, separating signal from the
    load_annotated_junctions confound.
    """
    if winner.n_annotated_junctions > loser.n_annotated_junctions:
        return 'annotated'
    if winner.n_canonical_junctions > loser.n_canonical_junctions:
        return 'canonical'
    if loser.has_false_3prime_junction and not winner.has_false_3prime_junction:
        return 'false_3prime'
    if loser.n_softclip > winner.n_softclip:
        return 'softclip'
    if (loser.n_insertions + loser.n_deletions) > (winner.n_insertions + winner.n_deletions):
        return 'indel'
    return 'residue'


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

    # ---- Loss/win-reason instrumentation (the measurement of WHY a segment
    # went the way it did).  Without this, a GMAP-vs-rest segment comparison
    # repeats the annotation-off confound: if `annotated_junctions` is empty
    # (the load_annotated_junctions exon-GTF bug), the annotated (+8) term is
    # dead and every junction decision collapses onto the canonical (+5/-3)
    # term — so we must record which term actually decided each contested
    # segment, not just who won.
    # `decisive_term`: across all CONTESTED (>1 distinct score) scored
    # segments, what term separated winner from runner-up.
    decisive_term: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    # Per-aligner: segments it competed in (was scored) and won.
    segments_competed_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    segments_won_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    # Per-aligner loss reasons: when this aligner LOST a contested segment,
    # which winner-favouring term beat it (annotated / canonical / false_3prime
    # / softclip / other).
    loss_reason_by_aligner: Dict[str, Dict[str, int]] = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    # Per-aligner win reasons: when this aligner WON a contested segment, what
    # the winning margin rode on (same term vocabulary).
    win_reason_by_aligner: Dict[str, Dict[str, int]] = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    # ---- Discovery accounting (guards against the annotated bonus silently
    # suppressing novel junctions). novel_canonical_won: novel canonical
    # junctions the WINNER carried into consensus, per aligner. The key number
    # for the GMAP question — annotated_loss_with_own_novel: segments an aligner
    # LOST on the annotated term while ITS OWN losing segment carried a
    # canonical-but-uncatalogued junction (i.e. a novel candidate the catalog
    # bonus outscored). A large value = the annotated win/loss ratio is
    # penalizing discovery, not measuring correctness.
    novel_canonical_won_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    annotated_loss_with_own_novel_by_aligner: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    # Support-gate diagnostic: of annotated junctions in the consensus winner,
    # how many earned the bonus (anchor ok) vs were withheld (short anchor). If
    # `unsupported` swamps `supported`, the gate is eating real junctions (likely
    # an anchor-truncation artifact), not filtering genuine short overhangs.
    annotated_supported: int = 0
    annotated_unsupported: int = 0

    def update(self, result: ChimericResult):
        """Update stats with one read's chimeric result."""
        self.total_reads += 1

        if result.is_fallback:
            self.fallback_reads += 1
        elif result.is_chimeric:
            self.chimeric_reads += 1
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

        # ---- Loss/win-reason: only contested scored segments carry this.
        # Agreement segments and single-candidate segments are skipped.
        for seg in result.all_segment_scores:
            scores = getattr(seg, 'scores', None)
            if not scores or len(scores) < 2:
                continue
            winner = seg.winning_aligner
            if winner not in scores:
                continue
            distinct = {round(s.score, 6) for s in scores.values()}
            contested = len(distinct) > 1
            for aligner in scores:
                self.segments_competed_by_aligner[aligner] += 1
            self.segments_won_by_aligner[winner] += 1
            if not contested:
                continue
            w = scores[winner]
            # Novel canonical junctions the winner carried into consensus.
            self.novel_canonical_won_by_aligner[winner] += getattr(
                w, 'n_novel_canonical_junctions', 0)
            # Support-gate diagnostic (consensus winner's annotated junctions).
            self.annotated_supported += getattr(w, 'n_annotated_junctions', 0)
            self.annotated_unsupported += getattr(w, 'n_annotated_unsupported', 0)
            for aligner, s in scores.items():
                if aligner == winner:
                    continue
                reason = _classify_term(w, s)
                self.loss_reason_by_aligner[aligner][reason] += 1
                self.win_reason_by_aligner[winner][reason] += 1
                self.decisive_term[reason] += 1
                # Discovery-suppression flag: lost on the annotated term while
                # ITS OWN losing segment carried a novel canonical junction.
                if reason == 'annotated' and getattr(
                        s, 'n_novel_canonical_junctions', 0) > 0:
                    self.annotated_loss_with_own_novel_by_aligner[aligner] += 1

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
            'decisive_term': dict(self.decisive_term),
            'segments_competed_by_aligner': dict(self.segments_competed_by_aligner),
            'segments_won_by_aligner': dict(self.segments_won_by_aligner),
            'loss_reason_by_aligner': {
                a: dict(d) for a, d in self.loss_reason_by_aligner.items()
            },
            'win_reason_by_aligner': {
                a: dict(d) for a, d in self.win_reason_by_aligner.items()
            },
            'novel_canonical_won_by_aligner': dict(self.novel_canonical_won_by_aligner),
            'annotated_loss_with_own_novel_by_aligner': dict(
                self.annotated_loss_with_own_novel_by_aligner),
            'annotated_supported': self.annotated_supported,
            'annotated_unsupported': self.annotated_unsupported,
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

        # ---- Decisive-term breakdown (the annotation-off confound check) ----
        if self.decisive_term:
            total_contested = sum(self.decisive_term.values())
            logger.info("")
            logger.info(f"  Contested-segment decisions (n={total_contested:,}) "
                         f"— which term separated winner from runner-up:")
            for term in ('annotated', 'canonical', 'false_3prime',
                         'softclip', 'indel', 'residue'):
                n = self.decisive_term.get(term, 0)
                if n:
                    logger.info(f"    {term:<14} {n:>10,} "
                                 f"({100*n/max(1,total_contested):.1f}%)")
            if self.decisive_term.get('annotated', 0) == 0:
                logger.info("    NOTE: 0 decisions on the annotated term — if "
                             "annotated_junctions is empty (exon-GTF bug), the "
                             "+8 reward is dead and junction calls ride on the "
                             "canonical term only.")

        # ---- Per-aligner win/loss reasons (GMAP signal vs gaming) ----
        reason_aligners = sorted(set(
            list(self.win_reason_by_aligner.keys()) +
            list(self.loss_reason_by_aligner.keys())
        ))
        if reason_aligners:
            logger.info("")
            logger.info("  Per-aligner contested win/loss reasons:")
            for a in reason_aligners:
                won = self.segments_won_by_aligner.get(a, 0)
                competed = self.segments_competed_by_aligner.get(a, 0)
                wr = dict(self.win_reason_by_aligner.get(a, {}))
                lr = dict(self.loss_reason_by_aligner.get(a, {}))
                logger.info(f"    {a:<12} won {won:,}/{competed:,} contested  "
                             f"win_terms={wr}  loss_terms={lr}")

        logger.info("=" * 60)
