"""
Best-alignment selection across multiple aligners for a single read.

Takes a ``Dict[aligner_name, AlignmentInfo]`` produced by
``extract.extract_alignment_info`` (one per aligner), scores each via
``scoring.score_alignment``, and returns a ``ConsensusResult`` naming
the winning aligner plus tiebreaker / confidence metadata.

The streaming orchestration that drives this per-read selection across a
full multi-aligner BAM set lives in ``consensus.run_consensus_selection``.
"""

import bisect
from typing import Dict, Optional, Set, Tuple

from .extract import AlignmentInfo, ConsensusResult
from .scoring import score_alignment


# Memoized per-chrom position index over an annotated-junction set. Lets the
# per-read 5' soft-clip rescue pool be built by locus window instead of a full
# chrom-wide scan, turning O(reads × junctions_on_chrom) — which hangs at
# genome-wide short-read scale — into O(reads × log n). Keyed by id() of the
# set (stable for a consensus run, the same object is threaded through every
# per-read call); the stored ref guards against id reuse after GC.
_ANNOT_INDEX_CACHE: dict = {}


def _annot_index(annotated_junctions):
    key = id(annotated_junctions)
    cached = _ANNOT_INDEX_CACHE.get(key)
    if cached is not None and cached[0] is annotated_junctions:
        return cached[1]
    by_chrom: Dict[str, Tuple[list, list]] = {}
    for j in annotated_junctions:
        chrom = j[0]
        starts, ends = by_chrom.setdefault(chrom, ([], []))
        starts.append((j[1], j))
        ends.append((j[2], j))
    index = {}
    for chrom, (by_start, by_end) in by_chrom.items():
        by_start.sort(key=lambda x: x[0])
        by_end.sort(key=lambda x: x[0])
        index[chrom] = (
            [k for k, _ in by_start], [v for _, v in by_start],
            [k for k, _ in by_end], [v for _, v in by_end],
        )
    _ANNOT_INDEX_CACHE[key] = (annotated_junctions, index)
    return index


def _query_annot_window(annotated_junctions, chrom, lo, hi):
    """Annotated junctions on ``chrom`` whose start OR end lies in [lo, hi].

    A superset of everything ``_rescue_5prime_softclip`` can match (it filters
    to ≤ ``search_window_bp`` of the read's 5' boundary), so restricting the
    annotated pool to this window leaves the scoring result byte-identical.
    """
    idx = _annot_index(annotated_junctions)
    entry = idx.get(chrom)
    if entry is None:
        return set()
    skeys, svals, ekeys, evals = entry
    out: Set[Tuple[str, int, int, str]] = set()
    i = bisect.bisect_left(skeys, lo)
    while i < len(skeys) and skeys[i] <= hi:
        out.add(svals[i]); i += 1
    i = bisect.bisect_left(ekeys, lo)
    while i < len(ekeys) and ekeys[i] <= hi:
        out.add(evals[i]); i += 1
    return out


def select_best_alignment(
    alignments: Dict[str, AlignmentInfo],
    genome: Dict[str, str],
    annotated_junctions: Optional[Set[Tuple[str, int, int, str]]] = None,
    tiebreak: str = 'rectify',
) -> ConsensusResult:
    """
    Select the best alignment from multiple aligners for a single read.

    Args:
        alignments: Dict mapping aligner name to AlignmentInfo
        genome: Dict mapping chrom to sequence
        annotated_junctions: Optional set of (chrom, start, end, strand) for annotated junctions
        tiebreak: Tiebreak ordering applied when multiple aligners share the top
            junction score. Two orders are supported:

            - ``'rectify'`` (default, long-read end-correction mode): prefer the
              alignment whose corrected 3' end agrees with the majority, then more
              annotated junctions, then more canonical (GT/AG) splice sites. This
              is the established long-read order and is kept byte-identical.
            - ``'compass'`` (short-read splice-junction mode): the COMPASS published
              order ``ungapped > gapped > annotated > shorter-intron`` — prefer the
              alignment with no introns, then more annotated junctions, then the
              shorter total intron length. On a score tie this conservatively favors
              the simplest explanation, which is what lets the short-read panel flag
              suspected false junctions as artifacts. The function default stays
              ``'rectify'`` so existing long-read callers are unaffected; the
              short-read entrypoint passes ``tiebreak='compass'`` explicitly.

    Returns:
        ConsensusResult with best alignment selected
    """
    if tiebreak not in ('rectify', 'compass'):
        raise ValueError(
            f"tiebreak must be 'rectify' or 'compass', got {tiebreak!r}"
        )
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
        # Locus-scoped annotated pool: a window-superset of the rescue's ≤300bp
        # reach around this read, so it is byte-identical to a full chrom scan
        # but O(log n) per read instead of O(n) (the chrom-wide scan hangs at
        # genome-wide short-read scale — ~tens of thousands of junctions/chrom
        # × ~50k reads).
        _pos = []
        for _a in alignments.values():
            _rs = getattr(_a, 'reference_start', None)
            _re = getattr(_a, 'reference_end', None)
            if _rs is not None:
                _pos.append(_rs)
            if _re is not None:
                _pos.append(_re)
        if _pos:
            candidate_junctions |= _query_annot_window(
                annotated_junctions, chrom_for_read, min(_pos) - 350, max(_pos) + 350)
        else:
            # No positional info (shouldn't happen) — fall back to chrom scan.
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
        # Number of this alignment's junctions that match an annotated intron.
        # Shared by both tiebreak orders.
        def _n_annotated(aligner_name):
            a = alignments[aligner_name]
            if not (annotated_junctions and a.junctions):
                return 0
            return sum(
                1 for junc in a.junctions
                if (a.chrom, junc[0], junc[1], a.strand) in annotated_junctions
            )

        if tiebreak == 'compass':
            # COMPASS published order: ungapped > gapped > annotated > shorter-intron.
            # Encoded as a tuple maximized by max():
            #   is_ungapped   1 if the alignment has no introns (prefer ungapped)
            #   n_annotated   more annotated junctions preferred
            #   -intron_len   shorter total intron length preferred (negate to maximize)
            #   canonical     more canonical GT/AG splice sites (scientific tiebreaker,
            #                 not part of the published order but harmless and improves
            #                 determinism vs. raw name fallback)
            #   aligner_name  deterministic final fallback
            def _tiebreak_key(aligner_name):
                a = alignments[aligner_name]
                is_ungapped = 0 if a.junctions else 1
                total_intron_len = sum(end - start for start, end in a.junctions)
                return (
                    is_ungapped,
                    _n_annotated(aligner_name),
                    -total_intron_len,
                    a.canonical_count,
                    aligner_name,
                )
        else:
            # RECTIFY long-read order (kept byte-identical):
            # Tiebreaker 1: prefer alignment whose corrected 3' end agrees with majority
            all_corrected = [a.corrected_3prime for a in alignments.values()
                             if a.corrected_3prime is not None]
            def _count_3prime_agreement(aligner_name):
                pos = alignments[aligner_name].corrected_3prime
                return sum(1 for p in all_corrected if p == pos) if pos is not None else 0

            # Tiebreaker 2: prefer alignment with more annotated junctions
            # Tiebreaker 3: prefer alignment with more canonical splice sites (GT/AG)
            def _tiebreak_key(aligner_name):
                return (
                    _count_3prime_agreement(aligner_name),
                    _n_annotated(aligner_name),
                    alignments[aligner_name].canonical_count,
                )

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
        tied_aligners=tied_aligners,
        confidence=confidence,
        was_5prime_rescued=was_rescued,
        false_junction_removed=False,  # Not tracked; walk back correction handles this
        all_alignments=alignments,
    )
