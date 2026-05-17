"""
Best-alignment selection across multiple aligners for a single read.

Takes a ``Dict[aligner_name, AlignmentInfo]`` produced by
``extract.extract_alignment_info`` (one per aligner), scores each via
``scoring.score_alignment``, and returns a ``ConsensusResult`` naming
the winning aligner plus tiebreaker / confidence metadata.

The streaming orchestration that drives this per-read selection across a
full multi-aligner BAM set lives in ``consensus.run_consensus_selection``.
"""

from typing import Dict, Optional, Set, Tuple

from .extract import AlignmentInfo, ConsensusResult
from .scoring import score_alignment


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
        tied_aligners=tied_aligners,
        confidence=confidence,
        was_5prime_rescued=was_rescued,
        false_junction_removed=False,  # Not tracked; walk back correction handles this
        all_alignments=alignments,
    )
