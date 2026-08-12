"""Refuse to REPORT a physically impossible intron in the consensus output.

THE DEFECT THIS FIXES (Chanfreau planning/684c, verified independently):
`uLTRA` and `deSALT` emit alignments whose N-ops run to hundreds of kilobases —
and the consensus scorer SELECTS them, because "spans more query" outranks "is
physically possible". Measured on 400,001 primary reads of one 684 sample's
``multialigned.bam``: 268 alignments (0.067 %) with an N-op > 10 kb, a maximum
N-op of **261,350 bp**, and 3 alignments whose reference span runs off the end
of the contig entirely. The longest annotated intron in *S. cerevisiae* is
**~1 kb**. The minimap2 arm produced ZERO of these — ``-G`` constrains it, and
nothing equivalent constrained the others or the selection step.

Two consequences, both bad: ``multialigned.bam`` is the input to
``cdna-analyze``, so these become candidate junctions in a discovery pipeline;
and a CIGAR whose reference span exceeds the contig is malformed by definition,
which is how it was found — it crashed ``overhang_informativeness.ambiguity_window``
with an IndexError while indexing the reference by alignment coordinates.

WHY SOFT-CLIP RATHER THAN REJECT THE READ (Kevin's call, 2026-08-12):
the 5' portion of these alignments is usually well supported — in the example
below, 800+ bp of it agrees with the other two aligners. Dropping the read
would discard that evidence; rejecting the whole ARM would change
aligner-selection semantics for every dataset. Truncating at the offending
junction keeps everything upstream of it and returns the rest to soft clip,
which is what an aligner should have emitted in the first place::

    r030_7056  1S67M384N245M1I12M1I1M1I181M 190957N 2M1D53M57S   ends 61.7 kb past chrIV
        ->     1S67M384N245M1I12M1I1M1I181M 112S                 in bounds, 5' evidence kept

STATION C IS A DIFFERENT GATE. ``pool_gate`` (Station C) governs which junctions
are admitted to the candidate POOL that 5' soft-clip rescue draws on. It never
sees, and cannot constrain, which aligner arm WINS consensus — which is why an
impossible junction reached ``multialigned.bam`` with Station C in place.

THRESHOLD. 10 kb for *S. cerevisiae*: an order of magnitude above the longest
annotated intron (~1 kb), so it is a physical-impossibility bound rather than a
biological one and cannot plausibly clip real biology. Override per-organism
with ``RECTIFY_MAX_REPORTABLE_INTRON`` (bp; 0 disables).
"""

import os
from typing import List, Optional, Tuple

__all__ = [
    'DEFAULT_MAX_REPORTABLE_INTRON',
    'max_reportable_intron_from_env',
    'truncate_impossible_introns',
]

# S. cerevisiae: ~10x the longest annotated intron. See module docstring.
DEFAULT_MAX_REPORTABLE_INTRON = 10_000

_ENV_VAR = 'RECTIFY_MAX_REPORTABLE_INTRON'

_QUERY_OPS = {0, 1, 4, 7, 8}   # M I S = X consume query
_REF_OPS = {0, 2, 3, 7, 8}     # M D N = X consume reference
_ALIGNED_OPS = {0, 7, 8}       # M = X — a real aligned block


def max_reportable_intron_from_env(default: int = DEFAULT_MAX_REPORTABLE_INTRON) -> int:
    """Threshold from ``RECTIFY_MAX_REPORTABLE_INTRON`` (0 disables); else default."""
    raw = os.environ.get(_ENV_VAR, '').strip()
    if not raw:
        return default
    try:
        val = int(raw)
    except ValueError:
        return default
    return max(0, val)


def truncate_impossible_introns(
    cigartuples: List[Tuple[int, int]],
    max_intron_bp: int = DEFAULT_MAX_REPORTABLE_INTRON,
    reference_start: int = 0,
    contig_len: Optional[int] = None,
) -> Tuple[Optional[List[Tuple[int, int]]], int]:
    """Soft-clip an alignment at its first impossible intron.

    Returns ``(new_cigartuples, offending_intron_bp)``. ``new_cigartuples`` is
    ``None`` when nothing needed changing, so callers can skip the rewrite.

    An intron is impossible when it exceeds ``max_intron_bp``, OR when the
    alignment's reference span would run past ``contig_len`` (a CIGAR that walks
    off the end of the chromosome is malformed regardless of op sizes).

    Everything from the offending N-op onward is dropped and the query bases it
    covered become a single trailing soft clip, so query length is preserved
    exactly — a changed query length would desync SEQ/QUAL and crash pysam on
    write.

    Returns ``(None, 0)`` when no aligned block would survive the truncation:
    emitting an all-soft-clip record is worse than leaving the read alone, and
    the counter still records it via the returned bp when that is non-zero.
    """
    if not cigartuples or max_intron_bp <= 0:
        return None, 0

    # Locate the first offending N-op: oversize, or one that pushes the
    # reference span past the contig end.
    ref_pos = reference_start
    bad_index = -1
    offending_bp = 0
    for idx, (op, length) in enumerate(cigartuples):
        if op == 3:  # N
            over_size = length > max_intron_bp
            over_end = contig_len is not None and (ref_pos + length) > contig_len
            if over_size or over_end:
                bad_index = idx
                offending_bp = length
                break
        if op in _REF_OPS:
            ref_pos += length
    if bad_index < 0:
        # No offending N-op. The span may still overrun the contig through
        # ordinary ops (a malformed record); report it without truncating here,
        # since there is no junction to cut at.
        if contig_len is not None and ref_pos > contig_len:
            return None, ref_pos - contig_len
        return None, 0

    head = cigartuples[:bad_index]
    if not any(op in _ALIGNED_OPS for op, _ in head):
        return None, offending_bp

    # Every query base at or after the offending op becomes soft clip.
    tail_query = sum(l for op, l in cigartuples[bad_index:] if op in _QUERY_OPS)

    # Trailing D/N in the head consume no query and would leave a dangling
    # reference-only op against a soft clip; drop them.
    while head and head[-1][0] in (2, 3):
        head = head[:-1]
    if not head:
        return None, offending_bp

    # Merge with an existing trailing soft clip rather than emitting `12S34S`.
    if head[-1][0] == 4:
        tail_query += head[-1][1]
        head = head[:-1]

    new_cigar = list(head)
    if tail_query:
        new_cigar.append((4, tail_query))
    return new_cigar, offending_bp
