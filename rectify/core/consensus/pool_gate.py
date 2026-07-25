"""Junction-pool admission gate.

The per-read candidate-junction pool (``select.py``) is what the 5' soft-clip rescue
draws on: any junction in the pool is a licence to convert a soft clip into an aligned
match, because clips cost 1.0/base in the HP-aware edit distance while intron N-ops are
free. Historically the pool admitted **every junction from every aligner with no filter
at all**, so a short spurious overhang slid to the nearest GT-AG entered the pool on
equal footing with a 30-nt-anchored annotated intron.

Measured on ONT DRS (Chanfreau ``planning/456``, minimap2 + gapmm2, 18,295 junction
observations):

===========================  ======  ==================  ==========
class                             n   median min anchor    1-read
===========================  ======  ==================  ==========
annotated (trusted)           1,329               29 nt      18.7%
novel canonical              12,448                3 nt      91.2%
non-canonical                   850               38 nt      91.1%
mega (> 3 kb)                 3,668                6 nt      89.4%
===========================  ======  ==================  ==========

68% of the pool is 3-nt-anchor "novel canonical" — the aligner sliding a tiny overhang
to the nearest GT-AG. The canonical motif is satisfied *because* the anchor is short, so
requiring a canonical motif is nearly useless here; requiring **anchor length** is what
works.

**Anchor, not read count.** 18.7% of annotated junction observations are single-read, so
an ``n_reads >= 2`` gate discards 100% of them while ``min_anchor >= 8`` keeps 99.6%, at
the same pool shrink (85.7% vs 85.6%). Rare junctions are the discovery target of this
work, so abundance must not be an admission criterion — see the workspace rule "read
count is not molecule count". Mismatch and canonicity are usable as soft down-weights,
never as hard gates (``mismatch == 0`` costs 18% of single-read trusted junctions on DRS).

Recommended setting: ``min_anchor_bp`` 8-12 with ``max_intron_len`` 3000 for
*S. cerevisiae*. **Annotated junctions bypass the gate unconditionally**, so the filter
can never damage known biology.

Both parameters default to 0 (disabled) — the gate must be validated against the cryptic
intron catalogue, not merely against annotated introns, before it governs discovery.
"""

from typing import Iterable, List, Optional, Set, Tuple

__all__ = ['junction_min_anchors', 'admit_junction', 'gated_alignment_junctions']

# CIGAR ops that make up a contiguous aligned run flanking an intron.
_ALIGNED_OPS = frozenset('M=X')


def _parse_cigar(cigar_string: str) -> List[Tuple[int, str]]:
    ops: List[Tuple[int, str]] = []
    num = ''
    for ch in cigar_string or '':
        if ch.isdigit():
            num += ch
        else:
            if not num:
                return []          # malformed — treat as unusable
            ops.append((int(num), ch))
            num = ''
    return [] if num else ops


def junction_min_anchors(cigar_string: str) -> List[int]:
    """Min flanking anchor length for each ``N`` op, in CIGAR order.

    The anchor on each side is the **contiguous** aligned run (``M``/``=``/``X``)
    immediately adjacent to the N op; an intervening ``I``/``D``/``S``/``H`` or a second
    ``N`` terminates it. Returns one value per N op, aligned positionally with
    ``AlignmentInfo.junctions`` (both are derived from the N ops in order).

    Returns an empty list if the CIGAR is missing or malformed, which callers must treat
    as "unknown" — i.e. admit, rather than silently drop real junctions.
    """
    ops = _parse_cigar(cigar_string)
    if not ops:
        return []

    def _run(indices: Iterable[int]) -> int:
        n = 0
        for i in indices:
            ln, op = ops[i]
            if op in _ALIGNED_OPS:
                n += ln
            else:
                break
        return n

    return [min(_run(range(i - 1, -1, -1)), _run(range(i + 1, len(ops))))
            for i, (_, op) in enumerate(ops) if op == 'N']


def admit_junction(min_anchor: Optional[int], intron_len: int,
                   min_anchor_bp: int = 0, max_intron_len: int = 0) -> bool:
    """Should this junction enter the per-read candidate pool?

    ``min_anchor`` of ``None`` means the anchor could not be computed (malformed CIGAR);
    such junctions are admitted, because dropping them would silently lose real biology
    on an input-format problem.
    """
    if max_intron_len > 0 and intron_len > max_intron_len:
        return False
    if min_anchor_bp > 0 and min_anchor is not None and min_anchor < min_anchor_bp:
        return False
    return True


def gated_alignment_junctions(
    alignment,
    annotated_junctions: Optional[Set[Tuple[str, int, int, str]]] = None,
    min_anchor_bp: int = 0,
    max_intron_len: int = 0,
) -> List[Tuple[str, int, int, str]]:
    """Junctions from one alignment that are allowed into the candidate pool.

    Annotated junctions bypass both criteria unconditionally.
    """
    juncs = list(getattr(alignment, 'junctions', ()) or ())
    if not juncs:
        return []
    if min_anchor_bp <= 0 and max_intron_len <= 0:
        return [(alignment.chrom, s, e, alignment.strand) for s, e in juncs]

    anchors = junction_min_anchors(getattr(alignment, 'cigar_string', '') or '')
    # Positional correspondence only holds when the counts agree; otherwise treat every
    # anchor as unknown rather than mis-pairing anchors with junctions.
    if len(anchors) != len(juncs):
        anchors = [None] * len(juncs)

    kept = []
    for (start, end), anchor in zip(juncs, anchors):
        key = (alignment.chrom, start, end, alignment.strand)
        if annotated_junctions and key in annotated_junctions:
            kept.append(key)
            continue
        if admit_junction(anchor, end - start, min_anchor_bp, max_intron_len):
            kept.append(key)
    return kept
