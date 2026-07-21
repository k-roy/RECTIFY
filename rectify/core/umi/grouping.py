"""Fragment-span grouping keys for paired-end UMI deduplication.

A UMI alone does not identify a molecule -- the same 12-mer recurs across a
library by chance (4^12 = 16.8M, and real UMIs are non-uniform, so effective
diversity is lower still). A molecule is identified by (fragment span x UMI).

WHY EXACT POSITIONS, NOT A TOLERANCE WINDOW
-------------------------------------------
PCR duplicates of one molecule share BOTH fragment ends *by construction* -- they
are copies. The real end-wobble comes from soft-clipping, which
:func:`five_prime_unclipped` corrects DETERMINISTICALLY (recovering the true 5' end
of the molecule) rather than fuzzily. This matches every mainstream tool:
umi_tools groups on exact start "corrected for soft-clipping"; fgbio requires "the
same genomic position"; gencore requires "identical mapping position". None
implements a positional tolerance.

A tolerance window is not a free safety margin -- it is actively costly:
paired-end grouping keys on BOTH ends, which shatters a locus's molecules into
many small (span x UMI) bins and is precisely why PE grouping is more specific
than single-end. A window merges those bins back together, raising the effective
number of molecules per bin -- the same term that drives the false-merge rate. So
tolerance and edit distance interact MULTIPLICATIVELY; raising both is the worst
case. A window also makes grouping non-transitive (a~b, b~c, but not a~c), which
only a graph discipline can resolve.

``tolerance`` is therefore exposed (default 0) so the hypothesis can be TESTED on
real data rather than assumed either way.

STRANDEDNESS
------------
``r1_sense`` is protocol-dependent and there is no safe default across
chemistries, so it is explicit:

  * **Lexogen CORALL (V1/V2): r1_sense=True.** Its User Guide flags this as
    unusual -- "In contrast to most other library preparation protocols, CORALL
    libraries generate reads in forward orientation"; "Read 1 reflects the RNA
    transcript sequence not the cDNA sequence" (App. E p.36). Strand specificity
    >99%.
  * dUTP-based protocols (TruSeq stranded) and QuantSeq REV are ANTISENSE:
    ``r1_sense=False``. RECTIFY's existing protocol wrappers
    (``core/correct/protocols/quantseq_rev.py``) encode that opposite convention --
    do not inherit it here by reflex.
"""
from __future__ import annotations

import enum
from typing import Any, List, NamedTuple, Optional, Tuple

# pysam CIGAR operation codes.
_CIGAR_SOFT_CLIP = 4
_CIGAR_HARD_CLIP = 5
_CIGAR_REF_CONSUMING = frozenset({0, 2, 3, 7, 8})  # M D N = X


class PairStatus(enum.Enum):
    """How completely a template's span is known.

    This is part of the grouping key ON PURPOSE. A degraded (single-end) key is
    LESS specific than a full-span key, so if the two were allowed to share a
    bucket, a degraded read could absorb molecules that the full span would have
    kept separate. Keeping the status in the key makes that impossible.

    This is not a rare edge case: reads spanning novel junctions are exactly where
    a mate fails to align or aligns discordantly, so the degraded path carries
    real signal and must be handled deliberately.
    """
    PROPER_PAIR = "proper_pair"      # both mates mapped -- full span known
    SINGLETON_R1 = "singleton_r1"    # only R1 usable
    SINGLETON_R2 = "singleton_r2"    # only R2 usable
    UNPAIRED = "unpaired"            # single-end library


class FragmentKey(NamedTuple):
    """The molecule-identifying key, excluding the UMI.

    Reads sharing a FragmentKey are candidates to be the same molecule; the UMI
    (clustered directionally) then separates true molecules within the bucket.
    """
    reference_name: str
    strand: str                 # '+' / '-' -- the RNA/gene strand (see r1_sense)
    five_prime_r1: Optional[int]
    five_prime_r2: Optional[int]
    pair_status: PairStatus
    spliced: bool               # spliced != unspliced at one position (umi_tools --spliced-is-unique)


def leading_soft_clip(read: Any) -> int:
    """Soft-clipped bases before the first aligned base (hard clips do not count)."""
    ct = read.cigartuples
    if not ct:
        return 0
    n = 0
    for op, length in ct:
        if op == _CIGAR_SOFT_CLIP:
            n += length
        elif op == _CIGAR_HARD_CLIP:
            continue
        else:
            break
    return n


def trailing_soft_clip(read: Any) -> int:
    """Soft-clipped bases after the last aligned base."""
    ct = read.cigartuples
    if not ct:
        return 0
    n = 0
    for op, length in reversed(ct):
        if op == _CIGAR_SOFT_CLIP:
            n += length
        elif op == _CIGAR_HARD_CLIP:
            continue
        else:
            break
    return n


def is_spliced(read: Any) -> bool:
    """True if the alignment contains a CIGAR N operation (a splice junction)."""
    ct = read.cigartuples
    if not ct:
        return False
    return any(op == 3 for op, _ in ct)


def five_prime_unclipped(read: Any) -> Optional[int]:
    """The read's 5' end in reference coordinates, corrected for soft-clipping.

    This is the deterministic wobble correction that makes a tolerance window
    unnecessary: it recovers where the MOLECULE actually started, even when the
    aligner soft-clipped the first bases (adapter remnant, low-quality 5' cycle,
    or a base that would not align across a junction).

    Forward read: reference_start - leading soft clip.
    Reverse read: reference_end   + trailing soft clip (its 5' end is at the
    higher coordinate).
    """
    if read.is_unmapped or read.reference_start is None:
        return None
    if read.is_reverse:
        end = read.reference_end
        if end is None:
            return None
        return end + trailing_soft_clip(read)
    return read.reference_start - leading_soft_clip(read)


def fragment_strand(read1: Optional[Any], read2: Optional[Any], r1_sense: bool) -> str:
    """The RNA/gene strand of the fragment.

    Derived from R1 where available (R1 defines the orientation), else inferred
    from R2 by inverting. See the module docstring on ``r1_sense``.
    """
    anchor, invert = (read1, False) if read1 is not None else (read2, True)
    if anchor is None:
        raise ValueError("fragment_strand requires at least one mate")
    read_is_forward = not anchor.is_reverse
    sense = read_is_forward if not invert else anchor.is_reverse
    if not r1_sense:
        sense = not sense
    return "+" if sense else "-"


def _bin_position(pos: Optional[int], tolerance: int) -> Optional[int]:
    """Quantise a position when a tolerance is requested.

    Binning (rather than fuzzy pairwise matching) keeps grouping TRANSITIVE and
    O(1): equal keys group, unequal keys do not. It is a coarser instrument than
    pairwise tolerance -- two positions 1 bp apart can still land in different
    bins -- but it cannot produce the a~b, b~c, not-a~c pathology that a pairwise
    window does. Given tolerance is off by default and exists to be measured, a
    predictable approximation beats an unpredictable exact one.
    """
    if pos is None or tolerance <= 0:
        return pos
    width = 2 * tolerance + 1
    return pos // width


def fragment_key(
    read1: Optional[Any],
    read2: Optional[Any] = None,
    r1_sense: bool = True,
    tolerance: int = 0,
) -> FragmentKey:
    """Build the molecule-identifying key for a template.

    Pass both mates when available. When a mate is missing or unmapped the key
    degrades to single-end AND records that fact in ``pair_status``, so degraded
    templates can never share a bucket with full-span ones.
    """
    r1_ok = read1 is not None and not read1.is_unmapped
    r2_ok = read2 is not None and not read2.is_unmapped
    if not r1_ok and not r2_ok:
        raise ValueError("fragment_key requires at least one mapped mate")

    anchor = read1 if r1_ok else read2
    paired = bool(getattr(anchor, "is_paired", False))

    if r1_ok and r2_ok:
        status = PairStatus.PROPER_PAIR
    elif not paired:
        status = PairStatus.UNPAIRED
    elif r1_ok:
        status = PairStatus.SINGLETON_R1
    else:
        status = PairStatus.SINGLETON_R2

    fp1 = five_prime_unclipped(read1) if r1_ok else None
    fp2 = five_prime_unclipped(read2) if r2_ok else None

    spliced = (r1_ok and is_spliced(read1)) or (r2_ok and is_spliced(read2))

    return FragmentKey(
        reference_name=anchor.reference_name,
        strand=fragment_strand(read1 if r1_ok else None, read2 if r2_ok else None, r1_sense),
        five_prime_r1=_bin_position(fp1, tolerance),
        five_prime_r2=_bin_position(fp2, tolerance),
        pair_status=status,
        spliced=spliced,
    )
