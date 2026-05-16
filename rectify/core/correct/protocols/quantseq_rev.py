"""QuantSeq REV protocol wrapper for the read-vs-reference walkback.

QuantSeq-3' REV chemistry produces cDNA antisense to mRNA. After alignment,
BAM strand is the *opposite* of the gene strand:

    is_reverse = False  ->  read aligns to the genome '+' strand  ->  gene is '-'
    is_reverse = True   ->  read aligns to the genome '-' strand  ->  gene is '+'

PolyA-side base orientation in BAM SEQ (symmetric to DRS / cDNA):

    is_reverse = False  ->  polyT (RC of basecalled polyA) at the LEFT side
    is_reverse = True   ->  polyA at the RIGHT side

Read polyT-at-LEFT for ``is_reverse=False`` is the reverse complement of the
basecalled polyA: pysam returns ``query_sequence`` in alignment orientation,
which for ``is_reverse=False`` is the same as the cDNA basecall. The cDNA
basecall starts (5') with the dT-primer's T's, so BAM SEQ at the LEFT
carries T's that need to be walked past with ``stop_base='T'``.
"""
from __future__ import annotations

from typing import Tuple

import pysam

from rectify.core.correct.walkback import (
    APPLIED_NONE,
    APPLIED_WALKBACK,
    THREE_PRIME_SIDE_LEFT,
    THREE_PRIME_SIDE_RIGHT,
    walkback_3prime_guarded,
)


def walkback_quantseq_rev(
    read: pysam.AlignedSegment, chrom_seq: str
) -> Tuple[int, int, str, str]:
    """Walkback for QuantSeq-3' REV.

    Returns ``(original_3prime, corrected_3prime, applied, gene_strand)``.

    Gene strand is the inverse of BAM strand. The polyA-side base in BAM
    SEQ depends on ``is_reverse``: A on the RIGHT for ``is_reverse=True``
    (BAM SEQ is the pysam RC of the basecall, so basecalled T's become A's
    on the right); T on the LEFT for ``is_reverse=False`` (BAM SEQ is the
    basecall in alignment orientation, so the dT-primer T's stay at the
    left).

    Delegates to :func:`walkback_3prime_guarded` so QSrev gets the same
    three artifact guards as DRS (large-deletion pre-scan, intron-bridging
    pre-scan, K-context tail false-stop). The homopolymer early-exit is
    disabled — QSrev V-primer-tip artifacts and internal-priming reads
    routinely terminate at the last genomic A of a tract, so the window
    that early-exit inspects (forward of ``original_3prime``) contains
    no A's even though a real artifact is present.
    """
    if read.is_reverse:
        gene_strand = "+"
        side = THREE_PRIME_SIDE_RIGHT
        stop_base = "A"
        original_3prime = read.reference_end - 1 if read.reference_end is not None else -1
    else:
        gene_strand = "-"
        side = THREE_PRIME_SIDE_LEFT
        stop_base = "T"
        original_3prime = read.reference_start if read.reference_start is not None else -1

    result = walkback_3prime_guarded(
        read,
        chrom_seq,
        side,
        stop_base=stop_base,
        early_exit_homopolymer_check=False,
        right_side_bridging_guard=True,
    )
    if result is None:
        return original_3prime, original_3prime, APPLIED_NONE, gene_strand
    return (
        result["original_pos"],
        result["corrected_pos"],
        APPLIED_WALKBACK,
        gene_strand,
    )
