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
    THREE_PRIME_SIDE_LEFT,
    THREE_PRIME_SIDE_RIGHT,
    walkback_3prime,
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
    """
    if read.is_reverse:
        gene_strand = "+"
        side = THREE_PRIME_SIDE_RIGHT
        stop_base = "A"
    else:
        gene_strand = "-"
        side = THREE_PRIME_SIDE_LEFT
        stop_base = "T"
    orig, corr, applied = walkback_3prime(read, chrom_seq, side, stop_base=stop_base)
    return orig, corr, applied, gene_strand
