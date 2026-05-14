"""QuantSeq REV protocol wrapper for the read-vs-reference walkback.

QuantSeq-3' REV chemistry produces cDNA antisense to mRNA. After alignment,
BAM strand is the *opposite* of the gene strand:

    is_reverse = False  ->  read aligns to the genome '+' strand  ->  gene is '-'
    is_reverse = True   ->  read aligns to the genome '-' strand  ->  gene is '+'

Poly-A tail (after pysam reverse-complementation) appears as A's at:
    is_reverse = False  ->  the LEFT side of the alignment (lowest ref coord)
    is_reverse = True   ->  the RIGHT side (highest ref coord)
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

    Gene strand is the inverse of BAM strand. The 3' poly-A side of the
    aligned read maps to whichever genomic end the antisense cDNA puts it.
    """
    if read.is_reverse:
        gene_strand = "+"
        side = THREE_PRIME_SIDE_RIGHT
    else:
        gene_strand = "-"
        side = THREE_PRIME_SIDE_LEFT
    orig, corr, applied = walkback_3prime(read, chrom_seq, side, stop_base="A")
    return orig, corr, applied, gene_strand
