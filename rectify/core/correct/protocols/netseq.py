"""NET-seq protocol wrapper for the read-vs-reference walkback.

NET-seq (Churchman-style yeast single-end NET-seq, and human mNET-seq after
paired-end reduction to read 2) is *antisense*, with the nascent RNA 3' end and
its non-templated poly(A) tail on the read 5' side — geometry IDENTICAL to
QuantSeq REV:

    is_reverse = False  ->  gene is '-' strand, 3' end / poly-A at the LEFT
    is_reverse = True   ->  gene is '+' strand, 3' end / poly-A at the RIGHT

For paired-end mNET-seq the nascent 3' end is the 5' end of read 2. Isolate
read 2 UPSTREAM (:mod:`rectify.core.netseq_cpa.pe_reduce`) so every record
reaching this wrapper is single-end. Read 2's 5'-end side then maps by
``is_reverse`` exactly as single-end QuantSeq REV (is_reverse=True -> 5' end at
``reference_end``, RIGHT; is_reverse=False -> 5' end at ``reference_start``,
LEFT), so the same mapping applies unchanged. (Whether read 2's 5'-side base is
basecalled A or T is chemistry-dependent; the QuantSeq REV default per-side
stop_base is used and should be confirmed on a real mNET-seq BAM — see the
``netseq-cpa`` docs / open questions.)

The per-read walkback is therefore identical to QuantSeq REV; this wrapper
delegates to :func:`walkback_quantseq_rev` to guarantee identical, already-
validated behaviour and to give NET-seq a distinct, labelled protocol entry
point (any future protocol-specific divergence belongs here).
"""
from __future__ import annotations

from typing import Tuple

import pysam

from rectify.core.correct.protocols.quantseq_rev import walkback_quantseq_rev


def walkback_netseq(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    *,
    artifact_analyses=None,
) -> Tuple[int, int, str, str]:
    """Walkback for NET-seq / mNET-seq (single-end, or PE after read-2 isolation).

    Returns ``(original_3prime, corrected_3prime, applied, gene_strand)``.

    NET-seq antisense geometry is identical to QuantSeq REV (poly-A on the read
    5' side), so this delegates to :func:`walkback_quantseq_rev`. Paired-end
    read-2 isolation is handled upstream by
    :mod:`rectify.core.netseq_cpa.pe_reduce`, so this wrapper always operates on
    single-end records.
    """
    return walkback_quantseq_rev(read, chrom_seq, artifact_analyses=artifact_analyses)
