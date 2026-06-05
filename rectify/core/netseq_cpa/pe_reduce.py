"""Paired-end reduction for human mNET-seq.

In paired-end mNET-seq (Nojima 2015 and successors) the nascent RNA 3' end is
the **5' end of read 2**. To reuse the single-end antisense walkback unchanged,
isolate read 2 UPSTREAM so every record reaching
:func:`...protocols.netseq.walkback_netseq` is single-end. Read 2's 5'-end side
then maps by ``is_reverse`` exactly as single-end NET-seq / QuantSeq REV
(is_reverse=True -> 5' end at ``reference_end``, RIGHT; is_reverse=False -> 5'
end at ``reference_start``, LEFT).

Two entry points depending on the orchestrator's input:

  * :func:`reduce_bam_to_read2` — a paired-end BAM -> single-end BAM of read-2
    primary records, with the paired/mate flags cleared so downstream tools and
    the walkback treat each as a standalone single-end read.
  * :func:`select_read2_fastq` — for a FASTQ pair, read 2 IS the input; this
    just validates and returns the R2 path (no transform needed pre-alignment).

NOTE (open question, see ``netseq-cpa`` docs): whether read 2's 5'-side poly-A
base is basecalled A or T is library-chemistry-dependent. A real mNET-seq BAM
showed the R2 5' end ~7x poly-T-enriched over poly-A, consistent with antisense;
confirm the per-side ``stop_base`` on aligned data before trusting human poly-A
lengths. The QuantSeq REV default is used until then.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict

import pysam


def reduce_bam_to_read2(in_bam: str | Path, out_bam: str | Path) -> Dict[str, int]:
    """Write a single-end BAM containing only read-2 primary alignments.

    Mate/paired flags (``is_paired``, ``is_proper_pair``, ``mate_is_*``,
    ``is_read1``, ``is_read2``) are cleared so each surviving record is a valid
    standalone single-end alignment. Secondary/supplementary records are
    dropped. Returns ``{seen, written}``.
    """
    bam = pysam.AlignmentFile(str(in_bam), "rb")
    out = pysam.AlignmentFile(str(out_bam), "wb", template=bam)
    n_seen = n_written = 0
    try:
        for read in bam.fetch(until_eof=True):
            n_seen += 1
            if read.is_secondary or read.is_supplementary:
                continue
            if not (read.is_paired and read.is_read2):
                continue
            # Collapse to single-end: clear all pairing-related flags.
            read.is_paired = False
            read.is_proper_pair = False
            read.mate_is_unmapped = False
            read.mate_is_reverse = False
            read.is_read1 = False
            read.is_read2 = False
            read.next_reference_id = -1
            read.next_reference_start = -1
            read.template_length = 0
            out.write(read)
            n_written += 1
    finally:
        out.close()
        bam.close()
    return {"seen": n_seen, "written": n_written}


def select_read2_fastq(r1_fastq: str | Path, r2_fastq: str | Path) -> Path:
    """For a paired FASTQ, read 2 is the nascent-3'-end read; return its path.

    Raises ``FileNotFoundError`` if R2 is missing. ``r1_fastq`` is accepted for
    a symmetric, self-documenting call site but is not read.
    """
    r2 = Path(r2_fastq)
    if not r2.exists():
        raise FileNotFoundError(f"mNET-seq read-2 FASTQ not found: {r2}")
    return r2
