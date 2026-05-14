"""Protocol-agnostic 3'-end walkback using read-vs-reference comparison.

Background
----------
Existing Module 2E ("A-tract walk-back") consults *only* the reference: it
flags an ambiguity window wherever the genome has a runs of A/T near the 3'
alignment end, but it doesn't validate that the *read* actually contains a
poly(A) extension. The result is correct ambiguity detection but no
correction — so for many protocols (QuantSeq REV in particular) ~92% of reads
emit ``correction_applied="none"`` while their alignments are silently
drifting 1–N bp into a genomic A-stretch.

The cDNA pipeline (``rectify/core/cdna_correct_command.py::walk_back_anchor_
and_tail``) already implements the correct algorithm — compare the read base
to the reference base at each aligned position, walk past stop-base (A/T)
matches, and anchor at the first non-stop match. The logic is protocol-
agnostic once a few wiring decisions are made (which side of the alignment is
the 3' end; what the "stop base" is on the read; how the gene strand maps to
``read.is_reverse``).

This module is the consolidated home for that algorithm. Each protocol-
specific entry point in ``rectify/core/correct/protocols/`` calls
:func:`walkback_3prime` with the appropriate wiring.

Author: Kevin R. Roy
Issue: ``.github/ISSUE_walkback_integration.md`` (NEW-075)
"""
from __future__ import annotations

from typing import Optional, Tuple

import pysam


THREE_PRIME_SIDE_RIGHT = "right"
THREE_PRIME_SIDE_LEFT = "left"

# Outcome string written to ``correction_applied`` column of corrected_3ends.tsv
APPLIED_WALKBACK = "polya_walkback_readgenome"
APPLIED_NONE = "none"


def walkback_3prime(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    three_prime_side: str,
    stop_base: str = "A",
) -> Tuple[int, int, str]:
    """Walk back from the 3' alignment end through a stop-base tract.

    The walk-back scans aligned (query, reference) pairs starting at the side
    of the alignment that holds the 3' end (``three_prime_side``) and moves
    inward toward the gene body. At each position the read base is compared
    to the reference base. The first position where they agree *and* the
    base is not the stop base is the canonical cleavage anchor.

    One terminal gate is applied upfront to skip the walk-back in the
    unambiguously correct case:

    1. Terminal is a non-stop base AND matches the reference → already
       anchored; no walk-back.

    All other cases (terminal is a stop-base, or terminal mismatches the
    reference) proceed to the walk-back scan. Critically, a terminal stop-base
    (A) that happens to match a genomic A is NOT skipped — non-genomically
    encoded poly-A tails frequently align over genomic A-tracts (the exact
    scenario RECTIFY was designed to correct), so the scan must walk inward to
    find the true cleavage anchor regardless of genomic A-content.

    Parameters
    ----------
    read
        ``pysam.AlignedSegment`` to correct. Must be mapped (no upstream
        check is performed).
    chrom_seq
        Reference sequence for ``read.reference_name``. 0-based, on the
        ``+`` genomic strand. The caller is responsible for selecting the
        right chromosome from the cache.
    three_prime_side
        Either :data:`THREE_PRIME_SIDE_RIGHT` (poly-A at the rightmost ref
        position, i.e. ``reference_end - 1``) or
        :data:`THREE_PRIME_SIDE_LEFT` (poly-A at ``reference_start``). The
        wrapping protocol module decides this from ``read.is_reverse`` and
        protocol chemistry.
    stop_base
        Nucleotide to walk past. Default ``"A"`` (correct for poly-A tails
        on every protocol after pysam reverse-complementation).

    Returns
    -------
    original_3prime
        Reference position of the original 3' end before any correction.
    corrected_3prime
        Reference position of the canonical cleavage anchor. Equal to
        ``original_3prime`` when no correction was applied.
    applied
        One of :data:`APPLIED_WALKBACK` or :data:`APPLIED_NONE`.

    Notes
    -----
    *Strand handling.* This function does not know or care about gene
    strand. It returns positions in reference (``+`` genomic) coordinates
    only. Mapping ``read.is_reverse`` ↔ gene strand belongs to the
    protocol wrapper, which then passes ``three_prime_side`` accordingly.

    *Intron handling.* ``get_aligned_pairs(matches_only=False)`` is used
    with a manual ``qp is not None and rp is not None`` filter — this
    drops insertions, deletions, and ``N``-cigar (splice) gaps so the scan
    naturally crosses introns without special-case code.

    *No walk-back cap.* The aligned-pair list is bounded by read length
    (~150 bp for QuantSeq, ~longer for ONT), which serves as the natural
    cap.

    *V-primer tip artifact.* QuantSeq libraries sometimes terminate with a
    spurious non-A base (e.g. ``...AAAAAAAG``). The terminal-mismatch gate
    falls into case 3 and the scan continues through the A-run until a
    valid genomic anchor is found.
    """
    if three_prime_side not in (THREE_PRIME_SIDE_RIGHT, THREE_PRIME_SIDE_LEFT):
        raise ValueError(
            f"three_prime_side must be 'right' or 'left', got {three_prime_side!r}"
        )
    stop_base = stop_base.upper()
    if len(stop_base) != 1 or stop_base not in "ACGT":
        raise ValueError(f"stop_base must be one of A/C/G/T, got {stop_base!r}")

    if three_prime_side == THREE_PRIME_SIDE_RIGHT:
        original_3prime = read.reference_end - 1
    else:
        original_3prime = read.reference_start

    pairs = [
        (qp, rp)
        for qp, rp in read.get_aligned_pairs(matches_only=False)
        if qp is not None and rp is not None
    ]
    if not pairs:
        return original_3prime, original_3prime, APPLIED_NONE

    if three_prime_side == THREE_PRIME_SIDE_RIGHT:
        scan_pairs = list(reversed(pairs))
    else:
        scan_pairs = pairs

    qs = read.query_sequence
    if qs is None:
        return original_3prime, original_3prime, APPLIED_NONE

    # --- Terminal-position gate: skip walkback in the obvious-correct cases. ---
    first_qp, first_rp = scan_pairs[0]
    first_rb = qs[first_qp].upper()
    first_gb = chrom_seq[first_rp].upper()

    # Gate: non-stop read base matching genome → already anchored.
    if first_rb != stop_base and first_rb == first_gb:
        return original_3prime, original_3prime, APPLIED_NONE

    # Walk back until a non-stop read-vs-ref match is found.
    # This covers terminal stop-base (A) regardless of genomic A-content:
    # poly-A tails routinely align over genomic A-tracts and must be walked
    # back to the true cleavage anchor.
    corrected = original_3prime
    changed = False
    for qp, rp in scan_pairs:
        read_base = qs[qp].upper()
        ref_base = chrom_seq[rp].upper()
        if read_base == ref_base and read_base != stop_base:
            corrected = rp
            changed = corrected != original_3prime
            break

    return original_3prime, corrected, (APPLIED_WALKBACK if changed else APPLIED_NONE)


def is_minus_strand_dRNA(read: pysam.AlignedSegment) -> bool:
    """Return ``True`` if a direct-RNA read is on the genome ``-`` strand.

    For ONT direct-RNA (dRNA-seq) the BAM strand passes through unchanged:
    ``is_reverse=False`` is a ``+`` strand RNA, ``is_reverse=True`` is a
    ``-`` strand RNA.
    """
    return read.is_reverse


def walkback_drs(
    read: pysam.AlignedSegment, chrom_seq: str
) -> Tuple[int, int, str, str]:
    """Walk-back wrapper for ONT direct-RNA sequencing (dRNA-seq).

    Strand convention: BAM strand == gene strand. So a ``+`` strand RNA
    has its 3' poly(A) tail aligning at the rightmost reference position;
    a ``-`` strand RNA has its 3' poly(A) tail (in cDNA coordinates → A's
    after pysam RC) at the leftmost reference position.

    Returns ``(original_3prime, corrected_3prime, applied, gene_strand)``.
    """
    if is_minus_strand_dRNA(read):
        gene_strand = "-"
        side = THREE_PRIME_SIDE_LEFT
    else:
        gene_strand = "+"
        side = THREE_PRIME_SIDE_RIGHT
    orig, corr, applied = walkback_3prime(read, chrom_seq, side, stop_base="A")
    return orig, corr, applied, gene_strand
