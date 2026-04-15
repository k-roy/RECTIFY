#!/usr/bin/env python3
"""
BAM output writers for RECTIFY.

Handles CIGAR surgery and BAM writing for all correction categories:

- ``write_dual_bam``        — single-pass writer for both hardclip + softclip BAMs (preferred)
- ``write_corrected_bam``  — hard-clips 3' end; extends 5' for Cat3
- ``write_softclipped_bam`` — soft-clips 3' end; extends 5' for Cat3
- ``write_polya_trimmed_bam`` — strips 3' poly(A) soft-clips
- ``write_netseq_assigned_bedgraph`` — Cat6 fractional position bedgraph

CIGAR helpers (in-place pysam surgery):
- ``clip_read_to_corrected_3prime``
- ``softclip_read_to_corrected_3prime``
- ``extend_read_5prime_for_junction_rescue``
- ``softclip_intronic_tail_5prime``   (available but not called from write functions)
- ``fix_homopolymer_mismatches``  — (deprecated) converts X at homopolymer positions to I/D
- ``realign_exon_blocks``          — global NW re-alignment of exon blocks with homopolymer scoring

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, Optional, Tuple
import logging
import pysam
from pathlib import Path

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# CIGAR op classification constants
# ---------------------------------------------------------------------------
# Ops that advance the reference coordinate
_REF_CONSUMING: frozenset = frozenset([0, 2, 3, 7, 8])   # M, D, N, =, X
# Ops that consume stored query bases (i.e. present in query_sequence)
_QUERY_CONSUMING: frozenset = frozenset([0, 1, 4, 7, 8])  # M, I, S, =, X
# Note: H (5) is excluded from both — hard-clipped bases are NOT in query_sequence.


def clip_read_to_corrected_3prime(
    read: pysam.AlignedSegment,
    corrected_3prime: int,
    strand: str,
) -> bool:
    """
    Hard-clip *read* in-place so its 3' alignment boundary matches *corrected_3prime*.

    For + strand reads the 3' boundary is ``reference_end - 1`` (rightmost ref base).
    For − strand reads it is ``reference_start`` (leftmost ref base, which is the 3'
    end for reverse-strand RNA).

    Bases removed from the CIGAR are converted to a single trailing / leading H
    (hard-clip) op and deleted from the stored query sequence and quality array.
    Any existing soft-clip at the clipped end is subsumed into the new hard clip.

    If *corrected_3prime* equals the current boundary the read is unchanged and
    False is returned.  Unmapped reads, reads without CIGAR, and reads where the
    corrected position would leave an empty alignment are also left unchanged.

    Args:
        read:             pysam AlignedSegment — modified in-place.
        corrected_3prime: Target 3' reference coordinate (0-based, inclusive).
        strand:           '+' or '-'.

    Returns:
        True if clipping was applied; False otherwise.
    """
    cigar = list(read.cigartuples or [])
    if not cigar or read.is_unmapped:
        return False

    seq   = read.query_sequence
    quals = read.query_qualities  # numpy uint8 array or None

    if strand == '+':
        current_end = read.reference_end - 1  # 0-based inclusive right boundary
        if corrected_3prime >= current_end:
            return False  # no clipping needed

        n_ref_clip = current_end - corrected_3prime  # reference bases to remove

        n_query_remove = 0

        # 1. Strip trailing soft-clips: their bases ARE in query_sequence.
        while cigar and cigar[-1][0] == 4:  # S
            n_query_remove += cigar.pop()[1]

        # 2. Strip trailing hard-clips: their bases are NOT in query_sequence.
        while cigar and cigar[-1][0] == 5:  # H
            cigar.pop()

        # 3. Walk CIGAR ops from the right, removing n_ref_clip reference bases.
        n_ref_removed = 0
        while cigar and n_ref_removed < n_ref_clip:
            op, length = cigar[-1]
            ref_in_op   = length if op in _REF_CONSUMING  else 0
            query_in_op = length if op in _QUERY_CONSUMING else 0
            need        = n_ref_clip - n_ref_removed

            if ref_in_op <= need:
                cigar.pop()
                n_ref_removed  += ref_in_op
                n_query_remove += query_in_op
            else:
                # Partial removal: trim `need` reference positions from this op.
                cigar[-1] = (op, length - need)
                n_ref_removed  += need
                n_query_remove += need if op in _QUERY_CONSUMING else 0
                break

        if not cigar:
            return False  # degenerate: entire alignment clipped

        # 4. Append hard-clip for removed query bases (if any).
        if n_query_remove > 0:
            cigar.append((5, n_query_remove))

        # 5. Apply to read.  Set sequence before qualities (setting seq resets quals).
        if n_query_remove > 0 and seq is not None:
            new_seq  = seq[:-n_query_remove]
            new_qual = quals[:-n_query_remove] if quals is not None else None
            read.query_sequence  = new_seq
            read.query_qualities = new_qual
        read.cigartuples = cigar
        return True

    else:  # minus strand — clip from the left (remove 3' bases at low coordinates)
        current_start = read.reference_start  # 0-based inclusive left boundary
        if corrected_3prime <= current_start:
            return False

        n_ref_clip = corrected_3prime - current_start

        n_query_remove = 0

        # 1. Strip leading soft-clips.
        while cigar and cigar[0][0] == 4:  # S
            n_query_remove += cigar.pop(0)[1]

        # 2. Strip leading hard-clips.
        while cigar and cigar[0][0] == 5:  # H
            cigar.pop(0)

        # 3. Walk CIGAR ops from the left.
        n_ref_removed = 0
        while cigar and n_ref_removed < n_ref_clip:
            op, length = cigar[0]
            ref_in_op   = length if op in _REF_CONSUMING  else 0
            query_in_op = length if op in _QUERY_CONSUMING else 0
            need        = n_ref_clip - n_ref_removed

            if ref_in_op <= need:
                cigar.pop(0)
                n_ref_removed  += ref_in_op
                n_query_remove += query_in_op
            else:
                cigar[0] = (op, length - need)
                n_ref_removed  += need
                n_query_remove += need if op in _QUERY_CONSUMING else 0
                break

        if not cigar:
            return False

        # 4. Prepend hard-clip for removed query bases.
        if n_query_remove > 0:
            cigar.insert(0, (5, n_query_remove))

        # 5. Apply to read; also shift reference_start for minus strand clips.
        if n_query_remove > 0 and seq is not None:
            new_seq  = seq[n_query_remove:]
            new_qual = quals[n_query_remove:] if quals is not None else None
            read.query_sequence  = new_seq
            read.query_qualities = new_qual
        read.cigartuples     = cigar
        read.reference_start = current_start + n_ref_removed
        return True


def softclip_read_to_corrected_3prime(
    read: pysam.AlignedSegment,
    corrected_3prime: int,
    strand: str,
) -> bool:
    """
    Soft-clip *read* in-place so its 3' alignment boundary matches *corrected_3prime*.

    Identical logic to :func:`clip_read_to_corrected_3prime` except that removed
    query bases are replaced with a ``S`` (soft-clip) op rather than an ``H``
    (hard-clip) op.  The bases remain in ``query_sequence`` and are visible in IGV
    when "Show soft-clipped bases" is enabled, making the poly(A) tail location
    apparent without affecting pileup or coverage tracks.

    ``reference_start`` is updated for minus-strand reads exactly as in the hard-clip
    path (soft-clips do not consume reference, so the first aligned base moves).

    Returns True if soft-clipping was applied; False otherwise.
    """
    cigar = list(read.cigartuples or [])
    if not cigar or read.is_unmapped:
        return False

    seq   = read.query_sequence
    quals = read.query_qualities

    if strand == '+':
        current_end = read.reference_end - 1
        if corrected_3prime >= current_end:
            return False

        n_ref_clip    = current_end - corrected_3prime
        n_query_remove = 0

        # 1. Absorb trailing soft-clips into the new soft-clip.
        while cigar and cigar[-1][0] == 4:   # S
            n_query_remove += cigar.pop()[1]
        # 2. Strip trailing hard-clips (not in sequence; discard silently).
        while cigar and cigar[-1][0] == 5:   # H
            cigar.pop()

        # 3. Walk CIGAR from right, accounting for n_ref_clip reference bases.
        n_ref_removed = 0
        while cigar and n_ref_removed < n_ref_clip:
            op, length = cigar[-1]
            ref_in_op   = length if op in _REF_CONSUMING  else 0
            query_in_op = length if op in _QUERY_CONSUMING else 0
            need        = n_ref_clip - n_ref_removed
            if ref_in_op <= need:
                cigar.pop()
                n_ref_removed  += ref_in_op
                n_query_remove += query_in_op
            else:
                cigar[-1] = (op, length - need)
                n_ref_removed  += need
                n_query_remove += need if op in _QUERY_CONSUMING else 0
                break

        if not cigar:
            return False

        # 4. Append soft-clip (sequence stays in read).
        if n_query_remove > 0:
            cigar.append((4, n_query_remove))  # S, not H
        read.cigartuples = cigar
        return True

    else:  # minus strand — clip from the left
        current_start = read.reference_start
        if corrected_3prime <= current_start:
            return False

        n_ref_clip    = corrected_3prime - current_start
        n_query_remove = 0

        # 1. Absorb leading soft-clips.
        while cigar and cigar[0][0] == 4:   # S
            n_query_remove += cigar.pop(0)[1]
        # 2. Strip leading hard-clips.
        while cigar and cigar[0][0] == 5:   # H
            cigar.pop(0)

        # 3. Walk CIGAR from left.
        n_ref_removed = 0
        while cigar and n_ref_removed < n_ref_clip:
            op, length = cigar[0]
            ref_in_op   = length if op in _REF_CONSUMING  else 0
            query_in_op = length if op in _QUERY_CONSUMING else 0
            need        = n_ref_clip - n_ref_removed
            if ref_in_op <= need:
                cigar.pop(0)
                n_ref_removed  += ref_in_op
                n_query_remove += query_in_op
            else:
                cigar[0] = (op, length - need)
                n_ref_removed  += need
                n_query_remove += need if op in _QUERY_CONSUMING else 0
                break

        if not cigar:
            return False

        # 4. Prepend soft-clip (sequence stays in read).
        if n_query_remove > 0:
            cigar.insert(0, (4, n_query_remove))  # S, not H
        read.cigartuples     = cigar
        read.reference_start = current_start + n_ref_removed
        return True


def fix_homopolymer_mismatches(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    min_homopolymer: int = 3,
) -> bool:
    """
    Convert mismatches at homopolymer positions to the more parsimonious indel form.

    Nanopore DRS systematically undercalls homopolymer lengths.  Wherever the
    aligner emits a mismatch (X) at a reference position that is within a
    homopolymer run ≥ min_homopolymer, the mismatch is better explained by a
    length error in the run than by a true sequence difference.

    Two patterns are handled for each X:

    **Deletion** (read is missing one copy of the homopolymer base):
      ref:   ...T T [T] A...   (T-run; the X is at the last T)
      read:  ...T T  A ...     (T-run undercalled by 1)
      Check: query[qp] == ref[rp+1]  →  replace  X  with  D

    **Insertion** (read has one extra copy of the homopolymer base):
      ref:   ... [T] G A...    (the X is at the first base AFTER the T-run)
      read:  ...T T  G A...    (extra T in read)
      Check: query[qp+1] == ref[rp]  →  replace  X  with  I

    Only length-1 X ops are converted (longer X runs are left unchanged).
    All aligned positions are checked, not just junction-adjacent ones.
    The CIGAR is modified in-place; query_sequence and reference_start are
    unchanged.

    Returns True if at least one X was converted, False otherwise.
    """
    if not read.cigartuples or not read.query_sequence:
        return False

    chrom = read.reference_name
    if not chrom or chrom not in genome:
        return False

    ref_seq = genome[chrom]
    query_seq = read.query_sequence

    # Build per-op (ref_start, query_start) lookup in a single pass.
    cigar = list(read.cigartuples)
    op_ref: list = []
    op_qry: list = []
    rp = read.reference_start
    qp = 0
    for op, length in cigar:
        op_ref.append(rp)
        op_qry.append(qp)
        if op in _REF_CONSUMING:
            rp += length
        if op in _QUERY_CONSUMING:
            qp += length

    def _homopolymer_run(pos: int) -> int:
        """Length of the homopolymer run that covers reference position pos (0 if none)."""
        if pos < 0 or pos >= len(ref_seq):
            return 0
        base = ref_seq[pos].upper()
        if base == 'N':
            return 0
        left = pos
        while left > 0 and ref_seq[left - 1].upper() == base:
            left -= 1
        right = pos + 1
        while right < len(ref_seq) and ref_seq[right].upper() == base:
            right += 1
        return right - left

    modified = False
    # Iterate in reverse so deletions/insertions don't invalidate forward indices.
    for i in range(len(cigar) - 1, -1, -1):
        op, length = cigar[i]
        if op != 8 or length != 1:  # only single-base X ops
            continue

        rp_x = op_ref[i]
        qp_x = op_qry[i]

        if rp_x >= len(ref_seq) or qp_x >= len(query_seq):
            continue

        ref_base = ref_seq[rp_x].upper()
        read_base = query_seq[qp_x].upper()

        # --- Deletion pattern ---
        # Read is missing one copy of a homopolymer base.  The reference has
        # the homopolymer base at rp_x, but the read already "ran out" and has
        # the next non-homopolymer base at qp_x.
        # Validity: ref[rp_x] is in a homopolymer AND query[qp_x] == ref[rp_x+1].
        if _homopolymer_run(rp_x) >= min_homopolymer:
            if rp_x + 1 < len(ref_seq) and read_base == ref_seq[rp_x + 1].upper():
                # Replace 1X with 1D (consumes ref, not query)
                cigar[i] = (2, 1)
                modified = True
                continue

        # --- Insertion pattern ---
        # Read has one extra copy of a homopolymer base before the current
        # reference position.  The reference at rp_x is some other base but the
        # preceding reference positions form a homopolymer run, and the read has
        # an extra copy of that run base here.
        # Validity: ref[rp_x - 1] is in a homopolymer run AND read_base == ref[rp_x-1]
        #           AND query[qp_x + 1] == ref[rp_x].
        if rp_x > 0 and _homopolymer_run(rp_x - 1) >= min_homopolymer:
            homo_base = ref_seq[rp_x - 1].upper()
            if (read_base == homo_base
                    and qp_x + 1 < len(query_seq)
                    and query_seq[qp_x + 1].upper() == ref_base):
                # Replace 1X with 1I (consumes query, not ref)
                cigar[i] = (1, 1)
                modified = True
                continue

    if modified:
        # Merge adjacent identical ops (e.g. two consecutive D → one D).
        merged: list = []
        for op, length in cigar:
            if merged and merged[-1][0] == op:
                merged[-1] = (op, merged[-1][1] + length)
            else:
                merged.append((op, length))
        read.cigartuples = merged

    return modified


def realign_exon_blocks(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    min_run: int = 3,
    max_block_bp: int = 120,
) -> bool:
    """
    Re-align exon blocks (between N intron-skip ops) using global affine-gap NW
    with homopolymer-aware mismatch scoring.

    Nanopore DRS systematically undercalls homopolymer lengths, causing mismatches
    (X ops) near splice junctions when the correct alignment is an indel. By re-aligning
    each exon block globally with a reduced penalty at homopolymer ref positions, the
    aligner can choose indels over mismatches where the evidence supports it.

    The global alignment guarantees that the query and reference span of each exon block
    are preserved, so ``query_sequence`` and ``reference_start`` are unchanged.

    Only reads that have at least one X op at a homopolymer ref position are
    processed; all others are returned unchanged (fast pre-check).  Exon blocks
    longer than *max_block_bp* are skipped because O(Q×R) Python DP is too slow
    for large blocks; a C extension can lift this limit in a future version.

    Args:
        read:         pysam AlignedSegment — CIGAR modified in-place if changed.
        genome:       Pre-loaded genome dict {chrom: sequence}.
        min_run:      Minimum homopolymer run length for the reduced mismatch penalty.
        max_block_bp: Skip blocks longer than this (ref bases). Default 120 bp.

    Returns:
        True if at least one CIGAR op was changed; False otherwise.
    """
    from .local_aligner import align_exon_block_global

    if not read.cigartuples or not read.query_sequence:
        return False

    chrom = read.reference_name
    if not chrom or chrom not in genome:
        return False

    chrom_ref = genome[chrom]
    query_seq = read.query_sequence
    cigar = list(read.cigartuples)

    # ── Fast pre-check: any X op at a homopolymer ref position? ──────────────
    def _homo_run(pos: int) -> int:
        if pos < 0 or pos >= len(chrom_ref):
            return 0
        base = chrom_ref[pos].upper()
        if base == 'N':
            return 0
        left = pos
        while left > 0 and chrom_ref[left - 1].upper() == base:
            left -= 1
        right = pos + 1
        while right < len(chrom_ref) and chrom_ref[right].upper() == base:
            right += 1
        return right - left

    has_homo_x = False
    rp = read.reference_start
    for op, length in cigar:
        if op == 8:  # X mismatch
            for k in range(length):
                if _homo_run(rp + k) >= min_run:
                    has_homo_x = True
                    break
        if has_homo_x:
            break
        if op in _REF_CONSUMING:
            rp += length

    if not has_homo_x:
        return False

    # ── Rebuild CIGAR, re-aligning each exon block ────────────────────────────
    new_cigar: list = []
    rp = read.reference_start
    qp = 0
    modified = False
    i = 0
    n = len(cigar)

    while i < n:
        op, length = cigar[i]

        if op == 4:   # S: soft-clip (query consumed, not ref)
            new_cigar.append((op, length))
            qp += length
            i += 1
            continue
        if op == 5:   # H: hard-clip (neither query nor ref consumed)
            new_cigar.append((op, length))
            i += 1
            continue
        if op == 3:   # N: intron skip (ref consumed, not query)
            new_cigar.append((op, length))
            rp += length
            i += 1
            continue

        # Exon block: collect consecutive M/I/D/X/= ops until N, S, or H
        block_start = i
        block_qp = qp
        block_rp = rp

        while i < n and cigar[i][0] not in (3, 4, 5):
            bop, blen = cigar[i]
            if bop in _QUERY_CONSUMING:
                qp += blen
            if bop in _REF_CONSUMING:
                rp += blen
            i += 1

        Q_block = qp - block_qp
        R_block = rp - block_rp

        if Q_block <= 0 or R_block <= 0:
            new_cigar.extend(cigar[block_start:i])
            continue

        # Skip blocks larger than max_block_bp — O(Q×R) DP too slow in Python
        if R_block > max_block_bp:
            new_cigar.extend(cigar[block_start:i])
            continue

        # Extract sequences for this block
        q_seq = query_seq[block_qp: block_qp + Q_block]
        r_seq = chrom_ref[block_rp: block_rp + R_block]

        if len(q_seq) != Q_block or len(r_seq) != R_block:
            new_cigar.extend(cigar[block_start:i])
            continue

        # Global realignment with homopolymer-aware scoring
        new_ops = align_exon_block_global(
            q_seq, r_seq,
            chrom_ref=chrom_ref,
            ref_offset=block_rp,
            min_run=min_run,
        )

        # Sanity-check: spans must be preserved
        new_q = sum(l for bop, l in new_ops if bop in (0, 1))  # M, I
        new_r = sum(l for bop, l in new_ops if bop in (0, 2))  # M, D
        if new_q != Q_block or new_r != R_block:
            logger.warning(
                "realign_exon_blocks: span mismatch for %s block at rp=%d "
                "(Q %d→%d, R %d→%d) — keeping original",
                read.query_name, block_rp, Q_block, new_q, R_block, new_r,
            )
            new_cigar.extend(cigar[block_start:i])
            continue

        # Compare normalized old vs new (=,X → M for comparison)
        old_norm = _normalize_cigar_ops(cigar[block_start:i])
        if new_ops != old_norm:
            modified = True

        new_cigar.extend(new_ops)

    if modified:
        # Merge adjacent same-type ops at block seams
        merged: list = []
        for op, length in new_cigar:
            if merged and merged[-1][0] == op:
                merged[-1] = (merged[-1][0], merged[-1][1] + length)
            else:
                merged.append([op, length])
        read.cigartuples = [tuple(x) for x in merged]

    return modified


def _normalize_cigar_ops(ops) -> list:
    """Return ops with = (7) and X (8) converted to M (0) and adjacent same-ops merged."""
    norm = [(0 if op in (7, 8) else op, length) for op, length in ops]
    merged: list = []
    for op, length in norm:
        if merged and merged[-1][0] == op:
            merged[-1] = (merged[-1][0], merged[-1][1] + length)
        else:
            merged.append([op, length])
    return [tuple(x) for x in merged]


def softclip_intronic_tail_5prime(
    read: pysam.AlignedSegment,
    clip_boundary: int,
    strand: str,
) -> bool:
    """Soft-clip intronic bases from the 5' end of an intronic-snap read.

    For reads where the aligner mapped into an intron without an N-op,
    ``rescue_3ss_truncation`` corrects ``five_prime_position`` in the TSV but
    leaves the BAM CIGAR unchanged.  This function converts the intronic portion
    to a soft-clip (S) so the read's visible alignment ends at the exon/intron
    boundary while preserving the sequence data in the BAM record.

    Hard-clips (H) are reserved for poly(A) tails.  5' intronic bases become S
    so they remain in ``query_sequence`` and are recoverable downstream.

    The SAM convention ``H? S? <alignment ops> S? H?`` is preserved: any
    pre-existing trailing/leading H is kept in place; the new S op is placed
    immediately inside it.

    Args:
        read:           pysam AlignedSegment — modified in-place.
        clip_boundary:  Exon-side intron boundary.
                        Minus strand: ``intron_start`` (first intron base).
                        Plus strand:  ``intron_end``   (exclusive end of intron).
        strand:         ``'+'`` or ``'-'``.

    Returns:
        True if the CIGAR was modified; False if preconditions were not met.
    """
    if read.is_unmapped or not read.cigartuples:
        return False
    cigar = list(read.cigartuples)

    if strand == '-':
        # The 5' end is at the RIGHT side of the CIGAR (minus strand).
        # reference_end is exclusive; last mapped base = reference_end - 1.
        # Goal: soft-clip until reference_end <= clip_boundary.
        if read.reference_end is None:
            return False

        soft_clip_bases = 0   # query bases that will become S

        # Step 1: extract any existing trailing H.
        # H bases are already absent from query_sequence; keep them as H and
        # re-append after the new S so the CIGAR stays valid (...S H).
        existing_h = 0
        while cigar and cigar[-1][0] == 5:
            existing_h += cigar[-1][1]
            cigar.pop()

        # Step 2: extract any existing trailing S.
        # These bases are already in query_sequence — absorb into the new S.
        existing_s = 0
        while cigar and cigar[-1][0] == 4:
            existing_s += cigar[-1][1]
            cigar.pop()

        # Early exit: nothing to trim and no existing S to re-attach.
        if existing_s == 0 and (
                _cigar_ref_end(read.reference_start, cigar) if cigar else 0
        ) <= clip_boundary:
            return False

        # Step 3: trim ref-consuming ops from the right until
        # reference_end <= clip_boundary.
        while cigar:
            current_ref_end = _cigar_ref_end(read.reference_start, cigar)
            if current_ref_end <= clip_boundary:
                break
            excess_ref = current_ref_end - clip_boundary
            op, length = cigar[-1]
            if op in _REF_CONSUMING:
                trim = min(length, excess_ref)
                if op in _QUERY_CONSUMING:
                    soft_clip_bases += trim
                cigar[-1] = (op, length - trim)
                if cigar[-1][1] == 0:
                    cigar.pop()
            elif op in (1, 4):  # I or S (safety — extracted above)
                soft_clip_bases += length
                cigar.pop()
            else:               # D / N: ref-only, no query bases
                cigar.pop()
            if not cigar:
                return False

        if not cigar:
            return False

        # Step 4: append S (existing + newly soft-clipped), then re-append H.
        # query_sequence is NOT modified — S bases are already present in it.
        total_s = existing_s + soft_clip_bases
        if total_s > 0:
            cigar.append((4, total_s))
        if existing_h > 0:
            cigar.append((5, existing_h))

        read.cigartuples = cigar
        return True

    else:  # plus strand
        # The 5' end is at the LEFT side of the CIGAR.
        # Goal: soft-clip until reference_start >= clip_boundary (= intron_end).
        soft_clip_bases = 0

        # Step 1: extract existing leading H (keep; re-insert before new S).
        existing_h = 0
        while cigar and cigar[0][0] == 5:
            existing_h += cigar[0][1]
            cigar.pop(0)

        # Step 2: extract existing leading S (absorb into new S).
        existing_s = 0
        while cigar and cigar[0][0] == 4:
            existing_s += cigar[0][1]
            cigar.pop(0)

        # Early exit: nothing to trim and no existing S to re-attach.
        if existing_s == 0 and read.reference_start >= clip_boundary:
            return False

        if not cigar:
            return False

        # Step 3: trim ref-consuming ops from the left until
        # reference_start >= clip_boundary.
        ref_pos = read.reference_start
        while cigar:
            op, length = cigar[0]
            if ref_pos >= clip_boundary:
                break
            if op in _REF_CONSUMING:
                deficit = clip_boundary - ref_pos
                trim = min(length, deficit)
                if op in _QUERY_CONSUMING:
                    soft_clip_bases += trim
                cigar[0] = (op, length - trim)
                if cigar[0][1] == 0:
                    cigar.pop(0)
                ref_pos += trim
            elif op in (1, 4):  # I or S (safety — extracted above)
                soft_clip_bases += length
                cigar.pop(0)
            else:               # H shouldn't appear here after step 1
                cigar.pop(0)
            if not cigar:
                return False

        # Step 4: insert H then S at the start; query_sequence unchanged.
        total_s = existing_s + soft_clip_bases
        if total_s > 0:
            cigar.insert(0, (4, total_s))
        if existing_h > 0:
            cigar.insert(0, (5, existing_h))

        if not cigar:
            return False

        read.cigartuples = cigar
        read.reference_start = clip_boundary
        return True


def reroute_intronic_tail_5prime_via_junction(
    read: pysam.AlignedSegment,
    clip_boundary: int,
    five_prime_position: int,
    exon_cigar_str: str,
    strand: str,
) -> bool:
    """Reroute intronic-mapped 5' bases through the splice junction to exon 1.

    For Case 1/2/2b reads where the aligner mapped into an intron without an
    N-op and without a 5' soft-clip, the TSV ``five_prime_position`` is already
    correct but the BAM CIGAR still shows the intronic mapping.  This function:

    1. Trims the CIGAR from the 5' end back to ``clip_boundary``
       (``intron_start`` for ``-``; ``intron_end`` for ``+``), discarding the
       intronic reference-consuming ops without touching ``query_sequence``.
    2. Appends (``-``) or prepends (``+``) N(*intron_len*) + exon ops so the
       same query bases are now assigned to exon-1 reference positions.

    ``query_sequence`` is **not modified** — the intronic query bases are
    simply reassigned to the exon-1 CIGAR ops.

    A sanity check verifies that the exon ops consume the same number of query
    bases as the trimmed intronic portion; the read is left unchanged on mismatch.

    Args:
        read:               pysam AlignedSegment — modified in-place.
        clip_boundary:      Exon-2–side intron boundary.
                            Minus strand: ``intron_start``.
                            Plus strand:  ``intron_end``.
        five_prime_position: Exon-1–side intron boundary.
                            Minus strand: ``intron_end``.
                            Plus strand:  ``intron_start``.
        exon_cigar_str:     SAM CIGAR string for the exon-1 segment
                            (from ``five_prime_exon_cigar`` in the TSV).
        strand:             ``'+'`` or ``'-'``.

    Returns:
        True if the CIGAR was modified; False if preconditions not met.
    """
    if read.is_unmapped or not read.cigartuples or not exon_cigar_str:
        return False
    if read.reference_end is None:
        return False

    try:
        from .local_aligner import cigar_str_to_ops
        exon_ops = cigar_str_to_ops(exon_cigar_str)
    except Exception:
        return False
    if not exon_ops:
        return False

    exon_q_bases = sum(l for op, l in exon_ops if op in _QUERY_CONSUMING)
    cigar = list(read.cigartuples)

    if strand == '-':
        if read.reference_end <= clip_boundary:
            return False  # already at or before boundary

        # Count query bases in the intronic tail (must match exon_ops query span).
        n_intronic_q = 0
        tmp = cigar[:]
        while tmp:
            cur_end = _cigar_ref_end(read.reference_start, tmp)
            if cur_end <= clip_boundary:
                break
            op, length = tmp[-1]
            excess = cur_end - clip_boundary
            if op in _REF_CONSUMING:
                trim = min(length, excess)
                if op in _QUERY_CONSUMING:
                    n_intronic_q += trim
                tmp[-1] = (op, length - trim)
                if tmp[-1][1] == 0:
                    tmp.pop()
            elif op in _QUERY_CONSUMING:
                n_intronic_q += length
                tmp.pop()
            else:
                tmp.pop()

        if n_intronic_q != exon_q_bases:
            return False  # query-length mismatch — don't risk corrupting the read

        # Trim the live CIGAR (no sequence change).
        while cigar:
            cur_end = _cigar_ref_end(read.reference_start, cigar)
            if cur_end <= clip_boundary:
                break
            op, length = cigar[-1]
            excess = cur_end - clip_boundary
            if op in _REF_CONSUMING:
                trim = min(length, excess)
                cigar[-1] = (op, length - trim)
                if cigar[-1][1] == 0:
                    cigar.pop()
            elif op in _QUERY_CONSUMING:
                cigar.pop()
            else:
                cigar.pop()
            if not cigar:
                return False

        if not cigar:
            return False

        intron_len = five_prime_position - _cigar_ref_end(read.reference_start, cigar)
        if intron_len <= 0:
            return False

        cigar.append((3, intron_len))   # N
        cigar.extend(exon_ops)
        read.cigartuples = cigar
        return True

    else:  # plus strand
        if read.reference_start >= clip_boundary:
            return False

        # Count query bases in the intronic head.
        n_intronic_q = 0
        tmp = cigar[:]
        ref_pos = read.reference_start
        while tmp:
            op, length = tmp[0]
            if ref_pos >= clip_boundary:
                break
            if op in _REF_CONSUMING:
                deficit = clip_boundary - ref_pos
                trim = min(length, deficit)
                if op in _QUERY_CONSUMING:
                    n_intronic_q += trim
                tmp[0] = (op, length - trim)
                if tmp[0][1] == 0:
                    tmp.pop(0)
                ref_pos += trim
            elif op in _QUERY_CONSUMING:
                n_intronic_q += length
                tmp.pop(0)
            else:
                tmp.pop(0)

        if n_intronic_q != exon_q_bases:
            return False

        # Trim the live CIGAR.
        ref_pos = read.reference_start
        while cigar:
            op, length = cigar[0]
            if ref_pos >= clip_boundary:
                break
            if op in _REF_CONSUMING:
                deficit = clip_boundary - ref_pos
                trim = min(length, deficit)
                cigar[0] = (op, length - trim)
                if cigar[0][1] == 0:
                    cigar.pop(0)
                ref_pos += trim
            elif op in _QUERY_CONSUMING:
                cigar.pop(0)
            else:
                cigar.pop(0)
            if not cigar:
                return False

        # intron_len = clip_boundary (intron_end) − five_prime_position (intron_start)
        intron_len = clip_boundary - five_prime_position
        if intron_len <= 0:
            return False

        exon_ref_span = sum(l for op, l in exon_ops if op in _REF_CONSUMING)
        new_ref_start = five_prime_position - exon_ref_span

        cigar.insert(0, (3, intron_len))        # N
        for op_tup in reversed(exon_ops):
            cigar.insert(0, op_tup)

        read.cigartuples = cigar
        read.reference_start = new_ref_start
        return True


def _cigar_ref_end(ref_start: int, cigar: list) -> int:
    """Compute the exclusive reference end from a CIGAR list."""
    _ref_consuming = frozenset([0, 2, 3, 7, 8])
    pos = ref_start
    for op, length in cigar:
        if op in _ref_consuming:
            pos += length
    return pos


def extend_read_5prime_for_junction_rescue(
    read: pysam.AlignedSegment,
    five_prime_position: int,
    soft_clip_len: int,
    strand: str,
    exon_cigar_str: str = '',
) -> bool:
    """
    Extend a read's 5' alignment to cover a rescued splice-junction exon.

    For Cat3 reads where :func:`rescue_3ss_truncation` rescued the 5' soft-clip
    to an upstream exon, this converts the soft-clip to exon CIGAR ops and
    inserts an ``N`` op for the intron — making the exon-intron-exon structure
    visible in IGV without any additional annotation track.

    When *exon_cigar_str* is provided (computed by :mod:`local_aligner`), the
    exon segment uses the actual M/I/D alignment ops.  If absent or empty, a
    flat ``nM`` is used as before.

    The intron length is derived from the gap between the current alignment
    boundary and *five_prime_position* (the junction boundary marker):

    * **Plus strand**: intron_len = ``reference_start − five_prime_position − 1``
      (``reference_start`` is the first base of the downstream exon;
      ``five_prime_position`` is the last base of the upstream exon =
      ``intron_start − 1``).
    * **Minus strand**: intron_len = ``five_prime_position − reference_end``
      (``reference_end`` is ``intron_start``; ``five_prime_position`` is
      ``intron_end``).

    Args:
        read:               pysam AlignedSegment — modified in-place.
        five_prime_position: Junction boundary marker from the TSV
                            (``intron_start − 1`` for ``+``; ``intron_end``
                            for ``−``).
        soft_clip_len:      Number of 5'-end soft-clip bases to convert.
                            If the actual leading/trailing ``S`` op length
                            differs, the actual length is used.
        strand:             ``'+'`` or ``'-'``.
        exon_cigar_str:     Optional SAM CIGAR string for the exon segment
                            (e.g. ``"8M1D3M"``).  When provided, its M/I/D ops
                            replace the flat ``nM``.  Falls back to flat ``nM``
                            if empty or unparseable.

    Returns:
        True if the CIGAR was modified; False if preconditions were not met.
    """
    cigar = list(read.cigartuples or [])
    if not cigar or read.is_unmapped or soft_clip_len <= 0:
        return False

    # Parse exon CIGAR if provided.
    exon_ops: Optional[list] = None
    if exon_cigar_str:
        try:
            from .local_aligner import cigar_str_to_ops
            _parsed = cigar_str_to_ops(exon_cigar_str)
            if _parsed:
                exon_ops = _parsed
        except Exception:
            pass  # fall back to flat M

    if strand == '+':
        # Leading soft-clip is the 5' end for plus-strand reads.
        if not cigar or cigar[0][0] != 4:  # S
            return False
        actual_sc = cigar.pop(0)[1]
        n = actual_sc

        intron_len = read.reference_start - five_prime_position - 1

        # Build exon CIGAR ops (M/I/D).  Reference bases consumed by exon ops
        # must equal (five_prime_position + 1) - exon_ref_start, but we can
        # derive the ref start from intron_len and exon ops ref span.
        if exon_ops is None:
            exon_ops = [(0, n)]  # flat M fallback

        # ref bases consumed by exon ops
        _ref_consuming_exon = frozenset([0, 2, 7, 8])  # M, D, =, X
        exon_ref_span = sum(l for op, l in exon_ops if op in _ref_consuming_exon)
        # New reference_start = five_prime_position - exon_ref_span + 1
        new_ref_start = five_prime_position - exon_ref_span + 1

        if intron_len <= 0:
            # No intron gap — just prepend exon ops.
            for op_tup in reversed(exon_ops):
                cigar.insert(0, op_tup)
            read.cigartuples = cigar
            read.reference_start = new_ref_start
            return True

        # Prepend: exon_ops + N(intron_len)
        cigar.insert(0, (3, intron_len))  # N
        for op_tup in reversed(exon_ops):
            cigar.insert(0, op_tup)
        read.cigartuples = cigar
        read.reference_start = new_ref_start
        return True

    else:  # minus strand — trailing soft-clip is the 5' end
        if not cigar or cigar[-1][0] != 4:  # S
            return False
        actual_sc = cigar.pop()[1]
        n = actual_sc

        intron_len = five_prime_position - read.reference_end

        if exon_ops is None:
            exon_ops = [(0, n)]  # flat M fallback

        if intron_len <= 0:
            for op_tup in exon_ops:
                cigar.append(op_tup)
            read.cigartuples = cigar
            return True

        # Append: N(intron_len) then exon_ops
        cigar.append((3, intron_len))  # N
        for op_tup in exon_ops:
            cigar.append(op_tup)
        read.cigartuples = cigar
        # reference_end is recomputed automatically by pysam from the new cigar
        return True


def extend_read_3prime_for_softclip_rescue(
    read: pysam.AlignedSegment,
    strand: str,
    homopolymer_extension: int,
    rescued_seq: str,
    original_softclip_len: int,
    hard_clip: bool = True,
) -> bool:
    """
    Extend a read's 3' alignment to cover a rescued soft-clip at a homopolymer boundary.

    For Cat2 reads where :func:`rescue_softclip_at_homopolymer` rescued the 3'
    soft-clip at an under-called homopolymer, this converts the soft-clip to a
    DELETION + MATCH structure, making the true RNA 3' end visible in IGV.

    The resulting CIGAR structure is:

    * **Plus strand** (3' is rightmost): append ``{D}D {M}M {poly_a}S|H``
      where ``D`` = under-called T's, ``M`` = len(rescued_seq), ``poly_a`` =
      remaining soft-clip after rescued bases.
    * **Minus strand** (3' is leftmost): prepend ``{poly_a}S|H {D}D {M}M``
      and shift ``reference_start`` leftward by ``D + M``.

    In *hard_clip* mode (the default, for ``rectified_pA_hardclip.bam``), the
    poly-A/T tail bases are removed from ``query_sequence`` and encoded as ``H``
    so they are invisible in IGV pileups.  In soft-clip mode (for
    ``rectified_pA_softclip.bam``), the tail remains in the sequence as ``S``.

    Args:
        read:                  pysam AlignedSegment — modified in-place.
        strand:                ``'+'`` or ``'-'``.
        homopolymer_extension: Number of under-called homopolymer bases
                               (reference DELETION between alignment end and
                               rescued bases).
        rescued_seq:           Sequence of non-poly-A/T bases from the soft-clip
                               that match the reference (the new M segment).
        original_softclip_len: Total length of the original 3' soft-clip,
                               used to compute the poly-A/T tail length
                               (``original_softclip_len - len(rescued_seq)``).
        hard_clip:             If True, poly-A/T bases are hard-clipped (H) and
                               removed from ``query_sequence``.  If False, they
                               are soft-clipped (S) and remain in the sequence.

    Returns:
        True if the CIGAR was modified; False if preconditions were not met.
    """
    rescued_len = len(rescued_seq)
    if rescued_len == 0 or original_softclip_len < rescued_len:
        return False

    cigar = list(read.cigartuples or [])
    if not cigar or read.is_unmapped:
        return False

    seq   = read.query_sequence
    quals = read.query_qualities

    poly_a_len = original_softclip_len - rescued_len

    if strand == '+':
        # 3' end is at the right (high coordinates).  The trailing soft-clip
        # encodes: [rescued_seq][poly-A tail].  Walk from right to left in
        # query_sequence: last `poly_a_len` bases are poly-A, preceding
        # `rescued_len` bases are the rescued genomic sequence.

        # 1. Strip trailing soft-clip (includes rescued + poly-A).
        if not cigar or cigar[-1][0] != 4:   # must end with S
            return False
        actual_sc = cigar.pop()[1]
        if actual_sc != original_softclip_len:
            # Mismatch between stored metadata and current CIGAR — skip surgery.
            cigar.append((4, actual_sc))
            return False

        # 2. Strip any trailing hard-clips (shouldn't exist, but be safe).
        while cigar and cigar[-1][0] == 5:
            cigar.pop()

        # 3. Append: D(homopolymer_extension) + M(rescued_len)
        if homopolymer_extension > 0:
            cigar.append((2, homopolymer_extension))   # D
        cigar.append((0, rescued_len))                 # M

        # 4. Append poly-A as H or S.
        if poly_a_len > 0:
            if hard_clip:
                cigar.append((5, poly_a_len))          # H
            else:
                cigar.append((4, poly_a_len))          # S

        # 5. For hard-clip: remove poly-A bases from query_sequence BEFORE
        #    setting cigartuples (pysam validates seq length against cigar at
        #    write time; set sequence first to keep them consistent).
        if hard_clip and poly_a_len > 0 and seq is not None:
            new_seq  = seq[:-poly_a_len]
            new_qual = quals[:-poly_a_len] if quals is not None else None
            read.query_sequence  = new_seq
            read.query_qualities = new_qual

        # 6. Apply CIGAR.
        read.cigartuples = cigar

        return True

    else:  # minus strand — 3' end is at the left (low coordinates).
        # The leading soft-clip encodes: [poly-T tail][rescued_seq] in BAM
        # orientation (poly-A in RNA = poly-T in BAM RC).

        # 1. Strip leading soft-clip.
        if not cigar or cigar[0][0] != 4:   # must start with S
            return False
        actual_sc = cigar.pop(0)[1]
        if actual_sc != original_softclip_len:
            cigar.insert(0, (4, actual_sc))
            return False

        # 2. Strip any leading hard-clips.
        while cigar and cigar[0][0] == 5:
            cigar.pop(0)

        # 3. Prepend ops so the final left-to-right CIGAR reads:
        #    [H/S(poly-T)] [M(rescued)] [D(T-tract)] [original alignment]
        #
        # Insert order (each insert goes to front):
        #   a) Insert D → [D, original_alignment]
        #   b) Insert M → [M, D, original_alignment]
        #   c) (step 4) Insert H/S → [H/S, M, D, original_alignment]
        #
        # Reference layout (low → high coords):
        #   corrected_pos … rescued_seq … T-tract (D) … original alignment start
        if homopolymer_extension > 0:
            cigar.insert(0, (2, homopolymer_extension))  # D
        cigar.insert(0, (0, rescued_len))              # M

        # 4. Prepend poly-T as H or S.
        if poly_a_len > 0:
            if hard_clip:
                cigar.insert(0, (5, poly_a_len))       # H
            else:
                cigar.insert(0, (4, poly_a_len))       # S

        # 5. For hard-clip: trim poly-T from front of query_sequence BEFORE
        #    setting cigartuples.
        old_ref_start = read.reference_start
        if hard_clip and poly_a_len > 0 and seq is not None:
            new_seq  = seq[poly_a_len:]
            new_qual = quals[poly_a_len:] if quals is not None else None
            read.query_sequence  = new_seq
            read.query_qualities = new_qual

        # 6. Apply CIGAR; update reference_start.
        read.cigartuples = cigar
        read.reference_start = old_ref_start - homopolymer_extension - rescued_len

        return True


def _load_corrections_from_tsv(corrected_tsv_path: str) -> Dict[str, dict]:
    """
    Parse a corrected TSV and return a per-read correction dict.

    Only the first row per read_id is used (dominant position for Cat6
    multi-peak reads).  Returns a dict mapping read_id → correction info:

        ``corrected_3prime``       int
        ``strand``                 str  ('+' or '-')
        ``five_prime_position``    int or None
        ``five_prime_rescued``     bool
        ``five_prime_soft_clip``   int
    """
    corrections: Dict[str, dict] = {}
    try:
        with open(corrected_tsv_path) as _f:
            hdr = _f.readline().strip().split('\t')
            try:
                i_id     = hdr.index('read_id')
                i_pos    = hdr.index('corrected_3prime')
                i_strand = hdr.index('strand')
            except ValueError as exc:
                raise ValueError(
                    f"Required column missing in {corrected_tsv_path}: {exc}"
                ) from exc

            # Optional columns — indices default to -1 if absent (old TSV).
            i_5p      = hdr.index('five_prime_position')         if 'five_prime_position'         in hdr else -1
            i_5p_res  = hdr.index('five_prime_rescued')          if 'five_prime_rescued'          in hdr else -1
            i_5p_clip = hdr.index('five_prime_soft_clip_length') if 'five_prime_soft_clip_length' in hdr else -1
            i_5p_cig  = hdr.index('five_prime_exon_cigar')       if 'five_prime_exon_cigar'       in hdr else -1
            # Cat2 soft-clip rescue columns (v2.9.1)
            i_sc_ext   = hdr.index('sc_homopolymer_extension')  if 'sc_homopolymer_extension'  in hdr else -1
            i_sc_seq   = hdr.index('sc_rescued_seq')             if 'sc_rescued_seq'             in hdr else -1
            i_sc_sclen = hdr.index('sc_original_softclip_len')  if 'sc_original_softclip_len'  in hdr else -1
            # Case 4 intronic-snap BAM hard-clip column (v2.9.8)
            i_5p_icp   = hdr.index('five_prime_intron_clip_pos') if 'five_prime_intron_clip_pos' in hdr else -1

            for line in _f:
                parts = line.rstrip('\n').split('\t')
                rid = parts[i_id] if len(parts) > i_id else None
                if rid is None or rid in corrections:
                    continue
                try:
                    corr_pos = int(parts[i_pos])
                    strand   = parts[i_strand]
                except (ValueError, IndexError):
                    continue

                five_prime_pos     = int(parts[i_5p])        if i_5p >= 0     and len(parts) > i_5p     and parts[i_5p]     else None
                five_prime_rescued = (parts[i_5p_res] == '1') if i_5p_res >= 0 and len(parts) > i_5p_res else False
                five_prime_sc      = int(parts[i_5p_clip])   if i_5p_clip >= 0 and len(parts) > i_5p_clip and parts[i_5p_clip] else 0
                five_prime_exon_cig = parts[i_5p_cig]        if i_5p_cig >= 0  and len(parts) > i_5p_cig  and parts[i_5p_cig]  else ''
                # Cat2 fields
                sc_ext   = int(parts[i_sc_ext])   if i_sc_ext   >= 0 and len(parts) > i_sc_ext   and parts[i_sc_ext]   else 0
                sc_seq   = parts[i_sc_seq]         if i_sc_seq   >= 0 and len(parts) > i_sc_seq   else ''
                sc_sclen = int(parts[i_sc_sclen])  if i_sc_sclen >= 0 and len(parts) > i_sc_sclen and parts[i_sc_sclen] else 0
                # Case 4 intronic-snap BAM clip position (-1 = not applicable)
                _icp_raw = parts[i_5p_icp] if i_5p_icp >= 0 and len(parts) > i_5p_icp and parts[i_5p_icp] else '-1'
                five_prime_icp = int(_icp_raw) if _icp_raw.lstrip('-').isdigit() else -1

                corrections[rid] = {
                    'corrected_3prime':           corr_pos,
                    'strand':                     strand,
                    'five_prime_position':        five_prime_pos,
                    'five_prime_rescued':         five_prime_rescued,
                    'five_prime_soft_clip':       five_prime_sc,
                    'five_prime_exon_cigar':      five_prime_exon_cig,
                    'five_prime_intron_clip_pos': five_prime_icp,
                    'sc_homopolymer_extension':   sc_ext,
                    'sc_rescued_seq':             sc_seq,
                    'sc_original_softclip_len':   sc_sclen,
                }
    except OSError as exc:
        raise OSError(
            f"Cannot read corrected TSV {corrected_tsv_path}: {exc}"
        ) from exc
    return corrections


def write_corrected_bam(
    input_bam_path: str,
    corrected_tsv_path: str,
    output_bam_path: str,
    genome: Optional[Dict[str, str]] = None,
) -> Dict[str, int]:
    """
    Write a new BAM with every read hard-clipped at its corrected 3' end and,
    for Cat3 junction-rescued reads, extended at its corrected 5' end.

    Reads per-read corrections from *corrected_tsv_path* (the TSV produced by
    ``rectify correct``) and applies:

    0. **Junction homopolymer surgery** (when *genome* provided):
       :func:`fix_homopolymer_mismatches` converts N-X patterns
       at homopolymer / splice-site boundaries to the more parsimonious N+1 form.
    1. **5' extension** (Cat3): :func:`extend_read_5prime_for_junction_rescue`
       converts the 5' soft-clip to M+N+M to show the exon-intron-exon structure.
    2. **3' hard-clip** (all other categories): :func:`clip_read_to_corrected_3prime`.

    Only the first row per read is used (covers Cat6 NET-seq multi-peak reads).
    Reads absent from the TSV are written unchanged.

    Args:
        input_bam_path:      Path to the original input BAM.
        corrected_tsv_path:  Path to the corrected_3ends.tsv from ``rectify correct``.
        output_bam_path:     Destination BAM path.  Overwritten if it exists.
        genome:              Optional pre-loaded genome dict for homopolymer surgery.

    Returns:
        Dict with summary counts: ``'total'``, ``'clipped'``, ``'unchanged'``.
    """
    corrections = _load_corrections_from_tsv(corrected_tsv_path)

    logger.info(
        "write_corrected_bam: loaded %d corrected positions from %s",
        len(corrections), corrected_tsv_path,
    )

    stats: Dict[str, int] = {'total': 0, 'clipped': 0, 'unchanged': 0}

    with pysam.AlignmentFile(input_bam_path, 'rb') as bam_in, \
         pysam.AlignmentFile(output_bam_path, 'wb', header=bam_in.header) as bam_out:

        for read in bam_in:
            stats['total'] += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                bam_out.write(read)
                stats['unchanged'] += 1
                continue

            correction = corrections.get(read.query_name)
            if correction is None:
                bam_out.write(read)
                stats['unchanged'] += 1
                continue

            modified = False

            # Homopolymer CIGAR surgery: re-align each exon block with reduced
            # mismatch penalty at homopolymer ref positions so that nanopore DRS
            # length errors are represented as indels rather than mismatches.
            # Global alignment preserves query and reference span by construction.
            if genome is not None:
                modified |= realign_exon_blocks(read, genome)

            # 5' junction rescue: extend soft-clip to exon 1 (Cat3).
            if correction['five_prime_rescued'] and correction['five_prime_position'] is not None:
                modified |= extend_read_5prime_for_junction_rescue(
                    read,
                    correction['five_prime_position'],
                    correction['five_prime_soft_clip'],
                    correction['strand'],
                    exon_cigar_str=correction.get('five_prime_exon_cigar', ''),
                )

            # 5' junction rescue: reroute intronic M ops to exon 1 (Cases 1/2/2b).
            # Fires for reads with no soft-clip but a non-empty five_prime_exon_cigar
            # (the aligner mapped into the intron using M/X/D ops rather than N).
            _icp = correction.get('five_prime_intron_clip_pos', -1)
            _exon_cig = correction.get('five_prime_exon_cigar', '')
            if (_icp >= 0 and _exon_cig and correction.get('five_prime_rescued')
                    and correction['five_prime_position'] is not None):
                modified |= reroute_intronic_tail_5prime_via_junction(
                    read,
                    clip_boundary=_icp,
                    five_prime_position=correction['five_prime_position'],
                    exon_cigar_str=_exon_cig,
                    strand=correction['strand'],
                )

            # Cat2 soft-clip rescue: extend 3' alignment outward into homopolymer.
            if correction.get('sc_rescued_seq'):
                modified |= extend_read_3prime_for_softclip_rescue(
                    read,
                    correction['strand'],
                    correction['sc_homopolymer_extension'],
                    correction['sc_rescued_seq'],
                    correction['sc_original_softclip_len'],
                    hard_clip=True,
                )

            # 3' correction: hard-clip to corrected position (Cat1/4/5/6).
            # Cat2 reads already have their 3' end extended outward; this call
            # will return False for them (corrected >= current_end).
            modified |= clip_read_to_corrected_3prime(
                read, correction['corrected_3prime'], correction['strand']
            )

            bam_out.write(read)
            if modified:
                stats['clipped'] += 1
            else:
                stats['unchanged'] += 1

    return stats


def write_softclipped_bam(
    input_bam_path: str,
    corrected_tsv_path: str,
    output_bam_path: str,
    genome: Optional[Dict[str, str]] = None,
) -> Dict[str, int]:
    """
    Write a new BAM with every read soft-clipped at its corrected 3' end and,
    for Cat3 junction-rescued reads, extended at its corrected 5' end.

    Identical to :func:`write_corrected_bam` except the 3' correction uses
    :func:`softclip_read_to_corrected_3prime` (``S`` ops, bases retained) rather
    than hard-clip.  The 5' extension always uses ``M`` ops in both variants
    since the rescued exon bases are aligned (not soft-clipped).
    """
    corrections = _load_corrections_from_tsv(corrected_tsv_path)

    logger.info(
        "write_softclipped_bam: loaded %d corrected positions from %s",
        len(corrections), corrected_tsv_path,
    )

    stats: Dict[str, int] = {'total': 0, 'clipped': 0, 'unchanged': 0}

    with pysam.AlignmentFile(input_bam_path, 'rb') as bam_in, \
         pysam.AlignmentFile(output_bam_path, 'wb', header=bam_in.header) as bam_out:

        for read in bam_in:
            stats['total'] += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                bam_out.write(read)
                stats['unchanged'] += 1
                continue

            correction = corrections.get(read.query_name)
            if correction is None:
                bam_out.write(read)
                stats['unchanged'] += 1
                continue

            modified = False

            # Homopolymer CIGAR surgery: re-align exon blocks.
            if genome is not None:
                modified |= realign_exon_blocks(read, genome)

            # 5' junction rescue: extend into upstream exon (Cat3).
            if correction['five_prime_rescued'] and correction['five_prime_position'] is not None:
                modified |= extend_read_5prime_for_junction_rescue(
                    read,
                    correction['five_prime_position'],
                    correction['five_prime_soft_clip'],
                    correction['strand'],
                    exon_cigar_str=correction.get('five_prime_exon_cigar', ''),
                )

            # 5' junction rescue: reroute intronic M ops to exon 1 (Cases 1/2/2b).
            _icp = correction.get('five_prime_intron_clip_pos', -1)
            _exon_cig = correction.get('five_prime_exon_cigar', '')
            if (_icp >= 0 and _exon_cig and correction.get('five_prime_rescued')
                    and correction['five_prime_position'] is not None):
                modified |= reroute_intronic_tail_5prime_via_junction(
                    read,
                    clip_boundary=_icp,
                    five_prime_position=correction['five_prime_position'],
                    exon_cigar_str=_exon_cig,
                    strand=correction['strand'],
                )

            # Cat2 soft-clip rescue: extend 3' alignment outward into homopolymer.
            if correction.get('sc_rescued_seq'):
                modified |= extend_read_3prime_for_softclip_rescue(
                    read,
                    correction['strand'],
                    correction['sc_homopolymer_extension'],
                    correction['sc_rescued_seq'],
                    correction['sc_original_softclip_len'],
                    hard_clip=False,
                )

            # 3' correction: soft-clip to corrected position (Cat1/4/5/6).
            # Cat2 reads are already extended outward; this call returns False.
            modified |= softclip_read_to_corrected_3prime(
                read, correction['corrected_3prime'], correction['strand']
            )

            bam_out.write(read)
            if modified:
                stats['clipped'] += 1
            else:
                stats['unchanged'] += 1

    return stats


def write_dual_bam(
    input_bam_path: str,
    corrected_tsv_path: str,
    output_hardclip_path: str,
    output_softclip_path: str,
    genome: Optional[Dict[str, str]] = None,
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Write hardclip and softclip corrected BAMs in a single pass over the input BAM.

    Equivalent to calling :func:`write_corrected_bam` and
    :func:`write_softclipped_bam` sequentially, but reads the input BAM only
    once — saving ~50% input I/O compared to the two-function approach.

    The two output BAMs differ only in how the 3' poly(A) region is encoded:

    * *hardclip* (``rectified_pA_hardclip.bam``): poly(A) bases removed from
      ``query_sequence`` and stored as ``H`` ops — invisible in IGV pileups.
    * *softclip* (``rectified_pA_softclip.bam``): poly(A) bases retained in
      ``query_sequence`` and stored as ``S`` ops — visible in IGV when
      "Show soft-clipped bases" is enabled.

    The shared pre-pass (homopolymer junction surgery + 5' rescue / Cat3) is
    applied before state diverges, so it is only computed once per read.
    For the diverging 3' operations the read's mutable state (CIGAR, sequence,
    qualities, reference_start) is saved after shared operations and restored
    before the softclip path.

    Args:
        input_bam_path:       Path to the original input BAM.
        corrected_tsv_path:   Path to the corrected_3ends.tsv from ``rectify correct``.
        output_hardclip_path: Destination BAM path for the hardclip variant.
        output_softclip_path: Destination BAM path for the softclip variant.
        genome:               Optional pre-loaded genome dict for homopolymer surgery.

    Returns:
        Tuple of (hardclip_stats, softclip_stats) dicts, each with keys
        ``'total'``, ``'clipped'``, ``'unchanged'``.
    """
    corrections = _load_corrections_from_tsv(corrected_tsv_path)

    logger.info(
        "write_dual_bam: loaded %d corrected positions from %s",
        len(corrections), corrected_tsv_path,
    )

    hc_stats: Dict[str, int] = {'total': 0, 'clipped': 0, 'unchanged': 0}
    sc_stats: Dict[str, int] = {'total': 0, 'clipped': 0, 'unchanged': 0}

    with pysam.AlignmentFile(input_bam_path, 'rb') as bam_in, \
         pysam.AlignmentFile(output_hardclip_path, 'wb', header=bam_in.header) as bam_hc, \
         pysam.AlignmentFile(output_softclip_path, 'wb', header=bam_in.header) as bam_sc:

        for read in bam_in:
            hc_stats['total'] += 1
            sc_stats['total'] += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                bam_hc.write(read)
                bam_sc.write(read)
                hc_stats['unchanged'] += 1
                sc_stats['unchanged'] += 1
                continue

            correction = corrections.get(read.query_name)
            if correction is None:
                bam_hc.write(read)
                bam_sc.write(read)
                hc_stats['unchanged'] += 1
                sc_stats['unchanged'] += 1
                continue

            # Apply shared pre-pass — identical for both BAMs.
            shared_modified = False

            # Homopolymer CIGAR surgery: re-align exon blocks.
            if genome is not None:
                shared_modified |= realign_exon_blocks(read, genome)

            # 5' rescue (Cat3) — identical for both BAMs.
            if correction['five_prime_rescued'] and correction['five_prime_position'] is not None:
                shared_modified |= extend_read_5prime_for_junction_rescue(
                    read,
                    correction['five_prime_position'],
                    correction['five_prime_soft_clip'],
                    correction['strand'],
                    exon_cigar_str=correction.get('five_prime_exon_cigar', ''),
                )

            # 5' junction rescue: reroute intronic M ops to exon 1 (Cases 1/2/2b).
            _icp = correction.get('five_prime_intron_clip_pos', -1)
            _exon_cig = correction.get('five_prime_exon_cigar', '')
            if (_icp >= 0 and _exon_cig and correction.get('five_prime_rescued')
                    and correction['five_prime_position'] is not None):
                shared_modified |= reroute_intronic_tail_5prime_via_junction(
                    read,
                    clip_boundary=_icp,
                    five_prime_position=correction['five_prime_position'],
                    exon_cigar_str=_exon_cig,
                    strand=correction['strand'],
                )

            # Save state at the divergence point (after shared ops, before 3' ops).
            # Only cigar, seq, quals, and reference_start are mutated by any path.
            saved_cigar     = list(read.cigartuples or [])
            saved_seq       = read.query_sequence
            saved_quals_arr = read.query_qualities          # numpy array or None
            saved_refstart  = read.reference_start

            # ── Hardclip path ────────────────────────────────────────────────
            hc_modified = shared_modified
            if correction.get('sc_rescued_seq'):
                hc_modified |= extend_read_3prime_for_softclip_rescue(
                    read,
                    correction['strand'],
                    correction['sc_homopolymer_extension'],
                    correction['sc_rescued_seq'],
                    correction['sc_original_softclip_len'],
                    hard_clip=True,
                )
            hc_modified |= clip_read_to_corrected_3prime(
                read, correction['corrected_3prime'], correction['strand']
            )
            bam_hc.write(read)
            if hc_modified:
                hc_stats['clipped'] += 1
            else:
                hc_stats['unchanged'] += 1

            # ── Restore state before softclip path ───────────────────────────
            read.cigartuples     = saved_cigar
            read.query_sequence  = saved_seq
            read.query_qualities = saved_quals_arr
            read.reference_start = saved_refstart

            # ── Softclip path ────────────────────────────────────────────────
            sc_modified = shared_modified
            if correction.get('sc_rescued_seq'):
                sc_modified |= extend_read_3prime_for_softclip_rescue(
                    read,
                    correction['strand'],
                    correction['sc_homopolymer_extension'],
                    correction['sc_rescued_seq'],
                    correction['sc_original_softclip_len'],
                    hard_clip=False,
                )
            sc_modified |= softclip_read_to_corrected_3prime(
                read, correction['corrected_3prime'], correction['strand']
            )
            bam_sc.write(read)
            if sc_modified:
                sc_stats['clipped'] += 1
            else:
                sc_stats['unchanged'] += 1

    return hc_stats, sc_stats


def write_polya_trimmed_bam(
    input_bam_path: str,
    output_bam_path: str,
    threshold: float = 0.8,
) -> Dict[str, int]:
    """
    Write a new BAM with 3' poly(A) soft-clips removed from each read.

    Iterates all reads (including secondary and supplementary) from
    *input_bam_path*, strips the RNA-3' poly(A) soft-clip from primary reads
    that pass the A-richness threshold, and writes every read to
    *output_bam_path*.  Header and all BAM tags are preserved unchanged.

    Secondary and supplementary reads are written as-is without trimming
    because their soft-clips may have different semantics and their 3' end
    is not independently defined.

    Args:
        input_bam_path:  Path to input BAM (sorted, indexed or not).
        output_bam_path: Destination BAM path.  Will be overwritten if it
                         exists.  Caller is responsible for sorting/indexing
                         the output if needed.
        threshold:       Minimum A (plus) or T (minus) fraction required to
                         consider a 3' soft-clip a poly(A) tail (default 0.8).

    Returns:
        Dict with summary counts:
            'total'        — total reads written
            'trimmed'      — reads whose poly(A) tail was removed
            'bases_trimmed'— total bases removed across all reads
    """
    from .polya_trimmer import trim_polya_from_bam_read

    stats = {'total': 0, 'trimmed': 0, 'bases_trimmed': 0}

    with pysam.AlignmentFile(input_bam_path, 'rb') as bam_in, \
         pysam.AlignmentFile(output_bam_path, 'wb', header=bam_in.header) as bam_out:

        for read in bam_in:
            stats['total'] += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                bam_out.write(read)
                continue

            strand = '-' if read.is_reverse else '+'
            read, n_trimmed = trim_polya_from_bam_read(read, strand, threshold=threshold)

            if n_trimmed > 0:
                stats['trimmed']       += 1
                stats['bases_trimmed'] += n_trimmed

            bam_out.write(read)

    return stats


def write_netseq_assigned_bedgraph(
    corrected_tsv_path: str,
    output_prefix: str,
) -> Dict[str, int]:
    """
    Write strand-specific bedGraph files for NET-seq-assigned 3' ends (Cat6).

    Reads all rows from *corrected_tsv_path* where ``correction_applied``
    contains ``'netseq_refinement'``, accumulates fractional counts per
    ``(chrom, strand, corrected_3prime)`` position, and writes two bedGraph
    files (plus and minus strand) named::

        {output_prefix}.plus.bedgraph
        {output_prefix}.minus.bedgraph

    Counts are **not** RPM-normalised (``normalize_rpm=False``) because the
    fractions already represent proportional signal and a single-sample RPM
    factor is meaningless for the subset of NET-seq-assigned reads.

    Args:
        corrected_tsv_path: Path to ``corrected_3ends.tsv`` produced by
            ``rectify correct``.
        output_prefix: Path prefix for output bedGraph files (no extension).

    Returns:
        Dict with counts per strand: ``{'plus': n_positions, 'minus': n_positions}``.
    """
    from .netseq_output import write_bedgraph as _write_bedgraph

    counts: Dict[Tuple[str, str, int], float] = {}
    n_netseq_rows = 0

    try:
        with open(corrected_tsv_path) as _f:
            hdr = _f.readline().strip().split('\t')
            try:
                i_chrom  = hdr.index('chrom')
                i_strand = hdr.index('strand')
                i_pos    = hdr.index('corrected_3prime')
                i_corr   = hdr.index('correction_applied')
            except ValueError as exc:
                raise ValueError(
                    f"Required column missing in {corrected_tsv_path}: {exc}"
                ) from exc

            i_frac = hdr.index('fraction') if 'fraction' in hdr else -1

            for line in _f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(i_chrom, i_strand, i_pos, i_corr):
                    continue
                corrections_str = parts[i_corr]
                if 'netseq_refinement' not in corrections_str:
                    continue
                try:
                    chrom  = parts[i_chrom]
                    strand = parts[i_strand]
                    pos    = int(parts[i_pos])
                    frac   = float(parts[i_frac]) if i_frac >= 0 else 1.0
                except (ValueError, IndexError):
                    continue
                key = (chrom, strand, pos)
                counts[key] = counts.get(key, 0.0) + frac
                n_netseq_rows += 1
    except OSError as exc:
        raise OSError(
            f"Cannot read corrected TSV {corrected_tsv_path}: {exc}"
        ) from exc

    if n_netseq_rows == 0:
        logger.info(
            "write_netseq_assigned_bedgraph: no netseq_refinement rows in %s; "
            "bedgraph files not written.",
            corrected_tsv_path,
        )
        return {'plus': 0, 'minus': 0}

    result_counts: Dict[str, int] = {}
    for strand_char, strand_name in [('+', 'plus'), ('-', 'minus')]:
        bg_path = Path(f"{output_prefix}.{strand_name}.bedgraph")
        _write_bedgraph(
            counts,
            bg_path,
            strand_char,
            total_reads=1,           # RPM normalisation disabled
            normalize_rpm=False,
            track_name=bg_path.stem,
        )
        n_pos = sum(1 for (_, s, _) in counts if s == strand_char)
        logger.info(
            "  Wrote %s (%d positions, %.1f fractional reads)",
            bg_path,
            n_pos,
            sum(v for (_, s, _), v in counts.items() if s == strand_char),
        )
        result_counts[strand_name] = n_pos

    return result_counts


def write_corrected_3ends_bedgraph(
    corrected_tsv_path: str,
    output_prefix: str,
) -> Dict[str, int]:
    """
    Write strand-specific bedGraph files for all corrected 3' ends (Cat1–6).

    Reads every row from *corrected_tsv_path*, accumulates fractional counts
    per ``(chrom, strand, corrected_3prime)`` position using the ``fraction``
    column (defaults to 1.0 when absent), and writes two bedGraph files::

        {output_prefix}.plus.bedgraph
        {output_prefix}.minus.bedgraph

    Cat6 multi-peak rows contribute their fractional values (< 1.0) so the
    bedgraph shows proportional signal at each peak.  All other reads
    contribute 1.0.

    Counts are **not** RPM-normalised.

    Returns:
        Dict with counts per strand: ``{'plus': n_positions, 'minus': n_positions}``.
    """
    from .netseq_output import write_bedgraph as _write_bedgraph

    counts: Dict[Tuple[str, str, int], float] = {}
    n_rows = 0

    try:
        with open(corrected_tsv_path) as _f:
            hdr = _f.readline().strip().split('\t')
            try:
                i_chrom  = hdr.index('chrom')
                i_strand = hdr.index('strand')
                i_pos    = hdr.index('corrected_3prime')
            except ValueError as exc:
                raise ValueError(
                    f"Required column missing in {corrected_tsv_path}: {exc}"
                ) from exc

            i_frac = hdr.index('fraction') if 'fraction' in hdr else -1

            for line in _f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(i_chrom, i_strand, i_pos):
                    continue
                try:
                    chrom  = parts[i_chrom]
                    strand = parts[i_strand]
                    pos    = int(parts[i_pos])
                    frac   = float(parts[i_frac]) if i_frac >= 0 else 1.0
                except (ValueError, IndexError):
                    continue
                key = (chrom, strand, pos)
                counts[key] = counts.get(key, 0.0) + frac
                n_rows += 1
    except OSError as exc:
        raise OSError(
            f"Cannot read corrected TSV {corrected_tsv_path}: {exc}"
        ) from exc

    if n_rows == 0:
        logger.info(
            "write_corrected_3ends_bedgraph: no rows in %s; "
            "bedgraph files not written.",
            corrected_tsv_path,
        )
        return {'plus': 0, 'minus': 0}

    result_counts: Dict[str, int] = {}
    for strand_char, strand_name in [('+', 'plus'), ('-', 'minus')]:
        bg_path = Path(f"{output_prefix}.{strand_name}.bedgraph")
        _write_bedgraph(
            counts,
            bg_path,
            strand_char,
            total_reads=1,
            normalize_rpm=False,
            track_name=bg_path.stem,
        )
        n_pos = sum(1 for (_, s, _) in counts if s == strand_char)
        logger.info(
            "  Wrote %s (%d positions, %.1f fractional reads)",
            bg_path,
            n_pos,
            sum(v for (_, s, _), v in counts.items() if s == strand_char),
        )
        result_counts[strand_name] = n_pos

    return result_counts
