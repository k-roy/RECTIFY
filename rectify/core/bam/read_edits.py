#!/usr/bin/env python3
"""
Per-read CIGAR/sequence editing primitives for RECTIFY BAM writers.

Pure functions of ``(pysam.AlignedSegment, args...) -> bool``: each modifies
the read in place and returns whether anything changed. Composable toolkit
used by the file-level writers in ``bam_writer.py``:

- ``clip_read_to_corrected_3prime``           — hard-clip 3' end
- ``softclip_read_to_corrected_3prime``       — soft-clip 3' end
- ``extend_read_5prime_for_junction_rescue``  — Cat3: convert 5' S to exon + N
- ``extend_read_3prime_for_softclip_rescue``  — Cat2: D + M + H/S at homopolymer
- ``softclip_intronic_tail_5prime``           — clip 5' intronic tail to S
- ``reroute_intronic_tail_5prime_via_junction`` — reroute via N to exon 1
- ``fix_homopolymer_mismatches``              — (deprecated) X → I/D at HP
- ``realign_exon_blocks``                     — global NW per exon block
- ``_hardclip_trailing_a_run``                — extra A-run trim for HC BAM
- ``_normalize_cigar_ops`` / ``_cigar_ref_end`` — small CIGAR helpers
"""

from typing import Dict, Optional
import logging
import pysam

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

        # 3.5. Strip trailing D/N/I ops left dangling by the walk loop.
        # D/N: corrected_3prime falls inside a deletion/intron span — D/N before H
        #      is invalid SAM.
        # I:   corrected_3prime falls immediately after an insertion — I before H is
        #      invalid (an alignment cannot end in an insertion; absorb those query
        #      bases into the hard-clip).
        while cigar and cigar[-1][0] in (1, 2, 3):  # I=1, D=2, N=3
            op, length = cigar.pop()
            if op == 1:  # I is query-consuming
                n_query_remove += length
        if not cigar:
            return False

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

        # 3.5. Strip leading D/N/I ops dangling after the walk (minus strand).
        # I: absorb into hard-clip; D/N: invalid at start of CIGAR.
        while cigar and cigar[0][0] in (1, 2, 3):  # I=1, D=2, N=3
            op, length = cigar.pop(0)
            if op == 1:
                n_query_remove += length
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

        # 3.5. Strip trailing D/N/I dangling after the walk (plus strand softclip).
        # I: absorb into soft-clip (query bases kept in sequence).
        # D/N: invalid immediately before a soft-clip.
        while cigar and cigar[-1][0] in (1, 2, 3):  # I=1, D=2, N=3
            op, length = cigar.pop()
            if op == 1:
                n_query_remove += length
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

        # 3.5. Strip leading D/N/I dangling after the walk (minus strand softclip).
        # I: absorb into soft-clip (query bases kept in sequence).
        # D/N: invalid immediately after a leading soft-clip.
        while cigar and cigar[0][0] in (1, 2, 3):  # I=1, D=2, N=3
            op, length = cigar.pop(0)
            if op == 1:
                n_query_remove += length
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

    .. deprecated::
        This function is no longer called by any write path.  It is
        superseded by ``realign_exon_blocks``, which performs global
        NW re-alignment with HP-aware scoring.  Retained for reference
        only; do not use in new code.
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
    from ..align.local_aligner import align_exon_block_global

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
        from ..align.local_aligner import cigar_str_to_ops
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
        # Merge consecutive N ops (can arise when an existing N is partially
        # trimmed and a new N is appended — e.g. Case 5 n_boundary_adjust).
        merged: list = []
        for _op, _l in cigar:
            if merged and merged[-1][0] == _op == 3:
                merged[-1] = (3, merged[-1][1] + _l)
            else:
                merged.append((_op, _l))
        read.cigartuples = merged
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

        # Merge consecutive N ops (Case 5 n_boundary_adjust: existing N may abut
        # the newly prepended N when the N's exon-2 start is used as clip_boundary).
        merged_p: list = []
        for _op, _l in cigar:
            if merged_p and merged_p[-1][0] == _op == 3:
                merged_p[-1] = (3, merged_p[-1][1] + _l)
            else:
                merged_p.append((_op, _l))
        cigar = merged_p

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
            from ..align.local_aligner import cigar_str_to_ops
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
        _query_consuming_exon = frozenset([0, 1, 4, 7, 8])  # M, I, S, =, X
        exon_ref_span = sum(l for op, l in exon_ops if op in _ref_consuming_exon)
        exon_query_span = sum(l for op, l in exon_ops if op in _query_consuming_exon)
        # Guard: exon_ops must not consume more query bases than the soft-clip.
        # If they do (local aligner artefact, e.g. one extra M from a D at the
        # boundary), fall back to a flat M so the BAM record stays valid.
        if exon_query_span != actual_sc:
            exon_ops = [(0, actual_sc)]
            exon_ref_span = actual_sc
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

        _query_consuming_exon = frozenset([0, 1, 4, 7, 8])  # M, I, S, =, X
        exon_query_span = sum(l for op, l in exon_ops if op in _query_consuming_exon)
        if exon_query_span != actual_sc:
            exon_ops = [(0, actual_sc)]

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

    In *hard_clip* mode (the default, for ``rectified_corrected_3end.bam``), the
    poly-A/T tail bases are removed from ``query_sequence`` and encoded as ``H``
    so they are invisible in IGV pileups.  In soft-clip mode (for
    ``rectified_pA_tail_trimmed.bam``), the tail remains in the sequence as ``S``.

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


def _hardclip_trailing_a_run(
    read: pysam.AlignedSegment,
    strand: str,
) -> bool:
    """
    Additional hard-clip step for the poly(A) hard-clip BAM.

    After the main ``clip_read_to_corrected_3prime`` call, genomic A-rich
    regions at the RNA 3' end (e.g. A-rich 3' UTRs) may remain in the
    alignment.  These are indistinguishable from poly(A) tail in the read
    sequence, so the hard-clip BAM removes them to show only the unambiguous
    exon body.

    For − strand: counts consecutive leading T's in the BAM sequence (= A's
    in RNA orientation), walks the CIGAR to compute the corresponding
    reference offset, and calls :func:`clip_read_to_corrected_3prime` for
    the extra clip.
    For + strand: counts trailing A's and clips symmetrically.

    This step is applied to the *hard-clip* BAM only.  The soft-clip BAM
    retains genomic A's as aligned bases because they are part of the
    transcript sequence (relevant for comparing sequence- vs. signal-based
    poly(A) length estimates).

    Returns True if additional clipping was applied; False otherwise.
    """
    if read.is_unmapped or not read.cigartuples or not read.query_sequence:
        return False

    seq   = read.query_sequence
    cigar = list(read.cigartuples)

    if strand == '-':
        # Count leading T's in BAM seq (= A's in RNA = trailing poly-A zone).
        n_t = 0
        for ch in seq:
            if ch.upper() == 'T':
                n_t += 1
            else:
                break
        if n_t == 0:
            return False

        # Walk CIGAR from left, consuming n_t query bases, tracking ref_pos.
        # H ops do not appear in query_sequence and are skipped.
        ref_pos    = read.reference_start
        q_consumed = 0
        for op, length in cigar:
            if op == 5:          # H — not in query_sequence
                continue
            if op in (0, 7, 8):  # M, =, X — both ref and query
                step = min(length, n_t - q_consumed)
                q_consumed += step
                ref_pos    += step
                if q_consumed >= n_t:
                    break
            elif op == 1:        # I — query only
                step = min(length, n_t - q_consumed)
                q_consumed += step
                if q_consumed >= n_t:
                    break
            elif op == 4:        # S — query only (unlikely in HC path but guard)
                step = min(length, n_t - q_consumed)
                q_consumed += step
                if q_consumed >= n_t:
                    break
            elif op in (2, 3):   # D, N — ref only, advance ref without query
                ref_pos += length

        new_3prime = ref_pos
        if new_3prime <= read.reference_start:
            return False
        return clip_read_to_corrected_3prime(read, new_3prime, strand)

    else:  # plus strand
        # Count trailing A's in BAM seq (= A's in RNA = trailing poly-A zone).
        n_a = 0
        for ch in reversed(seq):
            if ch.upper() == 'A':
                n_a += 1
            else:
                break
        if n_a == 0:
            return False

        # Walk CIGAR from right, consuming n_a query bases, tracking ref offset.
        ref_end    = read.reference_end - 1  # 0-based inclusive right boundary
        q_consumed = 0
        for op, length in reversed(cigar):
            if op == 5:          # H
                continue
            if op in (0, 7, 8):  # M, =, X
                step = min(length, n_a - q_consumed)
                q_consumed += step
                ref_end    -= step
                if q_consumed >= n_a:
                    break
            elif op == 1:        # I
                step = min(length, n_a - q_consumed)
                q_consumed += step
                if q_consumed >= n_a:
                    break
            elif op == 4:        # S
                step = min(length, n_a - q_consumed)
                q_consumed += step
                if q_consumed >= n_a:
                    break
            elif op in (2, 3):   # D, N
                ref_end -= length

        new_3prime = ref_end
        if new_3prime >= read.reference_end - 1:
            return False
        return clip_read_to_corrected_3prime(read, new_3prime, strand)
