#!/usr/bin/env python3
"""
Semi-global affine-gap aligner for RECTIFY Cat3 5' junction rescue.

Two alignment variants are provided:

  _align_right_anchored — free prefix in ref, right end of ref is fixed.
      Used for plus-strand reads where the clip must end at intron_start.

  _align_left_anchored — left end of ref is fixed, free suffix in ref.
      Used for minus-strand reads where the clip must start at intron_end.

Public API
----------
  align_clip_to_exon(clip_seq, genome_seq, intron_start, intron_end, strand)
    → (cigar_ops, exon_ref_start)

  cigar_ops_to_str(ops)  → SAM CIGAR string (e.g. "8M1I3M")
  cigar_str_to_ops(s)    → list of (op_code, length) tuples

Scoring (affine gap, Gotoh 1982):
  match       = +2
  mismatch    = -4
  gap_open    = -4  (paid once when opening a new gap, in addition to extend)
  gap_extend  = -1  (paid per base in the gap)

  Total cost for a gap of length k: gap_open + k * gap_extend

  Rationale: with linear gaps, isolated single-base deletions separated by
  matches can outscore a single consolidated deletion of the same total length
  because intermediate matches flip the score positive.  Affine gap makes
  multiple gap-open events expensive, so e.g. 3D scores -7 while three separate
  1D ops separated by matches score -9 or worse.

Empirical Calibration Notes
----------------------------
Two parameters in ``align_exon_block_global`` are currently fixed heuristics:

  homo_mismatch = -2.0  (reduced mismatch cost at HP positions, default)
  min_run       =  3    (HP run length threshold to apply homo_mismatch)

These could be tuned using the empirical substitution rates from
``empirical_cigar_error_profiler.py``.  The profiler's ``X`` (substitution)
rows give base-class-specific sub rates per HP length; the corresponding
mismatch score would be::

    sub_rate(base_class, hp) / sub_rate(base_class, hp=1) * _MISMATCH

In practice, empirical sub rates are nearly flat with HP length (the dominant
HP error is deletion, not substitution), so ``homo_mismatch=-2`` is a
reasonable approximation.  Re-evaluate if a future dataset shows systematic
mismatch enrichment at HP positions.

The gap scoring constants (_GAP_OPEN, _GAP_EXTEND) define a linear cost per
gap length.  The HP-context del costs from ``HpPenaltyTable`` use a different
model (per-position cost as a function of HP run length, not gap length).  These
are not directly interchangeable.  If HP-specific gap costs are needed here,
the DP would need to be restructured to look up del cost per ref position rather
than per gap-event.  See ``junction_refiner._score_hp_anchored`` for the
per-position HP del cost model.

Author: Kevin R. Roy
Date: 2026-04-11
"""

import re
from typing import List, Tuple
import logging

logger = logging.getLogger(__name__)

# Scoring constants
_MATCH      =  2
_MISMATCH   = -4
_GAP_OPEN   = -4   # penalty for opening a new gap (paid once per gap event)
_GAP_EXTEND = -1   # penalty per base within a gap

_NEG_INF = float('-inf')

# CIGAR op codes (pysam / SAM convention)
_OP_M = 0  # sequence match or mismatch
_OP_I = 1  # insertion in query (gap in reference)
_OP_D = 2  # deletion from query (gap in query)

# Traceback source codes for the H (match) matrix
_TBH_H = 1  # came from H (match/mismatch)
_TBH_D = 2  # came from D state (deletion was ending)
_TBH_I = 3  # came from I state (insertion was ending)

# Traceback codes for D and I matrices
_TBX_OPEN   = 1  # gap was opened from H
_TBX_EXTEND = 2  # gap was extended


def _compress(ops: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge adjacent same-type CIGAR ops into a single (op, length) tuple."""
    if not ops:
        return ops
    merged: List[Tuple[int, int]] = [ops[0]]
    for op, length in ops[1:]:
        if op == merged[-1][0]:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))
    return merged


def _align_right_anchored(
    query: str,
    ref: str,
) -> Tuple[List[Tuple[int, int]], int]:
    """
    Semi-global affine-gap NW: query fully consumed, ref has a free prefix.

    The alignment is anchored to the RIGHT end of *ref*.  The query must
    consume all of its bases; the alignment may begin at any position within
    ref without penalty (free-prefix initialisation).

    Uses the Gotoh three-matrix algorithm (H / D / I) for affine gap scoring,
    which prevents the "staircase" artifact where many isolated single-base
    deletions outscore a single consolidated deletion.

    Used for **plus-strand** reads where the 5' soft-clip must end exactly at
    ``intron_start`` (the last base of exon 1).

    Args:
        query: Query sequence (5'→3' alignment orientation).
        ref:   Reference sequence; ``ref[-1]`` is the last exon base before
               the intron donor.

    Returns:
        ``(cigar_ops, ref_skip)`` where *cigar_ops* is the compressed alignment
        CIGAR and *ref_skip* is the number of leading ref bases NOT consumed
        (so the alignment begins at ``ref[ref_skip]``).
    """
    Q, R = len(query), len(ref)
    if Q == 0:
        return [], R
    if R == 0:
        return [(_OP_I, Q)], 0

    # Three DP matrices.
    # H[i][j]: best score for query[0:i] aligned to ref[0:j], last op = M/X
    # D[i][j]: best score ending with a deletion from query (consuming ref, not query)
    # I_[i][j]: best score ending with an insertion in query (consuming query, not ref)
    H  = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]
    D  = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]
    I_ = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]

    tbH  = [[0] * (R + 1) for _ in range(Q + 1)]
    tbD  = [[0] * (R + 1) for _ in range(Q + 1)]
    tbI_ = [[0] * (R + 1) for _ in range(Q + 1)]

    # Free prefix in ref: alignment can start at any ref position for free.
    # H[0][j] = 0: the alignment "starts here" at ref column j with no prior cost.
    # D[0][j] = 0: being in deletion state at the start of alignment is free
    #              (models skipping the ref prefix without any gap cost).
    for j in range(R + 1):
        H[0][j] = 0.0
        D[0][j] = 0.0

    # Column 0 (no ref consumed): only insertions are possible.
    for i in range(1, Q + 1):
        I_[i][0] = _GAP_OPEN + i * _GAP_EXTEND
        tbI_[i][0] = _TBX_OPEN if i == 1 else _TBX_EXTEND

    # Main DP.
    for i in range(1, Q + 1):
        qi = query[i - 1].upper()
        for j in range(1, R + 1):
            s = _MATCH if qi == ref[j - 1].upper() else _MISMATCH

            # H: match/mismatch — diagonal move from any prior state.
            h_h = H[i-1][j-1]
            h_d = D[i-1][j-1]
            h_i = I_[i-1][j-1]
            best_prev = max(h_h, h_d, h_i)
            if best_prev == _NEG_INF:
                H[i][j] = _NEG_INF
            else:
                H[i][j] = best_prev + s
                if best_prev == h_h:
                    tbH[i][j] = _TBH_H
                elif best_prev == h_d:
                    tbH[i][j] = _TBH_D
                else:
                    tbH[i][j] = _TBH_I

            # D: deletion from query (move left: j decreases, i unchanged).
            # Opening: pay gap_open + gap_extend to start a new gap from H.
            # Extending: pay only gap_extend to continue an existing D gap.
            h_prev = H[i][j-1]
            d_prev = D[i][j-1]
            d_open   = (h_prev + _GAP_OPEN + _GAP_EXTEND) if h_prev != _NEG_INF else _NEG_INF
            d_extend = (d_prev + _GAP_EXTEND) if d_prev != _NEG_INF else _NEG_INF
            if d_open >= d_extend:
                D[i][j] = d_open
                tbD[i][j] = _TBX_OPEN
            else:
                D[i][j] = d_extend
                tbD[i][j] = _TBX_EXTEND

            # I: insertion in query (move up: i decreases, j unchanged).
            h_prev2 = H[i-1][j]
            i_prev  = I_[i-1][j]
            i_open   = (h_prev2 + _GAP_OPEN + _GAP_EXTEND) if h_prev2 != _NEG_INF else _NEG_INF
            i_extend = (i_prev + _GAP_EXTEND) if i_prev != _NEG_INF else _NEG_INF
            if i_open >= i_extend:
                I_[i][j] = i_open
                tbI_[i][j] = _TBX_OPEN
            else:
                I_[i][j] = i_extend
                tbI_[i][j] = _TBX_EXTEND

    # Choose the best ending state at (Q, R).
    end_h, end_d, end_i = H[Q][R], D[Q][R], I_[Q][R]
    best_end = max(end_h, end_d, end_i)
    if end_h == best_end:
        cur_state = 'H'
    elif end_d == best_end:
        cur_state = 'D'
    else:
        cur_state = 'I'

    # Traceback from (Q, R) until i == 0.
    ops: List[Tuple[int, int]] = []
    i, j = Q, R

    while i > 0:
        if j == 0:
            # No ref remaining — emit all remaining query bases as insertions.
            ops.append((_OP_I, i))
            i = 0
            break

        if cur_state == 'H':
            ops.append((_OP_M, 1))
            src = tbH[i][j]
            i -= 1
            j -= 1
            cur_state = 'H' if src == _TBH_H else ('D' if src == _TBH_D else 'I')

        elif cur_state == 'D':
            # Deletion: consume one ref base without consuming query.
            ops.append((_OP_D, 1))
            src = tbD[i][j]
            j -= 1
            cur_state = 'H' if src == _TBX_OPEN else 'D'

        else:  # 'I'
            # Insertion: consume one query base without consuming ref.
            ops.append((_OP_I, 1))
            src = tbI_[i][j]
            i -= 1
            cur_state = 'H' if src == _TBX_OPEN else 'I'

    ref_skip = j  # leading ref bases not consumed by the alignment
    ops.reverse()
    return _compress(ops), ref_skip


def _align_left_anchored(
    query: str,
    ref: str,
) -> Tuple[List[Tuple[int, int]], int]:
    """
    Semi-global affine-gap NW: query fully consumed, left end of ref is fixed,
    free suffix.

    The alignment must start at ``ref[0]`` (the exon–intron boundary); the best
    ending column in ``ref`` is chosen via free-suffix selection.

    Uses the Gotoh three-matrix algorithm for affine gap scoring.

    Used for **minus-strand** reads where the 5' soft-clip must begin exactly
    at ``intron_end`` (the first base of exon 1 in genomic coordinates).

    Args:
        query: Query sequence.
        ref:   Reference sequence; ``ref[0]`` is the first exon base after the
               intron acceptor.

    Returns:
        ``(cigar_ops, ref_consumed)`` where *cigar_ops* is the compressed
        alignment CIGAR and *ref_consumed* is the number of ref bases consumed.
    """
    Q, R = len(query), len(ref)
    if Q == 0:
        return [], 0
    if R == 0:
        return [(_OP_I, Q)], 0

    H  = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]
    D  = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]
    I_ = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]

    tbH  = [[0] * (R + 1) for _ in range(Q + 1)]
    tbD  = [[0] * (R + 1) for _ in range(Q + 1)]
    tbI_ = [[0] * (R + 1) for _ in range(Q + 1)]

    # Standard NW boundary (fixed left anchor — must start from position 0).
    H[0][0] = 0.0
    for j in range(1, R + 1):
        D[0][j] = _GAP_OPEN + j * _GAP_EXTEND   # j ref bases skipped as deletion
        tbD[0][j] = _TBX_OPEN if j == 1 else _TBX_EXTEND
    for i in range(1, Q + 1):
        I_[i][0] = _GAP_OPEN + i * _GAP_EXTEND  # i query bases inserted before any ref
        tbI_[i][0] = _TBX_OPEN if i == 1 else _TBX_EXTEND

    # Main DP.
    for i in range(1, Q + 1):
        qi = query[i - 1].upper()
        for j in range(1, R + 1):
            s = _MATCH if qi == ref[j - 1].upper() else _MISMATCH

            # H: match/mismatch.
            h_h = H[i-1][j-1]
            h_d = D[i-1][j-1]
            h_i = I_[i-1][j-1]
            best_prev = max(h_h, h_d, h_i)
            if best_prev == _NEG_INF:
                H[i][j] = _NEG_INF
            else:
                H[i][j] = best_prev + s
                if best_prev == h_h:
                    tbH[i][j] = _TBH_H
                elif best_prev == h_d:
                    tbH[i][j] = _TBH_D
                else:
                    tbH[i][j] = _TBH_I

            # D: deletion from query.
            h_prev = H[i][j-1]
            d_prev = D[i][j-1]
            d_open   = (h_prev + _GAP_OPEN + _GAP_EXTEND) if h_prev != _NEG_INF else _NEG_INF
            d_extend = (d_prev + _GAP_EXTEND) if d_prev != _NEG_INF else _NEG_INF
            if d_open >= d_extend:
                D[i][j] = d_open
                tbD[i][j] = _TBX_OPEN
            else:
                D[i][j] = d_extend
                tbD[i][j] = _TBX_EXTEND

            # I: insertion in query.
            h_prev2 = H[i-1][j]
            i_prev  = I_[i-1][j]
            i_open   = (h_prev2 + _GAP_OPEN + _GAP_EXTEND) if h_prev2 != _NEG_INF else _NEG_INF
            i_extend = (i_prev + _GAP_EXTEND) if i_prev != _NEG_INF else _NEG_INF
            if i_open >= i_extend:
                I_[i][j] = i_open
                tbI_[i][j] = _TBX_OPEN
            else:
                I_[i][j] = i_extend
                tbI_[i][j] = _TBX_EXTEND

    # Free suffix: pick the ending column with the best overall score.
    j_best = max(
        range(R + 1),
        key=lambda j: max(H[Q][j], D[Q][j], I_[Q][j])
    )
    best_end = max(H[Q][j_best], D[Q][j_best], I_[Q][j_best])
    if H[Q][j_best] == best_end:
        cur_state = 'H'
    elif D[Q][j_best] == best_end:
        cur_state = 'D'
    else:
        cur_state = 'I'

    # Traceback.
    ops: List[Tuple[int, int]] = []
    i, j = Q, j_best

    while i > 0 or j > 0:
        if i == 0:
            # Only ref left — emit deletions for remaining j.
            ops.append((_OP_D, j))
            j = 0
            break

        if j == 0:
            # Only query left — emit insertions for remaining i.
            ops.append((_OP_I, i))
            i = 0
            break

        if cur_state == 'H':
            ops.append((_OP_M, 1))
            src = tbH[i][j]
            i -= 1
            j -= 1
            cur_state = 'H' if src == _TBH_H else ('D' if src == _TBH_D else 'I')

        elif cur_state == 'D':
            ops.append((_OP_D, 1))
            src = tbD[i][j]
            j -= 1
            cur_state = 'H' if src == _TBX_OPEN else 'D'

        else:  # 'I'
            ops.append((_OP_I, 1))
            src = tbI_[i][j]
            i -= 1
            cur_state = 'H' if src == _TBX_OPEN else 'I'

    ref_consumed = j_best
    ops.reverse()
    return _compress(ops), ref_consumed


def score_left_anchored(query: str, ref: str) -> Tuple[float, int]:
    """Return the best affine-gap score and ref_consumed for a left-anchored
    semi-global alignment of *query* against *ref*.

    Runs the same DP as :func:`_align_left_anchored` but returns only the
    score and the number of reference bases consumed at the optimal end column
    (the soft-clip boundary), without tracing back the full CIGAR.

    This is used by the junction refiner's tier-1 scoring step: for each
    candidate 5'SS ``je``, align ``rescue[k:]`` against ``g[je:je+buffer]``
    to find the best (score, ref_consumed) pair.  ``ref_consumed`` tells us
    how many exon-1 bases were confidently matched before noise begins —
    bases beyond that are soft-clipped.

    Args:
        query: Query sequence (rescue suffix in BAM orientation).
        ref:   Reference sequence starting at intron_end / the first exon-1 base.

    Returns:
        ``(best_score, ref_consumed)`` where *best_score* is the Gotoh affine
        score at the optimal free-suffix end column and *ref_consumed* is the
        number of ref bases consumed (= soft-clip boundary).
    """
    Q, R = len(query), len(ref)
    if Q == 0:
        return 0.0, 0
    if R == 0:
        # All query bases become insertions; score = gap_open + Q * gap_extend
        return _GAP_OPEN + Q * _GAP_EXTEND, 0

    H  = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]
    D  = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]
    I_ = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]

    H[0][0] = 0.0
    for j in range(1, R + 1):
        D[0][j] = _GAP_OPEN + j * _GAP_EXTEND
    for i in range(1, Q + 1):
        I_[i][0] = _GAP_OPEN + i * _GAP_EXTEND

    for i in range(1, Q + 1):
        qi = query[i - 1].upper()
        for j in range(1, R + 1):
            s = _MATCH if qi == ref[j - 1].upper() else _MISMATCH

            best_h = max(H[i-1][j-1], D[i-1][j-1], I_[i-1][j-1])
            H[i][j] = best_h + s if best_h > _NEG_INF else _NEG_INF

            d_open = H[i][j-1] + _GAP_OPEN + _GAP_EXTEND if H[i][j-1] > _NEG_INF else _NEG_INF
            d_ext  = D[i][j-1] + _GAP_EXTEND if D[i][j-1] > _NEG_INF else _NEG_INF
            D[i][j] = max(d_open, d_ext)

            i_open = H[i-1][j] + _GAP_OPEN + _GAP_EXTEND if H[i-1][j] > _NEG_INF else _NEG_INF
            i_ext  = I_[i-1][j] + _GAP_EXTEND if I_[i-1][j] > _NEG_INF else _NEG_INF
            I_[i][j] = max(i_open, i_ext)

    # Free-suffix: best end column across all j in [0, R]
    j_best = max(range(R + 1), key=lambda j: max(H[Q][j], D[Q][j], I_[Q][j]))
    best_score = max(H[Q][j_best], D[Q][j_best], I_[Q][j_best])

    return best_score, j_best


def affine_score_to_edit_distance(score: float, query_len: int) -> float:
    """Convert an affine-gap alignment score to an approximate edit distance.

    The affine score uses match=+2, mismatch=-4, gap_open=-4, gap_extend=-1.
    A perfect alignment of *query_len* bases scores ``query_len * MATCH``.

    We normalize to ``[0, query_len]`` where 0 = perfect match and
    *query_len* = completely unrelated sequences:

        edit_dist = (perfect_score - score) / (MATCH - MISMATCH)
                  = (query_len * 2 - score) / 6

    This is comparable to the hp_edit_distance values used in tier-2 gap
    scoring, allowing the two scores to be added meaningfully.

    Args:
        score:      Output of :func:`score_left_anchored`.
        query_len:  Number of query bases used in the alignment.

    Returns:
        Float edit distance in [0, query_len].
    """
    if query_len == 0:
        return 0.0
    perfect = query_len * _MATCH
    # Clamp: score can't exceed perfect, but allow slightly better (no penalty)
    return max(0.0, (perfect - score) / (_MATCH - _MISMATCH))


# ---------------------------------------------------------------------------
# Homopolymer-aware global alignment (exon block realignment)
# ---------------------------------------------------------------------------

def _is_homopolymer_ref(seq: str, pos: int, min_run: int = 3) -> bool:
    """Return True if position *pos* in *seq* is within a homopolymer run >= *min_run*."""
    if pos < 0 or pos >= len(seq):
        return False
    base = seq[pos].upper()
    if base == 'N':
        return False
    left = pos
    while left > 0 and seq[left - 1].upper() == base:
        left -= 1
    right = pos + 1
    while right < len(seq) and seq[right].upper() == base:
        right += 1
    return (right - left) >= min_run


def align_exon_block_global(
    query: str,
    ref: str,
    chrom_ref: str = '',
    ref_offset: int = 0,
    homo_mismatch: float = -2.0,
    min_run: int = 3,
) -> List[Tuple[int, int]]:
    """
    Global (Needleman-Wunsch) affine-gap alignment with homopolymer-aware scoring.

    Both query and ref are fully consumed, so the result can safely replace exon-block
    CIGAR ops without changing the query or reference span:

      sum of lengths where op in {M=0, I=1} == len(query)
      sum of lengths where op in {M=0, D=2} == len(ref)

    Homopolymer-aware scoring: at ref positions within a homopolymer run >= *min_run*,
    mismatches receive *homo_mismatch* (default -2) instead of the standard -4 penalty.
    This allows the aligner to prefer an indel over a mismatch at homopolymer positions,
    correcting the nanopore DRS systematic length undercalling artifact.

    Scoring (Gotoh 1982 affine gap):
      match           = +2
      mismatch        = -4  (or homo_mismatch at homopolymer positions)
      gap_open        = -4  (paid once per gap event)
      gap_extend      = -1  (paid per base in the gap)

    Args:
        query:        Query sequence for the exon block.
        ref:          Reference sequence spanning the exon block.
        chrom_ref:    Full chromosome sequence for homopolymer detection.
                      If empty, standard mismatch scoring is used throughout.
        ref_offset:   0-based position of ref[0] within chrom_ref.
        homo_mismatch: Mismatch score at homopolymer ref positions (default -2.0).
        min_run:      Minimum homopolymer run length to apply the reduced penalty.

    Returns:
        Compressed list of (op_code, length) CIGAR tuples (M/I/D only).
    """
    Q, R = len(query), len(ref)
    if Q == 0 and R == 0:
        return []
    if Q == 0:
        return [(_OP_D, R)]
    if R == 0:
        return [(_OP_I, Q)]

    # Three Gotoh DP matrices
    H  = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]
    D  = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]
    I_ = [[_NEG_INF] * (R + 1) for _ in range(Q + 1)]

    tbH  = [[0] * (R + 1) for _ in range(Q + 1)]
    tbD  = [[0] * (R + 1) for _ in range(Q + 1)]
    tbI_ = [[0] * (R + 1) for _ in range(Q + 1)]

    # Standard global (NW) boundary: both ends anchored
    H[0][0] = 0.0
    for j in range(1, R + 1):
        D[0][j] = _GAP_OPEN + j * _GAP_EXTEND
        tbD[0][j] = _TBX_OPEN if j == 1 else _TBX_EXTEND
    for i in range(1, Q + 1):
        I_[i][0] = _GAP_OPEN + i * _GAP_EXTEND
        tbI_[i][0] = _TBX_OPEN if i == 1 else _TBX_EXTEND

    # Precompute homopolymer mask for ref positions
    if chrom_ref:
        homo_mask = [
            _is_homopolymer_ref(chrom_ref, ref_offset + j, min_run)
            for j in range(R)
        ]
    else:
        homo_mask = [False] * R

    # Main DP
    for i in range(1, Q + 1):
        qi = query[i - 1].upper()
        for j in range(1, R + 1):
            if qi == ref[j - 1].upper():
                s = _MATCH
            elif homo_mask[j - 1]:
                s = homo_mismatch
            else:
                s = _MISMATCH

            # H: match / mismatch (diagonal)
            h_h = H[i-1][j-1]
            h_d = D[i-1][j-1]
            h_i = I_[i-1][j-1]
            best_prev = max(h_h, h_d, h_i)
            if best_prev == _NEG_INF:
                H[i][j] = _NEG_INF
            else:
                H[i][j] = best_prev + s
                if best_prev == h_h:
                    tbH[i][j] = _TBH_H
                elif best_prev == h_d:
                    tbH[i][j] = _TBH_D
                else:
                    tbH[i][j] = _TBH_I

            # D: deletion (gap in query, reference consumed)
            h_prev = H[i][j-1]
            d_prev = D[i][j-1]
            d_open   = (h_prev + _GAP_OPEN + _GAP_EXTEND) if h_prev != _NEG_INF else _NEG_INF
            d_extend = (d_prev + _GAP_EXTEND)              if d_prev != _NEG_INF else _NEG_INF
            if d_open >= d_extend:
                D[i][j] = d_open
                tbD[i][j] = _TBX_OPEN
            else:
                D[i][j] = d_extend
                tbD[i][j] = _TBX_EXTEND

            # I: insertion (gap in reference, query consumed)
            h_prev2 = H[i-1][j]
            i_prev  = I_[i-1][j]
            i_open   = (h_prev2 + _GAP_OPEN + _GAP_EXTEND) if h_prev2 != _NEG_INF else _NEG_INF
            i_extend = (i_prev  + _GAP_EXTEND)              if i_prev  != _NEG_INF else _NEG_INF
            if i_open >= i_extend:
                I_[i][j] = i_open
                tbI_[i][j] = _TBX_OPEN
            else:
                I_[i][j] = i_extend
                tbI_[i][j] = _TBX_EXTEND

    # Endpoint at (Q, R) — global: must finish here
    end_h, end_d, end_i = H[Q][R], D[Q][R], I_[Q][R]
    best_end = max(end_h, end_d, end_i)
    if end_h == best_end:
        cur_state = 'H'
    elif end_d == best_end:
        cur_state = 'D'
    else:
        cur_state = 'I'

    # Traceback from (Q, R) to (0, 0)
    ops: List[Tuple[int, int]] = []
    i, j = Q, R

    while i > 0 or j > 0:
        if i == 0:
            ops.append((_OP_D, j))
            break
        if j == 0:
            ops.append((_OP_I, i))
            break

        if cur_state == 'H':
            ops.append((_OP_M, 1))
            src = tbH[i][j]
            i -= 1
            j -= 1
            cur_state = 'H' if src == _TBH_H else ('D' if src == _TBH_D else 'I')
        elif cur_state == 'D':
            ops.append((_OP_D, 1))
            src = tbD[i][j]
            j -= 1
            cur_state = 'H' if src == _TBX_OPEN else 'D'
        else:  # 'I'
            ops.append((_OP_I, 1))
            src = tbI_[i][j]
            i -= 1
            cur_state = 'H' if src == _TBX_OPEN else 'I'

    ops.reverse()
    return _compress(ops)


# ---------------------------------------------------------------------------
# CIGAR string utilities
# ---------------------------------------------------------------------------

_SAM_OP_CODES = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
_SAM_OP_CHARS = {v: k for k, v in _SAM_OP_CODES.items()}

_CIGAR_TOKEN_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def cigar_ops_to_str(ops: List[Tuple[int, int]]) -> str:
    """Convert a list of ``(op_code, length)`` tuples to a SAM CIGAR string."""
    return ''.join(f'{length}{_SAM_OP_CHARS.get(op, "?")}' for op, length in ops)


def cigar_str_to_ops(cigar_str: str) -> List[Tuple[int, int]]:
    """Parse a SAM CIGAR string into a list of ``(op_code, length)`` tuples."""
    return [(_SAM_OP_CODES[op], int(length))
            for length, op in _CIGAR_TOKEN_RE.findall(cigar_str)]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def align_clip_to_exon(
    clip_seq: str,
    genome_seq: str,
    intron_start: int,
    intron_end: int,
    strand: str,
    max_indel: int = 5,
) -> Tuple[List[Tuple[int, int]], int]:
    """
    Align a 5' soft-clip to the upstream exon using semi-global affine-gap NW.

    The alignment is anchored at the exon–intron boundary:

    * **Plus strand** (``'+'``): clip must end at ``intron_start``
      (right-anchored in exon).  Reference region used:
      ``genome_seq[intron_start - clip_len - max_indel : intron_start]``.

    * **Minus strand** (``'-'``): clip must start at ``intron_end``
      (left-anchored in exon).  Reference region used:
      ``genome_seq[intron_end : intron_end + clip_len + max_indel]``.

    The *clip_seq* is in BAM orientation — no reverse-complement is needed
    because BAM already stores the reverse-complemented minus-strand sequence,
    which cancels out with the genome complement.

    Args:
        clip_seq:     Soft-clip bases (BAM orientation).
        genome_seq:   Chromosome sequence as a plain Python string.  Must cover
                      at least the alignment region (see above).
        intron_start: 0-based intron start (first intron base; exclusive end of
                      exon 1 for ``'+'`` strand).
        intron_end:   0-based intron end (exclusive; first base of the
                      downstream exon for ``'+'`` strand).
        strand:       ``'+'`` or ``'-'``.
        max_indel:    Buffer (bp) added to each side of the expected exon
                      region to accommodate insertions/deletions.

    Returns:
        ``(cigar_ops, exon_ref_start)`` where *cigar_ops* is the list of
        ``(op_code, length)`` tuples for the exon segment (M/I/D only) and
        *exon_ref_start* is the 0-based reference coordinate where the exon
        alignment begins.
    """
    clip_len = len(clip_seq)
    if clip_len == 0:
        return [], intron_start if strand == '+' else intron_end

    if strand == '+':
        region_start = max(0, intron_start - clip_len - max_indel)
        ref_region = genome_seq[region_start:intron_start]
        if not ref_region:
            return [(_OP_M, clip_len)], max(0, intron_start - clip_len)

        cigar_ops, ref_skip = _align_right_anchored(clip_seq, ref_region)
        exon_ref_start = region_start + ref_skip
        return cigar_ops, exon_ref_start

    else:  # minus strand
        region_end = intron_end + clip_len + max_indel
        ref_region = genome_seq[intron_end:region_end]
        if not ref_region:
            return [(_OP_M, clip_len)], intron_end

        cigar_ops, _ref_consumed = _align_left_anchored(clip_seq, ref_region)
        return cigar_ops, intron_end
