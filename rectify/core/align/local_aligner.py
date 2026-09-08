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
from typing import List, NamedTuple, Tuple
import logging

logger = logging.getLogger(__name__)

# Scoring constants
_MATCH      =  2
_MISMATCH   = -4
_GAP_OPEN   = -4   # penalty for opening a new gap (paid once per gap event)
_GAP_EXTEND = -1   # penalty per base within a gap

_NEG_INF = float('-inf')

# Buffer (bp) added to the far side of the expected exon window by
# align_clip_to_exon — and, since ISSUE-020, by the 2F ranking scorers, which
# must see the SAME reference window as the placement.
ANCHOR_MAX_INDEL = 5

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


def score_right_anchored(query: str, ref: str) -> Tuple[float, str]:
    """Best affine-gap score of a RIGHT-anchored semi-global alignment of
    *query* against *ref* — the score :func:`_align_right_anchored` traces back,
    without the traceback.

    Query fully consumed; free prefix in *ref*; the alignment must end at
    ``ref[-1]`` (the last exon-1 base before the donor on the plus strand).
    Same recurrence, same boundary conditions, same four constants as
    ``_align_right_anchored`` — so for a given (query, ref) this score is
    exactly the score of the CIGAR ``align_clip_to_exon`` emits (the ISSUE-020
    consistency invariant; ``tests/test_2f_anchored_ranking.py``). Rolling rows:
    O(R) memory, no traceback matrices.

    Returns ``(best_score, end_state)`` with *end_state* in ``'H'/'D'/'I'`` —
    the state at the anchored end, ``'D'`` meaning the alignment closes with a
    deletion at the junction (a junction-side gap). Ties are broken H, D, I,
    as the traceback does.
    """
    Q, R = len(query), len(ref)
    if Q == 0:
        return 0.0, 'H'
    if R == 0:
        return _GAP_OPEN + Q * _GAP_EXTEND, 'I'
    qu = query.upper()
    ru = ref.upper()
    go = _GAP_OPEN + _GAP_EXTEND
    # Row 0: free prefix — the alignment may start at any ref column for free,
    # in the match state or "already skipping" in the deletion state.
    h_prev = [0.0] * (R + 1)
    d_prev = [0.0] * (R + 1)
    i_prev = [_NEG_INF] * (R + 1)
    for i in range(1, Q + 1):
        qi = qu[i - 1]
        h_cur = [_NEG_INF] * (R + 1)
        d_cur = [_NEG_INF] * (R + 1)
        i_cur = [_NEG_INF] * (R + 1)
        i_cur[0] = _GAP_OPEN + i * _GAP_EXTEND
        for j in range(1, R + 1):
            s = _MATCH if qi == ru[j - 1] else _MISMATCH
            best_h = max(h_prev[j - 1], d_prev[j - 1], i_prev[j - 1])
            if best_h > _NEG_INF:
                h_cur[j] = best_h + s
            hl = h_cur[j - 1]
            dl = d_cur[j - 1]
            d_open = hl + go if hl > _NEG_INF else _NEG_INF
            d_ext = dl + _GAP_EXTEND if dl > _NEG_INF else _NEG_INF
            d_cur[j] = d_open if d_open >= d_ext else d_ext
            hu = h_prev[j]
            iu = i_prev[j]
            i_open = hu + go if hu > _NEG_INF else _NEG_INF
            i_ext = iu + _GAP_EXTEND if iu > _NEG_INF else _NEG_INF
            i_cur[j] = i_open if i_open >= i_ext else i_ext
        h_prev, d_prev, i_prev = h_cur, d_cur, i_cur
    end_h, end_d, end_i = h_prev[R], d_prev[R], i_prev[R]
    best = max(end_h, end_d, end_i)
    state = 'H' if end_h == best else ('D' if end_d == best else 'I')
    return best, state


def affine_cigar_score(ops: List[Tuple[int, int]], query: str, ref: str) -> float:
    """Score an M/I/D CIGAR over *query* and *ref* (both fully consumed) with
    this module's four constants — the quantity the anchored DPs maximize.

    Each M/=/X base scores ``_MATCH`` or ``_MISMATCH``; each I or D run costs
    ``_GAP_OPEN + len * _GAP_EXTEND`` (one gap event per run — the DPs never
    place two gap runs of the same kind back to back, and ``_compress`` merges
    adjacent ops, so a run IS an event). Used by the ISSUE-020 consistency
    check: the emitted exon CIGAR scored this way must equal the ranking score.

    Raises ``ValueError`` when the ops do not consume exactly *query* and *ref*.
    """
    qi = ri = 0
    score = 0.0
    for op, ln in ops:
        if op in (_OP_M, 7, 8):
            for k in range(ln):
                score += _MATCH if query[qi + k].upper() == ref[ri + k].upper() else _MISMATCH
            qi += ln
            ri += ln
        elif op == _OP_I:
            score += _GAP_OPEN + ln * _GAP_EXTEND
            qi += ln
        elif op == _OP_D:
            score += _GAP_OPEN + ln * _GAP_EXTEND
            ri += ln
        else:
            raise ValueError(f"affine_cigar_score: unsupported op {op}")
    if qi != len(query) or ri != len(ref):
        raise ValueError(
            f"affine_cigar_score: ops consume {qi} query / {ri} ref bases, "
            f"sequences have {len(query)} / {len(ref)}")
    return score


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


def _homopolymer_run_len(seq: str, pos: int) -> Tuple[int, str]:
    """Length of the homopolymer run covering position *pos*, and its base.

    Returns ``(run_len, base)``; ``(0, '')`` if *pos* is out of range or the base
    is N. Used by the C1 length-law gap cost to look up the empirical per-(run
    length, base) deletion/insertion penalty. Shares the scan logic with
    ``_is_homopolymer_ref`` but returns the LENGTH (the load-bearing quantity for
    the penalty table) rather than a >= min_run boolean."""
    if pos < 0 or pos >= len(seq):
        return 0, ""
    base = seq[pos].upper()
    if base == 'N':
        return 0, ""
    left = pos
    while left > 0 and seq[left - 1].upper() == base:
        left -= 1
    right = pos + 1
    while right < len(seq) and seq[right].upper() == base:
        right += 1
    return right - left, base


def align_exon_block_global(
    query: str,
    ref: str,
    chrom_ref: str = '',
    ref_offset: int = 0,
    homo_mismatch: float = -2.0,
    min_run: int = 3,
    penalty_table=None,
    lam: float = 1.0,
    ins_lengthlaw: bool = False,
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

    # C1 length-law gap-OPEN deltas (added to gap_open at HP-mask ref positions
    # only; penalty_table is any object exposing del_open_delta(hp, base, lam) /
    # ins_open_delta(...). penalty_table=None => all-zero => byte-identical to the
    # legacy DP, the Cat3 / junction-rescue regression guard). The delta is
    # POSITIVE in longer runs (rate rises with hp) => the gap-open is less
    # negative => the DP prefers the in-run deletion over out-of-run misplacement.
    #
    # INSERTION discount is GATED OFF by default (ins_lengthlaw=False): the real-SIRV
    # over-call ablation (c1_real_sirv_ablation.py) showed it HALLUCINATES indels —
    # a cheap insertion lets the DP rewrite a length-preserving substitution as a
    # spurious D+I (3-7% over-call on sub-only windows, vs 0% for the deletion-only
    # law, replicated on LRGASP + SG-NEx real SIRV). The deletion law is safe (0%).
    # Re-enable ins only once it is independently validated (the injection-simulator
    # Claim B; see dev/C1_DESIGN.md). The generator injects deletions only anyway.
    del_open_d = [0.0] * R
    ins_open_d = [0.0] * R
    if penalty_table is not None and chrom_ref:
        rseq = chrom_ref
        for j in range(R):
            if not homo_mask[j]:
                continue
            run_len, base = _homopolymer_run_len(rseq, ref_offset + j)
            if run_len <= 0:
                continue
            del_open_d[j] = penalty_table.del_open_delta(run_len, base, lam)
            if ins_lengthlaw:
                ins_open_d[j] = penalty_table.ins_open_delta(run_len, base, lam)

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

            # D: deletion (gap in query, reference consumed). The length-law delta
            # adjusts the OPEN cost only (paid once per gap), gated on ref[j-1]'s
            # HP context; extend stays flat (continuation is a separate axis).
            h_prev = H[i][j-1]
            d_prev = D[i][j-1]
            d_open   = (h_prev + _GAP_OPEN + _GAP_EXTEND + del_open_d[j-1]) if h_prev != _NEG_INF else _NEG_INF
            d_extend = (d_prev + _GAP_EXTEND)              if d_prev != _NEG_INF else _NEG_INF
            if d_open >= d_extend:
                D[i][j] = d_open
                tbD[i][j] = _TBX_OPEN
            else:
                D[i][j] = d_extend
                tbD[i][j] = _TBX_EXTEND

            # I: insertion (gap in reference, query consumed). Length-law delta on
            # OPEN only, using the HP context of the preceding ref base (j-1).
            h_prev2 = H[i-1][j]
            i_prev  = I_[i-1][j]
            i_delta = ins_open_d[j-1] if j >= 1 else 0.0
            i_open   = (h_prev2 + _GAP_OPEN + _GAP_EXTEND + i_delta) if h_prev2 != _NEG_INF else _NEG_INF
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
# Evidence shape of a placed exon block (ISSUE-028 / invariant E, 2026-09-06)
# ---------------------------------------------------------------------------

_OP_S = 4  # soft clip (query consumed, no reference)

#: Genome homopolymer run length at (or above) which a mismatch or an indel
#: inside the run is charged HALF (DRS basecallers slip inside long runs). The
#: bundled empirical error table is iteration-4 material; this is the interim.
EVIDENCE_HP_FORGIVE_RUN = 5

# The evidence SCORE of a placed block, in bits (Kevin's ruling, 2026-09-06):
# a matched base against a four-letter alphabet is log2(4) = 2 bits of
# evidence; a mismatch and an affine gap cost the anchored aligner's own
# constants scaled by 1/2 (_MISMATCH -4 -> -2 bits, _GAP_OPEN -4 -> -2 bits per
# gap event, _GAP_EXTEND -1 -> -0.5 bits per gap base), so the costs keep the
# placement model's relative weights while the match keeps its information
# value (scaling the match too would read a matched base as 1 bit and no
# 14-nt anchor could reach the cutoffs the model is run at). Inside a genome
# homopolymer run >= EVIDENCE_HP_FORGIVE_RUN a mismatch or a gap is charged
# half (the -log2 of the bundled DRS error rates replaces the halving when the
# error table lands). Worked examples: `6=1X7=` = 26 - 2 = 24 bits; 22f609c6's
# `1I15M` (13=/2X) = 26 - 4 - 2.5 = 19.5 bits; `5M4D2M` = 14 - 4 = 10 bits;
# `9=/23X` = 18 - 46 = -28 bits; `11=/15X` = 22 - 30 = -8 bits.
_BITS_MATCH = 2.0                 # log2(4)
_BITS_MISMATCH = _MISMATCH / 2.0
_BITS_GAP_OPEN = _GAP_OPEN / 2.0
_BITS_GAP_EXTEND = _GAP_EXTEND / 2.0


class EvidenceShape(NamedTuple):
    """What the placed 5' block IS, walked base by base against the genome.

    ``matched`` / ``mismatches`` are raw counts of ``=`` / ``X`` columns;
    ``identity`` is ``matched / (matched + mismatches_eff)`` with a mismatch
    inside a genome homopolymer run >= ``EVIDENCE_HP_FORGIVE_RUN`` counting
    half; ``bits`` is the evidence score of the block (see the module
    comment above ``_BITS_MATCH``); ``junction_clean_run`` is the number of
    consecutive matched REFERENCE bases abutting the junction (a mismatch or a
    deletion breaks it; an insertion consumes no reference and is passed over,
    the same way the base-by-base review view does not show it) — reported,
    not gated on (a hard clean-run floor was the wrong primitive: `6=1X7=` is
    an excellent 14-nt anchor); ``leading_indel_len`` is the number of I / D /
    S bases at the 5' (free) end before the first aligned base — the S ops
    :func:`strip_leading_indel` left there, when it ran.
    """
    matched: int
    mismatches: int
    identity: float
    bits: float
    junction_clean_run: int
    leading_indel_len: int


def strip_leading_indel(ops: List[Tuple[int, int]], strand: str) -> Tuple[List[Tuple[int, int]], int]:
    """A LEADING insertion at the 5' (free) end of an exon CIGAR is a soft clip
    in disguise (``8I6M`` = six placed bases, not fourteen); a leading deletion
    consumes no query and is dropped. Returns ``(ops, unplaced_query_bases)``
    where the leading I bases are re-emitted as one ``S`` op at the 5' end (the
    writer prepends / appends the exon ops verbatim and counts S as
    query-consuming in both 5' surgeries, so the bases stay soft-clipped in the
    BAM). The 5' end is the LEFT end for the plus strand and the RIGHT end for
    the minus strand (BAM orientation). ``[]`` when nothing aligned survives.
    """
    if not ops:
        return list(ops or []), 0
    ops = list(ops)
    if strand != '+':
        ops.reverse()
    unplaced = 0
    while ops and ops[0][0] in (_OP_I, _OP_D, _OP_S):
        op, ln = ops.pop(0)
        if op in (_OP_I, _OP_S):
            unplaced += ln
    if not ops:
        return [], unplaced
    if unplaced:
        ops.insert(0, (_OP_S, unplaced))
    if strand != '+':
        ops.reverse()
    return ops, unplaced


def evidence_shape(
    cigar_ops: List[Tuple[int, int]],
    clip_seq: str,
    genome_seq: str,
    intron_start: int,
    intron_end: int,
    strand: str,
    hp_forgive_run: int = EVIDENCE_HP_FORGIVE_RUN,
) -> EvidenceShape:
    """Walk *cigar_ops* (the placement :func:`align_clip_to_exon` describes,
    anchored at the junction: right end at ``intron_start`` for ``+``, left end
    at ``intron_end`` for ``-``) over *clip_seq* and the genome and return the
    :class:`EvidenceShape`. S / I ops consume query only, D consumes reference
    only. Sequence comparison is case-insensitive; a query or reference base
    that is off the end of its sequence counts as a mismatch."""
    ops = list(cigar_ops or [])
    ref_len = sum(ln for op, ln in ops if op in (_OP_M, _OP_D, 7, 8))
    if strand == '+':
        ref = genome_seq[max(0, intron_start - ref_len):intron_start]
        ref_off = intron_start - len(ref)
    else:
        ref = genome_seq[intron_end:intron_end + ref_len]
        ref_off = intron_end
    ref_u = ref.upper()
    q_u = (clip_seq or '').upper()
    cols: List[str] = []       # one mark per column in BAM (left -> right) order
    matched = mism = 0
    forgiven = 0.0
    bits = 0.0
    qi = ri = 0

    def _in_hp_run(ref_pos: int) -> bool:
        return _homopolymer_run_len(genome_seq, ref_pos)[0] >= hp_forgive_run

    for op, ln in ops:
        if op in (_OP_M, 7, 8):
            for k in range(ln):
                qb = q_u[qi + k] if qi + k < len(q_u) else ''
                rb = ref_u[ri + k] if ri + k < len(ref_u) else ''
                if qb and rb and qb == rb:
                    matched += 1
                    bits += _BITS_MATCH
                    cols.append('=')
                else:
                    mism += 1
                    cols.append('X')
                    if rb and _in_hp_run(ref_off + ri + k):
                        forgiven += 0.5
                        bits += _BITS_MISMATCH / 2.0
                    else:
                        bits += _BITS_MISMATCH
            qi += ln
            ri += ln
        elif op == _OP_S:
            # the leading bases strip_leading_indel left unplaced: not part of the block
            cols.extend('I' * ln)
            qi += ln
        elif op in (_OP_I, _OP_D):
            cols.extend(('I' if op == _OP_I else 'D') * ln)
            # a gap event: charged where it sits in the genome (an insertion at
            # the reference position it interrupts, a deletion over its bases)
            gap_cost = _BITS_GAP_OPEN + ln * _BITS_GAP_EXTEND
            hp_pos = ref_off + ri if op == _OP_I else ref_off + ri
            if _in_hp_run(hp_pos) or (op == _OP_D and _in_hp_run(ref_off + ri + ln - 1)):
                gap_cost /= 2.0
            bits += gap_cost
            if op == _OP_I:
                qi += ln
            else:
                ri += ln
    # Junction side: the RIGHT end for '+', the LEFT end for '-'.
    junction_first = cols[::-1] if strand == '+' else cols
    run = 0
    for c in junction_first:
        if c == '=':
            run += 1
        elif c == 'I':
            continue
        else:
            break
    lead = 0
    for c in reversed(junction_first):
        if c in ('I', 'D'):
            lead += 1
        else:
            break
    aligned_eff = matched + mism - forgiven
    identity = (matched / aligned_eff) if aligned_eff > 0 else 0.0
    return EvidenceShape(matched, mism, identity, bits, run, lead)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def align_clip_to_exon(
    clip_seq: str,
    genome_seq: str,
    intron_start: int,
    intron_end: int,
    strand: str,
    max_indel: int = ANCHOR_MAX_INDEL,
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

    # Guard: O(clip_len²) Python Gotoh DP is prohibitive for long intronic reads.
    # Biological 5' clips at exon boundaries are <100 bp; clips >500 are
    # artifact/pathology (e.g. a read mapped deeply inside an intron).  Return an
    # all-M CIGAR identical to the existing `not ref_region` early-return semantics.
    _MAX_CLIP_ALIGN_LEN = 500
    if clip_len > _MAX_CLIP_ALIGN_LEN:
        if strand == '+':
            return [(_OP_M, clip_len)], max(0, intron_start - clip_len)
        else:
            return [(_OP_M, clip_len)], intron_end

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
