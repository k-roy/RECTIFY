"""Information-bounded splice-overhang search: the shared gate.

One module, three call sites (design: planning/641 SPEC §2 + the 2026-08-09
scope correction): the mapPacBio-role overhang resolver
(``core/align/overhang_resolver.py``), the deSALT-role reimplementation
(planning/638, future), and ``rescue_3ss_truncation``
(``core/splice/splice_aware_5prime.py``).

Principle: **never search a window larger than the overhang can distinguish
itself within.** For an overhang of effective information content ``I_eff``
bits, the expected number of chance matches in a window of width ``W`` is
``E = W * 2**(-I_eff)``. Requiring ``E <= alpha`` gives

    W_max = alpha * 2**I_eff

clamped to ``[0, max_window]``. ``W_max < 1`` means the overhang cannot
distinguish its true placement from chance anywhere — the correct outcome is
to REFUSE the sequence search (emit the soft clip unchanged), which is a
first-class result, not a failure. A poly(A) tail is the canonical refused
overhang: searching a splice acceptor for ``AAAAAAAAAAAATAAA`` is searching
for an intron in the poly(A) tail.

``I_eff`` is the min of two cheap O(L) estimators (both in bits):

- ``L_eff * H_comp``: base-composition entropy with homopolymer runs
  discounted — a run of length r contributes effective length
  ``1 + hp_discount*(r-1)`` (default 0.5, mirroring ``_hp_edit_distance``'s
  HP indel cost: ONT miscalls run length, so length inside a run carries
  less information than composition suggests).
- ``L * H_cond``: order-1 conditional plug-in entropy
  ``H(dinucleotides) - H(bases)``. This catches PERIODIC low complexity
  ((AG)n, (AAG)n ...) that composition entropy misses — the graded version
  of what ``repeat_expansion.py`` detects binarily.

Worked example (matches the 641 SPEC table): ``AAAAAAAAAAAATAAA`` scores
~4.6 bits => W_max ~ 0.24 bp at alpha=0.01 => refused. A random 16-mer
scores ~30 bits => W_max clamps to max_window. ``(AG)8`` scores ~0 via the
conditional term => refused.

Counters: module-level ``COUNTERS`` (per-process; workers under
``ProcessPoolExecutor`` each hold their own copy). The 641 acceptance test
T1 asserts on the counter, not the runtime.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

DEFAULT_ALPHA = 0.01
DEFAULT_HP_DISCOUNT = 0.5
# Exponent cap before 2**I_eff: beyond this W_max exceeds any genome and the
# clamp decides anyway; capping avoids float overflow on long complex clips.
_MAX_I_EFF_EXP = 64.0

_ACGT = frozenset('ACGT')

# Per-process instrumentation. Cheap ints, always on; reset in tests.
COUNTERS: Dict[str, int] = {
    'assessed': 0,          # overhangs assessed
    'refused': 0,           # W_max < 1 => sequence search refused
    'window_bounded': 0,    # W_max < the caller's default window (search shrunk)
    'candidates_evaluated': 0,   # candidate placements actually scored (resolver)
    'candidates_skipped': 0,     # candidates excluded by the W_max bound (rescue)
}


def reset_counters() -> None:
    for k in COUNTERS:
        COUNTERS[k] = 0


def _maybe_register_counter_dump() -> None:
    """When RECTIFY_OVERHANG_INFO_COUNTS names a file, append this process's
    counters as one JSON line at exit. Best-effort instrumentation for A/B
    gates (planning/644 T8): reliable in the main process; under spawn-start
    workers each re-import registers its own handler, but worker atexit is
    not guaranteed (planning/596) — run single-threaded when the counts must
    be complete."""
    path = os.environ.get('RECTIFY_OVERHANG_INFO_COUNTS')
    if not path:
        return
    import atexit
    import json

    def _dump() -> None:
        try:
            with open(path, 'a') as fh:
                fh.write(json.dumps({'pid': os.getpid(), **COUNTERS}) + '\n')
        except OSError:
            pass

    atexit.register(_dump)


_maybe_register_counter_dump()


def gate_enabled() -> bool:
    """Rescue-path wiring is DARK by default (policy change, not an exact
    optimization — see planning/596's T4 exclusion). Enable with
    ``RECTIFY_OVERHANG_INFO_GATE=1`` after the T8 identity/cost gate."""
    return os.environ.get('RECTIFY_OVERHANG_INFO_GATE', '0') == '1'


def gate_alpha() -> float:
    try:
        return float(os.environ.get('RECTIFY_OVERHANG_INFO_ALPHA', str(DEFAULT_ALPHA)))
    except ValueError:
        return DEFAULT_ALPHA


def _entropy_bits(counts: Dict[str, float], total: float) -> float:
    if total <= 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            h -= p * math.log2(p)
    return h


def effective_information_bits(seq: str, hp_discount: float = DEFAULT_HP_DISCOUNT) -> float:
    """I_eff of an overhang, in bits. Non-ACGT bases contribute zero."""
    s = [b for b in seq.upper() if b in _ACGT]
    n = len(s)
    if n < 2:
        return 0.0

    # --- Estimator 1: HP-discounted composition entropy ------------------
    # Weight each base by its effective (run-discounted) contribution, so a
    # 12-A run adds 6.5 effective A's, not 12.
    eff_counts: Dict[str, float] = {}
    l_eff = 0.0
    i = 0
    while i < n:
        j = i
        while j < n and s[j] == s[i]:
            j += 1
        run_eff = 1.0 + hp_discount * (j - i - 1)
        eff_counts[s[i]] = eff_counts.get(s[i], 0.0) + run_eff
        l_eff += run_eff
        i = j
    h_comp = _entropy_bits(eff_counts, l_eff)
    i_comp = l_eff * h_comp

    # --- Estimator 2: order-1 conditional entropy ------------------------
    base_counts: Dict[str, float] = {}
    for b in s:
        base_counts[b] = base_counts.get(b, 0.0) + 1.0
    di_counts: Dict[str, float] = {}
    for k in range(n - 1):
        d = s[k] + s[k + 1]
        di_counts[d] = di_counts.get(d, 0.0) + 1.0
    h1 = _entropy_bits(base_counts, float(n))
    h2 = _entropy_bits(di_counts, float(n - 1))
    h_cond = max(0.0, h2 - h1)
    i_cond = n * h_cond

    # Plug-in H2 is biased low at small n (only n-1 dinucleotide samples), so
    # i_cond under-credits very short overhangs. That is the conservative
    # direction: a <=6 nt overhang genuinely cannot distinguish its placement,
    # and refusing it is the correct outcome, so the bias is kept, not
    # corrected.
    return min(i_comp, i_cond)


def max_search_window_bp(
    seq: str,
    alpha: float = DEFAULT_ALPHA,
    max_window: Optional[int] = None,
    hp_discount: float = DEFAULT_HP_DISCOUNT,
) -> int:
    """W_max = alpha * 2**I_eff, clamped to [0, max_window]. 0 => refuse."""
    i_eff = effective_information_bits(seq, hp_discount=hp_discount)
    w = alpha * (2.0 ** min(i_eff, _MAX_I_EFF_EXP))
    if w < 1.0:
        return 0
    if max_window is not None and w > max_window:
        return int(max_window)
    return int(w)


@dataclass(frozen=True)
class OverhangAssessment:
    length: int            # ACGT length actually assessed
    i_eff_bits: float
    w_max_bp: int          # clamped; 0 == refused
    refused: bool
    alpha: float


def assess_overhang(
    seq: str,
    alpha: float = DEFAULT_ALPHA,
    max_window: Optional[int] = None,
    hp_discount: float = DEFAULT_HP_DISCOUNT,
) -> OverhangAssessment:
    """Assess an overhang; update COUNTERS. The one entry point call sites use."""
    s = [b for b in seq.upper() if b in _ACGT]
    i_eff = effective_information_bits(seq, hp_discount=hp_discount)
    w = alpha * (2.0 ** min(i_eff, _MAX_I_EFF_EXP))
    if w < 1.0:
        w_bp = 0
    elif max_window is not None and w > max_window:
        w_bp = int(max_window)
    else:
        w_bp = int(w)
    refused = w_bp == 0
    COUNTERS['assessed'] += 1
    if refused:
        COUNTERS['refused'] += 1
    elif max_window is not None and w_bp < max_window:
        COUNTERS['window_bounded'] += 1
    return OverhangAssessment(
        length=len(s), i_eff_bits=i_eff, w_max_bp=w_bp, refused=refused, alpha=alpha,
    )


# ---------------------------------------------------------------------------
# Junction ambiguity canonicalisation (planning/621).
#
# For an intron [s, e) (0-based half-open) the spliced product is invariant
# under shift +1 iff genome[s] == genome[e], and under -1 iff
# genome[s-1] == genome[e-1]. Canonical form = the LEFTMOST member of the
# ambiguity class, so (chrom, s - l_amb, e - l_amb) is a well-defined
# junction identity key (length is invariant under sliding).
# ---------------------------------------------------------------------------

def ambiguity_window(chrom_seq: str, start: int, end: int, max_shift: int = 60) -> Tuple[int, int]:
    """(l_amb, r_amb): how far the junction can slide left/right as a no-op."""
    g = chrom_seq
    n = len(g)
    r_amb = 0
    while (r_amb < max_shift and end + r_amb < n
           and g[start + r_amb].upper() == g[end + r_amb].upper()):
        r_amb += 1
    l_amb = 0
    while (l_amb < max_shift and start - 1 - l_amb >= 0
           and g[start - 1 - l_amb].upper() == g[end - 1 - l_amb].upper()):
        l_amb += 1
    return l_amb, r_amb


def canonicalize_junction(chrom_seq: str, start: int, end: int, max_shift: int = 60) -> Tuple[int, int]:
    """Leftmost member of the junction's ambiguity class (planning/621)."""
    l_amb, _ = ambiguity_window(chrom_seq, start, end, max_shift=max_shift)
    return start - l_amb, end - l_amb


def same_junction(chrom_seq: str, a: Tuple[int, int], b: Tuple[int, int], max_shift: int = 60) -> bool:
    """True iff two (start, end) introns are the same junction after
    ambiguity canonicalisation. Length must match (invariant under slide)."""
    if (a[1] - a[0]) != (b[1] - b[0]):
        return False
    return canonicalize_junction(chrom_seq, *a, max_shift=max_shift) == \
        canonicalize_junction(chrom_seq, *b, max_shift=max_shift)


# ---------------------------------------------------------------------------
# HP-aware global edit distance with an exact early-exit cutoff.
#
# Same recurrence as splice_aware_5prime._hp_edit_distance (HP indel cost
# 0.5, substitution 1.0, global corner value), implemented here so the
# resolver does not depend on that module's uncommitted planning/596 WIP.
# Cutoff exactness (596 T2 argument): all step costs are non-negative, so
# dp[n][m] >= min_j dp[i][j] for every row i; if the row minimum exceeds the
# cutoff the true value must too. Prune on ``>`` ONLY — a value equal to the
# cutoff must be computed exactly so caller tiebreakers still see it.
# ---------------------------------------------------------------------------

_HP_ED_MAX_LEN = 200


def hp_edit_distance_bounded(s1: str, s2: str, cutoff: float = -1.0) -> float:
    """HP-aware global edit distance; returns > cutoff early when provably
    beaten (cutoff < 0 disables pruning; result then exact always)."""
    n = len(s1)
    m = len(s2)
    if n > _HP_ED_MAX_LEN:
        n = _HP_ED_MAX_LEN
        s1 = s1[:n]
    if m > _HP_ED_MAX_LEN:
        m = _HP_ED_MAX_LEN
        s2 = s2[:m]
    if n == 0:
        return float(m)
    if m == 0:
        return float(n)

    del_cost = [
        0.5 if (i >= 2 and s1[i - 2] == s1[i - 1]) or (i < n and s1[i - 1] == s1[i])
        else 1.0
        for i in range(1, n + 1)
    ]
    ins_cost = [
        0.5 if (j >= 2 and s2[j - 2] == s2[j - 1]) or (j < m and s2[j - 1] == s2[j])
        else 1.0
        for j in range(1, m + 1)
    ]

    prune = cutoff >= 0.0
    prev = [0.0] * (m + 1)
    for j in range(1, m + 1):
        prev[j] = prev[j - 1] + ins_cost[j - 1]
    cur = [0.0] * (m + 1)
    for i in range(1, n + 1):
        dc = del_cost[i - 1]
        c1 = s1[i - 1]
        cur[0] = prev[0] + dc
        row_min = cur[0]
        for j in range(1, m + 1):
            if c1 == s2[j - 1]:
                v = prev[j - 1]
            else:
                v = prev[j - 1] + 1.0
                d = prev[j] + dc
                if d < v:
                    v = d
                ins = cur[j - 1] + ins_cost[j - 1]
                if ins < v:
                    v = ins
            cur[j] = v
            if v < row_min:
                row_min = v
        if prune and row_min > cutoff:
            return cutoff + 1.0
        prev, cur = cur, prev
    return prev[m]
