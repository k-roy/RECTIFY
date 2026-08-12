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
    'period_bounded': 0,    # W_max capped at (self-match period - 1)
    'candidates_evaluated': 0,   # candidate placements actually scored (resolver)
    'candidates_skipped': 0,     # candidates excluded by the W_max bound (rescue)
    'skipped_region': 0,         # reads bypassing rescue via RECTIFY_SKIP_REGIONS
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


# --- Tandem-period self-match (the 2026-08-11 reporter-construct fix) -------
#
# The two I_eff estimators are blind to periodicity with unit length > ~3:
# an MS2 stem-loop array (19-nt unit) scores 76-81 bits and passes the gate,
# yet it matches ITSELF at every 19-nt shift — placements differing by the
# period are indistinguishable, so W_max must not exceed the period. Measured
# consequence of the blind spot: >500x resolver slowdown on MS2 reporter
# constructs (depth x full candidate DP ending in rejected_ambiguous), and
# the same pathology at genomic tandem loci (rDNA precedent; clip legs on
# plain DRS). This term GENERALIZES the existing refusals: poly(A) is the
# period-1 case, (AG)n the period-2 case.

DEFAULT_PERIOD_MAX_SHIFT = 32
DEFAULT_PERIOD_MATCH_FRAC = 0.8   # tolerant of ~10-15% read error in arrays
DEFAULT_PERIOD_MIN_OVERLAP = 12


def min_self_match_period(
    seq: str,
    max_shift: int = DEFAULT_PERIOD_MAX_SHIFT,
    match_frac: float = DEFAULT_PERIOD_MATCH_FRAC,
    min_overlap: int = DEFAULT_PERIOD_MIN_OVERLAP,
) -> Optional[int]:
    """Smallest shift ``s <= max_shift`` at which the overhang matches itself
    over the full overlap at >= ``match_frac`` identity, or None.

    A hit means placements differing by ``s`` are (near-)indistinguishable,
    so the sequence search may not range wider than ``s - 1``. Uppercased
    raw sequence; non-ACGT characters simply fail to match (conservative:
    N-runs are not reported as periodic — their low I_eff already gates them).
    Chance identity is 0.25/base, so 0.8 over >= 12 bases is far off random.
    """
    s = seq.upper()
    n = len(s)
    top = min(max_shift, n - min_overlap)
    for shift in range(1, top + 1):
        ov = n - shift
        need = match_frac * ov
        matches = 0
        for i in range(ov):
            if s[i] == s[i + shift]:
                matches += 1
        if matches >= need:
            return shift
    return None


def max_search_window_bp(
    seq: str,
    alpha: float = DEFAULT_ALPHA,
    max_window: Optional[int] = None,
    hp_discount: float = DEFAULT_HP_DISCOUNT,
) -> int:
    """W_max = min(alpha * 2**I_eff, period - 1), clamped to [0, max_window].
    0 => refuse. ``period`` is the tandem self-match period (None => no cap)."""
    i_eff = effective_information_bits(seq, hp_discount=hp_discount)
    w = alpha * (2.0 ** min(i_eff, _MAX_I_EFF_EXP))
    period = min_self_match_period(seq)
    if period is not None:
        w = min(w, float(period - 1))
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
    period: Optional[int] = None   # tandem self-match period (None = aperiodic)


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
    period = min_self_match_period(seq)
    period_capped = False
    if period is not None and (period - 1) < w:
        w = float(period - 1)
        period_capped = True
    if w < 1.0:
        w_bp = 0
    elif max_window is not None and w > max_window:
        w_bp = int(max_window)
    else:
        w_bp = int(w)
    refused = w_bp == 0
    COUNTERS['assessed'] += 1
    if period_capped:
        COUNTERS['period_bounded'] += 1
    if refused:
        COUNTERS['refused'] += 1
    elif max_window is not None and w_bp < max_window:
        COUNTERS['window_bounded'] += 1
    return OverhangAssessment(
        length=len(s), i_eff_bits=i_eff, w_max_bp=w_bp, refused=refused,
        alpha=alpha, period=period,
    )


# ---------------------------------------------------------------------------
# Pool-admission complexity gate (planning/648 §1 form 2, planning/649).
#
# Same I_eff, different decision form from the search-radius gate above: the
# junction's intron span D is OBSERVED, so instead of deriving a window from
# the information we ask whether the observed placement beats chance —
#
#     E_chance = D * 2**(-I_eff)      admit iff E_chance <= alpha
#
# scored on the WORSE of the two genomic exonic flanks (the read-independent
# per-junction feature; planning/465b). The D-dependence is the point: a
# long-range junction on a low-complexity flank is refused while a short,
# annotated-length intron on the same flank is not. Measured on the dense
# PTC7+SUS1 pool (planning/649): alpha=0.01 refuses 19/109 novel junctions,
# all long-range/low-complexity (worst D=8,684 on TTCTTTTTTTTTTTT,
# E_chance=376), with 0/19 refused matching the planning/617 gold catalogue.
# ---------------------------------------------------------------------------

DEFAULT_POOL_FLANK_BP = 15


def flank_information_bits(
    chrom_seq: str,
    intron_start: int,
    intron_end: int,
    flank_bp: int = DEFAULT_POOL_FLANK_BP,
    hp_discount: float = DEFAULT_HP_DISCOUNT,
) -> float:
    """I_eff (bits) of the WORSE of the two exonic flanks of an intron.

    Coordinates are 0-based half-open (intron = ``[intron_start, intron_end)``),
    matching the pool convention. Flanks are exonic: ``flank_bp`` bases ending
    at ``intron_start`` (exon-1 side) and starting at ``intron_end`` (exon-2
    side). I_eff is strand-safe: composition entropy and the dinucleotide
    distribution are invariant under reverse complement.
    """
    left = chrom_seq[max(0, intron_start - flank_bp):intron_start]
    right = chrom_seq[intron_end:intron_end + flank_bp]
    return min(
        effective_information_bits(left, hp_discount=hp_discount),
        effective_information_bits(right, hp_discount=hp_discount),
    )


def junction_chance_expectation(
    chrom_seq: str,
    intron_start: int,
    intron_end: int,
    flank_bp: int = DEFAULT_POOL_FLANK_BP,
    hp_discount: float = DEFAULT_HP_DISCOUNT,
) -> float:
    """E_chance = D * 2**(-I_eff) for one junction (see module comment above)."""
    d = max(0, intron_end - intron_start)
    i_eff = flank_information_bits(
        chrom_seq, intron_start, intron_end,
        flank_bp=flank_bp, hp_discount=hp_discount,
    )
    return d * (2.0 ** -min(i_eff, _MAX_I_EFF_EXP))


def filter_pool_by_flank_complexity(
    all_junctions,
    annotated_set,
    genome: Dict[str, str],
    alpha: float,
    flank_bp: int = DEFAULT_POOL_FLANK_BP,
    hp_discount: float = DEFAULT_HP_DISCOUNT,
):
    """Drop NOVEL pool junctions whose span exceeds what their worse exonic
    flank can distinguish from chance (``E_chance > alpha``).

    Annotated junctions are always retained. Junctions on chromosomes missing
    from ``genome`` are retained (conservative: never drop for lack of
    evidence). Returns ``(filtered_set, n_refused)``; the input sets are not
    mutated. ``alpha=None`` is handled by callers as gate-off (byte parity).
    """
    refused = 0
    kept = set()
    for j in all_junctions:
        if j in annotated_set:
            kept.add(j)
            continue
        chrom_seq = genome.get(j[0])
        if chrom_seq is None:
            kept.add(j)
            continue
        e_chance = junction_chance_expectation(
            chrom_seq, int(j[1]), int(j[2]),
            flank_bp=flank_bp, hp_discount=hp_discount,
        )
        if e_chance <= alpha:
            kept.add(j)
        else:
            refused += 1
    return kept, refused


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


_CANONICAL_DINUCS = frozenset(
    # (donor-side, acceptor-side) on the FORWARD genome: GT/GC..AG for
    # plus-strand introns, CT..AC / CT..GC for minus-strand.
    [('GT', 'AG'), ('GC', 'AG'), ('CT', 'AC'), ('CT', 'GC')]
)


def is_canonical_junction(chrom_seq: str, start: int, end: int) -> bool:
    """Canonical splice grammar at the WRITTEN coordinates (either strand)."""
    if start < 0 or end > len(chrom_seq):
        return False
    return (chrom_seq[start:start + 2].upper(),
            chrom_seq[end - 2:end].upper()) in _CANONICAL_DINUCS


def canonical_in_class(chrom_seq: str, start: int, end: int, max_shift: int = 60) -> bool:
    """True iff ANY member of the junction's ambiguity class is canonical."""
    l, r = ambiguity_window(chrom_seq, start, end, max_shift=max_shift)
    return any(is_canonical_junction(chrom_seq, start + d, end + d)
               for d in range(-l, r + 1))


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

# --- Optional numba kernel (the 2026-08-12 cDNA-throughput fix) -------------
#
# Measured: the resolver on Path-A cDNA consensus molecules runs at ~1-4
# reads/s (vs ~1,600 r/s on trimmed DRS) because long high-complexity clips
# legitimately pass the gate and each one scores ~hundreds of candidates with
# this pure-Python DP (cProfile: 122 s of a 125 s / 50-read run inside this
# function). Same knob + guard pattern as splice_aware_5prime's kernel
# (RECTIFY_HP_ED_NUMBA; import-guarded; default OFF for the spawn-RSS reason
# documented there — production job scripts opt in explicitly). This kernel
# prunes EVERY row (not every 8th) and returns cutoff + 1.0 at the same
# detection point as the Python loop, so results are bit-identical in both
# the pruned and unpruned regimes.

_HP_ED_NUMBA_REQUESTED = os.environ.get('RECTIFY_HP_ED_NUMBA', '0').strip() not in (
    '', '0', 'false', 'False', 'no', 'off',
)
_hp_ed_bounded_numba = None
_HP_ED_NUMBA_MIN_CELLS = 400

if _HP_ED_NUMBA_REQUESTED:
    try:
        import numba as _numba_mod
        import numpy as _np_mod

        @_numba_mod.njit(cache=True)
        def _hp_ed_bounded_numba(a1, a2, cutoff):
            """Bit-identical twin of the Python loop below (uint8 inputs,
            already truncated to _HP_ED_MAX_LEN; every-row pruning)."""
            n = len(a1)
            m = len(a2)
            if n == 0:
                return float(m)
            if m == 0:
                return float(n)
            dcost = _np_mod.empty(n, dtype=_np_mod.float64)
            for i in range(1, n + 1):
                if (i >= 2 and a1[i - 2] == a1[i - 1]) or (i < n and a1[i - 1] == a1[i]):
                    dcost[i - 1] = 0.5
                else:
                    dcost[i - 1] = 1.0
            icost = _np_mod.empty(m, dtype=_np_mod.float64)
            for j in range(1, m + 1):
                if (j >= 2 and a2[j - 2] == a2[j - 1]) or (j < m and a2[j - 1] == a2[j]):
                    icost[j - 1] = 0.5
                else:
                    icost[j - 1] = 1.0
            prev = _np_mod.empty(m + 1, dtype=_np_mod.float64)
            curr = _np_mod.empty(m + 1, dtype=_np_mod.float64)
            prev[0] = 0.0
            for j in range(1, m + 1):
                prev[j] = prev[j - 1] + icost[j - 1]
            prune = cutoff >= 0.0
            for i in range(1, n + 1):
                dc = dcost[i - 1]
                c1 = a1[i - 1]
                curr[0] = prev[0] + dc
                rmin = curr[0]
                for j in range(1, m + 1):
                    if c1 == a2[j - 1]:
                        v = prev[j - 1]
                    else:
                        v = prev[j - 1] + 1.0
                        d = prev[j] + dc
                        if d < v:
                            v = d
                        ins = curr[j - 1] + icost[j - 1]
                        if ins < v:
                            v = ins
                    curr[j] = v
                    if v < rmin:
                        rmin = v
                if prune and rmin > cutoff:
                    return cutoff + 1.0
                for j in range(m + 1):
                    prev[j] = curr[j]
            return prev[m]
    except ImportError:  # pragma: no cover - numba genuinely absent
        _hp_ed_bounded_numba = None


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

    if (_hp_ed_bounded_numba is not None
            and n * m >= _HP_ED_NUMBA_MIN_CELLS):
        a1 = _np_mod.frombuffer(s1.encode('ascii', 'replace'), dtype=_np_mod.uint8)
        a2 = _np_mod.frombuffer(s2.encode('ascii', 'replace'), dtype=_np_mod.uint8)
        return float(_hp_ed_bounded_numba(a1, a2, float(cutoff)))

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
