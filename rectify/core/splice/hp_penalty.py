"""
Homopolymer-aware edit-distance primitives for splice junction refinement.

Low-level scoring kernels: HP/STR run detection, per-position deletion cost
precomputation, the ``HpPenaltyTable`` empirical lookup, the heuristic
``_hp_edit_distance`` DP, and the Numba-compiled rolling NW DP used by
``junction_scoring._score_hp_anchored``.

Designed to be imported by ``junction_scoring`` (and indirectly by the
refinement orchestrator in ``junction_refiner``). No dependencies on
either of those modules — keep this layer free of upward imports so the
``hp_penalty → scoring → refiner`` dependency direction stays clean.
"""

from __future__ import annotations

import csv
from typing import Dict, List, Optional, Tuple

import numpy as _np


# ---------------------------------------------------------------------------
# Optional Numba JIT for the DP inner loop
# ---------------------------------------------------------------------------
# When numba is installed, _score_hp_dp_numba() compiles the rolling NW DP
# to native machine code (~86x faster than pure Python).  cache=True writes
# the compiled bytecode to __pycache__/, so compilation (~1.2 s) only happens
# once per environment.  Falls back to the pure-Python DP when unavailable.

try:
    import numba as _numba
    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False

# _score_hp_dp_numba is importable in both cases (None when numba is absent) so
# callers can `from .hp_penalty import _score_hp_dp_numba` unconditionally. The
# actual JIT path is guarded by `if _NUMBA_AVAILABLE:` at every call site.
_score_hp_dp_numba = None

if _NUMBA_AVAILABLE:
    @_numba.njit(cache=True)
    def _score_hp_dp_numba(  # noqa: F811  (intentional override of None stub)
        q_arr: _np.ndarray,      # uint8 array of uppercased ASCII codes, length Q
        r_arr: _np.ndarray,      # uint8 array of uppercased ASCII codes, length R
        del_costs_arr: _np.ndarray,  # float64[R] per-position deletion costs
        ins_costs_arr: _np.ndarray,  # float64[Q] per-position insertion costs
    ) -> float:
        """Numba-compiled rolling NW DP for HP-aware semi-global alignment.

        Semantics are identical to the pure-Python loop in _score_hp_anchored:
        left end of ref is fixed (no free prefix), right suffix of ref is free.
        Returns the minimum total cost over all possible ref end positions.
        """
        Q = len(q_arr)
        R = len(r_arr)
        INF = 1e18
        prev = _np.full(R + 1, INF)
        prev[0] = 0.0
        for j in range(1, R + 1):
            prev[j] = prev[j - 1] + del_costs_arr[j - 1]
        for i in range(1, Q + 1):
            curr = _np.full(R + 1, INF)
            curr[0] = _np.float64(i) * ins_costs_arr[i - 1]
            qi = q_arr[i - 1]
            ic = ins_costs_arr[i - 1]
            for j in range(1, R + 1):
                cost_sub = 0.0 if qi == r_arr[j - 1] else 1.0
                d = prev[j - 1] + cost_sub
                a = prev[j]     + ic
                l = curr[j - 1] + del_costs_arr[j - 1]
                val = d
                if a < val:
                    val = a
                if l < val:
                    val = l
                curr[j] = val
            prev = curr
        return prev.min()


# ---------------------------------------------------------------------------
# Homopolymer / STR run detection
# ---------------------------------------------------------------------------

def _hp_run_length(seq: str, pos: int) -> int:
    """Return the length of the homopolymer run containing position *pos*."""
    if pos < 0 or pos >= len(seq):
        return 1  # position outside the contig — neutral HP length, never index-error
    c = seq[pos]
    left = pos
    while left > 0 and seq[left - 1] == c:
        left -= 1
    right = pos
    while right < len(seq) - 1 and seq[right + 1] == c:
        right += 1
    return right - left + 1


def _str_repeat_info(
    ref_seq: str,
    pos: int,
    max_unit: int = 3,
    min_copies: int = 3,
) -> Tuple[Optional[str], int]:
    """Return (canonical_unit, n_copies) if *pos* is in an STR of unit length 2-max_unit.

    Only considers units with ≥2 distinct characters (mononucleotide runs are
    handled by _hp_run_length).  Canonical unit = lexicographically smallest
    rotation (e.g. TATATA → 'AT').  Returns (None, 0) if not in a qualifying STR.
    """
    n = len(ref_seq)
    best_unit: Optional[str] = None
    best_copies = 0

    for unit_len in range(2, max_unit + 1):
        for phase in range(unit_len):
            unit_start = pos - phase
            if unit_start < 0 or unit_start + unit_len > n:
                continue
            unit = ref_seq[unit_start:unit_start + unit_len]
            if len(set(unit)) == 1:
                continue                             # mononucleotide — skip
            left = unit_start
            while left >= unit_len and ref_seq[left - unit_len:left] == unit:
                left -= unit_len
            copies, p = 0, left
            while p + unit_len <= n and ref_seq[p:p + unit_len] == unit:
                copies += 1
                p += unit_len
            if not (left <= pos < left + copies * unit_len):
                continue
            if copies >= min_copies and copies > best_copies:
                best_copies = copies
                best_unit = min(unit[i:] + unit[:i] for i in range(unit_len))

    return (best_unit, best_copies) if best_unit is not None else (None, 0)


# ---------------------------------------------------------------------------
# Empirical penalty table
# ---------------------------------------------------------------------------

class HpPenaltyTable:
    """Empirical HP-context and STR-context penalty table for Nanopore error scoring.

    Loaded from ``penalty_scores.tsv`` (HP penalties, split by base class AT/CG)
    and optionally ``str_penalty_scores.tsv`` (STR repeat penalties).

    Lookup methods:

    * ``del_cost(hp, ref_base)``        — deletion cost in a HP run of length *hp*,
                                          using the AT or CG sub-table based on *ref_base*.
    * ``ins_cost(hp, ref_base)``        — insertion cost likewise.
    * ``str_del_cost(unit, n_copies)``  — deletion cost in a STR context.
    * ``str_ins_cost(unit, n_copies)``  — insertion cost in a STR context.
    """

    def __init__(
        self,
        del_tables: Dict[str, Dict[int, float]],   # {'AT': {hp: pen}, 'CG': {hp: pen}}
        ins_tables: Dict[str, Dict[int, float]],
        str_del_table: Optional[Dict[Tuple[str, int], float]] = None,
        str_ins_table: Optional[Dict[Tuple[str, int], float]] = None,
        default_ins: float = 1.25,
        del_rate_tables: Optional[Dict[str, Dict[int, float]]] = None,
        ins_rate_tables: Optional[Dict[str, Dict[int, float]]] = None,
    ) -> None:
        self._del = del_tables                      # keyed by base_class
        self._ins = ins_tables
        self._str_del = str_del_table or {}
        self._str_ins = str_ins_table or {}
        self._default_ins = default_ins
        self._del_max_hp = {
            bc: max(t.keys()) for bc, t in del_tables.items() if t
        }
        self._ins_max_hp = {
            bc: max(t.keys()) for bc, t in ins_tables.items() if t
        }
        # rate_mean per (base_class, hp) — the CALIBRATED per-reference-column
        # error probability (rate_mean over {M,X,D,I} sums to ~1.0 per context).
        # The C1 length-law gap-OPEN delta is derived from THIS, NOT from the
        # penalty_score column returned by del_cost/ins_cost (which is a
        # reciprocal-rate heuristic, c/rate_mean, NOT -logP, and is incoherent to
        # sum in an additive DP). See dev/C1_DESIGN.md.
        self._del_rate = del_rate_tables or {}
        self._ins_rate = ins_rate_tables or {}

    @classmethod
    def from_tsv(
        cls,
        path: str,
        str_path: Optional[str] = None,
        min_count: int = 100,
        default_ins: float = 1.25,
    ) -> "HpPenaltyTable":
        """Load from ``penalty_scores.tsv`` and optionally ``str_penalty_scores.tsv``.

        The HP penalty file may be in one of two formats:
        - New format (base_class column): rows have 'base_class' ∈ {'AT','CG'}.
        - Legacy format (no base_class): a single set of penalties applied to both
          base classes.

        Args:
            path:        Path to ``penalty_scores.tsv``.
            str_path:    Path to ``str_penalty_scores.tsv`` (optional).
            min_count:   Skip rows with fewer agreed events than this.
            default_ins: Fallback insertion cost for sparse HP lengths.
        """
        del_tables: Dict[str, Dict[int, float]] = {'AT': {}, 'CG': {}}
        ins_tables: Dict[str, Dict[int, float]] = {'AT': {}, 'CG': {}}
        del_rate: Dict[str, Dict[int, float]] = {'AT': {}, 'CG': {}}
        ins_rate: Dict[str, Dict[int, float]] = {'AT': {}, 'CG': {}}

        with open(path, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                op = row.get('op_type', '').strip()
                try:
                    hp      = int(row['hp_length'])
                    penalty = float(row.get('penalty_score') or row.get('penalty') or 0)
                    count   = int(float(row.get('count_total') or row.get('count') or 0))
                    low_count = row.get('low_count', '').strip().lower() in ('true', '1', 'yes')
                except (KeyError, ValueError):
                    continue
                if count < min_count or low_count:
                    continue
                # rate_mean = the calibrated per-column error rate (for the C1
                # length-law log-odds delta); may be absent in legacy files.
                try:
                    rate = float(row.get('rate_mean') or 0) or None
                except ValueError:
                    rate = None
                # Determine which base classes this row applies to
                bc_raw = row.get('base_class', '').strip()
                classes = [bc_raw] if bc_raw in ('AT', 'CG') else ['AT', 'CG']
                for bc in classes:
                    if op in ('D', 'del'):
                        del_tables[bc][hp] = penalty
                        if rate is not None:
                            del_rate[bc][hp] = rate
                    elif op in ('I', 'ins'):
                        ins_tables[bc][hp] = penalty
                        if rate is not None:
                            ins_rate[bc][hp] = rate

        if not del_tables['AT'] and not del_tables['CG']:
            raise ValueError(
                f"HpPenaltyTable: no 'del' rows with count >= {min_count} found in {path}"
            )

        # --- Load STR penalty table (optional) ---
        str_del: Dict[Tuple[str, int], float] = {}
        str_ins: Dict[Tuple[str, int], float] = {}
        if str_path:
            try:
                with open(str_path, newline='') as fh:
                    reader = csv.DictReader(fh, delimiter='\t')
                    for row in reader:
                        op = row.get('op_type', '').strip()
                        try:
                            unit    = row['unit'].strip()
                            copies  = int(float(row['n_copies']))
                            penalty = float(row.get('penalty_score') or row.get('penalty') or 0)
                            count   = int(float(row.get('count') or row.get('count_total') or 0))
                            low_count = row.get('low_count', '').strip().lower() in ('true', '1', 'yes')
                        except (KeyError, ValueError):
                            continue
                        if count < min_count or low_count or penalty <= 0:
                            continue
                        if op in ('D', 'del'):
                            str_del[(unit, copies)] = penalty
                        elif op in ('I', 'ins'):
                            str_ins[(unit, copies)] = penalty
            except FileNotFoundError:
                import warnings
                warnings.warn(f'STR penalty table not found: {str_path}', stacklevel=2)

        return cls(del_tables, ins_tables, str_del, str_ins, default_ins=default_ins,
                   del_rate_tables=del_rate, ins_rate_tables=ins_rate)

    def _base_class(self, ref_base: str) -> str:
        return 'AT' if ref_base.upper() in ('A', 'T') else 'CG'

    def del_cost(self, hp: int, ref_base: str = 'A') -> float:
        """Deletion cost for a ref position in a HP run of *hp* bases."""
        bc = self._base_class(ref_base)
        table = self._del.get(bc) or self._del.get('AT') or {}
        if hp in table:
            return table[hp]
        max_hp = self._del_max_hp.get(bc, 1)
        return table.get(max_hp, table.get(1, 1.0))

    def ins_cost(self, hp: int, ref_base: str = 'A') -> float:
        """Insertion cost for a position in a HP run of *hp* bases."""
        bc = self._base_class(ref_base)
        table = self._ins.get(bc) or self._ins.get('AT') or {}
        if hp in table:
            return table[hp]
        max_hp = self._ins_max_hp.get(bc, 1)
        return table.get(max_hp, self._default_ins)

    # ------------------------------------------------------------------
    # C1 length-law gap-OPEN deltas (the de-novo aligner facet).
    # Derived from rate_mean as a BASELINE-ANCHORED log-odds, NOT from the
    # penalty_score column (which is a reciprocal-rate heuristic, not -logP, and
    # incoherent to sum in an additive DP). Δ(hp) = λ·ln(rate(hp)/rate(1)); zero
    # at hp=1 so ONLY the run-length DEPENDENCE is added on top of the legacy gap
    # cost (the partial that is provably coherent in a global full-consumption DP
    # — see dev/C1_DESIGN.md). rate rises with hp → Δ>0 → cheaper in-run gap.
    # ------------------------------------------------------------------
    def _open_delta(self, rate_tbl: Dict[str, Dict[int, float]],
                    hp: int, ref_base: str, lam: float) -> float:
        bc = self._base_class(ref_base)
        tbl = rate_tbl.get(bc) or rate_tbl.get('AT') or {}
        r1 = tbl.get(1)
        if not r1:                       # no rate_mean (legacy file) -> no delta
            return 0.0
        rh = tbl.get(hp)
        if rh is None:                   # clamp to the largest tabulated hp
            if not tbl:
                return 0.0
            rh = tbl[max(tbl.keys())] if hp > max(tbl.keys()) else r1
        if rh <= 0:
            return 0.0
        import math
        return lam * math.log(rh / r1)

    def del_open_delta(self, hp: int, ref_base: str = 'A', lam: float = 1.0) -> float:
        """Length-law deletion gap-OPEN delta (added to the legacy gap_open)."""
        return self._open_delta(self._del_rate, hp, ref_base, lam)

    def ins_open_delta(self, hp: int, ref_base: str = 'A', lam: float = 1.0) -> float:
        """Length-law insertion gap-OPEN delta. UNVALIDATED — the controlled
        generator injects deletions only (no insertion stratum yet)."""
        return self._open_delta(self._ins_rate, hp, ref_base, lam)

    def str_del_cost(self, unit: str, n_copies: int) -> float:
        """Deletion cost in an STR context (unit, n_copies).

        Falls back to HP=1 deletion cost if STR table is absent or has no
        matching entry.
        """
        if not self._str_del:
            return self.del_cost(1)
        key = (unit, n_copies)
        if key in self._str_del:
            return self._str_del[key]
        # Nearest n_copies for same unit
        same = {k: v for k, v in self._str_del.items() if k[0] == unit}
        if same:
            nearest = min(same, key=lambda k: abs(k[1] - n_copies))
            return same[nearest]
        return self.del_cost(1)

    def str_ins_cost(self, unit: str, n_copies: int) -> float:
        """Insertion cost in an STR context."""
        if not self._str_ins:
            return self.ins_cost(1)
        key = (unit, n_copies)
        if key in self._str_ins:
            return self._str_ins[key]
        same = {k: v for k, v in self._str_ins.items() if k[0] == unit}
        if same:
            nearest = min(same, key=lambda k: abs(k[1] - n_copies))
            return same[nearest]
        return self.ins_cost(1)


# ---------------------------------------------------------------------------
# Per-position deletion cost precomputation
# ---------------------------------------------------------------------------

def _precompute_del_costs(
    ref: str,
    genome_seq: Optional[str],
    ref_genome_start: int,
    ref_genome_rev: bool,
    penalty_table: Optional[HpPenaltyTable],
    hp_min_run: int = 4,
    del_normal: float = 1.0,
    del_hp: float = 0.5,
) -> List[float]:
    """Precompute per-position deletion costs for a reference sequence.

    Extracted from :func:`_score_hp_anchored` so callers that score the same
    reference window multiple times (e.g. the k-loop in :func:`_score_junction`)
    can compute these costs once and reuse them across all calls.

    Args:
        ref:              Reference sequence (already sliced/reversed as needed).
        genome_seq:       Full chromosome sequence, or None.
        ref_genome_start: Index in *genome_seq* corresponding to ref[0].
        ref_genome_rev:   True when *ref* is reversed; position j maps to
                          ``ref_genome_start - j`` in the genome.
        penalty_table:    Optional empirical penalty table.
        hp_min_run:       Minimum run length for HP discount.
        del_normal:       Normal (non-HP) deletion cost.
        del_hp:           HP deletion cost.

    Returns:
        List of length ``len(ref)`` with per-position deletion costs.
    """
    R = len(ref)
    del_costs: List[float] = []
    for j in range(R):
        base = ref[j]
        gp = (ref_genome_start - j) if ref_genome_rev else (ref_genome_start + j)
        if genome_seq is not None and 0 <= gp < len(genome_seq):
            run = _hp_run_length(genome_seq, gp)
            g_base = genome_seq[gp]
        else:
            run = 1
            kk = j - 1
            while kk >= 0 and ref[kk] == base:
                run += 1; kk -= 1
            kk = j + 1
            while kk < R and ref[kk] == base:
                run += 1; kk += 1
            g_base = base
        if penalty_table is not None:
            if run == 1 and genome_seq is not None and 0 <= gp < len(genome_seq):
                su, sc = _str_repeat_info(genome_seq, gp)
                if su is not None:
                    del_costs.append(penalty_table.str_del_cost(su, sc))
                    continue
            del_costs.append(penalty_table.del_cost(run, g_base))
        else:
            del_costs.append(del_hp if run >= hp_min_run else del_normal)
    return del_costs


# ---------------------------------------------------------------------------
# Heuristic HP-aware edit distance
# ---------------------------------------------------------------------------

def _hp_edit_distance(
    s1: str,
    s2: str,
    hp_pen: float = 0.25,   # retained for API compat; not used by new formula
    min_run: int = 4,
    base_pen: float = 0.5,
    del_scale: float = 1.0,
    ins_scale: float = 0.7,
    penalty_table: Optional[HpPenaltyTable] = None,
) -> float:
    """Edit distance calibrated to the Nanopore error profile.

    Nanopore sequencing has an asymmetric error distribution:
    - Deletions: most common (~5%+, especially in homopolymers)
    - Substitutions: moderate (~2-3%)
    - Insertions: least common (~2-3%, lower than deletions)

    **Penalty model:**

    Substitutions always cost 1.0 (no context dependence).

    Deletions (gap in s2 = removing a base from s1):
        Cost = del_scale × base_indel_cost(s1, pos)
        where base_indel_cost = base_pen × min_run / run_length for run ≥ min_run,
        else 1.0.  del_scale=1.0 by default (deletions are the most common error,
        so we don't further discount them beyond the hp-run adjustment).

    Insertions (gap in s1 = adding a base from s2 to query):
        Cost = ins_scale × base_indel_cost(s2, pos)
        ins_scale=0.7 reflects that insertions are less common than deletions
        but more common than substitutions.  A non-hp insertion costs 0.7;
        a hp insertion (run ≥ min_run) costs 0.7 × base_pen × min_run / run.

    **Homopolymer run discount:**

    - Run length < *min_run* (default 4): no discount — short runs are within
      the pore's resolving power (R9.4 k-mer ~5 bp, R10.4 ~10 bp).
    - Run length ≥ *min_run*: discounted by base_pen × min_run / run_length.
      Longer runs are harder to count, so errors cost less.

    Examples (del_scale=1.0, ins_scale=0.7, min_run=4, base_pen=0.5):

    ============  ======  ======  =======
    context       del     ins     sub
    ============  ======  ======  =======
    non-hp        1.00    0.70    1.00
    run=4         0.50    0.35    1.00
    run=5         0.40    0.28    1.00
    run=8         0.25    0.18    1.00
    ============  ======  ======  =======

    The ``hp_pen`` parameter is retained for API compatibility but ignored.

    Args:
        s1:            Query sequence.
        s2:            Reference sequence.
        hp_pen:        Ignored (API compatibility only).
        min_run:       Minimum homopolymer run length to activate the discount.
        base_pen:      Penalty at exactly *min_run* length; scales down for longer.
        del_scale:     Multiplier for deletion costs (default 1.0).
        ins_scale:     Multiplier for insertion costs (default 0.7).
        penalty_table: Optional empirical penalty table from
                       ``HpPenaltyTable.from_tsv()``.  When provided, overrides
                       the heuristic formula for del/ins costs.

    Returns:
        Float edit distance.
    """
    n, m = len(s1), len(s2)
    if n == 0:
        if penalty_table is not None:
            return float(m) * penalty_table.ins_cost(1, s2[0] if m > 0 else 'A')
        return float(m) * ins_scale   # all insertions
    if m == 0:
        if penalty_table is not None:
            return float(n) * penalty_table.del_cost(1, s1[0] if n > 0 else 'A')
        return float(n) * del_scale   # all deletions

    def _base_indel_cost(seq: str, pos: int) -> float:
        run = _hp_run_length(seq, pos)
        return base_pen * min_run / run if run >= min_run else 1.0

    if penalty_table is not None:
        def _del_cost(i: int) -> float:
            pos = i - 1
            run = _hp_run_length(s1, pos)
            return penalty_table.del_cost(run, s1[pos])

        def _ins_cost(j: int) -> float:
            pos = j - 1
            run = _hp_run_length(s2, pos)
            return penalty_table.ins_cost(run, s2[pos])
    else:
        def _del_cost(i: int) -> float:
            """Cost to delete s1[i-1] (gap in s2 = s1 base missing from reference)."""
            return del_scale * _base_indel_cost(s1, i - 1)

        def _ins_cost(j: int) -> float:
            """Cost to insert s2[j-1] (gap in s1 = extra base in reference)."""
            return ins_scale * _base_indel_cost(s2, j - 1)

    dp = [[0.0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + _del_cost(i)
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + _ins_cost(j)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = min(
                    dp[i - 1][j - 1] + 1.0,       # substitution
                    dp[i - 1][j]     + _del_cost(i),  # deletion from s1
                    dp[i][j - 1]     + _ins_cost(j),  # insertion into s1
                )
    return dp[n][m]
