"""
Per-junction scoring + candidate pool construction for splice refinement.

Defines:

- Splice-site canonicality (``_canonical_tier``, ``_3ss_tier_from_rna_trinucleotide``,
  ``_CANONICAL_5SS_*``, ``_3SS_CANONICAL_*``, ``_CANONICAL_HP_PRIOR``).
- Junction pool plumbing (``collect_junctions_from_bam``,
  ``build_junction_pool``, ``_build_junction_index``, ``_candidates_near``).
- Scoring kernels (``_score_hp_anchored``, ``_score_junction``,
  ``_score_junction_fallback``) that consume :class:`HpPenaltyTable` and
  the per-position cost helpers from ``hp_penalty``.

CIGAR-op constants and the ``Junction`` type alias also live here so the
scoring helpers and the refiner orchestrator (which re-imports them) share
a single source of truth.

Dependency direction: ``hp_penalty → junction_scoring → junction_refiner``.
Do not import from ``junction_refiner`` here.
"""

from __future__ import annotations

import logging
import os
import time
from bisect import bisect_left, bisect_right
from collections import Counter
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

import numpy as _np
import pysam

from ...utils.genome import standardize_chrom_name
from .hp_penalty import (
    HpPenaltyTable,
    _NUMBA_AVAILABLE,
    _hp_run_length,
    _precompute_del_costs,
    _score_hp_dp_numba,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
Junction = Tuple[str, int, int]   # (chrom, intron_start, intron_end)


# ---------------------------------------------------------------------------
# CIGAR op constants
# ---------------------------------------------------------------------------
_M = 0   # alignment match (M)
_I = 1   # insertion to reference
_D = 2   # deletion from reference
_N = 3   # intron skip
_S = 4   # soft clip
_H = 5   # hard clip
_EQ = 7  # sequence match (=)
_X = 8   # sequence mismatch (X)

_QUERY_CONSUMING: FrozenSet[int] = frozenset([_M, _I, _S, _EQ, _X])
_REF_CONSUMING:   FrozenSet[int] = frozenset([_M, _D, _EQ, _X])


# ---------------------------------------------------------------------------
# Splice-site canonicality
# ---------------------------------------------------------------------------

# Plus-strand: intron starts GT, ends AG
# Minus-strand (genomic orientation): intron starts CT (RC of AG), ends AC (RC of GT)
_CANONICAL_5SS_PLUS  = frozenset(['GT', 'GC'])   # donor dinucleotide (intron_start:+2)
_CANONICAL_5SS_MINUS = frozenset(['AC', 'GC'])   # genomic = RC of donor

# Canonical HP prior: edit-distance discount applied to canonical-tier junctions
# (tier < 4) when the current N-op is non-canonical (tier_beats_alt=True).  The
# value 0.5 equals one Nanopore HP-deletion equivalent — the expected noise floor
# for splice-site scoring.  This ensures canonical junctions win over non-canonical
# ones within the HP noise floor regardless of which penalty table is in use.
# Non-canonical candidates must beat canonical alternatives by more than this amount.
_CANONICAL_HP_PRIOR = 0.5

# 3'SS acceptor quality hierarchy based on yeast splicing observations.
# In RNA orientation (reading the intron 5'→3'):
#   Tier 0: YAG (Y = pyrimidine: C or T) — most common and highest-efficiency
#   Tier 1: RAG (R = purine: A or G)
#   Tier 2: NBG (N = any, B = C/G/T) — non-canonical but observed in yeast
#   Tier 3: NAT — very rare non-canonical
#   Tier 4: other
#
# Genomic representation (+ strand): RC of RNA acceptor
# + strand: last 2/3 bases of intron = intron_end-2:intron_end
# - strand: first 2/3 bases of intron = intron_start:intron_start+3
#
# The 3-base genomic context needed for full YAG/RAG/NBG classification:
#   Plus strand (last 3 of intron):  s[-3:] in RNA = g[intron_end-3:intron_end]
#   Minus strand (first 3 of intron): s[0:3] in RNA = RC(g[intron_start:intron_start+3])
_3SS_CANONICAL_PLUS  = frozenset(['AG'])   # YAG/RAG both end in AG; use dinucleotide for tier 0
_3SS_CANONICAL_MINUS = frozenset(['CT'])   # RC(AG) = CT

# 3'SS quality tiers for tiebreaking (lower = better).
# Keyed on the 3-base RNA sequence (last 3 of intron in RNA orientation).
# YAG: CAG, TAG → 0
# RAG: AAG, GAG → 1
# NBG: T/C/G/A + C/G/T + G → 2
# NAT: NAT patterns → 3
# other → 4

def _3ss_tier_from_rna_trinucleotide(tri: str) -> int:
    """Return 3'SS quality tier (0=best) from the RNA-orientation 3-nt intron tail."""
    if len(tri) < 3:
        return 4
    b1, b2, b3 = tri[0], tri[1], tri[2]
    if b3 == 'G':
        if b2 == 'A':
            if b1 in ('C', 'T'):   # YAG
                return 0
            if b1 in ('A', 'G'):   # RAG
                return 1
            return 2  # NAG
        if b2 in ('G', 'C', 'T'):   # NBG (B = not A)
            return 2
    if b3 == 'T' and b2 == 'A':     # NAT
        return 3
    return 4


def _canonical_tier(
    intron_start: int,
    intron_end: int,
    genome_seq: str,
    strand: str,
) -> int:
    """Return splice-site quality tier (0=best): considers 5'SS GT-AG and 3'SS YAG hierarchy.

    Returns the combined tier:
      0 — 5'SS GT/GC AND 3'SS YAG (canonical)
      1 — 5'SS GT/GC AND 3'SS RAG
      2 — 5'SS GT/GC AND 3'SS NBG, OR non-canonical 5'SS with YAG
      3 — various semi-canonical combinations
      4 — non-canonical at both sites
    """
    gs = len(genome_seq)
    if strand == '+':
        di5  = genome_seq[intron_start:intron_start + 2].upper() if intron_start + 2 <= gs else 'NN'
        tri3 = genome_seq[intron_end - 3:intron_end].upper() if intron_end >= 3 else 'NNN'
        ok5  = di5 in _CANONICAL_5SS_PLUS
        t3   = _3ss_tier_from_rna_trinucleotide(tri3)
    else:
        # Minus strand: 5'SS is at intron_end (genomic), 3'SS at intron_start
        di5  = genome_seq[intron_end - 2:intron_end].upper() if intron_end >= 2 else 'NN'
        # 3'SS trinucleotide in RNA = RC of g[intron_start:intron_start+3]
        tri3_genomic = genome_seq[intron_start:intron_start + 3].upper() if intron_start + 3 <= gs else 'NNN'
        tri3_rna     = tri3_genomic[::-1].translate(str.maketrans('ACGT', 'TGCA'))
        ok5 = di5 in _CANONICAL_5SS_MINUS
        t3  = _3ss_tier_from_rna_trinucleotide(tri3_rna)

    # Combined tier: 5'SS canonical and 3'SS quality
    if ok5:
        return t3          # 0=YAG, 1=RAG, 2=NBG, 3=NAT, 4=other
    else:
        return 4 + t3      # non-canonical 5'SS: always worse than any canonical 5'SS


# ---------------------------------------------------------------------------
# Junction pool construction
# ---------------------------------------------------------------------------

def collect_junctions_from_bam(
    bam_path: str,
    chrom_filter: Optional[str] = None,
    max_junction_size: Optional[int] = None,
) -> Set[Tuple[str, int, int]]:
    """Extract all N-op intervals from a BAM file.

    Args:
        bam_path:     Path to sorted, indexed BAM.
        chrom_filter: If given, only include junctions on this chromosome.

    Returns:
        Set of (chrom, intron_start, intron_end) tuples (0-based half-open).
    """
    junctions: Set[Tuple[str, int, int]] = set()
    try:
        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            for read in bam:
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or not read.cigartuples
                ):
                    continue
                chrom = standardize_chrom_name(read.reference_name)
                if chrom_filter and chrom != chrom_filter:
                    continue
                pos = read.reference_start
                for op, length in read.cigartuples:
                    if op == _N:
                        if max_junction_size is None or length <= max_junction_size:
                            junctions.add((chrom, pos, pos + length))
                    if op in _REF_CONSUMING:
                        pos += length
    except Exception as exc:
        logger.warning("collect_junctions_from_bam(%s): %s", bam_path, exc)
    return junctions


def collect_junction_counts_from_bam(
    bam_path: str,
    chrom_filter: Optional[str] = None,
    max_junction_size: Optional[int] = None,
) -> Counter[Junction]:
    """Count N-op intervals from a BAM file."""
    junctions: Counter[Junction] = Counter()
    try:
        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            for read in bam:
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or not read.cigartuples
                ):
                    continue
                chrom = standardize_chrom_name(read.reference_name)
                if chrom_filter and chrom != chrom_filter:
                    continue
                pos = read.reference_start
                for op, length in read.cigartuples:
                    if op == _N:
                        if max_junction_size is None or length <= max_junction_size:
                            junctions[(chrom, pos, pos + length)] += 1
                    if op in _REF_CONSUMING:
                        pos += length
    except Exception as exc:
        logger.warning("collect_junction_counts_from_bam(%s): %s", bam_path, exc)
    return junctions


def build_junction_pool(
    aligner_bams: List[str],
    annotated_junctions: Set[Tuple],
    chrom_filter: Optional[str] = None,
    min_observed_support: int = 1,
    max_junction_size: Optional[int] = None,
) -> Tuple[Set[Junction], Set[Junction]]:
    """Build union of annotated + per-aligner junctions.

    Args:
        aligner_bams:         Paths to per-aligner BAMs (pre-consensus).
        annotated_junctions:  Set of (chrom, start, end[, strand]) tuples from
                              ``load_annotated_junctions()``.
        chrom_filter:         Optional chromosome to restrict to.
        min_observed_support: Minimum total read support across all aligner BAMs
                              required for observed, non-annotated junctions.
                              Annotated junctions are always retained.
        max_junction_size:    Optional maximum observed N-op length to include.
                              Annotated junctions are always retained.

    Returns:
        (all_junctions, annotated_set) where annotated_set is the normalised
        3-tuple subset of annotated_junctions for fast membership testing.
    """
    # Normalise annotated junctions to 3-tuples (drop optional strand element)
    annot_3: Set[Junction] = set()
    for j in annotated_junctions:
        if len(j) >= 3:
            chrom = standardize_chrom_name(str(j[0]))
            if chrom_filter and chrom != chrom_filter:
                continue
            annot_3.add((chrom, int(j[1]), int(j[2])))

    all_j: Set[Junction] = set(annot_3)
    min_observed_support = max(1, int(min_observed_support))

    if not aligner_bams:
        # No aligner BAMs — pool is just the annotated set.
        pass
    elif len(aligner_bams) == 1:
        # Single aligner — avoid fork overhead.
        observed_counts = collect_junction_counts_from_bam(
            aligner_bams[0],
            chrom_filter=chrom_filter,
            max_junction_size=max_junction_size,
        )
        all_j.update(
            j for j, n in observed_counts.items()
            if n >= min_observed_support
        )
    else:
        # Multiple aligners — process in parallel (one process per BAM).
        observed_counts: Counter[Junction] = Counter()
        try:
            from concurrent.futures import ProcessPoolExecutor
            n_workers = min(len(aligner_bams), os.cpu_count() or 4)
            with ProcessPoolExecutor(max_workers=n_workers) as ex:
                futures = [
                    ex.submit(
                        collect_junction_counts_from_bam,
                        bp,
                        chrom_filter,
                        max_junction_size,
                    )
                    for bp in aligner_bams
                ]
                for fut in futures:
                    observed_counts.update(fut.result())
        except (OSError, PermissionError, NotImplementedError) as exc:
            logger.warning(
                "build_junction_pool: ProcessPoolExecutor unavailable (%s); "
                "falling back to sequential BAM scans",
                exc,
            )
            for bp in aligner_bams:
                observed_counts.update(
                    collect_junction_counts_from_bam(
                        bp,
                        chrom_filter,
                        max_junction_size=max_junction_size,
                    )
                )
        all_j.update(
            j for j, n in observed_counts.items()
            if n >= min_observed_support
        )

    logger.debug(
        "build_junction_pool: %d annotated, %d observed-support>=%d, %d total",
        len(annot_3), len(all_j - annot_3), min_observed_support, len(all_j),
    )
    return all_j, annot_3


# ---------------------------------------------------------------------------
# Interval index for fast radius lookup
# ---------------------------------------------------------------------------

def _build_junction_index(
    junctions: Set[Junction],
) -> Dict[str, List[Tuple[int, int]]]:
    """Build per-chromosome sorted list of (intron_start, intron_end) pairs."""
    idx: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, js, je in junctions:
        idx.setdefault(chrom, []).append((js, je))
    for chrom in idx:
        idx[chrom].sort()
    return idx


def _candidates_near(
    idx: Dict[str, List[Tuple[int, int]]],
    chrom: str,
    ns: int,
    ne: int,
    radius: int,
    start_radius: Optional[int] = None,
    end_radius: Optional[int] = None,
    max_junction_size: Optional[int] = None,
) -> List[Tuple[int, int]]:
    """Return all (js, je) near (ns, ne) on *chrom*.

    Args:
        idx:           Junction index from _build_junction_index.
        chrom:         Chromosome name.
        ns:            Current intron start.
        ne:            Current intron end.
        radius:        Discovery radius for intron_start (bp).  Used when
                       *start_radius* is not supplied.  Also serves as the
                       upper bound for scanning the sorted list.
        start_radius:  Maximum allowed shift of intron_start (bp).  When
                       supplied, only candidates with |js - ns| ≤ start_radius
                       are returned.  Defaults to *radius*.
        end_radius:    Maximum allowed shift of intron_end (bp).  When
                       supplied, only candidates with |je - ne| ≤ end_radius
                       are returned.  Defaults to *radius*.
        max_junction_size: Optional maximum candidate intron length.

    Separating *radius* from *start_radius*/*end_radius* allows a large
    discovery window (e.g. 5000 bp) for annotated junctions while keeping the
    boundary-shift constraints tight (e.g. 50 bp) to avoid pairing the current
    N-op with junctions from neighbouring genes that happen to score well due
    to sequence coincidence.
    """
    if start_radius is None:
        start_radius = radius
    if end_radius is None:
        end_radius = radius

    entries = idx.get(chrom, [])
    if not entries:
        return []

    # The index is sorted by intron_start.  Bound the scan to the same start
    # interval the legacy filter accepted instead of scanning from chromosome
    # start for every N-op.
    start_window = min(radius, start_radius)
    left = bisect_left(entries, (ns - start_window, -10**18))
    right = bisect_right(entries, (ns + start_window, 10**18))

    results = []
    for js, je in entries[left:right]:
        if max_junction_size is not None and (je - js) > max_junction_size:
            continue
        if abs(js - ns) <= start_radius and abs(je - ne) <= end_radius:
            results.append((js, je))
    return results


# ---------------------------------------------------------------------------
# Core scoring: HP-aware semi-global alignment helper
# ---------------------------------------------------------------------------

def _score_hp_anchored(
    query: str,
    ref: str,
    sub: float = 1.0,
    del_normal: float = 1.0,
    del_hp: float = 0.5,
    ins: float = 1.25,
    hp_min_run: int = 4,
    penalty_table: Optional[HpPenaltyTable] = None,
    genome_seq: Optional[str] = None,
    ref_genome_start: int = 0,
    ref_genome_rev: bool = False,
    precomputed_del_costs: Optional[List[float]] = None,
) -> float:
    """Left-anchored HP-aware semi-global alignment.

    Aligns *query* to *ref* with the left end of *ref* fixed (no free prefix)
    and the right suffix of *ref* free (any unmatched trailing ref bases cost
    nothing).  Returns the minimum total edit cost over all possible ref end
    positions.

    For intron_start-proximal alignment (anchor the last query base to
    ``intron_start - 1``), call with reversed sequences::

        _score_hp_anchored(query[::-1], ref[::-1])

    Costs (per base):

    * substitution : ``sub``         (default 1.0)
    * deletion     : ``del_normal``  (default 1.0), or ``del_hp`` (default 0.5)
                     if the ref base lies within a homopolymer run of length
                     >= ``hp_min_run`` in *ref* (reflects Nanopore HP dropout).
                     Overridden per-HP-length when *penalty_table* is provided.
    * insertion    : ``ins``         (default 1.25 — insertions are rarer than
                     HP deletions in Nanopore data).
                     Overridden per-HP-length when *penalty_table* is provided.

    Args:
        penalty_table: Optional empirical penalty table.  When provided, the
                       per-position del/ins costs are looked up from the table
                       rather than computed from the heuristic step function.
        genome_seq:        Full chromosome sequence.  When provided together with
                           *penalty_table*, enables base-class-aware (AT/CG) and
                           STR-context cost lookups.  Without it, costs fall back
                           to base-agnostic HP-length lookup.
        ref_genome_start:  Index in *genome_seq* corresponding to ref[0] (forward).
        ref_genome_rev:    Set True when *ref* is reversed (e.g. the tier-2 intron-end
                           reference).  Position j in ref then maps to genome position
                           ref_genome_start - j.
    """
    Q, R = len(query), len(ref)
    if Q == 0:
        return 0.0
    if R == 0:
        if penalty_table is not None:
            return Q * penalty_table.ins_cost(1, query[0] if Q > 0 else 'A')
        return Q * ins  # all insertions, no ref to match

    # Precompute per-position deletion cost for each ref base.
    # When genome_seq is provided: use full-genome HP run lengths and STR context.
    # Callers may supply precomputed_del_costs to skip this work when scoring the
    # same reference window multiple times (e.g. the k-loop in _score_junction).
    if precomputed_del_costs is not None:
        del_costs = precomputed_del_costs
    else:
        del_costs = _precompute_del_costs(
            ref, genome_seq, ref_genome_start, ref_genome_rev,
            penalty_table, hp_min_run, del_normal, del_hp,
        )

    # Insertion cost: use table lookup for first query base HP context when table provided.
    # For the DP, we need per-query-position ins cost; we approximate using the HP
    # context of each query base.
    if penalty_table is not None:
        ins_costs: List[float] = [
            penalty_table.ins_cost(_hp_run_length(query, i), query[i]) for i in range(Q)
        ]
    else:
        ins_costs = [ins] * Q

    # Fast path: use Numba JIT when available.
    if _NUMBA_AVAILABLE:
        q_arr = _np.frombuffer(query.upper().encode('ascii'), dtype=_np.uint8)
        r_arr = _np.frombuffer(ref.upper().encode('ascii'), dtype=_np.uint8)
        dc_arr = _np.array(del_costs, dtype=_np.float64)
        ic_arr = _np.array(ins_costs, dtype=_np.float64)
        return float(_score_hp_dp_numba(q_arr, r_arr, dc_arr, ic_arr))

    INF = float('inf')
    # One-row rolling DP (space O(R))
    prev = [INF] * (R + 1)
    prev[0] = 0.0
    for j in range(1, R + 1):
        prev[j] = prev[j - 1] + del_costs[j - 1]  # leading ref deletions

    for i in range(1, Q + 1):
        curr = [INF] * (R + 1)
        curr[0] = i * ins_costs[i - 1]  # i query insertions with no ref consumed
        qi = query[i - 1].upper()
        for j in range(1, R + 1):
            cost_sub = 0.0 if qi == ref[j - 1].upper() else sub
            diag  = prev[j - 1] + cost_sub        # match / mismatch
            above = prev[j]     + ins_costs[i - 1]  # insertion (query base, no ref)
            left  = curr[j - 1] + del_costs[j - 1]  # deletion (ref base, no query)
            curr[j] = min(diag, above, left)
        prev = curr

    return min(prev)  # free right suffix: best end column over all j in [0, R]


# ---------------------------------------------------------------------------
# Core scoring: bilateral anchor + rescue alignment
# ---------------------------------------------------------------------------

def _score_junction(
    query: str,
    q_split: int,
    intron_start: int,
    intron_end: int,
    genome_seq: str,
    hp_pen: float,
    W: int,            # kept for API compatibility; ignored by this implementation
    max_slide: int,    # kept for API compatibility; ignored by this implementation
    current_ns: int = -1,  # current N-op intron_start (anchor ref position for tier 2)
    profile: Optional[Any] = None,
    penalty_table: Optional[HpPenaltyTable] = None,
) -> Tuple[float, int]:
    """Score a candidate junction using HP-aware intron_end-anchored alignment.

    **Algorithm:**

    For each ``k`` in ``[0, L)`` (L = min(len(query[q_split:]), 30)):

    * **t1(k)**: ``rescue[k:]`` aligned left-anchored to ``g[intron_end : intron_end+buf]``.
      The free right suffix of the ref soft-clips terminal noise without penalty.
      Asks: "does the suffix start cleanly at intron_end?"

    * **t2(k)**: ``rescue[:k]`` reversed aligned left-anchored to
      ``g[intron_end-buf : intron_end]`` reversed (reversal trick: left-anchored DP
      on reversed sequences ≡ right-anchored).  Asks: "does the prefix end cleanly
      at intron_end - 1?"  When k=0, t2=0.

    * ``score(k) = t1(k) + t2(k)``; best ``k`` minimises this.

    Bilateral scoring ensures the aligner accounts for ALL rescue bases.  Without t2,
    a degenerate k near L-1 can give t1≈0 via a 1-char coincidental match even for
    the wrong junction.  t2 penalises the ignored prefix, preventing this artefact
    without any candidate pre-filtering.

    Args:
        query:        Full BAM query sequence.
        q_split:      Current query index of the start of the rescue window.
        intron_start: Candidate intron start (0-based, = js); unused by scorer.
        intron_end:   Candidate intron end (0-based, exclusive, = je).
        genome_seq:   Full chromosome sequence.
        hp_pen:       Unused (retained for API compatibility).
        W:            Unused (retained for API compatibility).
        max_slide:    Unused (retained for API compatibility).
        current_ns:   Unused (retained for API compatibility).

    Returns:
        ``(best_score, 0)`` where ``best_score`` is ``min_k (t1(k) + t2(k))``.
    """
    if profile is not None:
        profile.inc('score_junction_calls')
    _t_total = time.perf_counter() if profile is not None else 0.0
    gs = len(genome_seq)

    # Cap rescue length: only the first ~30 bp are needed to identify the junction.
    _t_setup = time.perf_counter() if profile is not None else 0.0
    _MAX_RESCUE = 30
    rescue = query[q_split : q_split + _MAX_RESCUE].upper()
    L = len(rescue)

    if L == 0:
        if profile is not None:
            profile.inc('score_junction_empty_rescue')
            profile.add_time('score_junction_total', time.perf_counter() - _t_total)
        return 0.0, 0

    # Tier-1 reference: sequence starting at intron_end (exon2 body after the intron).
    _BUF = max(L + 20, 50)
    ref_exon2_start = genome_seq[intron_end : min(gs, intron_end + _BUF)].upper()

    # Tier-2 reference: sequence immediately before intron_end (end of intron).
    # Used reversed so that the last base of rescue[:k] anchors to intron_end-1.
    ref_intron_end = genome_seq[max(0, intron_end - _BUF) : intron_end].upper()
    ref_intron_end_rev = ref_intron_end[::-1]
    if profile is not None:
        profile.add_time('score_junction_setup', time.perf_counter() - _t_setup)

    # Precompute per-position deletion costs for the two fixed reference windows.
    # del_costs depend only on the ref sequence and genome position, not on the
    # query, so they are identical across all k values in the loop below.
    # Computing them once here and passing them in avoids ~58 redundant scans of
    # _hp_run_length / _str_repeat_info per junction call.
    _t_delcosts = time.perf_counter() if profile is not None else 0.0
    del_costs_fwd = _precompute_del_costs(
        ref_exon2_start, genome_seq, intron_end, False, penalty_table,
    )
    del_costs_rev = _precompute_del_costs(
        ref_intron_end_rev, genome_seq, intron_end - 1, True, penalty_table,
    )
    if profile is not None:
        profile.add_time('score_junction_del_costs', time.perf_counter() - _t_delcosts)

    best_score = float("inf")

    for k in range(L):
        if profile is not None:
            profile.inc('score_junction_k_iterations')
        # Tier-1: rescue[k:] must align cleanly to exon2 starting at intron_end.
        q1 = rescue[k:]
        _t_t1 = time.perf_counter() if profile is not None else 0.0
        t1 = _score_hp_anchored(
            q1, ref_exon2_start, penalty_table=penalty_table,
            genome_seq=genome_seq, ref_genome_start=intron_end,
            precomputed_del_costs=del_costs_fwd,
        )
        if profile is not None:
            profile.inc('score_hp_anchored_t1_calls')
            profile.add_time('score_hp_anchored_t1', time.perf_counter() - _t_t1)

        # Tier-2: rescue[:k] must end cleanly at intron_end - 1 (last intronic base).
        # Reversal trick: left-anchored DP on reversed sequences ≡ right-anchored DP.
        if k > 0:
            _t_t2 = time.perf_counter() if profile is not None else 0.0
            t2 = _score_hp_anchored(
                rescue[:k][::-1], ref_intron_end_rev, penalty_table=penalty_table,
                genome_seq=genome_seq, ref_genome_start=intron_end - 1,
                ref_genome_rev=True,
                precomputed_del_costs=del_costs_rev,
            )
            if profile is not None:
                profile.inc('score_hp_anchored_t2_calls')
                profile.add_time('score_hp_anchored_t2', time.perf_counter() - _t_t2)
        else:
            t2 = 0.0

        score = t1 + t2
        if score < best_score:
            best_score = score
            if best_score == 0.0:
                if profile is not None:
                    profile.inc('score_junction_perfect_breaks')
                break  # perfect match; can't improve

    if profile is not None:
        profile.add_time('score_junction_total', time.perf_counter() - _t_total)
    return best_score, 0


def _score_junction_fallback(
    query: str, q_split: int, intron_start: int, intron_end: int,
    genome_seq: str, hp_pen: float,
    penalty_table: Optional[HpPenaltyTable] = None,
) -> float:
    """Simple edit-distance fallback when local_aligner is unavailable."""
    from .hp_penalty import _hp_edit_distance

    gs = len(genome_seq)
    rescue = query[q_split:].upper()
    L = len(rescue)
    if L == 0:
        return 0.0
    best = float('inf')
    for k in range(L + 1):
        ref_gap = genome_seq[max(0, intron_start - k): intron_start].upper()
        ref_ex1 = genome_seq[intron_end: min(gs, intron_end + L - k)].upper()
        ed2 = _hp_edit_distance(rescue[:k], ref_gap, hp_pen,
                                penalty_table=penalty_table) if rescue[:k] or ref_gap else 0.0
        ed1 = _hp_edit_distance(rescue[k:], ref_ex1, hp_pen,
                                penalty_table=penalty_table) if rescue[k:] or ref_ex1 else 0.0
        best = min(best, ed2 + ed1)
    return best
