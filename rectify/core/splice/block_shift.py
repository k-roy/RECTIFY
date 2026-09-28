"""Conservative, atomic whole-exon placement for junction refinement.

This supplements local boundary surgery. A gap-free block may shift only to a
unique exact sequence placement with strictly better identity, sufficient
sequence information for the complete searched endpoint window, and supported
canonical boundaries. Internal blocks couple both adjacent introns. The shared
information model is a conservative chance-match heuristic, not a calibrated
likelihood ratio; noisy/gapped blocks and annotated-to-novel moves are deferred.
"""
from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass
from typing import Optional, Tuple

from ...utils.genome import standardize_chrom_name
from .junction_scoring import _candidates_near, _canonical_tier
from .overhang_informativeness import (
    DEFAULT_ALPHA, effective_information_bits, min_self_match_period,
)

_ALIGNED = (0, 7, 8)
_QUERY = (0, 1, 4, 7, 8)
_REF = (0, 2, 3, 7, 8)


@dataclass(frozen=True)
class PlacementState:
    """Immutable original decoded input; includes SEQ/QUAL to reject stale plans."""
    name: str
    reference_id: int
    reference_name: str
    reference_start: int
    flag: int
    cigar: Tuple[Tuple[int, int], ...]
    sequence: str
    qualities: Optional[Tuple[int, ...]]


@dataclass(frozen=True)
class BlockPlacementEvidence:
    query_start: int
    query_end: int
    old_start: int
    new_start: int
    old_matches: int
    new_matches: int
    information_bits: float
    search_trials: int
    left_anchor: int
    right_anchor: int
    changed_junction_ordinals: Tuple[int, ...]


@dataclass(frozen=True)
class AtomicReadPlan:
    """One whole-read replacement; never intermixed with boundary-edit tuples."""
    expected: PlacementState
    strand: str
    reference_start: int
    cigar: Tuple[Tuple[int, int], ...]
    # All Ns, including unchanged Ns and legacy changes: (old_s, old_e, new_s, new_e).
    junction_map: Tuple[Tuple[int, int, int, int], ...]
    evidence: Tuple[BlockPlacementEvidence, ...]


def decoded_copy(read, genome_seq):
    """Decode SEQ '=' against ORIGINAL placement before any scoring or edit."""
    out = read.__copy__()
    if out.query_sequence and '=' in out.query_sequence:
        # Local import avoids loading the correction writer on the common path.
        from ..bam.bam_writer import _decode_eq_seq_inplace
        _decode_eq_seq_inplace(out, {out.reference_name: genome_seq})
    return out


def placement_state(read):
    return PlacementState(
        read.query_name, read.reference_id, read.reference_name, read.reference_start, read.flag,
        tuple(read.cigartuples or ()), read.query_sequence or '',
        None if read.query_qualities is None else tuple(read.query_qualities),
    )


def _geometry(read):
    q, r = 0, read.reference_start
    positions, introns = [], []
    for i, (op, n) in enumerate(read.cigartuples or ()):
        positions.append((q, r))
        if op == 3:
            introns.append((i, r, r + n))
        if op in _QUERY:
            q += n
        if op in _REF:
            r += n
    positions.append((q, r))
    return positions, introns


def _canonical(g, s, e, strand):
    if not (0 <= s < e <= len(g)) or e - s < 4:
        return False
    pair = (g[s:s + 2].upper(), g[e - 2:e].upper())
    return pair in ({('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')} if strand == '+'
                    else {('CT', 'AC'), ('CT', 'GC'), ('GT', 'AT')})


def _supported(entries, pair):
    i = bisect_left(entries, pair)
    return i < len(entries) and entries[i] == pair


def _gapfree_block(read, ordinal, strand):
    positions, introns = _geometry(read)
    ops = read.cigartuples
    ni = introns[ordinal][0]
    if strand == '+':
        partner = ordinal + 1 if ordinal + 1 < len(introns) else None
        lo, hi = ni + 1, introns[partner][0] if partner is not None else len(ops)
        if partner is None:
            while hi > lo and ops[hi - 1][0] in (4, 5):
                hi -= 1
    else:
        partner = ordinal - 1 if ordinal else None
        lo, hi = introns[partner][0] + 1 if partner is not None else 0, ni
        if partner is None:
            while lo < hi and ops[lo][0] in (4, 5):
                lo += 1
    if lo >= hi or any(op not in _ALIGNED or n <= 0 for op, n in ops[lo:hi]):
        return None
    q0, r0 = positions[lo]
    q1, r1 = positions[hi]
    return lo, hi, q0, q1, r0, r1, partner, introns


def plan_block_shifts(original, evolved, junctions_idx, annotated, genome_seq,
                      strand, *, max_boundary_shift=50, search_radius=5000,
                      max_junction_size=None, max_candidates_per_nop=None):
    """Return a whole-read atomic plan, or None when no block earns a move.

    ``original`` and ``evolved`` must already be literal decoded copies. The
    latter contains actual successful legacy edits. Changed legacy Ns and both
    boundaries of each accepted new block are reserved. Discovery always uses
    the current evolved geometry. An internal block's same bases support both
    ends, counted ONCE in the information test, with separate actual anchors.
    """
    # Xb owns recovered-exon coordinates. Until its provenance can be remapped
    # atomically, stand down even for unknown/invalid payloads; legacy boundary
    # refinement retains its existing behavior. Do not clear or reinterpret it.
    if (strand not in ('+', '-') or not original.query_sequence
            or original.is_supplementary or original.has_tag('SA')
            or original.has_tag('Xb')):
        return None
    _, old_ns = _geometry(original)
    _, evolved_ns = _geometry(evolved)
    if len(old_ns) != len(evolved_ns) or not old_ns:
        return None
    reserved = {i for i, (a, b) in enumerate(zip(old_ns, evolved_ns)) if a[1:] != b[1:]}
    trial = evolved.__copy__()
    chrom = standardize_chrom_name(original.reference_name)
    entries = junctions_idx.get(chrom, [])
    # Never expand the existing 50-bp safety guard; charge every possible
    # endpoint, not just annotated/pool survivors. Union bound across read Ns.
    radius = max(0, min(50, max_boundary_shift, search_radius))
    if radius == 0:
        return None
    search_trials = (2 * radius + 1) * len(old_ns)
    accepted = []
    for ordinal in range(len(old_ns)):
        block = _gapfree_block(trial, ordinal, strand)
        if block is None:
            continue
        lo, hi, q0, q1, r0, r1, partner, introns = block
        touched = (ordinal,) if partner is None else tuple(sorted((ordinal, partner)))
        if any(i in reserved for i in touched):
            continue
        query = trial.query_sequence[q0:q1].upper()
        if not query or any(b not in 'ACGT' for b in query) or not (0 <= r0 < r1 <= len(genome_seq)):
            continue
        old_reference = genome_seq[r0:r1].upper()
        if query == old_reference:
            continue
        ni, ns, ne = introns[ordinal]
        candidates = _candidates_near(
            junctions_idx, chrom, ns, ne, radius,
            start_radius=radius, end_radius=radius,
            max_junction_size=max_junction_size,
        )
        # A capped family cannot establish uniqueness; do not silently choose a
        # winner from the subset allowed to the legacy scorer.
        if max_candidates_per_nop and len(candidates) > max_candidates_per_nop:
            continue
        candidates = {(js, je) for js, je in candidates
                      if (js == ns if strand == '+' else je == ne)
                      and (js, je) != (ns, ne)}
        if not candidates:
            continue
        information = effective_information_bits(query)
        if search_trials * 2.0 ** (-information) > DEFAULT_ALPHA:
            continue
        period = min_self_match_period(query)
        # Uniqueness is over the complete endpoint window, including sequence
        # placements absent from the junction pool. Annotation cannot hide an
        # equally exact repeat and create a spurious positional preference.
        window_start = max(0, r0 - radius)
        window = genome_seq[window_start:min(len(genome_seq), r1 + radius)].upper()
        first = window.find(query)
        if first < 0 or window.find(query, first + 1) >= 0:
            continue
        unique_start = window_start + first
        exact = []
        for js, je in candidates:
            if (strand == '+' and js != ns) or (strand == '-' and je != ne):
                continue
            delta = je - ne if strand == '+' else js - ns
            if not delta or (period is not None and abs(delta) >= period):
                continue
            if not (0 <= r0 + delta < r1 + delta <= len(genome_seq)):
                continue
            if r0 + delta != unique_start:
                continue
            # Count sequence ties before motif/annotation gates, so those
            # priors cannot manufacture unique sequence placement evidence.
            exact.append((js, je, delta))
        if len(exact) != 1:
            continue
        js, je, delta = exact[0]
        changes = [(ordinal, (ns, ne), (js, je))]
        if partner is not None:
            _, ps, pe = introns[partner]
            new_partner = (ps + delta, pe) if strand == '+' else (ps, pe + delta)
            changes.append((partner, (ps, pe), new_partner))
        permitted = True
        for idx, old_pair, new_pair in changes:
            s, e = new_pair
            if (not _canonical(genome_seq, s, e, strand)
                    or not _supported(entries, new_pair)
                    or (max_junction_size is not None and e - s > max_junction_size)):
                permitted = False
                break
            # Strong initial hold: annotated canonical (including repo-specific
            # tiered acceptors) may only move to another annotated canonical.
            if ((chrom, *old_pair) in annotated
                    and _canonical_tier(*old_pair, genome_seq, strand) < 4
                    and (chrom, *new_pair) not in annotated):
                permitted = False
                break
        if not permitted:
            continue
        new_ops = list(trial.cigartuples)
        for idx, old_pair, new_pair in changes:
            new_ops[introns[idx][0]] = (3, new_pair[1] - new_pair[0])
        # '=' and 'X' op semantics describe the OLD placement. All moved bases
        # are exact at the new placement; use M without changing op indices.
        for idx in range(lo, hi):
            new_ops[idx] = (0, new_ops[idx][1])
        trial.cigartuples = new_ops
        if strand == '-' and partner is None:
            trial.reference_start += delta
        # Exact gap-free block supplies real anchors at both ends. They are
        # not independent evidence and are never multiplied into the score.
        old_matches = sum(a == b for a, b in zip(query, old_reference))
        accepted.append(BlockPlacementEvidence(
            q0, q1, r0, r0 + delta, old_matches, len(query), information,
            search_trials, len(query), len(query), touched,
        ))
        reserved.update(touched)
    if not accepted:
        return None
    _, final_ns = _geometry(trial)
    return AtomicReadPlan(
        placement_state(original), strand, trial.reference_start, tuple(trial.cigartuples),
        tuple((a[1], a[2], b[1], b[2]) for a, b in zip(old_ns, final_ns)),
        tuple(accepted),
    )


def apply_atomic_plan(read, plan, genome_seq):
    """Fail closed on stale input; preserve SEQ/QUAL/tags without partial edit."""
    if read.is_supplementary or read.has_tag('SA') or read.has_tag('Xb'):
        return read, False
    decoded = decoded_copy(read, genome_seq)
    if placement_state(decoded) != plan.expected:
        return read, False
    # Also fail closed if the caller supplies a different reference at a moved
    # block. This verifies the realized sequence placement without rescoring.
    for evidence in plan.evidence:
        # Revalidate motifs only for boundaries owned by the new block path;
        # legacy moves included in the atomic snapshot retain their own policy.
        for ordinal in evidence.changed_junction_ordinals:
            _, _, start, end = plan.junction_map[ordinal]
            if not _canonical(genome_seq, start, end, plan.strand):
                return read, False
        query = decoded.query_sequence[evidence.query_start:evidence.query_end]
        if query.upper() != genome_seq[evidence.new_start:evidence.new_start + len(query)].upper():
            return read, False
    decoded.reference_start = plan.reference_start
    decoded.cigartuples = list(plan.cigar)
    # Alignment-derived scores are stale after true base relocation. Keep all
    # provenance/signal/array tags; remove individual tags without rebuilding.
    for tag in ('MD', 'cs', 'NM', 'AS', 'ms', 'de', 'dv', 'UQ'):
        if decoded.has_tag(tag):
            decoded.set_tag(tag, None)
    return decoded, True
