"""Information-bounded splice-overhang resolver — the mapPacBio-role rewrite.

planning/641: mapPacBio contributes ZERO unique reads to the panel
(planning/583 §4a) — its actual role is *alternative placement of terminal
overhangs across splice junctions*. This module does that deliberately, with
an explicit false-discovery budget, instead of running a general-purpose
aligner with a mammalian constant (``maxindel=200000`` == ~12,500 chance AG
acceptors per window on both strands — simultaneously the compute bill and
the false-positive rate).

Per terminal soft clip:

1. **Assess** the clip with the shared informativeness gate
   (:mod:`rectify.core.splice.overhang_informativeness`): ``I_eff`` ->
   ``W_max = alpha * 2**I_eff``. A refused clip (poly(A) tail, repeat
   expansion, any low-information overhang) is emitted UNCHANGED — refusing
   to search is a first-class correct outcome, and no candidate is ever
   evaluated for it (641 acceptance test T1 asserts on the counter).
2. **Look up** candidate splice sites in the precomputed
   :class:`~rectify.core.splice.splice_site_index.SpliceSiteIndex` — a
   binary-search range query bounded by ``W_max``, never a scan.
3. **Score** each candidate placement with the HP-aware global edit
   distance using the exact ``>``-only early-exit cutoff (planning/596 T2).
4. **Accept** only an unambiguous winner: ``ed <= max_edit_frac * L`` and
   every competitor within ``min_margin`` is the SAME junction after
   ambiguity canonicalisation (planning/621) — never a fixed +/-k window.

Geometry (0-based half-open introns ``[start, end)``; BAM SEQ is stored in
forward-genome orientation, so clip bases compare literally — no revcomp):

=======  ==========  =================  ==================  ===============
clip     transcript  near site          far site (in W)     exon match
side     strand      (within slop of                        segment
                     the aligned edge)
=======  ==========  =================  ==================  ===============
LEFT     ``+``       acc_plus  ``e<=R``  don_plus ``f``     ``[f-m, f)``
LEFT     ``-``       don_minus ``e<=R``  acc_minus ``f``    ``[f-m, f)``
RIGHT    ``+``       don_plus  ``d>=E``  acc_plus ``e``     ``[e, e+m)``
RIGHT    ``-``       acc_minus ``d>=E``  don_minus ``e``    ``[e, e+m)``
=======  ==========  =================  ==================  ===============

v1 limits (documented, deliberate): the near site must sit at-or-outside the
aligned edge (aligned-into-the-intron cases remain ``rescue_3ss_truncation``'s
reanchor territory); one intron per clip; clips are matched on their
junction-proximal ``max_clip_match`` bases, the remainder stays soft-clipped.

The resolver consumes the (name-sorted) minimap2 arm BAM and emits a BAM in
the same order, so it drops into the panel exactly where mapPacBio's arm
did. Accepted placements carry ``XJ:Z:<intron_start>-<intron_end>:<ed>:<side>``;
MD/NM are dropped on rewritten records (stale after CIGAR surgery).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from ..splice.overhang_informativeness import (
    COUNTERS,
    assess_overhang,
    canonical_in_class,
    hp_edit_distance_bounded,
    same_junction,
)
from ..splice.repeat_expansion import is_repeat_expansion
from ..splice.splice_site_index import SpliceSiteIndex
from ...config import CHROM_TO_GENOME
from ...utils.genome import load_genome, standardize_chrom_name

logger = logging.getLogger(__name__)

_LEFT = 'L'
_RIGHT = 'R'


@dataclass
class ResolverConfig:
    alpha: float = 0.01          # false-discovery budget (641 T7: a published knob)
    max_intron: int = 5000       # per-reference intron cap (clamps W_max)
    min_intron: int = 40         # matches BBMap intronlen=40 D->N semantics
    min_clip: int = 8            # clips shorter than this are never assessed
    max_clip_match: int = 200    # junction-proximal bases used for matching
    edge_slop: int = 5           # near site within this many bp of the aligned edge
    max_edit_frac: float = 0.2   # accept threshold on best ED / matched length
    min_margin: float = 1.0      # runner-up (different junction) must be this much worse
    min_far_match: int = 8       # minimum exon bases beyond the far site
    # Near sites up to this many bp INSIDE the aligned block are also
    # candidates, re-assigning those aligned bases across the junction.
    # planning/644 T4c: minimap2 routinely mis-anchors the aligned start
    # 10-18 bp past a true acceptor into the intron tail (local homology),
    # which both hides the true near site and lets a chance site win
    # uncontested. Letting inside-edge sites compete fixes both.
    edge_inside_slop: int = 24

    # --- v2 junction re-arbitration (planning/644b) --------------------------
    # The peelback move generalized: treat the aligner's junction ASSIGNMENT
    # (and suspicious linear stretches) as hypotheses and re-score competing
    # spliced interpretations from the splice-site index, under the same
    # information budget and margin discipline. Targets the T3 residual:
    # junctions minimap2 mis-assigned to a nearby boundary (alt-SS class) or
    # missed entirely in a linear alignment.
    arb_enable: bool = True
    # Canonical-preference tiebreak at equal/better ED when the CURRENT
    # junction is non-canonical-class (the SRC1 adjudicator). Measured cost
    # (realigner noncanon control, 2026-08-09): on genuine non-canonical
    # junctions with a canonical decoy in range, this snap flattens real
    # signal (-6.0 pp cryptic recovery on the smoke mixture). Leave ON for
    # annotated/SMD accuracy work; turn OFF for non-canonical discovery
    # missions (upf1d/prp18d-type). Grammar-driven moves are marked ':g' in
    # the XB tag and counted as arb_grammar_tiebreak either way.
    arb_grammar: bool = True
    arb_window: int = 300        # max boundary shift searched (also clamped by W_max)
    arb_seg: int = 40            # query bases per junction side used for scoring
    arb_margin: float = 2.0      # an alternative must beat the current placement by this
    arb_dop_min: int = 20        # D op at least this long = putative intron (Case B)
    arb_dop_slop: int = 12       # boundary snap window around a D op
    arb_mm_win: int = 50         # mismatch-cluster trigger: window (query bases)
    arb_mm_frac: float = 0.30    # ...and the mismatch fraction that flags it


@dataclass
class ResolverStats:
    reads: int = 0
    passthrough_nonprimary: int = 0
    clips_seen: int = 0
    clips_assessed: int = 0
    refused_low_info: int = 0
    refused_repeat: int = 0
    no_candidates: int = 0
    rejected_edit: int = 0
    rejected_ambiguous: int = 0
    resolved: int = 0
    resolved_left: int = 0
    resolved_right: int = 0
    candidates_evaluated: int = 0
    extra: Dict[str, int] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, int]:
        d = {k: v for k, v in self.__dict__.items() if k != 'extra'}
        d.update(self.extra)
        return d


@dataclass
class _Placement:
    ed: float
    intron_start: int
    intron_end: int
    lead: int            # exon bases between the near site and the aligned edge
    m: int               # exon bases matched beyond the far site
    canonical_rank: int  # 0 = GT/AG-class donor, 1 = GC-class
    k_inside: int = 0    # aligned bases re-assigned across the junction
                         # (inside-edge near site; planning/644 T4c)


def _clip_lens(cigartuples) -> Tuple[int, int]:
    """(left_softclip, right_softclip) lengths; 0 when absent."""
    if not cigartuples:
        return 0, 0
    left = cigartuples[0][1] if cigartuples[0][0] == 4 else 0
    right = cigartuples[-1][1] if cigartuples[-1][0] == 4 else 0
    return left, right


def _site_kinds(side: str, strand: str) -> Tuple[str, str]:
    """(near_kind, far_kind) for a clip side x transcript strand."""
    if side == _LEFT:
        return ('acc_plus', 'don_plus') if strand == '+' else ('don_minus', 'acc_minus')
    return ('don_plus', 'acc_plus') if strand == '+' else ('acc_minus', 'don_minus')


def _donor_rank(chrom_seq: str, side: str, strand: str,
                intron_start: int, intron_end: int) -> int:
    """0 for the canonical GT-class donor, 1 for GC-class."""
    if strand == '+':
        dinuc = chrom_seq[intron_start:intron_start + 2].upper()
        return 0 if dinuc == 'GT' else 1
    dinuc = chrom_seq[intron_end - 2:intron_end].upper()
    return 0 if dinuc == 'AC' else 1


def resolve_clip(
    chrom_seq: str,
    index: SpliceSiteIndex,
    chrom_key: str,
    side: str,
    strand: str,
    clip_seq: str,
    edge: int,
    cfg: ResolverConfig,
    stats: Optional[ResolverStats] = None,
    inside_seq: str = '',
) -> Optional[_Placement]:
    """Resolve one terminal soft clip against the splice-site index.

    ``edge`` is the aligned boundary the clip attaches to: ``reference_start``
    for a LEFT clip, ``reference_end`` (exclusive) for a RIGHT clip. Returns
    an accepted :class:`_Placement` or None (refusal / no unambiguous winner).

    ``inside_seq``: the aligned query bases immediately adjacent to the clip
    (already '='-decoded by the caller), enabling inside-edge near sites —
    the mis-anchored-edge geometry from planning/644 T4c. Empty disables.
    """
    if stats is None:
        stats = ResolverStats()
    if len(clip_seq) < cfg.min_clip:
        return None
    stats.clips_seen += 1

    if '=' in clip_seq:
        # calmd never encodes soft-clip bases as '='; seeing one means the
        # query sequence is corrupt. Loud failure, per planning/638 §3.
        raise ValueError(
            "overhang_resolver: '=' byte inside a soft clip — query sequence "
            "is corrupt upstream (see the planning/581 eq-encoding note)"
        )

    # Junction-proximal portion: LEFT clips attach at their right end,
    # RIGHT clips at their left end.
    if len(clip_seq) > cfg.max_clip_match:
        clip_used = clip_seq[-cfg.max_clip_match:] if side == _LEFT else clip_seq[:cfg.max_clip_match]
    else:
        clip_used = clip_seq
    lc = len(clip_used)

    # --- The gate: refuse before any candidate is touched -----------------
    if is_repeat_expansion(clip_used):
        stats.refused_repeat += 1
        return None
    stats.clips_assessed += 1
    assessment = assess_overhang(clip_used, alpha=cfg.alpha, max_window=cfg.max_intron)
    if assessment.refused:
        stats.refused_low_info += 1
        return None
    w = assessment.w_max_bp

    near_kind, far_kind = _site_kinds(side, strand)

    # --- Candidate lookup (binary search, bounded by W) --------------------
    placements: List[_Placement] = []
    best_ed = float('inf')

    k_cap = min(cfg.edge_inside_slop, len(inside_seq))

    def _cutoff() -> float:
        # Exact pruning bound: a candidate with ed > max_edit_frac*len(cmp)
        # can never be accepted, and one with ed >= best + margin can neither
        # win nor block on ambiguity — so anything above
        # min(best, threshold) + margin is irrelevant. The threshold uses the
        # LONGEST comparison length (lc + k_cap) so the bound is valid for
        # every candidate; values AT the bound are still computed exactly
        # (the DP prunes on > only).
        return min(best_ed, cfg.max_edit_frac * (lc + k_cap)) + cfg.min_margin
    if side == _LEFT:
        # near site e in [edge - slop, edge + k_cap]; far f in [e - w, e - min_intron]
        near_sites = index.sites_in(chrom_key, near_kind,
                                    edge - cfg.edge_slop, edge + k_cap + 1)
        for e in near_sites:
            e = int(e)
            if e <= edge:
                lead, k = edge - e, 0
                cmp_seq = clip_used
            else:
                lead, k = 0, e - edge
                cmp_seq = clip_used + inside_seq[:k]
            m = len(cmp_seq) - lead
            if m < cfg.min_far_match:
                continue
            far_lo = max(0, e - w)
            far_hi = e - cfg.min_intron + 1
            for f in index.sites_in(chrom_key, far_kind, far_lo, far_hi):
                f = int(f)
                if f - m < 0:
                    continue
                ref = chrom_seq[f - m:f] + (chrom_seq[e:edge] if k == 0 else '')
                stats.candidates_evaluated += 1
                COUNTERS['candidates_evaluated'] += 1
                _c = _cutoff()
                ed = hp_edit_distance_bounded(cmp_seq, ref, cutoff=_c)
                if ed <= _c:
                    placements.append(_Placement(
                        ed=ed, intron_start=f, intron_end=e, lead=lead, m=m,
                        canonical_rank=_donor_rank(chrom_seq, side, strand, f, e),
                        k_inside=k,
                    ))
                    best_ed = min(best_ed, ed)
    else:
        # near site d in [edge - k_cap, edge + slop]; far e in [d + min_intron, d + w]
        near_sites = index.sites_in(chrom_key, near_kind,
                                    edge - k_cap, edge + cfg.edge_slop + 1)
        for d in near_sites:
            d = int(d)
            if d >= edge:
                lead, k = d - edge, 0
                cmp_seq = clip_used
            else:
                lead, k = 0, edge - d
                cmp_seq = inside_seq[len(inside_seq) - k:] + clip_used
            m = len(cmp_seq) - lead
            if m < cfg.min_far_match:
                continue
            far_lo = d + cfg.min_intron
            far_hi = d + w + 1
            for e in index.sites_in(chrom_key, far_kind, far_lo, far_hi):
                e = int(e)
                if e + m > len(chrom_seq):
                    continue
                ref = (chrom_seq[edge:d] if k == 0 else '') + chrom_seq[e:e + m]
                stats.candidates_evaluated += 1
                COUNTERS['candidates_evaluated'] += 1
                _c = _cutoff()
                ed = hp_edit_distance_bounded(cmp_seq, ref, cutoff=_c)
                if ed <= _c:
                    placements.append(_Placement(
                        ed=ed, intron_start=d, intron_end=e, lead=lead, m=m,
                        canonical_rank=_donor_rank(chrom_seq, side, strand, d, e),
                        k_inside=k,
                    ))
                    best_ed = min(best_ed, ed)

    if not placements:
        stats.no_candidates += 1
        return None

    # --- Selection: unambiguous winner only --------------------------------
    placements.sort(key=lambda p: (p.ed, p.canonical_rank, p.intron_start))
    best = placements[0]
    # comparison length = m + lead for k==0 (== lc) and m for k>0 (== lc + k);
    # both equal p.m + p.lead.
    if best.ed > cfg.max_edit_frac * (best.m + best.lead):
        stats.rejected_edit += 1
        return None
    for other in placements[1:]:
        if other.ed - best.ed >= cfg.min_margin:
            break  # sorted: everything after is at least as far
        if not same_junction(
            chrom_seq,
            (best.intron_start, best.intron_end),
            (other.intron_start, other.intron_end),
        ):
            stats.rejected_ambiguous += 1
            return None
    return best


def _rewrite_cigar(
    read: pysam.AlignedSegment,
    side: str,
    placement: _Placement,
    clip_len: int,
    clip_used_len: int,
) -> None:
    """Apply an accepted placement to ``read`` in place (CIGAR + position).

    The matched portion of the clip becomes ``M / N / M``; any clip bases
    beyond ``max_clip_match`` remain soft-clipped. MD/NM are dropped (stale).
    """
    ct = list(read.cigartuples)
    remainder = clip_len - clip_used_len
    intron_len = placement.intron_end - placement.intron_start
    k = placement.k_inside

    if side == _LEFT:
        new_ops: List[Tuple[int, int]] = []
        if remainder > 0:
            new_ops.append((4, remainder))
        new_ops.append((0, placement.m))
        new_ops.append((3, intron_len))
        rest = ct[1:]
        if k > 0:
            # k aligned bases moved across the junction into the exon-1 M;
            # the caller guarantees rest[0] is an M op longer than k.
            rest = [(0, rest[0][1] - k)] + rest[1:]
        if placement.lead > 0:
            if rest and rest[0][0] == 0:
                rest = [(0, placement.lead + rest[0][1])] + rest[1:]
            else:
                new_ops.append((0, placement.lead))
        read.cigartuples = new_ops + rest
        read.reference_start = placement.intron_start - placement.m
    else:
        rest = ct[:-1]
        if k > 0:
            rest = rest[:-1] + [(0, rest[-1][1] - k)]
        new_ops = []
        if placement.lead > 0:
            if rest and rest[-1][0] == 0:
                rest = rest[:-1] + [(0, rest[-1][1] + placement.lead)]
            else:
                new_ops.append((0, placement.lead))
        new_ops.append((3, intron_len))
        new_ops.append((0, placement.m))
        if remainder > 0:
            new_ops.append((4, remainder))
        read.cigartuples = rest + new_ops

    for tag in ('MD', 'NM'):
        if read.has_tag(tag):
            read.set_tag(tag, None)
    read.set_tag(
        'XJ',
        f'{placement.intron_start}-{placement.intron_end}:{placement.ed:.1f}:{side}',
    )


# ---------------------------------------------------------------------------
# v2 junction re-arbitration (planning/644b) — the peelback move generalized
# to fully-aligned reads and internal junctions.
# ---------------------------------------------------------------------------

def _bump(stats: ResolverStats, key: str, n: int = 1) -> None:
    stats.extra[key] = stats.extra.get(key, 0) + n


def _boundary_kinds(strand: str) -> Tuple[str, str]:
    """(left_kind, right_kind) of an intron's genomic boundaries."""
    if strand == '+':
        return 'don_plus', 'acc_plus'
    return 'acc_minus', 'don_minus'


def _decoded_query(read: pysam.AlignedSegment, chrom_seq: str) -> str:
    """Query sequence with calmd '=' bytes decoded from the reference.
    Soft clips and insertions are returned verbatim (never '='-encoded)."""
    q = read.query_sequence
    if '=' not in q:
        return q
    out = []
    qi = 0
    ref = read.reference_start
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):  # M/=/X consume both
            seg = q[qi:qi + ln]
            if '=' in seg:
                seg = ''.join(
                    chrom_seq[ref + i] if b == '=' else b
                    for i, b in enumerate(seg))
            out.append(seg)
            qi += ln
            ref += ln
        elif op in (1, 4):   # I/S consume query only
            out.append(q[qi:qi + ln])
            qi += ln
        elif op in (2, 3):   # D/N consume reference only
            ref += ln
    return ''.join(out)


def _iter_ref_ops(cigartuples):
    """Yield (op_index, op, length, ref_start, q_start) walking the CIGAR."""
    ref = 0
    qi = 0
    for i, (op, ln) in enumerate(cigartuples):
        yield i, op, ln, ref, qi
        if op in (0, 2, 3, 7, 8):
            ref += ln
        if op in (0, 1, 4, 7, 8):
            qi += ln


def _score_alts(qwin, cur_ref, alts, chrom_seq, cfg, cur_junc=None, ed_cur=None):
    """Score the current interpretation exactly, then alternatives under an
    exact cutoff. ``alts`` entries are (d, e, alt_ref, alt_qwin): a
    LEFT-boundary shift by delta moves the junction's query split by delta
    too, so each alternative is scored against its own query window (equal
    length; same locality).

    Acceptance: an alternative must beat the current placement by
    ``arb_margin`` — OR, when ``cfg.arb_grammar`` is on, the CURRENT
    junction is non-canonical-class and the alternative is canonical,
    match-or-beat it (splicing grammar as the tiebreaker, mirroring the
    rescue path's canonical-donor tiebreak; planning/644b: the SRC1 4-bp
    donor dispute ties the DP and only the GT/GC..AG grammar can adjudicate
    it). Returns (ed_cur, winner, via_grammar) with winner = (ed, d, e) or
    None; via_grammar is True when the winner was admitted ONLY by the
    grammar tiebreak (it did not clear the margin bound)."""
    if ed_cur is None:
        ed_cur = hp_edit_distance_bounded(qwin, cur_ref)
    cur_canon = (canonical_in_class(chrom_seq, *cur_junc)
                 if cur_junc is not None else True)
    bound = ed_cur - cfg.arb_margin
    grammar_bound = ed_cur if (cfg.arb_grammar and not cur_canon) else -1.0
    prune_at = max(bound, grammar_bound)
    if prune_at < 0:
        return ed_cur, None, False
    scored = []
    for d_alt, e_alt, alt_ref, alt_qwin in alts:
        ed = hp_edit_distance_bounded(alt_qwin, alt_ref, cutoff=prune_at)
        if ed > cfg.max_edit_frac * len(alt_qwin):
            continue
        if ed <= bound or (
                ed <= grammar_bound
                and canonical_in_class(chrom_seq, d_alt, e_alt)):
            scored.append((ed, d_alt, e_alt))
    if not scored:
        return ed_cur, None, False
    scored.sort()
    best = scored[0]
    for other in scored[1:]:
        if other[0] - best[0] >= cfg.min_margin:
            break
        if not same_junction(chrom_seq, (best[1], best[2]), (other[1], other[2])):
            return ed_cur, None, False   # ambiguous among alternatives
    return ed_cur, best, best[0] > bound


def _rearbitrate_read(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    index: SpliceSiteIndex,
    chrom_key: str,
    strand: str,
    cfg: ResolverConfig,
    stats: ResolverStats,
) -> bool:
    """Case A: re-arbitrate the FIRST/LAST N-op's boundaries against nearby
    index sites (mis-ASSIGNED junctions, the alt-SS class). Case B: convert
    intron-length D ops to index-snapped N (MISSED junction expressed as a
    deletion), and splice mismatch-flagged terminal linear blocks. All under
    the information budget + margin discipline; query length is conserved by
    every rewrite."""
    ct = list(read.cigartuples)
    q = read.query_sequence
    if not ct or not q:
        return False
    dq = None  # decoded lazily, once
    left_kind, right_kind = _boundary_kinds(strand)
    changed = False

    # ---- Case A0: boundary-deletion merge (planning/644c, the SRC1 smoking
    # gun). A D op abutting an N is the aligner encoding "these reference
    # bases have no query support" while asserting the wrong boundary —
    # minimap2 wrote M D4 N at SRC1, i.e. the GC-donor evidence AND the GT
    # claim in one CIGAR. Merging the D into the intron is ALIGNMENT-
    # IDENTICAL (same query, same reference consumption), strictly more
    # parsimonious, and accepted only when the merged junction is
    # canonical-in-class. The boundary-adjacent generalization of BBMap's
    # intronlen D->N semantics.
    merged = True
    while merged:
        merged = False
        for i, (op, ln) in enumerate(ct):
            if op != 3:
                continue
            walk = {j: (o, l, rs, qs) for j, o, l, rs, qs in _iter_ref_ops(ct)}
            _, nlen, n_rs, _ = walk[i]
            d = read.reference_start + n_rs
            e = d + nlen
            if i >= 1 and ct[i - 1][0] == 2:
                dl = ct[i - 1][1]
                if (e - (d - dl)) <= cfg.max_intron and \
                        canonical_in_class(chrom_seq, d - dl, e):
                    ct = ct[:i - 1] + [(3, nlen + dl)] + ct[i + 1:]
                    _bump(stats, 'arb_dmerge')
                    merged = True
                    changed = True
                    break
            if i + 1 < len(ct) and ct[i + 1][0] == 2:
                dl = ct[i + 1][1]
                if ((e + dl) - d) <= cfg.max_intron and \
                        canonical_in_class(chrom_seq, d, e + dl):
                    ct = ct[:i] + [(3, nlen + dl)] + ct[i + 2:]
                    _bump(stats, 'arb_dmerge')
                    merged = True
                    changed = True
                    break
    if changed:
        read.cigartuples = ct
        for tag in ('MD', 'NM'):
            if read.has_tag(tag):
                read.set_tag(tag, None)
        read.set_tag('XB', 'dmerge')

    n_idx = [i for i, (op, _) in enumerate(ct) if op == 3]

    # ---- Case A: boundary shifts on EVERY N op ----------------------------
    # planning/644 stage-1 vicinity triggers (Kevin, 2026-08-09): the
    # junction-proximal ed_cur below is a local mismatch/indel density meter,
    # and any junction above the margin enters arbitration — including
    # INTERIOR junctions. Frame safety per family: a left-boundary (genomic-
    # donor-side) shift changes only the flanking M swap and the N length
    # with the downstream anchor FIXED at e, so it is safe on ANY N;
    # right-boundary and diagonal shifts move everything downstream, so they
    # are offered only on the LAST N. All families for one junction compete
    # in ONE round (sequential greedy rounds lock in local optima).
    # Rewrites never change the op count, so indexes stay valid.
    targets = [('both' if i == n_idx[-1] else 'left_only', i) for i in n_idx]
    for which, i in targets:
        if i == 0 or i == len(ct) - 1:
            continue
        walk = {j: (op, ln, rs, qs) for j, op, ln, rs, qs in _iter_ref_ops(ct)}
        _, nlen, n_rs, _ = walk[i]
        d = read.reference_start + n_rs
        e = d + nlen
        # M flank lengths gate only the left/diagonal REWRITE bookkeeping;
        # SCORING windows come from the decoded query itself, so ONT indel
        # fragmentation at a mis-assigned boundary (planning/644c: the very
        # signature of mis-assignment) cannot exempt a junction from
        # arbitration. Windows stop at the soft clips.
        m_l = ct[i - 1][1] if ct[i - 1][0] in (0, 7, 8) else 0
        m_r = ct[i + 1][1] if ct[i + 1][0] in (0, 7, 8) else 0
        # Frame safety (v5.1 — the chrII 98b21cfd corruption, 2026-08-10).
        # (1) A d-moving rewrite swaps query between the FLANKING ops, so
        # both must be M-type: the m_l/m_r=0 sentinel let a delta<0 swap
        # overwrite a flanking 1I with a manufactured M (query-length
        # violation; htslib refuses the record — 11/17 chromosome BAMs).
        # (2) EVERY shift changes net reference consumption through the
        # junction, displacing all ref-consuming ops downstream of the
        # right flank (measured: a +2 donor shift moved a trailing
        # 4M4D55M block by +2; an interior-N shift displaces the NEXT
        # junction by -delta) — so shifts are offered only when nothing
        # beyond the right flank consumes reference. Interior junctions
        # keep the D-merge; their shift arbitration belongs to stage-2's
        # junction-local traceback realigner.
        flank_m = ct[i - 1][0] in (0, 7, 8) and ct[i + 1][0] in (0, 7, 8)
        if not all(o in (1, 4, 5) for o, _ in ct[i + 2:]):
            _bump(stats, 'arb_frame_unsafe_skip')
            continue
        _, _, _, q_split = walk[i + 1]      # query index where exon-2 resumes
        lclip = ct[0][1] if ct[0][0] == 4 else 0
        rclip = ct[-1][1] if ct[-1][0] == 4 else 0
        if dq is None:
            dq = _decoded_query(read, chrom_seq)
        L1 = min(cfg.arb_seg, q_split - lclip)
        L2 = min(cfg.arb_seg, len(dq) - rclip - q_split)
        if L1 < 10 or L2 < 10 or d - L1 < 0 or e + L2 > len(chrom_seq):
            continue
        qwin = dq[q_split - L1:q_split + L2]
        _bump(stats, 'arb_njunc_checked')
        # Fast path: a canonical-class junction whose window already scores
        # under the margin cannot be displaced (an alt would need ed < 0 and
        # the grammar tiebreak is off) — skip candidate enumeration entirely.
        # This is most junctions, and most of Case A's cost.
        cur_ref = chrom_seq[d - L1:d] + chrom_seq[e:e + L2]
        ed_cur = hp_edit_distance_bounded(qwin, cur_ref)
        if ed_cur < cfg.arb_margin and canonical_in_class(chrom_seq, d, e):
            _bump(stats, 'arb_clean_skip')
            continue
        a = assess_overhang(qwin, alpha=cfg.alpha, max_window=cfg.arb_window)
        if a.refused:
            _bump(stats, 'arb_refused')
            continue
        W = a.w_max_bp
        alts = []
        # right-boundary shifts (last N only): N length changes, nothing else
        if which == 'both':
            for e_alt in index.sites_in(chrom_key, right_kind, e - W, e + W + 1):
                e_alt = int(e_alt)
                if e_alt == e or not (cfg.min_intron <= e_alt - d <= cfg.max_intron):
                    continue
                if e_alt + L2 > len(chrom_seq):
                    continue
                alts.append((d, e_alt,
                             chrom_seq[d - L1:d] + chrom_seq[e_alt:e_alt + L2], qwin))
        # left-boundary shifts: query swaps between the flanking Ms (both
        # flanks must BE Ms — flank_m; see the frame-safety note above)
        if which in ('left_only', 'both') and flank_m:
            for d_alt in index.sites_in(chrom_key, left_kind, d - W, d + W + 1):
                d_alt = int(d_alt)
                delta = d_alt - d
                if delta == 0 or not (cfg.min_intron <= e - d_alt <= cfg.max_intron):
                    continue
                if m_l + delta < 1 or m_r - delta < 1 or d_alt - L1 < 0:
                    continue
                if q_split + delta - L1 < 0 or q_split + delta + L2 > len(dq):
                    continue
                alts.append((d_alt, e,
                             chrom_seq[d_alt - L1:d_alt] + chrom_seq[e:e + L2],
                             dq[q_split + delta - L1:q_split + delta + L2]))
        # diagonal (pseudo-slide) shifts: BOTH boundaries move by delta with
        # the intron length preserved — the beyond-legal-slide class where
        # partial homology lets two placements near-tie (the SRC1 4-bp donor
        # dispute). Needs the same M-flank swap as the left family.
        if which == 'both' and flank_m:
            span = e - d
            for d_alt in index.sites_in(chrom_key, left_kind, d - W, d + W + 1):
                d_alt = int(d_alt)
                delta = d_alt - d
                if delta == 0:
                    continue
                e_alt = d_alt + span
                if index.sites_in(chrom_key, right_kind, e_alt, e_alt + 1).size == 0:
                    continue
                if m_l + delta < 1 or m_r - delta < 1 or d_alt - L1 < 0 \
                        or e_alt + L2 > len(chrom_seq):
                    continue
                if q_split + delta - L1 < 0 or q_split + delta + L2 > len(dq):
                    continue
                alts.append((d_alt, e_alt,
                             chrom_seq[d_alt - L1:d_alt] + chrom_seq[e_alt:e_alt + L2],
                             dq[q_split + delta - L1:q_split + delta + L2]))
        if not alts:
            continue
        ed_cur, win, via_grammar = _score_alts(qwin, cur_ref, alts, chrom_seq,
                                               cfg, cur_junc=(d, e),
                                               ed_cur=ed_cur)
        if win is None:
            _bump(stats, 'arb_no_gain')
            continue
        _, d_new, e_new = win
        if d_new != d:                       # left shift: swap query between Ms
            delta = d_new - d
            ct[i - 1] = (0, m_l + delta)
            ct[i + 1] = (0, m_r - delta)
        ct[i] = (3, e_new - d_new)
        read.cigartuples = ct
        for tag in ('MD', 'NM'):
            if read.has_tag(tag):
                read.set_tag(tag, None)
        gmark = ':g' if via_grammar else ''
        read.set_tag('XB', f'shift:{d}-{e}>{d_new}-{e_new}'
                           f':{ed_cur:.1f}>{win[0]:.1f}{gmark}')
        _bump(stats, 'arb_shifted')
        if via_grammar:
            _bump(stats, 'arb_grammar_tiebreak')
        changed = True
        ct = list(read.cigartuples)
        n_idx = [j for j, (op, _) in enumerate(ct) if op == 3]

    # ---- Case B1: intron-length D ops -> index-snapped N ------------------
    for i, (op, ln) in enumerate(list(ct)):
        if op != 2 or ln < cfg.arb_dop_min:
            continue
        if i == 0 or i == len(ct) - 1 or ct[i - 1][0] not in (0, 7, 8) or ct[i + 1][0] not in (0, 7, 8):
            continue
        walk = {j: (o, l, rs, qs) for j, o, l, rs, qs in _iter_ref_ops(ct)}
        _, dlen, d_rs, _ = walk[i]
        p_start = read.reference_start + d_rs
        p_end = p_start + dlen
        m_l = ct[i - 1][1]
        m_r = ct[i + 1][1]
        _, _, _, q_split = walk[i + 1]
        L1 = min(cfg.arb_seg, m_l)
        L2 = min(cfg.arb_seg, m_r)
        if L1 < 10 or L2 < 10 or p_start - L1 < 0 or p_end + L2 > len(chrom_seq):
            continue
        if dq is None:
            dq = _decoded_query(read, chrom_seq)
        qwin = dq[q_split - L1:q_split + L2]
        _bump(stats, 'arb_dop_checked')
        s2 = cfg.arb_dop_slop
        # frame safety: an asymmetric snap (delta_e != delta_d) changes net
        # reference consumption and displaces every ref-consuming op beyond
        # the right flank — allow it only when there are none (v5.1)
        b1_tail_clear = all(o in (1, 4, 5) for o, _ in ct[i + 2:])
        alts = []
        for d_alt in index.sites_in(chrom_key, left_kind, p_start - s2, p_start + s2 + 1):
            d_alt = int(d_alt)
            delta = d_alt - p_start
            if m_l + delta < 1 or m_r - delta < 1 or d_alt - L1 < 0:
                continue
            if q_split + delta - L1 < 0 or q_split + delta + L2 > len(dq):
                continue
            for e_alt in index.sites_in(chrom_key, right_kind, p_end - s2, p_end + s2 + 1):
                e_alt = int(e_alt)
                if not (cfg.min_intron <= e_alt - d_alt <= cfg.max_intron):
                    continue
                if e_alt + L2 > len(chrom_seq):
                    continue
                if not b1_tail_clear and (e_alt - p_end) != (d_alt - p_start):
                    continue
                alts.append((d_alt, e_alt,
                             chrom_seq[d_alt - L1:d_alt] + chrom_seq[e_alt:e_alt + L2],
                             dq[q_split + delta - L1:q_split + delta + L2]))
        if not alts:
            continue
        cur_ref = chrom_seq[p_start - L1:p_start] + chrom_seq[p_end:p_end + L2]
        # A deletion and an intron at identical bounds are the same spliced
        # product, so no margin is demanded — canonical-class snapping within
        # the slop is the win; ties go to the index site.
        ed_cur = hp_edit_distance_bounded(qwin, cur_ref)
        scored = sorted(
            (hp_edit_distance_bounded(aq, ref, cutoff=ed_cur), da, ea)
            for da, ea, ref, aq in alts)
        best = scored[0]
        if best[0] > ed_cur or best[0] > cfg.max_edit_frac * len(qwin):
            _bump(stats, 'arb_no_gain')
            continue
        _, d_new, e_new = best
        delta = d_new - p_start
        if delta:
            ct[i - 1] = (0, m_l + delta)
            ct[i + 1] = (0, m_r - delta)
        ct[i] = (3, e_new - d_new)
        read.cigartuples = ct
        for tag in ('MD', 'NM'):
            if read.has_tag(tag):
                read.set_tag(tag, None)
        read.set_tag('XB', f'dop:{p_start}-{p_end}>{d_new}-{e_new}')
        _bump(stats, 'arb_dop_spliced')
        changed = True
        ct = list(read.cigartuples)

    # ---- Case B2: mismatch-flagged TERMINAL linear block -> spliced -------
    # Right-side hypothesis on the read's last M block (no N after it): a
    # spliced molecule aligned linearly shows a mismatch cluster downstream
    # of the missed donor. Left side is the mirror; v2 implements the right
    # side (the common missed-3'-portion case) and counts what it flags.
    last_m = None
    for j, (op, ln) in enumerate(ct):
        if op in (0, 7, 8):
            last_m = j
    if last_m is not None and ct[last_m][1] >= 3 * cfg.arb_mm_win \
            and not any(op == 3 for op, _ in ct[last_m:]) \
            and all(op in (1, 4, 5) for op, _ in ct[last_m + 1:]):
        walk = {j: (o, l, rs, qs) for j, o, l, rs, qs in _iter_ref_ops(ct)}
        _, mlen, m_rs, m_qs = walk[last_m]
        block_ref = read.reference_start + m_rs
        if dq is None:
            dq = _decoded_query(read, chrom_seq)
        # coarse mismatch scan. On calmd '='-encoded input (the production
        # case) a mismatch is simply a non-'=' byte in the M block, countable
        # at C speed; otherwise fall back to a per-base compare.
        eq_mode = '=' in q
        onset = None
        for off in range(cfg.arb_mm_win, mlen - cfg.arb_mm_win, cfg.arb_mm_win // 2):
            if eq_mode:
                seg = q[m_qs + off:m_qs + off + cfg.arb_mm_win]
                mm = len(seg) - seg.count('=')
            else:
                seg = dq[m_qs + off:m_qs + off + cfg.arb_mm_win]
                ref = chrom_seq[block_ref + off:block_ref + off + cfg.arb_mm_win]
                mm = sum(1 for a_, b_ in zip(seg, ref) if a_ != b_)
            if mm >= cfg.arb_mm_frac * cfg.arb_mm_win:
                onset = off
                break
        if onset is not None:
            _bump(stats, 'arb_mm_flagged')
            L1 = cfg.arb_seg
            # The true donor sits where matching stops — anywhere within the
            # flagged window (mismatches begin mid-window), so search a span
            # covering the window plus slack on both sides.
            donors = index.sites_in(
                chrom_key, left_kind,
                block_ref + max(0, onset - 20),
                block_ref + onset + cfg.arb_mm_win + 20)
            best_overall = None
            ed_cur_o = None
            for d_alt in donors:
                d_alt = int(d_alt)
                o = d_alt - block_ref
                if o < L1 or mlen - o < 10:
                    continue
                q_j = m_qs + o
                L2 = min(cfg.arb_seg, mlen - o)
                qwin = dq[q_j - L1:q_j + L2]
                a = assess_overhang(dq[q_j:q_j + L2], alpha=cfg.alpha,
                                    max_window=cfg.max_intron)
                if a.refused:
                    _bump(stats, 'arb_refused')
                    continue
                W = a.w_max_bp
                alts = []
                for e_alt in index.sites_in(chrom_key, right_kind,
                                            d_alt + cfg.min_intron, d_alt + W + 1):
                    e_alt = int(e_alt)
                    if e_alt + L2 > len(chrom_seq):
                        continue
                    alts.append((d_alt, e_alt,
                                 chrom_seq[d_alt - L1:d_alt] + chrom_seq[e_alt:e_alt + L2],
                                 qwin))
                if not alts:
                    continue
                cur_ref = chrom_seq[d_alt - L1:d_alt + L2]
                ed_cur, win, _ = _score_alts(qwin, cur_ref, alts, chrom_seq, cfg)
                if win is not None and (best_overall is None or win[0] < best_overall[0]):
                    best_overall = win
                    ed_cur_o = ed_cur
                    best_o = o
            if best_overall is not None:
                _, d_new, e_new = best_overall
                new_ct = ct[:last_m] + [(0, best_o), (3, e_new - d_new),
                                        (0, mlen - best_o)] + ct[last_m + 1:]
                read.cigartuples = new_ct
                for tag in ('MD', 'NM'):
                    if read.has_tag(tag):
                        read.set_tag(tag, None)
                read.set_tag('XB', f'mm:{d_new}-{e_new}:{ed_cur_o:.1f}>{best_overall[0]:.1f}')
                _bump(stats, 'arb_mm_spliced')
                changed = True
                ct = list(read.cigartuples)

    # ---- Case B3: mismatch storm at the HEAD of the first aligned block —
    # the 5' mirror of B2 (stage-1 vicinity triggers): a linear alignment
    # through an intron shows garbage from the block start until matching
    # resumes at the true acceptor. Re-anchoring the pre-junction bases
    # moves reference_start LEFT and leaves everything downstream in frame,
    # so this is safe regardless of later junctions.
    first_m = next((j for j, (op, _) in enumerate(ct) if op in (0, 7, 8)), None)
    if first_m is not None and ct[first_m][1] >= 3 * cfg.arb_mm_win \
            and all(op in (4, 5) for op, _ in ct[:first_m]):
        walk = {j: (o, l, rs, qs) for j, o, l, rs, qs in _iter_ref_ops(ct)}
        _, mlen, m_rs, m_qs = walk[first_m]
        block_ref = read.reference_start + m_rs
        if dq is None:
            dq = _decoded_query(read, chrom_seq)
        eq_mode = '=' in q
        resume = None
        saw_flag = False
        for off in range(0, mlen - cfg.arb_mm_win, cfg.arb_mm_win // 2):
            if eq_mode:
                seg = q[m_qs + off:m_qs + off + cfg.arb_mm_win]
                mm = len(seg) - seg.count('=')
            else:
                seg = dq[m_qs + off:m_qs + off + cfg.arb_mm_win]
                ref = chrom_seq[block_ref + off:block_ref + off + cfg.arb_mm_win]
                mm = sum(1 for a_, b_ in zip(seg, ref) if a_ != b_)
            flagged = mm >= cfg.arb_mm_frac * cfg.arb_mm_win
            if off == 0 and not flagged:
                break
            if flagged:
                saw_flag = True
            elif saw_flag:
                resume = off
                break
        if resume is not None:
            _bump(stats, 'arb_mm_flagged')
            L2c = cfg.arb_seg
            accs = index.sites_in(
                chrom_key, right_kind,
                block_ref + max(10, resume - cfg.arb_mm_win - 20),
                block_ref + resume + 20)
            best_overall = None
            ed_cur_o = None
            best_o = None
            for e_alt in accs:
                e_alt = int(e_alt)
                o = e_alt - block_ref
                if o < 10 or mlen - o < 10:
                    continue
                q_j = m_qs + o
                L1 = min(cfg.arb_seg, q_j)
                L2 = min(L2c, mlen - o)
                if L1 < 10:
                    continue
                qwin = dq[q_j - L1:q_j + L2]
                a = assess_overhang(dq[q_j - L1:q_j], alpha=cfg.alpha,
                                    max_window=cfg.max_intron)
                if a.refused:
                    _bump(stats, 'arb_refused')
                    continue
                W = a.w_max_bp
                alts = []
                for d_alt in index.sites_in(chrom_key, left_kind,
                                            max(0, e_alt - W), e_alt - cfg.min_intron + 1):
                    d_alt = int(d_alt)
                    if d_alt - L1 < 0:
                        continue
                    alts.append((d_alt, e_alt,
                                 chrom_seq[d_alt - L1:d_alt] + chrom_seq[e_alt:e_alt + L2],
                                 qwin))
                if not alts:
                    continue
                cur_ref = chrom_seq[e_alt - L1:e_alt + L2]
                ed_cur, win, _ = _score_alts(qwin, cur_ref, alts, chrom_seq, cfg)
                if win is not None and (best_overall is None or win[0] < best_overall[0]):
                    best_overall = win
                    ed_cur_o = ed_cur
                    best_o = o
            if best_overall is not None:
                _, d_new, e_new = best_overall
                new_ct = ct[:first_m] + [(0, best_o), (3, e_new - d_new),
                                         (0, mlen - best_o)] + ct[first_m + 1:]
                read.cigartuples = new_ct
                read.reference_start = d_new - best_o
                for tag in ('MD', 'NM'):
                    if read.has_tag(tag):
                        read.set_tag(tag, None)
                read.set_tag('XB', f'mmL:{d_new}-{e_new}:{ed_cur_o:.1f}>{best_overall[0]:.1f}')
                _bump(stats, 'arb_mm_spliced')
                changed = True

    return changed


def resolve_read(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    index: SpliceSiteIndex,
    cfg: ResolverConfig,
    stats: ResolverStats,
) -> bool:
    """Attempt overhang resolution on both terminal clips of a primary
    alignment. Mutates ``read`` in place; returns True if anything changed."""
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        stats.passthrough_nonprimary += 1
        return False
    if not read.cigartuples or not read.query_sequence:
        return False

    chrom = standardize_chrom_name(read.reference_name) if read.reference_name else None
    if not chrom:
        return False
    if chrom in genome:
        chrom_key = chrom
    else:
        chrom_key = CHROM_TO_GENOME.get(chrom, '')
        if chrom_key not in genome:
            return False
    chrom_seq = genome[chrom_key]

    strand = '-' if read.is_reverse else '+'
    changed = False

    left_len, _ = _clip_lens(read.cigartuples)
    if left_len >= cfg.min_clip:
        clip_seq = read.query_sequence[:left_len]
        clip_used_len = min(left_len, cfg.max_clip_match)
        # Aligned bases adjacent to the clip, '='-decoded from the reference
        # (calmd encodes matches as '='), for inside-edge near sites. Capped
        # at first-M-length - 1 so the CIGAR rewrite always leaves a valid op.
        ct = read.cigartuples
        inside = ''
        if len(ct) > 1 and ct[1][0] == 0 and ct[1][1] > 1:
            n_in = min(cfg.edge_inside_slop, ct[1][1] - 1)
            raw = read.query_sequence[left_len:left_len + n_in]
            rs = read.reference_start
            inside = ''.join(
                chrom_seq[rs + i] if b == '=' else b for i, b in enumerate(raw))
        placement = resolve_clip(
            chrom_seq, index, chrom_key, _LEFT, strand, clip_seq,
            edge=read.reference_start, cfg=cfg, stats=stats, inside_seq=inside,
        )
        if placement is not None:
            _rewrite_cigar(read, _LEFT, placement, left_len, clip_used_len)
            stats.resolved += 1
            stats.resolved_left += 1
            if placement.k_inside:
                stats.extra['resolved_inside_edge'] = stats.extra.get('resolved_inside_edge', 0) + 1
            changed = True

    _, right_len = _clip_lens(read.cigartuples)
    if right_len >= cfg.min_clip and read.reference_end is not None:
        clip_seq = read.query_sequence[-right_len:]
        clip_used_len = min(right_len, cfg.max_clip_match)
        ct = read.cigartuples
        inside = ''
        if len(ct) > 1 and ct[-2][0] == 0 and ct[-2][1] > 1:
            n_in = min(cfg.edge_inside_slop, ct[-2][1] - 1)
            qlen = len(read.query_sequence)
            raw = read.query_sequence[qlen - right_len - n_in:qlen - right_len]
            re_ = read.reference_end
            inside = ''.join(
                chrom_seq[re_ - n_in + i] if b == '=' else b
                for i, b in enumerate(raw))
        placement = resolve_clip(
            chrom_seq, index, chrom_key, _RIGHT, strand, clip_seq,
            edge=read.reference_end, cfg=cfg, stats=stats, inside_seq=inside,
        )
        if placement is not None:
            _rewrite_cigar(read, _RIGHT, placement, right_len, clip_used_len)
            stats.resolved += 1
            stats.resolved_right += 1
            if placement.k_inside:
                stats.extra['resolved_inside_edge'] = stats.extra.get('resolved_inside_edge', 0) + 1
            changed = True

    # v2: re-arbitrate junction assignments and suspicious linear structure
    # (runs after clip resolution so a freshly placed junction competes too).
    if cfg.arb_enable:
        if _rearbitrate_read(read, chrom_seq, index, chrom_key, strand, cfg, stats):
            changed = True

    return changed


def run_overhang_resolver(
    base_bam: str,
    genome_path: str,
    output_bam: str,
    threads: int = 1,
    max_intron: int = 5000,
    alpha: float = 0.01,
    config: Optional[ResolverConfig] = None,
) -> str:
    """Stream the (name-sorted) minimap2 arm BAM through the resolver.

    Every input record is written (resolved or passthrough) in input order,
    so the output keeps the panel's name-sort convention with no extra sort
    step. Returns ``output_bam``. Stats are logged and attached to the
    function as ``run_overhang_resolver.last_stats`` for tests/drivers.
    """
    cfg = config or ResolverConfig(alpha=alpha, max_intron=max_intron)
    genome = load_genome(Path(genome_path))
    index = SpliceSiteIndex.load_or_build(str(genome_path), genome)
    stats = ResolverStats()

    with pysam.AlignmentFile(base_bam, 'rb', check_sq=False) as inp:
        header = inp.header.to_dict()
        header.setdefault('PG', []).append({
            'ID': 'rectify-overhang-resolver',
            'PN': 'rectify-overhang-resolver',
            'CL': f'alpha={cfg.alpha} max_intron={cfg.max_intron} base={base_bam}',
        })
        with pysam.AlignmentFile(output_bam, 'wb', header=header) as out:
            for read in inp.fetch(until_eof=True):
                stats.reads += 1
                resolve_read(read, genome, index, cfg, stats)
                out.write(read)

    logger.info(
        'overhang_resolver: %d reads, %d clips seen, %d resolved '
        '(L=%d R=%d), refused low-info=%d repeat=%d, no-candidates=%d, '
        'rejected edit=%d ambiguous=%d, candidates evaluated=%d',
        stats.reads, stats.clips_seen, stats.resolved,
        stats.resolved_left, stats.resolved_right,
        stats.refused_low_info, stats.refused_repeat, stats.no_candidates,
        stats.rejected_edit, stats.rejected_ambiguous, stats.candidates_evaluated,
    )
    run_overhang_resolver.last_stats = stats
    return output_bam
