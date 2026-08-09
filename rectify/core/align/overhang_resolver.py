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
