"""Information-bounded splice-overhang resolver — the mapPacBio-role rewrite.

planning/641: mapPacBio contributes ZERO unique reads to the panel
(planning/583 §4a) — its actual role is *alternative placement of terminal
overhangs across splice junctions*. This module does that deliberately, with
an explicit false-discovery budget, instead of running a general-purpose
aligner with a mammalian constant (``maxindel=200000`` == ~12,500 chance AG
acceptors per window on both strands — simultaneously the compute bill and
the false-positive rate).

🔴 **SCOPE LIMIT — this module does NOT fully replace mapPacBio, and must not be
described as doing so** (measured, planning/721; do not re-litigate without new
data). On the SAME upf1Δ reads through minimap2 and mapPacBio, scoring the
junctions mapPacBio found that minimap2 missed:

===========================  =====  =========
class of mapPacBio-unique    n      recovered
===========================  =====  =========
annotated                    1,053  **99.1%**
novel NON-CANONICAL            417  **34.5%**
===========================  =====  =========

**Cause is structural, not a tuning knob.** Candidates come from
:class:`~rectify.core.splice.splice_site_index.SpliceSiteIndex`, which is built
on GT/AG-class dinucleotides only (``don_gt_plus`` / ``don_gc_plus`` /
``acc_plus`` / ``don_minus`` / ``acc_minus``). **A non-canonical junction has no
entry in that index and therefore cannot be enumerated as a candidate at all.**
mapPacBio aligns and lets the sequence decide; this module can only propose
placements the index already contains. ``arb_grammar`` (canonical-preference
tiebreak, measured −6.0 pp cryptic recovery) is a second, smaller motif
dependency layered on top.

Consequences for callers:
  * as a panel arm this is a net gain and belongs in the default panel
    (planning/720: +26.5k annotated junctions per 900k cDNA reads, 31 reads
    harmed, 0 impossible junctions introduced);
  * for **non-canonical discovery missions** (upf1Δ, prp18Δ, cryptic splicing)
    the panel is weaker without mapPacBio than with it, and this module does not
    close that gap. Set ``arb_grammar=False`` for such runs, and treat the
    remaining shortfall as a known, unfilled gap rather than an absence of
    signal.

**Partial mitigation (2026-08-17, planning/722b, opt-in):**
``ResolverConfig(acceptor_classes='prp18')`` extends the ACCEPTOR candidate
classes with the alternative-3'SS set measured in Roy et al. 2023 NAR
(gkad968): BG (TG/CG/GG) + non-G HAU (AT). Measured benefit on the published
junction table: utilized alt-3'SS become enumerable 47%→88% (prp18 mutants)
and 62%→85% (upf1Δ-only — NMD stabilizes the same isoforms, so this matters
for plain upf1Δ libraries). Measured price: acceptor candidate density
x4.82 (+2.3 bits of candidate space; ~x5 resolver CPU), which is why it is
an opt-in mission flag, not the default. It does NOT close the rest of the
721 gap: donors stay GT/GC (deliberate — 4/1,833 published alt-3'SS
junctions had non-canonical donors), and the flat-spectrum singleton noise
that dominates the mapPacBio-unique residual (planning/722) is not made
enumerable by any dinucleotide class. Everything below about a motif-FREE
search stands: extending further would change both the candidate space and
the false-discovery budget (``alpha`` is calibrated against the candidate
density), so it is a design change, not a flag.

**AT-AC introns (2026-09-04, opt-in, PAIRED):** ``ResolverConfig(atac=True)`` /
``--resolver-atac`` adds a second clip-enumeration pass over the index's
``don_at_*`` / ``acc_ac_*`` arrays with the pair enforced (AT donor ↔ AC
acceptor only), ranks such placements below GT..AG and GC..AG at equal ED,
and teaches the arbiter's canonical-class tests to accept AT..AC so a real
one is not grammar-snapped onto a chance GT..AG. Motivation: yeast splices
AT-AC through its major spliceosome (Talkish et al. 2019 PLoS Genet
15:e1008249 — an AT-AC junction in SUT635), and human AT-AC introns are the
U12-type class. Because the pair is enforced this is NOT the motif-free
extension refuted in planning/722; it adds one dinucleotide pair, not a
spectrum.

Every enumeration path is wired: the clip resolver (Station A's job), the
class-preserving boundary shifts, and — since 2026-09-05 — the arbiter's two
two-boundary DISCOVERY paths, the Case B1 intron-length-D -> N snap and both
mismatch-flagged linear rescues (B2 right, B3 left). In each, AT-AC runs as a
SEPARATE PAIRED pass over ``don_at_*``/``acc_ac_*``, never as a union with the
GT/GC..AG kinds — a union would admit AT..AG and GT..AC, which is exactly the
motif-free extension refuted above. The GT/GC pass always runs first and every
cross-pass comparison is strict, so an exact edit-distance tie keeps the
GT/GC..AG answer (the same ordering ``_donor_rank`` gives the clip resolver).

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
MD/NM are dropped on rewritten records (stale after CIGAR surgery), and a
calmd '='-compressed SEQ is decoded to real letters on any rewrite (a '=' only
means anything under the alignment it was written for — see
:func:`_set_decoded_sequence`). Records the resolver leaves alone keep their
'=' bytes, and their compression, untouched.

**Emission spells out the indels the score assumed** (:func:`_block_cigar`).
Scoring is indel-tolerant, so writing the placed exon block as one flat ``M``
produced alignments that could not express their own score — read r060_2601
(P06) was written ``...405N24M`` for a 24-bp block at 13/24 apparent identity.
The block is re-aligned against the reference at the accepted junction with the
affine-gap semi-global aligner from :mod:`~rectify.core.align.local_aligner`
(the same one ``rescue_3ss_truncation`` uses) and emitted as real M/I/D. This
is emission only: scoring, acceptance, the junction and ``XJ`` are untouched,
a gap-free block still emits exactly the old flat ``M``, and a degenerate
alignment falls back to it and is counted as ``emit_fallback_flat``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from .local_aligner import _align_left_anchored, _align_right_anchored
from ..splice.overhang_informativeness import (
    COUNTERS,
    assess_overhang,
    canonical_in_class,
    hp_edit_distance_bounded,
    same_junction,
)
from ..splice.region_skip import overlaps_skip_region, skip_regions_from_env
from ..splice.repeat_expansion import is_repeat_expansion
from ..splice.splice_site_index import SpliceSiteIndex
from ...config import CHROM_TO_GENOME
from ...utils.genome import load_genome, standardize_chrom_name

logger = logging.getLogger(__name__)

_LEFT = 'L'
_RIGHT = 'R'

# Emission allowance for the exon block a placement puts across the junction
# (mirrors ``local_aligner.align_clip_to_exon``'s ``max_indel``). The reference
# window handed to the affine-gap aligner is the block length plus this many
# bases on the FREE side only — the junction side is anchored — so a wildly
# wrong placement cannot manufacture a huge deletion. A block whose optimal
# alignment needs more indel than this is refused back to the flat M and
# counted as ``emit_fallback_flat``.
_EMIT_MAX_INDEL = 5


@dataclass
class ResolverConfig:
    alpha: float = 0.01          # false-discovery budget (641 T7: a published knob)
    max_intron: int = 5000       # per-reference intron cap (clamps W_max)
    # Acceptor candidate classes. 'canonical' = the AG-class index (default —
    # the planning/720 ADOPT verdict was measured here). 'prp18' additionally
    # enumerates the Roy et al. 2023 NAR (gkad968) alternative-3'SS classes:
    # BG (TG/CG/GG) + non-G HAU (AT) — see the SpliceSiteIndex docstring for
    # the measured price (acceptor density x4.82, +2.3 bits of candidate
    # space) and benefit (published utilized alt-3'SS enumerable 47%->88% in
    # prp18, 62%->85% in upf1Δ-only; planning/722b). OPT-IN for splicing
    # missions (upf1Δ / prp18Δ / prp18-AA); donors stay GT/GC — the published
    # set had 4/1,833 non-canonical donors. For such missions also consider
    # arb_grammar=False (the canonical-preference snap works against exactly
    # these junctions).
    acceptor_classes: str = 'canonical'
    # AT-AC introns — opt-in, PAIRED class (module docstring). True runs a
    # second clip-enumeration pass over don_at_*/acc_ac_* with the pair
    # enforced (never AT..AG / GT..AC), ranks AT-AC placements after GT..AG
    # and GC..AG at equal ED, and makes the arbiter's canonical-class tests
    # accept AT..AC. Default False = the planning/720-measured candidate
    # space, byte-identical.
    atac: bool = False
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

    # --- Pathological-contig circuit breaker --------------------------------
    # Candidate enumeration is near_sites x far_sites, and on a repetitive or
    # low-complexity contig that product explodes: measured 361 candidates per
    # clip on untrimmed cDNA (planning/681 CP6), and a reporter-construct contig
    # ran ~555x below baseline (planning/673) until an operator set
    # RECTIFY_SKIP_REGIONS by hand. Skip regions require knowing the bad contig
    # IN ADVANCE, so any unanticipated custom reference degrades to a silent
    # multi-hour stall. This bounds the per-clip DP work instead: past the
    # ceiling the clip is abandoned (the read passes through untouched) and
    # `refused_candidate_blowup` records it.
    #
    # CALIBRATION — read before changing this number.
    # Measured (planning/681 CP6, real resolver over 2,000 real cDNA molecules,
    # only the input sequence differing between arms):
    #     pre-trim-fix   17 clips ->  6,137 candidates  =  361 / clip
    #     post-trim-fix   9 clips ->     60 candidates  =  6.7 / clip
    # Production now lives in the ~6.7/clip regime, so 2000 is ~300x typical and
    # should never fire on healthy data.
    #
    # 🔴 The 361 figure is the PRE-FIX PATHOLOGY, not a healthy maximum, and
    # `tests/test_resolver_candidate_ceiling.py` pins the default above it as a
    # deliberately conservative floor. Do NOT "tighten toward the post-fix mean"
    # — a ceiling near 6.7 would refuse legitimate clips on repetitive or custom
    # references, which is precisely the case this guard exists to survive.
    # Re-deriving percentiles from pre-fix data would calibrate against a
    # distribution that no longer exists.
    #
    # The guard bounds COUNT, not time, and the two differ by the DP backend:
    # 2000 candidates is ~24 ms/clip with the numba kernel (RECTIFY_HP_ED_NUMBA)
    # but ~7 s/clip on the pure-Python path, so a numba-less run on a truly
    # pathological contig is still slow — bounded, but slow. If that case ever
    # shows up, prefer enabling numba over lowering this.
    #
    # A ceiling set too LOW silently refuses legitimate clips — the failure mode
    # this guard must not become. That is why every refusal is counted and
    # logged rather than silent.
    max_candidates_per_clip: int = 2000

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

    # Reads overlapping these reference regions bypass ALL junction rescue
    # (clip resolution + re-arbitration) and are written through untouched.
    # {chrom: [(start, end), ...]}, 0-based half-open. The canonical use is
    # the yeast rDNA repeat (region_skip.YEAST_RDNA_SPEC): rRNA is not a
    # spliceosomal substrate, and on the v5.1 run the locus was 47% of ALL
    # CPU (planning/644b). Populated from RECTIFY_SKIP_REGIONS by the driver.
    skip_regions: Dict[str, List[Tuple[int, int]]] = field(default_factory=dict)


@dataclass
class ResolverStats:
    reads: int = 0
    passthrough_nonprimary: int = 0
    clips_seen: int = 0
    clips_assessed: int = 0
    refused_low_info: int = 0
    refused_repeat: int = 0
    refused_candidate_blowup: int = 0
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
    canonical_rank: int  # 0 = GT/AG-class donor, 1 = GC-class, 2 = AT-AC (paired, opt-in)
    k_inside: int = 0    # aligned bases re-assigned across the junction
                         # (inside-edge near site; planning/644 T4c)
    qseg: str = ''       # the query bases this placement was SCORED on
                         # (``cmp_seq``: lead + m, already '='-decoded). Kept so
                         # the emitter can align them and write a real M/I/D
                         # CIGAR instead of the flat M the score never implied.


def _clip_lens(cigartuples) -> Tuple[int, int]:
    """(left_softclip, right_softclip) lengths; 0 when absent."""
    if not cigartuples:
        return 0, 0
    left = cigartuples[0][1] if cigartuples[0][0] == 4 else 0
    right = cigartuples[-1][1] if cigartuples[-1][0] == 4 else 0
    return left, right


def _acc_kind(base: str, acceptor_classes: str) -> str:
    """Map an acceptor kind to its union variant under 'prp18' classes."""
    if acceptor_classes == 'prp18' and base in ('acc_plus', 'acc_minus'):
        return base + '_all'
    return base


def _site_kinds(side: str, strand: str,
                acceptor_classes: str = 'canonical') -> Tuple[str, str]:
    """(near_kind, far_kind) for a clip side x transcript strand."""
    if side == _LEFT:
        near, far = (('acc_plus', 'don_plus') if strand == '+'
                     else ('don_minus', 'acc_minus'))
    else:
        near, far = (('don_plus', 'acc_plus') if strand == '+'
                     else ('acc_minus', 'don_minus'))
    return (_acc_kind(near, acceptor_classes), _acc_kind(far, acceptor_classes))


def _site_kinds_atac(side: str, strand: str) -> Tuple[str, str]:
    """(near_kind, far_kind) for the PAIRED AT-AC class (ResolverConfig.atac).

    Same geometry as :func:`_site_kinds` with the four AT-AC arrays
    substituted; the pair is enforced because the two kinds are only ever
    queried together."""
    if side == _LEFT:
        return (('acc_ac_plus', 'don_at_plus') if strand == '+'
                else ('don_at_minus', 'acc_ac_minus'))
    return (('don_at_plus', 'acc_ac_plus') if strand == '+'
            else ('acc_ac_minus', 'don_at_minus'))


def _boundary_kinds_atac(strand: str) -> Tuple[str, str]:
    """(left_kind, right_kind) of an AT-AC intron's genomic boundaries."""
    return (('don_at_plus', 'acc_ac_plus') if strand == '+'
            else ('acc_ac_minus', 'don_at_minus'))


def _is_atac(chrom_seq: str, strand: str, intron_start: int, intron_end: int) -> bool:
    """True iff the WRITTEN junction is an AT-AC intron on ``strand``."""
    if intron_start < 0 or intron_end > len(chrom_seq):
        return False
    pair = (chrom_seq[intron_start:intron_start + 2].upper(),
            chrom_seq[intron_end - 2:intron_end].upper())
    return pair == (('AT', 'AC') if strand == '+' else ('GT', 'AT'))


def _donor_rank(chrom_seq: str, side: str, strand: str,
                intron_start: int, intron_end: int) -> int:
    """0 for the canonical GT-class donor, 1 for GC-class, 2 for a (paired)
    AT-AC intron. AT-AC is only ever enumerated under ``ResolverConfig.atac``,
    so the default ranking is unchanged."""
    if _is_atac(chrom_seq, strand, intron_start, intron_end):
        return 2
    if strand == '+':
        dinuc = chrom_seq[intron_start:intron_start + 2].upper()
        return 0 if dinuc == 'GT' else 1
    dinuc = chrom_seq[intron_end - 2:intron_end].upper()
    return 0 if dinuc == 'AC' else 1


_BLOWUP_WARNED: set = set()


def _warn_blowup(chrom: str, ceiling: int) -> None:
    """Log a candidate blow-up ONCE per contig.

    A pathological contig produces one blow-up per clip — potentially millions
    — so an unthrottled warning would itself become the performance problem it
    is reporting. One line names the contig, which is the actionable part: the
    operator adds it to RECTIFY_SKIP_REGIONS. The per-clip tally stays exact in
    `stats.refused_candidate_blowup`.
    """
    if chrom in _BLOWUP_WARNED:
        return
    _BLOWUP_WARNED.add(chrom)
    logger.warning(
        "overhang_resolver: candidate blow-up on contig %r — a clip exceeded "
        "%d candidates and was ABANDONED (read passes through untouched). "
        "This is the repetitive/low-complexity-contig class; further "
        "occurrences on this contig are counted in refused_candidate_blowup "
        "but not logged. Consider adding %r to RECTIFY_SKIP_REGIONS.",
        chrom, ceiling, chrom,
    )


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

    # One (near, far) kind pair per enabled class. The AT-AC pair is a
    # separate pass, never a union, so the pair constraint holds by
    # construction (ResolverConfig.atac).
    kind_pairs = [_site_kinds(side, strand, cfg.acceptor_classes)]
    if cfg.atac:
        kind_pairs.append(_site_kinds_atac(side, strand))

    # --- Candidate lookup (binary search, bounded by W) --------------------
    placements: List[_Placement] = []
    best_ed = float('inf')
    # Per-clip candidate budget (see ResolverConfig.max_candidates_per_clip).
    # Checked inside both enumeration loops rather than pre-counted: the DP's
    # pruning cutoff tightens as `best_ed` improves, so the loops must run in
    # their existing order, and a pre-count would have to duplicate the range
    # arithmetic and could drift from it. The cost of bailing late is bounded
    # by the ceiling itself.
    n_cand = 0

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
    for near_kind, far_kind in kind_pairs:
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
                    n_cand += 1
                    if n_cand > cfg.max_candidates_per_clip:
                        stats.refused_candidate_blowup += 1
                        _warn_blowup(chrom_key, cfg.max_candidates_per_clip)
                        return None
                    _c = _cutoff()
                    ed = hp_edit_distance_bounded(cmp_seq, ref, cutoff=_c)
                    if ed <= _c:
                        placements.append(_Placement(
                            ed=ed, intron_start=f, intron_end=e, lead=lead, m=m,
                            canonical_rank=_donor_rank(chrom_seq, side, strand, f, e),
                            k_inside=k, qseg=cmp_seq,
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
                    n_cand += 1
                    if n_cand > cfg.max_candidates_per_clip:
                        stats.refused_candidate_blowup += 1
                        _warn_blowup(chrom_key, cfg.max_candidates_per_clip)
                        return None
                    _c = _cutoff()
                    ed = hp_edit_distance_bounded(cmp_seq, ref, cutoff=_c)
                    if ed <= _c:
                        placements.append(_Placement(
                            ed=ed, intron_start=d, intron_end=e, lead=lead, m=m,
                            canonical_rank=_donor_rank(chrom_seq, side, strand, d, e),
                            k_inside=k, qseg=cmp_seq,
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


def _block_cigar(
    side: str,
    qblock: str,
    chrom_seq: str,
    intron_start: int,
    intron_end: int,
    max_indel: int = _EMIT_MAX_INDEL,
) -> Optional[Tuple[List[Tuple[int, int]], int]]:
    """Real M/I/D CIGAR for the exon block a placement puts across the junction.

    Returns ``(ops, block_ref_start)``, or None when the aligner's answer is
    degenerate and the caller must fall back to the flat M.

    Why this exists: :func:`resolve_clip` SCORES a candidate with
    ``hp_edit_distance_bounded`` (indels allowed) but the emitter used to write
    the matched block as ONE flat M, so an accepted placement whose true fit
    needs a 1-bp shift was written as an ungapped block the alignment could not
    express — read r060_2601 (P06) got ``...405N24M`` for a 24-bp block at
    13/24 apparent identity. MD/NM are dropped in the same breath, so nothing
    downstream noticed the disagreement. The block is re-aligned here with the
    same affine-gap semi-global aligner ``rescue_3ss_truncation`` uses
    (``local_aligner.align_clip_to_exon`` is the precedent), anchored at the
    junction it must not move:

      RIGHT clip — the block starts AT ``intron_end``: left-anchored, free suffix.
      LEFT clip  — the block ends AT ``intron_start``: right-anchored, free prefix.

    EMISSION ONLY. Scoring, acceptance, the accepted junction and the ``XJ``
    edit distance are all untouched; the only thing that changes is how the
    already-accepted block is spelled.

    Degeneracy (fall back to the flat M rather than write something worse):
      * total indel above ``max_indel`` — the window allows no more anyway,
        and a placement needing that much was mis-scored;
      * a D touching the junction — that D silently redefines the intron
        boundary the caller is about to write as ``N``. An I touching the
        junction is NOT refused: it adds bases without moving the junction.
    """
    m = len(qblock)
    if m == 0:
        return None
    if side == _RIGHT:
        region_start = intron_end
        region_end = min(len(chrom_seq), intron_end + m + max_indel)
    else:
        region_start = max(0, intron_start - m - max_indel)
        region_end = intron_start
    ref_region = chrom_seq[region_start:region_end]
    if len(ref_region) < m:
        return None

    # Exact shortcut, not a heuristic: a block with zero ungapped mismatches
    # scores 2*m, and every gap costs at least gap_open+gap_extend = -5 against
    # a best-case +2, so no gapped alignment can tie it. Reproducing the flat M
    # here keeps the (common) clean case free of the O(m^2) Gotoh DP and makes
    # the "gap-free segment emits exactly the old CIGAR" invariant structural.
    ungapped = ref_region[:m] if side == _RIGHT else ref_region[-m:]
    if qblock.upper() == ungapped.upper():
        return ([(0, m)],
                intron_end if side == _RIGHT else intron_start - m)

    if side == _RIGHT:
        ops, _ref_consumed = _align_left_anchored(qblock, ref_region)
        block_ref_start = intron_end
    else:
        ops, ref_skip = _align_right_anchored(qblock, ref_region)
        block_ref_start = region_start + ref_skip
    if not ops:
        return None
    if sum(ln for op, ln in ops if op in (0, 1)) != m:
        return None          # defensive: the aligner must consume all of qblock
    if sum(ln for op, ln in ops if op in (1, 2)) > max_indel:
        return None
    junction_op = ops[0][0] if side == _RIGHT else ops[-1][0]
    if junction_op == 2:
        return None
    return ops, block_ref_start


def _rewrite_cigar(
    read: pysam.AlignedSegment,
    side: str,
    placement: _Placement,
    clip_len: int,
    clip_used_len: int,
    chrom_seq: Optional[str] = None,
    stats: Optional[ResolverStats] = None,
) -> None:
    """Apply an accepted placement to ``read`` in place (CIGAR + position).

    The matched portion of the clip becomes ``M[/I/D] / N / M[/I/D]``; any clip
    bases beyond ``max_clip_match`` remain soft-clipped. MD/NM are dropped
    (stale). ``chrom_seq`` enables the gapped emission of :func:`_block_cigar`;
    without it (or without ``placement.qseg``) the historical flat M is written.

    The ``lead`` bases stay a flat M merged into the neighbouring op on
    purpose: they sit between the near site and the untouched aligned edge, so
    their reference span is pinned at BOTH ends and any gap inside them would
    have to be a balanced I+D pair — churn, not information, and it would move
    the N. Only the ``m`` block, which is free at its outer end, is re-aligned.
    """
    ct = list(read.cigartuples)
    remainder = clip_len - clip_used_len
    intron_len = placement.intron_end - placement.intron_start
    k = placement.k_inside

    block_ops: Optional[List[Tuple[int, int]]] = None
    block_ref_start: Optional[int] = None
    if chrom_seq is not None and placement.qseg:
        # RIGHT: cmp_seq is lead-then-m; LEFT: m-then-lead (see resolve_clip).
        qblock = (placement.qseg[placement.lead:] if side == _RIGHT
                  else placement.qseg[:placement.m])
        got = _block_cigar(side, qblock, chrom_seq,
                           placement.intron_start, placement.intron_end)
        if got is None:
            if stats is not None:
                _bump(stats, 'emit_fallback_flat')
        else:
            block_ops, block_ref_start = got
    if block_ops is None:
        block_ops = [(0, placement.m)]
        block_ref_start = (placement.intron_end if side == _RIGHT
                           else placement.intron_start - placement.m)

    if side == _LEFT:
        new_ops: List[Tuple[int, int]] = []
        if remainder > 0:
            new_ops.append((4, remainder))
        new_ops.extend(block_ops)
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
        # The block's REFERENCE span is no longer placement.m once it carries
        # I/D ops, so the new start comes from the aligner's own consumption
        # (block_ref_start + ref span == intron_start by construction).
        read.reference_start = block_ref_start
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
        new_ops.extend(block_ops)
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


def _boundary_kinds(strand: str,
                    acceptor_classes: str = 'canonical') -> Tuple[str, str]:
    """(left_kind, right_kind) of an intron's genomic boundaries."""
    if strand == '+':
        left, right = 'don_plus', 'acc_plus'
    else:
        left, right = 'acc_minus', 'don_minus'
    return (_acc_kind(left, acceptor_classes), _acc_kind(right, acceptor_classes))


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


def _set_decoded_sequence(read: pysam.AlignedSegment, decoded: str) -> None:
    """Replace a calmd '='-compressed SEQ with the decoded letters.

    A '=' byte means "identical to the reference AT THIS POSITION", so it is
    only interpretable under the alignment it was written for. Every rewrite in
    this module moves query bases to new reference positions — the clip
    resolver re-assigns ``k_inside`` bases across the junction, and the
    arbiter's boundary shifts, D-op snaps and mismatch rescues move whole
    blocks by an intron length — so a surviving '=' would decode against the
    WRONG reference downstream. The resolver already decodes for SCORING
    (``_decoded_query`` / the ``inside_seq`` the caller builds); this writes the
    same answer back into the record so the two agree.

    ``decoded`` MUST have been computed against the INPUT alignment, before any
    CIGAR surgery. pysam clears ``query_qualities`` whenever ``query_sequence``
    is assigned, so they are saved across the assignment and restored — a
    silently qual-stripped BAM would be a worse bug than the one being fixed.

    Input BAMs that carry '=' are common, not exotic: ``samtools calmd -e``
    writes them and the 748 fixture's minimap2 BAM has '=' on 42,409/42,409
    reads.
    """
    if len(decoded) != len(read.query_sequence):
        # Only reachable from a CIGAR that does not cover the whole query, i.e.
        # a corrupt record. Loud, per planning/638 §3 — assigning a shorter SEQ
        # under the existing CIGAR would corrupt it further.
        raise ValueError(
            "overhang_resolver: '='-decoded SEQ is %d bases but the record "
            "carries %d — the CIGAR does not cover the query (corrupt input)"
            % (len(decoded), len(read.query_sequence)))
    quals = read.query_qualities
    read.query_sequence = decoded
    if quals is not None:
        read.query_qualities = quals


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
    cur_canon = (canonical_in_class(chrom_seq, *cur_junc, atac=cfg.atac)
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
                and canonical_in_class(chrom_seq, d_alt, e_alt, atac=cfg.atac)):
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
    left_kind, right_kind = _boundary_kinds(strand, cfg.acceptor_classes)
    # One (left, right) boundary-kind pair per enabled class, for the arbiter's
    # two-boundary DISCOVERY paths (Case B1's D-op snap and the B2/B3
    # mismatch-flagged linear rescues, which choose BOTH boundaries). The AT-AC
    # pair is a separate pass and never a union with the GT/GC..AG kinds — a
    # union would enumerate AT..AG and GT..AC, the motif-free extension
    # planning/722 refuted (module docstring; Chanfreau planning/722). The
    # GT/GC pass is first, so an exact ED tie keeps the GT/GC..AG answer,
    # matching _donor_rank's rank-2 placement of AT-AC.
    pair_kinds = [(left_kind, right_kind)]
    if cfg.atac:
        pair_kinds.append(_boundary_kinds_atac(strand))
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
                        canonical_in_class(chrom_seq, d - dl, e, atac=cfg.atac):
                    ct = ct[:i - 1] + [(3, nlen + dl)] + ct[i + 1:]
                    _bump(stats, 'arb_dmerge')
                    merged = True
                    changed = True
                    break
            if i + 1 < len(ct) and ct[i + 1][0] == 2:
                dl = ct[i + 1][1]
                if ((e + dl) - d) <= cfg.max_intron and \
                        canonical_in_class(chrom_seq, d, e + dl, atac=cfg.atac):
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
        # AT-AC (opt-in): boundary shifts preserve the junction's CLASS — an
        # AT..AC intron's alternatives are other AT/AC sites, never GT/AG
        # (which would let a chance GT..AG displace a real AT-AC site).
        if cfg.atac and _is_atac(chrom_seq, strand, d, e):
            lk, rk = _boundary_kinds_atac(strand)
        else:
            lk, rk = left_kind, right_kind
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
        if ed_cur < cfg.arb_margin and canonical_in_class(chrom_seq, d, e, atac=cfg.atac):
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
            for e_alt in index.sites_in(chrom_key, rk, e - W, e + W + 1):
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
            for d_alt in index.sites_in(chrom_key, lk, d - W, d + W + 1):
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
            for d_alt in index.sites_in(chrom_key, lk, d - W, d + W + 1):
                d_alt = int(d_alt)
                delta = d_alt - d
                if delta == 0:
                    continue
                e_alt = d_alt + span
                if index.sites_in(chrom_key, rk, e_alt, e_alt + 1).size == 0:
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
        # The pair loop encloses BOTH site loops: an AT donor may only be
        # offered with an AC acceptor (see `pair_kinds` above).
        for lk_b, rk_b in pair_kinds:
            for d_alt in index.sites_in(chrom_key, lk_b, p_start - s2, p_start + s2 + 1):
                d_alt = int(d_alt)
                delta = d_alt - p_start
                if m_l + delta < 1 or m_r - delta < 1 or d_alt - L1 < 0:
                    continue
                if q_split + delta - L1 < 0 or q_split + delta + L2 > len(dq):
                    continue
                for e_alt in index.sites_in(chrom_key, rk_b, p_end - s2, p_end + s2 + 1):
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
        #
        # The rank term keeps that tiebreak class-aware. Unlike _score_alts,
        # this selection takes scored[0] outright with no ambiguity check, so
        # without it an opt-in AT-AC alt could displace a GT..AG one at equal
        # ED purely on coordinate order — the opposite of the rank-2 ordering
        # _donor_rank gives the clip resolver.
        ed_cur = hp_edit_distance_bounded(qwin, cur_ref)
        scored = sorted(
            (hp_edit_distance_bounded(aq, ref, cutoff=ed_cur),
             2 if _is_atac(chrom_seq, strand, da, ea) else 0, da, ea)
            for da, ea, ref, aq in alts)
        best = scored[0]
        if best[0] > ed_cur or best[0] > cfg.max_edit_frac * len(qwin):
            _bump(stats, 'arb_no_gain')
            continue
        _, _, d_new, e_new = best
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
            best_overall = None
            ed_cur_o = None
            # The pair loop encloses the DONOR loop, not just the acceptor
            # lookup: pairing an AT donor with an AG acceptor is exactly the
            # union `pair_kinds` exists to prevent. The GT/GC pass runs first
            # and the comparison below is strict, so an exact ED tie keeps the
            # GT/GC..AG answer.
            for lk_b, rk_b in pair_kinds:
                # The true donor sits where matching stops — anywhere within
                # the flagged window (mismatches begin mid-window), so search a
                # span covering the window plus slack on both sides.
                donors = index.sites_in(
                    chrom_key, lk_b,
                    block_ref + max(0, onset - 20),
                    block_ref + onset + cfg.arb_mm_win + 20)
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
                    for e_alt in index.sites_in(chrom_key, rk_b,
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
            best_overall = None
            ed_cur_o = None
            best_o = None
            # Mirror of B2: the pair loop encloses the ACCEPTOR loop here, so
            # an AC acceptor is only ever offered with an AT donor.
            for lk_b, rk_b in pair_kinds:
                accs = index.sites_in(
                    chrom_key, rk_b,
                    block_ref + max(10, resume - cfg.arb_mm_win - 20),
                    block_ref + resume + 20)
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
                    for d_alt in index.sites_in(chrom_key, lk_b,
                                                max(0, e_alt - W),
                                                e_alt - cfg.min_intron + 1):
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
    sides: str = 'LR',
) -> bool:
    """Attempt overhang resolution on both terminal clips of a primary
    alignment. Mutates ``read`` in place; returns True if anything changed.

    ``sides`` restricts which genomic clip sides are attempted (``'L'``,
    ``'R'``, or the default ``'LR'``). The triage terminal-clip leg passes a
    single side so the 5' clip stays the Cat3/rescue leg's territory and every
    overhang is assessed exactly once (one refusal discipline across legs —
    dev/REALIGNER_LANDSCAPE_AND_PATH_20260809.md §4 step 2)."""
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        stats.passthrough_nonprimary += 1
        return False
    if not read.cigartuples or not read.query_sequence:
        return False

    chrom = standardize_chrom_name(read.reference_name) if read.reference_name else None
    if not chrom:
        return False
    if cfg.skip_regions and overlaps_skip_region(
            cfg.skip_regions, chrom, read.reference_start, read.reference_end):
        _bump(stats, 'skipped_region')
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

    # calmd '='-compressed SEQ: decode ONCE against the INPUT alignment (the
    # only one the '=' bytes are valid under) and write the letters back on the
    # first rewrite. It has to happen before the arbiter runs, not at the end:
    # the arbiter builds its own `dq = _decoded_query(read, chrom_seq)` from
    # whatever alignment it is handed, so if the clip resolver has already moved
    # `k_inside` bases across a junction, those '=' would decode against the new
    # reference positions and the arbiter would score garbage exactly at the
    # junction it is arbitrating. A read the resolver leaves alone keeps its '='
    # bytes untouched — including the arbiter's `eq_mode` mismatch-count fast
    # path, which is the common case.
    decoded_seq = (_decoded_query(read, chrom_seq)
                   if '=' in read.query_sequence else None)

    def _decode_seq_now():
        nonlocal decoded_seq
        if decoded_seq is not None:
            _set_decoded_sequence(read, decoded_seq)
            _bump(stats, 'seq_eq_decoded')
            decoded_seq = None

    left_len, _ = _clip_lens(read.cigartuples)
    if _LEFT in sides and left_len >= cfg.min_clip:
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
            _rewrite_cigar(read, _LEFT, placement, left_len, clip_used_len,
                           chrom_seq=chrom_seq, stats=stats)
            stats.resolved += 1
            stats.resolved_left += 1
            if placement.k_inside:
                stats.extra['resolved_inside_edge'] = stats.extra.get('resolved_inside_edge', 0) + 1
            _decode_seq_now()
            changed = True

    _, right_len = _clip_lens(read.cigartuples)
    if _RIGHT in sides and right_len >= cfg.min_clip and read.reference_end is not None:
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
            _rewrite_cigar(read, _RIGHT, placement, right_len, clip_used_len,
                           chrom_seq=chrom_seq, stats=stats)
            stats.resolved += 1
            stats.resolved_right += 1
            if placement.k_inside:
                stats.extra['resolved_inside_edge'] = stats.extra.get('resolved_inside_edge', 0) + 1
            _decode_seq_now()
            changed = True

    # v2: re-arbitrate junction assignments and suspicious linear structure
    # (runs after clip resolution so a freshly placed junction competes too).
    if cfg.arb_enable:
        if _rearbitrate_read(read, chrom_seq, index, chrom_key, strand, cfg, stats):
            # Every arbiter rewrite except the A0 D-merge moves query bases to
            # new reference positions (a boundary shift swaps |delta| bases
            # between the flanking Ms; a B1 snap with delta != 0 does the same;
            # B2/B3 split an M into M-N-M and displace everything past the split
            # by the intron length), so the surviving '=' bytes are stale here
            # too. The D-merge alone is alignment-identical, but distinguishing
            # it would buy nothing: decoding is idempotent and correct either way.
            _decode_seq_now()
            changed = True

    return changed


def run_overhang_resolver(
    base_bam: str,
    genome_path: str,
    output_bam: str,
    threads: int = 1,
    max_intron: int = 5000,
    alpha: float = 0.01,
    acceptor_classes: str = 'canonical',
    atac: bool = False,
    config: Optional[ResolverConfig] = None,
) -> str:
    """Stream the (name-sorted) minimap2 arm BAM through the resolver.

    Every input record is written (resolved or passthrough) in input order,
    so the output keeps the panel's name-sort convention with no extra sort
    step. Returns ``output_bam``. Stats are logged and attached to the
    function as ``run_overhang_resolver.last_stats`` for tests/drivers.

    .. warning::
       ``threads`` is **accepted but NOT implemented** — the body below is a
       single-threaded stream. The parameter is kept because callers already
       pass it (``align_command`` forwards ``args.threads``), but it buys
       nothing, and jobs have been sized against it: a DRS array requesting 8
       slots got one core's worth of resolver. It is logged rather than
       silently ignored so capacity planning stops inheriting the wrong
       number. Per-contig sharding would be the natural implementation and
       would also bound pathological contigs (the rDNA / reporter-construct
       class that needs ``RECTIFY_SKIP_REGIONS`` today).
    """
    if threads and threads > 1:
        logger.warning(
            "overhang_resolver: threads=%d requested but the resolver is "
            "SINGLE-THREADED (parameter accepted for API compatibility, not "
            "implemented). Size jobs for 1 core on this stage.",
            threads,
        )
    cfg = config or ResolverConfig(alpha=alpha, max_intron=max_intron,
                                   acceptor_classes=acceptor_classes, atac=atac)
    if not cfg.skip_regions:
        cfg.skip_regions = skip_regions_from_env()
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
        '(L=%d R=%d), refused low-info=%d repeat=%d blowup=%d, '
        'no-candidates=%d, rejected edit=%d ambiguous=%d, '
        'candidates evaluated=%d',
        stats.reads, stats.clips_seen, stats.resolved,
        stats.resolved_left, stats.resolved_right,
        stats.refused_low_info, stats.refused_repeat,
        stats.refused_candidate_blowup, stats.no_candidates,
        stats.rejected_edit, stats.rejected_ambiguous, stats.candidates_evaluated,
    )
    # The emission counters live in `extra` (which `as_dict()` folds into the
    # per-dataset stats JSON align_command writes), but the info line above
    # enumerates named fields only, so they would never reach an operator
    # reading a log. `emit_fallback_flat` in particular says the resolver wrote
    # a block it could not justify a gapped spelling for — surface it.
    _fallbacks = stats.extra.get('emit_fallback_flat', 0)
    _decoded = stats.extra.get('seq_eq_decoded', 0)
    if _fallbacks or _decoded:
        logger.info(
            'overhang_resolver: %d block(s) emitted as a flat M after a '
            'degenerate gapped alignment (emit_fallback_flat), %d record(s) '
            "had a calmd '=' SEQ decoded on rewrite",
            _fallbacks, _decoded,
        )
    # Escalate the blow-up count out of the info line: an acceptance gate that
    # skims the summary must not have to notice a non-zero field buried mid-row.
    if stats.refused_candidate_blowup:
        logger.warning(
            'overhang_resolver: %d clip(s) ABANDONED on the candidate ceiling '
            '(%d). Those reads passed through unresolved — real junctions may '
            'be unplaced on the affected contig(s). Investigate before trusting '
            'junction counts there.',
            stats.refused_candidate_blowup, cfg.max_candidates_per_clip,
        )
    run_overhang_resolver.last_stats = stats
    return output_bam
