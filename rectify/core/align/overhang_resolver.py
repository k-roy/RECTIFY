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
from ..splice.region_skip import overlaps_skip_region, skip_regions_from_env
from ..splice.clip_mappability import KmerIndex, best_out_of_window
from ..splice.repeat_expansion import is_repeat_expansion
from ..splice.splice_site_index import SpliceSiteIndex
from ...config import CHROM_TO_GENOME
from ...utils.genome import load_genome, standardize_chrom_name

logger = logging.getLogger(__name__)

# Arm 1 scores with PLAIN Levenshtein, never hp-ED: hp-ED charges 0.5 for a
# homopolymer indel, which is exactly the tolerance arm 1 exists to refuse
# (planning/770 -- the first pass was run with hp-ED and had to be discarded).
try:
    from rapidfuzz.distance import Levenshtein as _Lev
    # planning/619/651: a pip sdist fallback installs rapidfuzz pure-Python and
    # runs ~50x slower with no error. Check the BACKEND, not the version.
    if 'cpp' not in getattr(_Lev.distance, '__module__', ''):
        logger.warning(
            'rapidfuzz C++ backend unavailable (Levenshtein.distance resolves to %s) -- '
            'resolver arm 1 will run ~50x slower. Fix: pip install --force-reinstall '
            '--only-binary=:all: rapidfuzz',
            getattr(_Lev.distance, '__module__', '?'))
    _plain_ed = _Lev.distance
except ImportError:  # pragma: no cover - rapidfuzz is a core dependency
    _Lev = None

    def _plain_ed(a, b, score_cutoff=None):
        raise RuntimeError('resolver arm 1 requires rapidfuzz (core dependency)')

_LEFT = 'L'
_RIGHT = 'R'


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

    # --- v3 refusal clauses (planning/771 Brief A; measured in 769b-769d) ----
    # The 2026-08-09 `edge_inside_slop` feature (fc85930) let a near site 19 bp
    # INSIDE minimap2's aligned block create a candidate that did not previously
    # exist, and the P02 Ty1-delta-LTR artifact rode it through the gate at
    # hp-ED 35.0 vs threshold 36.0 — margin 1.0 (planning/769b). Three clauses
    # refuse that class without regressing the peelback the feature exists for.
    #
    # Clause A (Kevin's constraint, planning/769c): a peelback (k_inside > 0)
    # placement may win only if it beats the best no-peelback (k == 0)
    # placement by `min_margin` — today they compete on raw ED across DIFFERENT
    # comparison lengths, which is apples to oranges. If no k == 0 candidate
    # exists a k > 0 candidate may still win (clauses B/C then adjudicate).
    # 🔴 This is NOT a ban on re-assigning well-aligned bases: peelback is the
    # feature (`resolved_inside_edge` is 1,351 of 3,627 rescues in one
    # production sample), and the blanket prohibition floated in 769b rec 2 is
    # withdrawn. A k > 0 candidate that fails clause A is DEMOTED to the best
    # k == 0 candidate, not refused.
    clause_a: bool = True
    # Clause B (the within-window null, planning/769c/769d): before accepting,
    # score the placed query at EVERY position of the search window, not only
    # index-derived candidates; refuse when any position beats the accepted
    # placement by `clause_b_margin`. This is the clause that refuses P02 —
    # ~48 positions inside the resolver's own window score better than the one
    # it picked, and the resolver's ambiguity test could not see them because
    # it only compares indexed candidates (`rejected_ambiguous` fired 4 times
    # in 52,591 reads).
    clause_b: bool = True
    # MARGIN_B, hp-edits. FROZEN at 5.0 by planning/771 against 3,187 accepted
    # rescues over 6 fixture samples: refuses 26 (0.8%), ALL of them P02; zero
    # other panel rows, zero of 47 off-panel rescues. The separation is wide,
    # not marginal — the largest non-P02 improvement anywhere is 3.0 and P02's
    # are 13.0-18.5, so any margin in 4-12 is equivalent. Exposed as a config
    # field rather than a literal so the planning/772 PRODUCTION measurement
    # can move the default without a code change.
    clause_b_margin: float = 5.0
    # Clause C (genome-wide clip mappability, planning/769d): required for the
    # P07 class, which clause B structurally cannot see — P07 scores gain 0.0
    # in-window (its placement really is rank-1 within +/-5 kb) while the clip
    # places at hp-ED 2.0 sixty to sixty-seven kb away.
    # 🔴 MARGIN_C IS NOT CALIBRATED (planning/771 A3). Clause C therefore ships
    # REPORT-ONLY: with `clause_c_margin=None` it scores and records the
    # genome-wide competitor but never refuses. Set a float only after running
    # the planning/769d-style panel-wide impact measurement.
    # 🔴 Uniqueness is computed on the CLIP, never on the read's MAPQ — all four
    # artifact reads are MAPQ 60; the MAPQ 0 in 769b was the clip re-mapped
    # alone (769d correction 1).
    clause_c: bool = False
    clause_c_margin: Optional[float] = None

    # --- Arm 1: index-blind near-exact peel (planning/770, JUNCTION LEG ONLY) -
    # Peel-forward/peel-back around junctions minimap2 already found, with no
    # splice-site index and PLAIN Levenshtein (🔴 never hp-ED — hp-ED collapses
    # homopolymers, exactly the tolerance arm 1 exists to refuse).
    # 🔴 The 5' CLIP leg is deliberately absent: measured yield ZERO (173 clips
    # >= 24 nt, not one ed <= 2 20-mer match anywhere within +/-5 kb). Clips are
    # the error-burdened, motif-dependent case arm 2 exists for.
    # OFF by default: planning/770's discovery set is 11 reads at 3 loci from
    # one sample's first 3,000 reads and the whole-cohort number is unmeasured.
    arm1: bool = False
    arm1_x: int = 20             # query bases per side straddling the junction
    arm1_max_ed: float = 2.0     # plain-Levenshtein budget over the 2X-mer
    arm1_shift: int = 25         # max boundary shift searched (both directions)

    # Emit NM on every rewritten record. planning/769 defect 1: resolver
    # records carry NO NM (minimap2 NM=0, deSALT NM=19, mapPacBio NM=8, resolver
    # `-`), so a 70%-mismatched block is invisible to every consumer that
    # trusts NM — which is why nothing caught P02. Computed against the
    # '='-DECODED query, so it is correct even though calmd is deliberately not
    # re-run after the CIGAR surgery.
    emit_nm: bool = True

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
    canonical_rank: int  # 0 = GT/AG-class donor, 1 = GC-class
    k_inside: int = 0    # aligned bases re-assigned across the junction
                         # (inside-edge near site; planning/644 T4c)
    # The CONTIGUOUS half of the comparison — the query bases that become the
    # new M block and the genomic start they are placed at. Clause B scans this
    # against every window position; storing it avoids re-deriving the
    # side/lead geometry outside resolve_clip.
    qblk: str = ''
    blk_start: int = 0
    ed_blk: float = -1.0   # exact hp-ED of qblk at blk_start (clause B baseline)
    gain_b: float = 0.0    # best in-window improvement found by clause B
    gain_c: float = 0.0    # best out-of-window improvement found by clause C
    arm: str = 'arm2'      # which arm produced the call (planning/771 A5)


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


def _donor_rank(chrom_seq: str, side: str, strand: str,
                intron_start: int, intron_end: int) -> int:
    """0 for the canonical GT-class donor, 1 for GC-class."""
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


def _within_window_null(
    chrom_seq: str,
    qblk: str,
    blk_start: int,
    edge: int,
    ed_blk: float,
    margin: float,
    half_window: int,
) -> Tuple[float, int]:
    """Clause B: how much better does ``qblk`` score somewhere else in the window?

    Scans EVERY position in ``[edge - half_window, edge + half_window]`` -- not
    only the splice-site-indexed candidates -- with the DP cutoff seeded at
    ``ed_blk - margin``, and returns ``(gain, best_pos)`` where ``gain`` is
    ``ed_blk - best_alternative_ed`` (0.0 when nothing beats the bound).

    This is the hole planning/769c found: the resolver's ambiguity test compares
    only *indexed* candidates, so a dispersed repeat is unambiguous-by-
    construction inside W. For the P02 read ~48 non-site positions in the
    resolver's own window score strictly better than the placement it accepted;
    genuine rescues sampled the same way are rank 1 (two at rank 2).

    The scan runs ONCE per accepted rescue (722 of 52,591 reads in the fixture,
    3,627 of 913,703 in production), not once per candidate, and the cutoff is
    tight enough that almost every position prunes on the DP's first row.

    🔴 Window-local BY CONSTRUCTION. It cannot catch the P07 class, whose better
    placements sit 60-67 kb away and which scores gain 0.0 here -- that is what
    clause C is for (planning/769d correction 2).
    """
    cutoff = ed_blk - margin
    if cutoff < 0:
        return 0.0, -1
    n = len(qblk)
    lo = max(0, edge - half_window)
    hi = min(len(chrom_seq), edge + half_window) - n
    best_ed = None
    best_pos = -1
    for pos in range(lo, hi + 1):
        # The accepted placement is not its own competitor.
        if pos == blk_start:
            continue
        e2 = hp_edit_distance_bounded(qblk, chrom_seq[pos:pos + n], cutoff=cutoff)
        if e2 <= cutoff and (best_ed is None or e2 < best_ed):
            best_ed, best_pos = e2, pos
            if best_ed <= 0.0:
                break
    if best_ed is None:
        return 0.0, -1
    return ed_blk - best_ed, best_pos


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

    near_kind, far_kind = _site_kinds(side, strand, cfg.acceptor_classes)

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
                        k_inside=k,
                    ))
                    best_ed = min(best_ed, ed)

    if not placements:
        stats.no_candidates += 1
        return None

    # --- Selection: unambiguous winner only --------------------------------
    placements.sort(key=lambda p: (p.ed, p.canonical_rank, p.intron_start))
    best = placements[0]
    pool = placements

    # --- Clause A (planning/769c): peelback must BEAT no-peelback ----------
    # k>0 and k==0 candidates were competing on raw ED across DIFFERENT
    # comparison lengths (k>0 is scored over lc+k bases, k==0 over lc) -- apples
    # to oranges, and the accept threshold `max_edit_frac * m` RISES with
    # peelback length, so peeling was structurally cheap. Require the peelback
    # to win by `min_margin`; otherwise DEMOTE to the best no-peelback
    # placement rather than refuse. If no k==0 candidate exists at all (the P02
    # geometry -- `edge_slop=5` admitted nothing outside the edge) a k>0
    # candidate may still win, and clauses B/C adjudicate it.
    if cfg.clause_a and best.k_inside > 0:
        k0 = [p for p in placements if p.k_inside == 0]  # preserves the sort
        if k0 and not (best.ed <= k0[0].ed - cfg.min_margin):
            best, pool = k0[0], k0
            _bump(stats, 'clauseA_demoted')

    # comparison length = m + lead for k==0 (== lc) and m for k>0 (== lc + k);
    # both equal p.m + p.lead.
    if best.ed > cfg.max_edit_frac * (best.m + best.lead):
        stats.rejected_edit += 1
        return None
    for other in pool[1:]:
        if other.ed - best.ed >= cfg.min_margin:
            break  # sorted: everything after is at least as far
        if not same_junction(
            chrom_seq,
            (best.intron_start, best.intron_end),
            (other.intron_start, other.intron_end),
        ):
            stats.rejected_ambiguous += 1
            return None

    # --- The contiguous half of the comparison (clause B/C scan unit) ------
    # LEFT : cmp_seq[:m]  <-> chrom[f-m:f]   and cmp_seq[m:]  <-> chrom[e:edge]
    # RIGHT: cmp_seq[:lead] <-> chrom[edge:d] and cmp_seq[-m:] <-> chrom[e:e+m]
    # For k>0 lead is 0 and the whole comparison is contiguous.
    k = best.k_inside
    if side == _LEFT:
        cmp_seq = (clip_used + inside_seq[:k]) if k else clip_used
        best.qblk = cmp_seq[:best.m]
        best.blk_start = best.intron_start - best.m
    else:
        cmp_seq = (inside_seq[len(inside_seq) - k:] + clip_used) if k else clip_used
        best.qblk = cmp_seq[len(cmp_seq) - best.m:]
        best.blk_start = best.intron_end

    # --- Clause B (planning/769c/769d): the within-window null -------------
    # `best.ed` covers m+lead bases against a TWO-PIECE reference; the scan
    # compares m bases against contiguous windows, so the baseline is
    # recomputed on the contiguous block alone. (planning/769d scanned the
    # block but kept the whole-comparison ED as the baseline, which inflates
    # the gain slightly -- conservative in the refusing direction, and it
    # measured zero collateral at margin 5 anyway.)
    if cfg.clause_b:
        best.ed_blk = hp_edit_distance_bounded(
            best.qblk, chrom_seq[best.blk_start:best.blk_start + best.m])
        gain, _pos = _within_window_null(
            chrom_seq, best.qblk, best.blk_start, edge,
            best.ed_blk, cfg.clause_b_margin, cfg.max_intron)
        best.gain_b = gain
        if gain >= cfg.clause_b_margin:
            _bump(stats, 'refused_clause_b')
            _bump(stats, 'refused_clause_b_' + side)
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
    # APPEND, never overwrite: a read can be resolved on BOTH terminal clips,
    # and the previous unconditional set_tag silently dropped the first entry
    # (every consumer already splits XJ on ',' -- 769d/772 both do).
    entry = f'{placement.intron_start}-{placement.intron_end}:{placement.ed:.1f}:{side}'
    prev = str(read.get_tag('XJ')) if read.has_tag('XJ') else ''
    read.set_tag('XJ', f'{prev},{entry}' if prev else entry)

    # planning/771 A5 asks for the clause-B gain and the producing arm on XJ.
    # 🔴 They go on a SEPARATE tag instead, deliberately: XJ's three
    # colon-separated fields are load-bearing for existing consumers
    # (769d_clauseB_full.py and 772_clauseB_prod.py both do
    # `span, ed, side = fld.split(':')`), and a fourth field makes every one of
    # them raise ValueError and SILENTLY skip the record -- exactly the
    # invisible-failure class this whole change set exists to remove.
    #
    # 🔴 The tag is XW, NOT XQ. XQ:i is ALREADY TAKEN -- it is the cDNA 5'
    # pre-trim strip length written by core/cdna/io.py:248 and declared
    # READ_INTRINSIC in multialign/cma_schema.py:34. Every read of the 748
    # fixture's minimap2 arm carries one (5,000/5,000 sampled), so appending to
    # it produced `XQ:Z:140,arm=arm2;...` -- a type change (i -> Z) that
    # CLOBBERS a value correct-cdna depends on. Caught by reading a real
    # production record, not by any test. XW/XD/XE/XH/XP/XZ are all free.
    # Parallel to XJ: one comma-separated entry per XJ entry, same order.
    xw = f'arm={placement.arm};gB={placement.gain_b:.1f};gC={placement.gain_c:.1f};k={k}'
    prevw = str(read.get_tag('XW')) if read.has_tag('XW') else ''
    read.set_tag('XW', f'{prevw},{xw}' if prevw else xw)


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


def _compute_nm(cigartuples, ref_start: int, dq: str, chrom_seq: str) -> int:
    """SAM ``NM`` for a (possibly rewritten) record: mismatches + inserted +
    deleted bases; ``N`` is a skip, not a deletion, so it is not charged.

    planning/769 defect 1 -- the resolver emitted NO ``NM`` at all (minimap2
    NM=0, deSALT NM=19, mapPacBio NM=8, resolver ``-``), so a 70%-mismatched
    block was invisible to every consumer that trusts NM. That is the single
    largest reason the P02 artifact reached consensus unnoticed.

    ``dq`` must be the '='-DECODED query taken under the record's ORIGINAL
    placement: ``rectify align`` runs ``samtools calmd -e`` before this module
    and the resolver deliberately does not re-run it, so ``read.query_sequence``
    still carries '=' bytes that mean "matches the reference where minimap2 put
    me" -- decoding them at the NEW coordinates would score moved bases as
    perfect matches and defeat the whole point of the tag.
    """
    nm = 0
    qi = 0
    rp = ref_start
    for op, ln in cigartuples:
        if op in (0, 7, 8):          # M / = / X
            for t in range(ln):
                if dq[qi + t] != chrom_seq[rp + t]:
                    nm += 1
            qi += ln
            rp += ln
        elif op == 1:                # I
            nm += ln
            qi += ln
        elif op == 2:                # D
            nm += ln
            rp += ln
        elif op == 4:                # S
            qi += ln
        elif op == 3:                # N -- a skip, not an edit
            rp += ln
    return nm


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
    left_kind, right_kind = _boundary_kinds(strand, cfg.acceptor_classes)
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


def _apply_clause_c(
    placement: _Placement,
    genome: Dict[str, str],
    chrom_key: str,
    chrom_seq: str,
    cfg: ResolverConfig,
    stats: ResolverStats,
    kmer_index: Optional[KmerIndex],
) -> Optional[_Placement]:
    """Clause C: score the placed block genome-wide and refuse when an
    OUT-OF-WINDOW placement beats it by ``clause_c_margin``.

    Required for the P07 class, which clause B structurally cannot see: P07's
    placement is genuinely rank-1 within +/-5 kb (clause-B gain 0.0) while the
    clip places at hp-ED 2.0 sixty to sixty-seven kb away (planning/769d
    correction 2).

    🔴 REPORT-ONLY unless ``clause_c_margin`` is set. planning/771 A3 is explicit
    that MARGIN_C is **not calibrated** — derive it and report the panel-wide
    impact the way planning/769d did before adopting a value. With the margin
    unset this records ``gain_c`` on the XQ tag and a counter, and refuses
    nothing.
    """
    if kmer_index is None or not placement.qblk:
        return placement
    ed_base = placement.ed_blk
    if ed_base < 0:
        ed_base = hp_edit_distance_bounded(
            placement.qblk,
            chrom_seq[placement.blk_start:placement.blk_start + placement.m])
        placement.ed_blk = ed_base
    # Probe margin: with MARGIN_C unset we still want a usable REPORT, so probe
    # at 1.0 (anything better than the accepted placement is worth recording).
    probe = cfg.clause_c_margin if cfg.clause_c_margin is not None else 1.0
    gain, cn, cp = best_out_of_window(
        kmer_index, genome, placement.qblk, chrom_key,
        placement.blk_start - cfg.max_intron,
        placement.blk_start + cfg.max_intron,
        ed_base, probe, hp_edit_distance_bounded,
    )
    placement.gain_c = gain
    if gain > 0:
        _bump(stats, 'clauseC_better_elsewhere')
    if cfg.clause_c_margin is not None and gain >= cfg.clause_c_margin:
        _bump(stats, 'refused_clause_c')
        logger.debug('clause C refused %s:%d (gain %.1f, better at %s:%d)',
                     chrom_key, placement.blk_start, gain, cn, cp)
        return None
    return placement


# ---------------------------------------------------------------------------
# Arm 1 (planning/770): index-blind near-exact peel — JUNCTION LEG ONLY.
# ---------------------------------------------------------------------------

def _arm1_window_unique(chrom_seq: str, q: str, pos: int,
                        max_ed: int, half_window: int) -> bool:
    """Clause B for arm 1: no OTHER position in the window matches ``q`` within
    ``max_ed`` PLAIN edits.

    planning/771 A4: "an index-blind exact match inside a Ty1 LTR is the same
    repeat trap with a tighter threshold". The rescue path's clause B compares
    hp-EDs with a margin of 5; that is vacuous here because a near-exact 20-mer
    cannot be beaten by 5 edits. The equivalent test for a near-exact match is
    plain UNIQUENESS at the same tolerance."""
    n = len(q)
    lo = max(0, pos - half_window)
    hi = min(len(chrom_seq), pos + half_window) - n
    for p in range(lo, hi + 1):
        if p == pos:
            continue
        if _plain_ed(q, chrom_seq[p:p + n], score_cutoff=max_ed) <= max_ed:
            return False
    return True


def _arm1_junction_peel(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    cfg: ResolverConfig,
    stats: ResolverStats,
    dq0: Optional[str],
    ref_start0: int,
) -> bool:
    """Peel-forward/peel-back around a junction the aligner ALREADY found, with
    no splice-site index and PLAIN Levenshtein.

    Why it exists (planning/770): arm 2 can only propose placements the
    GT/AG-class :class:`SpliceSiteIndex` already contains, so a boundary the
    aligner mis-assigned by 1-4 bp to a non-indexed dinucleotide is invisible to
    it. Measured on one sample's first 3,000 reads: minimap2's current boundary
    already scores ed 0 on 1,447/1,813 N ops (80%); 11 junctions at 3 loci
    (chrII:407028-407122, chrII:406883-407183, chrII:604515-604927 = panel
    P03/P22/P36) have current ed > 2 AND exactly one ed <= 2 alternative — and
    arm 2 found NOTHING on any of them (`XJ` absent).

    🔴 PLAIN Levenshtein, never hp-ED. hp-ED charges 0.5 for a homopolymer indel;
    collapsing homopolymers is exactly the tolerance this arm exists to refuse.

    🔴 Ambiguity discipline is mandatory, not optional. Even at X=30 the median
    junction has 3-8 alternative ed<=2 placements within +/-25 (mean 4.43 at X=12
    -> 3.74 at X=30 — raising X barely helps). Requiring UNIQUENESS after a
    `same_junction` collapse is what took the discovery set from ~500 to 11.

    🔴 The 5' CLIP leg is deliberately not built: 173 clips >= 24 nt, ZERO with an
    ed<=2 20-mer match anywhere within +/-5 kb. Unsurprising in hindsight — if a
    clip matched near-exactly nearby, the aligner would have aligned it.

    Geometry note: this uses the EXACT formulation — a donor shift of ``a`` moves
    the query split by ``a`` too, so each candidate is scored against its own
    query window (the same convention :func:`_score_alts` uses). planning/770's
    measurement held the query window fixed, which is an approximation that lets
    the effective split float inside the DP; hit counts may therefore differ
    slightly from that table.
    """
    if _Lev is None:
        return False
    X = int(cfg.arm1_x)
    SH = int(cfg.arm1_shift)
    MAXED = int(cfg.arm1_max_ed)
    ct = list(read.cigartuples or [])
    if len(ct) < 3:
        return False
    dq = dq0 if dq0 is not None else _decoded_query(read, chrom_seq)
    L = len(chrom_seq)
    nops = len(ct)

    for i, op, ln, ref_off, q_off in list(_iter_ref_ops(ct)):
        if op != 3 or i == 0 or i + 1 >= nops:
            continue
        if ct[i - 1][0] not in (0, 7, 8) or ct[i + 1][0] not in (0, 7, 8):
            continue
        m_l, m_r = ct[i - 1][1], ct[i + 1][1]
        d = read.reference_start + ref_off
        e = d + ln
        qp = q_off
        if qp < X + SH or qp + X + SH > len(dq):
            continue
        if d - X - SH < 0 or e + X + SH > L:
            continue
        _bump(stats, 'arm1_junctions')

        # A boundary move that CHANGES the intron length displaces every
        # ref-consuming op downstream of the right flank (the same frame-safety
        # constraint Case A enforces). Rigid moves (b == a) are always safe.
        downstream_ref_free = all(o in (1, 4, 5) for o, _ in ct[i + 2:])

        cur = _plain_ed(dq[qp - X:qp + X], chrom_seq[d - X:d] + chrom_seq[e:e + X])
        if cur <= MAXED:
            _bump(stats, 'arm1_clean')
            continue

        hits: List[Tuple[int, int, int, int, int]] = []
        for a in range(-SH, SH + 1):
            if m_l + a < 1 or m_r - a < 1:
                continue
            # Left-half prefilter: ed<=MAXED over the 2X window requires the
            # left X bases alone to fit the budget AT THIS SPLIT. Every split is
            # scanned, so nothing reachable is lost — and this is what keeps the
            # scan ~200 DPs per junction instead of (2*SH+1)^2 = 2,601.
            if _plain_ed(dq[qp + a - X:qp + a], chrom_seq[d + a - X:d + a],
                         score_cutoff=MAXED) > MAXED:
                continue
            qwin = dq[qp + a - X:qp + a + X]
            for b in range(-SH, SH + 1):
                if a == 0 and b == 0:
                    continue
                if b != a and not downstream_ref_free:
                    continue
                d2, e2 = d + a, e + b
                if not (cfg.min_intron <= e2 - d2 <= cfg.max_intron):
                    continue
                if e2 + (m_r - a) > L:
                    continue
                v = _plain_ed(qwin, chrom_seq[d2 - X:d2] + chrom_seq[e2:e2 + X],
                              score_cutoff=MAXED)
                if v <= MAXED:
                    hits.append((v, a, b, d2, e2))
        if not hits:
            _bump(stats, 'arm1_no_hit')
            continue

        hits.sort()
        distinct: List[Tuple[int, int, int, int, int]] = []
        for h in hits:
            if not any(same_junction(chrom_seq, (h[3], h[4]), (g[3], g[4]))
                       for g in distinct):
                distinct.append(h)
        if len(distinct) != 1:
            _bump(stats, 'arm1_ambiguous')
            continue

        v, a, b, d2, e2 = distinct[0]
        if cfg.clause_b and not (
                _arm1_window_unique(chrom_seq, dq[qp + a - X:qp + a],
                                    d2 - X, MAXED, cfg.max_intron)
                and _arm1_window_unique(chrom_seq, dq[qp + a:qp + a + X],
                                        e2, MAXED, cfg.max_intron)):
            _bump(stats, 'arm1_refused_clause_b')
            continue

        new_ct = list(ct)
        new_ct[i - 1] = (ct[i - 1][0], m_l + a)
        new_ct[i] = (3, e2 - d2)
        new_ct[i + 1] = (ct[i + 1][0], m_r - a)
        read.cigartuples = new_ct
        for tag in ('MD', 'NM'):
            if read.has_tag(tag):
                read.set_tag(tag, None)
        entry = f'{d2}-{e2}:{float(v):.1f}:A1'
        prev = str(read.get_tag('XJ')) if read.has_tag('XJ') else ''
        read.set_tag('XJ', f'{prev},{entry}' if prev else entry)
        xw = f'arm=arm1;ed={v};cur={cur};from={d}-{e}'
        prevw = str(read.get_tag('XW')) if read.has_tag('XW') else ''
        read.set_tag('XW', f'{prevw},{xw}' if prevw else xw)
        _bump(stats, 'arm1_resolved')
        # One move per read: the CIGAR walk above is now stale.
        return True

    return False


def resolve_read(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    index: SpliceSiteIndex,
    cfg: ResolverConfig,
    stats: ResolverStats,
    sides: str = 'LR',
    kmer_index: Optional[KmerIndex] = None,
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

    # 🔴 Decode the query ONCE, HERE, before any rewrite. `rectify align` runs
    # `samtools calmd -e` upstream, so query_sequence carries '=' bytes meaning
    # "matches the reference AT MY CURRENT PLACEMENT". Every rewrite below moves
    # bases, after which those bytes are stale; decoding later would silently
    # score moved bases as perfect matches (planning/769 §2 -- the trap that
    # inflated the first mismatch survey to 89.4% and had to be discarded).
    dq0 = _decoded_query(read, chrom_seq) if cfg.emit_nm else None
    ref_start0 = read.reference_start

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
        if placement is not None and cfg.clause_c:
            placement = _apply_clause_c(placement, genome, chrom_key, chrom_seq,
                                        cfg, stats, kmer_index)
        if placement is not None:
            _rewrite_cigar(read, _LEFT, placement, left_len, clip_used_len)
            stats.resolved += 1
            stats.resolved_left += 1
            if placement.k_inside:
                stats.extra['resolved_inside_edge'] = stats.extra.get('resolved_inside_edge', 0) + 1
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
        if placement is not None and cfg.clause_c:
            placement = _apply_clause_c(placement, genome, chrom_key, chrom_seq,
                                        cfg, stats, kmer_index)
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

    # Arm 1 (planning/770): index-blind near-exact peel around junctions
    # minimap2 already found. Junction leg only -- the 5' clip leg measured
    # ZERO yield and is deliberately not built.
    if cfg.arm1:
        if _arm1_junction_peel(read, chrom_seq, cfg, stats, dq0, ref_start0):
            changed = True

    # planning/771 A5 / planning/769 defect 1: emit NM on every rewritten
    # record. Passthroughs keep the input aligner's NM untouched.
    if changed and cfg.emit_nm and dq0 is not None:
        try:
            read.set_tag('NM', _compute_nm(read.cigartuples, read.reference_start,
                                           dq0, chrom_seq), 'i')
        except IndexError:  # pragma: no cover - a rewrite ran off the contig
            _bump(stats, 'nm_compute_failed')

    return changed


def run_overhang_resolver(
    base_bam: str,
    genome_path: str,
    output_bam: str,
    threads: int = 1,
    max_intron: int = 5000,
    alpha: float = 0.01,
    acceptor_classes: str = 'canonical',
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
                                   acceptor_classes=acceptor_classes)
    if not cfg.skip_regions:
        cfg.skip_regions = skip_regions_from_env()
    genome = load_genome(Path(genome_path))
    index = SpliceSiteIndex.load_or_build(str(genome_path), genome)
    kmer_index = None
    if cfg.clause_c:
        kmer_index = KmerIndex.load_or_build(str(genome_path), genome)
        if cfg.clause_c_margin is None:
            logger.info('overhang_resolver: clause C is REPORT-ONLY '
                        '(clause_c_margin unset -- MARGIN_C is not calibrated, '
                        'planning/771 A3). Gains are recorded on XQ and counted '
                        'as clauseC_better_elsewhere; nothing is refused.')
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
                resolve_read(read, genome, index, cfg, stats,
                             kmer_index=kmer_index)
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
