"""
triage.py — read-evidence triage of consensus alignments + targeted re-alignment.

The consensus-triage layer of the Rectify Re-aligner
(design: dev/DENOVO_ALIGNER_REASSESSMENT_20260809.md §3):

  1. **Classify** each consensus/corrected alignment from READ EVIDENCE only —
     junction-proximal mismatch/indel weight (``_count_junction_proximity_errors``),
     unexplained 5'/3' soft clips, and junction annotation status. Clean
     alignments are labelled ``high_confidence`` and BYPASSED (the do-no-harm
     data say they are at ceiling: touching them is downside + compute).
  2. **Re-align** only the triaged minority, one leg per failure mode.
     JUNCTION leg: motif-blind refinement
     (``refine_bam_junctions(motif_blind=True)``, compensating-indel invariant
     always-on) over a junction pool built from the FULL BAM's evidence.
     CLIP legs (landscape §4 step 2; OFF by default —
     ``TriagePolicy.clip_legs_enable``): the TERMINAL-CLIP leg hands the 3'
     clip to the overhang resolver (``resolve_read``, shared
     ``SpliceSiteIndex``); the 5'-CLIP leg routes to the Cat3 rescue
     machinery (``rescue_3ss_truncation`` + the bam_writer application path).
     ONE refusal discipline across both: every overhang is assessed exactly
     once via ``assess_overhang`` (the resolver assesses internally; the 5'
     leg assesses before invoking the rescue) and a refused overhang is never
     sequence-searched — refusal is a counted, first-class outcome.
     ORIGINAL-ALIGNMENT leg (``original_bams``): the pre-correction record for
     each triaged read is offered as a further candidate. The junction leg's
     proposer IS Module 2H, so on a BAM 2H already refined it re-derives its
     own fixed point — the one placement it can never propose is the one it
     moved away from (measured: the arbiter accepted the original on 42/42
     damaged reads when handed it, while only 6/93 were actually reverted;
     dev/sumner_misplaced_panel_20260904 CP3). This adds a CANDIDATE, not an
     acceptance path.
  3. **Re-entry, never auto-accept**: a re-aligned candidate replaces the
     original ONLY if it strictly improves the HP-aware edit distance
     (``_cigar_hp_edit_distance``) — the validated consensus arbiter (the C3
     headroom result: given a truth candidate, hp_ed picks it; the panel's
     failures are missing-candidate failures). The re-aligner proposes; the
     arbiter disposes.

NOT implemented here (second gate, by design — see the reassessment §3
caution 1): the POOL-LEVEL discovery gate. A flattened non-canonical junction
looks CLEAN read-by-read (the aligner snaps it to a nearby canonical motif
with tidy overhangs), so read-level distress triage cannot see it. Discovery
triage operates on loci (minority-distress / recurrent alternate placements)
and is a separate, cluster-validated stage.

Classification is deliberately blind to any internal model score: every
signal is observable read evidence. (The score-gaming history is the reason
— see dev/ALIGNER_BENCH_STATE_AUDIT_20260721.md.)
"""

from __future__ import annotations

import logging
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pysam

from .extract import extract_junctions_from_cigar, get_softclip_lengths
from .scoring import _count_junction_proximity_errors

logger = logging.getLogger(__name__)

LABEL_HIGH_CONFIDENCE = 'high_confidence'
LABEL_TRIAGE = 'triage'

# Reasons a read is triaged (stable strings — they land in the output TSV).
REASON_JUNCTION_ERRORS = 'junction_proximal_errors'
REASON_CLIP_5P = 'softclip_5p'
REASON_CLIP_3P = 'softclip_3p'
REASON_UNANNOTATED_JUNCTION = 'unannotated_junction'

# Reasons the JUNCTION re-align leg acts on. Soft-clip-only triage waits for
# the clip legs (Cat3 5' rescue / overhang_resolver).
_JUNCTION_LEG_REASONS = frozenset({REASON_JUNCTION_ERRORS, REASON_UNANNOTATED_JUNCTION})

# Strict-improvement epsilon for the hp_ed re-entry comparison: a candidate
# must beat the incumbent by MORE than float noise. Ties keep the incumbent
# (the aligner's placement is the null hypothesis).
_REENTRY_EPS = 1e-9


@dataclass
class TriagePolicy:
    """Read-evidence thresholds for the high-confidence bypass.

    Defaults are deliberately conservative starting points, to be tuned on the
    validation corpora; they are policy knobs, not measured constants.
    """
    max_junction_proximal_errors: float = 1.0   # HP-weighted, summed over windows
    max_clip_5p: int = 30                       # bases of 5' soft clip tolerated
    max_clip_3p: int = 30                       # bases of 3' soft clip tolerated
    triage_unannotated_junctions: bool = True   # unannotated junction => re-check placement
    junction_window_bp: int = 5                 # window for the proximity-error signal
    # Clip legs (landscape §4 step 2). OFF by default: with False the triage
    # output is byte-identical to the junction-leg-only behavior.
    clip_legs_enable: bool = False


@dataclass
class TriageResult:
    read_id: str
    label: str
    reasons: List[str] = field(default_factory=list)
    junction_proximal_errors: float = 0.0
    clip_5p: int = 0
    clip_3p: int = 0
    n_junctions: int = 0
    n_unannotated: int = 0

    @property
    def junction_leg(self) -> bool:
        """True when the junction re-align leg should act on this read."""
        return any(r in _JUNCTION_LEG_REASONS for r in self.reasons)

    @property
    def clip5_leg(self) -> bool:
        """True when the 5'-clip (Cat3 rescue) leg should act on this read."""
        return REASON_CLIP_5P in self.reasons

    @property
    def clip3_leg(self) -> bool:
        """True when the terminal-clip (resolver) leg should act on this read."""
        return REASON_CLIP_3P in self.reasons


def _annotated_3tuples(
    annotated_junctions: Optional[Set[Tuple]],
) -> Optional[Set[Tuple[str, int, int]]]:
    """Normalise an annotated-junction set to (chrom, start, end) 3-tuples."""
    if annotated_junctions is None:
        return None
    out: Set[Tuple[str, int, int]] = set()
    for j in annotated_junctions:
        if len(j) >= 3:
            out.add((str(j[0]), int(j[1]), int(j[2])))
    return out


def classify_read(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    annotated_3: Optional[Set[Tuple[str, int, int]]] = None,
    policy: Optional[TriagePolicy] = None,
) -> TriageResult:
    """Classify one alignment as ``high_confidence`` or ``triage``.

    Annotation membership is exact-coordinate for now; ambiguity-window
    matching (the up_amb/down_amb enhancement) is a follow-up, so a junction
    one bp inside a repeat may read as unannotated — which errs toward
    re-checking, never toward silently trusting.
    """
    policy = policy or TriagePolicy()
    res = TriageResult(read_id=read.query_name or '', label=LABEL_HIGH_CONFIDENCE)

    res.junction_proximal_errors = _count_junction_proximity_errors(
        read, genome, junction_window_bp=policy.junction_window_bp)
    res.clip_5p, res.clip_3p = get_softclip_lengths(read)

    junctions = extract_junctions_from_cigar(read)
    res.n_junctions = len(junctions)
    if annotated_3 is not None and junctions:
        chrom = read.reference_name or ''
        res.n_unannotated = sum(
            1 for (js, je) in junctions if (chrom, js, je) not in annotated_3)

    if res.junction_proximal_errors > policy.max_junction_proximal_errors:
        res.reasons.append(REASON_JUNCTION_ERRORS)
    if res.clip_5p > policy.max_clip_5p:
        res.reasons.append(REASON_CLIP_5P)
    if res.clip_3p > policy.max_clip_3p:
        res.reasons.append(REASON_CLIP_3P)
    if (policy.triage_unannotated_junctions and annotated_3 is not None
            and res.n_unannotated > 0):
        res.reasons.append(REASON_UNANNOTATED_JUNCTION)

    if res.reasons:
        res.label = LABEL_TRIAGE
    return res


def _classifiable(read: pysam.AlignedSegment) -> bool:
    return (not read.is_unmapped and not read.is_secondary
            and not read.is_supplementary and read.cigartuples is not None)


def classify_bam(
    bam_path: str,
    genome: Dict[str, str],
    annotated_junctions: Optional[Set[Tuple]] = None,
    policy: Optional[TriagePolicy] = None,
) -> List[TriageResult]:
    """Classify every primary mapped read in a BAM. Order = BAM order."""
    annotated_3 = _annotated_3tuples(annotated_junctions)
    out: List[TriageResult] = []
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if _classifiable(read):
                out.append(classify_read(read, genome, annotated_3, policy))
    return out


def _read_signature(read: pysam.AlignedSegment) -> Tuple:
    return (read.reference_name, read.reference_start, tuple(read.cigartuples or ()))


def reentry_accept(
    original: pysam.AlignedSegment,
    candidate: pysam.AlignedSegment,
    genome: Dict[str, str],
    penalty_table=None,
) -> bool:
    """The re-entry gate: accept the candidate ONLY on strict hp_ed improvement.

    The candidate never auto-replaces the original — it competes under the
    same HP-aware edit distance the consensus winner-selection uses. Ties and
    regressions keep the original.
    """
    from .corrected_consensus import _cigar_hp_edit_distance
    if _read_signature(original) == _read_signature(candidate):
        return False
    ed_old = _cigar_hp_edit_distance(original, genome, penalty_table)
    ed_new = _cigar_hp_edit_distance(candidate, genome, penalty_table)
    return ed_new < ed_old - _REENTRY_EPS


# ---------------------------------------------------------------------------
# Clip legs (landscape §4 step 2). Both legs PROPOSE candidate rewrites on a
# COPY of the incumbent record; the existing arbiter (``reentry_accept``,
# strict hp_ed improvement) DECIDES. Refusals are counted first-class.
# ---------------------------------------------------------------------------

def two_sided_candidate_seam(read, index, cfg) -> List:
    """SEAM — two-sided junction-hypothesis enumeration (landscape §4b).

    Post-landing resolver-library work (owned by the resolver line, NOT this
    branch) will enumerate BOTH one-sided-anchor hypotheses per flagged
    junction — (A) donor held, acceptor scanned ±w; (B) acceptor held, donor
    scanned ±w (grammar-free on the scanned side) — and emit per-read evidence
    records for pool-recurrence adjudication (the chrXI:93365 11/11 off-by-1
    class; see dev/REALIGNER_LANDSCAPE_AND_PATH_20260809.md §4b and
    644h_realigner_integration.md §3). Until that lands this hook returns no
    candidates; anything it returns will enter the SAME adjudication path as
    the resolver's proposals (reentry_accept + pool aggregation).
    """
    return []


def _clip_copy(read: pysam.AlignedSegment) -> pysam.AlignedSegment:
    """Materialized deep copy of a record (safe to mutate / outlive files)."""
    return pysam.AlignedSegment.from_dict(read.to_dict(), read.header)


def _five_prime_clip_leg(
    cand: pysam.AlignedSegment,
    genome: Dict[str, str],
    candidate_junctions: Set[Tuple[str, int, int]],
    rcfg,
) -> Tuple[bool, bool]:
    """5'-clip leg: route the 5' soft clip through the Cat3 rescue machinery.

    Gate first (the ONE refusal discipline): the junction-proximal portion of
    the 5' clip is assessed exactly once via ``assess_overhang`` (plus the
    repeat-expansion short-circuit, mirroring the resolver's ``resolve_clip``);
    a refused overhang is NEVER sequence-searched. ``rescue_3ss_truncation``'s
    internal informativeness gate stays env-dark (``RECTIFY_OVERHANG_INFO_GATE``),
    so by default the assessment here is the only one. On sequence-confirmed
    rescue the BAM surgery is applied with the same read_edits machinery
    bam_writer uses (reanchor pre-pass + ``extend_read_5prime_for_junction_rescue``)
    — reuse, not reimplementation.

    Mutates ``cand`` in place on a proposal. Returns ``(proposed, refused)``.
    """
    from rectify.core.splice.overhang_informativeness import assess_overhang
    from rectify.core.splice.repeat_expansion import is_repeat_expansion

    ct = cand.cigartuples
    seq = cand.query_sequence
    if not ct or not seq:
        return False, False
    strand = '-' if cand.is_reverse else '+'
    if strand == '+':
        clip = seq[:ct[0][1]] if ct[0][0] == 4 else ''
        clip_used = clip[-rcfg.max_clip_match:]   # attaches at its right end
    else:
        clip = seq[-ct[-1][1]:] if ct[-1][0] == 4 else ''
        clip_used = clip[:rcfg.max_clip_match]    # attaches at its left end
    if not clip:
        return False, False

    # --- The shared gate: refuse before any sequence search ----------------
    if is_repeat_expansion(clip_used):
        return False, True
    assessment = assess_overhang(clip_used, alpha=rcfg.alpha,
                                 max_window=rcfg.max_intron)
    if assessment.refused:
        return False, True

    if not candidate_junctions:
        return False, False

    from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation
    res = rescue_3ss_truncation(cand, genome, candidate_junctions, strand)
    if not res.get('rescued') or res.get('rescue_type') not in (
            'softclip', 'mpb_mismatch'):
        # intronic_snap / proximity are TSV-level outcomes in the correct
        # stage; the triage leg proposes BAM rewrites only for the
        # sequence-confirmed classes.
        return False, False

    from rectify.core.bam.read_edits import (
        _apply_reanchor_from_clip_len,
        extend_read_5prime_for_junction_rescue,
    )
    changed = False
    rcl = int(res.get('reanchor_clip_len', 0) or 0)
    if rcl > 0:
        changed |= _apply_reanchor_from_clip_len(cand, rcl)
    sc_len = rcl if rcl > 0 else len(clip)
    changed |= extend_read_5prime_for_junction_rescue(
        cand,
        res['five_prime_corrected'],
        sc_len,
        strand,
        exon_cigar_str=res.get('five_prime_exon_cigar', '') or '',
        upstream_trim=int(res.get('five_prime_upstream_trim', 0) or 0),
    )
    return changed, False


def _terminal_clip_leg(
    cand: pysam.AlignedSegment,
    genome: Dict[str, str],
    index,
    rcfg,
    rstats,
) -> Tuple[bool, int]:
    """Terminal-clip leg: hand the 3'-end soft clip to the overhang resolver.

    The resolver is restricted to the 3' clip's genomic side (``sides=``) so
    the 5' clip stays the Cat3 leg's territory and each overhang is assessed
    exactly once (the resolver's internal ``assess_overhang`` call is the one
    assessment for this clip). Mutates ``cand`` in place on a proposal.
    Returns ``(proposed, n_refused_delta)``.
    """
    from rectify.core.align.overhang_resolver import resolve_read

    side = 'L' if cand.is_reverse else 'R'   # genomic side of the 3' clip
    refused_before = rstats.refused_low_info + rstats.refused_repeat
    changed = resolve_read(cand, genome, index, rcfg, rstats, sides=side)
    # §4b seam: two-sided enumeration candidates would compete here, through
    # the same adjudication as the resolver's rewrite. No-op until the
    # resolver-library work lands.
    for _cand in two_sided_candidate_seam(cand, index, rcfg):
        pass
    refused_delta = (rstats.refused_low_info + rstats.refused_repeat) - refused_before
    return changed, refused_delta


def triage_realign_bam(
    input_bam: str,
    output_bam: str,
    genome: Dict[str, str],
    annotated_junctions: Set[Tuple],
    policy: Optional[TriagePolicy] = None,
    pool_bams: Optional[List[str]] = None,
    original_bams: Optional[List[str]] = None,
    penalty_table=None,
    realign: bool = True,
    refine_kwargs: Optional[dict] = None,
    sort_and_index: bool = True,
    resolver_config=None,
    splice_site_index=None,
) -> Tuple[List[dict], Dict[str, int]]:
    """Classify a consensus BAM; motif-blind-re-align the triaged junction reads.

    Pipeline (see module docstring):
      classify → junction pool from the FULL evidence (``pool_bams`` or the
      input BAM itself) → ``refine_bam_junctions(motif_blind=True)`` on ONLY
      the triaged junction-leg reads → hp_ed re-entry per read → the
      PRE-CORRECTION alignment as a further candidate (``original_bams``) →
      clip legs (``policy.clip_legs_enable``, OFF by default) → final BAM.

    ``original_bams`` (ISSUE-009): the junction leg's only proposer is
    ``refine_bam_junctions(motif_blind=True)`` — Module 2H itself — so on a BAM
    2H has already refined it re-derives its own fixed point and can never
    offer the placement it moved away from. Measured on the Sumner RNA004 panel
    (dev/sumner_misplaced_panel_20260904 CP3): triage flagged the damage
    correctly (R1 16/16, R2 16/16 distressed) and ``reentry_accept`` accepted
    the ORIGINAL back on 42/42 reads when handed it directly, yet only 6/93
    damaged reads were reverted, because nothing proposed the original. Passing
    the pre-correction BAM(s) here adds that alignment to the candidate pool.

    This adds a CANDIDATE, never an acceptance path: each original record
    competes against the current incumbent under the same strict-hp_ed
    ``reentry_accept``, and only reads already triaged for the junction leg are
    offered one. Records are rebuilt against the INPUT BAM's header
    (``from_dict`` re-resolves ``ref_name`` by NAME), so a pre-correction BAM
    whose @SQ order differs cannot relocate a read to the wrong chromosome.

    Clip legs (landscape §4 step 2): when ``policy.clip_legs_enable`` is True,
    reads triaged for a 3' soft clip are handed per-read to the overhang
    resolver (``resolve_read`` with the shared ``splice_site_index`` — built
    from ``genome`` when not supplied — and ``resolver_config``, default
    ``ResolverConfig(arb_enable=False)``); reads triaged for a 5' soft clip
    route through the Cat3 rescue machinery behind the shared informativeness
    gate. Every proposal enters the SAME strict hp_ed re-entry arbiter;
    refusals are counted and the read passes through unchanged.

    Returns ``(rows, stats)`` where ``rows`` are per-read TSV dicts and
    ``stats`` counts {classified, high_confidence, triaged, junction_leg,
    realigned, accepted} plus the original-alignment counters (``orig_leg``,
    ``orig_proposed``, ``orig_accepted``, ``orig_skipped_unknown_chrom``) and
    the clip-leg counters (``clip5_leg``, ``clip3_leg``,
    ``clip{5,3}_refused/_proposed/_accepted`` — all zero when the legs are
    disabled).
    """
    from rectify.core.splice.junction_refiner import refine_bam_junctions
    from rectify.core.splice.junction_scoring import build_junction_pool

    policy = policy or TriagePolicy()
    annotated_3 = _annotated_3tuples(annotated_junctions)

    # ── Pass 1: classify ────────────────────────────────────────────────────
    results: Dict[str, TriageResult] = {}
    n_dup = 0
    with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam:
        # Kept for the original-alignment candidates: every candidate record is
        # rebuilt against THIS header so `ref_name` re-resolves by name.
        input_header = bam.header
        for read in bam.fetch(until_eof=True):
            if not _classifiable(read):
                continue
            r = classify_read(read, genome, annotated_3, policy)
            if r.read_id in results:
                n_dup += 1
                continue
            results[r.read_id] = r
    if n_dup:
        logger.warning("triage: %d duplicate primary read ids — first record kept", n_dup)

    junction_ids = {rid for rid, r in results.items() if r.junction_leg}
    stats = {
        'classified': len(results),
        'high_confidence': sum(1 for r in results.values() if r.label == LABEL_HIGH_CONFIDENCE),
        'triaged': sum(1 for r in results.values() if r.label == LABEL_TRIAGE),
        'junction_leg': len(junction_ids),
        'realigned': 0,
        'accepted': 0,
        # Original-alignment candidates (zero without `original_bams`)
        'orig_leg': 0,
        'orig_proposed': 0,
        'orig_accepted': 0,
        'orig_skipped_unknown_chrom': 0,
        # Clip-leg counters (always present; zero while clip_legs_enable=False)
        'clip5_leg': 0,
        'clip3_leg': 0,
        'clip5_refused': 0,
        'clip3_refused': 0,
        'clip5_proposed': 0,
        'clip3_proposed': 0,
        'clip5_accepted': 0,
        'clip3_accepted': 0,
    }

    accepted: Dict[str, pysam.AlignedSegment] = {}
    realigned_ids: Set[str] = set()
    reverted_ids: Set[str] = set()
    originals: Dict[str, pysam.AlignedSegment] = {}

    if realign and junction_ids:
        with tempfile.TemporaryDirectory(prefix='rectify_triage_') as tmp:
            triaged_bam = os.path.join(tmp, 'triaged.bam')
            refined_bam = os.path.join(tmp, 'refined.bam')

            with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam, \
                    pysam.AlignmentFile(triaged_bam, 'wb', header=bam.header) as out:
                for read in bam.fetch(until_eof=True):
                    if (_classifiable(read) and read.query_name in junction_ids
                            and read.query_name not in originals):
                        out.write(read)
                        originals[read.query_name] = read

            # Junction pool from the FULL evidence — never from the triaged
            # subset alone (a subset pool would weaken exactly the candidates
            # the triaged reads need).
            pool, annot3_std = build_junction_pool(
                pool_bams or [input_bam], annotated_junctions)

            kwargs = dict(
                motif_blind=True,
                sort_and_index=False,
                n_workers=1,
            )
            if refine_kwargs:
                kwargs.update(refine_kwargs)
            refine_bam_junctions(
                triaged_bam, refined_bam,
                aligner_bams=[],
                annotated_junctions=annotated_junctions,
                genome=genome,
                prebuilt_junction_pool=pool,
                prebuilt_annotated_set=annot3_std,
                **kwargs,
            )

            # ── Re-entry: strict hp_ed improvement or the original stays ────
            with pysam.AlignmentFile(refined_bam, 'rb', check_sq=False) as ref_bam:
                for cand in ref_bam.fetch(until_eof=True):
                    rid = cand.query_name
                    orig = originals.get(rid)
                    if orig is None or not _classifiable(cand):
                        continue
                    if _read_signature(cand) == _read_signature(orig):
                        continue
                    realigned_ids.add(rid)
                    if reentry_accept(orig, cand, genome, penalty_table):
                        # Materialise before the temp file vanishes.
                        accepted[rid] = pysam.AlignedSegment.from_dict(
                            cand.to_dict(), cand.header)
            stats['realigned'] = len(realigned_ids)
            stats['accepted'] = len(accepted)

    # ── The PRE-CORRECTION alignment as a candidate (ISSUE-009) ─────────────
    # The junction leg's proposer is Module 2H, so on a BAM 2H already refined
    # it re-derives its own fixed point: the placement it moved AWAY from is
    # the one candidate it structurally cannot offer. The pre-correction record
    # is enumerable evidence — an alignment that was actually observed — so it
    # belongs in the candidate pool, judged by the SAME arbiter as every other
    # proposal. Nothing here can accept anything the strict hp_ed gate would
    # not accept; if 2H's proposal genuinely wins, the original loses to it and
    # the read is not reverted.
    if realign and junction_ids and original_bams:
        input_names = set(input_header.references or ())
        for orig_path in original_bams:
            seen_orig: Set[str] = set()
            with pysam.AlignmentFile(orig_path, 'rb', check_sq=False) as ob:
                for rec in ob.fetch(until_eof=True):
                    rid = rec.query_name
                    if (rid not in junction_ids or rid in seen_orig
                            or not _classifiable(rec)):
                        continue
                    seen_orig.add(rid)
                    if rec.reference_name not in input_names:
                        # A reference the output header cannot express. Skipping
                        # is the only safe move: a reference_id carried across
                        # headers silently relocates the read.
                        stats['orig_skipped_unknown_chrom'] += 1
                        continue
                    # Rebuild against the INPUT header: to_dict() writes
                    # `ref_name` as a string, so from_dict re-resolves it by
                    # NAME and a differing @SQ order cannot move the read.
                    cand = pysam.AlignedSegment.from_dict(rec.to_dict(),
                                                          input_header)
                    # NB explicit None checks: pysam defines __len__ on a
                    # record (the query length), so `a or b` would silently
                    # skip a zero-length one.
                    incumbent = accepted.get(rid)
                    if incumbent is None:
                        incumbent = originals.get(rid)
                    if incumbent is None:
                        continue
                    stats['orig_leg'] += 1
                    if _read_signature(cand) == _read_signature(incumbent):
                        continue
                    stats['orig_proposed'] += 1
                    realigned_ids.add(rid)
                    if reentry_accept(incumbent, cand, genome, penalty_table):
                        accepted[rid] = cand
                        reverted_ids.add(rid)
                        stats['orig_accepted'] += 1
        stats['realigned'] = len(realigned_ids)
        stats['accepted'] = len(accepted)
        if stats['orig_skipped_unknown_chrom']:
            logger.warning(
                'triage: %d pre-correction record(s) skipped — their reference '
                'name is not in the input BAM header (chromosome naming '
                'mismatch between the aligner BAM and the corrected BAM)',
                stats['orig_skipped_unknown_chrom'])
        logger.info(
            'triage original-alignment candidates: offered=%d proposed=%d '
            'accepted (reverted to pre-correction)=%d',
            stats['orig_leg'], stats['orig_proposed'], stats['orig_accepted'])

    # ── Clip legs (OFF by default; landscape §4 step 2) ─────────────────────
    clip5_ids: Set[str] = set()
    clip3_ids: Set[str] = set()
    if realign and policy.clip_legs_enable:
        clip5_ids = {rid for rid, r in results.items() if r.clip5_leg}
        clip3_ids = {rid for rid, r in results.items() if r.clip3_leg}
    if clip5_ids or clip3_ids:
        from rectify.core.align.overhang_resolver import (
            ResolverConfig, ResolverStats)
        from rectify.core.splice.splice_site_index import SpliceSiteIndex

        # arb_enable=False: junction re-arbitration is the junction leg's /
        # Station-C's territory — the clip leg stays clip-scoped.
        rcfg = resolver_config or ResolverConfig(arb_enable=False)
        index = splice_site_index or SpliceSiteIndex.build(genome)
        rstats = ResolverStats()
        seen_clip: Set[str] = set()
        with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam:
            for read in bam.fetch(until_eof=True):
                rid = read.query_name
                if (not _classifiable(read) or rid in seen_clip
                        or (rid not in clip5_ids and rid not in clip3_ids)):
                    continue
                seen_clip.add(rid)
                incumbent = accepted.get(rid, read)

                if rid in clip5_ids:
                    stats['clip5_leg'] += 1
                    cand = _clip_copy(incumbent)
                    proposed, refused = _five_prime_clip_leg(
                        cand, genome, annotated_3 or set(), rcfg)
                    if refused:
                        stats['clip5_refused'] += 1
                    if proposed:
                        stats['clip5_proposed'] += 1
                        realigned_ids.add(rid)
                        if reentry_accept(incumbent, cand, genome, penalty_table):
                            accepted[rid] = cand
                            incumbent = cand
                            # the clip leg's rewrite is now the winner, so the
                            # record is no longer the pre-correction one
                            reverted_ids.discard(rid)
                            stats['clip5_accepted'] += 1

                if rid in clip3_ids:
                    stats['clip3_leg'] += 1
                    cand = _clip_copy(incumbent)
                    proposed, refused_delta = _terminal_clip_leg(
                        cand, genome, index, rcfg, rstats)
                    stats['clip3_refused'] += refused_delta
                    if proposed:
                        stats['clip3_proposed'] += 1
                        realigned_ids.add(rid)
                        if reentry_accept(incumbent, cand, genome, penalty_table):
                            accepted[rid] = cand
                            reverted_ids.discard(rid)
                            stats['clip3_accepted'] += 1
        logger.info(
            'triage clip legs: 5p leg=%d (refused=%d proposed=%d accepted=%d) '
            '3p leg=%d (refused=%d proposed=%d accepted=%d)',
            stats['clip5_leg'], stats['clip5_refused'],
            stats['clip5_proposed'], stats['clip5_accepted'],
            stats['clip3_leg'], stats['clip3_refused'],
            stats['clip3_proposed'], stats['clip3_accepted'],
        )

    # ── Final BAM: originals, with accepted replacements swapped in ─────────
    replaced: Set[str] = set()
    with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam, \
            pysam.AlignmentFile(output_bam, 'wb', header=bam.header) as out:
        for read in bam.fetch(until_eof=True):
            rid = read.query_name
            if (_classifiable(read) and rid in accepted and rid not in replaced):
                out.write(accepted[rid])
                replaced.add(rid)
            else:
                out.write(read)

    if sort_and_index:
        sorted_path = output_bam + '.sorted.tmp'
        pysam.sort('-o', sorted_path, output_bam)
        os.replace(sorted_path, output_bam)
        pysam.index(output_bam)

    rows = []
    for rid, r in results.items():
        rows.append({
            'read_id': rid,
            'label': r.label,
            'reasons': ','.join(r.reasons),
            'junction_proximal_errors': r.junction_proximal_errors,
            'clip_5p': r.clip_5p,
            'clip_3p': r.clip_3p,
            'n_junctions': r.n_junctions,
            'n_unannotated': r.n_unannotated,
            'realigned': rid in realigned_ids,
            'accepted': rid in accepted,
            # True when the record that won is the PRE-CORRECTION alignment
            # (ISSUE-009) — the metric for "how much 2F/2H damage did Station B
            # undo", which cannot be read off `accepted` alone.
            'reverted_to_original': rid in reverted_ids,
        })
    return rows, stats
