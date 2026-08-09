"""
triage.py — read-evidence triage of consensus alignments + targeted re-alignment.

The consensus-triage layer of the Rectify Re-aligner
(design: dev/DENOVO_ALIGNER_REASSESSMENT_20260809.md §3):

  1. **Classify** each consensus/corrected alignment from READ EVIDENCE only —
     junction-proximal mismatch/indel weight (``_count_junction_proximity_errors``),
     unexplained 5'/3' soft clips, and junction annotation status. Clean
     alignments are labelled ``high_confidence`` and BYPASSED (the do-no-harm
     data say they are at ceiling: touching them is downside + compute).
  2. **Re-align** only the triaged minority, one leg per failure mode. This
     module implements the JUNCTION leg: motif-blind refinement
     (``refine_bam_junctions(motif_blind=True)``, compensating-indel invariant
     always-on) over a junction pool built from the FULL BAM's evidence.
     The 5'/3' soft-clip legs route to the existing clip-rescue machinery
     (Cat3 / ``overhang_resolver``) and are NOT wired here yet — soft-clip
     triage reads keep their label for those legs.
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


def triage_realign_bam(
    input_bam: str,
    output_bam: str,
    genome: Dict[str, str],
    annotated_junctions: Set[Tuple],
    policy: Optional[TriagePolicy] = None,
    pool_bams: Optional[List[str]] = None,
    penalty_table=None,
    realign: bool = True,
    refine_kwargs: Optional[dict] = None,
    sort_and_index: bool = True,
) -> Tuple[List[dict], Dict[str, int]]:
    """Classify a consensus BAM; motif-blind-re-align the triaged junction reads.

    Pipeline (see module docstring):
      classify → junction pool from the FULL evidence (``pool_bams`` or the
      input BAM itself) → ``refine_bam_junctions(motif_blind=True)`` on ONLY
      the triaged junction-leg reads → hp_ed re-entry per read → final BAM.

    Returns ``(rows, stats)`` where ``rows`` are per-read TSV dicts and
    ``stats`` counts {classified, high_confidence, triaged, junction_leg,
    realigned, accepted}.
    """
    from rectify.core.splice.junction_refiner import refine_bam_junctions
    from rectify.core.splice.junction_scoring import build_junction_pool

    policy = policy or TriagePolicy()
    annotated_3 = _annotated_3tuples(annotated_junctions)

    # ── Pass 1: classify ────────────────────────────────────────────────────
    results: Dict[str, TriageResult] = {}
    n_dup = 0
    with pysam.AlignmentFile(input_bam, 'rb', check_sq=False) as bam:
        header = bam.header.to_dict()
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
    }

    accepted: Dict[str, pysam.AlignedSegment] = {}
    realigned_ids: Set[str] = set()

    if realign and junction_ids:
        with tempfile.TemporaryDirectory(prefix='rectify_triage_') as tmp:
            triaged_bam = os.path.join(tmp, 'triaged.bam')
            refined_bam = os.path.join(tmp, 'refined.bam')

            originals: Dict[str, pysam.AlignedSegment] = {}
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
        })
    return rows, stats
