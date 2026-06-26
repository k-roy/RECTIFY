#!/usr/bin/env python3
"""Ambiguity-aware per-junction + per-indel accuracy scorer (component 3).

Consumes an aligner/consensus BAM + the per-read truth table (component 2) and
emits TP/FP/FN per aligner, stratified by junction class (annotated/NIC/NNC) AND
read-class, plus the **framing metric**: EXACT INDEL-POSITION CONCORDANCE WITH
TRUTH, ambiguity-aware, NEVER edit distance.

Ambiguity-awareness reuses the in-tree primitives
(``chimeric_consensus.normalize_junction`` / ``junction_ambiguity_window`` /
``_canonical_within_window``, ``chimeric_consensus.py:59-155``) so that:

* a junction called one bp into a donor/acceptor repeat is NOT charged FP — both
  truth and call are left-normalized to the leftmost ambiguity-equivalent
  coordinate, so set membership IS the ambiguity-aware match (the exact trap that
  produced the GMAP 0.09 artifact);
* an HP/STR indel placed anywhere within its run is credited as a TP — the net
  (D-I) over the truth ``[eq_start, eq_end)`` span is compared to the truth net,
  not the literal edit position (edit distance ties by construction).

Mirrors ``validate_command.py``'s CPA-truth scoring pattern, lifted to junctions
+ indels (it does not score junctions today).

Author: Kevin R. Roy
"""
from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Iterable

import pysam

from ..consensus.chimeric_consensus import (
    normalize_junction,
    _canonical_within_window,
    junction_ambiguity_window,
)
from .truth_schema import (
    ReadTruth,
    JunctionTruth,
    IndelTruth,
    JunctionClass,
    IndelKind,
    read_truth_table,
)

logger = logging.getLogger(__name__)

_REF_CONSUMING = {0, 2, 3, 7, 8}   # M D N = X
_QUERY_CONSUMING = {0, 1, 4, 7, 8}  # M I S = X


# ---------------------------------------------------------------------------
# CIGAR feature extraction
# ---------------------------------------------------------------------------
def extract_junctions(ref_start: int, cigartuples) -> List[Tuple[int, int]]:
    """Return intron ``[start, end)`` spans (N-ops) in genome coords."""
    out = []
    pos = ref_start
    for op, ln in cigartuples or []:
        if op == 3:  # N
            out.append((pos, pos + ln))
            pos += ln
        elif op in _REF_CONSUMING:
            pos += ln
    return out


def net_indel_in_span(ref_start: int, cigartuples,
                      span_start: int, span_end: int) -> Tuple[int, int]:
    """Net ``(deletion - insertion)`` base count whose REF position lies within
    ``[span_start, span_end)``, and the indel base count OUTSIDE the span.

    Deletions consume ref (their bases overlap the span by genome position);
    insertions are charged to the span if the insertion point (the current ref
    position) lies within ``[span_start, span_end]``. Mirrors the vertical-slice
    ``net_indel_in_run`` ambiguity-aware match generalized to a named span.
    """
    pos = ref_start
    in_span = 0
    out_span = 0
    for op, ln in cigartuples or []:
        if op in (0, 7, 8):            # M/=/X
            pos += ln
        elif op == 2:                  # D consumes ref
            overlap = max(0, min(pos + ln, span_end) - max(pos, span_start))
            in_span += overlap
            out_span += (ln - overlap)
            pos += ln
        elif op == 1:                  # I consumes query, ref pos fixed
            if span_start <= pos <= span_end:
                in_span -= ln
            else:
                out_span += ln
        elif op == 3:                  # N intron — not an indel
            pos += ln
        # S/H clipping ignored
    return in_span, out_span


def all_indel_positions(ref_start: int, cigartuples) -> List[Tuple[int, int, int]]:
    """Return every indel as ``(ref_pos, length, kind)`` where kind 2=D, 1=I.
    ``ref_pos`` is the genome position at the start of the op."""
    pos = ref_start
    out = []
    for op, ln in cigartuples or []:
        if op == 2:
            out.append((pos, ln, 2))
            pos += ln
        elif op == 1:
            out.append((pos, ln, 1))
        elif op in _REF_CONSUMING:
            pos += ln
    return out


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------
@dataclass
class JunctionScore:
    tp: int = 0
    fp: int = 0
    fn: int = 0
    # stratified TP/FP/FN keyed by class label ("ANNOTATED"/"NIC"/"NNC") and by
    # canonicity ("canonical"/"noncanonical")
    by_class: Dict[str, Dict[str, int]] = field(default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    by_canon: Dict[str, Dict[str, int]] = field(default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    fp_variant_adjacent: int = 0   # C6: false junctions near a truth variant

    @property
    def recall(self) -> float:
        d = self.tp + self.fn
        return self.tp / d if d else 0.0

    @property
    def precision(self) -> float:
        d = self.tp + self.fp
        return self.tp / d if d else 0.0

    @property
    def fdr(self) -> float:
        d = self.tp + self.fp
        return self.fp / d if d else 0.0


@dataclass
class IndelScore:
    """Position-exact (ambiguity-aware) indel concordance — the framing metric."""
    correct: int = 0          # net-in-span matched truth net (TP)
    incorrect: int = 0        # truth indel present but net mismatched (FN)
    false_indels: int = 0     # reads/positions with an indel outside every truth span
    clean_reads: int = 0      # k=0 reads scored for the false-indel-rate control
    clean_reads_with_false_indel: int = 0
    # C1 cell accounting: keyed by (base_class, run_copies, context)
    by_cell: Dict[Tuple[str, int, str], Dict[str, int]] = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(int)))

    @property
    def position_exact_concordance(self) -> float:
        d = self.correct + self.incorrect
        return self.correct / d if d else 0.0

    @property
    def false_indel_rate(self) -> float:
        return (self.clean_reads_with_false_indel / self.clean_reads
                if self.clean_reads else 0.0)


@dataclass
class AlignerScore:
    aligner: str
    reads_scored: int = 0
    reads_placed: int = 0          # present + mapped in this BAM
    junction: JunctionScore = field(default_factory=JunctionScore)
    indel: IndelScore = field(default_factory=IndelScore)
    cpa_abs_errors: List[int] = field(default_factory=list)

    def summary(self) -> Dict:
        cpa = sorted(self.cpa_abs_errors)
        median_cpa = cpa[len(cpa) // 2] if cpa else None
        return {
            "aligner": self.aligner,
            "reads_scored": self.reads_scored,
            "reads_placed": self.reads_placed,
            "junction_recall": round(self.junction.recall, 4),
            "junction_precision": round(self.junction.precision, 4),
            "junction_fdr": round(self.junction.fdr, 4),
            "junction_tp": self.junction.tp,
            "junction_fp": self.junction.fp,
            "junction_fn": self.junction.fn,
            "junction_fp_variant_adjacent": self.junction.fp_variant_adjacent,
            "junction_by_class": {k: dict(v) for k, v in self.junction.by_class.items()},
            "junction_by_canonicity": {k: dict(v) for k, v in self.junction.by_canon.items()},
            "indel_position_exact_concordance": round(self.indel.position_exact_concordance, 4),
            "indel_correct": self.indel.correct,
            "indel_incorrect": self.indel.incorrect,
            "indel_false_indel_rate": round(self.indel.false_indel_rate, 4),
            "indel_clean_reads": self.indel.clean_reads,
            "median_abs_cpa_error": median_cpa,
        }


# ---------------------------------------------------------------------------
# Genome loading (light FASTA reader — pysam.FastaFile preferred)
# ---------------------------------------------------------------------------
def load_genome(fasta_path: Union[str, Path]) -> Dict[str, str]:
    fa = pysam.FastaFile(str(fasta_path))
    return {c: fa.fetch(c) for c in fa.references}


# ---------------------------------------------------------------------------
# Core scorer
# ---------------------------------------------------------------------------
def _normalized_truth_juncs(rt: ReadTruth, genome: Dict[str, str]) -> Dict[Tuple[int, int], JunctionTruth]:
    """Truth junctions keyed by (already-normalized) intron coords."""
    return {(j.intron_start, j.intron_end): j for j in rt.junctions}


def _variant_positions(rt: ReadTruth) -> List[int]:
    return [v.pos for v in rt.variants]


def score_bam(bam_path: Union[str, Path],
              truth: Union[Dict[str, ReadTruth], List[ReadTruth]],
              genome: Dict[str, str],
              aligner_name: Optional[str] = None,
              variant_adjacency_bp: int = 10) -> AlignerScore:
    """Score one BAM against the truth table. Returns an ``AlignerScore``.

    Junctions: each aligner N-op is left-normalized (``normalize_junction``) and
    matched against the (normalized) truth set — ambiguity-aware TP. Unmatched
    aligner junctions are FP (tagged variant-adjacent if within
    ``variant_adjacency_bp`` of a truth variant); unmatched truth junctions FN.

    Indels: per truth indel, ``net_indel_in_span`` over ``[eq_start, eq_end)`` is
    compared to the truth net (DEL=+L, INS=-L) → position-exact TP/FN. A
    clean/no-error read with any indel outside every truth span counts toward the
    false-indel rate.
    """
    if isinstance(truth, list):
        truth = {t.read_id: t for t in truth}
    name = aligner_name or Path(bam_path).stem
    score = AlignerScore(aligner=name)

    placed_ids = set()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            rt = truth.get(read.query_name)
            if rt is None:
                continue
            if read.is_unmapped:
                continue
            placed_ids.add(read.query_name)
            chrom = read.reference_name
            seq = genome.get(chrom, "")
            _score_read(score, read, rt, seq, variant_adjacency_bp)

    # FN bookkeeping: truth reads NOT placed in this BAM contribute all their
    # truth junctions as FN (the aligner failed to place => missed every junction)
    score.reads_placed = len(placed_ids)
    for rid, rt in truth.items():
        score.reads_scored += 1
        if rid not in placed_ids:
            for j in rt.junctions:
                score.junction.fn += 1
                score.junction.by_class[j.klass.value]["fn"] += 1
                score.junction.by_canon["canonical" if j.canonical else "noncanonical"]["fn"] += 1
            for ind in rt.indels:
                score.indel.incorrect += 1
    return score


def _score_read(score: AlignerScore, read, rt: ReadTruth, genome_seq: str,
                variant_adjacency_bp: int) -> None:
    cig = read.cigartuples
    rstart = read.reference_start

    # ---- Junctions -------------------------------------------------------
    truth_set = {(j.intron_start, j.intron_end) for j in rt.junctions}
    truth_by_coord = {(j.intron_start, j.intron_end): j for j in rt.junctions}
    called = extract_junctions(rstart, cig)
    matched_truth = set()
    var_pos = _variant_positions(rt)
    for (cs, ce) in called:
        ns, ne = normalize_junction(cs, ce, genome_seq) if genome_seq else (cs, ce)
        if (ns, ne) in truth_set:
            j = truth_by_coord[(ns, ne)]
            score.junction.tp += 1
            score.junction.by_class[j.klass.value]["tp"] += 1
            score.junction.by_canon["canonical" if j.canonical else "noncanonical"]["tp"] += 1
            matched_truth.add((ns, ne))
        else:
            score.junction.fp += 1
            # canonicity of the FALSE call (for the non-canonical FDR track)
            if genome_seq:
                l_amb, r_amb = junction_ambiguity_window(ns, ne, genome_seq)
                canon = _canonical_within_window(ns, ne, genome_seq, l_amb, r_amb)
            else:
                canon = False
            score.junction.by_canon["canonical" if canon else "noncanonical"]["fp"] += 1
            score.junction.by_class["FP_NOVEL"]["fp"] += 1
            if var_pos and min(abs(ns - vp) for vp in var_pos) <= variant_adjacency_bp:
                score.junction.fp_variant_adjacent += 1
    for (ts, te) in truth_set - matched_truth:
        j = truth_by_coord[(ts, te)]
        score.junction.fn += 1
        score.junction.by_class[j.klass.value]["fn"] += 1
        score.junction.by_canon["canonical" if j.canonical else "noncanonical"]["fn"] += 1

    # ---- Indels (the framing metric) ------------------------------------
    truth_spans: List[Tuple[int, int]] = []
    is_clean_read = True
    for ind in rt.indels:
        truth_net = ind.length if ind.kind == IndelKind.DEL else -ind.length
        if truth_net != 0:
            is_clean_read = False
        truth_spans.append((ind.eq_start, ind.eq_end))
        in_span, _out = net_indel_in_span(rstart, cig, ind.eq_start, ind.eq_end)
        cell = (ind.base_class, ind.run_copies, ind.context)
        if in_span == truth_net:
            score.indel.correct += 1
            score.indel.by_cell[cell]["correct"] += 1
        else:
            score.indel.incorrect += 1
            score.indel.by_cell[cell]["incorrect"] += 1

    # false-indel control: indel bases the aligner introduced OUTSIDE every truth
    # span. A read with NO truth indel (clean / k=0) is the FP control.
    if is_clean_read and not any(net != 0 for net in
                                 [(ind.length if ind.kind == IndelKind.DEL else -ind.length)
                                  for ind in rt.indels]):
        score.indel.clean_reads += 1
        false_found = False
        for (ipos, ilen, ikind) in all_indel_positions(rstart, cig):
            inside = any(s <= ipos < e for (s, e) in truth_spans)
            if not inside:
                false_found = True
                break
        if false_found:
            score.indel.clean_reads_with_false_indel += 1
            score.indel.false_indels += 1

    # ---- CPA -------------------------------------------------------------
    if rt.true_cpa is not None:
        # 3' end genome coord per RECTIFY convention
        if read.is_reverse:
            est = read.reference_start
        else:
            est = read.reference_end - 1
        score.cpa_abs_errors.append(abs(est - rt.true_cpa))


def score_panel(bam_paths: Dict[str, Union[str, Path]],
                truth_path: Union[str, Path],
                genome_path: Union[str, Path]) -> Dict[str, Dict]:
    """Score a panel of per-aligner BAMs against truth. Returns a dict of
    ``{aligner: summary}`` plus a ``_which_aligners_placed`` per-read map (so
    'panel-unplaced' is decidable, §9)."""
    truth = read_truth_table(truth_path)
    truth_map = {t.read_id: t for t in truth}
    genome = load_genome(genome_path)
    summaries = {}
    placed_by_aligner: Dict[str, set] = {}
    for aln, bp in bam_paths.items():
        s = score_bam(bp, truth_map, genome, aligner_name=aln)
        summaries[aln] = s.summary()
        with pysam.AlignmentFile(str(bp), "rb") as bam:
            placed = {r.query_name for r in bam
                      if not r.is_unmapped and not r.is_secondary and not r.is_supplementary}
        placed_by_aligner[aln] = placed & set(truth_map)

    which_placed = {}
    panel_unplaced = []
    for rid in truth_map:
        plc = [a for a, s in placed_by_aligner.items() if rid in s]
        which_placed[rid] = plc
        if not plc:
            panel_unplaced.append(rid)
    summaries["_which_aligners_placed"] = which_placed
    summaries["_panel_unplaced_reads"] = panel_unplaced
    summaries["_panel_unplaced_fraction"] = round(
        len(panel_unplaced) / len(truth_map), 4) if truth_map else 0.0
    return summaries
