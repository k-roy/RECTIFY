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
from collections import Counter, defaultdict
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

import numpy as _np
import pysam

from ...utils.genome import standardize_chrom_name
from .hp_penalty import (
    HpPenaltyTable,
    _NUMBA_AVAILABLE,
    _hp_run_length,
    _precompute_del_costs,
    _precompute_refcol_ins_costs,
    _score_hp_dp_numba,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# EXPERIMENTAL FLAG — full-run (cut-independent) insertion cost
# ---------------------------------------------------------------------------
# When True, _score_junction precomputes per-position insertion costs ONCE on the
# FULL rescue window (using each base's homopolymer run length measured over the
# whole rescue substring) and threads them, cut-INDEPENDENTLY, into every per-k
# _score_hp_anchored call — analogous to how del_costs_fwd/rev are already
# precomputed on the fixed reference windows.
#
# Default False -> byte-identical to the legacy per-cut behaviour (ins_costs are
# recomputed on the per-k TRUNCATED substring inside _score_hp_anchored). See
# dev/INSCOST_INVESTIGATION.md for the fabrication witness, call-change rate, and
# guard-interaction analysis. The env var lets a subprocess opt in with a fresh
# fork (so a reused worker pool cannot silently retain the old value).
_USE_FULL_RUN_INS: bool = os.environ.get("RECTIFY_FULL_RUN_INS", "") == "1"


# ---------------------------------------------------------------------------
# EXPERIMENTAL FLAG — reference-column (cut-independent) insertion cost
# ---------------------------------------------------------------------------
# When True, _score_junction precomputes per-GAP insertion costs ONCE on the two
# FIXED reference windows (exon2-start and intron-end) using the GENOME/REFERENCE
# homopolymer run length at each aligned column — exactly mirroring how
# _precompute_del_costs already indexes DELETION cost by the genome HP run. An
# insertion aligns to NO reference base, so its cost is charged against the
# reference context at the DP COLUMN (gap) where the insertion fires:
# hp_len = max(run(genome[gp-1]), run(genome[gp])) — the same axis the penalty
# table was CALIBRATED on (profiler Phase 5, empirical_cigar_error_profiler.py).
#
# This makes the ins-cost vector length R+1 (one per DP gap, not per query base),
# a function of the fixed reference window + absolute genome position ONLY — so it
# is computed ONCE outside the k-loop and passed UNCHANGED to every k (no slicing).
# That is the cut-independence proof (stronger than full-run, which had to slice a
# query-indexed vector). Like full-run it unlocks the single-pass concat DP, but it
# charges the over-call against the run the read aligns TO (the genome), not the
# read's own possibly-over-called run — the physically/calibration-correct axis.
#
# Default False -> byte-identical to the legacy per-cut behaviour. Refcol takes
# PRECEDENCE over _USE_FULL_RUN_INS when both are set (mutually exclusive semantics;
# they index different axes). See dev/INSCOST_REFCOL_BUILD.md and
# dev/INSCOST_AUDIT_model-correctness.md (KEY FINDING #3, the reference-column axis).
_USE_REFCOL_INS: bool = os.environ.get("RECTIFY_REFCOL_INS", "") == "1"


# ---------------------------------------------------------------------------
# PERF FLAG — vectorized (concat) DP for the table-free re-placer
# ---------------------------------------------------------------------------
# _score_junction's bilateral k-sweep does up to 2L (~60) _score_hp_anchored DP
# calls per candidate. When penalty_table is None the insertion cost is flat, so
# the per-cut ins truncation that blocks the general concat DP does NOT exist —
# ALL t1(k) and ALL t2(k) collapse into ONE vectorized DP pass each (the query-
# suffix reversal trick in _all_suffix_scores), reproducing min_k[t1(k)+t2(k)]
# EXACTLY. ~14x fewer cell-ops, verified byte-identical (8000/8000, table=None).
# GATED strictly on penalty_table is None (the flat-ins path) — the ONLY config
# where the fast path is provably identical (MECH2 per-cut ins vanishes; costs are
# exactly float-representable so MECH3 FP tolerance vanishes; MECH1 boundary is
# handled by k in [0,L)). The shipped native re-aligner (motif-blind + guard) runs
# table-free, so this covers it. Default False. See dev/PERF_PIVOT_TABLEFREE_CONCAT.md.
_USE_CONCAT_DP: bool = os.environ.get("RECTIFY_CONCAT_DP", "1") != "0"  # DEFAULT ON (2026-07-09); RECTIFY_CONCAT_DP=0 forces legacy (A/B)

# Flat (penalty_table=None) DP cost constants — the SINGLE source of truth shared by
# _score_hp_anchored and the vectorized _all_suffix_scores fast path. Byte-identity of
# the _USE_CONCAT_DP fast path depends on both using IDENTICAL flat costs; factoring
# them here prevents a future retune of one path from silently breaking identity
# (auditor recommendation, dev/CONCAT_DP_TABLEFREE_AUDIT.md). All dyadic (exactly
# float-representable) — the basis of the 0-ULP guarantee.
_FLAT_SUB: float = float(os.environ.get("RECTIFY_FLAT_SUB", "1.0"))
_FLAT_DEL_NORMAL: float = float(os.environ.get("RECTIFY_FLAT_DEL_NORMAL", "1.0"))
_FLAT_DEL_HP: float = float(os.environ.get("RECTIFY_FLAT_DEL_HP", "0.5"))
_FLAT_INS: float = float(os.environ.get("RECTIFY_FLAT_INS", "1.25"))


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
# EXON-flank op set: reference-consuming ops that belong to an ALIGNED block.
# N is deliberately EXCLUDED — the flank walkers in
# junction_refiner._apply_junction_replacement (find/consume the exon ops that
# abut an N-op) must stop at an intron, never walk through it.
_REF_CONSUMING:   FrozenSet[int] = frozenset([_M, _D, _EQ, _X])
# POSITION-tracking op set: every op that advances the genomic cursor per the SAM
# spec, N included. Any walker that maintains a reference position while stepping
# through a CIGAR MUST use this set. Using _REF_CONSUMING for a cursor reports
# every N-op after the first short by the summed length of all preceding introns
# (the ISSUE-004 defect: 81/81 multi-intron human panel reads wrong, 85 % of the
# observed junction pool phantom). Yeast hid it — ~95 % of intron-bearing yeast
# genes have a single intron.
_REF_CONSUMING_POS: FrozenSet[int] = _REF_CONSUMING | frozenset([_N])
_MATCH_OPS:       FrozenSet[int] = frozenset([_M, _EQ])  # consume query+ref, no indel

# Minimum clean exon overhang (bp) on BOTH flanks for an observed N-op to seed a
# pool junction. Aligners (notably gapmm2) occasionally emit a tiny-anchor
# spurious split — e.g. 4M4250N223M, a 4 bp anchor before a 4.25 kb "intron" —
# which then propagates to many cleanly-aligned reads during 3'SS rescue. A real
# junction is crossed by reads with a substantial clean anchor, so requiring one
# here keeps real junctions while dropping the artifacts. Kevin 2026-05-24.
DEFAULT_MIN_JUNCTION_ANCHOR = 10

# Aligner → independent algorithm family. The pool concordance relaxation
# (build_junction_pool) treats agreement ACROSS distinct families as independent
# evidence that a short-anchor junction is real, but agreement WITHIN a family is
# not: gapmm2 wraps minimap2 (mappy + edlib), so a systematic minimap2 mis-seed at
# a genomic A-tract is replicated by gapmm2 — counting them as two votes would
# re-admit the very artifacts the anchor floor drops. minimap2/gapmm2 = one family;
# uLTRA (annotation-graph collinear chaining), mapPacBio (BBMap), and deSALT (RdBG)
# are each algorithmically independent. Kevin 2026-05-25.
ALIGNER_FAMILY: Dict[str, str] = {
    "minimap2": "minimap2",
    "gapmm2": "minimap2",      # minimap2 wrapper (mappy) + edlib terminal refinement
    "overhang_resolver": "minimap2",  # re-placement of the minimap2 arm's clips
                                      # (planning/641) — NOT an independent vote
    "ultra": "ultra",          # annotation-guided collinear chaining (MEM)
    "mappacbio": "bbmap",      # BBMap long-read mode
    "desalt": "desalt",        # de Bruijn graph / deBGA (RdBG)
}
DEFAULT_RELAX_MIN_FAMILIES = 2


def _family_from_bam_path(path: str) -> str:
    """Infer the algorithm family from a per-aligner BAM filename.

    Matches a known aligner token anywhere in the basename (robust to suffixes
    like ``.namesorted``). An unrecognised name maps to the shared ``"?"``
    sentinel so it cannot satisfy the cross-family relaxation on its own — an
    unknown panel falls back to the pure hard floor rather than relaxing
    unsoundly.
    """
    base = os.path.basename(path).lower()
    for aln, fam in ALIGNER_FAMILY.items():
        if aln in base:
            return fam
    return "?"


def _resolve_families(
    aligner_bams: List[str],
    aligner_families: Optional[List[str]],
) -> List[str]:
    """Per-BAM family list: explicit ``aligner_families`` if given (mapped through
    ``ALIGNER_FAMILY``), else inferred from each BAM filename."""
    if aligner_families is not None and len(aligner_families) == len(aligner_bams):
        return [ALIGNER_FAMILY.get(a.lower(), a.lower()) for a in aligner_families]
    return [_family_from_bam_path(bp) for bp in aligner_bams]


def _merge_del_into_intron(
    cigartuples: List[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """Collapse any run of adjacent D/N ops containing an N into a single N op.

    D (deletion) and N (intron skip) both consume reference and not query, so an
    aligner's ``1D112N`` is geometrically identical to ``113N`` — the split is a
    representational quirk that leaves a phantom deletion abutting the junction.
    RECTIFY treats them as one intron: this normalises ``...M 1D 112N M...`` to
    ``...M 113N M...`` so junction extraction and the anchor check see the real
    exon match adjacent to the junction. Pure-deletion runs (no N) are untouched.
    """
    if not cigartuples:
        return cigartuples
    out: List[Tuple[int, int]] = []
    i, n = 0, len(cigartuples)
    while i < n:
        op = cigartuples[i][0]
        if op in (_D, _N):
            run_len = 0
            has_n = False
            j = i
            while j < n and cigartuples[j][0] in (_D, _N):
                run_len += cigartuples[j][1]
                has_n = has_n or cigartuples[j][0] == _N
                j += 1
            if has_n:
                out.append((_N, run_len))
            else:
                out.extend(cigartuples[i:j])  # pure deletion(s) — keep as-is
            i = j
        else:
            out.append(cigartuples[i])
            i += 1
    return out


def _is_real_nucleotide_seq(seq: str) -> bool:
    """True if seq is composed only of A/C/G/T/N (real bases, not '=' encoding)."""
    return bool(seq) and all(b in "ACGTNacgtn" for b in seq)


def _is_low_complexity_anchor(seq: str) -> bool:
    """True if seq is a period-1/2/3 repeat (homopolymer, di- or tri-nucleotide).

    A short anchor over low-complexity sequence matches spuriously in many
    places, so it should not be trusted to define a junction boundary.
    """
    if not seq:
        return True
    s = seq.upper()
    n = len(s)
    for period in (1, 2, 3):
        if n >= period and all(s[i] == s[i % period] for i in range(n)):
            return True
    return False


def _seq_periodicity(seq: str) -> float:
    """Graded low-complexity score: fraction of bases consistent with the best
    period-1/2/3 repeat. 1.0 = a perfect homopolymer / di- / tri-nucleotide
    repeat (matches ``_is_low_complexity_anchor``); diverse sequence sits near
    0.4–0.5. Motif-agnostic — measures repeat structure, never splice motifs.
    """
    s = seq.upper()
    n = len(s)
    if n == 0:
        return 1.0
    best = 0.0
    for period in (1, 2, 3):
        if n < period:
            continue
        m = sum(1 for i in range(n) if s[i] == s[i % period])
        f = m / n
        if f > best:
            best = f
    return best


def _junction_worst_flank_periodicity(
    genome: Dict[str, str],
    chrom: str,
    intron_start: int,
    intron_end: int,
    window: int = 10,
    reach: int = 30,
) -> float:
    """Worst-flank low-complexity score for a junction's genomic anchors.

    For each exon flank (``reach`` bp immediately abutting the intron on the
    upstream and downstream sides), slide a ``window``-bp window and take the
    MINIMUM periodicity over windows — i.e. the most-complex window that flank
    offers. Then return the MAX over the two flanks: a junction is well-anchored
    only if BOTH flanks contain at least one complex (uniquely-placeable) window,
    so the worst flank governs. A genomic A-tract flank yields ~1.0; a flank with
    any diverse 10-mer within ``reach`` yields its complex window's low score
    (this is what lets ``...AGCTGAGTC|AAAA...`` pass — the diverse prefix anchors
    it). Returns 0.0 when the genome/chrom is unavailable (cannot assess → caller
    must not flag on this dimension).
    """
    chrom_seq = genome.get(chrom) if genome else None
    if not chrom_seq:
        return 0.0
    L = len(chrom_seq)

    def _flank_min(seq: str) -> float:
        seq = seq.upper()
        if len(seq) < window:
            return _seq_periodicity(seq) if seq else 1.0
        return min(_seq_periodicity(seq[i:i + window])
                   for i in range(len(seq) - window + 1))

    left = chrom_seq[max(0, intron_start - reach):intron_start]
    right = chrom_seq[intron_end:min(L, intron_end + reach)]
    return max(_flank_min(left), _flank_min(right))


def _junction_anchor_ok(
    cigar: List[Tuple[int, int]],
    idx: int,
    qpos: int,
    query_seq: str,
    min_overhang: int,
) -> bool:
    """Whether the N-op at ``cigar[idx]`` has a trustworthy exon anchor on both
    flanks.

    ``qpos`` is the query coordinate at the N-op (N consumes no query, so it is
    both the end of the left op and the start of the right op). Requires the
    immediately-flanking ops to be match ops (M/=) of length >= ``min_overhang``
    — which guarantees no indel within the anchor, and for =/X aligners no
    mismatch either — and the ``min_overhang`` anchor bases on each side to be
    non-low-complexity.
    """
    if min_overhang <= 0:
        return True
    if idx == 0 or idx == len(cigar) - 1:
        return False  # junction at a read edge has no flanking exon
    lop, llen = cigar[idx - 1]
    rop, rlen = cigar[idx + 1]
    if lop not in _MATCH_OPS or llen < min_overhang:
        return False
    if rop not in _MATCH_OPS or rlen < min_overhang:
        return False
    if not query_seq:
        return True  # length/indel gate already passed; cannot complexity-check
    left_anchor = query_seq[qpos - min_overhang:qpos]
    right_anchor = query_seq[qpos:qpos + min_overhang]
    # Only judge complexity on real nucleotide sequence. Some BAMs store SEQ as
    # '=' (reference-match encoding); pysam returns literal '=' chars, which are
    # not assessable for complexity — fall back to the length/indel gate.
    for anchor in (left_anchor, right_anchor):
        if _is_real_nucleotide_seq(anchor) and _is_low_complexity_anchor(anchor):
            return False
    return True


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

    .. warning::
       Before the ISSUE-004 fix this walker advanced the reference cursor with
       ``_REF_CONSUMING`` (no N), so on a multi-intron read every junction after
       the first was reported short by the summed length of the preceding
       introns. **Any ``rectify prescan`` junction-pool pickle written before
       that fix is contaminated with those phantom coordinates and must be
       rebuilt** — it cannot be repaired in place.
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
                for op, length in _merge_del_into_intron(read.cigartuples):
                    if op == _N:
                        if max_junction_size is None or length <= max_junction_size:
                            junctions.add((chrom, pos, pos + length))
                    if op in _REF_CONSUMING_POS:
                        pos += length
    except Exception as exc:
        logger.warning("collect_junctions_from_bam(%s): %s", bam_path, exc)
    return junctions


def _collect_junction_counts_core(
    bam_path: str,
    chrom_filter: Optional[str] = None,
    max_junction_size: Optional[int] = None,
    min_anchor_overhang: int = DEFAULT_MIN_JUNCTION_ANCHOR,
) -> Tuple[Counter, Counter]:
    """One pass over a BAM → ``(anchor_pass_counts, raw_counts)``.

    ``anchor_pass_counts`` counts only N-ops with a trustworthy exon anchor on
    both flanks (see ``_junction_anchor_ok``). ``raw_counts`` counts EVERY
    observed (size-filtered) N-op junction regardless of anchor quality, so the
    pool builder can apply the cross-aligner concordance relaxation: a real
    junction whose every supporting read has a short anchor (e.g. a short
    exon-1) is kept when independently called by multiple aligners, while a
    single-aligner tiny-anchor artifact is dropped.
    """
    anchor: Counter = Counter()
    raw: Counter = Counter()
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
                cigar = _merge_del_into_intron(read.cigartuples)
                query_seq = read.query_sequence or ""
                pos = read.reference_start
                qpos = 0
                for idx, (op, length) in enumerate(cigar):
                    if op == _N and (max_junction_size is None
                                     or length <= max_junction_size):
                        j = (chrom, pos, pos + length)
                        raw[j] += 1
                        if _junction_anchor_ok(
                                cigar, idx, qpos, query_seq, min_anchor_overhang):
                            anchor[j] += 1
                    if op in _REF_CONSUMING_POS:
                        pos += length
                    if op in _QUERY_CONSUMING:
                        qpos += length
    except Exception as exc:
        logger.warning("_collect_junction_counts_core(%s): %s", bam_path, exc)
    return anchor, raw


def collect_junction_counts_from_bam(
    bam_path: str,
    chrom_filter: Optional[str] = None,
    max_junction_size: Optional[int] = None,
    min_anchor_overhang: int = DEFAULT_MIN_JUNCTION_ANCHOR,
) -> Counter[Junction]:
    """Count N-op intervals from a BAM file (anchor-passing junctions only).

    Only N-ops with a trustworthy exon anchor on both flanks (see
    ``_junction_anchor_ok``) are counted, so a tiny-anchor spurious split does
    not seed a pool junction. ``min_anchor_overhang=0`` disables the filter.
    """
    return _collect_junction_counts_core(
        bam_path, chrom_filter, max_junction_size, min_anchor_overhang)[0]


def build_junction_pool(
    aligner_bams: List[str],
    annotated_junctions: Set[Tuple],
    chrom_filter: Optional[str] = None,
    min_observed_support: int = 1,
    max_junction_size: Optional[int] = None,
    min_anchor_overhang: int = DEFAULT_MIN_JUNCTION_ANCHOR,
    relax_min_families: int = DEFAULT_RELAX_MIN_FAMILIES,
    aligner_families: Optional[List[str]] = None,
) -> Tuple[Set[Junction], Set[Junction]]:
    """Build union of annotated + per-aligner junctions.

    An observed junction enters the pool when EITHER:

    * it has well-anchored support — at least ``min_observed_support`` reads
      (summed across aligners) cross it with a clean ≥ ``min_anchor_overhang``
      exon anchor on both flanks (the hard low floor that drops tiny-anchor
      artifacts); OR
    * it has cross-FAMILY concordance — it is independently reported by
      ``>= relax_min_families`` distinct algorithm families (the relaxation
      escape hatch). Agreement is counted across INDEPENDENT families, not raw
      aligners: gapmm2 wraps minimap2, so a minimap2+gapmm2 pair is one vote, not
      two, and a systematic minimap2 mis-seed cannot relax itself in. This
      recovers real junctions whose every supporting read has a short anchor
      (e.g. a short exon-1) when genuinely orthogonal aligners agree, while still
      dropping single-family tiny-anchor splits.

    Args:
        aligner_bams:         Paths to per-aligner BAMs (pre-consensus).
        annotated_junctions:  Set of (chrom, start, end[, strand]) tuples from
                              ``load_annotated_junctions()``.
        chrom_filter:         Optional chromosome to restrict to.
        min_observed_support: Minimum well-anchored read support across all
                              aligner BAMs for an observed junction.
        max_junction_size:    Optional maximum observed N-op length to include.
        min_anchor_overhang:  Clean exon overhang required on both flanks for a
                              read to count toward anchored support.
        relax_min_families:   Number of DISTINCT independent algorithm families
                              that must report a junction for the concordance
                              relaxation to admit it without anchored support.
                              ``<= 0`` disables the relaxation (pure hard floor).
                              The unknown-family sentinel never counts, so an
                              unrecognised panel falls back to the hard floor.
        aligner_families:     Optional per-BAM family labels (same order/length as
                              ``aligner_bams``); mapped through ``ALIGNER_FAMILY``.
                              When omitted, family is inferred from each BAM
                              filename (the ``<sample>.<aligner>.bam`` convention).

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

    # Per-aligner (anchor_pass_counts, raw_counts).
    per_bam: List[Tuple[Counter, Counter]] = []
    if not aligner_bams:
        pass
    elif len(aligner_bams) == 1:
        per_bam.append(_collect_junction_counts_core(
            aligner_bams[0], chrom_filter, max_junction_size, min_anchor_overhang))
    else:
        try:
            from concurrent.futures import ProcessPoolExecutor
            n_workers = min(len(aligner_bams), os.cpu_count() or 4)
            with ProcessPoolExecutor(max_workers=n_workers) as ex:
                futures = [
                    ex.submit(_collect_junction_counts_core, bp, chrom_filter,
                              max_junction_size, min_anchor_overhang)
                    for bp in aligner_bams
                ]
                per_bam = [fut.result() for fut in futures]
        except (OSError, PermissionError, NotImplementedError) as exc:
            logger.warning(
                "build_junction_pool: ProcessPoolExecutor unavailable (%s); "
                "falling back to sequential BAM scans",
                exc,
            )
            per_bam = [
                _collect_junction_counts_core(
                    bp, chrom_filter, max_junction_size, min_anchor_overhang)
                for bp in aligner_bams
            ]

    # Accumulate anchored support (summed reads) and the set of DISTINCT
    # algorithm families reporting each junction (for the concordance relaxation).
    families = _resolve_families(aligner_bams, aligner_families)
    anchor_total: Counter = Counter()
    family_report: Dict[Junction, Set[str]] = defaultdict(set)
    for (anchor_c, raw_c), fam in zip(per_bam, families):
        anchor_total.update(anchor_c)
        for j in raw_c:
            family_report[j].add(fam)

    n_anchored = n_relaxed = 0
    for j in set(anchor_total) | set(family_report):
        if anchor_total[j] >= min_observed_support:
            all_j.add(j)
            n_anchored += 1
        elif relax_min_families > 0:
            # Count distinct REAL families; the "?" unknown sentinel never counts.
            n_fam = len(family_report.get(j, frozenset()) - {"?"})
            if n_fam >= relax_min_families:
                all_j.add(j)
                n_relaxed += 1

    logger.debug(
        "build_junction_pool: %d annotated, %d anchored-support>=%d, "
        "%d family-concordance-relaxed (>=%d families), %d total",
        len(annot_3), n_anchored, min_observed_support,
        n_relaxed, relax_min_families, len(all_j),
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
    sub: float = _FLAT_SUB,
    del_normal: float = _FLAT_DEL_NORMAL,
    del_hp: float = _FLAT_DEL_HP,
    ins: float = _FLAT_INS,
    hp_min_run: int = 4,
    penalty_table: Optional[HpPenaltyTable] = None,
    genome_seq: Optional[str] = None,
    ref_genome_start: int = 0,
    ref_genome_rev: bool = False,
    precomputed_del_costs: Optional[List[float]] = None,
    precomputed_ins_costs: Optional[List[float]] = None,
    precomputed_refcol_ins: Optional[List[float]] = None,
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
        # No ref to match: every query base is an insertion.
        # Reference-column (experimental): the single gap-0 cost (no ref context)
        # applies to every inserted base — cut-independent by construction.
        if precomputed_refcol_ins is not None:
            return float(Q * precomputed_refcol_ins[0])
        # Full-run (experimental): caller supplies cut-independent per-query costs.
        if precomputed_ins_costs is not None:
            return float(sum(precomputed_ins_costs))
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
    #
    # Callers may supply precomputed_ins_costs (already aligned to this call's query
    # positions) to make the ins cost cut-INDEPENDENT — i.e. measured on the full
    # rescue run rather than the per-k truncated substring. This is the experimental
    # full-run path (see _USE_FULL_RUN_INS); when None, fall back to the legacy
    # per-substring HP-run lookup.
    #
    # Reference-column (experimental) path: precomputed_refcol_ins is a length-(R+1)
    # vector indexed by DP GAP/column (genome HP context at that column), NOT by
    # query position. It changes the insertion transition to charge the cost of the
    # COLUMN where the insertion fires. This is the calibration-correct axis; it must
    # take the pure-Python branch because the numba kernel expects a length-Q,
    # query-indexed ins vector (feeding a length-(R+1) array would be mis-indexed).
    if precomputed_refcol_ins is not None:
        ins_col = precomputed_refcol_ins  # length R+1, indexed by gap/column j
        INF = float('inf')
        prev = [INF] * (R + 1)
        prev[0] = 0.0
        for j in range(1, R + 1):
            prev[j] = prev[j - 1] + del_costs[j - 1]  # leading ref deletions
        for i in range(1, Q + 1):
            curr = [INF] * (R + 1)
            curr[0] = i * ins_col[0]  # i insertions at gap 0 (before ref[0])
            qi = query[i - 1].upper()
            for j in range(1, R + 1):
                cost_sub = 0.0 if qi == ref[j - 1].upper() else sub
                diag  = prev[j - 1] + cost_sub          # match / mismatch
                above = prev[j]     + ins_col[j]        # insertion at column j
                left  = curr[j - 1] + del_costs[j - 1]  # deletion (ref base, no query)
                curr[j] = min(diag, above, left)
            prev = curr
        return min(prev)  # free right suffix: best end column over all j in [0, R]

    ins_costs: List[float]
    if precomputed_ins_costs is not None:
        ins_costs = precomputed_ins_costs
    elif penalty_table is not None:
        ins_costs = [
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
# Vectorized table-free DP: ALL query-suffix scores in one pass (perf fast path)
# ---------------------------------------------------------------------------

def _all_suffix_scores(
    query: str,
    ref: str,
    del_costs: List[float],
    ins: float = _FLAT_INS,
    sub: float = _FLAT_SUB,
) -> List[float]:
    """Return ``[score(query[k:], ref) for k in range(len(query)+1)]`` for the
    flat-insertion (penalty_table=None) DP of :func:`_score_hp_anchored`, computed
    in a SINGLE O(len(query)*len(ref)) pass instead of ``len(query)`` separate DPs.

    ``score(query[k:], ref)`` = query-global (every base of ``query[k:]`` consumed)
    + ref left end deletable at ``del_costs`` + ref right suffix free (``min`` over
    the end column) — identical to ``_score_hp_anchored(query[k:], ref, ...)`` when
    ``penalty_table`` is None so insertion cost is the flat ``ins``.

    Reversal trick: ``score(query[k:], ref)`` over the reversed sequences becomes a
    query-PREFIX-length-(L-k) alignment with the ref right end anchored (col R) and
    the ref left prefix free, so ONE forward DP over ``reverse(query)`` vs
    ``reverse(ref)`` exposes every ``k`` at row ``L-k``, column ``R``. Verified
    byte-identical to the per-``k`` reference (dev/concat_dp_prototype.py, 0/3000).
    """
    L, R = len(query), len(ref)
    out: List[float] = [None] * (L + 1)  # type: ignore[list-item]
    if L == 0:
        return [0.0]
    if R == 0:
        # all-insertion: score(query[k:], "") = (L-k)*ins
        return [(L - k) * ins for k in range(L + 1)]
    rq = query[::-1]
    rr = ref[::-1]
    dc = del_costs[::-1]  # del_costs[j] indexes ref[j]; reversed to align with rr
    INF = float("inf")
    # Row m=0: reverse(query)[:0] empty -> free ref left-prefix skip = 0 at every col.
    prev = [0.0] * (R + 1)
    out[L] = prev[R]  # k = L: empty query suffix -> score 0
    for m in range(1, L + 1):
        curr = [INF] * (R + 1)
        curr[0] = m * ins  # m query bases inserted, no ref consumed
        qm = rq[m - 1].upper()
        for j in range(1, R + 1):
            cost_sub = 0.0 if qm == rr[j - 1].upper() else sub
            diag = prev[j - 1] + cost_sub
            above = prev[j] + ins
            left = curr[j - 1] + dc[j - 1]
            curr[j] = diag if diag < above else above
            if left < curr[j]:
                curr[j] = left
        prev = curr
        out[L - m] = prev[R]  # score(query[(L-m):], ref)
    return out


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

    # PERF FAST PATH (table-free only): with penalty_table None the insertion cost is
    # flat, so ALL t1(k) and ALL t2(k) collapse into two vectorized DP passes that
    # reproduce min_k[t1(k)+t2(k)] EXACTLY (~14x fewer cell-ops). Strictly gated on
    # penalty_table is None — the only config where flat ins makes the fast path
    # provably byte-identical. See _USE_CONCAT_DP / dev/PERF_PIVOT_TABLEFREE_CONCAT.md.
    if _USE_CONCAT_DP and penalty_table is None:
        t1_vec = _all_suffix_scores(rescue, ref_exon2_start, del_costs_fwd)
        t2_vec = _all_suffix_scores(rescue[::-1], ref_intron_end_rev, del_costs_rev)
        best = float("inf")
        for k in range(L):
            # t2(k) = score(reverse(rescue[:k]), ref_intron_end_rev)
            #       = score(reverse(rescue)[L-k:], ...) = t2_vec[L-k]; t2(0)=0.
            s = t1_vec[k] + (0.0 if k == 0 else t2_vec[L - k])
            if s < best:
                best = s
                if best == 0.0:
                    break
        if profile is not None:
            profile.add_time('score_junction_total', time.perf_counter() - _t_total)
        return best, 0

    # Full-run (cut-independent) insertion costs — experimental. Precompute ONCE on
    # the FULL rescue window using each base's homopolymer run measured over the
    # whole rescue substring, then slice per-k (analogous to how the query itself is
    # sliced) so a homopolymer straddling the cut is charged at its FULL run length
    # regardless of where the k-cut falls. Only active with a penalty_table (the
    # flat-ins path is already cut-independent). See _USE_FULL_RUN_INS.
    full_ins_costs: Optional[List[float]] = None
    # Reference-column (cut-independent) insertion costs — experimental, PRECEDENCE
    # over full-run. Indexed by DP column (genome HP context at each aligned gap),
    # computed ONCE on the two fixed reference windows and reused UNCHANGED across
    # all k (never sliced) — the cut-independence proof. Mirrors del_costs_fwd/rev.
    refcol_ins_fwd: Optional[List[float]] = None   # length len(ref_exon2_start)+1
    refcol_ins_rev: Optional[List[float]] = None   # length len(ref_intron_end_rev)+1
    if _USE_REFCOL_INS and penalty_table is not None:
        refcol_ins_fwd = _precompute_refcol_ins_costs(
            ref_exon2_start, genome_seq, intron_end, False, penalty_table,
        )
        refcol_ins_rev = _precompute_refcol_ins_costs(
            ref_intron_end_rev, genome_seq, intron_end - 1, True, penalty_table,
        )
    elif _USE_FULL_RUN_INS and penalty_table is not None:
        full_ins_costs = [
            penalty_table.ins_cost(_hp_run_length(rescue, j), rescue[j]) for j in range(L)
        ]

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
            precomputed_ins_costs=(full_ins_costs[k:] if full_ins_costs is not None else None),
            precomputed_refcol_ins=refcol_ins_fwd,  # cut-independent: passed unchanged
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
                precomputed_ins_costs=(full_ins_costs[:k][::-1] if full_ins_costs is not None else None),
                precomputed_refcol_ins=refcol_ins_rev,  # cut-independent: passed unchanged
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
