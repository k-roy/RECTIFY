"""Station C v0 — the pool-level junction discovery gate (REPORT-FIRST).

The third station of the Re-aligner architecture
(dev/REALIGNER_LANDSCAPE_AND_PATH_20260809.md §3): Stations A (overhang
resolver) and B (consensus triage) act per READ; Station C judges per
JUNCTION, over the whole pool of supporting evidence.

v0 is deliberately **report-only**: it censuses junctions from a consensus /
corrected BAM, scores each one, and writes a per-junction admission table.
It annotates — it never deletes a read or a junction. The verdicts and their
thresholds come from the measured 644-series campaign:

- **Overhang quality** ``q`` (dev/PHASE2_OVERHANG_644H_20260811.md): per
  supporting read, the SHORT-exon-side overhang (the side of the N-op with
  the smaller contiguous aligned anchor) is scored
  ``max(0, I_eff(<=60 junction-adjacent bases) - 2.0 * errors_in_window)``,
  with ``I_eff`` from :mod:`overhang_informativeness` (the resolver's own
  currency). Junction ``q_max`` = max over supporting reads. Measured:
  canonical track support>=2 + q>=40 ~= 79% precision at 11/12 gold.
- **Repeat-context flags** (dev/STATIONC_REPEAT_FLAG_644I_20260811.md):
  annotation flag (rRNA/Ty/LTR/tRNA/telomere-class GFF features, zero
  measured gold cost) + genome self-homology track (bundled BED for yeast;
  ``minimap2 -DP -k19 -w19 -m200 G G`` merged +/-200 bp). Variant
  multiplicity was measured and REJECTED as a flag. A flag DEMOTES a
  candidate — to ``review`` on the canonical track, to
  needs-orthogonal-evidence on the non-canonical track; it never discards.
  **Flags are consulted on BOTH tracks** (planning/684c: the flag was
  originally unreachable when ``canonical_in_class=1``, which admitted a
  111 kb LTR-to-LTR artifact with both flags lit; ``canonical_in_class``
  is a weak shield exactly in repeat context, since one chance GT-AG
  anywhere in the ambiguity window satisfies it).
- **Length pre-gate** (planning/684c): junctions longer than
  ``cfg.max_intron`` (annotation-derived; 5,000 bp for yeast) demote
  BEFORE the verdict branches — never eligible for ``admit_candidate``.
  Deterministic aligner artifacts defeat every recurrence check, so
  length is the one term that stops the physically impossible.
- **Two-track admission** (dev/STATIONC_MAPPACBIO_HARVEST_20260810.md): the
  canonical-in-class and non-canonical tracks never share a threshold.

Verdicts (non-annotated junctions), applied in this order:

==============================  =============================================
``demote_orthogonal_evidence``  length > max_intron (the pre-gate, any
                                track); OR overlapping a known background-SV
                                region of the reference (bundled
                                ``R64_background_sv.bed`` — e.g. the chrIII
                                SRD1 flank-A Ty1 replacement, yKR888 T2T /
                                planning/730 W6 — any track); OR
                                non-canonical AND (repeat-flagged
                                OR q < q_noncanon): admissible only with
                                orthogonal evidence (short-read/COMPASS,
                                cross-sample recurrence, mm2-side distress)
``review``                      canonical AND repeat-flagged (flag caps the
                                canonical track at review); or evidence
                                below the admit bar but nothing against it
                                (low support and/or low q)
``admit_candidate``             canonical track: unflagged AND support>=2
                                AND q>=q_canon; non-canonical track:
                                unflagged AND support>=2 AND q>=q_noncanon
==============================  =============================================

Cross-sample recurrence and the orthogonal-evidence channels are recorded as
columns for downstream tooling but not computed in v0 (single-sample scope).

Two coverage properties the gate depends on (ISSUE-008, Sumner RNA004 panel,
dev/sumner_misplaced_panel_20260904):

- **The census must see the junctions RECTIFY creates.** Module 2H realizes a
  single-boundary move by planting a compensating I/D op immediately beside the
  N (``junction_refiner.py`` L1300-1309, deliberately exempt from the
  both-boundary refusal). Reading only the op next to the N scored those as
  anchor 0 and dropped them from the census entirely — neither admitted,
  reviewed nor demoted. ``_anchor_run`` sums the aligned run across a bounded
  number of intervening indel ops (``cfg.adj_indel_max_ops`` /
  ``cfg.adj_indel_max_bp``, both exposed on the CLI) and records what it
  stepped over (``adj_indel_l``/``adj_indel_r``/``n_adj_indel``). Measured:
  16/121 → 80/121 created N-ops censused at the defaults; ``_anchor_run``'s
  docstring carries the full budget sensitivity, including which of the two
  budgets binds. Thresholds are unchanged — this is coverage, not policy.
- **A missing track must say so.** The repeat, self-homology and background-SV
  tracks have no human equivalent, so all three read empty on a human run while
  the table looked like a clean bill of health. An absent track now writes
  ``TRACK_UNAVAILABLE`` in its column and lands in
  ``summary['tracks_unavailable']`` with a logged warning. The annotated-intron
  shield is likewise parsed from GFF3 *and* GTF dialects (:func:`_exon_parents`).
"""

from __future__ import annotations

import gzip
import json
import logging
import os
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from ..splice.overhang_informativeness import (
    ambiguity_window,
    canonical_in_class,
    effective_information_bits,
)

logger = logging.getLogger(__name__)

# GFF feature types whose neighbourhoods form the annotation repeat flag
# (644i: zero measured gold cost at every threshold).
REPEAT_FEATURE_TYPES = frozenset({
    'rRNA_gene', 'LTR_retrotransposon', 'long_terminal_repeat',
    'transposable_element_gene', 'telomere', 'telomeric_repeat',
    'X_element', 'X_element_combinatorial_repeat', 'Y_prime_element',
    'tRNA_gene',
})


#: Value written into a flag column when the TRACK IT COMES FROM does not
#: exist for this genome, as opposed to "the track was consulted and this
#: junction is clean". Human input has no bundled self-homology / background-SV
#: track and no REPEAT_FEATURE_TYPES equivalent, so all three columns used to
#: read 0/'' — indistinguishable from a clean bill of health (the Sumner
#: RNA004 panel, dev/sumner_misplaced_panel_20260904 CP2: three of Station C's
#: four demotion terms were silently dead). A sentinel here is deliberately
#: TRUTHY, so a downstream reader that tests the cell errs toward caution.
#: Station C's own verdict never sees it — the verdict reads the boolean flags.
TRACK_UNAVAILABLE = 'track_unavailable'


@dataclass
class PoolGateConfig:
    min_anchor: int = 8          # census: min adjacent exon anchor per N-op
    # Census anchor walk (see `_anchor_run`): how much intervening I/D the walk
    # may step over between the N-op and the aligned run that anchors it.
    # Module 2H realizes a single-boundary junction move by planting a
    # COMPENSATING indel right beside the N, so reading only ops[i±1] scored
    # 87% of the junctions RECTIFY creates as anchor 0 — never censused, never
    # gated. Both budgets are per SIDE of the N-op. The full measured
    # sensitivity of the pair is tabulated in `_anchor_run`; on the Sumner
    # panel the OP limit is the binding one, not the bp budget.
    adj_indel_max_ops: int = 2   # intervening I/D ops tolerated, per side
    adj_indel_max_bp: int = 30   # their summed length, per side
    overhang_cap: int = 60       # junction-adjacent query bases assessed
    err_bits: float = 2.0        # bits removed per alignment error in window
    q_canon: float = 40.0        # canonical-track admit threshold (644h)
    q_noncanon: float = 80.0     # non-canonical-track admit threshold (644h/i)
    min_support: int = 2         # recurrence gate within this sample
    repeat_margin: int = 500     # bp around an annotated repeat feature
    max_ambiguity_shift: int = 30
    # Length pre-gate (planning/684c): junctions longer than this are never
    # ELIGIBLE for admit_candidate — they demote before the verdict branches
    # run. Callers should pass the annotation-derived value
    # (multi_aligner.derive_max_intron; 5000 is exactly what the bundled
    # yeast annotation derives). A deterministic aligner artifact defeats
    # every recurrence check, so length is the one term that stops the
    # physically impossible (the measured case: a 111 kb LTR-to-LTR junction
    # with support=19, q=96.5).
    max_intron: int = 5000


# ---------------------------------------------------------------------------
# Repeat-context flags
# ---------------------------------------------------------------------------

class IntervalFlags:
    """Merged, sorted per-chrom intervals with an overlap query."""

    def __init__(self) -> None:
        self._iv: Dict[str, List[Tuple[int, int, str]]] = {}
        self.n_intervals = 0

    def add(self, chrom: str, start: int, end: int, label: str) -> None:
        self._iv.setdefault(chrom, []).append((max(0, start), end, label))

    def freeze(self) -> 'IntervalFlags':
        for chrom in self._iv:
            self._iv[chrom].sort()
        self.n_intervals = sum(len(v) for v in self._iv.values())
        return self

    def flag_of(self, chrom: str, start: int, end: int) -> Optional[str]:
        for a, b, label in self._iv.get(chrom, ()):
            if a >= end:
                break
            if start < b and end > a:
                return label
        return None


def _gff_open(path: Path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def load_repeat_flags(annotation_path: Path, margin: int = 500) -> IntervalFlags:
    """Annotation leg of the repeat flag: REPEAT_FEATURE_TYPES +/- margin."""
    flags = IntervalFlags()
    with _gff_open(annotation_path) as fh:
        for line in fh:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 5 or f[2] not in REPEAT_FEATURE_TYPES:
                continue
            flags.add(f[0], int(f[3]) - 1 - margin, int(f[4]) + margin, f[2])
    return flags.freeze()


def load_selfhom_bed(bed_path: Path) -> IntervalFlags:
    flags = IntervalFlags()
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.split('\t')
            flags.add(f[0], int(f[1]), int(f[2]), 'selfhom')
    return flags.freeze()


def find_bundled_selfhom_bed(genome_path: Path) -> Optional[Path]:
    """A ``*selfhomology.bed`` beside the genome, or the bundled yeast track."""
    candidates = list(Path(genome_path).parent.glob('*selfhomology.bed'))
    if candidates:
        return candidates[0]
    try:
        import rectify
        bundled = (Path(rectify.__file__).parent / 'data' / 'genomes'
                   / 'saccharomyces_cerevisiae' / 'R64_selfhomology.bed')
        if bundled.exists() and Path(genome_path).stem.startswith('S288C'):
            return bundled
    except Exception:  # pragma: no cover - import-shape guard
        pass
    return None


def load_background_sv_bed(bed_path: Path) -> IntervalFlags:
    """Known background-SV regions of the REFERENCE (BED col 4 = SV name).

    Unlike the repeat/self-homology tracks (priors), these are intervals where
    the reference assembly is KNOWN not to represent common strain genomes —
    e.g. the R64 chrIII SRD1 flank-A segment, deleted-and-replaced-by-Ty1 in
    real strains (yKR888 T2T; Kevin 2026-08-20). Reads bridging such a region
    produce canonical-looking junctions that are reference artifacts, so an
    overlapping junction demotes to ``demote_orthogonal_evidence`` on BOTH
    tracks (planning/730 W6).
    """
    flags = IntervalFlags()
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            label = f[3] if len(f) > 3 and f[3] else 'background_sv'
            flags.add(f[0], int(f[1]), int(f[2]), label)
    return flags.freeze()


def find_bundled_background_sv_bed(genome_path: Path) -> Optional[Path]:
    """A ``*background_sv.bed`` beside the genome, or the bundled R64 track."""
    candidates = list(Path(genome_path).parent.glob('*background_sv.bed'))
    if candidates:
        return candidates[0]
    try:
        import rectify
        bundled = (Path(rectify.__file__).parent / 'data' / 'genomes'
                   / 'saccharomyces_cerevisiae' / 'R64_background_sv.bed')
        if bundled.exists() and Path(genome_path).stem.startswith('S288C'):
            return bundled
    except Exception:  # pragma: no cover - import-shape guard
        pass
    return None


# ---------------------------------------------------------------------------
# Per-read short-exon-side overhang quality (the 644h scorer)
# ---------------------------------------------------------------------------

def _offsets(read: pysam.AlignedSegment) -> Tuple[List[int], List[int]]:
    ops = read.cigartuples
    qoff = [0]
    roff = [read.reference_start]
    for op, ln in ops:
        qoff.append(qoff[-1] + (ln if op in (0, 1, 4, 7, 8) else 0))
        roff.append(roff[-1] + (ln if op in (0, 2, 3, 7, 8) else 0))
    return qoff, roff


def _side_features(read, qoff, roff, i_n, chrom_seq, direction, cap):
    """One side of the N-op at cigar index ``i_n`` (see 644h).

    Handles SAM ``=`` SEQ compression (mapPacBio-style) by decoding matches
    from the reference — without this every base reads as a mismatch and all
    scores are silently 0 (the 644h smoke lesson).
    """
    ops = read.cigartuples
    qseq = read.query_sequence
    collected: List[str] = []
    errors = 0
    anchor_q = 0

    def note_err(n: int = 1) -> None:
        nonlocal errors
        errors += n

    rng = range(i_n + 1, len(ops)) if direction > 0 else range(i_n - 1, -1, -1)
    for k in rng:
        op, ln = ops[k]
        if op in (3, 4, 5):
            break
        if op in (0, 7, 8):
            for t in range(ln):
                if direction > 0:
                    qi, ri = qoff[k] + t, roff[k] + t
                else:
                    qi, ri = qoff[k + 1] - 1 - t, roff[k + 1] - 1 - t
                if len(collected) < cap:
                    b = qseq[qi].upper()
                    if b == '=':
                        collected.append(chrom_seq[ri] if ri < len(chrom_seq)
                                         else 'N')
                    else:
                        collected.append(b)
                        if ri >= len(chrom_seq) or b != chrom_seq[ri]:
                            note_err()
                else:
                    break
            anchor_q += ln
        elif op == 1:
            for t in range(ln):
                if len(collected) < cap:
                    qi = (qoff[k] + t) if direction > 0 else (qoff[k + 1] - 1 - t)
                    b = qseq[qi].upper()
                    collected.append('N' if b == '=' else b)
                    note_err()
                else:
                    break
            anchor_q += ln
        elif op == 2:
            if len(collected) < cap:
                note_err(ln)
    return {'anchor_q': anchor_q, 'seq': ''.join(collected), 'errors': errors}


def read_junction_q(read, qoff, roff, i_n, chrom_seq, cfg: PoolGateConfig) -> float:
    """Short-exon-side overhang quality of one junction on one read."""
    left = _side_features(read, qoff, roff, i_n, chrom_seq, -1, cfg.overhang_cap)
    right = _side_features(read, qoff, roff, i_n, chrom_seq, +1, cfg.overhang_cap)
    short = left if left['anchor_q'] <= right['anchor_q'] else right
    i_eff = effective_information_bits(short['seq'])
    return max(0.0, i_eff - cfg.err_bits * short['errors'])


# ---------------------------------------------------------------------------
# Census + verdicts
# ---------------------------------------------------------------------------

def _canonicalize(chrom_seq: str, start: int, end: int, max_shift: int) -> Tuple[int, int]:
    l_amb, _ = ambiguity_window(chrom_seq, start, end, max_shift=max_shift)
    return start - l_amb, end - l_amb


#: What ended an anchor walk, by the CIGAR op it stopped on (ISSUE-016 ledger).
_WALK_STOP = {3: 'n_op', 4: 'softclip', 5: 'hardclip', 6: 'pad'}


def _anchor_walk(ops, i_n: int, direction: int,
                 adj_indel_max_ops: int = 2,
                 adj_indel_max_bp: int = 30) -> Tuple[int, str, str]:
    """Aligned bases anchoring one side of the N-op at cigar index ``i_n``.

    Returns ``(anchor_bases, adjacent_indel_label, stop)``; :func:`_anchor_run`
    is the two-tuple view every earlier caller uses. ``stop`` names what ended
    the walk — ``n_op`` (the next intron: a short exon between two N-ops),
    ``softclip`` / ``hardclip`` / ``pad``, ``read_end``, or the budget it
    exceeded (``indel_ops`` / ``indel_bp``) — so a REFUSED N-op can say why its
    anchor was short instead of vanishing (ISSUE-016). The label uses the
    ``747_panel.py::junctions()`` encoding — ``'I1'``, ``'D10'``, or ``''`` for
    the FIRST indel op met walking toward the flank.

    The run sums the contiguous ``M``/``=``/``X`` ops toward the flank, stepping
    over at most ``adj_indel_max_ops`` intervening ``I``/``D`` ops whose lengths
    sum to at most ``adj_indel_max_bp``; it stops at ``N``/``S``/``H``/``P``, and
    at an indel run past either budget. Both budgets apply PER SIDE, and the bp
    budget is a TOTAL over the stepped-over ops (not per op).

    **Why walk at all.** Module 2H realizes a single-boundary junction move by
    re-labelling the boundary with a COMPENSATING indel placed immediately
    beside the N-op — deliberate behavior, explicitly exempted from the
    both-boundary refusal (``junction_refiner.py`` L1300-1309):
    ``69M 3699N 61M 468N`` becomes ``69M 3699N 60M 1I 469N``. Reading only
    ``ops[i±1]`` and requiring it to be M/=/X therefore scored exactly the
    junctions most in need of gating as anchor 0, so they were never enumerated
    — not admitted, not reviewed, not demoted.

    **Measured sensitivity** — created N-ops censused on the Sumner RNA004
    panel (dev/sumner_misplaced_panel_20260904, 100 reads, 121 junctions that
    ``correct`` created), replayed from the panel's own CIGARs. The single
    adjacent op (the pre-fix rule) censuses **16/121**:

    ==========  ====  ====  ====  ====  =====  =====
    max_ops \\   10    20    30    50    100    inf   (max_bp)
    ==========  ====  ====  ====  ====  =====  =====
    1             60    69    69    75     75     75
    **2**         69  **80**  **80**  87     87     87
    3             70    84    85    92     92     92
    4             70    85    90    97     97     97
    inf           70    85    90    97     97     97
    ==========  ====  ====  ====  ====  =====  =====

    Read the table before changing either default. On this panel the **op count
    is the binding constraint, not the bp budget**: at ``max_ops=2`` the 20 and
    30 bp budgets are identical (80/121), and the 97/121 ceiling needs
    ``max_ops>=4`` *and* ``max_bp>=50``. The extra junctions bought by relaxing
    the op limit are 2F's ``4M4I3M4I1M``-style rescued exons — shapes where no
    contiguous anchor really exists and the "anchor" is a sum over fragments —
    which is why the default keeps ``max_ops=2`` while allowing 30 bp of
    single-boundary indel. Both are exposed on ``rectify pool-gate``
    (``--adj-indel-max-ops`` / ``--adj-indel-max-bp``) so the trade-off can be
    re-measured rather than argued.

    This is a **coverage** fix, not a threshold change: ``min_anchor`` and both
    track thresholds are untouched, and a genuinely short anchor
    (``4M 60N 30M``) still yields 4 and is still refused. The returned anchor is
    monotone non-decreasing in both budgets and is never smaller than the
    pre-fix single-op rule, so widening a budget can only ADD junctions to the
    census, never drop one.

    The sibling ``_side_features`` (the q-scorer) has always walked I/D ops, so
    this also removes an internal inconsistency: the score was indel-tolerant
    while the census gate that fed it was not.
    """
    rng = (range(i_n + 1, len(ops)) if direction > 0
           else range(i_n - 1, -1, -1))
    total = 0
    n_ops = 0
    n_bp = 0
    label = ''
    stop = 'read_end'
    for k in rng:
        op, ln = ops[k]
        if op in (0, 7, 8):          # M / = / X
            total += ln
        elif op in (1, 2):           # I / D
            if not label:
                label = ('I' if op == 1 else 'D') + str(ln)
            n_ops += 1
            n_bp += ln
            if n_ops > adj_indel_max_ops:
                stop = 'indel_ops'
                break
            if n_bp > adj_indel_max_bp:
                stop = 'indel_bp'
                break
        else:                        # N / S / H / P
            stop = _WALK_STOP.get(op, 'other')
            break
    return total, label, stop


def _anchor_run(ops, i_n: int, direction: int,
                adj_indel_max_ops: int = 2,
                adj_indel_max_bp: int = 30) -> Tuple[int, str]:
    """``(anchor_bases, adjacent_indel_label)`` — see :func:`_anchor_walk`."""
    total, label, _stop = _anchor_walk(ops, i_n, direction,
                                       adj_indel_max_ops, adj_indel_max_bp)
    return total, label


def _refusal_reason(lf: int, rt: int, stop_l: str, stop_r: str,
                    min_anchor: int) -> str:
    """``L=<stop>`` / ``R=<stop>`` / ``L=<stop>+R=<stop>`` for the short side(s)."""
    parts = []
    if lf < min_anchor:
        parts.append(f'L={stop_l}')
    if rt < min_anchor:
        parts.append(f'R={stop_r}')
    return '+'.join(parts)


@dataclass
class CensusLedger:
    """Everything ``census_bam`` saw that the junction table cannot show.

    ISSUE-016: the tester matched the junctions RECTIFY *created* on a
    145k-read slice against ``<prefix>.pool_gate.tsv`` and found 87 % ABSENT at
    every op budget (0/2/4) and support gate (1/2). Neither knob was binding
    because the table cannot show two things: it lists NON-annotated junctions
    only (2F/2H create mostly annotated ones — that is their job), and it keys
    every junction at the LEFTMOST ambiguity-equivalent coordinate
    (``_canonicalize``) while the BAM's N-op — and any list derived from it —
    sits at the motif coordinate. On the 100-read panel: 10 created junctions,
    0 in the table, 271 censused of which 231 annotated.

    The ledger records every N-op outcome per RAW coordinate, so any junction
    can be attributed as censused (annotated or reported), refused (with the
    anchor-walk stop reason per side), or never seen.
    """
    n_reads: int = 0
    reads_skipped: Dict[str, int] = field(default_factory=dict)
    n_ops_seen: int = 0
    n_ops_censused: int = 0
    n_ops_refused: int = 0
    #: RAW (chrom, start, end) -> refusal reason -> N-op occurrences.
    refusals: Dict[Tuple[str, int, int], Dict[str, int]] = field(default_factory=dict)
    #: RAW key -> the best min(anchor) any occurrence reached (how close it came).
    best_anchor: Dict[Tuple[str, int, int], int] = field(default_factory=dict)
    #: canonical key -> RAW (start, end) -> occurrences censused under it.
    raw_keys: Dict[Tuple[str, int, int], Dict[Tuple[int, int], int]] = field(default_factory=dict)

    def note_skip(self, why: str) -> None:
        self.reads_skipped[why] = self.reads_skipped.get(why, 0) + 1

    def note_census(self, key: Tuple[str, int, int], raw: Tuple[int, int]) -> None:
        self.n_ops_censused += 1
        d = self.raw_keys.setdefault(key, {})
        d[raw] = d.get(raw, 0) + 1

    def note_refusal(self, raw_key: Tuple[str, int, int], reason: str,
                     anchor: int) -> None:
        self.n_ops_refused += 1
        d = self.refusals.setdefault(raw_key, {})
        d[reason] = d.get(reason, 0) + 1
        if anchor > self.best_anchor.get(raw_key, -1):
            self.best_anchor[raw_key] = anchor

    def reason_counts(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for reasons in self.refusals.values():
            for reason, n in reasons.items():
                out[reason] = out.get(reason, 0) + n
        return dict(sorted(out.items(), key=lambda kv: (-kv[1], kv[0])))

    def majority_raw(self, key: Tuple[str, int, int]) -> Tuple[int, int, int]:
        """``(start_raw, end_raw, n_raw_variants)`` for a censused key — the
        N-op coordinates most supporting reads actually carry; falls back to
        the key itself when nothing was recorded."""
        variants = self.raw_keys.get(key)
        if not variants:
            return key[1], key[2], 0
        (s, e), _n = sorted(variants.items(), key=lambda kv: (-kv[1], kv[0]))[0]
        return s, e, len(variants)


def _pick_label(counts: Dict[str, int]) -> str:
    """Most frequent adjacent-indel label; ties broken by label, not BAM order."""
    if not counts:
        return ''
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def census_bam(bam_path: str, genome: Dict[str, str], cfg: PoolGateConfig,
               max_q_reads_per_junction: int = 50,
               ledger: Optional[CensusLedger] = None) -> Dict[Tuple[str, int, int], dict]:
    """One streaming pass: per canonical junction, support + q statistics.

    ``max_q_reads_per_junction`` caps the per-read overhang scoring (q_max is
    a max — extra reads beyond the cap still count toward support).

    The per-side anchor is the aligned run reached by ``_anchor_run``, which
    steps over the compensating indel Module 2H plants beside every junction it
    moves. Each record also carries what was stepped over: ``adj_l`` / ``adj_r``
    (label -> count over supporting reads) and ``n_adj_indel`` (supporting reads
    with an adjacent indel on either side), so the flag survives into the
    report instead of being silently absorbed.

    ``ledger`` (optional): a :class:`CensusLedger` that receives every N-op
    outcome — censused under which canonical key from which RAW coordinates,
    or refused with the anchor-walk reason per side — and every skipped read.
    """
    J: Dict[Tuple[str, int, int], dict] = {}
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as fh:
        for read in fh.fetch(until_eof=True):
            if ledger is not None:
                ledger.n_reads += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                if ledger is not None:
                    ledger.note_skip('nonprimary')
                continue
            if not read.cigartuples:
                if ledger is not None:
                    ledger.note_skip('no_cigar')
                continue
            chrom = read.reference_name
            chrom_seq = genome.get(chrom)
            if chrom_seq is None:
                if ledger is not None:
                    ledger.note_skip('chrom_missing')
                continue
            offs = None
            ops = read.cigartuples
            ref = read.reference_start
            for i, (op, ln) in enumerate(ops):
                if op == 3:
                    lf, adj_l, stop_l = _anchor_walk(ops, i, -1,
                                                     cfg.adj_indel_max_ops,
                                                     cfg.adj_indel_max_bp)
                    rt, adj_r, stop_r = _anchor_walk(ops, i, +1,
                                                     cfg.adj_indel_max_ops,
                                                     cfg.adj_indel_max_bp)
                    if ledger is not None:
                        ledger.n_ops_seen += 1
                    if min(lf, rt) < cfg.min_anchor:
                        if ledger is not None:
                            ledger.note_refusal(
                                (chrom, ref, ref + ln),
                                _refusal_reason(lf, rt, stop_l, stop_r,
                                                cfg.min_anchor),
                                min(lf, rt))
                    else:
                        s, e = _canonicalize(chrom_seq, ref, ref + ln,
                                             cfg.max_ambiguity_shift)
                        key = (chrom, s, e)
                        if ledger is not None:
                            ledger.note_census(key, (ref, ref + ln))
                        rec = J.get(key)
                        if rec is None:
                            rec = J[key] = {'support': 0, 'q_max': 0.0,
                                            'q_2nd': 0.0, 'q_scored': 0,
                                            'adj_l': {}, 'adj_r': {},
                                            'n_adj_indel': 0}
                        rec['support'] += 1
                        if adj_l:
                            rec['adj_l'][adj_l] = rec['adj_l'].get(adj_l, 0) + 1
                        if adj_r:
                            rec['adj_r'][adj_r] = rec['adj_r'].get(adj_r, 0) + 1
                        if adj_l or adj_r:
                            rec['n_adj_indel'] += 1
                        if (rec['q_scored'] < max_q_reads_per_junction
                                and read.query_sequence is not None):
                            if offs is None:
                                offs = _offsets(read)
                            q = read_junction_q(read, offs[0], offs[1], i,
                                                chrom_seq, cfg)
                            rec['q_scored'] += 1
                            if q > rec['q_max']:
                                rec['q_2nd'] = rec['q_max']
                                rec['q_max'] = q
                            elif q > rec['q_2nd']:
                                rec['q_2nd'] = q
                if op in (0, 2, 3, 7, 8):
                    ref += ln
    return J


#: GFF3 exon -> transcript link. Matched FIRST, so the yeast path is unchanged.
_PARENT_RE = re.compile(r'Parent=([^;]+)')
#: GTF exon -> transcript link (GENCODE/Ensembl `transcript_id "X";`, and the
#: unquoted `transcript_id=X;` some GFF3 writers add). The leading
#: `(?:^|[;\s])` stops it matching a longer key that merely ends in
#: `transcript_id`.
_TRANSCRIPT_ID_RE = re.compile(r'(?:^|[;\s])transcript_id[ =]+"?([^";]+)"?')


def _exon_parents(attrs: str) -> List[str]:
    """Transcript key(s) an exon line belongs to, across annotation dialects.

    GFF3 (the bundled yeast R64 annotation) writes ``Parent=<id>[,<id>…]``;
    GTF (GENCODE, Ensembl, StringTie) writes ``transcript_id "<id>";`` and no
    ``Parent=`` at all. ``Parent=`` is matched FIRST, so every GFF3 line
    resolves exactly as it did before this function existed — the fallback is
    only reached by a line with no ``Parent=``.

    Station C used to look for ``Parent=`` only, so a GENCODE GTF (which also
    has no ``intron`` feature) yielded **zero** annotated introns and every real
    annotated junction was reported as a discovery candidate — measured on
    GENCODE v48 chr5, 222 junctions censused, ``n_annotated=0``
    (dev/sumner_misplaced_panel_20260904 CP2). The sibling
    ``multi_aligner.derive_max_intron`` has known the GTF dialect all along.
    """
    m = _PARENT_RE.search(attrs)
    if m:
        return m.group(1).split(',')
    m = _TRANSCRIPT_ID_RE.search(attrs)
    return [m.group(1).strip()] if m else []


def load_annotated_canonical(annotation_path: Path, genome: Dict[str, str],
                             cfg: PoolGateConfig) -> set:
    """Annotated introns (intron features + inferred from exon adjacency),
    ambiguity-canonicalised so census keys match.

    Both annotation dialects are handled: explicit ``intron`` features (yeast
    GFF3) and per-transcript exon adjacency keyed on ``Parent=`` (GFF3) or
    ``transcript_id`` (GTF — see :func:`_exon_parents`).
    """
    from collections import defaultdict
    exons = defaultdict(list)
    chrom_pool: Dict[str, str] = {}
    ann = set()
    with _gff_open(annotation_path) as fh:
        for line in fh:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            chrom, ftype, s1, e1, attrs = f[0], f[2], f[3], f[4], f[8]
            if ftype in ('intron', 'five_prime_UTR_intron'):
                ann.add((chrom, int(s1) - 1, int(e1)))
            elif ftype in ('exon', 'noncoding_exon'):
                # Pool the chrom string: a human GTF has ~1.6 M exon lines and
                # str.split makes a fresh object for each.
                chrom = chrom_pool.setdefault(chrom, chrom)
                for parent in _exon_parents(attrs):
                    exons[parent].append((chrom, int(s1) - 1, int(e1)))
    for parent, ex in exons.items():
        if len(ex) < 2:
            continue
        ex.sort(key=lambda x: x[1])
        for (c1, _s, e_prev), (c2, s_next, _e) in zip(ex, ex[1:]):
            if c1 == c2 and s_next > e_prev:
                ann.add((c1, e_prev, s_next))
    out = set()
    for chrom, s, e in ann:
        seq = genome.get(chrom)
        if seq is None:
            out.add((chrom, s, e))
        else:
            cs, ce = _canonicalize(seq, s, e, cfg.max_ambiguity_shift)
            out.add((chrom, cs, ce))
    return out


#: Same bounds ``multi_aligner.derive_max_intron`` clamps to, so the pre-gate
#: can never be looser than the aligner's own cap or absurdly tight.
_POOL_GATE_MAX_INTRON_BOUNDS = (1000, 500_000)


def _quantile(sorted_values: List[int], q: float) -> float:
    """Linear-interpolated quantile of a pre-sorted list (numpy-free)."""
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    pos = q * (len(sorted_values) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_values) - 1)
    return sorted_values[lo] + (pos - lo) * (sorted_values[hi] - sorted_values[lo])


def annotated_intron_lengths(annotation_path: Path) -> List[int]:
    """Every annotated intron length in an annotation, sorted ascending.

    Explicit ``intron`` features when the annotation has them (the R64 GFF3
    does); otherwise per-transcript exon gaps, across both dialects via
    :func:`_exon_parents`. Same precedence as
    ``multi_aligner.derive_max_intron``, which takes only the max of this.
    """
    from collections import defaultdict
    lens: List[int] = []
    exons = defaultdict(list)
    with _gff_open(annotation_path) as fh:
        for line in fh:
            if line.startswith('##FASTA') or line.startswith('>'):
                break
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            if f[2] == 'intron':
                lens.append(int(f[4]) - int(f[3]) + 1)
            elif f[2] in ('exon', 'noncoding_exon') and not lens:
                # Only pay the exon bookkeeping while no intron feature has
                # been seen — the derive_max_intron precedence.
                for parent in _exon_parents(f[8]):
                    exons[parent].append((int(f[3]), int(f[4])))
    if not lens:
        for ivs in exons.values():
            if len(ivs) < 2:
                continue
            ivs.sort()
            for (_s1, e1), (s2, _e2) in zip(ivs, ivs[1:]):
                gap = s2 - e1 - 1
                if gap > 0:
                    lens.append(gap)
    lens.sort()
    return lens


def derive_pool_gate_max_intron(
    annotation_path: Optional[Path],
    quantile: float = 0.995,
    multiplier: int = 2,
    fallback: int = 5000,
) -> int:
    """Station C's length pre-gate bound: ``multiplier x`` the ``quantile``
    of the annotated intron-length distribution, rounded up to 100 and clamped.

    **Why not the aligner's bound.** ``multi_aligner.derive_max_intron`` takes
    ``2 x the LONGEST annotated intron`` and clamps to 500,000. That is the
    right shape for an aligner's hard ``-G`` cap — it must not amputate real
    biology — but it makes Station C's pre-gate inert on any organism with a
    long tail. Measured on GENCODE v48 chr5: longest annotated intron 772,519
    bp, so the aligner rule saturates at the 500,000 clamp and NO junction can
    ever trip the pre-gate. The single term meant to stop the physically
    impossible was dead on every human run.

    A high quantile keeps the property that matters (real introns pass) while
    restoring a live ceiling:

    ===================  =========  =========  ============  ==============
    annotation                    n       max   p99.5 x 2    derive_max_intron
    ===================  =========  =========  ============  ==============
    R64 GFF3 (yeast)           378      2,483       **5,000**          5,000
    R64 GTF (yeast)            378      2,483       **5,000**          5,000
    GENCODE v48 chr5        38,849    772,519       310,100    500,000 (clamped)
    ===================  =========  =========  ============  ==============

    Yeast is unchanged — with 378 introns the top 0.5% IS the 2,483 bp maximum,
    so p99.5 x 2 reproduces the historical 5,000 exactly. This is asserted, not
    assumed (``test_station_c_max_intron.py``). Human gains a bound at 310 kb,
    below the clamp, so a 400 kb "intron" now demotes instead of admitting.

    ``p99.9 x 2`` was measured too and rejected: 524,000 on chr5, i.e. back
    above the clamp and inert again.

    The aligner's own ``max_intron`` is deliberately NOT changed — a pre-gate
    that annotates and a cap that amputates want different bounds.
    """
    if not annotation_path:
        return fallback
    try:
        lens = annotated_intron_lengths(Path(annotation_path))
    except Exception as exc:  # unreadable/malformed annotation: fall back
        logger.warning('station-c: could not derive the length pre-gate from '
                       '%s (%s); using %d', annotation_path, exc, fallback)
        return fallback
    if not lens:
        logger.warning('station-c: no annotated introns in %s; length pre-gate '
                       'falls back to %d', annotation_path, fallback)
        return fallback
    raw = multiplier * _quantile(lens, quantile)
    derived = -(-int(raw) // 100) * 100          # round UP to the nearest 100
    lo, hi = _POOL_GATE_MAX_INTRON_BOUNDS
    value = max(lo, min(hi, derived))
    # Log the PRE-clamp value too: when the clamp bites, the pre-gate is back
    # to being inert and the operator has to be able to see that it happened.
    logger.info('station-c length pre-gate: %d bp (%gx p%.4g of %d annotated '
                'introns; pre-clamp %d, max annotated %d)',
                value, multiplier, quantile * 100, len(lens), derived, lens[-1])
    if derived > hi:
        logger.warning('station-c: the derived length pre-gate (%d bp) exceeds '
                       'the %d bp ceiling and was clamped — the pre-gate is '
                       'inert for junctions between %d and %d bp.',
                       derived, hi, hi, derived)
    return value


def pool_gate(
    bam_path: str,
    genome: Dict[str, str],
    annotation_path: Path,
    cfg: Optional[PoolGateConfig] = None,
    selfhom_bed: Optional[Path] = None,
    background_sv_bed: Optional[Path] = None,
    attribute: Optional[List[Tuple[str, int, int]]] = None,
) -> Tuple[List[dict], dict]:
    """Run Station C v0 over one BAM. Returns (rows, summary).

    Rows cover NON-annotated junctions only (one dict per junction, sorted by
    verdict then descending support); annotated junctions are counted in the
    summary AND listed in ``summary['tables']['annotated']`` (ISSUE-016), the
    refused N-ops in ``summary['tables']['census_refusals']``, and — when
    ``attribute`` (RAW junction triples, see :func:`load_junction_list`) is
    given — one status per listed junction in ``summary['tables']['attribution']``.
    Report-only: nothing upstream is modified.

    A flag whose TRACK does not exist for this genome is reported as
    ``TRACK_UNAVAILABLE`` in its column, never as 0/'' — an absent track and a
    consulted-and-clean track must not read the same (Sumner RNA004: three of
    the four demotion terms were dead on human and the table looked clean).
    ``summary['tracks_unavailable']`` lists them and each is logged as a
    warning.
    """
    cfg = cfg or PoolGateConfig()
    ann_flags = load_repeat_flags(annotation_path, margin=cfg.repeat_margin)
    sh_flags = load_selfhom_bed(selfhom_bed) if selfhom_bed else None
    bg_flags = (load_background_sv_bed(background_sv_bed)
                if background_sv_bed else None)
    annotated = load_annotated_canonical(annotation_path, genome, cfg)

    # Which demotion terms actually have something behind them here. The
    # repeat track is annotation-derived, so "no REPEAT_FEATURE_TYPES feature
    # in this annotation" is the same practical condition as "no track": the
    # column would otherwise report a clean bill of health it never checked.
    have_repeat = ann_flags.n_intervals > 0
    have_selfhom = sh_flags is not None
    have_bg = bg_flags is not None
    tracks_unavailable = [name for name, ok in (
        ('repeat', have_repeat), ('selfhom', have_selfhom),
        ('background_sv', have_bg)) if not ok]
    for name in tracks_unavailable:
        logger.warning(
            'station-c: %s track UNAVAILABLE for this genome/annotation — the '
            '%s column reports %r, not a clean result. Junctions are gated on '
            'the remaining terms only.',
            name, f'{name}_flag', TRACK_UNAVAILABLE)
    if not annotated:
        logger.warning(
            'station-c: 0 annotated introns parsed from %s — every annotated '
            'junction will be reported as a discovery candidate.',
            annotation_path)

    ledger = CensusLedger()
    J = census_bam(bam_path, genome, cfg, ledger=ledger)

    rows: List[dict] = []
    annotated_rows: List[dict] = []
    n_annotated = 0
    verdict_counts: Dict[str, int] = {}
    for (chrom, s, e), rec in J.items():
        raw_s, raw_e, n_raw = ledger.majority_raw((chrom, s, e))
        if (chrom, s, e) in annotated:
            n_annotated += 1
            # ISSUE-016: annotated junctions were counted but never listed, so
            # a junction 2F/2H created ONTO the annotation read as ABSENT.
            annotated_rows.append({
                'chrom': chrom, 'start': s, 'end': e,
                'start_raw': raw_s, 'end_raw': raw_e, 'n_raw_variants': n_raw,
                'support': rec['support'],
                'q_max': round(rec['q_max'], 2), 'q_2nd': round(rec['q_2nd'], 2),
                'adj_indel_l': _pick_label(rec['adj_l']),
                'adj_indel_r': _pick_label(rec['adj_r']),
                'n_adj_indel': rec['n_adj_indel'],
            })
            continue
        chrom_seq = genome.get(chrom, '')
        canon = canonical_in_class(chrom_seq, s, e,
                                   max_shift=cfg.max_ambiguity_shift)
        ann_flag = ann_flags.flag_of(chrom, s, e)
        sh_flag = sh_flags.flag_of(chrom, s, e) if sh_flags else None
        bg_flag = bg_flags.flag_of(chrom, s, e) if bg_flags else None
        flagged = ann_flag or sh_flag
        q = rec['q_max']
        support = rec['support']

        # Length pre-gate (planning/684c fix 1, Kevin-approved 2026-08-17):
        # runs BEFORE the verdict branches so an implausible junction is
        # never eligible for admit_candidate. Recurrence cannot save it — a
        # deterministic aligner artifact reproduces by construction, so the
        # 111 kb LTR class satisfied support>=19 across ten libraries.
        over_len = (e - s) > cfg.max_intron

        if over_len or bg_flag:
            # Background-SV regions (planning/730 W6, Kevin 2026-08-20): the
            # reference is KNOWN wrong there (e.g. R64 chrIII SRD1 flank-A,
            # deleted-and-replaced-by-Ty1 in real strains — yKR888 T2T), so a
            # junction bridging one is a reference artifact regardless of
            # motif, support, or overhang quality. Demotes on BOTH tracks;
            # like every flag it annotates, never deletes.
            verdict = 'demote_orthogonal_evidence'
        elif canon:
            # planning/684c fix 2: `flagged` is consulted on BOTH branches.
            # On the canonical track a set repeat flag caps the verdict at
            # `review` (one-step demotion — the ambiguity class carries a
            # canonical member, a stronger prior than the non-canonical
            # track, so the flag demands review rather than orthogonal
            # evidence). Previously the flag was unreachable here, which is
            # how an LTR-to-LTR artifact with both flags lit was admitted;
            # canonical_in_class is a weak shield exactly there, since one
            # chance GT-AG anywhere in the ambiguity window satisfies it.
            if flagged:
                verdict = 'review'
            elif support >= cfg.min_support and q >= cfg.q_canon:
                verdict = 'admit_candidate'
            else:
                verdict = 'review'
        else:
            if flagged or q < cfg.q_noncanon:
                verdict = 'demote_orthogonal_evidence'
            elif support >= cfg.min_support:
                verdict = 'admit_candidate'
            else:
                verdict = 'review'

        verdict_counts[verdict] = verdict_counts.get(verdict, 0) + 1
        rows.append({
            'chrom': chrom, 'start': s, 'end': e,
            # The key above is the LEFTMOST ambiguity-equivalent coordinate;
            # these are the N-op coordinates most supporting reads carry
            # (ISSUE-016: what a BAM-derived list will be keyed on).
            'start_raw': raw_s, 'end_raw': raw_e, 'n_raw_variants': n_raw,
            'support': support,
            'q_max': round(q, 2), 'q_2nd': round(rec['q_2nd'], 2),
            'canonical_in_class': int(canon),
            'repeat_flag': (ann_flag or '') if have_repeat else TRACK_UNAVAILABLE,
            'selfhom_flag': (int(bool(sh_flag)) if have_selfhom
                             else TRACK_UNAVAILABLE),
            'background_sv_flag': ((bg_flag or '') if have_bg
                                   else TRACK_UNAVAILABLE),
            'over_max_intron': int(over_len),
            # What the census anchor walk stepped over on each side (ISSUE-008:
            # 2H's compensating indel is the signature of a junction it MOVED,
            # so it must survive into the report rather than be absorbed).
            'adj_indel_l': _pick_label(rec['adj_l']),
            'adj_indel_r': _pick_label(rec['adj_r']),
            'n_adj_indel': rec['n_adj_indel'],
            'verdict': verdict,
            'orthogonal_evidence': '',   # v0 placeholder columns for the
            'cross_sample_support': '',  # downstream evidence channels
        })

    order = {'admit_candidate': 0, 'review': 1, 'demote_orthogonal_evidence': 2}
    rows.sort(key=lambda r: (order[r['verdict']], -r['support'], -r['q_max']))
    annotated_rows.sort(key=lambda r: (-r['support'], r['chrom'], r['start']))

    # ISSUE-016: every refused N-op, per RAW coordinate, with why — and whether
    # the same junction was still censused through other reads.
    rows_by_key = {(r['chrom'], r['start'], r['end']): r for r in rows}
    refusal_rows: List[dict] = []
    for (chrom, rs, re_), reasons in ledger.refusals.items():
        chrom_seq = genome.get(chrom, '')
        cs, ce = _canonicalize(chrom_seq, rs, re_, cfg.max_ambiguity_shift)
        ckey = (chrom, cs, ce)
        refusal_rows.append({
            'chrom': chrom, 'start': rs, 'end': re_,
            'canonical_start': cs, 'canonical_end': ce,
            'n_ops': sum(reasons.values()),
            'reasons': ';'.join(f'{k}:{v}' for k, v in sorted(
                reasons.items(), key=lambda kv: (-kv[1], kv[0]))),
            'best_anchor': ledger.best_anchor.get((chrom, rs, re_), 0),
            'annotated': int(ckey in annotated),
            'censused_elsewhere': int(ckey in J),
            'in_table': int(ckey in rows_by_key),
        })
    refusal_rows.sort(key=lambda r: (-r['n_ops'], r['chrom'], r['start']))

    attribution_rows: List[dict] = []
    if attribute:
        attribution_rows = attribute_junctions(
            attribute, genome, cfg, J, rows_by_key, annotated, ledger)

    summary = {
        'bam': str(bam_path),
        'params': asdict(cfg),
        'selfhom_bed': str(selfhom_bed) if selfhom_bed else None,
        'background_sv_bed': (str(background_sv_bed)
                              if background_sv_bed else None),
        'n_junctions_censused': len(J),
        'n_annotated': n_annotated,
        'n_reported': len(rows),
        'verdicts': verdict_counts,
        'repeat_intervals': ann_flags.n_intervals,
        'selfhom_intervals': sh_flags.n_intervals if sh_flags else 0,
        'background_sv_intervals': bg_flags.n_intervals if bg_flags else 0,
        # Which demotion terms were live for this run. An empty list means all
        # three were consulted; anything listed did not run at all.
        'tracks_unavailable': tracks_unavailable,
        'n_annotated_introns_parsed': len(annotated),
        'n_junctions_with_adjacent_indel': sum(
            1 for r in rows if r['n_adj_indel']),
        # ISSUE-016: the census ledger — what the table cannot show.
        'census': {
            'n_reads': ledger.n_reads,
            'reads_skipped': dict(ledger.reads_skipped),
            'n_ops_seen': ledger.n_ops_seen,
            'n_ops_censused': ledger.n_ops_censused,
            'n_ops_refused': ledger.n_ops_refused,
            'n_refused_junctions': len(refusal_rows),
            'n_refused_junctions_censused_elsewhere': sum(
                1 for r in refusal_rows if r['censused_elsewhere']),
            'refusal_reasons': ledger.reason_counts(),
            'attribution': _attribution_counts(attribution_rows),
        },
        # Extra tables. write_pool_gate_outputs() writes each beside the main
        # TSV and replaces the list here with the file path, so the JSON stays
        # small; a caller that never writes gets the rows themselves.
        'tables': {
            'annotated': annotated_rows,
            'census_refusals': refusal_rows,
            'attribution': attribution_rows,
        },
    }
    logger.info('station-c: %d junctions censused (%d annotated); verdicts %s; '
                'N-ops seen %d, refused %d (%s)',
                len(J), n_annotated, verdict_counts, ledger.n_ops_seen,
                ledger.n_ops_refused, ledger.reason_counts() or 'none')
    return rows, summary


#: Attribution statuses, in the order the summary reports them.
ATTRIBUTION_STATUSES = ('reported', 'annotated', 'refused',
                        'annotated_not_seen', 'not_seen', 'chrom_missing')


def _attribution_counts(rows: List[dict]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for r in rows:
        out[r['status']] = out.get(r['status'], 0) + 1
    return {s: out[s] for s in ATTRIBUTION_STATUSES if s in out}


def attribute_junctions(listed: List[Tuple[str, int, int]], genome: Dict[str, str],
                        cfg: PoolGateConfig, J: Dict[Tuple[str, int, int], dict],
                        rows_by_key: Dict[Tuple[str, int, int], dict],
                        annotated: set, ledger: CensusLedger) -> List[dict]:
    """One status per listed junction (RAW N-op coordinates, 0-based half-open).

    Canonicalizes each junction exactly as the census does, then asks the
    census in order: ``reported`` (a row in the table, with its verdict),
    ``annotated`` (censused, annotated — listed in the annotated table),
    ``refused`` (every occurrence failed the anchor gate — reasons attached),
    ``annotated_not_seen`` (annotated, but no read carried an N-op there),
    ``not_seen`` (no N-op at those coordinates in the BAM at all — the list was
    built from a different BAM, or keyed differently), ``chrom_missing``.
    """
    refused_by_canon: Dict[Tuple[str, int, int], Dict[str, int]] = {}
    for (chrom, rs, re_), reasons in ledger.refusals.items():
        chrom_seq = genome.get(chrom, '')
        cs, ce = _canonicalize(chrom_seq, rs, re_, cfg.max_ambiguity_shift)
        d = refused_by_canon.setdefault((chrom, cs, ce), {})
        for k, v in reasons.items():
            d[k] = d.get(k, 0) + v
    out: List[dict] = []
    for chrom, s, e in listed:
        row = {'chrom': chrom, 'start': s, 'end': e,
               'canonical_start': s, 'canonical_end': e,
               'status': 'not_seen', 'verdict': '', 'support': 0,
               'reasons': ''}
        chrom_seq = genome.get(chrom)
        if chrom_seq is None:
            row['status'] = 'chrom_missing'
            out.append(row)
            continue
        cs, ce = _canonicalize(chrom_seq, s, e, cfg.max_ambiguity_shift)
        key = (chrom, cs, ce)
        row['canonical_start'], row['canonical_end'] = cs, ce
        reasons = refused_by_canon.get(key) or ledger.refusals.get((chrom, s, e))
        if key in rows_by_key:
            row['status'] = 'reported'
            row['verdict'] = rows_by_key[key]['verdict']
            row['support'] = rows_by_key[key]['support']
        elif key in J and key in annotated:
            row['status'] = 'annotated'
            row['support'] = J[key]['support']
        elif reasons:
            row['status'] = 'refused'
        elif key in annotated:
            row['status'] = 'annotated_not_seen'
        if reasons:
            row['reasons'] = ';'.join(f'{k}:{v}' for k, v in sorted(
                reasons.items(), key=lambda kv: (-kv[1], kv[0])))
        out.append(row)
    return out


def load_junction_list(path: Path) -> List[Tuple[str, int, int]]:
    """Junctions to attribute, as RAW ``(chrom, start, end)`` triples.

    Accepts a JSON file (a list of ``[chrom, start, end]`` /
    ``{"chrom","start","end"}`` / ``"chrom:start-end"`` items, or a dict whose
    keys or values are such items), a headed TSV with ``chrom``/``start``/``end``
    columns, an ``fpfn_events.tsv`` (``chrom`` + ``new_junction`` ``start-end``,
    rows with an empty ``new_junction`` skipped), or a headerless BED-like
    file (first three columns). Duplicates are dropped, order kept.
    """
    path = Path(path)
    items: List[Tuple[str, int, int]] = []

    def _add(chrom, s, e):
        try:
            items.append((str(chrom), int(s), int(e)))
        except (TypeError, ValueError):
            pass

    def _item(x):
        if isinstance(x, str):
            chrom, _, span = x.partition(':')
            s, _, e = span.partition('-')
            _add(chrom, s, e)
        elif isinstance(x, dict):
            if 'new_junction' in x and x.get('new_junction'):
                s, _, e = str(x['new_junction']).partition('-')
                _add(x.get('chrom'), s, e)
            else:
                _add(x.get('chrom'), x.get('start'), x.get('end'))
        elif isinstance(x, (list, tuple)) and len(x) >= 3:
            _add(x[0], x[1], x[2])

    if path.suffix.lower() == '.json':
        with open(path) as fh:
            data = json.load(fh)
        if isinstance(data, dict):
            for k, v in data.items():
                if isinstance(v, (list, dict)) and not (
                        isinstance(v, list) and len(v) >= 3
                        and not isinstance(v[0], (list, dict, str))):
                    if isinstance(v, list) and v and isinstance(v[0], (list, dict, str)):
                        for x in v:
                            _item(x)
                    else:
                        _item(v)
                else:
                    _item(k)
        else:
            for x in data:
                _item(x)
    else:
        with open(path) as fh:
            lines = [ln.rstrip('\n') for ln in fh if ln.strip() and not ln.startswith('#')]
        if not lines:
            return []
        head = lines[0].split('\t')
        if 'chrom' in head and ('new_junction' in head or ('start' in head and 'end' in head)):
            for ln in lines[1:]:
                f = dict(zip(head, ln.split('\t')))
                _item(f)
        else:
            for ln in lines:
                f = ln.split('\t')
                if len(f) >= 3:
                    _add(f[0], f[1], f[2])
    seen = set()
    out = []
    for it in items:
        if it not in seen:
            seen.add(it)
            out.append(it)
    return out


def write_pool_gate_outputs(rows: List[dict], summary: dict, out_prefix: Path) -> Tuple[Path, Path]:
    """Write ``<prefix>.pool_gate.tsv`` + ``<prefix>.pool_gate.json``.

    ``out_prefix`` is a PREFIX, not a directory: ``-o results/wt_rep1`` yields
    ``results/wt_rep1.pool_gate.tsv``.

    Uses string concatenation rather than ``Path.with_suffix``, which REPLACES
    an existing suffix and is silently lossy on a dotted prefix — ``sample.v2``
    became ``sample.pool_gate.tsv``, so ``sample.v1`` and ``sample.v2`` would
    overwrite each other with no warning. A directory-looking prefix is
    rejected outright rather than writing ``<dir>.pool_gate.tsv`` *beside* the
    directory, which cost two sessions a confused pass each.
    """
    raw = str(out_prefix)
    if raw.endswith(('/', os.sep)) or Path(raw).is_dir():
        raise ValueError(
            f"-o takes an output PREFIX, not a directory: got {raw!r}. "
            f"Include a filename stem, e.g. {raw.rstrip('/' + os.sep)}/<sample> "
            f"-> <sample>.pool_gate.tsv"
        )
    out_prefix = Path(raw)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv = Path(raw + '.pool_gate.tsv')
    js = Path(raw + '.pool_gate.json')
    cols = ['chrom', 'start', 'end', 'support', 'q_max', 'q_2nd',
            'canonical_in_class', 'repeat_flag', 'selfhom_flag',
            'background_sv_flag', 'over_max_intron',
            'adj_indel_l', 'adj_indel_r', 'n_adj_indel', 'verdict',
            'orthogonal_evidence', 'cross_sample_support',
            # ISSUE-016: appended, so column-name readers are unaffected.
            'start_raw', 'end_raw', 'n_raw_variants']
    with open(tsv, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')

    # ISSUE-016 side tables: written beside the main TSV; the JSON keeps the
    # path, not the rows. A summary that never had them (older callers) is
    # written unchanged.
    tables = summary.get('tables')
    if isinstance(tables, dict):
        spec = {
            'annotated': ('.pool_gate.annotated.tsv',
                          ['chrom', 'start', 'end', 'start_raw', 'end_raw',
                           'n_raw_variants', 'support', 'q_max', 'q_2nd',
                           'adj_indel_l', 'adj_indel_r', 'n_adj_indel']),
            'census_refusals': ('.census_refusals.tsv',
                                ['chrom', 'start', 'end', 'canonical_start',
                                 'canonical_end', 'n_ops', 'reasons',
                                 'best_anchor', 'annotated',
                                 'censused_elsewhere', 'in_table']),
            'attribution': ('.attribution.tsv',
                            ['chrom', 'start', 'end', 'canonical_start',
                             'canonical_end', 'status', 'verdict', 'support',
                             'reasons']),
        }
        written = {}
        for name, (suffix, tcols) in spec.items():
            table_rows = tables.get(name)
            if not isinstance(table_rows, list):
                continue
            if name == 'attribution' and not table_rows:
                written[name] = None
                continue
            tpath = Path(raw + suffix)
            with open(tpath, 'w') as fh:
                fh.write('\t'.join(tcols) + '\n')
                for r in table_rows:
                    fh.write('\t'.join(str(r.get(c, '')) for c in tcols) + '\n')
            written[name] = str(tpath)
        summary = dict(summary)
        summary['tables'] = written

    with open(js, 'w') as fh:
        json.dump(summary, fh, indent=1)
    return tsv, js
