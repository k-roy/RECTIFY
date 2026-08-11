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
  candidate to needs-orthogonal-evidence; it never discards.
- **Two-track admission** (dev/STATIONC_MAPPACBIO_HARVEST_20260810.md): the
  canonical-in-class and non-canonical tracks never share a threshold.

Verdicts (non-annotated junctions):

==============================  =============================================
``admit_candidate``             canonical track: support>=2 AND q>=q_canon;
                                non-canonical track: unflagged AND support>=2
                                AND q>=q_noncanon
``review``                      evidence below the admit bar but nothing
                                against it (low support and/or low q)
``demote_orthogonal_evidence``  non-canonical AND (repeat-flagged OR
                                q < q_noncanon): admissible only with
                                orthogonal evidence (short-read/COMPASS,
                                cross-sample recurrence, mm2-side distress)
==============================  =============================================

Cross-sample recurrence and the orthogonal-evidence channels are recorded as
columns for downstream tooling but not computed in v0 (single-sample scope).
"""

from __future__ import annotations

import gzip
import json
import logging
from dataclasses import asdict, dataclass
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


@dataclass
class PoolGateConfig:
    min_anchor: int = 8          # census: min adjacent exon anchor per N-op
    overhang_cap: int = 60       # junction-adjacent query bases assessed
    err_bits: float = 2.0        # bits removed per alignment error in window
    q_canon: float = 40.0        # canonical-track admit threshold (644h)
    q_noncanon: float = 80.0     # non-canonical-track admit threshold (644h/i)
    min_support: int = 2         # recurrence gate within this sample
    repeat_margin: int = 500     # bp around an annotated repeat feature
    max_ambiguity_shift: int = 30


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


def census_bam(bam_path: str, genome: Dict[str, str], cfg: PoolGateConfig,
               max_q_reads_per_junction: int = 50) -> Dict[Tuple[str, int, int], dict]:
    """One streaming pass: per canonical junction, support + q statistics.

    ``max_q_reads_per_junction`` caps the per-read overhang scoring (q_max is
    a max — extra reads beyond the cap still count toward support).
    """
    J: Dict[Tuple[str, int, int], dict] = {}
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as fh:
        for read in fh.fetch(until_eof=True):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or not read.cigartuples):
                continue
            chrom = read.reference_name
            chrom_seq = genome.get(chrom)
            if chrom_seq is None:
                continue
            offs = None
            ops = read.cigartuples
            ref = read.reference_start
            for i, (op, ln) in enumerate(ops):
                if op == 3:
                    lf = ops[i - 1][1] if i >= 1 and ops[i - 1][0] in (0, 7, 8) else 0
                    rt = (ops[i + 1][1]
                          if i + 1 < len(ops) and ops[i + 1][0] in (0, 7, 8) else 0)
                    if min(lf, rt) >= cfg.min_anchor:
                        s, e = _canonicalize(chrom_seq, ref, ref + ln,
                                             cfg.max_ambiguity_shift)
                        key = (chrom, s, e)
                        rec = J.get(key)
                        if rec is None:
                            rec = J[key] = {'support': 0, 'q_max': 0.0,
                                            'q_2nd': 0.0, 'q_scored': 0}
                        rec['support'] += 1
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


def load_annotated_canonical(annotation_path: Path, genome: Dict[str, str],
                             cfg: PoolGateConfig) -> set:
    """Annotated introns (intron features + inferred from exon adjacency),
    ambiguity-canonicalised so census keys match."""
    import re
    from collections import defaultdict
    exons = defaultdict(list)
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
                m = re.search(r'Parent=([^;]+)', attrs)
                if m:
                    for parent in m.group(1).split(','):
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


def pool_gate(
    bam_path: str,
    genome: Dict[str, str],
    annotation_path: Path,
    cfg: Optional[PoolGateConfig] = None,
    selfhom_bed: Optional[Path] = None,
) -> Tuple[List[dict], dict]:
    """Run Station C v0 over one BAM. Returns (rows, summary).

    Rows cover NON-annotated junctions only (one dict per junction, sorted by
    verdict then descending support); annotated junctions are counted in the
    summary. Report-only: nothing upstream is modified.
    """
    cfg = cfg or PoolGateConfig()
    ann_flags = load_repeat_flags(annotation_path, margin=cfg.repeat_margin)
    sh_flags = load_selfhom_bed(selfhom_bed) if selfhom_bed else None
    annotated = load_annotated_canonical(annotation_path, genome, cfg)

    J = census_bam(bam_path, genome, cfg)

    rows: List[dict] = []
    n_annotated = 0
    verdict_counts: Dict[str, int] = {}
    for (chrom, s, e), rec in J.items():
        if (chrom, s, e) in annotated:
            n_annotated += 1
            continue
        chrom_seq = genome.get(chrom, '')
        canon = canonical_in_class(chrom_seq, s, e,
                                   max_shift=cfg.max_ambiguity_shift)
        ann_flag = ann_flags.flag_of(chrom, s, e)
        sh_flag = sh_flags.flag_of(chrom, s, e) if sh_flags else None
        flagged = ann_flag or sh_flag
        q = rec['q_max']
        support = rec['support']

        if canon:
            if support >= cfg.min_support and q >= cfg.q_canon:
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
            'support': support,
            'q_max': round(q, 2), 'q_2nd': round(rec['q_2nd'], 2),
            'canonical_in_class': int(canon),
            'repeat_flag': ann_flag or '',
            'selfhom_flag': int(bool(sh_flag)),
            'verdict': verdict,
            'orthogonal_evidence': '',   # v0 placeholder columns for the
            'cross_sample_support': '',  # downstream evidence channels
        })

    order = {'admit_candidate': 0, 'review': 1, 'demote_orthogonal_evidence': 2}
    rows.sort(key=lambda r: (order[r['verdict']], -r['support'], -r['q_max']))

    summary = {
        'bam': str(bam_path),
        'params': asdict(cfg),
        'selfhom_bed': str(selfhom_bed) if selfhom_bed else None,
        'n_junctions_censused': len(J),
        'n_annotated': n_annotated,
        'n_reported': len(rows),
        'verdicts': verdict_counts,
        'repeat_intervals': ann_flags.n_intervals,
        'selfhom_intervals': sh_flags.n_intervals if sh_flags else 0,
    }
    logger.info('station-c: %d junctions censused (%d annotated); verdicts %s',
                len(J), n_annotated, verdict_counts)
    return rows, summary


def write_pool_gate_outputs(rows: List[dict], summary: dict, out_prefix: Path) -> Tuple[Path, Path]:
    """Write ``<prefix>.pool_gate.tsv`` + ``<prefix>.pool_gate.json``."""
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv = out_prefix.with_suffix('.pool_gate.tsv')
    js = out_prefix.with_suffix('.pool_gate.json')
    cols = ['chrom', 'start', 'end', 'support', 'q_max', 'q_2nd',
            'canonical_in_class', 'repeat_flag', 'selfhom_flag', 'verdict',
            'orthogonal_evidence', 'cross_sample_support']
    with open(tsv, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')
    with open(js, 'w') as fh:
        json.dump(summary, fh, indent=1)
    return tsv, js
