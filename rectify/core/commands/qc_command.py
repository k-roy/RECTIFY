#!/usr/bin/env python3
"""
RECTIFY QC Command — per-library sequencing QC.

Ported verbatim (logic-for-logic) from the validated prototype
``planning/630b_qc.py`` in the Chanfreau workspace, which was run over all nine
anchor-away DRS libraries.  This is the ONE genuinely new computation in that
unit: ``rectify analyze`` emits no N50, no read-length distribution, no base
quality and no error rate; nothing else in ``rectify/`` computes N50 at all.
Everything else the Analysis tab shows is packaging of ``analyze``'s own tables
(see :mod:`rectify.core.commands.browser_pack_command`).

────────────────────────────────────────────────────────────────────────────────
🔴 THE FIVE RULES THIS MODULE EXISTS TO OBEY — each cost someone real time
────────────────────────────────────────────────────────────────────────────────
 1. N50 is on READ LENGTH, and it says so. "N50" alone is ambiguous.
 2. Base quality is the mean ERROR PROBABILITY converted back to Phred — NOT the
    mean of the Phred values.  Measured: the naive mean overstates by ~8 Q units
    (24.95 vs 16.92 on one read).
 3. 🔴🔴 THE ERROR RATE.  ``NM`` is contaminated TWICE, and the bigger source is
    not the one you would guess:
      (a) INTRONS.  These are ``map-ont`` alignments — not splice-aware — so
          every intron is scored as a long DELETION and lands in NM.  Measured on
          ysh1_rep1 (n=20,000): long D (>=20 bp) is **71.6% of all NM**.  Naive
          NM/aligned-ref = 15.03%; excluding long D = 4.78%.  A 3.1x
          overstatement, and 15% would have been published as "the DRS error
          rate".
      (b) THE POLY(A) TAIL, which for DRS aligns *inside* the alignment over
          genomic A-tracts — it is not in the soft clip.  That is the whole
          reason RECTIFY exists (the 3' end moves on 60.3% of reads).  So
          "aligned blocks only" does NOT by itself remove it; you must exclude
          the walkback interval, original_3prime <-> corrected_3prime.
    ⇒ We report BOTH the naive and the corrected rate, so the contamination
      stays visible rather than being silently fixed.  The corrected rate is the
      headline; the naive one is the evidence.
 4. 5' and 3' soft clips are DIFFERENT QUANTITIES and are never pooled.  A 5'
    clip is predominantly a spliced read whose junction failed to align.  We take
    both from RECTIFY's own strand-aware TSV columns rather than re-deriving from
    CIGAR — re-deriving would risk disagreeing with the pipeline over the strand
    convention, which is exactly the class of second-set-of-numbers error this
    work is trying to avoid.
 5. TWO poly(A) measures, TWO distributions, NEVER merged.  dorado ``pt`` is
    signal-level; RECTIFY's ``polya_length`` is read off the basecall and
    SATURATES.  The gap between them IS the comparison.

And a sixth, from the modality rule: DRS is never UMI-deduplicated, so
``dup_rate`` is ``None`` — the tab must render "n/a", never "0%".

────────────────────────────────────────────────────────────────────────────────
SAMPLING — why a hash and not a reservoir
────────────────────────────────────────────────────────────────────────────────
Reads are selected by ``crc32(read_id) % 10000 < rate``.  This is:
  * deterministic and seed-free — reruns give the identical read set;
  * uniform over the file, so it is NOT chrI-biased the way a head-of-file
    ``head -n`` sample is on a coordinate-sorted BAM;
  * streaming — no reservoir to hold;
  * and crucially THE SAME PREDICATE IN THE BAM AND IN THE TSV, so the two
    passes select the same reads and the join is free.  The walkback-corrected
    error rate needs exactly that join.
``n_sampled`` is reported in the output and belongs on every panel caption.

Usage:
  rectify qc --bam raw.bam --tsv corrected_reads.tsv --sample WTAA_rep1 \\
      --modality DRS --genotype WT-AA --replicate 1 \\
      --annotation genes.gff --out qc/qc_WTAA_rep1.json [--target-reads 200000]

Author: Kevin R. Roy
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
import time
import zlib
from array import array
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple

# NOTE: numpy / pysam are imported INSIDE the functions that need them.  This
# module is imported by ``rectify.cli.create_parser`` on EVERY rectify
# invocation; a module-scope heavy import would tax every command and make the
# whole CLI unusable if pysam were missing.

# A deletion this long or longer is treated as an intron, not an error.  Yeast
# introns are almost all >=50 nt; the 20 bp floor is deliberately conservative
# (it excludes a little real indel signal rather than admitting any intron).  The
# count of these is reported so the choice stays auditable.
INTRON_D_MIN = 20

# Histogram bin edges.  Coarse enough to keep analysis.json small, fine enough to
# see the shapes that matter (the poly(A) saturation shoulder, the 5'-clip mode).
LEN_BINS = [0, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000,
            4000, 5000, 7500, 10000, 15000, 25000]
Q_BINS = list(range(0, 32, 1))
CLIP_BINS = [0, 1, 2, 3, 5, 8, 12, 20, 30, 50, 75, 100, 150, 250, 500, 1000, 5000]
TAIL_BINS = [0, 1, 5, 10, 14, 15, 20, 25, 30, 40, 50, 65, 80, 100, 125, 150, 200,
             300, 500]

MD_RE = re.compile(r'(\d+)|(\^[A-Za-z]+)|([A-Za-z])')

# The manifest header written by ``rectify correct`` (Commit B per-region form).
# Mirrors ``rectify.core.analyze.loaders._MANIFEST_HEADER_COLS``.
_MANIFEST_HEADER_COLS = ['region_id', 'chrom', 'start', 'end', 'tsv_path',
                         'n_rows', 'sha256']

_NULLISH = ('', 'NA', 'nan', 'NaN', 'None', 'null', '.')


# ──────────────────────────────────────────────────────────────────────────────
# Small numeric helpers
# ──────────────────────────────────────────────────────────────────────────────

def hist(values, bins: Sequence[float]) -> List[List[int]]:
    """``[[lower_edge, count], ...]`` — lower edge, so the last bin is an
    open-ended overflow."""
    import numpy as np
    if len(values) == 0:
        return []
    counts, _ = np.histogram(np.asarray(values, dtype=float),
                             bins=list(bins) + [float('inf')])
    return [[bins[i], int(counts[i])] for i in range(len(bins)) if counts[i] > 0]


def n50_of(lengths) -> Optional[int]:
    """N50 ON READ LENGTH: the length L such that reads >= L hold half the
    sequenced bases.  (Rule 1 — this is a read-length N50, not a contig N50.)"""
    import numpy as np
    if len(lengths) == 0:
        return None
    s = np.sort(np.asarray(lengths))[::-1]
    half = s.sum() / 2.0
    return int(s[np.searchsorted(np.cumsum(s), half)])


def mean_q_from_quals(quals) -> Tuple[Optional[float], Optional[float]]:
    """Rule 2: mean ERROR PROBABILITY -> Phred.

    Returns ``(mean_q, naive_mean_of_phred)``.  The mean of Phred values is not a
    mean quality; it is returned only so the overstatement stays visible.
    """
    import numpy as np
    a = np.frombuffer(bytes(quals), dtype=np.uint8).astype(np.float64)
    if a.size == 0:
        return None, None
    p = np.power(10.0, -a / 10.0)
    return float(-10.0 * np.log10(p.mean())), float(a.mean())


def md_mismatch_refpos(md: str, ref_start: int) -> List[int]:
    """Reference positions of every mismatch, from the MD tag.

    MD alternates match-run-lengths, deletions (``^SEQ``) and single mismatch ref
    bases.  Numbers and deletions consume reference; a mismatch consumes one
    reference base.  Insertions do not appear in MD at all, which is why
    insertions are counted from the CIGAR instead.
    """
    pos = ref_start
    out: List[int] = []
    for num, dele, mm in MD_RE.findall(md):
        if num:
            pos += int(num)
        elif dele:
            pos += len(dele) - 1
        elif mm:
            out.append(pos)
            pos += 1
    return out


def overlap(a0: float, a1: float, b0: float, b1: float) -> float:
    return max(0, min(a1, b1) - max(a0, b0))


def _rate_pct(n, d) -> Optional[float]:
    return round(100.0 * n / d, 4) if d else None


# ──────────────────────────────────────────────────────────────────────────────
# Inputs: the coding-gene set, and the corrected-reads TSV (flat or manifest)
# ──────────────────────────────────────────────────────────────────────────────

def coding_genes_from_gff(gff_path) -> Set[str]:
    """Protein-coding gene ids = genes with >=1 ``CDS`` feature.

    Ported from ``planning/630c_prep.py::coding_genes`` so the two agree by
    construction rather than by luck.

    🔴 In the SGD GFF a ``gene`` feature is the ISOFORM UNION and includes UTRs;
    ``CDS`` is the ORF.  CDS ``Parent=`` values arrive as ``<ORF>_mRNA`` /
    ``<ORF>_id001`` (and comma-joined multi-parent forms), so the parent is split
    on ``_`` and the first field taken.
    """
    coding: Set[str] = set()
    p = str(gff_path)
    opener = gzip.open if p.endswith('.gz') else open
    with opener(p, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            for kv in f[8].split(';'):
                if kv.startswith('Parent='):
                    for par in kv[7:].split(','):
                        coding.add(par.split('_')[0].strip())
    return coding


def load_coding_gene_ids(path) -> Set[str]:
    """One gene id per line."""
    with open(str(path)) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def _is_manifest(path) -> bool:
    """True if *path* is a ``corrected_reads.manifest.tsv`` (per-region form).

    🔴 Required columns PRESENT, not header equality — see the matching note in
    ``analyze/loaders._is_manifest``. Exact equality silently reclassified any manifest that
    gained a column as a flat reads TSV, which yields one bogus row and an empty QC block
    (no poly(A), no clips, no untailed rate) with no error anywhere.
    """
    try:
        with open(str(path)) as fh:
            first = fh.readline().rstrip('\n')
        cols = first.split('\t')
        return all(c in cols for c in _MANIFEST_HEADER_COLS)
    except OSError:
        return False


def resolve_corrected_tsvs(path) -> List[Path]:
    """Resolve ``--tsv`` to the list of flat TSVs to stream.

    ``rectify correct`` may emit either a single concatenated
    ``corrected_reads.tsv`` (legacy ``--emit-merged-tsv``) or a
    ``corrected_reads.manifest.tsv`` pointing at per-region shards.  Both are
    accepted; the shards are streamed one at a time (never concatenated in
    memory).
    """
    p = Path(str(path))
    if not _is_manifest(p):
        return [p]
    try:
        from ..bam.tsv_partition import load_manifest
        return [Path(str(e['tsv_path'])) for e in load_manifest(p)]
    except Exception:
        # Fall back to parsing the manifest directly; a manifest we cannot walk
        # must not silently become "zero rows".
        out: List[Path] = []
        with open(p, newline='') as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                tp = (row.get('tsv_path') or '').strip()
                if tp:
                    q = Path(tp)
                    out.append(q if q.is_absolute() else (p.parent / q))
        return out


def _iter_tsv_rows(paths: Sequence[Path]) -> Iterator[dict]:
    for tp in paths:
        with open(str(tp), newline='') as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                yield row


def _num(row: dict, key: str) -> Optional[float]:
    v = row.get(key, '')
    if v is None or v in _NULLISH:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# ──────────────────────────────────────────────────────────────────────────────
# The computation
# ──────────────────────────────────────────────────────────────────────────────

def compute_qc(
    bam_path,
    tsv_path,
    sample: str,
    modality: str = 'DRS',
    genotype: str = '',
    replicate: str = '',
    coding: Optional[Set[str]] = None,
    target_reads: int = 200000,
    max_scan: int = 0,
    verbose: bool = True,
) -> dict:
    """Compute the per-library QC block.  Returns the JSON-ready dict.

    Callable directly (``run-all`` uses this, bypassing argparse) or via
    :func:`run_qc`.
    """
    import numpy as np
    import pysam

    t0 = time.time()

    def log(msg: str) -> None:
        if verbose:
            print(msg, flush=True)

    bam_path = str(bam_path)

    # ── size the sample ───────────────────────────────────────────────────────
    # idxstats is index-only (instant) and counts mapped RECORDS, so it
    # over-counts relative to primary reads.  That is fine here: it is only used
    # to pick the sampling rate, and erring high makes the sample slightly
    # smaller than the target, never larger than memory allows.
    approx = 0
    precounted = False
    try:
        with pysam.AlignmentFile(bam_path, 'rb') as _probe:
            _has_index = _probe.has_index()
    except Exception:
        _has_index = False
    if _has_index:
        try:
            approx = sum(int(l.split('\t')[2])
                         for l in pysam.idxstats(bam_path).strip().split('\n') if l)
        except Exception:
            approx = 0
    if approx == 0 and not max_scan:
        # No index (or an unusable one).  Rather than fall back to "sample
        # everything" — which would be unbounded in memory AND, if capped, would
        # give a head-of-file (chrI-biased) sample, destroying the uniformity
        # that is the whole point of the hash predicate — do one cheap counting
        # pass.  Iterating without per-read work is decompression only.
        try:
            with pysam.AlignmentFile(bam_path, 'rb') as _cnt:
                for _r in _cnt:
                    if not (_r.is_unmapped or _r.is_secondary or _r.is_supplementary):
                        approx += 1
            precounted = True
            log(f'[{sample}] no BAM index: pre-counted {approx:,} primary reads')
        except Exception as exc:                       # pragma: no cover
            approx = 0
            log(f'[{sample}] WARNING: pre-count pass failed ({exc}); '
                f'sampling rate falls back to 10000/10000')

    rate = 10000 if approx <= target_reads or approx == 0 else \
        max(1, int(round(10000.0 * target_reads / approx)))

    def selected(rid: str) -> bool:
        return (zlib.crc32(rid.encode()) % 10000) < rate

    log(f'[{sample}] {"pre-counted" if precounted else "idxstats"} '
        f'~{approx:,} reads; sampling {rate}/10000')

    # Hard memory bound, independent of how ``rate`` was derived.  It should
    # never fire when the rate was sized from a real count; it exists so an
    # unexpectedly large library cannot OOM a 12-hour run-all at the last step.
    hard_cap = max(1, int(target_reads * 1.5))

    # ── pass 1: the BAM (geometry, quality, mismatches) ───────────────────────
    qlen: List[int] = []
    meanq: List[float] = []
    meanq_naive: List[float] = []
    nm_naive_num = nm_naive_den = 0        # the contaminated rate, kept as evidence
    mm_by_read: Dict[str, tuple] = {}      # read_id -> (mismatch refpos, ref_start, ref_end)
    short_indel: Dict[str, int] = {}       # read_id -> short I+D bases
    long_d: Dict[str, List[Tuple[int, int]]] = {}   # read_id -> intron-like deletions
    n_primary = n_sampled = 0
    n_intron_like = 0
    sampling_capped = False

    bam = pysam.AlignmentFile(bam_path, 'rb')
    try:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            n_primary += 1
            if max_scan and n_primary >= max_scan:
                break
            rid = r.query_name
            if not selected(rid):
                continue
            if n_sampled >= hard_cap:
                sampling_capped = True
                continue
            n_sampled += 1

            qlen.append(r.query_length or 0)
            if r.query_qualities is not None:
                mq, mqn = mean_q_from_quals(r.query_qualities)
                if mq is not None:
                    meanq.append(mq)
                    meanq_naive.append(mqn)

            reflen = r.reference_length or 0
            if r.has_tag('NM'):
                nm_naive_num += int(r.get_tag('NM'))
                nm_naive_den += reflen

            # CIGAR: split deletions into intron-like (excluded) and real short
            # indels (counted).
            pos = r.reference_start
            si = 0
            lds: List[Tuple[int, int]] = []
            for op, ln in (r.cigartuples or []):
                if op in (0, 7, 8):        # M / = / X  -- consume both
                    pos += ln
                elif op == 2:              # D
                    if ln >= INTRON_D_MIN:
                        lds.append((pos, pos + ln))
                        n_intron_like += 1
                    else:
                        si += ln
                    pos += ln
                elif op == 3:              # N -- a real splice op; excluded like an intron
                    lds.append((pos, pos + ln))
                    pos += ln
                elif op == 1:              # I -- query only
                    if ln < INTRON_D_MIN:
                        si += ln
            long_d[rid] = lds
            short_indel[rid] = si
            mm_by_read[rid] = (
                array('i', md_mismatch_refpos(r.get_tag('MD'), r.reference_start))
                if r.has_tag('MD') else array('i'),
                r.reference_start,
                r.reference_end if r.reference_end is not None else r.reference_start,
            )
    finally:
        bam.close()

    log(f'[{sample}] BAM pass: {n_primary:,} primary, {n_sampled:,} sampled '
        f'({time.time() - t0:.0f}s)')
    if sampling_capped:
        print(f'[{sample}] WARNING: sampling hit the {hard_cap:,}-read memory cap; '
              f'the sample is truncated and no longer uniform over the file.',
              file=sys.stderr)

    # ── pass 2: the RECTIFY TSV (tails, clips, walkback interval, gene) ───────
    tsv_paths = resolve_corrected_tsvs(tsv_path)
    tsv_seen: Set[str] = set()
    pt_vals: List[float] = []
    bc_vals: List[float] = []
    clip5: List[float] = []
    clip3: List[float] = []
    n_untailed = n_tail_rows = 0
    n_pt0_bc_pos = n_bc0_pt_pos = 0
    n_coding = n_gene_rows = 0
    walkback: Dict[str, Tuple[float, float]] = {}
    n_tsv_rows = n_tsv_sel = 0

    for row in _iter_tsv_rows(tsv_paths):
        n_tsv_rows += 1
        rid = row.get('read_id') or ''
        if not rid or not selected(rid):
            continue
        n_tsv_sel += 1
        tsv_seen.add(rid)

        pt, bc = _num(row, 'pt_tag'), _num(row, 'polya_length')
        if pt is not None:
            pt_vals.append(pt)
        if bc is not None:
            bc_vals.append(bc)
        # Rule 5: BOTH measures must be zero.  They fail independently and
        # disagree a lot, so pt==0 alone materially overstates the untailed rate.
        if pt is not None and bc is not None:
            n_tail_rows += 1
            if pt == 0 and bc == 0:
                n_untailed += 1
            elif pt == 0 and bc > 0:
                n_pt0_bc_pos += 1
            elif bc == 0 and pt > 0:
                n_bc0_pt_pos += 1

        c5 = _num(row, 'five_prime_soft_clip_length')
        c3 = _num(row, 'three_prime_soft_clip_length')
        if c5 is not None:
            clip5.append(c5)
        if c3 is not None:
            clip3.append(c3)

        g = (row.get('gene_id') or '').strip()
        if g and g not in ('NA', 'intergenic', 'None'):
            n_gene_rows += 1
            # A read can list SEVERAL overlapping genes, comma- or pipe-joined,
            # and CDS parents come as <ORF>_mRNA / <ORF>_id001.  A naive
            # ``g in coding`` silently scores every multi-gene read as
            # non-coding.
            if coding is not None and any(
                    x.split('_')[0].strip() in coding
                    for x in g.replace('|', ',').split(',')):
                n_coding += 1

        o3, c3p = _num(row, 'original_3prime'), _num(row, 'corrected_3prime')
        if o3 is not None and c3p is not None and rid in mm_by_read:
            lo, hi = (c3p, o3) if c3p <= o3 else (o3, c3p)
            if hi > lo:
                walkback[rid] = (lo, hi)

    log(f'[{sample}] TSV pass: {n_tsv_rows:,} rows, {n_tsv_sel:,} selected '
        f'({time.time() - t0:.0f}s)')

    # ── the corrected error rate (rule 3) ─────────────────────────────────────
    # Numerator   = mismatches + short indels, excluding anything inside an
    #               intron-like deletion or inside the poly(A) walkback interval.
    # Denominator = aligned reference span, minus intron-like deletions, minus
    #               the walkback interval.
    err_num = err_den = 0
    n_walkback = n_joined = 0
    wb_excluded_mm = 0
    for rid, (mms, rs, re_) in mm_by_read.items():
        if rid in tsv_seen:
            n_joined += 1
        span = max(0, re_ - rs)
        lds = long_d.get(rid, [])
        ld_bp = sum(e - s for s, e in lds)
        wb = walkback.get(rid)
        if wb is not None:
            n_walkback += 1
        wb_bp = overlap(rs, re_, wb[0], wb[1]) if wb else 0
        # a walkback interval can overlap an intron-like deletion; do not
        # subtract it twice
        if wb:
            wb_bp -= sum(overlap(s, e, wb[0], wb[1]) for s, e in lds)
            wb_bp = max(0, wb_bp)
        den = span - ld_bp - wb_bp
        if den <= 0:
            continue
        n_mm = 0
        for p in mms:
            if wb and wb[0] <= p < wb[1]:
                wb_excluded_mm += 1
                continue
            if any(s <= p < e for s, e in lds):
                continue
            n_mm += 1
        err_num += n_mm + short_indel.get(rid, 0)
        err_den += den

    ql = np.asarray(qlen, dtype=float)
    out = {
        'sample': sample,
        'modality': modality,
        'genotype': genotype,
        'replicate': replicate,
        'n_reads_primary': n_primary,
        'n_sampled': n_sampled,
        'sampling_rate_per_10000': rate,
        'sampling_capped': sampling_capped,
        'n_tsv_rows': n_tsv_rows,
        # Two DIFFERENT numbers that are easy to confuse.  The join rate should
        # be ~100%; the walkback rate should be ~60% (the 3' end moves on 60.3%
        # of reads).  Reporting only the second under a "joined" label reads as a
        # 40% join FAILURE.
        'n_joined_bam_tsv': n_joined,
        'n_with_walkback_interval': n_walkback,

        'read_len': {
            'n50': n50_of(ql),
            'median': float(np.median(ql)) if ql.size else None,
            'mean': float(ql.mean()) if ql.size else None,
            'total_bases_sampled': int(ql.sum()) if ql.size else 0,
            'hist': hist(ql, LEN_BINS),
            'note': 'N50 is on READ LENGTH: the length L such that reads >= L '
                    'hold half the sampled bases.',
        },
        'quality': {
            'mean_q': float(np.mean(meanq)) if meanq else None,
            'median_q': float(np.median(meanq)) if meanq else None,
            'hist': hist(meanq, Q_BINS),
            'mean_q_naive_phred_average': float(np.mean(meanq_naive)) if meanq_naive else None,
            'note': ('mean error probability converted back to Phred. The naive '
                     'mean-of-Phred is reported alongside only to show how far it '
                     'overstates.'),
        },
        'error': {
            'rate_pct': _rate_pct(err_num, err_den),
            'naive_nm_rate_pct': _rate_pct(nm_naive_num, nm_naive_den),
            'n_bases_scored': int(err_den),
            'n_intron_like_deletions': n_intron_like,
            'intron_d_min_bp': INTRON_D_MIN,
            'n_mismatches_excluded_by_walkback': wb_excluded_mm,
            'note': ('mismatches (MD) + indels <%d bp, over aligned reference span, '
                     'EXCLUDING deletions >=%d bp (introns -- these are map-ont '
                     'alignments, not splice-aware) and EXCLUDING the poly(A) '
                     'walkback interval original_3prime<->corrected_3prime. '
                     'naive_nm_rate_pct is NM/aligned-ref with neither exclusion, '
                     'kept as evidence.' % (INTRON_D_MIN, INTRON_D_MIN)),
        },
        'clip5': {'hist': hist(clip5, CLIP_BINS), 'n': len(clip5),
                  'median': float(np.median(clip5)) if clip5 else None,
                  'frac_zero': _rate_pct(sum(1 for c in clip5 if c == 0), len(clip5)),
                  'note': "5' soft clip is NOT the same quantity as 3' -- "
                          "predominantly a spliced read whose junction failed to "
                          "align. Never pooled with clip3."},
        'clip3': {'hist': hist(clip3, CLIP_BINS), 'n': len(clip3),
                  'median': float(np.median(clip3)) if clip3 else None,
                  'frac_zero': _rate_pct(sum(1 for c in clip3 if c == 0), len(clip3))},

        'tail_pt': {'hist': hist(pt_vals, TAIL_BINS), 'n': len(pt_vals),
                    'median': float(np.median(pt_vals)) if pt_vals else None,
                    'source': 'dorado pt:i (signal-level)'},
        'tail_rectify': {'hist': hist(bc_vals, TAIL_BINS), 'n': len(bc_vals),
                         'median': float(np.median(bc_vals)) if bc_vals else None,
                         'source': 'rectify polya_length (basecall-level; SATURATES ~14 nt)'},
        'tail_note': 'TWO measures, TWO distributions, never merged. The gap '
                     'between them IS the comparison.',
        'untailed': {
            'frac_pct': _rate_pct(n_untailed, n_tail_rows),
            'n': n_untailed, 'n_scored': n_tail_rows,
            'n_pt0_but_basecall_positive': n_pt0_bc_pos,
            'n_basecall0_but_pt_positive': n_bc0_pt_pos,
            'definition': 'pt_tag == 0 AND rectify polya_length == 0 (both, '
                          'conservatively). The two disagree often, so pt==0 alone '
                          'materially overstates it.',
        },
        'coding': {
            'frac_of_gene_assigned_pct': _rate_pct(n_coding, n_gene_rows) if coding else None,
            'n_gene_assigned': n_gene_rows,
            'caveat': ('gene_id is assigned by OVERLAP, so a read overlapping both '
                       'an ncRNA and a coding gene counts as coding: this is an '
                       'UPPER BOUND.') if coding else None,
        },
        # cDNA only; DRS is never UMI-deduped -> the tab must render n/a, not 0%.
        'dup_rate': None,
        'bam_path': str(bam_path),
        'tsv_path': str(tsv_path),
        'elapsed_s': round(time.time() - t0, 1),
    }
    return out


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def run_qc(args: argparse.Namespace) -> int:
    """Run the QC pass and write the JSON block.  Returns an exit code."""
    coding: Optional[Set[str]] = None
    coding_src = getattr(args, 'coding_genes', None)
    ann = getattr(args, 'annotation', None)
    if coding_src:
        coding = load_coding_gene_ids(coding_src)
        print(f'coding-gene set: {len(coding):,} ids from {coding_src}', flush=True)
    elif ann:
        coding = coding_genes_from_gff(ann)
        print(f'coding-gene set: {len(coding):,} genes with >=1 CDS in {ann}', flush=True)

    out = compute_qc(
        bam_path=args.bam,
        tsv_path=args.tsv,
        sample=args.sample,
        modality=getattr(args, 'modality', 'DRS'),
        genotype=getattr(args, 'genotype', '') or '',
        replicate=str(getattr(args, 'replicate', '') or ''),
        coding=coding,
        target_reads=getattr(args, 'target_reads', 200000),
        max_scan=getattr(args, 'max_scan', 0),
    )

    out_path = Path(str(args.out))
    if out_path.parent and str(out_path.parent):
        out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as fh:
        json.dump(out, fh, indent=1)

    print(f'[{args.sample}] wrote {out_path}  '
          f'err={out["error"]["rate_pct"]}% (naive {out["error"]["naive_nm_rate_pct"]}%)  '
          f'N50={out["read_len"]["n50"]}  untailed={out["untailed"]["frac_pct"]}%',
          flush=True)
    return 0


def create_qc_parser(subparsers) -> argparse.ArgumentParser:
    """Create argument parser for the ``qc`` command."""
    parser = subparsers.add_parser(
        'qc',
        help='Per-library sequencing QC (N50, base quality, error rate, clips, tails)',
        description=(
            "Per-library sequencing QC for one BAM + its RECTIFY corrected-reads TSV. "
            "Emits the JSON block the browser's Analysis tab consumes. "
            "Reports read-length N50, mean base quality as an error-probability "
            "average, an error rate with introns and the poly(A) walkback interval "
            "excluded (alongside the naive NM rate as evidence), 5' and 3' soft-clip "
            "distributions kept separate, and both poly(A) measures kept separate."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('--bam', required=True,
                        help='Alignment BAM for the library')
    parser.add_argument('--tsv', required=True,
                        help='RECTIFY corrected-reads TSV for the SAME library. '
                             'Either a flat corrected_reads.tsv or a '
                             'corrected_reads.manifest.tsv (per-region shards are '
                             'streamed one at a time).')
    parser.add_argument('--sample', required=True,
                        help='Sample name; also the key used in the packed browser JSON')
    parser.add_argument('--modality', default='DRS', choices=['DRS', 'cDNA'],
                        help='Sequencing modality. DRS is never UMI-deduped, so '
                             'dup_rate is reported as null (n/a), never 0%%.')
    parser.add_argument('--genotype', default='', help='Genotype label for the tab')
    parser.add_argument('--replicate', default='', help='Replicate label for the tab')

    parser.add_argument('--coding-genes', dest='coding_genes',
                        help='File with one protein-coding gene id per line. Overrides '
                             '--annotation when both are available.')
    parser.add_argument('--annotation',
                        help='Gene annotation GFF/GFF3 (optionally .gz). The '
                             'protein-coding set is derived as the genes with >=1 CDS '
                             'feature (CDS Parent=<ORF>_mRNA / <ORF>_id001 is split on '
                             '"_"). Used when --coding-genes is absent. NOTE: '
                             '--organism/--Scer fills this in from bundled data, so '
                             '`rectify qc --Scer` derives the coding set from the '
                             'bundled GFF. Omit all three to skip the coding fraction, '
                             'which then reports as not computed.')

    parser.add_argument('--target-reads', dest='target_reads', type=int, default=200000,
                        help='Approximate number of reads to sample. Selection is '
                             'crc32(read_id) %% 10000 < rate: deterministic, uniform '
                             'over the file, and the SAME predicate in the BAM and the '
                             'TSV so the two passes join for free.')
    parser.add_argument('--max-scan', dest='max_scan', type=int, default=0,
                        help='Debug only: stop the BAM pass after N primary reads.')

    parser.add_argument('-o', '--out', required=True,
                        help='Output JSON path. For the browser packer the '
                             'convention is <qc-dir>/qc_<sample>.json.')

    from ...data import add_organism_args
    add_organism_args(parser)

    return parser


def run(args: argparse.Namespace) -> int:
    """Alias so the module matches the ``<x>_command.run(args)`` convention."""
    return run_qc(args)
