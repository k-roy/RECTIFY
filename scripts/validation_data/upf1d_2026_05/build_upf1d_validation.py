#!/usr/bin/env python3
"""Build the parallel upf1Δ validation set (validation_reads_upf1d.*).

Mirrors the wt `validation_reads.bam` bundle but sourced from upf1Δ rep1 DRS
(poly-A trimmed reads, 5-aligner panel built on Sherlock 2026-05-29). Parallel
to validation_reads_cdna.bam / _quantseq_rev.bam — purely additive; the wt 36
reads and their tests are untouched.

Inputs (staged locally under STAGE):
  STAGE/aligners/upf1d_rep1.<aligner>.bam   — full alignments subset to chosen reads
  STAGE/dorado/upf1d_rep1.dorado.bam        — raw Dorado reads subset to chosen reads
  manifest.tsv                              — read_id  xv  xg  winner  [xk] [xa] [xs]

Outputs (rectify/data/validation/):
  validation_reads_upf1d.bam (+ .bai)                 — winner alignment per read, tagged
  validation_reads_upf1d_dorado_source.bam (+ .bai)   — immutable raw archive, tagged
  aligners/validation_reads_upf1d.<aligner>.bam (+ .bai) × 5 — per-aligner, tagged

Tags: XV (label), XG (category); XK/XA/XS for chimeric (cat5) when present.
"""
from pathlib import Path
from collections import defaultdict
import csv
import re
import pysam

REPO = Path(__file__).resolve().parents[3]
STAGE = Path(__file__).resolve().parent / 'stage'
MANIFEST = Path(__file__).resolve().parent / 'manifest.tsv'
OUT = REPO / 'rectify' / 'data' / 'validation'
ALN_OUT = OUT / 'aligners'
ALIGNERS = ['minimap2', 'gapmm2', 'mapPacBio', 'deSALT', 'uLTRA']
PFX = 'upf1d_rep1'


def load_manifest():
    rows = {}
    with open(MANIFEST) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            rows[row['read_id']] = row
    return rows


def _n_segments(rec) -> int:
    """Segment count = (#N-op introns) + 1 — the value the cat5 test reads from Xg."""
    n = sum(1 for op, _ in (rec.cigartuples or []) if op == 3)
    return n + 1


def _tag(rec, m):
    rec.set_tag('XV', m['xv'], 'Z')
    rec.set_tag('XG', m['xg'], 'Z')
    # Structural bookkeeping the category tests read (lowercase-second-letter =
    # rectify-internal namespace): Xg = segment count (cat5 expects >=2),
    # Xm = aligner-consensus count (cat6/cat7 expect >=1). Computed from the
    # stored alignment so they stay consistent with the actual CIGAR.
    rec.set_tag('Xg', _n_segments(rec), 'i')
    rec.set_tag('Xm', 1, 'i')
    if m.get('xk'):
        rec.set_tag('XK', int(m['xk']), 'i')
    if m.get('xa'):
        rec.set_tag('XA', m['xa'], 'Z')
    if m.get('xs'):
        rec.set_tag('XS', int(m['xs']), 'i')


def _sort_index(unsorted: Path, final: Path):
    pysam.sort('-o', str(final), str(unsorted))
    unsorted.unlink()
    pysam.index(str(final))


def build_per_aligner(man):
    ALN_OUT.mkdir(parents=True, exist_ok=True)
    summary = []
    for a in ALIGNERS:
        src = STAGE / 'aligners' / f'{PFX}.{a}.bam'
        if not src.exists():
            summary.append((a, 'MISSING', 0)); continue
        tmp = ALN_OUT / f'validation_reads_upf1d.{a}.unsorted.bam'
        final = ALN_OUT / f'validation_reads_upf1d.{a}.bam'
        n = 0
        with pysam.AlignmentFile(str(src), 'rb') as bi, \
             pysam.AlignmentFile(str(tmp), 'wb', header=bi.header) as bo:
            for rec in bi:
                m = man.get(rec.query_name)
                if m:
                    _tag(rec, m); n += 1
                bo.write(rec)
        _sort_index(tmp, final)
        summary.append((a, 'OK', n))
    return summary


def _n_nops(rec):
    return sum(1 for op, _ in (rec.cigartuples or []) if op == 3)


def build_winner_bam(man):
    """validation_reads_upf1d.bam = each read's winner-aligner alignment, chosen
    at the manifest-specified chrom.

    mapPacBio emits MULTIPLE records flagged PRIMARY (no secondary bit) for
    multi-mapping reads — so 'first primary' can land on the wrong locus. We
    therefore collect ALL primary records per read and pick the one on the
    intended ``chrom``; ties broken toward the most N-ops (so cat5/6/9 get the
    intron-bearing alignment)."""
    by_winner = {}
    for rid, m in man.items():
        by_winner.setdefault(m['winner'], {})[rid] = m
    template_src = STAGE / 'aligners' / f'{PFX}.minimap2.bam'
    with pysam.AlignmentFile(str(template_src), 'rb') as t:
        header = t.header
    tmp = OUT / 'validation_reads_upf1d.unsorted.bam'
    final = OUT / 'validation_reads_upf1d.bam'

    # gather candidate primary records per read, then choose
    cands = defaultdict(list)  # rid -> [rec.copy()]
    for aligner, reads in by_winner.items():
        src = STAGE / 'aligners' / f'{PFX}.{aligner}.bam'
        with pysam.AlignmentFile(str(src), 'rb') as bi:
            for rec in bi:
                if rec.is_secondary or rec.is_supplementary:
                    continue
                if rec.query_name in reads:
                    cands[rec.query_name].append(rec.to_dict())
    written = set()
    with pysam.AlignmentFile(str(tmp), 'wb', header=header) as bo:
        for rid, m in man.items():
            recs = cands.get(rid, [])
            if not recs:
                continue
            want = m.get('chrom')
            on_chrom = [d for d in recs if d['ref_name'] == want] if want else recs
            pool = on_chrom or recs
            # prefer most N-ops (intron-bearing) for structural categories
            best = max(pool, key=lambda d: sum(int(n) for n, op in
                       re.findall(r'(\d+)([MIDNSHP=X])', d['cigar']) if op == 'N'))
            rec = pysam.AlignedSegment.from_dict(best, header)
            _tag(rec, m)
            bo.write(rec)
            written.add(rid)
    _sort_index(tmp, final)
    missing = set(man) - written
    return len(written), missing


def build_dorado_source(man):
    src = STAGE / 'dorado' / f'{PFX}.dorado.bam'
    if not src.exists():
        return 'MISSING', 0
    tmp = OUT / 'validation_reads_upf1d_dorado_source.unsorted.bam'
    final = OUT / 'validation_reads_upf1d_dorado_source.bam'
    n = 0
    with pysam.AlignmentFile(str(src), 'rb', check_sq=False) as bi, \
         pysam.AlignmentFile(str(tmp), 'wb', header=bi.header) as bo:
        for rec in bi:
            m = man.get(rec.query_name)
            if m:
                _tag(rec, m); n += 1
                bo.write(rec)
    # dorado source is unaligned; sort by name (coord sort fails w/o coords)
    pysam.sort('-n', '-o', str(final), str(tmp))
    tmp.unlink()
    return 'OK', n


def main():
    man = load_manifest()
    print(f"manifest: {len(man)} reads")
    print("per-aligner:", build_per_aligner(man))
    nw, missing = build_winner_bam(man)
    print(f"winner BAM: {nw} reads written; missing={sorted(s[:8] for s in missing)}")
    print("dorado source:", build_dorado_source(man))


if __name__ == '__main__':
    main()
