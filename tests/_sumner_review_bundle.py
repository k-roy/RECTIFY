"""Loader + production-path replayer for the fixer's Sumner REVIEW bundles (STAGING.md §1 layout; ISSUE-028).

A review bundle is a directory of BAMs — ``stock.bam`` (the events-tier input record) plus one ``<sha>.bam`` per
arm — with a ``manifest.tsv`` (read, library, chrom, strand, class, ...), next to a ``<name>.slices.fa`` (hg38
slices, ``>chrom:start-end offset0based=start``) and a ``<name>.slices.gtf`` (the GENCODE v48 basic records inside
the slices). Collaborator data, kept OUTSIDE the repository
(``~/work/rectify/dev/sumner_misplaced_panel_20260904/holdout/events/<sha>/``); the tests that use this helper skip
when the bundle is absent.

Unlike the tester's d4-format bundles (``tests/_sumner_replay_bundle``) there is no pool file: the candidates are
the GENCODE introns of the slice plus the read's own N-ops (``correct_read_3prime`` adds those itself). Every read
of the 4993253 review that gained a rescue gained it on an ANNOTATED junction, so the annotated candidate set
reproduces the reviewed placements.
"""
import collections
import copy
import csv
import os
import re

import pysam

EVENTS = os.path.expanduser('~/work/rectify/dev/sumner_misplaced_panel_20260904/holdout/events')
BUNDLES = {
    '4993253': ('4993253', 'review_2f_4993253'),
}


def bundle_paths(name):
    d, p = BUNDLES[name]
    base = f"{EVENTS}/{d}"
    return f"{base}/{p}/stock.bam", f"{base}/{p}/manifest.tsv", f"{base}/{p}.slices.fa", f"{base}/{p}.slices.gtf"


def bundle_present(name):
    return all(os.path.exists(p) for p in bundle_paths(name))


def _load_slices(fa):
    slices = []
    for line in open(fa):
        if line.startswith('>'):
            m = re.match(r'>(\S+):(\d+)-(\d+) offset0based=(\d+)', line)
            slices.append(dict(chrom=m.group(1), off=int(m.group(4)), seq=[]))
        else:
            slices[-1]['seq'].append(line.strip())
    for s in slices:
        s['seq'] = ''.join(s['seq']).upper()
        s['end'] = s['off'] + len(s['seq'])
    return slices


def _load_introns(gtf):
    """{chrom: {(intron_start, intron_end, strand)}} from the GTF's exon records, per transcript (0-based half-open)."""
    ex = collections.defaultdict(list)
    for line in open(gtf):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        ex[(f[0], m.group(1), f[6])].append((int(f[3]) - 1, int(f[4])))
    introns = collections.defaultdict(set)
    for (chrom, _tid, strand), e in ex.items():
        e.sort()
        for (_s1, e1), (s2, _e2) in zip(e, e[1:]):
            if s2 > e1:
                introns[chrom].add((e1, s2, strand))
    return introns


def load_bundle(name):
    """{read8: entry} with the stock record (real coordinates), its slice, and the slice-local annotated introns."""
    bam, manifest, fa, gtf = bundle_paths(name)
    slices = _load_slices(fa)
    introns = _load_introns(gtf)
    rows = {}
    with open(manifest) as fh:
        lines = [ln for ln in fh if not ln.startswith('#')]
    for r in csv.DictReader(lines, delimiter='\t'):
        rows[r['read'][:8]] = r
    out = {}
    for a in pysam.AlignmentFile(bam):
        if a.is_secondary or a.is_supplementary or a.is_unmapped:
            continue
        read8 = a.query_name[:8]
        if read8 not in rows:
            continue
        cands = [s for s in slices if s['chrom'] == a.reference_name
                 and s['off'] <= a.reference_start and a.reference_end <= s['end']]
        if not cands:
            continue
        sl = min(cands, key=lambda s: s['end'] - s['off'])
        chrom, off, seq = sl['chrom'], sl['off'], sl['seq']
        ann4 = {(chrom, s - off, e - off, st) for s, e, st in introns[chrom]
                if s - off >= 0 and e - off <= len(seq)}
        out[read8] = dict(rows[read8], read=a.query_name, read8=read8, sam=a.to_string(), chrom=chrom, off=off,
                          seq=seq, ann4=ann4, strand=rows[read8]['strand'])
    return out


def stock_record(entry):
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'},
                                           'SQ': [{'SN': entry['chrom'], 'LN': len(entry['seq'])}]})
    a0 = pysam.AlignedSegment.fromstring(entry['sam'], hdr)
    a0.reference_start = a0.reference_start - entry['off']
    return a0


def nops(read):
    out, pos = [], read.reference_start
    for op, ln in read.cigartuples or []:
        if op == 3:
            out.append((pos, pos + ln))
        if op in (0, 2, 3, 7, 8):
            pos += ln
    return out


def real(js, off):
    return [(s + off, e + off) for s, e in js]


def replay(entry, monkeypatch=None, gate='report'):
    """Production path on the slice. Returns (row, raw 2F result, materialized record, stock record).

    ``monkeypatch`` may be None (a script caller): the environment is then set directly."""
    from rectify.core.bam import bam_processor as bp
    from rectify.core.bam import bam_writer as bw
    from rectify.utils.genome import register_genome_contigs
    from tests._sumner_replay_bundle import correction_from_row

    if monkeypatch is not None:
        monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', gate)
        monkeypatch.setenv('RECTIFY_2F_CHECK_CONSISTENCY', '1')
    else:
        os.environ['RECTIFY_2F_NOVEL_GATE'] = gate
        os.environ['RECTIFY_2F_CHECK_CONSISTENCY'] = '1'
    register_genome_contigs([entry['chrom']])
    genome = {entry['chrom']: entry['seq']}
    a0 = stock_record(entry)
    captured = {}
    orig = bp._rescue_3ss

    def wrap(read, genome_, cands, strand, *args, **kw):
        res = orig(read, genome_, cands, strand, *args, **kw)
        captured['res'] = dict(res)
        return res
    if monkeypatch is not None:
        monkeypatch.setattr(bp, '_rescue_3ss', wrap)
    else:
        bp._rescue_3ss = wrap
    try:
        rows = bp.correct_read_3prime(copy.deepcopy(a0), genome, annotated_junctions=entry['ann4'],
                                      pool_chrom_index=None)
    finally:
        if monkeypatch is None:
            bp._rescue_3ss = orig
    row = rows[0]
    rec = copy.deepcopy(a0)
    bw.apply_corrected_edits_to_read(rec, correction_from_row(row), genome)
    return row, captured.get('res', {}), rec, a0
