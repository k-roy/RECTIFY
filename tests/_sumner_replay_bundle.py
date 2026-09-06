"""Shared loader + production-path replayer for the tester's Sumner replay bundles (ISSUE-026).

The bundles are collaborator data kept OUTSIDE the repository
(``~/work/rectify/dev/sumner_misplaced_panel_20260904/holdout/events/<sha>/``); the test modules that use this
helper skip when their bundle is absent. Bundle format (the tester's d4 format):

  <prefix>.stock.sam      three records per read behind ``#<library> <arm>`` comment lines
                          (stock / baseline_6485226 / <sha>); a read listed under two classes appears twice
  <prefix>.slices.fa      one hg38 slice per listing, in order (``>chrom:start-end offset0based=start``)
  <prefix>.pool.tsv       junction-pool entries near the read (chrom, pool_start, pool_end, annotated flag)
  <prefix>.annotated.tsv  GENCODE v48 basic introns inside the slice (no strand column — inferred from the motif)
  <prefix>.reads.tsv      (optional) read -> class

``replay`` runs the PRODUCTION path on the slice: ``correct_read_3prime`` with the pool ∪ GENCODE candidates
(the exact candidate construction, 2F call, reanchor propagation and writer verdict), then the corrected-BAM
writer (``apply_corrected_edits_to_read``) on the stock record — so a test can compare what 2F REPORTED with
what the writer DREW.
"""
import collections
import copy
import csv
import os
import re

import pysam

EVENTS = os.path.expanduser('~/work/rectify/dev/sumner_misplaced_panel_20260904/holdout/events')
BUNDLES = {
    'e3906ce': ('e3906ce', 't0_e3906ce_regress'),
    'f53d770': ('f53d770', 'f53d770_replay31'),
}
_RCT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def bundle_files(name):
    d, p = BUNDLES[name]
    return [f"{EVENTS}/{d}/{p}.{ext}" for ext in ('stock.sam', 'slices.fa', 'pool.tsv', 'annotated.tsv')]


def bundle_present(name):
    return all(os.path.exists(f) for f in bundle_files(name))


def _rc(s):
    return s.translate(_RCT)[::-1]


def strand_of(seq, s, e):
    """Strand of an intron [s, e) from its splice motif in the slice-local sequence ('.' = unknown)."""
    d, a = seq[s:s + 2], seq[e - 2:e]
    plus = (d, a) in {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
    minus = (_rc(a), _rc(d)) in {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
    return '+' if plus and not minus else ('-' if minus and not plus else '.')


def load_bundle(name):
    """{key: entry}; key = read8, or read8#2 for a read the tester listed twice."""
    sam, fa, pool_p, ann_p = bundle_files(name)
    reads, cur = [], None
    for line in open(sam):
        if line.startswith('@') or not line.strip():
            continue
        if line.startswith('#'):
            cur = tuple(line[1:].split())
            continue
        f = line.rstrip('\n').split('\t')
        if cur[1] == 'stock' or not reads or reads[-1]['read'] != f[0]:
            n_prev = sum(1 for r in reads if r['read'] == f[0])
            reads.append(dict(read=f[0], read8=f[0][:8], key=f[0][:8] + ('' if n_prev == 0 else f'#{n_prev + 1}'),
                              library=cur[0], arms={}))
        reads[-1]['arms'][cur[1]] = line.rstrip('\n')
    slices = []
    for line in open(fa):
        if line.startswith('>'):
            m = re.match(r'>(\S+):(\d+)-(\d+) offset0based=(\d+)', line)
            slices.append(dict(chrom=m.group(1), off=int(m.group(4)), seq=[]))
        else:
            slices[-1]['seq'].append(line.strip())
    for s in slices:
        s['seq'] = ''.join(s['seq']).upper()
    assert len(reads) == len(slices), (len(reads), len(slices))
    pool = collections.defaultdict(list)
    for r in csv.DictReader(open(pool_p), delimiter='\t'):
        pool[r['read']].append((r['chrom'], int(r['pool_start']), int(r['pool_end'])))
    annot = collections.defaultdict(list)
    for r in csv.DictReader(open(ann_p), delimiter='\t'):
        annot[r['read']].append((r['chrom'], int(r['start']), int(r['end'])))
    classes = collections.defaultdict(list)
    reads_tsv = f"{EVENTS}/{BUNDLES[name][0]}/{BUNDLES[name][1]}.reads.tsv"
    if os.path.exists(reads_tsv):
        for r in csv.DictReader(open(reads_tsv), delimiter='\t'):
            classes[r['read'][:8]].append(r['class'])
    out = {}
    for r, sl in zip(reads, slices):
        chrom, off, seq = sl['chrom'], sl['off'], sl['seq']
        ann4 = {(chrom, s - off, e - off, strand_of(seq, s - off, e - off))
                for c, s, e in annot[r['read']] if 0 <= s - off and e - off <= len(seq)}
        pool3 = {(chrom, s - off, e - off) for c, s, e in pool[r['read']] if 0 <= s - off and e - off <= len(seq)}
        out[r['key']] = dict(r, chrom=chrom, off=off, seq=seq, ann4=ann4, pool3=pool3, classes=classes.get(r['read8'], []))
    return out


def nops(read):
    out, pos = [], read.reference_start
    for op, ln in read.cigartuples or []:
        if op == 3:
            out.append((pos, pos + ln))
        if op in (0, 2, 3, 7, 8):
            pos += ln
    return out


def stock_record(entry):
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': entry['chrom'], 'LN': len(entry['seq'])}]})
    a0 = pysam.AlignedSegment.fromstring(entry['arms']['stock'], hdr)
    a0.reference_start = a0.reference_start - entry['off']
    return a0


def correction_from_row(row):
    """bam_writer._load_corrections_from_single_tsv's field map applied to the in-memory row."""
    icp = row.get('five_prime_intron_clip_pos')
    return {
        'corrected_3prime': int(row['corrected_3prime']), 'strand': row['strand'],
        'five_prime_position': row.get('five_prime_position'),
        'five_prime_rescued': bool(row.get('five_prime_rescued')),
        'five_prime_soft_clip': int(row.get('five_prime_soft_clip_length') or 0),
        'five_prime_exon_cigar': row.get('five_prime_exon_cigar') or '',
        'five_prime_upstream_trim': int(row.get('five_prime_upstream_trim') or 0),
        'reanchor_clip_len': int(row.get('reanchor_clip_len') or 0),
        'five_prime_exon2_prefix': int(row.get('five_prime_exon2_prefix') or 0),
        'five_prime_intron_clip_pos': int(icp) if isinstance(icp, int) else -1,
        'sc_homopolymer_extension': 0, 'sc_rescued_seq': '', 'sc_original_softclip_len': 0,
        'oc_homopolymer_extension': 0, 'oc_overcall_count': 0, 'oc_terminal_base': '',
    }


def replay(entry, monkeypatch, gate='report'):
    """Production path on the slice. Returns (row, raw 2F result, materialized record, stock record)."""
    from rectify.core.bam import bam_processor as bp
    from rectify.core.bam import bam_writer as bw
    from rectify.utils.genome import register_genome_contigs

    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', gate)
    monkeypatch.setenv('RECTIFY_2F_CHECK_CONSISTENCY', '1')   # the ISSUE-020 debug invariant raises on violation
    register_genome_contigs([entry['chrom']])
    genome = {entry['chrom']: entry['seq']}
    a0 = stock_record(entry)
    captured = {}
    orig = bp._rescue_3ss

    def wrap(read, genome_, cands, strand, *args, **kw):
        res = orig(read, genome_, cands, strand, *args, **kw)
        captured['res'] = dict(res)
        return res
    monkeypatch.setattr(bp, '_rescue_3ss', wrap)
    rows = bp.correct_read_3prime(copy.deepcopy(a0), genome, annotated_junctions=entry['ann4'],
                                  pool_chrom_index=bp._build_pool_chrom_index(entry['pool3']))
    row = rows[0]
    rec = copy.deepcopy(a0)
    bw.apply_corrected_edits_to_read(rec, correction_from_row(row), genome)
    return row, captured.get('res', {}), rec, a0


def new_nops(rec, stock):
    return [n for n in nops(rec) if n not in nops(stock)]


def body_side_d(rec, new, strand):
    """A D op on the BODY side of a writer-created N-op = the writer's acceptor repair (ISSUE-023 shape)."""
    ct = rec.cigartuples or []
    pos = rec.reference_start
    for i, (op, ln) in enumerate(ct):
        if op == 3 and (pos, pos + ln) in new:
            if strand == '+' and i + 1 < len(ct) and ct[i + 1][0] == 2:
                return True
            if strand == '-' and i > 0 and ct[i - 1][0] == 2:
                return True
        if op in (0, 2, 3, 7, 8):
            pos += ln
    return False


def real(js, off):
    return [(s + off, e + off) for s, e in js]
