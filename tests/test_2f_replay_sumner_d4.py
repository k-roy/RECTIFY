"""T0 replay fixture — the tester's ISSUE-017 read bundle (Sumner human RNA004 DRS).

Replays real reads that landed 4 nt into the intron (the GTRAGT +5 GT decoy)
against their real annotated / pool candidates on genome slices. The records
are collaborator data and live OUTSIDE the repository
(``dev/sumner_misplaced_panel_20260904/holdout/events/269ebe7/``); the test
skips when they are absent. Expectations are the ISSUE-017 (slice fix + floor
+ annotated tie-break) outcomes; ISSUE-020 (anchored-score ranking) will
tighten them to "all annotated at the tester's extended EDs".
"""
import csv
import os
import re

import pysam
import pytest

from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation
from rectify.utils.genome import register_genome_contigs

E = os.path.expanduser('~/work/rectify/dev/sumner_misplaced_panel_20260904/holdout/events/269ebe7')
BUNDLE = [E + '/d4_replay10.tsv', E + '/d4_replay10.stock.sam', E + '/d4_replay10.slices.fa']

pytestmark = pytest.mark.skipif(not all(os.path.exists(p) for p in BUNDLE),
                                reason='Sumner d4 replay bundle not present (collaborator data, kept outside the repo)')


def _load():
    slices, name = {}, None
    for line in open(BUNDLE[2]):
        if line.startswith('>'):
            m = re.match(r'>(\S+):(\d+)-(\d+) offset0based=(\d+)', line)
            name = (m.group(1), int(m.group(4)))
            slices[name] = []
        else:
            slices[name].append(line.strip())
    slices = {k: ''.join(v).upper() for k, v in slices.items()}
    sam = {}
    for line in open(BUNDLE[1]):
        if line.startswith('#') or line.startswith('@') or not line.strip():
            continue
        sam[line.split('\t')[0]] = line.rstrip('\n')
    rows = list(csv.DictReader(open(BUNDLE[0]), delimiter='\t'))
    return rows, sam, slices


def _replay(row, sam, slices, monkeypatch):
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', 'report')
    chrom, strand = row['chrom'], row['strand']
    ns, ne, ba = int(row['novel_start']), int(row['novel_end']), int(row['best_annot'])
    key = [k for k in slices if k[0] == chrom and k[1] <= min(ns, ba)
           and k[1] + len(slices[k]) >= max(ne, ba)][0]
    off, seq = key[1], slices[key]
    register_genome_contigs([chrom])
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': chrom, 'LN': len(seq)}]})
    a = pysam.AlignedSegment.fromstring(sam[row['read']], hdr)
    a.reference_start = a.reference_start - off
    novel = (chrom, ns - off, ne - off)
    annotated = (chrom, ba - off, ne - off) if strand == '+' else (chrom, ns - off, ba - off)
    res = rescue_3ss_truncation(a, {chrom: seq}, {novel, annotated}, strand,
                                annotated_junctions={annotated})
    j = res.get('rescued_junction')
    return res, (j[0], j[1] + off, j[2] + off) if j else None, (ns, ne), \
        ((ba, ne) if strand == '+' else (ns, ba))


# What the 017 tree does per read: 'annotated' | 'within1' (≤ 1 nt of the annotated site) | 'novel'
EXPECT = {
    'c887bc16': 'novel',      # the terminal peel still reaches the +4 window — ISSUE-020's population
    '31fab950': 'within1',
    'c64bc988': 'annotated',
    '5b20c72a': 'within1',
    'c41c7314': 'annotated',
    'de84a10a': 'annotated',
    'a5a5a1bb': 'within1',
    '7f41e755': 'annotated',
    '6fc67f58': 'annotated',
    '885f6430': 'annotated',
}


@pytest.mark.parametrize('read8', sorted(EXPECT))
def test_bundle_read(read8, monkeypatch):
    rows, sam, slices = _load()
    row = [r for r in rows if r['read'].startswith(read8)][0]
    res, junction, novel, annotated = _replay(row, sam, slices, monkeypatch)
    assert res['rescued'], row['read']
    # Never the 1–2-mer artefact: the comparison covered the whole clip.
    assert res['query_bp'] >= int(row['clip_len']) or res['rescue_type'] == 'intronic_snap'
    if EXPECT[read8] == 'annotated':
        assert junction[1:] == annotated, (junction, annotated)
        assert res['landing_annotated'] is True
    elif EXPECT[read8] == 'within1':
        assert junction[1:] != novel
        assert abs(junction[1] - annotated[0]) + abs(junction[2] - annotated[1]) <= 1
    else:
        assert junction[1:] == novel


def test_exon_cigar_query_span_matches_the_writer_contract(monkeypatch):
    """The exon CIGAR must consume exactly the query bases the BAM writer will
    convert, or `extend` falls back to a flat M block and the evidence ops never
    reach the BAM. Contract (splice_aware_5prime, "Which sequence the exon CIGAR
    is sized from"): 5' end INSIDE the rescued intron -> the intron-mapped run
    (clip + bases mapped past the boundary; `reroute` demands exon_q ==
    n_intronic_q, no upstream trim); 5' end outside -> clip + upstream_trim.
    At 3834686 the equivalence extension borrowed the overshoot base a second
    time on 5b20c72a / a5a5a1bb / 31fab950 (span 22 for clip 20 + trim 1)."""
    from rectify.core.align.local_aligner import cigar_str_to_ops
    from rectify.core.splice.splice_aware_5prime import _get_intronic_query_bases

    rows, sam, slices = _load()
    checked = 0
    for row in rows:
        strand = row['strand']
        res, junction, *_ = _replay(row, sam, slices, monkeypatch)
        if not res.get('rescued') or res.get('rescue_type') != 'softclip' or not res.get('five_prime_exon_cigar'):
            continue
        ops = cigar_str_to_ops(res['five_prime_exon_cigar'])
        span = sum(l for op, l in ops if op in (0, 1, 4, 7, 8))
        # Rebuild the read the same way _replay did, to measure the contract.
        chrom = row['chrom']
        ns, ne, ba = int(row['novel_start']), int(row['novel_end']), int(row['best_annot'])
        key = [k for k in slices if k[0] == chrom and k[1] <= min(ns, ba)
               and k[1] + len(slices[k]) >= max(ne, ba)][0]
        off, seq = key[1], slices[key]
        hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': chrom, 'LN': len(seq)}]})
        a = pysam.AlignedSegment.fromstring(sam[row['read']], hdr)
        a.reference_start = a.reference_start - off
        clip = a.cigartuples[0][1] if (strand == '+' and a.cigartuples[0][0] == 4) else \
            (a.cigartuples[-1][1] if (strand == '-' and a.cigartuples[-1][0] == 4) else 0)
        j_start, j_end = junction[1] - off, junction[2] - off
        five = a.reference_start if strand == '+' else a.reference_end - 1
        inside = j_start <= five < j_end
        trim = int(res.get('five_prime_upstream_trim', 0) or 0)
        rcl = int(res.get('reanchor_clip_len', 0) or 0)
        if rcl > 0:
            # Terminal peel / re-anchor: the writer's pre-pass converts the
            # 5' bases through the re-anchor point into a soft clip of length
            # reanchor_clip_len, and `extend` converts exactly that clip.
            assert span == rcl + trim, (row['read'], res['five_prime_exon_cigar'], span, rcl, trim)
        elif inside:
            intronic = _get_intronic_query_bases(a, j_start if strand == '-' else j_end, strand)
            assert trim == 0, (row['read'], trim)
            assert span == len(intronic), (row['read'], res['five_prime_exon_cigar'], span, len(intronic))
        else:
            assert span == clip + trim, (row['read'], res['five_prime_exon_cigar'], span, clip, trim)
        checked += 1
    assert checked >= 3, checked


def test_no_bundle_read_lands_by_a_two_mer(monkeypatch):
    """Every placement scores its whole clip: no `edit_distance == 0.0` with a
    comparison shorter than the clip (the ISSUE-017 artefact)."""
    rows, sam, slices = _load()
    for row in rows:
        res, *_ = _replay(row, sam, slices, monkeypatch)
        if res['rescue_type'] != 'intronic_snap':
            assert not (res['edit_distance'] == 0.0 and res['query_bp'] < int(row['clip_len'])), row['read']
