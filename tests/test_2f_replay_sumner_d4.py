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


def test_no_bundle_read_lands_by_a_two_mer(monkeypatch):
    """Every placement scores its whole clip: no `edit_distance == 0.0` with a
    comparison shorter than the clip (the ISSUE-017 artefact)."""
    rows, sam, slices = _load()
    for row in rows:
        res, *_ = _replay(row, sam, slices, monkeypatch)
        if res['rescue_type'] != 'intronic_snap':
            assert not (res['edit_distance'] == 0.0 and res['query_bp'] < int(row['clip_len'])), row['read']
