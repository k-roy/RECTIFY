"""Stored NET query evidence survives H and never invents ambiguous bases."""
from dataclasses import asdict
import gzip
import csv

import pysam
import pytest

from rectify.core.bam.netseq_bam_processor import process_netseq_read, process_netseq_bam
from rectify.core.netseq import netseq_rescue as nr
from rectify.core.netseq_cpa.pileup import walkback_pileup


def edge_read(strand, kind='rescue', hard=0, encoded=False):
    bases = list(('CAGTCGGCATGC' * 100)[:1000])
    bases[300:302], bases[418:420], bases[498:500] = 'GT', 'AG', 'AG'
    prefix = 'NNNNNN' if kind == 'unknown' else 'GCTGCA'
    bases[420:426], bases[500:506] = 'CCCCCC', prefix
    if kind == 'tail':
        bases[294], bases[295:300] = 'C', 'AAAAA'
        prefix = 'AAA' + 'TCGGCT'
    forward = ''.join(bases)
    query = 'TGAC' + forward[260:300] + prefix
    ops = [(5, 5), (4, 4), (0, 40), (4, len(prefix))]
    if hard:
        ops += [(5, hard)]
    genome = forward
    if strand == '-':
        genome, query, ops = nr.revcomp(forward), nr.revcomp(query), ops[::-1]
    header = pysam.AlignmentHeader.from_references(['netedge'], [1000])
    read = pysam.AlignedSegment(header)
    read.query_name, read.reference_id = f'{strand}_{kind}', 0
    read.reference_start, read.flag = (260, 16) if strand == '+' else (700, 0)
    read.mapping_quality, read.cigartuples = 60, ops
    read.query_sequence, read.query_qualities = query, [31] * len(query)
    read.set_tag('ZZ', 'preserve')
    if encoded:
        seq = list(query)
        for q, p in read.get_aligned_pairs(matches_only=True):
            assert seq[q] == genome[p]
            seq[q] = '='
        read.query_sequence = ''.join(seq)
        read.query_qualities = [31] * len(seq)
    coords = [(300, 420), (300, 500)] if strand == '+' else [(580, 700), (500, 700)]
    pool = nr.JunctionPool([nr.make_junction('netedge', a, b, strand) for a, b in coords])
    return read, {'netedge': genome}, pool


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('kind', ['rescue', 'tail'])
@pytest.mark.parametrize('encoded', [False, True])
def test_hardclip_keeps_retained_clip_evidence_and_metadata(strand, kind, encoded):
    rows = []
    for hard in (0, 7):
        read, genome, pool = edge_read(strand, kind, hard, encoded)
        before = read.to_string()
        row = process_netseq_read(read, 'netedge', genome=genome, junction_pool=pool,
                                  umi_length=6 if kind == 'tail' else 0)
        assert read.to_string() == before
        assert row.five_prime_soft_clip_length == 4
        assert row.three_prime_soft_clip_length == (9 if kind == 'tail' else 6)
        if kind == 'rescue':
            assert row.rescue_status == 'spliced_rescued' and row.rescue_k == 6
            assert row.three_prime_corrected == (505 if strand == '+' else 494)
        else:
            assert row.tail_len == 8 and row.tail_walkback == 5
        record = asdict(row)
        record.pop('cigar_summary')
        rows.append(record)
    assert rows[0] == rows[1]


@pytest.mark.parametrize('strand', ['+', '-'])
def test_unknown_bases_do_not_support_an_acceptor(strand):
    read, genome, pool = edge_read(strand, 'unknown')
    before = read.to_string()
    row = process_netseq_read(read, 'netedge', genome=genome, junction_pool=pool)
    assert row.rescue_status == 'exon1_end' and row.rescue_k == 0
    assert row.three_prime_corrected == row.three_prime_raw
    assert read.to_string() == before
    assert nr.longest_common_prefix('ACN', 'ACN') == 2
    assert nr.longest_common_prefix('ACR', 'ACR') == 2
    assert nr.longest_common_prefix('acgt', 'ACGT') == 4


@pytest.mark.parametrize('strand', ['+', '-'])
def test_missing_sequence_is_counted_without_tail_or_rescue(tmp_path, strand):
    read, genome, pool = edge_read(strand)
    read.cigartuples, read.query_sequence = [(0, 42)], None
    if strand == '-':
        read.reference_start = 698
    before = read.to_string()
    row = process_netseq_read(read, 'netedge', genome=genome, junction_pool=pool,
                              walkback_requires_clip_a=False)
    assert row.rescue_status == 'none' and row.rescue_k == 0
    assert row.tail_len == 0 and row.three_prime_corrected == row.three_prime_raw
    assert read.to_string() == before
    path = tmp_path / 'missing.bam'
    with pysam.AlignmentFile(str(path), 'wb', header=read.header) as out:
        out.write(read)
    streamed = list(process_netseq_bam(path, genome=genome, junction_pool=pool,
                                      walkback_requires_clip_a=False, show_progress=False))
    assert len(streamed) == 1 and asdict(streamed[0]) == asdict(row)


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('encoded', [False, True])
def test_cpa_counts_retained_softclip_inside_hardclip(tmp_path, strand, encoded):
    outputs = []
    for hard in (0, 7):
        read, genome, _ = edge_read(strand, 'tail', hard, encoded)
        reference, path, result = [tmp_path / f'{label}{hard}' for label in ('ref.fa', 'read.bam', 'out.tsv.gz')]
        reference.write_text('>netedge\n' + genome['netedge'] + '\n')
        pysam.faidx(str(reference))
        with pysam.AlignmentFile(str(path), 'wb', header=read.header) as out:
            out.write(read)
        initial = path.read_bytes()
        stats = walkback_pileup(path, reference, result, sample='edge')
        assert path.read_bytes() == initial
        with gzip.open(result, 'rt') as inp:
            rows = list(csv.DictReader(inp, delimiter='\t'))
        assert stats['reads_seen'] == stats['reads_used'] == 1
        assert len(rows) == 1 and int(rows[0]['sum_oaNT']) == 3
        outputs.append((stats, rows))
    assert outputs[0] == outputs[1]
