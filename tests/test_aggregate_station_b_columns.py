"""ISSUE-040 — the junction table carries station B's provenance and its equally good alternatives.

Kevin, 2026-09-07: "keep the other one as a noted equally good alternative as a auxiliary Sam tag,
also perhaps in our junction parquet/tsv." The SAM tags are written by `bam_writer`
(`Xb` = versioned successful calls); this is the junction-level half — `rectify aggregate`
reads those tags off the BAM so a junction table says which junctions exist because a micro-exon was
drawn, and what the reads behind them could equally have been.

Credit belongs only to the new junctions in the successful call. Untagged reads
are not automatically independent stock evidence; other stages may have moved them.
"""
import json
import pysam

from rectify.core.aggregate.junctions import aggregate_junctions


def _bam(tmp_path, reads):
    p = tmp_path / 'x.bam'
    hdr = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrT', 'LN': 10_000}]}
    with pysam.AlignmentFile(str(p), 'wb', header=hdr) as out:
        for name, cigar, start, tags in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name = name
            a.reference_name = 'chrT'
            a.reference_start = start
            a.cigartuples = cigar
            a.mapping_quality = 60
            a.query_sequence = 'A' * sum(n for op, n in cigar if op in (0, 1, 4, 7, 8))
            for k, v in tags.items():
                a.set_tag(k, v)
            out.write(a)
    pysam.index(str(p))
    return str(p)


CIG = [(0, 50), (3, 50), (0, 6), (3, 144), (0, 50)]


def _tag(alternative):
    # Literal schema fixture: chosen exon [100,106) splits intron [50,250).
    return {'Xb': json.dumps({'v': 1, 'chrom': 'chrT', 'strand': '+', 'calls': [
        {'intron': [50, 250], 'exons': [[100, 106]], 'alternatives': alternative}]})}


def test_station_b_columns_are_emitted_even_when_the_station_never_ran(tmp_path):
    df = aggregate_junctions(_bam(tmp_path, [('r1', CIG, 0, {})]))
    assert list(df['station_b_reads']) == [0, 0]
    assert list(df['station_b_unverified_reads']) == [0, 0]
    assert list(df['station_b_alternatives']) == ['', '']


def test_a_junction_drawn_by_station_b_is_marked_with_its_alternatives(tmp_path):
    reads = [
        ('drawn', CIG, 0, _tag('chrT:150-156')),
        ('drawn2', CIG, 0, _tag('chrT:150-156')),
        ('untagged', CIG, 0, {}),
    ]
    df = aggregate_junctions(_bam(tmp_path, reads))
    assert list(df['full_junction_reads']) == [3, 3]
    assert list(df['station_b_reads']) == [2, 2]
    assert list(df['station_b_alternatives']) == ['chrT:150-156'] * 2


def test_alternatives_are_deduplicated_and_capped(tmp_path):
    reads = [(f'r{i}', CIG, 0, _tag(f'chrT:{200 + i}-{206 + i}'))
             for i in range(8)]
    df = aggregate_junctions(_bam(tmp_path, reads))
    alts = df.iloc[0]['station_b_alternatives'].split(';')
    assert len(alts) == 4                              # capped so one wild read cannot bloat a table
    assert len(set(alts)) == 4
