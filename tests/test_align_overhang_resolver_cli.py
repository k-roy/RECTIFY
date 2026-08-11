"""CLI exposure of the overhang resolver as a --junction-aligners panel member.

The resolver (planning/641/644) is a post-pass on the finished minimap2 arm,
not an external aligner: `rectify align --junction-aligners overhang_resolver`
must parse, must run the post-pass after the aligner loop, and must fail LOUD
when minimap2 is not part of the same invocation (a silent skip would feed
downstream junction-pool prescans an unresolved pool with exit 0). Its output
SUBSTITUTES the minimap2 arm downstream (passthrough-or-rewrite, same read
set) — carrying both arms would buy a duplicate `correct` pass over
~98%-identical records (planning/669 §1).

RN passthrough matters here: the consensus K-way merge keys on RN:i, and the
resolver arm inherits its records from the RN-injected minimap2 BAM.
"""

import argparse

import pysam

from rectify.core.commands.align_command import create_align_parser
from rectify.core.align.overhang_resolver import run_overhang_resolver


def _parse(argv):
    root = argparse.ArgumentParser()
    sub = root.add_subparsers(dest='command')
    create_align_parser(sub)
    return root.parse_args(argv)


def test_parser_accepts_overhang_resolver(tmp_path):
    args = _parse([
        'align', 'reads.fastq', '--genome', 'g.fa', '-o', str(tmp_path),
        '--aligners', 'minimap2',
        '--junction-aligners', 'uLTRA', 'deSALT', 'overhang_resolver',
    ])
    assert 'overhang_resolver' in args.junction_aligners


def test_resolver_passthrough_preserves_rn_and_count(tmp_path):
    """A no-clip BAM streams through unchanged: same records, RN intact."""
    genome_path = tmp_path / 'g.fa'
    genome_path.write_text('>chrI\n' + 'ACGTACGTGG' * 40 + '\n')

    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': 'chrI', 'LN': 400}],
    })
    bam_in = tmp_path / 'mm2.bam'
    with pysam.AlignmentFile(bam_in, 'wb', header=header) as out:
        for i, name in enumerate(['read-a', 'read-b']):
            r = pysam.AlignedSegment(header)
            r.query_name = name
            r.reference_id = 0
            r.reference_start = 10 + 50 * i
            r.mapping_quality = 60
            r.cigarstring = '20M'
            r.query_sequence = 'ACGTACGTGGACGTACGTGG'
            r.query_qualities = pysam.qualitystring_to_array('I' * 20)
            r.set_tag('RN', i + 1, value_type='i')
            out.write(r)

    bam_out = tmp_path / 'resolver.bam'
    run_overhang_resolver(str(bam_in), str(genome_path), str(bam_out))

    with pysam.AlignmentFile(bam_out, 'rb', check_sq=False) as f:
        recs = list(f)
    assert [r.query_name for r in recs] == ['read-a', 'read-b']
    assert [r.get_tag('RN') for r in recs] == [1, 2]
