"""ISSUE-094: consensus must not let an arm placement that leaves its contig compete.

deSALT indexes the genome as one concatenated sequence and reports a boundary-crossing alignment on the earlier
contig. The July Path-A consensus chose 875 such records (all `Xa=deSALT`, lone low-confidence winners) over the
minimap2/gapmm2 placements, and every tool downstream had to cope with invalid SAM. Consensus now drops any arm
record whose span runs past its contig (by the BAM header's length) before scoring, so a valid arm wins.
"""

import random

import pysam

from rectify.core.consensus.consensus import run_consensus_selection
from rectify.utils.genome import register_genome_contigs

LEN = 1000


def _genome():
    rng = random.Random(94)
    g = {'chrA': ''.join(rng.choice('ACGT') for _ in range(LEN)),
         'chrB': ''.join(rng.choice('ACGT') for _ in range(LEN))}
    register_genome_contigs(g)
    return g


def _arm(path, header, recs):
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for name, chrom, start, cigar, seq in recs:
            r = pysam.AlignedSegment(header)
            r.query_name = name
            r.reference_name = chrom
            r.reference_start = start
            r.flag = 0
            r.mapping_quality = 60
            r.cigarstring = cigar
            r.query_sequence = seq
            r.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
            out.write(r)


def test_a_valid_arm_wins_over_a_placement_past_the_contig_end(tmp_path):
    genome = _genome()
    header = pysam.AlignmentHeader.from_references(['chrA', 'chrB'], [LEN, LEN])
    # The molecule really is chrB 100..159 (60 nt). deSALT reports the first 40 nt on chrA past its end,
    # as a crossing of its concatenated index would; minimap2 places it correctly on chrB.
    seq = genome['chrB'][100:160]
    _arm(tmp_path / 'deSALT.bam', header, [('m1', 'chrA', 980, '60M', seq)])
    _arm(tmp_path / 'minimap2.bam', header, [('m1', 'chrB', 100, '60M', seq)])
    out = tmp_path / 'consensus.bam'
    stats = run_consensus_selection({'deSALT': str(tmp_path / 'deSALT.bam'),
                                     'minimap2': str(tmp_path / 'minimap2.bam')}, genome, str(out), n_workers=1)
    with pysam.AlignmentFile(str(out)) as bam:
        recs = list(bam)
    assert [(r.query_name, r.reference_name, r.reference_start) for r in recs] == [('m1', 'chrB', 100)]
    assert stats['arm_records_past_contig_dropped'] == 1
    assert stats['reads_all_arms_past_contig'] == 0


def test_a_read_with_no_valid_arm_is_not_written(tmp_path):
    genome = _genome()
    header = pysam.AlignmentHeader.from_references(['chrA', 'chrB'], [LEN, LEN])
    seq = genome['chrB'][100:160]
    _arm(tmp_path / 'deSALT.bam', header, [('m1', 'chrA', 980, '60M', seq)])
    out = tmp_path / 'consensus.bam'
    stats = run_consensus_selection({'deSALT': str(tmp_path / 'deSALT.bam')}, genome, str(out), n_workers=1)
    with pysam.AlignmentFile(str(out)) as bam:
        assert list(bam) == []
    assert stats['reads_all_arms_past_contig'] == 1
