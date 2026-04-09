"""
Test: no duplicate primary records in the rectified BAM output.

Bug 4: gapmm2 PAF records were written as primary alignments (FLAG=0/16)
even when another aligner already held the primary for the same read_id.
This created 2-3+ primary records per read, causing double-counting in
downstream tools.

Fix locations:
  - rectify/core/multi_aligner.py  _paf_to_bam(): check tp:A: tag; mark
    non-primary PAF records with FLAG |= 0x100 (secondary).
  - rectify/core/consensus.py  _process_and_write_batch(): clear FLAG bits
    0x900 on the winning record to guarantee it is written as primary.

Author: Kevin R. Roy
"""

import io
import os
import tempfile
from typing import Dict
from unittest.mock import patch

import pysam
import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_header():
    """Return a minimal pysam AlignmentHeader for test BAMs."""
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'unsorted'},
        'SQ': [{'SN': 'chrI', 'LN': 230218}],
    })


def _make_read(header, read_name: str, flag: int, pos: int = 100) -> pysam.AlignedSegment:
    """Return a minimal mapped pysam.AlignedSegment."""
    seg = pysam.AlignedSegment(header)
    seg.query_name = read_name
    seg.flag = flag
    seg.reference_id = 0        # chrI
    seg.reference_start = pos
    seg.mapping_quality = 60
    seg.cigarstring = '50M'
    seg.query_sequence = 'A' * 50
    seg.query_qualities = pysam.qualitystring_to_array('I' * 50)
    return seg


def _write_bam(path: str, reads, header):
    """Write a list of AlignedSegment objects to a name-sorted BAM."""
    with pysam.AlignmentFile(path, 'wb', header=header) as bam:
        for r in reads:
            bam.write(r)
    # Name-sort (needed by _iter_name_grouped_bams)
    sorted_path = path + '.ns.bam'
    pysam.sort('-n', '-o', sorted_path, path)
    os.replace(sorted_path, path)


# ---------------------------------------------------------------------------
# Unit test: _paf_to_bam correctly marks non-primary PAF records secondary
# ---------------------------------------------------------------------------

def test_paf_to_bam_marks_secondary_records():
    """
    When gapmm2 emits two PAF lines for the same read_id (one primary, one
    secondary as indicated by tp:A:S), _paf_to_bam must set FLAG|=0x100 on
    the secondary record so _filtered_read_iterator excludes it.
    """
    from rectify.core.multi_aligner import _paf_to_bam

    # Minimal PAF: two lines for the same read_id.
    # Fields: qname qlen qst qend strand tname tlen tst tend nmatch alen mapq [tags...]
    # tp:A:P = primary, tp:A:S = secondary
    paf_content = (
        "read1\t100\t0\t100\t+\tchrI\t230218\t1000\t1100\t98\t100\t60\ttp:A:P\tcs:Z::100\n"
        "read1\t100\t0\t100\t+\tchrI\t230218\t5000\t5100\t90\t100\t30\ttp:A:S\tcs:Z::100\n"
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        paf_path = os.path.join(tmpdir, 'test.paf')
        bam_path = os.path.join(tmpdir, 'test.bam')
        genome_path = os.path.join(tmpdir, 'genome.fa')

        with open(paf_path, 'w') as fh:
            fh.write(paf_content)

        # Create a minimal genome FASTA so pysam.FastaFile can read it
        with open(genome_path, 'w') as fh:
            fh.write('>chrI\n' + 'A' * 230218 + '\n')
        pysam.faidx(genome_path)

        from pathlib import Path
        _paf_to_bam(Path(paf_path), Path(bam_path), genome_path, threads=1)

        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            records = list(bam)

        assert len(records) == 2, f"Expected 2 records, got {len(records)}"

        primary_records = [r for r in records if not (r.flag & 0x900)]
        secondary_records = [r for r in records if r.flag & 0x100]

        assert len(primary_records) == 1, (
            f"Expected 1 primary record, got {len(primary_records)}. "
            f"Flags: {[r.flag for r in records]}"
        )
        assert len(secondary_records) == 1, (
            f"Expected 1 secondary record, got {len(secondary_records)}. "
            f"Flags: {[r.flag for r in records]}"
        )


# ---------------------------------------------------------------------------
# Integration test: rectified BAM writer produces at most one primary per read
# ---------------------------------------------------------------------------

def test_rectified_bam_no_duplicate_primaries():
    """
    Simulate two aligner BAMs for the same read_id: aligner_a provides a
    primary alignment and aligner_b also provides a primary alignment (the
    gapmm2 bug scenario).  After consensus selection only the winning record
    should be primary; any other written records must be secondary.

    This test exercises _process_and_write_batch() directly to confirm the
    FLAG &= ~0x900 guard works independently of the _paf_to_bam fix.
    """
    from rectify.core.consensus import _process_and_write_batch, extract_alignment_info

    header = _make_header()

    # Both aligners produce primary FLAG reads for the same read_id
    read_a = _make_read(header, 'read1', flag=0, pos=100)   # + strand primary
    read_b = _make_read(header, 'read1', flag=0, pos=100)   # + strand primary

    # Provide chrI in the genome dict so extract_alignment_info does not call
    # standardize_chrom_name() (which has an unrelated pre-existing signature
    # mismatch: it only accepts 1 arg but consensus.py passes 2).
    genome: Dict[str, str] = {'chrI': 'A' * 230218}

    aligner_reads = {'aligner_a': read_a, 'aligner_b': read_b}
    alignments = {}
    for aligner, read in aligner_reads.items():
        alignments[aligner] = extract_alignment_info(read, aligner, genome)

    read_batch = [('read1', alignments)]
    raw_read_batch = [('read1', aligner_reads)]

    stats = {
        'consensus_high': 0,
        'consensus_medium': 0,
        'consensus_low': 0,
        '5prime_rescued': 0,
        'tied_score': 0,
        'by_aligner': {},
        'by_aligner_combo': {},
    }
    from collections import defaultdict
    stats['by_aligner'] = defaultdict(int)
    stats['by_aligner_combo'] = defaultdict(int)

    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = os.path.join(tmpdir, 'rectified.bam')
        with pysam.AlignmentFile(out_path, 'wb', header=header) as out_bam:
            _process_and_write_batch(
                read_batch, raw_read_batch,
                genome, None, out_bam, stats
            )

        with pysam.AlignmentFile(out_path, 'rb') as bam:
            records = list(bam)

    # Exactly one record written (we only write the winner)
    assert len(records) == 1, f"Expected 1 record in output BAM, got {len(records)}"

    # That record must be primary (FLAG & 0x900 == 0)
    assert records[0].flag & 0x900 == 0, (
        f"Winning record is not primary: FLAG={records[0].flag:#06x}"
    )


# ---------------------------------------------------------------------------
# Test: winner already marked secondary gets promoted to primary
# ---------------------------------------------------------------------------

def test_consensus_winner_promoted_from_secondary():
    """
    If the winning aligner's record has FLAG=0x100 (secondary) set — the
    scenario where gapmm2 wins but its record arrived from _paf_to_bam with
    the wrong flag — _process_and_write_batch must clear 0x900 and write it
    as primary.
    """
    from rectify.core.consensus import _process_and_write_batch, extract_alignment_info

    header = _make_header()

    # Winning aligner's read is incorrectly flagged secondary (bug scenario)
    read_winner = _make_read(header, 'read1', flag=0x100, pos=100)

    genome: Dict[str, str] = {'chrI': 'A' * 230218}

    aligner_reads = {'gapmm2': read_winner}
    alignments = {'gapmm2': extract_alignment_info(read_winner, 'gapmm2', genome)}

    read_batch = [('read1', alignments)]
    raw_read_batch = [('read1', aligner_reads)]

    from collections import defaultdict
    stats = {
        'consensus_high': 0, 'consensus_medium': 0, 'consensus_low': 0,
        '5prime_rescued': 0, 'tied_score': 0,
        'by_aligner': defaultdict(int), 'by_aligner_combo': defaultdict(int),
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = os.path.join(tmpdir, 'rectified.bam')
        with pysam.AlignmentFile(out_path, 'wb', header=header) as out_bam:
            _process_and_write_batch(
                read_batch, raw_read_batch,
                genome, None, out_bam, stats
            )

        with pysam.AlignmentFile(out_path, 'rb') as bam:
            records = list(bam)

    assert len(records) == 1
    assert records[0].flag & 0x900 == 0, (
        f"Winner not promoted to primary: FLAG={records[0].flag:#06x}"
    )
