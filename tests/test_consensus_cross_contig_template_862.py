#!/usr/bin/env python3
"""862 — the chimeric-consensus template must come from the WINNER's contig.

`build_chimeric_read` takes ``reference_id`` from the template read and
``reference_start`` from the chimeric/fallback result.  ``select_best_chimeric``
guarantees a single contig only on the true-chimeric path; its cross-contig
FALLBACK (``_fallback_simple_selection``) returns the winner's ``ref_start``
while the arms are still on DIFFERENT contigs.  Picking the template by dict
order then emits <template's chrom> x <winner's position>.

Measured on rectify master @ fd2e2d2 (planning/862), a QuantSeq
`run-all --short-read --dT-primed-cDNA` smoke:
read ``D00689:118:C890GANXX:8:2204:16881:55011`` was placed by bwa at
chrIX:300228 (2S48M) and by bbmap at chrXIII:480619 (1=1X47=1X); the consensus
wrote **chrIX:480619**, past the end of chrIX (439,888 nt).  The run died only
because chrIX is shorter than the borrowed coordinate — with a longer template
contig the fabricated locus is written silently.
"""

import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pysam
import pytest

from rectify.core.consensus.consensus import _process_and_write_batch
from rectify.core.consensus.extract import extract_alignment_info


HEADER = pysam.AlignmentHeader.from_dict({
    'HD': {'VN': '1.6', 'SO': 'unsorted'},
    'SQ': [{'SN': 'chrA', 'LN': 1000}, {'SN': 'chrB', 'LN': 1000}],
})

GENOME = {'chrA': 'ACGT' * 250, 'chrB': 'TGCA' * 250}
SEQ = 'ACGTACGTAC' * 5          # 50 nt


class _FakeOutBam:
    """Minimal stand-in for the output ``pysam.AlignmentFile``."""

    def __init__(self, header):
        self.header = header
        self.written = []

    def get_reference_length(self, name):
        return self.header.get_reference_length(name)

    def write(self, read):
        self.written.append(read)


def _read(chrom, pos, cigartuples):
    r = pysam.AlignedSegment(HEADER)
    r.query_name = 'read1'
    r.query_sequence = SEQ
    r.query_qualities = pysam.qualitystring_to_array('I' * len(SEQ))
    r.flag = 0
    r.reference_id = HEADER.get_tid(chrom)
    r.reference_start = pos
    r.mapping_quality = 60
    r.cigartuples = cigartuples
    return r


def _fresh_stats():
    return {
        'consensus_high': 0, 'consensus_medium': 0, 'consensus_low': 0,
        'chimeric_reads': 0, '5prime_rescued': 0, 'tied_score': 0,
        'by_aligner': Counter(), 'by_aligner_combo': Counter(),
        'impossible_intron_truncated': 0, 'impossible_intron_bp': 0,
        'total_reads': 0,
    }


def _run(aligner_reads):
    """Drive one read through the chimeric branch; return the written records."""
    alignments = {
        name: extract_alignment_info(r, name, GENOME)
        for name, r in aligner_reads.items()
    }
    out = _FakeOutBam(HEADER)
    stats = defaultdict(int)
    stats.update(_fresh_stats())
    _process_and_write_batch(
        read_batch=[('read1', alignments)],
        raw_read_batch=[('read1', aligner_reads)],
        genome=GENOME,
        annotated_junctions=set(),
        out_bam=out,
        stats=stats,
        use_chimeric=True,
    )
    return out.written


def test_cross_contig_fallback_keeps_the_winners_contig():
    """bwa on chrA (5'-clipped, loses) + bbmap on chrB (clean, wins).

    The written record must be chrB:400 — never chrA:400, which is the
    fabricated coordinate the 862 smoke produced.
    """
    # dict order puts the LOSER first, reproducing the field failure exactly.
    aligner_reads = {
        'bwa': _read('chrA', 100, [(4, 10), (0, 40)]),   # 10S40M — 5' clip penalty
        'bbmap': _read('chrB', 400, [(0, 50)]),          # 50M — clean, wins
    }
    written = _run(aligner_reads)
    assert len(written) == 1, "the read must not be dropped"
    rec = written[0]
    assert rec.reference_start == 400
    assert rec.reference_name == 'chrB', (
        f"template contig leaked: got {rec.reference_name}:{rec.reference_start}, "
        "which is the winner's POSITION on the loser's CHROMOSOME"
    )


def test_same_contig_behaviour_is_unchanged():
    """Control: with both arms on one contig the template pool is untouched."""
    aligner_reads = {
        'bwa': _read('chrA', 100, [(4, 10), (0, 40)]),
        'bbmap': _read('chrA', 400, [(0, 50)]),
    }
    written = _run(aligner_reads)
    assert len(written) == 1
    assert written[0].reference_name == 'chrA'


def test_single_aligner_is_unaffected():
    aligner_reads = {'bbmap': _read('chrB', 400, [(0, 50)])}
    written = _run(aligner_reads)
    assert len(written) == 1
    assert written[0].reference_name == 'chrB'
    assert written[0].reference_start == 400
