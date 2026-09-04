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

SEQ = 'ACGTTGCAAGTCCTGAACGTTAGCCATGGTTACCAGGTTCAAGGTTCACA'   # 50 nt
assert len(SEQ) == 50

# Filler that does not resemble SEQ, so the only exact matches are the ones planted below.
_FILL = ('GGTTCC' * 200)[:1000]


def _with_insert(fill, seq, at):
    return fill[:at] + seq + fill[at + len(seq):]


# chrA carries SEQ[10:] at 100  -> bwa's 10S40M matches the reference EXACTLY
#                                  (40 aligned bases, but a 10-base 5' soft clip: -20)
# chrB carries the whole SEQ at 400 -> bbmap's 50M matches EXACTLY, no clip: 0
# => the fallback scorer must pick bbmap, which is NOT first in dict order.
GENOME = {
    'chrA': _with_insert(_FILL, SEQ[10:], 100),
    'chrB': _with_insert(_FILL, SEQ, 400),
}
assert len(GENOME['chrA']) == 1000 and len(GENOME['chrB']) == 1000


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


def _arm_positions(aligner_reads):
    return {(r.reference_name, r.reference_start) for r in aligner_reads.values()}


def test_cross_contig_record_is_self_consistent():
    """THE invariant: the written (chrom, pos) must be a pair that some ARM produced.

    On `master` the record came out as <loser's chrom> x <winner's pos>, a pair no
    aligner ever reported — and it is only caught downstream when that chromosome
    happens to be too short.  Scorer-independent, so it keeps holding if the
    scoring terms are ever retuned.
    """
    aligner_reads = {
        'bwa': _read('chrA', 100, [(4, 10), (0, 40)]),   # 10S40M, exact on chrA
        'bbmap': _read('chrB', 400, [(0, 50)]),          # 50M, exact on chrB
    }
    written = _run(aligner_reads)
    assert len(written) == 1, "the read must not be dropped"
    rec = written[0]
    assert (rec.reference_name, rec.reference_start) in _arm_positions(aligner_reads), (
        f"fabricated locus {rec.reference_name}:{rec.reference_start} — no arm reported it "
        f"(arms: {sorted(_arm_positions(aligner_reads))})"
    )


def test_cross_contig_fallback_keeps_the_winners_contig():
    """The clean 50M arm on chrB outscores the 10S40M arm on chrA (-20), and it is
    NOT first in dict order — exactly the shape that failed in the field."""
    aligner_reads = {
        'bwa': _read('chrA', 100, [(4, 10), (0, 40)]),
        'bbmap': _read('chrB', 400, [(0, 50)]),
    }
    # Guard the premise: if the scorer stops preferring bbmap this test is vacuous.
    infos = {n: extract_alignment_info(r, n, GENOME) for n, r in aligner_reads.items()}
    from rectify.core.consensus.scoring import score_alignment
    scores = {n: score_alignment(i, GENOME, set()) for n, i in infos.items()}
    assert scores['bbmap'] > scores['bwa'], f"premise broken: {scores}"

    written = _run(aligner_reads)
    assert len(written) == 1
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
