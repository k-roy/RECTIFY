#!/usr/bin/env python3
"""864 — a consensus record's RNAME, POS, CIGAR and STRAND must all come from
the SAME candidate.

`chimeric_consensus.build_chimeric_read` assembles a FRESH `AlignedSegment`.
Before this fix it took `reference_id`, the strand bits and `mapping_quality`
from a "template" read chosen by `_process_and_write_batch` as the first entry
of `aligner_reads` whose SEQ length matched — dict order, no reference to who
won — while `reference_start` and the CIGAR came from the winning candidate.

Whenever the candidates disagreed on contig or strand (exactly the case that
routes to `_fallback_simple_selection`) the emitted record was a splice of two
different alignments.  Measured on rectify master @ fd2e2d2:

  * TruSeq/COMPASS `run-all --short-read`, PE 7 arms (planning/861 §S8b):
    25 / 20,000 records = 0.125 % on the wrong chromosome at the right
    coordinate; SE 5 arms 115 / 157,475 = 0.073 %.  Only ~8 % overran the
    borrowed contig and raised — the rest were written silently.
    Real record, `SRR40431829.159` R1: bbmap won with chrI:165474 `49=1X25=`;
    `HISAT2_default` (first in the dict) supplied chrX, MAPQ 60, FLAG 89 and
    its entire aux-tag set.  `Xa` credited bbmap.
  * QuantSeq bbmap+bwa (planning/862): 2,532 / 179,062 reads took that path.

planning/864 supersedes the contig-only fix on `fix/runall-quantseq-862`
(whose tests are kept, in `test_consensus_cross_contig_template_862.py`).
"""

import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pysam
import pytest

from rectify.core.consensus.consensus import _process_and_write_batch
from rectify.core.consensus.extract import extract_alignment_info
from rectify.core.consensus.scoring import score_alignment


HEADER = pysam.AlignmentHeader.from_dict({
    'HD': {'VN': '1.6', 'SO': 'unsorted'},
    'SQ': [{'SN': 'chrA', 'LN': 1000}, {'SN': 'chrB', 'LN': 1000}],
})

SEQ = 'ACGTTGCAAGTCCTGAACGTTAGCCATGGTTACCAGGTTCAAGGTTCACA'   # 50 nt
assert len(SEQ) == 50

_FILL = ('GGTTCC' * 200)[:1000]


def _with_insert(fill, seq, at):
    return fill[:at] + seq + fill[at + len(seq):]


# chrA: SEQ[10:] planted at 100 -> a 10S40M arm matches EXACTLY but pays the
#       5' soft-clip penalty.
# chrB: the whole SEQ planted at 400 -> a 50M arm matches exactly, no clip.
#       SEQ[10:] therefore also sits at 410 on chrB.
GENOME = {
    'chrA': _with_insert(_FILL, SEQ[10:], 100),
    'chrB': _with_insert(_FILL, SEQ, 400),
}
assert len(GENOME['chrA']) == 1000 and len(GENOME['chrB']) == 1000


class _FakeOutBam:
    def __init__(self, header):
        self.header = header
        self.written = []

    def get_reference_length(self, name):
        return self.header.get_reference_length(name)

    def write(self, read):
        self.written.append(read)


def _read(chrom, pos, cigartuples, *, mapq=60, is_reverse=False, tags=()):
    r = pysam.AlignedSegment(HEADER)
    r.query_name = 'read1'
    r.query_sequence = SEQ
    r.query_qualities = pysam.qualitystring_to_array('I' * len(SEQ))
    r.flag = 16 if is_reverse else 0
    r.reference_id = HEADER.get_tid(chrom)
    r.reference_start = pos
    r.mapping_quality = mapq
    r.cigartuples = cigartuples
    for name, value, vtype in tags:
        r.set_tag(name, value, value_type=vtype)
    return r


def _run(aligner_reads):
    alignments = {
        name: extract_alignment_info(r, name, GENOME)
        for name, r in aligner_reads.items()
    }
    out = _FakeOutBam(HEADER)
    stats = defaultdict(int)
    stats.update({
        'consensus_high': 0, 'consensus_medium': 0, 'consensus_low': 0,
        'chimeric_reads': 0, '5prime_rescued': 0, 'tied_score': 0,
        'by_aligner': Counter(), 'by_aligner_combo': Counter(),
        'impossible_intron_truncated': 0, 'total_reads': 0,
        'placement_invariant_violations': 0,
    })
    _process_and_write_batch(
        read_batch=[('read1', alignments)],
        raw_read_batch=[('read1', aligner_reads)],
        genome=GENOME,
        annotated_junctions=set(),
        out_bam=out,
        stats=stats,
        use_chimeric=True,
    )
    return out.written, stats


# The loser is deliberately FIRST in insertion order and carries a DIFFERENT
# MAPQ, so a record that borrows anything from it is detectable field by field.
def _cross_contig_arms():
    return {
        'loser_first': _read('chrA', 100, [(4, 10), (0, 40)], mapq=60),
        'winner': _read('chrB', 400, [(0, 50)], mapq=7),
    }


def _premise(aligner_reads, winner='winner'):
    infos = {n: extract_alignment_info(r, n, GENOME) for n, r in aligner_reads.items()}
    scores = {n: score_alignment(i, GENOME, set()) for n, i in infos.items()}
    losers = [n for n in scores if n != winner]
    assert all(scores[winner] > scores[n] for n in losers), \
        f"premise broken, the intended winner does not win: {scores}"


def test_every_placement_field_comes_from_the_winning_candidate():
    """RNAME, POS, CIGAR, strand AND MAPQ must all be the winner's.

    On master this fails on `reference_name` (chrA leaked from the template)
    and on `mapping_quality` (60 leaked from the template).
    """
    arms = _cross_contig_arms()
    _premise(arms)
    written, stats = _run(arms)
    assert len(written) == 1, "the read must not be dropped"
    rec = written[0]
    win = arms['winner']

    assert rec.reference_name == win.reference_name, (
        f"contig leaked from the template: {rec.reference_name}:{rec.reference_start} "
        f"is the winner's POSITION on the loser's CHROMOSOME"
    )
    assert rec.reference_start == win.reference_start
    assert rec.cigartuples == win.cigartuples
    assert rec.is_reverse == win.is_reverse
    assert rec.mapping_quality == win.mapping_quality, (
        "MAPQ leaked from the template — the record claims the loser's "
        "confidence for the winner's placement"
    )
    assert rec.get_tag('Xa') == 'winner'
    assert stats['placement_invariant_violations'] == 0


def test_emitted_locus_was_reported_by_some_arm():
    """Scorer-independent: the (chrom, pos) pair must be one an arm produced."""
    arms = _cross_contig_arms()
    written, _ = _run(arms)
    reported = {(r.reference_name, r.reference_start) for r in arms.values()}
    rec = written[0]
    assert (rec.reference_name, rec.reference_start) in reported, (
        f"fabricated locus {rec.reference_name}:{rec.reference_start}; "
        f"arms reported {sorted(reported)}"
    )


def test_strand_comes_from_the_winner_not_the_template():
    """Same contig, DIFFERENT strands -> also the disagreement fallback.

    On master `out.flag = template_read.flag` copied the template's 0x10 bit, so
    the record carried the loser's orientation with the winner's coordinates —
    and because a BAM SEQ is stored reference-forward, the emitted SEQ then did
    not correspond to its own CIGAR.
    """
    arms = {
        'loser_first': _read('chrB', 410, [(4, 10), (0, 40)], is_reverse=False),
        'winner': _read('chrB', 400, [(0, 50)], is_reverse=True),
    }
    _premise(arms)
    written, stats = _run(arms)
    assert len(written) == 1
    rec = written[0]
    assert rec.is_reverse is True, "strand leaked from the template read"
    assert (rec.reference_name, rec.reference_start) == ('chrB', 400)
    assert rec.cigartuples == [(0, 50)]
    assert rec.query_sequence == SEQ, "SEQ must stay reference-forward"
    assert stats['placement_invariant_violations'] == 0


def test_template_tags_are_not_borrowed_from_a_foreign_placement():
    """The aux tags must come from the record that supplied the placement.

    On master the output carried the losing arm's tags while `Xa` credited the
    winner — the provenance lie that made the 861 forensics possible in the
    first place (`SRR40431829.159` shipped HISAT2's ZS/XN/XM/XO/XG/YT/NH).
    """
    arms = {
        'loser_first': _read('chrA', 100, [(4, 10), (0, 40)],
                             tags=(('ZZ', 'loser', 'Z'),)),
        'winner': _read('chrB', 400, [(0, 50)],
                        tags=(('ZZ', 'winner', 'Z'),)),
    }
    _premise(arms)
    written, _ = _run(arms)
    assert written[0].get_tag('ZZ') == 'winner'


def test_invariant_raises_if_the_writer_reintroduces_the_splice():
    """The emit-time assertion must catch a two-candidate record, not just the
    one code path that produced it.  Re-create the pre-fix writer and confirm
    the record is refused instead of written."""
    from rectify.core.consensus import chimeric_consensus as cc

    # Patch the SOURCE module, not consensus.py: `_process_and_write_batch`
    # does `from .chimeric_consensus import ... build_chimeric_read` INSIDE the
    # function, so the name is rebound locally on every call and patching
    # consensus.build_chimeric_read is a silent no-op.
    real = cc.build_chimeric_read

    def _pre_fix(template_read, ref_start, cigar_tuples, chimeric_result,
                 header, anchor_read=None):
        # Re-create master's output: the winner's POS/CIGAR on a bystander's
        # contig. Overwriting reference_id AFTER the real build is the only
        # faithful way -- passing anchor_read=template_read would be a no-op,
        # because the fixed template selection already picks the anchor first.
        out = real(template_read, ref_start, cigar_tuples, chimeric_result,
                   header, anchor_read=anchor_read)
        out.reference_id = header.get_tid('chrA')   # the LOSER's contig
        return out

    arms = _cross_contig_arms()
    monkey = pytest.MonkeyPatch()
    monkey.setattr(cc, 'build_chimeric_read', _pre_fix)
    try:
        with pytest.raises(RuntimeError, match='placement invariant violated'):
            _run(arms)
    finally:
        monkey.undo()


def test_single_contig_single_strand_is_unchanged():
    """Control: when every arm agrees there is nothing to disagree about."""
    arms = {
        'loser_first': _read('chrA', 100, [(4, 10), (0, 40)]),
        'winner': _read('chrA', 400, [(0, 50)]),
    }
    written, stats = _run(arms)
    assert len(written) == 1
    assert written[0].reference_name == 'chrA'
    assert stats['placement_invariant_violations'] == 0


def test_single_aligner_is_unaffected():
    arms = {'only': _read('chrB', 400, [(0, 50)])}
    written, stats = _run(arms)
    assert len(written) == 1
    assert (written[0].reference_name, written[0].reference_start) == ('chrB', 400)
    assert written[0].mapping_quality == 60
    assert stats['placement_invariant_violations'] == 0


def test_contig_is_resolved_by_name_not_by_tid():
    """The anchor's `reference_id` indexes ITS OWN header. Copying the integer
    is only right while every arm's @SQ block matches the output's in order.

    In the runs that exposed this bug all seven arm BAMs and the output shared a
    byte-identical 17-entry @SQ block, so a tid copy happened to be correct —
    but a panel whose indices order contigs differently would silently relabel
    every record, which is the same failure with a much larger blast radius.
    """
    from rectify.core.consensus.chimeric_consensus import (
        build_chimeric_read, _fallback_simple_selection,
    )

    # An arm BAM whose @SQ order is the REVERSE of the output's.
    arm_header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'unsorted'},
        'SQ': [{'SN': 'chrB', 'LN': 1000}, {'SN': 'chrA', 'LN': 1000}],
    })
    anchor = pysam.AlignedSegment(arm_header)
    anchor.query_name = 'read1'
    anchor.query_sequence = SEQ
    anchor.flag = 0
    anchor.reference_id = arm_header.get_tid('chrB')   # tid 0 here, tid 1 in HEADER
    anchor.reference_start = 400
    anchor.mapping_quality = 60
    anchor.cigartuples = [(0, 50)]
    assert anchor.reference_id != HEADER.get_tid('chrB'), "premise: the tids must differ"

    result = _fallback_simple_selection({'winner': anchor}, GENOME, set())
    out = build_chimeric_read(
        template_read=anchor, ref_start=400, cigar_tuples=[(0, 50)],
        chimeric_result=result, header=HEADER, anchor_read=anchor,
    )
    assert out.reference_name == 'chrB', (
        f"tid was copied blind across headers: got {out.reference_name}"
    )
