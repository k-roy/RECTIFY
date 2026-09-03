"""Module 2F must read the TRANSCRIPT 5' clip, which under an antisense protocol is not the
``read.is_reverse`` clip.

``rescue_3ss_truncation`` computes ``align_5prime`` from the GENE strand
(``reference_start`` on '+', ``reference_end - 1`` on '-') but used to take the clip from
``read.is_reverse``. For a sense protocol the two coincide; for ``--netseq`` / ``--dT-primed-cDNA``
they are OPPOSITE ends of the read, so the module compared the RNA-3'-end clip (which carries the
randomer and the poly(A) tail) against exon-1 sequence and could write a false ``five_prime_position``.
Verified by execution in planning 834 §7.1; this locks the fix.
"""
from __future__ import annotations

import pysam

from rectify.core.splice.splice_aware_5prime import (
    _extract_5prime_rescue_seq,
    _get_5prime_softclip_len,
    _transcript_5prime_is_right,
)

_M, _S = 0, 4
_HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrI", "LN": 10_000}]}
)


def _read(cigartuples, is_reverse, seq, reference_start=100):
    read = pysam.AlignedSegment(_HEADER)
    read.query_name = "r"
    read.reference_id = 0
    read.reference_start = reference_start
    read.mapping_quality = 60
    read.flag = 16 if is_reverse else 0
    read.cigartuples = list(cigartuples)
    read.query_sequence = seq
    read.query_qualities = [30] * len(seq)
    return read


def test_transcript_5prime_side_follows_the_gene_strand_when_given():
    fwd = _read([(_M, 30)], False, "C" * 30)
    rev = _read([(_M, 30)], True, "C" * 30)
    # no strand -> the pre-fix behaviour, byte-identical for every sense caller
    assert _transcript_5prime_is_right(fwd) is False
    assert _transcript_5prime_is_right(rev) is True
    # with a gene strand -> the gene strand decides, whatever the BAM flag says
    assert _transcript_5prime_is_right(rev, '+') is False
    assert _transcript_5prime_is_right(fwd, '-') is True


def test_antisense_read_reports_the_left_clip_as_the_transcript_5prime_clip():
    """--netseq: gene '+' <=> is_reverse True, so the transcript 5' end is the LEFT end."""
    left_clipped = _read([(_S, 5), (_M, 30)], True, "ACGTA" + "C" * 30)
    assert _get_5prime_softclip_len(left_clipped) == 0            # the bug: clip missed entirely
    assert _get_5prime_softclip_len(left_clipped, '+') == 5       # the fix


def test_antisense_read_does_not_mistake_the_rna_3prime_clip_for_the_5prime_clip():
    """The RIGHT clip on an antisense read holds the randomer / poly(A) tail, not exon 1."""
    right_clipped = _read([(_M, 30), (_S, 5)], True, "C" * 30 + "AAAAA")
    assert _get_5prime_softclip_len(right_clipped) == 5           # the bug: tail read as exon 1
    assert _get_5prime_softclip_len(right_clipped, '+') == 0      # the fix


def test_sense_protocol_is_unchanged():
    sense = _read([(_S, 5), (_M, 30)], False, "ACGTA" + "C" * 30)
    assert _get_5prime_softclip_len(sense) == 5
    assert _get_5prime_softclip_len(sense, '+') == 5
    sense_rev = _read([(_M, 30), (_S, 5)], True, "C" * 30 + "ACGTA")
    assert _get_5prime_softclip_len(sense_rev) == 5
    assert _get_5prime_softclip_len(sense_rev, '-') == 5


def test_rescue_seq_is_extracted_from_the_gene_strand_end():
    genome = "G" * 10_000
    left_clipped = _read([(_S, 5), (_M, 30)], True, "ACGTA" + "G" * 30)
    # is_reverse scan starts at the WRONG end: it walks the whole 30-nt M block before reaching the
    # clip, so the "rescue sequence" is all 35 query bases (planning 834 §7.1 measured exactly this)
    assert _extract_5prime_rescue_seq(left_clipped, genome_seq=genome) == "ACGTA" + "G" * 30
    # gene-strand '+' scan starts at the transcript 5' end and recovers just the soft clip
    assert _extract_5prime_rescue_seq(left_clipped, genome_seq=genome, strand='+') == "ACGTA"
