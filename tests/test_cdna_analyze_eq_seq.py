"""cdna-analyze read calmd `=` SEQ as sequence and skipped the poly(A) walkback (Chanfreau planning/578).

Every `rectify align` arm is `samtools calmd -e` encoded, and the multialigned BAM kept that SEQ until
#6 (83af3d1). The walkback anchors at the first read base that equals the reference and is not the
stop base; a stored `=` never equals a genome letter, so the scan never anchored and
``corrected_3prime`` stayed the raw alignment end. On the 2026-10-09 route smoke that left 20,160 of
78,898 molecules (25.6 %) 1-50 nt downstream of where the walkback puts them. ``correct`` decodes at
intake; cdna-analyze now does too, and the shared walkback core treats `=` as the reference base, so
any other caller is safe.
"""

import pysam

from rectify.core.commands.cdna_analyze_command import _read_info_from_bam_record
from rectify.core.correct.walkback import (
    THREE_PRIME_SIDE_LEFT,
    THREE_PRIME_SIDE_RIGHT,
    walkback_3prime,
    walkback_3prime_guarded,
    walkback_3prime_with_qpos,
)

CHROM = "chrT"
GLEN = 2000
BODY_LEN = 100
PLUS_START, MINUS_START = 500, 1000
TAIL = 20
CLIP = "G" * 30

_genome = list(("CGT" * (GLEN // 3 + 1))[:GLEN])
# Plus-strand molecule: the aligned body ends C-AAAA, so the walkback must step back over the A's.
_genome[PLUS_START + BODY_LEN - 5:PLUS_START + BODY_LEN] = "CAAAA"
# Minus-strand molecule: the poly(A) side is reference_start; the body starts TTTT-G.
_genome[MINUS_START:MINUS_START + 5] = "TTTTG"
GENOME = "".join(_genome)
PLUS_BODY = GENOME[PLUS_START:PLUS_START + BODY_LEN]
MINUS_BODY = GENOME[MINUS_START:MINUS_START + BODY_LEN]
PLUS_ANCHOR = PLUS_START + BODY_LEN - 5      # the C before the genomic AAAA
MINUS_ANCHOR = MINUS_START + 4               # the G after the genomic TTTT


def _record(strand, encoded, pretrimmed=False):
    """The same molecule, with the aligned body literal or `calmd -e` encoded (all matches -> `=`).

    ``pretrimmed``: correct-cdna cut the 3' end before alignment (no tail in SEQ, XK:i:60) and carried
    the measured tail as XA:i:33."""
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]})
    r = pysam.AlignedSegment(header)
    r.query_name = f"eq_{strand}_{int(encoded)}"
    r.reference_id = 0
    r.mapping_quality = 60
    if strand == "+":
        body = "=" * BODY_LEN if encoded else PLUS_BODY
        seq, cigar, flag, start, xo = CLIP + body + "A" * TAIL, [(4, len(CLIP)), (0, BODY_LEN), (4, TAIL)], 0, PLUS_START, "fwd"
        if pretrimmed:
            seq, cigar = CLIP + body, [(4, len(CLIP)), (0, BODY_LEN)]
    else:
        body = "=" * BODY_LEN if encoded else MINUS_BODY
        seq, cigar, flag, start, xo = "T" * TAIL + body + CLIP, [(4, TAIL), (0, BODY_LEN), (4, len(CLIP))], 16, MINUS_START, "rev"
    r.flag = flag
    r.reference_start = start
    r.cigartuples = cigar
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    for tag, val, typ in [("XU", "A" * 27, "Z"), ("XO", xo, "Z"), ("XT", 1, "i"), ("XY", f"umi_captured_{xo}", "Z"),
                          ("XC", 3, "i"), ("XF", 2, "i"), ("XA", 33 if pretrimmed else TAIL, "i"), ("XN", 1, "i")]:
        r.set_tag(tag, val, typ)
    if pretrimmed:
        r.set_tag("XK", 60, "i")
    return r


def test_cdna_analyze_walks_back_on_calmd_encoded_seq():
    for strand, expected in (("+", PLUS_ANCHOR), ("-", MINUS_ANCHOR)):
        literal, _ = _read_info_from_bam_record(_record(strand, False), GENOME, {})
        encoded, _ = _read_info_from_bam_record(_record(strand, True), GENOME, {})
        assert literal.anchor == expected, (strand, literal.anchor)
        assert encoded.anchor == expected, (strand, encoded.anchor)   # was the raw alignment end
        assert encoded.tail_len == literal.tail_len
        assert encoded.pos5_corrected == literal.pos5_corrected


def test_cdna_analyze_decodes_the_record_it_tags():
    rec = _record("+", True)
    _read_info_from_bam_record(rec, GENOME, {})
    assert "=" not in rec.query_sequence and rec.query_sequence == _record("+", False).query_sequence


def test_shared_walkback_core_treats_equals_as_the_reference_base():
    for strand, side, stop in (("+", THREE_PRIME_SIDE_RIGHT, "A"), ("-", THREE_PRIME_SIDE_LEFT, "T")):
        for fn in (walkback_3prime, walkback_3prime_with_qpos, walkback_3prime_guarded):
            assert fn(_record(strand, True), GENOME, side, stop_base=stop) == \
                fn(_record(strand, False), GENOME, side, stop_base=stop), (strand, fn.__name__)
    assert walkback_3prime_with_qpos(_record("+", True), GENOME, THREE_PRIME_SIDE_RIGHT)[1] == PLUS_ANCHOR


def test_pretrimmed_record_keeps_the_carried_tail_length():
    # correct-cdna cut the tail (XK > 0) after measuring it (XA). The four A's left over the genomic AAAA are
    # not the tail; while SEQ was '='-encoded the walk counted nothing and XA won by accident, and with
    # literal bases (#6) the walk's residual count replaced it (smoke: median 33 -> 5 for 20,301 molecules).
    for encoded in (False, True):
        stats = {}
        info, _ = _read_info_from_bam_record(_record("+", encoded, pretrimmed=True), GENOME, stats)
        assert info.anchor == PLUS_ANCHOR
        assert info.tail_len == 33, (encoded, info.tail_len)
        assert stats == {"carried": 1}


def test_equals_left_inside_a_clip_does_not_stop_the_run():
    # A later CIGAR edit can leave `=` inside a clip, where no reference base defines it, and the strict decoder
    # refuses the record (r000_24648 on the c287748 smoke BAM). cdna-analyze analyses it as stored instead of
    # aborting the whole library; the walkback still resolves the aligned `=`.
    rec = _record("+", True)
    rec.query_sequence = rec.query_sequence[:-1] + "="
    rec.query_qualities = pysam.qualitystring_to_array("I" * len(rec.query_sequence))
    stats = {}
    info, _ = _read_info_from_bam_record(rec, GENOME, stats)
    assert info.anchor == PLUS_ANCHOR
    assert stats.get("eq_unresolved") == 1


def _clip_equals_record():
    """A plus-strand molecule whose last soft-clipped base is `=` (no reference base defines it)."""
    rec = _record("+", True)
    rec.query_sequence = rec.query_sequence[:-1] + "="
    rec.query_qualities = pysam.qualitystring_to_array("I" * len(rec.query_sequence))
    return rec


def test_decoder_refuses_by_default_and_writes_n_on_request():
    import pytest
    from rectify.core.bam.alignment_tags import decode_original_sequence
    with pytest.raises(ValueError, match="soft clip"):
        decode_original_sequence(_clip_equals_record(), {CHROM: GENOME})
    rec = _clip_equals_record()
    assert decode_original_sequence(rec, {CHROM: GENOME}, unresolved="N")
    assert rec.query_sequence == CLIP + PLUS_BODY + "A" * (TAIL - 1) + "N"


def test_correct_intake_survives_equals_in_a_clip():
    # master's `correct` aborted the whole run on such a record from a pre-#6 multialigned BAM (2026-10-10:
    # read f2ee18b7 on the c287748 DRS route-smoke BAM); it now corrects it with the clip base written N.
    from rectify.core.bam import bam_processor as bp
    rec = _clip_equals_record()
    bp.correct_read_3prime(rec, {CHROM: GENOME})
    assert "=" not in rec.query_sequence and rec.query_sequence.endswith("N")
