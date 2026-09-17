"""Chanfreau planning/936: cdna-analyze put the poly(A) on the wrong end of a molecule.

``XO`` is decided on the PRE-alignment record and describes the molecule in that
record's BAM-SEQ frame.  ``correct-cdna`` emits every consensus RNA-sense
(``XN:i:1``), so after alignment the flag alone fixes the frame: flag 0 = poly(A)
at the RIGHT ('fwd'), flag 16 = poly(T) at the LEFT ('rev').  When the consensus
re-aligns in the other frame from its pre-alignment record — reporter molecules
parked on chrXV instead of the construct contig, 5'-truncated molecules whose
clipped 5' part decided the pre-alignment strand — the carried XO is stale, the
walkback runs on the wrong end and ``corrected_3prime`` lands on the 5' boundary
(measured: XO disagreed with the flag on 43 % of reverse-aligned RPL20B records,
2026-09-14).  On an oriented record the flag is authoritative; a record without
``XN`` (pre-3457ecc output) keeps XO, the only frame label it has.
"""

import pysam
import pytest

from rectify.core.commands.cdna_analyze_command import _read_info_from_bam_record

CHROM = "chrT"
GLEN = 2000
# No A anywhere, so the poly(A) walkback has no genomic A-tract to wander into and
# the anchor is the raw alignment end on the poly(A) side.
GENOME = ("CGT" * (GLEN // 3 + 1))[:GLEN]
BODY_START, BODY_LEN = 500, 100
CLIP = "G" * 30           # the 5' part the reference lacks (non-A, unmatched)
TAIL = 20


def _record(flag, xo, seq, cigar, xn=1):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]})
    r = pysam.AlignedSegment(header)
    r.query_name = "cluster_936"
    r.flag = flag
    r.reference_id = 0
    r.reference_start = BODY_START
    r.mapping_quality = 60
    r.cigartuples = cigar
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    tags = [("XU", "A" * 27, "Z"), ("XO", xo, "Z"), ("XT", 2, "i"),
            ("XY", "umi_not_captured", "Z"), ("XC", 3, "i"), ("XF", 1, "i"),
            ("XA", TAIL, "i")]
    if xn is not None:
        tags.append(("XN", xn, "i"))
    for tag, val, typ in tags:
        r.set_tag(tag, val, typ)
    return r


BODY = GENOME[BODY_START:BODY_START + BODY_LEN]


def _minus_strand_molecule(xo, xn=1):
    """The reporter's recipe: an is_reverse record whose 5' part is a 30-nt
    non-A soft clip at reference_end and whose A-tail sits at reference_start
    (poly(T) in BAM-SEQ frame)."""
    seq = "T" * TAIL + BODY + CLIP
    return _record(16, xo, seq, [(4, TAIL), (0, BODY_LEN), (4, len(CLIP))], xn)


def _plus_strand_molecule(xo, xn=1):
    seq = CLIP + BODY + "A" * TAIL
    return _record(0, xo, seq, [(4, len(CLIP)), (0, BODY_LEN), (4, TAIL)], xn)


ALN_START = BODY_START
ALN_END_1 = BODY_START + BODY_LEN - 1


def test_a_reverse_aligned_oriented_molecule_is_rev_whatever_xo_says():
    stats = {}
    info, xc = _read_info_from_bam_record(_minus_strand_molecule(xo="fwd"), GENOME, stats)
    assert info.orient == "rev"
    assert info.is_reverse
    # the 3' end is the poly(A) side = reference_start, never the clipped 5' boundary
    assert abs(info.anchor - ALN_START) <= 2, info.anchor
    assert abs(info.anchor - ALN_END_1) > 50
    assert abs(info.pos5_corrected - ALN_END_1) <= 2
    assert stats["xo_frame_corrected"] == 1
    assert xc == 3


def test_a_forward_aligned_oriented_molecule_is_fwd_whatever_xo_says():
    stats = {}
    info, _ = _read_info_from_bam_record(_plus_strand_molecule(xo="rev"), GENOME, stats)
    assert info.orient == "fwd"
    assert abs(info.anchor - ALN_END_1) <= 2, info.anchor
    assert abs(info.pos5_corrected - ALN_START) <= 2
    assert stats["xo_frame_corrected"] == 1


@pytest.mark.parametrize("make,xo", [(_minus_strand_molecule, "rev"),
                                     (_plus_strand_molecule, "fwd")])
def test_an_agreeing_xo_is_not_counted_as_corrected(make, xo):
    stats = {}
    info, _ = _read_info_from_bam_record(make(xo=xo), GENOME, stats)
    assert info.orient == xo
    assert "xo_frame_corrected" not in stats


def test_a_record_without_xn_keeps_the_carried_xo():
    """Pre-3457ecc output has no orientation guarantee: XO is all it has."""
    stats = {}
    info, _ = _read_info_from_bam_record(_minus_strand_molecule(xo="fwd", xn=None),
                                         GENOME, stats)
    assert info.orient == "fwd"
    assert "xo_frame_corrected" not in stats


def test_the_frame_counter_reaches_the_run_log(tmp_path, caplog):
    """End to end through `run()`: the correction is reported, not silent."""
    import logging
    from rectify.core.commands import cdna_analyze_command

    ref = tmp_path / "ref.fa"
    ref.write_text(f">{CHROM}\n{GENOME}\n")
    pysam.faidx(str(ref))
    gff = tmp_path / "genes.gff"
    gff.write_text(f"{CHROM}\tsrc\tgene\t{BODY_START + 1}\t{BODY_START + BODY_LEN}\t.\t-\t.\tID=g936;Name=g936\n")
    bam = tmp_path / "cons.bam"
    rec = _minus_strand_molecule(xo="fwd")
    with pysam.AlignmentFile(str(bam), "wb", header=rec.header) as f:
        f.write(rec)
    pysam.index(str(bam))

    caplog.set_level(logging.INFO)
    import argparse
    parser = argparse.ArgumentParser()
    cdna_analyze_command.create_cdna_analyze_parser(parser.add_subparsers(dest="cmd"))
    args = parser.parse_args(["cdna-analyze", str(bam), "-o", str(tmp_path / "out"),
                              "--gff", str(gff), "--reference", str(ref)])
    cdna_analyze_command.run(args)
    assert any("XO frame taken from the alignment flag on 1" in m for m in caplog.messages), caplog.messages
    # and the per-molecule table carries the flag's strand and the poly(A)-side 3' end
    rows = [l.split("\t") for l in (tmp_path / "out" / "corrected_reads.tsv").read_text().splitlines()]
    hdr, body = rows[0], rows[1:]
    assert len(body) == 1
    row = dict(zip(hdr, body[0]))
    assert row["strand"] == "-"
    assert abs(int(row["corrected_3prime"]) - ALN_START) <= 2, row
