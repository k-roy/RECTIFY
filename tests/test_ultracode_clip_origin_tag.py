"""Real 2F clip-origin evidence must coexist with cDNA orientation metadata."""
import copy

import pysam
import pytest

from rectify.core.bam import bam_writer as bw
from rectify.core.bam.bam_processor import correct_read_3prime
from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row
from rectify.core.cdna.read_info import revcomp
from rectify.core.commands.cdna_analyze_command import _read_info_from_bam_record


def case(reverse, xn):
    base = "T" * 20 + "A" + "ACGTTGCATGCAGTCCATG" + "GT" + "N" * 96 + "AG" + "C" * 100
    genome = revcomp(base) if reverse else base
    h = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chrT", "LN": len(genome)}]})
    r = pysam.AlignedSegment(h)
    r.query_name, r.reference_id, r.mapping_quality = "clip_origin", 0, 60
    r.reference_start, r.flag = (40, 16) if reverse else (140, 0)
    seq = base[128:140] + "C" * 60
    r.query_sequence = revcomp(seq) if reverse else seq
    r.query_qualities = [31] * len(seq)
    r.cigarstring = "60M12S" if reverse else "12S60M"
    orient = "rev" if reverse else "fwd"
    for tag, value in {"XU": "ACG" * 9, "XO": orient, "XT": 1,
                       "XY": "umi_captured_" + orient, "XC": 1, "XF": 1, "XA": 20}.items():
        r.set_tag(tag, value)
    if xn:
        r.set_tag("XN", 1)
    intron = (100, 200) if reverse else (40, 140)
    return r, genome, {("chrT", *intron)}


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("xn", [False, True])
def test_actual_origin_survives_every_writer_without_overwriting_xo(tmp_path, reverse, xn):
    read, genome, annotation = case(reverse, xn)
    before_info, _ = _read_info_from_bam_record(read, genome)
    row = correct_read_3prime(copy.deepcopy(read), {"chrT": genome},
                             annotated_junctions=annotation, ont_cDNA=True)[0]
    assert row["five_prime_clip_origin"] == "intron"
    bam, tsv = tmp_path / "stock.bam", tmp_path / "corrected.tsv"
    with pysam.AlignmentFile(str(bam), "wb", header=read.header) as handle:
        handle.write(read)
    tsv.write_text("\t".join(CORRECTION_TSV_HEADER) + "\n" +
                   "\t".join(correction_result_to_tsv_row(row)) + "\n")
    paths = {arm: str(tmp_path / (arm + ".bam")) for arm in ("hard", "soft", "dual_hard", "dual_soft")}
    bw.write_corrected_bam(str(bam), str(tsv), paths["hard"], {"chrT": genome})
    bw.write_softclipped_bam(str(bam), str(tsv), paths["soft"], {"chrT": genome})
    bw.write_dual_bam(str(bam), str(tsv), paths["dual_hard"], paths["dual_soft"], {"chrT": genome})
    for path in paths.values():
        with pysam.AlignmentFile(path, "rb") as handle:
            output = next(handle)
        assert output.get_tag("XO") == read.get_tag("XO")
        assert output.get_tag("Xo") == "intron"
        assert output.query_sequence == read.query_sequence
        assert output.query_qualities == read.query_qualities
        assert output.cigarstring == read.cigarstring
        info, _ = _read_info_from_bam_record(output, genome)
        assert (info.orient, info.anchor, info.pos5_corrected) == (
            before_info.orient, before_info.anchor, before_info.pos5_corrected)


@pytest.mark.parametrize("reverse", [False, True])
def test_invalid_legacy_orientation_refuses_but_sense_frame_remains_usable(reverse):
    read, genome, _ = case(reverse, False)
    read.set_tag("XO", "intron")
    assert _read_info_from_bam_record(read, genome) is None
    read.set_tag("XN", 1)
    info, _ = _read_info_from_bam_record(read, genome)
    assert info.orient == ("rev" if reverse else "fwd")
