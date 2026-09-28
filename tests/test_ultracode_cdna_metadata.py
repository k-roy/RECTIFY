"""RNA-sense frame survives correction and standalone sibling restoration."""
import csv
import copy

import pysam
import pytest

from rectify.core.bam.bam_processor import correct_read_3prime
from rectify.core.bam.bam_writer import (
    _load_corrections_from_single_tsv, apply_corrected_edits_to_read,
)
from rectify.core.commands.cdna_analyze_command import _read_info_from_bam_record
from rectify.core.consensus.consensus import _restore_comment_tags_from_siblings
from rectify.core.correct.protocols.ont_cdna import resolve_rna_strand


GENOME = "C" * 1000


def record(reverse, xo, xn=1, encoded=False):
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chrT", "LN": 1000}]})
    read = pysam.AlignedSegment(header)
    read.query_name = "frame_control"
    read.flag = 16 if reverse else 0
    read.reference_id = 0
    read.reference_start = 200
    read.cigarstring = "100M"
    read.mapping_quality = 60
    read.query_sequence = "=" * 100 if encoded else GENOME[200:300]
    read.query_qualities = list(range(20, 40)) * 5
    for tag, value in {"XU": "ACG" * 9, "XO": xo, "XT": 2,
                       "XY": "umi_not_captured", "XC": 3, "XF": 1,
                       "XA": 20, "XQ": 53, "XK": 37}.items():
        read.set_tag(tag, value)
    if xn is not None:
        read.set_tag("XN", xn)
    return read


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("xo", ["fwd", "rev"])
@pytest.mark.parametrize("xn", [None, 1])
@pytest.mark.parametrize("encoded", [False, True])
def test_correct_tsv_writer_agrees_with_current_frame(tmp_path, reverse, xo, xn, encoded):
    read = record(reverse, xo, xn, encoded)
    orient = ("rev" if reverse else "fwd") if xn == 1 else xo
    strand = "+" if orient == "fwd" else "-"
    end = 299 if strand == "+" else 200
    qualities = read.query_qualities[:]
    result = correct_read_3prime(read, {"chrT": GENOME}, ont_cDNA=True,
                                apply_3ss_rescue=False, apply_atract=False)[0]
    assert result["strand"] == strand
    assert result["original_3prime"] == result["corrected_3prime"] == end
    columns = ["read_id", "corrected_3prime", "strand", "five_prime_position",
               "five_prime_rescued", "five_prime_clip_origin"]
    path = tmp_path / "correction.tsv"
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerow(result)
    correction = _load_corrections_from_single_tsv(str(path))[read.query_name]
    apply_corrected_edits_to_read(read, correction, {"chrT": GENOME})
    info, _ = _read_info_from_bam_record(read, GENOME)
    assert info.orient == orient and info.anchor == end
    assert (read.reference_start, read.cigarstring, read.query_sequence) == (200, "100M", "C" * 100)
    assert read.query_qualities == qualities
    assert read.get_tag("XO") == xo


@pytest.mark.parametrize("xn", [None, 0, 2, "invalid"])
def test_non_sense_marker_keeps_legacy_orientation(xn):
    read = record(True, "fwd", xn)
    assert resolve_rna_strand(read)[0] == "+"


@pytest.mark.parametrize("reverse", [False, True])
def test_sense_frame_precedes_missing_or_stale_orientation(reverse):
    read = record(reverse, "invalid")
    read.set_tag("ro", "A", "A")
    expected = "-" if reverse else "+"
    assert resolve_rna_strand(read)[0] == expected
    read.set_tag("XO", None)
    assert resolve_rna_strand(read)[0] == expected


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("xo", ["fwd", "rev"])
def test_sibling_restores_frame_with_typed_comment_block(reverse, xo):
    source = record(reverse, xo)
    winner = copy.deepcopy(source)
    winner.set_tags([("XC", "NO_SPLICE"), ("XA", "ultra_splice")])
    before = winner.query_sequence, winner.query_qualities, winner.cigarstring, winner.reference_start
    _restore_comment_tags_from_siblings(winner, {"minimap2": source, "uLTRA": winner})
    for tag in ["XN", "XO", "XQ", "XK", "XC", "XA"]:
        assert winner.get_tag(tag, with_value_type=True) == source.get_tag(tag, with_value_type=True)
    assert before == (winner.query_sequence, winner.query_qualities,
                      winner.cigarstring, winner.reference_start)
    info, _ = _read_info_from_bam_record(winner, GENOME)
    assert info.anchor == (200 if reverse else 299)


def test_authoritative_winner_and_missing_source_are_unchanged():
    source = record(True, "fwd")
    winner = record(False, "rev", xn=None)
    before = winner.to_string()
    _restore_comment_tags_from_siblings(winner, {"minimap2": source, "uLTRA": winner})
    assert winner.to_string() == before
    winner.set_tags([])
    before = winner.to_string()
    _restore_comment_tags_from_siblings(winner, {"uLTRA": winner})
    assert winner.to_string() == before
