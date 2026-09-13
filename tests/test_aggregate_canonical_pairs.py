"""Canonical filtering must retain minor splice sites without inventing motif pairs."""
import pysam
import pytest

from rectify.core.aggregate.junctions import aggregate_junctions, filter_junctions


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("donor,acceptor,canonical", [
    ("GT", "AG", True), ("GC", "AG", True), ("AT", "AC", True),
    ("AT", "AG", False), ("GT", "AC", False), ("GC", "AC", False),
])
def test_canonical_motif_pairs_on_both_strands(tmp_path, strand, donor, acceptor, canonical):
    complement = str.maketrans("ACGT", "TGCA")
    if strand == "-":
        left, right = acceptor.translate(complement)[::-1], donor.translate(complement)[::-1]
    else:
        left, right = donor, acceptor
    genome = "C" * 20 + left + "C" * 96 + right + "C" * 20
    bam_path = tmp_path / "junction.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header={"SQ": [{"SN": "chrT", "LN": 140}]}) as bam:
        rec = pysam.AlignedSegment(bam.header)
        rec.query_name = "molecule"
        rec.reference_id = 0
        rec.reference_start = 0
        rec.cigarstring = "20M100N20M"
        rec.query_sequence = "C" * 40
        rec.flag = 16 if strand == "-" else 0
        bam.write(rec)
    table = aggregate_junctions(str(bam_path), genome={"chrT": genome})
    assert table.iloc[0]["five_ss_dinuc"] == donor
    assert table.iloc[0]["three_ss_dinuc"] == acceptor
    assert bool(table.iloc[0]["is_canonical"]) == canonical
    filtered = filter_junctions(table, min_reads=1, require_canonical=True)
    assert len(filtered) == int(canonical)
