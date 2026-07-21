"""Tests for UMI-aware deduplication: rectify/core/umi/dedup.py.

Two layers:
  * ``select_molecules`` -- the pure clustering core, tested with literal
    FragmentRecords (no BAM). This is where the molecule-identity logic lives.
  * ``dedup_bam`` -- the end-to-end BAM path on tiny synthetic BAMs, verifying
    PCR duplicates collapse, both mates of a survivor are kept, and MI/cD tags land.
"""
from __future__ import annotations

import pysam
import pytest

from rectify.core.umi.dedup import (
    FragmentRecord,
    UmiDedupStats,
    dedup_bam,
    extract_fragments,
    select_molecules,
)

_M, _N, _S = 0, 3, 4


def _frag(qname, umi, five_prime=100, span=300, strand="+", spliced=False, score=50,
          contig="chrI"):
    return FragmentRecord(qname=qname, contig=contig, strand=strand,
                          five_prime=five_prime, span=span, spliced=spliced,
                          umi=umi, score=score)


# ---------------------------------------------------------------------------
# select_molecules -- molecule identity
# ---------------------------------------------------------------------------

def test_pcr_duplicates_collapse_to_one_molecule():
    """Same span + same UMI = same molecule (the core dedup case)."""
    frags = [_frag("a", "ACGTACGTACGT"), _frag("b", "ACGTACGTACGT"),
             _frag("c", "ACGTACGTACGT")]
    keepers, family = select_molecules(frags, edit_distance=1)
    assert len(keepers) == 1
    assert list(family.values()) == [3]


def test_distinct_umis_same_position_are_distinct_molecules():
    frags = [_frag("a", "ACGTACGTACGT"), _frag("b", "TTTTGGGGCCCC")]
    keepers, _ = select_molecules(frags, edit_distance=1)
    assert len(keepers) == 2


def test_error_umi_absorbed_at_ed1():
    """A 1-substitution UMI with lower count is a sequencing error -> absorbed."""
    frags = [_frag(f"hi{i}", "ACGTACGTACGT") for i in range(5)]
    frags.append(_frag("err", "ACGTACGTACGA"))  # 1 sub, count 1
    keepers, family = select_molecules(frags, edit_distance=1)
    assert len(keepers) == 1
    assert list(family.values()) == [6]


def test_equal_count_neighbours_stay_separate():
    """Two real molecules (no count gradient) are NOT merged even at ED1."""
    frags = ([_frag(f"a{i}", "ACGTACGTACGT") for i in range(3)] +
             [_frag(f"b{i}", "ACGTACGTACGA") for i in range(3)])
    keepers, _ = select_molecules(frags, edit_distance=1)
    assert len(keepers) == 2


def test_different_span_separates_same_umi():
    """Random priming: same 5' + same UMI but different fragment length = 2 molecules."""
    frags = [_frag("a", "ACGTACGTACGT", span=300),
             _frag("b", "ACGTACGTACGT", span=500)]
    keepers, _ = select_molecules(frags, edit_distance=1)
    assert len(keepers) == 2


def test_spliced_and_unspliced_are_distinct():
    frags = [_frag("a", "ACGTACGTACGT", spliced=False),
             _frag("b", "ACGTACGTACGT", spliced=True)]
    keepers, _ = select_molecules(frags, edit_distance=1)
    assert len(keepers) == 2


def test_strand_separates_same_umi_and_position():
    frags = [_frag("a", "ACGTACGTACGT", strand="+"),
             _frag("b", "ACGTACGTACGT", strand="-")]
    keepers, _ = select_molecules(frags, edit_distance=1)
    assert len(keepers) == 2


def test_representative_is_highest_score():
    frags = [_frag("low", "ACGTACGTACGT", score=10),
             _frag("high", "ACGTACGTACGT", score=90)]
    keepers, _ = select_molecules(frags, edit_distance=1)
    assert "high" in keepers and "low" not in keepers


def test_stats_family_size_and_edit_distributions():
    frags = [_frag(f"hi{i}", "ACGTACGTACGT") for i in range(4)]
    frags.append(_frag("err", "ACGTACGTACGA"))  # 1 substitution
    stats = UmiDedupStats()
    select_molecules(frags, edit_distance=1, stats=stats)
    assert stats.n_molecules == 1
    assert stats.family_size_hist == {5: 1}
    # 4 members match canonical exactly (d0), 1 differs by 1 (d1)
    assert stats.within_family_edit_hist == {0: 4, 1: 1}


def test_duplication_rate():
    frags = [_frag("a", "ACGTACGTACGT"), _frag("b", "ACGTACGTACGT"),
             _frag("c", "TTTTGGGGCCCC")]
    stats = UmiDedupStats()
    stats.n_input_fragments = 3
    select_molecules(frags, edit_distance=1, stats=stats)
    assert stats.n_molecules == 2
    assert stats.duplication_rate == pytest.approx(1 - 2 / 3)


def test_components_clustering_mode_available():
    frags = [_frag("a", "ACGTACGTACGT"), _frag("b", "ACGTACGTACGT")]
    keepers, _ = select_molecules(frags, edit_distance=1, clustering="components")
    assert len(keepers) == 1


# ---------------------------------------------------------------------------
# dedup_bam -- end to end
# ---------------------------------------------------------------------------

def _header():
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrI", "LN": 5000}],
    })


def _pair(header, qname, umi, r1_pos=100, tlen=300, spliced=False):
    """Build an R1/R2 pair carrying RX on both, R1 forward at r1_pos."""
    r1 = pysam.AlignedSegment(header)
    r1.query_name = qname
    r1.flag = 0x1 | 0x2 | 0x40  # paired, proper, read1
    r1.reference_id = 0
    r1.reference_start = r1_pos
    r1.mapping_quality = 60
    r1.cigarstring = "20M100N30M" if spliced else "50M"
    r1.query_sequence = "A" * 50
    r1.query_qualities = pysam.qualitystring_to_array("I" * 50)
    r1.next_reference_id = 0
    r1.next_reference_start = r1_pos + tlen - 50
    r1.template_length = tlen
    r1.set_tag("RX", umi, value_type="Z")

    r2 = pysam.AlignedSegment(header)
    r2.query_name = qname
    r2.flag = 0x1 | 0x2 | 0x80 | 0x10  # paired, proper, read2, reverse
    r2.reference_id = 0
    r2.reference_start = r1_pos + tlen - 50
    r2.mapping_quality = 60
    r2.cigarstring = "50M"
    r2.query_sequence = "T" * 50
    r2.query_qualities = pysam.qualitystring_to_array("I" * 50)
    r2.next_reference_id = 0
    r2.next_reference_start = r1_pos
    r2.template_length = -tlen
    r2.set_tag("RX", umi, value_type="Z")
    return [r1, r2]


def _write_sorted_bam(path, pairs, header):
    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
        for reads in pairs:
            for r in reads:
                bam.write(r)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
    return path


def _read_all(path):
    with pysam.AlignmentFile(str(path), "rb") as bam:
        return list(bam)


def test_dedup_bam_collapses_pcr_duplicates_keeps_both_mates(tmp_path):
    hdr = _header()
    pairs = [
        _pair(hdr, "dup1", "ACGTACGTACGT"),   # 3 PCR duplicates of one molecule
        _pair(hdr, "dup2", "ACGTACGTACGT"),
        _pair(hdr, "dup3", "ACGTACGTACGT"),
        _pair(hdr, "other", "TTTTGGGGCCCC"),  # a distinct molecule, same position
    ]
    inp = _write_sorted_bam(tmp_path / "in.bam", pairs, hdr)
    stats = dedup_bam(str(inp), str(tmp_path / "out.bam"))

    assert stats.n_input_fragments == 4      # 4 R1 reads
    assert stats.n_molecules == 2            # 3 dups -> 1, plus 'other'
    out = _read_all(tmp_path / "out.bam")
    kept_names = {r.query_name for r in out}
    assert len(kept_names) == 2              # one survivor per molecule
    # both mates of each survivor are written
    assert len(out) == 4
    for r in out:
        assert r.has_tag("MI") and r.has_tag("cD")


def test_dedup_bam_cd_tag_is_family_size(tmp_path):
    hdr = _header()
    pairs = [_pair(hdr, f"d{i}", "ACGTACGTACGT") for i in range(5)]
    inp = _write_sorted_bam(tmp_path / "in.bam", pairs, hdr)
    dedup_bam(str(inp), str(tmp_path / "out.bam"))
    out = _read_all(tmp_path / "out.bam")
    assert {r.get_tag("cD") for r in out} == {5}


def test_dedup_bam_distinct_spans_not_collapsed(tmp_path):
    """Same UMI, same 5', different fragment length -> two molecules survive."""
    hdr = _header()
    pairs = [_pair(hdr, "short", "ACGTACGTACGT", tlen=300),
             _pair(hdr, "long", "ACGTACGTACGT", tlen=800)]
    inp = _write_sorted_bam(tmp_path / "in.bam", pairs, hdr)
    stats = dedup_bam(str(inp), str(tmp_path / "out.bam"))
    assert stats.n_molecules == 2


def test_dedup_bam_reads_without_umi_counted_and_skipped(tmp_path):
    hdr = _header()
    good = _pair(hdr, "good", "ACGTACGTACGT")
    noumi = _pair(hdr, "noumi", "ACGTACGTACGT")
    for r in noumi:
        # remove the RX tag
        r.set_tags([t for t in r.get_tags() if t[0] != "RX"])
    inp = _write_sorted_bam(tmp_path / "in.bam", [good, noumi], hdr)
    stats = dedup_bam(str(inp), str(tmp_path / "out.bam"))
    assert stats.n_no_umi == 1               # the R1 of the untagged pair
    assert stats.n_molecules == 1


def test_dedup_bam_consensus_mode_not_yet_implemented(tmp_path):
    hdr = _header()
    inp = _write_sorted_bam(tmp_path / "in.bam", [_pair(hdr, "a", "ACGTACGTACGT")], hdr)
    with pytest.raises(NotImplementedError):
        dedup_bam(str(inp), str(tmp_path / "out.bam"), collapse="consensus")


def test_extract_fragments_strand_and_span(tmp_path):
    hdr = _header()
    inp = _write_sorted_bam(tmp_path / "in.bam", [_pair(hdr, "a", "ACGTACGTACGT", tlen=300)], hdr)
    with pysam.AlignmentFile(str(inp), "rb") as bam:
        frags = extract_fragments(bam, r1_sense=True)
    assert len(frags) == 1              # only R1 yields a fragment
    assert frags[0].span == 300
    assert frags[0].strand == "+"       # R1 forward, CORALL sense
    assert frags[0].umi == "ACGTACGTACGT"


def test_extract_fragments_antisense_flips_strand(tmp_path):
    hdr = _header()
    inp = _write_sorted_bam(tmp_path / "in.bam", [_pair(hdr, "a", "ACGTACGTACGT")], hdr)
    with pysam.AlignmentFile(str(inp), "rb") as bam:
        frags = extract_fragments(bam, r1_sense=False)
    assert frags[0].strand == "-"       # R1 forward but antisense protocol -> minus
