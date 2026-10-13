"""Tests for the UMI consensus extension that uses the random-primer bases (rectify.core.umi.consensus).

Every rule has a test that would fail if the rule broke: the column call, the primer span on both strands and both
sources, the mate 5′ end from MC, the primer-based split of a UMI cluster (and the cases where it must NOT split),
masking on whichever mate covers the primer, the mismatch classes, and the span key that survives soft-clipped primers.
"""
from __future__ import annotations

import random

import pysam
import pytest

from rectify.core.umi.consensus import (
    call_column,
    consensus_bam,
    mate_five_prime_unclipped,
    primer_bases_read_orientation,
    primer_span,
    split_by_primer,
)

COMP = str.maketrans("ACGT", "TGCA")


# ------------------------------------------------------------------------------------------ primitives
def test_call_column_unanimous_high_quality():
    b, q, post = call_column("AAA", [30, 30, 30])
    assert b == "A" and q >= 30 and post > 0.999


def test_call_column_outvotes_a_single_error():
    b, q, _ = call_column("AAC", [30, 30, 30])
    assert b == "A"
    assert q < 40          # the dissenting base lowers confidence


def test_call_column_tie_becomes_n():
    b, q, post = call_column("AC", [30, 30])
    assert b == "N" and q == 2 and post < 0.9


def test_call_column_ignores_n_and_empty():
    assert call_column("NN", [30, 30])[0] == "N"
    assert call_column("", [])[0] == "N"
    assert call_column("NA", [30, 30])[0] == "A"


def test_primer_span_both_strands_and_sources():
    assert primer_span(100, False, 9, "seq") == (100, 108)
    assert primer_span(100, True, 9, "seq") == (92, 100)
    assert primer_span(100, False, 9, "tag") == (91, 99)
    assert primer_span(100, True, 9, "tag") == (101, 109)
    with pytest.raises(ValueError):
        primer_span(100, True, 9, "other")


def _hdr(L=3000):
    return pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrI", "LN": L}]})


def _read(hdr, name, flag, start, cigar, seq, mate_start=None, mate_cigar=None, tlen=0, umi=None):
    r = pysam.AlignedSegment(hdr)
    r.query_name = name
    r.flag = flag
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigarstring = cigar
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    if mate_start is not None:
        r.next_reference_id = 0
        r.next_reference_start = mate_start
        r.template_length = tlen
        r.set_tag("MC", mate_cigar, value_type="Z")
    if umi:
        r.set_tag("RX", umi, value_type="Z")
    return r


def test_primer_bases_read_orientation_forward_reverse_tag():
    h = _hdr()
    fwd = _read(h, "f", 0x1 | 0x80, 100, "12M", "ACGTACGTACGT")
    assert primer_bases_read_orientation(fwd, 4, "seq")[0] == "ACGT"
    rev = _read(h, "r", 0x1 | 0x80 | 0x10, 100, "12M", "AAAACCCCGGTT")
    # read orientation = reverse complement of SEQ: first 4 bases = revcomp of the last 4 of SEQ ("GGTT" -> "AACC")
    assert primer_bases_read_orientation(rev, 4, "seq")[0] == "AACC"
    tag = _read(h, "t", 0x1 | 0x80, 100, "12M", "ACGTACGTACGT")
    tag.set_tag("pR", "TTTTTTTTT")
    tag.set_tag("pQ", "IIIIIIIII")
    s, q = primer_bases_read_orientation(tag, 9, "tag")
    assert s == "TTTTTTTTT" and q == [40] * 9


def test_mate_five_prime_from_mc():
    h = _hdr()
    r1 = _read(h, "x", 0x1 | 0x40 | 0x20, 100, "50M", "A" * 50, mate_start=300, mate_cigar="41M9S", tlen=241)
    # reverse mate: 5′ = start + aligned length - 1 + trailing clip = 300 + 41 - 1 + 9
    assert mate_five_prime_unclipped(r1) == 349
    r1b = _read(h, "y", 0x1 | 0x40, 400, "50M", "A" * 50, mate_start=300, mate_cigar="5S45M", tlen=-145)
    assert mate_five_prime_unclipped(r1b) == 295   # forward mate: start - leading clip
    r1c = _read(h, "z", 0x1 | 0x40, 400, "50M", "A" * 50)
    assert mate_five_prime_unclipped(r1c) is None


def test_split_by_primer_two_raw_umis_two_primers():
    umis = ["AAAA", "AAAA", "AAAT", "AAAT"]
    prim = [("ACGTACGTA", [30] * 9)] * 2 + [("ACGTTCGTA", [30] * 9)] * 2
    groups, hetero = split_by_primer(umis, prim, min_support=2)
    assert sorted(map(sorted, groups)) == [[0, 1], [2, 3]] and not hetero


def test_split_by_primer_same_primer_does_not_split():
    umis = ["AAAA", "AAAA", "AAAT", "AAAT"]
    prim = [("ACGTACGTA", [30] * 9)] * 4
    groups, _ = split_by_primer(umis, prim, min_support=2)
    assert groups == [[0, 1, 2, 3]]


def test_split_by_primer_difference_within_one_raw_umi_is_flagged_not_split():
    umis = ["AAAA"] * 4
    prim = [("ACGTACGTA", [30] * 9)] * 2 + [("ACGTTCGTA", [30] * 9)] * 2
    groups, hetero = split_by_primer(umis, prim, min_support=2)
    assert groups == [[0, 1, 2, 3]] and hetero


def test_split_by_primer_low_quality_difference_does_not_split():
    umis = ["AAAA", "AAAA", "AAAT", "AAAT"]
    lowq = [30] * 4 + [5] + [30] * 4
    prim = [("ACGTACGTA", [30] * 9)] * 2 + [("ACGTTCGTA", lowq)] * 2
    groups, _ = split_by_primer(umis, prim, min_support=2)
    assert groups == [[0, 1, 2, 3]]


def test_split_by_primer_unsupported_raw_umi_joins_best_match():
    umis = ["AAAA", "AAAA", "AAAT", "AAAT", "AAAC"]
    prim = [("ACGTACGTA", [30] * 9)] * 2 + [("ACGTTCGTA", [30] * 9)] * 2 + [("ACGTTCGTA", [30] * 9)]
    groups, _ = split_by_primer(umis, prim, min_support=2)
    assert sorted(map(sorted, groups)) == [[0, 1], [2, 3, 4]]


# ------------------------------------------------------------------------------------------ end to end
@pytest.fixture()
def ref_fasta(tmp_path):
    rnd = random.Random(1136)
    seq = "".join(rnd.choice("ACGT") for _ in range(3000))
    p = tmp_path / "ref.fa"
    p.write_text(">chrI\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    pysam.faidx(str(p))
    return str(p), seq


def _alt(b):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[b]


def _pair(h, ref, name, umi, r1_start, tlen, rlen=50, mol_mut=None, r1_err=None, r2_err=None, r2_cigar=None):
    """R1 forward at r1_start, R2 reverse ending at r1_start + tlen - 1. mol_mut: {ref pos: base} on every read
    covering it (a molecule-level change); r1_err / r2_err: {ref pos: base} on that read only."""
    mol_mut, r1_err, r2_err = mol_mut or {}, r1_err or {}, r2_err or {}
    r2_start = r1_start + tlen - rlen

    def seq_for(start, errs):
        return "".join(errs.get(x, mol_mut.get(x, ref[x])) for x in range(start, start + rlen))
    r2_cig = r2_cigar or f"{rlen}M"
    r1 = _read(h, name, 0x1 | 0x2 | 0x40 | 0x20, r1_start, f"{rlen}M", seq_for(r1_start, r1_err),
               mate_start=r2_start, mate_cigar=r2_cig, tlen=tlen, umi=umi)
    r2 = _read(h, name, 0x1 | 0x2 | 0x80 | 0x10, r2_start, r2_cig, seq_for(r2_start, r2_err),
               mate_start=r1_start, mate_cigar=f"{rlen}M", tlen=-tlen, umi=umi)
    if r2_cigar and r2_cigar.endswith("S"):
        clip = int(r2_cigar.rstrip("S").split("M")[-1])
        r1.template_length = tlen - clip
        r2.template_length = -(tlen - clip)
    return [r1, r2]


def _write(path, pairs, h):
    tmp = str(path) + ".u.bam"
    with pysam.AlignmentFile(tmp, "wb", header=h) as out:
        for pr in pairs:
            for r in pr:
                out.write(r)
    pysam.sort("-o", str(path), tmp)
    pysam.index(str(path))
    return str(path)


def _reads(path):
    with pysam.AlignmentFile(path, "rb") as bam:
        return list(bam)


def test_consensus_corrects_minority_keeps_molecule_change_masks_primer(tmp_path, ref_fasta):
    fa, ref = ref_fasta
    h = _hdr()
    R1S, TLEN = 100, 300
    r2_end = R1S + TLEN - 1                         # R2 reverse: its 5′ end, primer span = r2_end-8 .. r2_end
    prim_pos = r2_end - 3
    mol = {130: _alt(ref[130]), prim_pos: _alt(ref[prim_pos])}
    pairs = [_pair(h, ref, f"m{i}", "ACGTACGTACGT", R1S, TLEN, mol_mut=mol,
                   r1_err={120: _alt(ref[120])} if i == 0 else None) for i in range(3)]
    pairs.append(_pair(h, ref, "solo", "TTTTGGGGCCCC", R1S, TLEN))
    inp = _write(tmp_path / "in.bam", pairs, h)
    out = str(tmp_path / "out.bam")
    st = consensus_bam(inp, out, reference=fa, primer_length=9, primer_mate="r2")
    assert st.n_molecules == 2 and st.n_input_fragments == 4
    recs = _reads(out)
    assert len(recs) == 4
    fam = [r for r in recs if r.get_tag("cD") == 3]
    r1 = next(r for r in fam if r.is_read1)
    r2 = next(r for r in fam if r.is_read2)
    assert r1.query_sequence[120 - R1S] == ref[120]          # sequencing error outvoted
    assert r1.query_sequence[130 - R1S] == mol[130]          # molecule-level change kept
    r2_start = R1S + TLEN - 50
    for x in range(r2_end - 8, r2_end + 1):                    # the primer span is masked on R2
        assert r2.query_sequence[x - r2_start] == "N"
    assert r2.get_tag("pK") == 9 and r1.get_tag("pK") == 0     # R1 does not reach the primer here
    ev = st.as_dict()["mismatch_events"]
    assert ev["minority"]["3"]["pass"] >= 1
    assert ev["molecule"]["3"]["pass"] >= 1
    assert ev["primer_template"]["3"]["pass"] >= 1


def test_short_insert_masks_primer_on_both_mates(tmp_path, ref_fasta):
    fa, ref = ref_fasta
    h = _hdr()
    R1S, TLEN = 500, 40
    pairs = [_pair(h, ref, f"s{i}", "GGGGAAAACCCC", R1S, TLEN, rlen=40) for i in range(2)]
    inp = _write(tmp_path / "in.bam", pairs, h)
    out = str(tmp_path / "out.bam")
    consensus_bam(inp, out, reference=fa, primer_length=9, primer_mate="r2")
    recs = _reads(out)
    assert {r.get_tag("pK") for r in recs} == {9}             # R1's 3′ end covers the same primer bases
    r1 = next(r for r in recs if r.is_read1)
    assert r1.query_sequence[-9:] == "N" * 9


def test_primer_split_separates_two_molecules_merged_by_umi_clustering(tmp_path, ref_fasta):
    fa, ref = ref_fasta
    h = _hdr()
    R1S, TLEN = 1000, 250
    r2_end = R1S + TLEN - 1
    p_a = {r2_end - 2: _alt(ref[r2_end - 2])}
    p_b = {r2_end - 5: _alt(ref[r2_end - 5])}
    pairs = [_pair(h, ref, f"a{i}", "ACGTACGTACGA", R1S, TLEN, mol_mut=p_a) for i in range(4)]
    pairs += [_pair(h, ref, f"b{i}", "ACGTACGTACGT", R1S, TLEN, mol_mut=p_b) for i in range(2)]   # 1 nt from the A UMI
    inp = _write(tmp_path / "in.bam", pairs, h)
    out = str(tmp_path / "out.bam")
    st = consensus_bam(inp, out, reference=fa, primer_length=9, primer_mate="r2")
    assert st.clusters_checked_for_primer == 1 and st.clusters_split_by_primer == 1
    assert st.n_molecules == 2
    assert sorted({r.get_tag("cD") for r in _reads(out)}) == [2, 4]


def test_no_split_when_primers_agree(tmp_path, ref_fasta):
    fa, ref = ref_fasta
    h = _hdr()
    pairs = [_pair(h, ref, f"a{i}", "ACGTACGTACGA", 1000, 250) for i in range(4)]
    pairs += [_pair(h, ref, f"b{i}", "ACGTACGTACGT", 1000, 250) for i in range(2)]
    inp = _write(tmp_path / "in.bam", pairs, h)
    st = consensus_bam(inp, str(tmp_path / "out.bam"), reference=fa, primer_length=9, primer_mate="r2")
    assert st.n_molecules == 1 and st.clusters_split_by_primer == 0


def test_unclipped_key_survives_soft_clipped_primer(tmp_path, ref_fasta):
    """Two copies of one molecule, one with its primer bases soft-clipped by the aligner: TLEN differs, the
    soft-clip-corrected ends do not."""
    fa, ref = ref_fasta
    h = _hdr()
    pairs = [_pair(h, ref, "c1", "CCCCAAAATTTT", 2000, 300),
             _pair(h, ref, "c2", "CCCCAAAATTTT", 2000, 300, r2_cigar="41M9S")]
    inp = _write(tmp_path / "in.bam", pairs, h)
    st = consensus_bam(inp, str(tmp_path / "out.bam"), reference=fa, primer_length=9, span_key="unclipped")
    assert st.n_molecules == 1 and st.n_molecules_other_key == 2
    st2 = consensus_bam(inp, str(tmp_path / "out2.bam"), reference=fa, primer_length=9, span_key="tlen")
    assert st2.n_molecules == 2 and st2.n_molecules_other_key == 1


def test_umi_from_qname(tmp_path, ref_fasta):
    fa, ref = ref_fasta
    h = _hdr()
    pairs = []
    for i in range(3):
        pr = _pair(h, ref, f"q{i}_ACGTACGTACGT", None, 100, 300)
        pairs.append(pr)
    inp = _write(tmp_path / "in.bam", pairs, h)
    st = consensus_bam(inp, str(tmp_path / "out.bam"), reference=fa, umi_from_qname=True)
    assert st.n_molecules == 1 and st.family_size_hist[3] == 1
    st0 = consensus_bam(inp, str(tmp_path / "out0.bam"), reference=fa, umi_from_qname=False)
    assert st0.n_molecules == 0 and st0.n_no_umi == 3


def test_primer_off_leaves_reads_unmasked(tmp_path, ref_fasta):
    fa, ref = ref_fasta
    h = _hdr()
    pairs = [_pair(h, ref, f"m{i}", "ACGTACGTACGT", 100, 300) for i in range(2)]
    inp = _write(tmp_path / "in.bam", pairs, h)
    out = str(tmp_path / "out.bam")
    consensus_bam(inp, out, reference=fa, primer_length=0)
    assert all("N" not in r.query_sequence for r in _reads(out))
