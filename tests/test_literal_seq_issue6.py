"""#6: the multialigned BAM must carry the input read's literal bases, oriented by each record's own flag."""
import gzip

import pysam
import pytest

from rectify.core.align.literal_seq import rebuild_record_sequence, restore_literal_sequences

READ = "ACGTTGCAAGGCTTACGATC"          # 20 nt, as given to the aligners (RNA 5'->3')
QUAL = "".join(chr(33 + q) for q in range(10, 30))


def _rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _header():
    return pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrT", "LN": 1000}]})


def _rec(header, name, seq, cigar, rev=False, pos=100, unmapped=False):
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.flag = (16 if rev else 0) | (4 if unmapped else 0)
    if not unmapped:
        r.reference_id = 0
        r.reference_start = pos
        r.cigarstring = cigar
        r.mapping_quality = 60
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    r.set_tag("Xa", "uLTRA")
    return r


def _eq(s, keep):
    """`s` with '=' everywhere except the positions in `keep` (as calmd -e leaves it)."""
    return "".join(c if i in keep else "=" for i, c in enumerate(s))


def test_plus_strand_eq_encoded():
    h = _header()
    r = _rec(h, "r1", _eq(READ, {3, 9}), "20M")
    assert rebuild_record_sequence(r, READ, QUAL) == "rebuilt"
    assert r.query_sequence == READ
    assert pysam.qualities_to_qualitystring(r.query_qualities) == QUAL
    assert r.get_tag("Xa") == "uLTRA" and r.cigarstring == "20M"


def test_minus_strand_is_reverse_complemented():
    h = _header()
    r = _rec(h, "r2", _eq(_rc(READ), {0, 5}), "20M", rev=True)
    assert rebuild_record_sequence(r, READ, QUAL) == "rebuilt"
    assert r.query_sequence == _rc(READ)
    assert pysam.qualities_to_qualitystring(r.query_qualities) == QUAL[::-1]


def test_seq_stored_opposite_to_flag_is_reoriented():
    h = _header()
    r = _rec(h, "r3", _eq(_rc(READ), {1, 2, 7, 11, 13}), "20M", rev=False)   # flag says +, SEQ stored as the reverse complement
    assert rebuild_record_sequence(r, READ, QUAL) == "rebuilt_flipped"
    assert r.query_sequence == READ


def test_eq_inside_soft_clip_is_filled():
    h = _header()
    r = _rec(h, "r4", _eq(READ, {6, 8}), "5S10M5S")
    assert rebuild_record_sequence(r, READ, QUAL) == "rebuilt"
    assert r.query_sequence == READ and "=" not in r.query_sequence


def test_hard_clips_take_the_spanned_part():
    h = _header()
    r = _rec(h, "r5", _eq(READ[3:15], {0}), "3H12M5H")
    assert rebuild_record_sequence(r, READ, QUAL) == "rebuilt"
    assert r.query_sequence == READ[3:15]
    assert pysam.qualities_to_qualitystring(r.query_qualities) == QUAL[3:15]


def test_inconsistent_and_length_mismatch_left_unchanged():
    h = _header()
    bad = "T" * 20                                      # spelled bases agree with neither orientation
    r = _rec(h, "r6", bad, "20M")
    assert rebuild_record_sequence(r, READ, QUAL) == "unchanged_inconsistent"
    assert r.query_sequence == bad
    r2 = _rec(h, "r7", "=" * 18, "18M")
    assert rebuild_record_sequence(r2, READ, QUAL) == "unchanged_length"
    assert r2.query_sequence == "=" * 18


def test_unmapped_record_gets_the_read():
    h = _header()
    r = _rec(h, "r8", "N" * 20, None, unmapped=True)
    assert rebuild_record_sequence(r, READ, QUAL) == "rebuilt"
    assert r.query_sequence == READ


def _write_inputs(tmp_path, names=("a", "b", "c", "d")):
    reads = {
        "a": READ,
        "b": "GGGTTTCCCAAAGGGTTTCCAA",
        "c": "ACACACGTGTGTTTAAACCC",
        "d": "TTTTGGGGCCCCAAAATTTT",
    }
    fq = tmp_path / "reads.fastq.gz"
    with gzip.open(fq, "wt") as fh:
        for n in names:
            fh.write(f"@{n} runid=x ch=1\n{reads[n]}\n+\n{'I' * len(reads[n])}\n")
    h = _header()
    recs = [
        _rec(h, "a", _eq(READ, {2}), "20M", pos=10),
        _rec(h, "b", _eq(_rc(reads["b"]), {4}), "4S18M", rev=True, pos=20),
        _rec(h, "c", _eq(_rc(reads["c"]), {0, 1}), "20M", rev=False, pos=30),     # flipped
        _rec(h, "d", _eq(reads["d"], {5}), "5S10M5S", pos=40),                    # '=' in clips
    ]
    bam = tmp_path / "x.multialigned.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=h) as out:
        for r in recs:
            out.write(r)
    pysam.index(str(bam))
    return bam, fq, reads


@pytest.mark.parametrize("max_reads", [1_000_000, 1])
def test_restore_rewrites_every_record_and_keeps_order(tmp_path, max_reads):
    bam, fq, reads = _write_inputs(tmp_path)
    st = restore_literal_sequences(str(bam), str(fq), max_reads_in_memory=max_reads)
    assert st["records"] == 4 and st["rebuilt"] == 3 and st["rebuilt_flipped"] == 1
    with pysam.AlignmentFile(str(bam)) as f:
        got = [(r.query_name, r.reference_start, r.query_sequence, r.is_reverse) for r in f]
    assert [g[1] for g in got] == [10, 20, 30, 40]
    for name, _, seq, rev in got:
        assert seq == (_rc(reads[name]) if rev else reads[name])
    assert (tmp_path / "x.multialigned.bam.bai").exists()


def test_restore_refuses_when_reads_are_missing(tmp_path):
    bam, fq, _ = _write_inputs(tmp_path)
    fq2 = tmp_path / "partial.fastq.gz"
    with gzip.open(fq2, "wt") as fh:
        fh.write(f"@a\n{READ}\n+\n{'I' * 20}\n")
    before = bam.read_bytes()
    with pytest.raises(RuntimeError):
        restore_literal_sequences(str(bam), str(fq2))
    assert bam.read_bytes() == before
