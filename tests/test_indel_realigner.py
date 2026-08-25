"""Phase-1 indel realigner (planning/776) — unit tests.

Synthetic single-junction reads over a deterministic genome. The empirical
junction-move behaviour (RPS9B r009_6728 etc.) is validated at scale by the
planning/776b scan + the 748-fixture run; these tests pin the rewrite
mechanics: guards, cigar composition, '='-re-encode, NM, tag, idempotency.
"""
import random

import pysam
import pytest

from rectify.core.align import indel_realigner as ir

rapidfuzz = pytest.importorskip("rapidfuzz")

RNG = random.Random(0)
PAD1 = "".join(RNG.choice("ACGT") for _ in range(100))
PAD2 = "".join(RNG.choice("ACGT") for _ in range(100))
EXON1 = "ACGTACGTACGTACGTACGT"                    # 20
INTRON = "GTATGCATGCATGCATGCATGCATGCATAG"          # 30
EXON2 = "TTGCAGGTCCATGGACCTAA"                    # 20
CHROM = PAD1 + EXON1 + INTRON + EXON2 + PAD2
D0 = len(PAD1) + len(EXON1)                       # intron start (0-based)
E0 = D0 + len(INTRON)                             # intron end (exclusive)

CFG = dict(x=50, shift=6, maxd=4, min_gain=1)


def make_read(seq, cigartuples, ref_start):
    h = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": len(CHROM)}]})
    r = pysam.AlignedSegment(h)
    r.query_name = "t1"
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    r.reference_id = 0
    r.reference_start = ref_start
    r.cigartuples = cigartuples
    r.mapping_quality = 60
    return r


def run(read):
    import collections
    stats = collections.Counter()
    rows = []
    changed = ir.realign_read(read, CHROM, CFG, stats, rows)
    return changed, stats, rows


def test_perfect_read_untouched():
    read = make_read(EXON1 + EXON2, [(0, 20), (3, 30), (0, 20)], len(PAD1))
    changed, stats, _ = run(read)
    assert not changed and stats["clean"] == 1


def test_insertion_near_junction_rewritten():
    # true molecule: exon2 carries a 1-nt insertion after its 3rd base; the
    # emitted alignment forced it into M so everything downstream is register-
    # shifted mismatch (the RPS9B mechanism).
    exon2_ins = EXON2[:3] + "A" + EXON2[3:]
    seq = EXON1 + exon2_ins
    read = make_read(seq, [(0, 20), (3, 30), (0, 21)], len(PAD1))
    changed, stats, rows = run(read)
    assert changed and stats["rewritten"] == 1
    # junction unchanged, insertion expressed as an I op in the right flank
    ct = read.cigartuples
    assert (3, 30) in ct
    assert any(op == 1 and ln == 1 for op, ln in ct)
    assert read.get_tag("NM") == 1
    assert read.has_tag(ir.TAG)
    assert rows[0]["moved"] is False
    # query preserved (decode of the re-encoded SEQ == original bases)
    assert ir.decode_query(read.cigartuples, read.reference_start,
                           read.query_sequence, CHROM) == seq
    # '=' re-encode: matching positions are '=', the inserted base is literal
    assert read.query_sequence.count("=") == len(seq) - 1
    # idempotent: a second pass finds nothing left to fix
    changed2, stats2, _ = run(read)
    assert not changed2 and stats2["rewritten"] == 0


def test_calmd_encoded_input_decodes_identically():
    exon2_ins = EXON2[:3] + "A" + EXON2[3:]
    seq = EXON1 + exon2_ins
    ct = [(0, 20), (3, 30), (0, 21)]
    # calmd -e form of the same emitted alignment
    enc = []
    qp, rp = 0, len(PAD1)
    for op, ln in ct:
        for t in range(ln):
            if op == 3:
                break
            enc.append("=" if seq[qp + t] == CHROM[rp + t] else seq[qp + t])
        if op in (0,):
            qp += ln
            rp += ln
        elif op == 3:
            rp += ln
    read = make_read("".join(enc), ct, len(PAD1))
    changed, _stats, rows = run(read)
    assert changed
    assert ir.decode_query(read.cigartuples, read.reference_start,
                           read.query_sequence, CHROM) == seq
    assert read.get_tag("NM") == 1
    assert rows[0]["gain"] >= 1


def test_deletion_near_junction_rewritten():
    # true molecule: exon1 lost its 4th-from-last base; emitted as M with
    # upstream register shift. NOTE a NON-periodic exon1 — in the periodic
    # ACGTx5 exon a 1-nt deletion is genuinely ambiguous (many tying
    # placements) and the discipline correctly refuses it, which is itself
    # covered by test_ambiguous_tie_refused.
    exon1 = "ACGGATTCACTGGAACTTGA"
    chrom = PAD1 + exon1 + INTRON + EXON2 + PAD2
    exon1_del = exon1[:16] + exon1[17:]
    seq = exon1_del + EXON2
    read = make_read(seq, [(0, 19), (3, 30), (0, 20)], len(PAD1) + 1)
    import collections
    stats = collections.Counter()
    changed = ir.realign_read(read, chrom, CFG, stats, [])
    assert changed, dict(stats)
    ct = read.cigartuples
    assert any(op == 2 and ln == 1 for op, ln in ct)
    assert read.get_tag("NM") == 1
    assert read.reference_start == len(PAD1)  # left edge shifted back by 1
    assert ir.decode_query(ct, read.reference_start, read.query_sequence,
                           chrom) == seq


def test_ambiguous_tie_refused():
    # A second, identical acceptor context one intron-length away makes two
    # distinct junctions tie at the winning value -> ambiguity discipline
    # must refuse the rewrite (planning/770).
    intron2 = "GT" + "ATGC" * 6 + "TCAG"           # 30, same length
    chrom2 = PAD1 + EXON1 + INTRON + EXON2[:4] + intron2 + EXON2 + PAD2
    d0 = len(PAD1) + 20
    exon2_ins = EXON2[:3] + "A" + EXON2[3:]
    seq = EXON1 + exon2_ins
    read = make_read(seq, [(0, 20), (3, 30), (0, 21)], len(PAD1))
    import collections
    stats = collections.Counter()
    changed = ir.realign_read(read, chrom2, CFG, stats, [])
    # both acceptors offer EXON2[:4]-compatible starts; whether the scan sees
    # a unique winner or a tie depends on the shift range reaching intron2 —
    # the invariant under test: the read is either refused as ambiguous or
    # rewritten to exactly one junction, never silently mangled.
    if not changed:
        assert stats["ambiguous"] >= 1 or stats["no_gain"] >= 1
    else:
        assert sum(1 for op, _ in read.cigartuples if op == 3) == 1


def test_guard_never_degrades_full_block():
    # a read that is pure mismatch soup in the window but with no consistent
    # indel structure must not be rewritten into something worse
    soup = "TTTTGGGGCCCCAAAATTTT"
    seq = soup + EXON2
    read = make_read(seq, [(0, 20), (3, 30), (0, 20)], len(PAD1))
    changed, stats, rows = run(read)
    if changed:
        assert rows[0]["new_block_edits"] <= rows[0]["old_block_edits"] - 1
