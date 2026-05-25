"""Tests for the min-anchor-overhang pool-junction filter.

Aligners (notably gapmm2) occasionally emit a tiny-anchor spurious split — e.g.
``4M4250N223M``, a 4 bp anchor before a 4.25 kb "intron" — which then propagates
to many cleanly-aligned reads during 3'SS rescue. ``collect_junction_counts_from_bam``
must drop such N-ops from the pool while keeping well-anchored real junctions.
"""

import pysam
import pytest

from rectify.core.splice.junction_scoring import (
    _is_low_complexity_anchor,
    _junction_anchor_ok,
    _junction_worst_flank_periodicity,
    _merge_del_into_intron,
    _seq_periodicity,
    build_junction_pool,
    collect_junction_counts_from_bam,
    DEFAULT_MIN_JUNCTION_ANCHOR,
)


# --- adjacent D/N collapse ---------------------------------------------------

@pytest.mark.parametrize("cigar, expected", [
    ([(0, 92), (2, 1), (3, 112), (0, 57)], [(0, 92), (3, 113), (0, 57)]),   # 1D before N
    ([(0, 20), (3, 96), (2, 2), (0, 20)], [(0, 20), (3, 98), (0, 20)]),     # D after N
    ([(0, 10), (2, 2), (3, 100), (2, 3), (0, 10)], [(0, 10), (3, 105), (0, 10)]),  # D on both sides
    ([(0, 10), (2, 3), (0, 10)], [(0, 10), (2, 3), (0, 10)]),               # pure deletion untouched
    ([(0, 10), (3, 100), (0, 10)], [(0, 10), (3, 100), (0, 10)]),           # already clean
    ([(4, 5), (0, 10), (1, 2), (3, 50), (0, 10)], [(4, 5), (0, 10), (1, 2), (3, 50), (0, 10)]),  # I not merged
])
def test_merge_del_into_intron(cigar, expected):
    assert _merge_del_into_intron(cigar) == expected


# --- low-complexity detection ------------------------------------------------

@pytest.mark.parametrize("seq, expected", [
    ("AAAAAAAAAA", True),     # homopolymer (period 1)
    ("ATATATATAT", True),     # dinucleotide (period 2)
    ("ATGATGATGA", True),     # trinucleotide (period 3)
    ("", True),               # empty = untrustworthy
    ("ACGTGCATGC", False),    # diverse
    ("ATGCATGCAT", False),    # period 4 — not flagged by the 1/2/3 rule
])
def test_is_low_complexity_anchor(seq, expected):
    assert _is_low_complexity_anchor(seq) is expected


# --- anchor gate -------------------------------------------------------------

# A 60-base diverse query is plenty for 10 bp anchors on both flanks.
_Q = "ACGTGCATGCttacggatccGGCATGCATGCacgtgcatgcTTACGGATCCAggcatgcatgc"


def test_tiny_left_anchor_rejected():
    # 4M 4250N 223M — the gapmm2 artifact. Left anchor (4) < 10.
    cigar = [(0, 4), (3, 4250), (0, 223)]
    assert _junction_anchor_ok(cigar, 1, 4, _Q, 10) is False


def test_both_anchors_long_and_diverse_accepted():
    # 20M 96N 20M with non-repetitive flanks.
    cigar = [(0, 20), (3, 96), (0, 20)]
    assert _junction_anchor_ok(cigar, 1, 20, _Q, 10) is True


def test_junction_at_read_edge_rejected():
    cigar = [(3, 100), (0, 30)]            # leading N — no upstream exon
    assert _junction_anchor_ok(cigar, 0, 0, _Q, 10) is False
    cigar2 = [(0, 30), (3, 100)]           # trailing N — no downstream exon
    assert _junction_anchor_ok(cigar2, 1, 30, _Q, 10) is False


def test_indel_or_mismatch_op_adjacent_rejected():
    # Insertion immediately before the N → not a clean match anchor.
    cigar = [(0, 20), (1, 5), (3, 96), (0, 20)]
    assert _junction_anchor_ok(cigar, 2, 25, _Q, 10) is False
    # Mismatch op (X) immediately before the N (=/X aligner).
    cigar_x = [(7, 20), (8, 3), (3, 96), (7, 20)]
    assert _junction_anchor_ok(cigar_x, 2, 23, _Q, 10) is False


def test_low_complexity_anchor_rejected():
    q = "A" * 15 + "C" * 15           # homopolymer flanks
    cigar = [(0, 15), (3, 96), (0, 15)]
    assert _junction_anchor_ok(cigar, 1, 15, q, 10) is False


def test_disabled_when_overhang_zero():
    cigar = [(0, 4), (3, 4250), (0, 223)]
    assert _junction_anchor_ok(cigar, 1, 4, _Q, 0) is True


# --- full BAM integration ----------------------------------------------------

def _write_bam(path, reads):
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrI", "LN": 230218}]}
    )
    with pysam.AlignmentFile(str(path), "wb", header=hdr) as out:
        for i, (start, cigar, seq) in enumerate(reads):
            r = pysam.AlignedSegment(hdr)
            r.query_name = f"r{i}"
            r.reference_id = 0
            r.reference_start = start
            r.mapping_quality = 60
            r.cigartuples = cigar
            r.query_sequence = seq
            out.write(r)
    pysam.index(str(path))


def test_collect_counts_drops_tiny_anchor_keeps_real(tmp_path):
    bam = tmp_path / "t.bam"
    real_q = _Q[:20] + _Q[:20]          # 40 bases for a 20M..20M read
    spurious_q = "ACGT" + _Q[:30]       # 4 + 30 bases for 4M..30M read
    _write_bam(bam, [
        # real, well-anchored junction at 1000: spans [1020, 1116)
        (1000, [(0, 20), (3, 96), (0, 20)], real_q),
        # gapmm2-style tiny-anchor split at 5000: 4M 4250N 30M
        (5000, [(0, 4), (3, 4250), (0, 30)], spurious_q),
    ])
    counts = collect_junction_counts_from_bam(str(bam))
    assert ("chrI", 1020, 1116) in counts          # real kept
    assert ("chrI", 5004, 9254) not in counts       # tiny-anchor dropped

    # With the filter disabled, both appear.
    counts_off = collect_junction_counts_from_bam(str(bam), min_anchor_overhang=0)
    assert ("chrI", 1020, 1116) in counts_off
    assert ("chrI", 5004, 9254) in counts_off


def test_default_overhang_is_ten():
    assert DEFAULT_MIN_JUNCTION_ANCHOR == 10


def test_boundary_deletion_absorbed_then_anchor_passes(tmp_path):
    # A nanopore read with a 1bp deletion at the splice boundary (92M 1D 112N
    # 57M) must, after D/N collapse, present a clean 92M anchor and seed the
    # (merged) junction — not be rejected by the anchor filter.
    bam = tmp_path / "bd.bam"
    real_q = (_Q * 3)[:149]             # 92 + 57 query bases (D/N consume no query)
    _write_bam(bam, [
        (2000, [(0, 92), (2, 1), (3, 112), (0, 57)], real_q),
    ])
    counts = collect_junction_counts_from_bam(str(bam))
    # merged intron spans [2000+92, 2000+92+113) = [2092, 2205)
    assert ("chrI", 2092, 2205) in counts


# --- graded periodicity helpers (consensus complexity dimension) -------------

@pytest.mark.parametrize("seq, lo, hi", [
    ("ACGTGCATGC", 0.0, 0.6),   # diverse → low
    ("AAAAAAAAAA", 1.0, 1.0),   # homopolymer → perfect
    ("ATATATATAT", 1.0, 1.0),   # dinucleotide → perfect
    ("ATGATGATGA", 1.0, 1.0),   # trinucleotide → perfect
    ("", 1.0, 1.0),             # empty → treated as max (untrustworthy)
])
def test_seq_periodicity_range(seq, lo, hi):
    assert lo <= _seq_periodicity(seq) <= hi


def test_worst_flank_periodicity_atract_flagged_but_diverse_prefix_rescues():
    intron = "N" * 100
    # Both flanks diverse → low (well below any sane T2).
    up = "GCATGCTAGCATCGATCGTAGC"
    dn = "TCGATCGATGCATTAGCATCGA"
    g = {"chrI": up + intron + dn}
    js, je = len(up), len(up) + len(intron)
    assert _junction_worst_flank_periodicity(g, "chrI", js, je) < 0.9

    # Both flanks pure A-tract → worst flank ~1.0 (low-complexity, flaggable).
    g2 = {"chrI": "A" * 22 + intron + "A" * 22}
    assert _junction_worst_flank_periodicity(g2, "chrI", 22, 22 + len(intron)) >= 0.9

    # Diverse prefix abutting an A-tract on one flank still rescues that flank
    # (a complex window exists within reach); the diverse other flank is fine.
    g3 = {"chrI": ("AGCTGAGTC" + "A" * 13) + intron + dn}
    assert _junction_worst_flank_periodicity(g3, "chrI", 22, 22 + len(intron)) < 0.9

    # Missing genome/chrom → 0.0 (cannot assess; caller must not flag).
    assert _junction_worst_flank_periodicity({}, "chrI", js, je) == 0.0


# --- pool concordance relaxation ---------------------------------------------

def test_family_concordance_relaxation(tmp_path):
    # Short anchors (4M 96N 4M, 4 bp < 10) fail the hard floor; whether the
    # junction is relaxed into the pool depends on cross-FAMILY concordance.
    cigar_short = [(0, 4), (3, 96), (0, 4)]
    q = "ACGTGCAT"                      # 8 query bases (4M + 4M)
    a = tmp_path / "a.bam"
    b = tmp_path / "b.bam"
    _write_bam(a, [(1000, cigar_short, q)])
    _write_bam(b, [(1000, cigar_short, q)])
    jx = ("chrI", 1004, 1100)

    # Two INDEPENDENT families (minimap2 + deSALT) → relaxed in.
    pool_indep, _ = build_junction_pool(
        [str(a), str(b)], set(), aligner_families=["minimap2", "deSALT"])
    assert jx in pool_indep

    # minimap2 + gapmm2 = ONE family (gapmm2 wraps minimap2) → NOT relaxed.
    pool_same, _ = build_junction_pool(
        [str(a), str(b)], set(), aligner_families=["minimap2", "gapmm2"])
    assert jx not in pool_same

    # Single aligner → no cross-family concordance → dropped.
    pool_solo, _ = build_junction_pool(
        [str(a)], set(), aligner_families=["minimap2"])
    assert jx not in pool_solo

    # Relaxation disabled → even independent families can't rescue a short anchor.
    pool_off, _ = build_junction_pool(
        [str(a), str(b)], set(),
        aligner_families=["minimap2", "deSALT"], relax_min_families=0)
    assert jx not in pool_off


def test_unknown_family_names_fall_back_to_hard_floor(tmp_path):
    # Unrecognised BAM names (no aligner token) → "?" sentinel, which never
    # counts toward concordance → pure hard floor even with 2 BAMs.
    cigar_short = [(0, 4), (3, 96), (0, 4)]
    q = "ACGTGCAT"
    a = tmp_path / "mystery_one.bam"
    b = tmp_path / "mystery_two.bam"
    _write_bam(a, [(1000, cigar_short, q)])
    _write_bam(b, [(1000, cigar_short, q)])
    pool, _ = build_junction_pool([str(a), str(b)], set())
    assert ("chrI", 1004, 1100) not in pool


def test_well_anchored_junction_kept_by_single_aligner(tmp_path):
    # A clean >= 10 bp anchored junction needs no concordance — one aligner is
    # enough (the hard-floor path, not the relaxation).
    cigar_ok = [(0, 20), (3, 96), (0, 20)]
    q = _Q[:20] + _Q[:20]
    a = tmp_path / "solo.minimap2.bam"
    _write_bam(a, [(3000, cigar_ok, q)])
    pool, _ = build_junction_pool([str(a)], set())
    assert ("chrI", 3020, 3116) in pool
