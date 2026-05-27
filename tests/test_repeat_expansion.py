"""Tests for read-side repeat-expansion artifact detection (repeat_expansion.py).

Anchored on real cases observed in human RNA004 DRS soft-clips (2026-05-26):
genuine (AAG)n/(CTT)n expansions must be flagged; poly-A/poly-T tails (incl. with
sparse basecall errors) must NOT be flagged, since those are the poly-A walkback's
job, not the rescue short-circuit's.
"""

from rectify.core.splice.repeat_expansion import (
    dominant_repeat_period,
    is_repeat_expansion,
    _canonical,
)


def test_canonical_rotation():
    assert _canonical("AAG") == "AAG"
    assert _canonical("AGA") == "AAG"
    assert _canonical("GAA") == "AAG"
    assert _canonical("CTT") == "CTT"


def test_flags_genuine_aag_expansion():
    # 48 bp (AAG)n — a real soft-clip example from SMA_GSB2394
    seq = "AAG" * 16
    assert is_repeat_expansion(seq)
    res = dominant_repeat_period(seq)
    assert res is not None and res[0] == 3 and res[1] == "AAG"


def test_flags_ctt_expansion():
    # (CTT)n as seen on reverse-strand reads (RNA-space AAG)
    seq = "CTT" * 25
    assert is_repeat_expansion(seq)
    assert dominant_repeat_period(seq)[1] == "CTT"


def test_flags_a_rich_aag():
    # A-rich (AAG)n with an A lead-in (maxbase ~0.69) — still genuine, must flag
    seq = "AAAAAG" + "AAG" * 14
    assert is_repeat_expansion(seq)


def test_rejects_pure_polya():
    # Pure poly-A passes a period-2 "AA" family but the motif is single-base -> reject
    seq = "A" * 60
    assert not is_repeat_expansion(seq)
    # the underlying period detector DOES see a homopolymer family
    res = dominant_repeat_period(seq)
    assert res is not None and len(set(res[1])) == 1


def test_rejects_polyt():
    assert not is_repeat_expansion("T" * 45)


def test_rejects_polya_with_sparse_errors():
    # poly-A with sparse G basecall errors: dominant period-3 family is "AAA" -> reject
    seq = "AAAAAGAAAAAAAGAAAAAAAAAGAAAAAAAGAAAAA"
    assert not is_repeat_expansion(seq)


def test_length_floor():
    # below MIN_LEN (30) even a clean (AAG)n is not flagged
    assert not is_repeat_expansion("AAG" * 9)   # 27 bp
    assert is_repeat_expansion("AAG" * 10)      # 30 bp


def test_flags_dinucleotide_expansion():
    # (AT)n di-nucleotide expansion is multi-base -> flagged
    assert is_repeat_expansion("AT" * 20)


def test_rejects_diverse_sequence():
    # genuine genomic exon sequence (diverse) must not be flagged
    seq = "ACGTACAGTCCGATTGCAACGTGACTGCATACGGATCAGT"
    assert not is_repeat_expansion(seq)


def test_empty_and_none():
    assert not is_repeat_expansion("")
    assert not is_repeat_expansion(None)


def test_impure_below_threshold_not_flagged():
    # an AAG run diluted with diverse sequence below 75% purity is not flagged
    seq = "AAGAAGAAG" + "CGTACGTACGTACGTACGTAC" + "AAGAAG"
    assert not is_repeat_expansion(seq)
