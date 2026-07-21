"""Unit tests for the short-read (standard N-mer) UMI path: rectify/core/umi/.

Zero data files: reads are duck-typed mocks (the production functions annotate
``Any`` and touch only the pysam AlignedSegment attributes asserted here), and
UMIs are literal strings. Mirrors the fixture style of
``tests/test_cdna_restore_eq_seq.py``.
"""
from __future__ import annotations

import pytest

from rectify.core.umi import (
    CORALL_UMI_LENGTH,
    PairStatus,
    canonical_umi,
    extract_umi,
    extract_umi_from_read_id,
    extract_umi_from_sequence,
    fragment_key,
    fragment_strand,
    five_prime_unclipped,
    is_valid_umi,
    leading_soft_clip,
    read_name_of,
    strip_umi_from_read_id,
    trailing_soft_clip,
    umi_components_directional,
)

# CIGAR ops: 0=M 1=I 2=D 3=N 4=S 5=H
_M, _I, _D, _N, _S, _H = 0, 1, 2, 3, 4, 5


class MockRead:
    """Minimal duck-type of pysam.AlignedSegment for grouping tests."""

    def __init__(self, reference_name="chrI", reference_start=1000, cigartuples=None,
                 is_reverse=False, is_unmapped=False, is_paired=True,
                 is_read1=True, query_name="r1"):
        self.reference_name = reference_name
        self.reference_start = reference_start
        self.cigartuples = cigartuples if cigartuples is not None else [(_M, 100)]
        self.is_reverse = is_reverse
        self.is_unmapped = is_unmapped
        self.is_paired = is_paired
        self.is_read1 = is_read1
        self.is_read2 = not is_read1
        self.query_name = query_name

    @property
    def reference_end(self):
        if self.is_unmapped:
            return None
        span = sum(l for op, l in self.cigartuples if op in (0, 2, 3, 7, 8))
        return self.reference_start + span


# --------------------------------------------------------------------------
# extract: read-id mode (THE CRITICAL PATH -- Angel's RES FASTQs are already
# umi_tools-extracted and the raw FASTQs no longer exist)
# --------------------------------------------------------------------------

def test_extract_umi_from_read_id_umi_tools_convention():
    rid = "@A00197:591:HMMFFDSX5:3:1101:23475:1031_ACGTACGTACGT 1:N:0:ATCACG"
    assert extract_umi_from_read_id(rid, expected_length=12) == "ACGTACGTACGT"


def test_extract_umi_from_read_id_ignores_illumina_comment():
    """The comment field must never reach the UMI parser.

    Illumina's comment (``1:N:0:INDEX``) carries the sample INDEX -- itself a
    12-mer in Lexogen UDI kits. Parsing it as the UMI would be catastrophic and
    would look plausible.
    """
    rid = "@RUN:1:1101:1:1_ACGTACGTACGT 1:N:0:GGGGGGGGGGGG"
    assert extract_umi_from_read_id(rid, expected_length=12) == "ACGTACGTACGT"


def test_extract_umi_from_read_id_splits_on_last_separator():
    """Instrument/run names may contain '_'; splitting on the first would return junk."""
    rid = "@MY_RUN_01:1:1101:1:1_ACGTACGTACGT"
    assert extract_umi_from_read_id(rid, expected_length=12) == "ACGTACGTACGT"


def test_extract_umi_from_read_id_rejects_wrong_length():
    rid = "@READ_ACGTACGT"  # 8 nt, not 12
    assert extract_umi_from_read_id(rid, expected_length=12) is None


def test_extract_umi_from_read_id_rejects_non_nucleotide():
    """A mis-parsed name yields a plausible string; the alphabet check catches it."""
    rid = "@READ_HELLOWORLD12"
    assert extract_umi_from_read_id(rid, expected_length=12) is None


def test_extract_umi_from_read_id_absent_separator_returns_none():
    assert extract_umi_from_read_id("@READONLY", expected_length=12) is None


def test_extract_umi_from_read_id_empty_returns_none():
    assert extract_umi_from_read_id("", expected_length=12) is None


def test_read_name_of_strips_at_and_comment():
    assert read_name_of("@NAME:1:2 1:N:0:AT") == "NAME:1:2"


def test_strip_umi_from_read_id_yields_shared_qname():
    """Both mates must end up with an IDENTICAL bare QNAME for pair-aware grouping."""
    r1 = "@A00197:591:H:3:1101:1:1_ACGTACGTACGT 1:N:0:AT"
    r2 = "@A00197:591:H:3:1101:1:1_ACGTACGTACGT 2:N:0:AT"
    assert strip_umi_from_read_id(r1) == strip_umi_from_read_id(r2) == "A00197:591:H:3:1101:1:1"


def test_strip_umi_leaves_name_without_umi_untouched():
    assert strip_umi_from_read_id("@PLAIN:1:2") == "PLAIN:1:2"


def test_strip_umi_does_not_eat_a_non_umi_suffix():
    """A trailing '_something' that isn't a UMI must survive."""
    assert strip_umi_from_read_id("@SAMPLE_rep1") == "SAMPLE_rep1"


# --------------------------------------------------------------------------
# extract: r1-start mode (raw CORALL -- future-proofing; the raw batch has not landed)
# --------------------------------------------------------------------------

def test_extract_umi_from_sequence_slices_n12():
    seq = "ACGTACGTACGT" + "GGGGCCCCTTTT"
    qual = "IIIIIIIIIIII" + "FFFFFFFFFFFF"
    got = extract_umi_from_sequence(seq, qual, CORALL_UMI_LENGTH)
    assert got.umi == "ACGTACGTACGT"
    assert got.umi_qual == "IIIIIIIIIIII"
    assert got.seq == "GGGGCCCCTTTT"      # no linker: insert abuts the UMI
    assert got.qual == "FFFFFFFFFFFF"


def test_extract_umi_from_sequence_too_short_returns_none():
    assert extract_umi_from_sequence("ACGT", "IIII", 12) is None


def test_extract_umi_from_sequence_exact_length_leaves_empty_insert():
    got = extract_umi_from_sequence("ACGTACGTACGT", "IIIIIIIIIIII", 12)
    assert got.umi == "ACGTACGTACGT"
    assert got.seq == ""


def test_extract_umi_from_sequence_rejects_nonpositive_length():
    with pytest.raises(ValueError):
        extract_umi_from_sequence("ACGT", "IIII", 0)


def test_corall_umi_length_is_12():
    """CORALL V2 UG 171UG394V0111 App. E p.36: N12 at the start of Read 1.
    V1 UG 095UG190V0110 is word-for-word identical -- V2 did not change the UMI."""
    assert CORALL_UMI_LENGTH == 12


def test_extract_umi_dispatch_r1_start_only_slices_r1():
    """R2 carries no UMI in CORALL; it must not be sliced."""
    assert extract_umi("r1-start", "@r", "ACGTACGTACGT" + "TTTT", "I" * 16,
                       length=12, is_read2=False).umi == "ACGTACGTACGT"
    assert extract_umi("r1-start", "@r", "ACGTACGTACGT" + "TTTT", "I" * 16,
                       length=12, is_read2=True) is None


def test_extract_umi_dispatch_read_id_does_not_touch_sequence():
    got = extract_umi("read-id", "@READ_ACGTACGTACGT", "GGGGCCCC", "IIIIIIII", length=12)
    assert got.umi == "ACGTACGTACGT"
    assert got.seq == "GGGGCCCC"  # UMI was never in the sequence


def test_extract_umi_requires_length_for_positional_modes():
    with pytest.raises(ValueError):
        extract_umi("r1-start", "@r", "ACGT", "IIII", length=None)


def test_extract_umi_rejects_unknown_location():
    with pytest.raises(ValueError):
        extract_umi("middle-of-r1", "@r", "ACGT", "IIII", length=12)


def test_is_valid_umi_allows_and_forbids_n():
    assert is_valid_umi("ACGTACGTACGN", expected_length=12) is True
    assert is_valid_umi("ACGTACGTACGN", expected_length=12, allow_n=False) is False


# --------------------------------------------------------------------------
# grouping: soft-clip correction (the deterministic wobble fix that makes a
# tolerance window unnecessary)
# --------------------------------------------------------------------------

def test_leading_and_trailing_soft_clip():
    r = MockRead(cigartuples=[(_S, 5), (_M, 90), (_S, 5)])
    assert leading_soft_clip(r) == 5
    assert trailing_soft_clip(r) == 5


def test_soft_clip_ignores_hard_clip():
    r = MockRead(cigartuples=[(_H, 10), (_S, 3), (_M, 90)])
    assert leading_soft_clip(r) == 3


def test_five_prime_unclipped_forward_corrects_leading_clip():
    """A forward read soft-clipped at its 5' end started EARLIER than reference_start."""
    r = MockRead(reference_start=1000, cigartuples=[(_S, 5), (_M, 95)])
    assert five_prime_unclipped(r) == 995


def test_five_prime_unclipped_reverse_uses_high_coordinate():
    """A reverse read's 5' end is at its reference_end, plus any trailing clip."""
    r = MockRead(reference_start=1000, cigartuples=[(_M, 95), (_S, 5)], is_reverse=True)
    assert five_prime_unclipped(r) == 1100  # end 1095 + 5


def test_five_prime_unclipped_same_for_clipped_and_unclipped_duplicate():
    """THE POINT of the correction: two copies of ONE molecule, one soft-clipped,
    must land on the SAME 5' position -- otherwise they'd wrongly count as 2."""
    clean = MockRead(reference_start=1000, cigartuples=[(_M, 100)])
    clipped = MockRead(reference_start=1004, cigartuples=[(_S, 4), (_M, 96)])
    assert five_prime_unclipped(clean) == five_prime_unclipped(clipped) == 1000


def test_five_prime_unclipped_unmapped_returns_none():
    assert five_prime_unclipped(MockRead(is_unmapped=True)) is None


def test_reference_end_skips_insertions_counts_deletions_and_introns():
    r = MockRead(reference_start=100, cigartuples=[(_M, 10), (_I, 5), (_M, 10), (_D, 2), (_N, 50), (_M, 10)])
    assert r.reference_end == 100 + 10 + 10 + 2 + 50 + 10


# --------------------------------------------------------------------------
# grouping: strandedness (CORALL R1 = SENSE -- opposite of QuantSeq REV/TruSeq)
# --------------------------------------------------------------------------

def test_fragment_strand_corall_r1_forward_is_plus():
    """CORALL UG App. E: 'CORALL libraries generate reads in forward orientation';
    'Read 1 reflects the RNA transcript sequence not the cDNA sequence'."""
    assert fragment_strand(MockRead(is_reverse=False), None, r1_sense=True) == "+"
    assert fragment_strand(MockRead(is_reverse=True), None, r1_sense=True) == "-"


def test_fragment_strand_antisense_protocol_inverts():
    """dUTP/QuantSeq-REV chemistries are antisense -- the opposite call."""
    assert fragment_strand(MockRead(is_reverse=False), None, r1_sense=False) == "-"
    assert fragment_strand(MockRead(is_reverse=True), None, r1_sense=False) == "+"


# --------------------------------------------------------------------------
# grouping: the fragment key + the pair_status safety property
# --------------------------------------------------------------------------

def test_fragment_key_proper_pair_uses_both_ends():
    r1 = MockRead(reference_start=1000, is_read1=True)
    r2 = MockRead(reference_start=1300, is_read1=False, is_reverse=True)
    k = fragment_key(r1, r2)
    assert k.pair_status is PairStatus.PROPER_PAIR
    assert k.five_prime_r1 == 1000
    assert k.five_prime_r2 == 1400
    assert k.reference_name == "chrI"


def test_fragment_key_pcr_duplicates_share_a_key():
    """Two copies of one molecule -> identical key -> same bucket."""
    a = fragment_key(MockRead(reference_start=1000), MockRead(reference_start=1300, is_reverse=True, is_read1=False))
    b = fragment_key(MockRead(reference_start=1000), MockRead(reference_start=1300, is_reverse=True, is_read1=False))
    assert a == b


def test_fragment_key_different_fragment_length_separates():
    """Random priming makes the SPAN discriminating: same R1 start, different R2
    end = two distinct molecules, even with the same UMI."""
    a = fragment_key(MockRead(reference_start=1000), MockRead(reference_start=1300, is_reverse=True, is_read1=False))
    b = fragment_key(MockRead(reference_start=1000), MockRead(reference_start=1500, is_reverse=True, is_read1=False))
    assert a != b


def test_fragment_key_degraded_never_collides_with_full_span():
    """THE SAFETY PROPERTY: a singleton key is less specific than a full span, so
    it must never share a bucket with one -- pair_status is in the key to prevent
    a degraded read absorbing molecules the full span would have separated."""
    full = fragment_key(MockRead(reference_start=1000), MockRead(reference_start=1300, is_reverse=True, is_read1=False))
    degraded = fragment_key(MockRead(reference_start=1000), MockRead(is_unmapped=True, is_read1=False))
    assert degraded.pair_status is PairStatus.SINGLETON_R1
    assert full != degraded
    assert degraded.five_prime_r2 is None


def test_fragment_key_r1_unmapped_falls_back_to_r2():
    """Reads spanning novel junctions are exactly where a mate fails to align."""
    k = fragment_key(MockRead(is_unmapped=True), MockRead(reference_start=1300, is_read1=False))
    assert k.pair_status is PairStatus.SINGLETON_R2
    assert k.five_prime_r1 is None
    assert k.five_prime_r2 == 1300


def test_fragment_key_unpaired_library():
    k = fragment_key(MockRead(is_paired=False))
    assert k.pair_status is PairStatus.UNPAIRED


def test_fragment_key_requires_a_mapped_mate():
    with pytest.raises(ValueError):
        fragment_key(MockRead(is_unmapped=True), MockRead(is_unmapped=True, is_read1=False))


def test_fragment_key_spliced_is_distinct_from_unspliced():
    """umi_tools --spliced-is-unique: a spliced and an unspliced read at one
    position are different molecules."""
    unspliced = fragment_key(MockRead(cigartuples=[(_M, 100)]))
    spliced = fragment_key(MockRead(cigartuples=[(_M, 50), (_N, 500), (_M, 50)]))
    assert unspliced.spliced is False and spliced.spliced is True
    assert unspliced != spliced


def test_fragment_key_tolerance_defaults_to_exact():
    """Default MUST be exact: 1 bp apart = different molecules with tolerance 0."""
    a = fragment_key(MockRead(reference_start=1000))
    b = fragment_key(MockRead(reference_start=1001))
    assert a != b


def test_fragment_key_tolerance_bins_are_transitive():
    """With tolerance on, grouping stays transitive (binning, not a fuzzy window):
    equal keys group, unequal keys do not -- no a~b, b~c, not-a~c pathology."""
    keys = [fragment_key(MockRead(reference_start=p), tolerance=2) for p in (1000, 1001, 1002)]
    assert keys[0] == keys[1] == keys[2]


# --------------------------------------------------------------------------
# clustering: re-export + the ED-1 fast path (the setting we default to)
# --------------------------------------------------------------------------

def test_clustering_reexport_is_importable_from_umi_package():
    assert callable(umi_components_directional) and callable(canonical_umi)


def test_directional_absorbs_single_error_child_at_ed1():
    """A sequencing-error UMI is RARER than its parent -> absorbed (2n-1 rule).
    Exercises the O(n*L) hash-neighbour fast path, which is gated on max_edit==1."""
    umis = ["ACGTACGTACGT"] * 10 + ["ACGTACGTACGA"]  # 1 substitution, count 1
    clusters = umi_components_directional(umis, max_edit=1)
    assert len(clusters) == 1
    assert len(clusters[0]) == 11


def test_directional_keeps_equal_count_neighbours_apart():
    """Two real molecules have NO count relationship -> the 2x gradient rejects the
    merge even though they are 1 edit apart. This is why ED-1 is safe."""
    umis = ["ACGTACGTACGT"] * 5 + ["ACGTACGTACGA"] * 5
    clusters = umi_components_directional(umis, max_edit=1)
    assert len(clusters) == 2


def test_directional_distant_umis_never_merge():
    umis = ["AAAAAAAAAAAA"] * 5 + ["TTTTTTTTTTTT"] * 5
    assert len(umi_components_directional(umis, max_edit=1)) == 2


def test_directional_empty_and_singleton():
    assert umi_components_directional([], max_edit=1) == []
    assert umi_components_directional(["ACGTACGTACGT"], max_edit=1) == [[0]]


def test_canonical_umi_picks_most_frequent():
    assert canonical_umi(["ACGT"] * 3 + ["ACGA"]) == "ACGT"
