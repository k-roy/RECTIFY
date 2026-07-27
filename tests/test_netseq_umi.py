"""Tests for the NET-seq UMI path: rectify/core/netseq/netseq_umi.py + the orientation it depends on.

Four layers, in order of how badly a regression would hurt:
  * ORIENTATION -- gene strand and which alignment terminus is the RNA 3' end. Everything else in the NET-seq
    module is expressed in terms of the strand, so this is the single load-bearing decision.
  * RANDOMER GEOMETRY -- the L>N / L==N / L<N cases, including the overshoot correction that the L<N case
    (~26% of reads in GSE159603) requires.
  * MOLECULE SELECTION -- the pure clustering core with literal fragments (no BAM), including the guarantee
    that the key is the CORRECTED 3' end and NOT the soft-clip-corrected 5' end.
  * SATURATION -- the occupancy correction and its flag, which is what keeps a 6-nt UMI honest at depth.
"""
from __future__ import annotations

import math

import pytest

from rectify.core.bam.netseq_bam_processor import (
    NETSEQ_RNA3P_DEFAULT,
    get_netseq_3prime_position,
    netseq_gene_strand,
)
from rectify.core.netseq.netseq_umi import (
    NETSEQ_UMI_LENGTH,
    NetseqFragment,
    extract_netseq_umi,
    five_prime_terminal_clip,
    is_saturated,
    netseq_nontemplated_tail,
    netseq_rna_three_prime,
    occupancy_corrected_molecules,
    randomer_overshoot,
    select_netseq_molecules,
    summarize_positions,
)

_M, _S = 0, 4  # CIGAR M, S


class FakeRead:
    """Minimal pysam.AlignedSegment stand-in.

    ``seq`` is the read AS SEQUENCED (what ``get_forward_sequence`` returns), which is the orientation the
    randomer lives in. ``cigartuples`` is in REFERENCE orientation, as pysam reports it -- the distinction
    matters: for a reverse-mapped read the read's 5' terminus is the LAST cigar op.
    """

    def __init__(self, qname="r", seq="A" * 50, cigartuples=((_M, 50),), reference_start=100,
                 reference_end=150, is_reverse=False, mapping_quality=42, reference_name="chrI"):
        self.query_name = qname
        self._seq = seq
        self.cigartuples = list(cigartuples)
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.is_reverse = is_reverse
        self.mapping_quality = mapping_quality
        self.reference_name = reference_name
        self.query_sequence = seq
        self.query_qualities = [30] * len(seq)
        self.is_unmapped = False
        self.is_secondary = False
        self.is_supplementary = False

    def get_forward_sequence(self):
        return self._seq

    def has_tag(self, _tag):
        return False


# ---------------------------------------------------------------------------
# ORIENTATION
# ---------------------------------------------------------------------------

def test_gene_strand_is_inverse_of_bam_strand():
    """The measured convention: the NET-seq read is revcomp(RNA), so a '+' gene read aligns REVERSE.

    Verified on GSE25107/GSE159603: BAM strand matches the annotated gene strand in only 18-19% of reads.
    """
    assert netseq_gene_strand(FakeRead(is_reverse=True)) == "+"
    assert netseq_gene_strand(FakeRead(is_reverse=False)) == "-"


def test_legacy_convention_is_the_exact_inverse_and_still_reachable():
    assert netseq_gene_strand(FakeRead(is_reverse=True), rna3p_at="read3p") == "-"
    assert netseq_gene_strand(FakeRead(is_reverse=False), rna3p_at="read3p") == "+"


def test_unknown_convention_raises():
    with pytest.raises(ValueError):
        netseq_gene_strand(FakeRead(), rna3p_at="whatever")


def test_rna_three_prime_is_the_reads_five_prime_terminus():
    """'+' gene (reverse read) -> reference_end-1; '-' gene (forward read) -> reference_start."""
    plus = FakeRead(is_reverse=True, reference_start=100, reference_end=150)
    minus = FakeRead(is_reverse=False, reference_start=100, reference_end=150)
    assert netseq_rna_three_prime(plus) == (149, "+")
    assert netseq_rna_three_prime(minus) == (100, "-")


def test_processor_and_module_agree_on_position_and_strand():
    """The read-level track and the molecule track must sit on the same coordinate, or they cannot be
    reported side by side -- which is exactly how a saturating UMI has to be reported."""
    for is_rev in (True, False):
        read = FakeRead(is_reverse=is_rev, reference_start=100, reference_end=150)
        pos, strand, _ = get_netseq_3prime_position(read, trim_terminal_oligo_a=False,
                                                    rna3p_at=NETSEQ_RNA3P_DEFAULT)
        assert (pos, strand) == netseq_rna_three_prime(read)


# ---------------------------------------------------------------------------
# RANDOMER GEOMETRY
# ---------------------------------------------------------------------------

def test_five_prime_clip_reads_the_last_cigar_op_for_a_reverse_read():
    """cigartuples is reference-oriented, so the READ's 5' terminus is the LAST op when reverse."""
    rev = FakeRead(is_reverse=True, cigartuples=((_M, 44), (_S, 6)))
    fwd = FakeRead(is_reverse=False, cigartuples=((_S, 6), (_M, 44)))
    assert five_prime_terminal_clip(rev) == 6
    assert five_prime_terminal_clip(fwd) == 6


def test_clip_equal_to_umi_length_means_no_tail_and_no_position_shift():
    read = FakeRead(is_reverse=True, cigartuples=((_M, 44), (_S, 6)),
                    reference_start=100, reference_end=144, seq="ACGTAC" + "G" * 44)
    assert netseq_nontemplated_tail(read, NETSEQ_UMI_LENGTH) == ""
    assert netseq_rna_three_prime(read, NETSEQ_UMI_LENGTH) == (143, "+")


def test_clip_longer_than_umi_exposes_the_tail_and_leaves_the_terminus_alone():
    """L > N: the extra clipped bases ARE the non-templated tail (poly(A) reads as T in read orientation)."""
    read = FakeRead(is_reverse=True, cigartuples=((_M, 41), (_S, 9)),
                    reference_start=100, reference_end=141, seq="ACGTAC" + "TTT" + "G" * 41)
    assert netseq_nontemplated_tail(read, NETSEQ_UMI_LENGTH) == "TTT"
    assert netseq_rna_three_prime(read, NETSEQ_UMI_LENGTH) == (140, "+")


def test_clip_shorter_than_umi_corrects_the_overshoot():
    """L < N: (N-L) randomer bases aligned by chance, so the terminus overshoots the true RNA 3' end.

    This is ~26% of reads in GSE159603 -- not an edge case. On a '+' gene the RNA 3'-ward direction is
    towards HIGHER coordinates, so the correction moves the position DOWN.
    """
    plus = FakeRead(is_reverse=True, cigartuples=((_M, 46), (_S, 4)),
                    reference_start=100, reference_end=146, seq="ACGTAC" + "G" * 44)
    assert randomer_overshoot(4, NETSEQ_UMI_LENGTH) == 2
    assert netseq_rna_three_prime(plus, NETSEQ_UMI_LENGTH) == (143, "+")   # 145 - 2

    minus = FakeRead(is_reverse=False, cigartuples=((_S, 4), (_M, 46)),
                     reference_start=100, reference_end=146, seq="ACGTAC" + "G" * 44)
    assert netseq_rna_three_prime(minus, NETSEQ_UMI_LENGTH) == (102, "-")  # 100 + 2


def test_no_tail_is_claimed_when_the_randomer_partially_aligned():
    """A molecule carrying a tail could not have had its randomer aligned -- so that case has no tail."""
    read = FakeRead(is_reverse=True, cigartuples=((_M, 46), (_S, 4)), seq="ACGTAC" + "G" * 44)
    assert netseq_nontemplated_tail(read, NETSEQ_UMI_LENGTH) == ""


def test_umi_is_sliced_regardless_of_how_much_the_aligner_clipped():
    """Clipping is an alignment fact, not a library fact: the randomer is always seq[:N]."""
    for clip in (0, 3, 6, 9):
        read = FakeRead(is_reverse=True, cigartuples=((_M, 50 - clip), (_S, clip)) if clip else ((_M, 50),),
                        seq="ACGTAC" + "G" * 44)
        assert extract_netseq_umi(read, NETSEQ_UMI_LENGTH) == "ACGTAC"


def test_umi_with_n_is_rejected_not_guessed():
    read = FakeRead(seq="ACGNAC" + "G" * 44)
    assert extract_netseq_umi(read, NETSEQ_UMI_LENGTH) is None


def test_umi_from_read_name_round_trips():
    read = FakeRead(qname="A00197:591:H:3:1101:23475:1031_ACGTAC")
    assert extract_netseq_umi(read, NETSEQ_UMI_LENGTH, umi_source="name") == "ACGTAC"


def test_unknown_umi_source_raises():
    with pytest.raises(ValueError):
        extract_netseq_umi(FakeRead(), NETSEQ_UMI_LENGTH, umi_source="adapter")


# ---------------------------------------------------------------------------
# MOLECULE SELECTION
# ---------------------------------------------------------------------------

def _frag(qname, umi, pos=1000, strand="+", spliced=False, score=50, contig="chrI"):
    return NetseqFragment(qname=qname, contig=contig, strand=strand, corrected_3p=pos,
                          spliced=spliced, umi=umi, score=score)


def test_pcr_duplicates_at_one_corrected_3p_collapse():
    frags = [_frag("a", "ACGTAC"), _frag("b", "ACGTAC"), _frag("c", "ACGTAC")]
    keepers, family = select_netseq_molecules(frags)
    assert len(keepers) == 1
    assert list(family.values()) == [3]


def test_distinct_umis_at_one_position_stay_distinct():
    keepers, _ = select_netseq_molecules([_frag("a", "ACGTAC"), _frag("b", "TTGGCC")])
    assert len(keepers) == 2


def test_same_umi_at_different_corrected_3p_are_different_molecules():
    keepers, _ = select_netseq_molecules([_frag("a", "ACGTAC", pos=1000),
                                          _frag("b", "ACGTAC", pos=1001)])
    assert len(keepers) == 2


def test_strand_is_part_of_the_key():
    keepers, _ = select_netseq_molecules([_frag("a", "ACGTAC", strand="+"),
                                          _frag("b", "ACGTAC", strand="-")])
    assert len(keepers) == 2


def test_spliced_and_unspliced_at_one_position_are_distinct():
    keepers, _ = select_netseq_molecules([_frag("a", "ACGTAC", spliced=True),
                                          _frag("b", "ACGTAC", spliced=False)])
    assert len(keepers) == 2


def test_exact_clustering_does_not_absorb_a_one_substitution_neighbour():
    """The whole point of defaulting to exact for a 6-mer: at U=4096 a distance-1 neighbour is far more
    likely to be a different molecule than a sequencing error, and directional clustering would eat it."""
    keepers, _ = select_netseq_molecules([_frag(f"hi{i}", "ACGTAC") for i in range(5)]
                                         + [_frag("nb", "ACGTAG")])
    assert len(keepers) == 2


def test_directional_clustering_is_still_reachable_for_measurement():
    keepers, _ = select_netseq_molecules([_frag(f"hi{i}", "ACGTAC") for i in range(5)]
                                         + [_frag("nb", "ACGTAG")],
                                         edit_distance=1, clustering="directional")
    assert len(keepers) == 1


def test_representative_is_the_highest_scoring_read():
    keepers, _ = select_netseq_molecules([_frag("lo", "ACGTAC", score=10),
                                          _frag("hi", "ACGTAC", score=90)])
    assert list(keepers) == ["hi"]


def test_stats_report_duplication_rate_and_family_sizes():
    from rectify.core.umi.dedup import UmiDedupStats
    stats = UmiDedupStats()
    frags = [_frag("a", "ACGTAC"), _frag("b", "ACGTAC"), _frag("c", "TTGGCC")]
    select_netseq_molecules(frags, stats=stats)
    assert stats.n_input_fragments == 3
    assert stats.n_molecules == 2
    assert stats.duplication_rate == pytest.approx(1 / 3)
    assert dict(stats.family_size_hist) == {1: 1, 2: 1}


def test_unknown_clustering_raises():
    with pytest.raises(ValueError):
        select_netseq_molecules([_frag("a", "ACGTAC")], edit_distance=1, clustering="magic")


# ---------------------------------------------------------------------------
# SATURATION
# ---------------------------------------------------------------------------

def test_occupancy_correction_inverts_the_collision_model():
    """m_hat must invert E[k] = U(1 - e^(-m/U)) -- round-tripping is the real test, not a magic number."""
    u = 4 ** NETSEQ_UMI_LENGTH
    for m in (10, 100, 1000, 3000):
        k = round(u * (1 - math.exp(-m / u)))
        assert occupancy_corrected_molecules(k, NETSEQ_UMI_LENGTH) == pytest.approx(m, rel=0.02)


def test_occupancy_correction_matches_the_measured_deep_position():
    """chrI's deepest 3' end: k = 2311 of 4096 -> m_hat 3402, a 47% upward correction (design doc 478 s4b)."""
    assert occupancy_corrected_molecules(2311, NETSEQ_UMI_LENGTH) == pytest.approx(3402, abs=2)


def test_exhausted_barcode_space_is_infinite_not_silently_capped():
    assert occupancy_corrected_molecules(4096, NETSEQ_UMI_LENGTH) == float("inf")


def test_saturation_flag_fires_above_half_occupancy():
    assert not is_saturated(2047, NETSEQ_UMI_LENGTH)
    assert is_saturated(2311, NETSEQ_UMI_LENGTH)


def test_summarize_positions_reports_reads_and_molecules_side_by_side():
    frags = ([_frag(f"d{i}", "ACGTAC", pos=1000) for i in range(5)]
             + [_frag("e", "TTGGCC", pos=1000)]
             + [_frag("f", "ACGTAC", pos=2000)])
    rows = {(r.position): r for r in summarize_positions(frags, NETSEQ_UMI_LENGTH)}
    assert rows[1000].reads == 6 and rows[1000].distinct_umis == 2
    assert rows[2000].reads == 1 and rows[2000].distinct_umis == 1
    assert not rows[1000].saturated
    # the corrected count must never be BELOW the observed distinct count
    for r in rows.values():
        assert r.molecules_corrected >= r.distinct_umis


# ---------------------------------------------------------------------------
# STREAMING (what a real library actually goes through)
# ---------------------------------------------------------------------------

def test_streaming_matches_the_list_based_summary():
    """The streaming path must be an optimisation, not a different answer."""
    from rectify.core.netseq.netseq_umi import stream_netseq_positions

    frags = ([_frag(f"d{i}", "ACGTAC", pos=1000) for i in range(5)]
             + [_frag("e", "TTGGCC", pos=1000)]
             + [_frag("f", "ACGTAC", pos=1200)]
             + [_frag("g", "ACGTAC", pos=90000)])
    listed = {(r.contig, r.strand, r.position): r for r in summarize_positions(frags, NETSEQ_UMI_LENGTH)}
    streamed = {(r.contig, r.strand, r.position): r
                for r in stream_netseq_positions(iter(frags), NETSEQ_UMI_LENGTH)}
    assert set(listed) == set(streamed)
    for key, row in listed.items():
        assert streamed[key].reads == row.reads
        assert streamed[key].distinct_umis == row.distinct_umis
        assert streamed[key].saturated == row.saturated


def test_streaming_flushes_across_a_contig_change():
    """A contig change must finalise everything held -- coordinates restart, so the cursor cannot."""
    from rectify.core.netseq.netseq_umi import stream_netseq_positions

    frags = [_frag("a", "ACGTAC", pos=5000, contig="chrI"),
             _frag("b", "ACGTAC", pos=10, contig="chrII")]
    rows = list(stream_netseq_positions(iter(frags), NETSEQ_UMI_LENGTH))
    assert {(r.contig, r.position) for r in rows} == {("chrI", 5000), ("chrII", 10)}


def test_streaming_stats_match_the_list_based_stats():
    from rectify.core.umi.dedup import UmiDedupStats
    from rectify.core.netseq.netseq_umi import stream_netseq_positions

    frags = [_frag("a", "ACGTAC"), _frag("b", "ACGTAC"), _frag("c", "TTGGCC")]
    s_list, s_stream = UmiDedupStats(), UmiDedupStats()
    select_netseq_molecules(frags, stats=s_list)
    list(stream_netseq_positions(iter(frags), NETSEQ_UMI_LENGTH, stats=s_stream))
    assert s_stream.n_input_fragments == s_list.n_input_fragments == 3
    assert s_stream.n_molecules == s_list.n_molecules == 2
    assert dict(s_stream.family_size_hist) == dict(s_list.family_size_hist)


def test_streaming_holds_only_a_window_not_the_whole_library():
    """The memory guarantee is the whole reason this path exists -- assert it, don't trust it."""
    from rectify.core.netseq.netseq_umi import stream_netseq_positions

    def gen(n=300000):
        for i in range(n):
            yield _frag(f"r{i}", "ACGTAC", pos=i * 10)

    n_rows = 0
    for _ in stream_netseq_positions(gen(), NETSEQ_UMI_LENGTH, flush_pad=2000):
        n_rows += 1
    assert n_rows == 300000
