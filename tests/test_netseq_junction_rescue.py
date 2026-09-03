"""Tests for the NET-seq donor-side junction rescue and non-templated tail call.

Layers, in order of how badly a regression would hurt:
  * GEOMETRY -- donor / exon1_last / exon2_first on both strands, and the RNA-orientation clip.
    Everything else is expressed in those terms, so an error here is silent and total.
  * RESCUE -- the k = 1..10 sweep Kevin asked for, the randomer remainder rule, the mis-extension
    case, and the exon1_end case that must NOT be moved (it is the real splicing-intermediate signal).
  * TAIL -- the randomer must not be counted as a tail, and the walkback must recover tail A's the
    aligner put on genomic A's (Chanfreau invariant 7).
  * SYMMETRY -- a synthetic read mirrored onto both strands must give mirrored answers. The verified
    minus-strand oligo(A) trim-sign bug (planning 834 Appendix) is exactly what this catches.
"""
from __future__ import annotations

import random

import pysam
import pytest

from rectify.core.bam.netseq_bam_processor import (
    get_netseq_3prime_position,
    process_netseq_read,
    rna_5prime_ward,
)
from rectify.core.netseq.netseq_rescue import (
    JunctionPool,
    NetseqCorrectionSummary,
    allowed_remainders,
    call_tail,
    load_junction_tsv,
    longest_common_prefix,
    make_junction,
    rescue_read,
    revcomp,
    rna_clip,
)

_M, _S = 0, 4

# ---------------------------------------------------------------------------------------------
# A synthetic chromosome: exon1 [0,100)  intron [100,160)  exon2 [160,260)
# ---------------------------------------------------------------------------------------------
INTRON_START, INTRON_END = 100, 160
CHROM = "chrT"


def _random_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _build_genome():
    exon1 = _random_seq(100, 1)
    intron = "GT" + _random_seq(56, 2) + "AG"
    exon2 = _random_seq(100, 3)
    return {CHROM: exon1 + intron + exon2}


GENOME = _build_genome()
SEQ = GENOME[CHROM]
EXON2_PLUS = SEQ[INTRON_END:]                    # gene-strand exon-2 for a '+' gene
EXON2_MINUS = revcomp(SEQ[:INTRON_START])        # gene-strand exon-2 for a '-' gene, from pos 99 down

_HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": len(SEQ)}]}
)


def make_read(strand, aligned_seq, reference_start, clip_rna="", qname="r", md=None):
    """One synthetic NET-seq read.

    ``strand`` is the GENE strand. Under ``--rna3p-at read5p`` a '+' gene read maps REVERSE, so its
    read-5' terminus (= the RNA 3' end + the clip) is the RIGHT end of the reference-forward record;
    a '-' gene read maps forward and it is the LEFT end, reverse-complemented.
    """
    read = pysam.AlignedSegment(_HEADER)
    read.query_name = qname
    read.reference_id = 0
    read.reference_start = reference_start
    read.mapping_quality = 255
    read.flag = 16 if strand == '+' else 0
    if strand == '+':
        read.query_sequence = aligned_seq + clip_rna
        read.cigartuples = [(_M, len(aligned_seq))] + ([(_S, len(clip_rna))] if clip_rna else [])
    else:
        clip_genomic = revcomp(clip_rna)
        read.query_sequence = clip_genomic + aligned_seq
        read.cigartuples = ([(_S, len(clip_genomic))] if clip_rna else []) + [(_M, len(aligned_seq))]
    read.query_qualities = [30] * len(read.query_sequence)
    if md is not None:
        read.set_tag("MD", md)
    return read


def non_matching_run(next_expected_base, length):
    """A run of one base guaranteed NOT to extend an exon-2 match past ``next_expected_base``.

    Test fixtures must not depend on the random genome happening to differ -- a randomer that
    accidentally continues the exon-2 sequence changes k and makes the test assert the wrong thing.
    """
    base = 'T' if next_expected_base != 'T' else 'G'
    return base * length


def plus_pool():
    return JunctionPool([make_junction(CHROM, INTRON_START, INTRON_END, '+')])


def minus_pool():
    return JunctionPool([make_junction(CHROM, INTRON_START, INTRON_END, '-')])


# ---------------------------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------------------------

def test_junction_coordinates_are_gene_strand_resolved():
    plus = make_junction(CHROM, 100, 160, '+')
    assert (plus.donor, plus.exon1_last, plus.exon2_first) == (100, 99, 160)
    minus = make_junction(CHROM, 100, 160, '-')
    assert (minus.donor, minus.exon1_last, minus.exon2_first) == (159, 160, 99)


def test_rna_clip_is_the_read5p_clip_in_rna_orientation():
    plus = make_read('+', SEQ[60:100], 60, clip_rna="ACGTA")
    minus = make_read('-', SEQ[160:200], 160, clip_rna="ACGTA")
    assert rna_clip(plus, '+') == "ACGTA"
    assert rna_clip(minus, '-') == "ACGTA"
    # ...and it really is the opposite physical end of the record on the two strands
    assert plus.query_sequence.endswith("ACGTA")
    assert minus.query_sequence.startswith(revcomp("ACGTA"))


def test_pool_lookup_finds_the_donor_and_reports_n_intronic():
    pool = plus_pool()
    assert pool.lookup(CHROM, '+', 99, 10) == (0, pool.lookup(CHROM, '+', 99, 10)[1])
    assert pool.lookup(CHROM, '+', 99, 10)[0] == 0     # aligned end exactly at exon1_last
    assert pool.lookup(CHROM, '+', 104, 10)[0] == 5    # 5 nt mis-extended into the intron
    assert pool.lookup(CHROM, '+', 110, 10) is None    # 11 nt: out of the window
    assert pool.lookup(CHROM, '+', 99, 10) is not None
    assert pool.lookup(CHROM, '-', 99, 10) is None     # wrong strand
    mpool = minus_pool()
    assert mpool.lookup(CHROM, '-', 160, 10)[0] == 0
    assert mpool.lookup(CHROM, '-', 155, 10)[0] == 5


def test_longest_common_prefix():
    assert longest_common_prefix("ACGT", "ACGA") == 3
    assert longest_common_prefix("", "ACGT") == 0
    assert longest_common_prefix("ACGT", "ACGT") == 4


# ---------------------------------------------------------------------------------------------
# RESCUE -- the k = 1..10 sweep
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("k", list(range(1, 11)))
def test_clean_overhang_of_k_nt_is_rescued_into_exon2_plus_strand(k):
    """A soft-clipped exon-2 overhang of 1-10 nt moves the 3' end k-1 nt into exon 2."""
    read = make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:k])
    call = rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, umi_length=0)
    assert call.status == "spliced_rescued"
    assert call.k == k and call.r == 0 and call.n_intronic == 0
    assert call.position == INTRON_END + k - 1


@pytest.mark.parametrize("k", list(range(1, 11)))
def test_clean_overhang_of_k_nt_is_rescued_into_exon2_minus_strand(k):
    """Mirror image: on a '-' gene the rescued end runs DOWN from exon2_first = intron_start - 1."""
    read = make_read('-', SEQ[160:200], 160, clip_rna=EXON2_MINUS[:k])
    call = rescue_read(read, '-', CHROM, 160, minus_pool(), SEQ, umi_length=0)
    assert call.status == "spliced_rescued"
    assert call.k == k and call.r == 0
    assert call.position == (INTRON_START - 1) - (k - 1)


def test_randomer_remainder_is_tolerated_only_when_a_randomer_is_declared():
    """The 6-nt randomer sits DISTAL to the overhang in the same clip (planning 478 §1 / 829 §4)."""
    randomer = non_matching_run(EXON2_PLUS[5], 6)      # must not extend the exon-2 match
    clip = EXON2_PLUS[:5] + randomer
    read = make_read('+', SEQ[60:100], 60, clip_rna=clip)

    with_umi = rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, umi_length=6)
    assert with_umi.status == "spliced_rescued"
    assert with_umi.k == 5 and with_umi.r == 6
    assert with_umi.position == INTRON_END + 4

    without_umi = rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, umi_length=0)
    assert without_umi.status == "ambiguous"   # k>=1 but a 6-nt remainder is not allowed
    assert without_umi.position is None


def test_allowed_remainders_covers_the_partially_aligned_randomer():
    """478 §1 measured a 6:5:4 clip ratio of 1 : 0.28 : 0.09 -- the randomer's terminal 1-2 nt align
    by chance, so a declared 6-nt randomer presents as 5, 6 or 7 leftover nt."""
    assert allowed_remainders(0) == (0,)
    assert allowed_remainders(6) == (0, 5, 6, 7)


def test_mis_extension_into_the_intron_is_rescued_and_counted():
    """The other aligner failure mode: exon-2 bases forced onto the first intronic bases."""
    n_intronic = 3
    aligned = SEQ[60:100] + EXON2_PLUS[:n_intronic]     # read bases are exon 2, placed on the intron
    read = make_read('+', aligned, 60, clip_rna=EXON2_PLUS[n_intronic:5])
    p = 100 + n_intronic - 1                            # reference_end - 1
    assert read.reference_end - 1 == p
    call = rescue_read(read, '+', CHROM, p, plus_pool(), SEQ, umi_length=0)
    assert call.status == "spliced_rescued"
    assert call.n_intronic == n_intronic
    assert call.k == 5 and call.position == INTRON_END + 4


def test_exon1_end_is_kept_where_it_is():
    """A read that stops exactly at the exon-1 3' end with nothing matching exon 2 IS the splicing
    intermediate -- moving it would erase the signal the rescue exists to protect."""
    read = make_read('+', SEQ[60:100], 60, clip_rna="")
    call = rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, umi_length=0)
    assert call.status == "exon1_end"
    assert call.k == 0 and call.position is None


def test_randomer_only_clip_that_does_not_match_exon2_is_an_exon1_end():
    randomer = non_matching_run(EXON2_PLUS[0], 6)
    read = make_read('+', SEQ[60:100], 60, clip_rna=randomer)
    call = rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, umi_length=6)
    assert call.status == "exon1_end"


def test_min_k_raises_the_floor():
    """A read that fails the floor keeps the class its GEOMETRY implies -- raising the floor must
    return chance-matched reads to `exon1_end` (the splicing-intermediate readout), not hide them
    in `ambiguous`."""
    read = make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:1])
    assert rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, min_k=1).status == "spliced_rescued"
    assert rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, min_k=2).status == "exon1_end"


def test_min_k_floor_is_remainder_aware():
    """THE CHANCE CHANNEL IS THE REMAINDER. With r == 0 the whole clip is exon-2 sequence, so k=1 is
    evidence; with r > 0 a randomer is being invoked and a randomer's first base matches exon 2 a
    quarter of the time, so the floor has to be higher."""
    pool, genome = plus_pool(), SEQ
    randomer = non_matching_run(EXON2_PLUS[1], 6)

    clean = make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:1])
    with_rem = make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:1] + randomer)

    # default floors: 1 with no remainder, 4 with one
    c = rescue_read(clean, '+', CHROM, 99, pool, genome, umi_length=6)
    assert (c.status, c.k, c.r) == ("spliced_rescued", 1, 0)
    w = rescue_read(with_rem, '+', CHROM, 99, pool, genome, umi_length=6)
    assert (w.status, w.k, w.r) == ("exon1_end", 1, 6)          # k=1 < 4 -> not moved

    # a 4-nt overhang ahead of the same randomer clears the higher floor
    randomer4 = non_matching_run(EXON2_PLUS[4], 6)
    deep = make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:4] + randomer4)
    d = rescue_read(deep, '+', CHROM, 99, pool, genome, umi_length=6)
    assert (d.status, d.k, d.r) == ("spliced_rescued", 4, 6)

    # setting the two equal restores the old flat rule
    flat = rescue_read(with_rem, '+', CHROM, 99, pool, genome, umi_length=6,
                       min_k_with_remainder=1)
    assert flat.status == "spliced_rescued" and flat.k == 1


def test_decoy_uses_the_same_remainder_aware_rule():
    """`decoy_would_rescue` is the chance floor only if it runs the IDENTICAL acceptance test."""
    randomer = non_matching_run(EXON2_PLUS[1], 6)
    read = make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:1] + randomer)
    loose = rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, umi_length=6,
                        min_k_with_remainder=1)
    strict = rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, umi_length=6,
                         min_k_with_remainder=9)
    # the decoy verdict must track the floor, not be frozen at min_k
    assert strict.decoy_would_rescue is False
    assert loose.decoy_would_rescue == (loose.decoy_k >= 1 and
                                        (len(loose.s_seq) - loose.decoy_k) in (0, 5, 6, 7))


def test_no_pooled_donor_means_no_call():
    read = make_read('+', SEQ[20:60], 20, clip_rna=EXON2_PLUS[:4])
    call = rescue_read(read, '+', CHROM, 59, plus_pool(), SEQ)
    assert call.status == "none" and call.position is None


def test_decoy_acceptor_null_is_computed():
    """The chance-match control travels WITH the call, so a k histogram can be read as evidence."""
    read = make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:6])
    call = rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, decoy_offset=50)
    assert call.k == 6
    assert call.decoy_k == longest_common_prefix(EXON2_PLUS[:6], SEQ[INTRON_END + 50:])


# ---------------------------------------------------------------------------------------------
# TAIL
# ---------------------------------------------------------------------------------------------

def test_plain_a_tail_in_the_clip():
    read = make_read('+', SEQ[60:100], 60, clip_rna="AAAA")
    tail = call_tail(read, '+', 99, SEQ, umi_length=0)
    assert tail.tail_len == 4 and tail.clip_a_run == 4 and tail.walkback == 0
    assert tail.tail_class == "A_only"


def test_randomer_is_stripped_before_the_tail_call():
    """The whole point of the randomer-aware rule: a 6-nt randomer-only clip is NOT a 6-nt tail, and
    the legacy 0.8-A-fraction test over the whole clip cannot express that."""
    read = make_read('+', SEQ[60:100], 60, clip_rna="AAAAAA")
    assert call_tail(read, '+', 99, SEQ, umi_length=6).tail_len == 0
    assert call_tail(read, '+', 99, SEQ, umi_length=0).tail_len == 6


def test_tail_plus_randomer_reports_only_the_tail():
    read = make_read('+', SEQ[60:100], 60, clip_rna="AAAA" + "CGTCGT")
    tail = call_tail(read, '+', 99, SEQ, umi_length=6)
    assert tail.tail_len == 4 and tail.randomer == "CGTCGT" and tail.tail_region == "AAAA"


def test_walkback_recovers_tail_a_bases_the_aligner_put_on_genomic_a():
    """Chanfreau invariant 7: walk past every A, stop at the first non-A."""
    stop = 'C' if SEQ[96] == 'A' else SEQ[96]
    genome = {CHROM: SEQ[:96] + stop + "AAA" + SEQ[100:]}
    aligned = genome[CHROM][60:100]                    # ends ...AAA over genomic AAA
    read = make_read('+', aligned, 60, clip_rna="")
    tail = call_tail(read, '+', 99, genome[CHROM], umi_length=0)
    assert tail.walkback == 3 and tail.tail_len == 3
    assert genome[CHROM][96] != "A"                    # the walk really did stop at a non-A


def test_walkback_is_gated_when_a_non_a_base_bounds_the_tail():
    """If the clip's tail region ends in a non-A, the RNA 3' terminus is bounded there and the
    aligned A's cannot be part of the same tail."""
    stop = 'C' if SEQ[96] == 'A' else SEQ[96]
    genome = {CHROM: SEQ[:96] + stop + "AAA" + SEQ[100:]}
    read = make_read('+', genome[CHROM][60:100], 60, clip_rna="AAC")
    tail = call_tail(read, '+', 99, genome[CHROM], umi_length=0)
    assert tail.clip_a_run == 2 and tail.walkback == 0 and tail.tail_len == 2
    assert tail.tail_class == "A_rich"


def test_walkback_mirrors_on_the_minus_strand():
    """On a '-' gene the RNA A is a genomic T and the walk runs the other way."""
    stop = 'G' if SEQ[163] == 'T' else SEQ[163]
    genome = {CHROM: SEQ[:160] + "TTT" + stop + SEQ[164:]}
    read = make_read('-', genome[CHROM][160:200], 160, clip_rna="")
    tail = call_tail(read, '-', 160, genome[CHROM], umi_length=0)
    assert tail.walkback == 3 and tail.tail_len == 3


def test_tail_class_other_for_a_non_a_clip():
    read = make_read('+', SEQ[60:100], 60, clip_rna="CGCGCG")
    tail = call_tail(read, '+', 99, SEQ, umi_length=0)
    assert tail.tail_len == 0 and tail.tail_class == "other"


# ---------------------------------------------------------------------------------------------
# SYMMETRY -- the minus-strand oligo(A) trim sign (planning 834 Appendix)
# ---------------------------------------------------------------------------------------------

def test_terminal_oligo_a_trim_moves_the_end_5prime_ward_on_both_strands():
    """Reproduces the 834 Appendix construction: three tail A's aligned as MD mismatches.

    '+' gene: reference_end-1 = 119, true 3' end 116.
    '-' gene: reference_start = 100, true 3' end 103 -- the old code returned 97 (2 x n the wrong way).
    """
    minus = make_read('-', "T" * 3 + SEQ[103:120], 100, md="0A0A0A17")
    pos, strand, n = get_netseq_3prime_position(minus, trim_terminal_oligo_a=True)
    assert (strand, n) == ('-', 3)
    assert pos == 103

    plus = make_read('+', SEQ[100:117] + "A" * 3, 100, md="17C0C0C")
    pos, strand, n = get_netseq_3prime_position(plus, trim_terminal_oligo_a=True)
    assert (strand, n) == ('+', 3)
    assert pos == 116


def test_rna_5prime_ward_is_the_single_source_of_the_sign():
    assert rna_5prime_ward('+') == -1
    assert rna_5prime_ward('-') == 1


def test_record_raw_position_survives_the_trim_on_both_strands():
    """three_prime_raw must be the UNCORRECTED terminus. It used to be derived with the buggy sign,
    so the .raw.* bedgraph would have broken on exactly the reads the sign fix touches."""
    minus = make_read('-', "T" * 3 + SEQ[103:120], 100, md="0A0A0A17")
    rec = process_netseq_read(minus, CHROM, detect_tail=False)
    assert rec.strand == '-' and rec.three_prime_raw == 100 and rec.three_prime_corrected == 103

    plus = make_read('+', SEQ[100:117] + "A" * 3, 100, md="17C0C0C")
    rec = process_netseq_read(plus, CHROM, detect_tail=False)
    assert rec.strand == '+' and rec.three_prime_raw == 119 and rec.three_prime_corrected == 116


def test_mirrored_synthetic_read_gives_mirrored_answers_end_to_end():
    """One read, both strands: rescue k, tail length and the signed end shift must mirror exactly."""
    plus_read = make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:4])
    minus_read = make_read('-', SEQ[160:200], 160, clip_rna=EXON2_MINUS[:4])
    p_rec = process_netseq_read(plus_read, CHROM, junction_pool=plus_pool(), genome=GENOME)
    m_rec = process_netseq_read(minus_read, CHROM, junction_pool=minus_pool(), genome=GENOME)
    assert p_rec.rescue_status == m_rec.rescue_status == "spliced_rescued"
    assert p_rec.rescue_k == m_rec.rescue_k == 4
    assert p_rec.three_prime_corrected == INTRON_END + 3
    assert m_rec.three_prime_corrected == (INTRON_START - 1) - 3
    # the shift is the mirror image: +64 on '+', -64 on '-'
    assert p_rec.three_prime_shift == -m_rec.three_prime_shift


# ---------------------------------------------------------------------------------------------
# PLUMBING
# ---------------------------------------------------------------------------------------------

def test_process_netseq_read_records_the_rescue_and_the_tail():
    randomer = non_matching_run(EXON2_PLUS[4], 6)   # k=4 clears the default remainder floor
    read = make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:4] + randomer)
    rec = process_netseq_read(read, CHROM, junction_pool=plus_pool(), genome=GENOME, umi_length=6)
    assert rec.rescue_status == "spliced_rescued" and rec.rescue_k == 4 and rec.rescue_r == 6
    assert rec.three_prime_corrected == INTRON_END + 3
    assert rec.clip_rna == EXON2_PLUS[:4] + randomer


def test_rescue_and_tail_are_no_ops_without_a_pool_or_a_genome():
    read = make_read('+', SEQ[60:100], 60, clip_rna="AAAA")
    rec = process_netseq_read(read, CHROM)
    assert rec.rescue_status == "none"
    assert rec.three_prime_corrected == rec.three_prime_raw == 99
    assert rec.tail_len == 4       # the clip A-run needs no genome; only the walkback does
    assert rec.tail_walkback == 0


def test_summary_counters():
    summary = NetseqCorrectionSummary()
    reads = [
        make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:k], qname=f"k{k}") for k in (1, 2, 5)
    ] + [make_read('+', SEQ[60:100], 60, clip_rna="", qname="e1")]
    for read in reads:
        summary.observe(process_netseq_read(read, CHROM, junction_pool=plus_pool(), genome=GENOME))
    assert summary.reads == 4
    assert summary.rescued == 3
    assert summary.exon1_end == 1
    assert summary.rescued_by_k == {1: 1, 2: 1, 5: 1}
    assert summary.rescued_by_k_clean == {1: 1, 2: 1, 5: 1}
    assert summary.end_moved == 3
    payload = summary.as_dict()
    assert payload['rescued_by_k'] == {'1': 1, '2': 1, '5': 1}


def test_junction_tsv_round_trip(tmp_path):
    """donor = FIRST intronic base, acceptor = LAST intronic base, both gene-strand."""
    path = tmp_path / "pool.tsv"
    path.write_text(
        "chrom\tdonor\tacceptor\tstrand\n"
        f"{CHROM}\t{INTRON_START}\t{INTRON_END - 1}\t+\n"
        f"{CHROM}\t{INTRON_END - 1}\t{INTRON_START}\t-\n"
    )
    junctions = load_junction_tsv(path)
    assert len(junctions) == 2
    plus = next(j for j in junctions if j.strand == '+')
    minus = next(j for j in junctions if j.strand == '-')
    assert (plus.intron_start, plus.intron_end) == (INTRON_START, INTRON_END)
    assert (minus.intron_start, minus.intron_end) == (INTRON_START, INTRON_END)
    assert plus.exon2_first == INTRON_END and minus.exon2_first == INTRON_START - 1


def test_junction_tsv_rejects_a_bad_header(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("chrom\tstart\tend\tstrand\nchrT\t1\t2\t+\n")
    with pytest.raises(ValueError, match="chrom/donor/acceptor/strand"):
        load_junction_tsv(path)


def test_walkback_can_be_gated_on_clip_evidence():
    """A nascent-RNA library's 3' ends mostly carry NO tail, so an unconditional walkback moves
    ends on no evidence -- measured on wt_rep3 chrI+chrII: 41,711 of 42,644 walkbacks had an empty
    or randomer-only clip, and at RPL32 (exon 1 ends ...AAAA) 24 of 33 reads on the exon-1 3' end
    were walked 4 nt off it. The gate is opt-in; the default stays invariant 7."""
    stop = 'C' if SEQ[96] == 'A' else SEQ[96]
    genome = {CHROM: SEQ[:96] + stop + "AAA" + SEQ[100:]}
    no_clip = make_read('+', genome[CHROM][60:100], 60, clip_rna="")
    assert call_tail(no_clip, '+', 99, genome[CHROM]).walkback == 3
    assert call_tail(no_clip, '+', 99, genome[CHROM], require_clip_evidence=True).walkback == 0
    # a read WITH a non-templated A next to the alignment still walks back
    with_clip = make_read('+', genome[CHROM][60:100], 60, clip_rna="AA")
    gated = call_tail(with_clip, '+', 99, genome[CHROM], require_clip_evidence=True)
    assert gated.clip_a_run == 2 and gated.walkback == 3 and gated.tail_len == 5


def test_molecule_track_uses_the_same_corrected_position_as_the_read_track(tmp_path):
    """planning 834 defect (c): the molecule track keyed on a DIFFERENT position from the read
    track. Now that the read track keys on trim + walkback + rescue, `iter_netseq_fragments` has to
    be handed the same inputs, or every walked-back / rescued read lands on its own coordinate."""
    from rectify.core.netseq.netseq_umi import iter_netseq_fragments

    stop = 'C' if SEQ[96] == 'A' else SEQ[96]
    genome = {CHROM: SEQ[:96] + stop + "AAA" + SEQ[100:]}
    # one walked-back read and one junction-rescued read, both with a usable 6-nt UMI
    walked = make_read('+', genome[CHROM][60:100], 60, clip_rna="ACGTAC", qname="wb")
    rescued = make_read('+', SEQ[60:100], 60, clip_rna=EXON2_PLUS[:4] + non_matching_run(EXON2_PLUS[4], 6),
                        qname="rescue")
    bam_path = tmp_path / "m.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=_HEADER) as out:
        for r in (walked, rescued):
            out.write(r)

    read_positions = {
        process_netseq_read(r, CHROM, junction_pool=plus_pool(), genome=genome, umi_length=6).three_prime_corrected
        for r in (walked, rescued)
    }
    frags = list(iter_netseq_fragments(
        str(bam_path), umi_length=6, genome=genome, junction_pool=plus_pool(),
    ))
    assert len(frags) == 2
    assert {f.corrected_3p for f in frags} == read_positions


def test_walkback_gate_is_the_netseq_pipeline_default():
    """`rectify netseq`'s default is now the GATED walkback (nascent 3' ends carry no tail);
    `call_tail` itself stays neutral (invariant 7) because it is a protocol-agnostic primitive."""
    import inspect

    from rectify.core.bam.netseq_bam_processor import process_netseq_bam, process_netseq_read

    for fn in (process_netseq_read, process_netseq_bam):
        assert inspect.signature(fn).parameters['walkback_requires_clip_a'].default is True
    stop = 'C' if SEQ[96] == 'A' else SEQ[96]
    genome = {CHROM: SEQ[:96] + stop + "AAA" + SEQ[100:]}
    read = make_read('+', genome[CHROM][60:100], 60, clip_rna="")
    assert process_netseq_read(read, CHROM, genome=genome).tail_walkback == 0
    assert process_netseq_read(read, CHROM, genome=genome,
                               walkback_requires_clip_a=False).tail_walkback == 3


# ---------------------------------------------------------------------------------------------
# CLASSIFICATION ORDER -- the floor is tested BEFORE the remainder
# ---------------------------------------------------------------------------------------------

def test_randomer_clip_matching_exon2_by_chance_is_an_exon1_end_not_ambiguous():
    """L=6 randomer clip, first 2 bases match exon 2 by chance -> k=2, r=4.

    4 is not an allowed remainder, so a remainder-FIRST test called this `ambiguous` and removed it
    from the splicing-intermediate tally. Measured on the full wt_rep3 library: 66 reads at k=2 and
    12 at k=3. Below the floor the match is not evidence of anything, so the read keeps the class
    its geometry implies.
    """
    clip = EXON2_PLUS[:2] + non_matching_run(EXON2_PLUS[2], 4)     # 6 nt total, k=2, r=4
    read = make_read('+', SEQ[60:100], 60, clip_rna=clip)
    call = rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, umi_length=6)
    assert (call.k, call.r) == (2, 4)
    assert call.r not in allowed_remainders(6)        # the remainder really is disallowed
    assert call.status == "exon1_end"                 # ...and it is STILL an exon-1 end
    assert call.position is None


def test_same_read_inside_the_intron_is_an_intronic_end():
    n_intronic = 3
    aligned = SEQ[60:100] + EXON2_PLUS[:n_intronic]
    clip = non_matching_run(EXON2_PLUS[n_intronic], 3)
    read = make_read('+', aligned, 60, clip_rna=clip)
    call = rescue_read(read, '+', CHROM, 102, plus_pool(), SEQ, umi_length=6)
    assert call.k == 3 and call.r == 3 and call.r not in allowed_remainders(6)
    assert call.status == "intronic_end"


def test_disallowed_remainder_ABOVE_the_floor_is_still_ambiguous():
    """`ambiguous` survives for the case it was meant for: the match clears the floor but the
    leftover is not a length any randomer can explain."""
    clip = EXON2_PLUS[:5] + non_matching_run(EXON2_PLUS[5], 3)     # k=5 >= 4, r=3 not allowed
    read = make_read('+', SEQ[60:100], 60, clip_rna=clip)
    call = rescue_read(read, '+', CHROM, 99, plus_pool(), SEQ, umi_length=6)
    assert (call.k, call.r) == (5, 3)
    assert call.status == "ambiguous" and call.position is None


# ---------------------------------------------------------------------------------------------
# POOL FILTER -- tRNA introns are not spliceosomal
# ---------------------------------------------------------------------------------------------

_SYNTHETIC_GFF = """\
##gff-version 3
chrT\tSGD\tgene\t101\t900\t.\t+\t.\tID=YAA001W;Name=MRNA1
chrT\tSGD\tmRNA\t101\t900\t.\t+\t.\tID=YAA001W_mRNA;Parent=YAA001W
chrT\tSGD\tintron\t201\t300\t.\t+\t.\tParent=YAA001W_mRNA;Name=YAA001W_intron
chrT\tSGD\ttRNA_gene\t1001\t1105\t.\t+\t.\tID=YNCT0001W;Name=YNCT0001W
chrT\tSGD\ttRNA\t1001\t1105\t.\t+\t.\tID=YNCT0001W_tRNA;Parent=YNCT0001W
chrT\tSGD\tintron\t1037\t1069\t.\t+\t.\tParent=YNCT0001W_tRNA;Name=YNCT0001W_intron
chrT\tSGD\tncRNA_gene\t2001\t2400\t.\t-\t.\tID=snR100;Name=snR100
chrT\tSGD\tsnoRNA\t2001\t2400\t.\t-\t.\tID=snR100_snoRNA;Parent=snR100
chrT\tSGD\tintron\t2101\t2200\t.\t-\t.\tParent=snR100_snoRNA;Name=snR100_intron
"""


def _write_gff(tmp_path, gzipped=False):
    if gzipped:
        import gzip
        path = tmp_path / "ann.gff.gz"
        with gzip.open(path, "wt") as fh:
            fh.write(_SYNTHETIC_GFF)
    else:
        path = tmp_path / "ann.gff"
        path.write_text(_SYNTHETIC_GFF)
    return path


def test_parse_gff_introns_resolves_the_parent_feature_type(tmp_path):
    from rectify.core.netseq.netseq_rescue import parse_gff_introns

    rows = parse_gff_introns(_write_gff(tmp_path))
    types = {r[4] for r in rows}
    assert types == {"mRNA", "tRNA", "snoRNA"}
    trna = next(r for r in rows if r[4] == "tRNA")
    assert trna[1:4] == (1036, 1069, '+') and trna[5] == "YNCT0001W_tRNA"


def test_trna_introns_are_dropped_from_the_pool_by_default(tmp_path):
    """A tRNA intron is excised by the tRNA endonuclease, not the spliceosome, at a Pol III locus
    this command already drops from its tracks -- so it can only manufacture rescues. Measured on
    wt_rep3: the single locus YNCO0031W_tRNA collected 46."""
    gff = _write_gff(tmp_path)
    default = JunctionPool.from_annotation(gff)
    assert len(default) == 2
    assert default.by_parent_type == {"mRNA": 1, "snoRNA": 1}
    assert default.dropped_by_parent_type == {"tRNA": 1}
    assert default.lookup(CHROM, '+', 1036, 10) is None      # the tRNA donor is not in the pool

    kept = JunctionPool.from_annotation(gff, include_trna=True)
    assert len(kept) == 3 and kept.by_parent_type["tRNA"] == 1
    assert kept.dropped_by_parent_type == {}


def test_pol2_introns_of_every_class_are_kept(tmp_path):
    pool = JunctionPool.from_annotation(_write_gff(tmp_path, gzipped=True))
    assert pool.lookup(CHROM, '+', 200, 10) is not None      # mRNA exon1_last = 200
    assert pool.lookup(CHROM, '-', 2200, 10) is not None     # snoRNA, minus strand


def test_is_trna_parent_matches_the_type_or_the_id_suffix():
    from rectify.core.netseq.netseq_rescue import is_trna_parent

    assert is_trna_parent("tRNA", "YNCO0031W_tRNA")
    assert is_trna_parent("tRNA_gene", "YNCO0031W")
    assert is_trna_parent("unknown", "YNCO0031W_tRNA")       # annotation omits the tRNA line
    assert not is_trna_parent("mRNA", "YAL003W_id001")
    assert not is_trna_parent("snoRNA", "snR100_snoRNA")


# ---------------------------------------------------------------------------------------------
# ORGANELLAR INTRONS -- group I/II self-splicing, not spliceosomal (planning 856, 2026-09-03)
# ---------------------------------------------------------------------------------------------

def _gff_with_a_mito_intron(tmp_path):
    """SGD-shaped GFF: one nuclear intron, one mito intron, one tRNA intron. All parents differ."""
    p = tmp_path / 'mix.gff'
    p.write_text(
        "##gff-version 3\n"
        "chrI\tSGD\tmRNA\t1\t1000\t.\t+\t.\tID=YAL001C_mRNA\n"
        "chrI\tSGD\tintron\t101\t160\t.\t+\t.\tID=i_nuc;Parent=YAL001C_mRNA\n"
        "chrmt\tSGD\tmRNA\t1\t9000\t.\t+\t.\tID=Q0045_mRNA\n"
        "chrmt\tSGD\tintron\t501\t1500\t.\t+\t.\tID=aI1;Parent=Q0045_mRNA\n"
        "chrIII\tSGD\ttRNA\t1\t100\t.\t+\t.\tID=YNCC0001W_tRNA\n"
        "chrIII\tSGD\tintron\t40\t60\t.\t+\t.\tID=i_trna;Parent=YNCC0001W_tRNA\n"
    )
    return p


def test_mito_introns_are_dropped_from_the_pool_by_default(tmp_path):
    """🔴 The wt_rep3 finding: mito group I/II introns have parent type `mRNA`, so the tRNA filter
    misses them, and 94 of 580 rescues were fabricated on chrMito."""
    from rectify.core.netseq.netseq_rescue import is_organellar_contig

    pool = JunctionPool.from_annotation(_gff_with_a_mito_intron(tmp_path))
    contigs = {c for c, _s in pool._by_exon1}
    assert contigs == {'chrI'}, f"pool should hold the nuclear intron only, got {contigs}"
    assert any('organellar' in k for k in pool.dropped_by_parent_type), pool.dropped_by_parent_type
    assert is_organellar_contig('chrMito') and is_organellar_contig('chrmt')
    assert is_organellar_contig('chrM') and is_organellar_contig('MT')
    assert not is_organellar_contig('chrI') and not is_organellar_contig('chrXIII')


def test_pool_include_organellar_restores_them(tmp_path):
    pool = JunctionPool.from_annotation(_gff_with_a_mito_intron(tmp_path),
                                        include_organellar=True)
    assert {c for c, _s in pool._by_exon1} == {'chrI', 'chrMito'}


def test_trna_and_organellar_filters_are_independent(tmp_path):
    gff = _gff_with_a_mito_intron(tmp_path)
    both_on = JunctionPool.from_annotation(gff, include_trna=True, include_organellar=True)
    assert {c for c, _s in both_on._by_exon1} == {'chrI', 'chrIII', 'chrMito'}
    trna_only = JunctionPool.from_annotation(gff, include_trna=True)
    assert {c for c, _s in trna_only._by_exon1} == {'chrI', 'chrIII'}


def test_exon_derived_fallback_also_drops_organellar(tmp_path):
    """The exon-only-GTF path must not reintroduce them."""
    p = tmp_path / 'exons.gff'
    p.write_text(
        "##gff-version 3\n"
        "chrI\tX\tmRNA\t1\t1000\t.\t+\t.\tID=m1\n"
        "chrI\tX\texon\t1\t100\t.\t+\t.\tParent=m1\n"
        "chrI\tX\texon\t161\t400\t.\t+\t.\tParent=m1\n"
        "chrmt\tX\tmRNA\t1\t9000\t.\t+\t.\tID=m2\n"
        "chrmt\tX\texon\t1\t500\t.\t+\t.\tParent=m2\n"
        "chrmt\tX\texon\t1501\t2000\t.\t+\t.\tParent=m2\n"
    )
    pool = JunctionPool.from_annotation(p)
    assert {c for c, _s in pool._by_exon1} == {'chrI'}, "exon-derived fallback kept a mito intron"
