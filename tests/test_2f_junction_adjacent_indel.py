"""ISSUE-038 — no I/D touching the N in a 2F exon block (Kevin 2026-09-07: "Fix the leading indel blocks too").
The 74 that reached the BAM at aa75dfb were ALL minus-strand: the block abuts the intron at its FIRST op in
reference order there, and reading it as if it were plus-strand let `N 1D`, `N 3I`, `N 6I` through."""
import pytest
from rectify.core.splice.splice_aware_5prime import _junction_adjacent_indel_refusal as R, JUNCTION_INDEL_REFUSAL as TOK

# genome: exon side of a plus-strand donor at 20 ends in AAAA; intron follows
G = "C" * 12 + "TTAAAA" + "GT" + "A" * 20          # exon ...CCTTAAAA | GT..., donor at index 18
GM = "C" * 20 + "CT" + "TTTT" + "G" * 18            # minus: acceptor CT at 20, exon side starts at 22


def test_deletion_touching_the_n_is_always_refused():
    assert R([('M', 10), ('D', 1)], '+', G, 18, 40) == TOK
    assert R([('D', 3), ('M', 10)], '-', GM, 0, 22) == TOK


def test_insertion_is_refused_unless_a_run_explains_it():
    # exon side reads ...TTAAAA| : a 2-nt insertion continues the A-run (run 4) -> allowed
    assert R([('M', 10), ('I', 2)], '+', G, 18, 40) == ''
    # 6 inserted bases is more than the run -> refused
    assert R([('M', 10), ('I', 6)], '+', G, 18, 40) == TOK


def test_minus_strand_reads_the_block_from_its_first_op():
    # the same block on the minus strand: the I abuts the N at the FIRST op
    assert R([('I', 2), ('M', 10)], '-', GM, 0, 22) == ''      # T-run of 4 on the exon side
    assert R([('I', 9), ('M', 10)], '-', GM, 0, 22) == TOK
    # and a plus-strand reading of a minus-strand block must NOT be what decides it
    assert R([('M', 10), ('I', 9)], '-', GM, 0, 22) == ''      # last op is irrelevant on minus


def test_a_clean_block_passes():
    assert R([('M', 14)], '+', G, 18, 40) == ''
    assert R([('S', 2), ('M', 14)], '-', GM, 0, 22) == ''


def test_a_homopolymer_does_not_get_a_second_chance_as_a_dinucleotide():
    """The exon flank AAAA satisfies "two copies of the same 2-mer" trivially. Without requiring the
    2-mer's bases to DIFFER, a 6-nt insertion at a 4-nt A-run was admitted by the dinucleotide branch
    after the homopolymer branch had correctly refused it."""
    assert R([('M', 10), ('I', 6)], '+', G, 18, 40) == TOK
    assert R([('M', 10), ('I', 4)], '+', G, 18, 40) == ''


def test_dinucleotide_repeat_is_bounded_by_the_repeat_length():
    GD = "C" * 10 + "GAGAGA" + "GT" + "A" * 10        # exon flank reads AGAGAG...
    assert R([('M', 10), ('I', 6)], '+', GD, 16, 40) == ''
    assert R([('M', 10), ('I', 8)], '+', GD, 16, 40) == TOK
