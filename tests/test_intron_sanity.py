"""Consensus must not REPORT a physically impossible intron.

Real defect (Chanfreau planning/684c): uLTRA/deSALT emit alignments with N-ops
of hundreds of kb and the consensus scorer selects them — 268 of 400,001
primary reads (0.067 %) with an N-op > 10 kb, max 261,350 bp, 3 running off the
end of the contig, in a genome whose longest annotated intron is ~1 kb. The
minimap2 arm produced zero.

The two CIGARs exercised below are REAL records I pulled from
`684_p1cdna_1M/WT_BY4742_rep1/align/*.multialigned.bam` while verifying the
report — not invented shapes. Per this repo's CLAUDE.md ("vet alignments at the
INDIVIDUAL-READ level"), the fix is pinned against the actual reads it exists
to repair.
"""

import pytest

from rectify.core.consensus.intron_sanity import (
    DEFAULT_MAX_REPORTABLE_INTRON,
    max_reportable_intron_from_env,
    truncate_impossible_introns,
)

# op codes: 0=M 1=I 2=D 3=N 4=S
_QUERY_OPS = {0, 1, 4, 7, 8}


def _qlen(cig):
    return sum(l for op, l in cig if op in _QUERY_OPS)


def _reflen(cig):
    return sum(l for op, l in cig if op in {0, 2, 3, 7, 8})


# r030_7056, chrIV (len 1,531,933): 1S67M384N245M1I12M1I1M1I181M190957N2M1D53M57S
_R030 = [(4, 1), (0, 67), (3, 384), (0, 245), (1, 1), (0, 12), (1, 1), (0, 1),
         (1, 1), (0, 181), (3, 190957), (0, 2), (2, 1), (0, 53), (4, 57)]

# r036_7091, chrV (len 576,874): 239M525N803M301N68M150801N1D29M99483N1M196I28M102507N26M66S
_R036 = [(0, 239), (3, 525), (0, 803), (3, 301), (0, 68), (3, 150801), (2, 1),
         (0, 29), (3, 99483), (0, 1), (1, 196), (0, 28), (3, 102507), (0, 26), (4, 66)]


def test_clean_alignment_is_untouched():
    """A normal spliced yeast read must not be rewritten at all."""
    cig = [(4, 5), (0, 100), (3, 412), (0, 250), (4, 10)]
    new, bp = truncate_impossible_introns(cig, contig_len=1_000_000)
    assert new is None and bp == 0


def test_longest_real_yeast_intron_is_safe():
    """~1 kb is the longest annotated S. cerevisiae intron; 10 kb is 10x that."""
    cig = [(0, 100), (3, 1_002), (0, 250)]
    new, bp = truncate_impossible_introns(cig, contig_len=1_000_000)
    assert new is None and bp == 0


def test_r030_real_read_is_truncated_at_the_impossible_intron():
    new, bp = truncate_impossible_introns(_R030, contig_len=1_531_933)
    assert bp == 190_957
    # keep everything before the bad N; the rest of the QUERY becomes soft clip
    assert new == [(4, 1), (0, 67), (3, 384), (0, 245), (1, 1), (0, 12), (1, 1),
                   (0, 1), (1, 1), (0, 181), (4, 112)]
    assert _qlen(new) == _qlen(_R030), "query length must be preserved exactly"


def test_r036_real_read_truncates_at_the_FIRST_impossible_intron():
    """Three oversize N-ops; the cut must be at the first, not the largest."""
    new, bp = truncate_impossible_introns(_R036, contig_len=576_874)
    assert bp == 150_801
    assert new == [(0, 239), (3, 525), (0, 803), (3, 301), (0, 68), (4, 346)]
    assert _qlen(new) == _qlen(_R036)


def test_truncation_brings_the_read_back_inside_the_contig():
    """The point of the exercise: r030 ended 61.7 kb past chrIV."""
    start = 1_360_000
    contig = 1_531_933
    assert start + _reflen(_R030) > contig          # the bug
    new, _ = truncate_impossible_introns(_R030, reference_start=start,
                                         contig_len=contig)
    assert start + _reflen(new) <= contig           # the fix


def test_overrun_past_contig_end_is_caught_even_when_the_intron_is_small():
    """A CIGAR walking off the chromosome is malformed regardless of op size."""
    cig = [(0, 100), (3, 500), (0, 200)]
    new, bp = truncate_impossible_introns(cig, reference_start=999_000,
                                          contig_len=999_400)
    assert bp == 500 and new is not None
    assert 999_000 + _reflen(new) <= 999_400


def test_refuses_when_no_aligned_block_would_survive():
    """An all-soft-clip record is worse than leaving the read alone."""
    cig = [(4, 20), (3, 500_000), (0, 80)]
    new, bp = truncate_impossible_introns(cig, contig_len=1_000_000)
    assert new is None
    assert bp == 500_000, "still reported, so the counter sees it"


def test_trailing_soft_clip_is_merged_not_duplicated():
    cig = [(0, 50), (3, 900_000), (0, 10), (4, 5)]
    new, _ = truncate_impossible_introns(cig, contig_len=2_000_000)
    assert new == [(0, 50), (4, 15)]
    assert sum(1 for op, _ in new if op == 4) == 1


def test_trailing_reference_only_op_is_dropped():
    """A head ending in D/N would dangle against the new soft clip."""
    cig = [(0, 50), (2, 3), (3, 900_000), (0, 10)]
    new, _ = truncate_impossible_introns(cig, contig_len=2_000_000)
    assert new == [(0, 50), (4, 10)]


def test_disabled_by_zero_threshold():
    new, bp = truncate_impossible_introns(_R030, max_intron_bp=0,
                                          contig_len=1_531_933)
    assert new is None and bp == 0


def test_env_override(monkeypatch):
    monkeypatch.delenv('RECTIFY_MAX_REPORTABLE_INTRON', raising=False)
    assert max_reportable_intron_from_env() == DEFAULT_MAX_REPORTABLE_INTRON
    monkeypatch.setenv('RECTIFY_MAX_REPORTABLE_INTRON', '3000')
    assert max_reportable_intron_from_env() == 3000
    monkeypatch.setenv('RECTIFY_MAX_REPORTABLE_INTRON', 'nonsense')
    assert max_reportable_intron_from_env() == DEFAULT_MAX_REPORTABLE_INTRON


@pytest.mark.parametrize('cig', [_R030, _R036])
def test_query_length_is_never_changed(cig):
    """A changed query length desyncs SEQ/QUAL and crashes pysam on write."""
    new, _ = truncate_impossible_introns(cig, contig_len=2_000_000)
    assert _qlen(new) == _qlen(cig)


# --- the crash this defect caused, and the wiring that prevents it ----------

def test_ambiguity_window_survives_an_out_of_contig_junction():
    """`rectify pool-gate` died here with IndexError on a real 684 BAM.

    The unguarded `g[end - 1 - l_amb]` indexed past the reference string when a
    consensus-selected N-op ran off the contig. The guard must REPORT (0, 0)
    rather than raise — while never being the only fix, since alone it would
    hide the impossible alignment.
    """
    from rectify.core.splice.overhang_informativeness import ambiguity_window

    g = 'ACGT' * 250  # 1,000 bp contig
    # junction end far past the end of the contig — the real shape
    assert ambiguity_window(g, 900, 261_350) == (0, 0)
    assert ambiguity_window(g, 999_999, 1_000_050) == (0, 0)
    assert ambiguity_window(g, -5, 100) == (0, 0)
    # a normal in-bounds junction still computes a real window
    l_amb, r_amb = ambiguity_window(g, 100, 200)
    assert (l_amb, r_amb) != (0, 0)


def test_consensus_write_path_is_wired_to_the_guard():
    """Guard the call sites: the helper must run before BOTH out_bam.write()s.

    The unit tests above prove the CIGAR surgery; this proves it is actually
    reached. Without it the module could be perfect and never invoked — the
    same vacuity trap as tests/test_resolver_candidate_ceiling.py.
    """
    import inspect
    from rectify.core.consensus import consensus as cons

    src = inspect.getsource(cons)
    assert 'truncate_impossible_introns' in src
    assert src.count('_enforce_intron_sanity(') >= 3  # 1 def + 2 call sites
    assert "'impossible_intron_truncated'" in src or \
           '"impossible_intron_truncated"' in src

    body = inspect.getsource(cons._enforce_intron_sanity)
    assert "set_tag('Xn'" in body, "reads must carry the Xn flag for rbrowse"
