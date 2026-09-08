"""One pipeline, two species — the yeast side of everything built for the human SMA panel.

Kevin, 2026-09-08: *"Let's make sure to include yeast regression tests on all these human dev
implementations too. I want our single pipeline to work well on all species, with the species knobs
set appropriately. In our case, both human and yeast can utilize GT-AG/GC-AG/AT-AC, but yeast also
have the non-canonical 3' SS arms BG and AT, etc."*

He is describing a hierarchy the tree already carries, in
``junction_scoring._3ss_tier_from_rna_trinucleotide``, derived from yeast splicing observations:

===== ======= ==============================================
tier  motif
===== ======= ==============================================
0     YAG     (C/T)AG — most common, highest efficiency
1     RAG     (A/G)AG
2     NBG     B = C/G/T — "non-canonical but observed in yeast"   <- his "BG"
3     NAT     very rare non-canonical                             <- his "AT"
4     other
===== ======= ==============================================

ISSUE-040 shipped with a hard-coded ``{GT-AG, GC-AG, AT-AC}`` pair set, which silently refuses every
yeast NBG and NAT acceptor. That is the defect these tests pin: the ceiling is a per-organism knob
read from ``--organism``, the model is the shared one, and there is not a second copy to drift.

What is covered here, for yeast, of the work built on human data 2026-09-07/08:
  ISSUE-039  station C site support        ISSUE-042  the conditioned gap-cap relaxation
  ISSUE-040  micro-exon recovery           ISSUE-043  every insertion, not only the first
  ISSUE-044  junction mismatch enrichment
"""
import collections

import pysam
import pytest

from rectify.core.splice import microexon as MX
from rectify.core.splice import junction_scoring as J
from rectify.core.splice import splice_aware_5prime as S

CHROM = 'chrIV'          # a real S. cerevisiae contig name; standardize_chrom_name leaves it alone


@pytest.fixture(autouse=True)
def _clean(monkeypatch):
    monkeypatch.delenv('RECTIFY_MICROEXON_MAX_3SS_TIER', raising=False)
    monkeypatch.delenv('RECTIFY_2F_GAP_RELAX', raising=False)
    MX.set_species(None)
    MX.set_microexon_index(None)
    S.set_site_support(None)
    S.set_junction_mismatch(None)
    yield
    MX.set_species(None)
    MX.set_microexon_index(None)
    S.set_site_support(None)
    S.set_junction_mismatch(None)


# --------------------------------------------------------------------------- the species knob
def test_the_3ss_ceiling_comes_from_the_organism():
    MX.set_species(None)
    assert MX.max_3ss_tier() == MX.MICROEXON_MAX_3SS_TIER_DEFAULT == 1     # YAG / RAG
    MX.set_species('homo_sapiens')
    assert MX.max_3ss_tier() == 1
    MX.set_species('Saccharomyces cerevisiae')                             # spaces and case tolerated
    assert MX.max_3ss_tier() == 3                                          # NBG and NAT reachable
    MX.set_species('saccharomyces_cerevisiae')
    assert MX.max_3ss_tier() == 3


def test_the_tier_model_is_the_shared_one_not_a_second_copy():
    """If these ever disagree, two motif models have drifted apart — which is the thing the
    ISSUE-040 hard-coded pair set was doing before this."""
    assert J._3ss_tier_from_rna_trinucleotide('CAG') == 0        # YAG
    assert J._3ss_tier_from_rna_trinucleotide('AAG') == 1        # RAG
    assert J._3ss_tier_from_rna_trinucleotide('ACG') == 2        # NBG  <- Kevin's "BG"
    assert J._3ss_tier_from_rna_trinucleotide('CAT') == 3        # NAT  <- Kevin's "AT"
    assert J._3ss_tier_from_rna_trinucleotide('CCC') == 4


# --------------------------------------------------------------------------- ISSUE-040 on yeast
MICRO = 'CAGCTC'


def _locus(acceptor3_upstream='CAG', acceptor3_downstream='CAG', donor_after='GT'):
    """exon1 [0,60) | intron [60,460) with a 6-nt micro-exon at [200,206) | exon2 [460,560).

    The two acceptors are settable so a yeast NBG/NAT arm can be planted where a human run would
    only ever see YAG.
    """
    g = list('A' * 60 + 'GT' + 'C' * 396 + 'CAG'[-2:] + 'T' * 100)
    g = list('A' * 60 + 'GT' + 'C' * 395 + 'CAG' + 'T' * 100)   # intron ends ...CAG at 457:460
    g[197:200] = list(acceptor3_upstream)      # the 3 bases closing intron 1, RNA orientation
    g[200:206] = list(MICRO)
    g[206:208] = list(donor_after)             # the 2 bases opening intron 2
    g[457:460] = list(acceptor3_downstream)    # the 3 bases closing intron 2
    return ''.join(g)


INDEX = {CHROM: [(200, 206)]}


def test_a_yeast_NBG_acceptor_microexon_is_refused_on_the_default_ceiling_and_kept_on_yeast():
    """The defect, pinned. `ACG` is tier 2 (NBG): real in yeast, refused everywhere else."""
    g = _locus(acceptor3_upstream='ACG')
    MX.set_species('homo_sapiens')
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', g, INDEX) is None
    MX.set_species('saccharomyces_cerevisiae')
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', g, INDEX) == [(200, 206)]


def test_a_yeast_NAT_acceptor_microexon_behaves_the_same_way():
    """`CAT` is tier 3 (NAT) — "very rare non-canonical", and still yeast."""
    g = _locus(acceptor3_upstream='CAT')
    MX.set_species('homo_sapiens')
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', g, INDEX) is None
    MX.set_species('saccharomyces_cerevisiae')
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', g, INDEX) == [(200, 206)]


def test_tier_4_is_refused_in_yeast_too():
    """The ceiling is a ceiling, not an amnesty: `CCC` closes nothing in any organism."""
    MX.set_species('saccharomyces_cerevisiae')
    g = _locus(acceptor3_upstream='CCC')
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', g, INDEX) is None


def test_the_canonical_classes_work_in_BOTH_species():
    """GT-AG, GC-AG and AT-AC are Kevin's shared set — neither organism may lose them."""
    for organism in ('homo_sapiens', 'saccharomyces_cerevisiae'):
        MX.set_species(organism)
        assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+',
                                       _locus(donor_after='GT'), INDEX) == [(200, 206)], organism
        assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+',
                                       _locus(donor_after='GC'), INDEX) == [(200, 206)], organism
        # AT-AC: the U12 pair, admitted in every organism rectify runs on (Talkish 2019)
        g = _locus(donor_after='AT', acceptor3_downstream='CAC')
        assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', g, INDEX) == [(200, 206)], organism


def test_a_mixed_pair_is_still_refused_in_yeast():
    """An AT donor with an AG acceptor is one U12 end and one U2 end. A wider 3'SS ceiling widens
    the U2 hierarchy; it does not dissolve the pairing."""
    MX.set_species('saccharomyces_cerevisiae')
    g = _locus(donor_after='AT', acceptor3_downstream='CAG')
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', g, INDEX) is None


def test_minus_strand_yeast_locus_reads_its_ends_in_transcript_orientation():
    """The ISSUE-038 lesson, on the species path: resolve orientation BEFORE judging ends."""
    MX.set_species('saccharomyces_cerevisiae')
    g = list(_locus())
    # Each intron is judged on its OWN ends in transcript orientation, so for a minus-strand read
    # the donor is the reverse complement of its RIGHT end and the acceptor of its LEFT three bases.
    # intron (60, 200):  donor = rc(g[198:200]) = GT ; acceptor = rc(g[60:63]) = CAG
    g[198:200] = list('AC'); g[60:63] = list('CTG')
    # intron (206, 460): donor = rc(g[458:460]) = GT ; acceptor = rc(g[206:209]) = CAG
    g[458:460] = list('AC'); g[206:209] = list('CTG')
    gg = ''.join(g)
    assert gg[200:206] == MICRO                       # the micro-exon itself is untouched
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '-', gg, INDEX) == [(200, 206)]
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', gg, INDEX) is None


# --------------------------------------------------------------------------- ISSUE-039 on yeast
def _spanning_bam(tmp_path, name, n_reads, anchor, junction, aperiodic):
    p = tmp_path / f'{name}.minimap2.bam'
    hdr = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': CHROM, 'LN': 1_531_933}]}   # real chrIV length
    _c, i_s, i_e = junction
    with pysam.AlignmentFile(str(p), 'wb', header=hdr) as out:
        for k in range(n_reads):
            a = pysam.AlignedSegment(out.header)
            a.query_name = f'{name}_{k}'
            a.reference_name = CHROM
            a.reference_start = i_s - anchor
            a.cigartuples = [(0, anchor), (3, i_e - i_s), (0, anchor)]
            a.mapping_quality = 60
            a.query_sequence = aperiodic[:anchor] + aperiodic[-anchor:]
            out.write(a)
    pysam.index(str(p))
    return str(p)


APERIODIC = 'ACGTTGCATGCAGTCCATGACGGATCTAGCATCG'
YEAST_J = (CHROM, 100_000, 100_400)      # a plausible yeast intron: short


def test_station_C_site_support_on_a_yeast_locus(tmp_path):
    """Yeast introns are short and yeast loci are dense; the support rule must behave identically —
    it is about anchors and arms, not about intron length."""
    arms = [_spanning_bam(tmp_path, f'arm{i}', 4, J.SITE_SUPPORT_ANCHOR, YEAST_J, APERIODIC)
            for i in range(3)]
    _all, _annot, sig = J.build_junction_pool(arms, set(), return_signal=True)
    assert sig['site_support'][YEAST_J] == 4          # max over arms, never 12
    S.set_site_support(sig['site_support'])
    assert S.site_established(YEAST_J)
    assert S.site_support_n(YEAST_J) == 4


def test_station_C_short_anchors_do_not_establish_a_yeast_site(tmp_path):
    """A yeast exon can be shorter than the 20-base strict anchor; such a junction is NOT
    established, and that is the correct conservative answer rather than a species exception."""
    arm = _spanning_bam(tmp_path, 'short', 6, 12, YEAST_J, APERIODIC)
    _all, _annot, sig = J.build_junction_pool([arm], set(), return_signal=True)
    assert sig['site_support'].get(YEAST_J, 0) == 0


# --------------------------------------------------------------------------- ISSUE-042 on yeast
def test_the_conditioned_gap_relaxation_is_species_neutral():
    """The two conditions are about the block's shape and the population, not about motifs, so they
    must read identically on a yeast junction."""
    from rectify.core.align.local_aligner import cigar_str_to_ops
    ops = cigar_str_to_ops('6M1I9M6D3M1D3M1I3M')
    matched = sum(n for o, n in ops if o in (0, 7, 8))
    S.set_site_support(None)
    assert S._gap_refusal(ops, junction=YEAST_J) == S.EXON_GAP_REFUSAL
    S.set_site_support({YEAST_J: 9})
    assert S._gap_relax_allowed(ops, matched, YEAST_J)
    assert S._gap_refusal(ops, junction=YEAST_J) == ''
    heavy = cigar_str_to_ops('7M1I4M1I3M1D12M12I1M2I9M8I6M13I10M')
    assert not S._gap_relax_allowed(heavy, sum(n for o, n in heavy if o in (0, 7, 8)), YEAST_J)


# --------------------------------------------------------------------------- ISSUE-044 on yeast
def test_junction_mismatch_enrichment_reads_the_same_on_a_yeast_junction():
    stats = {YEAST_J: [8, 100, 10, 1000, 12]}        # 8 % near vs 1 % far
    assert J.junction_mismatch_enrichment(stats, YEAST_J) == pytest.approx(8.0)
    S.set_junction_mismatch(stats)
    assert S.junction_mismatch_enrichment(YEAST_J) == pytest.approx(8.0)
    # too little evidence either side -> no number rather than a wrong one
    assert J.junction_mismatch_enrichment({YEAST_J: [1, 10, 5, 500, 2]}, YEAST_J) is None
    assert J.junction_mismatch_enrichment({YEAST_J: [5, 100, 0, 1000, 9]}, YEAST_J) is None
    assert J.junction_mismatch_enrichment({}, YEAST_J) is None


def test_mismatch_bins_split_a_yeast_read_at_its_own_junction(tmp_path):
    """The metric on a real record: a read whose mismatches sit AT the junction must score near-rate
    above far-rate, and the same read with the mismatches moved into its body must not."""
    genome = tmp_path / 'g.fa'
    body = (APERIODIC * 12)[:400]
    ref = 'A' * 40 + body + 'GT' + 'C' * 96 + 'AG' + body + 'T' * 40
    genome.write_text(f'>{CHROM}\n' + ref + '\n')
    pysam.faidx(str(genome))
    fa = pysam.FastaFile(str(genome))

    i_s, i_e = 440, 540                       # the intron in `ref`
    left, right = ref[340:440], ref[540:640]

    def read_with(muts):
        seq = list(left + right)
        for i in muts:
            seq[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[seq[i]]
        hdr = pysam.AlignmentHeader.from_dict(
            {'HD': {'VN': '1.6'}, 'SQ': [{'SN': CHROM, 'LN': len(ref)}]})
        a = pysam.AlignedSegment(hdr)
        a.query_name = 'r'
        a.reference_name = CHROM
        a.reference_start = 340
        a.cigartuples = [(0, 100), (3, i_e - i_s), (0, 100)]
        a.mapping_quality = 60
        a.query_sequence = ''.join(seq)
        return a

    edges = [i_s, i_e]
    at_junction = read_with([95, 96, 97, 98, 102, 103, 104, 105])       # within 12 nt of an edge
    n_mm, n_bp, f_mm, f_bp, prox = J._mismatch_bins(
        at_junction, fa, edges, J.JUNCTION_MM_WINDOW, J.JUNCTION_MM_VARIANT_REACH)
    assert n_mm == 8 and f_mm == 0
    assert (n_mm / n_bp) > (f_mm / f_bp if f_bp else 0)

    in_body = read_with([5, 6, 7, 8, 190, 191, 192, 193])               # far from both edges
    n_mm2, n_bp2, f_mm2, f_bp2, _ = J._mismatch_bins(
        in_body, fa, edges, J.JUNCTION_MM_WINDOW, J.JUNCTION_MM_VARIANT_REACH)
    assert n_mm2 == 0 and f_mm2 == 8
    # and the proximal tally only reaches as far as it claims to
    assert all(abs(p - i_s) <= J.JUNCTION_MM_VARIANT_REACH or
               abs(p - i_e) <= J.JUNCTION_MM_VARIANT_REACH for p, _ in prox)


def test_the_pool_carries_the_enrichment_blocks_for_yeast(tmp_path):
    """End to end on the pool builder: with a genome in hand the signal dict gains both ISSUE-044
    blocks, and without one it does not — an older cache stays valid."""
    genome = tmp_path / 'g.fa'
    ref = 'A' * 200 + 'GT' + 'C' * 396 + 'AG' + 'T' * 200
    genome.write_text(f'>{CHROM}\n' + ref + '\n')
    pysam.faidx(str(genome))
    junction = (CHROM, 200, 600)
    arm = _spanning_bam(tmp_path, 'arm', 5, 30, junction, APERIODIC)

    _a, _b, sig = J.build_junction_pool([arm], set(), return_signal=True, fasta_path=str(genome))
    assert junction in sig['junction_mismatch']
    near_mm, near_bp, far_mm, far_bp, n_reads = sig['junction_mismatch'][junction]
    assert n_reads == 5 and near_bp > 0
    assert 'mismatch_positions' in sig

    _a2, _b2, sig2 = J.build_junction_pool([arm], set(), return_signal=True)
    assert sig2['junction_mismatch'] == {} and sig2['mismatch_positions'] == {}


# --------------------------------------------------------------------------- ISSUE-043 on yeast
def test_two_yeast_microexons_on_one_read_are_both_recovered():
    """Yeast introns are short, so two micro-exons on one read sit close together — the iteration
    must still resolve both, and the second must not be silently dropped."""
    MX.set_species('saccharomyces_cerevisiae')
    unit = _locus()
    genome = unit + unit                        # two identical units, the second offset by len(unit)
    off = len(unit)
    index = {CHROM: [(200, 206), (off + 200, off + 206)]}
    seq = (genome[0:60] + MICRO + genome[460:560]
           + genome[off:off + 60] + MICRO + genome[off + 460:off + 560])
    hdr = pysam.AlignmentHeader.from_dict(
        {'HD': {'VN': '1.6'}, 'SQ': [{'SN': CHROM, 'LN': len(genome)}]})
    a = pysam.AlignedSegment(hdr)
    a.query_name = 'two_yeast'
    a.reference_name = CHROM
    a.reference_start = 0
    a.cigartuples = [(0, 60), (3, 400), (1, 6), (0, 100),
                     (0, 60), (3, 400), (1, 6), (0, 100)]
    a.mapping_quality = 60
    a.query_sequence = seq
    calls = MX.recover_all_microexons(a, genome, '+', index)
    assert [c.segments for c in calls] == [[(200, 206)], [(off + 200, off + 206)]]
    assert MX.recover_read_microexons(a, genome, '+', index).intron == (60, 460)
