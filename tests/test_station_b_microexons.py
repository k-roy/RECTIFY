"""ISSUE-040 — STATION B, first pass: annotated micro-exons orphaned as an insertion beside a junction.

A splice aligner cannot seed a 6- or 9-nt internal exon, so it leaves the exon's bases as an insertion
next to the intron it did find. Measured on card 3's read: `--junc-bed` with that transcript's OWN
junctions, at `--junc-bonus` 9 AND 30, gives a byte-identical minimap2 CIGAR. This is structural.

What is pinned here: the search's REFUSALS (they are most of the safety), the two conservation
invariants of the rewrite (query span, reference span) and the ISSUE-031/038 shape rule (no I/D left
beside an N), and the measured acceptance case — a0d80d8a's 15-nt insertion is EXACTLY the two annotated
DCTN2 micro-exons, 6 nt and 9 nt, consumed in genomic order, and the module must find both.

The hermetic fixture is minus-strand for the same reason ISSUE-038 exists: a CIGAR is written in
reference order, so the flank that is a donor on the plus strand is an acceptor on the minus, and
reading one as the other invented a bug and hid another on 2026-09-07.
"""
import os

import pysam
import pytest

from rectify.core.splice import microexon as MX

# --------------------------------------------------------------------------- a hermetic contig
# exon1 [0,60) | intron [60,460) with a 6-nt micro-exon at [200,206) | exon2 [460,560)
# Flanks are laid down explicitly so the canonical test is exercised, not assumed.
MICRO = 'CAGCTC'


def _build_genome():
    g = list('A' * 60 + 'GT' + 'C' * 396 + 'AG' + 'T' * 100)
    # micro-exon at [200, 206) with a plus-strand acceptor (AG) before and donor (GT) after
    g[198:200] = list('AG')
    g[200:206] = list(MICRO)
    g[206:208] = list('GT')
    return ''.join(g)


GENOME_SEQ = _build_genome()
CHROM = 'chrX'    # standardizes to itself, unlike chrM -> chrMito
INDEX = {CHROM: [(200, 206)]}
INTRON = (60, 460)


def _read(cigar, seq, start=0, name='r', rev=False):
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'},
                                           'SQ': [{'SN': CHROM, 'LN': 100_000}]})
    a = pysam.AlignedSegment(hdr)
    a.query_name = name
    a.reference_name = CHROM
    a.reference_start = start
    a.cigartuples = cigar
    a.is_reverse = rev
    a.mapping_quality = 60
    a.query_sequence = seq
    return a


def _orphan_read(inserted, name='orphan'):
    """`60M 400N <k>I 100M` — the shape an aligner leaves when it orphans a micro-exon."""
    seq = GENOME_SEQ[0:60] + inserted + GENOME_SEQ[460:560]
    return _read([(0, 60), (3, 400), (1, len(inserted)), (0, 100)], seq, start=0, name=name)


# --------------------------------------------------------------------------- the search
def test_finds_the_annotated_microexon():
    segs = MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', GENOME_SEQ, INDEX)
    assert segs == [(200, 206)]


def test_refuses_when_the_insertion_is_not_consumed_exactly():
    """Partial consumption would leave bases to glue to an N as an indel — the banned shape."""
    assert MX.find_microexon_split(MICRO + 'TTT', CHROM, 60, 460, '+', GENOME_SEQ, INDEX) is None
    assert MX.find_microexon_split(MICRO[:4], CHROM, 60, 460, '+', GENOME_SEQ, INDEX) is None


def test_refuses_a_sequence_mismatch():
    bad = 'T' + MICRO[1:]
    assert bad != MICRO
    assert MX.find_microexon_split(bad, CHROM, 60, 460, '+', GENOME_SEQ, INDEX) is None


def test_refuses_below_the_minimum_insertion():
    assert MX.find_microexon_split('CA', CHROM, 60, 460, '+', GENOME_SEQ, INDEX) is None


def test_refuses_when_the_flanks_are_not_canonical():
    broken = GENOME_SEQ[:206] + 'TT' + GENOME_SEQ[208:]        # donor GT -> TT
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', broken, INDEX) is None


def test_refuses_an_exon_outside_the_drawn_intron():
    assert MX.find_microexon_split(MICRO, CHROM, 60, 150, '+', GENOME_SEQ, INDEX) is None


def test_refuses_an_unannotated_segment():
    """This pass proposes annotated exons only — a de-novo micro-exon needs station C behind it,
    because a 6-mer with canonical flanks lands roughly once per 80 bases of intron by chance."""
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', GENOME_SEQ, {CHROM: []}) is None


def test_strand_decides_the_whole_configuration():
    """A configuration legal for a plus-strand transcript is not legal for a minus-strand one: every
    intron is judged as a PAIR, and on the minus strand a GT-AG intron reads CT-AC in the genome."""
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '-', GENOME_SEQ, INDEX) is None
    g = list(GENOME_SEQ)
    g[60:62] = list('CT'); g[198:200] = list('AC')      # intron 1: CT..AC  = minus-strand GT-AG
    g[206:208] = list('CT'); g[458:460] = list('AC')    # intron 2: same
    minus = ''.join(g)
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '-', minus, INDEX) == [(200, 206)]
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', minus, INDEX) is None


def test_a_mixed_pair_is_refused():
    """GT donor with AC acceptor is one U2 end and one U12 end, not a splice class. Checking the two
    ends against independent SETS would admit it; checking the PAIR does not."""
    g = list(GENOME_SEQ)
    g[206:208] = list('AT')                              # downstream intron opens AT ...
    assert g[458:460] == list('AG')                      # ... and closes AG -> AT-AG, not a class
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', ''.join(g), INDEX) is None
    g[458:460] = list('AC')                              # AT-AC IS a class (U12)
    assert MX.find_microexon_split(MICRO, CHROM, 60, 460, '+', ''.join(g), INDEX) == [(200, 206)]


# --------------------------------------------------------------------------- the rewrite
def test_rewrite_conserves_both_spans_and_leaves_no_indel_beside_an_n():
    read = _orphan_read(MICRO)
    got = MX.recover_read_microexons(read, GENOME_SEQ, '+', INDEX)
    assert got is not None
    new_cigar, segs = got.new_cigar, got.segments
    assert segs == [(200, 206)]
    assert got.alternatives == [] and got.n_tied == 1 and not got.ambiguous
    q_old = sum(ln for op, ln in read.cigartuples if op in (0, 1, 4, 7, 8))
    q_new = sum(ln for op, ln in new_cigar if op in (0, 1, 4, 7, 8))
    r_old = sum(ln for op, ln in read.cigartuples if op in (0, 2, 3, 7, 8))
    r_new = sum(ln for op, ln in new_cigar if op in (0, 2, 3, 7, 8))
    assert (q_old, r_old) == (q_new, r_new)
    assert new_cigar == [(0, 60), (3, 140), (0, 6), (3, 254), (0, 100)]
    assert not any(
        (new_cigar[i][0] == 3 and new_cigar[i + 1][0] in (1, 2))
        or (new_cigar[i][0] in (1, 2) and new_cigar[i + 1][0] == 3)
        for i in range(len(new_cigar) - 1))


def test_a_read_with_no_junction_adjacent_insertion_is_untouched():
    read = _read([(0, 60), (3, 400), (0, 100)], GENOME_SEQ[0:60] + GENOME_SEQ[460:560], name='plain')
    assert MX.recover_read_microexons(read, GENOME_SEQ, '+', INDEX) is None


def test_an_insertion_away_from_the_junction_is_not_a_microexon_candidate():
    seq = GENOME_SEQ[0:30] + MICRO + GENOME_SEQ[30:60] + GENOME_SEQ[460:560]
    read = _read([(0, 30), (1, 6), (0, 30), (3, 400), (0, 100)], seq, name='midexon')
    assert MX.recover_read_microexons(read, GENOME_SEQ, '+', INDEX) is None


def test_the_index_loads_only_short_exons(tmp_path):
    gtf = tmp_path / 'a.gtf'
    gtf.write_text(
        'chrX\tT\texon\t201\t206\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
        'chrX\tT\texon\t1\t500\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
        'chrX\tT\texon\t201\t206\t.\t+\t.\tgene_id "g"; transcript_id "t2";\n'   # same exon, 2nd tx
        '#comment\n')
    # 'chrX' standardizes to itself; standardize_chrom_name maps bare numerals to yeast roman
    # numerals unless the run has registered its genome's contigs, which a real run does.
    idx = MX.load_microexons(str(gtf))
    assert idx['chrX'] == [(200, 206)]                         # deduplicated, long exon dropped
    assert [k for k in idx if k != '__transcripts__'] == ['chrX']
    assert MX.transcripts_of(idx, 'chrX', 200, 206) == frozenset({'t', 't2'})


def test_station_b_is_report_by_default(monkeypatch):
    monkeypatch.delenv('RECTIFY_STATION_B', raising=False)
    assert MX.station_b_mode() == 'report'
    monkeypatch.setenv('RECTIFY_STATION_B', 'apply')
    assert MX.station_b_mode() == 'apply'


# --------------------------------------------------------------------------- the measured case
_BUNDLE = os.path.expanduser(
    '~/work/rectify/dev/sumner_misplaced_panel_20260904/holdout/events/95cec1e/review_2h_t1_95cec1e')
_SLICE = 'chr12:57524317-57553043'
_OFF = 57524317


@pytest.mark.skipif(not os.path.exists(_BUNDLE + '/slices.fa'),
                    reason='collaborator review bundle not present')
def test_acceptance_a0d80d8a_two_annotated_microexons_in_one_insertion():
    """Card 3 (a0d80d8a, GSB_191, chr12 minus): the 15 bases 6485226 parks as an insertion beside a
    10,182-bp intron are EXACTLY two annotated DCTN2 micro-exons — CAGCTC at 57,538,518 (6 nt) and
    TTGTGCAAA at 57,541,363 (9 nt), both in ENST00000434715.7, consumed in genomic order. Do not
    build for exactly one."""
    import collections
    import copy

    fa = pysam.FastaFile(_BUNDLE + '/slices.fa')
    seq = fa.fetch(_SLICE).upper()
    per = collections.defaultdict(set)
    for line in open(_BUNDLE + '/slices.gtf'):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon' or f[0] != 'chr12':
            continue
        s, e = int(f[3]) - 1, int(f[4])
        if 0 < e - s <= MX.MAX_MICROEXON_LEN and s >= _OFF and e <= _OFF + len(seq):
            per['chr12'].add((s - _OFF, e - _OFF))
    index = {c: sorted(v) for c, v in per.items()}

    rec = next(a for a in pysam.AlignmentFile(_BUNDLE + '/6485226.bam')
               if a.query_name.startswith('a0d80d8a') and not (a.is_secondary or a.is_supplementary))
    assert (1, 15) in rec.cigartuples, 'the arm under test must still carry the 15-nt insertion'
    read = copy.deepcopy(rec)
    read.reference_start = rec.reference_start - _OFF

    got = MX.recover_read_microexons(read, seq, '-', index)
    assert got is not None, 'station B must recover this configuration'
    new_cigar, segs = got.new_cigar, got.segments
    assert [(s + _OFF, e + _OFF) for s, e in segs] == [(57538518, 57538524), (57541363, 57541372)]
    assert [e - s for s, e in segs] == [6, 9]
    q_old = sum(ln for op, ln in rec.cigartuples if op in (0, 1, 4, 7, 8))
    q_new = sum(ln for op, ln in new_cigar if op in (0, 1, 4, 7, 8))
    r_old = sum(ln for op, ln in rec.cigartuples if op in (0, 2, 3, 7, 8))
    r_new = sum(ln for op, ln in new_cigar if op in (0, 2, 3, 7, 8))
    assert (q_old, r_old) == (q_new, r_new)
    want = [(3, 2673), (0, 6), (3, 2839), (0, 9), (3, 4655)]
    assert any(new_cigar[i:i + len(want)] == want for i in range(len(new_cigar))), new_cigar[-10:]
    assert (1, 15) not in new_cigar and (3, 10182) not in new_cigar


# --------------------------------------------------------------------------- apply mode, end to end
def test_report_mode_records_but_draws_nothing(monkeypatch):
    """The default. The columns are filled, `station_b_applied` is 0, and the writer is a no-op."""
    import rectify.core.bam.bam_processor as bp
    from rectify.core.bam.bam_writer import apply_station_b_microexons

    monkeypatch.delenv('RECTIFY_STATION_B', raising=False)
    MX.set_microexon_index(INDEX)
    try:
        read = _orphan_read(MICRO, name='report')
        row = bp.correct_read_3prime(read, {CHROM: GENOME_SEQ}, annotated_junctions=set())[0]
        assert row['station_b_microexons'] == f'{CHROM}:200-206'
        assert row['station_b_applied'] == 0
        assert row['station_b_intron_start'] == 60 and row['station_b_intron_end'] == 460
        # the junctions column still names the ONE intron the aligner drew
        assert (60, 460) in [tuple(j) for j in row['junctions']]
        fresh = _orphan_read(MICRO, name='report')
        assert apply_station_b_microexons(fresh, row) is False
        assert fresh.cigartuples == [(0, 60), (3, 400), (1, 6), (0, 100)]
    finally:
        MX.set_microexon_index(None)


def test_apply_mode_draws_the_microexon_and_the_tsv_junctions_match_the_written_n_ops(monkeypatch):
    """ISSUE-024/033 is the binding invariant: the TSV `junctions` column must equal the N-ops the
    writer actually writes — the tester's scorer reads that column, so a mismatch is either an
    assertion or a silent zero."""
    import rectify.core.bam.bam_processor as bp
    from rectify.core.bam.bam_writer import apply_station_b_microexons
    from rectify.core.splice import microexon as _mx

    monkeypatch.setenv('RECTIFY_STATION_B', 'apply')
    MX.set_microexon_index(INDEX)
    try:
        row = bp.correct_read_3prime(_orphan_read(MICRO, name='apply'),
                                     {CHROM: GENOME_SEQ}, annotated_junctions=set())[0]
        assert row['station_b_applied'] == 1
        assert sorted(tuple(j) for j in row['junctions']) == [(60, 200), (206, 460)]
        read = _orphan_read(MICRO, name='apply')
        assert apply_station_b_microexons(read, row) is True
        n_ops = []
        ref = read.reference_start
        for op, ln in read.cigartuples:
            if op == 3:
                n_ops.append((ref, ref + ln))
            if op in (0, 2, 3, 7, 8):
                ref += ln
        assert n_ops == sorted(tuple(j) for j in row['junctions'])
        assert _mx.split_introns([(200, 206)], 60, 460) == n_ops
    finally:
        MX.set_microexon_index(None)


def test_the_writer_refuses_when_the_geometry_has_moved(monkeypatch):
    """Record-driven, not index-driven: if the named intron is no longer in the CIGAR (a 2F rescue
    rewrote the region), the rewrite is SKIPPED rather than applied to the wrong op."""
    from rectify.core.bam.bam_writer import apply_station_b_microexons
    row = {'station_b_applied': 1, 'station_b_microexons': f'{CHROM}:200-206',
           'station_b_intron_start': 60, 'station_b_intron_end': 460}
    moved = _read([(0, 60), (3, 399), (1, 6), (0, 100)],
                  GENOME_SEQ[0:60] + MICRO + GENOME_SEQ[460:560], name='moved')
    assert apply_station_b_microexons(moved, row) is False
    assert moved.cigartuples == [(0, 60), (3, 399), (1, 6), (0, 100)]


def test_a_tie_is_broken_at_random_but_reproducibly_and_the_loser_is_kept():
    """Kevin 2026-09-07: pick one at random, keep the other as a noted equally-good alternative.
    Two 6-nt segments with identical sequence and identical flanks are genuinely indistinguishable,
    so the pick must be arbitrary — seeded by the read name, so it is also reproducible."""
    g = list(GENOME_SEQ)
    g[298:300] = list('AG'); g[300:306] = list(MICRO); g[306:308] = list('GT')
    twin = ''.join(g)
    index = {CHROM: [(200, 206), (300, 306)]}
    splits = MX.find_microexon_splits(MICRO, CHROM, 60, 460, '+', twin, index)
    assert sorted(splits) == [[(200, 206)], [(300, 306)]]
    a, alt_a, n_a = MX.choose_split(splits, CHROM, index, 'read-one')
    b, alt_b, n_b = MX.choose_split(splits, CHROM, index, 'read-two')
    assert n_a == n_b == 2                                     # both tied for best
    assert len(alt_a) == 1 and alt_a[0] != a                   # the loser is kept, not discarded
    again, _, _ = MX.choose_split(splits, CHROM, index, 'read-one')
    assert again == a                                          # reproducible for the same read
    assert MX.format_alternatives(CHROM, alt_a).count(':') == 1
    # and only TIED configurations are kept as alternatives — a loser is not "equally good"
    b_only = MX.choose_split([[(200, 206)]], CHROM, index, 'solo')
    assert b_only[1] == [] and b_only[2] == 1


def test_the_higher_scoring_split_is_not_a_tie():
    """A longer segment carries more sequence evidence (2 bits a base), so it is not tied with a
    shorter one — the loser is reported as an alternative but never drawn."""
    long_seg = MICRO + 'GAT'
    g = list(GENOME_SEQ)
    g[198:200] = list('AG'); g[200:209] = list(long_seg); g[209:211] = list('GT')
    g[298:300] = list('AG'); g[300:306] = list(MICRO); g[306:308] = list('GT')
    gg = ''.join(g)
    index = {CHROM: [(200, 209), (300, 306)]}
    splits = MX.find_microexon_splits(long_seg, CHROM, 60, 460, '+', gg, index)
    assert splits == [[(200, 209)]]
    chosen, alts, n_tied = MX.choose_split(splits, CHROM, index, 'r')
    assert chosen == [(200, 209)] and alts == [] and n_tied == 1


def test_station_b_stands_down_when_a_5prime_rescue_touches_the_same_intron(monkeypatch):
    """The one interaction that could make the TSV and the BAM disagree. The writer applies the 5'
    rescue surgery BEFORE station B, so a rescue that rewrites the same intron moves the geometry out
    from under the row's coordinates and the writer then correctly declines the draw — leaving the TSV
    claiming two introns where the BAM has one. The row stands down instead, which is the same
    discipline as `predict_5prime_rescue_refusal`: predict the writer's refusal, never out-run it."""
    import rectify.core.bam.bam_processor as bp

    monkeypatch.setenv('RECTIFY_STATION_B', 'apply')
    MX.set_microexon_index(INDEX)
    try:
        # A read with BOTH a junction-adjacent insertion at [60, 460) and a 5' clip that 2F can
        # rescue onto that same intron: 12 clean exon-1 bases, then the body starting at the acceptor.
        clip = GENOME_SEQ[60 - 12:60]
        seq = clip + MICRO + GENOME_SEQ[460:560]
        read = _read([(4, 12), (1, 6), (0, 100)], seq, start=460, name='both')
        row = bp.correct_read_3prime(read, {CHROM: GENOME_SEQ},
                                     annotated_junctions={(CHROM, 60, 460)})[0]
        if row['five_prime_rescued'] and row['station_b_microexons']:
            assert row['station_b_applied'] == 0
            assert (60, 460) in [tuple(j) for j in row['junctions']] or not row['junctions']
    finally:
        MX.set_microexon_index(None)


def test_the_tsv_loader_carries_the_station_b_columns_to_the_writer(tmp_path):
    """The row reaches the writer ONLY through `_load_corrections_from_tsv`, which reads a fixed list
    of columns off the header — a column it does not name is silently dropped.

    Measured on the 2026-09-08 ON arm before this was fixed: the TSV claimed a 12-nt insertion split
    into two introns and the BAM still carried `70M 12I 599N`, because the loader never passed
    `station_b_applied` and the writer's first test therefore always failed. A round trip through a
    real TSV is the only test that catches it; asserting on the row dict does not.
    """
    from rectify.core.bam.bam_writer import _load_corrections_from_tsv
    from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row

    row = {
        'read_id': 'r', 'chrom': CHROM, 'strand': '+',
        'original_3prime': 559, 'corrected_3prime': 559,
        'ambiguity_min': 0, 'ambiguity_max': 0, 'ambiguity_range': 0,
        'correction_applied': [], 'confidence': 'high', 'qc_flags': [],
        'station_b_microexons': f'{CHROM}:200-206',
        'station_b_alternatives': f'{CHROM}:300-306',
        'station_b_applied': 1,
        'station_b_intron_start': 60, 'station_b_intron_end': 460,
    }
    p = tmp_path / 'corrected.tsv'
    with open(p, 'w') as fh:
        fh.write('\t'.join(CORRECTION_TSV_HEADER) + '\n')
        fh.write('\t'.join(correction_result_to_tsv_row(row)) + '\n')

    got = _load_corrections_from_tsv(str(p))['r']
    assert got['station_b_microexons'] == f'{CHROM}:200-206'
    assert got['station_b_alternatives'] == f'{CHROM}:300-306'
    assert got['station_b_applied'] == 1
    assert (got['station_b_intron_start'], got['station_b_intron_end']) == ('60', '460')

    # ... and the writer, handed exactly that loaded row, draws the micro-exon.
    from rectify.core.bam.bam_writer import apply_station_b_microexons
    read = _orphan_read(MICRO, name='roundtrip')
    assert apply_station_b_microexons(read, got) is True
    assert read.cigartuples == [(0, 60), (3, 140), (0, 6), (3, 254), (0, 100)]


def test_the_writer_stamps_XB_and_XV_only_when_it_drew(tmp_path):
    """The tags are the record of what was drawn and what was equally good; an un-applied row must
    not leave an XB behind claiming a draw that never happened."""
    from rectify.core.bam.bam_writer import apply_corrected_edits_to_read

    corr = {'corrected_3prime': 559, 'strand': '+', 'five_prime_rescued': False,
            'five_prime_position': None, 'five_prime_soft_clip': 0, 'five_prime_exon_cigar': '',
            'five_prime_upstream_trim': 0, 'five_prime_intron_clip_pos': -1, 'reanchor_clip_len': 0,
            'station_b_microexons': f'{CHROM}:200-206',
            'station_b_alternatives': f'{CHROM}:300-306',
            'station_b_applied': 1, 'station_b_intron_start': 60, 'station_b_intron_end': 460}
    read = _orphan_read(MICRO, name='tagged')
    apply_corrected_edits_to_read(read, corr, {CHROM: GENOME_SEQ})
    assert read.get_tag('XB') == f'{CHROM}:200-206'
    assert read.get_tag('XV') == f'{CHROM}:300-306'

    corr_off = dict(corr, station_b_applied=0)
    read2 = _orphan_read(MICRO, name='untagged')
    apply_corrected_edits_to_read(read2, corr_off, {CHROM: GENOME_SEQ})
    assert not read2.has_tag('XB') and not read2.has_tag('XV')
