"""Station C on HUMAN-shaped input: a GENCODE exon-GTF and a 2H-refined N-op.

Adopted verbatim in substance from the tester's hermetic regression
(``dev/sumner_misplaced_panel_20260904/PROPOSED_test_station_c_human_gtf.py``,
which failed 2/3 at base 4eefd1f). Both regressions were measured on real ONT
RNA004 DRS (Sumner GSB_2394 chr5, ``FINDINGS_agent4_stations.md``):

  A. ``load_annotated_canonical`` parsed GFF3 only (``Parent=`` + an ``intron``
     feature). GENCODE has neither, so ``n_annotated`` came back 0 on 222 real
     junctions and every annotated junction was reported as a discovery
     candidate.  ``derive_max_intron`` in the same repo already handles
     ``transcript_id "`` — this asserts Station C's copy does too.

  B. ``census_bam``'s ``min_anchor`` read the SINGLE adjacent CIGAR op and only
     when it was M/=/X.  Module 2H realizes a single-boundary junction move
     with a COMPENSATING I/D op placed right beside the N-op
     (``junction_refiner.py`` L1300-1309), which zeroed the anchor — so 87% of
     the junctions RECTIFY created (16/121 censused on the panel) were never
     censused at all: not admitted, not reviewed, not demoted. Station C's own
     q-scorer ``_side_features`` already walked I/D ops; only the census gate
     did not.

  C. On human, three of the four demotion terms have no track behind them.
     They must report ``track_unavailable``, not a clean 0/''.

The genome below is hermetic and 139 bp: 40 bp exon, a 59 bp GT..AG intron at
[40, 99) (0-based half-open), 40 bp exon.
"""
import pysam
import pytest

from rectify.core.consensus.station_c import (
    TRACK_UNAVAILABLE, PoolGateConfig, census_bam, load_annotated_canonical,
    pool_gate, write_pool_gate_outputs)

# ── a genome with one GT..AG intron at [40, 99) and 40 bp exons ─────────────
EX1 = 'ACGTTGCATTACGGCATTGCACGTTACGGTACCATGCATA'          # 40
INTRON = 'GT' + 'CATGACTGACTTGCACGATTGCAAGTACCTGATGCACGTTGCACGTACCTGACTG'[:56] + 'AG'  # 59
EX2 = 'TTGCACGTAACCGGTTACGATCGGATCCATTACGGCATGA'          # 40
CHRSEQ = EX1 + INTRON + EX2
GENOME = {'chr5': CHRSEQ}
J = (len(EX1), len(EX1) + len(INTRON))                     # (40, 99), 0-based half-open


def test_the_fixture_is_what_the_tests_assume():
    """Anti-vacuity: the planted intron really is a 59 bp GT..AG at [40, 99)."""
    assert len(CHRSEQ) == 139 and J == (40, 99)
    assert CHRSEQ[J[0]:J[0] + 2] == 'GT' and CHRSEQ[J[1] - 2:J[1]] == 'AG'


@pytest.fixture
def gencode_gtf(tmp_path):
    """A GENCODE-shaped exon GTF: transcript_id attributes, NO intron feature,
    NO Parent=. 1-based inclusive, exactly as GENCODE writes it."""
    p = tmp_path / 'mini.gencode.gtf'
    attrs = 'gene_id "ENSG00000000001.1"; transcript_id "ENST00000000001.1"; ' \
            'gene_type "protein_coding"; exon_number {}; tag "basic";'
    p.write_text(
        '##description: evidence-based annotation (test fixture)\n'
        f'chr5\tHAVANA\tgene\t1\t{len(CHRSEQ)}\t.\t+\t.\t{attrs.format(0)}\n'
        f'chr5\tHAVANA\ttranscript\t1\t{len(CHRSEQ)}\t.\t+\t.\t{attrs.format(0)}\n'
        f'chr5\tHAVANA\texon\t1\t{J[0]}\t.\t+\t.\t{attrs.format(1)}\n'
        f'chr5\tHAVANA\texon\t{J[1] + 1}\t{len(CHRSEQ)}\t.\t+\t.\t{attrs.format(2)}\n'
    )
    return p


def test_annotated_set_is_built_from_a_gencode_exon_gtf(gencode_gtf):
    """REGRESSION A — GENCODE GTF must yield the annotated intron, not {}."""
    ann = load_annotated_canonical(gencode_gtf, GENOME, PoolGateConfig())
    assert ('chr5', J[0], J[1]) in ann, (
        f"Station C read 0 annotated introns from a GENCODE exon-GTF "
        f"(got {sorted(ann)!r}); every annotated junction will be reported as "
        f"a discovery candidate. The parser needs the transcript_id \" branch "
        f"that multi_aligner.derive_max_intron already has.")


def test_gff3_parent_still_wins_over_transcript_id(tmp_path):
    """The GTF branch is a FALLBACK: a GFF3 line with Parent= must group on
    Parent=, even when a transcript_id attribute is also present (NCBI GFF3
    carries both). This is what keeps the yeast path byte-identical."""
    p = tmp_path / 'both.gff3'
    p.write_text(
        '##gff-version 3\n'
        f'chr5\ttest\texon\t1\t{J[0]}\t.\t+\t.\tParent=rna-A;transcript_id=shared\n'
        f'chr5\ttest\texon\t{J[1] + 1}\t{len(CHRSEQ)}\t.\t+\t.\t'
        f'Parent=rna-B;transcript_id=shared\n'
    )
    ann = load_annotated_canonical(p, GENOME, PoolGateConfig())
    # Two DIFFERENT Parent transcripts, one exon each -> no intron. If the
    # transcript_id fallback had run too, both would land in "shared" and a
    # phantom intron would appear.
    assert ann == set()


def _bam(tmp_path, name, cigar_reads):
    path = tmp_path / name
    header = {'HD': {'VN': '1.6'},
              'SQ': [{'SN': 'chr5', 'LN': len(CHRSEQ)}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for i, (start, cig, seq) in enumerate(cigar_reads):
            a = pysam.AlignedSegment()
            a.query_name = f'r{i}'
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigarstring = cig
            a.query_sequence = seq
            out.write(a)
    return str(path)


def test_census_sees_a_junction_whose_flank_carries_a_compensating_indel(tmp_path):
    """REGRESSION B — a 2H-refined N-op (compensating 1I beside it) must still
    be censused. Two reads so support clears min_support."""
    # 30M 1I 60N 30M : the exact shape Module 2H writes for a 1-nt donor slide.
    # 30 M + 1 I of query before the intron, 30 M after.
    reads = [(10, '30M1I60N30M', CHRSEQ[10:40] + 'A' + CHRSEQ[100:130])
             for _ in range(2)]
    bam = _bam(tmp_path, 'refined.bam', reads)
    J_census = census_bam(bam, GENOME, PoolGateConfig())
    assert J_census, (
        "census_bam enumerated NO junction: min_anchor reads only the single "
        "adjacent CIGAR op, so 2H's compensating indel (junction_refiner.py "
        "L1300-1309) makes every refined junction invisible to Station C — it "
        "is neither admitted, reviewed nor demoted. Sum the contiguous "
        "M/=/X run across intervening I/D ops, as _side_features already does.")


def test_census_still_rejects_a_genuinely_short_anchor(tmp_path):
    """Guard on the fix: a real 4 bp anchor must STILL be refused."""
    reads = [(36, '4M60N30M', CHRSEQ[36:40] + CHRSEQ[100:130]) for _ in range(2)]
    bam = _bam(tmp_path, 'shortanchor.bam', reads)
    assert not census_bam(bam, GENOME, PoolGateConfig())


def test_census_records_the_indel_it_stepped_over(tmp_path):
    """The compensating indel is the CIGAR signature of a MOVED boundary, so it
    must survive the walk into the census record — not be silently absorbed."""
    reads = [(10, '30M1I60N30M', CHRSEQ[10:40] + 'A' + CHRSEQ[100:130])
             for _ in range(2)]
    bam = _bam(tmp_path, 'refined.bam', reads)
    rec = next(iter(census_bam(bam, GENOME, PoolGateConfig()).values()))
    assert rec['support'] == 2
    assert rec['n_adj_indel'] == 2
    assert rec['adj_l'] == {'I1': 2}      # the 1I sits on the left of the N
    assert rec['adj_r'] == {}


def test_the_walk_is_bounded_in_ops_and_in_bases(tmp_path):
    """The walk is a coverage fix, not an unbounded anchor search: an indel run
    past either budget still terminates the anchor, so a rescued exon made of
    tiny match runs (2F's `4M4I3M4I1M` shape) stays refused."""
    cfg = PoolGateConfig()
    assert (cfg.max_adj_indel_ops, cfg.max_adj_indel_bp) == (2, 20)

    # 30 bp of intervening deletion > max_adj_indel_bp -> left anchor stops at 0
    seq = CHRSEQ[10:30] + CHRSEQ[100:130]
    over_bp = [(10, '20M30D40N30M', seq) for _ in range(2)]
    assert not census_bam(_bam(tmp_path, 'over_bp.bam', over_bp), GENOME, cfg)

    # three separate indel ops > max_adj_indel_ops, each small: the walk gets
    # 1M + 5M = 6 aligned bases before the third indel stops it -> < min_anchor
    q = (CHRSEQ[10:30] + 'AAA' + CHRSEQ[30:35] + 'CC' + CHRSEQ[36:37]
         + CHRSEQ[97:127])
    over_ops = [(10, '20M3I5M2I1D1M60N30M', q) for _ in range(2)]
    assert not census_bam(_bam(tmp_path, 'over_ops.bam', over_ops), GENOME, cfg)

    # ... and within budget it IS censused (anti-vacuity for the two above)
    ok = [(10, '20M3I10M60N30M', CHRSEQ[10:30] + 'AAA' + CHRSEQ[30:40]
           + CHRSEQ[100:130]) for _ in range(2)]
    assert census_bam(_bam(tmp_path, 'ok.bam', ok), GENOME, cfg)


def test_unavailable_flag_tracks_are_reported_not_silently_clean(tmp_path, gencode_gtf):
    """REGRESSION C — with no self-homology / background-SV BED and a GENCODE
    GTF (no REPEAT_FEATURE_TYPES features), all three flag columns must say
    'track_unavailable'. A 0/'' there reads as 'consulted and clean', which is
    how a human run looked identical to a fully-gated yeast one."""
    reads = [(10, '30M1I60N30M', CHRSEQ[10:40] + 'A' + CHRSEQ[100:130])
             for _ in range(2)]
    bam = _bam(tmp_path, 'refined.bam', reads)
    rows, summary = pool_gate(bam, GENOME, gencode_gtf)

    assert sorted(summary['tracks_unavailable']) == [
        'background_sv', 'repeat', 'selfhom']
    assert rows, "the fixture must produce a reportable junction"
    for r in rows:
        assert r['repeat_flag'] == TRACK_UNAVAILABLE
        assert r['selfhom_flag'] == TRACK_UNAVAILABLE
        assert r['background_sv_flag'] == TRACK_UNAVAILABLE

    # ... and the annotated shield is alive on this GTF, so the summary can be
    # trusted about what it did NOT report.
    assert summary['n_annotated_introns_parsed'] == 1

    # the new census columns reach the TSV (over_max_intron was computed and
    # silently dropped by the writer's fixed column list once already)
    tsv, _ = write_pool_gate_outputs(rows, summary, tmp_path / 'human')
    header = tsv.read_text().splitlines()[0].split('\t')
    for col in ('adj_indel_l', 'adj_indel_r', 'n_adj_indel'):
        assert col in header
