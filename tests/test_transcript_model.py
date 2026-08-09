#!/usr/bin/env python3
"""
Tests for the transcript-model-aware CPA classifier (planning/167 + 169).

Oracle coordinates (planning/170d) are verified against the bundled SGD GFF; the
model is 0-based half-open (loader-native) so oracle 1-based coords are converted
in-test. Categorical assertions (region_class, is_premature, intron_evidence,
antisense_context) are the load-bearing oracle; continuous coords are asserted by
sign + recomputation from the model's own coords (denominator conventions differ
by ±1 vs the oracle's hand arithmetic — see 170d conventions).

Author: Kevin R. Roy
"""

import os
import textwrap

import pandas as pd
import pytest

from rectify.core.analyze.loaders import load_annotation
from rectify.core.analyze.clustering import annotate_clusters_with_genes
from rectify.core.analyze.transcript_model import (
    TranscriptModel,
    NcrnaTrack,
    CLASSIFIER_COLUMNS,
    containment_attributions_from_clusters,
    annotate_clusters_with_transcript_model,
    collect_cryptic_junctions_from_tsvs,
    load_ncrna_atlas,
    build_ncrna_tracks_from_args,
    parse_ncrna_annotation_spec,
    _assert_build_compatible,
)
from rectify.core.commands.analyze_command import _make_cluster_lookup


_BUNDLED_GFF = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "rectify", "data", "genomes", "saccharomyces_cerevisiae",
    "saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz",
)


# =============================================================================
# Fixtures
# =============================================================================
@pytest.fixture(scope="module")
def sgd_annotation():
    if not os.path.exists(_BUNDLED_GFF):
        pytest.skip("bundled SGD GFF not available")
    return load_annotation(_BUNDLED_GFF, normalize_chroms=False)


@pytest.fixture(scope="module")
def model(sgd_annotation):
    return TranscriptModel(sgd_annotation)


def _cluster_row(cid, chrom, strand, pos, n_reads=50, radius=5):
    return {
        "cluster_id": cid, "chrom": chrom, "strand": strand,
        "start": pos - radius, "end": pos + radius,
        "modal_position": pos, "n_reads": n_reads,
    }


def _write_gff(tmp_path, lines):
    p = tmp_path / "mini.gff"
    p.write_text("##gff-version 3\n" + textwrap.dedent(lines).strip() + "\n")
    return str(p)


# =============================================================================
# (b) RPL7B / YPL198W — clean multi-intron gene
# =============================================================================
def test_multi_intron_annotated(model):
    g = model.genes["YPL198W"]
    assert g.strand == "+"
    # introns (0-based): [173162,173571) rank1, [173665,174072) rank2
    assert g.introns[0][2] == 1 and g.introns[1][2] == 2

    r1 = model.classify_site("chrXVI", "+", 173400)   # inside intron1
    assert r1["region_class"] == "intronic"
    assert r1["within_intron_id"] == "intron1"
    assert r1["intron_rank"] == 1
    assert r1["intron_evidence"] == "annotated"
    assert r1["is_premature"] is True           # upstream of the stop codon
    assert r1["attributed_gene"] == "YPL198W"   # clean systematic id, not a token

    r2 = model.classify_site("chrXVI", "+", 173900)   # inside intron2
    assert r2["region_class"] == "intronic"
    assert r2["intron_rank"] == 2

    # 3'UTR proximal (78 nt past the stop)
    r3 = model.classify_site("chrXVI", "+", g.cds_end - 1 + 78)
    assert r3["region_class"] == "3primeUTR_proximal"
    assert r3["is_premature"] is False
    assert r3["distance_to_stop_codon"] == 78


def test_utr3_proximal_distal_split(model, sgd_annotation):
    g = model.genes["YPL198W"]
    stop = g.cds_end - 1
    prox = model.classify_site("chrXVI", "+", stop + 50)
    dist = model.classify_site("chrXVI", "+", stop + 200)  # >150 split
    assert prox["region_class"] == "3primeUTR_proximal"
    assert dist["region_class"] == "3primeUTR_distal"
    # split is configurable: with a 300-nt split the same site is proximal
    m2 = TranscriptModel(sgd_annotation, utr3_proximal_distal_split=300)
    assert m2.classify_site("chrXVI", "+", stop + 200)["region_class"] == "3primeUTR_proximal"


# =============================================================================
# (c) RAD3 / YER171W — deep CDS-internal (MS2/NSD axis)
# =============================================================================
def test_cds_internal_deep(model):
    g = model.genes["YER171W"]
    assert g.strand == "+" and g.cds_start is not None
    r = model.classify_site("chrV", "+", 528000)
    assert r["region_class"] == "CDS_internal"
    assert r["is_premature"] is True
    assert r["distance_to_stop_codon"] < 0        # upstream of the stop
    assert 0.0 < r["orf_fraction"] < 1.0
    assert abs(r["orf_fraction"] - 0.393) < 0.01  # oracle 170d (c1)
    assert r["attributed_gene"] == "YER171W"


def test_containment_rescues_deep_cds(model):
    """Legacy 3'-window loses a deep-CDS cluster; containment rescues it."""
    clusters = pd.DataFrame([_cluster_row("c1", "chrV", "+", 528000)])
    legacy = annotate_clusters_with_genes(
        clusters.copy(), load_annotation(_BUNDLED_GFF, normalize_chroms=False)
    )
    assert legacy.loc[0, "gene_id"] in (None, "") or pd.isna(legacy.loc[0, "gene_id"])
    cont = containment_attributions_from_clusters(clusters, model)
    assert set(cont["gene_id"]) == {"YER171W"}


# =============================================================================
# (d) PGK1 / YCR012W — canonical 3'UTR; attribution MUST stay stable
# =============================================================================
def test_canonical_attribution_stable_both_modes(model):
    clusters = pd.DataFrame([_cluster_row("d1", "chrIII", "+", 139100)])
    ann = load_annotation(_BUNDLED_GFF, normalize_chroms=False)
    legacy = annotate_clusters_with_genes(clusters.copy(), ann)
    cont = containment_attributions_from_clusters(clusters, model)
    assert legacy.loc[0, "gene_id"] == "YCR012W"          # legacy window
    assert set(cont["gene_id"]) == {"YCR012W"}            # containment
    r = model.classify_site("chrIII", "+", 139100)
    assert r["region_class"] == "3primeUTR_proximal"
    assert r["is_premature"] is False


def test_readthrough_zone_keeps_gene_id(model):
    """A canonical cluster just past an annotated 3' end is NOT body-contained;
    the containment default must keep its gene_id (proximity fallback) and the
    classifier must call it downstream_readthrough — coherent, not intergenic
    (advisor items 1 & 2)."""
    g = model.genes["YCR012W"]
    pos = g.canonical_cpa + 40                              # +40 downstream, readthrough zone
    # not body-contained
    assert model.primary_sense_feature("chrIII", "+", pos) is None
    # classifier: coherent with a legacy fallback gene_id
    r = model.classify_site("chrIII", "+", pos, fallback_gene="YCR012W")
    assert r["region_class"] == "downstream_readthrough"
    assert r["attributed_gene"] == "YCR012W"
    # beyond the window -> intergenic
    r2 = model.classify_site("chrIII", "+", g.canonical_cpa + 500, fallback_gene="YCR012W")
    assert r2["region_class"] == "intergenic"


# =============================================================================
# (a) MSL5 / YLR116W — cryptic 3'UTR intron + convergent antisense CLF1
# =============================================================================
def test_msl5_antisense_axis(model):
    # antisense flips CLF1 3'UTR -> CDS at CLF1 CDS-start (382470 0-based)
    r_utr = model.classify_site("chrXII", "+", 382350)   # < CLF1 cds-start
    r_cds = model.classify_site("chrXII", "+", 382720)   # >= CLF1 cds-start (canonical CPA)
    assert r_utr["antisense_context"] == "3primeUTR:YLR117C"
    assert r_cds["antisense_context"] == "CDS:YLR117C"
    # antisense is NEVER a region_class (169 D4)
    assert not str(r_utr["region_class"]).startswith("antisense")
    assert r_cds["distance_to_canonical_CPA"] == 0        # canonical CPA


def test_msl5_cryptic_intron_registration(model):
    """Without the cryptic layer a1 is 3'UTR; after registering the population
    cryptic junction it is cryptic_intronic (167 §B.2 — the reference-only gap)."""
    # fresh model so registration does not leak across tests
    m = TranscriptModel(load_annotation(_BUNDLED_GFF, normalize_chroms=False))
    before = m.classify_site("chrXII", "+", 382350)
    assert before["region_class"].startswith("3primeUTR")   # reference-only miss
    # cryptic intron 0-based [382300, 382412)
    assert m.is_annotated_junction("chrXII", "+", 382300, 382412) is False
    n = m.register_cryptic_introns([("chrXII", "+", 382300, 382412)])
    assert n == 1
    after = m.classify_site("chrXII", "+", 382350)
    assert after["region_class"] == "cryptic_intronic"
    assert after["intron_evidence"] == "cryptic-read-junction"
    assert after["within_intron_id"] == "cryptic:chrXII:382300-382412"
    assert after["is_premature"] is False                   # downstream of the stop


def test_cryptic_noise_guards(tmp_path):
    """DRS junction-noise guards: min_support (drop 1-read junctions) + require a
    junction to lie within a transcribed body (drop genome-wide intergenic noise)."""
    gff = _write_gff(tmp_path, """
        chr1\ttest\tgene\t1000\t3000\t.\t+\t.\tID=G1;gene_id=G1;gene_name=GENE1
        chr1\ttest\tmRNA\t1000\t3000\t.\t+\t.\tID=G1_id1;Parent=G1;Name=G1_id1
        chr1\ttest\tCDS\t1100\t2900\t.\t+\t.\tID=c1;Parent=G1_id1;Name=G1_CDS
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    # within-body, support 2 -> registered; within-body support 1 -> dropped;
    # out-of-body support 2 -> dropped
    n = m.register_cryptic_introns([
        ("chr1", "+", 1500, 1700), ("chr1", "+", 1500, 1700),   # support 2, in body
        ("chr1", "+", 1800, 1900),                               # support 1, in body
        ("chr1", "+", 5000, 5200), ("chr1", "+", 5000, 5200),   # support 2, NO body
    ], min_support=2)
    assert n == 1
    assert m.classify_site("chr1", "+", 1600)["region_class"] == "cryptic_intronic"
    assert m.classify_site("chr1", "+", 1850)["region_class"] != "cryptic_intronic"


def test_annotated_junction_fuzzy_absorbs_ont_imprecision(model):
    """F2 (SUPERSEDES 170c §4's exact-match intent): ONT/DRS splice-site imprecision
    on a REAL annotated intron must NOT register a false cryptic intron. A junction
    within ±tol on BOTH donor and acceptor is treated as the annotated intron; a
    genuinely-novel junction (beyond tol) is still cryptic. ``tol=0`` forces exact.

    (170c §4 deliberately made a 1-bp shift read as cryptic; 172c F2 reverses that as
    a bug for the RECTIFY target data type — this test encodes the corrected policy.)"""
    m = TranscriptModel(load_annotation(_BUNDLED_GFF, normalize_chroms=False))
    # RPL7B intron1 = 0-based [173162, 173571); default junction_match_tol = 8
    assert m.is_annotated_junction("chrXVI", "+", 173162, 173571) is True          # exact
    # a 3-bp acceptor shift is WITHIN tol -> treated as the annotated intron (fuzzy)
    assert m.is_annotated_junction("chrXVI", "+", 173162, 173574) is True
    # ... but with exact matching (tol=0) the shift is NOT annotated
    assert m.is_annotated_junction("chrXVI", "+", 173162, 173574, tol=0) is False
    # a 50-bp shift is a genuinely-distinct junction -> NOT annotated (still cryptic)
    assert m.is_annotated_junction("chrXVI", "+", 173162, 173621) is False
    # registration: exact + the 3-bp-imprecise reads collapse onto the annotated
    # intron (0 cryptic there); only the 50-bp-distinct junction registers as cryptic.
    n = m.register_cryptic_introns([
        ("chrXVI", "+", 173162, 173571),   # annotated exact -> ignored
        ("chrXVI", "+", 173162, 173574),   # 3bp imprecise   -> absorbed (not cryptic)
        ("chrXVI", "+", 173162, 173574),   #                    (2nd supporting read)
        ("chrXVI", "+", 173162, 173621),   # 50bp distinct   -> cryptic
        ("chrXVI", "+", 173162, 173621),   #                    (2nd supporting read)
    ], min_support=2)
    assert n == 1


# =============================================================================
# Cryptic-intron read pass end-to-end (corrected TSV → cluster classification)
# =============================================================================
def test_cryptic_read_pass_end_to_end(model, tmp_path):
    m = TranscriptModel(load_annotation(_BUNDLED_GFF, normalize_chroms=False))
    clusters = pd.DataFrame([
        _cluster_row("a1", "chrXII", "+", 382350),   # inside cryptic span (no own junction)
        _cluster_row("a2", "chrXII", "+", 382451),   # intron-proximal pA (carries junction)
        _cluster_row("a3", "chrXII", "+", 382720),   # distal pA (carries junction)
    ])
    # synthetic corrected TSV: the downstream pAs carry the cryptic junction
    tsv = tmp_path / "corrected.tsv"
    rows = [
        "chrom\tstrand\tcorrected_3prime\tjunctions\tn_junctions",
        "chrXII\t+\t382451\t382300-382412\t1",
        "chrXII\t+\t382451\t382300-382412\t1",
        "chrXII\t+\t382720\t382300-382412\t1",
        "chrXII\t+\t382720\t\t0",                     # a read w/o junction
    ]
    tsv.write_text("\n".join(rows) + "\n")
    lookup = _make_cluster_lookup(clusters)
    samples = [{"sample_id": "s", "path": str(tsv)}]

    cryptic, spliced = collect_cryptic_junctions_from_tsvs(samples, m, lookup)
    assert ("chrXII", "+", 382300, 382412) in cryptic
    assert "a2" in spliced and "a3" in spliced and "a1" not in spliced

    out = annotate_clusters_with_transcript_model(
        clusters, m, samples=samples, cluster_lookup=lookup,
    )
    by = out.set_index("cluster_id")
    assert by.loc["a1", "region_class"] == "cryptic_intronic"    # rescued via population
    assert bool(by.loc["a2", "spliced_support"]) is True
    assert by.loc["a1", "region_class"] == "cryptic_intronic"


def test_read_pass_skips_tsv_without_junctions(model, tmp_path):
    """Older corrected TSVs lack junction columns → skip gracefully, no crash."""
    tsv = tmp_path / "old.tsv"
    tsv.write_text("chrom\tstrand\tcorrected_3prime\nchrXII\t+\t382720\n")
    clusters = pd.DataFrame([_cluster_row("a3", "chrXII", "+", 382720)])
    lookup = _make_cluster_lookup(clusters)
    cryptic, spliced = collect_cryptic_junctions_from_tsvs(
        [{"sample_id": "s", "path": str(tsv)}], model, lookup
    )
    assert cryptic == [] and spliced == set()


# =============================================================================
# (e) ncRNA atlas — CUT antisense to a CDS
# =============================================================================
def test_ncrna_cut_antisense_to_cds(model, sgd_annotation):
    cut_df = pd.DataFrame([{
        "chrom": "chrXV", "start": 624206, "end": 624338, "strand": "-",
        "gene_id": "CUT310", "gene_name": "CUT310", "feature_type": "CUT",
    }])
    m = TranscriptModel(
        sgd_annotation,
        ncrna_tracks=[NcrnaTrack(df=cut_df, source="VeraDowell2016", ncrna_class="CUT")],
    )
    r = m.classify_site("chrXV", "-", 624206)         # CUT 3' end (- strand, low coord)
    assert r["region_class"] == "ncRNA:CUT"
    assert r["attributed_feature_id"] == "CUT310"
    assert r["attribution_source"] == "VeraDowell2016"
    assert r["antisense_context"] == "CDS:YOR153W"    # PDR5 opposite strand
    assert r["distance_to_ncrna_3prime"] == 0
    assert abs(r["pos_within_ncrna_frac"] - 1.0) < 1e-9

    # SGD-core-only fallback: intergenic on the - strand, antisense STILL populated
    r0 = model.classify_site("chrXV", "-", 624206)
    assert r0["region_class"] == "intergenic"
    assert r0["antisense_context"] == "CDS:YOR153W"


def test_ncrna_atlas_registry_default():
    tracks = load_ncrna_atlas(genome_build="R64-5-1")   # default = demo_synthetic
    assert tracks and tracks[0].ncrna_class == "CUT"


def test_ncrna_atlas_ad_hoc_spec(tmp_path):
    assert parse_ncrna_annotation_spec("/a/b.gff:Xu2009:CUT") == ("/a/b.gff", "Xu2009", "CUT")
    with pytest.raises(ValueError):
        parse_ncrna_annotation_spec("bad_spec")
    gff = _write_gff(tmp_path, """
        chrXV\tsyn\tSUT\t624207\t624338\t.\t-\t.\tID=s1;Name=s1
    """)
    tracks = build_ncrna_tracks_from_args(None, [f"{gff}:MyStudy:SUT"], genome_build="R64-5-1")
    assert tracks[0].source == "MyStudy" and tracks[0].ncrna_class == "SUT"


def test_build_assertion_tolerates_label_lag():
    _assert_build_compatible("R64-4-1", "R64-5-1")   # must not raise
    _assert_build_compatible("sacCer3", "R64-1-1")
    with pytest.raises(ValueError):
        _assert_build_compatible("hg38", "R64-5-1")


# =============================================================================
# Edge cases: minus strand, nested sense genes, single-exon, degraded input
# =============================================================================
def test_minus_strand_intron(model):
    """ACT1 / YFL039C is - strand with one intron [54377,54686)."""
    g = model.genes["YFL039C"]
    assert g.strand == "-"
    r = model.classify_site("chrVI", "-", 54500)      # inside the intron
    assert r["region_class"] == "intronic"
    assert r["intron_evidence"] == "annotated"
    assert r["is_premature"] is True
    assert r["attributed_gene"] == "YFL039C"


def test_nested_sense_genes_upstream_tiebreak(tmp_path):
    """Two overlapping sense-strand coding genes: containment picks the
    transcription-upstream one (advisor item 6)."""
    # mRNA rows carry Name=<transcript-token> exactly as SGD does (170a §2);
    # that is the token the CDS Parent join resolves through.
    gff = _write_gff(tmp_path, """
        chr1\ttest\tgene\t1000\t3000\t.\t+\t.\tID=GBIG
        chr1\ttest\tmRNA\t1000\t3000\t.\t+\t.\tID=GBIG_id1;Parent=GBIG;Name=GBIG_id1
        chr1\ttest\tCDS\t1100\t2900\t.\t+\t.\tID=c1;Parent=GBIG_id1;Name=GBIG_CDS
        chr1\ttest\tgene\t2000\t2600\t.\t+\t.\tID=GSMALL
        chr1\ttest\tmRNA\t2000\t2600\t.\t+\t.\tID=GSMALL_id1;Parent=GSMALL;Name=GSMALL_id1
        chr1\ttest\tCDS\t2100\t2500\t.\t+\t.\tID=c2;Parent=GSMALL_id1;Name=GSMALL_CDS
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    pos = 2300   # inside both bodies (0-based)
    primary = m.primary_sense_feature("chr1", "+", pos)
    assert primary["gene_id"] == "GBIG"           # upstream start wins
    r = m.classify_site("chr1", "+", pos)
    assert r["attributed_gene"] == "GBIG"
    assert "GSMALL" in (r["also_within"] or "")   # loser recorded, nothing dropped


def test_single_exon_gene_utr_derivation(tmp_path):
    gff = _write_gff(tmp_path, """
        chr1\ttest\tgene\t1000\t2200\t.\t+\t.\tID=G1;gene_id=G1;gene_name=GENE1
        chr1\ttest\tmRNA\t1000\t2200\t.\t+\t.\tID=G1_id1;Parent=G1;Name=G1_id1
        chr1\ttest\tCDS\t1100\t2000\t.\t+\t.\tID=cds1;Parent=G1_id1;Name=G1_CDS
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    g = m.genes["G1"]
    assert g.is_coding and g.cds_start == 1099 and g.cds_end == 2000
    # 5'UTR
    assert m.classify_site("chr1", "+", 1050)["region_class"] == "5primeUTR"
    # CDS
    assert m.classify_site("chr1", "+", 1500)["region_class"] == "CDS_internal"
    # 3'UTR (2000..2200)
    assert m.classify_site("chr1", "+", 2100)["region_class"].startswith("3primeUTR")


def test_cds_gap_intron_fallback(tmp_path):
    """Annotation with NO explicit intron rows -> derive introns from CDS gaps."""
    gff = _write_gff(tmp_path, """
        chr1\ttest\tgene\t1000\t3000\t.\t+\t.\tID=G1;gene_id=G1;gene_name=GENE1
        chr1\ttest\tmRNA\t1000\t3000\t.\t+\t.\tID=G1_id1;Parent=G1;Name=G1_id1
        chr1\ttest\tCDS\t1100\t1500\t.\t+\t.\tID=c1;Parent=G1_id1;Name=G1_CDS
        chr1\ttest\tCDS\t2000\t2800\t.\t+\t.\tID=c2;Parent=G1_id1;Name=G1_CDS
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    g = m.genes["G1"]
    assert len(g.introns) == 1                      # gap [1500, 1999) derived
    r = m.classify_site("chr1", "+", 1700)          # inside the derived intron
    assert r["region_class"] == "intronic"


def test_degraded_annotation_warns_not_crash(tmp_path):
    """Exon-only GTF (no CDS/mRNA) -> warns, no coding classes, no crash."""
    gff = _write_gff(tmp_path, """
        chr1\ttest\texon\t1000\t2000\t.\t+\t.\tID=e1;gene_id=G1
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    assert any("coding" in w for w in m.warnings)
    assert m.classify_site("chr1", "+", 1500)["region_class"] == "intergenic"


def test_empty_annotation():
    m = TranscriptModel(pd.DataFrame(
        columns=["chrom", "start", "end", "strand", "feature_type"]
    ))
    assert m.classify_site("chr1", "+", 100)["region_class"] == "intergenic"


def test_classifier_columns_present(model):
    clusters = pd.DataFrame([
        _cluster_row("x1", "chrXVI", "+", 173400),
        _cluster_row("x2", "chrIII", "+", 139100),
    ])
    out = annotate_clusters_with_transcript_model(clusters, model)
    for col in CLASSIFIER_COLUMNS:
        assert col in out.columns
    # the classifier never sets gene_id (that is the attribution dispatch's job)
    assert "gene_id" not in out.columns
    assert out.set_index("cluster_id").loc["x1", "region_class"] == "intronic"


# =============================================================================
# 172c adversarial-review fixes (planning/177): F1, F3, F7, F8
# =============================================================================
def test_f1_spliced_orf_coords_not_genomic(tmp_path):
    """F1: orf_fraction / distance_to_stop_codon are SPLICED, not genomic. Two-exon
    gene, 1799-bp intron: a pA in exon2 (intron 5' of it) and a pA in exon1 (intron
    between it and the stop) must both get spliced transcript coordinates."""
    # CDS exon1 [1099,1200) 101nt, intron [1200,2999) 1799nt, exon2 [2999,3100) 101nt
    # spliced ORF = 202 nt (orf_fraction denominator = 201)
    gff = _write_gff(tmp_path, """
        chr1\ttest\tgene\t1100\t3100\t.\t+\t.\tID=G1;gene_id=G1;gene_name=GENE1
        chr1\ttest\tmRNA\t1100\t3100\t.\t+\t.\tID=G1_id1;Parent=G1;Name=G1_id1
        chr1\ttest\tCDS\t1100\t1200\t.\t+\t.\tID=c1;Parent=G1_id1;Name=G1_id1
        chr1\ttest\tCDS\t3000\t3100\t.\t+\t.\tID=c2;Parent=G1_id1;Name=G1_id1
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    g = m.genes["G1"]
    assert len(g.introns) == 1 and g.introns[0][:2] == (1200, 2999)
    assert g.exon_blocks == ((1099, 1200), (2999, 3100))
    # pA in exon2 at 3050: spliced offset-from-start = 101 + 51 = 152 -> 152/201
    r2 = m.classify_site("chr1", "+", 3050)
    assert r2["region_class"] == "CDS_internal"
    assert abs(r2["orf_fraction"] - (152 / 201)) < 1e-9
    assert r2["orf_fraction"] < 0.80            # genomic would be 1951/1999 = 0.976
    assert r2["distance_to_stop_codon"] == -49  # no intron between 3050 and stop 3099
    # pA in exon1 at 1150: the 1799-bp intron lies BETWEEN it and the stop
    r1 = m.classify_site("chr1", "+", 1150)
    assert r1["region_class"] == "CDS_internal"
    assert r1["distance_to_stop_codon"] == -150     # spliced; genomic would be -1949
    assert abs(r1["orf_fraction"] - (51 / 201)) < 1e-9


def test_f1_spliced_coords_minus_strand_symmetric(tmp_path):
    """F1: the spliced projection is symmetric on the minus strand (the reviewer
    confirmed minus was already off-by-one-free; keep it so)."""
    gff = _write_gff(tmp_path, """
        chr1\ttest\tgene\t1100\t3100\t.\t-\t.\tID=GM;gene_id=GM;gene_name=GM
        chr1\ttest\tmRNA\t1100\t3100\t.\t-\t.\tID=GM_id1;Parent=GM;Name=GM_id1
        chr1\ttest\tCDS\t1100\t1200\t.\t-\t.\tID=c1;Parent=GM_id1;Name=GM_id1
        chr1\ttest\tCDS\t3000\t3100\t.\t-\t.\tID=c2;Parent=GM_id1;Name=GM_id1
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    # start codon at high coord (3099); stop at low coord (1099). spliced ORF = 202.
    r_near_start = m.classify_site("chr1", "-", 3050)    # near start (49 nt in)
    r_near_stop = m.classify_site("chr1", "-", 1150)     # near stop
    assert abs(r_near_start["orf_fraction"] - (49 / 201)) < 1e-9
    assert r_near_start["distance_to_stop_codon"] == -152
    assert abs(r_near_stop["orf_fraction"] - (150 / 201)) < 1e-9
    assert r_near_stop["distance_to_stop_codon"] == -51


def test_f3_nested_gene_intron_attributed_to_owner(tmp_path):
    """F3: a pA in an INNER sense gene's intron, while the primary (widest) gene is
    exonic there, is classified intronic and attributed to the intron's OWNER — the
    within_intron_id is no longer lost."""
    gff = _write_gff(tmp_path, """
        chr1\ttest\tgene\t1000\t5000\t.\t+\t.\tID=GBIG;gene_id=GBIG;gene_name=BIG
        chr1\ttest\tmRNA\t1000\t5000\t.\t+\t.\tID=GBIG_id1;Parent=GBIG;Name=GBIG_id1
        chr1\ttest\tCDS\t1100\t4900\t.\t+\t.\tID=cb;Parent=GBIG_id1;Name=GBIG_id1
        chr1\ttest\tgene\t2000\t3000\t.\t+\t.\tID=GSMALL;gene_id=GSMALL;gene_name=SMALL
        chr1\ttest\tmRNA\t2000\t3000\t.\t+\t.\tID=GSMALL_id1;Parent=GSMALL;Name=GSMALL_id1
        chr1\ttest\tCDS\t2100\t2300\t.\t+\t.\tID=cs1;Parent=GSMALL_id1;Name=GSMALL_id1
        chr1\ttest\tCDS\t2700\t2900\t.\t+\t.\tID=cs2;Parent=GSMALL_id1;Name=GSMALL_id1
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    assert m.genes["GSMALL"].introns[0][:2] == (2300, 2699)      # GSMALL's derived intron
    assert m.genes["GBIG"].introns == ()                        # GBIG spans it exonically
    assert m.primary_sense_feature("chr1", "+", 2500)["gene_id"] == "GBIG"   # widest wins
    r = m.classify_site("chr1", "+", 2500)                      # inside GSMALL's intron
    assert r["region_class"] == "intronic"      # was CDS_internal (bug)
    assert r["within_intron_id"] == "intron1"   # was None (bug)
    assert r["intron_rank"] == 1
    assert r["intron_evidence"] == "annotated"
    assert r["attributed_gene"] == "GSMALL"     # intron owner, not the primary GBIG
    assert "GBIG" in (r["also_within"] or "")   # primary still recorded, nothing dropped


def test_f7_cds_gap_fallback_is_per_gene(tmp_path):
    """F7: one gene's explicit intron row must NOT disable CDS-gap intron derivation
    for a different gene that has none of its own (was a silent global flag)."""
    gff = _write_gff(tmp_path, """
        chr1\ttest\tgene\t1000\t3000\t.\t+\t.\tID=GA;gene_id=GA;gene_name=GA
        chr1\ttest\tmRNA\t1000\t3000\t.\t+\t.\tID=GA_id1;Parent=GA;Name=GA_id1
        chr1\ttest\tCDS\t1100\t1500\t.\t+\t.\tID=ca1;Parent=GA_id1;Name=GA_id1
        chr1\ttest\tintron\t1500\t2000\t.\t+\t.\tID=iA;Parent=GA_id1;Name=GA_id1
        chr1\ttest\tCDS\t2000\t2900\t.\t+\t.\tID=ca2;Parent=GA_id1;Name=GA_id1
        chr2\ttest\tgene\t1000\t3000\t.\t+\t.\tID=GB;gene_id=GB;gene_name=GB
        chr2\ttest\tmRNA\t1000\t3000\t.\t+\t.\tID=GB_id1;Parent=GB;Name=GB_id1
        chr2\ttest\tCDS\t1100\t1500\t.\t+\t.\tID=cb1;Parent=GB_id1;Name=GB_id1
        chr2\ttest\tCDS\t2000\t2900\t.\t+\t.\tID=cb2;Parent=GB_id1;Name=GB_id1
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    assert m.genes["GA"].introns                             # explicit intron row honored
    assert m.genes["GB"].introns[0][:2] == (1500, 1999)      # gap-derived despite GA's row
    assert m.classify_site("chr2", "+", 1700)["region_class"] == "intronic"


def test_f7_short_gap_is_frameshift_not_intron(tmp_path):
    """F7 guard: a sub-``min_gap_intron_len`` CDS gap (a +1 frameshift site — 47 SGD
    Ty ORFs have 1-bp gaps) must NOT be derived as an intron."""
    gff = _write_gff(tmp_path, """
        chr1\ttest\tgene\t1000\t3000\t.\t+\t.\tID=GF;gene_id=GF;gene_name=GF
        chr1\ttest\tmRNA\t1000\t3000\t.\t+\t.\tID=GF_id1;Parent=GF;Name=GF_id1
        chr1\ttest\tCDS\t1100\t2000\t.\t+\t.\tID=cf1;Parent=GF_id1;Name=GF_id1
        chr1\ttest\tCDS\t2002\t2900\t.\t+\t.\tID=cf2;Parent=GF_id1;Name=GF_id1
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    assert m.genes["GF"].introns == ()                          # 1-bp gap [2000,2001) skipped
    assert m.classify_site("chr1", "+", 2000)["region_class"] == "CDS_internal"


def test_f8_multi_gene_overlap_emits_all_weights(tmp_path):
    """F8: a cluster contained in multiple sense genes emits ALL of them with weights
    (primary STRICTLY highest), not just the primary; a single-gene cluster is
    unchanged (one row at weight 1.0)."""
    gff = _write_gff(tmp_path, """
        chr1\ttest\tgene\t1000\t5000\t.\t+\t.\tID=GBIG;gene_id=GBIG;gene_name=BIG
        chr1\ttest\tmRNA\t1000\t5000\t.\t+\t.\tID=GBIG_id1;Parent=GBIG;Name=GBIG_id1
        chr1\ttest\tCDS\t1100\t4900\t.\t+\t.\tID=cb;Parent=GBIG_id1;Name=GBIG_id1
        chr1\ttest\tgene\t2000\t3000\t.\t+\t.\tID=GSMALL;gene_id=GSMALL;gene_name=SMALL
        chr1\ttest\tmRNA\t2000\t3000\t.\t+\t.\tID=GSMALL_id1;Parent=GSMALL;Name=GSMALL_id1
        chr1\ttest\tCDS\t2100\t2900\t.\t+\t.\tID=cs;Parent=GSMALL_id1;Name=GSMALL_id1
    """)
    m = TranscriptModel(load_annotation(gff, normalize_chroms=False))
    clusters = pd.DataFrame([_cluster_row("c1", "chr1", "+", 2500, n_reads=50)])
    attr = containment_attributions_from_clusters(clusters, m)
    by_gene = attr.set_index("gene_id")["attribution_weight"]
    assert set(attr["gene_id"]) == {"GBIG", "GSMALL"}    # both emitted, GSMALL not dropped
    assert by_gene["GBIG"] > by_gene["GSMALL"]           # primary strictly highest
    assert abs(by_gene.sum() - 1.0) < 1e-9               # normalized within cluster
    # single-gene cluster (the common case) still emits exactly one row at weight 1.0
    solo = pd.DataFrame([_cluster_row("c2", "chr1", "+", 4000, n_reads=30)])  # only GBIG
    solo_attr = containment_attributions_from_clusters(solo, m)
    assert list(solo_attr["gene_id"]) == ["GBIG"]
    assert abs(float(solo_attr["attribution_weight"].iloc[0]) - 1.0) < 1e-9
