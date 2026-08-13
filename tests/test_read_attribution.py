"""Tests for the per-read gene-attribution sidecar (`702` spec, `703` design).

The behaviours pinned here are the ones that were measured on real data and would
silently regress: the rescue-proof 5' source, the termination boundary that must
not invent readthrough in the wild type, and the per-read key that survives
re-clustering.
"""

import csv
import textwrap

import pytest

from rectify.core.analyze import read_attribution as ra


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

# Two tandem plus-strand genes with a gap, mirroring the PLB3/RCL1 geometry the
# spec was written from (PLB3 chrXV:304,592-307,689; RCL1 307,882-309,265).
GFF = """\
##gff-version 3
chrXV\tSGD\tgene\t304593\t307689\t.\t+\t.\tID=GENEA;Name=GENEA
chrXV\tSGD\tCDS\t304593\t307000\t.\t+\t0\tParent=GENEA_id001;Name=GENEA_CDS
chrXV\tSGD\tgene\t307883\t309265\t.\t+\t.\tID=GENEB;Name=GENEB
chrXV\tSGD\tCDS\t307938\t309041\t.\t+\t0\tParent=GENEB_id001;Name=GENEB_CDS
chrXV\tSGD\tgene\t320001\t322000\t.\t-\t.\tID=GENEC;Name=GENEC
chrXV\tSGD\tCDS\t320501\t322000\t.\t-\t0\tParent=GENEC_id001;Name=GENEC_CDS
"""

TSV_HEADER = [
    "read_id", "chrom", "strand", "original_3prime", "corrected_3prime",
    "five_prime_position", "five_prime_rescued",
    "alignment_start", "alignment_end", "strand_evidence",
]


@pytest.fixture
def gff(tmp_path):
    p = tmp_path / "test.gff"
    p.write_text(GFF)
    return str(p)


@pytest.fixture
def genes(gff):
    return ra.GeneIndex.from_gff(gff)


def write_tsv(tmp_path, rows, name="corrected_reads.tsv"):
    p = tmp_path / name
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(TSV_HEADER)
        w.writerows(rows)
    return str(p)


def read(rid, start, end, three_p, strand="+", rescued="0", five_p=None,
         chrom="chrXV", strand_ev=""):
    return [rid, chrom, strand, str(three_p), str(three_p),
            str(five_p if five_p is not None else start), rescued,
            str(start), str(end), strand_ev]


# --------------------------------------------------------------------------
# Annotation
# --------------------------------------------------------------------------

def test_gene_features_are_utr_inclusive_and_cds_is_separate(genes):
    """The GFF `gene` feature spans the transcript; CDS marks the stop codon.

    This is the fact that invalidated the original "gene intervals end at the
    CDS" premise, so it is pinned explicitly.
    """
    chrom, s, e, strand = genes.spans["GENEB"]
    assert (s, e) == (307882, 309265)          # transcript extent, UTR-inclusive
    assert genes.cds_end["GENEB"] == 309040    # stop codon, 224 bp inside it
    assert genes.annotated_end_offset_base("GENEB") == 309264


def test_minus_strand_cds_end_is_the_leftmost_base(genes):
    assert genes.cds_end["GENEC"] == 320500


def test_overlapping_orders_by_overlap_size(genes):
    # spans all of GENEA (3096 bp) and the first 100 bp of GENEB
    hits = genes.overlapping("chrXV", "+", 304592, 307982)
    assert hits == ["GENEA", "GENEB"]


# --------------------------------------------------------------------------
# The 5' coordinate trap
# --------------------------------------------------------------------------

def test_observed_5prime_uses_raw_alignment_not_the_rescued_column(tmp_path, genes):
    """A rescued read's 5' position must NOT be what origin5 is derived from.

    `extend_read_5prime_for_junction_rescue` moves `five_prime_position`
    upstream; deriving "began upstream of X" from it would manufacture the
    evidence. The raw alignment_start is used instead.
    """
    # five_prime_position claims 300000 (far upstream); alignment says 305000
    path = write_tsv(tmp_path, [
        read("r1", 305000, 307000, 306999, rescued="1", five_p=300000),
    ])
    rec = next(ra._iter_corrected(path))
    assert rec["five_p"] == 305000
    assert rec["rescued"] is True


def test_observed_5prime_minus_strand():
    assert ra.observed_5prime("-", 100, 200) == 199
    assert ra.observed_5prime("+", 100, 200) == 100


def test_signed_offset_is_downstream_positive_on_both_strands():
    assert ra.signed_offset_past("+", 150, 100) == 50
    assert ra.signed_offset_past("-", 50, 100) == 50


# --------------------------------------------------------------------------
# Termination boundary
# --------------------------------------------------------------------------

def test_boundary_falls_back_to_annotated_end_plus_measured_tail(genes):
    """A gene with no CPA profile still gets the measured genome-wide tail."""
    assert ra.termination_boundary("GENEA", genes, {}) == ra.GENOME_TAIL_BP


def test_boundary_extends_to_a_downstream_observed_cpa(genes):
    """A gene whose real CPA sits far past its annotation keeps its own CPA.

    This is the `YGL007C-A` case (n=2,136, 46% of reads +3,249 bp past the
    annotated end) that refutes any fixed window.
    """
    ref = {"GENEA": {"n": 2136, "mode": 3249, "q50": 3249, "q95": 3253}}
    assert ra.termination_boundary("GENEA", genes, ref) == 3249 + ra.GENOME_TAIL_BP


def test_boundary_never_pulled_upstream_by_an_internal_cpa(genes):
    """A gene whose modal CPA is INSIDE its body must not shrink the boundary.

    71.6% of genes have a negative modal offset; clamping at 0 stops those from
    reclassifying their own 3'UTR reads as escaping.
    """
    ref = {"GENEA": {"n": 500, "mode": -2321, "q50": -100, "q95": 7}}
    assert ra.termination_boundary("GENEA", genes, ref) == ra.GENOME_TAIL_BP


# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------

def _classify(tmp_path, genes, row, ref=None):
    path = write_tsv(tmp_path, [row])
    rec = next(ra._iter_corrected(path))
    return ra.classify_read(rec, genes, ref or {})


def test_read_confined_to_its_gene_does_not_escape(tmp_path, genes):
    out = _classify(tmp_path, genes, read("r", 305000, 307600, 307599))
    assert out["attr_primary"] == "GENEA"
    assert out["origin5"] == "internal"
    assert out["escapes_gene"] == ""
    assert out["escapes_gene_cpa"] == ""


def test_read_ending_just_past_the_annotation_does_not_escape(tmp_path, genes):
    """The wild-type population that containment-only turns into false readthrough.

    Measured: containment-only takes PLB3's WT escape from 0/845 to ~204/845.
    A 3' end in the intergenic gap is the gene's own 3'UTR, not an escape.
    """
    out = _classify(tmp_path, genes, read("r", 305000, 307800, 307799))
    assert out["escapes_gene"] == ""
    assert out["escapes_gene_cpa"] == ""


def test_readthrough_into_the_next_gene_escapes(tmp_path, genes):
    """The AT1 case: initiates in GENEA, terminates deep inside GENEB."""
    out = _classify(tmp_path, genes, read("r", 305000, 309276, 309275))
    assert out["escapes_gene"] == "GENEA"
    assert out["escapes_gene_cpa"] == "GENEA"
    assert out["attr_genes"] == "GENEA|GENEB"
    assert out["origin5"] == "internal"      # 5' end is inside GENEA


def test_escape_respects_a_long_observed_utr(tmp_path, genes):
    """With a distant real CPA, the same read is the gene's own 3'UTR, not escape."""
    ref = {"GENEA": {"n": 2136, "mode": 3249, "q50": 3249, "q95": 3253}}
    out = _classify(tmp_path, genes, read("r", 305000, 309276, 309275), ref=ref)
    assert out["escapes_gene_cpa"] == ""      # per-gene rule keeps it
    assert out["escapes_gene"] == "GENEA"     # fixed-window contrast column does not


def test_upstream_origin_is_recorded_with_its_evidence(tmp_path, genes):
    """5' end upstream of the gene it overlaps — provable, not inferred."""
    out = _classify(tmp_path, genes, read("r", 300000, 307000, 306999))
    assert out["origin5"] == "upstream_origin"
    assert out["origin5_evidence"] == "drs_5p_upstream"


def test_rescued_read_is_labelled_not_laundered(tmp_path, genes):
    """A moved 5' coordinate must be visible to the consumer."""
    out = _classify(tmp_path, genes, read("r", 300000, 307000, 306999, rescued="1"))
    assert out["origin5"] == "upstream_origin"
    assert out["origin5_evidence"] == "five_prime_rescued"


def test_initiating_is_unreachable_without_cdna_evidence(tmp_path, genes):
    """AT2: no DRS read may ever be called `initiating`."""
    for start in (304592, 305000, 307688):
        out = _classify(tmp_path, genes, read("r", start, 307689, 307688))
        assert out["origin5"] != "initiating"


def test_intergenic_read_is_unknown_not_forced(tmp_path, genes):
    out = _classify(tmp_path, genes, read("r", 312000, 313000, 312999))
    assert out["attr_primary"] == ""
    assert out["attr_rule"] == "none"
    assert out["origin5"] == "unknown"


def test_dist_to_stop_is_measured_from_the_cds_not_the_transcript_end(tmp_path, genes):
    out = _classify(tmp_path, genes, read("r", 308000, 309100, 309099))
    assert out["attr_primary"] == "GENEB"
    assert int(out["dist_to_stop"]) == 309099 - 309040


def test_region_class_is_empty_in_v1(tmp_path, genes):
    """Deliberate: the vocabulary is cluster-keyed and the spec forbids that join."""
    out = _classify(tmp_path, genes, read("r", 305000, 307000, 306999))
    assert out["region_class"] == ""


# --------------------------------------------------------------------------
# Schema / join
# --------------------------------------------------------------------------

def test_missing_alignment_columns_is_a_hard_error(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("read_id\tchrom\tstrand\tcorrected_3prime\nr\tchrXV\t+\t100\n")
    with pytest.raises(ValueError, match="alignment_start"):
        list(ra._iter_corrected(str(p)))


def test_umikey_is_accepted_as_the_join_key(tmp_path):
    """The UMI-collapsed cDNA arm has no read_id at all."""
    p = tmp_path / "umi.tsv"
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["umikey", "chrom", "strand", "corrected_3prime",
                    "alignment_start", "alignment_end"])
        w.writerow(["MOL1", "chrXV", "+", "307000", "305000", "307001"])
    rec = next(ra._iter_corrected(str(p)))
    assert rec["key"] == "MOL1"


def test_sidecar_header_records_unit_and_rule(tmp_path, genes):
    path = write_tsv(tmp_path, [read("r", 305000, 307000, 306999)])
    out = tmp_path / "out.tsv"
    n = ra.write_attribution_sidecar(path, str(out), genes, {},
                                     unit="molecules", control_labels=["wt1"])
    assert n == 1
    text = out.read_text()
    assert "# unit: molecules" in text
    assert "termination_rule:" in text
    assert "cpa_reference_from: wt1" in text
    # comment lines must precede the header so a naive reader can skip them
    body = [ln for ln in text.splitlines() if not ln.startswith("#")]
    assert body[0].split("\t") == ra.SIDECAR_HEADER


def test_ambiguity_range_passes_through_unmodified(tmp_path, genes):
    """Passed through for the consumer, never used for attribution.

    Its median over the reads in view is the assay's 3'-end resolution, which is
    the defensible way to pick a clustering window.
    """
    p = tmp_path / "amb.tsv"
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["read_id", "chrom", "strand", "corrected_3prime",
                    "alignment_start", "alignment_end", "ambiguity_range"])
        w.writerow(["r1", "chrXV", "+", "306999", "305000", "307000", "12"])
    rec = next(ra._iter_corrected(str(p)))
    assert rec["ambiguity_range"] == "12"
    assert ra.classify_read(rec, genes, {})["ambiguity_range"] == "12"


def test_ambiguity_range_absent_is_empty_not_an_error(tmp_path, genes):
    """The 11-column UMI-collapsed cDNA schema does not carry it."""
    path = write_tsv(tmp_path, [read("r", 305000, 307000, 306999)])
    rec = next(ra._iter_corrected(path))          # TSV_HEADER has no ambiguity_range
    assert rec["ambiguity_range"] == ""


def test_reference_ignores_genes_below_the_read_floor(tmp_path, genes):
    rows = [read(f"r{i}", 305000, 307000, 306999) for i in range(5)]
    path = write_tsv(tmp_path, rows)
    assert ra.build_cpa_reference([path], genes, min_reads=20) == {}
    assert "GENEA" in ra.build_cpa_reference([path], genes, min_reads=5)
