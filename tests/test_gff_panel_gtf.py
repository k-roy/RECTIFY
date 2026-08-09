"""Exon-feature GTF loader tests for the Tier-2 panel (``gff_panel``).

SIRV-Set 4 / Sequins / ERCC / LRGASP-GENCODE are exon-feature GTFs (one
transcript = a set of ``exon`` rows sharing ``transcript_id``; NO ``intron``
feature). ``build_panel_from_gtf`` must produce the SAME
``(models, pairs, donors, acceptors)`` tuple the GFF path produces.

The load-bearing invariant (advisor): build models from a GTF, classify every
model's junctions against the catalogue derived from THAT SAME GTF -> all come
back ANNOTATED (non-ANNOTATED == 0). That single check validates the coordinate
convention AND the intron-from-exon-gaps derivation simultaneously.

Hermetic: a synthetic exon-GTF is written inline (no network / no real SIRV
fetch) so the regression suite has no external dependency.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scripts.benchmark.sim.gff_panel import (  # noqa: E402
    _gtf_attrs, parse_gtf_exons, build_panel_from_gtf,
)
from rectify.core.benchmark.truth_schema import JunctionClass  # noqa: E402


# --- synthetic genome + exon-GTF ----------------------------------------
# One contig "chrT". Intron motifs are placed explicitly so canonicity is
# deterministic: plus-strand introns are GT..AG; the minus-strand intron's
# GENOME motif is CT..AC (== revcomp GT-AG), canonical when strand-aware.
def _build_genome():
    seq = list("C" * 400)

    def put(i, s):
        for k, ch in enumerate(s):
            seq[i + k] = ch

    # plus tx introns: [20,40) and [60,90)
    put(20, "GT"); put(38, "AG")
    put(60, "GT"); put(88, "AG")
    # minus tx intron: [170,200) genome motif CT..AC
    put(170, "CT"); put(198, "AC")
    return {"chrT": "".join(seq)}


# exon rows are deliberately NOT pre-sorted and use 1-based inclusive coords.
GTF = """\
chrT\ttest\texon\t41\t60\t.\t+\t.\tgene_id "G1"; transcript_id "tx_plus.1";
chrT\ttest\texon\t1\t20\t.\t+\t.\tgene_id "G1"; transcript_id "tx_plus.1";
chrT\ttest\texon\t91\t110\t.\t+\t.\tgene_id "G1"; transcript_id "tx_plus.1";
chrT\ttest\texon\t151\t170\t.\t-\t.\tgene_id "G2"; transcript_id "tx_minus.2";
chrT\ttest\texon\t201\t220\t.\t-\t.\tgene_id "G2"; transcript_id "tx_minus.2";
chrT\ttest\texon\t301\t360\t.\t+\t.\tgene_id "G3"; transcript_id "ercc_mono.3";
"""


def _write_gtf(tmp_path):
    p = os.path.join(str(tmp_path), "synth.gtf")
    with open(p, "w") as fh:
        fh.write(GTF)
    return p


# --- attribute parser ---------------------------------------------------
def test_gtf_attrs_quoted_space_separated():
    a = _gtf_attrs('gene_id "G1"; transcript_id "tx_plus.1";')
    assert a["gene_id"] == "G1"
    assert a["transcript_id"] == "tx_plus.1"  # version retained verbatim


# --- parse_gtf_exons: coords, intron derivation, span -------------------
def test_parse_gtf_exons_introns_and_span(tmp_path):
    gtf = _write_gtf(tmp_path)
    mrnas, introns = parse_gtf_exons(gtf)
    # span = min exon start (0-based) .. max exon end
    assert mrnas["tx_plus.1"] == ("chrT", "+", 0, 110)
    assert mrnas["tx_minus.2"] == ("chrT", "-", 150, 220)
    # introns derived from sorted adjacent-exon gaps, 0-based half-open
    assert introns["tx_plus.1"] == [(20, 40), (60, 90)]
    assert introns["tx_minus.2"] == [(170, 200)]
    # monoexonic transcript has no intron entry
    assert "ercc_mono.3" not in introns


# --- the load-bearing self-consistency invariant ------------------------
def test_all_model_junctions_classify_annotated(tmp_path):
    gtf = _write_gtf(tmp_path)
    genome = _build_genome()
    models, pairs, donors, acceptors = build_panel_from_gtf(gtf, genome)
    # spliced_only default drops the monoexonic ERCC
    assert sorted(m.name for m in models) == ["tx_minus.2", "tx_plus.1"]
    n_anno = n_other = 0
    for m in models:
        for j in m.junction_truths(pairs, donors, acceptors):
            if j.klass is JunctionClass.ANNOTATED:
                n_anno += 1
            else:
                n_other += 1
    assert n_other == 0, "every model junction must be ANNOTATED against its own GTF"
    assert n_anno == 3  # 2 introns (plus) + 1 (minus)


def test_strand_and_canonicity_round_trip(tmp_path):
    gtf = _write_gtf(tmp_path)
    genome = _build_genome()
    models, pairs, donors, acceptors = build_panel_from_gtf(gtf, genome)
    by_name = {m.name: m for m in models}
    plus = by_name["tx_plus.1"]
    assert plus.strand == "+"
    pj = plus.junction_truths(pairs, donors, acceptors)
    assert [(j.intron_start, j.intron_end) for j in pj] == [(20, 40), (60, 90)]
    assert all(j.canonical for j in pj), [j.motif for j in pj]  # GT..AG
    # minus-strand intron canonical via strand-aware motif (genome CT..AC)
    minus = by_name["tx_minus.2"]
    assert minus.strand == "-"
    mj = minus.junction_truths(pairs, donors, acceptors)
    assert len(mj) == 1 and mj[0].canonical, mj[0].motif


def test_include_intronless_keeps_monoexonic(tmp_path):
    gtf = _write_gtf(tmp_path)
    genome = _build_genome()
    models, *_ = build_panel_from_gtf(gtf, genome, spliced_only=False)
    assert "ercc_mono.3" in {m.name for m in models}


def test_missing_contig_skipped_not_crash(tmp_path):
    gtf = _write_gtf(tmp_path)
    # genome lacking chrT -> every transcript skipped, no crash, empty panel
    models, *_ = build_panel_from_gtf(gtf, {"other": "ACGT" * 50})
    assert models == []


if __name__ == "__main__":
    import tempfile

    class _TP:
        def __init__(self, d):
            self._d = d

        def __str__(self):
            return self._d

    with tempfile.TemporaryDirectory() as d:
        tp = _TP(d)
        test_gtf_attrs_quoted_space_separated()
        test_parse_gtf_exons_introns_and_span(tp)
        test_all_model_junctions_classify_annotated(tp)
        test_strand_and_canonicity_round_trip(tp)
        test_include_intronless_keeps_monoexonic(tp)
        test_missing_contig_skipped_not_crash(tp)
    print("all gff_panel GTF tests passed")
