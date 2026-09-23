"""`rectify cdna-length-correct` end to end, and the bundled S. cerevisiae calibration.

The CLI test walks genes -> count -> panel -> fit -> apply -> ratio on a synthetic annotation and clusters.tsv
files with a planted length tilt; the bundled-data tests check the shipped files against their PROVENANCE and
recover a tilt planted on the bundled direct RNA reference.
"""
import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from rectify.cli import main
from rectify.core.cdna import length_bias as lb
from rectify.data import get_bundled_cdna_length_bias


def _run(argv):
    with pytest.raises(SystemExit) as e:
        main(["cdna-length-correct", *map(str, argv)])
    return e.value.code


# ── synthetic annotation and libraries ───────────────────────────────────────

N = 240
SPACING = 12_000


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = tmp_path_factory.mktemp("clc")
    rng = np.random.default_rng(21)
    cds_len = np.round(10 ** rng.uniform(np.log10(200), np.log10(6000), N)).astype(int)
    genes = [f"YSYN{i:03d}W" for i in range(N)]
    start0 = np.arange(N) * SPACING + 1000
    lines = ["##gff-version 3"]
    for g, s0, L in zip(genes, start0, cds_len):
        lines.append(f"chrI\tSGD\tCDS\t{s0 + 1}\t{s0 + L}\t.\t+\t0\tParent={g}_mRNA;Name={g}_CDS")
    lines.append("chrmt\tSGD\tCDS\t101\t400\t.\t+\t0\tParent=Q0001_mRNA;Name=Q0001_CDS")
    (d / "a.gff").write_text("\n".join(lines) + "\n")
    tm = pd.DataFrame({"gene": genes, "tx_lo": start0 + 1 - 60, "tx_hi": start0 + cds_len + 140})
    tm.to_csv(d / "models.tsv", sep="\t", index=False)
    tx = pd.Series(cds_len + 200.0, index=genes)
    ref = pd.Series(np.exp(rng.normal(5, 1.0, N)), index=genes)
    ref = ref / ref.sum() * 1e6
    tilts = {"libA": lambda x: -2.0 * (x - 3.0), "libB": lambda x: 1.0 * (x - 3.0)}
    for lib, t in tilts.items():
        lam = ref * 2.0 ** t(np.log10(tx))
        k = rng.poisson((lam / lam.sum() * 60_000).values)
        rows = []
        for g, s0, L, n in zip(genes, start0, cds_len, k):
            stop = s0 + L - 1
            rows += [("chrI", stop + 120, "fwd", 1)] * int(n)
        rows += [("chrI", start0[0] + 200, "fwd", 0)] * 500      # XF 0: dropped
        rows += [("chrI", 5, "fwd", 1)] * 50                       # before any stop: unassigned
        c = pd.DataFrame(rows, columns=["chrom", "anchor", "orient", "xf"])
        c.insert(0, "cluster_id", range(len(c)))
        c["n_reads"] = 3
        (d / lib / "analyze").mkdir(parents=True)
        c.to_csv(d / lib / "analyze" / "clusters.tsv.gz", sep="\t", index=False)
    panel = list(rng.choice(genes, 120, replace=False))
    (d / "panel.txt").write_text("\n".join(panel) + "\n")
    pd.DataFrame({"gene": genes, "cpm": ref.values}).to_csv(d / "ref.tsv", sep="\t", index=False)
    return d, genes, tx, ref, panel, tilts


def test_cli_genes_count_fit_apply_ratio(world):
    d, genes, tx, ref, panel, tilts = world
    assert _run(["genes", "--gff", d / "a.gff", "--transcript-models", d / "models.tsv", "-o", d / "genes.tsv"]) == 0
    gt = pd.read_csv(d / "genes.tsv", sep="\t").set_index("gene")
    assert "Q0001" not in gt.index and len(gt) == N              # mito dropped
    assert np.allclose(gt.loc[genes, "tx_len"], tx.loc[genes])

    assert _run(["count", d / "libA/analyze/clusters.tsv.gz", d / "libB/analyze/clusters.tsv.gz",
                 "--gene-table", d / "genes.tsv", "-o", d / "counts"]) == 0
    counts = pd.read_csv(d / "counts/counts.tsv", sep="\t", index_col=0)
    assert list(counts.columns) == ["libA", "libB"]              # names from <lib>/analyze/
    stats = pd.read_csv(d / "counts/count_stats.tsv", sep="\t").set_index("library")
    assert (stats.n_molecules - stats.n_xf_pass == 500).all()
    assert (stats.n_xf_pass - stats.n_assigned == 50).all()

    meta = pd.DataFrame({"library": ["libA", "libB"], "background": ["W303-AA", "BY4742"],
                         "medium": ["YPD", "YPD"], "condition": ["rapamycin", "rapamycin"]})
    meta.to_csv(d / "meta.tsv", sep="\t", index=False)
    (d / "scope.json").write_text(json.dumps({"background": ["W303-AA"], "medium": ["YPD"],
                                              "condition": ["rapamycin"]}))
    assert _run(["fit", "--counts", d / "counts/counts.tsv", "--gene-table", d / "genes.tsv",
                 "--panel", d / "panel.txt", "--reference", d / "ref.tsv", "--library-meta", d / "meta.tsv",
                 "--panel-scope", d / "scope.json", "-o", d / "fit"]) == 0
    qc = pd.read_csv(d / "fit/qc.tsv", sep="\t").set_index("library")
    assert qc.loc["libA", "slope_equiv"] == pytest.approx(-2.0, abs=0.25)
    assert qc.loc["libB", "slope_equiv"] == pytest.approx(1.0, abs=0.25)
    assert qc.certified.to_dict() == {"libA": True, "libB": False}
    assert "background 'BY4742' outside scope" in qc.loc["libB", "certified_reason"]
    prov = json.loads((d / "fit/fit_provenance.json").read_text())
    assert prov["scale_c"] == 1.0 and prov["cohort_curve"]["n_panel_used"] == 120
    nf = pd.read_csv(d / "fit/deseq2_normalization_factors.tsv.gz", sep="\t", index_col=0)
    np.testing.assert_allclose(np.exp(np.log(nf).mean(axis=1)).values, 1.0, atol=1e-6)

    # apply the stored curves to the same table: identical factors (as written, 8 significant digits)
    assert _run(["apply", "--params", d / "fit/params.tsv", "--counts", d / "counts/counts.tsv",
                 "--gene-table", d / "genes.tsv", "-o", d / "apply"]) == 0
    b_fit = pd.read_csv(d / "fit/bias_factors.tsv.gz", sep="\t", index_col=0)
    b_app = pd.read_csv(d / "apply/bias_factors.tsv.gz", sep="\t", index_col=0)
    pd.testing.assert_frame_equal(b_fit, b_app)

    # another gene-length definition is refused
    gt2 = pd.read_csv(d / "genes.tsv", sep="\t")
    gt2.loc[0, "tx_len"] += 5
    gt2.to_csv(d / "genes2.tsv", sep="\t", index=False)
    assert _run(["apply", "--params", d / "fit/params.tsv", "--counts", d / "counts/counts.tsv",
                 "--gene-table", d / "genes2.tsv", "-o", d / "apply2"]) != 0

    assert _run(["ratio", "--params", d / "fit/params.tsv", "--numerator-length", 2604,
                 "--denominator-length", 443, "--scale", 0.5, "-o", d / "ratio.tsv"]) == 0
    r = pd.read_csv(d / "ratio.tsv", sep="\t").set_index("library")
    assert (r.scale_c == 0.5).all()
    assert r.loc["libA", "log2_correction"] > 0 > r.loc["libB", "log2_correction"]


def test_cli_ratio_refuses_without_scale(world):
    d = world[0]
    assert _run(["ratio", "--params", d / "whatever.tsv", "--numerator-length", 2604,
                 "--denominator-length", 443]) == 2


def test_cli_count_contig_mismatch_is_an_error(world, tmp_path):
    d = world[0]
    c = pd.read_csv(d / "libA/analyze/clusters.tsv.gz", sep="\t")
    c["chrom"] = "NC_001133.9"
    c.to_csv(tmp_path / "clusters.tsv", sep="\t", index=False)
    assert _run(["count", tmp_path / "clusters.tsv", "--names", "x", "--gene-table", d / "genes.tsv",
                 "-o", tmp_path / "out"]) == 1


# ── bundled S. cerevisiae calibration ────────────────────────────────────────

def test_bundled_calibration_matches_its_provenance():
    paths = get_bundled_cdna_length_bias("yeast")
    assert paths is not None
    prov = json.loads(paths["provenance"].read_text())
    for key in ("gene_table", "panel", "reference"):
        p = paths[key]
        assert hashlib.sha256(p.read_bytes()).hexdigest() == prov["files"][p.name]["sha256"], p.name
    gt = lb.load_gene_table(paths["gene_table"])
    panel = pd.read_csv(paths["panel"], sep="\t")
    ref = pd.read_csv(paths["reference"], sep="\t").set_index("gene").cpm
    assert len(gt) == 5899 and len(panel) == 405 and int(panel.strict.sum()) == 127 and len(ref) == 5875
    assert set(panel.gene) <= set(gt.gene) and set(panel.gene) <= set(ref.index)
    assert ref.sum() == pytest.approx(1e6, rel=1e-6)
    L, _ = lb.resolve_lengths(gt, "tx")
    assert (L[panel.gene].min(), L[panel.gene].max()) == (314.0, 7346.0)
    assert prov["certified_scope"] == {"background": ["W303-AA"], "medium": ["YPD"], "condition": ["rapamycin"]}
    # the reporter-host loci are outside the gene universe by construction
    assert not {"YOL086C", "YLR044C", "YML017W"} & set(gt.gene)


def test_bundled_panel_and_reference_recover_a_planted_tilt(tmp_path):
    paths = get_bundled_cdna_length_bias("saccharomyces_cerevisiae")
    gt = lb.load_gene_table(paths["gene_table"])
    ref = pd.read_csv(paths["reference"], sep="\t").set_index("gene").cpm
    L, _ = lb.resolve_lengths(gt, "tx")
    rng = np.random.default_rng(5)
    genes = ref.index
    lam = ref * 2.0 ** (-3.0 * (np.log10(L[genes]) - 3.0))
    counts = pd.DataFrame({"short_lib": rng.poisson((lam / lam.sum() * 2e6).values)}, index=genes)
    counts.index.name = "gene"
    counts.to_csv(tmp_path / "counts.tsv", sep="\t")
    pd.DataFrame({"library": ["short_lib"], "background": ["W303-AA"], "medium": ["SD"],
                  "condition": ["rapamycin"]}).to_csv(tmp_path / "meta.tsv", sep="\t", index=False)
    assert _run(["fit", "--Scer", "--counts", tmp_path / "counts.tsv", "--library-meta", tmp_path / "meta.tsv",
                 "-o", tmp_path / "fit"]) == 0
    qc = pd.read_csv(tmp_path / "fit/qc.tsv", sep="\t").iloc[0]
    assert qc.slope_equiv == pytest.approx(-3.0, abs=0.2)
    assert qc.n_panel_used == 405
    assert not qc.certified and "medium 'SD' outside scope" in qc.certified_reason
    params = pd.read_csv(tmp_path / "fit/params.tsv", sep="\t")
    assert params.len_lo_nt.iloc[0] == pytest.approx(314.0) and params.len_hi_nt.iloc[0] == pytest.approx(7346.0)
