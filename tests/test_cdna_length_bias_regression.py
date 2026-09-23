"""Regression on the calibration cohort: 51 PCR-cDNA and 33 direct RNA libraries.

The count tables are unpublished and are NOT in the repository. Point ``RECTIFY_CDNA_LENGTH_FIXTURES`` at a
directory holding:

  counts_cdna.tsv.gz, counts_drs.tsv.gz  gene x library molecule counts, one counting rule for both assays
  gene_lengths.tsv                       the gene table the counts were made with (tx_len, cds_len)
  libraries.tsv                          assay, lib, set, geno for every library
  panel_candidates.tsv                   considered genes: cdna_mean, drs_mean, sd_sets, panel (the strict flag)
  matched_drs.tsv                        assay, lib, drs_ref: the genotype-matched direct RNA set of each cDNA library
  reference_params.tsv                   (optional) the reference analysis's per-library coefficients and knots

Without it every test here skips. With it, the test reproduces the held-out acceptance numbers of the correction,
checking "before" (which proves the gene split and groupings were rebuilt) separately from "after" (which proves
the fit), so a failure says which half broke.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from rectify.core.cdna import length_bias as lb

FIX = os.environ.get("RECTIFY_CDNA_LENGTH_FIXTURES")
NEEDED = ("counts_cdna.tsv.gz", "counts_drs.tsv.gz", "gene_lengths.tsv", "libraries.tsv",
          "panel_candidates.tsv", "matched_drs.tsv")
pytestmark = pytest.mark.skipif(
    not FIX or not all((Path(FIX) / f).exists() for f in NEEDED),
    reason="RECTIFY_CDNA_LENGTH_FIXTURES unset or incomplete (unpublished calibration tables)")

HELD_OUT_SEED = 1021     # the seed the reference analysis split its held-out genes with


@pytest.fixture(scope="module")
def cohort():
    F = Path(FIX)
    Mc = pd.read_csv(F / "counts_cdna.tsv.gz", sep="\t", index_col=0)
    Md = pd.read_csv(F / "counts_drs.tsv.gz", sep="\t", index_col=0)
    gt = pd.read_csv(F / "gene_lengths.tsv", sep="\t")
    GL = gt.set_index("gene")
    libs = pd.read_csv(F / "libraries.tsv", sep="\t")
    PV = pd.read_csv(F / "panel_candidates.tsv", sep="\t", index_col=0)
    matched = pd.read_csv(F / "matched_drs.tsv", sep="\t")
    assert Mc.shape[1] == 51 and Md.shape[1] == 33
    genes = Mc.index.intersection(GL.index)
    L_tx, _ = lb.resolve_lengths(gt, "tx")
    relax = PV.index[(PV.cdna_mean >= 100) & (PV.drs_mean >= 100) & (PV.sd_sets <= 0.5)]
    ref = lb.reference_from_counts(Md, genes)
    fr = lb.fit_length_bias(Mc, ref, sorted(relax), L_tx)

    # held-out genes: the considered genes outside the panel, split at random within CDS-length deciles
    outside = PV.index.difference(relax)
    rng = np.random.default_rng(HELD_OUT_SEED)
    bins = pd.qcut(np.log10(GL.loc[outside, "cds_len"].astype(float)), 10, labels=False)
    sel = pd.Index([g for b in range(10)
                    for g in rng.permutation(outside[bins.values == b])[: int((bins.values == b).sum()) // 2]])
    held = outside.difference(sel)
    return dict(Mc=Mc, Md=Md, GL=GL, libs=libs, PV=PV, matched=matched, genes=genes, L_tx=L_tx,
                relax=relax, ref=ref, fr=fr, held=held)


def _metrics(c, corrected):
    Mc, Md, GL, genes, held = c["Mc"], c["Md"], c["GL"], c["genes"], c["held"]
    meta = c["libs"][c["libs"].assay == "cDNA"].set_index("lib")
    dmeta = c["libs"][c["libs"].assay == "DRS"].set_index("lib")
    group = (meta.set.astype(str) + "|" + meta.geno.astype(str)).to_dict()
    lc = np.log2((lb.cpm(Mc.loc[genes]) if corrected is None else corrected) + 0.5)
    x = np.log10(GL.loc[genes, "cds_len"].astype(float))
    long_ = [g for g in held if GL.loc[g, "cds_len"] > 2000]
    cd = lb.cpm(Md.reindex(genes).fillna(0))

    def drs(key):
        if key.startswith("prjna:"):
            ls = dmeta.index[(dmeta.set == "prjna") & (dmeta.geno == key.split(":")[1])]
        else:
            ls = dmeta.index[dmeta.set == key]
        assert len(ls) == 3, key
        return np.log2(cd.loc[:, ls].mean(axis=1) + 0.5)

    e1 = lb.replicate_slopes(lc, group, x, held).abs().median()
    e3 = lb.replicate_spread(lc, group, long_).median()
    m = c["matched"][c["matched"].assay == "cDNA"].set_index("lib").drs_ref.dropna()
    e4 = np.median([abs(np.polyfit(x.loc[held], (lc.loc[held, l] - drs(k).loc[held]).values, 1)[0])
                    for l, k in m.items()])
    l1 = meta.index[(meta.set == "aug") & (meta.geno == "Rrp6-AA")]
    l0 = meta.index[(meta.set == "aug") & (meta.geno == "WT-AA")]
    fc_c = lc.loc[held, l1].mean(axis=1) - lc.loc[held, l0].mean(axis=1)
    fc_d = (drs("rrp6aa") - drs("wtaa")).loc[held]
    e5 = stats.pearsonr(fc_c, fc_d)[0]
    return dict(e1=e1, e3=e3, e4=e4, e5=e5)


def test_panel_is_rebuilt_from_the_counts(cohort):
    c = cohort
    d = c["libs"][c["libs"].assay == "DRS"]
    groups = dict(zip(d.lib, np.where(d.set == "prjna", "prjna:" + d.geno, d.set)))
    P = lb.build_panel(c["Mc"], c["Md"], groups, genes=c["GL"].index)
    assert set(P.index) == set(c["PV"].index) and len(P) == 3376
    assert set(P.index[P.panel]) == set(c["relax"]) and int(P.panel.sum()) == 405
    assert set(P.index[P.strict]) == set(c["PV"].index[c["PV"].panel]) and int(P.strict.sum()) == 127


def test_fit_reproduces_the_reference_coefficients(cohort):
    p = Path(FIX) / "reference_params.tsv"
    if not p.exists():
        pytest.skip("reference_params.tsv not supplied")
    ref = pd.read_csv(p, sep="\t", dtype=str).set_index("lib")
    for lib, curve in cohort["fr"].curves.items():
        want = np.array([float(v) for v in ref.loc[lib, "coefs"].split(";")])
        knots = np.array([float(v) for v in ref.loc[lib, "knots"].split(";")])
        np.testing.assert_allclose(curve.coef, want, rtol=1e-8, atol=1e-10)
        np.testing.assert_allclose(curve.knots, knots, rtol=0, atol=1e-9)
        assert (round(curve.len_lo_nt), round(curve.len_hi_nt)) == (314, 7346)
    assert (cohort["fr"].qc.n_genes_clamped == 36).all()


def test_held_out_numbers_before_correction(cohort):
    assert len(cohort["held"]) == 1488
    m = _metrics(cohort, None)
    assert m["e1"] == pytest.approx(1.20, abs=0.01)     # replicate gene-length slope, median |.|
    assert m["e3"] == pytest.approx(0.75, abs=0.01)     # replicate SD, genes > 2 kb
    assert m["e4"] == pytest.approx(2.28, abs=0.02)     # slope against genotype-matched direct RNA, median |.|
    assert m["e5"] == pytest.approx(-0.16, abs=0.005)   # Rrp6-AA / WT-AA fold changes vs direct RNA, r


def test_held_out_numbers_after_correction(cohort):
    m = _metrics(cohort, cohort["fr"].corrected_cpm)
    assert m["e1"] == pytest.approx(0.10, abs=0.01)
    assert m["e3"] == pytest.approx(0.26, abs=0.01)
    assert m["e4"] == pytest.approx(0.32, abs=0.01)
    assert m["e5"] == pytest.approx(0.47, abs=0.005)


def test_a_shared_reference_cancels_in_fold_changes(cohort):
    """Fitting against the cDNA cohort mean instead of direct RNA gives the same fold changes."""
    c = cohort
    own = lb.reference_from_counts(c["Mc"].loc[c["genes"]])
    fr2 = lb.fit_length_bias(c["Mc"], own, sorted(c["relax"]), c["L_tx"])
    a = _metrics(c, c["fr"].corrected_cpm)["e5"]
    b = _metrics(c, fr2.corrected_cpm)["e5"]
    assert b == pytest.approx(a, abs=0.002)


def test_only_the_certified_scope_is_certified(cohort):
    meta = cohort["libs"][cohort["libs"].assay == "cDNA"].copy()
    reporter = meta.set == "psp2"
    meta["background"] = np.where(reporter, "reporter strain", "W303-AA")
    meta["medium"] = np.where(reporter, "SD", "YPD")
    meta["condition"] = "rapamycin"
    scope = {"background": ["W303-AA"], "medium": ["YPD"], "condition": ["rapamycin"]}
    cert = lb.certify_libraries(list(meta.lib), meta.rename(columns={"lib": "library"}), scope)
    by = cert.set_index("library").certified
    assert int((~by).sum()) == 18 and not by[meta.lib[reporter]].any()
    assert by[meta.lib[~reporter]].all()
