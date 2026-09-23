"""rectify cdna-length-correct: per-library length-bias correction of ONT PCR-cDNA gene counts.

Each test is written against a defect it must catch:
  * a library's known tilt (linear or sigmoid) must be recovered and removed from held-out genes;
  * the curve is never extrapolated: outside the panel's length range a gene gets the EDGE factor and is counted;
  * log2(B) is the fitted curve at each gene's length, and DESeq2 factors have per-gene geometric mean 1;
  * stored params reproduce the factors bit for bit, knots are fixed at fit time, and a different gene-length
    definition is refused;
  * the ratio helper has NO default scale c;
  * the counting rule's window edges (-30 .. +500 nt around the stop), strand, XF filter and the contig-naming
    guard behave exactly as specified.
"""
import gzip
import json
import math

import numpy as np
import pandas as pd
import pytest

from rectify.core.cdna import length_bias as lb


# ── synthetic cohort ──────────────────────────────────────────────────────────

N_GENES = 3000
DEPTH = 4e6


def _genes(seed=7):
    rng = np.random.default_rng(seed)
    genes = [f"G{i:05d}" for i in range(N_GENES)]
    L = pd.Series(10 ** rng.uniform(np.log10(150), np.log10(9000), N_GENES), index=genes)
    ref = pd.Series(np.exp(rng.normal(4.0, 1.3, N_GENES)), index=genes)
    ref = ref / ref.sum() * 1e6
    # panel: well-expressed genes between 300 and 6,000 nt, so genes on both sides are clamped
    ok = (L > 300) & (L < 6000) & (ref > np.quantile(ref, 0.5))
    panel = list(rng.choice(np.array(genes)[ok.values], 400, replace=False))
    return genes, L, ref, panel


def _library(ref, L, tilt, seed):
    rng = np.random.default_rng(seed)
    lam = ref * 2.0 ** tilt(np.log10(L))
    lam = lam / lam.sum() * DEPTH
    return pd.Series(rng.poisson(lam.values), index=ref.index)


def _linear(x):
    return -2.5 * (x - 3.0)


def _sigmoid(x):
    return 3.0 / (1 + np.exp(-6 * (x - 3.2)))


@pytest.fixture(scope="module")
def cohort():
    genes, L, ref, panel = _genes()
    C = pd.DataFrame({"lin": _library(ref, L, _linear, 1), "sig": _library(ref, L, _sigmoid, 2),
                      "flat": _library(ref, L, lambda x: 0 * x, 3)})
    fr = lb.fit_length_bias(C, ref, panel, L)
    return genes, L, ref, panel, C, fr


def _slope(x, y):
    return float(np.polyfit(x, y, 1)[0])


# ── recovery ──────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("lib,tilt", [("lin", _linear), ("sig", _sigmoid)])
def test_known_tilt_is_recovered_and_flattened_on_held_out_genes(cohort, lib, tilt):
    genes, L, ref, panel, C, fr = cohort
    curve = fr.curves[lib]
    inside = (L >= curve.len_lo_nt) & (L <= curve.len_hi_nt)
    x = np.log10(L[inside])
    f = curve.log2_bias(L[inside])
    t = tilt(x.values)
    # the intercept is arbitrary; compare shapes
    assert np.max(np.abs((f - f.mean()) - (t - t.mean()))) < 0.12

    held = [g for g in genes if g not in set(panel) and inside[g] and C.loc[g, lib] >= 50]
    assert len(held) > 300
    xh = np.log10(L[held])
    before = np.log2(lb.cpm(C)[lib][held] + 0.5) - np.log2(ref[held] + 0.5)
    after = np.log2(fr.corrected_cpm[lib][held] + 0.5) - np.log2(ref[held] + 0.5)
    b0, b1 = _slope(xh, before), _slope(xh, after)
    assert abs(b0) > 1.0, "the planted tilt must be visible before correction"
    assert abs(b1) < 0.1, f"held-out genes still tilted after correction: {b1:+.3f} log2/decade"


def test_untilted_library_is_left_alone(cohort):
    genes, L, ref, panel, C, fr = cohort
    f = fr.curves["flat"].log2_bias(L)
    assert np.ptp(f) < 0.25


# ── never extrapolated ────────────────────────────────────────────────────────

def test_curve_is_flat_outside_the_panel_range_and_clamped_genes_are_counted(cohort):
    genes, L, ref, panel, C, fr = cohort
    lo, hi = L[panel].min(), L[panel].max()
    for lib, curve in fr.curves.items():
        assert curve.len_lo_nt == pytest.approx(lo, rel=1e-12)
        assert curve.len_hi_nt == pytest.approx(hi, rel=1e-12)
        edge_lo, edge_hi = curve.log2_bias([lo, hi])
        below, above = L[L < lo], L[L > hi]
        assert len(below) and len(above)
        # flat outside the range: every clamped gene gets the edge value (to the last few ulps, which differ
        # between BLAS/SIMD code paths for different array sizes; an extrapolated curve differs by ~0.4 log2)
        np.testing.assert_allclose(curve.log2_bias(below), edge_lo, rtol=0, atol=1e-12)
        np.testing.assert_allclose(curve.log2_bias(above), edge_hi, rtol=0, atol=1e-12)
        # factors of out-of-range genes are exactly the edge factors
        B = fr.factors[lib]
        np.testing.assert_allclose(np.log2(B[below.index].values), edge_lo, rtol=0, atol=1e-12)
        np.testing.assert_allclose(np.log2(B[above.index].values), edge_hi, rtol=0, atol=1e-12)
        row = fr.qc.set_index("library").loc[lib]
        assert row.n_clamped_below == len(below)
        assert row.n_clamped_above == len(above)
        assert row.n_genes_clamped == len(below) + len(above)


# ── factors and DESeq2 ───────────────────────────────────────────────────────

def test_log2_factor_is_the_curve_at_each_gene_length(cohort):
    genes, L, ref, panel, C, fr = cohort
    for lib, curve in fr.curves.items():
        np.testing.assert_allclose(np.log2(fr.factors[lib].values), curve.log2_bias(L[fr.factors.index]),
                                   rtol=0, atol=1e-12)


def test_corrected_cpm_is_cpm_over_b_renormalized(cohort):
    genes, L, ref, panel, C, fr = cohort
    raw = lb.cpm(C) / fr.factors
    np.testing.assert_allclose(fr.corrected_cpm.values, (raw / raw.sum() * 1e6).values, rtol=1e-12)
    np.testing.assert_allclose(fr.corrected_cpm.sum().values, 1e6, rtol=1e-12)


def test_deseq2_factors_have_unit_geometric_mean_per_gene(cohort):
    genes, L, ref, panel, C, fr = cohort
    nf = lb.deseq2_normalization_factors(C, fr.factors)
    np.testing.assert_allclose(np.exp(np.log(nf).mean(axis=1)).values, 1.0, rtol=0, atol=1e-12)
    # within a library the factors vary with B alone (the size factor is one number per library)
    ratio = nf["lin"] / fr.factors["lin"]
    assert ratio.std() / ratio.mean() > 0  # rescaled per gene ...
    r = (nf["lin"] / nf["flat"]) / (fr.factors["lin"] / fr.factors["flat"])
    np.testing.assert_allclose(r.values, r.values[0], rtol=1e-10)  # ... but the library contrast is B's


def test_deseq2_factors_reduce_to_median_of_ratios_without_bias():
    counts = pd.DataFrame({"a": [10, 20, 30, 40, 0], "b": [20, 40, 60, 80, 5], "c": [5, 10, 15, 20, 1]},
                          index=list("vwxyz"))
    B = pd.DataFrame(1.0, index=counts.index, columns=counts.columns)
    nf = lb.deseq2_normalization_factors(counts, B)
    # DESeq2 median-of-ratios size factors: 1, 2, 0.5 up to a common scale
    sf = nf.iloc[0]
    np.testing.assert_allclose((sf / sf["a"]).values, [1.0, 2.0, 0.5], rtol=1e-12)
    np.testing.assert_allclose(np.exp(np.log(nf).mean(axis=1)).values, 1.0, atol=1e-12)


# ── stored parameters ─────────────────────────────────────────────────────────

def test_stored_params_reproduce_factors_exactly(cohort, tmp_path):
    genes, L, ref, panel, C, fr = cohort
    lsha = lb.lengths_fingerprint(L)
    params = lb.params_to_frame(fr.curves, lengths_sha256=lsha)
    p = tmp_path / "params.tsv"
    params.to_csv(p, sep="\t", index=False)
    back = lb.read_params(p)
    B, Ccor, qc = lb.apply_length_curves(C, back, L)
    assert np.array_equal(B.loc[fr.factors.index, fr.factors.columns].values, fr.factors.values)
    np.testing.assert_array_equal(Ccor.values, fr.corrected_cpm.loc[Ccor.index, Ccor.columns].values)


def test_knots_and_range_are_fixed_at_fit_time(cohort):
    """Applying to a table of only long genes must not re-derive knots or the clamp range from it."""
    genes, L, ref, panel, C, fr = cohort
    params = lb.params_to_frame(fr.curves, lengths_sha256=lb.lengths_fingerprint(L))
    long_genes = L.index[L > 2000]
    B, _, _ = lb.apply_length_curves(C.loc[long_genes], params, L)
    # a subset table goes through a different matmul shape, so compare to ulp level, not bit for bit
    np.testing.assert_allclose(B.values, fr.factors.loc[B.index, B.columns].values, rtol=1e-13, atol=0)


def test_apply_refuses_a_different_gene_length_definition(cohort):
    genes, L, ref, panel, C, fr = cohort
    params = lb.params_to_frame(fr.curves, lengths_sha256=lb.lengths_fingerprint(L))
    L2 = L.copy()
    L2.iloc[0] += 1.0
    with pytest.raises(ValueError, match="gene-length definition differs"):
        lb.apply_length_curves(C, params, L2)


# ── ratio helper ──────────────────────────────────────────────────────────────

def test_ratio_helper_has_no_default_scale(cohort):
    curve = cohort[5].curves["lin"]
    with pytest.raises(TypeError):
        lb.ratio_log2_correction(curve, 2604, 443)
    with pytest.raises(TypeError):
        lb.ratio_log2_correction(curve, 2604, 443, c=None)
    for bad in (float("nan"), float("inf"), -0.5):
        with pytest.raises(ValueError):
            lb.ratio_log2_correction(curve, 2604, 443, c=bad)
    with pytest.raises(TypeError):
        lb.ratio_corrections({"lin": curve}, 2604, 443)


def test_ratio_correction_scales_with_c_and_flags_short_molecules(cohort):
    genes, L, ref, panel, C, fr = cohort
    curve = fr.curves["lin"]
    full = float(lb.ratio_log2_correction(curve, 2604, 443, c=1.0)[0])
    half = float(lb.ratio_log2_correction(curve, 2604, 443, c=0.5)[0])
    f = curve.log2_bias([2604, 443])
    assert full == pytest.approx(-(f[0] - f[1]), abs=1e-12)
    assert half == pytest.approx(0.5 * full, abs=1e-12)
    tab = lb.ratio_corrections(fr.curves, 443, 66, c=0.5)
    assert (tab.scale_c == 0.5).all()
    assert tab.denominator_clamped.all() and not tab.numerator_clamped.any()


# ── counting rule ─────────────────────────────────────────────────────────────

def _gene_table(chrom="chrI"):
    # PLUS: CDS 1000..1999 (stop 1999). MINUS: CDS 5000..5999 (stop 5000). PLUS2: stop 2999, 1 kb downstream of PLUS.
    return pd.DataFrame([
        ("PLUS", chrom, "+", 1000, 2000, 1000, 1000, 1200.0),
        ("PLUS2", chrom, "+", 2500, 3000, 500, 500, 700.0),
        ("MINUS", chrom, "-", 5000, 6000, 1000, 1000, 1300.0),
    ], columns=["gene", "chrom", "strand", "cds_start0", "cds_end", "cds_span", "cds_len", "tx_len"])


def _clusters(rows):
    return pd.DataFrame(rows, columns=["chrom", "orient", "anchor", "xf"])


def test_assignment_window_edges_and_strand():
    idx = lb.StopIndex(_gene_table())
    # -31 / -30 / 0 / +500 / +501 around PLUS's stop; 969 nt past it (beyond its window, before PLUS2's);
    # then -30 / -29 around PLUS2's stop, which is the NEAREST upstream stop from there on
    plus = idx.assign("chrI", "+", [1968, 1969, 1999, 2499, 2500, 2968, 2969, 2970])
    assert list(plus) == [None, "PLUS", "PLUS", "PLUS", None, None, "PLUS2", "PLUS2"]
    minus = idx.assign("chrI", "-", [5031, 5030, 5000, 4500, 4499])
    assert list(minus) == [None, "MINUS", "MINUS", "MINUS", None]
    # a + gene never takes a - strand end at the same place
    assert list(idx.assign("chrI", "-", [2100])) == [None]


def test_count_keeps_both_read_types_drops_xf0_and_counts_rows_not_reads():
    idx = lb.StopIndex(_gene_table())
    d = _clusters([
        ("chrI", "fwd", 2100, 1), ("chrI", "fwd", 2100, 2), ("chrI", "fwd", 2100, 0),
        ("chrI", "rev", 4900, 1), ("chrI", "rev", 4900, 3), ("chrI", "fwd", 9000, 1),
    ])
    d["xt"] = [1, 2, 1, 2, 1, 1]        # read type: Type 2 (no UMI) must be counted like Type 1
    d["n_reads"] = [50, 1, 7, 9, 2, 1]  # cluster sizes must NOT weight the count
    counts, st = lb.count_cdna_clusters(d, idx, min_assigned_share=None)
    assert counts.to_dict() == {"PLUS": 2, "MINUS": 2}
    assert st["n_molecules"] == 6 and st["n_xf_pass"] == 5 and st["n_assigned"] == 4


def test_contig_naming_is_normalized_and_a_mismatch_is_loud():
    idx = lb.StopIndex(_gene_table(chrom="I"))
    d = _clusters([("chrI", "fwd", 2100, 1)] * 5)
    counts, _ = lb.count_cdna_clusters(d, idx)
    assert counts.to_dict() == {"PLUS": 5}
    bad = _clusters([("ref|NC_001133|", "fwd", 2100, 1)] * 5)
    with pytest.raises(ValueError, match="contig-naming mismatch"):
        lb.count_cdna_clusters(bad, idx)


def test_excluded_gene_is_removed_from_the_candidate_stops(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chrI\tSGD\tCDS\t1001\t2000\t.\t+\t0\tParent=PLUS_mRNA;Name=PLUS_CDS\n"
        "chrI\tSGD\tCDS\t2201\t2400\t.\t+\t0\tParent=NEAR_mRNA;Name=NEAR_CDS\n")
    full = lb.StopIndex(lb.gene_table_from_gff(gff))
    excl = lb.StopIndex(lb.gene_table_from_gff(gff, exclude_genes=["NEAR"]))
    assert list(full.assign("chrI", "+", [2450])) == ["NEAR"]
    assert list(excl.assign("chrI", "+", [2450])) == ["PLUS"]  # 451 nt past PLUS's stop, inside its window


# ── gene table from a GFF ─────────────────────────────────────────────────────

GFF_TEXT = (
    "##gff-version 3\n"
    "chrI\tSGD\tgene\t101\t1000\t.\t+\t.\tID=YA;Name=YA\n"
    "chrI\tSGD\tCDS\t101\t300\t.\t+\t0\tParent=YA_mRNA;Name=YA_CDS\n"
    "chrI\tSGD\tintron\t301\t400\t.\t+\t.\tParent=YA_mRNA\n"
    "chrI\tSGD\tCDS\t401\t1000\t.\t+\t1\tParent=YA_mRNA;Name=YA_CDS\n"
    "chrI\tSGD\tCDS\t401\t1000\t.\t+\t1\tParent=YA_id001;Name=YA_CDS\n"
    "chrI\tSGD\tCDS\t2001\t2600\t.\t-\t0\tParent=YB_mRNA;Name=YB_CDS\n"
    "chrXII\tSGD\tCDS\t455001\t456000\t.\t+\t0\tParent=RDN_mRNA;Name=RDN_CDS\n"
    "chrmt\tSGD\tCDS\t101\t400\t.\t+\t0\tParent=Q1_mRNA;Name=Q1_CDS\n"
    "##FASTA\n>chrI\nACGTACGT\nchrI\tSGD\tCDS\t1\t9\t.\t+\t0\tName=GHOST_CDS\n"
)


@pytest.mark.parametrize("gz", [False, True])
def test_gene_table_from_gff(tmp_path, gz):
    p = tmp_path / ("a.gff.gz" if gz else "a.gff")
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(GFF_TEXT)
    else:
        p.write_text(GFF_TEXT)
    tm = pd.DataFrame({"gene": ["YA", "YB"], "tx_lo": [51, 1901], "tx_hi": [1100, 2650]})
    gt = lb.gene_table_from_gff(p, transcript_models=tm, exclude_regions=[("XII", 450000, 470000)]).set_index("gene")
    assert sorted(gt.index) == ["YA", "YB"]  # rDNA region, mito and the post-##FASTA row are gone
    ya = gt.loc["YA"]
    assert (ya.cds_start0, ya.cds_end, ya.cds_span, ya.cds_len) == (100, 1000, 900, 800)
    assert ya.tx_len == 1050 - 100  # 51..1100 minus the 100-nt intron
    assert gt.loc["YB", "strand"] == "-" and gt.loc["YB", "tx_len"] == 750


def test_resolve_lengths_fills_missing_transcripts_from_cds_plus_median_utr():
    gt = pd.DataFrame({"gene": list("abcd"), "cds_len": [1000, 2000, 500, 800],
                       "tx_len": [1200, 2300, np.nan, 1000]})
    L, info = lb.resolve_lengths(gt, "tx")
    assert L["c"] == 500 + 200 and info["n_filled"] == 1 and info["fill_utr_nt"] == 200
    assert lb.resolve_lengths(gt, "cds")[0]["c"] == 500
    with pytest.raises(ValueError, match="no transcript lengths"):
        lb.resolve_lengths(gt.assign(tx_len=np.nan), "tx")


# ── panel, certification, covariates ─────────────────────────────────────────

def test_build_panel_keeps_stable_high_count_genes_only():
    rng = np.random.default_rng(3)
    genes = [f"g{i}" for i in range(200)]
    base = pd.Series(rng.uniform(200, 2000, 200), index=genes)
    ref = {}
    for grp in range(4):
        shift = np.ones(200)
        shift[:20] = 2.0 ** (grp - 1.5)  # the first 20 genes change between condition groups
        for rep in range(3):
            ref[f"d{grp}_{rep}"] = rng.poisson(base * shift)
    R = pd.DataFrame(ref, index=genes)
    C = pd.DataFrame({f"c{i}": rng.poisson(base) for i in range(3)}, index=genes)
    groups = {f"d{g}_{r}": f"grp{g}" for g in range(4) for r in range(3)}
    P = lb.build_panel(C, R, groups, sd_max=0.5)
    assert not P.loc[genes[:20], "panel"].any()
    assert P.loc[genes[20:], "panel"].mean() > 0.95


def test_certification_is_scope_driven_and_defaults_to_flagged():
    scope = {"background": ["W303-AA"], "medium": ["YPD"], "condition": ["rapamycin"]}
    meta = pd.DataFrame({"library": ["ok", "sd", "blank"], "background": ["w303-aa", "W303-AA", ""],
                         "medium": ["YPD", "SD", "YPD"], "condition": ["Rapamycin", "rapamycin", "rapamycin"]})
    c = lb.certify_libraries(["ok", "sd", "blank", "absent"], meta, scope).set_index("library")
    assert c.certified.to_dict() == {"ok": True, "sd": False, "blank": False, "absent": False}
    assert "medium 'SD' outside scope" in c.loc["sd", "certified_reason"]
    assert "background not declared" in c.loc["blank", "certified_reason"]
    none = lb.certify_libraries(["ok"], meta, None)
    assert not none.certified.iloc[0] and "no certification scope" in none.certified_reason.iloc[0]


def test_optional_covariate_is_recovered_and_round_trips(tmp_path):
    genes, L, ref, panel, C, _ = (*_genes(), None, None)
    rng = np.random.default_rng(11)
    gc = pd.Series(rng.normal(0.4, 0.05, len(genes)), index=genes)
    z = (gc - gc[panel].mean()) / gc[panel].std(ddof=1)
    lib = _library(ref * 2.0 ** (0.3 * z), L, _linear, 5)
    cov = gc.to_frame("gc")
    fr = lb.fit_length_bias(lib.to_frame("x"), ref, panel, L, covariates=cov)
    curve = fr.curves["x"]
    assert curve.cov_names == ("gc",) and curve.coef[-1] == pytest.approx(0.3, abs=0.05)
    params = lb.params_to_frame(fr.curves, lengths_sha256=lb.lengths_fingerprint(L))
    B, _, _ = lb.apply_length_curves(lib.to_frame("x"), params, L, covariates=cov)
    assert np.array_equal(B.values, fr.factors.values)
    with pytest.raises(ValueError, match="uses covariates"):
        lb.apply_length_curves(lib.to_frame("x"), params, L)
    # off by default: no covariate columns in a plain fit
    plain = lb.fit_length_bias(lib.to_frame("x"), ref, panel, L)
    assert plain.curves["x"].cov_names == () and len(plain.curves["x"].coef) == 4
