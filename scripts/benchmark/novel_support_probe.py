#!/usr/bin/env python3
"""Novel-feature-support discovery-FDR guardrail — RUNNABLE PROTOTYPE.

Operationalizes the "novel-feature support" principle (ALIGNER_BENCHMARK_OVERVIEW.md
§"Read-quality structure & the 'novel-feature support' principle"; PI design note
#2 in dev/HANDOFF_ALIGNER_BENCHMARK.md):

    A REAL novel isoform should be supported by reads sampled from the dataset's
    OVERALL read-quality spectrum. If a dataset is ~p0 low-Q, a real novel call is
    ~p0-supported by low-Q reads. A putative novel call that is TAIL-ENRICHED
    (over-represented in low-Q reads — errors concentrated in the read segments
    that dictate the call) is SUSPECT: a per-read-reliability-weighted FDR control
    on discovery. This is a SOFT prior, not a hard gate.

This prototype reuses the EXISTING junction-discovery machinery —
``controlled.py``-style locus construction, ``error_injector.py`` per-read
over-dispersion (the ground-truth "hotness"/reliability label), minimap2 -ax
splice, and the ``scorer.py`` junction primitives (``extract_junctions`` /
``normalize_junction``) — and adds the missing piece: aggregate per-read junction
calls into per-COORDINATE novel CALLS, stratify each call's support by read
reliability, and flag tail-enriched calls with a support-count-aware binomial
statistic. Mirrors the ``c1_ablation.py`` / ``smoke_roundtrip.py`` harness style.

WHAT IS DEMONSTRATED (and the honest framing — pre-committed):
  * The metric + harness run end-to-end and SEPARATE genuine vs spurious novel
    calls (ROC/AUC) on BOTH a ground-truth reliability axis (the injector hot/cold
    label) AND a DEPLOYABLE reference-free proxy (each read's own exonic error
    density, estimated from its alignment WITHOUT junction truth).
  * The binomial tail-enrichment statistic BEATS naive frac_lowQ thresholding in
    the LOW-SUPPORT regime (where a genuine low-prevalence novel can swing high by
    chance) — the non-tautological content.
  * The genuine-vs-spurious separation is PARTLY BY CONSTRUCTION (we place spurious
    junctions in hot reads). What this proves is the METRIC + HARNESS and the FDR
    lift in the overlap zone — NOT that the principle holds on emergent data. The
    calibrated production FDR number is MAGNITUDE-SENSITIVE and SIRV-gated (the
    injector params are placeholders); it is OUT OF SCOPE here.

TOP OPEN RISK (handoff line 376, surfaced not waved through): on REAL data the
reference-free proxy conflates per-MOLECULE hotness with per-TRANSCRIPT
alignability, so down-weighting hot reads could suppress genuine novels from
hard-to-align transcripts (a discovery bias). In THIS sim the reliability label is
per-molecule by construction, so the axis is clean — real deployment must isolate
per-molecule hotness (within- vs across-transcript autocorr) before trusting the
lift.

Usage:  python scripts/benchmark/novel_support_probe.py --out /tmp/novel_support
Deterministic, M1-light (a few thousand reads through minimap2). Needs minimap2 +
samtools on PATH.

Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import math
import os
import random
import statistics
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import (  # noqa: E402
    extract_junctions, load_genome,
)
from rectify.core.consensus.chimeric_consensus import normalize_junction  # noqa: E402
from scripts.benchmark.sim.controlled import _rand_unique  # noqa: E402
from scripts.benchmark.sim.error_injector import (  # noqa: E402
    InjectorParams, inject, events_from_bam_read,
)
from scripts.benchmark.sim.transcript_model import TranscriptModel, Exon  # noqa: E402

EXON = 300
SPLICE_MIMIC_LEN = 60          # >=40bp GT..AG block -> minimap2 fabricates an intron
COORD_TOL = 8                  # bp; cluster nearby normalized junctions into one CALL


# ---------------------------------------------------------------------------
# Read-quality spectrum: the injector's per-read over-dispersion IS the hotness
# ---------------------------------------------------------------------------
# HOT (low-Q) reads carry genuine GLOBAL over-dispersion + bursts, so the full
# background error load lands on them — this is what makes the reference-free
# proxy (axis b) able to SEE them (exonic error density elevated). The mechanism
# probe established the coupling: over-dispersion -> exonic predicts junction-region
# error (r=0.955); burst-only does not. We therefore use gamma over-dispersion as
# the primary driver. PLACEHOLDER magnitudes (SIRV-gated).
def hot_params(base_rate: float) -> InjectorParams:
    return InjectorParams(
        base_rate=base_rate, read_rate_dist="gamma", gamma_shape=0.5,
        burst_on=True, hot_factor=8.0, p_hot_to_cold=0.2, p_cold_to_hot=0.05,
        indel_run_pmf=((1, 0.6), (2, 0.25), (3, 0.15)))


def cold_params(base_rate: float) -> InjectorParams:
    p = InjectorParams.null()
    p.base_rate = base_rate
    return p


# ---------------------------------------------------------------------------
# Corpus generation
# ---------------------------------------------------------------------------
@dataclass
class Locus:
    chrom: str
    kind: str                       # "genuine" | "spurious"
    truth_junc: Optional[Tuple[int, int]] = None   # normalized (start,end) for genuine


def build_corpus(out_dir: str, rng: random.Random,
                 n_genuine: int, n_spurious: int,
                 p_spectrum: float, hot_rate: float, cold_rate: float,
                 spurious_hot_purity: float,
                 support_choices: Tuple[int, ...]
                 ) -> Tuple[str, str, Dict[str, Locus], Dict[str, str]]:
    """Construct genuine (real canonical NNC novel intron) and spurious
    (splice-mimic deletion fabricated only/mostly in hot reads) loci.

    Returns (ref_fa, reads_fastq, loci_by_chrom, read_class) where read_class maps
    read_id -> "hot"/"cold" (the GROUND-TRUTH reliability label, axis a).
    """
    os.makedirs(out_dir, exist_ok=True)
    refs: Dict[str, str] = {}
    reads: List[Tuple[str, str]] = []
    loci: Dict[str, Locus] = {}
    read_class: Dict[str, str] = {}

    def emit(rid: str, clean: str, is_hot: bool):
        params = hot_params(hot_rate) if is_hot else cold_params(cold_rate)
        dirty, _ev = inject(clean, params, rng)
        reads.append((rid, dirty))
        read_class[rid] = "hot" if is_hot else "cold"

    # ---- GENUINE loci: a REAL canonical GT-AG novel (NNC) intron -----------
    # Canonical so minimap2 reliably CALLS it (non-canonical would be motif-snapped
    # away — stratum G — entangling this axis with the snapping axis). No annotation
    # supplied -> NNC (a genuine DISCOVERY). All reads carry the true junction; each
    # read's quality is drawn from the spectrum, so support ~ Binom(n, p_spectrum).
    for i in range(n_genuine):
        chrom = f"gen_{i:03d}"
        e1 = _rand_unique(EXON, rng)
        intron = "GT" + _rand_unique(196, rng) + "AG"
        e2 = _rand_unique(EXON, rng)
        contig = e1 + intron + e2
        i_start, i_end = EXON, EXON + len(intron)
        model = TranscriptModel(name=chrom, chrom=chrom, strand="+",
                                exons=[Exon(0, EXON), Exon(i_end, len(contig))],
                                genome_seq=contig)
        jt = model.junction_truths()          # NNC (no annotation)
        ns, ne = normalize_junction(jt[0].intron_start, jt[0].intron_end, contig)
        refs[chrom] = contig
        loci[chrom] = Locus(chrom, "genuine", (ns, ne))
        spliced = model.spliced_transcript()
        n = rng.choice(support_choices)
        for r in range(n):
            emit(f"{chrom}_r{r:03d}", spliced, rng.random() < p_spectrum)

    # ---- SPURIOUS loci: a splice-mimic deletion fabricated by HOT reads -----
    # Single-exon contig with a GT..AG block (>=40bp). A read that DELETES the block
    # makes minimap2 fabricate an intron at a FIXED, recurrent coordinate (verified:
    # 19/20 reads -> one normalized coord). The deletion is carried PREFERENTIALLY by
    # hot reads (spurious_hot_purity), with a small COLD leak so the call is not
    # 100% low-Q (realistic overlap). Hot carriers ALSO get the full background error
    # load, so axis (b) sees them as low-Q (the deletion rides on top of genuine
    # global hotness). Cold reads carry the FULL contig -> no junction.
    for i in range(n_spurious):
        chrom = f"spur_{i:03d}"
        left = _rand_unique(EXON, rng)
        block = "GT" + _rand_unique(SPLICE_MIMIC_LEN - 4, rng) + "AG"
        right = _rand_unique(EXON, rng)
        contig = left + block + right
        refs[chrom] = contig
        loci[chrom] = Locus(chrom, "spurious", None)
        deleted = left + right            # block removed -> fabricated intron
        full = contig                     # intact -> no junction
        n = rng.choice(support_choices)   # number of deletion-CARRIER reads
        # carriers (the support population for the spurious call)
        for r in range(n):
            is_hot = rng.random() < spurious_hot_purity
            emit(f"{chrom}_carry_r{r:03d}", deleted, is_hot)
        # a few intact reads from the spectrum (do not support the spurious call;
        # contribute to the locus's normal coverage). Kept small + same molecule.
        for r in range(max(2, n // 4)):
            emit(f"{chrom}_full_r{r:03d}", full, rng.random() < p_spectrum)

    ref_fa = os.path.join(out_dir, "ref.fa")
    reads_fq = os.path.join(out_dir, "reads.fastq")
    with open(ref_fa, "w") as fh:
        for nm, s in refs.items():
            fh.write(f">{nm}\n{s}\n")
    with open(reads_fq, "w") as fh:
        for rid, s in reads:
            fh.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")
    return ref_fa, reads_fq, loci, read_class


# ---------------------------------------------------------------------------
# Alignment (reuse the smoke pattern; --MD so the axis-b proxy can read mismatches)
# ---------------------------------------------------------------------------
def run_minimap2_splice(ref_fa: str, reads_fq: str, out_bam: str) -> None:
    p = subprocess.run(
        ["minimap2", "-ax", "splice", "-uf", "--eqx", "--MD", "-k", "14", "-t", "2",
         ref_fa, reads_fq], capture_output=True)
    if p.returncode:
        raise RuntimeError("minimap2 failed: " + p.stderr.decode()[:500])
    sort = subprocess.run(["samtools", "sort", "-o", out_bam], input=p.stdout,
                          capture_output=True)
    if sort.returncode:
        raise RuntimeError("samtools sort failed: " + sort.stderr.decode()[:300])
    subprocess.run(["samtools", "index", out_bam], check=True)


# ---------------------------------------------------------------------------
# Call aggregation: per-read junctions -> per-COORDINATE novel CALLS
# ---------------------------------------------------------------------------
@dataclass
class NovelCall:
    chrom: str
    rep: Tuple[int, int]                 # representative normalized (start,end)
    kind: str                            # "genuine" | "spurious" (truth label)
    support_reads: List[str] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.support_reads)


def aggregate_calls(bam_path: str, genome: Dict[str, str], loci: Dict[str, Locus],
                    coord_tol: int = COORD_TOL
                    ) -> Tuple[List[NovelCall], Dict[str, float], List[str]]:
    """Cluster placed-read junction calls into per-coordinate CALLS (tolerance
    window) and measure each PLACED read's exonic error density (axis-b proxy).

    Returns (calls, exonic_rate_by_read, placed_read_ids).
    """
    # per-chrom list of (norm_start, norm_end, read_id)
    by_chrom: Dict[str, List[Tuple[int, int, str]]] = {}
    exonic_rate: Dict[str, float] = {}
    placed: List[str] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            placed.append(r.query_name)
            seq = genome.get(r.reference_name, "")
            # axis-b reference-free reliability proxy: the read's OWN exonic error
            # density (no junction truth used) — estimable on real data.
            try:
                ev, span = events_from_bam_read(r, exonic_coords=True)
                exonic_rate[r.query_name] = (len(ev) / span) if span else 0.0
            except Exception:
                pass
            for cs, ce in extract_junctions(r.reference_start, r.cigartuples):
                if not seq:
                    continue
                ns, ne = normalize_junction(cs, ce, seq)
                by_chrom.setdefault(r.reference_name, []).append((ns, ne, r.query_name))

    calls: List[NovelCall] = []
    for chrom, items in by_chrom.items():
        loc = loci.get(chrom)
        items.sort()
        clusters: List[List[Tuple[int, int, str]]] = []
        for (s, e, rid) in items:
            placed_in = False
            for cl in clusters:
                rs, re, _ = cl[0]
                if abs(s - rs) <= coord_tol and abs(e - re) <= coord_tol:
                    cl.append((s, e, rid))
                    placed_in = True
                    break
            if not placed_in:
                clusters.append([(s, e, rid)])
        for cl in clusters:
            rep = (cl[0][0], cl[0][1])
            # truth label: genuine ONLY if on a genuine locus AND coord matches the
            # truth junction within tolerance; else spurious (incl. stray mis-placed
            # junctions on a genuine locus).
            kind = "spurious"
            if loc and loc.kind == "genuine" and loc.truth_junc:
                ts, te = loc.truth_junc
                if abs(rep[0] - ts) <= coord_tol and abs(rep[1] - te) <= coord_tol:
                    kind = "genuine"
            calls.append(NovelCall(chrom, rep, kind,
                                   support_reads=[rid for _, _, rid in cl]))
    return calls, exonic_rate, placed


# ---------------------------------------------------------------------------
# The tail-enrichment statistic
# ---------------------------------------------------------------------------
def binom_upper_tail(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binom(n, p). Exact (n is small here). Returns 1.0 for
    k<=0, clamps p to (eps, 1-eps)."""
    if k <= 0:
        return 1.0
    if n <= 0:
        return 1.0
    p = min(max(p, 1e-9), 1 - 1e-9)
    total = 0.0
    for i in range(k, n + 1):
        total += math.comb(n, i) * (p ** i) * ((1 - p) ** (n - i))
    return min(1.0, max(0.0, total))


def tail_enrichment(k: int, n: int, p0: float) -> float:
    """-log10 P(Binom(n,p0) >= k). Higher = more tail-enriched = more SUSPECT.
    This is the per-call read-reliability-weighted FDR prior."""
    return -math.log10(binom_upper_tail(k, n, p0))


def auc(scores_pos: List[float], scores_neg: List[float]) -> float:
    """AUC = P(score(positive) > score(negative)) (Mann-Whitney), ties=0.5.
    positive = SPURIOUS (the class the score should rank high)."""
    if not scores_pos or not scores_neg:
        return float("nan")
    wins = 0.0
    for sp in scores_pos:
        for sn in scores_neg:
            if sp > sn:
                wins += 1.0
            elif sp == sn:
                wins += 0.5
    return wins / (len(scores_pos) * len(scores_neg))


# ---------------------------------------------------------------------------
# Scoring an axis (a label assignment -> per-call k_j, TE, ROC, FDR lift)
# ---------------------------------------------------------------------------
@dataclass
class AxisResult:
    name: str
    p0: float
    n_genuine: int
    n_spurious: int
    auc_te: float
    auc_naive: float
    precision_before: float
    precision_after: float
    fdr_before: float
    fdr_after: float
    recall_retained: float
    te_threshold: float
    soft_precision: float
    # advisor #3: matched-operating-point false-alarm comparison in the low-support
    # regime (where genuine frac_lowQ swings high by chance). At the SAME number of
    # spurious calls caught, how many GENUINE calls does each guardrail wrongly flag?
    ls_n_genuine: int
    ls_n_spurious: int
    ls_spur_caught: int
    ls_gen_falsealarm_te: int
    ls_gen_falsealarm_naive: int
    ls_gen_fa_naive_lowsupp: int


def score_axis(name: str, calls: List[NovelCall], is_low_q: Dict[str, bool],
               p0: float, min_support: int, alpha: float,
               low_support_max: int) -> AxisResult:
    """Compute TE per call under one reliability axis, then ROC + FDR lift."""
    rows = []  # (call, n, k, frac_lowQ, TE)
    for c in calls:
        if c.n < min_support:
            continue
        k = sum(1 for rid in c.support_reads if is_low_q.get(rid, False))
        frac = k / c.n
        te = tail_enrichment(k, c.n, p0)
        rows.append((c, c.n, k, frac, te))

    gen = [r for r in rows if r[0].kind == "genuine"]
    spu = [r for r in rows if r[0].kind == "spurious"]

    te_thr = -math.log10(alpha)
    # accept a call unless flagged tail-enriched (TE > threshold)
    kept = [r for r in rows if r[4] <= te_thr]
    n_all = len(rows)
    n_gen = len(gen)
    prec_before = n_gen / n_all if n_all else float("nan")
    fdr_before = (n_all - n_gen) / n_all if n_all else float("nan")
    kept_gen = sum(1 for r in kept if r[0].kind == "genuine")
    prec_after = kept_gen / len(kept) if kept else float("nan")
    fdr_after = (len(kept) - kept_gen) / len(kept) if kept else float("nan")
    recall_retained = kept_gen / n_gen if n_gen else float("nan")

    # SOFT down-weight (honor "soft prior not hard gate"): weight each call by its
    # binomial p-value w = P(Binom>=k) in [0,1] (tail-enriched -> ~0). Soft
    # precision = sum(w over genuine) / sum(w over all).
    w = [(r[0].kind, binom_upper_tail(r[2], r[1], p0)) for r in rows]
    sw_all = sum(wi for _, wi in w)
    sw_gen = sum(wi for k_, wi in w if k_ == "genuine")
    soft_prec = sw_gen / sw_all if sw_all else float("nan")

    a_te = auc([r[4] for r in spu], [r[4] for r in gen])
    a_naive = auc([r[3] for r in spu], [r[3] for r in gen])

    # ---- advisor #3: matched-operating-point false-alarm test ----
    # The non-tautological content. A naive frac_lowQ threshold ignores support
    # count, so to catch MODERATELY tail-enriched spurious calls (frac ~ 2-3x p0,
    # the realistic enrichment with cold-read leakage) it must drop its cutoff —
    # which then ALSO sweeps up GENUINE low-prevalence novels whose frac swung high
    # by chance on few reads. The binomial TE is support-count-aware: the SAME frac
    # is significant on a high-support spurious call but NOT on a low-support genuine
    # one. Compare at the SAME spurious-catch count: TE should raise FEWER genuine
    # false alarms, and those it avoids are concentrated in the low-support regime.
    spur_caught_te = sum(1 for r in spu if r[4] > te_thr)
    gen_fa_te = sum(1 for r in gen if r[4] > te_thr)
    spu_fracs = sorted((r[3] for r in spu), reverse=True)
    if spur_caught_te > 0 and spu_fracs:
        idx = min(spur_caught_te, len(spu_fracs)) - 1
        naive_thr = spu_fracs[idx]
        gen_fa_naive = sum(1 for r in gen if r[3] >= naive_thr)
        gen_fa_naive_lowsupp = sum(1 for r in gen
                                   if r[3] >= naive_thr and r[1] <= low_support_max)
    else:
        gen_fa_naive = gen_fa_naive_lowsupp = 0

    return AxisResult(
        name=name, p0=p0, n_genuine=n_gen, n_spurious=len(spu),
        auc_te=a_te, auc_naive=a_naive,
        precision_before=prec_before, precision_after=prec_after,
        fdr_before=fdr_before, fdr_after=fdr_after,
        recall_retained=recall_retained, te_threshold=te_thr,
        soft_precision=soft_prec,
        ls_n_genuine=len(gen), ls_n_spurious=len(spu),
        ls_spur_caught=spur_caught_te,
        ls_gen_falsealarm_te=gen_fa_te, ls_gen_falsealarm_naive=gen_fa_naive,
        ls_gen_fa_naive_lowsupp=gen_fa_naive_lowsupp)


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="/tmp/novel_support")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--n-genuine", type=int, default=100)
    ap.add_argument("--n-spurious", type=int, default=100)
    ap.add_argument("--p-spectrum", type=float, default=0.15,
                    help="dataset low-Q fraction (the quality spectrum)")
    ap.add_argument("--spurious-hot-purity", type=float, default=0.85,
                    help="P(a spurious deletion-carrier is hot); <1 = cold leakage")
    ap.add_argument("--hot-rate", type=float, default=0.05)
    ap.add_argument("--cold-rate", type=float, default=0.008)
    ap.add_argument("--min-support", type=int, default=4)
    ap.add_argument("--alpha", type=float, default=0.01,
                    help="binomial tail alpha; TE threshold = -log10(alpha)")
    ap.add_argument("--low-support-max", type=int, default=15,
                    help="n_support <= this = the low-support regime (advisor #3)")
    args = ap.parse_args()
    rng = random.Random(args.seed)
    os.makedirs(args.out, exist_ok=True)
    # varied support depth, WEIGHTED toward low-prevalence (the overlap regime where
    # a genuine novel's frac_lowQ swings high by chance and the support-count-aware
    # TE earns its keep over naive frac thresholding).
    support_choices = (40, 30, 20, 12, 8, 6, 5, 5, 4, 4, 4)

    print("[novel-support] building corpus ...", file=sys.stderr)
    ref_fa, reads_fq, loci, read_class = build_corpus(
        args.out, rng, args.n_genuine, args.n_spurious, args.p_spectrum,
        args.hot_rate, args.cold_rate, args.spurious_hot_purity, support_choices)
    genome = load_genome(ref_fa)
    bam = os.path.join(args.out, "mm2.bam")
    print("[novel-support] minimap2 -ax splice ...", file=sys.stderr)
    run_minimap2_splice(ref_fa, reads_fq, bam)

    calls, exonic_rate, placed = aggregate_calls(bam, genome, loci)
    placed_set = set(placed)

    # ----- axis (a): GROUND-TRUTH injector hot/cold label -----
    label_a = {rid: (read_class.get(rid) == "hot") for rid in placed}
    # p0 = low-Q fraction among PLACED GENUINE-LOCUS reads (the dataset spectrum;
    # the spurious deletion-carriers are the anomaly tested AGAINST this null).
    gen_reads = [rid for rid in placed if rid.startswith("gen_")]
    p0_a = (sum(1 for rid in gen_reads if label_a.get(rid)) / len(gen_reads)
            if gen_reads else args.p_spectrum)

    # ----- axis (b): DEPLOYABLE reference-free proxy (exonic error density) -----
    # threshold so the flagged fraction among GENUINE-locus reads == p0_a (quantile
    # rule) -> the proxy's null matches the spectrum; this is the self-calibrating
    # per-read reliability covariate (no junction truth used).
    gen_rates = sorted(exonic_rate.get(rid, 0.0) for rid in gen_reads)
    if gen_rates:
        q_idx = max(0, min(len(gen_rates) - 1, int((1 - p0_a) * len(gen_rates))))
        tau = gen_rates[q_idx]
    else:
        tau = 0.02
    label_b = {rid: (exonic_rate.get(rid, 0.0) > tau) for rid in placed}
    p0_b = (sum(1 for rid in gen_reads if label_b.get(rid)) / len(gen_reads)
            if gen_reads else p0_a)

    res_a = score_axis("ground-truth hotness (axis a)", calls, label_a, p0_a,
                       args.min_support, args.alpha, args.low_support_max)
    res_b = score_axis(f"reference-free exonic-rate proxy (axis b, tau={tau:.4f})",
                       calls, label_b, p0_b, args.min_support, args.alpha,
                       args.low_support_max)

    # diagnostic: hot vs cold exonic rates (confirms axis-b sees the hotness)
    hot_r = [exonic_rate[r] for r in placed if read_class.get(r) == "hot" and r in exonic_rate]
    cold_r = [exonic_rate[r] for r in placed if read_class.get(r) == "cold" and r in exonic_rate]

    # ---------------- report ----------------
    n_calls = sum(1 for c in calls if c.n >= args.min_support)
    print("\n" + "=" * 78)
    print("NOVEL-FEATURE-SUPPORT DISCOVERY-FDR GUARDRAIL — prototype results")
    print("=" * 78)
    print(f"corpus: genuine_loci={args.n_genuine} spurious_loci={args.n_spurious} "
          f"reads_placed={len(placed_set)}  novel_calls(min_support>={args.min_support})={n_calls}")
    print(f"reliability proxy check: exonic error rate  hot median="
          f"{statistics.median(hot_r) if hot_r else float('nan'):.4f}  "
          f"cold median={statistics.median(cold_r) if cold_r else float('nan'):.4f}  "
          f"(axis-b threshold tau={tau:.4f})")

    def block(r: AxisResult):
        print(f"\n--- {r.name} ---")
        print(f"  p0 (dataset low-Q fraction, the binomial null) = {r.p0:.4f}")
        print(f"  calls scored: genuine={r.n_genuine}  spurious={r.n_spurious}")
        print(f"  AUC (tail-enrichment TE separates spurious from genuine) = {r.auc_te:.4f}")
        print(f"  AUC (naive frac_lowQ threshold)                           = {r.auc_naive:.4f}")
        print(f"  matched-operating-point false alarms (advisor #3): at equal spurious-catch "
              f"({r.ls_spur_caught}/{r.ls_n_spurious})")
        print(f"     genuine novels WRONGLY flagged:  TE={r.ls_gen_falsealarm_te}  "
              f"vs naive-frac={r.ls_gen_falsealarm_naive} "
              f"({r.ls_gen_fa_naive_lowsupp} of them low-support)  (TE should be <= naive)")
        print(f"  FDR lift (hard flag at TE> {r.te_threshold:.2f} = alpha {args.alpha}):")
        print(f"     precision  before={r.precision_before:.4f} -> after={r.precision_after:.4f}")
        print(f"     FDR        before={r.fdr_before:.4f} -> after={r.fdr_after:.4f}")
        print(f"     genuine recall retained = {r.recall_retained:.4f}")
        print(f"  SOFT down-weight precision (binomial p-value weights) = {r.soft_precision:.4f}")

    block(res_a)
    block(res_b)

    # advisor #3 concrete illustration: a low-support genuine vs a high-support spurious
    print("\n--- advisor #3 illustration (support-count-aware TE vs naive frac) ---")
    examples = []
    for c in calls:
        if c.n < args.min_support:
            continue
        k = sum(1 for rid in c.support_reads if label_a.get(rid, False))
        examples.append((c.kind, c.n, k, k / c.n, tail_enrichment(k, c.n, p0_a)))
    gen_lo = [e for e in examples if e[0] == "genuine" and e[1] <= args.low_support_max]
    spu_hi = [e for e in examples if e[0] == "spurious" and e[1] > args.low_support_max]
    gen_lo.sort(key=lambda e: -e[3]); spu_hi.sort(key=lambda e: -e[1])
    te_thr = -math.log10(args.alpha)
    if gen_lo:
        k_, n_, kk, fr, te = gen_lo[0]
        print(f"  genuine low-support: n={n_} k={kk} frac_lowQ={fr:.2f} TE={te:.2f} "
              f"-> {'FLAGGED' if te > te_thr else 'kept (correct)'}  "
              f"(naive frac>{p0_a:.2f} would {'FLAG (wrong)' if fr > p0_a else 'keep'})")
    if spu_hi:
        k_, n_, kk, fr, te = spu_hi[0]
        print(f"  spurious high-support: n={n_} k={kk} frac_lowQ={fr:.2f} TE={te:.2f} "
              f"-> {'FLAGGED (correct)' if te > te_thr else 'kept (miss)'}")

    print("\n" + "=" * 78)
    print("HONEST FRAMING (pre-committed):")
    print("  * SOFT prior, not a hard gate. Magnitude is PLACEHOLDER (injector params")
    print("    SIRV-gated); the separation/lift is a MECHANISM demo, not a calibrated")
    print("    production FDR number.")
    print("  * Separation is PARTLY BY CONSTRUCTION (spurious junctions placed in hot")
    print("    reads). Demonstrated value = the metric + harness + the lift in the")
    print("    overlap (low-support) zone, where binomial TE beats naive frac.")
    print("  * TOP OPEN RISK (real data): the reference-free proxy conflates per-MOLECULE")
    print("    hotness with per-TRANSCRIPT alignability; down-weighting could suppress")
    print("    genuine novels from hard-to-align transcripts (discovery bias). In sim the")
    print("    label is per-molecule by construction (clean); real use must isolate it.")
    print("=" * 78)

    # exit non-zero only if the harness is broken (no calls of either class)
    ok = (res_a.n_genuine > 0 and res_a.n_spurious > 0
          and not math.isnan(res_a.auc_te))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
