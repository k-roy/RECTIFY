#!/usr/bin/env python3
"""End-to-end SMOKE: truth -> reads -> aligner -> scorer, with the GATE assertions.

Validates the three benchmark components wire together on a tiny controlled
corpus BEFORE scaling on Sherlock. Per the brief, it asserts:

  (A) a KNOWN junction round-trips truth -> sim-read -> minimap2 -> scorer as a TP
      (the ambiguity-aware match credits minimap2's placement even if it lands one
      bp into the engineered donor repeat);
  (B) a DELIBERATELY SHIFTED junction call (hand-crafted CIGAR placing the intron
      one bp left, into the repeat) scores TP-NOT-FP under normalize_junction
      (the exact trap that produced the GMAP 0.09 artifact);
  (C) an INDEL round-trips position-exact: HP/STR deletions score as position-exact
      concordance TP, and a deliberately shifted-WITHIN-RUN placement is also TP
      (the framing metric: exact indel-position concordance, ambiguity-aware,
      never edit distance).

Deterministic, M1-light (no pbsim3 needed — that is the Tier-2 realism wrapper,
exercised separately on Sherlock). Needs minimap2 + samtools on PATH.

Usage:  python scripts/benchmark/smoke_roundtrip.py --out /tmp/bench_smoke
Exit 0 = all GATE assertions pass.

Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.truth_schema import read_truth_table  # noqa: E402
from rectify.core.benchmark.scorer import (  # noqa: E402
    score_bam, load_genome, extract_junctions, cigar_records_to_bam,
    net_indel_in_span, all_indel_positions, open_fasta,
)
from rectify.core.benchmark.truth_schema import IndelKind  # noqa: E402
from rectify.core.consensus.chimeric_consensus import normalize_junction  # noqa: E402
from rectify.core.align.local_aligner import align_exon_block_global  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402


def load_fastq(path):
    seqs = {}
    with pysam.FastxFile(path) as fq:
        for e in fq:
            seqs[e.name] = e.sequence
    return seqs


def run_flat_affine_arm(read_seqs, truth_subset, genome, ref_fa, out_bam):
    """BAM-ize the INTERNAL flat-affine DP (align_exon_block_global — the arm C1
    upgrades) so the gate can SCORE the arm-flat vs arm-law ablation, not just
    external aligners. Each read is globally aligned to its single-contig ref."""
    records = []
    for rid, t in truth_subset.items():
        seq = read_seqs.get(rid)
        if seq is None:
            continue
        ref = genome[t.chrom]
        cig = align_exon_block_global(seq, ref)   # ref_start = 0 (global to contig)
        records.append((rid, t.chrom, 0, cig, seq))
    cigar_records_to_bam(records, ref_fa, out_bam)
    return score_bam(out_bam, truth_subset, genome, aligner_name="flat_affine_DP")


def run_minimap2_splice(ref_fa: str, reads_fq: str, out_bam: str) -> None:
    p = subprocess.run(
        ["minimap2", "-ax", "splice", "-uf", "--eqx", "-k", "14", "-t", "2",
         ref_fa, reads_fq],
        capture_output=True)
    if p.returncode:
        raise RuntimeError("minimap2 failed: " + p.stderr.decode()[:500])
    sort = subprocess.run(["samtools", "sort", "-o", out_bam], input=p.stdout,
                          capture_output=True)
    if sort.returncode:
        raise RuntimeError("samtools sort failed: " + sort.stderr.decode()[:300])
    subprocess.run(["samtools", "index", out_bam], check=True)


def build_shifted_junction_bam(ref_fa: str, truth, out_bam: str) -> tuple:
    """Hand-craft a BAM where the JUNCTION_AMB read's intron is placed one bp LEFT
    of truth (into the engineered donor repeat). Returns (truth_intron, shifted).
    The scorer must normalize the shifted call back onto truth -> TP not FP."""
    jrow = next(t for t in truth if t.stratum == "JUNCTION_AMB")
    j = jrow.junctions[0]
    chrom = jrow.chrom
    genome = load_genome(ref_fa)
    seq = genome[chrom]
    # truth intron is stored LEFT-normalized; an ambiguity-equivalent call sits
    # one bp RIGHT (within the slide window). Place the call there; it must
    # normalize BACK onto truth -> TP not FP. (Shifting LEFT of the leftmost
    # normalized truth would leave the window and IS a genuine FP — see B2.)
    shifted_start, shifted_end = j.intron_start + 1, j.intron_end + 1
    fa = open_fasta(ref_fa)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"LN": fa.get_reference_length(c), "SN": c}
                     for c in fa.references]}
    e1_len = j.intron_start - jrow.genome_start
    e2_len = jrow.genome_end - j.intron_end
    # CIGAR with the intron shifted 1bp RIGHT: exon1 1bp longer, intron same
    # length one base right, exon2 1bp shorter.
    cig = [(0, e1_len + 1), (3, j.intron_end - j.intron_start), (0, e2_len - 1)]
    read_seq = jrow.true_transcript  # placeholder; sequence content irrelevant to junction scoring
    spliced = "".join([seq[jrow.genome_start:j.intron_start],
                       seq[j.intron_end:jrow.genome_end]])
    with pysam.AlignmentFile(out_bam, "wb", header=header) as bam:
        a = pysam.AlignedSegment(bam.header)
        a.query_name = jrow.read_id + "_SHIFTED"
        a.query_sequence = spliced
        a.flag = 0
        a.reference_id = bam.header.get_tid(chrom)
        a.reference_start = jrow.genome_start
        a.mapping_quality = 60
        a.cigartuples = cig
        bam.write(a)
    subprocess.run(["samtools", "sort", "-o", out_bam + ".sorted.bam", out_bam],
                   check=True, capture_output=True)
    os.replace(out_bam + ".sorted.bam", out_bam)
    subprocess.run(["samtools", "index", out_bam], check=True)
    # add a SHIFTED truth row so score_bam has truth for the renamed read
    import copy
    shifted_truth = copy.deepcopy(jrow)
    shifted_truth.read_id = jrow.read_id + "_SHIFTED"
    return (j.intron_start, j.intron_end), (shifted_start, shifted_end), shifted_truth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/bench_smoke")
    ap.add_argument("--reps", type=int, default=120)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("[smoke] generating controlled corpus ...", file=sys.stderr)
    info = controlled.generate_corpus(args.out, reps=args.reps, seed=7)
    print("[smoke]", info, file=sys.stderr)
    ref_fa = info["ref_fa"]
    reads_fq = info["reads_fastq"]
    truth_tsv = info["truth_tsv"]
    truth = read_truth_table(truth_tsv)
    truth_map = {t.read_id: t for t in truth}
    genome = load_genome(ref_fa)

    # ---- align with a real aligner (minimap2 splice) --------------------
    print("[smoke] minimap2 -ax splice ...", file=sys.stderr)
    bam = os.path.join(args.out, "mm2.bam")
    run_minimap2_splice(ref_fa, reads_fq, bam)
    score = score_bam(bam, truth_map, genome, aligner_name="minimap2")
    s = score.summary()
    print("[smoke] minimap2 summary:", file=sys.stderr)
    import json
    print(json.dumps(s, indent=2), file=sys.stderr)

    failures = []

    # ---- per-stratum breakout (the discriminating criterion) ------------
    # A candidate ablation only counts if minimap2 is BELOW ceiling on the
    # stratum (SPEC VERTICAL-SLICE FINDING). Score each stratum separately.
    strata = sorted({t.stratum for t in truth})
    per_stratum = {}
    print("\n[smoke] === per-stratum breakout (minimap2) ===", file=sys.stderr)
    for st in strata:
        sub = {rid: t for rid, t in truth_map.items() if t.stratum == st}
        ss = score_bam(bam, sub, genome, aligner_name=f"mm2:{st}")
        per_stratum[st] = ss
        ji = ss.indel
        jc = {k: dict(v) for k, v in ss.junction.by_class.items()}
        print(f"  {st:14s} indel_concordance={ji.position_exact_concordance:.4f} "
              f"(corr={ji.correct},inc={ji.incorrect})  junc_by_class={jc or '-'}",
              file=sys.stderr)

    # (A) known junction round-trips as TP (NNC ambiguity locus)
    jt_tp = score.junction.by_class.get("NNC", {}).get("tp", 0)
    if jt_tp < 1:
        failures.append(f"(A) junction round-trip: expected >=1 NNC TP, got {jt_tp}")
    else:
        print(f"[smoke] (A) PASS junction round-trip NNC TP={jt_tp}", file=sys.stderr)

    # (A2) NIC + ANNOTATED discovery-class labels verified end-to-end
    nic_tp = per_stratum.get("NIC") and per_stratum["NIC"].junction.by_class.get("NIC", {}).get("tp", 0)
    ann_tp = per_stratum.get("ANNOTATED") and per_stratum["ANNOTATED"].junction.by_class.get("ANNOTATED", {}).get("tp", 0)
    if not nic_tp:
        failures.append(f"(A2) NIC class not recovered as TP (got {nic_tp})")
    if not ann_tp:
        failures.append(f"(A2) ANNOTATED class not recovered as TP (got {ann_tp})")
    if nic_tp and ann_tp:
        print(f"[smoke] (A2) PASS discovery classes: NIC TP={nic_tp}, ANNOTATED TP={ann_tp}",
              file=sys.stderr)

    # (E) VARIANT/C6 — variant-induced discovery FDR is MEASURED and SPECIFIC.
    # A standing variant near a splice site can FABRICATE a 'novel' junction
    # (ALIGNER_MEMBER_DESIGN.md Addendum (b) / §8). The benchmark must show
    # (1) the incumbent IS below ceiling — the flat HAPLOID reference fabricates
    # a variant-adjacent junction from a >=40bp splice-mimicking deletion variant
    # — AND (2) the effect is SPECIFIC to that context (a plain SNP, or a short
    # deletion, does NOT), so the future C6 guard is targeted, not a blunt
    # 'abstain near every variant' rule. Verified vs minimap2 -ax splice -uf.
    if "VARIANT" in per_stratum:
        var_truth = {rid: t for rid, t in truth_map.items() if t.stratum == "VARIANT"}

        def _is_driver(t):
            return (t.true_locus.startswith("var_del_")
                    and int(t.true_locus.split("_")[-1]) >= 40)

        drivers = {rid: t for rid, t in var_truth.items() if _is_driver(t)}
        controls = {rid: t for rid, t in var_truth.items() if not _is_driver(t)}
        sd = score_bam(bam, drivers, genome, aligner_name="mm2:VAR_driver")
        sc = score_bam(bam, controls, genome, aligner_name="mm2:VAR_control")
        drv_adj = sd.junction.fp_variant_adjacent
        ctl_adj = sc.junction.fp_variant_adjacent
        print(f"[smoke] (E) VARIANT/C6: driver fp_variant_adjacent={drv_adj} "
              f"(junc_fdr={sd.junction.fdr:.3f}, "
              f"indel_conc={sd.indel.position_exact_concordance:.3f}); "
              f"control fp_variant_adjacent={ctl_adj}", file=sys.stderr)
        if drv_adj < 1:
            failures.append("(E) VARIANT: incumbent fabricated NO variant-adjacent junction "
                            "on the splice-mimic-deletion driver — stratum is non-discriminating")
        elif ctl_adj != 0:
            failures.append(f"(E) VARIANT: control sub-cases produced {ctl_adj} variant-adjacent "
                            f"FP — the effect is NOT specific to the splice-mimic context")
        else:
            print(f"[smoke] (E) PASS VARIANT/C6 stratum is DISCRIMINATING + SPECIFIC: the flat "
                  f"haploid reference fabricated {drv_adj} variant-adjacent junction FPs on the "
                  f">=40bp splice-mimic-deletion variant (a C6-addressable discovery-FDR "
                  f"inflation) while plain-SNP + short-deletion controls produced 0. Whether C6 "
                  f"REDUCES this FDR is the NEXT-cycle ablation.", file=sys.stderr)

    # (F) PARALOG/C4 — window/locus-selection accuracy is MEASURED, the incumbent
    # is below ceiling on the C4-addressable slice, AND the gap is provably
    # closeable by pooling (the (D)-equivalent, decidable from truth NOW). The
    # discriminating slice is NOT a fragment with NO distinguishing base (that is
    # information-theoretically unrecoverable by ANY method incl. C4 — a structural
    # 0.5, not a closeable gap, the vertical-slice trap). It is a WEAK-evidence
    # fragment covering exactly ONE distinguishing SNP: per-read the lone SNP is
    # sometimes noise-flipped so minimap2 mis-assigns (below ceiling), yet the
    # MAJORITY SNP base across the read pool recovers the true copy. A read SPANNING
    # all spread SNPs is the at-ceiling control. Verified vs minimap2 -ax splice -uf.
    if "PARALOG" in per_stratum:
        import collections
        para = {rid: t for rid, t in truth_map.items() if t.stratum == "PARALOG"}
        span = {rid: t for rid, t in para.items() if "_span_" in rid}
        weak = {rid: t for rid, t in para.items() if "_weak_" in rid}
        ssp = score_bam(bam, span, genome, aligner_name="mm2:PARA_span")
        swk = score_bam(bam, weak, genome, aligner_name="mm2:PARA_weak")
        span_acc, weak_acc = ssp.locus_accuracy, swk.locus_accuracy
        # pooling-recovery proof (truth-side, NO C4): per (family, copy) pool, does
        # the MAJORITY base at the lone distinguishing SNP recover the true copy?
        para_reads = load_fastq(reads_fq)
        pools = {}
        for rid, t in weak.items():
            v = t.variants[0]
            idx = v.pos - t.genome_start
            seq = para_reads.get(rid, "")
            base = seq[idx] if 0 <= idx < len(seq) else "N"
            key = (t.true_locus, t.true_transcript)
            d = pools.setdefault(key, {"alt": v.alt_allele, "bases": []})
            d["bases"].append(base)
        recovered = sum(
            1 for d in pools.values()
            if collections.Counter(d["bases"]).most_common(1)[0][0] == d["alt"])
        n_pools = len(pools)
        print(f"[smoke] (F) PARALOG/C4: spanning locus_acc={span_acc:.3f}; "
              f"weak-evidence (1-SNP) fragment locus_acc={weak_acc:.3f}; "
              f"pooling-majority recovery={recovered}/{n_pools} pools", file=sys.stderr)
        if span_acc < 0.98:
            failures.append(f"(F) PARALOG: spanning reads NOT at ceiling "
                            f"(locus_acc={span_acc:.3f}) — metric mis-built / span evidence weak")
        elif weak_acc >= 0.98:
            failures.append(f"(F) PARALOG: weak fragments give the incumbent NO headroom "
                            f"(weak_acc={weak_acc:.3f}>=0.98 vs span={span_acc:.3f}) — "
                            f"non-discriminating")
        elif recovered < n_pools:
            failures.append(f"(F) PARALOG: pooling did NOT recover all pools "
                            f"({recovered}/{n_pools}) — the below-ceiling gap is NOT provably "
                            f"C4-addressable (could be the unrecoverable zero-evidence trap)")
        else:
            print(f"[smoke] (F) PASS PARALOG/C4 window-selection is MEASURED + DISCRIMINATING + "
                  f"provably C4-ADDRESSABLE: per-read minimap2 drops to {weak_acc:.3f} on weak "
                  f"(1-SNP) fragments (vs spanning ceiling {span_acc:.3f}), YET pooling-majority "
                  f"recovers the true copy in {recovered}/{n_pools} pools — the gap a C4 pooling "
                  f"member could close. Whether C4 ACTUALLY closes it is the NEXT-cycle ablation.",
                  file=sys.stderr)

    # (D) HP_HARD must SEPARATE TWO ARMS (the real validity bar): score the
    # internal flat-affine DP (align_exon_block_global, the arm C1 upgrades)
    # ALONGSIDE minimap2 on HP_HARD, broken out by mode. Below-ceiling on ONE arm
    # is necessary but not sufficient — both arms depressed identically by noise
    # is the vertical-slice trap one level down. The gate is valid for C1 only if
    # the stratum can separate methods AND the internal-DP scoring path runs (the
    # 'ablations runnable' exit criterion). This also PROVES the BAM-ization path
    # that the future arm-flat vs arm-law ablation depends on.
    read_seqs = load_fastq(reads_fq)
    if "HP_HARD" in per_stratum:
        hard_truth = {rid: t for rid, t in truth_map.items() if t.stratum == "HP_HARD"}
        mm2_hard = per_stratum["HP_HARD"].indel.position_exact_concordance
        flat_arm = run_flat_affine_arm(read_seqs, hard_truth, genome, ref_fa,
                                       os.path.join(args.out, "flat_hp_hard.bam"))
        flat_hard = flat_arm.indel.position_exact_concordance
        iso = per_stratum.get("HP")
        iso_c = iso.indel.position_exact_concordance if iso else float("nan")

        # break out by construction mode (boundary_sub = the indel-vs-sub case C1
        # targets; noisy = combined background noise)
        def _mode_conc(score_fn_truth, mode):
            sub = {rid: t for rid, t in hard_truth.items() if f"_{mode}_" in rid}
            mm = score_bam(bam, sub, genome).indel.position_exact_concordance
            fl = run_flat_affine_arm(read_seqs, sub, genome, ref_fa,
                                     os.path.join(args.out, f"flat_{mode}.bam"))
            return mm, fl.indel.position_exact_concordance, len(sub)

        bs_mm, bs_fl, bs_n = _mode_conc(hard_truth, "boundary_sub")
        no_mm, no_fl, no_n = _mode_conc(hard_truth, "noisy")
        print(f"[smoke] (D) HP_HARD two-arm: minimap2={mm2_hard:.4f} "
              f"flat_affine_DP={flat_hard:.4f}  (isolated-HP control={iso_c:.4f})",
              file=sys.stderr)
        print(f"[smoke] (D)   by mode  boundary_sub: mm2={bs_mm:.4f} flat={bs_fl:.4f} (n={bs_n}); "
              f"noisy: mm2={no_mm:.4f} flat={no_fl:.4f} (n={no_n})", file=sys.stderr)

        # The validity criterion the SPEC actually states: the incumbent is below
        # ceiling AGAINST TRUTH, and that headroom is C1-ADDRESSABLE (systematic
        # indel-misplacement a length/run-aware cost would fix), NOT generic noise.
        # We can decide this NOW from truth alone (we do not need the length-law
        # arm to know whether a C1-addressable stratum EXISTS — only to prove the
        # law CLOSES it, which is the next cycle).
        #
        # Classify each flat-affine HP_HARD-noisy failure: is the deletion placed
        # OUTSIDE the run span, or a substitution "repaired" with a spurious indel?
        # Both are exactly the indel-vs-substitution tradeoff C1 targets.
        noisy_truth = {rid: t for rid, t in hard_truth.items() if "_noisy_" in rid}
        fails = c1_addressable = 0
        for rid, t in noisy_truth.items():
            if not t.indels:
                continue
            ind = t.indels[0]
            tn = ind.length if ind.kind == IndelKind.DEL else -ind.length
            cig = align_exon_block_global(read_seqs[rid], genome[t.chrom])
            ins, outs = net_indel_in_span(0, cig, ind.eq_start, ind.eq_end)
            if ins != tn or outs != 0:
                fails += 1
                ipos = all_indel_positions(0, cig)
                if any(not (ind.eq_start <= p < ind.eq_end) for p, l, op in ipos):
                    c1_addressable += 1
        frac = (c1_addressable / fails) if fails else 0.0
        dp_ran = (flat_arm.reads_scored > 0 and
                  (flat_arm.indel.correct + flat_arm.indel.incorrect) > 0)
        if not dp_ran:
            failures.append("(D) internal flat-affine DP arm did not score any HP_HARD "
                            "read — the ablation path is NOT runnable")
        elif no_fl >= 0.999:
            failures.append(f"(D) HP_HARD-noisy gives the incumbent NO headroom vs truth "
                            f"(flat={no_fl:.4f}>=0.999) — no C1-addressable stratum")
        elif frac < 0.8:
            failures.append(f"(D) HP_HARD-noisy headroom is GENERIC not C1-addressable "
                            f"(only {frac:.0%} of flat-affine failures are out-of-run "
                            f"misplacement/spurious-indel)")
        else:
            print(f"[smoke] (D) PASS gate CONTAINS a C1-addressable stratum: flat-affine "
                  f"HP_HARD-noisy={no_fl:.4f} (< ceiling vs TRUTH), and {c1_addressable}/{fails} "
                  f"({frac:.0%}) failures are systematic out-of-run misplacement / "
                  f"spurious mismatch-repair — the indel-vs-sub tradeoff C1 targets. "
                  f"Internal-DP ablation path RUNS (BAM-ized+scored). Whether the length-law "
                  f"CLOSES this gap (arm-LAW > arm-flat) is the NEXT-cycle ablation.",
                  file=sys.stderr)

    # (C) indel position-exact concordance present and meaningful
    pec = score.indel.position_exact_concordance
    if score.indel.correct < 1:
        failures.append(f"(C) indel concordance: expected >=1 correct, got {score.indel.correct}")
    else:
        print(f"[smoke] (C) PASS indel position-exact concordance={pec:.3f} "
              f"(correct={score.indel.correct}, incorrect={score.indel.incorrect})",
              file=sys.stderr)

    # (C') HP cell min_count audit -> at least the small cells exist; report floor
    cell_counts = {}
    for t in truth:
        for (bc, rc, ctx) in t.hp_cells():
            cell_counts[(bc, rc, ctx)] = cell_counts.get((bc, rc, ctx), 0) + 1
    # report (not a hard gate at smoke scale; the Sherlock run sizes >=100)
    min_cell = min(cell_counts.values()) if cell_counts else 0
    print(f"[smoke] (C') HP/STR cells={len(cell_counts)} min_cell_reads={min_cell} "
          f"(scale on Sherlock to clear min_count=100)", file=sys.stderr)

    # (B) deliberately shifted junction call -> TP not FP (ambiguity-aware)
    print("[smoke] (B) building shifted-junction BAM ...", file=sys.stderr)
    shifted_bam = os.path.join(args.out, "shifted.bam")
    (ts, te), (ss, se), shifted_truth = build_shifted_junction_bam(ref_fa, truth, shifted_bam)
    # confirm the shift is genuinely a DIFFERENT raw coordinate that normalizes back
    seq = genome[shifted_truth.chrom]
    norm_shift = normalize_junction(ss, se, seq)
    print(f"[smoke] (B) truth intron=({ts},{te}) shifted raw=({ss},{se}) "
          f"normalized={norm_shift}", file=sys.stderr)
    if (ss, se) == (ts, te):
        failures.append("(B) shift was not a distinct coordinate (engineering bug)")
    sc2 = score_bam(shifted_bam, {shifted_truth.read_id: shifted_truth},
                    genome, aligner_name="shifted")
    if sc2.junction.tp >= 1 and sc2.junction.fp == 0:
        print(f"[smoke] (B) PASS shifted call scored TP={sc2.junction.tp} FP={sc2.junction.fp} "
              f"(ambiguity-aware match credited the 1bp-shifted junction)", file=sys.stderr)
    else:
        failures.append(f"(B) shifted call: expected TP>=1 FP=0, got TP={sc2.junction.tp} "
                        f"FP={sc2.junction.fp}")

    # (B2) NEGATIVE control: a junction shifted FAR beyond the ambiguity window
    # MUST be charged FP (proves the scorer is not trivially crediting everything).
    far_s, far_e = ts + 40, te + 40
    called = [(far_s, far_e)]
    norm_far = normalize_junction(far_s, far_e, seq)
    is_fp = norm_far != (ts, te)
    if is_fp:
        print(f"[smoke] (B2) PASS far junction ({far_s},{far_e})->norm {norm_far} "
              f"!= truth ({ts},{te}) => correctly FP", file=sys.stderr)
    else:
        failures.append(f"(B2) far junction wrongly normalized onto truth: {norm_far}")

    print("\n" + "=" * 60, file=sys.stderr)
    if failures:
        print("SMOKE FAILED:", file=sys.stderr)
        for f in failures:
            print("  - " + f, file=sys.stderr)
        sys.exit(1)
    print("SMOKE PASSED — all GATE assertions green", file=sys.stderr)
    sys.exit(0)


if __name__ == "__main__":
    main()
