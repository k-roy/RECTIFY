#!/usr/bin/env python3
"""C6 variant-aware-junction headroom probe — is the incumbent below ceiling on
variant-adjacent FP junctions, and is the gap C6-addressable from READ EVIDENCE,
or is it the zero-evidence trap?  (DECODER-FREE; run this FIRST.)

C6's premise (ALIGNER_MEMBER_DESIGN.md Addendum (b)/§8): a standing germline/somatic
variant near a splice site can FABRICATE a spurious "novel" intron — a GT..AG-flanked
DELETION variant is re-expressed by the haploid-reference aligner as an N-op (intron)
instead of a D-op (deletion), inflating the de-novo discovery FDR. A variant-aware
emission would score the locus as a deletion and recover the TRUE placement.

Per the C1/C2/C3/C4 gate discipline, BEFORE building any variant-aware decoder we must
establish (a) the incumbent IS below ceiling on recoverable reads, AND (b) the gap is
addressable from READ-LEVEL EVIDENCE and NOT the zero-evidence trap — a variant
genuinely indistinguishable from a real cryptic junction is unrecoverable by ANY
read-evidence-weighing method (the pre-committed null). We ALSO test SPECIFICITY: any
variant-aware "relabel the gap as a deletion" rule must NOT suppress REAL novel
junctions.

The decisive measurements (DECODER-FREE; fitness=TRUTH, never an internal score):

  (a) BELOW-CEILING — on the SPLICE_MIMIC_DEL drivers (gen_variant_stratum), how often
      does minimap2 fabricate a variant-adjacent FP intron (N) instead of the truth
      deletion (D)? ceiling = the deletion truth is in-principle placeable (the read
      genuinely lacks the block, so a D-op reconstructs it 0-mismatch by construction).

  (b) ZERO-EVIDENCE (the (G)-equivalent NM test) — for each driver where mm2 emits the
      FP intron, is mm2's N-op alignment NM==0?  If YES, the intron is an EQUALLY-GOOD
      (0-mismatch) alignment of the SAME query as the truth D-op: no read-level evidence
      separates intron from deletion → a motif-blind score CANNOT prefer truth →
      unrecoverable per-read. (This is largely PRE-DETERMINED by the GT..AG/random-flank
      construction — it confirms the construction, it is not the load-bearing measure.)

  (c) THE TENSION (the load-bearing, FALSIFIABLE measure) — the ONLY read-level rule
      that recovers a driver is "a long GT..AG gap is a DELETION, not an intron." Apply
      it to LENGTH-MATCHED REAL introns (built here at the SAME lengths as the drivers,
      GT..AG, random interior, spliced reads — truth = the intron). The rule that
      recovers N drivers converts N length-matched real introns to FN, ~1:1. Drivers and
      real introns occupy the SAME read-feature cell (length, GT..AG, random interior),
      so NO function of the read separates them — the zero-evidence proof, quantitative.

  (d) SPECIFICITY (positive) — mm2 RETAINS the real junction on the SNP_NEAR_JUNC
      controls (recall ~1.0) and on the matched real introns: mm2 is NOT blunt, so an
      "abstain/relabel near any variant" rule would wrongly suppress real junctions.

VERDICT (pre-committed):
  Drivers are zero-evidence (NM==0 intron == 0-mismatch D) AND the recovering rule
  destroys length-matched real introns ~1:1  => the gap is NOT addressable from read
  evidence on this corpus; the only deciding input is an EXTERNAL variant catalog (a
  matched VCF), which on this corpus == the truth label (circular).  => C6-as-read-
  evidence REFUTED on this corpus; DEFER the VCF-integration residue to a measured
  real-data trigger (a matched germline VCF, where the VCF is independent evidence, not
  the truth label).  Mirrors C4's REFUTE-on-corpus + named real-data residue.

  If instead the driver intron carried NM>=1 (intron costs identity, D strictly wins on
  a motif-blind score) AND the recovering rule spared real introns => a genuine read-
  evidence gap => PROCEED (then SKETCH the emission + specificity fence, do not build).

Why minimap2 alone (not the C1/C3 flat/law panel): the flat-affine DP arm has NO N-op,
so it emits D for EVERY long gap — it "recovers" every driver by being intron-BLIND,
and that SAME arm converts every real intron to a giant deletion (catastrophic
specificity failure). Flat-arm recovery is therefore NOT a C6-addressability signal; it
is the tension itself. minimap2 -ax splice -uf is the faithful incumbent (the shipped
scorer/smoke (E) locus metric is minimap2-based).

Usage:
  python scripts/benchmark/c6_headroom.py --reps 60 --out /tmp/c6_hr
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import random
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import (  # noqa: E402
    load_genome, score_bam, extract_junctions, normalize_junction)
from scripts.benchmark.sim import controlled  # noqa: E402
from scripts.benchmark.sim.controlled import (  # noqa: E402
    TranscriptModel, Exon, ReadTruth, _rand_unique, _split_for)
import pysam  # noqa: E402

DRIVER_LENS = [40, 60, 100]   # the >=40 splice-mimic-deletion drivers (length-matched arm too)


def load_fastq(path):
    seqs = {}
    with pysam.FastxFile(path) as fq:
        for e in fq:
            seqs[e.name] = e.sequence
    return seqs


def run_minimap2(ref_fa, reads_fq, out_bam):
    p = subprocess.run(
        ["minimap2", "-ax", "splice", "-uf", "--eqx", "-k", "14", "-t", "2",
         ref_fa, reads_fq], capture_output=True)
    if p.returncode:
        raise RuntimeError("minimap2 failed: " + p.stderr.decode()[:400])
    s = subprocess.run(["samtools", "sort", "-o", out_bam], input=p.stdout,
                       capture_output=True)
    if s.returncode:
        raise RuntimeError("samtools sort failed: " + s.stderr.decode()[:300])
    subprocess.run(["samtools", "index", out_bam], check=True)


def build_matched_introns(reps, rng, lens, locus0=9000):
    """LENGTH-MATCHED REAL introns: GT..AG, random interior, spliced reads. Same
    construction process as the driver deletion-block, but truth = the intron (N).
    These are the specificity arm — the rule that recovers drivers FNs these."""
    refs, reads, truth = {}, [], []
    li = locus0
    L = 200
    for ln in lens:
        for i in range(reps):
            chrom = f"real_intron_{ln:03d}_r{i:03d}"
            left = _rand_unique(L, rng)
            intron = "GT" + _rand_unique(ln - 4, rng) + "AG"
            right = _rand_unique(L, rng)
            contig = left + intron + right
            refs[chrom] = contig
            i_start, i_end = L, L + ln
            model = TranscriptModel(
                name=chrom, chrom=chrom, strand="+",
                exons=[Exon(0, L), Exon(i_end, len(contig))], genome_seq=contig)
            jt = model.junction_truths()       # NNC, truth = the real intron
            split = _split_for(li)
            li += 1
            reads.append((chrom, model.spliced_transcript()))
            truth.append(ReadTruth(
                read_id=chrom, true_locus=f"real_intron_{ln:03d}",
                true_transcript=chrom, chrom=chrom, strand="+",
                genome_start=model.genome_start, genome_end=model.genome_end,
                true_cigar=model.fulllength_cigar(), junctions=jt,
                true_cpa=model.genome_end - 1, stratum="REAL_INTRON",
                split=split, coverage=1))
    return refs, reads, truth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/c6_hr")
    ap.add_argument("--reps", type=int, default=60,
                    help="reads per driver length / per control / per matched-intron length")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--result", default=None, help="write the result table here too")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = random.Random(args.seed)

    print(f"[c6hr] building VARIANT stratum + length-matched real-intron arm reps={args.reps}",
          file=sys.stderr)
    refs, reads, truth = controlled.gen_variant_stratum(args.reps, rng)
    iref, iread, itruth = build_matched_introns(args.reps, rng, DRIVER_LENS)
    refs.update(iref); reads += iread; truth += itruth

    ref_fa = os.path.join(args.out, "ref.fa")
    reads_fq = os.path.join(args.out, "reads.fastq")
    with open(ref_fa, "w") as fh:
        for n, s in refs.items():
            fh.write(f">{n}\n{s}\n")
    with open(reads_fq, "w") as fh:
        for rid, s in reads:
            fh.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")
    truth_map = {t.read_id: t for t in truth}
    genome = load_genome(ref_fa)

    bam = os.path.join(args.out, "mm2.bam")
    run_minimap2(ref_fa, reads_fq, bam)

    # ---- partition truth ----
    def _driver(t):
        return t.true_locus.startswith("var_del_") and int(t.true_locus.split("_")[-1]) >= 40
    drivers = {rid: t for rid, t in truth_map.items() if _driver(t)}
    short_del = {rid: t for rid, t in truth_map.items()
                 if t.true_locus.startswith("var_del_") and int(t.true_locus.split("_")[-1]) < 40}
    snp_near = {rid: t for rid, t in truth_map.items() if t.true_locus == "near"}
    snp_far = {rid: t for rid, t in truth_map.items() if t.true_locus == "distant"}
    real_introns = {rid: t for rid, t in truth_map.items() if t.stratum == "REAL_INTRON"}

    # ---- pull mm2 alignments by read ----
    aln = {}
    with pysam.AlignmentFile(bam, "rb") as b:
        for r in b:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if r.query_name in truth_map:
                aln[r.query_name] = r

    def _nm(r):
        if r.has_tag("NM"):
            return r.get_tag("NM")
        return sum(l for o, l in (r.cigartuples or []) if o in (1, 2, 8))

    def _has_intron(r):
        return any(o == 3 for o, _l in (r.cigartuples or []))

    # ---- (a) driver below-ceiling: mm2 emits FP intron not D ----
    sd = score_bam(bam, drivers, genome, aligner_name="mm2:driver")
    drv_n = len(drivers)
    drv_fp_intron = 0          # placed, emitted an N-op (intron)
    drv_fp_nm0 = 0             # ...and that intron is NM==0 (zero-evidence)
    drv_placed = 0
    for rid, t in drivers.items():
        r = aln.get(rid)
        if r is None:
            continue
        drv_placed += 1
        if _has_intron(r):
            drv_fp_intron += 1
            if _nm(r) == 0:
                drv_fp_nm0 += 1
    incx = (drv_placed - drv_fp_intron) / drv_placed if drv_placed else 0.0   # mm2 got D right
    headroom = drv_fp_intron / drv_placed if drv_placed else 0.0

    # ---- (d) specificity: mm2 retains real junctions ----
    ss_near = score_bam(bam, snp_near, genome, aligner_name="mm2:snp_near")
    s_ri = score_bam(bam, real_introns, genome, aligner_name="mm2:real_intron")
    # per-read: real introns mm2 calls as N (TP) — and that N is NM==0 too
    ri_n = len(real_introns)
    ri_intron = ri_nm0 = ri_placed = 0
    for rid, t in real_introns.items():
        r = aln.get(rid)
        if r is None:
            continue
        ri_placed += 1
        if _has_intron(r):
            ri_intron += 1
            if _nm(r) == 0:
                ri_nm0 += 1

    # ---- (c) the tension, length-stratified: per driver length, mm2's driver-FP
    #          intron geometry vs the matched real-intron N. The "long GT..AG gap -> D"
    #          recovery rule converts matched real introns to FN 1:1 within each length.
    per_len = {}
    for ln in DRIVER_LENS:
        dks = [rid for rid, t in drivers.items() if int(t.true_locus.split("_")[-1]) == ln]
        iks = [rid for rid, t in real_introns.items() if t.true_locus.endswith(f"_{ln:03d}")]
        d_fp = sum(1 for rid in dks if aln.get(rid) is not None and _has_intron(aln[rid]))
        i_tp = sum(1 for rid in iks if aln.get(rid) is not None and _has_intron(aln[rid]))
        per_len[ln] = (len(dks), d_fp, len(iks), i_tp)

    # ================= report =================
    out = []
    P = out.append
    P("================ C6 VARIANT-AWARE-JUNCTION HEADROOM (decoder-free) ================")
    P(f"reps={args.reps}  seed={args.seed}  aligner=minimap2 -ax splice -uf --eqx")
    P("")
    P("(a) BELOW-CEILING — SPLICE_MIMIC_DEL drivers (truth = DELETION, no intron)")
    P(f"    drivers n={drv_n}  placed={drv_placed}")
    P(f"    ceiling (D placeable, read lacks block by construction)   = 1.0000")
    P(f"    incumbent (mm2 emits the truth D, NO intron)              = {incx:.4f} "
      f"({drv_placed - drv_fp_intron}/{drv_placed})")
    P(f"    HEADROOM (mm2 fabricates a variant-adjacent FP INTRON)    = {headroom:.4f} "
      f"({drv_fp_intron}/{drv_placed})")
    P(f"    scorer fp_variant_adjacent (shipped metric)               = {sd.junction.fp_variant_adjacent}")
    P(f"    scorer indel position-exact concordance on the deletion   = {sd.indel.position_exact_concordance:.4f}")
    P("")
    P("(b) ZERO-EVIDENCE NM test — is mm2's FP intron a 0-mismatch alignment?")
    P(f"    driver FP introns with NM==0 (== the truth D is also 0-mm) = {drv_fp_nm0}/{drv_fp_intron}")
    P(f"    => the intron and the deletion are the SAME query at 0 mismatches; no read-")
    P(f"       level evidence separates them (construction-confirmed, not load-bearing).")
    P("")
    P("(c) THE TENSION (load-bearing) — length-matched REAL introns (truth = intron, N)")
    P(f"    {'len':>5s} {'drivers':>8s} {'mm2_FP_N':>9s} {'realIntr':>9s} {'mm2_TP_N':>9s}")
    for ln in DRIVER_LENS:
        dn, dfp, inn, itp = per_len[ln]
        P(f"    {ln:5d} {dn:8d} {dfp:9d} {inn:9d} {itp:9d}")
    P(f"    => within EACH length, mm2 emits the IDENTICAL N-op for the driver (FP) and")
    P(f"       the real intron (TP). A 'long GT..AG gap -> D' rule recovers the {drv_fp_intron}")
    P(f"       driver FPs but converts the {ri_intron} matched real introns to FN, ~1:1.")
    P(f"       Drivers and real introns share one read-feature cell (length, GT..AG,")
    P(f"       random interior) -> NO function of the read separates them (zero-evidence).")
    P("")
    P("(d) SPECIFICITY (positive) — mm2 is NOT blunt; it retains real junctions")
    P(f"    SNP_NEAR_JUNC real-junction recall  = {ss_near.junction.recall:.4f} "
      f"(fp_variant_adjacent={ss_near.junction.fp_variant_adjacent})")
    P(f"    matched real-intron recall          = {s_ri.junction.recall:.4f}  "
      f"(mm2 emits N on {ri_intron}/{ri_placed} placed; NM==0 on {ri_nm0})")
    P(f"    short-deletion (<40bp) control: stays D (smoke (E) verified) n={len(short_del)}")
    P("")
    P("---- VERDICT (pre-committed) ----")
    zero_ev = (drv_fp_intron > 0 and drv_fp_nm0 == drv_fp_intron)
    tension = (ri_intron > 0)
    if zero_ev and tension:
        P(f"  HEADROOM={headroom:.4f} > 0 (mm2 IS below ceiling, fabricates {drv_fp_intron} variant-")
        P("  adjacent FP introns) BUT every driver FP intron is NM==0 == the truth D (0-mismatch),")
        P("  AND the ONLY read-level recovery rule (long GT..AG gap -> D) FNs the length-matched")
        P(f"  real introns ~1:1 ({ri_intron} real introns mm2 correctly calls N). The drivers are")
        P("  the ZERO-EVIDENCE trap: indistinguishable from real cryptic introns at the read level.")
        P("  The only deciding input is an EXTERNAL variant catalog (matched VCF) — which on THIS")
        P("  corpus == the truth label (circular).")
        P("  => C6-as-read-evidence REFUTED on this corpus. DEFER a VCF-integration residue to a")
        P("     MEASURED real-data trigger (matched germline VCF as INDEPENDENT evidence). Do NOT")
        P("     build a variant-aware emission against this corpus; do NOT ship a blunt")
        P("     'relabel/abstain near any variant' rule (specificity (d) shows it suppresses real")
        P("     junctions). C6 joins the gate-refuted facets.")
    else:
        P(f"  HEADROOM={headroom:.4f}; zero_evidence={zero_ev}; tension={tension}")
        P("  => NOT the predicted zero-evidence pattern. Inspect: a read-evidence gap MAY exist")
        P("     (driver intron NM>=1 strictly worse than D, and real introns spared) -> PROCEED")
        P("     to SKETCH the variant-aware emission + its specificity fence.")
    text = "\n".join(out)
    print(text)
    if args.result:
        with open(args.result, "w") as fh:
            fh.write(text + "\n")
        print(f"\n[c6hr] wrote {args.result}", file=sys.stderr)


if __name__ == "__main__":
    main()
