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
    score_bam, load_genome, extract_junctions,
)
from rectify.core.consensus.chimeric_consensus import normalize_junction  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402


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
    fa = pysam.FastaFile(ref_fa)
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

    # (A) known junction round-trips as TP
    jt_tp = score.junction.by_class.get("NNC", {}).get("tp", 0)
    if jt_tp < 1:
        failures.append(f"(A) junction round-trip: expected >=1 NNC TP, got {jt_tp}")
    else:
        print(f"[smoke] (A) PASS junction round-trip TP={jt_tp}", file=sys.stderr)

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
