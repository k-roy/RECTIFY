#!/usr/bin/env python3
"""P0 benchmark — THIN VERTICAL SLICE (one stratum: homopolymer-length).

Purpose (advisor's critical-path steer): build the smallest end-to-end slice that
EXPOSES the real spec gaps — truth-tag format, the min_count=100 floor actually
biting, and the ambiguity-aware position-exact match wiring — by running BOTH
minimap2 and the current flat-affine DP (align_exon_block_global, the one C1 will
upgrade) through a position-exact scorer on homopolymer-length truth.

Framing metric (per SIMULATION_BENCHMARK_SPEC.md): EXACT INDEL-POSITION
concordance, ambiguity-aware (a homopolymer deletion may sit anywhere in the run),
NOT edit distance (tied by construction).

Model (Tier-1, controlled, truth by construction): a true RNA molecule has an HP
run of length L; the basecaller miscalls it as L-k (del-dominant length error).
Aligning the read (run L-k) to the reference (run L) must yield a k-base DELETION
located WITHIN the run span. An aligner is CORRECT iff net (D-I) within the run ==
k AND it introduces no indel outside the run. k=0 = the clean-run FP control.

Usage:  python hp_vertical_slice.py --out /tmp/hp_slice [--reps 120]
Runs on M1 (light). Needs minimap2 + samtools on PATH and rectify importable.

Author: Kevin R. Roy
"""
import argparse, os, random, subprocess, sys
from collections import Counter, defaultdict

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from rectify.core.align.local_aligner import align_exon_block_global

BASES = "ACGT"
FLANK = 100                      # unique flank each side (confident anchoring)
LENS = list(range(1, 13))        # HP run lengths 1..12
# del-dominant length-error distribution P(k) for k bases UNDER-called
K_DIST = [(0, 0.50), (1, 0.30), (2, 0.14), (3, 0.06)]


def _rand_flank(n, avoid, rng):
    s = [rng.choice([b for b in BASES if b != avoid]) if i in (0, n - 1)
         else rng.choice(BASES) for i in range(n)]
    return "".join(s)


def draw_k(L, rng):
    r, acc = rng.random(), 0.0
    for k, p in K_DIST:
        acc += p
        if r <= acc:
            return min(k, L - 1)        # keep run length >= 1
    return 0


def generate(out, reps, rng):
    refs, reads, truth = {}, [], []
    for b in BASES:
        for L in LENS:
            name = f"hp_{b}_{L:02d}"
            lf, rf = _rand_flank(FLANK, b, rng), _rand_flank(FLANK, b, rng)
            refs[name] = lf + b * L + rf
            run_start, run_end = FLANK, FLANK + L     # ref coords of the run [start,end)
            for i in range(reps):
                k = draw_k(L, rng)
                rid = f"{name}_r{i:03d}_k{k}"
                reads.append((rid, lf + b * (L - k) + rf))
                truth.append((rid, name, run_start, run_end, k, b, L))
    # write reference, reads, truth
    with open(f"{out}/ref.fa", "w") as fh:
        for n, s in refs.items():
            fh.write(f">{n}\n{s}\n")
    with open(f"{out}/reads.fastq", "w") as fh:
        for rid, s in reads:
            fh.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")
    with open(f"{out}/truth.tsv", "w") as fh:
        fh.write("read_id\tref\trun_start\trun_end\ttrue_k\tbase\tL\n")
        for r in truth:
            fh.write("\t".join(map(str, r)) + "\n")
    return refs, {t[0]: t for t in truth}


def net_indel_in_run(cigar, ref_start, run_start, run_end):
    """net (deletion - insertion) bases whose REF position lies within [run_start,
    run_end); also flag any indel OUTSIDE the run (spurious). cigar = [(op,len)]
    with op 0=M/=/X consume both, 1=I consume query, 2=D consume ref."""
    pos = ref_start
    in_run, out_run = 0, 0
    for op, ln in cigar:
        if op in (0, 7, 8):            # M/=/X
            pos += ln
        elif op == 2:                  # D (consumes ref)
            overlap = max(0, min(pos + ln, run_end) - max(pos, run_start))
            in_run += overlap
            out_run += (ln - overlap)
            pos += ln
        elif op == 1:                  # I (consumes query, ref pos fixed)
            if run_start <= pos <= run_end:
                in_run -= ln
            else:
                out_run += ln
        # S/H/N ignored for this stratum
    return in_run, out_run


def score(method, per_read_cigar, truth):
    """per_read_cigar: {read_id: (ref_start, cigar)}. Returns metrics dict."""
    by_L_total, by_L_ok = defaultdict(int), defaultdict(int)
    pos_exact = run_len_ok = total = 0
    for rid, (rs, cig) in per_read_cigar.items():
        if rid not in truth:
            continue
        _, ref, run_start, run_end, k, base, L = truth[rid]
        in_run, out_run = net_indel_in_run(cig, rs, run_start, run_end)
        total += 1
        by_L_total[L] += 1
        size_ok = (in_run == k)
        if size_ok:
            run_len_ok += 1
            by_L_ok[L] += 1
        if size_ok and out_run == 0:
            pos_exact += 1
    return {
        "method": method, "scored": total,
        "run_length_accuracy": round(run_len_ok / total, 4) if total else 0.0,
        "position_exact_concordance": round(pos_exact / total, 4) if total else 0.0,
        "by_L_run_length_accuracy": {L: round(by_L_ok[L] / by_L_total[L], 3)
                                     for L in sorted(by_L_total) if by_L_total[L]},
        "min_cell_count": min(by_L_total.values()) if by_L_total else 0,
    }


def run_minimap2(out):
    bam = f"{out}/mm2.bam"
    p1 = subprocess.run(["minimap2", "-ax", "map-ont", "--eqx", "-t", "2",
                         f"{out}/ref.fa", f"{out}/reads.fastq"],
                        capture_output=True)
    if p1.returncode:
        raise RuntimeError("minimap2 failed: " + p1.stderr.decode()[:500])
    sort = subprocess.run(["samtools", "sort", "-o", bam], input=p1.stdout,
                          capture_output=True)
    if sort.returncode:
        raise RuntimeError("samtools sort failed: " + sort.stderr.decode()[:300])
    subprocess.run(["samtools", "index", bam], check=True)
    res = {}
    with pysam.AlignmentFile(bam, "rb") as b:
        for r in b:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            res[r.query_name] = (r.reference_start, r.cigartuples)
    return res


def run_flat_dp(out, refs):
    """Align each read to its reference with the live flat-affine DP (the one C1
    upgrades). penalty_table is implicit (homo_mismatch=-2 default = production)."""
    res = {}
    reads = {}
    with pysam.FastxFile(f"{out}/reads.fastq") as fq:
        for e in fq:
            reads[e.name] = e.sequence
    # ref name is encoded in the read id (hp_<b>_<LL>_r...)
    for rid, seq in reads.items():
        refname = "_".join(rid.split("_")[:3])
        ref = refs[refname]
        cig = align_exon_block_global(seq, ref)        # [(op,len)], ref_start=0
        res[rid] = (0, cig)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=120)   # >=100 to clear the floor
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = random.Random(args.seed)

    refs, truth = generate(args.out, args.reps, rng)
    print(f"[slice] generated {len(refs)} refs, {len(truth)} reads "
          f"({args.reps}/cell, floor=100)", file=sys.stderr)

    import json
    results = []
    results.append(score("minimap2", run_minimap2(args.out), truth))
    results.append(score("flat_affine_DP", run_flat_dp(args.out, refs), truth))

    print(json.dumps(results, indent=2))
    with open(f"{args.out}/slice_results.json", "w") as fh:
        json.dump({"results": results, "reps": args.reps, "seed": args.seed}, fh, indent=2)


if __name__ == "__main__":
    main()
