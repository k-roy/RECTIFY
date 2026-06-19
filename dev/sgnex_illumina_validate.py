#!/usr/bin/env python3
"""Orthogonal validation of A549 novel junctions against SG-NEx Illumina short reads.

Short-read split-read evidence is the gold standard for confirming a spliced junction
(annotation-independent). We extract junctions from the 3 SG-NEx A549 Illumina genome
BAMs (STAR-aligned, so spliced reads carry N-ops) with the SAME ambiguity-aware extractor
used to derive the candidates, then ask how many of each candidate set are corroborated.

CALIBRATED: three sets scored together —
  - the_111   : GMAP-only recurrent GT-AG novels (the target; expect SOME validation)
  - corrob_609: novel-canonical junctions already multi-ONT-aligner corroborated (POSITIVE
                control; expect HIGH validation)
  - noncanon  : a random sample of GMAP's non-canonical junctions (NEGATIVE control;
                expect ~0 — this is the noise the fences suppress, measured vs real truth)

A junction is "Illumina-validated" if >= MIN_READS anchored split reads support it in
>= MIN_REPS of the 3 replicates (after leftmost ambiguity-normalization).

Chrom-naming: SG-NEx Illumina is Ensembl-named ('5'); the candidates are 'chr5'. We register
each BAM's own contig names (verbatim-trust) and compare on the bare coordinate after
normalization (single chromosome, so chrom label is dropped for the match).

Run ON Sherlock. Author: Kevin R. Roy.
"""
import argparse, json, random, sys
from collections import defaultdict

import pysam
from rectify.core.splice.junction_scoring import collect_junction_counts_from_bam
from rectify.core.consensus.chimeric_consensus import normalize_junction
from rectify.utils.genome import register_genome_contigs

MIN_READS = 2      # anchored split reads per replicate
MIN_REPS = 2       # replicates needed
MIN_ANCHOR = 6     # short-read overhang (looser than the 10bp long-read default)
NONCANON_SAMPLE = 3000


def chr5_name(refs):
    for n in ("chr5", "5"):
        if n in refs:
            return n
    raise RuntimeError(f"no chr5/5 contig in BAM refs: {list(refs)[:5]}...")


def norm_set(pairs, seq):
    """leftmost-normalize a list of (start,end) -> set."""
    return {normalize_junction(s, e, seq) for (s, e) in pairs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates-json", required=True)   # gmap_corroboration.json
    ap.add_argument("--genome", required=True)            # for chr5 seq (normalization)
    ap.add_argument("--illumina-bams", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fa = pysam.FastaFile(args.genome)
    chrname = "chr5" if "chr5" in fa.references else "5"
    seq = fa.fetch(chrname)
    fa.close()

    d = json.load(open(args.candidates_json))["buckets"]
    rng = random.Random(13)
    sets = {
        "the_111": [tuple(r["junction"][1:3]) for r in d["gmap_only_recurrent"]],
        "corrob_609": [tuple(r["junction"][1:3]) for r in d["independently_corroborated"]],
        "noncanon": [tuple(r["junction"][1:3]) for r in
                     rng.sample(d["gmap_noncanonical"],
                                min(NONCANON_SAMPLE, len(d["gmap_noncanonical"])))],
    }
    norm = {k: norm_set(v, seq) for k, v in sets.items()}

    # per-replicate Illumina anchored junction support, normalized
    rep_support = []
    for bam in args.illumina_bams:
        with pysam.AlignmentFile(bam, "rb") as b:
            register_genome_contigs(b.references)
            cn = chr5_name(b.references)
        counts = collect_junction_counts_from_bam(bam, chrom_filter=cn,
                                                  min_anchor_overhang=MIN_ANCHOR)
        norm_counts = defaultdict(int)
        for (c, s, e), k in counts.items():
            ns, ne = normalize_junction(s, e, seq)
            norm_counts[(ns, ne)] += k
        rep_support.append(norm_counts)
        print(f"[val] {bam.split('/')[-1]}: {len(norm_counts):,} normalized chr5 junctions",
              file=sys.stderr)

    def validated(j):
        reps_ok = sum(1 for rc in rep_support if rc.get(j, 0) >= MIN_READS)
        total = sum(rc.get(j, 0) for rc in rep_support)
        return reps_ok >= MIN_REPS, reps_ok, total

    report = {}
    for name, jset in norm.items():
        n_val = 0
        detail = []
        for j in jset:
            ok, reps_ok, total = validated(j)
            if ok:
                n_val += 1
            detail.append({"junction": list(j), "validated": ok,
                           "reps_with_support": reps_ok, "total_illumina_reads": total})
        rate = round(n_val / len(jset), 4) if jset else 0.0
        report[name] = {"n_candidates": len(jset), "n_validated": n_val,
                        "validation_rate": rate, "detail": detail}
        print(f"[val] {name}: {n_val}/{len(jset)} validated ({rate:.1%})", file=sys.stderr)

    with open(args.out, "w") as fh:
        json.dump({"params": {"min_reads": MIN_READS, "min_reps": MIN_REPS,
                              "min_anchor": MIN_ANCHOR, "noncanon_sample": NONCANON_SAMPLE},
                   "report": report}, fh)
    print(f"[val] wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
