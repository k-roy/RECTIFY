#!/usr/bin/env python3
"""Adjudicate the 111 GMAP-only recurrent novel chr5 junctions against the
COMPASS short-read panel (RECTIFY short-read mode, P5 step 4).

Question: do INDEPENDENT short-read aligners (the COMPASS-7 panel: bbmap +
STAR x2 + HISAT2 x2 + magicblast + gsnap, adjudicated to one consensus call per
read) corroborate the 111 junctions that ONLY GMAP placed among the long-read
panel? If 111 INTERSECT COMPASS ~= 0 *and the positive control passes*, the 111
are confirmed GMAP artifacts.

Reuses the live rectify library (same functions + ambiguity normalization as the
long-read corroboration in dev/deliverable_b_gmap_corroboration.py) so the
coordinate convention matches the 111 TSV (which that pipeline produced):
  - collect_junction_counts_from_bam  -> anchor-passing junction Counter
  - normalize_junction                -> leftmost 3-tuple (dissolve ambiguity)
  - load_annotated_junctions          -> annotated catalogue (positive control)
  - register_genome_contigs_from_fasta-> guardrail #6 (trust 'chr5', not chr5->chrV)

Report the THREE numbers together: a near-zero 111-intersection only means
"artifact" if the positive control scored HIGH (else the panel just found nothing).

Run ON Sherlock in the `rectify` env (pysam + the COMPASS contigs registered).

Example:
  W=/scratch/users/kevinroy/compass_a549
  python dev/compass_shortread_adjudicate_111.py \
    --bam   $W/rectify_sr_full/final/A549_rep1.merged.consensus.bam \
    --genome $W/COMPASS/genome_references/GRCh38_gencode_v44.fasta \
    --gtf    $W/COMPASS/genome_references/GRCh38_gencode_v44.gtf \
    --novels /scratch/users/kevinroy/deliverable_b/rectify_src/dev/gmap_only_recurrent_novels_chr5.tsv \
    --out    $W/rectify_sr_full/adjudication_111.json

Author: Kevin R. Roy (+ Claude)
"""
import argparse
import csv
import json
import sys

import pysam  # noqa: F401  (ensures the rectify env)
from rectify.core.splice.junction_scoring import collect_junction_counts_from_bam
from rectify.core.consensus.chimeric_consensus import normalize_junction
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.utils.genome import register_genome_contigs_from_fasta

MIN_ANCHOR = 10        # clean exon overhang on both flanks (DEFAULT_MIN_JUNCTION_ANCHOR)
SUPPORT_DEPTH = 1      # >=1 anchored read = "COMPASS supports this junction"
DECOY_OFFSET = 1_000_000   # negative control: shift each novel 1 Mb downstream


def _load_chrom_seq(fa_path, chrom):
    fa = pysam.FastaFile(fa_path)
    name = chrom if chrom in fa.references else (
        "chr" + chrom if "chr" + chrom in fa.references else chrom.replace("chr", ""))
    seq = fa.fetch(name)
    fa.close()
    return seq


def _load_novels(tsv_path, chrom):
    """The 111: header chrom,start,end,gmap_reads,canonical."""
    out = []
    with open(tsv_path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row['chrom'] != chrom:
                continue
            out.append((int(row['start']), int(row['end']),
                        int(row.get('gmap_reads', 0) or 0)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True, help='merged COMPASS consensus BAM')
    ap.add_argument('--genome', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--novels', required=True, help='gmap_only_recurrent_novels_chr5.tsv')
    ap.add_argument('--chrom', default='chr5')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    # GUARDRAIL: register human contigs so standardize_chrom_name trusts 'chr5'
    # verbatim instead of the empty-registry roman fallback (chr5->chrV), which
    # would silently skip every read. Must precede any junction extraction.
    register_genome_contigs_from_fasta(args.genome)
    seq = _load_chrom_seq(args.genome, args.chrom)
    print(f"[adj] {args.chrom} loaded: {len(seq):,} bp", file=sys.stderr)

    # COMPASS junction set with anchored depth (normalized to leftmost).
    raw = collect_junction_counts_from_bam(
        args.bam, chrom_filter=args.chrom, min_anchor_overhang=MIN_ANCHOR)
    compass = {}
    for (c, s, e), k in raw.items():
        ns, ne = normalize_junction(s, e, seq)
        compass[(c, ns, ne)] = compass.get((c, ns, ne), 0) + k
    print(f"[adj] COMPASS {args.chrom} junctions (anchored): {len(compass):,}", file=sys.stderr)

    # Fail LOUD on a silent-zero (chrom-standardization mismatch / empty BAM):
    # a real chr5 short-read panel has thousands of anchored junctions.
    if len(compass) < 100:
        raise RuntimeError(
            f"only {len(compass)} COMPASS junctions on {args.chrom} — likely a "
            "chrom-standardization mismatch or empty/wrong BAM. Aborting rather "
            "than reporting a fake-null intersection.")

    # Annotated chr5 set (positive control), normalized identically.
    annot = set()
    for t in load_annotated_junctions(args.gtf):
        if t[0] != args.chrom:
            continue
        ns, ne = normalize_junction(t[1], t[2], seq)
        annot.add((args.chrom, ns, ne))
    annot_supported = [j for j in annot if compass.get(j, 0) >= SUPPORT_DEPTH]
    # Of the COMPASS junctions, how many are annotated (sanity: should be high).
    compass_annotated = [j for j in compass if j in annot]

    # The 111, normalized.
    novels_raw = _load_novels(args.novels, args.chrom)
    novels = set()
    novel_meta = {}
    for s, e, gmap_reads in novels_raw:
        ns, ne = normalize_junction(s, e, seq)
        novels.add((args.chrom, ns, ne))
        novel_meta[(args.chrom, ns, ne)] = gmap_reads
    novels_in_compass = {j: compass[j] for j in novels if j in compass}

    # Negative control: decoy junctions = each novel shifted DECOY_OFFSET bp.
    decoys = set()
    for s, e, _ in novels_raw:
        s2, e2 = s + DECOY_OFFSET, e + DECOY_OFFSET
        if e2 < len(seq):
            ns, ne = normalize_junction(s2, e2, seq)
            decoys.add((args.chrom, ns, ne))
    decoys_in_compass = {j: compass[j] for j in decoys if j in compass}

    summary = {
        'chrom': args.chrom,
        'min_anchor': MIN_ANCHOR,
        'support_depth': SUPPORT_DEPTH,
        'compass_junctions_total': len(compass),
        'compass_junctions_annotated': len(compass_annotated),
        'annotated_total': len(annot),
        'POSITIVE_control_annotated_supported_by_compass': len(annot_supported),
        'POSITIVE_control_fraction': round(len(annot_supported) / max(1, len(annot)), 4),
        'NEGATIVE_control_decoys_total': len(decoys),
        'NEGATIVE_control_decoys_in_compass': len(decoys_in_compass),
        'novels_111_total': len(novels),
        'INTERSECTION_111_in_compass': len(novels_in_compass),
    }
    verdict = (
        "ARTIFACTS_CONFIRMED" if (len(novels_in_compass) <= max(1, len(decoys_in_compass))
                                  and len(annot_supported) >= 100)
        else "INCONCLUSIVE_or_CORROBORATED")
    summary['verdict'] = verdict

    detail = {
        'novels_in_compass': {f"{c}:{s}-{e}": {'compass_depth': d,
                                               'gmap_reads': novel_meta.get((c, s, e))}
                              for (c, s, e), d in novels_in_compass.items()},
        'decoys_in_compass': {f"{c}:{s}-{e}": d for (c, s, e), d in decoys_in_compass.items()},
    }
    print(f"[adj] SUMMARY:\n{json.dumps(summary, indent=2)}", file=sys.stderr)
    with open(args.out, 'w') as fh:
        json.dump({'summary': summary, 'detail': detail}, fh, indent=2)
    print(f"[adj] wrote {args.out}", file=sys.stderr)
    print(f"[adj] VERDICT: {verdict} "
          f"(111∩COMPASS={len(novels_in_compass)}, decoys∩COMPASS={len(decoys_in_compass)}, "
          f"pos-control annotated supported={len(annot_supported)}/{len(annot)})", file=sys.stderr)


if __name__ == '__main__':
    main()
