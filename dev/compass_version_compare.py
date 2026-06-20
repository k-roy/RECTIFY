#!/usr/bin/env python3
"""Aligner-version sensitivity check: compare the chr5 junction sets produced by
the COMPASS-pinned panel vs the latest-version panel on the SAME reads.

If the two panels call largely the same junctions at correlated depths (high
Jaccard + high depth correlation), aligner version plays a minimal role and the
111-adjudication is robust to it. Both BAMs must come from the identical rectify
pipeline on the identical input, differing ONLY in aligner binaries/indices.

Run ON Sherlock in the `rectify` env.

Example:
  python dev/compass_version_compare.py \
    --bam-a $W/ver_cmp/compass/cmpver.rectified.bam --label-a compass_pinned \
    --bam-b $W/ver_cmp/latest/latestver.rectified.bam --label-b latest \
    --genome $W/COMPASS/genome_references/GRCh38_gencode_v44.fasta \
    --gtf    $W/COMPASS/genome_references/GRCh38_gencode_v44.gtf \
    --novels /scratch/users/kevinroy/deliverable_b/rectify_src/dev/gmap_only_recurrent_novels_chr5.tsv \
    --out    $W/ver_cmp/version_compare.json
"""
import argparse
import csv
import json
import sys

import pysam  # noqa: F401
from rectify.core.splice.junction_scoring import collect_junction_counts_from_bam
from rectify.core.consensus.chimeric_consensus import normalize_junction
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.utils.genome import register_genome_contigs_from_fasta

MIN_ANCHOR = 10


def _chrom_seq(fa_path, chrom):
    fa = pysam.FastaFile(fa_path)
    name = chrom if chrom in fa.references else "chr" + chrom
    seq = fa.fetch(name)
    fa.close()
    return seq


def _junctions(bam, chrom, seq):
    raw = collect_junction_counts_from_bam(bam, chrom_filter=chrom, min_anchor_overhang=MIN_ANCHOR)
    out = {}
    for (c, s, e), k in raw.items():
        ns, ne = normalize_junction(s, e, seq)
        out[(ns, ne)] = out.get((ns, ne), 0) + k
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam-a', required=True); ap.add_argument('--label-a', default='A')
    ap.add_argument('--bam-b', required=True); ap.add_argument('--label-b', default='B')
    ap.add_argument('--genome', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--novels', required=True)
    ap.add_argument('--chrom', default='chr5')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    register_genome_contigs_from_fasta(args.genome)
    seq = _chrom_seq(args.genome, args.chrom)

    A = _junctions(args.bam_a, args.chrom, seq)
    B = _junctions(args.bam_b, args.chrom, seq)
    print(f"[cmp] {args.label_a}: {len(A):,} junctions; {args.label_b}: {len(B):,}", file=sys.stderr)
    if len(A) < 100 or len(B) < 100:
        raise RuntimeError(f"too few junctions (A={len(A)}, B={len(B)}) — bad BAM/chrom mismatch")

    sa, sb = set(A), set(B)
    shared = sa & sb
    only_a, only_b = sa - sb, sb - sa
    union = sa | sb
    jaccard = len(shared) / len(union) if union else 0.0

    # Depth correlation on shared junctions.
    depth_corr = None
    if shared:
        xa = [A[j] for j in shared]
        xb = [B[j] for j in shared]
        try:
            from scipy.stats import spearmanr
            depth_corr = round(float(spearmanr(xa, xb).correlation), 4)
        except Exception:
            depth_corr = None

    # Annotated vs novel split (novel = not annotated). Version effects matter
    # most for NOVEL junctions (the 111 are novel).
    annot = set()
    for t in load_annotated_junctions(args.gtf):
        if t[0] == args.chrom:
            ns, ne = normalize_junction(t[1], t[2], seq)
            annot.add((ns, ne))
    novel_a = sa - annot
    novel_b = sb - annot
    novel_shared = novel_a & novel_b
    novel_union = novel_a | novel_b
    novel_jaccard = len(novel_shared) / len(novel_union) if novel_union else 0.0

    # The 111 on each panel.
    novels111 = set()
    with open(args.novels) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row['chrom'] != args.chrom:
                continue
            ns, ne = normalize_junction(int(row['start']), int(row['end']), seq)
            novels111.add((ns, ne))
    in_a = len(novels111 & sa)
    in_b = len(novels111 & sb)

    summary = {
        'chrom': args.chrom, 'min_anchor': MIN_ANCHOR,
        f'junctions_{args.label_a}': len(A),
        f'junctions_{args.label_b}': len(B),
        'shared': len(shared), f'only_{args.label_a}': len(only_a),
        f'only_{args.label_b}': len(only_b),
        'jaccard_all': round(jaccard, 4),
        'depth_spearman_on_shared': depth_corr,
        'novel_jaccard': round(novel_jaccard, 4),
        f'novel_{args.label_a}': len(novel_a), f'novel_{args.label_b}': len(novel_b),
        'novel_shared': len(novel_shared),
        'novels_111_total': len(novels111),
        f'111_in_{args.label_a}': in_a, f'111_in_{args.label_b}': in_b,
    }
    # Heuristic readout: high overall + novel Jaccard => version plays minimal role.
    summary['verdict'] = (
        "VERSION_MINIMAL_ROLE" if (jaccard >= 0.9 and novel_jaccard >= 0.75
                                   and in_a == in_b)
        else "VERSION_MATTERS_REVIEW")
    print(f"[cmp] SUMMARY:\n{json.dumps(summary, indent=2)}", file=sys.stderr)
    with open(args.out, 'w') as fh:
        json.dump(summary, fh, indent=2)
    print(f"[cmp] wrote {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
