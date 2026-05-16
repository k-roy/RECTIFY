#!/usr/bin/env python3
"""rectify correct-cdna — UMI-aware per-molecule consensus for PCR-cDNA.

This module is now a thin orchestration layer over the
:mod:`rectify.core.cdna` subpackage; algorithm code (read extraction, UMI
clustering, walkback, consensus, isoforms, GFF, I/O) lives there. Per-version
algorithm history and PCB114.24 chemistry details are documented in
``docs/algorithms/cdna_correct.md``.

Output files (in --out directory):
  - stage1_consensus.fastq.gz    one consensus record per molecule. Per-cluster
                                  SAM-format tags are appended to each read's
                                  comment line for `rectify align -y` to pass
                                  through to the final aligned BAM.

Downstream: run `rectify align` on the FASTQ, then `rectify cdna-analyze` on
the resulting BAM to produce clusters.tsv, isoforms.tsv, and t1t2_pairs.tsv
(on post-align coordinates).

Tag glossary:
  XU  canonical UMI                          XA  polyA tail length
  XO  orient (fwd / rev)                     XG  primary gene (XS join key)
  XC  cluster size                           XI  isoform_id
  XR  input read IDs                         XT  read type (1=SSP+UMI, 2=SSP-less)
  XM  consensus method                       XB  strand-split (n_top/n_bottom)
  XS  sense/antisense/unannotated            XL  partner cluster_id (T1↔T2)
  XF  full-length tier (0/1/2)              XY  read subtype (umi_captured_fwd/rev, umi_not_captured)
  XQ  5' pre-trim bases stripped             XK  3' pre-trim bases stripped
       (SSP+UMI+GGG for T1 / polyT for rev)       (polyA for fwd / SSP_RC suffix for rev)

Usage (via rectify CLI):
    rectify correct-cdna INPUT.bam --out OUTDIR [options]
"""

from __future__ import annotations

import argparse
import logging
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rectify.core.cdna.cluster import cluster_reads
from rectify.core.cdna.consensus import HAS_POA
from rectify.core.cdna.gff import load_rdna_intervals
from rectify.core.cdna.io import stream_reads, write_stage1_fastq
from rectify.core.cdna.umi import canonical_umi


def run(args) -> int:
    """v3 RECTIFY integration entry point. The argparse Namespace is built by
    rectify/cli.py's `correct-cdna` subparser; we coerce string paths to Path
    objects here (the rectify CLI uses bare `default=None` strings, while the
    algorithm code below expects pathlib.Path)."""
    args.bam = Path(args.bam) if not isinstance(args.bam, Path) else args.bam
    args.out = Path(args.out) if not isinstance(args.out, Path) else args.out
    if getattr(args, "reference", None) is not None:
        args.reference = Path(args.reference) if not isinstance(args.reference, Path) else args.reference
    for attr, default in [
        ("umi_clustering", "directional"),
        ("strand_aware_consensus", False),
        ("no_mask_rdna", False),
    ]:
        if not hasattr(args, attr): setattr(args, attr, default)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger("cdna_correct")

    if not args.bam.exists():
        log.error("BAM not found: %s", args.bam)
        return 2
    args.out.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    log.info("Streaming reads from %s region=%s ...", args.bam, args.region or "all")
    if args.reference is None:
        log.warning("No --reference provided: walk-back anchor canonicalization + "
                    "polyA tail-length measurement DISABLED (legacy aln-end anchor used)")

    # rDNA masking (default ON): skip reads overlapping rRNA_gene intervals to
    # avoid O(n²) UMI clustering on chrXII rDNA tandem repeats.
    rdna_intervals: Optional[Dict[str, List[Tuple[int, int]]]] = None
    if args.gff and not args.no_mask_rdna:
        rdna_intervals = load_rdna_intervals(args.gff)
        n_rdna_loci = sum(len(v) for v in rdna_intervals.values())
        if n_rdna_loci:
            log.info("rDNA masking ON: %d rRNA_gene loci across %d chrom(s) "
                     "(--no-mask-rdna to disable)", n_rdna_loci, len(rdna_intervals))

    reads, n_rdna_masked = stream_reads(args.bam, args.region,
                                        reference=args.reference,
                                        rdna_intervals=rdna_intervals)
    log.info("  %d reads with extractable UMIs (%.1fs)", len(reads), time.time() - t0)
    if n_rdna_masked:
        log.info("  rDNA masked: %d reads skipped (rRNA_gene overlap)", n_rdna_masked)

    t1 = time.time()
    log.info("Stage 1: clustering (anchor window %d, UMI Lev≤%d, method=%s)",
             args.anchor_window_bp, args.umi_edit_distance, args.umi_clustering)
    clusters, stats = cluster_reads(reads, args.anchor_window_bp,
                                     args.umi_edit_distance, args.per_cluster_cap,
                                     clustering_method=args.umi_clustering)
    log.info("  %d molecule clusters (%.1fs)  biggest bucket=%d reads",
             stats["molecule_clusters"], time.time() - t1, stats["biggest_bucket_size"])

    # Canonical UMI per cluster
    umi_canon = {cid: canonical_umi([r.umi for r in c]) for cid, c in enumerate(clusters)}

    # Per-cluster XF tier: max of constituent reads (any one full-length is enough)
    cluster_xf_tier = {cid: max((r.xf_tier for r in c), default=0) for cid, c in enumerate(clusters)}

    # Per-cluster polyA tail length: median across reads (sequence-level, biased
    # short for long tails — see walkback module). Singletons report the lone
    # read's per-read tail length.
    def _median(xs: List[int]) -> int:
        if not xs: return 0
        s = sorted(xs)
        return s[len(s) // 2]
    cluster_tail_len = {cid: _median([r.tail_len for r in c])
                        for cid, c in enumerate(clusters)}

    # Gene assignment (XG/XS), isoform clustering (XI), and Type-1↔Type-2
    # pairing (XL) moved to `rectify cdna-analyze`. They run on post-align
    # coordinates from `rectify align`, which are more accurate than the
    # pre-align positions available here. correct-cdna now emits ONLY the
    # alignment-independent stage-1 outputs: the consensus FASTQ.

    # Stage 1 FASTQ — per-cluster consensus sequences for `rectify align` to
    # remap with the triple aligner. Per-cluster tags ride along on the FASTQ
    # comment line so the downstream BAM picks them up via `minimap2 -y`.
    t2 = time.time()
    out_fastq = args.out / "stage1_consensus.fastq.gz"
    log.info("Writing Stage-1 consensus FASTQ → %s", out_fastq)
    # POA still uses pyabpoa for sequence-level consensus; no aligner needed here.
    use_poa = HAS_POA and not args.no_poa

    if use_poa:
        log.info("Multi-read cluster consensus: POA (pyabpoa)")
    else:
        log.info("Multi-read cluster consensus: pileup-style (substitution-corrective only)")

    fastq_stats = write_stage1_fastq(args.bam, out_fastq, clusters, umi_canon,
                                      cluster_xf_tier, cluster_tail_len,
                                      use_poa=use_poa,
                                      strand_aware_consensus=args.strand_aware_consensus)
    log.info("  wrote %d records (%d singletons, %d pileup consensus, %d rep fallback) in %.1fs",
             fastq_stats["written"], fastq_stats["from_singletons"],
             fastq_stats["from_multi_pileup"], fastq_stats["from_multi_fallback"],
             time.time() - t2)

    # Final report
    print()
    print("=" * 70)
    print("cdna_correct v1 — Stage-1 dedup complete")
    print("=" * 70)
    print(f"input reads (UMI-extractable): {len(reads):>8d}")
    print(f"output records (one per molecule): {fastq_stats['written']:>8d}"
          f"  ({100 * fastq_stats['written'] / max(1, len(reads)):.1f}% of input)")
    print(f"  singletons (passed through):     {fastq_stats['from_singletons']:>8d}")
    print(f"  multi-read (pileup consensus):   {fastq_stats['from_multi_pileup']:>8d}")
    print(f"  multi-read (rep fallback):       {fastq_stats['from_multi_fallback']:>8d}")
    print(f"polyA-pileup buckets dropped:    {stats['buckets_dropped_polyA_pileup']:>8d}"
          f"  ({stats['reads_in_dropped_buckets']} reads — these need POA + position-aware handling)")
    n_t1 = stats.get('type1_reads', 0); n_t2 = stats.get('type2_reads', 0)
    n_total_typed = max(1, n_t1 + n_t2)
    print(f"Read-type breakdown:")
    print(f"  Type 1 (SSP+UMI captured):       {n_t1:>8d}  ({100*n_t1/n_total_typed:.1f}%)  → {stats.get('type1_clusters', 0)} clusters")
    print(f"  Type 2 (SSP-less, 5'-truncated): {n_t2:>8d}  ({100*n_t2/n_total_typed:.1f}%)  → {stats.get('type2_clusters', 0)} clusters")
    print()
    print("XF tier breakdown (full-length confidence):")
    tier_counts: defaultdict[int, int] = defaultdict(int)
    for v in cluster_xf_tier.values():
        tier_counts[v] += 1
    n_clust = max(1, len(clusters))
    print(f"  XF=2 (anchored, HIGH confidence):     {tier_counts[2]:>8d} clusters"
          f"  ({100*tier_counts[2]/n_clust:.1f}%)")
    print(f"  XF=1 (unanchored, MEDIUM confidence): {tier_counts[1]:>8d} clusters"
          f"  ({100*tier_counts[1]/n_clust:.1f}%)")
    print(f"  XF=0 (not detected):                  {tier_counts[0]:>8d} clusters"
          f"  ({100*tier_counts[0]/n_clust:.1f}%)")
    n_full = tier_counts[1] + tier_counts[2]
    print(f"  XF≥1 (any full-length):               {n_full:>8d} clusters"
          f"  ({100*n_full/n_clust:.1f}%)")
    # PolyA tail length distribution (sequence-level; bias-uniform across conditions)
    print()
    print("PolyA tail-length distribution (per-cluster median, sequence-level, XF=2 only):")
    bins = [(0, 1), (1, 10), (10, 20), (20, 30), (30, 50), (50, 75),
            (75, 100), (100, 150), (150, 250), (250, 10_000)]
    bin_counts = [0] * len(bins)
    for cid in range(len(clusters)):
        if cluster_xf_tier[cid] < 2: continue
        tl = cluster_tail_len[cid]
        for i, (lo, hi) in enumerate(bins):
            if lo <= tl < hi:
                bin_counts[i] += 1; break
    n_with_tl = sum(bin_counts)
    print(f"  (restricted to XF=2 anchored clusters: N={n_with_tl})")
    print(f"  {'range':<12} {'count':>8} {'pct':>6}")
    for (lo, hi), n in zip(bins, bin_counts):
        label = f"{lo}-{hi-1}" if hi < 10_000 else f"≥{lo}"
        pct = 100 * n / max(1, n_with_tl)
        print(f"  {label:<12} {n:>8d} {pct:>5.1f}%")
    print()
    print(f"Output FASTQ: {out_fastq}")
    print("Next step: `rectify align` on this FASTQ → `rectify cdna-analyze` "
          "for clusters.tsv / isoforms.tsv / t1t2_pairs.tsv.")
    log.info("Total runtime: %.1fs", time.time() - t0)
    return 0


# Module is invoked via rectify CLI: `rectify correct-cdna ...` calls run(args).
# No __main__ block — the standalone entry point lives in the ont_cdna staging
# repo at /Users/kevinroy/work/ont_cdna/src/cdna_correct.py.


def create_correct_cdna_parser(subparsers):
    """Wire the `correct-cdna` subcommand into the given subparsers group."""
    # =========================================================================
    # correct-cdna command (UMI-aware per-molecule consensus from cDNA BAM)
    # =========================================================================
    correct_cdna_parser = subparsers.add_parser(
        'correct-cdna',
        help='UMI-aware per-molecule consensus for PCR-cDNA BAMs (PCB114.24 chemistry)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            'Cluster aligned PCR-cDNA reads (SQK-PCB114.24 SSP + 27-nt UMI architecture) '
            'by (locus, orientation, UMI) and emit one consensus record per starting RNA '
            'molecule.\n\n'
            'Operates AFTER alignment, complementing trim-cdna-polya which runs on the '
            'FASTQ before alignment. Emits a representative-read or pileup-based '
            'consensus per cluster (POA if pyabpoa is available).\n\n'
            'Output: stage1_consensus.fastq.gz — one consensus sequence per UMI cluster, '
            'with alignment-independent SAM-tag comments (XU/XO/XC/XR/XM/XF/XA/XT/XY/XQ/XK/XB) '
            'for `rectify align -y` to propagate into the post-align BAM. Gene assignment, '
            'isoform clustering, and Type-1↔Type-2 pairing run downstream in '
            '`rectify cdna-analyze`.'
        ),
    )
    correct_cdna_parser.add_argument(
        'bam',
        help='Input BAM (aligned PCR-cDNA reads, indexed)',
    )
    correct_cdna_parser.add_argument(
        '-o', '--out',
        dest='out',
        required=True,
        help='Output directory (will contain stage1_consensus.fastq.gz)',
    )
    correct_cdna_parser.add_argument(
        '--umi-edit-distance',
        dest='umi_edit_distance',
        type=int,
        default=3,
        help='Max Levenshtein distance between UMIs in the same cluster',
    )
    correct_cdna_parser.add_argument(
        '--anchor-window-bp',
        dest='anchor_window_bp',
        type=int,
        default=25,
        help='Window around alignment-end anchor for same-locus clustering (bp)',
    )
    correct_cdna_parser.add_argument(
        '--per-cluster-cap',
        dest='per_cluster_cap',
        type=int,
        default=200,
        help='Max reads per cluster (larger = potential PCR-jackpot signal, slower POA)',
    )
    correct_cdna_parser.add_argument(
        '--gff',
        default=None,
        help='Genome annotation GFF3 (gzip OK). Used for rDNA masking (reads '
             'overlapping rRNA_gene loci are excluded by default to prevent the '
             'O(n²) UMI clustering blow-up on chrXII tandem repeats). Sense/antisense '
             'XS classification, XG gene-name tagging, and isoform clustering now '
             'run downstream in `rectify cdna-analyze`.',
    )
    correct_cdna_parser.add_argument(
        '--reference',
        default=None,
        help='Reference FASTA (gzip OK). Required for walk-back anchor '
             'canonicalization + poly-A tail-length measurement during UMI '
             'extraction; without it, the legacy aln-end anchor is used.',
    )
    correct_cdna_parser.add_argument(
        '--no-poa',
        action='store_true',
        dest='no_poa',
        default=False,
        help='Force pileup-only consensus even if abPOA is available',
    )
    correct_cdna_parser.add_argument(
        '--region',
        default=None,
        help='Restrict to a single BAM region (e.g. chrI) — useful for testing',
    )
    correct_cdna_parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        default=False,
        help='Verbose DEBUG-level logging',
    )
    correct_cdna_parser.add_argument(
        '--umi-clustering',
        dest='umi_clustering',
        choices=('directional', 'components'),
        default='directional',
        help='UMI clustering method (default: directional, umi_tools-style 2× rule)',
    )
    correct_cdna_parser.add_argument(
        '--strand-aware-consensus',
        dest='strand_aware_consensus',
        action='store_true',
        default=False,
        help='v1.18 strand-aware POA: split reads by is_reverse, build a per-strand '
             'sub-consensus, then merge. Cancels strand-specific systematic errors.',
    )
    correct_cdna_parser.add_argument(
        '--no-mask-rdna',
        dest='no_mask_rdna',
        action='store_true',
        default=False,
        help='Disable rDNA masking. By default, reads overlapping rRNA_gene loci in '
             '--gff are excluded before clustering to prevent the O(n²) UMI '
             'bottleneck on chrXII tandem repeats (observed: 261k reads → 4h42m).',
    )
