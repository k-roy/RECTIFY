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
import gzip
import logging
import subprocess
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rectify.core.cdna.cluster import cluster_reads
from rectify.core.cdna.consensus import HAS_POA
from rectify.core.cdna.gff import load_rdna_intervals
from rectify.core.cdna.io import stream_reads, write_stage1_fastq
from rectify.core.cdna.umi import canonical_umi


log = logging.getLogger("cdna_correct")


# ---------------------------------------------------------------------------
# Region-parallel helpers
# ---------------------------------------------------------------------------


def _cdna_region_task(
    bam_path: str,
    region_str: str,
    region_name: str,
    out_fastq_path: str,
    reference_path: Optional[str],
    rdna_intervals: Optional[Dict],
    anchor_window_bp: int,
    umi_edit_distance: int,
    per_cluster_cap: int,
    umi_clustering: str,
    use_poa: bool,
    strand_aware_consensus: bool,
) -> Dict:
    """Process one BAM region into a FASTQ chunk.

    Workers open their own fresh pysam handle (no parent handle inherited).
    Cluster IDs are prefixed with region_name for global uniqueness.
    Returns stats dict.
    """
    from pathlib import Path as _P
    from rectify.core.cdna.cluster import cluster_reads as _cluster
    from rectify.core.cdna.io import stream_reads as _stream, write_stage1_fastq as _write
    from rectify.core.cdna.umi import canonical_umi as _canon_umi
    from rectify.core.cdna.consensus import HAS_POA as _HAS_POA

    reference = _P(reference_path) if reference_path else None
    reads, n_rdna_masked = _stream(
        _P(bam_path), region_str, reference=reference, rdna_intervals=rdna_intervals
    )
    if not reads:
        # Write empty gzip FASTQ so the merge step has a file to cat
        import gzip as _gz
        with _gz.open(out_fastq_path, "wt"):
            pass
        return {"input_reads": 0, "written": 0, "from_singletons": 0,
                "from_multi_pileup": 0, "from_multi_fallback": 0, "n_rdna_masked": 0}

    clusters, stats = _cluster(reads, anchor_window_bp, umi_edit_distance, per_cluster_cap,
                                clustering_method=umi_clustering)

    def _median(xs: List[int]) -> int:
        if not xs: return 0
        s = sorted(xs)
        return s[len(s) // 2]

    umi_canon = {cid: _canon_umi([r.umi for r in c]) for cid, c in enumerate(clusters)}
    cluster_xf_tier = {cid: max((r.xf_tier for r in c), default=0)
                       for cid, c in enumerate(clusters)}
    cluster_tail_len = {cid: _median([r.tail_len for r in c])
                        for cid, c in enumerate(clusters)}

    fq_stats = _write(
        _P(bam_path), _P(out_fastq_path), clusters, umi_canon,
        cluster_xf_tier, cluster_tail_len,
        use_poa=use_poa and _HAS_POA,
        strand_aware_consensus=strand_aware_consensus,
        cluster_name_prefix=region_name,
    )
    fq_stats["n_rdna_masked"] = n_rdna_masked
    return fq_stats


def _run_cdna_correct_parallel(
    bam_path: Path,
    out_fastq: Path,
    reference: Optional[Path],
    rdna_intervals: Optional[Dict],
    anchor_window_bp: int,
    umi_edit_distance: int,
    per_cluster_cap: int,
    umi_clustering: str,
    use_poa: bool,
    strand_aware_consensus: bool,
    workers: int,
) -> Dict:
    """Region-parallel stage-1 consensus FASTQ.

    Axis-2 safety: plan_regions opens+closes pysam before the pool is created;
    workers each open a fresh pysam handle.
    """
    import tempfile
    from rectify.core.bam.regions import plan_regions

    bam_path_str = str(bam_path)
    bai = bam_path_str + ".bai"
    alt_bai = bam_path_str.replace(".bam", ".bai")
    if not Path(bai).exists() and not Path(alt_bai).exists():
        log.info("No BAM index — indexing %s", bam_path_str)
        subprocess.run(["samtools", "index", bam_path_str], check=True)

    with tempfile.TemporaryDirectory(prefix="cdna_correct_par_") as tmp_str:
        tmp = Path(tmp_str)
        plans = plan_regions(bam_path_str, tmp)

        futures: Dict = {}
        per_region_fastqs: List[str] = []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for plan in plans:
                region_str = f"{plan.chrom}:{plan.start + 1}-{plan.end}"
                region_name = plan.region_id
                out_fq = str(tmp / f"{plan.region_id}.fastq.gz")
                per_region_fastqs.append(out_fq)
                fut = pool.submit(
                    _cdna_region_task,
                    bam_path_str, region_str, region_name,
                    out_fq,
                    str(reference) if reference else None,
                    rdna_intervals,
                    anchor_window_bp, umi_edit_distance, per_cluster_cap,
                    umi_clustering, use_poa, strand_aware_consensus,
                )
                futures[fut] = plan.region_id

            total_stats: Dict = {"input_reads": 0, "written": 0, "from_singletons": 0,
                                  "from_multi_pileup": 0, "from_multi_fallback": 0}
            for fut in as_completed(futures):
                s = fut.result()
                for k in total_stats:
                    total_stats[k] += s.get(k, 0)

        # Concatenate gzipped FASTQ shards into the final output
        existing_fastqs = [f for f in per_region_fastqs if Path(f).exists()]
        if existing_fastqs:
            with gzip.open(str(out_fastq), "wt") as fout:
                for fq_path in existing_fastqs:
                    with gzip.open(fq_path, "rt") as fin:
                        for line in fin:
                            fout.write(line)
        else:
            with gzip.open(str(out_fastq), "wt"):
                pass  # empty output

    return total_stats


def run(args) -> int:
    """v3 RECTIFY integration entry point. The argparse Namespace is built by
    rectify/cli.py's `correct-cdna` subparser; we coerce string paths to Path
    objects here (the rectify CLI uses bare `default=None` strings, while the
    algorithm code below expects pathlib.Path)."""
    import sys as _sys
    from datetime import datetime as _dt, timezone as _tz
    from time import perf_counter as _perf

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

    if not args.bam.exists():
        log.error("BAM not found: %s", args.bam)
        return 2
    args.out.mkdir(parents=True, exist_ok=True)

    # --- Resume skip-check ---
    _current_inputs = {"input_bam": str(args.bam)}
    from rectify.core.provenance import can_skip_stage
    from rectify.utils.version import get_rectify_git_sha as _get_sha
    _cdna_sha = _get_sha()
    _sample_id = args.out.name
    _skip = can_skip_stage(
        stage="cdna_stage1",
        sample_output=args.out,
        sample_id=_sample_id,
        current_inputs=_current_inputs,
        current_argv=_sys.argv,
        rectify_git_sha=_cdna_sha,
        force=getattr(args, "force_all", False),
        force_stages=set(getattr(args, "force_stage", "").split(","))
                     if getattr(args, "force_stage", None) else None,
        accept_prior_provenance=getattr(args, "accept_prior_provenance", False),
    )
    log.info("[RESUME] sample=%s stage=cdna_stage1 decision=%s reason=%s",
             _sample_id, "SKIP" if _skip.skip else "RUN", _skip.reason)
    if _skip.skip:
        return 0
    if getattr(args, "dry_run_resume", False):
        print(f"[dry-run-resume] stage=cdna_stage1 decision=RUN reason={_skip.reason}")
        return 0

    _stage_started_at = _dt.now(_tz.utc).isoformat()
    _t_start = _perf()

    t0 = time.time()
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

    use_poa = HAS_POA and not args.no_poa
    out_fastq = args.out / "stage1_consensus.fastq.gz"
    _workers = max(1, getattr(args, "workers", 1) or 1)

    if _workers > 1 and not args.region:
        log.info("Parallel stage-1 (%d workers) → %s", _workers, out_fastq)
        fastq_stats = _run_cdna_correct_parallel(
            args.bam, out_fastq, args.reference, rdna_intervals,
            args.anchor_window_bp, args.umi_edit_distance, args.per_cluster_cap,
            args.umi_clustering, use_poa, args.strand_aware_consensus, _workers,
        )
        print()
        print("=" * 70)
        print(f"cdna_correct v1 — Stage-1 dedup complete (parallel {_workers} workers)")
        print("=" * 70)
        print(f"output records (one per molecule): {fastq_stats['written']:>8d}")
        print(f"  singletons (passed through):     {fastq_stats['from_singletons']:>8d}")
        print(f"  multi-read (pileup/POA):         {fastq_stats['from_multi_pileup']:>8d}")
        print(f"  multi-read (rep fallback):       {fastq_stats['from_multi_fallback']:>8d}")
        print()
        print(f"Output FASTQ: {out_fastq}")
        print("Next step: `rectify align` on this FASTQ → `rectify cdna-analyze`")
    else:
        log.info("Streaming reads from %s region=%s ...", args.bam, args.region or "all")
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

        umi_canon = {cid: canonical_umi([r.umi for r in c]) for cid, c in enumerate(clusters)}
        cluster_xf_tier = {cid: max((r.xf_tier for r in c), default=0)
                           for cid, c in enumerate(clusters)}

        def _median(xs: List[int]) -> int:
            if not xs: return 0
            s = sorted(xs)
            return s[len(s) // 2]
        cluster_tail_len = {cid: _median([r.tail_len for r in c])
                            for cid, c in enumerate(clusters)}

        t2 = time.time()
        log.info("Writing Stage-1 consensus FASTQ → %s", out_fastq)
        if use_poa:
            log.info("Multi-read cluster consensus: POA (pyabpoa)")
        else:
            log.info("Multi-read cluster consensus: pileup-style (substitution-corrective only)")

        fastq_stats = write_stage1_fastq(args.bam, out_fastq, clusters, umi_canon,
                                          cluster_xf_tier, cluster_tail_len,
                                          use_poa=use_poa,
                                          strand_aware_consensus=args.strand_aware_consensus)
        log.info("  wrote %d records (%d singletons, %d pileup, %d rep fallback) in %.1fs",
                 fastq_stats["written"], fastq_stats["from_singletons"],
                 fastq_stats["from_multi_pileup"], fastq_stats["from_multi_fallback"],
                 time.time() - t2)

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

    # Write stage sidecar (non-fatal on failure)
    try:
        from rectify.core.provenance import ProvenanceRecord, write_stage_sidecar
        _sc_record = ProvenanceRecord.from_components(
            stage="cdna_stage1",
            sample_id=_sample_id,
            sample_output_dir=args.out,
            started_at=_stage_started_at,
            completed_at=_dt.now(_tz.utc).isoformat(),
            exit_status=0,
            inputs=_current_inputs,
            outputs={"stage1_fastq": str(out_fastq)},
            stats={
                "wall_seconds": _perf() - _t_start,
                "written": fastq_stats.get("written", 0),
                "from_singletons": fastq_stats.get("from_singletons", 0),
            },
            argv=_sys.argv,
            rectify_git_sha=_cdna_sha,
        )
        write_stage_sidecar(_sc_record, sample_output=args.out)
    except Exception as _sc_exc:
        log.warning("Failed to write cdna_stage1 sidecar: %s", _sc_exc)

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
    correct_cdna_parser.add_argument(
        '--workers',
        type=int,
        default=1,
        help='Number of parallel worker processes for stage-1 consensus (default: 1 = sequential). '
             'Parallelizes over BAM regions; requires a sorted, indexed BAM. '
             'Ignored when --region is set. Auto-indexes if .bai is absent.',
    )
    from rectify.core.provenance.cli import add_resume_args
    add_resume_args(correct_cdna_parser)
