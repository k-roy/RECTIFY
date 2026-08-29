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
  XN  oriented (always 1): the consensus is emitted RNA-sense, so after alignment
       is_reverse IS the gene strand and minimap2 -uf is valid (planning/730, 2026-08-21)

Usage (via rectify CLI):
    rectify correct-cdna INPUT.bam --out OUTDIR [options]
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
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


def _print_pretrim_health(fastq_stats: Dict) -> None:
    """Report adapter-pretrim no-op rates (planning/681).

    A silently no-op'd trim sends ~124 nt of adapter + random 27-nt UMI + poly-A into the
    aligner as a soft clip; the junction resolver then enumerates hundreds of GT/AG
    candidates against sequence that can never align, and — because a random UMI is
    maximally high-complexity — no complexity gate rejects it. That cost months as an
    unexplained "1-4 reads/s cDNA" pathology. Printing the rate makes a regression a number
    instead of a stall.

    Expected floor for ``5' trim NO-OP``: Type-2 molecules have no SSP by design, and
    pileup/POA consensuses are built from aligned positions only and so carry no adapter to
    strip at all. ``frame vs XO label mismatch`` should be ~0 once the caller propagates the
    frame correctly — a nonzero rate means the label and the sequence have drifted apart
    again, which is the original defect.
    """
    n = max(1, fastq_stats.get("written", 0))
    print("Adapter pretrim health:")
    for label, key in (("frame flipped (basecalled seq)", "trim_frame_flipped"),
                       ("frame vs XO label mismatch", "trim_frame_mismatch"),
                       ("5' trim NO-OP (Type-1 only)", "trim_noop_5p"),
                       ("3' trim NO-OP", "trim_noop_3p")):
        v = fastq_stats.get(key, 0)
        print(f"  {label:<32s} {v:>8d}  ({100 * v / n:.1f}%)")


def _region_cluster_prefix(region: Optional[str]) -> str:
    """Stage-1 cluster-name prefix that makes per-region output globally unique.

    The bare ``cluster_<cid>`` counter restarts at 0 in every ``correct-cdna``
    invocation, so per-region (per-chromosome) runs whose FASTQs are later
    concatenated for align2 collide (``cluster_0`` from chrI vs chrII vs …). The
    align2 K-way consensus merge keys on the name/RN and would collapse those
    distinct molecules to one record (the ~87% cDNA align2 loss; see
    planning/251, /250a-c). Prefixing with the region makes concatenated names
    globally unique. Sanitized to a normalize-safe token so ``chrI:1-1000`` ->
    ``cluster_chrI_1_1000`` (not ``cluster_chrI:1-1000``), which
    ``_normalize_bam_read_name`` leaves intact.
    """
    if not region:
        return "cluster"
    tok = re.sub(r"[^A-Za-z0-9]+", "_", str(region)).strip("_")
    return f"cluster_{tok}" if tok else "cluster"


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
    adaptive_threshold: int = 0,
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

    # ---- 624: ANCHOR-IN-REGION FILTER — makes --workers N EXACTLY equal --workers 1 ----
    # `stream_reads` fetches via pysam, which returns every read OVERLAPPING the interval,
    # and regions are contiguous. Two failure modes follow, and they are NOT the same bug:
    #
    #   (a) DUPLICATION — a read overlapping two regions is clustered in both and emitted
    #       twice. Fixed by any partition-by-position rule.
    #   (b) SPLITTING   — two reads of the SAME molecule land in different regions, are never
    #       compared, and the molecule is emitted as two. Measured on a 440,009-read subset:
    #       partitioning on `aln_start` conserved reads exactly (436,013 both arms) but still
    #       produced +293 molecules (+0.096%), all of this kind.
    #
    # planning/584b attributed its "+7 molecules on 63,508" to double-counting; the read-
    # conservation check above shows that diagnosis was wrong — nothing is duplicated, molecules
    # are fragmented.
    #
    # THE FIX: partition on the CLUSTERING KEY, not on the alignment start. `cluster_reads`
    # buckets on (chrom, anchor // anchor_window, orient, read_type) — so at anchor_window == 1
    # two reads co-cluster ONLY IF their anchors are IDENTICAL. Assigning each read to the region
    # containing its ANCHOR therefore keeps every member of a bucket in the same region by
    # construction: no duplication AND no splitting.
    #
    # Safety: the anchor is derived by walking back from the read's 3' end, so it always lies
    # inside the aligned span — a read whose anchor is in this region necessarily overlaps it and
    # is therefore returned by this fetch. And `get_processing_regions` tiles every chromosome
    # disjointly over [0, chrom_len), so each anchor falls in exactly one region.
    #
    # ⚠️ At anchor_window > 1 a bucket spans `anchor_window` positions, so a bucket straddling a
    # region boundary can still split. Exactness is guaranteed only for anchor_window == 1 (the
    # lab's production setting, planning/584d). Residual splitting is reported below.
    #
    # `region_str` is f"{chrom}:{start+1}-{end}" (1-based inclusive); rsplit(":", 1) so
    # chromosome names containing ':' are handled.
    try:
        _span = region_str.rsplit(":", 1)[1]
        _r_start0, _r_end = _span.split("-")
        _r_start0 = int(_r_start0) - 1          # back to 0-based, matching ref coords
        _r_end = int(_r_end)
        _n_before = len(reads)
        reads = [r for r in reads if _r_start0 <= r.anchor < _r_end]
        _n_boundary = _n_before - len(reads)
    except (IndexError, ValueError):
        # Unparseable region (e.g. a whole-chromosome string with no coords): no filtering.
        # Serial mode passes region=None and never reaches this function, so this is only a
        # guard against a future caller shape, not an expected path.
        _n_boundary = 0

    if not reads:
        # Write empty gzip FASTQ so the merge step has a file to cat
        import gzip as _gz
        with _gz.open(out_fastq_path, "wt"):
            pass
        return {"input_reads": 0, "written": 0, "from_singletons": 0,
                "from_multi_pileup": 0, "from_multi_fallback": 0, "n_rdna_masked": 0}

    clusters, stats = _cluster(reads, anchor_window_bp, umi_edit_distance, per_cluster_cap,
                                clustering_method=umi_clustering,
                                adaptive_threshold=adaptive_threshold)

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
    # 624: report what the start-in-region filter removed, so the correction is visible in the
    # run log rather than silent (workspace rule: never a silent cap).
    fq_stats["n_boundary_dropped"] = _n_boundary
    # 658: surface the per-region adaptive-ed decisions so the parallel path is not blind to them.
    fq_stats["adaptive_deep_buckets"] = stats.get("adaptive_deep_buckets", 0)
    fq_stats["adaptive_deep_reads"] = stats.get("adaptive_deep_reads", 0)
    # QC: these were computed here and then DISCARDED, so --workers>1 (the path any real
    # run takes) shipped without the read-type / XF / tail QC the serial path prints.
    # Fold them into the return so the parent can sum them across regions.
    from rectify.core.cdna.qc import TAIL_BINS as _TB
    fq_stats["input_reads"] = len(reads)
    for _k in ("type1_reads", "type2_reads", "type1_clusters", "type2_clusters",
               "buckets_dropped_polyA_pileup", "reads_in_dropped_buckets"):
        fq_stats[_k] = stats.get(_k, 0)
    fq_stats["subtype_reads"] = dict(stats.get("subtype_reads", {}) or {})
    _tiers: Dict[int, int] = {0: 0, 1: 0, 2: 0}
    for _v in cluster_xf_tier.values():
        _tiers[int(_v)] = _tiers.get(int(_v), 0) + 1
    fq_stats["xf_tier_counts"] = _tiers
    _hist = [0] * len(_TB)
    for _cid in range(len(clusters)):
        if cluster_xf_tier.get(_cid, 0) < 2:
            continue
        _tl = cluster_tail_len.get(_cid, 0)
        for _i, (_lo, _hi) in enumerate(_TB):
            if _lo <= _tl < _hi:
                _hist[_i] += 1
                break
    fq_stats["tail_hist"] = _hist
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
    adaptive_threshold: int = 0,
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
                    adaptive_threshold,
                )
                futures[fut] = plan.region_id

            # 624: n_boundary_dropped accumulates the reads removed by the start-in-region
            # filter — i.e. exactly the double-counting that --workers > 1 used to introduce.
            total_stats: Dict = {"input_reads": 0, "written": 0, "from_singletons": 0,
                                  "from_multi_pileup": 0, "from_multi_fallback": 0,
                                  "n_boundary_dropped": 0,
                                  "adaptive_deep_buckets": 0, "adaptive_deep_reads": 0,
                                  # planning/681 adapter-pretrim health counters
                                  "trim_frame_flipped": 0, "trim_frame_mismatch": 0,
                                  "trim_noop_5p": 0, "trim_noop_3p": 0}
            total_stats.update({"input_reads": 0, "type1_reads": 0, "type2_reads": 0,
                                "type1_clusters": 0, "type2_clusters": 0,
                                "buckets_dropped_polyA_pileup": 0,
                                "reads_in_dropped_buckets": 0, "n_rdna_masked": 0})
            agg_tiers: Dict[int, int] = {0: 0, 1: 0, 2: 0}
            agg_tail: List[int] = []
            agg_sub: Dict[str, int] = {}
            for fut in as_completed(futures):
                s = fut.result()
                for k in total_stats:
                    total_stats[k] += s.get(k, 0)
                for _t, _n in (s.get("xf_tier_counts") or {}).items():
                    agg_tiers[int(_t)] = agg_tiers.get(int(_t), 0) + _n
                _h = s.get("tail_hist") or []
                if not agg_tail:
                    agg_tail = list(_h)
                elif _h:
                    agg_tail = [a + b for a, b in zip(agg_tail, _h)]
                for _k, _n in (s.get("subtype_reads") or {}).items():
                    agg_sub[_k] = agg_sub.get(_k, 0) + _n
            total_stats["xf_tier_counts"] = agg_tiers
            total_stats["tail_hist"] = agg_tail
            total_stats["subtype_reads"] = agg_sub
            if adaptive_threshold > 0:
                log.info("  adaptive ed: %d deep buckets (%d reads) clustered at ed=1 "
                         "(threshold %d reads/bucket)",
                         total_stats["adaptive_deep_buckets"],
                         total_stats["adaptive_deep_reads"], adaptive_threshold)
            if total_stats["n_boundary_dropped"]:
                log.info(
                    "  region-boundary de-duplication: %d overlapping read placements "
                    "excluded (start-in-region filter; makes workers=%d identical to workers=1)",
                    total_stats["n_boundary_dropped"], workers,
                )

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
        ("umi_adaptive_threshold", 0),
        ("umi_profile", False),
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

    # 658: per-library UMI error-profile calibration (spacer ground truth; capped
    # pre-pass, ~1-2 min). Written as a sidecar for provenance/QC and downstream
    # calibrated merge rules; does not change clustering behaviour by itself.
    if args.umi_profile:
        import json as _json
        from rectify.core.cdna.profile import derive_umi_error_profile
        _t_prof = time.time()
        _prof = derive_umi_error_profile(str(args.bam))
        _prof_path = args.out / "umi_error_profile.json"
        with open(_prof_path, "w") as _fh:
            _json.dump(_prof, _fh, indent=1)
        for _reg, _r in _prof.get("regimes", {}).items():
            log.info("UMI profile [%s]: n=%d clean=%.1f%% frameshift=%.2f%% "
                     "sub/base=%.3f%% meanQ=%.1f", _reg, _r["n"], 100 * _r["clean"],
                     100 * _r["frameshift"], 100 * _r["inframe_sub_per_base"],
                     _r["mean_umi_Q"])
        log.info("UMI error profile -> %s (%.0fs)", _prof_path, time.time() - _t_prof)

    if args.umi_adaptive_threshold > 0:
        log.info("Adaptive UMI edit distance ON: ed=%d below %d reads/bucket, ed=1 at/above "
                 "(planning/654 policy)", args.umi_edit_distance, args.umi_adaptive_threshold)

    use_poa = HAS_POA and not args.no_poa
    out_fastq = args.out / "stage1_consensus.fastq.gz"
    _workers = max(1, getattr(args, "workers", 1) or 1)

    if _workers > 1 and not args.region:
        log.info("Parallel stage-1 (%d workers) → %s", _workers, out_fastq)
        fastq_stats = _run_cdna_correct_parallel(
            args.bam, out_fastq, args.reference, rdna_intervals,
            args.anchor_window_bp, args.umi_edit_distance, args.per_cluster_cap,
            args.umi_clustering, use_poa, args.strand_aware_consensus, _workers,
            adaptive_threshold=args.umi_adaptive_threshold,
        )
        from rectify.core.cdna.qc import collect_qc, render_qc, write_qc_json
        _qc = collect_qc(
            fastq_stats=fastq_stats, stats=fastq_stats,
            n_clusters=fastq_stats.get("written", 0),
            tier_counts=fastq_stats.get("xf_tier_counts"),
            tail_hist=fastq_stats.get("tail_hist"),
            n_input_reads=fastq_stats.get("input_reads"),
            workers=_workers, sample=args.out.name,
        )
        print()
        print(render_qc(_qc))
        _qp = write_qc_json(_qc, args.out)
        print()
        print(f"Output FASTQ: {out_fastq}")
        print(f"Stage-1 QC:   {_qp}")
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
                                         clustering_method=args.umi_clustering,
                                         adaptive_threshold=args.umi_adaptive_threshold)
        log.info("  %d molecule clusters (%.1fs)  biggest bucket=%d reads",
                 stats["molecule_clusters"], time.time() - t1, stats["biggest_bucket_size"])
        if stats.get("adaptive_deep_buckets"):
            log.info("  adaptive ed: %d deep buckets (%d reads) clustered at ed=1; "
                     "the rest at ed=%d", stats["adaptive_deep_buckets"],
                     stats["adaptive_deep_reads"], args.umi_edit_distance)

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

        # When restricted to a single region, namespace cluster names with the
        # region so that per-region invocations (the production per-chromosome
        # `correct-cdna --region chrX` pattern) produce GLOBALLY-UNIQUE read
        # names after the per-region FASTQs are concatenated for align2. Without
        # this the bare `cluster_<cid>` counter restarts at 0 in every region and
        # the concatenated FASTQ collides ~7.5-way; the align2 K-way consensus
        # merge then keys on the (colliding) name/RN and collapses ~7.5 distinct
        # molecules to one, silently dropping ~87% of the reads (see planning/251,
        # /250a-c). The internal multi-region parallel path already prefixes with
        # region_id (_cdna_region_task); this makes the external single-region
        # path consistent.
        _cluster_prefix = _region_cluster_prefix(getattr(args, "region", None))
        fastq_stats = write_stage1_fastq(args.bam, out_fastq, clusters, umi_canon,
                                          cluster_xf_tier, cluster_tail_len,
                                          use_poa=use_poa,
                                          cluster_name_prefix=_cluster_prefix,
                                          strand_aware_consensus=args.strand_aware_consensus,
                                          reference=args.reference)
        log.info("  wrote %d records (%d singletons, %d pileup, %d rep fallback) in %.1fs",
                 fastq_stats["written"], fastq_stats["from_singletons"],
                 fastq_stats["from_multi_pileup"], fastq_stats["from_multi_fallback"],
                 time.time() - t2)

        from rectify.core.cdna.qc import collect_qc, render_qc, write_qc_json
        _qc = collect_qc(
            fastq_stats=fastq_stats, stats=stats,
            n_clusters=len(clusters),
            cluster_xf_tier=cluster_xf_tier, cluster_tail_len=cluster_tail_len,
            n_input_reads=len(reads), workers=1, sample=args.out.name,
        )
        print()
        print(render_qc(_qc))
        _qp = write_qc_json(_qc, args.out)
        print()
        print(f"Output FASTQ: {out_fastq}")
        print(f"Stage-1 QC:   {_qp}")
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
# repo at ~/work/ont_cdna/src/cdna_correct.py.


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
            'with alignment-independent SAM-tag comments (XU/XO/XC/XR/XM/XF/XA/XT/XY/XQ/XK/XB/XN) '
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
        '--umi-adaptive-threshold',
        dest='umi_adaptive_threshold',
        type=int,
        default=0,
        help='Depth-adaptive UMI edit distance (0 = off). Type-1 anchor buckets with >= N '
             'reads cluster at ed=1; buckets below N at --umi-edit-distance. Ground-truth '
             'pricing (Chanfreau planning/652+654): over-merge at ed>=2 grows with bucket '
             'depth without bound while under-merge stays ~+2%% at any depth. Production '
             'policy: --umi-edit-distance 2 --umi-adaptive-threshold 5000.',
    )
    correct_cdna_parser.add_argument(
        '--umi-profile',
        dest='umi_profile',
        action='store_true',
        default=False,
        help='Derive the per-library UMI error profile from the fixed-T spacers before '
             'clustering (capped pre-pass, ~1-2 min) and write umi_error_profile.json '
             'beside the stage-1 output. Per-regime (pore-entry vs pore-exit) clean/'
             'frameshift/substitution rates — provenance + QC + input for calibrated '
             'merge rules (planning/650/655).',
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
