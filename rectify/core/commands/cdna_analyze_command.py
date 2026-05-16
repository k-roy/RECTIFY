#!/usr/bin/env python3
"""rectify cdna-analyze — post-alignment isoform clustering for cDNA pipeline.

Reads the BAM produced by `rectify align` on the per-cluster consensus FASTQ
output by `rectify correct-cdna`. Each BAM record is one UMI cluster's
consensus aligned via the multi-aligner (minimap2 + mapPacBio + gapmm2) and
chimeric-consensus integration.

Recomputes tail_len (XA), 5' TSS correction (pos5_corrected), gene assignment
(XG), sense/antisense classification (XS), isoform clustering (XI), and
Type-1↔Type-2 same-orient pairing (XL) using the post-align coordinates —
which are more accurate than the pre-align values that were used for UMI
bucketing in `correct-cdna`.

Outputs (in --out directory):
  - clusters.tsv          per-cluster manifest with corrected positions
  - isoforms.tsv          isoform-level aggregation
  - t1t2_pairs.tsv        Type-1↔Type-2 reconciliation pairs
  - consensus_tagged.bam  the input BAM rewritten with the newly computed
                          XA/XG/XS/XI/XL tags added to each record. Indexed.
"""
from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from rectify.core.cdna.gff import (
    classify_sense_antisense,
    load_gff_genes,
)
from rectify.core.cdna.isoforms import (
    assign_isoforms,
    reconcile_t1_t2_pairs,
)
from rectify.core.cdna.read_info import ReadInfo
from rectify.core.cdna.walkback import (
    walk_back_anchor_and_tail,
    walk_forward_tss,
)


def _read_info_from_bam_record(rec: pysam.AlignedSegment,
                                chrom_seq: str,
                                ) -> Optional[Tuple[ReadInfo, int]]:
    """Build a synthetic ReadInfo from one post-align consensus BAM record.

    Returns (read_info, cluster_size) or None if the record is unmapped or
    missing required tags.  cluster_size comes from the XC tag (count of
    original input reads that contributed to this consensus).
    """
    if rec.is_unmapped or not rec.cigartuples:
        return None
    try:
        umi = rec.get_tag("XU")
        orient = rec.get_tag("XO")
        read_type = int(rec.get_tag("XT"))
        read_subtype = rec.get_tag("XY")
        xc = int(rec.get_tag("XC"))
        xf = int(rec.get_tag("XF"))
    except KeyError:
        return None

    chrom = rec.reference_name
    aln_start = rec.reference_start
    aln_end = rec.reference_end or aln_start

    # Recompute tail_len + canonical cleavage anchor via walkback on post-align coords
    anchor, tail_len = walk_back_anchor_and_tail(rec, chrom_seq, orient)

    # Recompute pos5_corrected via bridge-G walk-forward
    pos5_corrected = walk_forward_tss(rec, orient, chrom_seq)

    info = ReadInfo(
        read_id=rec.query_name,
        chrom=chrom,
        anchor=anchor,
        orient=orient,
        umi=umi,
        is_reverse=rec.is_reverse,
        xf_tier=xf,
        tail_len=tail_len,
        aln_start=aln_start,
        aln_end=aln_end,
        read_type=read_type,
        pos5_corrected=pos5_corrected,
        read_subtype=read_subtype,
    )
    return (info, xc)


def run(args) -> int:
    """rectify cdna-analyze entry point."""
    args.bam = Path(args.bam) if not isinstance(args.bam, Path) else args.bam
    args.out = Path(args.out) if not isinstance(args.out, Path) else args.out
    if getattr(args, "gff", None) is not None:
        args.gff = Path(args.gff) if not isinstance(args.gff, Path) else args.gff
    if getattr(args, "reference", None) is not None:
        args.reference = Path(args.reference) if not isinstance(args.reference, Path) else args.reference

    for attr, default in [
        ("min_gene_frac", 0.3),
        ("min_read_frac", 0.8),
        ("isoform_tol_5", 5),
        ("isoform_tol_3", 5),
        ("t1t2_tol_5", 5),
        ("t1t2_tol_3", 5),
    ]:
        if not hasattr(args, attr):
            setattr(args, attr, default)

    args.out.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger("cdna-analyze")
    t0 = time.time()

    if args.reference is None:
        log.error("--reference is required (walkback needs the genome FASTA)")
        return 1
    if args.gff is None:
        log.error("--gff is required (gene assignment needs the annotation)")
        return 1

    log.info("Loading reference FASTA: %s", args.reference)
    chrom_cache: Dict[str, str] = {}
    with pysam.FastaFile(str(args.reference)) as fa:
        for c in fa.references:
            chrom_cache[c] = fa.fetch(c).upper()

    log.info("Loading GFF: %s", args.gff)
    gene_trees = load_gff_genes(args.gff)

    log.info("Streaming consensus BAM: %s", args.bam)
    clusters: List[List[ReadInfo]] = []
    cluster_umi: Dict[int, str] = {}
    cluster_xf: Dict[int, int] = {}
    cluster_xm: Dict[int, str] = {}
    cluster_xb: Dict[int, str] = {}
    cluster_xr: Dict[int, str] = {}
    cluster_tail_len: Dict[int, int] = {}
    # cid keyed by query_name so we can re-emit the tagged BAM on a second pass.
    qname_to_cid: Dict[str, int] = {}

    n_in = n_unmapped = n_missing_tags = 0

    with pysam.AlignmentFile(str(args.bam), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            n_in += 1
            if rec.is_unmapped:
                n_unmapped += 1
                continue
            chrom_seq = chrom_cache.get(rec.reference_name, "")
            result = _read_info_from_bam_record(rec, chrom_seq)
            if result is None:
                n_missing_tags += 1
                continue
            info, xc = result
            cid = len(clusters)
            # Inflate to cluster_size so existing weight=len(clusters[cid]) is correct
            clusters.append([info] * xc)
            cluster_umi[cid] = info.umi
            cluster_xf[cid] = info.xf_tier
            cluster_tail_len[cid] = info.tail_len
            qname_to_cid[rec.query_name] = cid
            for tag, store in (("XM", cluster_xm), ("XB", cluster_xb), ("XR", cluster_xr)):
                try:
                    store[cid] = rec.get_tag(tag)
                except KeyError:
                    store[cid] = ""

    log.info("  %d records read (%d unmapped, %d missing required tags)",
             n_in, n_unmapped, n_missing_tags)
    log.info("  %d usable consensus clusters", len(clusters))

    # Gene / sense-antisense classification on post-align positions
    t1 = time.time()
    cluster_xs: Dict[int, str] = {}
    cluster_xg: Dict[int, Optional[str]] = {}
    for cid, c in enumerate(clusters):
        r0 = c[0]
        xs, gene = classify_sense_antisense(
            r0.chrom, r0.aln_start, r0.aln_end, r0.orient, gene_trees,
            min_gene_frac=args.min_gene_frac, min_read_frac=args.min_read_frac)
        cluster_xs[cid] = xs
        cluster_xg[cid] = gene
    log.info("Gene/strand classification: %d clusters in %.1fs",
             len(clusters), time.time() - t1)

    # Isoform clustering
    t2 = time.time()
    cluster_xi, iso_records = assign_isoforms(
        clusters, cluster_xg, args.isoform_tol_5, args.isoform_tol_3)
    log.info("Isoform clustering (tol5=%d, tol3=%d): %d isoforms from %d annotated clusters (%.1fs)",
             args.isoform_tol_5, args.isoform_tol_3, len(iso_records),
             len(cluster_xi), time.time() - t2)

    # Type-1 ↔ Type-2 same-orient pairing
    t3 = time.time()
    cluster_xl, t1t2_pairs = reconcile_t1_t2_pairs(
        clusters, cluster_xg, args.t1t2_tol_5, args.t1t2_tol_3)
    log.info("Type-1↔Type-2 reconciliation (tol5=%d, tol3=%d): %d pairs (%.1fs)",
             args.t1t2_tol_5, args.t1t2_tol_3, len(t1t2_pairs), time.time() - t3)

    # ── Write manifests ──────────────────────────────────────────────────────
    manifest = args.out / "clusters.tsv"
    with manifest.open("w") as f:
        f.write("cluster_id\tchrom\tanchor\torient\tn_reads\tumi_canonical"
                "\txs\txf\txt\tread_subtype\ttail_len\tgene\tisoform_id"
                "\taln_start\taln_end\tread_ids\n")
        for cid, c in enumerate(clusters):
            gene_str = cluster_xg.get(cid) or ""
            iso_str = cluster_xi.get(cid) or ""
            r0 = c[0]
            f.write(f"{cid}\t{r0.chrom}\t{r0.anchor}\t{r0.orient}"
                    f"\t{len(c)}\t{cluster_umi[cid]}\t{cluster_xs[cid]}"
                    f"\t{cluster_xf[cid]}\t{r0.read_type}\t{r0.read_subtype}"
                    f"\t{cluster_tail_len[cid]}"
                    f"\t{gene_str}\t{iso_str}\t{r0.aln_start}\t{r0.aln_end}"
                    f"\t{cluster_xr[cid]}\n")
    log.info("Wrote cluster manifest → %s", manifest)

    iso_path = args.out / "isoforms.tsv"
    with iso_path.open("w") as f:
        f.write("isoform_id\tgene\tchrom\torient\tread_type\tread_subtype\tn_clusters\tn_reads_total"
                "\tpos5_modal\tpos3_modal\ttail_len_median\tcluster_ids\n")
        for r in iso_records:
            f.write(f"{r['isoform_id']}\t{r['gene']}\t{r['chrom']}\t{r['orient']}"
                    f"\t{r['read_type']}\t{r['read_subtype']}\t{r['n_clusters']}\t{r['n_reads_total']}"
                    f"\t{r['pos5_modal']}\t{r['pos3_modal']}\t{r['tail_len_median']}"
                    f"\t{r['cluster_ids']}\n")
    log.info("Wrote isoform manifest → %s (%d isoforms)", iso_path, len(iso_records))

    pairs_tsv = args.out / "t1t2_pairs.tsv"
    with pairs_tsv.open("w") as f:
        f.write("t1_cid\tt2_cid\tgene\torient\tt1_pos5\tt1_pos3\tt2_pos5\tt2_pos3"
                "\td5\td3\tt1_n_reads\tt2_n_reads\tt1_umi\n")
        for p in t1t2_pairs:
            f.write(f"{p['t1_cid']}\t{p['t2_cid']}\t{p['gene']}\t{p['orient']}"
                    f"\t{p['t1_pos5']}\t{p['t1_pos3']}\t{p['t2_pos5']}\t{p['t2_pos3']}"
                    f"\t{p['d5']}\t{p['d3']}\t{p['t1_n_reads']}\t{p['t2_n_reads']}"
                    f"\t{p['t1_umi']}\n")
    log.info("Wrote T1↔T2 pairs → %s (%d pairs)", pairs_tsv, len(t1t2_pairs))

    # ── Tagged BAM ───────────────────────────────────────────────────────────
    # Rewrite the input BAM with the new XA/XG/XS/XI/XL/XA' (corrected tail
    # length) tags added per record. Useful for downstream consumers that
    # prefer SAM-tag access over TSV joins. Records that didn't survive the
    # tag/walkback gates (unmapped or missing required tags) pass through
    # unmodified.
    t4 = time.time()
    tagged_bam = args.out / "consensus_tagged.bam"
    n_tagged = 0
    with pysam.AlignmentFile(str(args.bam), "rb") as src, \
         pysam.AlignmentFile(str(tagged_bam), "wb", template=src) as dst:
        for rec in src.fetch(until_eof=True):
            cid = qname_to_cid.get(rec.query_name)
            if cid is not None:
                rec.set_tag("XA", int(cluster_tail_len[cid]), value_type="i")
                xs = cluster_xs.get(cid, "unannotated")
                rec.set_tag("XS", xs, value_type="Z")
                gene = cluster_xg.get(cid)
                if gene:
                    rec.set_tag("XG", gene, value_type="Z")
                iso = cluster_xi.get(cid)
                if iso:
                    rec.set_tag("XI", iso, value_type="Z")
                partner = cluster_xl.get(cid)
                if partner is not None:
                    rec.set_tag("XL", int(partner), value_type="i")
                n_tagged += 1
            dst.write(rec)
    pysam.index(str(tagged_bam))
    log.info("Wrote tagged BAM → %s (%d records tagged, %.1fs)",
             tagged_bam, n_tagged, time.time() - t4)

    log.info("Total runtime: %.1fs", time.time() - t0)
    return 0


def create_cdna_analyze_parser(subparsers):
    """Wire the `cdna-analyze` subcommand into the given subparsers group."""
    import argparse
    # =========================================================================
    # cdna-analyze command (post-align isoform clustering for cDNA)
    # =========================================================================
    cdna_analyze_parser = subparsers.add_parser(
        'cdna-analyze',
        help='Post-alignment isoform clustering for the cDNA pipeline (consumes `rectify align` output of `correct-cdna` consensus FASTQ)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            'Reads the multi-aligner consensus BAM produced by `rectify align` on the '
            'per-cluster consensus FASTQ from `rectify correct-cdna`. Each BAM record '
            'is one UMI-deduplicated molecule; tags XU/XC/XR/XO/XT/XY/XF/XB ride along '
            'from the FASTQ comments via `minimap2 -y`.\n\n'
            'Recomputes tail_len (XA) and pos5_corrected via walkback / walk-forward on '
            'the post-align CIGAR, then runs gene assignment, isoform clustering '
            '(directional, tol5/tol3), and Type-1↔Type-2 same-orient pairing on the '
            'post-align coordinates.\n\n'
            'Output: clusters.tsv (per-cluster manifest), isoforms.tsv (isoform-level '
            'aggregation), t1t2_pairs.tsv (T1↔T2 reconciliation), consensus_tagged.bam '
            '(input BAM rewritten with the new XA/XG/XS/XI/XL tags).'
        ),
    )
    cdna_analyze_parser.add_argument(
        'bam',
        help='Consensus BAM from `rectify align` on the correct-cdna FASTQ',
    )
    cdna_analyze_parser.add_argument(
        '-o', '--out',
        dest='out',
        required=True,
        help='Output directory (clusters.tsv, isoforms.tsv, t1t2_pairs.tsv, consensus_tagged.bam)',
    )
    cdna_analyze_parser.add_argument(
        '--gff',
        required=True,
        help='Genome annotation (GFF3) for gene tree / XS classification',
    )
    cdna_analyze_parser.add_argument(
        '--reference',
        required=True,
        help='Genome FASTA — required for walkback (XA tail_len recomputation)',
    )
    cdna_analyze_parser.add_argument(
        '--min-gene-frac', dest='min_gene_frac', type=float, default=0.3,
        help='XS classifier: gene_overlap / gene_length threshold',
    )
    cdna_analyze_parser.add_argument(
        '--min-read-frac', dest='min_read_frac', type=float, default=0.8,
        help='XS classifier: gene_overlap / read_aln_length threshold',
    )
    cdna_analyze_parser.add_argument(
        '--isoform-tol-5', dest='isoform_tol_5', type=int, default=5,
        help='Isoform clustering bp tolerance on 5\' axis (Type-1 only)',
    )
    cdna_analyze_parser.add_argument(
        '--isoform-tol-3', dest='isoform_tol_3', type=int, default=5,
        help='Isoform clustering bp tolerance on 3\' axis',
    )
    cdna_analyze_parser.add_argument(
        '--t1t2-tol-5', dest='t1t2_tol_5', type=int, default=5,
        help='T1↔T2 reconciliation bp tolerance on 5\' axis',
    )
    cdna_analyze_parser.add_argument(
        '--t1t2-tol-3', dest='t1t2_tol_3', type=int, default=5,
        help='T1↔T2 reconciliation bp tolerance on 3\' axis',
    )
