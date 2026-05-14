#!/usr/bin/env python3
"""NMD-sensitive alternative-splicing audit.

For each BAM (wt + upf1Δ replicates of DRS, plus the wt_rep1 cDNA consensus
BAM), classify reads aligning to a panel of known NMD-AS yeast genes into
unspliced / annotated_spliced / alternative_spliced / novel_spliced.

The expected biological signal:
  - upf1Δ should ACCUMULATE NMD-sensitive isoforms vs wt.
  - For intron-retention NMD substrates (YRA1, PTC7): higher %unspliced in upf1Δ.
  - For alternative splicing NMD substrates: higher %alternative or %novel
    in upf1Δ for the NMD-targeted isoforms.

Usage (on H2, after `git pull`):

    /u/home/k/kevinroy/.conda/envs/rectify/bin/python3 \\
        /u/home/k/kevinroy/software/rectify/tools/nmd_as_splice_audit.py \\
        --out /u/scratch/k/kevinroy/nmd_as_audit/results.tsv
"""
from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import List, Optional, Tuple

from rectify.core.analyze.splice_summary import (
    GeneIntronSet,
    classify_bam_for_genes,
)


logger = logging.getLogger("nmd_as_audit")


# Gene panel — well-curated NMD-AS targets in S. cerevisiae.
# Columns: gene_name, sys_id, chrom, strand, [(intron_start, intron_end), ...]
GENE_PANEL = [
    ("YRA1",   "YDR381W",     "chrIV",  "+", [(1236843, 1237608)]),
    ("SUS1",   "YBR111W-A",   "chrII",  "+", [(462210, 462289), (462430, 462499)]),
    ("PTC7",   "YHR076W",     "chrVIII","+", [(251156, 251248)]),
    ("DBP2",   "YNL112W",     "chrXIV", "+", [(414912, 415913)]),
    ("RPL30",  "YGL030W",     "chrVII", "+", [(439094, 439323)]),
    ("RPL7B",  "YPL198W",     "chrXVI", "+", [(173163, 173571), (173666, 174072)]),
    ("GCR1",   "YPL075W",     "chrXVI", "+", [(412262, 413012)]),
    ("SRC1",   "YML034W",     "chrXIII","+", [(211445, 211570)]),
    ("RPL22B", "YFL034C-A",   "chrVI",  "-", [(64600, 64920)]),
    ("TFC3",   "YAL001C",     "chrI",   "-", [(151007, 151096)]),
]


# BAMs to audit. {sample_label: full_bam_path_on_H2}.
BAMS_H2 = {
    # DRS — 3 wt + 3 upf1Δ
    "drs_wt_rep1":    "/u/project/guillom/shared/raw/by4742-wt-upf1D_polya_drs_2025/wt_polya_rep1.sorted.bam",
    "drs_wt_rep2":    "/u/project/guillom/shared/raw/by4742-wt-upf1D_polya_drs_2025/wt_polya_rep2.sorted.bam",
    "drs_wt_rep3":    "/u/project/guillom/shared/raw/by4742-wt-upf1D_polya_drs_2025/wt_polya_rep3.sorted.bam",
    "drs_upf1d_rep1": "/u/project/guillom/shared/raw/by4742-wt-upf1D_polya_drs_2025/upf1d_polya_rep1.sorted.bam",
    "drs_upf1d_rep2": "/u/project/guillom/shared/raw/by4742-wt-upf1D_polya_drs_2025/upf1d_polya_rep2.sorted.bam",
    "drs_upf1d_rep3": "/u/project/guillom/shared/raw/by4742-wt-upf1D_polya_drs_2025/upf1d_polya_rep3.sorted.bam",
    # cDNA — the M1-produced v3 consensus BAM (one record per UMI cluster)
    "cdna_wt_rep1":   "/u/project/guillom/kevinroy/projects/ont_cdna/analyses/wt_rep1_v3_20260513/stage1_consensus.bam",
}


def resolve_bams(bams_dir: Optional[str]) -> List[Tuple[str, str]]:
    """Return [(label, path), ...].

    If ``bams_dir`` is given, look for ``<bams_dir>/<label>.nmdas.bam`` (the
    subset BAM naming used by the H2 extraction step). Otherwise return the
    full H2 paths in :data:`BAMS_H2` — useful when running on H2 itself.
    """
    if bams_dir:
        base = Path(bams_dir)
        return [(label, str(base / f"{label}.nmdas.bam")) for label in BAMS_H2]
    return [(label, path) for label, path in BAMS_H2.items()]


def build_gene_intron_sets() -> List[GeneIntronSet]:
    """Convert GFF (1-based inclusive) intron coords in GENE_PANEL into BAM
    (0-based half-open) coords used by N-cigars. The conversion is
    ``(gff_start, gff_end) -> (gff_start - 1, gff_end)`` — verified empirically
    against cDNA N-cigars at YRA1 (GFF 1236843,1237608 ↔ BAM 1236842,1237608).
    """
    out = []
    for name, sys_id, chrom, strand, gff_introns in GENE_PANEL:
        bam_introns = frozenset((s - 1, e) for s, e in gff_introns)
        out.append(
            GeneIntronSet(
                gene_id=f"{name} ({sys_id})",
                chrom=chrom,
                strand=strand,
                introns=bam_introns,
            )
        )
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="Output TSV path (long format)")
    ap.add_argument(
        "--bams-dir",
        default=None,
        help="Look for <label>.nmdas.bam in this dir instead of the full H2 paths. "
        "Use for local M1 iteration after pulling subset BAMs.",
    )
    ap.add_argument(
        "--snap-tolerance",
        type=int,
        default=0,
        help="Snap N-ops within N bp of an annotated intron (both ends). 0 = no snapping (strict).",
    )
    ap.add_argument("--min-mapq", type=int, default=0, help="Skip reads below this mapping quality")
    ap.add_argument(
        "--span-margin",
        type=int,
        default=5,
        help="Read must overlap intron by this many flanking bp to count as spanning",
    )
    args = ap.parse_args()
    bams = resolve_bams(args.bams_dir)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    genes = build_gene_intron_sets()
    rows = []

    for sample_label, bam_path in bams:
        if not Path(bam_path).exists():
            logger.warning("Missing BAM, skipping: %s", bam_path)
            continue
        logger.info("Processing %s -> %s", sample_label, bam_path)
        per_gene = classify_bam_for_genes(
            bam_path,
            genes,
            sample_name=sample_label,
            snap_tolerance=args.snap_tolerance,
            span_margin=args.span_margin,
            min_mapq=args.min_mapq,
        )
        for st in per_gene:
            rows.append(st.as_row())
        # Quick console summary per sample
        n_total = sum(st.total_spanning for st in per_gene)
        logger.info("  %s: %d total spanning reads across %d genes", sample_label, n_total, len(genes))

    # Write long-format TSV
    if not rows:
        logger.error("No rows to write — all BAMs missing?")
        sys.exit(1)
    with out_path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d rows to %s", len(rows), out_path)

    # Console pivot: per-gene pct_unspliced wt vs upf1Δ for quick eyeballing
    print("\n=== Quick console pivot: pct_unspliced (intron retention signal) ===")
    print(f"{'gene':<25} ", end="")
    for label, _ in bams:
        print(f"{label:>16} ", end="")
    print()
    for g in genes:
        print(f"{g.gene_id:<25} ", end="")
        for label, _ in bams:
            match = next((r for r in rows if r["gene_id"] == g.gene_id and r["sample"] == label), None)
            if match:
                print(f"{match['pct_unspliced']:>14.1f}%  ", end="")
            else:
                print(f"{'--':>16} ", end="")
        print()


if __name__ == "__main__":
    main()
