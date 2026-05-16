#!/usr/bin/env python3
"""Characterize the current walkback behavior on the bundled cDNA and
QuantSeq REV validation BAMs.

Calls the protocol-specific walkback wrappers (``walk_back_anchor_and_tail``
for cDNA, ``walkback_quantseq_rev`` for QSrev) on every primary read in the
bundled BAM and prints a TSV. The script is used to (a) lock in baseline
expected values for the per-arm regression tests, and (b) diff the output
between the current simple-core wiring and the planned guarded-core wiring
to see what enabling the artifact guards actually changes.

Usage::

    python scripts/validation_data/characterize_baseline.py --arm cdna
    python scripts/validation_data/characterize_baseline.py --arm quantseq_rev
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pysam

from rectify.utils.genome import load_genome
from rectify.config import CHROM_TO_GENOME


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "rectify" / "data"

CDNA_BAM = DATA / "validation" / "validation_reads_cdna.bam"
QSREV_BAM = DATA / "validation" / "validation_reads_quantseq_rev.bam"

# Use the bundled S288C reference; the cDNA + QSrev source BAMs all use it.
GENOME_FA = (
    DATA / "genomes" / "saccharomyces_cerevisiae"
    / "S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
)


def resolve_chrom_seq(genome, chrom):
    return genome.get(chrom) or genome.get(CHROM_TO_GENOME.get(chrom, ""), "")


def characterize_cdna(bam_path: Path, genome: dict) -> None:
    """Walk every primary read in the cDNA BAM through walk_back_anchor_and_tail.
    The cDNA pipeline uses ``orient`` (not is_reverse) — we derive it from the
    XU tag if present, else from is_reverse (proxy).
    """
    from rectify.core.cdna_correct_command import walk_back_anchor_and_tail

    print(
        "\t".join([
            "xv_label", "xg_category", "qname", "chrom", "is_reverse",
            "orient", "aln_start", "aln_end", "anchor_ref_pos", "tail_len",
        ])
    )
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            try:
                xv = r.get_tag("XV")
                xg = r.get_tag("XG")
            except KeyError:
                continue
            # Approximation: cDNA "orient" defaults to fwd if is_reverse=False,
            # rev if is_reverse=True. The real cdna pipeline derives orient
            # from SSP/UMI bridge detection, but for the simple-core test we
            # mirror DRS strand convention.
            orient = "fwd" if not r.is_reverse else "rev"
            chrom_seq = resolve_chrom_seq(genome, r.reference_name)
            if not chrom_seq:
                anchor, tail_len = -1, -1
            else:
                anchor, tail_len = walk_back_anchor_and_tail(r, chrom_seq, orient)
            print(
                "\t".join(map(str, [
                    xv, xg, r.query_name, r.reference_name,
                    int(r.is_reverse), orient,
                    r.reference_start, r.reference_end,
                    anchor, tail_len,
                ]))
            )


def characterize_quantseq_rev(bam_path: Path, genome: dict) -> None:
    """Walk every primary read in the QSrev BAM through walkback_quantseq_rev."""
    from rectify.core.correct.protocols.quantseq_rev import walkback_quantseq_rev

    print(
        "\t".join([
            "xv_label", "xg_category", "qname", "chrom", "is_reverse",
            "gene_strand", "aln_start", "aln_end", "original_3prime",
            "corrected_3prime", "applied", "correction_bp",
        ])
    )
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            try:
                xv = r.get_tag("XV")
                xg = r.get_tag("XG")
            except KeyError:
                continue
            chrom_seq = resolve_chrom_seq(genome, r.reference_name)
            if not chrom_seq:
                orig = corr = -1
                applied = "no_genome"
                gene_strand = "?"
            else:
                orig, corr, applied, gene_strand = walkback_quantseq_rev(r, chrom_seq)
            correction_bp = abs(orig - corr) if applied != "no_genome" else -1
            print(
                "\t".join(map(str, [
                    xv, xg, r.query_name, r.reference_name,
                    int(r.is_reverse), gene_strand,
                    r.reference_start, r.reference_end,
                    orig, corr, applied, correction_bp,
                ]))
            )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--arm", choices=["cdna", "quantseq_rev"], required=True,
        help="Which arm to characterize.",
    )
    args = p.parse_args()

    if not GENOME_FA.exists():
        print(f"Missing genome FASTA: {GENOME_FA}", file=sys.stderr)
        return 1
    genome = load_genome(str(GENOME_FA))

    if args.arm == "cdna":
        if not CDNA_BAM.exists():
            print(f"Missing cDNA BAM: {CDNA_BAM}", file=sys.stderr)
            return 1
        characterize_cdna(CDNA_BAM, genome)
    else:
        if not QSREV_BAM.exists():
            print(f"Missing QSrev BAM: {QSREV_BAM}", file=sys.stderr)
            return 1
        characterize_quantseq_rev(QSREV_BAM, genome)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
