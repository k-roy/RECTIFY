"""GFF parsing + gene-strand classification for the cdna pipeline.

Two responsibilities:
  - rDNA interval loading (used by stream_reads to mask rRNA_gene reads before
    UMI clustering, preventing the O(n²) Levenshtein blow-up on chrXII tandem
    repeats)
  - Gene-tree loading + overlap-based sense/antisense classification (XS) with
    first-TSS sense tiebreaker (v1.14).
"""
from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from intervaltree import IntervalTree


def _parse_gff_gene_name(attrs: str) -> Optional[str]:
    """Extract a stable gene identifier from a GFF attributes column.

    Preference order (matches SGD GFF conventions + cross-project handoff needs):
      1. ID=        — systematic name (e.g. YAL038W); the canonical join key
      2. Name=      — same as ID for SGD
      3. gene=      — common name (e.g. CDC19); fallback only
    Returns None if no identifier found.
    """
    for kv in attrs.split(";"):
        if kv.startswith("ID="):
            return kv[3:].strip()
    for kv in attrs.split(";"):
        if kv.startswith("Name="):
            return kv[5:].strip()
    for kv in attrs.split(";"):
        if kv.startswith("gene="):
            return kv[5:].strip()
    return None


def load_rdna_intervals(gff_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """Parse rRNA_gene features from GFF → {chrom: [(start, end), ...]} (0-based half-open).

    Used for rDNA masking: reads overlapping these intervals are dropped in
    stream_reads() before UMI extraction and clustering, preventing the O(n²)
    Levenshtein bottleneck that occurs when thousands of rDNA reads share
    identical anchor sequences (observed: 261k chrXII reads → 4h42m clustering).
    """
    result: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "rRNA_gene":
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1  # GFF is 1-based inclusive → 0-based half-open
            end = int(parts[4])
            result[chrom].append((start, end))
    return dict(result)


def load_gff_genes(gff_path: Path) -> Dict[str, IntervalTree]:
    """Load gene-like features → chrom → IntervalTree.

    Each interval's data is a tuple `(strand, gene_name)` for both XS
    classification and downstream join-key emission (XG tag, v1.14).
    """
    trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    # Gene-level types only (no mRNA/transcript — those expose transcript IDs
    # like "YAL038W_id001" which break the systematic-name join key used by
    # cross-project taxonomies. SGD GFF reliably emits a 'gene' feature at
    # every protein-coding locus.)
    GENE_TYPES = {"gene", "pseudogene", "ncRNA_gene",
                   "tRNA_gene", "snoRNA_gene", "snRNA_gene", "rRNA_gene"}
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] not in GENE_TYPES:
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            if start < end and strand in ("+", "-"):
                name = _parse_gff_gene_name(parts[8]) or f"{chrom}:{start}-{end}({strand})"
                trees[chrom][start:end] = (strand, name)
    return dict(trees)


def classify_sense_antisense(chrom: str, aln_start: int, aln_end: int, orient: str,
                              gene_trees: Dict[str, IntervalTree],
                              min_gene_frac: float = 0.3,
                              min_read_frac: float = 0.8
                              ) -> Tuple[str, Optional[str]]:
    """v1.14: Overlap-based gene attribution with sense-strand preference and first-TSS tiebreaker.

    For each gene whose body overlaps the cluster's alignment:
      - gene_frac = overlap / gene_length  (how much of the gene is covered)
      - read_frac = overlap / read_aln_length  (how much of the read is in the gene)

    Keep candidates passing `gene_frac >= min_gene_frac` OR `read_frac >= min_read_frac`.
    Among these, classify by sense-strand match: the implied mRNA strand is "+"
    for orient=fwd, "-" for orient=rev.

      - If any candidate gene's strand matches mRNA strand → SENSE. Tie-break
        by picking the gene whose TSS is FIRST in the transcription direction
        (lowest TSS for + reads, highest TSS for - reads). This naturally
        attributes dicistronic readthroughs to the upstream-promoter gene.
      - Else if any candidate gene exists on the opposite strand → ANTISENSE.
        Tie-break by maximum gene_frac.
      - Else → unannotated.

    Replaces the legacy anchor-proximity classifier (pre-v1.14), which mis-
    classified polyA-drift reads as antisense of downstream genes whose CAP
    start happened to lie just past the read's drift cleavage site.
    """
    tree = gene_trees.get(chrom)
    if tree is None:
        return ("unannotated", None)
    aln_len = aln_end - aln_start
    if aln_len <= 0:
        return ("unannotated", None)

    overlaps = list(tree.overlap(aln_start, aln_end))
    if not overlaps:
        return ("unannotated", None)

    mRNA_strand = "+" if orient == "fwd" else "-"

    # (gene_frac, read_frac, tss, gene_strand, gene_name)
    candidates: List[Tuple[float, float, int, str, str]] = []
    for iv in overlaps:
        ov_len = max(0, min(iv.end, aln_end) - max(iv.begin, aln_start))
        if ov_len <= 0:
            continue
        gene_len = iv.end - iv.begin
        gene_frac = ov_len / gene_len
        read_frac = ov_len / aln_len
        if gene_frac < min_gene_frac and read_frac < min_read_frac:
            continue
        gene_strand, gene_name = iv.data
        tss = iv.begin if gene_strand == "+" else iv.end
        candidates.append((gene_frac, read_frac, tss, gene_strand, gene_name))

    if not candidates:
        return ("unannotated", None)

    # Sense-preferred: pick first by transcription direction TSS.
    sense_candidates = [c for c in candidates if c[3] == mRNA_strand]
    if sense_candidates:
        # First TSS in transcription direction = lowest for +, highest for -.
        if mRNA_strand == "+":
            sense_candidates.sort(key=lambda c: c[2])
        else:
            sense_candidates.sort(key=lambda c: -c[2])
        return ("sense", sense_candidates[0][4])

    # Antisense fallback: pick the gene whose body the read overlaps most.
    anti_candidates = [c for c in candidates if c[3] != mRNA_strand]
    if anti_candidates:
        anti_candidates.sort(key=lambda c: -c[0])  # max gene_frac
        return ("antisense", anti_candidates[0][4])

    return ("unannotated", None)
