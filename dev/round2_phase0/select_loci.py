#!/usr/bin/env python
"""Phase-0 Step B — data-driven locus selection (coverage ∩ short-internal-exon structure).

Picks chr5/6 genes whose **MANE Select** transcript has a short internal exon (≤30 bp) that is
(a) actually INCLUDED in the A549 reads (not skipped) and (b) supported by enough spliced reads
to score — the stress cases Phase 0 needs. Replaces the a-priori gene list that came up
uncovered in A549 (`apriori_stress_loci_NOTES.md`).

Run on Sherlock in the rectify env (has pysam):
    python select_loci.py
Outputs <P0>/coverage/selected_loci.tsv (ranked) and prints the top rows.
"""
import gzip
import re
import sys
from collections import defaultdict

import pysam

P0 = "/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0"
GTF = "/oak/stanford/groups/larsms/Projects/split_tap/newvolume/genomes/Homo_sapiens.GRCh38.109.gtf.gz"
BAM = f"{P0}/data/a549_pooled_chr5_6.bam"
OUT = f"{P0}/coverage/selected_loci.tsv"

CHROMS = {"5", "6"}
MAX_SHORT = 30          # internal exon length (bp) to count as a micro/short-exon stressor
MIN_SPLICED = 50        # gene-level spliced-read floor (reads with an N-op over the gene)
MIN_INCL_DEPTH = 10     # the short exon must be INCLUDED with at least this mean depth
MIN_INCL_RATIO = 0.10   # short-exon depth / flank-exon depth (reject near-always-skipped exons)

_tx_re = re.compile(r'transcript_id "([^"]+)"')
_gn_re = re.compile(r'gene_name "([^"]+)"')


def parse_gtf():
    """Return MANE-Select transcripts on chr5/6: txid -> (gene, chrom, strand, [(s,e),...])."""
    exons = defaultdict(list)
    meta = {}
    mane = set()
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in CHROMS or f[2] != "exon":
                continue
            attrs = f[8]
            m = _tx_re.search(attrs)
            if not m:
                continue
            txid = m.group(1)
            exons[txid].append((int(f[3]), int(f[4])))     # 1-based inclusive
            if txid not in meta:
                gn = _gn_re.search(attrs)
                meta[txid] = (gn.group(1) if gn else txid, f[0], f[6])
            if "MANE_Select" in attrs:
                mane.add(txid)
    return {t: (meta[t][0], meta[t][1], meta[t][2], sorted(exons[t])) for t in mane}


def short_exon_candidates(tx):
    """Yield (gene, chrom, strand, txid, n_junctions, (s,e,len)_short, (fs,fe)_flank)."""
    for txid, (gene, chrom, strand, exons) in tx.items():
        if len(exons) < 3:
            continue
        for i in range(1, len(exons) - 1):            # internal exons only
            s, e = exons[i]
            L = e - s + 1
            if L <= MAX_SHORT:
                # constitutive-proxy flank = the longer of the two neighbours
                flank = max(exons[i - 1], exons[i + 1], key=lambda x: x[1] - x[0])
                yield (gene, chrom, strand, txid, len(exons) - 1, (s, e, L), flank)


def mean_depth(bam, chrom, start1, end1):
    """Mean inclusion depth over a 1-based inclusive interval (count_coverage excludes N skips)."""
    cov = bam.count_coverage(chrom, start1 - 1, end1)   # 4 base arrays, 0-based half-open
    n = end1 - start1 + 1
    if n <= 0:
        return 0.0
    return sum(cov[0][k] + cov[1][k] + cov[2][k] + cov[3][k] for k in range(n)) / n


def main():
    tx = parse_gtf()
    cands = list(short_exon_candidates(tx))
    sys.stderr.write(f"MANE tx on chr5/6: {len(tx)};  short-internal-exon candidates: {len(cands)}\n")

    bam = pysam.AlignmentFile(BAM)
    rows = []
    seen_gene = set()
    for gene, chrom, strand, txid, njunc, (s, e, L), (fs, fe) in cands:
        if gene in seen_gene:        # one short exon per gene is enough to select it
            continue
        exons = tx[txid][3]
        gs, ge = exons[0][0], exons[-1][1]
        spliced = total = 0
        for r in bam.fetch(chrom, gs - 1, ge):
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            total += 1
            if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
                spliced += 1
        if spliced < MIN_SPLICED:
            continue
        sd = mean_depth(bam, chrom, s, e)
        fd = mean_depth(bam, chrom, fs, fe)
        incl = (sd / fd) if fd > 0 else 0.0
        if sd < MIN_INCL_DEPTH or incl < MIN_INCL_RATIO:
            continue
        seen_gene.add(gene)
        rows.append((gene, chrom, strand, txid, njunc, s, e, L, gs, ge,
                     total, spliced, round(sd, 1), round(fd, 1), round(incl, 3)))

    rows.sort(key=lambda r: r[11], reverse=True)   # by spliced-read support
    header = ["gene", "chrom", "strand", "mane_txid", "n_junctions",
              "short_exon_start", "short_exon_end", "short_exon_len",
              "gene_start", "gene_end", "gene_reads", "gene_spliced_reads",
              "short_exon_meandepth", "flank_exon_meandepth", "inclusion_ratio"]
    with open(OUT, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
    sys.stderr.write(f"selected {len(rows)} loci -> {OUT}\n\n")
    print("\t".join(header))
    for r in rows[:20]:
        print("\t".join(map(str, r)))


if __name__ == "__main__":
    main()
