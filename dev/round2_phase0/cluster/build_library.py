#!/usr/bin/env python
"""Round-2 cDNA Phase-0 — Step C: build the padded-cDNA library + block maps.

For each selected locus emit one cDNA per distinct intron chain:
  - the MANE/canonical exon chain (Ensembl GTF.109), AND
  - any Round-1 read-supported intron chain whose every junction PASSES the k>=10
    perfect-match anchor gate (advisor: MANE-only manufactures the novelty-suppression
    false-GO; the de-novo chains are what the production feature ships).

Each cDNA = [5' genomic pad] + exon1 + ... + exonN + [3' genomic pad] (introns removed),
padded with REAL genomic sequence, intergenic-clipped (50 bp guard before the nearest
other-transcript exon). The CdnaModel block map ties transcript<->genome 1:1; we assert
per-base sequence identity (the real validation) on top of CdnaModel.validate's structural
checks. Outputs: library.fa (aligned to in Round 2) + blockmaps.json (consumed by lift-over).

Run in the rectify env on Sherlock.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
from collections import defaultdict

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from liftover import Block, CdnaModel, revcomp  # noqa: E402

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def comp(b): return b.translate(_COMP)


# ───────────────────────── GTF parsing ─────────────────────────
def parse_gtf_exons(gtf_path, chroms):
    """Return (tx_exons, tx_meta, all_exon_intervals).
    tx_exons[txid] = sorted list of (g_start0, g_end) half-open plus-strand.
    tx_meta[txid]  = (chrom, strand, gene_name).
    all_exon_intervals[chrom] = sorted list of (g_start0, g_end, txid) for neighbor clipping.
    """
    tx_exons = defaultdict(list)
    tx_meta = {}
    all_iv = defaultdict(list)
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            chrom = f[0]
            if chrom not in chroms:
                continue
            g0 = int(f[3]) - 1
            g1 = int(f[4])
            strand = f[6]
            attr = f[8]
            txid = _attr(attr, "transcript_id")
            if txid is None:
                continue
            txid = txid.split(".")[0]   # strip GENCODE version (ENST...8); Ensembl is already unversioned
            gname = _attr(attr, "gene_name") or _attr(attr, "gene_id") or "NA"
            tx_exons[txid].append((g0, g1))
            tx_meta[txid] = (chrom, strand, gname)
            all_iv[chrom].append((g0, g1, txid))
    for t in tx_exons:
        tx_exons[t].sort()
    for c in all_iv:
        all_iv[c].sort()
    return tx_exons, tx_meta, all_iv


def _attr(attr, key):
    k = key + ' "'
    i = attr.find(k)
    if i < 0:
        return None
    i += len(k)
    j = attr.find('"', i)
    return attr[i:j]


# ───────────────────────── pad clipping ─────────────────────────
def clip_pad(chrom, lo, hi, all_iv, exclude_txid, guard=50):
    """Shrink genomic pad interval [lo,hi) so it stops `guard` bp before the nearest
    exon of any OTHER transcript that intrudes into the window. Returns (lo,hi) clipped."""
    for (e0, e1, txid) in all_iv.get(chrom, []):
        if txid == exclude_txid:
            continue
        if e1 <= lo or e0 >= hi:
            continue
        # overlapping other-gene exon -> clip the pad away from it
        # if the exon is on the high side, pull hi down; on the low side, pull lo up
        mid = (lo + hi) // 2
        if e0 >= mid:
            hi = min(hi, max(lo, e0 - guard))
        else:
            lo = max(lo, min(hi, e1 + guard))
    return lo, hi


# ───────────────────────── model construction ─────────────────────────
def build_model(cdna_id, chrom, strand, exons, genome, all_iv, exclude_txid,
                pad=1500, chrom_len=None):
    """exons: sorted ascending list of (g0,g1) half-open plus-strand. Returns (CdnaModel, seq, junctions)."""
    chrom_len = chrom_len if chrom_len is not None else len(genome[chrom])
    # 5'/3' pad in genome plus-strand coords
    lo_exon = exons[0][0]
    hi_exon = exons[-1][1]
    pad_low = clip_pad(chrom, max(0, lo_exon - pad), lo_exon, all_iv, exclude_txid)
    pad_high = clip_pad(chrom, hi_exon, min(chrom_len, hi_exon + pad), all_iv, exclude_txid)

    # transcript-ordered blocks (5'->3'): + ascending genome, - descending genome
    # On '+': low pad, exons asc, high pad.  On '-': high pad, exons desc, low pad.
    ordered = []  # (g0,g1,is_pad,kind)
    if strand == "+":
        if pad_low[1] > pad_low[0]:
            ordered.append((pad_low[0], pad_low[1], True, "pad5"))
        for (g0, g1) in exons:
            ordered.append((g0, g1, False, "exon"))
        if pad_high[1] > pad_high[0]:
            ordered.append((pad_high[0], pad_high[1], True, "pad3"))
    else:
        if pad_high[1] > pad_high[0]:
            ordered.append((pad_high[0], pad_high[1], True, "pad5"))
        for (g0, g1) in reversed(exons):
            ordered.append((g0, g1, False, "exon"))
        if pad_low[1] > pad_low[0]:
            ordered.append((pad_low[0], pad_low[1], True, "pad3"))

    blocks = []
    seq_parts = []
    t = 0
    for (g0, g1, is_pad, kind) in ordered:
        blocks.append(Block(t_start=t, g_start=g0, g_end=g1, is_pad=is_pad, kind=kind))
        gsub = genome[chrom][g0:g1]
        seq_parts.append(revcomp(gsub) if strand == "-" else gsub)
        t += (g1 - g0)
    seq = "".join(seq_parts)
    model = CdnaModel(cdna_id=cdna_id, chrom=chrom, strand=strand, blocks=blocks, length=t)

    # structural validation
    model.validate(genome)
    # per-base sequence identity (the real 1:1 check)
    for i in range(0, len(seq), 1):
        g = model.t2g(i)
        exp = genome[chrom][g] if strand == "+" else comp(genome[chrom][g])
        if seq[i].upper() != exp.upper():
            raise AssertionError(f"{cdna_id}: seq mismatch at t={i} (cdna {seq[i]} vs genome {exp})")

    # junction list (genome donor/acceptor between consecutive exon blocks)
    junctions = []
    exon_blocks = [b for b in blocks if not b.is_pad]
    for a, b in zip(exon_blocks, exon_blocks[1:]):
        if strand == "+":
            junctions.append((a.g_end, b.g_start))   # donor, acceptor (plus coords)
        else:
            junctions.append((b.g_end, a.g_start))
    return model, seq, junctions


# ───────────────────────── de-novo chains from Round-1 ─────────────────────────
def read_chain(read):
    """Return (chrom, strand_unused, tuple of (donor,acceptor) introns) from a spliced read,
    plus per-junction left/right perfect-anchor lengths (in genome=match space)."""
    if read.is_unmapped or read.cigartuples is None:
        return None
    introns = []
    ref = read.reference_start
    # track match-run lengths around each N for the anchor gate (approximate: flanking M length)
    cig = read.cigartuples
    pos = read.reference_start
    # collect (intron_start, intron_end) and flanking aligned-match lengths
    blocks = []  # consecutive ref spans of M/=/X
    cur_start = None
    spans = []
    for op, ln in cig:
        if op in (0, 7, 8):  # M/=/X
            if cur_start is None:
                cur_start = pos
            pos += ln
        elif op == 2:  # D
            pos += ln
        elif op == 3:  # N
            if cur_start is not None:
                spans.append((cur_start, pos))
            introns.append((pos, None))
            cur_start = None
            pos += ln
        # I/S/H don't consume ref
    if cur_start is not None:
        spans.append((cur_start, pos))
    if not introns:
        return None
    # rebuild intron donor/acceptor from spans boundaries
    juncs = []
    for k in range(len(spans) - 1):
        donor = spans[k][1]
        acceptor = spans[k + 1][0]
        left_anchor = spans[k][1] - spans[k][0]
        right_anchor = spans[k + 1][1] - spans[k + 1][0]
        juncs.append((donor, acceptor, left_anchor, right_anchor))
    return read.reference_name, juncs, spans


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--loci", required=True, help="selected_loci.tsv (header row)")
    ap.add_argument("--genes", nargs="+", required=True, help="gene names to build")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pad", type=int, default=1500)
    ap.add_argument("--round1-bams", nargs="*", default=[],
                    help="corrected genome BAMs; mine read-supported gate-passed chains")
    ap.add_argument("--min-anchor", type=int, default=10)
    ap.add_argument("--min-chain-support", type=int, default=3,
                    help="min reads supporting a de-novo intron chain to admit it")
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    chroms = set()
    loci = {}
    with open(args.loci) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if not r or len(r) < len(header):
                continue
            loci[r[ix["gene"]]] = dict(chrom=r[ix["chrom"]], strand=r[ix["strand"]],
                                       txid=r[ix["mane_txid"]])
            chroms.add(r[ix["chrom"]])

    print(f"[lib] loading genome chroms {sorted(chroms)} ...", flush=True)
    fa = pysam.FastaFile(args.genome)
    genome = {c: fa.fetch(c) for c in chroms}
    chrom_lens = {c: len(genome[c]) for c in chroms}

    print(f"[lib] parsing GTF ...", flush=True)
    tx_exons, tx_meta, all_iv = parse_gtf_exons(args.gtf, chroms)

    lib_fa = open(out / "library.fa", "w")
    blockmaps = {}
    summary = []

    for gene in args.genes:
        if gene not in loci:
            print(f"[lib] WARN gene {gene} not in loci tsv; skip"); continue
        info = loci[gene]
        chrom, strand, txid = info["chrom"], info["strand"], info["txid"]
        if txid not in tx_exons:
            print(f"[lib] WARN {gene} MANE {txid} not in GTF; skip"); continue
        exons = tx_exons[txid]
        cdna_id = f"{gene}__MANE_{txid}"
        model, seq, juncs = build_model(cdna_id, chrom, strand, exons, genome, all_iv, txid,
                                        pad=args.pad, chrom_len=chrom_lens[chrom])
        lib_fa.write(f">{cdna_id}\n{seq}\n")
        blockmaps[cdna_id] = _serialize(model, juncs, source="MANE")
        summary.append((cdna_id, len(exons), len(seq), len(juncs)))
        print(f"[lib] {cdna_id}: {len(exons)} exons, {len(juncs)} junc, cdna {len(seq)} bp  OK")

        # de-novo gate-passed chains for this gene
        if args.round1_bams:
            denovo = mine_denovo_chains(gene, chrom, strand, exons, args.round1_bams,
                                        args.min_anchor, args.min_chain_support, tx_exons[txid])
            for ci, (chain_exons, support) in enumerate(denovo):
                did = f"{gene}__denovo{ci}_s{support}"
                try:
                    m2, s2, j2 = build_model(did, chrom, strand, chain_exons, genome, all_iv, txid,
                                             pad=args.pad, chrom_len=chrom_lens[chrom])
                except AssertionError as e:
                    print(f"[lib] denovo {did} failed build: {e}"); continue
                lib_fa.write(f">{did}\n{s2}\n")
                blockmaps[did] = _serialize(m2, j2, source=f"denovo_support{support}")
                summary.append((did, len(chain_exons), len(s2), len(j2)))
                print(f"[lib] {did}: {len(chain_exons)} exons, {len(j2)} junc (support {support})  OK [novel-chain admitted]")

    lib_fa.close()
    with open(out / "blockmaps.json", "w") as fh:
        json.dump(blockmaps, fh, indent=1)
    with open(out / "library_summary.tsv", "w") as fh:
        fh.write("cdna_id\tn_exons\tcdna_len\tn_junctions\n")
        for s in summary:
            fh.write("\t".join(map(str, s)) + "\n")
    print(f"[lib] wrote {len(blockmaps)} cDNAs to {out}/library.fa")


def _serialize(model, juncs, source):
    return dict(
        cdna_id=model.cdna_id, chrom=model.chrom, strand=model.strand, length=model.length,
        source=source,
        blocks=[dict(t_start=b.t_start, g_start=b.g_start, g_end=b.g_end,
                     is_pad=b.is_pad, kind=b.kind) for b in model.blocks],
        junctions=[list(j) for j in juncs],
    )


def mine_denovo_chains(gene, chrom, strand, mane_exons, bams, min_anchor, min_support, mane_exon_list):
    """Collect read-supported intron chains at this gene whose every junction passes the
    k>=min_anchor perfect-anchor gate, are NOT identical to the MANE chain, and have
    >=min_support reads. Returns list of (exon_chain, support). Exon chain derived by
    splicing the MANE exon span at the observed donor/acceptor set."""
    glo = mane_exons[0][0] - 2000
    ghi = mane_exons[-1][1] + 2000
    mane_introns = tuple((mane_exons[k][1], mane_exons[k + 1][0]) for k in range(len(mane_exons) - 1))
    chain_support = defaultdict(int)
    chain_juncs = {}
    for bam in bams:
        af = pysam.AlignmentFile(bam, "rb")
        try:
            it = af.fetch(chrom, max(0, glo), ghi)
        except (ValueError, KeyError):
            continue
        for read in it:
            rc = read_chain(read)
            if rc is None:
                continue
            rname, juncs, spans = rc
            # gate: every junction has both anchors >= min_anchor (approx via flanking match length)
            if any(min(la, ra) < min_anchor for (_, _, la, ra) in juncs):
                continue
            introns = tuple((d, a) for (d, a, _, _) in juncs)
            if introns == mane_introns:
                continue  # identical to MANE
            # require the chain to lie within the padded gene span
            if introns[0][0] < glo or introns[-1][1] > ghi:
                continue
            chain_support[introns] += 1
            chain_juncs[introns] = juncs
    out = []
    for introns, sup in sorted(chain_support.items(), key=lambda kv: -kv[1]):
        if sup < min_support:
            continue
        # build an exon chain spanning the MANE 5'/3' bounds but with the observed introns
        exon_chain = introns_to_exons(introns, mane_exons[0][0], mane_exons[-1][1])
        if exon_chain is None:
            continue
        out.append((exon_chain, sup))
    return out[:5]  # cap per gene


def introns_to_exons(introns, gene_lo, gene_hi):
    """Turn a sorted intron set into an exon chain spanning [gene_lo,gene_hi)."""
    introns = sorted(introns)
    exons = []
    cur = gene_lo
    for (d, a) in introns:
        if d <= cur or a <= d:
            return None
        exons.append((cur, d))
        cur = a
    if gene_hi <= cur:
        return None
    exons.append((cur, gene_hi))
    return exons


if __name__ == "__main__":
    main()
