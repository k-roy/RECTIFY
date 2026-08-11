#!/usr/bin/env python3
"""Genome-wide RECTIFY discovery on ONE Sumner WGS DRS sample (downsampled first pass).
Downsample -> refine (motif-blind + guard, FULL GRCh38 so chrom names register verbatim) ->
classify raw vs refined junctions (annotated/novel x canonical/non-canonical) -> per-sample TSV.

Discovery = junctions the motif-blind re-aligner places at NON-CANONICAL / UNANNOTATED positions
that raw minimap2 flattened (canonical/annotated). Reported genome-wide, per sample; aggregate SMA vs WT after.
Usage: python sumner_gw_discover.py --bam <wgs.bam> --sample <name> --frac 0.05 --genome <GRCh38.fa> --gff <gencode.gtf> --outdir <dir>
"""
import argparse, os, subprocess, sys
from collections import defaultdict
import pysam

RG = "/scratch/users/kevinroy/rectify_guard"
sys.path.insert(0, RG)
from real_drs_hp_drift import load_genome, read_junctions
from rectify.core.splice.junction_refiner import refine_bam_junctions, build_junction_pool
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.core.consensus.chimeric_consensus import normalize_junction
from rectify.utils.genome import register_genome_contigs_from_fasta, standardize_chrom_name


def canonical(chrom, donor, acceptor, g):
    """canonical if GT..AG (+ strand) or CT..AC (- strand) at the intron ends."""
    seq = g.get(chrom, "")
    if not seq or acceptor > len(seq) or donor < 0:
        return None
    d = seq[donor:donor + 2]; a = seq[acceptor - 2:acceptor]
    return (d == "GT" and a == "AG") or (d == "CT" and a == "AC")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True); ap.add_argument("--sample", required=True)
    ap.add_argument("--frac", type=float, default=0.05)
    ap.add_argument("--genome", required=True); ap.add_argument("--gff", required=True)
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    S = a.sample

    # 1. downsample (primary only), keep coord-sorted, index
    sub = os.path.join(a.outdir, f"{S}.sub.bam")
    subprocess.run(f"samtools view -s {a.frac} -b -F 0x900 {a.bam} -o {sub} && samtools index {sub}",
                   shell=True, check=True)
    nsub = int(subprocess.run(f"samtools view -c {sub}", shell=True, capture_output=True, text=True).stdout or 0)
    sys.stderr.write(f"[{S}] downsampled to {nsub} reads (frac {a.frac})\n"); sys.stderr.flush()

    # 2. genome (FULL -> chrom names register verbatim, no chr10->chrX collision) + comprehensive annot
    register_genome_contigs_from_fasta(a.genome)
    g = load_genome(a.genome)
    annot = list(load_annotated_junctions(a.gff))
    annot_norm = {(standardize_chrom_name(str(t[0])),) + normalize_junction(int(t[1]), int(t[2]),
                  g.get(standardize_chrom_name(str(t[0])), "")) for t in annot
                  if g.get(standardize_chrom_name(str(t[0])), "")}
    sys.stderr.write(f"[{S}] genome {len(g)} contigs; {len(annot_norm)} annotated junctions\n"); sys.stderr.flush()

    # 3. refine (motif-blind + guard)
    pool, aset = build_junction_pool([sub], annot)
    ref = os.path.join(a.outdir, f"{S}.refined.bam")
    refine_bam_junctions(input_bam=sub, output_bam=ref, aligner_bams=[sub], annotated_junctions=annot,
                         genome=g, penalty_table_path=None, prebuilt_junction_pool=pool,
                         prebuilt_annotated_set=aset, sort_and_index=True, n_workers=a.workers,
                         motif_blind=True, hp_drift_margin=3.0)

    # 4. classify raw vs refined junctions
    def classify(bam):
        rj = read_junctions(bam)
        stats = defaultdict(int); noncanon_novel = set()
        for (c, _rid), introns in rj.items():
            for (s, e) in introns:
                n = (c,) + normalize_junction(s, e, g.get(c, ""))
                ann = n in annot_norm
                can = canonical(c, s, e, g)
                key = ("annot" if ann else "novel") + ("_canon" if can else "_noncanon")
                stats[key] += 1
                if (not ann) and (can is False):
                    noncanon_novel.add((c, s, e))
        return stats, noncanon_novel

    raw_stats, raw_nn = classify(sub)
    ref_stats, ref_nn = classify(ref)
    revealed = ref_nn - raw_nn   # novel non-canonical the refine REVEALED (not in raw)

    # 5. write per-sample summary + revealed junctions
    with open(os.path.join(a.outdir, f"{S}.summary.tsv"), "w") as fh:
        fh.write(f"sample\t{S}\nsub_reads\t{nsub}\nfrac\t{a.frac}\n")
        for k in ("annot_canon", "annot_noncanon", "novel_canon", "novel_noncanon"):
            fh.write(f"raw_{k}\t{raw_stats.get(k,0)}\tref_{k}\t{ref_stats.get(k,0)}\n")
        fh.write(f"revealed_novel_noncanon\t{len(revealed)}\n")
    with open(os.path.join(a.outdir, f"{S}.revealed_noncanon.tsv"), "w") as fh:
        for (c, s, e) in sorted(revealed):
            fh.write(f"{c}\t{s}\t{e}\t{S}\n")
    print(f"[{S}] raw novel_noncanon {raw_stats.get('novel_noncanon',0)} -> refined {ref_stats.get('novel_noncanon',0)}; "
          f"REVEALED by refine: {len(revealed)}")


if __name__ == "__main__":
    main()
