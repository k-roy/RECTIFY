#!/usr/bin/env python3
"""Q1 positive control (the load-bearing go/no-go from the scope vet): can Illumina resolve
HP-abutting junctions BASE-EXACT? If it can, Q1 (guard fix/harm on novels) is answerable and
this fraction is the decidable floor to report every Q1 rate against. If not, Q1 is unanswerable
for the HP subset — learned cheaply before any effort on the 944.

Reports base-exact split-read concordance at (a) the ~3007 annotated HP-ABUTTING chr5 junctions
and (b) ALL annotated chr5 junctions (the canonical floor for comparison), across the 3 SG-NEx
Illumina reps. Position-only (strand-agnostic): an intron's (donor,acceptor) coords are
strand-independent. Illumina BAMs are Ensembl-named ('5'); annotated junctions are chr5 — single
chrom, so compare (donor,acceptor) integer tuples directly.
"""
import sys
from collections import defaultdict
import pysam

RG = "/scratch/users/kevinroy/rectify_guard"
sys.path.insert(0, RG)
from real_drs_hp_drift import load_genome, hp_abutting
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.utils.genome import register_genome_contigs_from_fasta

G = "/scratch/users/kevinroy/compass_a549/COMPASS/genome_references/GRCh38_gencode_v44_chr5.fasta"
F = "/scratch/users/kevinroy/sumner_lab/references/gencode.v44.basic.chr5.gtf"
ILL = "/scratch/users/kevinroy/sgnex_a549_illumina/replicate{r}/SGNex_A549_Illumina_replicate{r}_run1.bam"
REPS = (1, 3, 5)
MIN_ANCHOR = 8

register_genome_contigs_from_fasta(G)
g = load_genome(G)
annot = list(load_annotated_junctions(F))
annot_coords = set((int(t[1]), int(t[2])) for t in annot)                 # all chr5 annotated
hp = hp_abutting(g, annot, 4, 3)
hp_coords = set((s, e) for (c, s, e) in hp)                                # HP-abutting subset
print(f"annotated chr5 junctions: {len(annot_coords)}; HP-abutting: {len(hp_coords)}", flush=True)


def extract_sr_junctions(bam, chrom):
    """(donor,acceptor) -> read count, for split reads with >=MIN_ANCHOR flanking match on both sides."""
    b = pysam.AlignmentFile(bam)
    counts = defaultdict(int)
    n = 0
    for r in b.fetch(chrom):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        n += 1
        cig = r.cigartuples or []
        ref = r.reference_start
        for i, (op, ln) in enumerate(cig):
            if op == 3:  # N (intron)
                pre = cig[i - 1][1] if i > 0 and cig[i - 1][0] in (0, 7, 8) else 0
                post = cig[i + 1][1] if i + 1 < len(cig) and cig[i + 1][0] in (0, 7, 8) else 0
                if pre >= MIN_ANCHOR and post >= MIN_ANCHOR:
                    counts[(ref, ref + ln)] += 1
            if op in (0, 2, 3, 7, 8):
                ref += ln
    b.close()
    return counts, n


sr = {}
for r in REPS:
    bam = ILL.format(r=r)
    sr[r], nreads = extract_sr_junctions(bam, "5")   # Illumina chrom name is '5'
    print(f"rep{r}: {nreads} chr5 primary reads; {len(sr[r])} distinct anchored split junctions", flush=True)


def concordance(target, label):
    print(f"\n== {label} (n={len(target)}) base-exact split-read concordance ==")
    for K in (1, 2, 3):
        for min_reps in (1, 2):
            resolved = sum(1 for j in target
                           if sum(1 for r in REPS if sr[r].get(j, 0) >= K) >= min_reps)
            print(f"  K>={K}, >={min_reps}/3 reps: {resolved}/{len(target)} = {resolved/len(target):.3f}")


concordance(annot_coords, "ALL annotated chr5 (canonical floor)")
concordance(hp_coords, "HP-abutting annotated chr5 (Q1's decidable floor)")
print("\nGO/NO-GO: if HP-abutting base-exact concordance is non-trivial (say >0.3 at K>=1, >=2/3), "
      "Illumina can adjudicate HP junctions base-exact -> Q1 answerable, report Q1 vs this floor. "
      "If ~0 while the canonical floor is high -> Illumina can't resolve HP -> Q1 needs the sensitive truth.")
