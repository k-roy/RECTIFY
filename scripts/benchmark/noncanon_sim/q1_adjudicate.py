#!/usr/bin/env python3
"""Q1 adjudication (vetted plan): are the guard's changes at NON-annotated chr5 junctions FIXES or HARMS,
judged base-exact against Illumina split-read support? Reports over the DECIDABLE subset, relative to the
~0.52 HP-abutting positive-control floor.

Must-fixes folded in:
 - PER-JUNCTION (collapse the changed reads to distinct moves), not per-read.
 - Re-partition against COMPREHENSIVE gencode v44 first (basic-novel that is comprehensive-annotated =
   adjudicated free); carry only the truly-novel residual into Illumina adjudication.
 - BASE-EXACT (+/-0) Illumina support; the guard move must cross OUT of the ambiguity window (arm-B != guard
   after normalize) to be resolvable.
 - Support-RATIO margin (guard-pos reads >> arm-B-pos reads) so an Illumina HP-artifact can't fake HARM.
 - HARM = investigation set (one-directional), not a clean error rate.
"""
import sys
from collections import defaultdict
import pysam

RG = "/scratch/users/kevinroy/rectify_guard"
sys.path.insert(0, RG)
from real_drs_hp_drift import load_genome, read_junctions
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.core.consensus.chimeric_consensus import normalize_junction
from rectify.utils.genome import register_genome_contigs_from_fasta

G = "/scratch/users/kevinroy/compass_a549/COMPASS/genome_references/GRCh38_gencode_v44_chr5.fasta"
F_BASIC = "/scratch/users/kevinroy/sumner_lab/references/gencode.v44.basic.chr5.gtf"
F_COMP = "/scratch/users/kevinroy/rectify_human_validation/error_model_gm12878/refs/gencode.v44.annotation.gtf"
ARM_B = "/scratch/users/kevinroy/human_drs_out/real_arm_B.bam"
ARM_G = "/scratch/users/kevinroy/human_drs_out/real_arm_Bguard.bam"
ILL = "/scratch/users/kevinroy/sgnex_a549_illumina/replicate{r}/SGNex_A549_Illumina_replicate{r}_run1.bam"
REPS = (1, 3, 5)
MIN_ANCHOR, K, RATIO = 8, 2, 3.0

register_genome_contigs_from_fasta(G)
g = load_genome(G)
seq5 = next(iter(g.values()))


def annot_norm_set(gtf, chrom_filter="chr5"):
    out = set()
    for t in load_annotated_junctions(gtf):
        c = str(t[0])
        if chrom_filter and c not in (chrom_filter, chrom_filter.replace("chr", "")):
            continue
        out.add(normalize_junction(int(t[1]), int(t[2]), seq5))
    return out


basic = annot_norm_set(F_BASIC)
try:
    comp = annot_norm_set(F_COMP)
except Exception as e:
    print(f"[warn] comprehensive GTF load failed ({type(e).__name__}: {e}); using basic only", flush=True)
    comp = basic
print(f"annotated (basic): {len(basic)}; (comprehensive): {len(comp)}", flush=True)


def extract_sr(bam, chrom):
    b = pysam.AlignmentFile(bam)
    counts = defaultdict(int)
    for r in b.fetch(chrom):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        cig = r.cigartuples or []
        ref = r.reference_start
        for i, (op, ln) in enumerate(cig):
            if op == 3:
                pre = cig[i - 1][1] if i > 0 and cig[i - 1][0] in (0, 7, 8) else 0
                post = cig[i + 1][1] if i + 1 < len(cig) and cig[i + 1][0] in (0, 7, 8) else 0
                if pre >= MIN_ANCHOR and post >= MIN_ANCHOR:
                    counts[(ref, ref + ln)] += 1
            if op in (0, 2, 3, 7, 8):
                ref += ln
    b.close()
    return counts


sr = {r: extract_sr(ILL.format(r=r), "5") for r in REPS}
print("Illumina split-junction dicts:", {r: len(sr[r]) for r in REPS}, flush=True)


def sr_reads(intron):
    """max per-rep base-exact split-read count for this exact (donor,acceptor)."""
    return [sr[r].get(intron, 0) for r in REPS]


def supported(intron):
    return sum(1 for c in sr_reads(intron) if c >= K) >= 2   # base-exact, >=K in >=2/3 reps


def total_sr(intron):
    return sum(sr_reads(intron))


# per-read arm junctions -> distinct MOVES (arm-B intron -> guard intron), DRS-read weighted
B = read_junctions(ARM_B)
Gj = read_junctions(ARM_G)
moves = defaultdict(int)  # (b_intron, g_intron) -> DRS reads
for k in set(B) & set(Gj):
    bset, gset = set(B[k]), set(Gj[k])
    if bset == gset:
        continue
    b_only, g_only = bset - gset, gset - bset
    for bi in b_only:
        if not g_only:
            continue
        gi = min(g_only, key=lambda x: abs(x[0] - bi[0]) + abs(x[1] - bi[1]))
        moves[(bi, gi)] += 1
print(f"distinct guard MOVES (per-junction): {len(moves)}; DRS reads: {sum(moves.values())}", flush=True)


def is_annot(intron, aset):
    return normalize_junction(intron[0], intron[1], seq5) in aset


# classify each distinct move
stats = defaultdict(int)
fixes, harms = [], []
for (bi, gi), nreads in moves.items():
    # re-partition: if guard target is annotated in comprehensive -> adjudicated free (annotated)
    ga, ba = is_annot(gi, comp), is_annot(bi, comp)
    if ga or ba:
        stats["annotated_in_comprehensive"] += 1
        if ga and not ba:
            stats["comp_FIX_guard_to_annotated"] += 1
        elif ba and not ga:
            stats["comp_HARM_guard_off_annotated"] += 1; harms.append((bi, gi, total_sr(gi), total_sr(bi), nreads))
        else:
            stats["comp_NEUTRAL_both_annotated"] += 1
        continue
    stats["truly_novel_residual"] += 1
    # ambiguity gate: must cross OUT of the ambiguity window (normalized bi != gi)
    if normalize_junction(*bi, seq5) == normalize_junction(*gi, seq5):
        stats["ambiguity_equivalent_unresolvable"] += 1
        continue
    stats["resolvable"] += 1
    gs, bs = supported(gi), supported(bi)
    gt, bt = total_sr(gi), total_sr(bi)
    if not (gs or bs):
        stats["inconclusive_no_sr_coverage"] += 1
    elif gs and gt > RATIO * bt:
        stats["FIX_guard_supported"] += 1; fixes.append((bi, gi, gt, bt, nreads))
    elif bs and bt > RATIO * gt:
        stats["HARM_armB_supported"] += 1; harms.append((bi, gi, gt, bt, nreads))
    else:
        stats["inconclusive_ambiguous_support"] += 1

print("\n=== Q1 ADJUDICATION (guard changes at non-annotated chr5 junctions) ===")
for k in ("annotated_in_comprehensive", "comp_FIX_guard_to_annotated", "comp_HARM_guard_off_annotated",
          "comp_NEUTRAL_both_annotated", "truly_novel_residual",
          "ambiguity_equivalent_unresolvable", "resolvable", "inconclusive_no_sr_coverage",
          "inconclusive_ambiguous_support", "FIX_guard_supported", "HARM_armB_supported"):
    print(f"  {k}: {stats[k]}")
dec = stats["FIX_guard_supported"] + stats["HARM_armB_supported"]
print(f"\nDECIDABLE (FIX+HARM): {dec}; FIX:HARM = {stats['FIX_guard_supported']}:{stats['HARM_armB_supported']}")
print("HARM entries are an INVESTIGATION SET (one-directional: a real in-run acceptor the guard rightly holds")
print("presents identically). Sample HARM moves (bi, gi, guard_sr, armB_sr, drs_reads):")
for h in harms[:10]:
    print("   ", h)
print("Sample FIX moves:")
for f in fixes[:6]:
    print("   ", f)
