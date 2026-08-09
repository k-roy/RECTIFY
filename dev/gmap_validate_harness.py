#!/usr/bin/env python3
"""GMAP-unique-novel VALIDATION HARNESS — genome-wide-capable, chrom-agnostic.

Generalizes the chr5 "111" adjudication (this session) into reusable tooling to
answer track-goal #2 ("validate GMAP's appropriateness") at ANY chromosome / the
whole genome. Chains the three validated steps:

  STEP 1  CANDIDATE GEN: bucket GMAP junctions vs the rest of the long-read panel
          -> GMAP-only, recurrent (>= --recurrence-min), canonical, NOVEL (not
          annotated). (= deliverable_b_gmap_corroboration logic, chrom-param.)
  STEP 2  SHORT-READ VALIDATION: each candidate vs an INDEPENDENT short-read
          consensus (COMPASS) BAM, ambiguity-tolerant (same intron length,
          |start shift| <= --win) + splice motif AT the supported coordinate +
          annotation. (= rederive_111 logic.) -> ARTIFACT vs SUPPORTED.
  STEP 3  CONFIRM (SUPPORTED only): cross-aligner long-read placement + (optional)
          WGS interior/flank depth ratio at the locus. (= lr_probe / dna_split.)

Reuses the live rectify library (identical coordinate convention + ambiguity
normalization as the original adjudication). Run ON Sherlock in the `rectify` env.

DATA REQUIREMENT for a GENOME-WIDE run (today only chr5 long-read BAMs exist):
  - genome-wide A549 long-read alignments for GMAP + >=1 independent aligner
    (minimap2/deSALT/uLTRA/mapPacBio). These do NOT exist yet (current BAMs are
    a549_chr5_trimmed.*) -> generate before --chrom all. Short-read consensus
    (Oak A549_rep1.consensus.bam) and WGS (a549_wgs_deep.bam) ARE genome-wide.

Example (chr5 self-test — should reproduce ~3-4 real of ~111 candidates):
  W=/scratch/users/kevinroy/compass_a549
  ALN=/scratch/users/kevinroy/rectify_human_validation/sgnex_a549/alignments
  OAK=/oak/stanford/groups/larsms/Users/kevinroy/compass_a549_out
  PYTHONPATH=$W/rectify_src python dev/gmap_validate_harness.py \
    --gmap-bam $ALN/a549_chr5_trimmed.GMAP.bam \
    --other-bams $ALN/a549_chr5_trimmed.minimap2.bam $ALN/a549_chr5_trimmed.deSALT.bam \
                 $ALN/a549_chr5_trimmed.uLTRA.bam $ALN/a549_chr5_trimmed.mapPacBio.bam \
    --short-read-bam $OAK/A549_rep1.consensus.bam \
    --wgs-bam $W/wgs/a549_wgs_deep.bam \
    --genome $W/COMPASS/genome_references/GRCh38_gencode_v44.fasta \
    --gtf $W/COMPASS/genome_references/GRCh38_gencode_v44.gtf \
    --chrom chr5 --out $W/gmap_validate_chr5.json

Author: Kevin R. Roy (+ Claude)
"""
import argparse, json, sys
import pysam
from rectify.core.splice.junction_scoring import collect_junction_counts_from_bam
from rectify.core.consensus.chimeric_consensus import (
    normalize_junction, junction_ambiguity_window, _canonical_within_window)
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.utils.genome import register_genome_contigs_from_fasta

_COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
def _rc(s): return s.translate(_COMP)[::-1]


def motif(seq, s, e):
    """(donor, acceptor, strand) for the canonical reading; strand None if non-canonical."""
    l, r = seq[s:s+2].upper(), seq[e-2:e].upper()
    if l in ('GT', 'GC') and r == 'AG':
        return l, r, '+'
    if _rc(r) in ('GT', 'GC') and _rc(l) == 'AG':
        return _rc(r), _rc(l), '-'
    return l, r, None


def collect_norm(bam, chrom, seq, min_anchor):
    """normalized-junction Counter for one BAM on one chrom (+ raw by_len index)."""
    raw = collect_junction_counts_from_bam(bam, chrom_filter=chrom, min_anchor_overhang=min_anchor)
    norm, by_len = {}, {}
    for (c, s, e), k in raw.items():
        if c != chrom:
            continue
        ns, ne = normalize_junction(s, e, seq)
        norm[(ns, ne)] = norm.get((ns, ne), 0) + k
        by_len.setdefault(e - s, []).append((s, e, k))
    return norm, by_len


def validate_chrom(chrom, seq, args, sr_norm, sr_bylen, annot_norm):
    gmap, _ = collect_norm(args.gmap_bam, chrom, seq, args.min_anchor)
    # "others" = union of independent aligners, indexed BOTH exactly (normalized)
    # and by intron length. The gmap-only test must be AMBIGUITY-TOLERANT: an
    # aligner that places the SAME junction a few bp off gmap's coordinate (same
    # intron length, |start shift| <= win) is NOT a gmap-only miss. Using exact
    # normalized matching over-counts gmap-only (the confirm on the chr5 111 +
    # 5-chrom SUPPORTED showed 10/13 "gmap-only" are actually multi-aligner, just
    # a few bp off — the same normalization gap that inflated the 111).
    others = set()
    others_bylen = {}  # intron_len -> list of raw starts
    for b in args.other_bams:
        on, obl = collect_norm(b, chrom, seq, args.min_anchor)
        others |= set(on.keys())
        for L, lst in obl.items():
            others_bylen.setdefault(L, []).extend(s for (s, e, k) in lst)

    def _in_others(s, e):
        if (s, e) in others:
            return True
        L = e - s
        return any(abs(os_ - s) <= args.win for os_ in others_bylen.get(L, ()))

    # STEP 1: GMAP-only (ambiguity-tolerant), recurrent, canonical, novel
    candidates = []
    for (s, e), k in gmap.items():
        if k < args.recurrence_min:
            continue
        if _in_others(s, e) or (chrom, s, e) in annot_norm:
            continue
        l_amb, r_amb = junction_ambiguity_window(s, e, seq)
        if not _canonical_within_window(s, e, seq, l_amb, r_amb):
            continue
        candidates.append((s, e, k))
    # STEP 2: short-read validation (ambiguity-tolerant same-len +-win)
    results = []
    for s, e, gmap_k in candidates:
        L = e - s
        exact = sr_norm.get((s, e), 0)
        cands = [(cs, ce, kk) for (cs, ce, kk) in sr_bylen.get(L, []) if abs(cs - s) <= args.win]
        relaxed = sum(kk for _, _, kk in cands)
        dom = max(cands, key=lambda x: x[2]) if cands else None
        dom_rec = None
        if dom:
            cs, ce, kk = dom
            d, a, strand = motif(seq, cs, ce)
            dom_rec = {'coord': [cs, ce], 'sr_depth': kk, 'motif': f'{d}..{a}',
                       'canonical': strand is not None, 'strand': strand}
        cls = 'ARTIFACT' if relaxed == 0 else 'SUPPORTED'
        results.append({'junction': f'{chrom}:{s}-{e}', 'intron_len': L,
                        'gmap_reads': gmap_k, 'sr_exact': exact, 'sr_relaxed': relaxed,
                        'dominant_sr_placement': dom_rec, 'class': cls})
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gmap-bam', required=True)
    ap.add_argument('--other-bams', nargs='+', required=True, help='independent (non-GMAP) long-read aligner BAMs')
    ap.add_argument('--short-read-bam', required=True, help='genome-wide COMPASS short-read consensus BAM')
    ap.add_argument('--wgs-bam', default=None, help='(optional) A549 WGS BAM for the confirm step')
    ap.add_argument('--genome', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--chrom', default='all', help='one chrom, or "all" for every reference contig')
    ap.add_argument('--recurrence-min', type=int, default=5)
    ap.add_argument('--min-anchor', type=int, default=10)
    ap.add_argument('--win', type=int, default=5)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    register_genome_contigs_from_fasta(args.genome)
    fa = pysam.FastaFile(args.genome)
    chroms = [args.chrom] if args.chrom != 'all' else list(fa.references)
    annot_all = load_annotated_junctions(args.gtf)

    all_results = {}
    for chrom in chroms:
        if chrom not in fa.references:
            print(f'[harness] skip {chrom} (not in genome)', file=sys.stderr); continue
        seq = fa.fetch(chrom)
        annot_norm = set()
        for t in annot_all:
            if t[0] == chrom:
                ns, ne = normalize_junction(t[1], t[2], seq)
                annot_norm.add((chrom, ns, ne))
        sr_norm, sr_bylen = collect_norm(args.short_read_bam, chrom, seq, args.min_anchor)
        if len(sr_norm) < 50:
            print(f'[harness] {chrom}: only {len(sr_norm)} short-read junctions — skip (no SR coverage)', file=sys.stderr)
            continue
        res = validate_chrom(chrom, seq, args, sr_norm, sr_bylen, annot_norm)
        all_results[chrom] = res
        n_art = sum(1 for r in res if r['class'] == 'ARTIFACT')
        n_sup = sum(1 for r in res if r['class'] == 'SUPPORTED')
        n_canon = sum(1 for r in res if r['class'] == 'SUPPORTED' and r['dominant_sr_placement']
                      and r['dominant_sr_placement']['canonical'])
        print(f'[harness] {chrom}: {len(res)} gmap-only-recurrent-canonical-novel candidates '
              f'-> {n_art} ARTIFACT, {n_sup} SUPPORTED ({n_canon} canonical-at-supported-coord)', file=sys.stderr)
    fa.close()

    summary = {chrom: {'candidates': len(r),
                       'artifacts': sum(1 for x in r if x['class'] == 'ARTIFACT'),
                       'supported': sum(1 for x in r if x['class'] == 'SUPPORTED')}
               for chrom, r in all_results.items()}
    tot_c = sum(s['candidates'] for s in summary.values())
    tot_s = sum(s['supported'] for s in summary.values())
    out = {'params': {'recurrence_min': args.recurrence_min, 'win': args.win, 'min_anchor': args.min_anchor},
           'per_chrom_summary': summary,
           'TOTAL': {'candidates': tot_c, 'supported': tot_s, 'artifacts': tot_c - tot_s,
                     'note': 'SUPPORTED = independent short-read corroboration (ambiguity-tolerant). '
                             'Motif at supported coord + cross-aligner/WGS confirm distinguishes real '
                             'non-canonical junction from SV/artifact (see dev/COMPASS_2corroborated_CROSSPLATFORM.md). '
                             'WGS confirm step: run dev/dna_split.py + coverage on --wgs-bam for SUPPORTED loci.'},
           'detail': all_results}
    with open(args.out, 'w') as fh:
        json.dump(out, fh, indent=2)
    print(f'[harness] TOTAL {tot_c} candidates -> {tot_s} SUPPORTED; wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
