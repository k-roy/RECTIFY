#!/usr/bin/env python3
"""Re-derive the COMPASS short-read verdict on the 111 GMAP-only recurrent novels.

The locked adjudication_111.json reported "109 artifacts / 2 corroborated" using
EXACT ambiguity-normalized coordinate matching. That is non-robust: the strict
base-equality normalization (`normalize_junction`) leaves same-intron-length
placements a few bp apart UNMERGED, so COMPASS support sitting at a shifted
coordinate is false-negatived (the recount found 2 "artifacts" with 2959 and 323
short-read reads at a same-length <=5bp-shifted coord). And the original verdict
never checked the splice motif AT the supported coordinate (J1's bulk support is
non-canonical AG..CA, not GT-AG).

This re-derivation, for each of the 111, reports the EVIDENCE and a classification:
  - compass_exact   : COMPASS anchored depth at the EXACT normalized gmap coord.
  - compass_relaxed : total COMPASS depth over same-intron-length placements with
                      |start shift| <= WIN (the ambiguity-tolerant support).
  - dominant supported coordinate + its motif (canonical GT-AG/GC-AG on +/- strand
                      vs non-canonical) + annotated?  -> distinguishes a real
                      canonical novel from a non-canonical (SV/misalignment) signal.
  - alt-SS-of-annotated: shares one boundary (<=3bp) with an annotated junction
                      (different other boundary) -> likely an alt-splice-site of a
                      catalogued junction (e.g. J2 = +2 acceptor off annotated 804),
                      not an independent novel.
CLASS in {ARTIFACT, ALT_SS_OF_ANNOTATED, SUPPORTED_NONCANONICAL,
          SUPPORTED_CANONICAL_NOVEL, ANNOTATED}.

Same-intron-length is REQUIRED for the relaxed match (|dstart|<=WIN alone would
false-match a different-acceptor alt-SS). Run ON Sherlock in the `rectify` env.

Example:
  W=/scratch/users/kevinroy/compass_a549
  OAK=/oak/stanford/groups/larsms/Users/kevinroy/compass_a549_out
  PYTHONPATH=$W/rectify_src python dev/rederive_111.py \
    --bam $OAK/A549_rep1.consensus.bam \
    --genome $W/COMPASS/genome_references/GRCh38_gencode_v44.fasta \
    --gtf $W/COMPASS/genome_references/GRCh38_gencode_v44.gtf \
    --novels /scratch/users/kevinroy/deliverable_b/rectify_src/dev/gmap_only_recurrent_novels_chr5.tsv \
    --out $W/rederive_111.json

Author: Kevin R. Roy (+ Claude)
"""
import argparse, csv, json, sys
import pysam
from rectify.core.splice.junction_scoring import collect_junction_counts_from_bam
from rectify.core.consensus.chimeric_consensus import normalize_junction
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.utils.genome import register_genome_contigs_from_fasta

MIN_ANCHOR = 10
WIN = 5            # same-length start-shift window for the relaxed match
BND = 3           # boundary-sharing tolerance for alt-SS-of-annotated
CHROM = 'chr5'
_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(_COMP)[::-1]


def motif(seq, s, e):
    """Return (donor, acceptor, strand) for the canonical reading, or
    (donor+, acceptor+, None) if non-canonical on both strands.
    + : seq[s:s+2] in {GT,GC} & seq[e-2:e]==AG.
    - : revcomp(seq[e-2:e]) in {GT,GC} & revcomp(seq[s:s+2])==AG."""
    left = seq[s:s+2].upper(); right = seq[e-2:e].upper()
    if left in ('GT', 'GC') and right == 'AG':
        return left, right, '+'
    md, ma = rc(right), rc(left)
    if md in ('GT', 'GC') and ma == 'AG':
        return md, ma, '-'
    return left, right, None   # non-canonical (report the + reading)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True)
    ap.add_argument('--genome', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--novels', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    register_genome_contigs_from_fasta(args.genome)
    fa = pysam.FastaFile(args.genome); seq = fa.fetch(CHROM); fa.close()
    print(f'[rd] {CHROM} {len(seq):,} bp', file=sys.stderr)

    raw = collect_junction_counts_from_bam(args.bam, chrom_filter=CHROM, min_anchor_overhang=MIN_ANCHOR)
    norm_compass = {}
    by_len = {}                         # intron_len -> list[(start, end, count)]
    for (c, s, e), k in raw.items():
        if c != CHROM: continue
        ns, ne = normalize_junction(s, e, seq)
        norm_compass[(ns, ne)] = norm_compass.get((ns, ne), 0) + k
        by_len.setdefault(e - s, []).append((s, e, k))
    print(f'[rd] COMPASS chr5 raw junctions {len(raw)}', file=sys.stderr)

    annot_norm = set(); annot_starts = set(); annot_ends = set()
    for t in load_annotated_junctions(args.gtf):
        if t[0] != CHROM: continue
        ns, ne = normalize_junction(t[1], t[2], seq)
        annot_norm.add((ns, ne)); annot_starts.add(int(t[1])); annot_ends.add(int(t[2]))

    novels = []
    with open(args.novels) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row['chrom'] != CHROM: continue
            novels.append((int(row['start']), int(row['end']), int(row.get('gmap_reads', 0) or 0)))

    results = []
    for s_raw, e_raw, gmap in novels:
        ns, ne = normalize_junction(s_raw, e_raw, seq); L = ne - ns
        compass_exact = norm_compass.get((ns, ne), 0)
        # relaxed: same intron length, |start shift| <= WIN
        cands = [(cs, ce, k) for (cs, ce, k) in by_len.get(L, []) if abs(cs - ns) <= WIN]
        relaxed_total = sum(k for _, _, k in cands)
        # dominant supported placement + its motif/annotation
        dom = max(cands, key=lambda x: x[2]) if cands else None
        dom_rec = None
        if dom:
            cs, ce, k = dom
            d, a, strand = motif(seq, cs, ce)
            dom_rec = {'coord': [cs, ce], 'depth': k, 'motif': f'{d}..{a}',
                       'canonical': strand is not None, 'strand': strand,
                       'annotated': (normalize_junction(cs, ce, seq)) in annot_norm}
        ed, ea, estr = motif(seq, ns, ne)   # motif at the exact gmap coord (sanity: should be canonical)
        # alt-SS of an annotated junction: shares one boundary within BND (not exact-annotated)
        shares_start = any(abs(ns - a) <= BND for a in annot_starts)
        shares_end = any(abs(ne - a) <= BND for a in annot_ends)
        exact_annot = (ns, ne) in annot_norm
        alt_ss = (shares_start ^ shares_end) and not exact_annot

        if exact_annot:
            cls = 'ANNOTATED'
        elif relaxed_total == 0:
            cls = 'ARTIFACT'
        elif alt_ss:
            cls = 'ALT_SS_OF_ANNOTATED'
        elif dom_rec and not dom_rec['canonical']:
            cls = 'SUPPORTED_NONCANONICAL'
        else:
            cls = 'SUPPORTED_CANONICAL_NOVEL'

        results.append({
            'gmap_coord': f'{CHROM}:{ns}-{ne}', 'intron_len': L, 'gmap_reads': gmap,
            'gmap_coord_motif': f'{ed}..{ea}', 'gmap_coord_canonical': estr is not None,
            'compass_exact': compass_exact, 'compass_relaxed_total': relaxed_total,
            'dominant_compass_placement': dom_rec,
            'shares_annotated_boundary': bool(shares_start or shares_end),
            'class': cls,
        })

    # summary
    from collections import Counter
    cls_counts = Counter(r['class'] for r in results)
    summary = {
        'n_novels': len(results), 'window_bp': WIN,
        'EXACT_compass_supported': sum(1 for r in results if r['compass_exact'] > 0),
        'RELAXED_compass_supported': sum(1 for r in results if r['compass_relaxed_total'] > 0),
        'class_counts': dict(cls_counts),
        'note': ("ARTIFACT = no independent short-read support even relaxed. "
                 "SUPPORTED_CANONICAL_NOVEL = the real candidate novel junctions. "
                 "SUPPORTED_NONCANONICAL = supported but at a non-canonical coord (SV/misalignment suspect, e.g. J1). "
                 "ALT_SS_OF_ANNOTATED = shares one boundary with a catalogued junction (e.g. J2)."),
    }
    out = {'summary': summary,
           'junctions': sorted(results, key=lambda r: -r['compass_relaxed_total'])}
    with open(args.out, 'w') as fh:
        json.dump(out, fh, indent=2)
    print(f'[rd] SUMMARY: {json.dumps(summary, indent=2)}', file=sys.stderr)
    print(f'[rd] wrote {args.out}', file=sys.stderr)
    # human table: the supported ones
    print('\n[rd] SUPPORTED junctions (relaxed>0), by class:', file=sys.stderr)
    for r in out['junctions']:
        if r['compass_relaxed_total'] == 0: continue
        d = r['dominant_compass_placement']
        print(f"   {r['class']:26} {r['gmap_coord']} len{r['intron_len']} gmap{r['gmap_reads']} "
              f"exact{r['compass_exact']} relaxed{r['compass_relaxed_total']} "
              f"dom={d['coord'] if d else None} {d['motif'] if d else ''} "
              f"canon={d['canonical'] if d else None} annot={d['annotated'] if d else None}", file=sys.stderr)


if __name__ == '__main__':
    main()
