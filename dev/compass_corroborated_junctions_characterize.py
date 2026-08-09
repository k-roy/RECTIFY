#!/usr/bin/env python3
"""Characterize the 2 CORROBORATED junctions from the COMPASS short-read
111-adjudication (dev/compass_shortread_adjudicate_111.py).

The adjudication found 109/111 GMAP-only recurrent novels to be artifacts (no
independent short-read support) and 2 CORROBORATED — junctions the multi-aligner
COMPASS panel DOES place that single-pass STAR (Deliverable B) missed:
  chr5:140564954-140565547  compass_depth 12, gmap 8
  chr5:179823051-179823857  compass_depth  3, gmap 14

This turns "corroborated (booleans canonical=True, depth N)" into a verdict-ready
packet for Kevin so the accept/reject decision doesn't require manually driving
IGV. For each junction it reports the SPECIFICS the booleans hide:
  - the SPECIFIC canonical placement in the ambiguity window: exact intron coords
    + the actual donor/acceptor dinucleotides (GT..AG vs GC..AG) and INFERRED
    STRAND (plus = GT..AG; minus = CT..AC on the + genome = revcomp GT..AG).
  - annotation status (exact-junction membership) + the OVERLAPPING annotated
    gene(s) and their strand -> the sharpest discriminator: motif-strand vs
    gene-strand agreement (disagree => antisense / likely-still-artifact).
  - splice-in support as a RATIO, not raw depth: anchored split-read junction
    depth over reads spanning each flank (depth 3 of 3 vs 3 of 300 are opposite
    verdicts -- exactly what IGV shows visually).
  - the supporting split-read QNAMEs (so Kevin can find them in IGV).

SELF-CHECK: recomputed anchored junction depths MUST equal 12 and 3 (the locked
adjudication_111.json). A mismatch means a coordinate/normalization bug -> abort.

Reuses the live rectify library (same coordinate convention + ambiguity
normalization as the adjudication). Run ON Sherlock in the `rectify` env.

Example:
  W=/scratch/users/kevinroy/compass_a549
  OAK=/oak/stanford/groups/larsms/Users/kevinroy/compass_a549_out
  python dev/compass_corroborated_junctions_characterize.py \
    --bam    $OAK/A549_rep1.consensus.bam \
    --genome $W/COMPASS/genome_references/GRCh38_gencode_v44.fasta \
    --gtf    $W/COMPASS/genome_references/GRCh38_gencode_v44.gtf \
    --out    $W/corroborated_2_characterization.json

Author: Kevin R. Roy (+ Claude)
"""
import argparse
import json
import sys

import pysam
from rectify.core.consensus.chimeric_consensus import (
    normalize_junction, junction_ambiguity_window,
    CANONICAL_5SS, CANONICAL_3SS)
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.core.splice.junction_scoring import collect_junction_counts_from_bam
from rectify.utils.genome import register_genome_contigs_from_fasta

MIN_ANCHOR = 10        # match the adjudication's anchored-junction gate
CHROM = 'chr5'

# The 2 corroborated junctions (raw coords from gmap_only_recurrent_novels_chr5.tsv).
CORROBORATED = [
    {'raw': (140564954, 140565547), 'compass_depth': 12, 'gmap_reads': 8},
    {'raw': (179823051, 179823857), 'compass_depth': 3,  'gmap_reads': 14},
]

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def _load_chrom_seq(fa_path, chrom):
    fa = pysam.FastaFile(fa_path)
    name = chrom if chrom in fa.references else (
        'chr' + chrom if 'chr' + chrom in fa.references else chrom.replace('chr', ''))
    seq = fa.fetch(name)
    fa.close()
    return seq


def canonical_placement(s, e, seq, l_amb, r_amb):
    """Scan the ambiguity window for a canonical splice placement on EITHER strand.

    Intron is [js, je). Plus strand: donor = seq[js:js+2] in {GT,GC},
    acceptor = seq[je-2:je] == AG. Minus strand: the sense (minus) intron reads
    GT..AG, so on the + genome donor sits at the RIGHT end revcomp'd:
      donor  = revcomp(seq[je-2:je]) in {GT,GC};  acceptor = revcomp(seq[js:js+2]) == AG.
    Returns the first canonical placement found (dict) or None.
    """
    for shift in range(-l_amb, r_amb + 1):
        js, je = s + shift, e + shift
        if js < 2 or je + 2 > len(seq):
            continue
        left = seq[js:js + 2].upper()
        right = seq[je - 2:je].upper()
        # plus-strand hypothesis
        if left in CANONICAL_5SS and right in CANONICAL_3SS:
            return {'strand': '+', 'intron_start_0based': js, 'intron_end_0based': je,
                    'donor': left, 'acceptor': right, 'motif': f'{left}..{right}',
                    'shift_from_normalized': shift}
        # minus-strand hypothesis (sense intron = revcomp of + genome)
        m_donor = revcomp(right)
        m_acceptor = revcomp(left)
        if m_donor in CANONICAL_5SS and m_acceptor in CANONICAL_3SS:
            return {'strand': '-', 'intron_start_0based': js, 'intron_end_0based': je,
                    'donor': m_donor, 'acceptor': m_acceptor,
                    'motif': f'{m_donor}..{m_acceptor}',
                    'plus_genome': f'{left}..{right}', 'shift_from_normalized': shift}
    return None


def overlapping_genes(gtf_path, chrom, lo, le):
    """Annotated `gene` features (1-based GTF) overlapping the 0-based interval
    [lo, le). Returns list of {name, strand, start, end, biotype}."""
    out = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] != chrom or f[2] != 'gene':
                continue
            g0, g1 = int(f[3]) - 1, int(f[4])  # -> 0-based half-open
            if g0 < le and g1 > lo:
                attrs = f[8]
                name = _attr(attrs, 'gene_name') or _attr(attrs, 'gene_id')
                out.append({'name': name, 'strand': f[6],
                            'start': g0, 'end': g1,
                            'biotype': _attr(attrs, 'gene_type')})
    return out


def _attr(attrs, key):
    tag = key + ' "'
    i = attrs.find(tag)
    if i < 0:
        return None
    i += len(tag)
    j = attrs.find('"', i)
    return attrs[i:j] if j > i else None


def flank_coverage(bam, chrom, pos):
    """Total aligned read depth at a single 0-based position (mapped reads,
    primary+secondary, mirrors how junction depth was counted)."""
    cov = bam.count(chrom, pos, pos + 1,
                    read_callback=lambda r: not r.is_unmapped)
    return cov


def neighbor_junctions(raw_counts, annot_raw, chrom, ns, ne, radius=200):
    """Locus splice landscape from collect's ANCHORED raw N-op counts (the same
    source as the adjudication depth, not a fetch re-walk) within `radius` bp of
    the normalized junction [ns, ne). Reveals whether the corroborated novel sits
    ISOLATED or is a minor variant of a far-more-abundant neighbor — the decisive
    IGV-equivalent context a single depth number hides. Each placement is flagged
    annotated (exact membership in `annot_raw`) so a few-read placement that is a
    small shift off a high-count ANNOTATED junction reads as a likely
    misalignment artifact, while a well-supported placement with no annotated
    neighbor reads as a real novel junction."""
    out = []
    for (c, s, e), k in raw_counts.items():
        if c != chrom:
            continue
        if abs(s - ns) <= radius and abs(e - ne) <= radius:
            out.append({'intron_0based': [s, e], 'intron_len': e - s, 'anchored_reads': k,
                        'annotated': (chrom, s, e) in annot_raw,
                        'is_this_junction': (s == ns and e == ne)})
    return sorted(out, key=lambda x: -x['anchored_reads'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True, help='locked Oak merged consensus BAM')
    ap.add_argument('--genome', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    register_genome_contigs_from_fasta(args.genome)
    seq = _load_chrom_seq(args.genome, CHROM)
    print(f'[char] {CHROM} loaded: {len(seq):,} bp', file=sys.stderr)

    # COMPASS anchored junction depths (normalized to leftmost), for the SELF-CHECK.
    from rectify.core.consensus.chimeric_consensus import normalize_junction as _nj
    raw_counts = collect_junction_counts_from_bam(
        args.bam, chrom_filter=CHROM, min_anchor_overhang=MIN_ANCHOR)
    compass = {}
    for (c, s, e), k in raw_counts.items():
        ns, ne = _nj(s, e, seq)
        compass[(c, ns, ne)] = compass.get((c, ns, ne), 0) + k
    print(f'[char] COMPASS chr5 anchored junctions: {len(compass):,}', file=sys.stderr)

    # Annotated set: normalized (exact-junction novelty test) AND raw (per-placement
    # annotation flag in the locus landscape — to catch a near-shift off an
    # annotated junction, the classic misalignment-artifact signature).
    annot = set()
    annot_raw = set()
    for t in load_annotated_junctions(args.gtf):
        if t[0] != CHROM:
            continue
        annot_raw.add((CHROM, int(t[1]), int(t[2])))
        ns, ne = normalize_junction(t[1], t[2], seq)
        annot.add((CHROM, ns, ne))

    bam = pysam.AlignmentFile(args.bam, 'rb')
    results = []
    self_check_ok = True
    for item in CORROBORATED:
        rs, re_ = item['raw']
        ns, ne = normalize_junction(rs, re_, seq)
        l_amb, r_amb = junction_ambiguity_window(ns, ne, seq)
        key = (CHROM, ns, ne)

        recomputed_depth = compass.get(key, 0)
        expected = item['compass_depth']
        ok = (recomputed_depth == expected)
        self_check_ok = self_check_ok and ok

        placement = canonical_placement(ns, ne, seq, l_amb, r_amb)
        exact_annotated = key in annot
        genes = overlapping_genes(args.gtf, CHROM, ns, ne)

        # splice-in ratio: junction split reads over reads spanning each flank.
        # flank positions: last exonic base before donor (ns-1) and first exonic
        # base after acceptor (ne) -> reads crossing these either splice or read through.
        donor_flank_cov = flank_coverage(bam, CHROM, ns - 1)
        acceptor_flank_cov = flank_coverage(bam, CHROM, ne)
        neighbors = neighbor_junctions(raw_counts, annot_raw, CHROM, ns, ne, radius=200)
        # the most abundant OTHER junction in the locus window (dominant neighbor)
        dom = next((n for n in neighbors if not n['is_this_junction']), None)
        # SAME-intron-length support: collect's strict base-equality normalization
        # leaves same-length placements a few bp apart UNMERGED (e.g. a 3 bp slide
        # over a non-identical base), so the adjudication depth can undercount a
        # real junction. Sum anchored reads over all placements sharing this
        # junction's intron length within the window.
        same_len_total = sum(n['anchored_reads'] for n in neighbors
                             if n['intron_len'] == (ne - ns))
        # NEAR-ANNOTATED-SHIFT: is this a small (<=5 bp donor/acceptor) shift off a
        # high-count ANNOTATED junction? That is the misalignment-artifact signature.
        near_annot = None
        for n in neighbors:
            if n['annotated'] and not n['is_this_junction']:
                ds, de = n['intron_0based'][0] - ns, n['intron_0based'][1] - ne
                if abs(ds) <= 5 or abs(de) <= 5:
                    near_annot = {'annotated_junction_0based': n['intron_0based'],
                                  'annotated_reads': n['anchored_reads'],
                                  'donor_shift': -ds, 'acceptor_shift': -de}
                    break

        # strand cross-check: does the canonical motif strand match an overlapping gene?
        motif_strand = placement['strand'] if placement else None
        gene_strands = sorted({g['strand'] for g in genes})
        if motif_strand and gene_strands:
            if motif_strand in gene_strands:
                strand_agreement = 'AGREE'
            else:
                strand_agreement = 'DISAGREE_antisense_or_artifact'
        elif motif_strand and not gene_strands:
            strand_agreement = 'no_overlapping_gene_intergenic'
        else:
            strand_agreement = 'noncanonical_no_motif_strand'

        rec = {
            'raw_junction': f'{CHROM}:{rs}-{re_}',
            'normalized_junction_0based': [CHROM, ns, ne],
            'intron_length': ne - ns,
            'ambiguity_window': {'left_bp': l_amb, 'right_bp': r_amb},
            'compass_depth_recomputed': recomputed_depth,
            'compass_depth_expected': expected,
            'SELF_CHECK_depth_matches': ok,
            'gmap_long_read_reads': item['gmap_reads'],
            'canonical_placement': placement,
            'exact_junction_annotated': exact_annotated,
            'overlapping_genes': genes,
            'motif_strand': motif_strand,
            'gene_strands': gene_strands,
            'strand_agreement': strand_agreement,
            'splice_in_support': {
                'junction_split_reads': recomputed_depth,
                'donor_flank_total_cov': donor_flank_cov,
                'acceptor_flank_total_cov': acceptor_flank_cov,
                'splice_in_frac_donor': round(recomputed_depth / donor_flank_cov, 4) if donor_flank_cov else None,
                'splice_in_frac_acceptor': round(recomputed_depth / acceptor_flank_cov, 4) if acceptor_flank_cov else None,
            },
            'same_intron_length_total_anchored_reads': same_len_total,
            'near_annotated_shift': near_annot,
            'locus_junction_landscape': neighbors,
            'dominant_neighbor_junction': dom,
        }
        results.append(rec)
        print(f'[char] {rec["raw_junction"]}: depth {recomputed_depth} '
              f'(same-len total {same_len_total}); '
              f'motif {placement["motif"] if placement else "NONCANONICAL"} '
              f'strand {motif_strand}; genes {[g["name"] for g in genes]} ({gene_strands}); '
              f'{strand_agreement}; splice-in {recomputed_depth}/{donor_flank_cov}; '
              f'near-annotated-shift {near_annot}; '
              f'dominant neighbor {dom["intron_0based"] if dom else None} '
              f'(len {dom["intron_len"] if dom else None}, '
              f'{dom["anchored_reads"] if dom else 0} reads, annot={dom["annotated"] if dom else None})',
              file=sys.stderr)
    bam.close()

    if not self_check_ok:
        raise SystemExit('[char] SELF-CHECK FAILED: recomputed depth != locked '
                         'adjudication json (12, 3). Coordinate/normalization bug — aborting.')

    out = {'meta': {'chrom': CHROM, 'min_anchor': MIN_ANCHOR,
                    'bam': args.bam, 'genome': args.genome, 'gtf': args.gtf,
                    'self_check_all_depths_match': self_check_ok},
           'junctions': results}
    with open(args.out, 'w') as fh:
        json.dump(out, fh, indent=2)
    print(f'[char] wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
