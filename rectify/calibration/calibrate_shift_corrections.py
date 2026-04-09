#!/usr/bin/env python3
"""
Calibrate SHIFT_CORRECTIONS_BY_ACOUNT from NET-seq data.

This script derives the empirical shift correction table used by
calculate_atract_ambiguity() in rectify/core/atract_detector.py.

BACKGROUND
----------
Direct RNA-seq aligners (minimap2) place the called 3' end at the position
where the genomic sequence and poly(A) tail stop matching. When the terminal
mRNA sequence is A-rich (e.g., upstream of polyadenylation signals), the
observed tail length is shorter because some alignment positions are taken up
by genomic A's. This creates an apparent upstream shift in the called 3' end.

The calibration measures this shift empirically by comparing mean oligo-A tail
lengths at CPA sites stratified by upstream A-richness (total A's in the last
10 bp of the transcript body). Sites with more upstream A's show shorter tails,
meaning the aligner absorbed some of those A's into the alignment — shifting
the apparent 3' end upstream.

CALIBRATION METRIC
------------------
Metric: total A's in the upstream 10 bp window (gene-strand oriented).
This is the COUNT_UPSTREAM_A metric, not the downstream-contiguous count
used at runtime. The relationship between upstream A-richness and shift
magnitude is monotonically increasing and empirically robust across
S. cerevisiae NET-seq replicates.

WHAT YOU NEED TO RUN THIS
--------------------------
1. NET-seq BAM files aligned with soft-clip retention (e.g., bbmap with
   softclip=t), so oligo-A tails appear as soft-clipped bases.
2. A CPA cluster file (TSV with columns: chrom, strand, modal_position,
   wt_counts) — typically generated from NET-seq peak calling.
3. A reference genome FASTA.

HOW TO UPDATE config.py AFTER RUNNING
---------------------------------------
After running, the output file 'shift_correction_by_acount.tsv' contains
the new shift values. Copy the 'shift_correction_bp' column values into
SHIFT_CORRECTIONS_BY_ACOUNT in rectify/config.py.

USAGE
-----
    python -m rectify.calibration.calibrate_shift_corrections \\
        --bam-dir /scratch/netseq/aligned_softclip \\
        --samples netseq_wt1.bam,netseq_wt2.bam,netseq_wt3.bam \\
        --clusters /path/to/pa_clusters.tsv \\
        --genome /path/to/genome.fsa \\
        --output-dir /path/to/output \\
        --top-n 5000

    # Or run with the yeast roadblocks defaults:
    python -m rectify.calibration.calibrate_shift_corrections --use-yeast-defaults

Author: Kevin R. Roy
Date: 2026-04-02
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO

from rectify.utils.genome import reverse_complement  # noqa: F401 (re-used below)

# Default S. cerevisiae yeast config (roadblocks project)
YEAST_DEFAULTS = {
    'bam_dir': '/scratch/users/kevinroy/roadblocks_reprocess/aligned_softclip',
    'samples': {
        'WT_2022_rep1': 'netseq_2022_wt1_bbmap.bam',
        'WT_2022_rep2': 'netseq_2022_wt2_bbmap.bam',
        'WT_2022_rep3': 'netseq_2022_wt3_trimmed_bbmap.bam',
        'WT_2022_rep4': 'netseq_2022_wt4_trimmed_bbmap.bam',
        'WT_2011':      'netseq_2011_WT_bbmap_correct.bam',
    },
    'clusters': '/oak/stanford/groups/larsms/Users/kevinroy/projects/roadblocks/processed_data/clusters_adaptive/pa_clusters_adaptive.tsv',
    'genome':   '/oak/stanford/groups/larsms/Users/kevinroy/common/reference_genomes/S288C_reference_genome_R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529.chrnames.fsa',
    'output_dir': '/oak/stanford/groups/larsms/Users/kevinroy/projects/roadblocks/processed_data/softclip_analysis',
    'top_n': 5000,
    # Chromosome name mapping: standard name → BAM contig name
    'chrom_to_bam': {
        'chrI':    'ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]',
        'chrII':   'ref|NC_001134| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=II]',
        'chrIII':  'ref|NC_001135| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=III]',
        'chrIV':   'ref|NC_001136| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IV]',
        'chrV':    'ref|NC_001137| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=V]',
        'chrVI':   'ref|NC_001138| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VI]',
        'chrVII':  'ref|NC_001139| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VII]',
        'chrVIII': 'ref|NC_001140| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VIII]',
        'chrIX':   'ref|NC_001141| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IX]',
        'chrX':    'ref|NC_001142| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=X]',
        'chrXI':   'ref|NC_001143| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XI]',
        'chrXII':  'ref|NC_001144| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XII]',
        'chrXIII': 'ref|NC_001145| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIII]',
        'chrXIV':  'ref|NC_001146| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIV]',
        'chrXV':   'ref|NC_001147| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XV]',
        'chrXVI':  'ref|NC_001148| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XVI]',
    },
    'chrom_to_genome': {
        'chrI':    'ref|NC_001133|',
        'chrII':   'ref|NC_001134|',
        'chrIII':  'ref|NC_001135|',
        'chrIV':   'ref|NC_001136|',
        'chrV':    'ref|NC_001137|',
        'chrVI':   'ref|NC_001138|',
        'chrVII':  'ref|NC_001139|',
        'chrVIII': 'ref|NC_001140|',
        'chrIX':   'ref|NC_001141|',
        'chrX':    'ref|NC_001142|',
        'chrXI':   'ref|NC_001143|',
        'chrXII':  'ref|NC_001144|',
        'chrXIII': 'ref|NC_001145|',
        'chrXIV':  'ref|NC_001146|',
        'chrXV':   'ref|NC_001147|',
        'chrXVI':  'ref|NC_001148|',
    },
}


def load_genome(genome_file):
    print("Loading genome...")
    genome = {}
    for record in SeqIO.parse(str(genome_file), 'fasta'):
        genome[record.id] = str(record.seq).upper()
    print(f"  Loaded {len(genome)} chromosomes")
    return genome


def get_upstream_sequence(genome, chrom, pos, gene_strand, upstream_bp=10):
    """
    Get the last upstream_bp bases of transcript body ending at the 3' end pos.

    For + strand: genome[pos-10 : pos]  (gene-strand = forward)
    For - strand: reverse complement of genome[pos+1 : pos+11]

    This is the 'upstream 10bp window' used to stratify CPA sites by A-richness.
    More A's in this window → shorter observed oligo-A tail (alignment absorbs them).
    """
    seq = genome.get(chrom)
    if seq is None:
        return None
    if gene_strand == '+':
        raw = seq[max(0, pos - upstream_bp):pos]
    else:
        raw = seq[pos + 1:min(len(seq), pos + 1 + upstream_bp)]
        raw = reverse_complement(raw)
    return raw


def parse_cigar_stats(cigar):
    stats = {'matches': 0, 'insertions': 0, 'deletions': 0}
    for op, length in cigar:
        if op in (0, 7, 8):
            stats['matches'] += length
        elif op == 1:
            stats['insertions'] += length
        elif op == 2:
            stats['deletions'] += length
    return stats


def passes_alignment_quality(read, cigar_stats, max_mismatches=3, min_aligned=20):
    if cigar_stats['deletions'] > 0 or cigar_stats['insertions'] > 0:
        return False
    nm = read.get_tag('NM') if read.has_tag('NM') else 0
    if nm > max_mismatches:
        return False
    if cigar_stats['matches'] < min_aligned:
        return False
    return True


def find_oligo_a_length(clip_seq, read_strand, t_rich_threshold=0.6):
    """
    Detect oligo-A (or oligo-T) tail length in soft-clipped sequence.

    NET-seq reads have oligo-A added by the oligo-dA tailing protocol.
    For + gene strand, reads are on the - read strand (antisense), so the
    oligo-A appears as a leading A-run in the soft clip.
    For - gene strand, reads are on the + read strand, so it's a trailing T-run.
    """
    if len(clip_seq) < 1:
        return 0
    if read_strand == '-':
        # Leading A-run
        count = 0
        for base in clip_seq:
            if base == 'A':
                count += 1
            else:
                break
        if count >= 1:
            return count
        if clip_seq.count('A') / len(clip_seq) >= t_rich_threshold:
            return int(clip_seq.count('A'))
    else:
        # Trailing T-run
        count = 0
        for base in reversed(clip_seq):
            if base == 'T':
                count += 1
            else:
                break
        if count >= 1:
            return count
        if clip_seq.count('T') / len(clip_seq) >= t_rich_threshold:
            return int(clip_seq.count('T'))
    return 0


def analyze_sample(sample_name, bam_path, clusters_df, genome, chrom_to_bam, chrom_to_genome):
    """
    Analyze one NET-seq BAM, collecting oligo-A tail lengths stratified by
    upstream A-count at each CPA cluster site.

    Returns:
        results: dict {a_count: [tail_lengths]}
        read_counts: dict {a_count: n_reads}
        cpa_site_counts: dict {a_count: n_sites}
    """
    print(f"\nAnalyzing {sample_name}...")
    bam = pysam.AlignmentFile(str(bam_path), "rb")

    results = defaultdict(list)
    read_counts = defaultdict(int)
    cpa_sites_per_acount = defaultdict(set)

    for _, cluster in clusters_df.iterrows():
        chrom = cluster['chrom']
        bam_chrom = chrom_to_bam.get(chrom)
        genome_chrom = chrom_to_genome.get(chrom)
        if not bam_chrom or not genome_chrom:
            continue

        gene_strand = cluster['strand']
        modal_pos = int(cluster['modal_position'])
        cluster_id = f"{chrom}:{modal_pos}:{gene_strand}"
        expected_read_strand = '-' if gene_strand == '+' else '+'

        region_start = max(0, modal_pos - 10)
        region_end = modal_pos + 10

        try:
            for read in bam.fetch(bam_chrom, region_start, region_end):
                if read.is_unmapped:
                    continue
                cigar = read.cigartuples
                seq = read.query_sequence
                if not cigar or not seq:
                    continue

                read_strand = '-' if read.is_reverse else '+'
                if read_strand != expected_read_strand:
                    continue

                rna_3prime = read.reference_end - 1 if read_strand == '-' else read.reference_start
                if abs(rna_3prime - modal_pos) > 5:
                    continue

                cigar_stats = parse_cigar_stats(cigar)
                if not passes_alignment_quality(read, cigar_stats):
                    continue

                upstream_seq = get_upstream_sequence(genome, genome_chrom, rna_3prime, gene_strand)
                if not upstream_seq or len(upstream_seq) < 10:
                    continue

                a_count = upstream_seq.count('A')
                read_counts[a_count] += 1
                cpa_sites_per_acount[a_count].add(cluster_id)

                tail_len = 0
                if read_strand == '-':
                    if cigar[-1][0] == 4 and cigar[-1][1] >= 1:
                        clip_seq = seq[-cigar[-1][1]:]
                        tail_len = find_oligo_a_length(clip_seq, read_strand)
                else:
                    if cigar[0][0] == 4 and cigar[0][1] >= 1:
                        clip_seq = seq[:cigar[0][1]]
                        tail_len = find_oligo_a_length(clip_seq, read_strand)

                if tail_len >= 1:
                    results[a_count].append(tail_len)

        except Exception:
            continue

    bam.close()

    cpa_site_counts = {a: len(sites) for a, sites in cpa_sites_per_acount.items()}

    print(f"  Counts by upstream A count:")
    for a_count in range(11):
        n_reads = read_counts[a_count]
        n_cpa = cpa_site_counts.get(a_count, 0)
        n_tails = len(results[a_count])
        mean_tail = np.mean(results[a_count]) if results[a_count] else 0
        print(f"    {a_count}A: {n_cpa:,} CPA sites, {n_reads:,} reads, {n_tails:,} tails, mean={mean_tail:.2f}bp")

    return results, read_counts, cpa_site_counts


def compute_shift_corrections(sample_stats_by_name):
    """
    Compute shift corrections from per-sample tail length statistics.

    The shift for a_count=N is: baseline_tail_length - mean_tail_at_N.
    Baseline = 0A (no A-richness in upstream window).
    Negative shifts are clamped to 0.

    Args:
        sample_stats_by_name: dict {sample_name: {a_count: {'mean', 'std', 'n'}}}

    Returns:
        mean_tails_by_acount: dict {a_count: mean_tail}
        shift_corrections: dict {a_count: shift}
    """
    mean_tails_by_acount = {}
    for a_count in range(11):
        rep_means = [
            stats[a_count]['mean']
            for stats in sample_stats_by_name.values()
            if not np.isnan(stats[a_count]['mean'])
        ]
        if rep_means:
            mean_tails_by_acount[a_count] = np.mean(rep_means)

    if 0 not in mean_tails_by_acount:
        raise ValueError("No data for 0A sites — cannot compute baseline")

    baseline = mean_tails_by_acount[0]
    shift_corrections = {}
    for a_count in range(11):
        if a_count in mean_tails_by_acount:
            shift = baseline - mean_tails_by_acount[a_count]
            shift_corrections[a_count] = max(0.0, shift)
    shift_corrections[11] = shift_corrections.get(10, 0.0)

    return mean_tails_by_acount, shift_corrections


def print_config_snippet(shift_corrections):
    """Print the Python dict to paste into config.py."""
    print("\n" + "=" * 60)
    print("PASTE INTO rectify/config.py as SHIFT_CORRECTIONS_BY_ACOUNT:")
    print("=" * 60)
    print("SHIFT_CORRECTIONS_BY_ACOUNT = {")
    for a in range(12):
        val = shift_corrections.get(a, shift_corrections.get(10, 0.0))
        label = f"{a}+" if a == 11 else str(a)
        print(f"    {a}: {val:.4f},  # {label}A upstream")
    print("}")


def run_calibration(
    bam_dir,
    samples,            # dict {name: filename} or list of filenames
    cluster_file,
    genome_file,
    output_dir,
    top_n=5000,
    chrom_to_bam=None,
    chrom_to_genome=None,
):
    """
    Main calibration routine. Returns shift_corrections dict.

    Args:
        bam_dir: Path to directory containing NET-seq BAM files
        samples: dict {sample_name: bam_filename} or list of bam filenames
        cluster_file: Path to CPA cluster TSV (needs: chrom, strand, modal_position, wt_counts)
        genome_file: Path to reference genome FASTA
        output_dir: Where to write output TSV and plots
        top_n: Number of top CPA clusters (by wt_counts) to use
        chrom_to_bam: dict mapping standard chrom names to BAM contig names
        chrom_to_genome: dict mapping standard chrom names to genome FASTA keys
    """
    bam_dir = Path(bam_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(samples, list):
        samples = {Path(f).stem: f for f in samples}

    if chrom_to_bam is None:
        chrom_to_bam = {}
    if chrom_to_genome is None:
        chrom_to_genome = {}

    genome = load_genome(genome_file)

    print("\nLoading CPA clusters...")
    clusters_df = pd.read_csv(cluster_file, sep='\t')
    if chrom_to_bam:
        clusters_df = clusters_df[clusters_df['chrom'].isin(chrom_to_bam.keys())]
    clusters_df = clusters_df.sort_values('wt_counts', ascending=False).head(top_n)
    print(f"  Using top {len(clusters_df)} clusters by WT signal")

    all_results = {}
    all_read_counts = {}
    all_cpa_site_counts = {}
    for sample_name, bam_filename in samples.items():
        bam_path = bam_dir / bam_filename
        if not bam_path.exists():
            print(f"  WARNING: BAM not found: {bam_path} — skipping")
            continue
        results, read_counts, cpa_site_counts = analyze_sample(
            sample_name, bam_path, clusters_df, genome, chrom_to_bam, chrom_to_genome
        )
        all_results[sample_name] = results
        all_read_counts[sample_name] = read_counts
        all_cpa_site_counts[sample_name] = cpa_site_counts

    if not all_results:
        raise RuntimeError("No samples successfully analyzed")

    # Per-sample statistics
    sample_stats = {}
    for sample_name, results in all_results.items():
        sample_stats[sample_name] = {}
        for a_count in range(11):
            tails = results[a_count]
            if tails:
                sample_stats[sample_name][a_count] = {
                    'mean': np.mean(tails), 'std': np.std(tails), 'n': len(tails)
                }
            else:
                sample_stats[sample_name][a_count] = {'mean': np.nan, 'std': 0, 'n': 0}

    mean_tails_by_acount, shift_corrections = compute_shift_corrections(sample_stats)

    # Save detailed per-sample TSV
    rows = []
    for sample_name in all_results:
        for a_count in range(11):
            stats = sample_stats[sample_name][a_count]
            rows.append({
                'sample': sample_name,
                'upstream_a_count': a_count,
                'mean_tail_length': stats['mean'],
                'std_tail_length': stats['std'],
                'n_tails': stats['n'],
                'n_reads': all_read_counts[sample_name].get(a_count, 0),
            })
    pd.DataFrame(rows).to_csv(output_dir / 'tail_length_by_acount_all_samples.tsv', sep='\t', index=False)

    # Save shift correction table
    combined_cpa_counts = {
        a: max((all_cpa_site_counts[s].get(a, 0) for s in all_cpa_site_counts), default=0)
        for a in range(11)
    }
    shift_rows = [
        {
            'upstream_a_count': a if a < 11 else '11+',
            'mean_tail_length_bp': mean_tails_by_acount.get(min(a, 10), 0),
            'shift_correction_bp': shift_corrections.get(min(a, 10), 0),
            'n_cpa_sites': combined_cpa_counts.get(min(a, 10), 0),
        }
        for a in range(12)
    ]
    shift_file = output_dir / 'shift_correction_by_acount.tsv'
    pd.DataFrame(shift_rows).to_csv(shift_file, sep='\t', index=False)
    print(f"\nSaved: {shift_file}")

    print_config_snippet(shift_corrections)

    return shift_corrections


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--use-yeast-defaults', action='store_true',
                        help='Run with hardcoded yeast (S. cerevisiae) roadblocks defaults')
    parser.add_argument('--bam-dir', help='Directory containing NET-seq BAM files')
    parser.add_argument('--samples', help='Comma-separated list of BAM filenames (in --bam-dir)')
    parser.add_argument('--clusters', help='CPA cluster TSV file')
    parser.add_argument('--genome', help='Reference genome FASTA')
    parser.add_argument('--output-dir', help='Output directory')
    parser.add_argument('--top-n', type=int, default=5000,
                        help='Number of top CPA clusters to use (default: 5000)')
    args = parser.parse_args()

    if args.use_yeast_defaults:
        cfg = YEAST_DEFAULTS
        run_calibration(
            bam_dir=cfg['bam_dir'],
            samples=cfg['samples'],
            cluster_file=cfg['clusters'],
            genome_file=cfg['genome'],
            output_dir=cfg['output_dir'],
            top_n=cfg['top_n'],
            chrom_to_bam=cfg['chrom_to_bam'],
            chrom_to_genome=cfg['chrom_to_genome'],
        )
    else:
        if not all([args.bam_dir, args.samples, args.clusters, args.genome, args.output_dir]):
            parser.error("Must provide --bam-dir, --samples, --clusters, --genome, and --output-dir "
                         "(or use --use-yeast-defaults)")
        samples = args.samples.split(',')
        run_calibration(
            bam_dir=args.bam_dir,
            samples=samples,
            cluster_file=args.clusters,
            genome_file=args.genome,
            output_dir=args.output_dir,
            top_n=args.top_n,
        )


if __name__ == '__main__':
    main()
