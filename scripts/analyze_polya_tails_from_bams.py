#!/usr/bin/env python3
"""
Analyze poly(A) tail lengths from soft-clipped reads at 0A control sites.

Strategy:
1. Use WT BY4742 Nanopore reads to identify top 5000 CPA sites at 0A positions
   (where there are no downstream genomic A's)
2. For each NET-seq BAM file, extract soft-clip information at these sites
3. Aggregate to get the oligo(A) tail length distribution

Author: Kevin R. Roy
Date: 2026-03-19
"""

import pysam
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import json

# Paths
BASE_DIR = Path("/oak/stanford/groups/larsms/Users/kevinroy")
RECTIFY_DIR = BASE_DIR / "software/rectify"
GENOME_PATH = BASE_DIR / "scripts_and_keyfiles/by_project/common_reference_files/MAGESTIC_background_strain.fasta"

# WT BY4742 Nanopore BAMs (to define CPA sites)
NANOPORE_BAMS = [
    BASE_DIR / "projects/roadblocks/intermediate_data/nanopore/inhouse_by4742_dst1_4nqo/wt_by4742_rep1.sorted.bam",
    BASE_DIR / "projects/roadblocks/intermediate_data/nanopore/inhouse_by4742_dst1_4nqo/wt_by4742_rep2.sorted.bam",
    BASE_DIR / "projects/roadblocks/intermediate_data/nanopore/inhouse_by4742_dst1_4nqo/wt_by4742_rep3.sorted.bam",
]

# NET-seq BAMs (to analyze soft-clips)
# NOTE: Use aligned_hexamer_trimmed for proper trimming:
#   - First 6bp (random hexamer) removed
#   - Adapter INCLUDING linker A removed
# This preserves oligo-A tails at the read 5' end
NETSEQ_BAM_DIRS = [
    BASE_DIR / "projects/roadblocks/intermediate_data/netseq/GSE159603_couvillion_2022/aligned_hexamer_trimmed",
    BASE_DIR / "projects/roadblocks/intermediate_data/netseq/SRA054029_legacy",
]

OUTPUT_DIR = RECTIFY_DIR / "scripts/polya_analysis_results"

# Parameters
DOWNSTREAM_WINDOW = 10  # bp to check for genomic A's
TOP_N_SITES = 5000      # Number of top CPA sites to use
MIN_SOFTCLIP_LEN = 1    # Minimum soft-clip length to consider
MIN_A_RICHNESS = 0.7    # Minimum A-richness to be considered poly(A)
MAX_WORKERS = 8         # Parallel workers

# Chromosome name mappings (standard chrI -> BAM format)
CHROM_TO_BAM = {
    'chrI': 'ref|NC_001133|', 'chrII': 'ref|NC_001134|', 'chrIII': 'ref|NC_001135|',
    'chrIV': 'ref|NC_001136|', 'chrV': 'ref|NC_001137|', 'chrVI': 'ref|NC_001138|',
    'chrVII': 'ref|NC_001139|', 'chrVIII': 'ref|NC_001140|', 'chrIX': 'ref|NC_001141|',
    'chrX': 'ref|NC_001142|', 'chrXI': 'ref|NC_001143|', 'chrXII': 'ref|NC_001144|',
    'chrXIII': 'ref|NC_001145|', 'chrXIV': 'ref|NC_001146|', 'chrXV': 'ref|NC_001147|',
    'chrXVI': 'ref|NC_001148|', 'chrM': 'ref|NC_001224|',
}

BAM_TO_CHROM = {v: k for k, v in CHROM_TO_BAM.items()}


def load_genome(fasta_path: Path) -> Dict[str, str]:
    """Load genome from FASTA file."""
    print(f"Loading genome from {fasta_path}...")
    genome = {}
    current_chrom = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom:
                    genome[current_chrom] = ''.join(current_seq).upper()
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_chrom:
            genome[current_chrom] = ''.join(current_seq).upper()

    print(f"  Loaded {len(genome)} chromosomes")
    return genome


def count_downstream_as(genome: Dict[str, str], chrom: str, position: int,
                         strand: str, window: int = 10) -> int:
    """Count A's (or T's for minus strand) in downstream window."""
    seq = genome.get(chrom)
    if seq is None:
        return -1

    if strand == '+':
        start = position + 1
        end = min(position + 1 + window, len(seq))
        if start >= len(seq):
            return 0
        downstream_seq = seq[start:end]
        return downstream_seq.count('A')
    else:
        start = max(0, position - window)
        end = position
        if end <= 0:
            return 0
        downstream_seq = seq[start:end]
        return downstream_seq.count('T')


def get_3prime_end(read) -> Optional[Tuple[str, str, int]]:
    """Get the 3' end position of a read."""
    if read.is_unmapped:
        return None

    chrom = read.reference_name
    strand = '-' if read.is_reverse else '+'

    if strand == '+':
        # 3' end is at the right (reference_end is 0-based exclusive)
        pos_3prime = read.reference_end - 1
    else:
        # 3' end is at the left (reference_start is 0-based)
        pos_3prime = read.reference_start

    return (chrom, strand, pos_3prime)


def identify_cpa_sites_from_nanopore(
    nanopore_bams: List[Path],
    genome: Dict[str, str],
    downstream_window: int = 10,
    top_n: int = 5000
) -> List[Tuple[str, str, int, int]]:
    """
    Identify top CPA sites from Nanopore reads at 0A positions.

    Returns list of (chrom, strand, position, read_count) sorted by read_count.
    """
    print(f"\nIdentifying CPA sites from {len(nanopore_bams)} Nanopore BAMs...")

    # Count 3' ends across all BAMs
    end_counts = Counter()

    for bam_path in nanopore_bams:
        print(f"  Processing {bam_path.name}...")
        try:
            bam = pysam.AlignmentFile(bam_path, "rb")
            for read in bam.fetch():
                end = get_3prime_end(read)
                if end:
                    end_counts[end] += 1
            bam.close()
        except Exception as e:
            print(f"    Error: {e}")

    print(f"  Total unique 3' end positions: {len(end_counts):,}")

    # Filter to 0A sites
    print(f"  Filtering to 0A sites (no downstream A's in {downstream_window}bp)...")
    zero_a_sites = []

    for (chrom, strand, position), count in end_counts.items():
        # Standardize chromosome name
        std_chrom = BAM_TO_CHROM.get(chrom, chrom)
        a_count = count_downstream_as(genome, std_chrom, position, strand, downstream_window)
        if a_count == 0:
            zero_a_sites.append((chrom, strand, position, count))

    print(f"  0A sites: {len(zero_a_sites):,}")

    # Sort by count and take top N
    zero_a_sites.sort(key=lambda x: x[3], reverse=True)
    top_sites = zero_a_sites[:top_n]

    print(f"  Selected top {len(top_sites):,} sites")
    if len(top_sites) > 0:
        print(f"    Min reads: {top_sites[-1][3]}")
        print(f"    Max reads: {top_sites[0][3]}")
        print(f"    Median reads: {top_sites[len(top_sites)//2][3]}")

    return top_sites


def get_softclip_info(read, strand: str) -> Optional[Dict]:
    """
    Extract soft-clip information from a read.

    For plus strand: poly(A) tail is on the RIGHT (3' end)
    For minus strand: poly(A) tail is on the LEFT (appears as T's)
    """
    if read.is_unmapped or read.cigartuples is None:
        return None

    cigar = read.cigartuples
    seq = read.query_sequence

    if seq is None:
        return None

    # For plus strand reads, poly(A) is at the 3' end (right soft-clip)
    if strand == '+':
        if cigar[-1][0] == 4:  # Soft-clip at end
            clip_len = cigar[-1][1]
            clip_seq = seq[-clip_len:]
            a_count = clip_seq.count('A')
            a_richness = a_count / clip_len if clip_len > 0 else 0
            return {
                'length': clip_len,
                'sequence': clip_seq,
                'a_richness': a_richness,
            }
    else:
        # For minus strand, check left soft-clip (T-rich = poly(A) on RNA)
        if cigar[0][0] == 4:  # Soft-clip at start
            clip_len = cigar[0][1]
            clip_seq = seq[:clip_len]
            t_count = clip_seq.count('T')
            t_richness = t_count / clip_len if clip_len > 0 else 0
            return {
                'length': clip_len,
                'sequence': clip_seq,
                'a_richness': t_richness,
            }

    return None


def analyze_single_bam(args) -> Dict:
    """
    Analyze soft-clips at control sites in a single BAM file.

    Efficient single-pass approach: iterate through all reads once and check
    if each read's 3' position is in the control site set.
    """
    bam_path, control_sites, min_softclip_len, min_a_richness = args

    results = {
        'bam': str(bam_path.name),
        'tail_lengths': [],
        'a_richness_values': [],
        'site_ids': [],  # Track which site each read comes from
        'n_reads_analyzed': 0,
        'n_polya_reads': 0,
    }

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return results

    # Build dict mapping positions to site IDs for O(1) lookup
    # Key: (chrom, strand, position), Value: site_id (index in control_sites)
    position_to_site = {}
    for site_id, (chrom, strand, position, count) in enumerate(control_sites):
        # Control sites use Nanopore format (chrI) - NET-seq reads will be converted
        # Map the exact position and neighboring positions (±5bp window)
        for offset in range(-5, 6):
            pos = position + offset
            # Only map to closest site (don't overwrite if already mapped)
            if (chrom, strand, pos) not in position_to_site:
                position_to_site[(chrom, strand, pos)] = site_id

    # Single pass through BAM
    for read in bam.fetch():
        if read.is_unmapped:
            continue

        results['n_reads_analyzed'] += 1

        chrom = read.reference_name
        strand = '-' if read.is_reverse else '+'

        # Get 3' position
        if strand == '+':
            pos_3prime = read.reference_end - 1 if read.reference_end else None
        else:
            pos_3prime = read.reference_start

        if pos_3prime is None:
            continue

        # Check if near any control site
        # Convert NET-seq chrom name (ref|NC_001133|) to standard (chrI) for lookup
        std_chrom = BAM_TO_CHROM.get(chrom, chrom)
        site_id = position_to_site.get((std_chrom, strand, pos_3prime))
        if site_id is None:
            continue

        # Get soft-clip info
        clip_info = get_softclip_info(read, strand)
        if clip_info is None:
            continue

        if clip_info['length'] < min_softclip_len:
            continue

        if clip_info['a_richness'] >= min_a_richness:
            results['tail_lengths'].append(clip_info['length'])
            results['a_richness_values'].append(clip_info['a_richness'])
            results['site_ids'].append(site_id)
            results['n_polya_reads'] += 1

    bam.close()
    return results


def main():
    print("=" * 70)
    print("Analyzing poly(A) tail lengths from BAM soft-clips at 0A control sites")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load genome
    genome = load_genome(GENOME_PATH)

    # Identify top 5k CPA sites from Nanopore data at 0A positions
    control_sites = identify_cpa_sites_from_nanopore(
        NANOPORE_BAMS, genome, DOWNSTREAM_WINDOW, TOP_N_SITES
    )

    if len(control_sites) == 0:
        print("ERROR: No control sites found!")
        return

    # Find all NET-seq BAM files
    netseq_bams = []
    for bam_dir in NETSEQ_BAM_DIRS:
        dir_bams = sorted(bam_dir.glob("*.bam"))
        print(f"  {bam_dir.name}: {len(dir_bams)} BAM files")
        netseq_bams.extend(dir_bams)
    print(f"\nTotal NET-seq BAMs: {len(netseq_bams)}")

    # Process BAM files in parallel
    print(f"\nProcessing BAM files with {MAX_WORKERS} workers...")

    all_tail_lengths = []
    all_a_richness = []
    all_site_ids = []
    total_reads = 0
    total_polya = 0

    args_list = [
        (bam_path, control_sites, MIN_SOFTCLIP_LEN, MIN_A_RICHNESS)
        for bam_path in netseq_bams
    ]

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(analyze_single_bam, args): args[0] for args in args_list}

        for i, future in enumerate(as_completed(futures)):
            try:
                result = future.result()
                all_tail_lengths.extend(result['tail_lengths'])
                all_a_richness.extend(result['a_richness_values'])
                all_site_ids.extend(result['site_ids'])
                total_reads += result['n_reads_analyzed']
                total_polya += result['n_polya_reads']

                if (i + 1) % 20 == 0:
                    print(f"  Processed {i + 1}/{len(netseq_bams)} BAMs...")
            except Exception as e:
                print(f"  Error: {e}")

    print(f"\nProcessing complete!")
    print(f"  Total reads analyzed: {total_reads:,}")
    print(f"  Reads with poly(A) tails: {total_polya:,}")

    # Compute and display statistics
    if len(all_tail_lengths) > 0:
        tail_lengths = np.array(all_tail_lengths)
        a_richness = np.array(all_a_richness)

        print("\n" + "=" * 70)
        print("RESULTS: Poly(A) Tail Length Distribution at 0A Control Sites")
        print("=" * 70)

        print(f"\nTail Length Statistics (n={len(tail_lengths):,}):")
        print(f"  Mean:   {np.mean(tail_lengths):.1f} bp")
        print(f"  Median: {np.median(tail_lengths):.1f} bp")
        print(f"  Std:    {np.std(tail_lengths):.1f} bp")
        print(f"  Min:    {np.min(tail_lengths)} bp")
        print(f"  Max:    {np.max(tail_lengths)} bp")

        print(f"\nA-richness Statistics:")
        print(f"  Mean:   {np.mean(a_richness):.3f}")
        print(f"  Median: {np.median(a_richness):.3f}")

        print(f"\nLength Distribution:")
        for length in range(1, 21):
            count = np.sum(tail_lengths == length)
            pct = count / len(tail_lengths) * 100
            bar = '#' * int(pct * 2)
            print(f"  {length:2d} bp: {count:6d} ({pct:5.1f}%) {bar}")

        # Per-site analysis
        site_ids = np.array(all_site_ids)
        unique_sites, site_counts = np.unique(site_ids, return_counts=True)

        print(f"\n" + "-" * 70)
        print(f"Per-Site Distribution (checking for concentration at single loci):")
        print(f"-" * 70)
        print(f"  Sites with poly(A) reads: {len(unique_sites):,} / {len(control_sites):,} ({100*len(unique_sites)/len(control_sites):.1f}%)")
        print(f"  Reads per site - Mean: {np.mean(site_counts):.1f}, Median: {np.median(site_counts):.1f}")
        print(f"  Reads per site - Min: {np.min(site_counts)}, Max: {np.max(site_counts)}")

        # Top 20 sites by read count
        top_indices = np.argsort(site_counts)[-20:][::-1]
        print(f"\n  Top 20 sites by poly(A) read count:")
        cumulative = 0
        for rank, idx in enumerate(top_indices, 1):
            site_id = unique_sites[idx]
            count = site_counts[idx]
            cumulative += count
            cum_pct = 100 * cumulative / len(tail_lengths)
            chrom, strand, position, nanopore_count = control_sites[site_id]
            print(f"    {rank:2d}. {chrom}:{position} ({strand}) - {count:,} reads ({100*count/len(tail_lengths):.1f}%, cumulative: {cum_pct:.1f}%)")

        # Concentration summary
        top10_reads = sum(site_counts[np.argsort(site_counts)[-10:][::-1]])
        top50_reads = sum(site_counts[np.argsort(site_counts)[-50:][::-1]])
        print(f"\n  Concentration:")
        print(f"    Top 10 sites: {top10_reads:,} reads ({100*top10_reads/len(tail_lengths):.1f}%)")
        print(f"    Top 50 sites: {top50_reads:,} reads ({100*top50_reads/len(tail_lengths):.1f}%)")

        # Save results
        output_file = OUTPUT_DIR / "polya_tail_analysis.json"
        results = {
            'n_control_sites': len(control_sites),
            'n_bam_files': len(netseq_bams),
            'total_reads_analyzed': total_reads,
            'total_polya_reads': total_polya,
            'tail_length_mean': float(np.mean(tail_lengths)),
            'tail_length_median': float(np.median(tail_lengths)),
            'tail_length_std': float(np.std(tail_lengths)),
            'length_distribution': {
                str(i): int(np.sum(tail_lengths == i))
                for i in range(int(np.min(tail_lengths)), int(np.max(tail_lengths)) + 1)
            }
        }

        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to {output_file}")

        np.save(OUTPUT_DIR / "tail_lengths.npy", tail_lengths)
        np.save(OUTPUT_DIR / "a_richness.npy", a_richness)

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
