#!/usr/bin/env python3
"""
RECTIFY 'train-polya' command implementation.

Trains poly(A) tail model from control data (sites with 0 downstream A's).
At control sites, any A-rich soft-clipped sequence must genuinely be a
poly(A) tail, not genomic sequence misaligned.

Author: Kevin R. Roy
Date: 2026-03-16
"""

from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
import sys

import numpy as np
import pysam
from tqdm import tqdm

from ..slurm import set_thread_limits, get_available_cpus
from ..utils.genome import load_genome, reverse_complement
from ..utils.alignment import extract_soft_clips
from .polya_model import PolyAModel, PolyAModelParameters, PolyAModelStats
from .atract_detector import calculate_atract_ambiguity

logger = logging.getLogger(__name__)


@dataclass
class ControlSite:
    """A control CPA site with 0 downstream A's."""
    chrom: str
    position: int
    strand: str
    gene_name: Optional[str] = None


@dataclass
class TrainingData:
    """Collected training data from control sites."""
    soft_clip_sequences: List[str]
    soft_clip_lengths: List[int]
    a_richness_values: List[float]
    position_base_counts: Dict[int, Counter]  # position -> Counter({'A': n, 'G': n, ...})
    n_sites_used: int = 0
    n_sites_filtered: int = 0
    n_reads_total: int = 0


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for command execution."""
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def load_control_sites(filepath: Path) -> List[ControlSite]:
    """
    Load control CPA sites from TSV file.

    Expected format (tab-separated, header required):
        chrom    position    strand    [gene_name]    [downstream_a_count]

    Sites with downstream_a_count > 0 will be filtered out.

    Args:
        filepath: Path to TSV file

    Returns:
        List of ControlSite objects
    """
    sites = []
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Control sites file not found: {filepath}")

    with open(filepath, 'r') as f:
        # Read header
        header = f.readline().strip().split('\t')
        header_lower = [h.lower() for h in header]

        # Find column indices
        try:
            chrom_idx = header_lower.index('chrom')
        except ValueError:
            chrom_idx = 0

        try:
            pos_idx = header_lower.index('position')
        except ValueError:
            pos_idx = 1

        try:
            strand_idx = header_lower.index('strand')
        except ValueError:
            strand_idx = 2

        gene_idx = None
        for name in ['gene_name', 'gene', 'name']:
            if name in header_lower:
                gene_idx = header_lower.index(name)
                break

        acount_idx = None
        for name in ['downstream_a_count', 'a_count', 'acount']:
            if name in header_lower:
                acount_idx = header_lower.index(name)
                break

        # Read data lines
        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 3:
                logger.warning(f"Line {line_num}: insufficient fields, skipping")
                continue

            try:
                chrom = fields[chrom_idx]
                position = int(fields[pos_idx])
                strand = fields[strand_idx]

                # Filter out sites with downstream A's
                if acount_idx is not None and len(fields) > acount_idx:
                    a_count = int(fields[acount_idx])
                    if a_count > 0:
                        continue

                gene_name = None
                if gene_idx is not None and len(fields) > gene_idx:
                    gene_name = fields[gene_idx]

                sites.append(ControlSite(
                    chrom=chrom,
                    position=position,
                    strand=strand,
                    gene_name=gene_name
                ))

            except (ValueError, IndexError) as e:
                logger.warning(f"Line {line_num}: parse error ({e}), skipping")
                continue

    logger.info(f"Loaded {len(sites)} control sites from {filepath}")
    return sites


def verify_control_site(
    genome: Dict[str, str],
    site: ControlSite,
    window_size: int = 10
) -> bool:
    """
    Verify that a control site truly has 0 downstream A's.

    Args:
        genome: Genome dict from load_genome()
        site: Control site to verify
        window_size: Window to check for A's

    Returns:
        True if site has 0 downstream A's
    """
    ambiguity = calculate_atract_ambiguity(
        genome, site.chrom, site.position, site.strand,
        downstream_bp=window_size
    )
    return ambiguity['downstream_a_count'] == 0


def extract_soft_clips_at_site(
    bam: pysam.AlignmentFile,
    site: ControlSite,
    tolerance: int = 5,
    min_mapq: int = 10
) -> List[str]:
    """
    Extract 3' soft-clip sequences from reads at a control site.

    Args:
        bam: Opened BAM file
        site: Control CPA site
        tolerance: Position tolerance (±bp)
        min_mapq: Minimum mapping quality

    Returns:
        List of soft-clip sequences (in RNA orientation, i.e., A-rich)
    """
    soft_clips = []

    # Fetch reads overlapping the site
    start = max(0, site.position - tolerance)
    end = site.position + tolerance + 1

    try:
        for read in bam.fetch(site.chrom, start, end):
            # Filter reads
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            # Check if read's 3' end is near the site
            if site.strand == '+':
                read_3prime = read.reference_end - 1 if read.reference_end else None
            else:
                read_3prime = read.reference_start

            if read_3prime is None:
                continue

            if abs(read_3prime - site.position) > tolerance:
                continue

            # Extract soft-clips
            clips = extract_soft_clips(read)

            # Get 3' soft-clip
            if site.strand == '+':
                clip = next((c for c in clips if c['side'] == 'right'), None)
            else:
                clip = next((c for c in clips if c['side'] == 'left'), None)

            if clip and clip['length'] > 0:
                seq = clip['seq']
                # Convert to RNA orientation (A-rich)
                if site.strand == '-':
                    seq = reverse_complement(seq)
                soft_clips.append(seq)

    except ValueError:
        # Chromosome not in BAM
        pass

    return soft_clips


def collect_training_data(
    bam_path: Path,
    genome: Dict[str, str],
    control_sites: List[ControlSite],
    min_reads: int = 10,
    verify_sites: bool = True
) -> TrainingData:
    """
    Collect training data from control sites.

    Args:
        bam_path: Path to BAM file
        genome: Genome dict
        control_sites: List of control sites
        min_reads: Minimum reads per site to include
        verify_sites: Whether to verify 0A status against genome

    Returns:
        TrainingData object with collected sequences and statistics
    """
    training_data = TrainingData(
        soft_clip_sequences=[],
        soft_clip_lengths=[],
        a_richness_values=[],
        position_base_counts={},
    )

    bam = pysam.AlignmentFile(str(bam_path), 'rb')

    logger.info(f"Collecting training data from {len(control_sites)} control sites...")

    for site in tqdm(control_sites, desc="Processing sites"):
        # Optionally verify site
        if verify_sites and not verify_control_site(genome, site):
            training_data.n_sites_filtered += 1
            continue

        # Extract soft-clips
        clips = extract_soft_clips_at_site(bam, site)
        training_data.n_reads_total += len(clips)

        if len(clips) < min_reads:
            training_data.n_sites_filtered += 1
            continue

        training_data.n_sites_used += 1

        for clip_seq in clips:
            training_data.soft_clip_sequences.append(clip_seq)
            training_data.soft_clip_lengths.append(len(clip_seq))

            # Calculate A-richness
            if len(clip_seq) > 0:
                a_richness = clip_seq.upper().count('A') / len(clip_seq)
                training_data.a_richness_values.append(a_richness)

            # Track base composition by position
            for i, base in enumerate(clip_seq.upper()):
                if i not in training_data.position_base_counts:
                    training_data.position_base_counts[i] = Counter()
                training_data.position_base_counts[i][base] += 1

    bam.close()

    logger.info(f"Collected {len(training_data.soft_clip_sequences)} soft-clips "
                f"from {training_data.n_sites_used} sites")
    logger.info(f"Filtered {training_data.n_sites_filtered} sites "
                f"(failed verification or insufficient reads)")

    return training_data


def calculate_a_richness_distribution(values: List[float]) -> Dict[str, float]:
    """Calculate A-richness distribution statistics."""
    if not values:
        return {
            'mean': 0.91,
            'std': 0.08,
            'median': 0.92,
            'min': 0.0,
            'max': 1.0,
        }

    arr = np.array(values)
    return {
        'mean': float(np.mean(arr)),
        'std': float(np.std(arr)),
        'median': float(np.median(arr)),
        'min': float(np.min(arr)),
        'max': float(np.max(arr)),
    }


def calculate_length_distribution(lengths: List[int]) -> Dict[str, float]:
    """Calculate length distribution statistics."""
    if not lengths:
        return {
            'mean': 28.5,
            'std': 15.2,
            'median': 25.0,
            'min': 0,
            'max': 200,
            'p10': 10.0,
            'p90': 50.0,
        }

    arr = np.array(lengths)
    return {
        'mean': float(np.mean(arr)),
        'std': float(np.std(arr)),
        'median': float(np.median(arr)),
        'min': int(np.min(arr)),
        'max': int(np.max(arr)),
        'p10': float(np.percentile(arr, 10)),
        'p90': float(np.percentile(arr, 90)),
    }


def calculate_position_profile(
    position_counts: Dict[int, Counter],
    max_position: int = 50
) -> Dict[int, float]:
    """
    Calculate position-dependent A-frequency profile.

    Args:
        position_counts: Dict of position -> Counter({'A': n, 'G': n, ...})
        max_position: Maximum position to include

    Returns:
        Dict of position -> A-frequency
    """
    profile = {}

    for pos in range(max_position):
        if pos not in position_counts:
            # Assume high A-frequency for unobserved positions
            profile[pos] = 0.95
            continue

        counts = position_counts[pos]
        total = sum(counts.values())

        if total > 0:
            a_freq = counts.get('A', 0) / total
            profile[pos] = a_freq
        else:
            profile[pos] = 0.95

    return profile


def calculate_non_a_frequencies(
    position_counts: Dict[int, Counter]
) -> Dict[str, float]:
    """Calculate overall frequency of non-A bases."""
    total_counts = Counter()

    for counts in position_counts.values():
        total_counts.update(counts)

    total = sum(total_counts.values())
    if total == 0:
        return {'G': 0.04, 'C': 0.02, 'T': 0.01, 'N': 0.0}

    result = {}
    for base in ['G', 'C', 'T', 'N']:
        result[base] = total_counts.get(base, 0) / total

    return result


def determine_optimal_threshold(a_richness_values: List[float]) -> float:
    """
    Determine optimal A-richness threshold from training data.

    Uses the 10th percentile of observed A-richness values, but
    clamps to reasonable range.
    """
    if not a_richness_values:
        return 0.8

    p10 = np.percentile(a_richness_values, 10)

    # Clamp to reasonable range
    threshold = max(0.6, min(0.95, p10 - 0.05))  # Subtract small margin

    return float(threshold)


def fit_model(
    training_data: TrainingData,
    technology: str = "unknown",
    genome_name: str = "unknown"
) -> PolyAModel:
    """
    Fit poly(A) model from training data.

    Args:
        training_data: Collected training data
        technology: Sequencing technology name
        genome_name: Genome/organism name

    Returns:
        Fitted PolyAModel
    """
    logger.info("Fitting poly(A) model from training data...")

    # Calculate distributions
    a_richness_dist = calculate_a_richness_distribution(training_data.a_richness_values)
    length_dist = calculate_length_distribution(training_data.soft_clip_lengths)

    # Calculate profiles
    position_profile = calculate_position_profile(training_data.position_base_counts)
    non_a_freq = calculate_non_a_frequencies(training_data.position_base_counts)

    # Determine threshold
    threshold = determine_optimal_threshold(training_data.a_richness_values)

    # Create model
    model = PolyAModel(
        version="1.0",
        created=datetime.now().isoformat(),
        technology=technology,
        genome=genome_name,
        parameters=PolyAModelParameters(
            a_richness_threshold=threshold,
            a_richness_mean=a_richness_dist['mean'],
            a_richness_std=a_richness_dist['std'],
            min_tail_length=max(5, int(length_dist['p10'])),
            mean_tail_length=length_dist['mean'],
            median_tail_length=length_dist['median'],
            std_tail_length=length_dist['std'],
            max_tail_length=int(length_dist['max']),
        ),
        training_stats=PolyAModelStats(
            n_control_sites=training_data.n_sites_used,
            n_reads_used=training_data.n_reads_total,
            n_soft_clips_analyzed=len(training_data.soft_clip_sequences),
            n_sites_filtered=training_data.n_sites_filtered,
        ),
        position_profile=position_profile,
        non_a_frequencies=non_a_freq,
    )

    logger.info(f"Fitted model with threshold={threshold:.2f}, "
                f"mean_length={length_dist['mean']:.1f}")

    return model


def validate_inputs(args) -> Dict:
    """Validate command arguments and return config."""
    errors = []

    # Check BAM file
    if not args.bam.exists():
        errors.append(f"BAM file not found: {args.bam}")

    # Check genome file
    if not args.genome.exists():
        errors.append(f"Genome file not found: {args.genome}")

    # Check control sites file
    if not args.control_sites.exists():
        errors.append(f"Control sites file not found: {args.control_sites}")

    # Check output directory exists
    output_dir = args.output.parent
    if not output_dir.exists():
        try:
            output_dir.mkdir(parents=True)
        except Exception as e:
            errors.append(f"Cannot create output directory {output_dir}: {e}")

    if errors:
        for error in errors:
            logger.error(error)
        sys.exit(1)

    return {
        'min_reads': args.min_reads,
    }


def run(args) -> None:
    """
    Execute the 'train-polya' command.

    Args:
        args: Parsed command-line arguments from argparse
    """
    # Setup
    set_thread_limits(getattr(args, 'threads', None))
    setup_logging(getattr(args, 'verbose', False))

    logger.info("=" * 70)
    logger.info("RECTIFY - Poly(A) Tail Model Training")
    logger.info("=" * 70)
    logger.info("")

    # Validate inputs
    config = validate_inputs(args)

    logger.info("Inputs:")
    logger.info(f"  BAM file:         {args.bam}")
    logger.info(f"  Genome:           {args.genome}")
    logger.info(f"  Control sites:    {args.control_sites}")
    logger.info(f"  Output model:     {args.output}")
    logger.info(f"  Min reads/site:   {config['min_reads']}")
    logger.info("")

    # Load genome
    logger.info("Loading genome...")
    genome = load_genome(args.genome)
    logger.info(f"  Loaded {len(genome)} chromosomes")

    # Determine genome name from filename
    genome_name = args.genome.stem

    # Determine technology from BAM filename (heuristic)
    bam_name = args.bam.name.lower()
    if 'nanopore' in bam_name or 'ont' in bam_name:
        technology = 'nanopore'
    elif 'quantseq' in bam_name:
        technology = 'quantseq'
    elif 'helicos' in bam_name:
        technology = 'helicos'
    else:
        technology = 'unknown'

    # Load control sites
    logger.info("Loading control sites...")
    control_sites = load_control_sites(args.control_sites)

    if len(control_sites) == 0:
        logger.error("No control sites loaded. Check file format.")
        sys.exit(1)

    # Collect training data
    training_data = collect_training_data(
        args.bam,
        genome,
        control_sites,
        min_reads=config['min_reads'],
        verify_sites=True
    )

    if len(training_data.soft_clip_sequences) == 0:
        logger.error("No training data collected. Check that BAM has reads at control sites.")
        sys.exit(1)

    # Fit model
    model = fit_model(training_data, technology=technology, genome_name=genome_name)

    # Save model
    logger.info(f"Saving model to {args.output}...")
    model.to_json(args.output)

    # Print summary
    logger.info("")
    logger.info(model.summary())

    logger.info("")
    logger.info("=" * 70)
    logger.info("Training complete!")
    logger.info("=" * 70)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('bam', type=Path)
    parser.add_argument('--genome', type=Path, required=True)
    parser.add_argument('--control-sites', type=Path, required=True)
    parser.add_argument('-o', '--output', type=Path, required=True)
    parser.add_argument('--min-reads', type=int, default=10)
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()
    run(args)
