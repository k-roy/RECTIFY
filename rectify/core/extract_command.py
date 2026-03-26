"""
Extract Command for RECTIFY.

Extracts per-read information from BAM files into TSV format.

Output includes:
- Read coordinates (5' end, 3' end, aligned span)
- Junction information (spliced reads)
- Soft-clip information
- COMPASS/RECTIFY correction columns (if available)

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam

logger = logging.getLogger(__name__)


def create_extract_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Create extract subcommand parser."""
    parser = subparsers.add_parser(
        'extract',
        help='Extract per-read information from BAM to TSV',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument(
        'bam',
        type=Path,
        help='Input BAM file'
    )

    parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output TSV file'
    )

    # Reference files
    ref_group = parser.add_argument_group('Reference files')
    ref_group.add_argument(
        '--genome',
        type=Path,
        help='Reference genome FASTA (for sequence context)'
    )

    ref_group.add_argument(
        '--annotation',
        type=Path,
        help='Gene annotation GFF/GTF (for gene assignment)'
    )

    # Output columns
    col_group = parser.add_argument_group('Output columns')
    col_group.add_argument(
        '--include-sequence',
        action='store_true',
        help='Include read sequence in output'
    )

    col_group.add_argument(
        '--include-quality',
        action='store_true',
        help='Include quality string in output'
    )

    col_group.add_argument(
        '--include-junctions',
        action='store_true',
        default=True,
        help='Include junction coordinates for spliced reads'
    )

    col_group.add_argument(
        '--include-softclips',
        action='store_true',
        default=True,
        help='Include soft-clip information'
    )

    col_group.add_argument(
        '--include-context',
        action='store_true',
        help='Include genomic context around 3\' end (requires --genome)'
    )

    col_group.add_argument(
        '--context-length',
        type=int,
        default=20,
        help='Length of genomic context to extract'
    )

    # Filter options
    filter_group = parser.add_argument_group('Filter options')
    filter_group.add_argument(
        '--chrom',
        help='Process only this chromosome'
    )

    filter_group.add_argument(
        '--region',
        help='Process only this region (chr:start-end)'
    )

    filter_group.add_argument(
        '--strand',
        choices=['+', '-'],
        help='Process only this strand'
    )

    filter_group.add_argument(
        '--min-mapq',
        type=int,
        default=0,
        help='Minimum mapping quality'
    )

    filter_group.add_argument(
        '--min-length',
        type=int,
        default=0,
        help='Minimum aligned length'
    )

    filter_group.add_argument(
        '--spliced-only',
        action='store_true',
        help='Output only spliced reads'
    )

    # Performance options
    perf_group = parser.add_argument_group('Performance options')
    perf_group.add_argument(
        '--streaming',
        action='store_true',
        help='Stream output to file (lower memory)'
    )

    perf_group.add_argument(
        '--chunk-size',
        type=int,
        default=100000,
        help='Chunk size for streaming mode'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose logging'
    )

    return parser


def extract_read_info(
    read: pysam.AlignedSegment,
    genome: Optional[dict] = None,
    context_length: int = 20,
    include_sequence: bool = False,
    include_quality: bool = False,
    include_junctions: bool = True,
    include_softclips: bool = True,
    include_context: bool = False,
) -> dict:
    """Extract information from a single read.

    Args:
        read: pysam AlignedSegment
        genome: Optional genome dict for context extraction
        context_length: Length of context to extract
        include_sequence: Include read sequence
        include_quality: Include quality string
        include_junctions: Include junction info
        include_softclips: Include soft-clip info
        include_context: Include genomic context

    Returns:
        Dict with read information
    """
    chrom = read.reference_name
    strand = '-' if read.is_reverse else '+'

    # Basic coordinates
    if strand == '+':
        five_prime = read.reference_start
        three_prime = read.reference_end - 1
    else:
        five_prime = read.reference_end - 1
        three_prime = read.reference_start

    aligned_length = read.reference_end - read.reference_start

    info = {
        'read_id': read.query_name,
        'chrom': chrom,
        'strand': strand,
        'reference_start': read.reference_start,
        'reference_end': read.reference_end,
        'five_prime': five_prime,
        'three_prime': three_prime,
        'aligned_length': aligned_length,
        'mapq': read.mapping_quality,
    }

    # Junction information
    if include_junctions and read.cigartuples:
        junctions = []
        ref_pos = read.reference_start
        for op, length in read.cigartuples:
            if op == 3:  # N = skipped region (intron)
                junctions.append(f"{ref_pos}-{ref_pos + length}")
                ref_pos += length
            elif op in (0, 2, 7, 8):  # M, D, =, X consume reference
                ref_pos += length

        info['n_junctions'] = len(junctions)
        info['junctions'] = ';'.join(junctions) if junctions else ''
        info['is_spliced'] = len(junctions) > 0

    # Soft-clip information
    if include_softclips and read.cigartuples:
        left_clip = 0
        right_clip = 0

        if read.cigartuples[0][0] == 4:  # S
            left_clip = read.cigartuples[0][1]
        if read.cigartuples[-1][0] == 4:  # S
            right_clip = read.cigartuples[-1][1]

        info['left_softclip'] = left_clip
        info['right_softclip'] = right_clip

        # 5' soft-clip (depends on strand)
        if strand == '+':
            info['five_prime_softclip'] = left_clip
            info['three_prime_softclip'] = right_clip
        else:
            info['five_prime_softclip'] = right_clip
            info['three_prime_softclip'] = left_clip

    # Genomic context
    if include_context and genome and chrom in genome:
        genome_seq = genome[chrom]
        if strand == '+':
            # Context downstream of 3' end
            ctx_start = three_prime + 1
            ctx_end = min(ctx_start + context_length, len(genome_seq))
            context = genome_seq[ctx_start:ctx_end]
        else:
            # Context downstream of 3' end (upstream in genome)
            ctx_end = three_prime
            ctx_start = max(0, ctx_end - context_length)
            context = reverse_complement(genome_seq[ctx_start:ctx_end])

        info['downstream_context'] = context

    # Read sequence
    if include_sequence:
        info['sequence'] = read.query_sequence

    # Quality
    if include_quality and read.query_qualities is not None:
        info['quality'] = ''.join(chr(q + 33) for q in read.query_qualities)

    return info


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def run_extract(args: argparse.Namespace) -> int:
    """Run extract command."""
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Validate input
    if not args.bam.exists():
        logger.error(f"BAM file not found: {args.bam}")
        return 1

    # Load genome if needed
    genome = None
    if args.genome and args.genome.exists():
        if args.include_context:
            logger.info(f"Loading genome from {args.genome}")
            genome = {}
            fasta = pysam.FastaFile(str(args.genome))
            for chrom in fasta.references:
                genome[chrom] = fasta.fetch(chrom)
            fasta.close()

    # Create output directory
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Parse region if provided
    region_chrom = args.chrom
    region_start = None
    region_end = None

    if args.region:
        parts = args.region.replace('-', ':').split(':')
        region_chrom = parts[0]
        if len(parts) >= 2:
            region_start = int(parts[1])
        if len(parts) >= 3:
            region_end = int(parts[2])

    # Process BAM
    bam = pysam.AlignmentFile(str(args.bam), 'rb')

    if args.streaming:
        # Streaming mode - write chunks directly
        run_streaming_extract(
            bam, args, genome, region_chrom, region_start, region_end
        )
    else:
        # Batch mode - collect all then write
        run_batch_extract(
            bam, args, genome, region_chrom, region_start, region_end
        )

    bam.close()

    logger.info(f"Extraction complete: {args.output}")
    return 0


def run_streaming_extract(bam, args, genome, region_chrom, region_start, region_end):
    """Run extraction in streaming mode."""
    n_written = 0
    chunk = []
    first_chunk = True

    for read in bam.fetch(region_chrom, region_start, region_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Apply filters
        if args.min_mapq > 0 and read.mapping_quality < args.min_mapq:
            continue

        aligned_length = read.reference_end - read.reference_start
        if aligned_length < args.min_length:
            continue

        strand = '-' if read.is_reverse else '+'
        if args.strand and strand != args.strand:
            continue

        has_junction = any(op == 3 for op, _ in (read.cigartuples or []))
        if args.spliced_only and not has_junction:
            continue

        # Extract info
        info = extract_read_info(
            read,
            genome=genome,
            context_length=args.context_length,
            include_sequence=args.include_sequence,
            include_quality=args.include_quality,
            include_junctions=args.include_junctions,
            include_softclips=args.include_softclips,
            include_context=args.include_context,
        )
        chunk.append(info)
        n_written += 1

        # Write chunk
        if len(chunk) >= args.chunk_size:
            df = pd.DataFrame(chunk)
            df.to_csv(
                args.output,
                sep='\t',
                index=False,
                mode='w' if first_chunk else 'a',
                header=first_chunk
            )
            chunk = []
            first_chunk = False
            logger.debug(f"Written {n_written:,} reads")

    # Write remaining
    if chunk:
        df = pd.DataFrame(chunk)
        df.to_csv(
            args.output,
            sep='\t',
            index=False,
            mode='w' if first_chunk else 'a',
            header=first_chunk
        )

    logger.info(f"Extracted {n_written:,} reads")


def run_batch_extract(bam, args, genome, region_chrom, region_start, region_end):
    """Run extraction in batch mode."""
    records = []

    for read in bam.fetch(region_chrom, region_start, region_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Apply filters
        if args.min_mapq > 0 and read.mapping_quality < args.min_mapq:
            continue

        aligned_length = read.reference_end - read.reference_start
        if aligned_length < args.min_length:
            continue

        strand = '-' if read.is_reverse else '+'
        if args.strand and strand != args.strand:
            continue

        has_junction = any(op == 3 for op, _ in (read.cigartuples or []))
        if args.spliced_only and not has_junction:
            continue

        # Extract info
        info = extract_read_info(
            read,
            genome=genome,
            context_length=args.context_length,
            include_sequence=args.include_sequence,
            include_quality=args.include_quality,
            include_junctions=args.include_junctions,
            include_softclips=args.include_softclips,
            include_context=args.include_context,
        )
        records.append(info)

    # Write output
    df = pd.DataFrame(records)
    df.to_csv(args.output, sep='\t', index=False)

    logger.info(f"Extracted {len(df):,} reads")


def run(args: argparse.Namespace):
    """Entry point for CLI."""
    sys.exit(run_extract(args))
