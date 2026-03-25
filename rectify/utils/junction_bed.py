"""
Generate minimap2 junction BED files from GFF/GTF annotations.

minimap2's --junc-bed option improves alignment quality by preferring
annotated splice junctions, while still allowing novel junction detection.

This module provides utilities to:
1. Parse intron features from GFF files
2. Generate BED files for minimap2 --junc-bed

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import List, Tuple, Dict, Optional
import gzip
import logging

logger = logging.getLogger(__name__)


def parse_gff_introns(gff_path: str) -> List[Tuple[str, int, int, str, str]]:
    """Extract intron features from GFF file.

    Args:
        gff_path: Path to GFF file (can be gzipped)

    Returns:
        List of (chrom, start, end, name, strand) tuples
        Coordinates are 0-based, half-open (BED format)
    """
    introns = []

    # Handle gzipped files
    opener = gzip.open if gff_path.endswith('.gz') else open
    mode = 'rt' if gff_path.endswith('.gz') else 'r'

    with opener(gff_path, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs = fields

            # Match intron features (including UTR introns)
            if 'intron' in feature.lower():
                # Convert GFF 1-based to BED 0-based
                start_0based = int(start) - 1
                end_0based = int(end)  # BED is half-open, so end stays same

                # Extract name from attributes
                name = '.'
                for attr in attrs.split(';'):
                    attr = attr.strip()
                    if attr.startswith('Name='):
                        name = attr.split('=')[1]
                        break
                    elif attr.startswith('ID='):
                        name = attr.split('=')[1]
                        break

                introns.append((chrom, start_0based, end_0based, name, strand))

    return introns


def parse_gtf_introns(gtf_path: str) -> List[Tuple[str, int, int, str, str]]:
    """Infer intron coordinates from GTF exon features.

    GTF files typically don't have explicit intron features, so we infer
    introns from gaps between consecutive exons in the same transcript.

    Args:
        gtf_path: Path to GTF file (can be gzipped)

    Returns:
        List of (chrom, start, end, name, strand) tuples
        Coordinates are 0-based, half-open (BED format)
    """
    # Collect exons by transcript
    transcripts: Dict[str, List[Tuple[int, int, str, str]]] = {}

    opener = gzip.open if gtf_path.endswith('.gz') else open
    mode = 'rt' if gtf_path.endswith('.gz') else 'r'

    with opener(gtf_path, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs = fields

            if feature != 'exon':
                continue

            # Parse transcript_id from GTF attributes
            transcript_id = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    # GTF format: transcript_id "ENST00000..."
                    transcript_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    break

            if transcript_id:
                # Convert GTF 1-based to 0-based
                start_0based = int(start) - 1
                end_0based = int(end)

                if transcript_id not in transcripts:
                    transcripts[transcript_id] = []
                transcripts[transcript_id].append((start_0based, end_0based, chrom, strand))

    # Infer introns from exon gaps
    introns = []
    for transcript_id, exons in transcripts.items():
        if len(exons) < 2:
            continue

        # Sort by start position
        exons.sort(key=lambda x: x[0])
        chrom = exons[0][2]
        strand = exons[0][3]

        for i in range(len(exons) - 1):
            intron_start = exons[i][1]  # End of current exon
            intron_end = exons[i + 1][0]  # Start of next exon

            if intron_end > intron_start:  # Valid intron
                name = f"{transcript_id}_intron_{i + 1}"
                introns.append((chrom, intron_start, intron_end, name, strand))

    # Remove duplicates (same intron from multiple transcripts)
    unique_introns = list(set(introns))
    unique_introns.sort(key=lambda x: (x[0], x[1]))

    return unique_introns


def write_junction_bed(
    introns: List[Tuple[str, int, int, str, str]],
    output_path: str,
    compress: bool = False
):
    """Write introns to BED format for minimap2 --junc-bed.

    Args:
        introns: List of (chrom, start, end, name, strand) tuples
        output_path: Output BED file path
        compress: Whether to gzip the output
    """
    opener = gzip.open if compress else open
    mode = 'wt' if compress else 'w'

    with opener(output_path, mode) as f:
        for chrom, start, end, name, strand in introns:
            # 6-column BED format for minimap2
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")


def generate_junction_bed(
    annotation_path: str,
    output_path: Optional[str] = None,
    force: bool = False
) -> str:
    """Generate junction BED file from GFF or GTF annotation.

    If output_path is not specified, creates a .junc.bed file alongside
    the annotation file (or in a cache directory).

    Args:
        annotation_path: Path to GFF or GTF annotation file
        output_path: Optional output path for BED file
        force: Regenerate even if output exists

    Returns:
        Path to the generated BED file
    """
    annotation_path = Path(annotation_path)

    # Determine output path
    if output_path is None:
        output_path = annotation_path.with_suffix('.junc.bed')
    output_path = Path(output_path)

    # Check if we need to regenerate
    if output_path.exists() and not force:
        logger.info(f"Junction BED already exists: {output_path}")
        return str(output_path)

    # Determine file type and parse
    suffix = annotation_path.suffix.lower()
    if suffix == '.gz':
        # Check the second-to-last suffix
        inner_suffix = annotation_path.stem.split('.')[-1].lower()
    else:
        inner_suffix = suffix.lstrip('.')

    logger.info(f"Generating junction BED from {annotation_path}...")

    if inner_suffix in ('gff', 'gff3'):
        introns = parse_gff_introns(str(annotation_path))
    elif inner_suffix in ('gtf', 'gtf2'):
        introns = parse_gtf_introns(str(annotation_path))
    else:
        raise ValueError(f"Unknown annotation format: {inner_suffix}")

    logger.info(f"Extracted {len(introns)} introns")

    # Write BED file
    write_junction_bed(introns, str(output_path))
    logger.info(f"Written to {output_path}")

    return str(output_path)


def get_minimap2_junc_args(
    annotation_path: Optional[str] = None,
    junc_bed_path: Optional[str] = None,
    junc_bonus: int = 9,
    cache_dir: Optional[str] = None
) -> List[str]:
    """Get minimap2 arguments for junction annotation.

    Either annotation_path or junc_bed_path must be provided.

    Args:
        annotation_path: Path to GFF/GTF annotation file
        junc_bed_path: Path to pre-generated junction BED file
        junc_bonus: Bonus score for annotated junctions (default: 9)
        cache_dir: Directory to cache generated BED file

    Returns:
        List of minimap2 arguments: ['--junc-bed', path, '--junc-bonus', '9']
        Returns empty list if no annotation available
    """
    if junc_bed_path and Path(junc_bed_path).exists():
        return ['--junc-bed', junc_bed_path, '--junc-bonus', str(junc_bonus)]

    if annotation_path:
        # Generate or use cached BED file
        if cache_dir:
            output_path = Path(cache_dir) / 'annotation.junc.bed'
            Path(cache_dir).mkdir(parents=True, exist_ok=True)
        else:
            output_path = None

        try:
            bed_path = generate_junction_bed(annotation_path, output_path)
            return ['--junc-bed', bed_path, '--junc-bonus', str(junc_bonus)]
        except Exception as e:
            logger.warning(f"Could not generate junction BED: {e}")
            return []

    return []
