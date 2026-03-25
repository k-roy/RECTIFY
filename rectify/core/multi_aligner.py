"""
Multi-aligner alignment pipeline for RECTIFY.

Runs multiple aligners (minimap2, mapPacBio, gapmm2) on the same reads
and provides utilities for consensus analysis.

Strategy:
- minimap2: Fast seed-and-chain baseline with junction annotation support
- mapPacBio: BBTools long-read aligner with splice-aware mode
- gapmm2: minimap2 wrapper with terminal exon refinement

Junction annotations are used to IMPROVE alignment quality but scoring
remains BLIND to annotations (novel junctions can still be detected).

Author: Kevin R. Roy
"""

import subprocess
import logging
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field

from ..utils.junction_bed import get_minimap2_junc_args, generate_junction_bed

logger = logging.getLogger(__name__)


@dataclass
class AlignerConfig:
    """Configuration for a single aligner."""
    name: str
    enabled: bool = True
    path: str = ""  # Path to executable (empty = use PATH)
    threads: int = 8
    extra_args: List[str] = field(default_factory=list)


@dataclass
class MultiAlignerConfig:
    """Configuration for multi-aligner pipeline."""
    minimap2: AlignerConfig = field(default_factory=lambda: AlignerConfig(
        name="minimap2",
        enabled=True,
        path="minimap2"
    ))
    map_pacbio: AlignerConfig = field(default_factory=lambda: AlignerConfig(
        name="mapPacBio",
        enabled=True,
        path="mapPacBio.sh"
    ))
    gapmm2: AlignerConfig = field(default_factory=lambda: AlignerConfig(
        name="gapmm2",
        enabled=True,
        path="gapmm2"
    ))

    # Junction annotation options
    use_junction_annotation: bool = True
    junc_bonus: int = 9

    # Output options
    keep_sam: bool = False


def check_aligner_available(aligner: str) -> bool:
    """Check if an aligner is available in PATH."""
    return shutil.which(aligner) is not None


def run_minimap2(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    annotation_path: Optional[str] = None,
    junc_bonus: int = 9,
    cache_dir: Optional[str] = None,
    extra_args: Optional[List[str]] = None
) -> str:
    """Run minimap2 alignment with optional junction annotation.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_bam: Path for output BAM file
        threads: Number of threads
        annotation_path: Optional GFF/GTF for junction annotation
        junc_bonus: Bonus score for annotated junctions
        cache_dir: Directory to cache junction BED
        extra_args: Additional minimap2 arguments

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # Base minimap2 command for direct RNA
    cmd = [
        'minimap2',
        '-ax', 'splice',
        '-uf',  # Stranded (forward)
        '-k14',  # Smaller k-mer for sensitivity
        '-G', '5000',  # Max intron size
        '--splice-flank=no',  # Disable for compatibility
        '--secondary=no',
        '-t', str(threads),
    ]

    # Add junction annotation if available
    if annotation_path:
        junc_args = get_minimap2_junc_args(
            annotation_path=annotation_path,
            junc_bonus=junc_bonus,
            cache_dir=cache_dir
        )
        cmd.extend(junc_args)

    # Add extra arguments
    if extra_args:
        cmd.extend(extra_args)

    # Add reference and reads
    cmd.extend([genome_path, reads_path])

    logger.info(f"Running minimap2: {' '.join(cmd[:10])}...")

    # Run minimap2, pipe to samtools
    sam_output = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    # Sort and convert to BAM
    sort_cmd = [
        'samtools', 'sort',
        '-@', str(threads),
        '-o', str(output_bam)
    ]

    sort_proc = subprocess.Popen(
        sort_cmd,
        stdin=sam_output.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    # Wait for completion
    sam_output.stdout.close()
    _, stderr = sort_proc.communicate()

    if sort_proc.returncode != 0:
        raise RuntimeError(f"minimap2/samtools failed: {stderr.decode()}")

    # Index BAM
    subprocess.run(['samtools', 'index', str(output_bam)], check=True)

    logger.info(f"minimap2 complete: {output_bam}")
    return str(output_bam)


def run_map_pacbio(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    extra_args: Optional[List[str]] = None
) -> str:
    """Run BBTools mapPacBio alignment.

    mapPacBio is splice-aware for long reads up to 6kb.
    intronlen=50 converts deletions >=50bp to N (intron) CIGAR operations.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_bam: Path for output BAM file
        threads: Number of threads
        extra_args: Additional mapPacBio arguments

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # Check if mapPacBio is available
    map_pacbio_path = shutil.which('mapPacBio.sh')
    if not map_pacbio_path:
        # Try common anaconda location
        map_pacbio_path = '/home/groups/larsms/users/kevinroy/anaconda3/bin/mapPacBio.sh'
        if not Path(map_pacbio_path).exists():
            raise FileNotFoundError("mapPacBio.sh not found")

    # mapPacBio outputs SAM to stdout
    sam_path = output_bam.with_suffix('.sam')

    cmd = [
        map_pacbio_path,
        f'ref={genome_path}',
        f'in={reads_path}',
        f'out={sam_path}',
        f'threads={threads}',
        'intronlen=50',  # Convert large deletions to introns
        'nodisk',
        'minratio=0.4',
    ]

    if extra_args:
        cmd.extend(extra_args)

    logger.info(f"Running mapPacBio: {' '.join(cmd[:5])}...")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"mapPacBio failed: {result.stderr}")

    # Convert to sorted BAM
    subprocess.run([
        'samtools', 'view', '-bS', str(sam_path),
        '|', 'samtools', 'sort', '-@', str(threads), '-o', str(output_bam)
    ], shell=True, check=True)

    # Actually use proper piping
    view_proc = subprocess.Popen(
        ['samtools', 'view', '-bS', str(sam_path)],
        stdout=subprocess.PIPE
    )
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-@', str(threads), '-o', str(output_bam)],
        stdin=view_proc.stdout
    )
    view_proc.stdout.close()
    sort_proc.communicate()

    # Index
    subprocess.run(['samtools', 'index', str(output_bam)], check=True)

    # Clean up SAM
    sam_path.unlink(missing_ok=True)

    logger.info(f"mapPacBio complete: {output_bam}")
    return str(output_bam)


def run_gapmm2(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    extra_args: Optional[List[str]] = None
) -> str:
    """Run gapmm2 alignment with terminal exon refinement.

    gapmm2 is a minimap2 wrapper that uses edlib to refine
    alignments at terminal exons, improving 5' and 3' end accuracy.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_bam: Path for output BAM file
        threads: Number of threads
        extra_args: Additional gapmm2 arguments

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # Check if gapmm2 is available
    gapmm2_path = shutil.which('gapmm2')
    if not gapmm2_path:
        gapmm2_path = '/home/groups/larsms/users/kevinroy/anaconda3/bin/gapmm2'
        if not Path(gapmm2_path).exists():
            raise FileNotFoundError("gapmm2 not found")

    sam_path = output_bam.with_suffix('.sam')

    # gapmm2 arguments
    cmd = [
        gapmm2_path,
        '-t', str(threads),
        '-m', '1',  # Mode 1: direct RNA
        '-i', '5000',  # Max intron size
        '-o', str(sam_path),
        genome_path,
        reads_path,
    ]

    if extra_args:
        cmd.extend(extra_args)

    logger.info(f"Running gapmm2: {' '.join(cmd[:5])}...")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"gapmm2 failed: {result.stderr}")

    # Convert to sorted BAM
    view_proc = subprocess.Popen(
        ['samtools', 'view', '-bS', str(sam_path)],
        stdout=subprocess.PIPE
    )
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-@', str(threads), '-o', str(output_bam)],
        stdin=view_proc.stdout
    )
    view_proc.stdout.close()
    sort_proc.communicate()

    # Index
    subprocess.run(['samtools', 'index', str(output_bam)], check=True)

    # Clean up SAM
    sam_path.unlink(missing_ok=True)

    logger.info(f"gapmm2 complete: {output_bam}")
    return str(output_bam)


def run_multi_aligner(
    reads_path: str,
    genome_path: str,
    output_dir: str,
    sample_name: str,
    annotation_path: Optional[str] = None,
    config: Optional[MultiAlignerConfig] = None,
    threads: int = 8,
    aligners: Optional[List[str]] = None
) -> Dict[str, str]:
    """Run multiple aligners on the same reads.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_dir: Output directory for BAM files
        sample_name: Sample name for output files
        annotation_path: Optional GFF/GTF for junction annotation
        config: Optional MultiAlignerConfig
        threads: Number of threads per aligner
        aligners: List of aligners to run (default: all enabled)

    Returns:
        Dict mapping aligner name to output BAM path
    """
    if config is None:
        config = MultiAlignerConfig()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = output_dir / '.cache'

    results = {}

    # Determine which aligners to run
    if aligners is None:
        aligners = []
        if config.minimap2.enabled:
            aligners.append('minimap2')
        if config.map_pacbio.enabled:
            aligners.append('mapPacBio')
        if config.gapmm2.enabled:
            aligners.append('gapmm2')

    for aligner in aligners:
        output_bam = output_dir / f"{sample_name}.{aligner}.sorted.bam"

        try:
            if aligner == 'minimap2':
                results['minimap2'] = run_minimap2(
                    reads_path=reads_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    annotation_path=annotation_path if config.use_junction_annotation else None,
                    junc_bonus=config.junc_bonus,
                    cache_dir=str(cache_dir),
                    extra_args=config.minimap2.extra_args
                )
            elif aligner == 'mapPacBio':
                results['mapPacBio'] = run_map_pacbio(
                    reads_path=reads_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    extra_args=config.map_pacbio.extra_args
                )
            elif aligner == 'gapmm2':
                results['gapmm2'] = run_gapmm2(
                    reads_path=reads_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    extra_args=config.gapmm2.extra_args
                )
            else:
                logger.warning(f"Unknown aligner: {aligner}")

        except Exception as e:
            logger.error(f"Aligner {aligner} failed: {e}")

    logger.info(f"Multi-aligner pipeline complete: {len(results)} aligners succeeded")
    return results
