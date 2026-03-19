#!/usr/bin/env python3
"""
Preprocessing module for RECTIFY.

Handles input data preparation including:
- Detecting input file format (FASTQ, BAM)
- Running minimap2 alignment for FASTQ files
- Decompressing gzipped reference files
- Index creation

Author: Kevin R. Roy
Date: 2026-03-18
"""

import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple


def detect_input_type(input_path: Path) -> str:
    """
    Detect input file type based on extension.

    Args:
        input_path: Path to input file

    Returns:
        'bam', 'fastq', 'fastq.gz', or 'unknown'
    """
    name = input_path.name.lower()

    if name.endswith('.bam'):
        return 'bam'
    elif name.endswith('.fastq.gz') or name.endswith('.fq.gz'):
        return 'fastq.gz'
    elif name.endswith('.fastq') or name.endswith('.fq'):
        return 'fastq'
    elif name.endswith('.sam'):
        return 'sam'
    else:
        return 'unknown'


def check_minimap2_installed() -> bool:
    """Check if minimap2 is available in PATH."""
    try:
        result = subprocess.run(
            ['minimap2', '--version'],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def check_samtools_installed() -> bool:
    """Check if samtools is available in PATH."""
    try:
        result = subprocess.run(
            ['samtools', '--version'],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def decompress_file(gzipped_path: Path, output_path: Path) -> Path:
    """
    Decompress a gzipped file.

    Args:
        gzipped_path: Path to gzipped file
        output_path: Path for decompressed output

    Returns:
        Path to decompressed file
    """
    with gzip.open(gzipped_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return output_path


def index_genome(genome_path: Path, verbose: bool = True) -> Path:
    """
    Create minimap2 and samtools indices for genome.

    Args:
        genome_path: Path to genome FASTA
        verbose: Print progress

    Returns:
        Path to genome (with indices created)
    """
    # Create samtools index (.fai) if not exists
    fai_path = Path(str(genome_path) + '.fai')
    if not fai_path.exists():
        if verbose:
            print(f"Creating samtools index: {fai_path}")
        subprocess.run(
            ['samtools', 'faidx', str(genome_path)],
            check=True
        )

    return genome_path


def run_minimap2_alignment(
    fastq_path: Path,
    genome_path: Path,
    output_bam: Path,
    preset: str = 'splice:hq',
    threads: int = 4,
    secondary: bool = False,
    extra_args: Optional[List[str]] = None,
    verbose: bool = True,
) -> Path:
    """
    Align FASTQ reads to genome using minimap2.

    Uses settings optimized for Nanopore direct RNA-seq.

    Args:
        fastq_path: Path to FASTQ file (can be gzipped)
        genome_path: Path to reference genome FASTA
        output_bam: Path for output BAM file
        preset: minimap2 preset ('splice:hq' for nanopore RNA, 'splice' for others)
        threads: Number of threads for minimap2
        secondary: Keep secondary alignments
        extra_args: Additional minimap2 arguments
        verbose: Print progress

    Returns:
        Path to output BAM file

    Raises:
        RuntimeError: If minimap2 or samtools not available
        subprocess.CalledProcessError: If alignment fails
    """
    if not check_minimap2_installed():
        raise RuntimeError(
            "minimap2 not found. Install with: conda install -c bioconda minimap2"
        )

    if not check_samtools_installed():
        raise RuntimeError(
            "samtools not found. Install with: conda install -c bioconda samtools"
        )

    # Build minimap2 command
    # Settings optimized for Nanopore direct RNA-seq:
    # -ax splice:hq = High-quality splice-aware alignment for Nanopore
    # -uf = Force use of f-strand (for direct RNA-seq, reads are reverse complement)
    # -k14 = k-mer size 14 (good for noisy long reads)
    # --secondary=no = Don't output secondary alignments
    # --MD = Include MD tag for SNV/indel calling
    cmd = [
        'minimap2',
        '-ax', preset,
        '-uf',  # Forward strand only (direct RNA is antisense)
        '-k14',
        '--MD',
        '-t', str(threads),
    ]

    if not secondary:
        cmd.append('--secondary=no')

    if extra_args:
        cmd.extend(extra_args)

    cmd.extend([str(genome_path), str(fastq_path)])

    if verbose:
        print(f"Running minimap2 alignment...")
        print(f"  Input: {fastq_path}")
        print(f"  Reference: {genome_path}")
        print(f"  Preset: {preset}")
        print(f"  Threads: {threads}")

    # Run minimap2, pipe to samtools for BAM conversion and sorting
    # minimap2 | samtools view -bh | samtools sort -o output.bam
    try:
        # Create temp SAM file
        temp_sam = output_bam.parent / f".{output_bam.stem}.temp.sam"

        # Run minimap2
        with open(temp_sam, 'w') as sam_out:
            proc_mm2 = subprocess.run(
                cmd,
                stdout=sam_out,
                stderr=subprocess.PIPE if not verbose else None,
                check=True
            )

        # Convert to BAM and sort
        if verbose:
            print("Converting to sorted BAM...")

        subprocess.run(
            ['samtools', 'view', '-bh', '-@', str(threads), str(temp_sam)],
            stdout=subprocess.PIPE,
            check=True
        )

        # Sort and output
        subprocess.run(
            [
                'samtools', 'sort',
                '-@', str(threads),
                '-o', str(output_bam),
                str(temp_sam)
            ],
            check=True
        )

        # Index the BAM
        subprocess.run(
            ['samtools', 'index', '-@', str(threads), str(output_bam)],
            check=True
        )

        # Clean up temp file
        temp_sam.unlink()

        if verbose:
            print(f"Output BAM: {output_bam}")

        return output_bam

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Alignment failed: {e}")


def prepare_input(
    input_path: Path,
    genome_path: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    threads: int = 4,
    verbose: bool = True,
) -> Tuple[Path, str]:
    """
    Prepare input for RECTIFY processing.

    Detects input type and runs alignment if needed.

    Args:
        input_path: Path to input file (BAM or FASTQ)
        genome_path: Path to reference genome (required for FASTQ)
        output_dir: Output directory for generated files
        threads: Number of threads for alignment
        verbose: Print progress

    Returns:
        Tuple of (bam_path, input_type)

    Raises:
        ValueError: If input type unknown or genome missing for FASTQ
    """
    input_type = detect_input_type(input_path)

    if input_type == 'unknown':
        raise ValueError(
            f"Unknown input file type: {input_path}\n"
            "Supported formats: .bam, .fastq, .fastq.gz, .fq, .fq.gz"
        )

    if input_type == 'bam':
        if verbose:
            print(f"Input is BAM file: {input_path}")
        return input_path, 'bam'

    # FASTQ input - need to align
    if input_type in ('fastq', 'fastq.gz'):
        if genome_path is None:
            raise ValueError(
                "Genome reference required for FASTQ input. "
                "Use --genome or specify organism with bundled genome."
            )

        if output_dir is None:
            output_dir = input_path.parent

        output_bam = output_dir / f"{input_path.stem.replace('.fastq', '').replace('.fq', '')}.bam"

        if verbose:
            print(f"Input is FASTQ file: {input_path}")
            print(f"Will align to genome and output to: {output_bam}")

        # Run alignment
        bam_path = run_minimap2_alignment(
            fastq_path=input_path,
            genome_path=genome_path,
            output_bam=output_bam,
            threads=threads,
            verbose=verbose,
        )

        return bam_path, 'fastq'

    raise ValueError(f"Unsupported input type: {input_type}")


def prepare_bundled_genome(
    organism: str,
    output_dir: Optional[Path] = None,
    verbose: bool = True,
) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Prepare bundled genome for use (decompress if needed).

    Args:
        organism: Organism name (e.g., 'yeast')
        output_dir: Directory to decompress to (default: temp)
        verbose: Print progress

    Returns:
        Tuple of (genome_path, annotation_path)
    """
    from ..data import (
        get_bundled_genome_path,
        get_bundled_annotation_path,
        normalize_organism,
    )

    org = normalize_organism(organism)
    genome_gz = get_bundled_genome_path(org)
    ann_gz = get_bundled_annotation_path(org)

    if genome_gz is None:
        return None, None

    # Create output directory
    if output_dir is None:
        output_dir = Path(tempfile.gettempdir()) / 'rectify_genomes' / org
    output_dir.mkdir(parents=True, exist_ok=True)

    # Decompress genome if needed
    genome_path = None
    if genome_gz:
        genome_name = genome_gz.name.replace('.gz', '')
        genome_path = output_dir / genome_name

        if not genome_path.exists():
            if verbose:
                print(f"Decompressing bundled genome for {org}...")
            decompress_file(genome_gz, genome_path)
            index_genome(genome_path, verbose=verbose)
        else:
            if verbose:
                print(f"Using cached genome: {genome_path}")

    # Decompress annotation if needed
    ann_path = None
    if ann_gz:
        ann_name = ann_gz.name.replace('.gz', '')
        ann_path = output_dir / ann_name

        if not ann_path.exists():
            if verbose:
                print(f"Decompressing bundled annotation...")
            decompress_file(ann_gz, ann_path)
        else:
            if verbose:
                print(f"Using cached annotation: {ann_path}")

    return genome_path, ann_path
