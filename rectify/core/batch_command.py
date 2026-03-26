#!/usr/bin/env python3
"""
RECTIFY Batch Processing Command

Enables processing multiple samples with optional SLURM array job submission.

Usage:
    # Process samples listed in manifest
    rectify batch manifest.tsv --genome genome.fa --annotation genes.gtf -o results/

    # Generate and submit SLURM array job
    rectify batch manifest.tsv --genome genome.fa --annotation genes.gtf -o results/ \
        --slurm --partition larsms,owners --time 4:00:00

    # Or specify BAM files directly
    rectify batch *.bam --genome genome.fa --annotation genes.gtf -o results/

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-25
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import List, Optional, Dict
from datetime import datetime


SLURM_TEMPLATE = '''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --time={time}
#SBATCH --mem={mem}
#SBATCH --cpus-per-task={cpus}
#SBATCH --array=0-{max_idx}{array_limit}
#SBATCH --output={log_dir}/rectify_%A_%a.log

# =============================================================================
# RECTIFY Batch Processing - SLURM Array Job
# Generated: {timestamp}
# =============================================================================

set -euo pipefail

# CRITICAL: Set thread limits to avoid oversubscription
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK

# Python environment
PYTHON="{python}"

# Sample list
SAMPLES=(
{samples_array}
)

# Get current sample
SAMPLE="${{SAMPLES[$SLURM_ARRAY_TASK_ID]}}"
SAMPLE_NAME=$(basename "$SAMPLE" .bam | sed 's/.sorted$//')

echo "=============================================="
echo "RECTIFY Batch Processing"
echo "=============================================="
echo "Sample: ${{SAMPLE_NAME}}"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo ""

# Create output directory for this sample
OUTPUT_DIR="{output_dir}/${{SAMPLE_NAME}}"
mkdir -p "$OUTPUT_DIR"

# Run rectify
$PYTHON -m rectify run \\
    "$SAMPLE" \\
    --genome "{genome}" \\
    --annotation "{annotation}" \\
    --output-dir "$OUTPUT_DIR" \\
    --threads $SLURM_CPUS_PER_TASK \\
    {extra_args}

echo ""
echo "Complete: ${{SAMPLE_NAME}}"
echo "Output: $OUTPUT_DIR"
echo "Finished: $(date)"
'''


def parse_manifest(manifest_path: Path) -> List[Dict[str, str]]:
    """
    Parse sample manifest TSV file.

    Expected format (tab-separated):
        sample_id    bam_path    [condition]
        wt_rep1      /path/to/wt_rep1.bam    WT
        ko_rep1      /path/to/ko_rep1.bam    KO

    If no header, assumes columns are: sample_id, bam_path
    """
    samples = []

    with open(manifest_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    if not lines:
        raise ValueError(f"Empty manifest file: {manifest_path}")

    # Check for header
    header = lines[0].lower().split('\t')
    has_header = 'bam' in header or 'bam_path' in header or 'path' in header

    if has_header:
        # Parse with header
        col_names = lines[0].split('\t')
        col_map = {name.lower(): idx for idx, name in enumerate(col_names)}

        # Find bam column
        bam_col = None
        for name in ['bam', 'bam_path', 'path', 'file']:
            if name in col_map:
                bam_col = col_map[name]
                break

        if bam_col is None:
            raise ValueError(f"Could not find BAM column in manifest. Found columns: {col_names}")

        # Find sample column
        sample_col = col_map.get('sample', col_map.get('sample_id', 0))

        # Find condition column
        condition_col = col_map.get('condition', col_map.get('group', None))

        for line in lines[1:]:
            parts = line.split('\t')
            sample = {
                'sample_id': parts[sample_col],
                'bam_path': parts[bam_col],
            }
            if condition_col is not None and condition_col < len(parts):
                sample['condition'] = parts[condition_col]
            samples.append(sample)
    else:
        # No header - assume sample_id, bam_path format
        for line in lines:
            parts = line.split('\t')
            if len(parts) >= 2:
                sample = {
                    'sample_id': parts[0],
                    'bam_path': parts[1],
                }
                if len(parts) >= 3:
                    sample['condition'] = parts[2]
                samples.append(sample)
            elif len(parts) == 1:
                # Just BAM path
                bam_path = parts[0]
                sample_id = Path(bam_path).stem.replace('.sorted', '')
                samples.append({
                    'sample_id': sample_id,
                    'bam_path': bam_path,
                })

    return samples


def validate_inputs(samples: List[Dict[str, str]], genome: Path, annotation: Path) -> bool:
    """Validate that all input files exist."""
    valid = True

    # Check genome
    if not genome.exists():
        print(f"ERROR: Genome file not found: {genome}", file=sys.stderr)
        valid = False

    # Check annotation
    if not annotation.exists():
        print(f"ERROR: Annotation file not found: {annotation}", file=sys.stderr)
        valid = False

    # Check BAM files
    for sample in samples:
        bam_path = Path(sample['bam_path'])
        if not bam_path.exists():
            print(f"ERROR: BAM file not found: {bam_path}", file=sys.stderr)
            valid = False

    return valid


def generate_slurm_script(
    samples: List[Dict[str, str]],
    output_dir: Path,
    genome: Path,
    annotation: Path,
    args,
) -> str:
    """Generate SLURM array job script."""

    # Build samples array
    bam_paths = [sample['bam_path'] for sample in samples]
    samples_array = '\n'.join(f'    "{bam}"' for bam in bam_paths)

    # Build extra args
    extra_args = []
    if args.organism:
        extra_args.append(f'--organism "{args.organism}"')
    if args.polya_sequenced:
        extra_args.append('--polya-sequenced')
    if args.aligner:
        extra_args.append(f'--aligner {args.aligner}')
    if args.netseq_dir:
        extra_args.append(f'--netseq-dir "{args.netseq_dir}"')

    extra_args_str = ' \\\n    '.join(extra_args) if extra_args else ''

    # Array limit
    array_limit = f'%{args.max_concurrent}' if args.max_concurrent else ''

    # Get Python path
    python_path = sys.executable

    # Log directory
    log_dir = output_dir / 'slurm_logs'

    script = SLURM_TEMPLATE.format(
        job_name=args.job_name or 'rectify_batch',
        partition=args.partition,
        time=args.time,
        mem=args.mem,
        cpus=args.cpus,
        max_idx=len(samples) - 1,
        array_limit=array_limit,
        log_dir=log_dir,
        timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        python=python_path,
        samples_array=samples_array,
        output_dir=output_dir,
        genome=genome,
        annotation=annotation,
        extra_args=extra_args_str,
    )

    return script


def run(args) -> int:
    """Run batch command."""

    # Parse input
    if args.manifest:
        # Parse manifest file
        samples = parse_manifest(args.manifest)
    else:
        # Direct BAM file list
        samples = []
        for bam_path in args.bams:
            sample_id = Path(bam_path).stem.replace('.sorted', '')
            samples.append({
                'sample_id': sample_id,
                'bam_path': str(bam_path),
            })

    if not samples:
        print("ERROR: No samples found.", file=sys.stderr)
        return 1

    print(f"Found {len(samples)} samples to process:")
    for sample in samples:
        print(f"  - {sample['sample_id']}: {sample['bam_path']}")
    print()

    # Validate inputs
    if not validate_inputs(samples, args.genome, args.annotation):
        return 1

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.slurm:
        # Generate and optionally submit SLURM array job
        log_dir = args.output_dir / 'slurm_logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        script = generate_slurm_script(samples, args.output_dir, args.genome, args.annotation, args)

        # Write script
        script_path = args.output_dir / 'rectify_batch.sh'
        with open(script_path, 'w') as f:
            f.write(script)
        os.chmod(script_path, 0o755)

        print(f"SLURM script written to: {script_path}")

        if args.submit:
            # Submit job
            print("Submitting SLURM array job...")
            result = subprocess.run(
                ['sbatch', str(script_path)],
                capture_output=True,
                text=True,
            )

            if result.returncode == 0:
                print(f"Job submitted: {result.stdout.strip()}")
                return 0
            else:
                print(f"ERROR submitting job: {result.stderr}", file=sys.stderr)
                return 1
        else:
            print("\nTo submit the job, run:")
            print(f"  sbatch {script_path}")
            return 0

    else:
        # Run locally (sequential)
        print("Running samples locally (sequential)...")
        print()

        for i, sample in enumerate(samples):
            print(f"[{i+1}/{len(samples)}] Processing {sample['sample_id']}...")

            sample_output = args.output_dir / sample['sample_id']
            sample_output.mkdir(parents=True, exist_ok=True)

            # Build command
            cmd = [
                sys.executable, '-m', 'rectify', 'run',
                sample['bam_path'],
                '--genome', str(args.genome),
                '--annotation', str(args.annotation),
                '--output-dir', str(sample_output),
                '--threads', str(args.threads),
            ]

            if args.organism:
                cmd.extend(['--organism', args.organism])
            if args.polya_sequenced:
                cmd.append('--polya-sequenced')
            if args.aligner:
                cmd.extend(['--aligner', args.aligner])
            if args.netseq_dir:
                cmd.extend(['--netseq-dir', str(args.netseq_dir)])

            # Run
            result = subprocess.run(cmd)

            if result.returncode != 0:
                print(f"ERROR processing {sample['sample_id']}", file=sys.stderr)
                if not args.continue_on_error:
                    return 1

        print()
        print(f"All {len(samples)} samples processed.")
        print(f"Results in: {args.output_dir}")
        return 0


def create_batch_parser(subparsers) -> None:
    """Add batch command to argument parser."""

    batch_parser = subparsers.add_parser(
        'batch',
        help='Process multiple samples with optional SLURM array job submission',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Examples:
  # Process samples from manifest
  rectify batch manifest.tsv --genome genome.fa --annotation genes.gtf -o results/

  # Generate SLURM script without submitting
  rectify batch manifest.tsv --genome genome.fa --annotation genes.gtf -o results/ --slurm

  # Generate and submit SLURM array job
  rectify batch manifest.tsv --genome genome.fa --annotation genes.gtf -o results/ \\
      --slurm --submit --partition larsms,owners --time 4:00:00

  # Process BAM files directly (no manifest)
  rectify batch sample1.bam sample2.bam sample3.bam \\
      --genome genome.fa --annotation genes.gtf -o results/

Manifest format (TSV):
  sample_id    bam_path           condition
  wt_rep1      /path/wt_rep1.bam  WT
  ko_rep1      /path/ko_rep1.bam  KO
        """
    )

    # Input (manifest or BAMs)
    input_group = batch_parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--manifest', '-m',
        type=Path,
        help='Sample manifest TSV file (sample_id, bam_path, [condition])'
    )
    input_group.add_argument(
        '--bams',
        nargs='+',
        type=Path,
        help='BAM files to process (alternative to manifest)'
    )

    # Required
    batch_parser.add_argument(
        '--genome',
        type=Path,
        required=True,
        help='Reference genome FASTA file'
    )

    batch_parser.add_argument(
        '--annotation',
        type=Path,
        required=True,
        help='Gene annotation file (GTF/GFF)'
    )

    batch_parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory (subdirectories created per sample)'
    )

    # SLURM options
    slurm_group = batch_parser.add_argument_group('SLURM options')
    slurm_group.add_argument(
        '--slurm',
        action='store_true',
        help='Generate SLURM array job script'
    )

    slurm_group.add_argument(
        '--submit',
        action='store_true',
        help='Submit SLURM job after generating script'
    )

    slurm_group.add_argument(
        '--partition',
        default='larsms,owners',
        help='SLURM partition'
    )

    slurm_group.add_argument(
        '--time',
        default='4:00:00',
        help='SLURM time limit per sample'
    )

    slurm_group.add_argument(
        '--mem',
        default='32G',
        help='SLURM memory per sample'
    )

    slurm_group.add_argument(
        '--cpus',
        type=int,
        default=8,
        help='SLURM CPUs per sample'
    )

    slurm_group.add_argument(
        '--max-concurrent',
        type=int,
        default=None,
        help='Maximum concurrent array jobs (default: no limit)'
    )

    slurm_group.add_argument(
        '--job-name',
        default='rectify_batch',
        help='SLURM job name'
    )

    # Rectify options
    rectify_group = batch_parser.add_argument_group('RECTIFY options')
    rectify_group.add_argument(
        '--organism',
        default='yeast',
        help='Organism name for bundled NET-seq'
    )

    rectify_group.add_argument(
        '--aligner',
        choices=['minimap2', 'star', 'bowtie2', 'bwa'],
        default='minimap2',
        help='Aligner used for BAM files'
    )

    rectify_group.add_argument(
        '--polya-sequenced',
        action='store_true',
        help='Poly(A) tail was sequenced'
    )

    rectify_group.add_argument(
        '--netseq-dir',
        type=Path,
        help='Custom NET-seq directory'
    )

    rectify_group.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Threads per sample (for local execution)'
    )

    rectify_group.add_argument(
        '--continue-on-error',
        action='store_true',
        help='Continue processing other samples if one fails (local mode)'
    )


# For direct import
import argparse
