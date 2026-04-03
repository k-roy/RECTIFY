#!/usr/bin/env python3
"""
RECTIFY Batch Processing Command

Default mode: parallel interactive execution, auto-sizing to available CPUs.
SLURM mode: pass a profile YAML with cluster settings.

Usage:
    # Parallel interactive (default — uses all available cores)
    rectify batch --manifest manifest.tsv --genome genome.fa --annotation genes.gtf -o results/

    # SLURM array job via profile
    rectify batch --manifest manifest.tsv --genome genome.fa --annotation genes.gtf -o results/ \
        --profile slurm_profiles/sherlock_larsms.yaml

    # Or specify BAM files directly
    rectify batch --bams *.bam --genome genome.fa --annotation genes.gtf -o results/

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-25
"""

import os
import sys
import json
import subprocess
import concurrent.futures
import multiprocessing
from pathlib import Path
from typing import List, Optional, Dict, Tuple, Any
from datetime import datetime

from ..utils.provenance import init_provenance


SLURM_CORRECT_TEMPLATE = '''#!/bin/bash
#SBATCH --job-name={job_name}_correct
#SBATCH --partition={partition}
#SBATCH --time={time}
#SBATCH --mem={mem}
#SBATCH --cpus-per-task={cpus}
#SBATCH --array=0-{max_idx}{array_limit}
#SBATCH --output={log_dir}/rectify_correct_%A_%a.log

# =============================================================================
# RECTIFY Stage 1: Per-sample correction (SLURM array job)
# Generated: {timestamp}
# After all tasks complete, submit: rectify_batch_analyze.sh
# =============================================================================

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK

PYTHON="{python}"

# Parallel arrays: sample IDs and input paths
SAMPLE_IDS=(
{sample_ids_array}
)
INPUT_PATHS=(
{input_paths_array}
)

SAMPLE_ID="${{SAMPLE_IDS[$SLURM_ARRAY_TASK_ID]}}"
INPUT="${{INPUT_PATHS[$SLURM_ARRAY_TASK_ID]}}"

echo "=============================="
echo "RECTIFY Correction"
echo "=============================="
echo "Sample:   ${{SAMPLE_ID}}"
echo "Input:    ${{INPUT}}"
echo "Task ID:  $SLURM_ARRAY_TASK_ID"
echo "Date:     $(date)"
echo "Host:     $(hostname)"
echo ""

OAK_OUTPUT_DIR="{output_dir}/${{SAMPLE_ID}}"
mkdir -p "$OAK_OUTPUT_DIR"

{scratch_setup}

$PYTHON -m rectify correct \\
    "${{WORK_INPUT}}" \\
    --genome "{genome}" \\
    --annotation "{annotation}" \\
    -o "${{WORK_OUTPUT_DIR}}/corrected_3ends.tsv" \\
    --threads $SLURM_CPUS_PER_TASK \\
    {extra_args}

{scratch_teardown}

echo ""
echo "Complete: ${{SAMPLE_ID}}"
echo "Finished: $(date)"
'''

# Scratch staging block — inserted when use_scratch=True in the profile.
# $SCRATCH on Sherlock has ~75 GB/s vs Oak NFS; staging BAMs here reduces
# wall time and eliminates I/O contention across concurrent array tasks.
# BAMs are staged (copied) because they are accessed non-sequentially during
# correction; FASTQ inputs are left on Oak (read once, sequentially).
_SCRATCH_SETUP_BLOCK = '''\
# --- Scratch staging (use_scratch=true in profile) ---
# Stage BAM to $SCRATCH for high-bandwidth local I/O.
SCRATCH_DIR="$SCRATCH/rectify_${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"
mkdir -p "$SCRATCH_DIR"
trap 'rm -rf "$SCRATCH_DIR"' EXIT

if [[ "${{INPUT}}" == *.bam ]]; then
    echo "Staging BAM to scratch: ${{SCRATCH_DIR}}"
    cp "${{INPUT}}" "${{SCRATCH_DIR}}/"
    cp "${{INPUT}}.bai" "${{SCRATCH_DIR}}/" 2>/dev/null || true
    WORK_INPUT="${{SCRATCH_DIR}}/$(basename "${{INPUT}}")"
else
    # FASTQ is read sequentially once — no staging benefit
    WORK_INPUT="${{INPUT}}"
fi
WORK_OUTPUT_DIR="${{SCRATCH_DIR}}"
'''

_SCRATCH_TEARDOWN_BLOCK = '''\
# Copy all outputs from scratch back to Oak (exclude BAMs — already on Oak)
echo "Copying outputs to Oak: $OAK_OUTPUT_DIR"
rsync -a --exclude="*.bam" --exclude="*.bai" "${{SCRATCH_DIR}}/" "$OAK_OUTPUT_DIR/"
echo "Outputs copied: $(date)"
'''

_NO_SCRATCH_SETUP_BLOCK = '''\
WORK_INPUT="${{INPUT}}"
WORK_OUTPUT_DIR="${{OAK_OUTPUT_DIR}}"
'''

_NO_SCRATCH_TEARDOWN_BLOCK = ''  # outputs already in OAK_OUTPUT_DIR


SLURM_ANALYZE_TEMPLATE = '''#!/bin/bash
#SBATCH --job-name={job_name}_analyze
#SBATCH --partition={partition}
#SBATCH --time={analyze_time}
#SBATCH --mem={analyze_mem}
#SBATCH --cpus-per-task={analyze_cpus}
#SBATCH --output={log_dir}/rectify_analyze_%j.log

# =============================================================================
# RECTIFY Stage 2: Combine TSVs + combined analysis (runs after array job)
# Generated: {timestamp}
# =============================================================================

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK

PYTHON="{python}"

echo "=============================="
echo "RECTIFY Combined Analysis"
echo "=============================="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo ""

OUTPUT_DIR="{output_dir}"
COMBINED_DIR="$OUTPUT_DIR/combined"
mkdir -p "$COMBINED_DIR"

# ── Step 1: Add sample column and combine corrected TSVs ──────────────────
echo "Combining corrected TSVs..."
$PYTHON - <<'PYEOF'
import sys, pandas as pd
from pathlib import Path

samples = {samples_python_list}
output_dir = Path("{output_dir}")

dfs = []
for sample_id, input_path in samples:
    tsv = output_dir / sample_id / "corrected_3ends.tsv"
    if not tsv.exists():
        print(f"WARNING: {{tsv}} not found, skipping", file=sys.stderr)
        continue
    df = pd.read_csv(tsv, sep="\\t")
    df["sample"] = sample_id
    dfs.append(df)
    print(f"  {{sample_id}}: {{len(df):,}} reads")

if not dfs:
    print("ERROR: No corrected TSVs found", file=sys.stderr)
    sys.exit(1)

combined = pd.concat(dfs, ignore_index=True)
combined_tsv = output_dir / "combined" / "corrected_3ends_combined.tsv"
combined.to_csv(combined_tsv, sep="\\t", index=False)
print(f"Combined: {{len(dfs)}} samples, {{len(combined):,}} reads -> {{combined_tsv}}")
PYEOF

# ── Step 2: Run combined analysis ─────────────────────────────────────────
echo "Running combined analysis..."
$PYTHON -m rectify analyze \\
    "$COMBINED_DIR/corrected_3ends_combined.tsv" \\
    --annotation "{annotation}" \\
    --genome "{genome}" \\
    -o "$COMBINED_DIR" \\
    --threads $SLURM_CPUS_PER_TASK \\
    --run-deseq2 \\
    {extra_analyze_args}

echo ""
echo "Combined analysis complete: $COMBINED_DIR"
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


def validate_inputs(samples: List[Dict[str, str]], genome: Path, annotation: Optional[Path]) -> bool:
    """Validate that all input files exist."""
    valid = True

    # Check genome
    if genome and not genome.exists():
        print(f"ERROR: Genome file not found: {genome}", file=sys.stderr)
        valid = False

    # Check annotation (optional)
    if annotation and not annotation.exists():
        print(f"ERROR: Annotation file not found: {annotation}", file=sys.stderr)
        valid = False

    # Check BAM files
    for sample in samples:
        bam_path = Path(sample['bam_path'])
        if not bam_path.exists():
            print(f"ERROR: BAM file not found: {bam_path}", file=sys.stderr)
            valid = False

    return valid


def _get_available_cpus() -> int:
    """Get usable CPUs: respects SLURM allocation, falls back to system count."""
    slurm_cpus = os.environ.get('SLURM_CPUS_PER_TASK')
    if slurm_cpus:
        return int(slurm_cpus)
    return multiprocessing.cpu_count()


def _build_sample_cmd(sample: Dict[str, str], sample_output: Path, args) -> List[str]:
    """
    Build the rectify correct command for a single sample.

    Batch runs correct-only per sample (not run-all), then combines TSVs
    and runs a single shared analyze step for cross-sample DESeq2.
    """
    input_path = sample.get('path', sample.get('bam_path'))
    corrected_tsv = sample_output / 'corrected_3ends.tsv'

    cmd = [
        sys.executable, '-m', 'rectify', 'correct',
        str(input_path),
        '--genome', str(args.genome),
        '-o', str(corrected_tsv),
        '--threads', str(args.threads),
    ]
    if getattr(args, 'annotation', None):
        cmd.extend(['--annotation', str(args.annotation)])
    if getattr(args, 'organism', None):
        cmd.extend(['--organism', args.organism])
    if getattr(args, 'polya_sequenced', False):
        cmd.append('--polya-sequenced')
    if getattr(args, 'aligner', None):
        cmd.extend(['--aligner', args.aligner])
    if getattr(args, 'netseq_dir', None):
        cmd.extend(['--netseq-dir', str(args.netseq_dir)])
    if getattr(args, 'filter_spikein', None):
        cmd.extend(['--filter-spikein'] + args.filter_spikein)
    return cmd


def _run_sample_task(sample_id: str, cmd: List[str], sample_output: Path) -> Tuple[str, int, Path]:
    """
    Run one sample command, capturing output to a per-sample log.
    Returns (sample_id, returncode, log_path).
    Called from threads — no stdout interleaving.
    """
    log_file = sample_output / 'rectify_run.log'
    with open(log_file, 'w') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    return sample_id, result.returncode, log_file


def load_slurm_profile(profile_path: Path) -> Dict:
    """
    Load SLURM settings from a YAML or JSON profile file.

    YAML profile example (sherlock_larsms.yaml):
        partition: larsms,owners
        time: "4:00:00"    # MUST quote HH:MM:SS — bare colons parse as sexagesimal
        mem: 32G
        cpus: 8
        max_concurrent: 20
        job_name: rectify_batch
        submit: true          # auto-submit after generating script
    """
    suffix = profile_path.suffix.lower()
    with open(profile_path) as f:
        if suffix in ('.yaml', '.yml'):
            try:
                import yaml
                return yaml.safe_load(f) or {}
            except ImportError:
                raise ImportError(
                    "PyYAML is required for YAML profiles.\n"
                    "Install with:  pip install pyyaml\n"
                    "Or use a JSON profile file (.json) instead."
                )
        elif suffix == '.json':
            return json.load(f)
        else:
            raise ValueError(
                f"Profile must be .yaml, .yml, or .json (got '{suffix}')"
            )


def _apply_profile(args, profile: Dict) -> None:
    """Merge profile values into args, only where args has not been set by the user."""
    profile_to_arg = {
        'partition': 'partition',
        'time': 'time',
        'mem': 'mem',
        'cpus': 'cpus',
        'max_concurrent': 'max_concurrent',
        'job_name': 'job_name',
        'submit': 'submit',
    }
    defaults = {
        'partition': 'larsms,owners',
        'time': '4:00:00',
        'mem': '32G',
        'cpus': 8,
        'max_concurrent': None,
        'job_name': 'rectify_batch',
        'submit': False,
    }
    for profile_key, arg_key in profile_to_arg.items():
        if profile_key in profile:
            # Apply if user left arg at its argparse default
            current = getattr(args, arg_key, None)
            if current == defaults.get(arg_key) or current is None:
                setattr(args, arg_key, profile[profile_key])


def generate_slurm_scripts(
    samples: List[Dict[str, str]],
    output_dir: Path,
    genome: Path,
    annotation: Path,
    args,
) -> Tuple[str, str]:
    """
    Generate two SLURM scripts:
      1. rectify_batch_correct.sh  — array job, one task per sample
      2. rectify_batch_analyze.sh  — single job, combine + DESeq2

    Returns (correct_script, analyze_script) as strings.
    """
    # Parallel arrays of sample IDs and input paths
    sample_ids_array = '\n'.join(
        f'    "{s["sample_id"]}"' for s in samples
    )
    input_paths_array = '\n'.join(
        f'    "{s.get("path", s.get("bam_path", ""))}"' for s in samples
    )

    # Python list literal for the analyze script
    pairs = [(s['sample_id'], s.get('path', s.get('bam_path', ''))) for s in samples]
    samples_python_list = repr(pairs)

    # Scratch staging
    use_scratch = getattr(args, 'use_scratch', False)
    scratch_setup    = _SCRATCH_SETUP_BLOCK    if use_scratch else _NO_SCRATCH_SETUP_BLOCK
    scratch_teardown = _SCRATCH_TEARDOWN_BLOCK if use_scratch else _NO_SCRATCH_TEARDOWN_BLOCK

    # Correction extra args
    extra_args = []
    if getattr(args, 'organism', None):
        extra_args.append(f'--organism "{args.organism}"')
    if getattr(args, 'polya_sequenced', False):
        extra_args.append('--polya-sequenced')
    if getattr(args, 'aligner', None):
        extra_args.append(f'--aligner {args.aligner}')
    if getattr(args, 'netseq_dir', None):
        extra_args.append(f'--netseq-dir "{args.netseq_dir}"')
    if getattr(args, 'filter_spikein', None):
        extra_args.append('--filter-spikein ' + ' '.join(args.filter_spikein))
    if getattr(args, 'streaming', False):
        extra_args.append('--streaming')
    extra_args_str = ' \\\n    '.join(extra_args) if extra_args else ''

    # Analyze extra args (go annotations etc.)
    extra_analyze = []
    if getattr(args, 'go_annotations', None):
        extra_analyze.append(f'--go-annotations "{args.go_annotations}"')
    extra_analyze_str = ' \\\n    '.join(extra_analyze) if extra_analyze else ''

    max_concurrent = getattr(args, 'max_concurrent', None)
    array_limit = f'%{max_concurrent}' if max_concurrent else ''
    log_dir = output_dir / 'slurm_logs'
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    job_name = getattr(args, 'job_name', None) or 'rectify_batch'
    partition = getattr(args, 'partition', 'larsms,owners')

    correct_script = SLURM_CORRECT_TEMPLATE.format(
        job_name=job_name,
        partition=partition,
        time=getattr(args, 'time', '4:00:00'),
        mem=getattr(args, 'mem', '32G'),
        cpus=getattr(args, 'cpus', 8),
        max_idx=len(samples) - 1,
        array_limit=array_limit,
        log_dir=log_dir,
        timestamp=timestamp,
        python=sys.executable,
        sample_ids_array=sample_ids_array,
        input_paths_array=input_paths_array,
        output_dir=output_dir,
        genome=genome,
        annotation=annotation,
        extra_args=extra_args_str,
        scratch_setup=scratch_setup,
        scratch_teardown=scratch_teardown,
    )

    analyze_script = SLURM_ANALYZE_TEMPLATE.format(
        job_name=job_name,
        partition=partition,
        # Analyze may need more time/memory (DESeq2 is RAM-intensive)
        analyze_time=getattr(args, 'analyze_time', '8:00:00'),
        analyze_mem=getattr(args, 'analyze_mem', '64G'),
        analyze_cpus=getattr(args, 'analyze_cpus', getattr(args, 'cpus', 8)),
        log_dir=log_dir,
        timestamp=timestamp,
        python=sys.executable,
        samples_python_list=samples_python_list,
        output_dir=output_dir,
        annotation=annotation,
        genome=genome,
        extra_analyze_args=extra_analyze_str,
    )

    return correct_script, analyze_script


def generate_slurm_script(
    samples: List[Dict[str, str]],
    output_dir: Path,
    genome: Path,
    annotation: Path,
    args,
) -> str:
    """Legacy single-script generator — returns just the correction array script."""
    correct_script, _ = generate_slurm_scripts(samples, output_dir, genome, annotation, args)
    return correct_script


def _run_slurm_mode(
    samples: List[Dict[str, str]],
    args,
    profile: Optional[Dict] = None,
) -> int:
    """
    Generate two SLURM scripts with dependency chaining:
      Stage 1: rectify_batch_correct.sh  (array job, per-sample correction)
      Stage 2: rectify_batch_analyze.sh  (single job, combine + DESeq2)

    When --submit (or profile submit: true), submits both automatically:
      CORRECT_JOB=$(sbatch --parsable rectify_batch_correct.sh)
      sbatch --dependency=afterok:$CORRECT_JOB rectify_batch_analyze.sh
    """
    if profile:
        _apply_profile(args, profile)

    log_dir = args.output_dir / 'slurm_logs'
    log_dir.mkdir(parents=True, exist_ok=True)

    correct_script, analyze_script = generate_slurm_scripts(
        samples, args.output_dir, args.genome, args.annotation, args
    )

    correct_path = args.output_dir / 'rectify_batch_correct.sh'
    analyze_path = args.output_dir / 'rectify_batch_analyze.sh'

    with open(correct_path, 'w') as f:
        f.write(correct_script)
    os.chmod(correct_path, 0o755)

    with open(analyze_path, 'w') as f:
        f.write(analyze_script)
    os.chmod(analyze_path, 0o755)

    partition = getattr(args, 'partition', 'larsms,owners')
    print(f"SLURM scripts generated:")
    print(f"  Stage 1 (correction array):  {correct_path}")
    print(f"    Partition: {partition}  |  CPUs: {getattr(args, 'cpus', 8)}  "
          f"|  Mem: {getattr(args, 'mem', '32G')}  |  Time: {getattr(args, 'time', '4:00:00')}")
    print(f"    Array: 0-{len(samples) - 1}")
    print(f"  Stage 2 (combined analysis): {analyze_path}")
    print(f"    Partition: {partition}  |  CPUs: {getattr(args, 'analyze_cpus', 8)}  "
          f"|  Mem: {getattr(args, 'analyze_mem', '64G')}  "
          f"|  Time: {getattr(args, 'analyze_time', '8:00:00')}")

    should_submit = getattr(args, 'submit', False)
    if should_submit:
        print("\nSubmitting Stage 1 (correction array)...")
        result1 = subprocess.run(
            ['sbatch', '--parsable', str(correct_path)],
            capture_output=True, text=True,
        )
        if result1.returncode != 0:
            print(f"ERROR submitting correction job: {result1.stderr}", file=sys.stderr)
            return 1
        correct_job_id = result1.stdout.strip()
        print(f"  Correction array job submitted: {correct_job_id}")

        print(f"Submitting Stage 2 (combined analysis, depends on {correct_job_id})...")
        result2 = subprocess.run(
            ['sbatch', f'--dependency=afterok:{correct_job_id}', str(analyze_path)],
            capture_output=True, text=True,
        )
        if result2.returncode != 0:
            print(f"ERROR submitting analysis job: {result2.stderr}", file=sys.stderr)
            return 1
        analyze_job_id = result2.stdout.strip()
        print(f"  Analysis job submitted: {analyze_job_id}")
        print(f"\nJobs queued. Analysis will run automatically after all corrections complete.")
        return 0
    else:
        print("\nTo submit both stages with automatic dependency:")
        print(f"  CORRECT_JOB=$(sbatch --parsable {correct_path})")
        print(f"  sbatch --dependency=afterok:$CORRECT_JOB {analyze_path}")
        return 0


def _run_interactive_mode(samples: List[Dict[str, str]], args) -> int:
    """
    Three-stage batch pipeline (interactive parallel mode):

    Stage 1 (parallel):  rectify correct per sample
                         Workers = floor(available_cpus / threads_per_sample)
    Stage 2 (sequential): Combine corrected TSVs → add sample column
    Stage 3 (sequential): Combined rectify analyze (DESeq2, GO, motifs)
    """
    from .run_command import _combine_corrected_tsvs, _run_analysis

    n_cpus = _get_available_cpus()
    threads_per_sample = getattr(args, 'threads', 4)
    n_workers = getattr(args, 'workers', None) or max(1, n_cpus // threads_per_sample)

    print(f"[Stage 1/3] Correcting samples in parallel")
    print(f"  Available CPUs:  {n_cpus}")
    print(f"  Threads/sample:  {threads_per_sample}")
    print(f"  Workers:         {n_workers}")
    print()

    # Pre-build commands
    tasks = []
    for sample in samples:
        sample_output = args.output_dir / sample['sample_id']
        sample_output.mkdir(parents=True, exist_ok=True)
        cmd = _build_sample_cmd(sample, sample_output, args)
        tasks.append((sample['sample_id'], cmd, sample_output))

    failed = []
    completed = 0
    w = len(str(len(samples)))

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        future_to_id = {
            executor.submit(_run_sample_task, sid, cmd, out): sid
            for sid, cmd, out in tasks
        }
        for future in concurrent.futures.as_completed(future_to_id):
            sample_id = future_to_id[future]
            completed += 1
            try:
                sid, rc, log_file = future.result()
                status = "OK  " if rc == 0 else "FAIL"
                print(f"  [{completed:>{w}}/{len(samples)}] [{status}] {sid}")
                if rc != 0:
                    print(f"         Log: {log_file}", file=sys.stderr)
                    failed.append(sid)
                    if not getattr(args, 'continue_on_error', False):
                        print(
                            f"\nERROR: {sid} failed. Use --continue-on-error to keep going.",
                            file=sys.stderr,
                        )
                        return 1
            except Exception as e:
                print(f"  [{completed}/{len(samples)}] [ERROR] {sample_id}: {e}", file=sys.stderr)
                failed.append(sample_id)
                if not getattr(args, 'continue_on_error', False):
                    return 1

    successful = [s for s in samples if s['sample_id'] not in failed]
    print()
    if failed:
        print(f"  {len(failed)} failed: {', '.join(failed)}", file=sys.stderr)
    print(f"  {len(successful)}/{len(samples)} samples corrected.")

    if not successful:
        return 1

    # ── Stage 2: combine TSVs ────────────────────────────────────────────────
    print(f"\n[Stage 2/3] Combining corrected TSVs...")
    try:
        combined_tsv = _combine_corrected_tsvs(successful, args.output_dir)
    except Exception as e:
        print(f"ERROR combining TSVs: {e}", file=sys.stderr)
        return 1

    # ── Stage 3: combined analysis ───────────────────────────────────────────
    print(f"\n[Stage 3/3] Running combined analysis (DESeq2, GO, motifs)...")
    combined_dir = args.output_dir / 'combined'
    try:
        _run_analysis(
            corrected_tsv=combined_tsv,
            output_dir=combined_dir,
            genome_path=args.genome,
            annotation_path=args.annotation,
            args=args,
            n_samples=len(successful),
        )
    except Exception as e:
        print(f"ERROR in combined analysis: {e}", file=sys.stderr)
        return 1

    # ── Provenance ───────────────────────────────────────────────────────────
    provenance = init_provenance(
        args.output_dir,
        description=f"RECTIFY batch processing ({len(successful)} samples)",
        config=vars(args),
    )
    for sample in successful:
        sample_output = args.output_dir / sample['sample_id']
        input_path = Path(sample.get('path', sample.get('bam_path', '')))
        for tsv_file in sample_output.glob('*.tsv'):
            provenance.add_output_file(tsv_file, source_files=[input_path])
    provenance.save()

    print(f"\nBatch complete.")
    print(f"  Per-sample outputs:   {args.output_dir}/<sample_id>/")
    print(f"  Combined analysis:    {combined_dir}/")
    print(f"  DESeq2 results:       {combined_dir}/tables/deseq2_genes_*.tsv")
    print(f"  Batch provenance:     {args.output_dir / 'PROVENANCE.json'}")

    return 0


def _resolve_reference_paths(args) -> None:
    """Resolve genome/annotation/GO paths from explicit args or bundled organism data."""
    from rectify.data import ensure_reference_data, get_bundled_go_annotations_path, normalize_organism
    import sys

    organism = getattr(args, 'organism', None)
    genome_arg = getattr(args, 'genome', None)
    annotation_arg = getattr(args, 'annotation', None)

    if organism is None and genome_arg is None:
        print(
            "ERROR: No reference genome provided.\n"
            "  Supply --genome /path/to/genome.fa,\n"
            "  or use --Scer (S. cerevisiae bundled data),\n"
            "  or use --organism <name> for another supported organism.",
            file=sys.stderr,
        )
        sys.exit(1)

    if organism is not None:
        genome_path, annotation_path, data_source = ensure_reference_data(
            organism=organism,
            custom_genome=genome_arg,
            custom_annotation=annotation_arg,
            verbose=True,
        )
    else:
        genome_path = Path(genome_arg) if genome_arg else None
        annotation_path = Path(annotation_arg) if annotation_arg else None

    if genome_path is None:
        print(
            f"ERROR: No bundled genome available for organism '{organism}'. "
            "Use --genome to provide a custom reference.",
            file=sys.stderr,
        )
        sys.exit(1)

    args.genome = genome_path
    args.annotation = annotation_path

    if not getattr(args, 'go_annotations', None) and organism:
        go_path = get_bundled_go_annotations_path(normalize_organism(organism))
        if go_path:
            args.go_annotations = go_path


def run(args) -> int:
    """Run batch command."""

    _resolve_reference_paths(args)

    # Parse input
    if args.manifest:
        samples = parse_manifest(args.manifest)
    else:
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
        cond = f"  [{sample['condition']}]" if 'condition' in sample else ''
        print(f"  - {sample['sample_id']}: {sample['bam_path']}{cond}")
    print()

    # Validate inputs
    if not validate_inputs(samples, args.genome, args.annotation):
        return 1

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Determine execution mode
    profile = None
    if getattr(args, 'profile', None):
        try:
            profile = load_slurm_profile(args.profile)
            print(f"SLURM profile loaded: {args.profile}")
        except Exception as e:
            print(f"ERROR loading SLURM profile: {e}", file=sys.stderr)
            return 1

    if profile or getattr(args, 'slurm', False):
        return _run_slurm_mode(samples, args, profile)
    else:
        return _run_interactive_mode(samples, args)


def create_batch_parser(subparsers) -> None:
    """Add batch command to argument parser."""

    batch_parser = subparsers.add_parser(
        'batch',
        help='Process multiple samples (parallel interactive or SLURM array job)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Examples:
  # S. cerevisiae — no need to specify genome/annotation (bundled)
  rectify batch --manifest manifest.tsv --Scer -o results/

  # SLURM array job via profile file
  rectify batch --manifest manifest.tsv --Scer -o results/ \\
      --profile slurm_profiles/sherlock_larsms.yaml

  # With explicit references
  rectify batch --manifest manifest.tsv --genome genome.fa --annotation genes.gtf -o results/

  # Process BAM files directly (no manifest)
  rectify batch --bams sample1.bam sample2.bam sample3.bam --Scer -o results/

  # Spike-in filtering
  rectify batch --manifest manifest.tsv --Scer -o results/ --filter-spikein ENO2

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

    # Reference files (optional when --Scer or --organism is set)
    batch_parser.add_argument(
        '--genome',
        type=Path,
        default=None,
        help='Reference genome FASTA file. Not required when --Scer or --organism is set.'
    )

    batch_parser.add_argument(
        '--annotation',
        type=Path,
        default=None,
        help='Gene annotation file (GTF/GFF). Not required when --Scer or --organism is set.'
    )

    batch_parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory (subdirectories created per sample)'
    )

    # Execution mode
    mode_group = batch_parser.add_argument_group('Execution mode')
    mode_group.add_argument(
        '--profile',
        type=Path,
        metavar='PROFILE_YAML',
        help='SLURM profile YAML/JSON file. When provided, generates and optionally submits '
             'a SLURM array job. See rectify/slurm_profiles/ for examples.'
    )

    mode_group.add_argument(
        '--workers',
        type=int,
        default=None,
        help='Number of parallel workers for interactive mode '
             '(default: floor(available_cpus / --threads))'
    )

    # SLURM options (legacy / manual override)
    slurm_group = batch_parser.add_argument_group(
        'SLURM options (manual — prefer --profile)'
    )
    slurm_group.add_argument(
        '--slurm',
        action='store_true',
        help='Generate SLURM array job script (uses defaults below)'
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

    slurm_group.add_argument(
        '--analyze-time',
        default='8:00:00',
        help='SLURM time limit for combined analysis job (DESeq2 can be slow)'
    )

    slurm_group.add_argument(
        '--analyze-mem',
        default='64G',
        help='SLURM memory for combined analysis job'
    )

    slurm_group.add_argument(
        '--analyze-cpus',
        type=int,
        default=None,
        help='SLURM CPUs for analysis job (default: same as --cpus)'
    )

    # Rectify options
    rectify_group = batch_parser.add_argument_group('RECTIFY options')
    rectify_group.add_argument(
        '--organism',
        default=None,
        help='Organism name for bundled genome/annotation/NET-seq '
             '(e.g., yeast, saccharomyces_cerevisiae). '
             'Required when --genome and --annotation are not specified.'
    )
    rectify_group.add_argument(
        '--Scer',
        dest='organism',
        action='store_const',
        const='saccharomyces_cerevisiae',
        help='Shorthand for --organism saccharomyces_cerevisiae. '
             'Uses all bundled S. cerevisiae reference data.'
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
        help='Threads per sample'
    )

    rectify_group.add_argument(
        '--continue-on-error',
        action='store_true',
        help='Continue processing other samples if one fails'
    )

    rectify_group.add_argument(
        '--filter-spikein',
        nargs='+',
        metavar='GENE',
        help='Remove spike-in reads by gene name before processing (e.g., --filter-spikein ENO2)'
    )


# For direct import
import argparse
