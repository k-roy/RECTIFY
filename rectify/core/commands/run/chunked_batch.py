"""
Scheduler-array shell-script generator for ``rectify run-all --chunked-alignment``.

Single 900-line function (``_generate_chunked_pipeline``) that emits a
dependency-chained set of SLURM / UGE / PBS scripts: ``00_split.sh`` (FASTQ
split), per-aligner array jobs (mapPacBio + others), merge / consensus, and
``run_correct_analyze.sh``. Both single-sample (positional input) and
multi-sample (``--manifest``) paths are handled.

Self-contained — does not import from peer ``run`` modules. Heavy lifting
(scheduler header generation, alignment-array script emission, chunk-count
heuristics) is delegated to ``rectify.core.commands.split_command``.
"""

import stat
import sys
from pathlib import Path


def _generate_chunked_pipeline(args) -> int:
    """
    Generate a dependency-chained set of scheduler scripts for chunked alignment.

    For FASTQ input: generates split → alignment arrays → merge → correct+analyze chain.
    For BAM input: prints a warning and returns 0 (inline correction will proceed).
    For manifest: generates per-sample chunked alignment chains + correction array + analysis.

    Returns 0 on success, 1 on error.
    """
    from .. import split_command

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    scheduler    = getattr(args, 'scheduler', 'slurm')
    python_path  = getattr(args, 'python_path', None) or sys.executable
    rectify_src  = getattr(args, 'rectify_src', None)
    if rectify_src is None:
        # Auto-detect: walk up from this file to the package root
        rectify_src = str(Path(__file__).resolve().parent.parent.parent)

    target_reads = getattr(args, 'target_reads_per_chunk', split_command.DEFAULT_TARGET_READS)
    slurm_partition = getattr(args, 'slurm_partition', None)
    slurm_account   = getattr(args, 'slurm_account', None)
    uge_queue       = getattr(args, 'uge_queue', 'long.q')
    uge_pe          = getattr(args, 'uge_pe', 'smp')
    pbs_queue       = getattr(args, 'pbs_queue', 'workq')

    genome_path     = getattr(args, 'genome', None)
    annotation_path = getattr(args, 'annotation', None)

    sched_kwargs = dict(
        scheduler=scheduler,
        slurm_partition=slurm_partition,
        slurm_account=slurm_account,
        uge_queue=uge_queue,
        uge_pe=uge_pe,
        pbs_queue=pbs_queue,
    )

    def _genome_flag() -> str:
        return f'--genome "{genome_path}"' if genome_path else ''

    def _annot_flag() -> str:
        return f'--annotation "{annotation_path}"' if annotation_path else ''

    def _organism_flag() -> str:
        org = getattr(args, 'organism', None)
        return f'--organism "{org}"' if org else ''

    def _ref_flags() -> str:
        parts = [_genome_flag(), _annot_flag(), _organism_flag()]
        return ' '.join(p for p in parts if p)

    # ── Helper: correction+analysis script header (non-array) ────────────────
    def _make_correct_analyze_header(
        job_name: str, cores: int, mem_gb: int, time: str, log_dir: Path,
    ) -> str:
        log_pat = str(log_dir / '%j')
        return split_command._scheduler_headers(
            scheduler, job_name, 1, cores, mem_gb, time,
            str(log_dir), log_pat,
            is_array=False, **{k: v for k, v in sched_kwargs.items() if k != 'scheduler'},
        )

    # ── Helper: write an executable script ──────────────────────────────────
    def _write_script(path: Path, content: str) -> None:
        path.write_text(content)
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)

    # ── Helper: submission dependency syntax per scheduler ───────────────────
    def _submit_with_dep(script: Path, dep_var: str, dep_type: str = 'afterok') -> str:
        """Return a shell fragment that submits `script` after `dep_var` job IDs complete."""
        if scheduler == 'slurm':
            return (
                f'$(sbatch --dependency={dep_type}:{dep_var} {script} | awk \'{{print $4}}\')'
            )
        elif scheduler == 'uge':
            return f'$(qsub -hold_jid {dep_var} {script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')'
        else:  # pbs
            return f'$(qsub -W depend=afterok:{dep_var} {script})'

    def _submit_no_dep(script: Path) -> str:
        if scheduler == 'slurm':
            return f'$(sbatch {script} | awk \'{{print $4}}\')'
        elif scheduler == 'uge':
            return f'$(qsub {script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')'
        else:
            return f'$(qsub {script})'

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # Single-sample path
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    if getattr(args, 'input', None):
        input_path = Path(args.input)
        from ...align.preprocess import detect_input_type
        input_type = detect_input_type(input_path)

        # (BAM inputs are handled in run() before calling this function)
        if input_type not in ('fastq', 'fastq.gz'):
            # Shouldn't normally reach here; run() screens for this.
            return 0

        # Derive sample prefix from FASTQ name
        name = input_path.name
        sample_prefix = name
        for sfx in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
            if name.endswith(sfx):
                sample_prefix = name[:-len(sfx)]
                break

        # Count reads and compute n_chunks
        print(f"Counting reads in {input_path} ...")
        n_reads = split_command.count_reads(input_path)
        n_chunks = split_command.compute_n_chunks(n_reads, target_reads)
        print(f"  {n_reads:,} reads → {n_chunks} chunks (~{n_reads // n_chunks:,} reads/chunk)")

        chunks_dir = output_dir / 'chunks'
        chunks_dir.mkdir(parents=True, exist_ok=True)
        log_dir = output_dir / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        # ── 00_split.sh ────────────────────────────────────────────────────
        split_header = split_command._scheduler_headers(
            scheduler, f'{sample_prefix}_split', 1, 4, 16, '2:00:00',
            str(log_dir), str(log_dir / '%j'),
            is_array=False, **{k: v for k, v in sched_kwargs.items() if k != 'scheduler'},
        )
        split_script_content = f"""#!/bin/bash
# RECTIFY — split FASTQ into chunks for parallel alignment
# Generated by: rectify run-all --chunked-alignment

{split_header}

set -euo pipefail

export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

SPLIT_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-4}}}}}}
export OMP_NUM_THREADS=$SPLIT_CPUS
export OPENBLAS_NUM_THREADS=$SPLIT_CPUS
export MKL_NUM_THREADS=$SPLIT_CPUS
export LOKY_MAX_CPU_COUNT=$SPLIT_CPUS

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
INPUT="{input_path}"
CHUNKS_DIR="{chunks_dir}"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Start:   $(date)"
echo "Input:   $INPUT"
echo "Chunks:  $CHUNKS_DIR"

mkdir -p "$CHUNKS_DIR"
cd "$RECTIFY_SRC"

$PYTHON -m rectify split \\
    "$INPUT" \\
    -n {n_chunks} \\
    -o "$CHUNKS_DIR" \\
    --prefix "{sample_prefix}" \\
    --generate-slurm \\
    {_genome_flag()} \\
    {_annot_flag()} \\
    --python-path "{python_path}" \\
    --rectify-src "{rectify_src}" \\
    --scheduler {scheduler} \\
    --slurm-partition {slurm_partition} \\
    --slurm-account {slurm_account} \\
    --uge-queue {uge_queue} \\
    --uge-pe {uge_pe} \\
    --pbs-queue {pbs_queue}

echo "Split complete: $(date)"
ls -lh "$CHUNKS_DIR/"*chunk*.fastq.gz 2>/dev/null | head -5 || true
"""
        split_script = output_dir / '00_split.sh'
        _write_script(split_script, split_script_content)

        # ── Alignment scripts (in chunks_dir) ─────────────────────────────
        # Note: _generate_scripts writes into output_dir (= chunks_dir here)
        # Alignment scripts are generated now with n_chunks already known.
        align_scripts = split_command.generate_alignment_scripts(
            n_chunks=n_chunks,
            sample_prefix=sample_prefix,
            output_dir=chunks_dir,
            genome=genome_path,
            annotation=annotation_path,
            python_path=python_path,
            rectify_src=rectify_src,
            other_aligners=None,  # use defaults
            skip_map_pacbio=False,
            **{k: v for k, v in sched_kwargs.items() if k != 'scheduler'},
            scheduler=scheduler,
        )

        # ── run_correct_analyze.sh ─────────────────────────────────────────
        correct_header = _make_correct_analyze_header(
            f'{sample_prefix}_correct', 8, 32, '4:00:00', log_dir,
        )
        consensus_bam = output_dir / 'consensus' / f'{sample_prefix}.consensus.bam'
        ref_flags = _ref_flags()
        correct_analyze_content = f"""#!/bin/bash
# RECTIFY — correction + analysis after chunked alignment
# Generated by: rectify run-all --chunked-alignment
# Run after run_merge_consensus.sh completes.

{correct_header}

set -euo pipefail

export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

RECTIFY_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-8}}}}}}
{split_command._thread_limits_block()}

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
CONSENSUS_BAM="{consensus_bam}"
OUTDIR="{output_dir}"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Start:   $(date)"
echo "Input:   $CONSENSUS_BAM"

[ -f "$CONSENSUS_BAM" ] || {{ echo "ERROR: consensus BAM not found: $CONSENSUS_BAM" >&2; exit 1; }}

cd "$RECTIFY_SRC"
$PYTHON -m rectify run-all \\
    "$CONSENSUS_BAM" \\
    {ref_flags} \\
    --streaming \\
    --threads $RECTIFY_CPUS \\
    -o "$OUTDIR"

echo "Done: $(date)"
"""
        correct_analyze_script = output_dir / 'run_correct_analyze.sh'
        _write_script(correct_analyze_script, correct_analyze_content)

        # ── submit_pipeline.sh ─────────────────────────────────────────────
        mpb_script  = align_scripts['mpb']
        others_script = align_scripts['others']
        merge_script  = align_scripts['merge']

        if scheduler == 'slurm':
            submit_lines = [
                '#!/bin/bash',
                '# Submit the full chunked-alignment pipeline (SLURM)',
                '# Generated by: rectify run-all --chunked-alignment',
                '',
                f'SPLIT_JOB={_submit_no_dep(split_script)}',
                f'echo "Split job: $SPLIT_JOB"',
            ]
            if mpb_script:
                submit_lines += [
                    f'MPB_JOB=$(sbatch --dependency=afterok:$SPLIT_JOB {mpb_script} | awk \'{{print $4}}\')',
                    'echo "mapPacBio array: $MPB_JOB"',
                ]
            submit_lines += [
                f'OTHERS_JOB=$(sbatch --dependency=afterok:$SPLIT_JOB {others_script} | awk \'{{print $4}}\')',
                'echo "Others array: $OTHERS_JOB"',
            ]
            if mpb_script:
                merge_dep = '--dependency=afterok:$MPB_JOB:$OTHERS_JOB'
            else:
                merge_dep = '--dependency=afterok:$OTHERS_JOB'
            submit_lines += [
                f'MERGE_JOB=$(sbatch {merge_dep} {merge_script} | awk \'{{print $4}}\')',
                'echo "Merge+consensus: $MERGE_JOB"',
                f'CA_JOB=$(sbatch --dependency=afterok:$MERGE_JOB {correct_analyze_script} | awk \'{{print $4}}\')',
                'echo "Correct+analyze: $CA_JOB"',
                'echo ""',
                'echo "Pipeline submitted. Monitor with: squeue -u $USER"',
            ]
        elif scheduler == 'uge':
            submit_lines = [
                '#!/bin/bash',
                '# Submit the full chunked-alignment pipeline (UGE/SGE)',
                '',
                f'SPLIT_JOB={_submit_no_dep(split_script)}',
                f'echo "Split job: $SPLIT_JOB"',
            ]
            if mpb_script:
                submit_lines += [
                    f'MPB_JOB=$(qsub -hold_jid $SPLIT_JOB {mpb_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                    'echo "mapPacBio array: $MPB_JOB"',
                ]
            submit_lines += [
                f'OTHERS_JOB=$(qsub -hold_jid $SPLIT_JOB {others_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                'echo "Others array: $OTHERS_JOB"',
            ]
            if mpb_script:
                merge_hold = '-hold_jid $MPB_JOB,$OTHERS_JOB'
            else:
                merge_hold = '-hold_jid $OTHERS_JOB'
            submit_lines += [
                f'MERGE_JOB=$(qsub {merge_hold} {merge_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                'echo "Merge+consensus: $MERGE_JOB"',
                f'CA_JOB=$(qsub -hold_jid $MERGE_JOB {correct_analyze_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                'echo "Correct+analyze: $CA_JOB"',
                'echo ""',
                'echo "Pipeline submitted. Monitor with: qstat -u $USER"',
            ]
        else:  # pbs
            submit_lines = [
                '#!/bin/bash',
                '# Submit the full chunked-alignment pipeline (PBS/Torque)',
                '',
                f'SPLIT_JOB={_submit_no_dep(split_script)}',
                f'echo "Split job: $SPLIT_JOB"',
            ]
            if mpb_script:
                submit_lines += [
                    f'MPB_JOB=$(qsub -W depend=afterok:$SPLIT_JOB {mpb_script})',
                    'echo "mapPacBio array: $MPB_JOB"',
                ]
            submit_lines += [
                f'OTHERS_JOB=$(qsub -W depend=afterok:$SPLIT_JOB {others_script})',
                'echo "Others array: $OTHERS_JOB"',
            ]
            if mpb_script:
                merge_dep = '-W depend=afterok:$MPB_JOB:$OTHERS_JOB'
            else:
                merge_dep = '-W depend=afterok:$OTHERS_JOB'
            submit_lines += [
                f'MERGE_JOB=$(qsub {merge_dep} {merge_script})',
                'echo "Merge+consensus: $MERGE_JOB"',
                f'CA_JOB=$(qsub -W depend=afterok:$MERGE_JOB {correct_analyze_script})',
                'echo "Correct+analyze: $CA_JOB"',
                'echo ""',
                'echo "Pipeline submitted. Monitor with: qstat -u $USER"',
            ]

        submit_pipeline = output_dir / 'submit_pipeline.sh'
        _write_script(submit_pipeline, '\n'.join(submit_lines) + '\n')

        print("\nGenerated chunked-alignment pipeline scripts:")
        print(f"  Split:           {split_script}")
        print(f"  mapPacBio array: {align_scripts['mpb'] or '(skipped)'}")
        print(f"  Others array:    {align_scripts['others']}")
        print(f"  Merge+consensus: {align_scripts['merge']}")
        print(f"  Correct+analyze: {correct_analyze_script}")
        print(f"\nTo launch the pipeline:")
        print(f"  bash {submit_pipeline}")
        return 0

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # Multi-sample path (manifest)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    if getattr(args, 'manifest', None):
        from ..batch_command import parse_manifest
        manifest_path = Path(args.manifest)
        samples = parse_manifest(manifest_path)
        if not samples:
            print("ERROR: No samples found in manifest.", file=sys.stderr)
            return 1

        ref_flags = _ref_flags()
        log_dir = output_dir / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        fastq_samples = []
        bam_samples   = []

        for s in samples:
            p = Path(s.get('path', s.get('bam_path', '')))
            from ...align.preprocess import detect_input_type
            itype = detect_input_type(p)
            if itype in ('fastq', 'fastq.gz'):
                fastq_samples.append(s)
            else:
                bam_samples.append(s)

        # Per-sample split job IDs (for submit_all.sh)
        # We'll track them as shell variable names: SPLIT_<sample_id>
        sample_split_scripts    = {}
        sample_submit_pipelines = {}

        for s in fastq_samples:
            sample_id  = s['sample_id']
            input_path = Path(s.get('path', s.get('bam_path', '')))
            sample_dir = output_dir / sample_id
            sample_dir.mkdir(parents=True, exist_ok=True)
            chunks_dir = sample_dir / 'chunks'
            chunks_dir.mkdir(parents=True, exist_ok=True)
            slog_dir = sample_dir / 'logs'
            slog_dir.mkdir(parents=True, exist_ok=True)

            # Derive prefix
            name = input_path.name
            prefix = name
            for sfx in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
                if name.endswith(sfx):
                    prefix = name[:-len(sfx)]
                    break

            # Count reads
            print(f"  [{sample_id}] Counting reads ...")
            n_reads  = split_command.count_reads(input_path)
            n_chunks = split_command.compute_n_chunks(n_reads, target_reads)
            print(f"  [{sample_id}] {n_reads:,} reads → {n_chunks} chunks")

            # 00_split.sh for this sample
            split_header = split_command._scheduler_headers(
                scheduler, f'{sample_id}_split', 1, 4, 16, '2:00:00',
                str(slog_dir), str(slog_dir / '%j'),
                is_array=False, **{k: v for k, v in sched_kwargs.items() if k != 'scheduler'},
            )
            split_content = f"""#!/bin/bash
# RECTIFY — split FASTQ for {sample_id}

{split_header}

set -euo pipefail
export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

SPLIT_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-4}}}}}}
export OMP_NUM_THREADS=$SPLIT_CPUS
export OPENBLAS_NUM_THREADS=$SPLIT_CPUS
export MKL_NUM_THREADS=$SPLIT_CPUS
export LOKY_MAX_CPU_COUNT=$SPLIT_CPUS

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"

mkdir -p "{chunks_dir}"
cd "$RECTIFY_SRC"

echo "Start: $(date)  Host: $(hostname)"

$PYTHON -m rectify split \\
    "{input_path}" \\
    -n {n_chunks} \\
    -o "{chunks_dir}" \\
    --prefix "{prefix}" \\
    --generate-slurm \\
    {_genome_flag()} \\
    {_annot_flag()} \\
    --python-path "{python_path}" \\
    --rectify-src "{rectify_src}" \\
    --scheduler {scheduler} \\
    --slurm-partition {slurm_partition} \\
    --slurm-account {slurm_account} \\
    --uge-queue {uge_queue} \\
    --uge-pe {uge_pe} \\
    --pbs-queue {pbs_queue}

echo "Done: $(date)"
"""
            split_script = sample_dir / '00_split.sh'
            _write_script(split_script, split_content)
            sample_split_scripts[sample_id] = split_script

            # Alignment + correction scripts in chunks_dir (correct-first pipeline)
            _jpt = getattr(args, 'junction_penalty_table', None) or ''
            _spt = getattr(args, 'str_penalty_table', None) or ''
            align_scripts = split_command.generate_alignment_scripts(
                n_chunks=n_chunks,
                sample_prefix=prefix,
                output_dir=chunks_dir,
                genome=genome_path,
                annotation=annotation_path,
                python_path=python_path,
                rectify_src=rectify_src,
                scheduler=scheduler,
                other_aligners=None,
                skip_map_pacbio=False,
                slurm_partition=slurm_partition,
                slurm_account=slurm_account,
                uge_queue=uge_queue,
                uge_pe=uge_pe,
                pbs_queue=pbs_queue,
                junction_penalty_table=str(_jpt) if _jpt else '',
                str_penalty_table=str(_spt) if _spt else '',
            )

            # Per-sample submit_sample_alignment.sh — correct-first pipeline:
            # split → align (mpb + others) → merge_aligners → prescan
            #   → correct_* (per aligner) → chunk_merge → final_merge
            # Echoes the final_merge job ID as the last line for submit_all.sh.
            corrected_bam = chunks_dir / 'final' / f'{prefix}.corrected_reads.tsv'
            mpb_script          = align_scripts['mpb']
            others_script       = align_scripts['others']
            merge_aligners_script = align_scripts['merge_aligners']
            prescan_script      = align_scripts['prescan']
            correct_scripts     = align_scripts['correct']   # dict {aligner: Path}
            chunk_merge_script  = align_scripts['chunk_merge']
            final_merge_script  = align_scripts['final_merge']

            if scheduler == 'slurm':
                sub_lines = [
                    '#!/bin/bash',
                    f'# Submit correct-first pipeline for {sample_id}',
                    f'# Order: split → align → merge_aligners → prescan → correct_* → chunk_merge → final_merge',
                    f'# Diagnostic output goes to stderr; only the final SLURM job ID is written to stdout.',
                    '',
                    f'SPLIT_JOB={_submit_no_dep(split_script)}',
                    f'echo "  [{sample_id}] Split: $SPLIT_JOB" >&2',
                ]
                if mpb_script:
                    sub_lines += [
                        f'MPB_JOB=$(sbatch --dependency=afterok:$SPLIT_JOB {mpb_script} | awk \'{{print $4}}\')',
                        f'echo "  [{sample_id}] mapPacBio: $MPB_JOB" >&2',
                    ]
                sub_lines += [
                    f'OTHERS_JOB=$(sbatch --dependency=afterok:$SPLIT_JOB {others_script} | awk \'{{print $4}}\')',
                    f'echo "  [{sample_id}] Others: $OTHERS_JOB" >&2',
                ]
                align_dep = '--dependency=afterok:$MPB_JOB:$OTHERS_JOB' if mpb_script else '--dependency=afterok:$OTHERS_JOB'
                sub_lines += [
                    f'MERGE_ALIGNERS_JOB=$(sbatch {align_dep} {merge_aligners_script} | awk \'{{print $4}}\')',
                    f'echo "  [{sample_id}] MergeAligners: $MERGE_ALIGNERS_JOB" >&2',
                    f'PRESCAN_JOB=$(sbatch --dependency=afterok:$MERGE_ALIGNERS_JOB {prescan_script} | awk \'{{print $4}}\')',
                    f'echo "  [{sample_id}] Prescan: $PRESCAN_JOB" >&2',
                ]
                correct_vars = []
                for aligner, cscript in correct_scripts.items():
                    var = f'CORRECT_{aligner.upper()}_JOB'
                    sub_lines += [
                        f'{var}=$(sbatch --dependency=afterok:$PRESCAN_JOB {cscript} | awk \'{{print $4}}\')',
                        f'echo "  [{sample_id}] Correct_{aligner}: ${{{var}}}" >&2',
                    ]
                    correct_vars.append(f'${{{var}}}')
                correct_dep = ':'.join(correct_vars)
                sub_lines += [
                    f'CHUNK_MERGE_JOB=$(sbatch --dependency=afterok:{correct_dep} {chunk_merge_script} | awk \'{{print $4}}\')',
                    f'echo "  [{sample_id}] ChunkMerge: $CHUNK_MERGE_JOB" >&2',
                    f'FINAL_MERGE_JOB=$(sbatch --dependency=afterok:$CHUNK_MERGE_JOB {final_merge_script} | awk \'{{print $4}}\')',
                    f'echo "  [{sample_id}] FinalMerge: $FINAL_MERGE_JOB" >&2',
                    'echo $FINAL_MERGE_JOB',  # only stdout: bare job ID captured by submit_all.sh
                ]
            elif scheduler == 'uge':
                sub_lines = [
                    '#!/bin/bash',
                    f'# Submit correct-first pipeline for {sample_id}',
                    '',
                    f'SPLIT_JOB={_submit_no_dep(split_script)}',
                ]
                if mpb_script:
                    sub_lines += [f'MPB_JOB=$(qsub -hold_jid $SPLIT_JOB {mpb_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')']
                sub_lines += [f'OTHERS_JOB=$(qsub -hold_jid $SPLIT_JOB {others_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')']
                align_hold = '$MPB_JOB,$OTHERS_JOB' if mpb_script else '$OTHERS_JOB'
                sub_lines += [
                    f'MERGE_ALIGNERS_JOB=$(qsub -hold_jid {align_hold} {merge_aligners_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                    f'PRESCAN_JOB=$(qsub -hold_jid $MERGE_ALIGNERS_JOB {prescan_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                ]
                correct_vars = []
                for aligner, cscript in correct_scripts.items():
                    var = f'CORRECT_{aligner.upper()}_JOB'
                    sub_lines += [f'{var}=$(qsub -hold_jid $PRESCAN_JOB {cscript} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')']
                    correct_vars.append(f'${{{var}}}')
                correct_hold = ','.join(correct_vars)
                sub_lines += [
                    f'CHUNK_MERGE_JOB=$(qsub -hold_jid {correct_hold} {chunk_merge_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                    f'FINAL_MERGE_JOB=$(qsub -hold_jid $CHUNK_MERGE_JOB {final_merge_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                    'echo $FINAL_MERGE_JOB',
                ]
            else:  # pbs
                sub_lines = [
                    '#!/bin/bash',
                    f'# Submit correct-first pipeline for {sample_id}',
                    '',
                    f'SPLIT_JOB={_submit_no_dep(split_script)}',
                ]
                if mpb_script:
                    sub_lines += [f'MPB_JOB=$(qsub -W depend=afterok:$SPLIT_JOB {mpb_script})']
                sub_lines += [f'OTHERS_JOB=$(qsub -W depend=afterok:$SPLIT_JOB {others_script})']
                align_dep_pbs = '$MPB_JOB:$OTHERS_JOB' if mpb_script else '$OTHERS_JOB'
                sub_lines += [
                    f'MERGE_ALIGNERS_JOB=$(qsub -W depend=afterok:{align_dep_pbs} {merge_aligners_script})',
                    f'PRESCAN_JOB=$(qsub -W depend=afterok:$MERGE_ALIGNERS_JOB {prescan_script})',
                ]
                correct_vars = []
                for aligner, cscript in correct_scripts.items():
                    var = f'CORRECT_{aligner.upper()}_JOB'
                    sub_lines += [f'{var}=$(qsub -W depend=afterok:$PRESCAN_JOB {cscript})']
                    correct_vars.append(f'${{{var}}}')
                correct_dep_pbs = ':'.join(correct_vars)
                sub_lines += [
                    f'CHUNK_MERGE_JOB=$(qsub -W depend=afterok:{correct_dep_pbs} {chunk_merge_script})',
                    f'FINAL_MERGE_JOB=$(qsub -W depend=afterok:$CHUNK_MERGE_JOB {final_merge_script})',
                    'echo $FINAL_MERGE_JOB',
                ]

            sample_submit = sample_dir / 'submit_sample_alignment.sh'
            _write_script(sample_submit, '\n'.join(sub_lines) + '\n')
            sample_submit_pipelines[sample_id] = (sample_submit, corrected_bam)

        # ── Correction array script (BAM samples only) ────────────────────
        # FASTQ samples already include correction in their per-sample correct-first pipeline.
        # BAM samples bypass alignment entirely and go straight to correction here.
        all_samples_for_correct = []
        for s in bam_samples:
            all_samples_for_correct.append(
                (s['sample_id'], Path(s.get('path', s.get('bam_path', ''))))
            )

        n_correct_tasks = len(all_samples_for_correct)
        correct_header = split_command._scheduler_headers(
            scheduler, 'rectify_correct', n_correct_tasks, 8, 32, '4:00:00',
            str(log_dir), str(log_dir / '%A_%a'),
            is_array=True, **{k: v for k, v in sched_kwargs.items() if k != 'scheduler'},
        )

        # Build arrays for bash
        sample_ids_arr   = ' '.join(f'"{sid}"' for sid, _ in all_samples_for_correct)
        consensus_bams_arr = ' '.join(f'"{bam}"' for _, bam in all_samples_for_correct)

        # Penalty table flags for generated rectify correct calls
        import shlex as _shlex
        _jpt = getattr(args, 'junction_penalty_table', None)
        _spt = getattr(args, 'str_penalty_table', None)
        _penalty_flags = ''
        if _jpt:
            _penalty_flags += f'    --junction-penalty-table {_shlex.quote(str(_jpt))} \\\\\n'
        if _spt:
            _penalty_flags += f'    --str-penalty-table {_shlex.quote(str(_spt))} \\\\\n'

        _softclip_correct_flag = (
            '    --write-softclipped-bam "$SAMPLE_OUTDIR/${SAMPLE_ID}.rectified_pA_tail_trimmed.bam" \\\\\n'
            if getattr(args, 'write_softclip_bam', False) else ''
        )

        correct_array_content = f"""#!/bin/bash
# RECTIFY — correction array (one task per sample)
# Generated by: rectify run-all --chunked-alignment --manifest

{correct_header}

set -euo pipefail
export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

# ── Scheduler-agnostic task ID ────────────────────────────────────────────
if   [ -n "${{SLURM_ARRAY_TASK_ID:-}}" ]; then
    RECTIFY_TASK_ID=$SLURM_ARRAY_TASK_ID
    RECTIFY_CPUS=${{SLURM_CPUS_PER_TASK:-8}}
elif [ -n "${{SGE_TASK_ID:-}}" ]; then
    RECTIFY_TASK_ID=$(( SGE_TASK_ID - 1 ))
    RECTIFY_CPUS=${{NSLOTS:-8}}
elif [ -n "${{PBS_ARRAY_INDEX:-}}" ]; then
    RECTIFY_TASK_ID=$(( PBS_ARRAY_INDEX - 1 ))
    RECTIFY_CPUS=${{PBS_NUM_PPN:-8}}
else
    echo "ERROR: no scheduler array task variable found" >&2; exit 1
fi
# ─────────────────────────────────────────────────────────────────────────

{split_command._thread_limits_block()}

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"

SAMPLE_IDS=({sample_ids_arr})
CONSENSUS_BAMS=({consensus_bams_arr})

SAMPLE_ID="${{SAMPLE_IDS[$RECTIFY_TASK_ID]}}"
BAM="${{CONSENSUS_BAMS[$RECTIFY_TASK_ID]}}"
SAMPLE_OUTDIR="$OUTDIR/$SAMPLE_ID"
mkdir -p "$SAMPLE_OUTDIR"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Task:    $RECTIFY_TASK_ID  Sample: $SAMPLE_ID  CPUs: $RECTIFY_CPUS"
echo "Input:   $BAM"
echo "Start:   $(date)"

[ -f "$BAM" ] || {{ echo "ERROR: BAM not found: $BAM" >&2; exit 1; }}

cd "$RECTIFY_SRC"
$PYTHON -m rectify correct \\
    "$BAM" \\
    {_genome_flag()} \\
    {_annot_flag()} \\
{_penalty_flags}{_softclip_correct_flag}    -o "$SAMPLE_OUTDIR/corrected_reads.tsv" \\
    --streaming \\
    --threads $RECTIFY_CPUS

echo "Done: $(date)"
"""
        correct_array_script = output_dir / 'run_correct_array.sh'
        _write_script(correct_array_script, correct_array_content)

        # ── Analysis script ────────────────────────────────────────────────
        analyze_header = _make_correct_analyze_header(
            'rectify_analyze', 8, 64, '8:00:00', log_dir,
        )
        analyze_manifest_path = output_dir / 'corrected_manifest.tsv'
        analyze_content = f"""#!/bin/bash
# RECTIFY — combined analysis (DESeq2, GO, motifs)
# Generated by: rectify run-all --chunked-alignment --manifest
# Run after run_correct_array.sh completes.

{analyze_header}

set -euo pipefail
export PATH="$HOME/bin:$HOME/.rectify/bin:$PATH"

ANALYZE_CPUS=${{SLURM_CPUS_PER_TASK:-${{NSLOTS:-${{PBS_NUM_PPN:-8}}}}}}
{split_command._thread_limits_block('$ANALYZE_CPUS')}

PYTHON="{python_path}"
RECTIFY_SRC="{rectify_src}"
OUTDIR="{output_dir}"
MANIFEST="{analyze_manifest_path}"
COMBINED_DIR="$OUTDIR/combined"

echo "Python:  $($PYTHON --version)"
echo "Host:    $(hostname)"
echo "Start:   $(date)"

mkdir -p "$COMBINED_DIR"

# Write corrected manifest (sample_id, path, condition columns from original manifest)
# This manifest is pre-written by the script-generation step; just verify it exists.
[ -f "$MANIFEST" ] || {{ echo "ERROR: corrected manifest not found: $MANIFEST" >&2; exit 1; }}

cd "$RECTIFY_SRC"
$PYTHON -m rectify analyze /dev/null \\
    --manifest "$MANIFEST" \\
    {_genome_flag()} \\
    {_annot_flag()} \\
    --threads $ANALYZE_CPUS \\
    -o "$COMBINED_DIR"

echo "Done: $(date)"
"""
        analyze_script = output_dir / 'run_analyze.sh'
        _write_script(analyze_script, analyze_content)

        # Write the corrected manifest (paths to per-sample corrected TSVs)
        import csv
        manifest_rows = []
        for s in samples:
            row = {
                'sample_id': s['sample_id'],
                'path': str(output_dir / s['sample_id'] / 'corrected_reads.tsv'),
            }
            if 'condition' in s:
                row['condition'] = s['condition']
            manifest_rows.append(row)

        fieldnames = ['sample_id', 'path']
        if any('condition' in r for r in manifest_rows):
            fieldnames.append('condition')
        with open(analyze_manifest_path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t',
                                    extrasaction='ignore')
            writer.writeheader()
            writer.writerows(manifest_rows)

        # ── submit_all.sh ──────────────────────────────────────────────────
        if scheduler == 'slurm':
            submit_all_lines = [
                '#!/bin/bash',
                '# Submit all per-sample correct-first pipelines, then analysis',
                '# Generated by: rectify run-all --chunked-alignment --manifest',
                '# Each per-sample script echoes its final SLURM job ID to stdout.',
                '',
                '# ── Per-sample pipelines (all submitted in parallel) ──────────────',
                'FINAL_JOBS=""',
            ]
            for sid, (sample_submit, _) in sample_submit_pipelines.items():
                var = f'FINAL_{sid.replace("-", "_").replace(".", "_")}'
                submit_all_lines += [
                    f'{var}=$(bash {sample_submit})',
                    f'echo "  [{sid}] final_merge job: ${{{var}}}"',
                    f'FINAL_JOBS="${{FINAL_JOBS}}:${{{var}}}"',
                ]
            # BAM-only samples still use a global correction array
            if bam_samples and not fastq_samples:
                submit_all_lines += [
                    '',
                    '# ── Correction array (BAM-only manifest, no alignment dependency) ─',
                    f'CORRECT_JOB=$(sbatch {correct_array_script} | awk \'{{print $4}}\')',
                    'echo "Correction array: $CORRECT_JOB"',
                    '',
                    '# ── Analysis (after correction) ──────────────────────────────────',
                    f'ANALYZE_JOB=$(sbatch --dependency=afterok:$CORRECT_JOB {analyze_script} | awk \'{{print $4}}\')',
                ]
            elif bam_samples:
                # Mixed FASTQ+BAM: wait for both per-sample finals and BAM correction
                submit_all_lines += [
                    '',
                    '# ── Correction array for BAM samples ─────────────────────────────',
                    'FINAL_JOBS="${FINAL_JOBS#:}"  # strip leading colon',
                    f'CORRECT_JOB=$(sbatch --dependency=afterok:$FINAL_JOBS {correct_array_script} | awk \'{{print $4}}\')',
                    'echo "Correction array: $CORRECT_JOB"',
                    '',
                    '# ── Analysis (after all jobs complete) ───────────────────────────',
                    f'ANALYZE_JOB=$(sbatch --dependency=afterok:$CORRECT_JOB {analyze_script} | awk \'{{print $4}}\')',
                ]
            else:
                # FASTQ-only: no global correction needed — analysis after per-sample final_merge
                submit_all_lines += [
                    '',
                    '# ── Analysis (after all per-sample pipelines complete) ────────────',
                    'FINAL_JOBS="${FINAL_JOBS#:}"  # strip leading colon',
                    f'ANALYZE_JOB=$(sbatch --dependency=afterok:$FINAL_JOBS {analyze_script} | awk \'{{print $4}}\')',
                ]
            submit_all_lines += [
                'echo "Analysis: $ANALYZE_JOB"',
                'echo ""',
                'echo "Full pipeline submitted. Monitor with: squeue -u $USER"',
            ]
        elif scheduler == 'uge':
            submit_all_lines = [
                '#!/bin/bash',
                '# Submit all per-sample correct-first pipelines, then analysis',
                '',
                'FINAL_JOBS=""',
            ]
            for sid, (sample_submit, _) in sample_submit_pipelines.items():
                var = f'FINAL_{sid.replace("-", "_").replace(".", "_")}'
                submit_all_lines += [
                    f'{var}=$(bash {sample_submit})',
                    f'FINAL_JOBS="${{FINAL_JOBS}},${{{var}}}"',
                ]
            if bam_samples and not fastq_samples:
                submit_all_lines += [
                    f'CORRECT_JOB=$(qsub {correct_array_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                    f'ANALYZE_JOB=$(qsub -hold_jid $CORRECT_JOB {analyze_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                ]
            elif bam_samples:
                submit_all_lines += [
                    'FINAL_JOBS="${FINAL_JOBS#,}"',
                    f'CORRECT_JOB=$(qsub -hold_jid $FINAL_JOBS {correct_array_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                    f'ANALYZE_JOB=$(qsub -hold_jid $CORRECT_JOB {analyze_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                ]
            else:
                submit_all_lines += [
                    'FINAL_JOBS="${FINAL_JOBS#,}"',
                    f'ANALYZE_JOB=$(qsub -hold_jid $FINAL_JOBS {analyze_script} | awk \'{{print $3}}\' | grep -oE \'^[0-9]+\')',
                ]
            submit_all_lines += [
                'echo "Analysis: $ANALYZE_JOB"',
                'echo "Pipeline submitted. Monitor with: qstat -u $USER"',
            ]
        else:  # pbs
            submit_all_lines = [
                '#!/bin/bash',
                '# Submit all per-sample correct-first pipelines, then analysis',
                '',
                'FINAL_JOBS=""',
            ]
            for sid, (sample_submit, _) in sample_submit_pipelines.items():
                var = f'FINAL_{sid.replace("-", "_").replace(".", "_")}'
                submit_all_lines += [
                    f'{var}=$(bash {sample_submit})',
                    f'FINAL_JOBS="${{FINAL_JOBS}}:${{{var}}}"',
                ]
            if bam_samples and not fastq_samples:
                submit_all_lines += [
                    f'CORRECT_JOB=$(qsub {correct_array_script})',
                    f'ANALYZE_JOB=$(qsub -W depend=afterok:$CORRECT_JOB {analyze_script})',
                ]
            elif bam_samples:
                submit_all_lines += [
                    'FINAL_JOBS="${FINAL_JOBS#:}"',
                    f'CORRECT_JOB=$(qsub -W depend=afterok:$FINAL_JOBS {correct_array_script})',
                    f'ANALYZE_JOB=$(qsub -W depend=afterok:$CORRECT_JOB {analyze_script})',
                ]
            else:
                submit_all_lines += [
                    'FINAL_JOBS="${FINAL_JOBS#:}"',
                    f'ANALYZE_JOB=$(qsub -W depend=afterok:$FINAL_JOBS {analyze_script})',
                ]
            submit_all_lines += [
                'echo "Analysis: $ANALYZE_JOB"',
                'echo "Pipeline submitted. Monitor with: qstat -u $USER"',
            ]

        submit_all = output_dir / 'submit_all.sh'
        _write_script(submit_all, '\n'.join(submit_all_lines) + '\n')

        print("\nGenerated multi-sample correct-first pipeline scripts:")
        for sid in [s['sample_id'] for s in fastq_samples]:
            print(f"  [{sid}] {sample_split_scripts[sid]}")
            print(f"  [{sid}] {sample_submit_pipelines[sid][0]}")
        if bam_samples:
            print(f"  {len(bam_samples)} BAM sample(s): direct to correction")
            print(f"  Correction array: {correct_array_script}")
        print(f"  Analysis:         {analyze_script}")
        print(f"\nTo launch the full pipeline:")
        print(f"  bash {submit_all}")
        return 0

    print("ERROR: --chunked-alignment requires either a FASTQ input or --manifest.", file=sys.stderr)
    return 1


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
