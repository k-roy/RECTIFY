"""
Manifest-driven multi-sample pipeline for ``rectify run-all --manifest``.

- ``_run_multi_sample``: stage 1 fans samples out over a ``ThreadPoolExecutor``
  (one call to ``_process_one_sample`` per sample), stage 2 writes a manifest
  of per-sample corrected TSVs into ``output_dir/combined/``, stage 3 runs
  ``_run_analysis_manifest`` for the combined DESeq2 + GO + motif step.
- ``_run_analysis_manifest``: the manifest-mode equivalent of
  ``stages._run_analysis`` — invokes ``analyze`` with ``manifest=...`` so the
  analyzer streams each per-sample TSV in turn (no combined-frame load).
"""

import argparse
import concurrent.futures
import sys
from pathlib import Path
from typing import Optional

from .helpers import _resolve_reference_paths
from .single_sample import _process_one_sample
from .stages import (
    _readable_corrected_tsv,
)


def _run_analysis_manifest(
    manifest_path: Path,
    output_dir: Path,
    genome_path: Optional[Path],
    annotation_path: Optional[Path],
    args,
    n_samples: int = 1,
) -> None:
    """Run the analyze command in manifest mode (memory-efficient multi-sample path)."""
    from ..analyze_command import run_analyze

    run_deseq2 = n_samples > 1

    analyze_args = argparse.Namespace(
        # In manifest mode, 'input' is unused but run_analyze checks it after
        # the manifest dispatch, so set it to the manifest path as a fallback.
        input=str(manifest_path),
        manifest=str(manifest_path),
        output=output_dir,
        annotation=annotation_path,
        genome=genome_path,
        reference=getattr(args, 'reference', None),
        go_annotations=getattr(args, 'go_annotations', None),
        threads=getattr(args, 'threads', 4),
        # Clustering
        sample_column='sample',
        count_column=None,
        cluster_distance=25,
        min_reads=5,
        # Analysis flags
        run_deseq2=run_deseq2,
        run_motif=run_deseq2,
        sample_sets=None,
        # Filtering
        exclude_mito=True,
        include_mito=False,
        exclude_rdna=True,
        include_rdna=False,
        # Bedgraph and genomic distribution now run per-condition in manifest mode
        no_bedgraph=False,
        bedgraph_dir=None,
        no_genomic_distribution=False,
        # TSS clustering window (75 bp for DRS — 5' ends are noisier than 3' ends)
        tss_cluster_distance=75,
        # Motif windows
        motif_upstream=100,
        motif_downstream=50,
    )

    exit_code = run_analyze(analyze_args)
    if exit_code != 0:
        print(f"\nAnalysis completed with warnings (exit code: {exit_code})")
def _run_multi_sample(args) -> int:
    """
    Multi-sample pipeline via manifest:
      Stage 1 (parallel): align + correct per sample
      Stage 2:            combine corrected TSVs → add sample column
      Stage 3:            combined analyze (full DESeq2, GO, motifs)
    """
    from ..batch_command import parse_manifest, _get_available_cpus

    _resolve_reference_paths(args)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_path = args.genome
    annotation_path = args.annotation

    samples = parse_manifest(Path(args.manifest))
    if not samples:
        print("ERROR: No samples found in manifest.", file=sys.stderr)
        return 1

    # Validate 'path' column exists (manifest uses 'bam_path' or 'path')
    for s in samples:
        if 'path' not in s and 'bam_path' not in s:
            print(
                f"ERROR: manifest row for '{s['sample_id']}' has no path column.",
                file=sys.stderr,
            )
            return 1

    print(f"Found {len(samples)} samples:")
    for s in samples:
        p = s.get('path', s.get('bam_path', '?'))
        cond = f"  [{s['condition']}]" if 'condition' in s else ''
        print(f"  {s['sample_id']}: {p}{cond}")
    print()

    # ── Stage 1: align + correct in parallel ────────────────────────────────
    n_cpus = _get_available_cpus()
    _t = getattr(args, 'threads', 0)
    threads_per_sample = _t if _t > 0 else 4  # 0 = auto-detect, default 4
    n_workers = max(1, n_cpus // threads_per_sample)

    print(f"[Stage 1/3] Aligning and correcting ({n_workers} parallel workers)...")

    failed = []
    completed = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        future_to_sample = {
            executor.submit(
                _process_one_sample, s, output_dir, genome_path, annotation_path, args
            ): s
            for s in samples
        }

        for future in concurrent.futures.as_completed(future_to_sample):
            s = future_to_sample[future]
            completed += 1
            try:
                sample_id, rc = future.result()
                status = "OK  " if rc == 0 else "FAIL"
                n = len(str(len(samples)))
                print(f"  [{completed:>{n}}/{len(samples)}] [{status}] {sample_id}")
                if rc != 0:
                    failed.append(sample_id)
                    if not getattr(args, 'continue_on_error', False):
                        print(
                            f"\nERROR: {sample_id} failed. "
                            f"Use --continue-on-error to proceed anyway.",
                            file=sys.stderr,
                        )
                        return 1
            except Exception as e:
                completed += 1
                print(f"  [ERROR] {s['sample_id']}: {e}", file=sys.stderr)
                failed.append(s['sample_id'])
                if not getattr(args, 'continue_on_error', False):
                    return 1

    if failed:
        print(f"\n{len(failed)} sample(s) failed: {', '.join(failed)}", file=sys.stderr)
        if len(failed) == len(samples):
            return 1

    # Only combine successfully corrected samples
    successful_samples = [s for s in samples if s['sample_id'] not in failed]

    # ── Stage 2: write manifest pointing to per-sample corrected TSVs ────────
    # No combine step needed — manifest mode in analyze streams each file separately
    print(f"\n[Stage 2/3] Writing sample manifest for combined analysis...")
    import pandas as pd
    combined_dir = output_dir / 'combined'
    combined_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    unresolved = []
    for s in successful_samples:
        # 🔴 This used to hardcode `corrected_reads.tsv`, which has not existed since Commit B made
        # the per-region manifest the default artifact (the merged TSV is renamed to
        # `corrected_reads.region_000.tsv`). Every path written here pointed at a missing file:
        # Pass 1 happened to survive on `corrected_reads_index.bed.gz`, then Pass 1b died with
        # FileNotFoundError. Resolve to something that is actually readable as a reads table.
        resolved = _readable_corrected_tsv(output_dir / s['sample_id'])
        if resolved is None:
            unresolved.append(s['sample_id'])
            continue
        row = {
            'sample_id': s['sample_id'],
            'path': str(resolved),
        }
        if 'condition' in s:
            row['condition'] = s.get('condition', '')
        manifest_rows.append(row)

    if unresolved:
        # Fail loudly rather than analysing a silently smaller cohort.
        print(f"ERROR: no readable corrected TSV for: {', '.join(unresolved)}. "
              f"Expected corrected_reads.tsv or a corrected_reads.manifest.tsv naming one region "
              f"TSV under {output_dir}/<sample_id>/.", file=sys.stderr)
        return 1

    manifest_for_analyze = combined_dir / 'corrected_manifest.tsv'
    pd.DataFrame(manifest_rows).to_csv(manifest_for_analyze, sep='\t', index=False)
    print(f"  Wrote manifest: {manifest_for_analyze} ({len(manifest_rows)} samples)")

    # ── Stage 3: combined analysis (manifest mode) ───────────────────────────
    print(f"\n[Stage 3/3] Running combined analysis (DESeq2, GO, motifs)...")
    try:
        _run_analysis_manifest(
            manifest_path=manifest_for_analyze,
            output_dir=combined_dir,
            genome_path=genome_path,
            annotation_path=annotation_path,
            args=args,
            n_samples=len(successful_samples),
        )
    except Exception as e:
        print(f"ERROR: Combined analysis failed: {e}", file=sys.stderr)
        return 1

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("Pipeline Complete!")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}")
    print("\nPer-sample outputs:")
    for s in successful_samples:
        print(f"  {output_dir}/{s['sample_id']}/corrected_reads.tsv")
    print(f"\nCombined analysis: {combined_dir}/")
    print(f"  DESeq2 results:  {combined_dir}/tables/deseq2_genes_*.tsv")
    print(f"  HTML report:     {combined_dir}/report.html")
    print(f"  Provenance:      {combined_dir}/PROVENANCE.json")

    return 0
