"""
Per-pipeline-stage runners for ``rectify run-all``.

Each function is a thin wrapper around one stage of the pipeline:

- ``_run_alignment``               — multi-aligner consensus (or re-use of an
                                     existing rectified BAM).
- ``_run_correction``              — single ``rectify correct`` invocation.
- ``_run_correction_per_aligner``  — per-aligner correction (one call per BAM)
                                     into ``output_dir/per_aligner_corrected/``.
- ``_combine_corrected_tsvs``      — concat per-sample corrected TSVs.
- ``_run_analysis``                — ``rectify analyze`` on a (combined) TSV.
- ``_run_junction_aggregation``    — junction aggregation with partial rescue
                                     (no-op when annotation has no introns).

The resume-on-restart logic (skip-if-output-exists) is sample-level and stays
in ``single_sample`` / ``multi_sample``; these runners just execute.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .helpers import (
    _bam_has_md_tags,
    _collect_per_aligner_bams,
    _rectified_bam_path,
    _validate_bam_integrity,
)


def _run_alignment(
    input_path: Path,
    sample_id: str,
    sample_output_dir: Path,
    genome_path: Path,
    annotation_path: Optional[Path],
    threads: int,
    parallel_aligners: bool = False,
    base_aligners: Optional[List[str]] = None,
    junction_aligners: Optional[List[str]] = None,
    chimeric_consensus: bool = True,
    ultra_path: str = 'uLTRA',
    desalt_path: str = 'deSALT',
    mapPacBio_chunks: int = 1,
    checkpoint_dir: Optional[str] = None,
    short_read: bool = False,
    trust_existing_bams: bool = False,
) -> Tuple[Dict[str, Path], Path]:
    """
    Run multi-aligner alignment and selection, or return existing rectified.bam.

    Default aligners: minimap2 + mapPacBio + gapmm2 (long-read Tier 1).
    Pass base_aligners to restrict or change the set (e.g. ['mapPacBio']).
    Pass short_read=True to use bbmap + bwa instead of the long-read panel.
    Pass junction_aligners=[] to disable uLTRA + deSALT.

    Skips automatically if the rectified.bam already exists — safe to re-run.

    Returns
    -------
    Tuple of (per_aligner_bams, rectified_bam) where per_aligner_bams maps
    aligner name → BAM path for each per-aligner BAM found on disk.
    """
    rectified_bam = _rectified_bam_path(sample_id, sample_output_dir)

    # Backward compatibility: accept old consensus.bam filename from prior runs
    _legacy_bam = sample_output_dir / f"{sample_id}.consensus.bam"
    if not rectified_bam.exists() and _legacy_bam.exists():
        rectified_bam = _legacy_bam

    # Compute run provenance once — used by both the rectified.bam reuse gate
    # and the per-aligner BAM provenance filter.
    import sys as _sys
    from rectify.utils.bam_provenance import (
        compute_run_provenance,
        expected_provenance_for_aligner,
        matches_strict,
        read_sidecar,
    )
    _run_provenance = compute_run_provenance(command=_sys.argv)

    if _validate_bam_integrity(rectified_bam):
        if trust_existing_bams:
            print(f"    Skipping alignment — rectified.bam exists (--trust-existing-bams): {rectified_bam}")
            per_aligner_bams = _collect_per_aligner_bams(
                sample_id, sample_output_dir,
                run_provenance=None, trust_existing_bams=True,
            )
            return per_aligner_bams, rectified_bam
        _expected = expected_provenance_for_aligner(_run_provenance, "consensus")
        _stored = read_sidecar(rectified_bam)
        _ok, _reason = matches_strict(_stored, _expected)
        if _ok:
            print(f"    Skipping alignment — rectified.bam exists (provenance match): {rectified_bam}")
            per_aligner_bams = _collect_per_aligner_bams(
                sample_id, sample_output_dir,
                run_provenance=_run_provenance, trust_existing_bams=False,
            )
            return per_aligner_bams, rectified_bam
        print(
            f"    rectified.bam provenance mismatch ({_reason}); re-running alignment. "
            "Pass --trust-existing-bams to override."
        )

    if base_aligners is not None:
        _base_aligners = base_aligners
    elif short_read:
        _base_aligners = ['bbmap', 'bwa']
    else:
        _base_aligners = ['minimap2', 'mapPacBio', 'gapmm2']
    # Junction-aligner default depends on --short-read: BBMap's intronlen=20
    # already covers splicing for short reads, so uLTRA/deSALT aren't useful
    # and shouldn't be the default. Long-read protocols still default to
    # [uLTRA, deSALT]. Explicit user list (including []) overrides.
    if junction_aligners is None:
        _junction_aligners = [] if short_read else ['uLTRA', 'deSALT']
    else:
        _junction_aligners = junction_aligners
    all_aligners = _base_aligners + _junction_aligners
    aligner_desc = ' + '.join(all_aligners)
    print(f"    Running {len(all_aligners)}-aligner consensus ({aligner_desc})...")
    from ..align_command import run_align

    align_args = argparse.Namespace(
        reads=input_path,
        genome=genome_path,
        output_dir=sample_output_dir,
        annotation=annotation_path,
        threads=threads,
        aligners=_base_aligners,
        short_read=short_read,
        junction_aligners=_junction_aligners,
        no_consensus=False,
        chimeric_consensus=chimeric_consensus,
        junc_bonus=9,
        junc_bed=None,
        parallel_aligners=parallel_aligners,
        minimap2_path='minimap2',
        mapPacBio_path='mapPacBio.sh',
        gapmm2_path='gapmm2',
        ultra_path=ultra_path,
        desalt_path=desalt_path,
        mapPacBio_chunks=mapPacBio_chunks,
        mapPacBio_chunk_idx=None,  # merge mode: look for existing chunk BAMs
        # Use the canonical sample_id as prefix so the rectified BAM lands at
        # <sample_id>.rectified.bam (matching _rectified_bam_path). With an
        # empty prefix, align_command falls back to args.reads.stem — which
        # for DRS inputs is "<sample>_trimmed" after Step 0, breaking the
        # post-align rectified.bam lookup at stages.py:178.
        prefix=sample_id,
        keep_sam=False,
        sort=True,
        index=True,
        verbose=False,
        checkpoint_dir=checkpoint_dir,
        trust_existing_bams=trust_existing_bams,
        _run_provenance=_run_provenance,
    )

    # Temporarily hide scheduler array env vars so consensus.py runs in
    # single-sample mode. In run-all, array indices are for sample-level
    # parallelism; we never want within-sample read-level partitioning here.
    import os as _os
    _array_vars = {k: _os.environ.pop(k) for k in (
        # SLURM
        'SLURM_ARRAY_TASK_ID', 'SLURM_ARRAY_TASK_COUNT',
        'SLURM_ARRAY_TASK_MAX', 'SLURM_ARRAY_TASK_MIN',
        'SLURM_ARRAY_TASK_STEP',
        # UGE/SGE
        'SGE_TASK_ID', 'SGE_TASK_FIRST', 'SGE_TASK_LAST', 'SGE_TASK_STEPSIZE',
        # PBS/Torque
        'PBS_ARRAY_INDEX', 'PBS_ARRAYID',
    ) if k in _os.environ}
    try:
        rc = run_align(align_args)
    finally:
        _os.environ.update(_array_vars)
    if rc != 0:
        raise RuntimeError(f"Triple-aligner failed for {input_path}")

    if not rectified_bam.exists():
        raise RuntimeError(
            f"Alignment completed but rectified.bam not found: {rectified_bam}"
        )

    per_aligner_bams = _collect_per_aligner_bams(
        sample_id, sample_output_dir,
        run_provenance=_run_provenance, trust_existing_bams=trust_existing_bams,
    )
    return per_aligner_bams, rectified_bam
def _run_correction(
    bam_path: Path,
    output_dir: Path,
    genome_path: Path,
    annotation_path: Optional[Path],
    args,
    aligner_bams: Optional[List[Path]] = None,
) -> Path:
    """
    Run rectify correct on a BAM, writing corrected_reads.tsv into output_dir.
    Returns path to the corrected TSV.

    aligner_bams: Optional list of per-aligner BAM paths to feed Module 2H
        (junction refinement). Required to match the CLAUDE.md "validated"
        correct-first pipeline order. When provided, every per-aligner
        correction gets the full cross-aligner junction candidate pool.
    """
    from .. import correct_command

    corrected_tsv = output_dir / 'corrected_reads.tsv'

    # Indel correction and variant-aware rescue require MD tags in the BAM.
    # Disable them gracefully when MD tags are absent (e.g. rectified BAMs
    # generated by older runs without --MD).
    has_md = _bam_has_md_tags(bam_path)
    # --dT-primed-cDNA enables poly-A trimming / indel-correction modules.
    # Old --no-polya-sequenced is a deprecated alias for NOT passing --dT-primed-cDNA.
    # Default (no flag): direct RNA mode — poly-A IS in the read.
    polya_seq = getattr(args, 'dT_primed_cDNA', False)
    skip_indel = not has_md  # Can't correct indels without MD tags
    skip_variant = not has_md

    if not has_md and polya_seq:
        print(
            "    Note: BAM lacks MD tags — indel correction disabled. "
            "Re-run alignment with current RECTIFY to enable it."
        )

    # Auto-generate report path alongside the corrected TSV
    report_path = corrected_tsv.parent / (corrected_tsv.stem + '_report.html')

    # Derive sample prefix from input BAM stem (strip .rectified / .consensus suffixes)
    _stem = bam_path.stem
    for _sfx in ('.rectified', '.consensus'):
        if _stem.endswith(_sfx):
            _stem = _stem[:-len(_sfx)]
            break
    corrected_bam_path = (
        output_dir / f"{_stem}.rectified_corrected_3end.bam"
        if getattr(args, 'write_corrected_bam', False) else None
    )
    softclipped_bam_path = (
        output_dir / f"{_stem}.rectified_pA_tail_trimmed.bam"
        if getattr(args, 'write_softclip_bam', False) else None
    )

    # Module 2H (junction refinement) requires the full per-aligner BAM pool to
    # build the candidate junction set. Falling back to args.aligner_bams keeps
    # backwards-compat for callers that set the attribute directly.
    _aligner_bams = aligner_bams or getattr(args, 'aligner_bams', None) or []

    correct_args = argparse.Namespace(
        input=bam_path,
        genome=genome_path,
        annotation=annotation_path,
        output=corrected_tsv,
        corrected_bam=corrected_bam_path,
        softclipped_bam=softclipped_bam_path,
        netseq_dir=getattr(args, 'netseq_dir', None),
        organism=getattr(args, 'organism', None),
        aligner=getattr(args, 'aligner', 'minimap2'),
        aligner_bams=[str(p) for p in _aligner_bams],
        dT_primed_cDNA=polya_seq,
        polya_sequenced=polya_seq,  # deprecated attr kept for compat
        threads=getattr(args, 'threads', 4),
        filter_spikein=getattr(args, 'filter_spikein', None),
        streaming=getattr(args, 'streaming', False),
        # Defaults
        min_mapq=10,
        skip_secondary=True,
        skip_supplementary=True,
        skip_ag_check=False,
        skip_atract_check=False,
        skip_polya_trim=False,
        skip_indel_correction=skip_indel,
        skip_variant_aware=skip_variant,
        polya_model=None,
        report=report_path,
        max_downstream_a=20,
        chunk_size=getattr(args, 'chunk_size', 10000),
        junction_penalty_table=getattr(args, 'junction_penalty_table', None),
        str_penalty_table=getattr(args, 'str_penalty_table', None),
        junction_overhang_table=getattr(args, 'junction_overhang_table', None),
        junction_pool_cache=getattr(args, 'junction_pool_cache', None),
        debug=False,
        verbose=False,
    )

    correct_command.run(correct_args)

    if not corrected_tsv.exists():
        raise RuntimeError(f"Correction did not produce output: {corrected_tsv}")

    return corrected_tsv
def _run_correction_per_aligner(
    per_aligner_bams: Dict[str, Path],
    output_dir: Path,
    genome_path: Path,
    annotation_path: Optional[Path],
    args,
) -> Tuple[Dict[str, Path], Dict[str, Path]]:
    """
    Run ``rectify correct`` on each per-aligner BAM independently.

    Writes per-aligner outputs to ``output_dir/per_aligner_corrected/{aligner}/``.
    Skips any aligner whose TSV already exists (safe to resume).

    Returns
    -------
    Tuple ``(per_aligner_tsvs, per_aligner_corrected_bams)``:
        - ``per_aligner_tsvs``: Dict mapping aligner name → corrected_reads.tsv
          path for each aligner whose correction succeeded.
        - ``per_aligner_corrected_bams``: Dict mapping aligner name →
          rectified_corrected_3end.bam path for each aligner whose correction
          succeeded AND produced a corrected BAM. Empty when no corrected BAMs
          were materialized (e.g. ``--write-corrected-bam`` omitted upstream).
          When non-empty, this dict is the input ``merge_corrected_tsvs`` needs
          to activate HP-edit-distance mode for winner selection.
    """
    per_aligner_dir = output_dir / 'per_aligner_corrected'
    per_aligner_dir.mkdir(parents=True, exist_ok=True)

    per_aligner_tsvs: Dict[str, Path] = {}
    per_aligner_corrected_bams: Dict[str, Path] = {}
    for aligner_name, bam_path in per_aligner_bams.items():
        aligner_output_dir = per_aligner_dir / aligner_name
        aligner_output_dir.mkdir(exist_ok=True)
        tsv_path = aligner_output_dir / 'corrected_reads.tsv'

        # Validate TSV file exists and is readable (not truncated/corrupt)
        _tsv_valid = False
        if tsv_path.exists():
            try:
                import pandas as pd
                _df = pd.read_csv(tsv_path, sep='\t', nrows=1)
                _tsv_valid = len(_df) > 0
            except Exception:
                _tsv_valid = False

        if _tsv_valid:
            print(f"    [{aligner_name}] Skipping — TSV exists: {tsv_path}")
            per_aligner_tsvs[aligner_name] = tsv_path
            continue

        print(f"    [{aligner_name}] Correcting {bam_path.name}...")
        # rectify correct uses pysam.fetch() which requires a coordinate-sorted,
        # indexed BAM.  Aligners that sort internally (uLTRA, deSALT) have a .bai
        # placed by align_command.py; aligners that emit unsorted output
        # (mapPacBio, minimap2, gapmm2) do not.  Sort + index to a temp file when
        # the .bai is absent, then clean up after correction.
        sorted_bam: Optional[Path] = None
        correction_bam = bam_path
        bai_path = Path(str(bam_path) + '.bai')
        if not bai_path.exists():
            try:
                import pysam as _pysam
                sorted_bam = bam_path.with_suffix('.coord_sorted.bam')
                print(f"    [{aligner_name}] Coordinate-sorting for correction...")
                _pysam.sort('-m', '1G', '-o', str(sorted_bam), str(bam_path))
                _pysam.index(str(sorted_bam))
                correction_bam = sorted_bam
            except Exception as _sort_exc:
                print(
                    f"    WARNING: [{aligner_name}] sort failed ({_sort_exc}); "
                    "skipping correction",
                    file=sys.stderr,
                )
                if sorted_bam is not None and sorted_bam.exists():
                    sorted_bam.unlink(missing_ok=True)
                continue
        try:
            # Pass ALL per-aligner BAMs so Module 2H (junction refinement) can
            # build the cross-aligner candidate junction pool — the documented
            # "correct-first" pipeline order in CLAUDE.md. Without this every
            # per-aligner correction would run with Module 2H silently disabled.
            _all_aligner_bams = list(per_aligner_bams.values())
            _run_correction(
                bam_path=correction_bam,
                output_dir=aligner_output_dir,
                genome_path=genome_path,
                annotation_path=annotation_path,
                args=args,
                aligner_bams=_all_aligner_bams,
            )
        except Exception as exc:
            print(
                f"    WARNING: [{aligner_name}] correction failed: {exc}",
                file=sys.stderr,
            )
        finally:
            if sorted_bam is not None:
                sorted_bam.unlink(missing_ok=True)
                sorted_bai = Path(str(sorted_bam) + '.bai')
                sorted_bai.unlink(missing_ok=True)

        if tsv_path.exists():
            per_aligner_tsvs[aligner_name] = tsv_path
            # Find the corrected BAM emitted by _run_correction.  Naming
            # convention: ``{input_bam.stem stripped of .rectified/.consensus}
            # .rectified_corrected_3end.bam``.  When ``--write-corrected-bam``
            # was disabled upstream this file won't exist; that's fine —
            # absence triggers legacy 5-key sort in merge_corrected_tsvs.
            _stem = (sorted_bam.stem if sorted_bam is not None else bam_path.stem)
            for _sfx in ('.rectified', '.consensus', '.coord_sorted'):
                if _stem.endswith(_sfx):
                    _stem = _stem[:-len(_sfx)]
                    break
            _corr_bam = aligner_output_dir / f"{_stem}.rectified_corrected_3end.bam"
            if _corr_bam.exists():
                per_aligner_corrected_bams[aligner_name] = _corr_bam
        else:
            print(
                f"    WARNING: [{aligner_name}] correction produced no output",
                file=sys.stderr,
            )

    return per_aligner_tsvs, per_aligner_corrected_bams
def _combine_corrected_tsvs(
    samples: List[Dict[str, str]],
    output_dir: Path,
) -> Path:
    """
    Concatenate per-sample corrected TSVs into a single multi-sample TSV.

    Adds a `sample` column equal to sample_id so the analyze command can
    distinguish samples for DESeq2.
    """
    import pandas as pd

    combined_dir = output_dir / 'combined'
    combined_dir.mkdir(parents=True, exist_ok=True)
    combined_tsv = combined_dir / 'corrected_reads_combined.tsv'

    dfs = []
    missing = []
    for sample in samples:
        tsv_path = output_dir / sample['sample_id'] / 'corrected_reads.tsv'
        if not tsv_path.exists():
            missing.append(sample['sample_id'])
            print(f"  WARNING: {tsv_path} not found, skipping", file=sys.stderr)
            continue
        df = pd.read_csv(tsv_path, sep='\t')
        df['sample'] = sample['sample_id']
        dfs.append(df)
        print(f"  {sample['sample_id']}: {len(df):,} reads")

    if not dfs:
        raise RuntimeError(
            "No corrected TSVs found to combine. "
            "Ensure the correction step completed successfully."
        )

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(combined_tsv, sep='\t', index=False)
    print(f"  Combined → {combined_tsv}  ({len(dfs)} samples, {len(combined):,} reads)")

    if missing:
        print(
            f"  WARNING: {len(missing)} samples missing corrected TSV: {', '.join(missing)}",
            file=sys.stderr,
        )

    return combined_tsv
def _run_analysis(
    corrected_tsv: Path,
    output_dir: Path,
    genome_path: Optional[Path],
    annotation_path: Optional[Path],
    args,
    n_samples: int = 1,
) -> None:
    """Run the analyze command on a (possibly combined) corrected TSV."""
    from ..analyze_command import run_analyze

    # Only run DESeq2 when there are multiple samples
    run_deseq2 = n_samples > 1

    analyze_args = argparse.Namespace(
        input=corrected_tsv,
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
        run_motif=run_deseq2,   # motif discovery only makes sense with DESeq2 results
        sample_sets=None,
        # Filtering
        exclude_mito=True,
        include_mito=False,
        exclude_rdna=True,
        include_rdna=False,
        # Output options
        no_bedgraph=False,
        bedgraph_dir=None,
        no_genomic_distribution=False,
        # TSS clustering window (75 bp for DRS — 5' ends noisier than 3' ends)
        tss_cluster_distance=75,
        # Manifest mode (not used in single-file path)
        manifest=None,
        # Motif windows
        motif_upstream=100,
        motif_downstream=50,
    )

    exit_code = run_analyze(analyze_args)
    if exit_code != 0:
        print(f"\nAnalysis completed with warnings (exit code: {exit_code})")

def _run_junction_aggregation(
    bam_path: Path,
    genome_path: Path,
    annotation_path: Path,
    output_dir: Path,
    config: Optional[dict] = None,
) -> Optional[Path]:
    """
    Run junction aggregation with partial rescue.
    Only runs if annotation is GFF/GFF3 (has intron features).
    Returns path to junctions TSV, or None if skipped/failed.
    """
    # Handle compressed files: .gff.gz → stem is .gff, suffixes[-1] is .gz
    suffixes = [s.lower() for s in annotation_path.suffixes]
    has_gff = any(s in ('.gff', '.gff3') for s in suffixes)
    if not has_gff:
        print(
            "\n[Junctions] Skipping — requires GFF/GFF3 annotation "
            "(current annotation has no intron features)"
        )
        return None

    print("\n[Junctions] Aggregating splice junctions with partial rescue...")

    from ...aggregate.junctions import aggregate_junctions, merge_with_partial_evidence, export_junctions
    from ...splice.terminal_exon_refiner import load_splice_sites_from_gff, detect_partial_junction_crossings
    import pysam

    junctions_dir = output_dir / 'junctions'
    junctions_dir.mkdir(parents=True, exist_ok=True)
    junctions_tsv = junctions_dir / 'junctions.tsv'

    try:
        genome = {}
        fasta = pysam.FastaFile(str(genome_path))
        for chrom in fasta.references:
            genome[chrom] = fasta.fetch(chrom)
        fasta.close()

        junction_df = aggregate_junctions(
            bam_path=str(bam_path), genome=genome, min_reads=1,
        )
        print(f"  Found {len(junction_df)} junctions from CIGAR")

        splice_index = load_splice_sites_from_gff(str(annotation_path))
        partial_results = detect_partial_junction_crossings(
            bam_path=str(bam_path),
            genome=genome,
            splice_index=splice_index,
            min_clip_length=1,
            ambiguous_mode='proportional',
        )

        n_rescued = partial_results['stats']['reads_rescued_as_spliced']
        n_ambiguous = partial_results['stats']['reads_ambiguous']
        print(f"  Rescued {n_rescued} partial crossings ({n_ambiguous} ambiguous)")

        junction_df = merge_with_partial_evidence(
            junction_df, partial_results, ambiguous_mode='proportional',
        )
        export_junctions(junction_df, str(junctions_tsv), format='tsv')
        print(f"  Junctions written to {junctions_tsv}")

        # Provenance
        if junctions_tsv.exists():
            from ....utils.provenance import init_provenance
            prov = init_provenance(output_dir, description="RECTIFY junction aggregation", config=config)
            prov.add_output_file(
                junctions_tsv,
                source_files=[bam_path],
                metadata={
                    'n_junctions': len(junction_df),
                    'n_rescued': n_rescued,
                    'n_ambiguous': n_ambiguous,
                },
            )
            prov.save()

        return junctions_tsv

    except Exception as e:
        print(f"\nWarning: Junction aggregation failed: {e}", file=sys.stderr)
        print("Continuing without junction output...")
        return None
