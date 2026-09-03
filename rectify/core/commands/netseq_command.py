"""
NET-seq command for RECTIFY CLI.

Usage:
    rectify netseq input.bam --genome sacCer3.fa --gff genes.gff -o output/

    # With deconvolution disabled (raw 3' ends only)
    rectify netseq input.bam --genome sacCer3.fa --no-deconvolution -o output/

    # Process multiple samples
    rectify netseq sample1.bam sample2.bam --genome sacCer3.fa -o output/

Author: Kevin R. Roy
"""

import json
import sys
from pathlib import Path
from typing import List, Optional

from ..bam.netseq_bam_processor import (
    NETSEQ_RNA3P_CHOICES,
    NETSEQ_RNA3P_DEFAULT,
    process_netseq_bam,
    aggregate_netseq_stream,
    aggregate_positions,
)
from ..netseq.netseq_rescue import JunctionPool, allowed_remainders, load_junction_tsv
from ..netseq.netseq_umi import (
    NETSEQ_UMI_LENGTH,
    NETSEQ_UMI_SOURCES,
    SATURATION_KU_THRESHOLD,
    iter_netseq_fragments,
    stream_netseq_positions,
)
from ..umi.dedup import UmiDedupStats
from ..netseq.netseq_deconvolution import deconvolve_all_regions
from ..netseq.netseq_output import export_netseq_results, write_exclusion_stats
from ..exclusion_regions import ExclusionRegionDetector
from ...utils.genome import load_genome


def run_netseq(args) -> int:
    """
    Execute NET-seq processing pipeline.

    Args:
        args: Parsed command line arguments

    Returns:
        Exit code (0 for success)
    """
    print("=" * 60)
    print("RECTIFY NET-seq Processing Pipeline")
    print("=" * 60)

    rna3p_at = getattr(args, 'rna3p_at', NETSEQ_RNA3P_DEFAULT)
    umi_length = int(getattr(args, 'umi_length', 0) or 0)
    do_dedup = bool(getattr(args, 'dedup', False))
    if rna3p_at == 'read3p':
        print("\n  !! rna3p-at=read3p: the read is treated as SENSE. That is the LEGACY assumption and it is")
        print("     wrong for Churchman-style NET-seq (measured: gene strand is the inverse of the BAM strand,")
        print("     and the RNA 3' end is the read's 5' terminus). Use it only to reproduce old output.")
    if do_dedup and umi_length <= 0:
        print("Error: --dedup requires --umi-length (NET-seq randomer is "
              f"{NETSEQ_UMI_LENGTH} nt for GSE159603/GSE61332-style libraries)")
        return 1

    # Validate inputs
    input_paths = args.input
    output_dir = Path(args.output_dir)
    if args.genome is None:
        print("Error: --genome is required (or use --Scer / --organism for bundled data)")
        return 1
    genome_path = Path(args.genome)

    if not genome_path.exists():
        print(f"Error: Genome file not found: {genome_path}")
        return 1

    for input_path in input_paths:
        if not Path(input_path).exists():
            print(f"Error: Input file not found: {input_path}")
            return 1

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- output-format preflight ---------------------------------------------------------------
    # Fail (or warn) HERE, not after a multi-hour BAM pass. Which formats keep the streaming path:
    #   bedgraph / bigwig -> STREAMING (built from the position dicts; records are never held)
    #   parquet / tsv     -> the read-level table, so every record is materialised (~23 GB / 50 M)
    _formats = getattr(args, 'output_format', ['parquet', 'bedgraph'])
    if 'bigwig' in _formats:
        from ..netseq.netseq_output import HAS_PYBIGWIG
        if not HAS_PYBIGWIG:
            print("Error: --output-format bigwig requires pyBigWig, which is not importable.\n"
                  "       Install it (pip install pyBigWig) or drop 'bigwig' from --output-format.\n"
                  "       Refusing to start: the .bw files would be silently missing at the end of "
                  "the run.")
            return 1
    if 'parquet' in _formats:
        from ..unified_record import HAS_PYARROW
        if not HAS_PYARROW:
            print("  WARNING: pyarrow is not importable, so --output-format parquet will fall back "
                  "to TSV.\n"
                  "           That fallback MATERIALISES every read record (~23 GB per 50M reads). "
                  "For a large\n"
                  "           library either install pyarrow or use --output-format bedgraph bigwig,"
                  " which stream.")

    # ---- annotation resolution -----------------------------------------------------------------
    # resolve_reference_paths() fills args.annotation, never args.gff, and this parser has no
    # --annotation slot -- so under --Scer the GFF was silently absent (exclusion regions fell back
    # to the hard-coded rDNA box, and a --gff-gated junction rescue would never fire). Resolve here.
    annotation_path = _resolve_annotation(args)

    # ---- genome --------------------------------------------------------------------------------
    # Needed by deconvolution AND by the poly(A) walkback and the donor-side junction rescue, so it
    # is loaded unless every consumer is off.
    detect_tail = not getattr(args, 'no_tail_detection', False)
    want_rescue = _want_rescue(args, annotation_path)
    genome = None
    if (not args.no_deconvolution) or detect_tail or want_rescue:
        print(f"\nLoading genome: {genome_path}")
        genome = load_genome(str(genome_path))
        print(f"  Loaded {len(genome)} chromosomes")

    # ---- junction pool -------------------------------------------------------------------------
    junction_pool = None
    if want_rescue:
        pool_tsv = getattr(args, 'junction_pool', None)
        junction_pool = JunctionPool()
        if annotation_path:
            include_trna = bool(getattr(args, 'pool_include_trna', False))
            include_organellar = bool(getattr(args, 'pool_include_organellar', False))
            junction_pool = JunctionPool.from_annotation(
                annotation_path, include_trna=include_trna,
                include_organellar=include_organellar)
            kept = ", ".join(f"{t} {n}" for t, n in
                             sorted(junction_pool.by_parent_type.items(), key=lambda kv: -kv[1]))
            print(f"\nJunction pool: {len(junction_pool):,} annotated introns from "
                  f"{Path(annotation_path).name}")
            print(f"  by parent feature type: {kept or 'unknown'}")
            if junction_pool.dropped_by_parent_type:
                dropped = ", ".join(f"{t} {n}" for t, n in
                                    sorted(junction_pool.dropped_by_parent_type.items()))
                print(f"  DROPPED (not spliceosomal -- tRNA endonuclease introns and "
                      f"mitochondrial group I/II self-splicing introns; --pool-include-trna / "
                      f"--pool-include-organellar to keep): {dropped}")
        if pool_tsv:
            n_before = len(junction_pool)
            for j in load_junction_tsv(pool_tsv):
                junction_pool.add(j)
            print(f"  + {len(junction_pool) - n_before:,} junctions from {Path(pool_tsv).name}")
        if not len(junction_pool):
            print("  WARNING: junction rescue requested but the pool is EMPTY -- rescue disabled")
            junction_pool = None
        else:
            print(f"  rescue window: aligned end within {args.rescue_max_intronic} nt past the donor; "
                  f"min k = {args.rescue_min_k} (r == 0) / "
                  f"{args.rescue_min_k_with_remainder} (r > 0); allowed remainder = "
                  f"{list(allowed_remainders(umi_length))}")

    # Set up exclusion regions
    exclusion_detector = None
    exclude_rdna = not getattr(args, 'include_rdna', False)
    exclude_pol3 = not getattr(args, 'include_pol3', False)

    if exclude_rdna or exclude_pol3:
        print("\nSetting up exclusion regions...")
        exclusion_detector = ExclusionRegionDetector(flanking_bp=args.pol3_flanking)

        # 🔴 EXCLUSION uses the EXPLICIT --gff only, never the organism fallback.
        # `_resolve_annotation` exists so the junction pool can be built under --Scer, but feeding
        # it here too would silently turn `rectify netseq --Scer` from "exclude the rDNA box" into
        # "exclude 610 regions incl. every tRNA and snRNA" -- and the 829 QC gate's primary
        # observable is the SNR7 (U5 snRNA) 3' end, which that would delete. Pass --gff to opt in.
        gff_path = getattr(args, 'gff', None)
        if gff_path and Path(gff_path).exists():
            n_regions = exclusion_detector.load_from_gff(
                Path(gff_path),
                exclude_tRNA=exclude_pol3,
                exclude_snRNA=exclude_pol3,
                exclude_rDNA=exclude_rdna,
                exclude_mito=args.exclude_mito,
            )
            print(f"  Loaded {n_regions} exclusion regions from GFF")
        else:
            # Use default yeast rDNA coordinates
            if exclude_rdna:
                exclusion_detector.add_rdna_region()
                print("  Added default yeast rDNA region (chrXII:450,000-490,000)")
            if annotation_path:
                print(f"  NOTE: {Path(annotation_path).name} was resolved for the JUNCTION POOL "
                      "only. Pass --gff to use it for exclusion regions too (that would also "
                      "exclude every tRNA and snRNA unless --include-pol3 is set).")

        # Show summary
        stats = exclusion_detector.get_stats_by_reason()
        for reason, count in stats.items():
            print(f"  {reason}: {count} regions")

    # Process each input BAM
    for input_path in input_paths:
        input_path = Path(input_path)
        sample_name = input_path.stem.replace('.bam', '').replace('_bbmap_correct', '')

        print(f"\n{'='*60}")
        print(f"Processing: {input_path.name}")
        print(f"Sample name: {sample_name}")
        print(f"{'='*60}")

        output_formats = args.output_format if hasattr(args, 'output_format') else ['parquet', 'bedgraph']
        want_tsv = 'tsv' in output_formats
        # Both read-level formats need the records materialised; the position tracks never do.
        want_parquet = 'parquet' in output_formats
        want_records = want_parquet or want_tsv

        # Process BAM
        print("\n1. Processing BAM file...")
        record_stream = process_netseq_bam(
            input_path,
            exclusion_detector=exclusion_detector,
            min_mapq=args.min_mapq if hasattr(args, 'min_mapq') else 0,
            min_a_fraction=0.8,
            max_reads=args.max_reads if hasattr(args, 'max_reads') else None,
            show_progress=True,
            rna3p_at=rna3p_at,
            junction_pool=junction_pool,
            genome=genome,
            umi_length=umi_length,
            rescue_max_intronic=args.rescue_max_intronic,
            rescue_min_k=args.rescue_min_k,
            rescue_min_k_with_remainder=args.rescue_min_k_with_remainder,
            detect_tail=detect_tail,
            walkback_requires_clip_a=getattr(args, 'walkback_requires_clip_a', True),
        )

        # Only the read-level parquet needs every record in memory (~23 GB per 50M reads).
        # Position tracks AND the correction counters are accumulated in ONE streaming pass when
        # parquet is not requested, so a 100-180M-read library fits in a few GB.
        print("\n2. Aggregating positions (raw + corrected, single pass)...")
        raw_counts, corrected_counts, summary, records = aggregate_netseq_stream(
            record_stream, collect=want_records,
        )
        total_reads = summary.reads
        if not total_reads:
            print("  Warning: No records generated!")
            continue
        summary.junction_pool_size = len(junction_pool) if junction_pool else 0
        print(f"  Unique positions: raw {len(raw_counts):,}  corrected {len(corrected_counts):,}  "
              f"(reads: {total_reads:,})")
        _print_correction_summary(summary)

        # Which track the deconvolution and the downstream analyses key on.
        track_position = getattr(args, 'track_position', 'corrected')
        primary_counts = corrected_counts if track_position == 'corrected' else raw_counts
        print(f"  Primary track (deconvolution input): {track_position}")

        # Deconvolution
        deconv_counts = None
        if not args.no_deconvolution and genome:
            print("\n3. Applying deconvolution...")
            deconv_counts, deconv_results = deconvolve_all_regions(
                primary_counts,
                genome,
                min_downstream_a=args.min_atract_length if hasattr(args, 'min_atract_length') else 3,
                show_progress=True,
            )
        else:
            print("\n3. Skipping deconvolution (--no-deconvolution)")

        # Export outputs
        print("\n4. Exporting results...")
        normalize_rpm = not getattr(args, 'no_rpm_normalize', False)

        outputs = export_netseq_results(
            records=records,
            raw_counts=raw_counts,
            deconv_counts=deconv_counts,
            output_dir=output_dir,
            sample_name=sample_name,
            total_reads=total_reads,
            export_parquet=want_parquet,
            export_bedgraph='bedgraph' in output_formats,
            export_bigwig='bigwig' in output_formats,
            normalize_rpm=normalize_rpm,
            corrected_counts=corrected_counts,
            export_tsv=want_tsv,
        )

        summary_path = output_dir / f"{sample_name}.netseq_summary.json"
        payload = summary.as_dict()
        payload.update({
            'sample': sample_name, 'input_bam': str(input_path), 'rna3p_at': rna3p_at,
            'umi_length': umi_length, 'track_position': track_position,
            'junction_rescue': junction_pool is not None,
            'junction_pool_size': len(junction_pool) if junction_pool else 0,
            'junction_pool_by_parent_type': dict(junction_pool.by_parent_type) if junction_pool else {},
            'junction_pool_dropped_by_parent_type':
                dict(junction_pool.dropped_by_parent_type) if junction_pool else {},
            'pool_include_trna': bool(getattr(args, 'pool_include_trna', False)),
            'pool_include_organellar': bool(getattr(args, 'pool_include_organellar', False)),
            'rescue_max_intronic': args.rescue_max_intronic,
            'rescue_min_k': args.rescue_min_k,
            'rescue_min_k_with_remainder': args.rescue_min_k_with_remainder,
            'allowed_remainders': list(allowed_remainders(umi_length)),
            'tail_detection': detect_tail,
            'walkback_requires_clip_a': getattr(args, 'walkback_requires_clip_a', True),
            'annotation': str(annotation_path) if annotation_path else None,
        })
        summary_path.write_text(json.dumps(payload, indent=1))
        outputs['netseq_summary'] = summary_path

        print(f"\n  Generated {len(outputs)} output files")
        for name, path in outputs.items():
            print(f"    {name}: {path.name}")

        # ---- molecule-level (UMI) track ------------------------------------------------------------
        # Emitted ALONGSIDE the read-level track, never instead of it: a 6-nt UMI saturates at deep
        # positions, so molecules alone would be silently compressed exactly where signal is highest.
        if do_dedup:
            _emit_molecule_track(
                input_path=input_path, sample_name=sample_name, output_dir=output_dir,
                umi_length=umi_length, umi_source=getattr(args, 'umi_source', 'read5p'),
                clustering=getattr(args, 'umi_clustering', 'exact'),
                edit_distance=int(getattr(args, 'umi_edit_distance', 0) or 0),
                exclusion_detector=exclusion_detector,
                min_mapq=args.min_mapq if hasattr(args, 'min_mapq') else 0,
                rna3p_at=rna3p_at,
                max_reads=args.max_reads if hasattr(args, 'max_reads') else None,
                # SAME correction inputs as the read pass, or the two tracks sit on different
                # coordinates for every walked-back / rescued read (planning 834 defect (c)).
                genome=genome, junction_pool=junction_pool,
                rescue_max_intronic=args.rescue_max_intronic,
                rescue_min_k=args.rescue_min_k,
                rescue_min_k_with_remainder=args.rescue_min_k_with_remainder,
                detect_tail=detect_tail,
                walkback_requires_clip_a=getattr(args, 'walkback_requires_clip_a', True),
            )

    print(f"\n{'='*60}")
    print("NET-seq processing complete!")
    print(f"Output directory: {output_dir}")
    print(f"{'='*60}")

    return 0



def _resolve_annotation(args):
    """Annotation path for exclusion regions AND the junction pool.

    Priority: explicit ``--gff`` -> ``args.annotation`` (filled by ``resolve_reference_paths`` when
    the parser exposes it) -> the bundled annotation for ``--organism``/``--Scer``. The netseq parser
    has no ``--annotation`` slot, so without this a bundled-organism run silently had NO annotation:
    exclusion fell back to the hard-coded rDNA box and the junction rescue would never fire.
    """
    gff = getattr(args, 'gff', None)
    if gff and Path(gff).exists():
        return Path(gff)
    annotation = getattr(args, 'annotation', None)
    if annotation and Path(annotation).exists():
        return Path(annotation)
    organism = getattr(args, 'organism', None)
    if organism:
        from rectify.data import get_bundled_annotation_path, normalize_organism
        bundled = get_bundled_annotation_path(normalize_organism(organism))
        if bundled and Path(bundled).exists():
            return Path(bundled)
    return None


def _want_rescue(args, annotation_path) -> bool:
    """Whether the donor-side junction rescue runs.

    Default ON as soon as a junction source exists (an annotation or ``--junction-pool``);
    ``--no-junction-rescue`` forces it off, ``--junction-rescue`` makes its absence an explicit
    request that warns rather than silently doing nothing.
    """
    if getattr(args, 'junction_rescue', None) is False:
        return False
    have_source = bool(annotation_path) or bool(getattr(args, 'junction_pool', None))
    if getattr(args, 'junction_rescue', None) is True and not have_source:
        print("  WARNING: --junction-rescue given but no --gff/--organism annotation and no "
              "--junction-pool -- there is nothing to rescue against")
        return False
    return have_source


def _print_correction_summary(summary) -> None:
    """One compact block per sample; the full detail lands in <sample>.netseq_summary.json."""
    print(f"  near a pooled donor: {summary.reads_near_donor:,}"
          f"   rescued: {summary.rescued:,}"
          f" (mis-extended into intron: {summary.rescued_mis_extended:,})"
          f"   exon1_end: {summary.exon1_end:,}"
          f"   ambiguous: {summary.ambiguous:,}"
          f"   intronic_end: {summary.intronic_end:,}")
    if summary.rescued:
        by_k = ", ".join(f"k={k}:{n}" for k, n in sorted(summary.rescued_by_k.items()))
        print(f"    rescued by k: {by_k}")
        decoy = ", ".join(f"k={k}:{n}" for k, n in sorted(summary.decoy_k_at_donor.items()) if k)
        print(f"    decoy-acceptor null (same reads, exon2 + 50 nt): {decoy or 'k>=1: 0'}")
        print(f"    the SAME rule applied to the decoy acceptor would have rescued "
              f"{summary.decoy_rescued:,} of these reads -- that is the chance-match floor")
    pct = (lambda n: 100.0 * n / summary.reads if summary.reads else 0.0)
    print(f"  tails: >=1 nt {summary.tailed_ge1:,} ({pct(summary.tailed_ge1):.2f}%)"
          f"   >=3 nt {summary.tailed_ge3:,} ({pct(summary.tailed_ge3):.2f}%)"
          f"   partly aligned (walkback>0) {summary.tail_walkback_gt0:,}")
    print(f"    tail evidence: {summary.tail_clip_evidence:,} reads carry a non-templated A next to "
          f"the alignment; {summary.tail_walkback_only:,} were called from aligned A's alone")
    print(f"  corrected end differs from raw for {summary.end_moved:,} reads "
          f"({pct(summary.end_moved):.2f}%)")


def _emit_molecule_track(
    input_path,
    sample_name: str,
    output_dir: Path,
    umi_length: int,
    umi_source: str,
    clustering: str,
    edit_distance: int,
    exclusion_detector,
    min_mapq: int,
    rna3p_at: str,
    max_reads,
    genome=None,
    junction_pool=None,
    rescue_max_intronic: int = 10,
    rescue_min_k: int = 1,
    rescue_min_k_with_remainder: int = 4,
    detect_tail: bool = True,
    walkback_requires_clip_a: bool = True,
) -> None:
    """Write the UMI molecule-level outputs for one sample.

    Two files:
      ``<sample>_umi_stats.json``      duplication rate + family-size histogram -- the calibration outputs,
                                       not decoration.
      ``<sample>_molecules.tsv``       per 3'-end position: reads, distinct UMIs (k), occupancy-corrected
                                       molecules (m_hat) and the SATURATED flag. Reads and molecules sit
                                       side by side ON PURPOSE -- a 6-nt UMI saturates, so a molecule count
                                       alone is silently compressed exactly where signal is highest.

    Deliberately NOT a per-molecule read list: that would be one row per molecule (~20M for a real library)
    and nothing downstream consumes it. Use ``core.netseq.netseq_umi.select_netseq_molecules`` directly if
    you need representative QNAMEs for a subset.
    """
    import json

    print("\n5. Molecule-level (UMI) counting...")
    print("   ⚠️  molecule positions carry the randomer-OVERSHOOT correction, which the read track")
    print("       does not: it assumes EVERY read carries the randomer, so on a MIXED library it")
    print("       shifts the randomer-free class by the full --umi-length (planning 829 §4). Verify")
    print("       the randomer is universal before using this track -- 5'-clip histogram, not the kit.")
    print(f"   umi-length={umi_length} umi-source={umi_source} clustering={clustering} "
          f"edit-distance={edit_distance}")
    stats = UmiDedupStats()
    mol_path = output_dir / f"{sample_name}_molecules.tsv"
    n_positions = n_sat = 0
    sum_reads = sum_k = 0
    sum_mhat = 0.0
    # Streamed, not materialised: a real library is ~50M reads over ~4M positions, and this pass runs on top
    # of the read-level pass's ~23 GB. See stream_netseq_positions().
    fragments = iter_netseq_fragments(
        str(input_path), umi_length=umi_length, umi_source=umi_source,
        exclusion_detector=exclusion_detector, min_mapq=min_mapq,
        rna3p_at=rna3p_at, max_reads=max_reads, stats=stats,
        genome=genome, junction_pool=junction_pool,
        rescue_max_intronic=rescue_max_intronic, rescue_min_k=rescue_min_k,
        rescue_min_k_with_remainder=rescue_min_k_with_remainder,
        detect_tail=detect_tail, walkback_requires_clip_a=walkback_requires_clip_a,
    )
    tmp = mol_path.with_suffix(".tsv.tmp")
    with open(tmp, "w") as fh:
        fh.write("chrom\tstrand\tposition\treads\tdistinct_umis\tmolecules_corrected\tsaturated\n")
        for p in stream_netseq_positions(fragments, umi_length=umi_length,
                                         edit_distance=edit_distance, clustering=clustering,
                                         stats=stats):
            n_positions += 1
            sum_reads += p.reads
            sum_k += p.distinct_umis
            if p.saturated:
                n_sat += 1
            mh = "inf" if p.molecules_corrected == float('inf') else f"{p.molecules_corrected:.2f}"
            if p.molecules_corrected != float('inf'):
                sum_mhat += p.molecules_corrected
            fh.write(f"{p.contig}\t{p.strand}\t{p.position}\t{p.reads}\t{p.distinct_umis}\t"
                     f"{mh}\t{int(p.saturated)}\n")
    tmp.replace(mol_path)

    if not n_positions:
        print("   Warning: no reads carried a usable UMI -- molecule track not written")
        return

    stats_path = output_dir / f"{sample_name}_umi_stats.json"
    payload = stats.as_dict()
    payload.update({
        "umi_length": umi_length, "umi_source": umi_source, "clustering": clustering,
        "edit_distance": edit_distance, "rna3p_at": rna3p_at,
        "n_positions": n_positions, "n_saturated_positions": n_sat,
        "saturation_threshold_k_over_U": SATURATION_KU_THRESHOLD,
        "sum_reads": sum_reads, "sum_distinct_umis": sum_k,
        "sum_occupancy_corrected_molecules": round(sum_mhat, 1),
        "occupancy_inflation": round(sum_mhat / sum_k, 4) if sum_k else None,
    })
    stats_path.write_text(json.dumps(payload, indent=1))

    print(f"   reads with UMI: {stats.n_input_fragments:,}   molecules: {stats.n_molecules:,}   "
          f"duplication rate: {100 * stats.duplication_rate:.1f}%")
    print(f"   positions: {n_positions:,}   SATURATED (k/U>{SATURATION_KU_THRESHOLD}, fall back to reads "
          f"there): {n_sat:,}")
    print("   ⚠️  duplication is depth-dependent: at positions deeper than ~1000 reads a 6-nt UMI is")
    print("       exhausted, so the aggregate rate above is an UPPER bound. Read _molecules.tsv, not this.")
    for name, path in (("umi_stats", stats_path), ("molecules", mol_path)):
        print(f"    {name}: {path.name}")


def add_netseq_parser(subparsers) -> None:
    """Add the netseq subcommand parser."""
    parser = subparsers.add_parser(
        'netseq',
        help='Process NET-seq BAM files to extract 3\' end positions',
        description=(
            'Extract and deconvolve 3\' end positions from NET-seq data. '
            'Outputs parquet files with per-read records and RPM-normalized '
            'bedgraph files with position-level signal.'
        ),
    )

    # Required arguments
    parser.add_argument(
        'input',
        nargs='+',
        type=Path,
        help='Input BAM file(s)',
    )

    parser.add_argument(
        '--genome', '-g',
        type=Path,
        default=None,
        help='Reference genome FASTA (required for deconvolution). '
             'Optional when --Scer or --organism is set.',
    )

    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory',
    )

    # Annotation for exclusion + the junction pool
    parser.add_argument(
        '--gff',
        type=Path,
        help='GFF/GTF annotation for exclusion regions AND the donor-side junction rescue pool. '
             'Defaults to the bundled annotation under --Scer/--organism.',
    )

    rescue_group = parser.add_argument_group(
        'Donor-side junction rescue',
        "A nascent RNA whose 3' end sits 1-10 nt into exon 2 cannot be anchored across the junction "
        "by a short-read aligner: the overhang is soft-clipped (STAR Local, bbmap) or the alignment "
        "is mis-extended a few nt into the intron, and the read is called at the 5' SPLICE SITE -- "
        "manufacturing a false splicing intermediate at exactly the coordinate the real signal "
        "occupies. This pass matches that clip to the start of exon 2 of a pooled intron and moves "
        "the 3' end there. ON by default whenever an annotation or --junction-pool is available.",
    )
    rescue_group.add_argument(
        '--junction-rescue',
        dest='junction_rescue', action='store_true', default=None,
        help='Force the rescue on (warns if no junction source is available).',
    )
    rescue_group.add_argument(
        '--no-junction-rescue',
        dest='junction_rescue', action='store_false',
        help='Disable the rescue.',
    )
    rescue_group.add_argument(
        '--pool-include-trna',
        action='store_true',
        help="Keep tRNA introns in the junction pool. OFF by default: a tRNA intron is excised by "
             "the tRNA endonuclease/ligase pathway, not the spliceosome, at a Pol III locus that "
             "this command already drops from its tracks -- so it can only manufacture rescues "
             "(measured on wt_rep3: the single locus YNCO0031W_tRNA collected 46). mRNA, snoRNA, "
             "snRNA and other Pol II introns are always kept; the pool composition by parent "
             "feature type is printed at startup.",
    )
    rescue_group.add_argument(
        '--pool-include-organellar',
        action='store_true',
        help="Keep mitochondrial/organellar introns in the junction pool. OFF by default, for the "
             "same reason as tRNA introns and with the same kind of measurement behind it: mito "
             "introns are group I/II SELF-SPLICING introns on a genome Pol II does not transcribe, "
             "so a nascent Pol II 3' end cannot legitimately sit at one. Their parent feature type "
             "is `mRNA`, so the tRNA filter does not catch them -- the CONTIG is the signal. "
             "Measured on wt_rep3 with chrMito reads included: 94 of 580 rescues (16%%) were on "
             "chrMito.",
    )
    rescue_group.add_argument(
        '--junction-pool',
        type=Path,
        default=None,
        help="External junction table (TSV, header chrom/donor/acceptor/strand), merged with the "
             "annotated introns -- e.g. a long-read-derived pool. BOTH sites are INTRONIC and on "
             "the GENE strand: donor = first intronic base, acceptor = last intronic base, 0-based "
             "(so donor > acceptor on a '-' gene).",
    )
    rescue_group.add_argument(
        '--rescue-max-intronic',
        type=int, default=10,
        help="How far past the donor the aligned RNA 3' end may sit and still be considered "
             "(aligner mis-extension into the intron). Default 10.",
    )
    rescue_group.add_argument(
        '--rescue-min-k',
        type=int, default=1,
        help="Minimum recovered exon-2 length for a rescue with NO non-templated remainder "
             "(r == 0, i.e. the clip is exon-2 sequence and nothing else). Default 1: a 1-nt "
             "overhang with nothing beyond it is legitimate evidence.",
    )
    rescue_group.add_argument(
        '--rescue-min-k-with-remainder',
        type=int, default=4,
        help="Minimum recovered exon-2 length when a randomer remainder is invoked to explain the "
             "rest of the clip (r > 0). Default 4, because THE CHANCE CHANNEL IS THE REMAINDER: a "
             "randomer's first base matches exon 2 a quarter of the time. Measured on a 504-read "
             "candidate set, the decoy acceptor produced k=1 70 times against 61 observed, k=2 24 "
             "against 4, k=3 10 against 3 -- indistinguishable from chance -- while k>=4 observed "
             "11 against 2. Set equal to --rescue-min-k for the old flat behaviour.",
    )
    rescue_group.add_argument(
        '--walkback-requires-clip-a',
        dest='walkback_requires_clip_a', action='store_true', default=True,
        help="ON BY DEFAULT for `rectify netseq`. Only run the poly(A) walkback on reads whose clip "
             "carries a non-templated A next to the alignment. Nascent 3' ends have no tail by "
             "default, and ~25%% of reads end on an A over a genomic A by chance, so the "
             "unconditional walkback moved 22%% of ends 1-5 nt on no evidence and ERASED the RPL32 "
             "exon-1 3' end peak (33 reads -> 1) -- the splicing-intermediate signal this command "
             "exists to measure. Use --walkback-unconditional for a poly(A)+ input.",
    )
    rescue_group.add_argument(
        '--walkback-unconditional',
        dest='walkback_requires_clip_a', action='store_false',
        help="Restore invariant-7 behaviour: walk back through terminal A's over genomic A's even "
             "with no non-templated A in the clip. Correct for a poly(A)-SELECTED input where every "
             "read really does carry a tail (and it is what rectify's own walkback_3prime does); "
             "wrong for nascent RNA.",
    )
    rescue_group.add_argument(
        '--no-tail-detection',
        action='store_true',
        help="Disable the non-templated 3'-tail call (clip A-run + genome-aware walkback).",
    )

    from rectify.data import add_organism_args
    add_organism_args(parser)

    # Exclusion options
    excl_group = parser.add_argument_group('Exclusion regions')
    excl_group.add_argument(
        '--include-rdna',
        action='store_true',
        help='Include rDNA locus (default: exclude)',
    )
    excl_group.add_argument(
        '--include-pol3',
        action='store_true',
        help='Include Pol III genes (tRNAs, SNR6, etc.) (default: exclude)',
    )
    excl_group.add_argument(
        '--pol3-flanking',
        type=int,
        default=100,
        help='Flanking bp around Pol III genes to exclude (default: 100)',
    )
    excl_group.add_argument(
        '--exclude-mito',
        action='store_true',
        default=False,
        help='Exclude mitochondrial genome (default: False)',
    )

    # Deconvolution options
    deconv_group = parser.add_argument_group('Deconvolution')
    deconv_group.add_argument(
        '--no-deconvolution',
        action='store_true',
        help='Disable NNLS deconvolution (output raw 3\' ends only)',
    )
    deconv_group.add_argument(
        '--min-atract-length',
        type=int,
        default=3,
        help='Minimum downstream A\'s to trigger deconvolution (default: 3)',
    )

    # Output options
    output_group = parser.add_argument_group('Output options')
    output_group.add_argument(
        '--output-format',
        nargs='+',
        choices=['parquet', 'bedgraph', 'bigwig', 'tsv'],
        default=['parquet', 'bedgraph'],
        help="Output formats to generate (default: parquet bedgraph). 'parquet' and 'tsv' are the "
             "read-level table (per-read rescue and tail fields); parquet degrades to TSV with a "
             "warning when pyarrow is not installed. 'bedgraph'/'bigwig' are the position tracks.",
    )
    output_group.add_argument(
        '--track-position',
        choices=['raw', 'corrected'],
        default='corrected',
        help="Which 3' end drives the PRIMARY track (deconvolution input): 'corrected' (default; "
             "after terminal oligo(A) trim, poly(A) walkback and junction rescue) or 'raw' (the "
             "bare alignment terminus, the pre-2026-09 behaviour). Both the .raw.* and .corrected.* "
             "bedgraphs are written either way, so before/after is always diffable.",
    )
    output_group.add_argument(
        '--no-rpm-normalize',
        action='store_true',
        help='Disable RPM normalization for bedgraph/bigwig',
    )

    # Processing options
    parser.add_argument(
        '--min-mapq',
        type=int,
        default=0,
        help='Minimum mapping quality (default: 0)',
    )

    parser.add_argument(
        '--max-reads',
        type=int,
        help='Maximum reads to process (for testing)',
    )

    parser.add_argument(
        '--rna3p-at',
        choices=list(NETSEQ_RNA3P_CHOICES),
        default=NETSEQ_RNA3P_DEFAULT,
        help=(
            "Which terminus of the SEQUENCED read is the RNA 3' end. "
            "'read5p' (default) = Churchman-style NET-seq: the read is revcomp(RNA), so the gene strand is "
            "the INVERSE of the BAM strand -- this is the MEASURED convention for GSE25107 / GSE159603 / "
            "GSE61332. 'read3p' reproduces the legacy sense-read assumption and is wrong for those data."
        ),
    )

    umi_group = parser.add_argument_group(
        'UMI / molecule counting',
        'OFF by default -- opt in per chemistry. The NET-seq randomer sits at the read 5\' end, i.e. ON TOP '
        'of the RNA 3\' end, so the dedup key is the CORRECTED 3\' end x strand x UMI (never the '
        'soft-clip-corrected 5\' end, which here is non-genomic).'
    )
    umi_group.add_argument(
        '--umi-length',
        type=int,
        default=0,
        help=f'UMI (randomer) length; {NETSEQ_UMI_LENGTH} for Churchman/Couvillion NET-seq. 0 = no UMI.',
    )
    umi_group.add_argument(
        '--umi-source',
        choices=list(NETSEQ_UMI_SOURCES),
        default='read5p',
        help=("Where to find the UMI. 'read5p' slices it off the read as sequenced -- requires the UNTRIMMED "
              "build (on a hexamer-trimmed BAM the randomer is gone from SEQ). 'name' parses a "
              "umi_tools-style suffix on the read name, which requires re-extraction from FASTQ."),
    )
    umi_group.add_argument(
        '--dedup',
        action='store_true',
        help='Emit molecule-level output alongside the read-level output (requires --umi-length).',
    )
    umi_group.add_argument(
        '--umi-clustering',
        choices=['exact', 'directional', 'components'],
        default='exact',
        help=("How to group UMIs within a position. 'exact' (default) is the correct choice for a 6-nt UMI: "
              "with U=4096 the distance-1 graph is dense at deep positions, so directional clustering "
              "collapses giant components instead of correcting errors. The others are for measurement."),
    )
    umi_group.add_argument(
        '--umi-edit-distance',
        type=int,
        default=0,
        help='Max edit distance for non-exact clustering (default 0).',
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output',
    )

    parser.set_defaults(func=run_netseq)
