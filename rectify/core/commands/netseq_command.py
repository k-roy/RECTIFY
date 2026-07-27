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

import sys
from pathlib import Path
from typing import List, Optional

from ..bam.netseq_bam_processor import (
    NETSEQ_RNA3P_CHOICES,
    NETSEQ_RNA3P_DEFAULT,
    process_netseq_bam,
    aggregate_positions,
)
from ..netseq.netseq_umi import (
    NETSEQ_UMI_LENGTH,
    NETSEQ_UMI_SOURCES,
    iter_netseq_fragments,
    select_netseq_molecules,
    summarize_positions,
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

    # Load genome (needed for deconvolution)
    genome = None
    if not args.no_deconvolution:
        print(f"\nLoading genome: {genome_path}")
        genome = load_genome(str(genome_path))
        print(f"  Loaded {len(genome)} chromosomes")

    # Set up exclusion regions
    exclusion_detector = None
    exclude_rdna = not getattr(args, 'include_rdna', False)
    exclude_pol3 = not getattr(args, 'include_pol3', False)

    if exclude_rdna or exclude_pol3:
        print("\nSetting up exclusion regions...")
        exclusion_detector = ExclusionRegionDetector(flanking_bp=args.pol3_flanking)

        # Load from GFF if provided
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

        # Process BAM
        print("\n1. Processing BAM file...")
        records = list(process_netseq_bam(
            input_path,
            exclusion_detector=exclusion_detector,
            min_mapq=args.min_mapq if hasattr(args, 'min_mapq') else 0,
            min_a_fraction=0.8,
            max_reads=args.max_reads if hasattr(args, 'max_reads') else None,
            show_progress=True,
            rna3p_at=rna3p_at,
        ))

        if not records:
            print("  Warning: No records generated!")
            continue

        # Aggregate positions
        print("\n2. Aggregating positions...")
        raw_counts = aggregate_positions(iter(records))
        print(f"  Unique positions: {len(raw_counts):,}")

        # Deconvolution
        deconv_counts = None
        if not args.no_deconvolution and genome:
            print("\n3. Applying deconvolution...")
            deconv_counts, deconv_results = deconvolve_all_regions(
                raw_counts,
                genome,
                min_downstream_a=args.min_atract_length if hasattr(args, 'min_atract_length') else 3,
                show_progress=True,
            )
        else:
            print("\n3. Skipping deconvolution (--no-deconvolution)")

        # Export outputs
        print("\n4. Exporting results...")
        output_formats = args.output_format if hasattr(args, 'output_format') else ['parquet', 'bedgraph']
        normalize_rpm = not getattr(args, 'no_rpm_normalize', False)

        outputs = export_netseq_results(
            records=records,
            raw_counts=raw_counts,
            deconv_counts=deconv_counts,
            output_dir=output_dir,
            sample_name=sample_name,
            export_parquet='parquet' in output_formats,
            export_bedgraph='bedgraph' in output_formats,
            export_bigwig='bigwig' in output_formats,
            normalize_rpm=normalize_rpm,
        )

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
            )

    print(f"\n{'='*60}")
    print("NET-seq processing complete!")
    print(f"Output directory: {output_dir}")
    print(f"{'='*60}")

    return 0


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
) -> None:
    """Write the UMI molecule-level outputs for one sample.

    Three files, and all three matter:
      ``<sample>_umi_stats.json``      duplication rate + family-size + within-family edit histograms --
                                       the calibration outputs, not decoration.
      ``<sample>_molecules.tsv``       per 3'-end position: reads, distinct UMIs (k), occupancy-corrected
                                       molecules (m_hat) and the SATURATED flag.
      ``<sample>_molecule_ids.tsv``    one row per molecule (representative read -> molecule id, family size).
    """
    import json

    print("\n5. Molecule-level (UMI) counting...")
    print(f"   umi-length={umi_length} umi-source={umi_source} clustering={clustering} "
          f"edit-distance={edit_distance}")
    stats = UmiDedupStats()
    fragments = list(iter_netseq_fragments(
        str(input_path), umi_length=umi_length, umi_source=umi_source,
        exclusion_detector=exclusion_detector, min_mapq=min_mapq,
        rna3p_at=rna3p_at, max_reads=max_reads, stats=stats,
    ))
    if not fragments:
        print("   Warning: no reads carried a usable UMI -- molecule track not written")
        return
    keepers, family_size = select_netseq_molecules(
        fragments, edit_distance=edit_distance, clustering=clustering, stats=stats,
    )
    positions = summarize_positions(fragments, umi_length)
    n_sat = sum(1 for p in positions if p.saturated)
    sum_reads = sum(p.reads for p in positions)
    sum_k = sum(p.distinct_umis for p in positions)
    sum_mhat = sum(p.molecules_corrected for p in positions
                   if p.molecules_corrected != float('inf'))

    stats_path = output_dir / f"{sample_name}_umi_stats.json"
    payload = stats.as_dict()
    payload.update({
        "umi_length": umi_length, "umi_source": umi_source, "clustering": clustering,
        "edit_distance": edit_distance, "rna3p_at": rna3p_at,
        "n_positions": len(positions), "n_saturated_positions": n_sat,
        "saturation_threshold_k_over_U": 0.5,
        "sum_reads": sum_reads, "sum_distinct_umis": sum_k,
        "sum_occupancy_corrected_molecules": round(sum_mhat, 1),
        "occupancy_inflation": round(sum_mhat / sum_k, 4) if sum_k else None,
    })
    stats_path.write_text(json.dumps(payload, indent=1))

    mol_path = output_dir / f"{sample_name}_molecules.tsv"
    with open(mol_path, "w") as fh:
        fh.write("chrom\tstrand\tposition\treads\tdistinct_umis\tmolecules_corrected\tsaturated\n")
        for p in positions:
            mh = "inf" if p.molecules_corrected == float('inf') else f"{p.molecules_corrected:.2f}"
            fh.write(f"{p.contig}\t{p.strand}\t{p.position}\t{p.reads}\t{p.distinct_umis}\t"
                     f"{mh}\t{int(p.saturated)}\n")

    ids_path = output_dir / f"{sample_name}_molecule_ids.tsv"
    with open(ids_path, "w") as fh:
        fh.write("representative_read\tmolecule_id\tfamily_size\n")
        for qname, mol_id in keepers.items():
            fh.write(f"{qname}\t{mol_id}\t{family_size.get(mol_id, 1)}\n")

    print(f"   reads with UMI: {stats.n_input_fragments:,}   molecules: {stats.n_molecules:,}   "
          f"duplication rate: {100 * stats.duplication_rate:.1f}%")
    print(f"   positions: {len(positions):,}   SATURATED (k/U>0.5, fall back to reads there): {n_sat:,}")
    for name, path in (("umi_stats", stats_path), ("molecules", mol_path), ("molecule_ids", ids_path)):
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

    # Annotation for exclusion
    parser.add_argument(
        '--gff',
        type=Path,
        help='GFF annotation file for exclusion region detection',
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
        help='Output formats to generate (default: parquet bedgraph)',
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
