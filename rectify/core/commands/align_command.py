"""
Align Command for RECTIFY.

Runs multi-aligner alignment pipeline with consensus selection.

Supported aligners:
- minimap2: Fast seed-and-chain baseline
- mapPacBio: BBTools long-read aligner with splice-aware mode
- gapmm2: minimap2 wrapper with terminal exon refinement

Consensus selection:
- Runs all enabled aligners in parallel
- Compares alignments per read
- Selects best based on junction quality (canonical sites, annotation) and 5' rescue
- Outputs single rectified BAM with best alignment per read
- Note: 3' false junctions are handled by walk back correction, not consensus

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# The paired-end short-read COMPASS panel (the 111-adjudication set): the aligner
# set that `--short-read --read2` expands to. STAR and HISAT2 each appear twice —
# once in their canonical default mode and once in a non-canonical mode
# (--scoreGapNoncan 0 / --pen-noncansplice 0) — so the consensus pool contains
# both a canonical-biased and a non-canonical-tolerant call per read.
#
# Single source of truth: `run-all` (run/stages.py) imports this so its panel can
# never drift from `align`'s. Order is significant only for logging/reporting.
COMPASS_PE_ALIGNERS = [
    'bbmap',
    'STAR_default',
    'STAR_noncanonical',
    'HISAT2_default',
    'HISAT2_noncanonical',
    'magicblast',
    'gsnap',
]

# Single-end COMPASS subset: the default for TruSeq-style short reads without
# --read2. STAR/HISAT2 run single-end; magicblast/gsnap are PE-only and join
# only with --read2. QuantSeq-class 3'-end libraries (--dT-primed-cDNA) use
# bbmap + bwa instead (splice recall matters less than 3'-end placement there).
COMPASS_SE_ALIGNERS = [
    'bbmap',
    'STAR_default',
    'STAR_noncanonical',
    'HISAT2_default',
    'HISAT2_noncanonical',
]


def _commit_indexed_bam(temp_bam: Path, final_bam: Path, index_runner) -> None:
    """Index a temporary BAM before atomically replacing the final BAM and BAI."""
    temp_bai = Path(str(temp_bam) + '.bai')
    final_bai = Path(str(final_bam) + '.bai')
    index_runner(['samtools', 'index', str(temp_bam)], check=True)
    temp_bam.replace(final_bam)
    temp_bai.replace(final_bai)


def create_align_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Create align subcommand parser."""
    parser = subparsers.add_parser(
        'align',
        help='Align reads using multiple aligners with junction annotation support',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument(
        'reads',
        type=Path,
        help='Input FASTQ file (or FASTQ.GZ)'
    )

    parser.add_argument(
        '--genome',
        type=Path,
        default=None,
        help='Reference genome FASTA file. Optional when --Scer or --organism is set.'
    )

    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for BAM files'
    )

    from rectify.data import add_organism_args
    add_organism_args(parser)

    # Protocol flag
    parser.add_argument(
        '--short-read',
        dest='short_read',
        action='store_true',
        default=False,
        help=(
            'Input is short-read data (Illumina/Aviti ≤150 bp). When set, "all" '
            'expands by protocol: TruSeq-style RNA-seq (no --dT-primed-cDNA) '
            'uses the COMPASS splice-aware panel — bbmap + STAR×2 + HISAT2×2 '
            'single-end, plus magicblast + gsnap when paired (--read2); '
            "QuantSeq-class 3'-end libraries (--dT-primed-cDNA) use "
            'bbmap + bwa. Ignored if --aligners is specified explicitly.'
        ),
    )
    parser.add_argument(
        '--dT-primed-cDNA',
        dest='dT_primed_cDNA',
        action='store_true',
        default=False,
        help=(
            "Short-read input is a dT-primed 3'-end library (QuantSeq REV "
            'etc.). With --short-read, "all" selects the 3\'-end panel '
            '(bbmap + bwa) instead of the TruSeq COMPASS panel.'
        ),
    )
    parser.add_argument(
        '-2', '--read2',
        dest='read2',
        type=Path,
        default=None,
        help=(
            'Mate-2 (R2) FASTQ for paired-end short-read alignment. The COMPASS '
            'panel aligners (STAR/HISAT2/magicblast/gsnap) and bbmap run paired '
            'when this is set. Use the paired chunk FASTQs from '
            '`rectify split --read2`.'
        ),
    )
    parser.add_argument(
        '--read-length',
        dest='read_length',
        type=int,
        default=150,
        help='Read length for STAR sjdbOverhang / index selection (default: 150).',
    )

    # Aligner selection
    aligner_group = parser.add_argument_group('Aligner selection')
    aligner_group.add_argument(
        '--aligners',
        nargs='+',
        choices=['minimap2', 'mapPacBio', 'gapmm2', 'bbmap', 'bwa',
                 'winnowmap2', 'minisplice_mm2',
                 'STAR_default', 'STAR_noncanonical', 'HISAT2_default',
                 'HISAT2_noncanonical', 'magicblast', 'gsnap',
                 'all', 'none'],
        default=['all'],
        help=(
            'Aligners to run. "all" = minimap2 (long-read, default; the overhang resolver '
            'is added as the default junction aligner — see --junction-aligners); '
            'mapPacBio and gapmm2 are opt-in arms, listed explicitly. '
            'with --short-read "all" = the COMPASS panel (single-end subset; '
            'full panel with --read2), or bbmap + bwa when --dT-primed-cDNA '
            'is set. '
            '"winnowmap2" and "minisplice_mm2" are opt-in extras (not in "all"); '
            'winnowmap2 requires meryl on PATH; minisplice_mm2 requires --minisplice-model. '
            'Use "none" to run only --junction-aligners (deSALT/uLTRA). (default: all)'
        )
    )

    aligner_group.add_argument(
        '--junction-aligners',
        nargs='+',
        choices=['uLTRA', 'deSALT', 'gmap', 'overhang_resolver'],
        default=None,
        metavar='ALIGNER',
        help=(
            'Splice-aware aligners to add to the consensus pool '
            '(choices: uLTRA, deSALT, gmap, overhang_resolver). '
            'DEFAULT when omitted: uLTRA + deSALT + overhang_resolver on '
            'long-read runs (matches run-all; resolver per the planning/720 '
            'ADOPT verdict), nothing under --short-read. Pass '
            'an explicit list to override; pass only uLTRA/deSALT/gmap to '
            'drop the resolver. uLTRA requires '
            '--annotation; gmap requires a pre-built db (see --gmap-db). '
            'overhang_resolver is not an external aligner: it re-places '
            'terminal soft clips of the minimap2 arm across canonical '
            'junctions under an information bound (planning/641/644), and its '
            'output SUBSTITUTES the minimap2 arm downstream (one correct arm, '
            'not two — planning/669). Requires minimap2 in --aligners. '
            'RECOMMENDED for general use (measured: +26.5k annotated junctions '
            'recovered per 900k cDNA reads, 31 reads harmed, 0 impossible '
            'junctions — planning/720). '
            # NOTE: argparse runs every help string through %-formatting, so a literal
            # percent sign MUST be written '%%'. A bare '%' here raised
            # "TypeError: %o format: an integer is required, not dict" on
            # `rectify align --help` -- the help was unreachable, not merely ugly.
            # Independently found and fixed on both sides of this merge (planning
            # 830 G2 / 831 / 832, and ISSUE-014) — same bug, same fix.
            # tests/test_cli_help_all_subcommands.py pins it for every subcommand.
            'LIMITATION: it recovers 99%% of the annotated junctions mapPacBio '
            'finds but only ~35%% of the NON-CANONICAL ones, because its '
            'candidates come from a GT/AG-class splice-site index and a '
            'non-canonical junction has no entry there (planning/721). For '
            'non-canonical discovery work (upf1D, prp18D, cryptic splicing) it '
            'does NOT substitute for mapPacBio.'
        )
    )

    aligner_group.add_argument(
        '--no-junction-aligners',
        dest='junction_aligners',
        action='store_const',
        const=[],
        help=(
            'Disable all junction aligners, including the default '
            'overhang_resolver post-pass (run only --aligners).'
        )
    )

    aligner_group.add_argument(
        '--resolver-acceptor-classes',
        choices=['canonical', 'prp18'],
        default='canonical',
        help=(
            "Acceptor candidate classes for the overhang_resolver. "
            "'canonical' (default) = AG-class only, the planning/720-measured "
            "configuration. 'prp18' additionally enumerates the alternative "
            "3'SS classes measured in Roy et al. 2023 NAR (gkad968): BG "
            "(TG/CG/GG) + non-G HAU (AT). Opt-in for splicing missions "
            "(upf1D, prp18D/prp18-AA): published utilized alt-3'SS become "
            "enumerable 47%%->88%% (prp18) / 62%%->85%% (upf1D-only) at "
            "~x4.8 acceptor candidate density (planning/722b). Pair with "
            "non-canonical discovery settings (arb_grammar off) where "
            "appropriate."
        )
    )

    aligner_group.add_argument(
        '--resolver-atac',
        dest='resolver_atac',
        action='store_true',
        default=False,
        help=(
            "Also enumerate AT-AC introns in the overhang_resolver, as a PAIRED "
            "class (AT donor <-> AC acceptor only; never AT..AG or GT..AC). Yeast "
            "splices AT-AC through its major spliceosome (Talkish et al. 2019 "
            "PLoS Genet 15:e1008249, SUT635); in human, AT-AC is the U12-type "
            "minor-spliceosome class. Ranked below GT..AG / GC..AG at equal "
            "score. Default off = the planning/720-measured candidate space."
        )
    )

    aligner_group.add_argument(
        '--require-aligners',
        dest='require_aligners',
        action='store_true',
        help=(
            'Exit non-zero if any requested aligner produced no BAM. By '
            'default a missing binary (or an unmet precondition such as uLTRA '
            'without --annotation) is logged as a warning and the run '
            'continues, so a requested 3-aligner panel can silently become a '
            '2-aligner one and still exit 0 — which passes any acceptance gate '
            'that only checks the exit code. Use this in production/array '
            'jobs. The end-of-run DROPPED-ALIGNER summary is emitted either '
            'way; this flag only decides whether it is fatal.'
        )
    )

    aligner_group.add_argument(
        '--trust-existing-bams',
        dest='trust_existing_bams',
        action='store_true',
        help=(
            'Reuse existing per-aligner BAMs even when their provenance '
            'sidecar (rectify SHA / aligner version) does not match the '
            'current run. The checkpoint logic already consumed this flag '
            '(see bam_provenance.py) but only run-all exposed it; needed on '
            'align for reruns that add a panel arm (e.g. overhang_resolver) '
            'across a code update without re-paying the alignment.'
        )
    )

    aligner_group.add_argument(
        '--no-consensus',
        action='store_true',
        help='Skip consensus selection, output separate BAMs per aligner'
    )

    aligner_group.add_argument(
        '--chimeric-consensus',
        action='store_true',
        default=False,
        help=(
            'Use chimeric consensus selection: independently pick the best aligner '
            'for each read segment, then assemble a chimeric CIGAR from the winners. '
            'Experimental — requires further validation before enabling by default.'
        )
    )

    aligner_group.add_argument(
        '--junction-pool-max-intron-len', type=int, default=0, metavar='BP',
        help=(
            'Maximum intron length (nt) for a non-annotated junction to enter the '
            "candidate-junction pool used by consensus selection's 5' soft-clip "
            'rescue. 0 = no limit (default); 3000 suits S. cerevisiae. Mirrors '
            '`rectify consensus --junction-pool-max-intron-len`. No effect with '
            '--no-consensus.'
        ),
    )

    aligner_group.add_argument(
        '--junction-pool-min-anchor-bp', type=int, default=0, metavar='BP',
        help=(
            'Minimum flanking anchor (nt) for a non-annotated junction to enter the '
            "same candidate-junction pool used by consensus selection's 5' soft-clip "
            'rescue. 0 = off (default); 8 is the validated value. Mirrors '
            '`rectify consensus --junction-pool-min-anchor-bp`. No effect with '
            '--no-consensus.'
        ),
    )

    aligner_group.add_argument(
        '--minimap2-path',
        default='minimap2',
        help='Path to minimap2 executable'
    )

    aligner_group.add_argument(
        '--mapPacBio-path',
        default='mapPacBio.sh',
        help='Path to mapPacBio.sh script'
    )

    aligner_group.add_argument(
        '--bbmap-path',
        default='bbmap.sh',
        help='Path to bbmap.sh script (for short-read alignment)'
    )

    aligner_group.add_argument(
        '--bwa-path',
        default='bwa',
        help='Path to bwa executable (for short-read alignment)'
    )

    aligner_group.add_argument(
        '--mapPacBio-chunks',
        type=int,
        default=1,
        metavar='N',
        help=(
            'Split mapPacBio alignment into N independent chunks for parallel '
            'SLURM array execution. Each chunk processes 1/N of the reads '
            '(interleaved, so read-length distribution is even). '
            'Use with --mapPacBio-chunk-idx K to run chunk K, or without '
            '--mapPacBio-chunk-idx to merge existing chunk BAMs.'
        )
    )

    aligner_group.add_argument(
        '--mapPacBio-chunk-idx',
        type=int,
        default=None,
        metavar='K',
        help=(
            '0-based index of the mapPacBio chunk to run (0 to N-1). '
            'Requires --mapPacBio-chunks N. Omit to trigger chunk-merge mode.'
        )
    )

    aligner_group.add_argument(
        '--gapmm2-path',
        default='gapmm2',
        help='Path to gapmm2 executable'
    )

    aligner_group.add_argument(
        '--ultra-path',
        default='uLTRA',
        help='Path to uLTRA executable'
    )

    aligner_group.add_argument(
        '--desalt-path',
        default='deSALT',
        help='Path to deSALT executable'
    )

    aligner_group.add_argument(
        '--gmap-path',
        default='gmap',
        help='Path to gmap executable'
    )

    aligner_group.add_argument(
        '--gmap-db',
        default=None,
        metavar='DIR',
        help=(
            'Pre-built GMAP database directory (<-D dir>/<-d name>). If omitted, '
            'a gmap_db/<genome_stem> dir adjacent to the genome is used. '
            'Build once with: gmap_build -D <dir> -d <name> <genome.fa>'
        )
    )

    aligner_group.add_argument(
        '--winnowmap-repetitive-kmers',
        default=None,
        metavar='FILE',
        help=(
            'Pre-computed repetitive k-mers file for winnowmap2 (output of '
            '"meryl print greater-than distinct=0.9998 <merylDB>"). '
            'If omitted, meryl is invoked automatically and the result cached '
            'adjacent to the genome.'
        )
    )

    aligner_group.add_argument(
        '--minisplice-model',
        default=None,
        metavar='FILE',
        help=(
            'Path to minisplice model file (e.g. vi2-7k.kan). '
            'Required for minisplice_mm2 unless --minisplice-scores is provided.'
        )
    )

    aligner_group.add_argument(
        '--minisplice-model-cali',
        default=None,
        metavar='FILE',
        help='Optional minisplice calibration file (-c flag for minisplice predict).'
    )

    aligner_group.add_argument(
        '--minisplice-scores',
        default=None,
        metavar='FILE',
        help=(
            'Pre-computed splice site scores TSV from "minisplice predict". '
            'If provided, the predict step is skipped.'
        )
    )

    # Junction annotation
    junction_group = parser.add_argument_group('Junction annotation')
    junction_group.add_argument(
        '--annotation',
        type=Path,
        help='Gene annotation GFF/GTF for junction hints'
    )

    junction_group.add_argument(
        '--junc-bed',
        type=Path,
        help='Pre-computed junction BED file (overrides --annotation)'
    )

    junction_group.add_argument(
        '--junc-bonus',
        type=int,
        default=9,
        help='Bonus score for annotated junctions (minimap2)'
    )

    # Performance options
    perf_group = parser.add_argument_group('Performance options')
    perf_group.add_argument(
        '-t', '--threads',
        type=int,
        default=8,
        help='Number of threads per aligner'
    )

    perf_group.add_argument(
        '--parallel-aligners',
        action='store_true',
        help='Run aligners in parallel (uses more memory)'
    )

    perf_group.add_argument(
        '--keep-checkpoints',
        action='store_true',
        default=False,
        help='When consensus selection runs with a checkpoint dir, retain the '
             'per-batch checkpoint files (consensus_batch_*.bam/.done and '
             'consensus_checkpoint.json) after a SUCCESSFUL run. By default these '
             'dead resume-state files are deleted once the final BAM is written. '
             'Failed/interrupted runs always keep checkpoints for resume.'
    )

    perf_group.add_argument(
        '--emit-cma',
        action='store_true',
        default=False,
        help='Additionally write a compressed-multialign <prefix>.cma.bam next to '
             'the multialigned BAM: the pre-correct per-aligner placements with '
             'read SEQ/QUAL stored ONCE (deduplicated across aligners). Opt-in, '
             'non-destructive (nothing is deleted); a lossless store enabling '
             'add-an-aligner / resume-select at ~1/3 the retained footprint. See '
             'planning/254 and `rectify cma`.'
    )

    perf_group.add_argument(
        '--max-intron',
        type=int,
        default=None,
        metavar='BP',
        help=(
            'Maximum intron size passed to the aligner panel (minimap2 -G, '
            'gapmm2 -i, deSALT -I, uLTRA --max_intron, BBMap, STAR, GMAP). '
            'Default: derived from --annotation as 2x the longest annotated '
            'intron, rounded up to 100 and clamped to [1000, 500000] — for '
            'the bundled S. cerevisiae annotation this derives exactly the '
            'historical 5000. Without an annotation the fallback is 5000. '
            'Pass a value to override (e.g. 500000 or larger for human data '
            'with no annotation).'
        )
    )

    # Output options
    output_group = parser.add_argument_group('Output options')
    output_group.add_argument(
        '--prefix',
        default='',
        help='Output file prefix'
    )

    output_group.add_argument(
        '--keep-sam',
        action='store_true',
        help='Keep intermediate SAM files'
    )

    output_group.add_argument(
        '--sort',
        action='store_true',
        default=True,
        help='Sort output BAM files by coordinate'
    )

    output_group.add_argument(
        '--index',
        action='store_true',
        default=True,
        help='Index output BAM files'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose logging'
    )

    return parser


def run_align(args: argparse.Namespace) -> int:
    """Run align command."""
    from datetime import datetime as _dt_al, timezone as _tz_al
    from time import perf_counter as _perf_al
    _align_started_at = _dt_al.now(_tz_al.utc).isoformat()
    _t_align = _perf_al()
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Validate input
    if not args.reads.exists():
        logger.error(f"Reads file not found: {args.reads}")
        return 1

    if not args.genome.exists():
        logger.error(f"Genome file not found: {args.genome}")
        return 1

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Determine prefix
    prefix = args.prefix if args.prefix else args.reads.stem.replace('.fastq', '').replace('.fq', '')

    # Max intron: when not given explicitly (None = auto), derive it from the
    # annotation — 2x the longest annotated intron, rounded up to 100
    # (derive_max_intron docstring has the rule + bounds). For the bundled
    # yeast annotation this derives exactly the historical constant 5000, so
    # existing cohorts are unchanged; other organisms get an honest cap
    # instead of a yeast constant. No annotation -> fallback 5000.
    if getattr(args, 'max_intron', None) is None:
        from ..align.multi_aligner import derive_max_intron
        args.max_intron = derive_max_intron(
            str(args.annotation) if getattr(args, 'annotation', None) else None)
        logger.info(
            f"max_intron: derived {args.max_intron} bp from "
            f"{'annotation' if getattr(args, 'annotation', None) else 'fallback (no annotation)'}"
        )

    # Expand 'all' to list of aligners, then append any junction-mode aligners
    aligners = list(args.aligners)
    if 'none' in aligners:
        aligners = []
    elif 'all' in aligners:
        if getattr(args, 'short_read', False):
            if getattr(args, 'read2', None):
                # Paired short-read → full COMPASS panel (the 111-adjudication set)
                aligners = list(COMPASS_PE_ALIGNERS)
            elif getattr(args, 'dT_primed_cDNA', False):
                # QuantSeq-class 3'-end library → the dT-primed panel
                aligners = ['bbmap', 'bwa']
            else:
                # TruSeq-style RNA-seq, single-end → COMPASS SE subset
                aligners = list(COMPASS_SE_ALIGNERS)
        else:
            # De-paneled 2026-08-17: mapPacBio/gapmm2 are opt-in (--aligners),
            # the overhang resolver (default junction aligner) fills their role.
            aligners = ['minimap2']

    # Junction aligners. None = not given -> the default long-read set
    # [uLTRA, deSALT, overhang_resolver] (Kevin 2026-08-17: bare align matches
    # run-all — uLTRA/deSALT stay wired in unless measured to add little; the
    # resolver per the planning/720 ADOPT verdict). uLTRA gracefully skips
    # without --annotation; deSALT skips without its binary. Nothing on
    # short-read (no minimap2 arm to resolve). An explicit list (including
    # []) always wins.
    _ja_raw = getattr(args, 'junction_aligners', None)
    _ja_explicit = _ja_raw is not None
    if _ja_explicit:
        junction_aligners = list(_ja_raw)
    else:
        junction_aligners = (
            [] if getattr(args, 'short_read', False)
            else ['uLTRA', 'deSALT', 'overhang_resolver'])
    # overhang_resolver is a post-pass on the finished minimap2 arm, not a
    # dispatchable aligner — pulled out here and run after the aligner loop.
    want_resolver = 'overhang_resolver' in junction_aligners
    for ja in junction_aligners:
        if ja == 'overhang_resolver':
            continue
        if ja not in aligners:
            aligners.append(ja)

    # Import multi-aligner functions
    from ..align.multi_aligner import (
        run_minimap2,
        run_map_pacbio,
        run_gapmm2,
        run_ultra,
        run_desalt,
        run_gmap,
        run_bbmap,
        run_bwa_mem,
        run_winnowmap2,
        run_minisplice_mm2,
        run_star,
        run_hisat2,
        run_magicblast,
        run_gsnap,
        check_aligner_available,
    )

    # Paired-end (short-read COMPASS panel) mate-2 FASTQ + read length for STAR.
    _reads2 = getattr(args, 'read2', None)
    _reads2 = str(_reads2) if _reads2 else None
    _read_length = getattr(args, 'read_length', 150) or 150

    # Generate junction BED if needed
    junc_bed_path = None
    if args.junc_bed:
        junc_bed_path = str(args.junc_bed)
    elif args.annotation:
        from ...utils.junction_bed import generate_junction_bed
        junc_bed_path = str(args.output_dir / f"{prefix}_junctions.bed")
        generate_junction_bed(str(args.annotation), junc_bed_path)
        logger.info(f"Generated junction BED: {junc_bed_path}")

    # Run aligners
    results = {}

    parallel = getattr(args, 'parallel_aligners', False)

    # Two-phase scheduler: mapPacBio is ~10x slower than minimap2/gapmm2.
    # Phase 1: mapPacBio runs alone with all threads.
    # Phase 2: remaining aligners run in parallel with equal thread shares.
    aligner_thread_counts = {a: args.threads for a in aligners}

    def _run_one_aligner(aligner):
        """Run a single aligner and return (aligner, bam_path_or_None)."""
        if aligner == 'minimap2':
            exec_path = args.minimap2_path
        elif aligner == 'mapPacBio':
            exec_path = args.mapPacBio_path
        elif aligner == 'gapmm2':
            exec_path = args.gapmm2_path
        elif aligner == 'uLTRA':
            exec_path = getattr(args, 'ultra_path', 'uLTRA')
        elif aligner == 'deSALT':
            exec_path = getattr(args, 'desalt_path', 'deSALT')
        elif aligner == 'gmap':
            exec_path = getattr(args, 'gmap_path', 'gmap')
        elif aligner == 'bbmap':
            exec_path = getattr(args, 'bbmap_path', 'bbmap.sh')
        elif aligner == 'bwa':
            exec_path = getattr(args, 'bwa_path', 'bwa')
        elif aligner == 'winnowmap2':
            exec_path = 'winnowmap'  # binary is named 'winnowmap'; check_aligner_available handles fallback
        elif aligner == 'minisplice_mm2':
            exec_path = 'minimap2'  # uses system minimap2 with --spsc flag
        elif aligner in ('STAR_default', 'STAR_noncanonical', 'HISAT2_default',
                         'HISAT2_noncanonical', 'magicblast', 'gsnap'):
            # COMPASS short-read panel: binary resolved inside the wrapper
            # (_require_binary). STAR/HISAT2 run single- or paired-end;
            # magicblast/gsnap are paired-end only.
            exec_path = None
            if aligner in ('magicblast', 'gsnap') and not _reads2:
                logger.error(
                    "Aligner %s needs paired reads; pass --read2 R2.fastq.gz "
                    "(or paired chunk FASTQs from `rectify split --read2`). Skipping.",
                    aligner,
                )
                return aligner, None
        else:
            logger.warning(f"Unknown aligner: {aligner}")
            return aligner, None

        if aligner == 'winnowmap2':
            import shutil as _shutil
            if not (_shutil.which('winnowmap') or _shutil.which('winnowmap2')):
                logger.warning("winnowmap not found on PATH, skipping winnowmap2")
                return aligner, None
        elif exec_path is not None and not check_aligner_available(exec_path):
            logger.warning(f"{aligner} not found at {exec_path}, skipping")
            return aligner, None

        import time as _time
        n_threads = aligner_thread_counts[aligner]
        output_bam = args.output_dir / f"{prefix}.{aligner}.bam"

        # Per-aligner BAM checkpoint: skip if final output already exists and
        # is plausibly complete.  deSALT may intentionally produce empty BAMs
        # (crash fallback in run_desalt) — always honour those.  For all other
        # aligners, a file < 2 kB means a prior crash wrote a partial/empty
        # file; re-run rather than propagating a corrupt result.
        if output_bam.exists():
            size = output_bam.stat().st_size
            size_ok = (aligner == 'deSALT' or size > 2000)
            if size_ok:
                # Provenance gate: refuse to reuse an on-disk BAM whose
                # sidecar (rectify SHA / aligner version) doesn't match the
                # current run. --trust-existing-bams bypasses this; see
                # rectify/utils/bam_provenance.py.
                _trust = bool(getattr(args, 'trust_existing_bams', False))
                if _trust:
                    logger.info(f"{aligner} BAM already exists ({size}B), reusing "
                                f"(--trust-existing-bams): {output_bam}")
                    return aligner, str(output_bam)
                from ...utils.bam_provenance import (
                    read_sidecar, matches_strict,
                    expected_provenance_for_aligner, compute_run_provenance,
                )
                _run_prov = getattr(args, '_run_provenance', None)
                if _run_prov is None:
                    import sys as _sys
                    _run_prov = compute_run_provenance(command=_sys.argv)
                expected = expected_provenance_for_aligner(_run_prov, aligner)
                stored = read_sidecar(output_bam)
                ok, reason = matches_strict(stored, expected)
                if ok:
                    logger.info(f"{aligner} BAM already exists ({size}B); "
                                f"provenance matches — reusing: {output_bam}")
                    return aligner, str(output_bam)
                logger.warning(
                    f"{aligner} BAM exists at {output_bam} but provenance "
                    f"check failed ({reason}); re-running alignment. "
                    "Pass --trust-existing-bams to override."
                )
                # Fall through to re-run; the aligner wrapper will overwrite.
            else:
                logger.warning(
                    f"{aligner} BAM exists but is only {size}B — likely from a prior crash; "
                    "re-running alignment"
                )

        logger.info(f"Running {aligner} (threads={n_threads})...")
        _t_aligner = _time.perf_counter()

        try:
            if aligner == 'minimap2':
                run_minimap2(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    annotation_path=str(args.annotation) if args.annotation else None,
                    junc_bonus=args.junc_bonus,
                    cache_dir=str(args.output_dir),
                    max_intron=getattr(args, 'max_intron', 5000),
                )
            elif aligner == 'mapPacBio':
                _n_chunks = getattr(args, 'mapPacBio_chunks', 1) or 1
                _chunk_idx = getattr(args, 'mapPacBio_chunk_idx', None)
                _mpb_out = run_map_pacbio(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    chunk_idx=_chunk_idx,
                    n_chunks=_n_chunks if _n_chunks > 1 else None,
                    max_intron=getattr(args, 'max_intron', 5000),
                )
                # Chunk mode redirects the output to *.chunk_K_of_N.bam and
                # returns that path. Rebind so the coord-sort/index/summary
                # below operate on the artifact that exists — otherwise a
                # SUCCESSFUL chunk ends in "No aligners succeeded" against
                # the merged name (planning/644, chunk-0 smoke).
                if _mpb_out:
                    output_bam = Path(_mpb_out)
            elif aligner == 'gapmm2':
                run_gapmm2(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    max_intron=getattr(args, 'max_intron', 5000),
                )
            elif aligner == 'uLTRA':
                if not args.annotation:
                    logger.warning("uLTRA requires --annotation; skipping")
                    return aligner, None
                run_ultra(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    annotation_path=str(args.annotation),
                    threads=n_threads,
                    max_intron=getattr(args, 'max_intron', 5000),
                )
            elif aligner == 'deSALT':
                run_desalt(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    annotation_path=str(args.annotation) if args.annotation else None,
                    threads=n_threads,
                    max_intron=getattr(args, 'max_intron', 5000),
                )
            elif aligner == 'gmap':
                run_gmap(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    annotation_path=str(args.annotation) if args.annotation else None,
                    threads=n_threads,
                    gmap_db=getattr(args, 'gmap_db', None),
                    gmap_path=exec_path,
                    max_intron=getattr(args, 'max_intron', 5000),
                )
            elif aligner == 'bbmap':
                run_bbmap(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    bbmap_path=exec_path,
                    reads2_path=_reads2,
                )
            elif aligner == 'bwa':
                run_bwa_mem(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    bwa_path=exec_path,
                )
            elif aligner in ('STAR_default', 'STAR_noncanonical'):
                # planning/833 G-8: the COMPASS arms must receive rectify's
                # DERIVED --max-intron, not their own defaults (STAR ~589 kb,
                # HISAT2 200 kb, Magic-BLAST 500 kb, GSNAP 200 kb). Otherwise
                # they emit junctions rectify itself calls impossible, and those
                # candidates enter the consensus pool.
                run_star(
                    reads_path=str(args.reads),
                    reads2_path=_reads2,
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    read_length=_read_length,
                    noncanonical=(aligner == 'STAR_noncanonical'),
                    max_intron=getattr(args, 'max_intron', 5000),
                )
            elif aligner in ('HISAT2_default', 'HISAT2_noncanonical'):
                run_hisat2(
                    reads_path=str(args.reads),
                    reads2_path=_reads2,
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    noncanonical=(aligner == 'HISAT2_noncanonical'),
                    max_intron=getattr(args, 'max_intron', 5000),
                )
            elif aligner == 'magicblast':
                run_magicblast(
                    reads_path=str(args.reads),
                    reads2_path=_reads2,
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=min(n_threads, 12),
                    max_intron=getattr(args, 'max_intron', 5000),
                )
            elif aligner == 'gsnap':
                run_gsnap(
                    reads_path=str(args.reads),
                    reads2_path=_reads2,
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    max_intron=getattr(args, 'max_intron', 5000),
                )
            elif aligner == 'winnowmap2':
                run_winnowmap2(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    repetitive_kmers=getattr(args, 'winnowmap_repetitive_kmers', None),
                    cache_dir=str(args.output_dir),
                    max_intron=getattr(args, 'max_intron', 5000),
                )
            elif aligner == 'minisplice_mm2':
                run_minisplice_mm2(
                    reads_path=str(args.reads),
                    genome_path=str(args.genome),
                    output_bam=str(output_bam),
                    threads=n_threads,
                    model_path=getattr(args, 'minisplice_model', None),
                    model_cali_path=getattr(args, 'minisplice_model_cali', None),
                    splice_scores=getattr(args, 'minisplice_scores', None),
                    cache_dir=str(args.output_dir),
                    annotation_path=str(args.annotation) if args.annotation else None,
                    junc_bonus=args.junc_bonus,
                    max_intron=getattr(args, 'max_intron', 5000),
                )

            _elapsed = _time.perf_counter() - _t_aligner
            logger.info(f"{aligner} complete: {output_bam} [{_elapsed:.1f}s]")
            # Aligner wrappers emit name-sorted BAMs (so consensus can stream
            # across aligners).  Coordinate-sort + index here so the on-disk
            # BAM is directly usable by IGV, pysam.fetch(), and downstream
            # per-aligner-review tooling.  Consensus selection re-name-sorts
            # via _ensure_name_sorted() — coord-sort here is non-destructive.
            try:
                import subprocess as _sp
                _sort_tmp = str(output_bam) + '.coord_sort_tmp'
                _sp.run(
                    ['samtools', 'sort', '-@', str(n_threads),
                     '-o', _sort_tmp, str(output_bam)],
                    check=True, capture_output=True,
                )
                Path(_sort_tmp).replace(output_bam)
                _sp.run(['samtools', 'index', str(output_bam)], check=True,
                        capture_output=True)
            except Exception as _idx_err:
                logger.warning(f"{aligner}: coord-sort/index failed ({_idx_err}); "
                                "downstream correction may be skipped")
            # Drop a provenance sidecar (rectify SHA, aligner name+version,
            # timestamp) so a future run-all reuse gate can verify the BAM was
            # produced by the same tool versions before reusing it. See
            # rectify/utils/bam_provenance.py.
            try:
                from ...utils.bam_provenance import write_sidecar
                _run_prov = getattr(args, '_run_provenance', None)
                if _run_prov is None:
                    from ...utils.bam_provenance import compute_run_provenance
                    import sys as _sys
                    _run_prov = compute_run_provenance(command=_sys.argv)
                write_sidecar(output_bam, _run_prov, aligner_name=aligner)
            except Exception as _prov_err:
                logger.warning(f"{aligner}: failed to write provenance sidecar "
                               f"({_prov_err}); BAM emitted without sidecar")
            return aligner, str(output_bam)

        except Exception as e:
            _elapsed = _time.perf_counter() - _t_aligner
            logger.error(f"{aligner} failed after {_elapsed:.1f}s: {e}")
            return aligner, None

    if parallel:
        from concurrent.futures import ThreadPoolExecutor, as_completed
        # Phase 1: mapPacBio alone with all threads (~10x slower than others).
        if 'mapPacBio' in aligners:
            logger.info(
                f"Running mapPacBio first with all {args.threads} threads "
                f"(two-phase: mapPacBio → rest in parallel)"
            )
            aligner_name, bam_path = _run_one_aligner('mapPacBio')
            results['mapPacBio'] = bam_path
        # Phase 2: remaining aligners with equal thread shares.
        # deSALT runs sequentially after the parallel pool — it crashes with
        # "double free or corruption" when forked inside a multithreaded process.
        remaining = [a for a in aligners if a != 'mapPacBio']
        if remaining:
            per_thread = max(1, args.threads // len(remaining))
            for a in remaining:
                aligner_thread_counts[a] = per_thread

            alloc_summary = ', '.join(f"{a}={per_thread}" for a in remaining)
            logger.info(
                f"Running {len(remaining)} aligners in parallel "
                f"(threads: {alloc_summary}, total≤{args.threads})"
            )

            # deSALT and gmap run sequentially after the parallel pool: deSALT
            # crashes when forked inside a multithreaded process, and gmap is slow
            # + spawns its own worker threads, so a dedicated full-thread pass is
            # both safer and faster than contending in the parallel batch.
            _seq = ('deSALT', 'gmap')
            parallel_batch = [a for a in remaining if a not in _seq]
            sequential_batch = [a for a in remaining if a in _seq]

            with ThreadPoolExecutor(max_workers=max(1, len(parallel_batch))) as pool:
                futures = {pool.submit(_run_one_aligner, a): a for a in parallel_batch}
                for future in as_completed(futures):
                    aligner, bam_path = future.result()
                    results[aligner] = bam_path

            for aligner in sequential_batch:
                aligner_name, bam_path = _run_one_aligner(aligner)
                results[aligner_name] = bam_path
    else:
        for aligner in aligners:
            aligner_name, bam_path = _run_one_aligner(aligner)
            results[aligner_name] = bam_path

    # Dropped-aligner gate. Every `return aligner, None` path above (missing
    # binary, uLTRA without --annotation, a PE aligner without --read2, an
    # unknown name) logs a warning and lets the run continue to exit 0, so a
    # requested 3-aligner panel silently becomes a 2-aligner one and still
    # passes any acceptance gate that checks only the exit code (observed on
    # SCG: deSALT absent from the env → green "3-aligner" run with 2 arms).
    # Emit ONE greppable summary line regardless, and make it fatal under
    # --require-aligners. Same fail-loud principle as the resolver check below.
    dropped = [a for a in aligners if not results.get(a)]
    if dropped:
        logger.error(
            "DROPPED-ALIGNER: %d of %d requested aligners produced no BAM: %s "
            "(produced: %s). Scroll up for the per-aligner reason.",
            len(dropped), len(aligners), ', '.join(dropped),
            ', '.join(a for a in aligners if results.get(a)) or 'none',
        )
        if getattr(args, 'require_aligners', False):
            logger.error(
                "--require-aligners was passed; failing rather than emitting a "
                "partial panel with exit 0."
            )
            return 1

    # Overhang-resolver post-pass (planning/641/644): re-places terminal soft
    # clips of the minimap2 arm across canonical junctions under an
    # information bound. The resolved BAM SUBSTITUTES the minimap2 arm
    # downstream (passthrough-or-rewrite: same read set, strictly refined
    # placements) — it is NOT an additional arm. Carrying both would (a) buy a
    # duplicate `correct` arm over ~98%-identical records (+~1.0× the minimap2
    # correct bill, planning/669 §1) and (b) add near-duplicate tie noise to
    # consensus. The raw minimap2 BAM stays on disk for the delta census
    # (rewritten records carry XB tags). Runs after the aligner loop because
    # it consumes the finished (RN-injected) minimap2 BAM. Fails LOUD when
    # requested but unrunnable — a silent skip would feed downstream
    # junction-pool prescans an unresolved pool with exit 0.
    if want_resolver:
        if not results.get('minimap2'):
            if not _ja_explicit:
                # The resolver arrived by DEFAULT, not by request — a panel
                # deliberately run without minimap2 should not fail on a
                # default it never asked for. Degrade to a loud warning.
                logger.warning(
                    "overhang_resolver (default junction aligner) skipped: "
                    "no successful minimap2 arm in this invocation. Pass "
                    "--junction-aligners overhang_resolver to make this "
                    "fatal instead."
                )
                want_resolver = False
            else:
                logger.error(
                    "--junction-aligners overhang_resolver requires a successful "
                    "minimap2 arm in the same invocation (add minimap2 to "
                    "--aligners)."
                )
                return 1
    if want_resolver:
        from ..align.overhang_resolver import run_overhang_resolver
        resolver_bam = args.output_dir / f"{prefix}.overhang_resolver.bam"
        _base_bam = Path(results['minimap2'])
        if (resolver_bam.exists() and resolver_bam.stat().st_size > 2000
                and resolver_bam.stat().st_mtime > _base_bam.stat().st_mtime):
            # Same idempotent-resume spirit as the per-aligner BAM checkpoint:
            # reuse only an output newer than its input arm.
            logger.info(f"overhang_resolver: reusing existing {resolver_bam}")
        else:
            import time as _time
            _t_res = _time.perf_counter()
            run_overhang_resolver(
                base_bam=str(_base_bam),
                genome_path=str(args.genome),
                output_bam=str(resolver_bam),
                threads=args.threads,
                max_intron=getattr(args, 'max_intron', 5000),
                acceptor_classes=getattr(
                    args, 'resolver_acceptor_classes', 'canonical'),
                atac=getattr(args, 'resolver_atac', False),
            )
            logger.info(
                f"[TIMING] overhang_resolver: {_time.perf_counter() - _t_res:.1f}s"
            )
            # Beta-ledger record (realigner runbook 2026-08-11): per-dataset
            # resolver stats next to the BAM, for attributing user-reported
            # alignment oddities to a stage.
            _stats = getattr(run_overhang_resolver, 'last_stats', None)
            if _stats is not None:
                import json as _json
                _stats_path = args.output_dir / f"{prefix}.overhang_resolver.stats.json"
                try:
                    _stats_path.write_text(_json.dumps(
                        _stats.as_dict() if hasattr(_stats, 'as_dict')
                        else vars(_stats), indent=2))
                except (TypeError, OSError) as _e:
                    logger.warning(f"could not write resolver stats JSON: {_e}")
        logger.info(
            "overhang_resolver: SUBSTITUTING the minimap2 arm with the "
            f"resolved BAM ({resolver_bam}); raw minimap2 BAM retained at "
            f"{_base_bam} for the XB delta census (not passed downstream)."
        )
        results['minimap2'] = str(resolver_bam)

    # Summary of alignment step
    logger.info(f"\nAlignment summary:")
    for aligner, bam_path in results.items():
        status = "SUCCESS" if bam_path else "FAILED"
        logger.info(f"  {aligner}: {status}")
        if bam_path:
            logger.info(f"    Output: {bam_path}")

    # Validate outputs: check each reported BAM actually exists and is non-empty.
    # Aligners (especially deSALT) can fail silently, leaving 0-byte files.
    for aligner, bam_path in list(results.items()):
        if bam_path is None:
            continue
        p = Path(bam_path)
        if not p.exists() or p.stat().st_size == 0:
            logger.warning(
                f"{aligner}: output BAM is missing or empty ({bam_path}); "
                "treating as failed"
            )
            results[aligner] = None

    # Check if we should run consensus selection
    successful_aligners = {k: v for k, v in results.items() if v}

    if not successful_aligners:
        logger.error("No aligners succeeded")
        return 1

    # Skip consensus if only one aligner or --no-consensus
    if len(successful_aligners) == 1 or getattr(args, 'no_consensus', False):
        logger.info("Skipping consensus selection (single aligner or --no-consensus)")
        # Copy single BAM to consensus output path, sort and index
        if len(successful_aligners) == 1:
            single_bam = list(successful_aligners.values())[0]
            multialigned_bam = args.output_dir / f"{prefix}.multialigned.bam"
            import shutil
            import subprocess as _sp
            threads = getattr(args, 'threads', 1)
            sorted_tmp = str(multialigned_bam) + '.sorting_tmp'
            _sp.run(
                ['samtools', 'sort', '-@', str(threads), '-o', str(multialigned_bam), str(single_bam)],
                check=True,
            )
            _sp.run(['samtools', 'index', str(multialigned_bam)], check=True)
            logger.info(f"Single-aligner output (sorted+indexed): {multialigned_bam}")
        _write_align_sidecar(args, prefix, multialigned_bam, _align_started_at, _t_align)
        return 0

    # Run consensus selection
    import time as _time
    _t_consensus_start = _time.perf_counter()
    logger.info(f"\nRunning consensus selection across {len(successful_aligners)} aligners...")

    from ..consensus.consensus import run_consensus_selection, load_annotated_junctions

    # Load genome for consensus scoring
    _t_genome_load = _time.perf_counter()
    genome = {}
    import pysam as pysam_lib
    genome_path = str(args.genome)
    try:
        fasta = pysam_lib.FastaFile(genome_path)
    except (OSError, IOError) as e:
        # Handle gzip (not bgzip) compressed genome — auto-convert
        if genome_path.endswith('.gz'):
            import gzip as _gzip
            import subprocess as _sp
            from shutil import which as _which
            logger.warning(f"pysam cannot open {genome_path}: {e}")
            logger.warning("Attempting gzip→bgzip conversion...")
            raw_path = genome_path[:-3]  # strip .gz
            with _gzip.open(genome_path, 'rb') as _fin, open(raw_path, 'wb') as _fout:
                _fout.write(_fin.read())
            import os as _os
            _os.rename(genome_path, genome_path + '.gzip_bak')
            if _which('bgzip'):
                _sp.run(['bgzip', raw_path], check=True)
                _sp.run(['samtools', 'faidx', genome_path], check=True)
                fasta = pysam_lib.FastaFile(genome_path)
                logger.info(f"Successfully converted genome to bgzip: {genome_path}")
            else:
                # No bgzip — use uncompressed
                fasta = pysam_lib.FastaFile(raw_path)
                logger.info(f"Using uncompressed genome: {raw_path}")
        else:
            raise
    for chrom in fasta.references:
        genome[chrom] = fasta.fetch(chrom)
    fasta.close()
    logger.info(f"[TIMING] Genome load: {_time.perf_counter() - _t_genome_load:.1f}s")

    # Load annotated junctions if annotation provided
    annotated_junctions = None
    if args.annotation:
        _t_junc = _time.perf_counter()
        annotated_junctions = load_annotated_junctions(str(args.annotation))
        logger.info(f"[TIMING] Junction load: {_time.perf_counter() - _t_junc:.1f}s")

    # Run aligner selection → merged multi-aligner BAM (pre-correction)
    multialigned_bam = args.output_dir / f"{prefix}.multialigned.bam"

    try:
        _t_sel = _time.perf_counter()
        use_chimeric = getattr(args, 'chimeric_consensus', False)
        if use_chimeric:
            logger.info("Chimeric consensus enabled (experimental) — segments scored independently")
        # Short-read mode adjudicates novel junctions with the COMPASS published
        # tiebreak order (ungapped > gapped > annotated > shorter-intron); an
        # explicit --tiebreak overrides. Long-read keeps the rectify order.
        _tiebreak = getattr(args, 'tiebreak', None)
        if not _tiebreak:
            _tiebreak = 'compass' if getattr(args, 'short_read', False) else 'rectify'
        stats = run_consensus_selection(
            bam_paths=successful_aligners,
            genome=genome,
            output_bam=str(multialigned_bam),
            annotated_junctions=annotated_junctions,
            use_chimeric=use_chimeric,
            checkpoint_dir=getattr(args, 'checkpoint_dir', None),
            keep_checkpoints=getattr(args, 'keep_checkpoints', False),
            tiebreak=_tiebreak,
            pool_max_intron_len=getattr(args, 'junction_pool_max_intron_len', 0),
            pool_min_anchor_bp=getattr(args, 'junction_pool_min_anchor_bp', 0),
        )
        logger.info(f"[TIMING] Aligner selection: {_time.perf_counter() - _t_sel:.1f}s")

        logger.info(f"\nMulti-aligned (pre-correction) BAM: {multialigned_bam}")
        logger.info(f"  High confidence: {stats['consensus_high']} reads")
        logger.info(f"  5' rescued: {stats['5prime_rescued']} reads")
        logger.info(f"[TIMING] Aligner selection total (incl. genome/junctions): {_time.perf_counter() - _t_consensus_start:.1f}s")

        # Write aligner stats TSV and HTML report alongside the multi-aligned BAM
        try:
            from ..bam.processing_stats import write_consensus_stats_tsv
            _stats_tsv = args.output_dir / f"{prefix}.consensus_aligner_stats.tsv"
            write_consensus_stats_tsv(stats, str(_stats_tsv))
        except Exception as _e:
            logger.warning(f"Could not write consensus stats TSV: {_e}")

        try:
            from ..analyze.summary import generate_consensus_html_report
            _report_html = args.output_dir / f"{prefix}.consensus_report.html"
            generate_consensus_html_report(stats, str(_report_html), sample_name=prefix)
            logger.info(f"Consensus report: {_report_html}")
        except Exception as _e:
            logger.warning(f"Could not write consensus HTML report: {_e}")

    except Exception as e:
        logger.error(f"Consensus selection failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Opt-in compressed-multialign artifact (non-destructive; planning/254).
    # Wrapped so a CMA failure can never fail an otherwise-successful alignment.
    if getattr(args, 'emit_cma', False):
        try:
            from ..multialign import build_cma_from_bams, validate_cma
            _cma_path = args.output_dir / f"{prefix}.cma.bam"
            _cma_stats = build_cma_from_bams(
                successful_aligners, str(_cma_path),
                panel=list(successful_aligners), genome=genome,
            )
            _cma_probs = validate_cma(str(_cma_path))
            logger.info(
                f"CMA emitted: {_cma_path} ({_cma_stats['reads']} reads, "
                f"{_cma_stats['records']} records; "
                f"validate={'OK' if not _cma_probs else str(len(_cma_probs)) + ' problems'})"
            )
        except Exception as _e:
            logger.warning(f"--emit-cma failed (non-fatal): {_e}")

    # Add MD tags via samtools calmd (required for indel correction and
    # alignment identity calculation downstream).
    logger.info("Adding MD tags with samtools calmd...")
    try:
        import subprocess as _sp
        calmd_bam = args.output_dir / f"{prefix}.multialigned.md.bam"
        calmd_cmd = [
            'samtools', 'calmd', '-b',
            str(multialigned_bam),
            str(args.genome),
        ]
        with open(str(calmd_bam), 'wb') as fh_out:
            result = _sp.run(calmd_cmd, stdout=fh_out, stderr=_sp.PIPE)
        if result.returncode == 0 and calmd_bam.stat().st_size > 0:
            _commit_indexed_bam(calmd_bam, multialigned_bam, _sp.run)
            logger.info("  MD tags added successfully")
        else:
            logger.warning(
                f"  samtools calmd failed (rc={result.returncode}); "
                "proceeding without MD tags"
            )
            if calmd_bam.exists():
                calmd_bam.unlink()
    except Exception as e:
        logger.warning(f"  samtools calmd error: {e}; proceeding without MD tags")

    _write_align_sidecar(args, prefix, multialigned_bam, _align_started_at, _t_align)
    return 0


def run(args: argparse.Namespace):
    """Entry point for CLI."""
    sys.exit(run_align(args))


def _write_align_sidecar(args, prefix, out_bam, started_at, t_start, aligners=None):
    """Record which rectify (and which aligners) produced this BAM.

    `align` REWRITES coordinates -- every downstream 3'-end call inherits them --
    yet until 2026-08-02 it wrote no provenance at all. Non-fatal on failure,
    matching the convention in the other stages.
    """
    try:
        from datetime import datetime as _dt, timezone as _tz
        from time import perf_counter as _perf
        from rectify.core.provenance import ProvenanceRecord, write_stage_sidecar
        from rectify.utils.version import get_rectify_git_sha as _get_sha
        _stats = {"wall_seconds": _perf() - t_start}
        if aligners:
            _stats["aligners"] = list(aligners)
            try:
                from rectify.utils.bam_provenance import get_aligner_version
                _stats["aligner_versions"] = {
                    a: (get_aligner_version(a) or "unknown") for a in aligners
                }
            except Exception:
                pass
        _inputs = {"reads": args.reads, "genome": args.genome}
        _rec = ProvenanceRecord.from_components(
            stage="align",
            sample_id=prefix,
            sample_output_dir=args.output_dir,
            started_at=started_at,
            completed_at=_dt.now(_tz.utc).isoformat(),
            exit_status=0,
            inputs=_inputs,
            outputs={"bam": out_bam},
            stats=_stats,
            argv=sys.argv,
            rectify_git_sha=_get_sha(),
        )
        write_stage_sidecar(_rec, sample_output=args.output_dir)
    except Exception as exc:
        logger.warning("Failed to write align sidecar: %s", exc)
