#!/usr/bin/env python3
"""
Per-read 3' end correction for RECTIFY.

This module owns the ``correct_read_3prime`` workhorse and the small helpers it
needs (3'/5' position extraction, NET-seq loader factory, pool junction index).
Higher-level orchestration lives in peer modules:

  - ``regions``      — coverage-gap detection and region planning
  - ``parallel``     — worker dispatch, parallel/streaming pipelines
  - ``output``       — TSV writer, summary reports
  - ``variant_scan`` — variant-aware first-pass scan
  - ``bam_writer``   — CIGAR surgery and BAM writers
  - ``processing_stats`` — ProcessingStats accumulator

Per-read correction module order (within ``correct_read_3prime``):

  2E pre-pass  — filter poly(A)-artifact junctions before 5' rescue
  2F           — 3'SS truncation rescue (Cat3: extend 5' end to nearest splice donor)
  2A           — AG mispriming filter (dT-primed cDNA only; pass ``--dT-primed-cDNA``)
  2B           — poly(A) tail detection / trimming (all protocols with sequenced poly-A)
  2C           — indel artifact correction (genome-vs-query HP comparison)
  2D           — variant-aware homopolymer rescue (enabled by default)
  2G           — soft-clip rescue at 3' homopolymer boundaries (all protocols)
  2E           — poly-A walk-back: genome-aware 3' end correction (all protocols)
  NET-seq      — fractional position refinement (optional, NET-seq experiments only)

Note: DRS poly(A) pre-trimming (Step 0, ``drs_trim_command.trim_drs_bam_polya``) and
poly(A) tail restoration (Step 4, ``restore_polya_command.restore_polya_softclips``)
are orchestrated by ``run_command._run_single_sample``, NOT by this module.

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, List, Optional, Tuple
import logging

import pysam

from ..polya import atract_detector
from ..polya import ag_mispriming
from ..polya import polya_trimmer
from ..correct import indel_corrector
from ..netseq import netseq_refiner
from ..correct.indel_corrector import VariantAwareHomopolymerRescue
from ..correct.protocols.quantseq_rev import walkback_quantseq_rev
from ..correct.protocols.ont_cdna import (
    POLYA_SOURCE_TRIM as _POLYA_SOURCE_TRIM,
    resolve_rna_strand as _resolve_ont_cdna_strand,
    three_prime_position as _ont_cdna_3prime,
    trim_stage_tail_length as _trim_stage_tail_length,
    walkback_ont_cdna,
)
from ..correct.walkback import (
    APPLIED_WALKBACK as _APPLIED_WALKBACK_READGENOME,
    walkback_drs_full,
)
from ...utils.genome import (
    get_chrom_sequence,
    load_genome,
    reverse_complement,
    standardize_chrom_name,
)
from ..polya.polya_model import PolyAModel
from ...utils.alignment import (
    extract_junctions_simple,
    extract_soft_clips,
    format_junctions_string,
)
from ..splice.splice_aware_5prime import (
    rescue_3ss_truncation as _rescue_3ss,
    _get_5prime_softclip_len,
    MAX_SS_SHIFT as _MAX_SS_SHIFT,
    DEFAULT_JUNCTION_PROXIMITY_BP as _JUNCTION_PROXIMITY_BP,
    DEFAULT_PEEL_MAX_BP as _PEEL_MAX_BP,
)
from ..splice import false_junction_filter as _fjf

# Backwards-compat re-exports for callers that still reach in via this module.
from .processing_stats import ProcessingStats, write_stats_tsv, generate_stats_report  # noqa: F401
from .bam_writer import (  # noqa: F401  (re-exported)
    _load_corrections_from_tsv,
    _decode_eq_seq_inplace,
    write_dual_bam,
    write_corrected_bam,
    write_softclipped_bam,
    write_polya_trimmed_bam,
)
from .read_edits import (  # noqa: F401  (re-exported)
    clip_read_to_corrected_3prime,
    softclip_read_to_corrected_3prime,
    extend_read_5prime_for_junction_rescue,
    fix_homopolymer_mismatches,
    realign_exon_blocks,
)
from .bedgraph_writers import (  # noqa: F401  (re-exported)
    write_netseq_assigned_bedgraph,
    write_corrected_reads_bedgraph,
)
from ..position_index import write_position_index  # noqa: F401  (re-exported)

# ``process_bam_file`` (defined below) writes its result list via the shared
# TSV emitter in ``output``.  All other peer-module dispatch entry points
# (region planning, parallel/streaming workers, variant scan) are imported
# directly from their home modules by their callers.
from .output import consensus_tag_fields, write_output_tsv

logger = logging.getLogger(__name__)

# Half-window (bp) for the strand-agnostic pool fetch, EXCLUDING the per-read 5'
# soft-clip (added at fetch time). The pool is indexed as a per-chrom interval
# tree keyed on BOTH splice sites, so the fetch asks "which junction intervals
# come within this window of the read's 5' boundary" — an interval whose NEAR
# site sits by the read overlaps the window regardless of how far its OTHER site
# is (human introns can be >100 kb). That makes the fetch intron-length
# independent, unlike the old intron_start-only bisect (which needed a radius >=
# the intron length to surface a plus-strand junction by its acceptor).
#
# Sized to the rescue's actual reach so no winner is dropped: the baseline
# sequence rescue accepts a junction edge up to junction_proximity_bp + five_clip
# away, and the terminal peel reaches junction_proximity_bp + peel_max_bp; both
# then slide +/-MAX_SS_SHIFT. Derived from the SAME constants the rescue uses
# (imported above) so it cannot drift the way the old magic 10000 did.
_POOL_FETCH_HALF_WINDOW = _JUNCTION_PROXIMITY_BP + _PEEL_MAX_BP + _MAX_SS_SHIFT + 5


def _build_pool_chrom_index(
    pool_junctions: Optional[set],
) -> Optional[Dict[str, 'IntervalTree']]:
    """Build a per-chromosome interval tree for either-site pool lookup.

    Returns a dict mapping each standardised chromosome name to an IntervalTree
    of ``[intron_start, intron_end)`` intervals (each carrying its
    ``(intron_start, intron_end)`` tuple as data).  Returns None if
    pool_junctions is None or empty.

    Keying on the whole intron interval (rather than one coordinate) lets the
    per-read fetch surface a junction whenever EITHER splice site is near the
    read's 5' boundary OR the read's 5' end sits inside the intron — without a
    radius that has to scale with intron length.  Call once before BAM
    processing; pass the result as ``pool_chrom_index``.
    """
    if not pool_junctions:
        return None
    from intervaltree import IntervalTree
    index: Dict[str, 'IntervalTree'] = {}
    for entry in pool_junctions:
        chrom_raw, intron_start, intron_end = entry
        if intron_end <= intron_start:
            continue  # IntervalTree requires begin < end; skip degenerate entries
        chrom = standardize_chrom_name(chrom_raw)
        if chrom not in index:
            index[chrom] = IntervalTree()
        index[chrom][intron_start:intron_end] = (intron_start, intron_end)
    return index


def _load_netseq(netseq_dir: str) -> 'netseq_refiner.NetseqLoader':
    """
    Load NET-seq data into a NetseqLoader.

    Handles two cases:
    - 'bundled:<organism>': loads pre-deconvolved TSV bundled with the package
    - '<path>': loads BigWig files from a directory

    Args:
        netseq_dir: Either 'bundled:<organism>' or a filesystem path

    Returns:
        Populated NetseqLoader
    """
    loader = netseq_refiner.NetseqLoader()
    if netseq_dir.startswith('bundled:'):
        organism = netseq_dir[len('bundled:'):]
        loader.load_bundled(organism)
    else:
        loader.load_directory(netseq_dir, pattern="*.bw")
    return loader


def get_read_3prime_position(read: pysam.AlignedSegment) -> Tuple[int, str]:
    """
    Get 3' end position and strand from read.

    The 3' end is the cleavage and polyadenylation (CPA) site.

    Args:
        read: pysam AlignedSegment

    Returns:
        Tuple of (position, strand)

    Coordinate Details (0-based):
        Plus strand (+): 3' end = reference_end - 1 (rightmost aligned base)
        Minus strand (-): 3' end = reference_start (leftmost aligned base)
    """
    # Determine strand.
    # NOTE: DRS (Direct RNA Sequencing) reads are sequenced 3'→5', so minimap2
    # aligns them reverse-complemented; is_reverse=True for sense-strand reads.
    # This is corrected upstream by aligning with -uf (forward-strand-only) so
    # strand orientation here matches the RNA, not the sequencing direction.
    strand = '-' if read.is_reverse else '+'

    # Get 3' end position
    if strand == '+':
        ref_end = read.reference_end
        if ref_end is None:
            return None, strand  # Unmapped read — no valid position
        position = ref_end - 1  # 0-based inclusive
    else:
        ref_start = read.reference_start
        if ref_start is None or ref_start < 0:
            return None, strand  # Unmapped read — no valid position
        position = ref_start  # 0-based

    return position, strand


def get_read_5prime_position(read: pysam.AlignedSegment, strand: Optional[str] = None) -> int:
    """
    Get 5' end position from read (transcription start site end).

    The 5' end is the opposite end from the 3' end (CPA site).
    This represents where transcription started for this RNA molecule.

    Args:
        read: pysam AlignedSegment
        strand: Optional strand override. If None, uses read.is_reverse.

    Returns:
        5' end position (0-based genomic coordinate)

    Coordinate Details:
        Plus strand (+): RNA 5'→3' matches genomic left→right
            5' end = reference_start (leftmost = transcription start)

        Minus strand (-): RNA 5'→3' matches genomic right→left
            5' end = reference_end - 1 (rightmost = transcription start)

    Visual Example (plus strand):
        Genomic:  100       110       120       130
                  |---------|---------|---------|
        Read:     [====================]
                  5'                  3'
                  ↑                    ↑
                  start=100           end-1=129
        5' end = 100

    Visual Example (minus strand):
        Genomic:  100       110       120       130
                  |---------|---------|---------|
        Read:     [====================]
                  3'                  5'
                  ↑                    ↑
                  start=100           end-1=129
        5' end = 129
    """
    if strand is None:
        strand = '-' if read.is_reverse else '+'

    if strand == '+':
        return read.reference_start  # Leftmost = 5' for plus strand
    else:
        ref_end = read.reference_end
        if ref_end is None:
            return None  # Unmapped read — no valid position
        return ref_end - 1  # Rightmost = 5' for minus strand


def correct_read_3prime(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    ag_threshold: float = 0.65,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_loader: Optional[netseq_refiner.NetseqLoader] = None,
    variant_aware_rescue: Optional[VariantAwareHomopolymerRescue] = None,
    annotated_junctions: Optional[set] = None,
    pool_chrom_index: Optional[Dict] = None,
    apply_3ss_rescue: bool = True,
    gene_interval_trees: Optional[Dict] = None,
    polya_model: Optional[PolyAModel] = None,
    dt_primed_cDNA: bool = False,
    ont_cDNA: bool = False,
    exclusion_detector: Optional['ExclusionRegionDetector'] = None,
    use_dorado_polya: bool = False,
) -> List[Dict]:
    """
    Apply all corrections to a single read.

    Args:
        read: pysam AlignedSegment
        genome: Genome dict from load_genome()
        apply_atract: Apply A-tract ambiguity detection
        apply_ag_mispriming: Apply AG mispriming screening
        apply_polya_trim: Apply poly(A) tail trimming
        apply_indel_correction: Apply indel artifact correction
        netseq_loader: Optional NET-seq loader for refinement
        variant_aware_rescue: Optional VariantAwareHomopolymerRescue object for
            filtering out likely true variants during homopolymer rescue

    Returns:
        Dict with correction results
    """
    # Decode SAM-spec ``=`` shorthand in SEQ before any module reads
    # query_sequence. Aligners with ``=``-CIGAR emission (minimap2, gapmm2,
    # deSALT, uLTRA) propagate ``=`` chars into SEQ at match positions; any
    # downstream consumer comparing SEQ bytes to genome bytes (walkback,
    # softclip rescue, overcall rescue, indel corrector) would see ``=`` as
    # a literal mismatch and short-circuit. Decoding once at intake makes
    # every consumer ``=``-safe without per-module patches.
    if genome is not None and not read.is_unmapped:
        _decode_eq_seq_inplace(read, genome)

    # Multi-aligner consensus tags (Xa/Xc/Xn/Xt) carried through to the
    # per-read TSV as four separate columns.  Read once, here, ABOVE the
    # chimeric ``try/except KeyError`` below — a KeyError raised inside that
    # block would be swallowed and silently reroute the read.
    _consensus_tags = consensus_tag_fields(read)

    # Skip chimeric reads — they have already been reconstructed by multi_aligner
    # and must not be re-corrected (Xz=1 flag set by chimeric_consensus.py).
    try:
        if read.get_tag('Xz') == 1:
            original_position, strand = get_read_3prime_position(read)
            if original_position is None:
                return []
            if dt_primed_cDNA:
                # Antisense protocol: gene strand is opposite of read strand
                if read.is_reverse:
                    strand = '+'
                    original_position = (read.reference_end - 1) if read.reference_end else original_position
                else:
                    strand = '-'
                    original_position = read.reference_start
            chrom = read.reference_name
            chrom_std = standardize_chrom_name(chrom)
            _strand_evidence = ''
            if ont_cDNA:
                # PCR-cDNA reads occur in both orientations; resolve per read.
                _rna_strand, _strand_evidence = _resolve_ont_cdna_strand(
                    read, gene_interval_trees, chrom=chrom_std
                )
                if _rna_strand is not None:
                    strand = _rna_strand
                    _p = _ont_cdna_3prime(read, strand)
                    if _p is not None:
                        original_position = _p
            five_prime_position = get_read_5prime_position(read, strand)
            if five_prime_position is None:
                return []
            # Extract splice junctions from CIGAR even for chimeric reads —
            # we don't re-correct them, but the N-ops in the multi-aligner
            # stitch are still real and downstream consumers (validation
            # tests, junction aggregation, BED export) need them. Previously
            # this dict omitted junctions/n_junctions entirely, which
            # silently dropped them from the corrected_reads.tsv.
            _chimeric_junctions = extract_junctions_simple(read)
            _chimeric_junctions_str = format_junctions_string(_chimeric_junctions)
            return [{
                'read_id': read.query_name,
                'chrom': chrom_std,
                'strand': strand,
                'original_3prime': original_position,
                'corrected_3prime': original_position,
                'five_prime_position': five_prime_position,
                'five_prime_rescued': False,
                'alignment_start': read.reference_start,
                'alignment_end': read.reference_end,
                'ambiguity_min': original_position,
                'ambiguity_max': original_position,
                'ambiguity_range': 0,
                'correction_category': None,
                'correction_applied': [],
                'polya_length': 0,
                'aligned_a_length': 0,
                'soft_clip_a_length': 0,
                'atract_ambiguity': None,
                'fraction': 1.0,
                'mapq': read.mapping_quality,
                'confidence': 'chimeric',
                'qc_flags': '',
                'pt_tag': None,
                'polya_score': None,
                'polya_source': 'none',
                'junctions': _chimeric_junctions,
                'junctions_str': _chimeric_junctions_str,
                'n_junctions': len(_chimeric_junctions),
                'strand_evidence': _strand_evidence,
                **_consensus_tags,
            }]
    except KeyError:
        pass

    # Get original position and strand
    original_position, strand = get_read_3prime_position(read)
    if original_position is None:
        logger.warning(f"Could not compute 3' position for read {read.query_name}, skipping")
        return []
    if dt_primed_cDNA:
        # Antisense protocol (e.g. QuantSeq REV): reads are antisense to mRNA.
        # Gene strand is the opposite of read strand; CPA is at the opposite end.
        if read.is_reverse:
            # is_reverse=True → gene is plus-strand → CPA at rightmost aligned base
            strand = '+'
            original_position = (read.reference_end - 1) if read.reference_end else original_position
        else:
            # is_reverse=False → gene is minus-strand → CPA at leftmost aligned base
            strand = '-'
            original_position = read.reference_start
    chrom = read.reference_name

    # Standardize chromosome name
    chrom_std = standardize_chrom_name(chrom)

    # ONT PCR-cDNA: the library is double-stranded, so reads arrive in BOTH
    # orientations and `is_reverse` alone does not give the RNA strand (gene
    # strand is '+' iff read_is_sense XOR is_reverse).  Resolve per read from
    # the trim-stage `ro` tail-evidence tag, falling back to the maximally-
    # overlapping annotated gene, and take the terminus the resolved strand
    # implies.  Unresolved reads keep the DRS-rule default and are labelled
    # `unassigned` in strand_evidence so downstream can drop them — never
    # guessed silently.  See rectify.core.correct.protocols.ont_cdna.
    strand_evidence = ''
    if ont_cDNA:
        rna_strand, strand_evidence = _resolve_ont_cdna_strand(
            read, gene_interval_trees, chrom=chrom_std
        )
        if rna_strand is not None:
            strand = rna_strand
            _p = _ont_cdna_3prime(read, strand)
            if _p is not None:
                original_position = _p

    # Get 5' end position (transcription start site)
    five_prime_position = get_read_5prime_position(read, strand)
    if five_prime_position is None:
        logger.warning(f"Could not compute 5' position for read {read.query_name}, skipping")
        return []

    # Read Xj BAM tag (set by consensus.py when 5' soft-clip was rescued during consensus)
    try:
        five_prime_rescued = read.get_tag('Xj') == 1
    except KeyError:
        five_prime_rescued = False
    _five_prime_exon_cigar = ''      # set by 3'SS rescue when local alignment succeeds
    _five_prime_intron_clip_pos = -1 # set for intronic-snap reads (Case 4) only
    _five_prime_upstream_trim = 0    # set by 3'SS rescue equivalence-extension (cat3 - strand)
    _reanchor_clip_len = 0           # set by 3'SS rescue reanchor pre-pass (mpb 5'-edge cluster)

    # Module 2E (pre-pass): filter poly(A)-artifact junctions before 5' rescue
    # so they are never used as 3'SS rescue candidates.
    _genome_ref = _fjf._GenomeDictReference(genome)
    _real_junctions, _artifact_analyses = _fjf.filter_polya_artifact_junctions(read, _genome_ref, strand)

    # Module 2F: 3'SS truncation rescue (post-consensus).
    # Corrects five_prime_position for reads truncated or soft-clipped at the
    # exon 2 / 3' splice site boundary.
    #
    # rRNA / Pol III exclusion: rDNA, tRNAs, SNR6, RDN5, etc. are dropped
    # downstream by `analyze`, yet they dominate correct-stage wall time because
    # these are the highest-coverage loci and 3'SS rescue runs its O(n*m) HP-edit
    # distance per read there. Skip the rescue for reads mapping into an excluded
    # region — they are not alt-spliced mRNA and need no 3'SS rescue. (The read
    # still gets all other corrections + a normal corrected_3end row; only the
    # expensive, pointless rescue is skipped.)
    _excluded_region = bool(
        exclusion_detector is not None
        and exclusion_detector.is_excluded(chrom_std, five_prime_position)
    )
    if apply_3ss_rescue and not _excluded_region and (annotated_junctions or _real_junctions or pool_chrom_index):
        _ss_junctions: set = set()
        if annotated_junctions:
            _ss_junctions.update(annotated_junctions)
        for _js, _je in _real_junctions:
            _ss_junctions.add((chrom_std, _js, _je))
        # Add pool junctions (novel + annotated from prescan) whose intron comes
        # within the rescue's reach of this read's 5' boundary. Either-site /
        # containment lookup via the per-chrom interval tree: a junction is
        # surfaced whenever its NEAR splice site sits by the read (human introns
        # can put the far site >100 kb away, so a one-coordinate radius would miss
        # it) OR the 5' end falls inside the intron (Case-4 snap). The window adds
        # the read's 5' soft-clip because the baseline sequence rescue reaches a
        # junction edge up to junction_proximity_bp + five_clip away. The
        # downstream rescue applies the exact per-loop gate; this is only a safe
        # superset that keeps the per-read candidate count small.
        if pool_chrom_index:
            _tree = pool_chrom_index.get(chrom_std)
            if _tree:
                _w = _POOL_FETCH_HALF_WINDOW + _get_5prime_softclip_len(read)
                for _iv in _tree.overlap(five_prime_position - _w,
                                         five_prime_position + _w + 1):
                    _ss_junctions.add((chrom_std, _iv.data[0], _iv.data[1]))
        if _ss_junctions:
            _3ss_result = _rescue_3ss(read, genome, _ss_junctions, strand)
            if _3ss_result['rescued']:
                five_prime_rescued = True
                five_prime_position = _3ss_result['five_prime_corrected']
                _five_prime_exon_cigar = _3ss_result.get('five_prime_exon_cigar', '')
                _five_prime_upstream_trim = int(_3ss_result.get('five_prime_upstream_trim', 0) or 0)
                _reanchor_clip_len = int(_3ss_result.get('reanchor_clip_len', 0) or 0)
                # five_prime_intron_clip_pos ("icp") = the exon-2-side intron
                # boundary, recorded whenever the alignment's 5' end sits inside
                # the rescued intron — Case 4 (intronic_snap) and Cases 1/2 alike,
                # with or without a 5' soft clip.
                #
                # This field DRIVES BAM surgery; it is not a diagnostic. In
                # bam_writer.apply_corrected_edits_to_read it selects which of the
                # three mutually exclusive 5' helpers runs: `extend` only when the
                # N-op it would draw lands on icp, otherwise
                # `reroute_intronic_tail_5prime_via_junction` (which uses icp as
                # its clip_boundary), otherwise `softclip_intronic_tail_5prime`.
                # The rescue sizes `five_prime_exon_cigar` from the intronic query
                # run for exactly these reads so the reroute's query-span check
                # can pass (ISSUE-002). An earlier comment here claimed "no BAM
                # surgery is performed for these reads"; that stopped being true
                # when the reroute/softclip helpers landed, and believing it is
                # how the writer came to draw introns the rescue never chose.
                #
                # NOTE (latent, not fixed here): this is computed from the
                # RESTORED, pre-reanchor read, while the rescue body's
                # align_5prime is post-reanchor. For a read the reanchor pre-pass
                # materially changed, icp and the exon CIGAR can therefore be
                # sized against different geometries.
                _rj = _3ss_result.get('rescued_junction')
                if _rj:
                    _, _intron_start, _intron_end = _rj
                    if strand == '-':
                        _align_5p = (read.reference_end - 1) if read.reference_end else -1
                        _in_intron = _intron_start <= _align_5p < _intron_end
                        _five_prime_intron_clip_pos = _intron_start if _in_intron else -1
                    else:
                        _align_5p = read.reference_start
                        _in_intron = _intron_start <= _align_5p < _intron_end
                        _five_prime_intron_clip_pos = _intron_end if _in_intron else -1
                else:
                    _five_prime_intron_clip_pos = -1

    # Extract splice junctions from CIGAR
    junctions = extract_junctions_simple(read)
    # Junctions the aligner itself called. The writer-verdict block below may
    # retract a junction the RESCUE added; it must never drop one of these.
    _cigar_junctions = set(junctions)
    # When 5' rescue fired and located a real annotated junction, the read's
    # raw CIGAR doesn't yet carry an N-op (the surgery happens later in
    # bam_writer's extend_read_5prime_for_junction_rescue). For TSV reporting
    # purposes, include the rescued junction so n_junctions reflects the
    # post-rescue state — this is what downstream consumers (validation tests,
    # junction aggregation, BED export) expect.
    if five_prime_rescued:
        _rj_post = _3ss_result.get('rescued_junction') if '_3ss_result' in locals() else None
        if _rj_post:
            _, _r_start, _r_end = _rj_post
            _junc_tuple = (_r_start, _r_end)
            if _junc_tuple not in junctions:
                junctions = list(junctions) + [_junc_tuple]
    junctions_str = format_junctions_string(junctions)

    # Extract soft clips (returns list of dicts with 'side' and 'length' keys)
    soft_clips_list = extract_soft_clips(read)
    left_clip_length = 0
    right_clip_length = 0
    for clip in soft_clips_list:
        if clip['side'] == 'left':
            left_clip_length = clip['length']
        elif clip['side'] == 'right':
            right_clip_length = clip['length']

    # Determine which soft clip is 5' vs 3' based on strand
    if strand == '+':
        five_prime_soft_clip_len = left_clip_length
        three_prime_soft_clip_len = right_clip_length
    else:
        five_prime_soft_clip_len = right_clip_length
        three_prime_soft_clip_len = left_clip_length

    # Reanchor pre-pass: when rescue_3ss_truncation's reanchor materially
    # modified the CIGAR (e.g. mapPacBio's `1X 2= 7I` 5' edge collapsed into
    # a `10S`), the rescue's exon_cigar was computed against the post-reanchor
    # geometry. The TSV's five_prime_soft_clip_length must reflect the POST-
    # reanchor leading-S length so bam_writer's extend_read_5prime_for_junction_rescue
    # gate (soft_clip_len > 0) fires and the rescue's exon_cigar actually
    # replaces the soft-clip in the final BAM. Without this propagation, mpb-
    # style reads pay the per-base soft-clip penalty (1.0/bp) in HP-ED scoring
    # for a clip that the rescue was supposed to eliminate.
    if _reanchor_clip_len > 0:
        five_prime_soft_clip_len = _reanchor_clip_len

    # --- Writer verdict: does the rescue this row found actually reach the BAM?
    # The corrected TSV is written a whole pipeline stage BEFORE any BAM writer
    # runs, so a row could advertise `five_prime_rescued=1` with an exon CIGAR
    # and a junction that the writer then refuses to draw — the ISSUE-006
    # canonical-destination guard, the ISSUE-007 cut-N refusal, or the icp gate
    # sending the read to a reroute that declines. A live TSV consumer (the
    # browser) would draw a junction the BAM does not contain.
    #
    # Ask the writer itself, on a throwaway copy of the read, and downgrade the
    # row when it says no. This runs the same function on the same read after the
    # same pre-passes, so it cannot drift from the writer; and it is
    # self-consistent, because a downgraded row makes the writer skip the surgery
    # entirely — which is the outcome the refusal already produced.
    #
    # Placed after five_prime_soft_clip_len is final (the reanchor propagation
    # above), because the writer reads that value as `five_prime_soft_clip`.
    _five_prime_rescue_refused = ''
    if (not five_prime_rescued and '_3ss_result' in locals()
            and _3ss_result.get('displaced_canonical_refused')):
        # Not a writer verdict: the RESCUE declined every candidate that would
        # have destroyed a canonical junction the aligner already called
        # (splice_aware_5prime, hold-out read 34625d8e). Recorded in the same
        # column so one field answers "why is this read not rescued".
        _five_prime_rescue_refused = 'would_displace_canonical'
    if five_prime_rescued and genome is not None:
        from .bam_writer import predict_5prime_rescue_refusal as _predict_5p_refusal
        _five_prime_rescue_refused = _predict_5p_refusal(
            read,
            {
                'five_prime_rescued': True,
                'five_prime_position': five_prime_position,
                'five_prime_soft_clip': five_prime_soft_clip_len,
                'five_prime_exon_cigar': _five_prime_exon_cigar,
                'five_prime_upstream_trim': _five_prime_upstream_trim,
                'five_prime_intron_clip_pos': _five_prime_intron_clip_pos,
                'reanchor_clip_len': _reanchor_clip_len,
                'strand': strand,
            },
            genome,
        )
        if _five_prime_rescue_refused:
            # Always retract the JUNCTION — that is the claim the BAM does not
            # back. Never touch a junction the aligner's own CIGAR carries.
            # five_prime_position is deliberately KEPT: the rescue's estimate of
            # where the 5' end lies is still the module's best answer.
            from .bam_writer import REFUSAL_SOFTCLIP_ONLY as _REFUSAL_SC_ONLY
            _rj_ref = _3ss_result.get('rescued_junction') if '_3ss_result' in locals() else None
            if _rj_ref:
                _junc_drop = (_rj_ref[1], _rj_ref[2])
                if _junc_drop in junctions and _junc_drop not in _cigar_junctions:
                    junctions = [_j for _j in junctions if _j != _junc_drop]
                    junctions_str = format_junctions_string(junctions)
            _five_prime_exon_cigar = ''
            if _five_prime_rescue_refused != _REFUSAL_SC_ONLY:
                # No surgery at all: clear the rest so the writer skips it, which
                # reproduces exactly what the refusal already produced.
                five_prime_rescued = False
                _five_prime_intron_clip_pos = -1
            # For softclipped_no_junction the intronic bases ARE hidden at the
            # true acceptor, so five_prime_rescued and the icp must survive or the
            # writer would skip that surgery and leave the bases mapped inside the
            # intron — a worse BAM than the one this verdict describes.

    # Extract 3' soft-clip sequence for poly(A) model scoring.
    # Iterate the clip list rather than recomputing — seq is already present.
    _three_prime_clip_seq = None
    for _clip in soft_clips_list:
        if (strand == '+' and _clip['side'] == 'right') or \
                (strand == '-' and _clip['side'] == 'left'):
            _three_prime_clip_seq = _clip.get('seq')
            break

    # Read pt:i BAM tag (dorado signal-level poly(A) length estimate, if present).
    # Never absent from DRS reads processed by dorado >= 0.5; None otherwise.
    try:
        _pt_tag = read.get_tag('pt')
    except KeyError:
        _pt_tag = None

    # Score 3' soft-clip against the poly(A) model when one is provided.
    # For minus-strand reads the soft-clip is in genomic (anti-sense) orientation;
    # reverse-complement to RNA (5'→3') orientation before scoring.
    _polya_score = None
    if polya_model is not None and _three_prime_clip_seq:
        _seq_to_score = _three_prime_clip_seq
        if strand == '-':
            _seq_to_score = reverse_complement(_three_prime_clip_seq)
        _polya_score = round(
            polya_model.score_sequence(_seq_to_score.upper())['confidence'], 4
        )

    # Initialize result (use standardized chromosome name for output)
    result = {
        'read_id': read.query_name,
        'chrom': chrom_std,  # Use Roman numeral format (chrI, chrII, etc.)
        'strand': strand,
        # How `strand` was determined.  Only ONT PCR-cDNA populates this (the
        # other protocols have a fixed, protocol-level strand rule); one of
        # polyA_3p / polyT_5p / gene_overlap / unassigned.  `unassigned` marks a
        # read whose orientation could not be established — it keeps the DRS
        # default so the row is still emitted, but it must be filtered out of
        # any strand-sensitive analysis.
        'strand_evidence': strand_evidence,
        'original_3prime': original_position,
        'corrected_3prime': original_position,
        'five_prime_position': five_prime_position,  # TSS end of the read
        'alignment_start': read.reference_start,  # Leftmost aligned position
        'alignment_end': read.reference_end,  # Rightmost + 1 (exclusive)
        'ambiguity_min': original_position,
        'ambiguity_max': original_position,
        'ambiguity_range': 0,
        'correction_applied': [],
        'qc_flags': [],
        'confidence': 'high',
        # New fields for unified record
        'junctions': junctions,  # List of (start, end) tuples
        'junctions_str': junctions_str,  # Semicolon-separated string
        'n_junctions': len(junctions),
        'five_prime_soft_clip_length': five_prime_soft_clip_len,
        'three_prime_soft_clip_length': three_prime_soft_clip_len,
        'mapq': read.mapping_quality,
        'five_prime_rescued': five_prime_rescued,
        'five_prime_exon_cigar': _five_prime_exon_cigar,
        # -1 = not applicable; ≥0 = genomic boundary to hard-clip intronic tail to
        'five_prime_intron_clip_pos': _five_prime_intron_clip_pos,
        # Cat3 equivalence-extension: k bases the BAM writer should trim from
        # the end of the upstream M before applying the rescue surgery
        'five_prime_upstream_trim': _five_prime_upstream_trim,
        # Length of the leading soft-clip produced by the reanchor pre-pass
        # (rectify/core/bam/read_edits.py:reanchor_5prime_for_rescue); 0 when
        # reanchor did not materially modify the CIGAR. > 0 signals bam_writer
        # to apply the same reanchor before realign + 5'-rescue surgery.
        'reanchor_clip_len': _reanchor_clip_len,
        # '' when the 5' rescue reached the corrected BAM (or none was found);
        # otherwise the bam_writer REFUSAL_* token saying why it did not, with
        # five_prime_rescued/exon_cigar/intron_clip_pos already downgraded and the
        # rescued junction dropped from `junctions`.
        'five_prime_rescue_refused': _five_prime_rescue_refused,
        # Cat2 soft-clip rescue fields (v2.9.1) — populated if Module 2G fires
        'sc_homopolymer_extension': 0,   # under-called homopolymer bases → D op
        'sc_rescued_seq': '',            # non-poly-A bases matched to ref → M op
        'sc_original_softclip_len': 0,   # original 3' soft-clip length
        # Poly(A) evidence (v2.9.0)
        'pt_tag': _pt_tag,
        'polya_score': _polya_score,
        'polya_source': (
            'pt_tag' if _pt_tag is not None else
            ('model' if _polya_score is not None else 'none')
        ),
        # consensus_aligner / consensus_confidence / consensus_n_agree /
        # consensus_tied — '' when the input BAM carried no Xa/Xc/Xn/Xt.
        **_consensus_tags,
    }

    # Per-read gene attribution via read body overlap
    gene_id = None
    if gene_interval_trees is not None:
        try:
            from ..analyze.gene_attribution import compute_read_gene_attribution
            gene_ids = compute_read_gene_attribution(
                read, gene_interval_trees, chrom=chrom_std,
                # ONT PCR-cDNA: use the RESOLVED RNA strand, not is_reverse.
                # Without this the antisense half is looked up on the wrong
                # strand -- half get no gene_id at all and the rest are
                # attributed to a gene on the opposite strand.
                rna_strand=(strand if ont_cDNA else None))
            gene_id = gene_ids[0] if gene_ids else None
        except Exception as _e:
            logger.warning("Gene attribution failed for read %s: %s", read.query_name, _e)
    result['gene_id'] = gene_id

    current_position = original_position

    # Module 1: A-tract ambiguity (always applied by default)
    if apply_atract:
        atract_result = atract_detector.calculate_atract_ambiguity(
            genome, chrom_std, current_position, strand, downstream_bp=10
        )

        result['ambiguity_min'] = atract_result['ambiguity_min']
        result['ambiguity_max'] = atract_result['ambiguity_max']
        result['ambiguity_range'] = atract_result['ambiguity_range']
        result['downstream_a_count'] = atract_result.get('downstream_a_count')

        if atract_result['has_ambiguity']:
            result['correction_applied'].append('atract_ambiguity')

    # Module 2A: AG mispriming (for oligo-dT technologies)
    if apply_ag_mispriming:
        ag_result = ag_mispriming.screen_ag_mispriming(
            genome, chrom_std, current_position, strand,
            threshold=ag_threshold,
        )

        if ag_result['is_likely_misprimed']:
            qc_flag = ag_mispriming.get_ag_qc_flag(ag_result)
            result['qc_flags'].append(qc_flag)
            result['ag_content'] = ag_result['ag_content']
            result['ag_score'] = ag_result['ag_score']
            if result['confidence'] == 'high':
                result['confidence'] = 'low'

    # Module 2B: Poly(A) tail detection (when poly(A) is sequenced)
    # NOTE: This does NOT correct positions - it only measures poly(A) tail length.
    # Position correction for A-tract ambiguity is already handled above.
    polya_shift = 0
    if apply_polya_trim:
        # Pass atract_result so polya_trimmer can include aligned A's in the count
        atract_for_polya = {
            'tract_length': result.get('ambiguity_range', 0),  # Approximate aligned A's
        }
        polya_result = polya_trimmer.trim_polya_from_read(
            read, strand, atract_result=atract_for_polya,
            use_dorado_polya=use_dorado_polya
        )

        # Always record poly(A) length (even if 0) for completeness
        result['polya_length'] = polya_result['polya_length']
        result['aligned_a_length'] = polya_result['aligned_a_length']
        result['soft_clip_a_length'] = polya_result['soft_clip_a_length']

        # ONT PCR-cDNA: prefer the TRIM-STAGE tail measurement when present.
        # `trim-cdna-polya` removes the tail before alignment, so the read
        # measured just above no longer contains it and reports ~0 for 96% of
        # sense and 95% of antisense reads (planning/550). The `pl` tag carries
        # the pre-trim length, orientation-aware. Absent tag => fall back to the
        # post-alignment value rather than assume a zero-length tail.
        if ont_cDNA:
            _trim_tail = _trim_stage_tail_length(read)
            if _trim_tail is not None:
                result['polya_length'] = _trim_tail
                result['polya_source'] = _POLYA_SOURCE_TRIM
        # Dorado pt:i estimate (None when absent) — recorded for comparison and,
        # under --use-dorado-polya, already folded into polya_length above.
        result['dorado_polya_length'] = polya_result['dorado_polya_length']

        # NOTE: We no longer add 'polya_trim' to correction_applied because
        # soft-clips don't affect genomic position. The only position correction
        # comes from A-tract ambiguity detection (atract_ambiguity).

    # Module 2C: Indel artifact correction (when poly(A) is sequenced)
    indel_shift = 0
    if apply_indel_correction:
        indel_result = indel_corrector.correct_indels_from_read(read, strand, genome)

        if indel_result['has_artifacts']:
            result['correction_applied'].append('indel_correction')
            result['n_indel_artifacts'] = len(indel_result['artifacts'])
            indel_shift = indel_result['correction_bp']

            # Update current position
            current_position = indel_result['corrected_3prime']

    # Module 2D: Variant-aware homopolymer rescue (optional)
    # This filters out positions where high read frequency suggests true variant
    if variant_aware_rescue is not None:
        var_rescue_result = variant_aware_rescue.rescue_with_variant_filter(
            read, strand, genome, end='3prime'
        )
        if var_rescue_result:
            if var_rescue_result['variant_check'] == 'RESCUED':
                result['correction_applied'].append('homopolymer_rescue')
                result['homopolymer_rescue_bases'] = var_rescue_result['rescued_bases']
                current_position = var_rescue_result['corrected_pos']
            elif var_rescue_result['variant_check'] == 'SKIPPED_LIKELY_VARIANT':
                result['qc_flags'].append('LIKELY_VARIANT')
                result['variant_confidence'] = var_rescue_result['variant_confidence']

    # NEW-063 / NEW-068: Detect hard-clipped 3' end before Modules 2G and 2E.
    # Hard-clip means the 3' sequence is absent; both softclip rescue and poly-A
    # walkback would explore an unanchored region and must be skipped.
    _has_3prime_hardclip = False
    if read.cigartuples:
        _cigar_ops = read.cigartuples
        if strand == '+' and _cigar_ops[-1][0] == 5:
            _has_3prime_hardclip = True
        elif strand == '-' and _cigar_ops[0][0] == 5:
            _has_3prime_hardclip = True

    # Module 2G: Soft-clip rescue at homopolymer boundaries.
    # When the basecaller under-calls a homopolymer the aligner may end the
    # alignment inside the run and leave downstream matching bases as a 3' soft-
    # clip.  rescue_softclip_at_homopolymer() detects this pattern (homopolymer
    # at the alignment boundary + soft-clipped bases that match the reference
    # downstream) and extends the 3' end OUTWARD.
    #
    # This correction is mutually exclusive with poly-A walk-back (Module 2E),
    # which moves the 3' end INWARD.  If this module fires we skip walk-back.
    #
    # Minimum soft-clip length guard: single-base soft-clips are common
    # basecalling artifacts that should not override polya_walkback.  Only
    # attempt rescue when the 3' soft-clip is ≥ 3 bases.
    _3prime_sc_len = 0
    if read.cigartuples:
        if strand == '+' and read.cigartuples[-1][0] == 4:
            _3prime_sc_len = read.cigartuples[-1][1]
        elif strand == '-' and read.cigartuples[0][0] == 4:
            _3prime_sc_len = read.cigartuples[0][1]

    softclip_rescue_applied = False
    if genome and not _has_3prime_hardclip and _3prime_sc_len >= 3:
        _sc_result = indel_corrector.rescue_softclip_at_homopolymer(
            read, strand, genome, end='3prime'
        )
        if _sc_result is not None:
            result['correction_applied'].append('softclip_rescue')
            current_position = _sc_result['corrected_pos']
            softclip_rescue_applied = True
            # Store CIGAR surgery metadata so bam_writer can extend the alignment.
            result['sc_homopolymer_extension'] = _sc_result['homopolymer_extension']
            result['sc_rescued_seq'] = _sc_result['rescued_seq']
            result['sc_original_softclip_len'] = _3prime_sc_len
            # The soft-clip bases matched the reference exactly: the corrected
            # position is definitive.  Reset ambiguity to zero so that NET-seq
            # refinement (which handles positional ambiguity, not sequence
            # evidence) does not override the rescue.
            result['ambiguity_min'] = current_position
            result['ambiguity_max'] = current_position
            result['ambiguity_range'] = 0

    # Module 2G.5: Over-call terminal-match rescue (cat1_plus_1/2 pattern).
    # When the basecaller OVER-called the poly-A tail length, the 3' soft-clip's
    # inner bases are all pA-stop-base and the OUTERMOST base is a non-stop-base
    # that matches genome at the first non-stop-base position past the genomic
    # HP. Convert the soft-clip into `{hp_ext}D {overcall_count}I 1=` and advance
    # the corrected 3' end by hp_ext+1.
    #
    # Mutually exclusive with softclip_rescue (Module 2G) and polya_walkback
    # (Module 2E): all three operate on the same 3' soft-clip evidence. Fire only
    # when softclip_rescue did not already claim the soft-clip.
    overcall_rescue_applied = False
    if (
        genome
        and not _has_3prime_hardclip
        and not softclip_rescue_applied
        and _3prime_sc_len >= 1
    ):
        _oc_result = indel_corrector.rescue_overcall_terminal_match(
            read, strand, genome,
        )
        if _oc_result is not None:
            result['correction_applied'].append('overcall_rescue')
            current_position = _oc_result['corrected_pos']
            overcall_rescue_applied = True
            # Bam_writer reads these to apply the CIGAR surgery.
            result['oc_overcall_count']         = _oc_result['overcall_count']
            result['oc_homopolymer_extension']  = _oc_result['homopolymer_extension']
            result['oc_terminal_base']          = _oc_result['terminal_match']
            # Sequence evidence is definitive: collapse ambiguity to the rescue
            # target so NET-seq refinement does not override.
            result['ambiguity_min'] = current_position
            result['ambiguity_max'] = current_position
            result['ambiguity_range'] = 0

    # Module 2E: Poly-A walk-back — genome-aware correction for reads where the
    # poly-A tail aligned to a short genomic A-run, shifting the apparent 3' end
    # into the run. find_polya_boundary() walks backwards from the mapped 3' end
    # to the first position where genome and read agree on a non-A (+ strand) or
    # non-T (- strand) base — the true CPA site. Always applied when a genome is
    # available; no minimum poly-A length threshold.
    #
    # The walk-back span [corrected, original] is genuine positional ambiguity:
    # NET-seq refinement will use it to select the most likely CPA within that
    # window. ambiguity_range is updated to the max of the atract window and the
    # walk-back distance so the NET-seq refinement guard (range > 0) fires.
    # NEW-068: Skip poly-A walkback when softclip rescue (Module 2G) already
    # corrected this end — the two modules move in opposite directions.
    #
    # NEW-075: For QuantSeq REV (dt_primed_cDNA=True) the legacy genome-only
    # find_polya_boundary fires for only ~8% of reads because its 4-A homopolymer
    # gate filters out the dominant artifact (internal-priming over short genomic
    # A-runs). Route those reads through the read-vs-reference walkback
    # (rectify.core.correct.protocols.quantseq_rev). The new walkback also
    # handles the V-primer tip artifact (terminal G over a genomic A-run).
    polya_walkback_applied = False
    if (
        genome
        and not _has_3prime_hardclip
        and not softclip_rescue_applied
        and not overcall_rescue_applied
    ):
        if ont_cDNA:
            # ONT PCR-cDNA: same walkback core as DRS, but side/stop_base come
            # from the per-read resolved RNA strand rather than from is_reverse.
            _chrom_seq, _chrom_seq_key = get_chrom_sequence(genome, chrom)
            wb = walkback_ont_cdna(
                read, _chrom_seq, strand,
                artifact_analyses=_artifact_analyses,
            ) if _chrom_seq else None
            if wb is not None:
                result['correction_applied'].append('polya_walkback')
                polya_walkback_applied = True
                wb_bp = wb['correction_bp']
                current_position = wb['corrected_pos']
                if strand == '+':
                    result['ambiguity_min'] = current_position
                    result['ambiguity_max'] = original_position
                else:
                    result['ambiguity_min'] = original_position
                    result['ambiguity_max'] = current_position
                result['ambiguity_range'] = max(result['ambiguity_range'], abs(wb_bp))
        elif dt_primed_cDNA:
            _chrom_seq, _chrom_seq_key = get_chrom_sequence(genome, chrom)
            if _chrom_seq:
                _orig_wb, _corr_wb, _applied_wb, _gene_strand_wb = walkback_quantseq_rev(
                    read, _chrom_seq,
                    artifact_analyses=_artifact_analyses,
                )
                if _applied_wb == _APPLIED_WALKBACK_READGENOME:
                    result['correction_applied'].append('polya_walkback_readgenome')
                    polya_walkback_applied = True
                    current_position = _corr_wb
                    if strand == '+':
                        result['ambiguity_min'] = current_position
                        result['ambiguity_max'] = original_position
                    else:
                        result['ambiguity_min'] = original_position
                        result['ambiguity_max'] = current_position
                    wb_bp = abs(original_position - current_position)
                    result['ambiguity_range'] = max(result['ambiguity_range'], wb_bp)
        else:
            # DRS production walkback. The shared protocol-agnostic core is
            # walkback_3prime_guarded; walkback_drs_full applies the DRS
            # strand→side/stop_base mapping and returns the same dict shape
            # as the legacy find_polya_boundary (which now delegates to it).
            _chrom_seq, _chrom_seq_key = get_chrom_sequence(genome, chrom)
            if _chrom_seq:
                wb = walkback_drs_full(
                    read, _chrom_seq,
                    artifact_analyses=_artifact_analyses,
                )
            else:
                wb = None
            if wb is not None:
                result['correction_applied'].append('polya_walkback')
                polya_walkback_applied = True
                wb_bp = wb['correction_bp']
                current_position = wb['corrected_pos']
                # Set ambiguity window to the full walk-back span
                if strand == '+':
                    result['ambiguity_min'] = current_position
                    result['ambiguity_max'] = original_position
                else:
                    result['ambiguity_min'] = original_position
                    result['ambiguity_max'] = current_position
                result['ambiguity_range'] = max(result['ambiguity_range'], wb_bp)

    # NEW-061: Prevent corrected positions from landing inside artifact N op spans.
    # Two cases handled:
    #
    # Case A — walkback landed STRICTLY inside the N: snap to N boundary and
    # collapse ambiguity to zero so NET-seq has nothing to refine.
    # Strict < on left: current_position == junction_start means the read abuts
    # the N correctly and no snap is needed.
    #
    # Case B — walkback landed before the N but the ambiguity window extends into
    # or beyond it: clip the window so NET-seq cannot place signal inside the N.
    # For + strand, cap ambiguity_max at junction_start - 1.
    # For − strand, cap ambiguity_min at junction_end.
    if polya_walkback_applied and _artifact_analyses:
        for _art in _artifact_analyses:
            if _art.junction_start < current_position < _art.junction_end:
                # Case A: walkback landed strictly inside the N — snap to the N's
                # near edge so all poly-A-aligned M bases are consumed.
                # Plus strand: 3' end is on the right; N is to the left; snap to
                #   junction_start-1 (last M before N).
                # Minus strand: 3' end is on the left; N is to the right; snap to
                #   junction_end (first aligned base after N). junction_start is
                #   the first skipped base and would put corrected_3prime inside
                #   the artifact gap.
                if strand == '+':
                    current_position = _art.junction_start - 1
                else:
                    current_position = _art.junction_end
                result['ambiguity_min'] = current_position
                result['ambiguity_max'] = current_position
                result['ambiguity_range'] = 0
                break
            elif strand == '+' and result['ambiguity_max'] >= _art.junction_start:
                # Case B: walkback landed before the N but the ambiguity window
                # extends into or past it — clip max to junction_start-1.
                result['ambiguity_max'] = _art.junction_start - 1
                result['ambiguity_range'] = result['ambiguity_max'] - result['ambiguity_min']
                break
            elif strand == '-' and result['ambiguity_min'] <= _art.junction_end:
                # Case B (minus): ambiguity window extends into the N — clip min
                # to junction_end so NET-seq stays on the 3' side of the N.
                result['ambiguity_min'] = _art.junction_end
                result['ambiguity_range'] = result['ambiguity_max'] - result['ambiguity_min']
                break

    # NEW-073: Prevent poly-A walkback from crossing real (annotated) splice
    # junctions.  The existing guard (lines above) only covers artifact N-ops;
    # it misses the case where the 3' end lies in the *last exon* and the
    # walk-back retreats past the intron boundary into the preceding intron
    # (which can have long genomic A-runs on the reference).
    #
    # For + strand: the nearest intron whose right edge (intron_end) is <=
    #   original_position is the boundary.  Clip ambiguity_min there.
    # For − strand: the nearest intron whose left edge (intron_start) is >=
    #   original_position is the boundary.  Clip ambiguity_max there.
    if polya_walkback_applied and _real_junctions:
        if strand == '+':
            # All introns that are completely to the left of original_position
            left_introns = [(s, e) for s, e in _real_junctions if e <= original_position]
            if left_introns:
                nearest_intron_end = max(e for s, e in left_introns)
                if current_position < nearest_intron_end:
                    # Walkback crossed into the intron — snap to intron end,
                    # preserve ambiguity only within the valid exonic range.
                    current_position = nearest_intron_end
                    result['ambiguity_min'] = nearest_intron_end
                    result['ambiguity_max'] = original_position
                    result['ambiguity_range'] = original_position - nearest_intron_end
                elif result['ambiguity_min'] < nearest_intron_end:
                    # Ambiguity window extends into the intron — clip min.
                    result['ambiguity_min'] = nearest_intron_end
                    result['ambiguity_range'] = result['ambiguity_max'] - result['ambiguity_min']
        else:  # minus strand
            # All introns that are completely to the right of original_position
            right_introns = [(s, e) for s, e in _real_junctions if s >= original_position]
            if right_introns:
                nearest_intron_start = min(s for s, e in right_introns)
                if current_position > nearest_intron_start:
                    # Walkback crossed into the intron — snap to intron start.
                    current_position = nearest_intron_start
                    result['ambiguity_min'] = original_position
                    result['ambiguity_max'] = nearest_intron_start
                    result['ambiguity_range'] = nearest_intron_start - original_position
                elif result['ambiguity_max'] > nearest_intron_start:
                    # Ambiguity window extends into the intron — clip max.
                    result['ambiguity_max'] = nearest_intron_start
                    result['ambiguity_range'] = result['ambiguity_max'] - result['ambiguity_min']

    # Update corrected position (after all position-moving corrections)
    # Universal 100%-non-A enforcement (AGENT_FIXES [2026-06-14]).
    #
    # RECTIFY's DRS policy left-shifts EVERY corrected 3' end to the first
    # non-A: the pre-trim removes the poly-A tail, so a corrected end parked on
    # a gene-strand stop base (genomic A on +, genomic T on −) is always wrong.
    # The walkback's own _enforce_non_stop_anchor only covers the "walkback
    # returned None" path; this generalizes the guarantee to "ANY module
    # (homopolymer_rescue / indel_correction / atract / raw) left the end on a
    # stop base". Rather than a fresh genomic scan (which would re-cross real
    # introns / large dels), re-anchor via the guard-respecting walkback; if
    # that can't anchor off the stop base, revert the offending extension to
    # the raw 3' end when IT is non-A. This is the single chokepoint that
    # drives the corrected-on-stop-base rate to ~0 (see TestCorrectedEndsAreNonA
    # and the Sumner human-DRS verification).
    if genome:
        _cs_enf, _ = get_chrom_sequence(genome, chrom)
        _stop_enf = 'A' if strand == '+' else 'T'
        if (
            _cs_enf
            and 0 <= current_position < len(_cs_enf)
            and _cs_enf[current_position].upper() == _stop_enf
        ):
            _target = None
            _wb_enf = walkback_drs_full(
                read, _cs_enf, artifact_analyses=_artifact_analyses
            )
            if (
                _wb_enf is not None
                and 0 <= _wb_enf['corrected_pos'] < len(_cs_enf)
                and _cs_enf[_wb_enf['corrected_pos']].upper() != _stop_enf
            ):
                _target = _wb_enf['corrected_pos']
            elif (
                0 <= original_position < len(_cs_enf)
                and _cs_enf[original_position].upper() != _stop_enf
            ):
                # walkback could not anchor off the stop base — the offending
                # module extended past the raw 3' end onto a stop base; revert
                # to the (non-A) raw 3' end.
                _target = original_position
            if _target is not None:
                current_position = _target
                # Revert rescue CIGAR-surgery metadata so the output BAM CIGAR
                # stays coherent with the re-anchored position (the walkback /
                # raw path uses trimming, not CIGAR surgery).
                for _k in (
                    'sc_homopolymer_extension', 'sc_rescued_seq',
                    'sc_original_softclip_len', 'oc_homopolymer_extension',
                    'oc_overcall_count', 'oc_terminal_base',
                    'homopolymer_rescue_bases',
                ):
                    result.pop(_k, None)
                for _t in ('softclip_rescue', 'overcall_rescue',
                           'homopolymer_rescue'):
                    while _t in result['correction_applied']:
                        result['correction_applied'].remove(_t)
                if 'polya_walkback' not in result['correction_applied']:
                    result['correction_applied'].append('polya_walkback')
                # Mark as walkback so the NET-seq poly-A filter (NEW-065) treats
                # the re-anchored position as the safe fallback and collapse the
                # ambiguity window so NET-seq cannot move it back onto a stop base.
                polya_walkback_applied = True
                result['ambiguity_min'] = current_position
                result['ambiguity_max'] = current_position
                result['ambiguity_range'] = 0

    result['corrected_3prime'] = current_position

    # Update ambiguity window for indel/rescue corrections that moved the position
    # but did NOT set the window themselves (i.e. polya_walkback was not applied).
    if current_position != original_position and not polya_walkback_applied:
        # Mirror the walkback window logic (lines 579-584): the NET-seq window
        # spans [corrected_anchor, uncorrected_position] for each strand.
        # Plus strand:  correction moves LEFT  → ambiguity_min = current (left), max = original (right)
        # Minus strand: correction moves RIGHT → ambiguity_min = original (left), max = current (right)
        if strand == '+':
            result['ambiguity_min'] = current_position
            result['ambiguity_max'] = original_position
        else:
            result['ambiguity_min'] = original_position
            result['ambiguity_max'] = current_position
        if result['ambiguity_min'] > result['ambiguity_max']:
            result['ambiguity_min'], result['ambiguity_max'] = (
                result['ambiguity_max'],
                result['ambiguity_min'],
            )
        result['ambiguity_range'] = result['ambiguity_max'] - result['ambiguity_min']

    # NEW-062: record 5' rescue in correction_applied before any early return.
    if result.get('five_prime_rescued'):
        result['correction_applied'].append('five_prime_rescued')

    # Module 3: NET-seq refinement (optional)
    # Bundled TSV data is already deconvolved — skip re-deconvolution.
    # Custom BigWig data has not been deconvolved — apply NNLS.
    use_deconvolution = not getattr(netseq_loader, '_bundled_loaded', False)

    if netseq_loader is not None and result['ambiguity_range'] > 0:
        assignments = netseq_refiner.refine_with_netseq(
            netseq_loader,
            chrom_std,
            result['ambiguity_min'],
            result['ambiguity_max'],
            strand,
            current_position,
            use_deconvolution=use_deconvolution,
            proportional_split=True,
        )

        if result['ambiguity_range'] > 0:
            result['correction_applied'].append('netseq_refinement')

        # NEW-065: When poly-A walkback was applied, NET-seq must not assign
        # the corrected 3' end to a poly-A base (plus strand: 'A';
        # minus strand reference: 'T').  The walkback anchor (current_position)
        # is the first non-A/T base and is always a safe position.
        #
        # Filter out any assignment at a poly-A base.  If at least one non-A/T
        # assignment remains, renormalise those fractions and use them.  If all
        # assignments are at poly-A bases, fall back to the walkback anchor
        # (fraction=1.0).
        if polya_walkback_applied and genome and assignments:
            chrom_seq, _chrom_seq_key = get_chrom_sequence(genome, chrom)
            if not chrom_seq:
                chrom_seq = None
        else:
            chrom_seq = None
        if chrom_seq is not None:
            polya_base = 'A' if strand == '+' else 'T'
            filtered_assignments = [
                a for a in assignments
                if a['assigned_position'] < len(chrom_seq)
                and chrom_seq[a['assigned_position']].upper() != polya_base
            ]
            if not filtered_assignments:
                # All peaks were at poly-A bases — fall back to walkback anchor.
                assignments = [{
                    'assigned_position': current_position,
                    'fraction': 1.0,
                    'confidence': 'low',
                    'method': 'polya_walkback_fallback',
                    'peak_signal': 0.0,
                }]
            else:
                # Re-normalise fractions after removing poly-A positions.
                total_frac = sum(a['fraction'] for a in filtered_assignments)
                for a in filtered_assignments:
                    a['fraction'] = a['fraction'] / total_frac
                assignments = filtered_assignments

        # Build one output row per proportional assignment.
        # For a single dominant peak the list has length 1 (fraction=1.0).
        # Only the first row is the primary result — subsequent rows are
        # split-read duplicates and must not be double-counted in per-read stats.
        output_rows = []
        for i, assignment in enumerate(assignments):
            row = dict(result)  # shallow copy — all scalar fields
            row['correction_applied'] = list(result['correction_applied'])
            row['qc_flags'] = list(result['qc_flags'])
            row['corrected_3prime'] = assignment['assigned_position']
            row['fraction'] = assignment['fraction']
            row['netseq_confidence'] = assignment['confidence']
            row['netseq_method'] = assignment['method']
            row['netseq_peak_signal'] = assignment['peak_signal']
            row['is_primary_result'] = (i == 0)
            # Propagate NET-seq confidence into overall confidence
            if assignment['confidence'] == 'low':
                row['confidence'] = 'low'
            elif assignment['confidence'] == 'medium' and row['confidence'] == 'high':
                row['confidence'] = 'medium'
            if not row['qc_flags']:
                row['qc_flags'].append('PASS')
            output_rows.append(row)

        return output_rows

    # No NET-seq refinement — single row, full weight
    result['fraction'] = 1.0
    result['is_primary_result'] = True
    if not result['qc_flags']:
        result['qc_flags'].append('PASS')
    return [result]


def process_bam_file(
    bam_path: str,
    genome_path: str,
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    output_path: Optional[str] = None,
    max_reads: Optional[int] = None,
    apply_3ss_rescue: bool = True,
    dt_primed_cDNA: bool = False,
    ont_cDNA: bool = False,
    use_dorado_polya: bool = False,
) -> List[Dict]:
    """
    Process BAM file and apply all corrections.

    Args:
        bam_path: Path to BAM file
        genome_path: Path to genome FASTA
        apply_atract: Apply A-tract ambiguity detection
        apply_ag_mispriming: Apply AG mispriming screening
        apply_polya_trim: Apply poly(A) tail trimming
        apply_indel_correction: Apply indel artifact correction
        netseq_dir: Optional directory with NET-seq BigWig files
        output_path: Optional output TSV path
        max_reads: Optional maximum number of reads to process

    Returns:
        List of correction result dicts
    """
    # Load genome
    print(f"Loading genome from {genome_path}...")
    genome = load_genome(genome_path)

    # Load NET-seq data if provided
    netseq_loader = None
    if netseq_dir:
        print(f"Loading NET-seq data from {netseq_dir}...")
        netseq_loader = _load_netseq(netseq_dir)
        n_bw = len(netseq_loader.bigwigs)
        src = f"{n_bw} BigWig file(s)" if n_bw else "bundled TSV"
        print(f"  Loaded {src}")

    # Open BAM file
    print(f"Processing BAM file: {bam_path}")
    results = []
    n_processed = 0

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Skip secondary/supplementary alignments
            if read.is_secondary or read.is_supplementary:
                continue

            # Apply corrections
            read_results = correct_read_3prime(
                read,
                genome,
                apply_atract=apply_atract,
                apply_ag_mispriming=apply_ag_mispriming,
                apply_polya_trim=apply_polya_trim,
                apply_indel_correction=apply_indel_correction,
                netseq_loader=netseq_loader,
                apply_3ss_rescue=apply_3ss_rescue,
                dt_primed_cDNA=dt_primed_cDNA,
                ont_cDNA=ont_cDNA,
                use_dorado_polya=use_dorado_polya,
            )

            results.extend(read_results)
            n_processed += 1

            # Progress reporting
            if n_processed % 10000 == 0:
                print(f"  Processed {n_processed:,} reads...")

            # Limit for testing
            if max_reads and n_processed >= max_reads:
                print(f"  Reached max_reads limit ({max_reads})")
                break

    print(f"Completed processing {n_processed:,} reads")

    # Write output if requested
    if output_path:
        write_output_tsv(results, output_path)

    return results
