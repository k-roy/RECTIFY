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
import bisect as _bisect
import logging

import pysam

from ..polya import atract_detector
from ..polya import ag_mispriming
from ..polya import polya_trimmer
from ..correct import indel_corrector
from ..netseq import netseq_refiner
from ..correct.indel_corrector import VariantAwareHomopolymerRescue
from ..correct.protocols.quantseq_rev import walkback_quantseq_rev
from ..correct.walkback import (
    APPLIED_WALKBACK as _APPLIED_WALKBACK_READGENOME,
    walkback_drs_full,
)
from ...utils.genome import load_genome, standardize_chrom_name, reverse_complement
from ..polya.polya_model import PolyAModel
from ...utils.alignment import (
    extract_junctions_simple,
    extract_soft_clips,
    format_junctions_string,
)
from ..splice.splice_aware_5prime import rescue_3ss_truncation as _rescue_3ss
from ..splice import false_junction_filter as _fjf

# Backwards-compat re-exports for callers that still reach in via this module.
from .processing_stats import ProcessingStats, write_stats_tsv, generate_stats_report  # noqa: F401
from .bam_writer import (  # noqa: F401  (re-exported)
    clip_read_to_corrected_3prime,
    softclip_read_to_corrected_3prime,
    extend_read_5prime_for_junction_rescue,
    fix_homopolymer_mismatches,
    realign_exon_blocks,
    _load_corrections_from_tsv,
    write_dual_bam,
    write_corrected_bam,
    write_softclipped_bam,
    write_polya_trimmed_bam,
    write_netseq_assigned_bedgraph,
    write_corrected_reads_bedgraph,
)
from ..position_index import write_position_index  # noqa: F401  (re-exported)

# ``process_bam_file`` (defined below) writes its result list via the shared
# TSV emitter in ``output``.  All other peer-module dispatch entry points
# (region planning, parallel/streaming workers, variant scan) are imported
# directly from their home modules by their callers.
from .output import write_output_tsv

logger = logging.getLogger(__name__)

# Search radius (bp) used when pre-filtering pool junctions near a read's 5' end
# before passing them to rescue_3ss_truncation.  Must be ≥ the junction_proximity_bp
# default (5000 bp) used inside rescue_3ss_truncation, plus the maximum intron length
# (yeast: ~1 kb) to also capture minus-strand junctions sorted by intron_start.
_POOL_SEARCH_RADIUS = 10000


def _build_pool_chrom_index(
    pool_junctions: Optional[set],
) -> Optional[Dict[str, List[Tuple[int, int]]]]:
    """Build a per-chromosome sorted index for fast bisect-based pool lookup.

    Returns a dict mapping each standardised chromosome name to a list of
    (intron_start, intron_end) tuples sorted by intron_start.  Returns None
    if pool_junctions is None or empty.

    Call once before starting BAM processing; pass the result as
    ``pool_chrom_index`` to the process_bam_* and correct_read_3prime functions.
    """
    if not pool_junctions:
        return None
    index: Dict[str, List[Tuple[int, int]]] = {}
    for entry in pool_junctions:
        chrom_raw, intron_start, intron_end = entry
        chrom = standardize_chrom_name(chrom_raw)
        if chrom not in index:
            index[chrom] = []
        index[chrom].append((intron_start, intron_end))
    for chrom in index:
        index[chrom].sort()
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
    gene_interval_trees: Optional[Dict] = None,
    polya_model: Optional[PolyAModel] = None,
    dt_primed_cDNA: bool = False,
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
            five_prime_position = get_read_5prime_position(read, strand)
            if five_prime_position is None:
                return []
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

    # Module 2E (pre-pass): filter poly(A)-artifact junctions before 5' rescue
    # so they are never used as 3'SS rescue candidates.
    _genome_ref = _fjf._GenomeDictReference(genome)
    _real_junctions, _artifact_analyses = _fjf.filter_polya_artifact_junctions(read, _genome_ref, strand)

    # Module 2F: 3'SS truncation rescue (post-consensus).
    # Corrects five_prime_position for reads truncated or soft-clipped at the
    # exon 2 / 3' splice site boundary.
    if annotated_junctions or _real_junctions or pool_chrom_index:
        _ss_junctions: set = set()
        if annotated_junctions:
            _ss_junctions.update(annotated_junctions)
        for _js, _je in _real_junctions:
            _ss_junctions.add((chrom_std, _js, _je))
        # Add pool junctions (novel + annotated from prescan) near this read's 5' end.
        # Pool may contain novel junctions absent from the GFF annotation.
        # Use a bisect-indexed lookup to avoid iterating all 200k+ pool entries per read.
        if pool_chrom_index:
            _pool_entries = pool_chrom_index.get(chrom_std)
            if _pool_entries:
                _lo = _bisect.bisect_left(_pool_entries,
                                          (five_prime_position - _POOL_SEARCH_RADIUS,))
                _hi = _bisect.bisect_right(_pool_entries,
                                           (five_prime_position + _POOL_SEARCH_RADIUS,
                                            10 ** 9))
                for _pjs, _pje in _pool_entries[_lo:_hi]:
                    _ss_junctions.add((chrom_std, _pjs, _pje))
        if _ss_junctions:
            _3ss_result = _rescue_3ss(read, genome, _ss_junctions, strand)
            if _3ss_result['rescued']:
                five_prime_rescued = True
                five_prime_position = _3ss_result['five_prime_corrected']
                _five_prime_exon_cigar = _3ss_result.get('five_prime_exon_cigar', '')
                # Record the exon-side intron boundary for BAM hard-clip when:
                #   (a) the alignment's 5' end sits inside the rescued intron, AND
                #   (b) there is no 5' soft-clip to extend via extend_read_5prime.
                # This covers both Case 4 (intronic_snap, no soft clip) and Case 1/2
                # Diagnostic: flag reads whose 5' alignment end falls inside the intron.
                # five_prime_intron_clip_pos records the exon-side boundary for future
                # BAM rerouting surgery (converting intronic M ops → N + upstream exon M).
                # Currently no BAM surgery is performed for these reads — their TSV
                # five_prime_position is correct and no sequence data is hidden.
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
    }

    # Per-read gene attribution via read body overlap
    gene_id = None
    if gene_interval_trees is not None:
        try:
            from ..analyze.gene_attribution import compute_read_gene_attribution
            gene_ids = compute_read_gene_attribution(read, gene_interval_trees, chrom=chrom_std)
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
            read, strand, atract_result=atract_for_polya
        )

        # Always record poly(A) length (even if 0) for completeness
        result['polya_length'] = polya_result['polya_length']
        result['aligned_a_length'] = polya_result['aligned_a_length']
        result['soft_clip_a_length'] = polya_result['soft_clip_a_length']

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
    if genome and not _has_3prime_hardclip and not softclip_rescue_applied:
        if dt_primed_cDNA:
            _chrom_seq = genome.get(chrom_std) or genome.get(chrom)
            if _chrom_seq:
                _orig_wb, _corr_wb, _applied_wb, _gene_strand_wb = walkback_quantseq_rev(
                    read, _chrom_seq
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
            _chrom_seq = genome.get(chrom_std) or genome.get(chrom)
            if _chrom_seq:
                wb = walkback_drs_full(read, _chrom_seq)
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
                #   junction_start (N's left edge; all leading M bases consumed).
                if strand == '+':
                    current_position = _art.junction_start - 1
                else:
                    current_position = _art.junction_start
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
            elif strand == '-' and result['ambiguity_min'] <= _art.junction_start:
                # Case B (minus): ambiguity window extends into the N — clip min
                # to junction_start so NET-seq stays on the 3' side of the N.
                result['ambiguity_min'] = _art.junction_start
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
        if polya_walkback_applied and genome and assignments and chrom_std in genome:
            polya_base = 'A' if strand == '+' else 'T'
            chrom_seq = genome[chrom_std]
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
    dt_primed_cDNA: bool = False,
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
                dt_primed_cDNA=dt_primed_cDNA,
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


