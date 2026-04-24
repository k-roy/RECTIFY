#!/usr/bin/env python3
"""
BAM processing pipeline for RECTIFY.

This module orchestrates the per-read correction workflow (``correct_read_3prime``).
Modules are labelled by their internal ID; execution order is as listed.

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

Supports parallel processing via:
- Per-chromosome parallelization for small genomes (e.g., yeast)
- Region-based parallelization with gap splitting for large chromosomes

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, List, Optional, Tuple, Union
from multiprocessing import Pool
from functools import partial
import gzip
import logging
import pysam
from pathlib import Path

from . import atract_detector
from . import ag_mispriming
from . import polya_trimmer
from . import indel_corrector
from . import netseq_refiner
from .indel_corrector import VariantAwareHomopolymerRescue
from .processing_stats import ProcessingStats, write_stats_tsv, generate_stats_report
from ..utils.genome import load_genome, standardize_chrom_name, reverse_complement
from .polya_model import PolyAModel, load_model as load_polya_model
from ..utils.alignment import (
    extract_junctions_simple,
    extract_soft_clips,
    format_junctions_string,
)
from .splice_aware_5prime import rescue_3ss_truncation as _rescue_3ss
from . import false_junction_filter as _fjf
from ..slurm import get_available_cpus

# BAM writing and CIGAR surgery — imported here so callers that import from
# bam_processor continue to work unchanged.
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
    write_corrected_3ends_bedgraph,
)

# Position index — imported here for the same backwards-compat reason.
from .position_index import write_position_index  # noqa: F401  (re-exported)

logger = logging.getLogger(__name__)


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
    gene_interval_trees: Optional[Dict] = None,
    polya_model: Optional[PolyAModel] = None,
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
    # and must not be re-corrected (XK=1 flag set by consensus.py).
    try:
        if read.get_tag('XK') == 1:
            original_position, strand = get_read_3prime_position(read)
            if original_position is None:
                return []
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
    chrom = read.reference_name

    # Standardize chromosome name
    chrom_std = standardize_chrom_name(chrom)

    # Get 5' end position (transcription start site)
    five_prime_position = get_read_5prime_position(read, strand)
    if five_prime_position is None:
        logger.warning(f"Could not compute 5' position for read {read.query_name}, skipping")
        return []

    # Read XR BAM tag (set by consensus.py when 5' soft-clip was rescued during consensus)
    try:
        five_prime_rescued = read.get_tag('XR') == 1
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
    if annotated_junctions or _real_junctions:
        _ss_junctions: set = set()
        if annotated_junctions:
            _ss_junctions.update(annotated_junctions)
        for _js, _je in _real_junctions:
            _ss_junctions.add((chrom_std, _js, _je))
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
            from .analyze.gene_attribution import compute_read_gene_attribution
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
            result['ag_confidence'] = ag_result['confidence']

            # Reduce confidence for misprimed reads
            if result['confidence'] == 'high':
                result['confidence'] = 'medium' if ag_result['confidence'] == 'medium' else 'low'

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
    polya_walkback_applied = False
    if genome and not _has_3prime_hardclip and not softclip_rescue_applied:
        wb = indel_corrector.find_polya_boundary(read, strand, genome)
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
    max_reads: Optional[int] = None
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
                netseq_loader=netseq_loader
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


def write_output_tsv(results: List[Dict], output_path: str):
    """
    Write correction results to TSV file.

    Args:
        results: List of correction result dicts
        output_path: Path to output TSV file
    """
    print(f"Writing output to {output_path}...")

    with open(output_path, 'w') as f:
        # Write header
        header = [
            'read_id', 'chrom', 'strand',
            'original_3prime', 'corrected_3prime',
            'five_prime_position',  # TSS end of the read (v2.6.0)
            'five_prime_rescued',   # 1 if 5' end was corrected by junction rescue (v2.7.9)
            'five_prime_exon_cigar',  # SAM CIGAR for exon segment of Cat3 rescue (v2.8.0)
            'alignment_start', 'alignment_end',  # Full read body interval (v2.6.0)
            'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
            'polya_length',  # Total observed poly(A) tail length
            'aligned_a_length', 'soft_clip_a_length',  # Breakdown of poly(A)
            'junctions', 'n_junctions',  # Splice junctions (v2.7.0)
            'five_prime_soft_clip_length', 'three_prime_soft_clip_length',  # Soft clips (v2.7.0)
            'mapq',  # Mapping quality (v2.7.0)
            'correction_applied', 'confidence', 'qc_flags', 'fraction',
            'gene_id',  # Per-read gene attribution (optional)
            'pt_tag',      # dorado pt:i signal-level poly(A) length (blank if absent) (v2.9.0)
            'polya_score', # poly(A) model confidence 0-1 (blank if no model) (v2.9.0)
            'polya_source',  # 'pt_tag' | 'model' | 'none' (v2.9.0)
            'sc_homopolymer_extension',  # Cat2: under-called homopolymer bases (v2.9.1)
            'sc_rescued_seq',            # Cat2: non-poly-A bases matched to ref (v2.9.1)
            'sc_original_softclip_len',  # Cat2: original 3' soft-clip length (v2.9.1)
            'five_prime_intron_clip_pos',  # Case 4: exon-side intron boundary for BAM hard-clip (-1 if N/A)
        ]
        f.write('\t'.join(header) + '\n')

        # Write results
        for result in results:
            _pt = result.get('pt_tag')
            _ps = result.get('polya_score')
            row = [
                result['read_id'],
                result['chrom'],
                result['strand'],
                str(result['original_3prime']),
                str(result['corrected_3prime']),
                str(result.get('five_prime_position', '')),  # 5' end (TSS)
                '1' if result.get('five_prime_rescued') else '0',  # 5' rescue flag
                result.get('five_prime_exon_cigar') or '',  # exon CIGAR for Cat3
                str(result.get('alignment_start', '')),  # Read body start
                str(result.get('alignment_end', '')),  # Read body end (exclusive)
                str(result['ambiguity_min']),
                str(result['ambiguity_max']),
                str(result['ambiguity_range']),
                str(result.get('polya_length', 0)),  # poly(A) length, default 0 if not computed
                str(result.get('aligned_a_length', 0)),  # Aligned A's
                str(result.get('soft_clip_a_length', 0)),  # Soft-clipped A's
                result.get('junctions_str', ''),  # Junctions as semicolon-separated string
                str(result.get('n_junctions', 0)),  # Number of junctions
                str(result.get('five_prime_soft_clip_length', 0)),  # 5' soft clip
                str(result.get('three_prime_soft_clip_length', 0)),  # 3' soft clip
                str(result.get('mapq', 0)),  # Mapping quality
                ','.join(result['correction_applied']) if result['correction_applied'] else 'none',
                result['confidence'],
                ','.join(result['qc_flags']),
                f"{result.get('fraction', 1.0):.6f}",
                result.get('gene_id') or '',  # Per-read gene attribution (empty if not computed)
                str(_pt) if _pt is not None else '',
                f'{_ps:.4f}' if _ps is not None else '',
                result.get('polya_source', 'none'),
                str(result.get('sc_homopolymer_extension', 0)),  # Cat2 CIGAR surgery
                result.get('sc_rescued_seq', ''),
                str(result.get('sc_original_softclip_len', 0)),
                str(result.get('five_prime_intron_clip_pos', -1)),  # Case 4 BAM clip
            ]
            f.write('\t'.join(row) + '\n')

    print(f"  Wrote {len(results):,} corrected positions")


def generate_summary_report(results: List[Dict]) -> str:
    """
    Generate summary report from correction results.

    Uses single-pass accumulation for efficiency with large datasets.

    Args:
        results: List of correction result dicts

    Returns:
        Formatted report string
    """
    import numpy as np

    n_total = len(results)
    if n_total == 0:
        return "No reads processed."

    # Single-pass accumulation of all metrics
    n_with_ambiguity = 0
    n_corrected = 0
    ambiguity_ranges = []
    corrected_shifts = []
    by_confidence = {'high': 0, 'medium': 0, 'low': 0}

    for r in results:
        # Ambiguity check
        ambig_range = r['ambiguity_range']
        ambiguity_ranges.append(ambig_range)
        if ambig_range > 0:
            n_with_ambiguity += 1

        # Correction check
        shift = abs(r['corrected_3prime'] - r['original_3prime'])
        if shift > 0:
            n_corrected += 1
            corrected_shifts.append(shift)

        # Confidence
        conf = r['confidence']
        if conf in by_confidence:
            by_confidence[conf] += 1

    # Build report
    report = []
    report.append("=" * 60)
    report.append("RECTIFY Correction Summary")
    report.append("=" * 60)
    report.append("")
    report.append("Overall:")
    report.append(f"  Total reads:            {n_total:,}")
    report.append(f"  With ambiguity:         {n_with_ambiguity:,} ({100*n_with_ambiguity/n_total:.1f}%)")
    report.append(f"  Position corrected:     {n_corrected:,} ({100*n_corrected/n_total:.1f}%)")
    report.append("")

    report.append("Ambiguity Ranges:")
    report.append(f"  Mean:                   {np.mean(ambiguity_ranges):.2f} bp")
    report.append(f"  Median:                 {np.median(ambiguity_ranges):.1f} bp")
    report.append(f"  Maximum:                {max(ambiguity_ranges)} bp")
    report.append("")

    if n_corrected > 0:
        report.append("Position Corrections:")
        report.append(f"  Mean shift:             {np.mean(corrected_shifts):.2f} bp")
        report.append(f"  Median shift:           {np.median(corrected_shifts):.1f} bp")
        report.append(f"  Maximum shift:          {max(corrected_shifts)} bp")
        report.append("")

    report.append("Confidence Levels:")
    for conf in ['high', 'medium', 'low']:
        count = by_confidence.get(conf, 0)
        pct = 100.0 * count / n_total
        report.append(f"  {conf:10s}          {count:7,} ({pct:5.1f}%)")

    report.append("")
    report.append("=" * 60)

    return "\n".join(report)


# =============================================================================
# Parallel Processing Functions
# =============================================================================


def find_coverage_gaps(
    bam_path: str,
    chrom: str,
    min_gap_size: int = 10000
) -> List[Tuple[int, int]]:
    """
    Find gaps in read coverage that can be used to split chromosome.

    Large chromosomes can be split at coverage deserts where no reads exist,
    allowing independent parallel processing of regions.

    Args:
        bam_path: Path to BAM file
        chrom: Chromosome name
        min_gap_size: Minimum gap size in bp to consider a split point

    Returns:
        List of (start, end) tuples representing gap regions
    """
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        # Get chromosome length
        chrom_length = None
        for idx, ref in enumerate(bam.references):
            if ref == chrom:
                chrom_length = bam.lengths[idx]
                break

        if chrom_length is None:
            return []

        # Collect read intervals and merge to find covered regions
        # Uses O(reads) memory instead of O(genome_positions)
        intervals = []
        for read in bam.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            intervals.append((read.reference_start, read.reference_end))

    # Use index to efficiently find gaps
    # Strategy: sample positions and find where no reads map
    gaps = []
    last_covered = 0

    if not intervals:
        return [(0, chrom_length)]

    # Sort and merge overlapping intervals
    intervals.sort()
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])

    # Find gaps between merged intervals
    prev_end = 0
    for seg_start, seg_end in merged:
        if seg_start - prev_end >= min_gap_size:
            gaps.append((prev_end, seg_start))
        prev_end = seg_end

    # Check final gap
    if chrom_length - prev_end >= min_gap_size:
        gaps.append((prev_end, chrom_length))

    return gaps


def get_processing_regions(
    bam_path: str,
    min_gap_size: int = 10000,
    max_region_size: int = 100000
) -> List[Tuple[str, int, int]]:
    """
    Get list of regions for parallel processing.

    Splits large chromosomes at coverage gaps to create balanced work units.
    Each region is capped at max_region_size reference bases to limit per-worker
    result list size and prevent OOM when many workers run in parallel.
    The 100k default limits peak worker memory to ~150 MB on a 1552x-coverage
    nanopore DRS yeast dataset (vs ~1 GB at 500k). Workers with very high
    local coverage (e.g. RPL genes at 50x) stay under ~300 MB.

    Args:
        bam_path: Path to BAM file
        min_gap_size: Minimum gap to use as split point
        max_region_size: Target maximum region size

    Returns:
        List of (chrom, start, end) tuples representing processing regions
    """
    bam = pysam.AlignmentFile(bam_path, 'rb')
    chromosomes = list(bam.references)
    chrom_lengths = dict(zip(bam.references, bam.lengths))
    bam.close()

    regions = []

    for chrom in chromosomes:
        chrom_len = chrom_lengths[chrom]

        # For small chromosomes, process as single region
        if chrom_len <= max_region_size:
            regions.append((chrom, 0, chrom_len))
            continue

        # For large chromosomes, split at gaps
        gaps = find_coverage_gaps(bam_path, chrom, min_gap_size)

        if not gaps:
            # No gaps found — split uniformly by coordinate so workers don't
            # accumulate millions of reads in one list.
            n_sub = max(1, (chrom_len + max_region_size - 1) // max_region_size)
            sub_size = chrom_len // n_sub
            for i in range(n_sub):
                sub_start = i * sub_size
                sub_end = chrom_len if i == n_sub - 1 else (i + 1) * sub_size
                regions.append((chrom, sub_start, sub_end))
            continue

        # Create regions between gaps
        current_start = 0
        for gap_start, gap_end in gaps:
            if gap_start > current_start:
                regions.append((chrom, current_start, gap_start))
            current_start = gap_end

        # Final region after last gap
        if current_start < chrom_len:
            regions.append((chrom, current_start, chrom_len))

    return regions


def _process_region_worker(
    region: Tuple[str, int, int],
    bam_path: str,
    genome: Dict[str, str],
    apply_atract: bool,
    apply_ag_mispriming: bool,
    ag_threshold: float = 0.65,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    variant_aware_rescue: Optional[VariantAwareHomopolymerRescue] = None,
    annotated_junctions: Optional[set] = None,
    gene_interval_trees: Optional[Dict] = None,
    polya_model: Optional[PolyAModel] = None,
) -> List[Dict]:
    """
    Worker function to process a single region.

    Opens its own BAM handle to avoid multiprocessing issues.

    Args:
        region: Tuple of (chrom, start, end)
        bam_path: Path to BAM file
        genome: Pre-loaded genome dict (shared via fork)
        apply_*: Correction module flags
        netseq_dir: Optional NET-seq BigWig directory
        variant_aware_rescue: Optional variant-aware rescue object (from first pass)

    Returns:
        List of correction result dicts for reads in region
    """
    chrom, start, end = region
    results = []

    # Open BAM for this region
    bam = pysam.AlignmentFile(bam_path, 'rb')

    # Load NET-seq if needed (per-worker to avoid file handle issues)
    netseq_loader = None
    if netseq_dir:
        netseq_loader = _load_netseq(netseq_dir)

    try:
        for read in bam.fetch(chrom, start, end):
            # Skip unmapped/secondary/supplementary
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            # Boundary deduplication: when a chromosome is split by coordinate
            # (not by a coverage gap), pysam.fetch returns reads that *overlap*
            # [start, end), so a read spanning a sub-region boundary would be
            # processed by both workers.  Only process reads that START in this
            # region — reads starting before 'start' belong to the previous worker.
            if read.reference_start < start:
                continue

            # Apply corrections
            read_results = correct_read_3prime(
                read,
                genome,
                apply_atract=apply_atract,
                apply_ag_mispriming=apply_ag_mispriming,
                ag_threshold=ag_threshold,
                apply_polya_trim=apply_polya_trim,
                apply_indel_correction=apply_indel_correction,
                netseq_loader=netseq_loader,
                variant_aware_rescue=variant_aware_rescue,
                annotated_junctions=annotated_junctions,
                gene_interval_trees=gene_interval_trees,
                polya_model=polya_model,
            )
            results.extend(read_results)

    finally:
        bam.close()
        if netseq_loader:
            netseq_loader.close()

    return results


def process_bam_file_parallel(
    bam_path: str,
    genome_path: str,
    n_threads: int = 0,
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    ag_threshold: float = 0.65,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    output_path: Optional[str] = None,
    max_reads: Optional[int] = None,
    min_gap_size: int = 10000,
    show_progress: bool = True,
    return_stats: bool = False,
    variant_aware: bool = False,
    variant_output_path: Optional[str] = None,
    annotated_junctions: Optional[set] = None,
    gene_interval_trees: Optional[Dict] = None,
    polya_model_path: Optional[str] = None,
    variant_scan_cache: Optional[str] = None,
) -> Union[List[Dict], Tuple[List[Dict], ProcessingStats]]:
    """
    Process BAM file with parallel region-based processing.

    Splits chromosomes into independent regions at coverage gaps,
    then processes regions in parallel.

    Args:
        bam_path: Path to BAM file
        genome_path: Path to genome FASTA
        n_threads: Number of threads (0 = auto-detect from SLURM/system)
        apply_atract: Apply A-tract ambiguity detection
        apply_ag_mispriming: Apply AG mispriming screening
        apply_polya_trim: Apply poly(A) tail trimming
        apply_indel_correction: Apply indel artifact correction
        netseq_dir: Optional directory with NET-seq BigWig files
        output_path: Optional output TSV path
        max_reads: Optional maximum number of reads (for testing)
        min_gap_size: Minimum gap size for region splitting
        show_progress: Show progress information
        return_stats: Return ProcessingStats alongside results
        variant_aware: Enable variant-aware homopolymer rescue (two-pass)
        variant_output_path: Optional path to write potential variants TSV

    Returns:
        List of correction result dicts, or tuple of (results, stats) if return_stats=True
    """
    # Auto-detect threads if not specified
    if n_threads <= 0:
        n_threads = get_available_cpus()

    logger.info(f"Using {n_threads} thread(s) for parallel processing")

    # Load genome (shared across workers via fork)
    logger.info(f"Loading genome from {genome_path}...")
    genome = load_genome(genome_path)

    # Load poly(A) model if provided (once, shared across all workers via fork)
    polya_model = load_polya_model(Path(polya_model_path) if polya_model_path else None)
    if polya_model is not None:
        logger.info(f"Loaded poly(A) model from {polya_model_path}")

    # Run variant-aware scan if enabled (first pass).
    # When variant_scan_cache points to a pre-built pickle (from `rectify prescan`),
    # load it directly and skip the scan — used for chunked correction pipelines.
    import pickle as _pickle
    variant_aware_rescue = None
    if variant_aware:
        if variant_scan_cache and Path(variant_scan_cache).exists():
            logger.info(f"Loading pre-computed variant scan from cache: {variant_scan_cache}")
            with open(variant_scan_cache, 'rb') as _f:
                variant_aware_rescue = _pickle.load(_f)
        else:
            variant_aware_rescue = run_variant_aware_scan(
                bam_path=bam_path,
                genome=genome,
                min_variant_fraction=0.8,
                min_reads_for_variant_call=5,
                output_variants_path=variant_output_path,
            )

    # Count total reads and filtered reads for stats
    stats = ProcessingStats()
    logger.info("Counting reads in BAM file...")
    bam = pysam.AlignmentFile(bam_path, 'rb')
    for read in bam:
        stats.total_reads_in_bam += 1
        if read.is_unmapped:
            stats.reads_unmapped += 1
        elif read.is_secondary:
            stats.reads_secondary += 1
        elif read.is_supplementary:
            stats.reads_supplementary += 1
    bam.close()
    logger.info(f"  Total reads: {stats.total_reads_in_bam:,}")
    logger.info(f"  Unmapped: {stats.reads_unmapped:,}")
    logger.info(f"  Secondary: {stats.reads_secondary:,}")
    logger.info(f"  Supplementary: {stats.reads_supplementary:,}")

    # Get processing regions
    logger.info("Identifying processing regions...")
    regions = get_processing_regions(bam_path, min_gap_size=min_gap_size)
    logger.info(f"  Found {len(regions)} processing regions")

    # Single-threaded fallback
    if n_threads == 1:
        logger.info("Processing regions sequentially...")
        all_results = []
        for i, region in enumerate(regions):
            chrom, start, end = region
            if show_progress:
                logger.info(f"  Region {i+1}/{len(regions)}: {chrom}:{start:,}-{end:,}")

            results = _process_region_worker(
                region, bam_path, genome,
                apply_atract, apply_ag_mispriming,
                ag_threshold,
                apply_polya_trim, apply_indel_correction,
                netseq_dir, variant_aware_rescue,
                annotated_junctions,
                gene_interval_trees,
                polya_model,
            )
            all_results.extend(results)

            if max_reads and len(all_results) >= max_reads:
                logger.info(f"  Reached max_reads limit ({max_reads})")
                break

        logger.info(f"Completed processing {len(all_results):,} reads")

        # Update stats from results
        for result in all_results:
            stats.update_from_result(result)

        if output_path:
            write_output_tsv(all_results, output_path)
            write_position_index(all_results, output_path)

        if return_stats:
            return all_results, stats
        return all_results

    # Multi-threaded processing
    logger.info(f"Processing {len(regions)} regions across {n_threads} workers...")

    # Create partial function with fixed arguments
    worker_func = partial(
        _process_region_worker,
        bam_path=bam_path,
        genome=genome,
        apply_atract=apply_atract,
        apply_ag_mispriming=apply_ag_mispriming,
        ag_threshold=ag_threshold,
        apply_polya_trim=apply_polya_trim,
        apply_indel_correction=apply_indel_correction,
        netseq_dir=netseq_dir,
        variant_aware_rescue=variant_aware_rescue,
        annotated_junctions=annotated_junctions,
        gene_interval_trees=gene_interval_trees,
        polya_model=polya_model,
    )

    all_results = []

    # Process with pool
    with Pool(n_threads) as pool:
        try:
            # Try to use tqdm for progress if available
            from tqdm import tqdm
            results_iter = tqdm(
                pool.imap(worker_func, regions),
                total=len(regions),
                desc="Processing regions"
            )
        except ImportError:
            results_iter = pool.imap(worker_func, regions)

        for region_results in results_iter:
            all_results.extend(region_results)

            if max_reads and len(all_results) >= max_reads:
                logger.info(f"Reached max_reads limit ({max_reads})")
                break

    logger.info(f"Completed processing {len(all_results):,} reads")

    # Update stats from results
    for result in all_results:
        stats.update_from_result(result)

    # Write output if requested
    if output_path:
        write_output_tsv(all_results, output_path)
        write_position_index(all_results, output_path)

    if return_stats:
        return all_results, stats
    return all_results


# =============================================================================
# Streaming Output Mode (for memory efficiency on large BAMs)
# =============================================================================


def process_bam_streaming(
    bam_path: str,
    genome_path: str,
    output_path: str,
    chunk_size: int = 10000,
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    ag_threshold: float = 0.65,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    show_progress: bool = True,
    annotated_junctions: Optional[set] = None,
    gene_interval_trees: Optional[Dict] = None,
    polya_model_path: Optional[str] = None,
) -> ProcessingStats:
    """
    Process BAM file with streaming output to minimize memory usage.

    Writes results directly to output file instead of accumulating in memory.
    Recommended for BAM files > 10GB.

    Args:
        bam_path: Input BAM path
        genome_path: Genome FASTA path
        output_path: Output TSV path (supports .gz compression)
        chunk_size: Number of reads to process before writing
        apply_*: Correction module flags
        netseq_dir: Optional NET-seq BigWig directory
        show_progress: Show progress information

    Returns:
        ProcessingStats object with comprehensive statistics
    """
    # Load genome
    logger.info(f"Loading genome from {genome_path}...")
    genome = load_genome(genome_path)

    # Load NET-seq if provided
    netseq_loader = None
    if netseq_dir:
        logger.info(f"Loading NET-seq data from {netseq_dir}...")
        netseq_loader = _load_netseq(netseq_dir)

    # Load poly(A) model if provided
    polya_model = load_polya_model(Path(polya_model_path) if polya_model_path else None)
    if polya_model is not None:
        logger.info(f"Loaded poly(A) model from {polya_model_path}")

    # Initialize comprehensive stats
    stats = ProcessingStats()

    # Position count accumulator for the compact index
    from collections import defaultdict as _defaultdict
    _pos_counts = _defaultdict(float)

    # Open output file (support gzip)
    if output_path.endswith('.gz'):
        out_fh = gzip.open(output_path, 'wt')
    else:
        out_fh = open(output_path, 'w')

    _failed = False
    try:
        # Write header
        header = [
            'read_id', 'chrom', 'strand',
            'original_3prime', 'corrected_3prime',
            'five_prime_position',  # TSS end of the read
            'five_prime_rescued',   # 1 if 5' end was corrected by junction rescue (v2.7.9)
            'five_prime_exon_cigar',  # SAM CIGAR for exon segment of Cat3 rescue (v2.8.0)
            'alignment_start', 'alignment_end',  # Full read body interval
            'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
            'polya_length',  # Total observed poly(A) tail length
            'aligned_a_length', 'soft_clip_a_length',  # Breakdown of poly(A)
            'junctions', 'n_junctions',  # Splice junctions
            'five_prime_soft_clip_length', 'three_prime_soft_clip_length',  # Soft clips
            'mapq',  # Mapping quality
            'correction_applied', 'confidence', 'qc_flags', 'fraction',
            'gene_id',  # Per-read gene attribution (optional)
            'pt_tag',      # dorado pt:i signal-level poly(A) length (blank if absent) (v2.9.0)
            'polya_score', # poly(A) model confidence 0-1 (blank if no model) (v2.9.0)
            'polya_source',  # 'pt_tag' | 'model' | 'none' (v2.9.0)
            'sc_homopolymer_extension',  # Cat2: under-called homopolymer bases (v2.9.1)
            'sc_rescued_seq',            # Cat2: non-poly-A bases matched to ref (v2.9.1)
            'sc_original_softclip_len',  # Cat2: original 3' soft-clip length (v2.9.1)
        ]
        out_fh.write('\t'.join(header) + '\n')

        # Open BAM
        bam = pysam.AlignmentFile(bam_path, 'rb')
        chunk = []
        _progress_interval = 100000  # Log every N reads
        _next_progress = _progress_interval
        import time as _time
        _t_start = _time.monotonic()
        _t_last = _t_start

        try:
            for read in bam:
                stats.total_reads_in_bam += 1

                # Track filtered reads
                if read.is_unmapped:
                    stats.reads_unmapped += 1
                    continue
                if read.is_secondary:
                    stats.reads_secondary += 1
                    continue
                if read.is_supplementary:
                    stats.reads_supplementary += 1
                    continue

                read_results = correct_read_3prime(
                    read, genome,
                    apply_atract=apply_atract,
                    apply_ag_mispriming=apply_ag_mispriming,
                    ag_threshold=ag_threshold,
                    apply_polya_trim=apply_polya_trim,
                    apply_indel_correction=apply_indel_correction,
                    netseq_loader=netseq_loader,
                    annotated_junctions=annotated_junctions,
                    gene_interval_trees=gene_interval_trees,
                    polya_model=polya_model,
                )
                chunk.extend(read_results)

                # Update comprehensive stats and position index accumulator
                for result in read_results:
                    stats.update_from_result(result)
                    _pos_counts[(result['chrom'], result['corrected_3prime'], result['strand'])] += float(result.get('fraction', 1.0))

                # Write chunk when full
                if len(chunk) >= chunk_size:
                    _write_results_chunk(out_fh, chunk)
                    chunk = []

                # Progress log: fires every _progress_interval BAM reads (not rows)
                if show_progress and stats.total_reads_in_bam >= _next_progress:
                    _t_now = _time.monotonic()
                    _elapsed = _t_now - _t_start
                    _rate_per_min = (stats.total_reads_in_bam / _elapsed * 60) if _elapsed > 0 else 0
                    logger.info(
                        f"  Processed {stats.total_reads_in_bam:,} reads  "
                        f"({_rate_per_min / 1000:.0f}k reads/min, {_elapsed / 60:.1f} min elapsed)"
                    )
                    _next_progress += _progress_interval

            # Write remaining
            if chunk:
                _write_results_chunk(out_fh, chunk)

        except Exception:
            _failed = True
            bam.close()
            raise
        else:
            bam.close()

    finally:
        out_fh.close()
        if netseq_loader:
            netseq_loader.close()
        # Remove partial output on failure so callers don't see a truncated file
        if _failed:
            _partial = Path(output_path)
            if _partial.exists():
                _partial.unlink()
                logger.warning(f"Removed partial output file after error: {output_path}")

    logger.info(f"Completed processing {stats.reads_processed:,} reads")
    logger.info(f"  Output written to {output_path}")

    # Write compact position index
    _base = str(output_path)
    if _base.endswith('.tsv.gz'):
        _index_path = _base[:-len('.tsv.gz')] + '_index.bed.gz'
    elif _base.endswith('.tsv'):
        _index_path = _base[:-len('.tsv')] + '_index.bed.gz'
    else:
        _index_path = _base + '_index.bed.gz'

    with gzip.open(_index_path, 'wt') as _idx_fh:
        _idx_fh.write('chrom\tcorrected_3prime\tstrand\tcount\n')
        for (chrom, pos, strand), count in sorted(_pos_counts.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
            _idx_fh.write(f'{chrom}\t{pos}\t{strand}\t{count:.6f}\n')
    logger.info(f"Position index written to {_index_path}")

    return stats


def _rebuild_pos_counts_from_partial(output_path: str, pos_counts: dict) -> int:
    """Re-read an existing partial output TSV to rebuild pos_counts. Returns read count."""
    n = 0
    try:
        with open(output_path, 'r') as f:
            header = f.readline().strip().split('\t')
            try:
                ci = header.index('chrom')
                pi = header.index('corrected_3prime')
                si = header.index('strand')
                fi = header.index('fraction')
            except ValueError:
                return 0
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(ci, pi, si, fi):
                    continue
                try:
                    key = (parts[ci], int(parts[pi]), parts[si])
                    frac = float(parts[fi]) if parts[fi].strip() else 1.0
                    pos_counts[key] = pos_counts.get(key, 0.0) + frac
                    n += 1
                except (ValueError, IndexError):
                    pass
    except OSError:
        pass
    return n


def process_bam_streaming_parallel(
    bam_path: str,
    genome_path: str,
    output_path: str,
    n_threads: int = 0,
    apply_atract: bool = True,
    apply_ag_mispriming: bool = False,
    ag_threshold: float = 0.65,
    apply_polya_trim: bool = False,
    apply_indel_correction: bool = False,
    netseq_dir: Optional[str] = None,
    show_progress: bool = True,
    variant_aware: bool = False,
    variant_output_path: Optional[str] = None,
    annotated_junctions: Optional[set] = None,
    gene_interval_trees: Optional[Dict] = None,
    polya_model_path: Optional[str] = None,
    min_gap_size: int = 10000,
    checkpoint_dir: Optional[str] = None,
    variant_scan_cache: Optional[str] = None,
) -> ProcessingStats:
    """
    Process BAM file with parallel region workers and streaming output.

    Combines the memory efficiency of process_bam_streaming (results written
    immediately, not accumulated in RAM) with the throughput of
    process_bam_file_parallel (multiple worker processes).

    Regions are computed from genomic coverage gaps (same as
    process_bam_file_parallel). Each worker processes one region and returns
    its results; the main thread writes each batch to disk as it arrives.

    Args:
        bam_path: Input BAM path
        genome_path: Genome FASTA path
        output_path: Output TSV path (supports .gz compression)
        n_threads: Worker processes (0 = auto-detect from SLURM/system)
        apply_*: Correction module flags
        netseq_dir: Optional NET-seq BigWig directory
        show_progress: Log reads/min progress
        variant_aware: Enable variant-aware homopolymer rescue (two-pass)
        variant_output_path: Optional path to write potential variants TSV

    Returns:
        ProcessingStats object
    """
    import time as _time

    if n_threads <= 0:
        n_threads = get_available_cpus()

    logger.info(f"Using {n_threads} worker(s) for parallel streaming processing")

    # Load genome once — shared across workers via fork
    logger.info(f"Loading genome from {genome_path}...")
    genome = load_genome(genome_path)

    # Load poly(A) model once (shared via fork)
    polya_model = load_polya_model(Path(polya_model_path) if polya_model_path else None)
    if polya_model is not None:
        logger.info(f"Loaded poly(A) model from {polya_model_path}")

    # Checkpoint directory setup
    import pickle as _pickle
    _chk_dir = Path(checkpoint_dir) if checkpoint_dir else None
    if _chk_dir:
        _chk_dir.mkdir(parents=True, exist_ok=True)
    _rescue_pkl = (_chk_dir / 'rescue_scan.pkl') if _chk_dir else None

    # Variant-aware first pass (single-threaded pre-scan).
    # Priority: (1) variant_scan_cache from `rectify prescan` (chunked pipelines),
    #           (2) checkpoint pkl from a prior streaming run (resume),
    #           (3) run the scan fresh.
    variant_aware_rescue = None
    if variant_aware:
        if variant_scan_cache and Path(variant_scan_cache).exists():
            logger.info(f"Loading pre-computed variant scan from cache: {variant_scan_cache}")
            with open(variant_scan_cache, 'rb') as _f:
                variant_aware_rescue = _pickle.load(_f)
        elif _rescue_pkl and _rescue_pkl.exists():
            logger.info(f"Loading pre-computed variant scan from checkpoint: {_rescue_pkl}")
            with open(_rescue_pkl, 'rb') as _f:
                variant_aware_rescue = _pickle.load(_f)
        else:
            variant_aware_rescue = run_variant_aware_scan(
                bam_path=bam_path,
                genome=genome,
                min_variant_fraction=0.8,
                min_reads_for_variant_call=5,
                output_variants_path=variant_output_path,
            )
            if _rescue_pkl:
                with open(_rescue_pkl, 'wb') as _f:
                    _pickle.dump(variant_aware_rescue, _f, protocol=_pickle.HIGHEST_PROTOCOL)
                logger.info(f"Scan checkpoint saved: {_rescue_pkl}")

    # Count total reads and filtered reads for stats (pre-scan, single-threaded)
    stats = ProcessingStats()
    logger.info("Counting reads in BAM file...")
    _bam_prescan = pysam.AlignmentFile(bam_path, 'rb')
    for _r in _bam_prescan:
        stats.total_reads_in_bam += 1
        if _r.is_unmapped:
            stats.reads_unmapped += 1
        elif _r.is_secondary:
            stats.reads_secondary += 1
        elif _r.is_supplementary:
            stats.reads_supplementary += 1
    _bam_prescan.close()
    logger.info(f"  Total reads: {stats.total_reads_in_bam:,}")

    # Compute genomic regions
    logger.info("Identifying processing regions...")
    regions = get_processing_regions(bam_path, min_gap_size=min_gap_size)
    logger.info(f"  {len(regions)} regions across {n_threads} workers")
    _pos_counts: dict = {}
    _t_start = _time.monotonic()
    _next_progress = 100000

    # Checkpoint: find already-completed regions
    _done_region_idxs: set = set()
    if _chk_dir:
        for _sf in _chk_dir.glob('region_*.done'):
            try:
                _done_region_idxs.add(int(_sf.stem.split('_')[1]))
            except (ValueError, IndexError):
                pass
        if _done_region_idxs:
            logger.info(
                f"Checkpoint: {len(_done_region_idxs)}/{len(regions)} regions already done, resuming..."
            )
            _n_partial = _rebuild_pos_counts_from_partial(output_path, _pos_counts)
            logger.info(f"  Rebuilt pos_counts from {_n_partial:,} rows in partial output")

    _regions_to_run = [r for i, r in enumerate(regions) if i not in _done_region_idxs]
    _orig_idxs = [i for i, _r in enumerate(regions) if i not in _done_region_idxs]

    worker_func = partial(
        _process_region_worker,
        bam_path=bam_path,
        genome=genome,
        apply_atract=apply_atract,
        apply_ag_mispriming=apply_ag_mispriming,
        ag_threshold=ag_threshold,
        apply_polya_trim=apply_polya_trim,
        apply_indel_correction=apply_indel_correction,
        netseq_dir=netseq_dir,
        variant_aware_rescue=variant_aware_rescue,
        annotated_junctions=annotated_junctions,
        gene_interval_trees=gene_interval_trees,
        polya_model=polya_model,
    )

    # Build the TSV header (identical to process_bam_streaming)
    _header = [
        'read_id', 'chrom', 'strand',
        'original_3prime', 'corrected_3prime',
        'five_prime_position', 'five_prime_rescued', 'five_prime_exon_cigar',
        'alignment_start', 'alignment_end',
        'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
        'polya_length', 'aligned_a_length', 'soft_clip_a_length',
        'junctions', 'n_junctions',
        'five_prime_soft_clip_length', 'three_prime_soft_clip_length',
        'mapq', 'correction_applied', 'confidence', 'qc_flags', 'fraction',
        'gene_id', 'pt_tag', 'polya_score', 'polya_source',
        'sc_homopolymer_extension', 'sc_rescued_seq', 'sc_original_softclip_len',
        'five_prime_intron_clip_pos',
    ]

    _resuming = bool(_done_region_idxs) and Path(output_path).exists()
    if output_path.endswith('.gz'):
        out_fh = gzip.open(output_path, 'at' if _resuming else 'wt')
    else:
        out_fh = open(output_path, 'a' if _resuming else 'w')

    _failed = False
    try:
        if not _resuming:
            out_fh.write('\t'.join(_header) + '\n')

        _map_fn = 'imap' if _chk_dir else 'imap_unordered'
        with Pool(n_threads) as pool:
            _iter = getattr(pool, _map_fn)(worker_func, _regions_to_run)
            for _batch_num, region_results in enumerate(_iter):
                _write_results_chunk(out_fh, region_results)

                for result in region_results:
                    stats.update_from_result(result)
                    _pos_counts.setdefault(
                        (result['chrom'], result['corrected_3prime'], result['strand']), 0.0
                    )
                    _pos_counts[(result['chrom'], result['corrected_3prime'], result['strand'])] += float(result.get('fraction', 1.0))

                if _chk_dir:
                    orig_idx = _orig_idxs[_batch_num]
                    (_chk_dir / f'region_{orig_idx:04d}.done').touch()

                if show_progress and stats.reads_processed >= _next_progress:
                    _elapsed = _time.monotonic() - _t_start
                    _rate = (stats.reads_processed / _elapsed * 60) if _elapsed > 0 else 0
                    logger.info(
                        f"  Processed {stats.reads_processed:,} reads  "
                        f"({_rate / 1000:.0f}k reads/min, {_elapsed / 60:.1f} min elapsed)"
                    )
                    _next_progress += 100000

    except Exception:
        _failed = True
        raise
    finally:
        out_fh.close()
        if _failed and not _chk_dir:
            _partial = Path(output_path)
            if _partial.exists():
                _partial.unlink()
                logger.warning(f"Removed partial output after error: {output_path}")
        elif _failed and _chk_dir:
            logger.warning(
                f"Processing failed — partial output preserved for checkpoint resume: {output_path}"
            )

    logger.info(f"Completed processing {stats.reads_processed:,} reads")
    logger.info(f"  Output written to {output_path}")

    # Write compact position index
    _base = str(output_path)
    if _base.endswith('.tsv.gz'):
        _index_path = _base[:-len('.tsv.gz')] + '_index.bed.gz'
    elif _base.endswith('.tsv'):
        _index_path = _base[:-len('.tsv')] + '_index.bed.gz'
    else:
        _index_path = _base + '_index.bed.gz'

    with gzip.open(_index_path, 'wt') as _idx_fh:
        _idx_fh.write('chrom\tcorrected_3prime\tstrand\tcount\n')
        for (chrom, pos, strand), count in sorted(_pos_counts.items(), key=lambda x: (x[0][0], x[0][1])):
            _idx_fh.write(f'{chrom}\t{pos}\t{strand}\t{count:.6f}\n')
    logger.info(f"Position index written to {_index_path}")

    return stats


def _write_results_chunk(fh, results: List[Dict]):
    """Write a chunk of results to file handle."""
    for result in results:
        _pt = result.get('pt_tag')
        _ps = result.get('polya_score')
        row = [
            result['read_id'],
            result['chrom'],
            result['strand'],
            str(result['original_3prime']),
            str(result['corrected_3prime']),
            str(result.get('five_prime_position', '')),  # 5' end (TSS)
            '1' if result.get('five_prime_rescued') else '0',  # 5' rescue flag
            result.get('five_prime_exon_cigar') or '',  # exon CIGAR for Cat3
            str(result.get('alignment_start', '')),  # Read body start
            str(result.get('alignment_end', '')),  # Read body end (exclusive)
            str(result['ambiguity_min']),
            str(result['ambiguity_max']),
            str(result['ambiguity_range']),
            str(result.get('polya_length', 0)),  # poly(A) length, default 0 if not computed
            str(result.get('aligned_a_length', 0)),  # Aligned A's
            str(result.get('soft_clip_a_length', 0)),  # Soft-clipped A's
            result.get('junctions_str', ''),  # Junctions as semicolon-separated string
            str(result.get('n_junctions', 0)),  # Number of junctions
            str(result.get('five_prime_soft_clip_length', 0)),  # 5' soft clip
            str(result.get('three_prime_soft_clip_length', 0)),  # 3' soft clip
            str(result.get('mapq', 0)),  # Mapping quality
            ','.join(result['correction_applied']) if result['correction_applied'] else 'none',
            result['confidence'],
            ','.join(result['qc_flags']),
            f"{result.get('fraction', 1.0):.6f}",
            result.get('gene_id') or '',  # Per-read gene attribution (empty if not computed)
            str(_pt) if _pt is not None else '',
            f'{_ps:.4f}' if _ps is not None else '',
            result.get('polya_source', 'none'),
            str(result.get('sc_homopolymer_extension', 0)),  # Cat2 CIGAR surgery
            result.get('sc_rescued_seq', ''),
            str(result.get('sc_original_softclip_len', 0)),
            str(result.get('five_prime_intron_clip_pos', -1)),  # Case 4 BAM clip
        ]
        fh.write('\t'.join(row) + '\n')


def run_variant_aware_scan(
    bam_path: str,
    genome: Dict[str, str],
    min_variant_fraction: float = 0.8,
    min_reads_for_variant_call: int = 5,
    output_variants_path: Optional[str] = None,
) -> VariantAwareHomopolymerRescue:
    """
    Run first pass to scan reads and build variant frequency map.

    This enables variant-aware homopolymer rescue by identifying positions
    where high mismatch frequency indicates a true variant (not basecalling error).

    Args:
        bam_path: Path to BAM file
        genome: Pre-loaded genome dict
        min_variant_fraction: Threshold for variant call (default 0.8 = 80%)
        min_reads_for_variant_call: Minimum reads at position (default 5)
        output_variants_path: Optional path to write potential variants TSV

    Returns:
        VariantAwareHomopolymerRescue object ready for second pass
    """
    logger.info("Running variant-aware scan (first pass)...")

    rescue = VariantAwareHomopolymerRescue(
        min_variant_fraction=min_variant_fraction,
        min_reads_for_variant_call=min_reads_for_variant_call,
        min_homopolymer_len=4,
        max_rescue_bases=3,
    )

    n_scanned = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            strand = '-' if read.is_reverse else '+'
            rescue.scan_read(read, strand, genome, end='3prime')
            n_scanned += 1

            if n_scanned % 500000 == 0:
                logger.info(f"  Scanned {n_scanned:,} reads...")

    logger.info(f"  Total scanned: {n_scanned:,} reads")

    # Finalize scan
    rescue.finalize_scan()

    # Get statistics
    stats = rescue.get_statistics()
    logger.info(f"  Positions with mismatches: {stats['total_positions_scanned']:,}")
    logger.info(f"  With sufficient coverage: {stats['positions_with_sufficient_coverage']:,}")
    logger.info(f"  Potential variants: {stats['potential_variants_detected']:,}")

    # Write potential variants if requested
    if output_variants_path:
        variants = rescue.get_potential_variants()
        if variants:
            logger.info(f"  Writing {len(variants)} potential variants to {output_variants_path}")
            with open(output_variants_path, 'w') as f:
                header = ['chrom', 'position', 'ref_base', 'homopolymer_base',
                          'total_reads', 'mismatch_fraction',
                          'dominant_mismatch_base', 'dominant_mismatch_fraction']
                f.write('\t'.join(header) + '\n')
                for v in variants:
                    row = [
                        v['chrom'],
                        str(v['position']),
                        v['ref_base'],
                        v['homopolymer_base'],
                        str(v['total_reads']),
                        f"{v['mismatch_fraction']:.3f}",
                        v['dominant_mismatch_base'],
                        f"{v['dominant_mismatch_fraction']:.3f}",
                    ]
                    f.write('\t'.join(row) + '\n')

    return rescue


def generate_summary_from_stats(stats) -> str:
    """
    Generate summary report from ProcessingStats or legacy dict.

    Args:
        stats: ProcessingStats object or legacy Dict from process_bam_streaming

    Returns:
        Formatted report string
    """
    # Handle both ProcessingStats and legacy dict format
    if isinstance(stats, ProcessingStats):
        return generate_stats_report(stats)

    # Legacy dict format (backwards compatibility)
    n_total = stats.get('total_reads', 0)
    if n_total == 0:
        return "No reads processed."

    report = []
    report.append("=" * 60)
    report.append("RECTIFY Correction Summary (Streaming Mode)")
    report.append("=" * 60)
    report.append("")
    report.append("Overall:")
    report.append(f"  Total reads:            {n_total:,}")
    report.append(f"  With ambiguity:         {stats['with_ambiguity']:,} ({100*stats['with_ambiguity']/n_total:.1f}%)")
    report.append(f"  Position corrected:     {stats['position_corrected']:,} ({100*stats['position_corrected']/n_total:.1f}%)")
    report.append("")
    report.append("Confidence Levels:")
    for conf in ['high', 'medium', 'low']:
        count = stats['by_confidence'].get(conf, 0)
        pct = 100.0 * count / n_total
        report.append(f"  {conf:10s}          {count:7,} ({pct:5.1f}%)")
    report.append("")
    report.append("=" * 60)

    return "\n".join(report)
