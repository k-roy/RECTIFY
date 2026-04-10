#!/usr/bin/env python3
"""
BAM processing pipeline for RECTIFY.

This module orchestrates the full correction workflow:
1. A-tract ambiguity detection (universal)
2. AG mispriming screening (for oligo-dT)
3. Poly(A) tail trimming (when poly(A) sequenced)
4. Indel artifact correction (when poly(A) sequenced)
5. NET-seq refinement (optional)

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
from ..utils.genome import load_genome, standardize_chrom_name
from ..utils.alignment import (
    extract_junctions_simple,
    extract_soft_clips,
    format_junctions_string,
)
from .splice_aware_5prime import rescue_3ss_truncation as _rescue_3ss
from . import false_junction_filter as _fjf
from ..slurm import get_available_cpus

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
    # NEW-063: Skip poly-A walkback for reads with hard-clipped 3' ends.
    # The 3' sequence is absent; walkback would explore an unanchored ambiguity window.
    _has_3prime_hardclip = False
    if read.cigartuples:
        _cigar_ops = read.cigartuples
        if strand == '+' and _cigar_ops[-1][0] == 5:
            _has_3prime_hardclip = True
        elif strand == '-' and _cigar_ops[0][0] == 5:
            _has_3prime_hardclip = True

    polya_walkback_applied = False
    if genome and not _has_3prime_hardclip:
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

    # NEW-061: If walkback landed inside an artifact N op, snap to the N boundary
    # so the artifact is fully eliminated from the rectified CIGAR.
    if polya_walkback_applied and _artifact_analyses:
        for _art in _artifact_analyses:
            if _art.junction_start <= current_position < _art.junction_end:
                if strand == '+':
                    current_position = _art.junction_start - 1
                    result['ambiguity_min'] = current_position
                else:
                    current_position = _art.junction_end
                    result['ambiguity_max'] = current_position
                result['ambiguity_range'] = abs(current_position - original_position)
                break

    # Update corrected position (after all position-moving corrections)
    result['corrected_3prime'] = current_position

    # Update ambiguity window for indel/rescue corrections that moved the position
    # but did NOT set the window themselves (i.e. polya_walkback was not applied).
    if current_position != original_position and not polya_walkback_applied:
        if strand == '+':
            result['ambiguity_min'] = max(0, current_position - result['ambiguity_range'])
            result['ambiguity_max'] = current_position
        else:
            # Minus-strand corrections move toward lower coords
            result['ambiguity_min'] = max(0, current_position - result['ambiguity_range'])
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

        if result['ambiguity_range'] > 1:
            result['correction_applied'].append('netseq_refinement')

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
            'alignment_start', 'alignment_end',  # Full read body interval (v2.6.0)
            'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
            'polya_length',  # Total observed poly(A) tail length
            'aligned_a_length', 'soft_clip_a_length',  # Breakdown of poly(A)
            'junctions', 'n_junctions',  # Splice junctions (v2.7.0)
            'five_prime_soft_clip_length', 'three_prime_soft_clip_length',  # Soft clips (v2.7.0)
            'mapq',  # Mapping quality (v2.7.0)
            'correction_applied', 'confidence', 'qc_flags', 'fraction',
            'gene_id',  # Per-read gene attribution (optional)
        ]
        f.write('\t'.join(header) + '\n')

        # Write results
        for result in results:
            row = [
                result['read_id'],
                result['chrom'],
                result['strand'],
                str(result['original_3prime']),
                str(result['corrected_3prime']),
                str(result.get('five_prime_position', '')),  # 5' end (TSS)
                '1' if result.get('five_prime_rescued') else '0',  # 5' rescue flag
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
            ]
            f.write('\t'.join(row) + '\n')

    print(f"  Wrote {len(results):,} corrected positions")


# ---------------------------------------------------------------------------
# CIGAR op classification constants (used by clip_read_to_corrected_3prime)
# ---------------------------------------------------------------------------
# Ops that advance the reference coordinate
_REF_CONSUMING: frozenset = frozenset([0, 2, 3, 7, 8])   # M, D, N, =, X
# Ops that consume stored query bases (i.e. present in query_sequence)
_QUERY_CONSUMING: frozenset = frozenset([0, 1, 4, 7, 8])  # M, I, S, =, X
# Note: H (5) is excluded from both — hard-clipped bases are NOT in query_sequence.


def clip_read_to_corrected_3prime(
    read: pysam.AlignedSegment,
    corrected_3prime: int,
    strand: str,
) -> bool:
    """
    Hard-clip *read* in-place so its 3' alignment boundary matches *corrected_3prime*.

    For + strand reads the 3' boundary is ``reference_end - 1`` (rightmost ref base).
    For − strand reads it is ``reference_start`` (leftmost ref base, which is the 3'
    end for reverse-strand RNA).

    Bases removed from the CIGAR are converted to a single trailing / leading H
    (hard-clip) op and deleted from the stored query sequence and quality array.
    Any existing soft-clip at the clipped end is subsumed into the new hard clip.

    If *corrected_3prime* equals the current boundary the read is unchanged and
    False is returned.  Unmapped reads, reads without CIGAR, and reads where the
    corrected position would leave an empty alignment are also left unchanged.

    Args:
        read:             pysam AlignedSegment — modified in-place.
        corrected_3prime: Target 3' reference coordinate (0-based, inclusive).
        strand:           '+' or '-'.

    Returns:
        True if clipping was applied; False otherwise.
    """
    cigar = list(read.cigartuples or [])
    if not cigar or read.is_unmapped:
        return False

    seq   = read.query_sequence
    quals = read.query_qualities  # numpy uint8 array or None

    if strand == '+':
        current_end = read.reference_end - 1  # 0-based inclusive right boundary
        if corrected_3prime >= current_end:
            return False  # no clipping needed

        n_ref_clip = current_end - corrected_3prime  # reference bases to remove

        n_query_remove = 0

        # 1. Strip trailing soft-clips: their bases ARE in query_sequence.
        while cigar and cigar[-1][0] == 4:  # S
            n_query_remove += cigar.pop()[1]

        # 2. Strip trailing hard-clips: their bases are NOT in query_sequence.
        while cigar and cigar[-1][0] == 5:  # H
            cigar.pop()

        # 3. Walk CIGAR ops from the right, removing n_ref_clip reference bases.
        n_ref_removed = 0
        while cigar and n_ref_removed < n_ref_clip:
            op, length = cigar[-1]
            ref_in_op   = length if op in _REF_CONSUMING  else 0
            query_in_op = length if op in _QUERY_CONSUMING else 0
            need        = n_ref_clip - n_ref_removed

            if ref_in_op <= need:
                cigar.pop()
                n_ref_removed  += ref_in_op
                n_query_remove += query_in_op
            else:
                # Partial removal: trim `need` reference positions from this op.
                cigar[-1] = (op, length - need)
                n_ref_removed  += need
                n_query_remove += need if op in _QUERY_CONSUMING else 0
                break

        if not cigar:
            return False  # degenerate: entire alignment clipped

        # 4. Append hard-clip for removed query bases (if any).
        if n_query_remove > 0:
            cigar.append((5, n_query_remove))

        # 5. Apply to read.  Set sequence before qualities (setting seq resets quals).
        if n_query_remove > 0 and seq is not None:
            new_seq  = seq[:-n_query_remove]
            new_qual = quals[:-n_query_remove] if quals is not None else None
            read.query_sequence  = new_seq
            read.query_qualities = new_qual
        read.cigartuples = cigar
        return True

    else:  # minus strand — clip from the left (remove 3' bases at low coordinates)
        current_start = read.reference_start  # 0-based inclusive left boundary
        if corrected_3prime <= current_start:
            return False

        n_ref_clip = corrected_3prime - current_start

        n_query_remove = 0

        # 1. Strip leading soft-clips.
        while cigar and cigar[0][0] == 4:  # S
            n_query_remove += cigar.pop(0)[1]

        # 2. Strip leading hard-clips.
        while cigar and cigar[0][0] == 5:  # H
            cigar.pop(0)

        # 3. Walk CIGAR ops from the left.
        n_ref_removed = 0
        while cigar and n_ref_removed < n_ref_clip:
            op, length = cigar[0]
            ref_in_op   = length if op in _REF_CONSUMING  else 0
            query_in_op = length if op in _QUERY_CONSUMING else 0
            need        = n_ref_clip - n_ref_removed

            if ref_in_op <= need:
                cigar.pop(0)
                n_ref_removed  += ref_in_op
                n_query_remove += query_in_op
            else:
                cigar[0] = (op, length - need)
                n_ref_removed  += need
                n_query_remove += need if op in _QUERY_CONSUMING else 0
                break

        if not cigar:
            return False

        # 4. Prepend hard-clip for removed query bases.
        if n_query_remove > 0:
            cigar.insert(0, (5, n_query_remove))

        # 5. Apply to read; also shift reference_start for minus strand clips.
        if n_query_remove > 0 and seq is not None:
            new_seq  = seq[n_query_remove:]
            new_qual = quals[n_query_remove:] if quals is not None else None
            read.query_sequence  = new_seq
            read.query_qualities = new_qual
        read.cigartuples     = cigar
        read.reference_start = current_start + n_ref_removed
        return True


def write_corrected_bam(
    input_bam_path: str,
    corrected_tsv_path: str,
    output_bam_path: str,
) -> Dict[str, int]:
    """
    Write a new BAM with every read hard-clipped at its corrected 3' end.

    Reads per-read corrections from *corrected_tsv_path* (the TSV produced by
    ``rectify correct``) and applies CIGAR-level hard-clipping via
    :func:`clip_read_to_corrected_3prime`.  The corrected position is the
    ``corrected_3prime`` column; only the first row per read is used (covers
    NET-seq multi-peak Cat6 reads where all rows share the same dominant position).

    Reads absent from the TSV (unmapped, secondary, supplementary) are written
    unchanged.  The output BAM preserves the same header as the input.

    Caller is responsible for sorting and indexing the output BAM if needed
    (e.g. ``pysam.sort()`` / ``pysam.index()``).

    Args:
        input_bam_path:      Path to the original input BAM.
        corrected_tsv_path:  Path to the corrected_3ends.tsv from ``rectify correct``.
        output_bam_path:     Destination BAM path.  Overwritten if it exists.

    Returns:
        Dict with summary counts:
            ``'total'``     — total reads written
            ``'clipped'``   — reads whose CIGAR was modified
            ``'unchanged'`` — reads written without modification
    """
    # Load corrections: read_id -> (corrected_3prime, strand)
    corrections: Dict[str, tuple] = {}
    try:
        with open(corrected_tsv_path) as _f:
            hdr = _f.readline().strip().split('\t')
            try:
                i_id     = hdr.index('read_id')
                i_pos    = hdr.index('corrected_3prime')
                i_strand = hdr.index('strand')
            except ValueError as exc:
                raise ValueError(
                    f"Required column missing in {corrected_tsv_path}: {exc}"
                ) from exc
            for line in _f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(i_id, i_pos, i_strand):
                    continue
                rid = parts[i_id]
                if rid in corrections:
                    continue  # keep first (dominant) row only
                try:
                    corrections[rid] = (int(parts[i_pos]), parts[i_strand])
                except (ValueError, IndexError):
                    pass
    except OSError as exc:
        raise OSError(
            f"Cannot read corrected TSV {corrected_tsv_path}: {exc}"
        ) from exc

    logger.info(
        "write_corrected_bam: loaded %d corrected positions from %s",
        len(corrections), corrected_tsv_path,
    )

    stats: Dict[str, int] = {'total': 0, 'clipped': 0, 'unchanged': 0}

    with pysam.AlignmentFile(input_bam_path, 'rb') as bam_in, \
         pysam.AlignmentFile(output_bam_path, 'wb', header=bam_in.header) as bam_out:

        for read in bam_in:
            stats['total'] += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                bam_out.write(read)
                stats['unchanged'] += 1
                continue

            correction = corrections.get(read.query_name)
            if correction is None:
                bam_out.write(read)
                stats['unchanged'] += 1
                continue

            corrected_pos, strand = correction
            clipped = clip_read_to_corrected_3prime(read, corrected_pos, strand)

            bam_out.write(read)
            if clipped:
                stats['clipped'] += 1
            else:
                stats['unchanged'] += 1

    return stats


def softclip_read_to_corrected_3prime(
    read: pysam.AlignedSegment,
    corrected_3prime: int,
    strand: str,
) -> bool:
    """
    Soft-clip *read* in-place so its 3' alignment boundary matches *corrected_3prime*.

    Identical logic to :func:`clip_read_to_corrected_3prime` except that removed
    query bases are replaced with a ``S`` (soft-clip) op rather than an ``H``
    (hard-clip) op.  The bases remain in ``query_sequence`` and are visible in IGV
    when "Show soft-clipped bases" is enabled, making the poly(A) tail location
    apparent without affecting pileup or coverage tracks.

    ``reference_start`` is updated for minus-strand reads exactly as in the hard-clip
    path (soft-clips do not consume reference, so the first aligned base moves).

    Returns True if soft-clipping was applied; False otherwise.
    """
    cigar = list(read.cigartuples or [])
    if not cigar or read.is_unmapped:
        return False

    seq   = read.query_sequence
    quals = read.query_qualities

    if strand == '+':
        current_end = read.reference_end - 1
        if corrected_3prime >= current_end:
            return False

        n_ref_clip    = current_end - corrected_3prime
        n_query_remove = 0

        # 1. Absorb trailing soft-clips into the new soft-clip.
        while cigar and cigar[-1][0] == 4:   # S
            n_query_remove += cigar.pop()[1]
        # 2. Strip trailing hard-clips (not in sequence; discard silently).
        while cigar and cigar[-1][0] == 5:   # H
            cigar.pop()

        # 3. Walk CIGAR from right, accounting for n_ref_clip reference bases.
        n_ref_removed = 0
        while cigar and n_ref_removed < n_ref_clip:
            op, length = cigar[-1]
            ref_in_op   = length if op in _REF_CONSUMING  else 0
            query_in_op = length if op in _QUERY_CONSUMING else 0
            need        = n_ref_clip - n_ref_removed
            if ref_in_op <= need:
                cigar.pop()
                n_ref_removed  += ref_in_op
                n_query_remove += query_in_op
            else:
                cigar[-1] = (op, length - need)
                n_ref_removed  += need
                n_query_remove += need if op in _QUERY_CONSUMING else 0
                break

        if not cigar:
            return False

        # 4. Append soft-clip (sequence stays in read).
        if n_query_remove > 0:
            cigar.append((4, n_query_remove))  # S, not H
        read.cigartuples = cigar
        return True

    else:  # minus strand — clip from the left
        current_start = read.reference_start
        if corrected_3prime <= current_start:
            return False

        n_ref_clip    = corrected_3prime - current_start
        n_query_remove = 0

        # 1. Absorb leading soft-clips.
        while cigar and cigar[0][0] == 4:   # S
            n_query_remove += cigar.pop(0)[1]
        # 2. Strip leading hard-clips.
        while cigar and cigar[0][0] == 5:   # H
            cigar.pop(0)

        # 3. Walk CIGAR from left.
        n_ref_removed = 0
        while cigar and n_ref_removed < n_ref_clip:
            op, length = cigar[0]
            ref_in_op   = length if op in _REF_CONSUMING  else 0
            query_in_op = length if op in _QUERY_CONSUMING else 0
            need        = n_ref_clip - n_ref_removed
            if ref_in_op <= need:
                cigar.pop(0)
                n_ref_removed  += ref_in_op
                n_query_remove += query_in_op
            else:
                cigar[0] = (op, length - need)
                n_ref_removed  += need
                n_query_remove += need if op in _QUERY_CONSUMING else 0
                break

        if not cigar:
            return False

        # 4. Prepend soft-clip (sequence stays in read).
        if n_query_remove > 0:
            cigar.insert(0, (4, n_query_remove))  # S, not H
        read.cigartuples     = cigar
        read.reference_start = current_start + n_ref_removed
        return True


def write_softclipped_bam(
    input_bam_path: str,
    corrected_tsv_path: str,
    output_bam_path: str,
) -> Dict[str, int]:
    """
    Write a new BAM with every read soft-clipped at its corrected 3' end.

    Identical to :func:`write_corrected_bam` except clipping uses
    :func:`softclip_read_to_corrected_3prime` (``S`` ops, bases retained in
    query sequence) rather than hard-clip (``H`` ops, bases removed).  The
    poly(A) tail bases remain visible in IGV when soft-clip display is enabled.
    """
    corrections: Dict[str, tuple] = {}
    try:
        with open(corrected_tsv_path) as _f:
            hdr = _f.readline().strip().split('\t')
            try:
                i_id     = hdr.index('read_id')
                i_pos    = hdr.index('corrected_3prime')
                i_strand = hdr.index('strand')
            except ValueError as exc:
                raise ValueError(
                    f"Required column missing in {corrected_tsv_path}: {exc}"
                ) from exc
            for line in _f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(i_id, i_pos, i_strand):
                    continue
                rid = parts[i_id]
                if rid in corrections:
                    continue
                try:
                    corrections[rid] = (int(parts[i_pos]), parts[i_strand])
                except (ValueError, IndexError):
                    pass
    except OSError as exc:
        raise OSError(
            f"Cannot read corrected TSV {corrected_tsv_path}: {exc}"
        ) from exc

    logger.info(
        "write_softclipped_bam: loaded %d corrected positions from %s",
        len(corrections), corrected_tsv_path,
    )

    stats: Dict[str, int] = {'total': 0, 'clipped': 0, 'unchanged': 0}

    with pysam.AlignmentFile(input_bam_path, 'rb') as bam_in, \
         pysam.AlignmentFile(output_bam_path, 'wb', header=bam_in.header) as bam_out:

        for read in bam_in:
            stats['total'] += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                bam_out.write(read)
                stats['unchanged'] += 1
                continue

            correction = corrections.get(read.query_name)
            if correction is None:
                bam_out.write(read)
                stats['unchanged'] += 1
                continue

            corrected_pos, strand = correction
            clipped = softclip_read_to_corrected_3prime(read, corrected_pos, strand)

            bam_out.write(read)
            if clipped:
                stats['clipped'] += 1
            else:
                stats['unchanged'] += 1

    return stats


def write_polya_trimmed_bam(
    input_bam_path: str,
    output_bam_path: str,
    threshold: float = 0.8,
) -> Dict[str, int]:
    """
    Write a new BAM with 3' poly(A) soft-clips removed from each read.

    Iterates all reads (including secondary and supplementary) from
    *input_bam_path*, strips the RNA-3' poly(A) soft-clip from primary reads
    that pass the A-richness threshold, and writes every read to
    *output_bam_path*.  Header and all BAM tags are preserved unchanged.

    Secondary and supplementary reads are written as-is without trimming
    because their soft-clips may have different semantics and their 3' end
    is not independently defined.

    Args:
        input_bam_path:  Path to input BAM (sorted, indexed or not).
        output_bam_path: Destination BAM path.  Will be overwritten if it
                         exists.  Caller is responsible for sorting/indexing
                         the output if needed.
        threshold:       Minimum A (plus) or T (minus) fraction required to
                         consider a 3' soft-clip a poly(A) tail (default 0.8).

    Returns:
        Dict with summary counts:
            'total'        — total reads written
            'trimmed'      — reads whose poly(A) tail was removed
            'bases_trimmed'— total bases removed across all reads
    """
    from .polya_trimmer import trim_polya_from_bam_read

    stats = {'total': 0, 'trimmed': 0, 'bases_trimmed': 0}

    with pysam.AlignmentFile(input_bam_path, 'rb') as bam_in, \
         pysam.AlignmentFile(output_bam_path, 'wb', header=bam_in.header) as bam_out:

        for read in bam_in:
            stats['total'] += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                bam_out.write(read)
                continue

            strand = '-' if read.is_reverse else '+'
            read, n_trimmed = trim_polya_from_bam_read(read, strand, threshold=threshold)

            if n_trimmed > 0:
                stats['trimmed']       += 1
                stats['bases_trimmed'] += n_trimmed

            bam_out.write(read)

    return stats


def write_position_index(results_or_path, output_tsv_path: str):
    """
    Write a compact position index (chrom, corrected_3prime, strand -> count).

    Args:
        results_or_path: List of result dicts, OR a path to an existing TSV file
        output_tsv_path: Path to the main corrected TSV; index path is derived from it
    """
    from collections import defaultdict
    import gzip as _gzip

    # Derive index path: replace .tsv or .tsv.gz with _index.bed.gz
    _base = str(output_tsv_path)
    if _base.endswith('.tsv.gz'):
        _index_path = _base[:-len('.tsv.gz')] + '_index.bed.gz'
    elif _base.endswith('.tsv'):
        _index_path = _base[:-len('.tsv')] + '_index.bed.gz'
    else:
        _index_path = _base + '_index.bed.gz'

    pos_counts = defaultdict(float)

    if isinstance(results_or_path, (str, Path)):
        # Load from existing TSV
        import pandas as pd
        # Try with fraction; fall back if column missing
        try:
            _df2 = pd.read_csv(str(results_or_path), sep='\t',
                               usecols=['chrom', 'corrected_3prime', 'strand', 'fraction'])
            for _, row in _df2.iterrows():
                pos_counts[(row['chrom'], row['corrected_3prime'], row['strand'])] += float(row.get('fraction', 1.0))
        except (ValueError, KeyError):
            _df2 = pd.read_csv(str(results_or_path), sep='\t',
                               usecols=['chrom', 'corrected_3prime', 'strand'])
            for _, row in _df2.iterrows():
                pos_counts[(row['chrom'], row['corrected_3prime'], row['strand'])] += 1.0
    else:
        # List of result dicts
        for r in results_or_path:
            key = (r['chrom'], r['corrected_3prime'], r['strand'])
            pos_counts[key] += float(r.get('fraction', 1.0))

    # Write gzip-compressed TSV sorted by (chrom, corrected_3prime, strand)
    with _gzip.open(_index_path, 'wt') as f:
        f.write('chrom\tcorrected_3prime\tstrand\tcount\n')
        for (chrom, pos, strand), count in sorted(pos_counts.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
            f.write(f'{chrom}\t{pos}\t{strand}\t{count:.6f}\n')

    logger.info(f"Position index written to {_index_path}")


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
    max_region_size: int = 500000
) -> List[Tuple[str, int, int]]:
    """
    Get list of regions for parallel processing.

    Splits large chromosomes at coverage gaps to create balanced work units.

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
            # No gaps found, process as single region
            regions.append((chrom, 0, chrom_len))
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

    # Run variant-aware scan if enabled (first pass)
    variant_aware_rescue = None
    if variant_aware:
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
            'alignment_start', 'alignment_end',  # Full read body interval
            'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
            'polya_length',  # Total observed poly(A) tail length
            'aligned_a_length', 'soft_clip_a_length',  # Breakdown of poly(A)
            'junctions', 'n_junctions',  # Splice junctions
            'five_prime_soft_clip_length', 'three_prime_soft_clip_length',  # Soft clips
            'mapq',  # Mapping quality
            'correction_applied', 'confidence', 'qc_flags', 'fraction',
            'gene_id',  # Per-read gene attribution (optional)
        ]
        out_fh.write('\t'.join(header) + '\n')

        # Open BAM
        bam = pysam.AlignmentFile(bam_path, 'rb')
        chunk = []

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

                    if show_progress and stats.reads_processed % 100000 == 0:
                        logger.info(f"  Processed {stats.reads_processed:,} reads...")

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


def _write_results_chunk(fh, results: List[Dict]):
    """Write a chunk of results to file handle."""
    for result in results:
        row = [
            result['read_id'],
            result['chrom'],
            result['strand'],
            str(result['original_3prime']),
            str(result['corrected_3prime']),
            str(result.get('five_prime_position', '')),  # 5' end (TSS)
            '1' if result.get('five_prime_rescued') else '0',  # 5' rescue flag
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
