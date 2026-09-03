"""
NET-seq BAM processing for RECTIFY.

Extracts 3' ends from NET-seq BAM files with:
- Strand-aware soft-clip oligo(A) detection
- CIGAR-based alignment parsing
- Chromosome name standardization

Unlike nanopore processing, NET-seq reads are short (~40-76bp) and represent
the 3' end of nascent RNA. The aligned 3' end of the read IS the 3' end
of the nascent transcript.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Generator
import re

import pysam

from ...config import NCBI_TO_CHROM
from ...utils.alignment import extract_soft_clips, parse_cigar
from ..unified_record import UnifiedReadRecord
from ..exclusion_regions import ExclusionRegionDetector


# Extended NCBI mapping for full contig names with metadata
def extract_ncbi_accession(contig_name: str) -> str:
    """
    Extract NCBI accession from full contig name.

    Handles formats like:
    - 'ref|NC_001133| [org=Saccharomyces cerevisiae] ...'
    - 'ref|NC_001133|'
    - 'NC_001133'

    Returns:
        NCBI accession in 'ref|NC_XXXXXX|' format, or original name if not recognized
    """
    # Match ref|NC_XXXXXX| pattern
    match = re.match(r'(ref\|NC_\d+\|)', contig_name)
    if match:
        return match.group(1)

    # Match bare NC_XXXXXX
    match = re.match(r'(NC_\d+)', contig_name)
    if match:
        return f"ref|{match.group(1)}|"

    return contig_name


def standardize_netseq_chrom(contig_name: str) -> str:
    """
    Standardize NET-seq BAM chromosome name to chrI, chrII, etc. format.

    Args:
        contig_name: Full contig name from BAM header

    Returns:
        Standardized chromosome name (chrI, chrII, etc.)
    """
    # Extract NCBI accession
    ncbi_acc = extract_ncbi_accession(contig_name)

    # Look up in NCBI_TO_CHROM mapping
    if ncbi_acc in NCBI_TO_CHROM:
        return NCBI_TO_CHROM[ncbi_acc]

    # Try without the pipes
    ncbi_bare = ncbi_acc.replace('ref|', '').replace('|', '')
    for key, value in NCBI_TO_CHROM.items():
        if ncbi_bare in key:
            return value

    # Return as-is if unrecognized
    return contig_name


def count_terminal_oligo_a_mismatches(
    read: pysam.AlignedSegment,
    strand: str,
) -> int:
    """
    Count terminal A/T mismatches at the biological 3' end.

    BBMap sometimes aligns oligo(A) bases as mismatches (X in CIGAR)
    rather than soft-clipping them. This shifts the apparent 3' position
    and corrupts PSF/deconvolution if not corrected.

    For minus strand: oligo(A) in RNA = T in genomic orientation (at left)
    For plus strand: oligo(A) in RNA = A in genomic orientation (at right)

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'

    Returns:
        Number of terminal oligo(A) bases aligned as mismatches
    """
    if not read.cigartuples or not read.query_sequence:
        return 0

    query_seq = read.query_sequence.upper()
    cigar = read.cigartuples

    n_oligo_a_mismatches = 0

    if strand == '-':
        # Minus strand: 3' end at left side
        # Oligo(A) in RNA = T in genomic BAM orientation
        pos_in_read = 0

        for op, length in cigar:
            if op == 4:  # Soft clip at start - skip past it
                pos_in_read += length
            elif op == 8:  # X = explicit mismatch
                # Count consecutive T's (oligo A in RNA orientation)
                for i in range(length):
                    if pos_in_read + i < len(query_seq):
                        base = query_seq[pos_in_read + i]
                        if base == 'T':
                            n_oligo_a_mismatches += 1
                        else:
                            # Non-T mismatch, stop counting
                            return n_oligo_a_mismatches
                pos_in_read += length
            elif op == 0:  # M = match/mismatch (need MD tag to distinguish)
                # For M operations, check if there's an MD tag
                # If MD tag starts with a non-digit, it's a mismatch
                md_tag = read.get_tag('MD') if read.has_tag('MD') else None
                if md_tag:
                    md_str = str(md_tag)
                    # Count leading mismatches in MD tag (minus strand: oligo-A appears as T)
                    # MD format: "3A5T2" means 3 matches, A mismatch, 5 matches, T mismatch, 2 matches
                    # "0A0T" means consecutive leading mismatches (0-length match runs between them)
                    i = 0
                    while i < len(md_str):
                        if md_str[i].isdigit():
                            # Accumulate full multi-digit number
                            j = i
                            while j < len(md_str) and md_str[j].isdigit():
                                j += 1
                            run_len = int(md_str[i:j])
                            if run_len > 0:
                                break  # Non-zero matches precede mismatches — stop
                            i = j  # Zero-length run: mismatches are consecutive — continue
                        elif md_str[i].isalpha() and md_str[i] != '^':
                            # Mismatch - check if read base is T
                            if pos_in_read + n_oligo_a_mismatches < len(query_seq):
                                base = query_seq[pos_in_read + n_oligo_a_mismatches]
                                if base == 'T':
                                    n_oligo_a_mismatches += 1
                                    i += 1
                                else:
                                    break
                            else:
                                break
                        else:
                            break
                break  # Stop after first alignment block
            else:
                # Other operation (insertion, deletion), stop
                break

        return n_oligo_a_mismatches

    else:
        # Plus strand: 3' end at right side
        # Oligo(A) in RNA = A in genomic BAM orientation
        pos_in_read = len(query_seq)

        for op, length in reversed(cigar):
            if op == 4:  # Soft clip at end - skip past it
                pos_in_read -= length
            elif op == 8:  # X = explicit mismatch
                # Count consecutive A's from the end
                for i in range(length):
                    idx = pos_in_read - 1 - i
                    if idx >= 0:
                        base = query_seq[idx]
                        if base == 'A':
                            n_oligo_a_mismatches += 1
                        else:
                            return n_oligo_a_mismatches
                pos_in_read -= length
            elif op == 0:  # M = match/mismatch
                # Check MD tag for trailing mismatches
                md_tag = read.get_tag('MD') if read.has_tag('MD') else None
                if md_tag:
                    md_str = str(md_tag)
                    # Count trailing mismatches (plus strand: oligo-A appears as A)
                    # Traverse MD backwards; consecutive terminal mismatches are separated by "0"
                    i = len(md_str) - 1
                    while i >= 0:
                        if md_str[i].isdigit():
                            # Accumulate full multi-digit number (reading rightmost digits first)
                            j = i
                            while j >= 0 and md_str[j].isdigit():
                                j -= 1
                            run_len = int(md_str[j+1:i+1])
                            if run_len > 0:
                                break  # Non-zero matches separate mismatches — stop
                            i = j  # Zero-length run: mismatches are consecutive — continue
                        elif md_str[i].isalpha() and md_str[i] != '^':
                            idx = pos_in_read - 1 - n_oligo_a_mismatches
                            if idx >= 0:
                                base = query_seq[idx]
                                if base == 'A':
                                    n_oligo_a_mismatches += 1
                                    i -= 1
                                else:
                                    break
                            else:
                                break
                        else:
                            break
                break
            else:
                break

        return n_oligo_a_mismatches


#: Which terminus of the SEQUENCED read is the RNA 3' end.
#:
#: ``read5p`` -- Churchman-style NET-seq (GSE25107, GSE159603, GSE61332, Burke C. neoformans): the 3' linker
#: is ligated to the RNA 3' end and the read is revcomp(RNA), so the read's 5' terminus IS the RNA 3' end and
#: the **gene strand is the INVERSE of the BAM strand**. This is the MEASURED convention -- verified three
#: independent ways (strand concordance vs annotated genes 18-19%; the randomer-free 2011 library's read-5'
#: soft clip is 85.6% T = a poly(A) tail through the reverse complement; the read-3' A-richness is the
#: linker's constant A) and it reproduces ``408c_clipdump.py`` + the 697/697 intron self-test.
#: See planning ``478a`` CP1.
#:
#: ``read3p`` -- LEGACY behaviour of this module (read treated as sense: ``strand = '-' if is_reverse``).
#: ⚠️ It is WRONG for every NET-seq dataset in this project: it emits the RNA's 5'-ward terminus on the
#: ANTISENSE strand, and makes ``detect_oligo_a_in_softclip`` read the wrong soft clip. Retained only so
#: pre-existing output can be reproduced/diffed.
NETSEQ_RNA3P_CHOICES = ("read5p", "read3p")
NETSEQ_RNA3P_DEFAULT = "read5p"


def rna_5prime_ward(strand: str) -> int:
    """Reference-coordinate step that moves 5'-ward along the RNA: ``-1`` on ``+``, ``+1`` on ``-``.

    Every 3'-end correction in this module (terminal oligo(A) trim, poly(A) walkback) shortens the
    RNA from its 3' end, i.e. moves the called end 5'-ward. Expressing that as one signed step keeps
    the two strands symmetric by construction -- the inverted minus-strand sign fixed in
    :func:`get_netseq_3prime_position` was exactly a hand-written copy of this going the wrong way.
    """
    return -1 if strand == '+' else 1


def netseq_gene_strand(read: pysam.AlignedSegment, rna3p_at: str = NETSEQ_RNA3P_DEFAULT) -> str:
    """Gene/RNA strand of a NET-seq read under the requested convention.

    Everything downstream in this module is expressed in terms of the returned STRAND (which soft clip holds
    the tail, which alignment terminus is the 3' end, which base is oligo(A)), so this single function is
    what makes the whole module correct or incorrect for a given chemistry.
    """
    if rna3p_at not in NETSEQ_RNA3P_CHOICES:
        raise ValueError(f"unknown rna3p_at {rna3p_at!r}; expected one of {NETSEQ_RNA3P_CHOICES}")
    if rna3p_at == "read5p":
        return '+' if read.is_reverse else '-'
    return '-' if read.is_reverse else '+'


def get_netseq_3prime_position(
    read: pysam.AlignedSegment,
    trim_terminal_oligo_a: bool = True,
    rna3p_at: str = NETSEQ_RNA3P_DEFAULT,
) -> Tuple[int, str, int]:
    """
    Get 3' end position and gene strand from a NET-seq read.

    Once the gene strand is known the position follows without a second convention:
    - '+' gene: RNA 3' end = rightmost aligned base = reference_end - 1
    - '-' gene: RNA 3' end = leftmost  aligned base = reference_start

    ⚠️ The strand itself IS the convention. For Churchman-style NET-seq the read is revcomp(RNA), so the gene
    strand is the INVERSE of the BAM strand (``rna3p_at='read5p'``, the default and the measured truth).
    ``rna3p_at='read3p'`` reproduces this module's legacy sense-read assumption.

    If trim_terminal_oligo_a is True, detects and trims terminal A/T
    mismatches that should have been soft-clipped (oligo(A) bases).

    Args:
        read: pysam AlignedSegment
        trim_terminal_oligo_a: Whether to trim terminal oligo(A) mismatches
        rna3p_at: 'read5p' (NET-seq, default) or 'read3p' (legacy sense-read assumption)

    Returns:
        Tuple of (position, strand, n_trimmed) where:
        - position is 0-based, adjusted for terminal mismatches if requested
        - strand is the GENE strand, '+' or '-'
        - n_trimmed is number of terminal oligo(A) bases trimmed
    """
    if read.is_unmapped or read.reference_end is None:
        raise ValueError(
            f"get_netseq_3prime_position() called on unmapped read "
            f"{read.query_name!r}"
        )

    strand = netseq_gene_strand(read, rna3p_at)

    if strand == '+':
        # Plus strand: 3' end is rightmost aligned base
        position = read.reference_end - 1
    else:
        # Minus strand: 3' end is leftmost aligned base
        position = read.reference_start

    n_trimmed = 0
    if trim_terminal_oligo_a:
        n_trimmed = count_terminal_oligo_a_mismatches(read, strand)
        if n_trimmed > 0:
            # Trimming always moves the 3' end 5'-ward IN RNA, and RNA 5'-ward is the DIRECTION
            # SET BY THE GENE STRAND, not a fixed sign:
            #   '+' gene: the tail extends rightward, so reference_end-1 is too HIGH  -> subtract
            #   '-' gene: the tail extends leftward,  so reference_start  is too LOW  -> ADD
            # The old code subtracted on both strands, which on a '-' gene moved the position a
            # further n_trimmed nt AWAY from the truth (measured 2 x n error, planning 834 Appendix).
            position += rna_5prime_ward(strand) * n_trimmed

    return position, strand, n_trimmed


def detect_oligo_a_in_softclip(
    read: pysam.AlignedSegment,
    strand: str,
    min_a_fraction: float = 0.8,
) -> Dict:
    """
    Detect oligo(A) content in 3' soft-clipped sequence (strand-aware).

    NET-seq reads can have oligo(A) tails in soft-clips due to nascent
    transcript adenylation. This is strand-aware:
    - Plus strand: 3' end is RIGHT side, oligo(A) in right soft-clip
    - Minus strand: 3' end is LEFT side, oligo(A) in left soft-clip
      (soft-clip sequence contains T's which = A's in RNA)

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'
        min_a_fraction: Minimum A (or T for minus) content to classify as oligo(A)

    Returns:
        Dict with:
            - soft_clip_a_length: Number of A's (RNA orientation) in 3' soft-clip
            - three_prime_soft_clip_seq: Full soft-clip sequence (genomic orientation)
            - three_prime_soft_clip_length: Total soft-clip length
            - has_oligo_a: Boolean indicating if clip is A-rich
    """
    clips = extract_soft_clips(read)

    # Find the soft-clip at the biological 3' end
    if strand == '+':
        clip_side = 'right'  # 3' is right for plus strand
    else:
        clip_side = 'left'   # 3' is left for minus strand

    three_prime_clip = next((c for c in clips if c['side'] == clip_side), None)

    if three_prime_clip is None or not three_prime_clip['seq']:
        return {
            'soft_clip_a_length': 0,
            'three_prime_soft_clip_seq': '',
            'three_prime_soft_clip_length': 0,
            'has_oligo_a': False,
        }

    seq = three_prime_clip['seq']
    clip_length = len(seq)

    # Count A's in RNA orientation
    # For minus strand, T in genomic = A in RNA
    if strand == '-':
        a_count = seq.upper().count('T')
    else:
        a_count = seq.upper().count('A')

    a_fraction = a_count / clip_length if clip_length > 0 else 0
    is_a_rich = a_fraction >= min_a_fraction

    return {
        'soft_clip_a_length': a_count if is_a_rich else 0,
        'three_prime_soft_clip_seq': seq,
        'three_prime_soft_clip_length': clip_length,
        'has_oligo_a': is_a_rich,
    }


def get_5prime_softclip_info(
    read: pysam.AlignedSegment,
    strand: str,
) -> Dict:
    """
    Get 5' soft-clip information (for completeness).

    Args:
        read: pysam AlignedSegment
        strand: '+' or '-'

    Returns:
        Dict with 5' soft-clip length and sequence
    """
    clips = extract_soft_clips(read)

    # 5' end is opposite of 3' end
    if strand == '+':
        clip_side = 'left'   # 5' is left for plus strand
    else:
        clip_side = 'right'  # 5' is right for minus strand

    five_prime_clip = next((c for c in clips if c['side'] == clip_side), None)

    if five_prime_clip is None:
        return {
            'five_prime_soft_clip_length': 0,
            'five_prime_soft_clip_seq': '',
        }

    return {
        'five_prime_soft_clip_length': five_prime_clip['length'],
        'five_prime_soft_clip_seq': five_prime_clip['seq'] or '',
    }


def process_netseq_read(
    read: pysam.AlignedSegment,
    chrom_std: str,
    min_a_fraction: float = 0.8,
    trim_terminal_oligo_a: bool = True,
    rna3p_at: str = NETSEQ_RNA3P_DEFAULT,
    junction_pool=None,
    genome=None,
    umi_length: int = 0,
    rescue_max_intronic: int = 10,
    rescue_min_k: int = 1,
    detect_tail: bool = True,
) -> UnifiedReadRecord:
    """
    Process single NET-seq read into a UnifiedReadRecord (a NetseqReadRecord when the rescue or
    tail passes run).

    Creates a record with:
    - three_prime_raw: Raw 3' position from alignment (before any correction)
    - three_prime_corrected: terminal oligo(A) trim -> poly(A) walkback -> donor-side junction rescue
    - soft_clip_a_length: legacy 0.8-A-fraction oligo(A) call over the whole clip (kept for
      back-compat; the authoritative tail call is ``tail_len``)
    - aligned_a_length: Oligo(A) bases aligned as mismatches (trimmed)
    - rescue_* / tail_* : see :mod:`rectify.core.netseq.netseq_rescue`
    - All other standard fields

    Args:
        read: pysam AlignedSegment
        chrom_std: Standardized chromosome name
        min_a_fraction: Minimum A fraction for the legacy whole-clip oligo(A) call
        trim_terminal_oligo_a: Whether to trim terminal oligo(A) mismatches
        rna3p_at: which read terminus is the RNA 3' end -- 'read5p' (NET-seq, default) or 'read3p' (legacy)
        junction_pool: optional :class:`~rectify.core.netseq.netseq_rescue.JunctionPool`
        genome: optional ``{chrom: sequence}`` -- required for the walkback and the rescue
        umi_length: declared randomer length (0 = none)
        rescue_max_intronic: how far past the donor the aligned end may sit (nt)
        rescue_min_k: minimum recovered exon-2 length for a rescue
        detect_tail: run the randomer-aware tail call

    Returns:
        UnifiedReadRecord (NetseqReadRecord) with all fields populated
    """
    from ..netseq.netseq_rescue import (
        NetseqReadRecord, RescueCall, TailCall, _ref_to_query, call_tail, rescue_read, rna_clip,
    )

    # Get 3' position and strand (with optional trimming)
    three_prime_corrected, strand, n_trimmed = get_netseq_3prime_position(
        read, trim_terminal_oligo_a=trim_terminal_oligo_a, rna3p_at=rna3p_at
    )

    # The trim moved the end 5'-ward IN RNA, so the raw end is n_trimmed back the OTHER way.
    three_prime_raw = three_prime_corrected - rna_5prime_ward(strand) * n_trimmed

    # Get 5' position
    if strand == '+':
        five_prime_pos = read.reference_start
    else:
        five_prime_pos = read.reference_end - 1

    # Detect oligo(A) in 3' soft-clip
    oligo_a_info = detect_oligo_a_in_softclip(read, strand, min_a_fraction)

    # Get 5' soft-clip info
    five_prime_info = get_5prime_softclip_info(read, strand)

    # Parse CIGAR
    cigar_str = parse_cigar(read.cigartuples) if read.cigartuples else ''

    # ---- 3'-end correction: walkback, then donor-side junction rescue ---------------------------
    # Both need the chromosome sequence; both are no-ops without it.
    genome_seq = None
    if genome is not None:
        genome_seq = genome.get(chrom_std) or genome.get(read.reference_name)
    clip_rna = rna_clip(read, strand)
    ref_to_query = _ref_to_query(read) if (genome_seq is not None or junction_pool is not None) else None

    tail = TailCall(clip_rna=clip_rna)
    if detect_tail:
        tail = call_tail(
            read, strand, three_prime_corrected, genome_seq,
            umi_length=umi_length, clip_rna=clip_rna, ref_to_query=ref_to_query,
        )
        three_prime_corrected += rna_5prime_ward(strand) * tail.walkback

    rescue = rescue_read(
        read, strand, chrom_std, three_prime_raw, junction_pool, genome_seq,
        umi_length=umi_length, max_intronic=rescue_max_intronic, min_k=rescue_min_k,
        clip_rna=clip_rna, ref_to_query=ref_to_query,
    )
    if rescue.position is not None:
        # A rescued read's clip is exon-2 sequence, not a tail -- the walkback result is discarded.
        three_prime_corrected = rescue.position

    # Total oligo(A) = soft-clipped + aligned-as-mismatch
    total_oligo_a = oligo_a_info['soft_clip_a_length'] + n_trimmed

    # Create the record
    record = NetseqReadRecord(
        read_id=read.query_name,
        chrom=chrom_std,
        strand=strand,

        # 5' end positions
        five_prime_raw=five_prime_pos,
        five_prime_corrected=five_prime_pos,  # No correction for NET-seq
        first_exon_start=None,
        starts_in_intron=False,

        # 3' end positions
        three_prime_raw=three_prime_raw,
        three_prime_corrected=three_prime_corrected,  # trim -> walkback -> junction rescue

        # Alignment span
        alignment_start=read.reference_start,
        alignment_end=read.reference_end,

        # Splice junctions (NET-seq reads are short, rarely spliced)
        junctions=[],
        n_junctions=0,
        junctions_filtered=0,

        # Soft clips
        five_prime_soft_clip_length=five_prime_info['five_prime_soft_clip_length'],
        three_prime_soft_clip_length=oligo_a_info['three_prime_soft_clip_length'],
        five_prime_soft_clip_seq=five_prime_info['five_prime_soft_clip_seq'],
        three_prime_soft_clip_seq=oligo_a_info['three_prime_soft_clip_seq'],

        # Poly(A) / Oligo(A) information
        polya_length=total_oligo_a,  # Total oligo(A) detected
        aligned_a_length=n_trimmed,  # Oligo(A) bases aligned as mismatches (trimmed)
        soft_clip_a_length=oligo_a_info['soft_clip_a_length'],

        # Alignment quality
        mapq=read.mapping_quality,
        cigar_summary=cigar_str,
        alignment_identity=0.0,  # Not computed for NET-seq

        # Consensus info (not applicable for NET-seq)
        best_aligner='bbmap',
        consensus_confidence='',
        n_aligners_support=0,

        # Count
        count=1,

        # NET-seq rescue / tail
        rescue_status=rescue.status,
        rescue_k=rescue.k,
        rescue_r=rescue.r,
        rescue_n_intronic=rescue.n_intronic,
        rescue_decoy_k=rescue.decoy_k,
        rescue_decoy_would_rescue=rescue.decoy_would_rescue,
        rescue_intron_start=rescue.junction.intron_start if rescue.junction else -1,
        rescue_intron_end=rescue.junction.intron_end if rescue.junction else -1,
        tail_len=tail.tail_len,
        tail_clip_a_run=tail.clip_a_run,
        tail_walkback=tail.walkback,
        tail_class=tail.tail_class,
        clip_rna=tail.clip_rna,
        three_prime_shift=three_prime_corrected - three_prime_raw,
    )

    return record


def process_netseq_bam(
    bam_path: Path,
    exclusion_detector: Optional[ExclusionRegionDetector] = None,
    min_mapq: int = 0,
    min_a_fraction: float = 0.8,
    trim_terminal_oligo_a: bool = True,
    max_reads: Optional[int] = None,
    show_progress: bool = True,
    progress_interval: int = 1_000_000,
    rna3p_at: str = NETSEQ_RNA3P_DEFAULT,
    junction_pool=None,
    genome=None,
    umi_length: int = 0,
    rescue_max_intronic: int = 10,
    rescue_min_k: int = 1,
    detect_tail: bool = True,
) -> Generator[UnifiedReadRecord, None, None]:
    """
    Stream process NET-seq BAM file.

    Yields UnifiedReadRecord for each valid read, excluding:
    - Unmapped reads
    - Secondary/supplementary alignments
    - Reads in exclusion regions (if detector provided)
    - Reads below min_mapq

    Args:
        bam_path: Path to BAM file
        exclusion_detector: Optional ExclusionRegionDetector
        min_mapq: Minimum mapping quality (default 0)
        min_a_fraction: Minimum A fraction for oligo(A) detection
        trim_terminal_oligo_a: Trim terminal A/T mismatches (oligo(A) aligned as mismatch)
        max_reads: Limit for testing (None = process all)
        show_progress: Print progress messages
        progress_interval: Print progress every N reads
        rna3p_at: which read terminus is the RNA 3' end -- 'read5p' (NET-seq, default) or 'read3p' (legacy)
        junction_pool: optional JunctionPool enabling the donor-side rescue
        genome: optional {chrom: sequence}; required for the walkback and the rescue
        umi_length: declared randomer length (0 = none)
        rescue_max_intronic: how far past the donor the aligned end may sit (nt)
        rescue_min_k: minimum recovered exon-2 length for a rescue
        detect_tail: run the randomer-aware tail call

    Yields:
        UnifiedReadRecord for each valid read
    """
    bam_path = Path(bam_path)

    if show_progress:
        print(f"Processing NET-seq BAM: {bam_path.name}")

    # Build chromosome name mapping from BAM header
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    chrom_map = {}
    for sq in bam.header['SQ']:
        contig_name = sq['SN']
        chrom_map[contig_name] = standardize_netseq_chrom(contig_name)

    if show_progress:
        print(f"  Chromosome mapping: {len(chrom_map)} contigs")
        # Show a few examples
        for i, (orig, std) in enumerate(chrom_map.items()):
            if i < 3:
                print(f"    {orig[:50]}... -> {std}")

    # Process reads
    n_total = 0
    n_skipped_unmapped = 0
    n_skipped_secondary = 0
    n_skipped_mapq = 0
    n_skipped_excluded = 0
    n_yielded = 0
    n_trimmed_oligo_a = 0  # Reads with terminal oligo(A) mismatches trimmed

    try:
        for read in bam:
            n_total += 1

            # Progress reporting
            if show_progress and n_total % progress_interval == 0:
                print(f"  Processed {n_total:,} reads, yielded {n_yielded:,}...")

            # Skip invalid reads
            if read.is_unmapped:
                n_skipped_unmapped += 1
                continue
            if read.is_secondary or read.is_supplementary:
                n_skipped_secondary += 1
                continue
            if read.mapping_quality < min_mapq:
                n_skipped_mapq += 1
                continue

            # Get standardized chromosome name
            contig_name = read.reference_name
            chrom_std = chrom_map.get(contig_name, contig_name)

            # Get 3' position for exclusion check (with trimming)
            three_prime_pos, strand, _ = get_netseq_3prime_position(
                read, trim_terminal_oligo_a=trim_terminal_oligo_a, rna3p_at=rna3p_at
            )

            # Check exclusion regions
            if exclusion_detector and exclusion_detector.is_excluded(chrom_std, three_prime_pos):
                n_skipped_excluded += 1
                continue

            # Process read
            record = process_netseq_read(
                read, chrom_std, min_a_fraction,
                trim_terminal_oligo_a=trim_terminal_oligo_a,
                rna3p_at=rna3p_at,
                junction_pool=junction_pool,
                genome=genome,
                umi_length=umi_length,
                rescue_max_intronic=rescue_max_intronic,
                rescue_min_k=rescue_min_k,
                detect_tail=detect_tail,
            )
            n_yielded += 1

            # Track terminal oligo(A) trimming
            if record.aligned_a_length > 0:
                n_trimmed_oligo_a += 1

            yield record

            # Check max_reads limit
            if max_reads and n_yielded >= max_reads:
                if show_progress:
                    print(f"  Reached max_reads limit ({max_reads})")
                break
    finally:
        bam.close()

    if show_progress:
        print(f"  Completed: {n_total:,} total reads")
        print(f"    Yielded: {n_yielded:,}")
        print(f"    Skipped unmapped: {n_skipped_unmapped:,}")
        print(f"    Skipped secondary/supp: {n_skipped_secondary:,}")
        print(f"    Skipped low MAPQ: {n_skipped_mapq:,}")
        print(f"    Skipped excluded regions: {n_skipped_excluded:,}")
        if trim_terminal_oligo_a:
            pct = 100.0 * n_trimmed_oligo_a / n_yielded if n_yielded > 0 else 0
            print(f"    Terminal oligo(A) trimmed: {n_trimmed_oligo_a:,} ({pct:.1f}%)")


def aggregate_positions(
    records: Generator[UnifiedReadRecord, None, None],
    position_field: str = 'three_prime_raw',
) -> Dict[Tuple[str, str, int], int]:
    """
    Aggregate records into position counts.

    Args:
        records: Generator of UnifiedReadRecord
        position_field: Field to use for position ('three_prime_raw' or 'three_prime_corrected')

    Returns:
        Dict mapping (chrom, strand, position) -> count
    """
    counts: Dict[Tuple[str, str, int], int] = {}

    for record in records:
        position = getattr(record, position_field)
        key = (record.chrom, record.strand, position)
        counts[key] = counts.get(key, 0) + record.count

    return counts


def aggregate_netseq_stream(
    records,
    collect: bool = False,
):
    """One streaming pass over NET-seq records -> (raw_counts, corrected_counts, summary, records).

    ``aggregate_positions`` can only build ONE dict per generator, and the streaming path added in
    429f1c9 deliberately never materialises records (~23 GB per 50 M reads). Emitting BOTH a raw and
    a corrected track therefore has to happen in a single pass, together with the per-sample
    correction counters -- which is what this does.

    Args:
        records: iterable of records from :func:`process_netseq_bam`.
        collect: also return the records in a list (only when parquet is requested).

    Returns:
        ``(raw_counts, corrected_counts, summary, records_list)`` where the two count dicts are
        ``{(chrom, strand, position): count}`` and ``summary`` is a
        :class:`~rectify.core.netseq.netseq_rescue.NetseqCorrectionSummary`.
    """
    from ..netseq.netseq_rescue import NetseqCorrectionSummary

    raw_counts: Dict[Tuple[str, str, int], int] = {}
    corrected_counts: Dict[Tuple[str, str, int], int] = {}
    summary = NetseqCorrectionSummary()
    kept: List = [] if collect else []

    for record in records:
        raw_key = (record.chrom, record.strand, record.three_prime_raw)
        raw_counts[raw_key] = raw_counts.get(raw_key, 0) + record.count
        cor_key = (record.chrom, record.strand, record.three_prime_corrected)
        corrected_counts[cor_key] = corrected_counts.get(cor_key, 0) + record.count
        summary.observe(record)
        if collect:
            kept.append(record)

    return raw_counts, corrected_counts, summary, kept
