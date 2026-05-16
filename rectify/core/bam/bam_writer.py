#!/usr/bin/env python3
"""
File-level BAM writers for RECTIFY.

These writers orchestrate the per-read CIGAR/sequence editing primitives in
``read_edits.py`` over a whole input BAM. Each accepts ``corrected_tsv_path``
(the output of ``rectify correct``) and produces one or two output BAMs.

Write functions (preferred entry points):

- ``write_dual_bam``           — single-pass writer for both hardclip + softclip BAMs
- ``write_corrected_bam``      — hard-clips 3' end; extends 5' for Cat3
- ``write_softclipped_bam``    — soft-clips 3' end; extends 5' for Cat3
- ``write_polya_trimmed_bam``  — strips 3' poly(A) soft-clips

For per-read CIGAR helpers see ``rectify.core.bam.read_edits``.
For bedgraph emission see ``rectify.core.bam.bedgraph_writers``.

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, Optional, Tuple
import logging
import pysam

from .read_edits import (
    clip_read_to_corrected_3prime,
    softclip_read_to_corrected_3prime,
    extend_read_5prime_for_junction_rescue,
    extend_read_3prime_for_softclip_rescue,
    softclip_intronic_tail_5prime,
    reroute_intronic_tail_5prime_via_junction,
    realign_exon_blocks,
    _hardclip_trailing_a_run,
)

logger = logging.getLogger(__name__)


def _load_corrections_from_tsv(corrected_tsv_path: str) -> Dict[str, dict]:
    """
    Parse a corrected TSV and return a per-read correction dict.

    Only the first row per read_id is used (dominant position for Cat6
    multi-peak reads).  Returns a dict mapping read_id → correction info:

        ``corrected_3prime``       int
        ``strand``                 str  ('+' or '-')
        ``five_prime_position``    int or None
        ``five_prime_rescued``     bool
        ``five_prime_soft_clip``   int
    """
    corrections: Dict[str, dict] = {}
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

            # Optional columns — indices default to -1 if absent (old TSV).
            i_5p      = hdr.index('five_prime_position')         if 'five_prime_position'         in hdr else -1
            i_5p_res  = hdr.index('five_prime_rescued')          if 'five_prime_rescued'          in hdr else -1
            i_5p_clip = hdr.index('five_prime_soft_clip_length') if 'five_prime_soft_clip_length' in hdr else -1
            i_5p_cig  = hdr.index('five_prime_exon_cigar')       if 'five_prime_exon_cigar'       in hdr else -1
            # Cat2 soft-clip rescue columns (v2.9.1)
            i_sc_ext   = hdr.index('sc_homopolymer_extension')  if 'sc_homopolymer_extension'  in hdr else -1
            i_sc_seq   = hdr.index('sc_rescued_seq')             if 'sc_rescued_seq'             in hdr else -1
            i_sc_sclen = hdr.index('sc_original_softclip_len')  if 'sc_original_softclip_len'  in hdr else -1
            # Case 4 intronic-snap BAM hard-clip column (v2.9.8)
            i_5p_icp   = hdr.index('five_prime_intron_clip_pos') if 'five_prime_intron_clip_pos' in hdr else -1

            for line in _f:
                parts = line.rstrip('\n').split('\t')
                rid = parts[i_id] if len(parts) > i_id else None
                if rid is None or rid in corrections:
                    continue
                try:
                    corr_pos = int(parts[i_pos])
                    strand   = parts[i_strand]
                except (ValueError, IndexError):
                    continue

                five_prime_pos     = int(parts[i_5p])        if i_5p >= 0     and len(parts) > i_5p     and parts[i_5p]     else None
                five_prime_rescued = (parts[i_5p_res] == '1') if i_5p_res >= 0 and len(parts) > i_5p_res else False
                five_prime_sc      = int(parts[i_5p_clip])   if i_5p_clip >= 0 and len(parts) > i_5p_clip and parts[i_5p_clip] else 0
                five_prime_exon_cig = parts[i_5p_cig]        if i_5p_cig >= 0  and len(parts) > i_5p_cig  and parts[i_5p_cig]  else ''
                # Cat2 fields
                sc_ext   = int(parts[i_sc_ext])   if i_sc_ext   >= 0 and len(parts) > i_sc_ext   and parts[i_sc_ext]   else 0
                sc_seq   = parts[i_sc_seq]         if i_sc_seq   >= 0 and len(parts) > i_sc_seq   else ''
                sc_sclen = int(parts[i_sc_sclen])  if i_sc_sclen >= 0 and len(parts) > i_sc_sclen and parts[i_sc_sclen] else 0
                # Case 4 intronic-snap BAM clip position (-1 = not applicable)
                _icp_raw = parts[i_5p_icp] if i_5p_icp >= 0 and len(parts) > i_5p_icp and parts[i_5p_icp] else '-1'
                five_prime_icp = int(_icp_raw) if _icp_raw.lstrip('-').isdigit() else -1

                corrections[rid] = {
                    'corrected_3prime':           corr_pos,
                    'strand':                     strand,
                    'five_prime_position':        five_prime_pos,
                    'five_prime_rescued':         five_prime_rescued,
                    'five_prime_soft_clip':       five_prime_sc,
                    'five_prime_exon_cigar':      five_prime_exon_cig,
                    'five_prime_intron_clip_pos': five_prime_icp,
                    'sc_homopolymer_extension':   sc_ext,
                    'sc_rescued_seq':             sc_seq,
                    'sc_original_softclip_len':   sc_sclen,
                }
    except OSError as exc:
        raise OSError(
            f"Cannot read corrected TSV {corrected_tsv_path}: {exc}"
        ) from exc
    return corrections


def write_corrected_bam(
    input_bam_path: str,
    corrected_tsv_path: str,
    output_bam_path: str,
    genome: Optional[Dict[str, str]] = None,
) -> Dict[str, int]:
    """
    Write a new BAM with every read hard-clipped at its corrected 3' end and,
    for Cat3 junction-rescued reads, extended at its corrected 5' end.

    Reads per-read corrections from *corrected_tsv_path* (the TSV produced by
    ``rectify correct``) and applies:

    0. **Junction homopolymer surgery** (when *genome* provided):
       :func:`fix_homopolymer_mismatches` converts N-X patterns
       at homopolymer / splice-site boundaries to the more parsimonious N+1 form.
    1. **5' extension** (Cat3): :func:`extend_read_5prime_for_junction_rescue`
       converts the 5' soft-clip to M+N+M to show the exon-intron-exon structure.
    2. **3' hard-clip** (all other categories): :func:`clip_read_to_corrected_3prime`.

    Only the first row per read is used (covers Cat6 NET-seq multi-peak reads).
    Reads absent from the TSV are written unchanged.

    Args:
        input_bam_path:      Path to the original input BAM.
        corrected_tsv_path:  Path to the corrected_reads.tsv from ``rectify correct``.
        output_bam_path:     Destination BAM path.  Overwritten if it exists.
        genome:              Optional pre-loaded genome dict for homopolymer surgery.

    Returns:
        Dict with summary counts: ``'total'``, ``'clipped'``, ``'unchanged'``.
    """
    corrections = _load_corrections_from_tsv(corrected_tsv_path)

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

            modified = False

            # Homopolymer CIGAR surgery: re-align each exon block with reduced
            # mismatch penalty at homopolymer ref positions so that nanopore DRS
            # length errors are represented as indels rather than mismatches.
            # Global alignment preserves query and reference span by construction.
            if genome is not None:
                modified |= realign_exon_blocks(read, genome)

            # 5' junction rescue: extend soft-clip to exon 1 (Cat3).
            if correction['five_prime_rescued'] and correction['five_prime_position'] is not None:
                modified |= extend_read_5prime_for_junction_rescue(
                    read,
                    correction['five_prime_position'],
                    correction['five_prime_soft_clip'],
                    correction['strand'],
                    exon_cigar_str=correction.get('five_prime_exon_cigar', ''),
                )

            # 5' junction rescue: reroute intronic M ops to exon 1 (Cases 1/2/2b/4).
            # Fires for reads with no soft-clip but a non-empty five_prime_exon_cigar
            # (the aligner mapped into the intron using M/X/D ops rather than N).
            _icp = correction.get('five_prime_intron_clip_pos', -1)
            _exon_cig = correction.get('five_prime_exon_cigar', '')
            if (_icp >= 0 and _exon_cig and correction.get('five_prime_rescued')
                    and correction['five_prime_position'] is not None):
                modified |= reroute_intronic_tail_5prime_via_junction(
                    read,
                    clip_boundary=_icp,
                    five_prime_position=correction['five_prime_position'],
                    exon_cigar_str=_exon_cig,
                    strand=correction['strand'],
                )
            elif _icp >= 0 and not _exon_cig and correction.get('five_prime_rescued'):
                # Fallback: no exon CIGAR available (local alignment failed or
                # intronic sequence too short).  Soft-clip the intronic tail at
                # clip_boundary so the alignment ends at the exon/intron boundary.
                modified |= softclip_intronic_tail_5prime(
                    read,
                    clip_boundary=_icp,
                    strand=correction['strand'],
                )

            # Cat2 soft-clip rescue: extend 3' alignment outward into homopolymer.
            if correction.get('sc_rescued_seq'):
                modified |= extend_read_3prime_for_softclip_rescue(
                    read,
                    correction['strand'],
                    correction['sc_homopolymer_extension'],
                    correction['sc_rescued_seq'],
                    correction['sc_original_softclip_len'],
                    hard_clip=True,
                )

            # 3' correction: hard-clip to corrected position (Cat1/4/5/6).
            # Cat2 reads already have their 3' end extended outward; this call
            # will return False for them (corrected >= current_end).
            modified |= clip_read_to_corrected_3prime(
                read, correction['corrected_3prime'], correction['strand']
            )

            # Additional hard-clip: remove any trailing genomic A-run at the 3' end.
            # After the main correction, A-rich 3' UTR regions are indistinguishable
            # from poly(A) tail in the read sequence and are clipped away here.
            # Not applied to the soft-clip BAM, which retains them as aligned bases.
            modified |= _hardclip_trailing_a_run(read, correction['strand'])

            # Tag the final corrected 3' end so it is readable in IGV / samtools view.
            # cp:i: = corrected 3' end position, 0-based inclusive reference coordinate.
            # This equals reference_start (minus strand) or reference_end-1 (plus strand)
            # after clipping, but is written explicitly so it is visible even for reads
            # where no correction was applied (corrected == original).
            read.set_tag('cp', correction['corrected_3prime'])

            bam_out.write(read)
            if modified:
                stats['clipped'] += 1
            else:
                stats['unchanged'] += 1

    return stats


def write_softclipped_bam(
    input_bam_path: str,
    corrected_tsv_path: str,
    output_bam_path: str,
    genome: Optional[Dict[str, str]] = None,
) -> Dict[str, int]:
    """
    Write a new BAM with every read soft-clipped at its corrected 3' end and,
    for Cat3 junction-rescued reads, extended at its corrected 5' end.

    Identical to :func:`write_corrected_bam` except the 3' correction uses
    :func:`softclip_read_to_corrected_3prime` (``S`` ops, bases retained) rather
    than hard-clip.  The 5' extension always uses ``M`` ops in both variants
    since the rescued exon bases are aligned (not soft-clipped).
    """
    corrections = _load_corrections_from_tsv(corrected_tsv_path)

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

            modified = False

            # Homopolymer CIGAR surgery: re-align exon blocks.
            if genome is not None:
                modified |= realign_exon_blocks(read, genome)

            # 5' junction rescue: extend into upstream exon (Cat3).
            if correction['five_prime_rescued'] and correction['five_prime_position'] is not None:
                modified |= extend_read_5prime_for_junction_rescue(
                    read,
                    correction['five_prime_position'],
                    correction['five_prime_soft_clip'],
                    correction['strand'],
                    exon_cigar_str=correction.get('five_prime_exon_cigar', ''),
                )

            # 5' junction rescue: reroute intronic M ops to exon 1 (Cases 1/2/2b/4).
            _icp = correction.get('five_prime_intron_clip_pos', -1)
            _exon_cig = correction.get('five_prime_exon_cigar', '')
            if (_icp >= 0 and _exon_cig and correction.get('five_prime_rescued')
                    and correction['five_prime_position'] is not None):
                modified |= reroute_intronic_tail_5prime_via_junction(
                    read,
                    clip_boundary=_icp,
                    five_prime_position=correction['five_prime_position'],
                    exon_cigar_str=_exon_cig,
                    strand=correction['strand'],
                )
            elif _icp >= 0 and not _exon_cig and correction.get('five_prime_rescued'):
                modified |= softclip_intronic_tail_5prime(
                    read,
                    clip_boundary=_icp,
                    strand=correction['strand'],
                )

            # Cat2 soft-clip rescue: extend 3' alignment outward into homopolymer.
            if correction.get('sc_rescued_seq'):
                modified |= extend_read_3prime_for_softclip_rescue(
                    read,
                    correction['strand'],
                    correction['sc_homopolymer_extension'],
                    correction['sc_rescued_seq'],
                    correction['sc_original_softclip_len'],
                    hard_clip=False,
                )

            # 3' correction: soft-clip to corrected position (Cat1/4/5/6).
            # Cat2 reads are already extended outward; this call returns False.
            modified |= softclip_read_to_corrected_3prime(
                read, correction['corrected_3prime'], correction['strand']
            )

            read.set_tag('cp', correction['corrected_3prime'])

            bam_out.write(read)
            if modified:
                stats['clipped'] += 1
            else:
                stats['unchanged'] += 1

    return stats


def write_dual_bam(
    input_bam_path: str,
    corrected_tsv_path: str,
    output_hardclip_path: str,
    output_softclip_path: str,
    genome: Optional[Dict[str, str]] = None,
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Write hardclip and softclip corrected BAMs in a single pass over the input BAM.

    Equivalent to calling :func:`write_corrected_bam` and
    :func:`write_softclipped_bam` sequentially, but reads the input BAM only
    once — saving ~50% input I/O compared to the two-function approach.

    The two output BAMs differ only in how the 3' poly(A) region is encoded:

    * *hardclip* (``rectified_corrected_3end.bam``): poly(A) bases removed from
      ``query_sequence`` and stored as ``H`` ops — invisible in IGV pileups.
    * *softclip* (``rectified_pA_tail_trimmed.bam``): poly(A) bases retained in
      ``query_sequence`` and stored as ``S`` ops — visible in IGV when
      "Show soft-clipped bases" is enabled.

    The shared pre-pass (homopolymer junction surgery + 5' rescue / Cat3) is
    applied before state diverges, so it is only computed once per read.
    For the diverging 3' operations the read's mutable state (CIGAR, sequence,
    qualities, reference_start) is saved after shared operations and restored
    before the softclip path.

    Args:
        input_bam_path:       Path to the original input BAM.
        corrected_tsv_path:   Path to the corrected_reads.tsv from ``rectify correct``.
        output_hardclip_path: Destination BAM path for the hardclip variant.
        output_softclip_path: Destination BAM path for the softclip variant.
        genome:               Optional pre-loaded genome dict for homopolymer surgery.

    Returns:
        Tuple of (hardclip_stats, softclip_stats) dicts, each with keys
        ``'total'``, ``'clipped'``, ``'unchanged'``.
    """
    corrections = _load_corrections_from_tsv(corrected_tsv_path)

    logger.info(
        "write_dual_bam: loaded %d corrected positions from %s",
        len(corrections), corrected_tsv_path,
    )

    hc_stats: Dict[str, int] = {'total': 0, 'clipped': 0, 'unchanged': 0}
    sc_stats: Dict[str, int] = {'total': 0, 'clipped': 0, 'unchanged': 0}

    with pysam.AlignmentFile(input_bam_path, 'rb') as bam_in, \
         pysam.AlignmentFile(output_hardclip_path, 'wb', header=bam_in.header) as bam_hc, \
         pysam.AlignmentFile(output_softclip_path, 'wb', header=bam_in.header) as bam_sc:

        for read in bam_in:
            hc_stats['total'] += 1
            sc_stats['total'] += 1

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                bam_hc.write(read)
                bam_sc.write(read)
                hc_stats['unchanged'] += 1
                sc_stats['unchanged'] += 1
                continue

            correction = corrections.get(read.query_name)
            if correction is None:
                bam_hc.write(read)
                bam_sc.write(read)
                hc_stats['unchanged'] += 1
                sc_stats['unchanged'] += 1
                continue

            # Apply shared pre-pass — identical for both BAMs.
            shared_modified = False

            # Homopolymer CIGAR surgery: re-align exon blocks.
            if genome is not None:
                shared_modified |= realign_exon_blocks(read, genome)

            # 5' rescue (Cat3) — identical for both BAMs.
            if correction['five_prime_rescued'] and correction['five_prime_position'] is not None:
                shared_modified |= extend_read_5prime_for_junction_rescue(
                    read,
                    correction['five_prime_position'],
                    correction['five_prime_soft_clip'],
                    correction['strand'],
                    exon_cigar_str=correction.get('five_prime_exon_cigar', ''),
                )

            # 5' junction rescue: reroute intronic M ops to exon 1 (Cases 1/2/2b/4).
            _icp = correction.get('five_prime_intron_clip_pos', -1)
            _exon_cig = correction.get('five_prime_exon_cigar', '')
            if (_icp >= 0 and _exon_cig and correction.get('five_prime_rescued')
                    and correction['five_prime_position'] is not None):
                shared_modified |= reroute_intronic_tail_5prime_via_junction(
                    read,
                    clip_boundary=_icp,
                    five_prime_position=correction['five_prime_position'],
                    exon_cigar_str=_exon_cig,
                    strand=correction['strand'],
                )
            elif _icp >= 0 and not _exon_cig and correction.get('five_prime_rescued'):
                shared_modified |= softclip_intronic_tail_5prime(
                    read,
                    clip_boundary=_icp,
                    strand=correction['strand'],
                )

            # Save state at the divergence point (after shared ops, before 3' ops).
            # Only cigar, seq, quals, and reference_start are mutated by any path.
            saved_cigar     = list(read.cigartuples or [])
            saved_seq       = read.query_sequence
            saved_quals_arr = read.query_qualities          # numpy array or None
            saved_refstart  = read.reference_start

            # ── Hardclip path ────────────────────────────────────────────────
            hc_modified = shared_modified
            if correction.get('sc_rescued_seq'):
                hc_modified |= extend_read_3prime_for_softclip_rescue(
                    read,
                    correction['strand'],
                    correction['sc_homopolymer_extension'],
                    correction['sc_rescued_seq'],
                    correction['sc_original_softclip_len'],
                    hard_clip=True,
                )
            hc_modified |= clip_read_to_corrected_3prime(
                read, correction['corrected_3prime'], correction['strand']
            )
            # Additional hard-clip: remove any trailing genomic A-run at the 3' end.
            # Genomic A-rich 3' UTR regions are indistinguishable from poly(A) tail
            # in the read sequence; hard-clip removes them for an unambiguous view.
            # The soft-clip path below does NOT apply this — it retains them as aligned.
            hc_modified |= _hardclip_trailing_a_run(read, correction['strand'])
            read.set_tag('cp', correction['corrected_3prime'])
            bam_hc.write(read)
            if hc_modified:
                hc_stats['clipped'] += 1
            else:
                hc_stats['unchanged'] += 1

            # ── Restore state before softclip path ───────────────────────────
            read.cigartuples     = saved_cigar
            read.query_sequence  = saved_seq
            read.query_qualities = saved_quals_arr
            read.reference_start = saved_refstart

            # ── Softclip path ────────────────────────────────────────────────
            sc_modified = shared_modified
            if correction.get('sc_rescued_seq'):
                sc_modified |= extend_read_3prime_for_softclip_rescue(
                    read,
                    correction['strand'],
                    correction['sc_homopolymer_extension'],
                    correction['sc_rescued_seq'],
                    correction['sc_original_softclip_len'],
                    hard_clip=False,
                )
            sc_modified |= softclip_read_to_corrected_3prime(
                read, correction['corrected_3prime'], correction['strand']
            )
            read.set_tag('cp', correction['corrected_3prime'])
            bam_sc.write(read)
            if sc_modified:
                sc_stats['clipped'] += 1
            else:
                sc_stats['unchanged'] += 1

    return hc_stats, sc_stats


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
    from ..polya.polya_trimmer import trim_polya_from_bam_read

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
