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

from ...utils.genome import get_chrom_sequence
from .read_edits import (
    clip_read_to_corrected_3prime,
    softclip_read_to_corrected_3prime,
    extend_read_5prime_for_junction_rescue,
    extend_read_3prime_for_overcall_rescue,
    extend_read_3prime_for_softclip_rescue,
    softclip_intronic_tail_5prime,
    reroute_intronic_tail_5prime_via_junction,
    realign_exon_blocks,
    reanchor_5prime_for_rescue,
    _hardclip_trailing_a_run,
)

logger = logging.getLogger(__name__)


def _decode_eq_seq_inplace(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
) -> bool:
    """Replace ``=`` characters in ``read.query_sequence`` with the genome
    base at the corresponding reference position.

    BAM/SAM spec allows ``SEQ`` to contain ``=`` as a shorthand for "matches
    reference at this position". Aligners that emit this compressed form
    (minimap2, gapmm2, deSALT, uLTRA on long ``=``-CIGAR runs) propagate
    the ``=`` chars through to downstream BAMs unless explicitly decoded.
    Downstream consumers that read ``query_sequence`` and compare to the
    reference (notably ``_cigar_hp_edit_distance`` in winner-selection)
    misinterpret ``=`` as a literal mismatch and score every M-block base
    as wrong — silently biasing winner-selection against any aligner that
    used the compressed form.

    Decoding at BAM-write time produces explicit-SEQ BAMs that all consumers
    can read uniformly. Bloat is typically +20–40 % on minimap2-family BAMs
    (perfect-match runs compress poorly once expanded) and zero on aligners
    that already emit explicit bases. Returns True if any ``=`` was decoded.

    Only M/=/X CIGAR ops can carry ``=`` in SEQ per spec (those positions
    have a defined ref base); I, S positions are passed through untouched.
    """
    if read.is_unmapped:
        return False
    seq = read.query_sequence
    if seq is None or '=' not in seq:
        return False
    chrom_seq, _chrom_key = get_chrom_sequence(genome, read.reference_name) if genome else (None, None)
    if chrom_seq is None:
        return False
    new_chars = list(seq)
    ref_pos = read.reference_start
    q_pos = 0
    decoded_any = False
    for op, length in read.cigartuples or []:
        if op in (0, 7, 8):                # M, =, X — consume ref + query
            for i in range(length):
                if new_chars[q_pos + i] == '=' and ref_pos + i < len(chrom_seq):
                    new_chars[q_pos + i] = chrom_seq[ref_pos + i].upper()
                    decoded_any = True
            ref_pos += length
            q_pos += length
        elif op in (1, 4):                 # I, S — consume query only
            q_pos += length
        elif op in (2, 3):                 # D, N — consume ref only
            ref_pos += length
        # H (5), P (6) — consume neither
    if not decoded_any:
        return False
    # Re-assigning query_sequence clobbers query_qualities; save and restore.
    saved_qual = read.query_qualities
    read.query_sequence = ''.join(new_chars)
    if saved_qual is not None:
        read.query_qualities = saved_qual
    return True


def _load_corrections_from_single_tsv(corrected_tsv_path: str) -> Dict[str, dict]:
    """Load corrections from a single (non-manifest) corrected_reads TSV.

    Internal helper — callers should use :func:`_load_corrections_from_tsv`
    which auto-detects manifest vs single-TSV format.
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
            i_5p_trim = hdr.index('five_prime_upstream_trim')    if 'five_prime_upstream_trim'    in hdr else -1
            i_5p_reanc = hdr.index('reanchor_clip_len')          if 'reanchor_clip_len'           in hdr else -1
            # Cat2 soft-clip rescue columns (v2.9.1)
            i_sc_ext   = hdr.index('sc_homopolymer_extension')  if 'sc_homopolymer_extension'  in hdr else -1
            i_sc_seq   = hdr.index('sc_rescued_seq')             if 'sc_rescued_seq'             in hdr else -1
            i_sc_sclen = hdr.index('sc_original_softclip_len')  if 'sc_original_softclip_len'  in hdr else -1
            # Case 4 intronic-snap BAM hard-clip column (v2.9.8)
            i_5p_icp   = hdr.index('five_prime_intron_clip_pos') if 'five_prime_intron_clip_pos' in hdr else -1
            # Over-call rescue columns
            i_oc_ext   = hdr.index('oc_homopolymer_extension')   if 'oc_homopolymer_extension'   in hdr else -1
            i_oc_cnt   = hdr.index('oc_overcall_count')          if 'oc_overcall_count'          in hdr else -1
            i_oc_term  = hdr.index('oc_terminal_base')           if 'oc_terminal_base'           in hdr else -1

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
                five_prime_trim    = int(parts[i_5p_trim])   if i_5p_trim >= 0 and len(parts) > i_5p_trim and parts[i_5p_trim] else 0
                five_prime_reanc   = int(parts[i_5p_reanc])  if i_5p_reanc >= 0 and len(parts) > i_5p_reanc and parts[i_5p_reanc] else 0
                # Cat2 fields
                sc_ext   = int(parts[i_sc_ext])   if i_sc_ext   >= 0 and len(parts) > i_sc_ext   and parts[i_sc_ext]   else 0
                sc_seq   = parts[i_sc_seq]         if i_sc_seq   >= 0 and len(parts) > i_sc_seq   else ''
                sc_sclen = int(parts[i_sc_sclen])  if i_sc_sclen >= 0 and len(parts) > i_sc_sclen and parts[i_sc_sclen] else 0
                # Case 4 intronic-snap BAM clip position (-1 = not applicable)
                _icp_raw = parts[i_5p_icp] if i_5p_icp >= 0 and len(parts) > i_5p_icp and parts[i_5p_icp] else '-1'
                five_prime_icp = int(_icp_raw) if _icp_raw.lstrip('-').isdigit() else -1
                # Over-call rescue fields
                oc_ext  = int(parts[i_oc_ext])  if i_oc_ext  >= 0 and len(parts) > i_oc_ext  and parts[i_oc_ext]  else 0
                oc_cnt  = int(parts[i_oc_cnt])  if i_oc_cnt  >= 0 and len(parts) > i_oc_cnt  and parts[i_oc_cnt]  else 0
                oc_term = parts[i_oc_term]      if i_oc_term >= 0 and len(parts) > i_oc_term else ''

                corrections[rid] = {
                    'corrected_3prime':           corr_pos,
                    'strand':                     strand,
                    'five_prime_position':        five_prime_pos,
                    'five_prime_rescued':         five_prime_rescued,
                    'five_prime_soft_clip':       five_prime_sc,
                    'five_prime_exon_cigar':      five_prime_exon_cig,
                    'five_prime_upstream_trim':   five_prime_trim,
                    'reanchor_clip_len':          five_prime_reanc,
                    'five_prime_intron_clip_pos': five_prime_icp,
                    'sc_homopolymer_extension':   sc_ext,
                    'sc_rescued_seq':             sc_seq,
                    'sc_original_softclip_len':   sc_sclen,
                    'oc_homopolymer_extension':   oc_ext,
                    'oc_overcall_count':          oc_cnt,
                    'oc_terminal_base':           oc_term,
                }
    except OSError as exc:
        raise OSError(
            f"Cannot read corrected TSV {corrected_tsv_path}: {exc}"
        ) from exc
    return corrections


# ---------------------------------------------------------------------------
# Manifest-aware public loader (Commit B, 2026-05-20)
# ---------------------------------------------------------------------------

_MANIFEST_HEADER_COLS = ['region_id', 'chrom', 'start', 'end', 'tsv_path', 'n_rows', 'sha256']

def _load_corrections_from_tsv(corrected_tsv_path: str) -> Dict[str, dict]:
    """Load corrections from a corrected_reads.tsv OR a corrected_reads.manifest.tsv.

    Auto-detects which format the input is by reading the first header line:

    * **Manifest** header: ``region_id\\tchrom\\tstart\\tend\\ttsv_path\\tn_rows\\tsha256``
    * **Single TSV** header: rectify-correct row header (``read_id\\tchrom\\t...``)

    Both code paths ultimately return the same ``{read_id: correction_dict}``
    mapping so ALL callers (``write_corrected_bam``, ``write_softclipped_bam``,
    ``write_dual_bam``, ``bam_writer_parallel``) work unchanged.

    Args:
        corrected_tsv_path:  Path to either a ``corrected_reads.tsv`` or a
            ``corrected_reads.manifest.tsv`` produced by ``rectify correct``.

    Returns:
        Dict mapping read_id → correction info dict.
    """
    from pathlib import Path as _Path
    p = _Path(corrected_tsv_path)
    with p.open() as fh:
        header = fh.readline().rstrip('\n').split('\t')

    if header == _MANIFEST_HEADER_COLS:
        # Manifest format — load per-region TSVs and merge.
        from .tsv_partition import load_manifest
        entries = load_manifest(p)
        merged: Dict[str, dict] = {}
        for entry in entries:
            region_tsv = entry['tsv_path']  # already absolute Path from load_manifest
            merged.update(_load_corrections_from_single_tsv(str(region_tsv)))
        return merged
    else:
        # Legacy single-TSV format.
        return _load_corrections_from_single_tsv(corrected_tsv_path)


def apply_corrected_edits_to_read(
    read: pysam.AlignedSegment,
    correction: Optional[Dict],
    genome: Optional[Dict[str, str]] = None,
) -> bool:
    """Apply the canonical hard-clipped corrected-BAM edits to one read.

    This mutates *read* in place but does not write it.  It is the shared
    implementation used by the sequential writer, parallel writer, lazy
    HP-edit-distance scoring, and final corrected-consensus BAM generation.

    Returns True when correction surgery changed the alignment/CIGAR.  Decoding
    SAM-spec ``=`` shorthand in SEQ is intentionally not counted as "modified",
    preserving the legacy writer's clipped/unchanged statistics.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False

    # Decode '='-compressed SEQ so every downstream scorer/writer sees explicit
    # bases.  This mirrors the legacy write_corrected_bam pre-pass.
    if genome is not None:
        _decode_eq_seq_inplace(read, genome)

    if correction is None:
        return False

    modified = False

    # 5'-edge reanchor pre-pass: when bam_processor's 3'SS rescue used a
    # reanchored copy of the read (TSV reanchor_clip_len > 0), apply the same
    # deterministic reanchor here BEFORE realign_exon_blocks so the live CIGAR
    # matches the geometry that exon_cigar was sized for.
    if genome is not None and correction.get('reanchor_clip_len', 0) > 0:
        modified |= reanchor_5prime_for_rescue(read, genome, anchor_min_run=10)

    # Homopolymer CIGAR surgery: re-align exon blocks with X ops at homopolymer
    # positions so under-called DRS homopolymers are represented as indels.
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
            upstream_trim=correction.get('five_prime_upstream_trim', 0),
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
            hard_clip=True,
        )

    # Over-call rescue (cat1_plus_1/2): convert terminal overcall geometry.
    if correction.get('oc_terminal_base'):
        modified |= extend_read_3prime_for_overcall_rescue(
            read,
            correction['strand'],
            correction['oc_homopolymer_extension'],
            correction['oc_overcall_count'],
            correction['oc_terminal_base'],
        )

    # 3' correction: hard-clip to corrected position (Cat1/4/5/6).
    modified |= clip_read_to_corrected_3prime(
        read, correction['corrected_3prime'], correction['strand']
    )

    # Additional hard-clip: remove any trailing genomic A-run at the 3' end.
    modified |= _hardclip_trailing_a_run(read, correction['strand'])

    # Tag the final corrected 3' end so it is visible in IGV / samtools view.
    read.set_tag('cp', correction['corrected_3prime'])

    return modified


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

            correction = None
            if not (read.is_unmapped or read.is_secondary or read.is_supplementary):
                correction = corrections.get(read.query_name)
            modified = apply_corrected_edits_to_read(read, correction, genome)
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

            # Decode '='-compressed SEQ so every emitted read has explicit
            # bases (see _decode_eq_seq_inplace docstring).
            if genome is not None:
                _decode_eq_seq_inplace(read, genome)

            correction = corrections.get(read.query_name)
            if correction is None:
                bam_out.write(read)
                stats['unchanged'] += 1
                continue

            modified = False

            # 5'-edge reanchor pre-pass (see write_corrected_bam for rationale).
            if genome is not None and correction.get('reanchor_clip_len', 0) > 0:
                modified |= reanchor_5prime_for_rescue(read, genome, anchor_min_run=10)

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
                    upstream_trim=correction.get('five_prime_upstream_trim', 0),
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

            # Over-call rescue (cat1_plus_1/2): consume the soft-clip into
            # `{hp_ext}D {overcall_count}I 1=` extending the alignment by 1
            # ref position to capture the terminal body match.
            if correction.get('oc_terminal_base'):
                modified |= extend_read_3prime_for_overcall_rescue(
                    read,
                    correction['strand'],
                    correction['oc_homopolymer_extension'],
                    correction['oc_overcall_count'],
                    correction['oc_terminal_base'],
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

            # Decode '='-compressed SEQ so both emitted BAMs carry explicit
            # bases (see _decode_eq_seq_inplace docstring).
            if genome is not None:
                _decode_eq_seq_inplace(read, genome)

            correction = corrections.get(read.query_name)
            if correction is None:
                bam_hc.write(read)
                bam_sc.write(read)
                hc_stats['unchanged'] += 1
                sc_stats['unchanged'] += 1
                continue

            # Apply shared pre-pass — identical for both BAMs.
            shared_modified = False

            # 5'-edge reanchor pre-pass (see write_corrected_bam for rationale).
            if genome is not None and correction.get('reanchor_clip_len', 0) > 0:
                shared_modified |= reanchor_5prime_for_rescue(read, genome, anchor_min_run=10)

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
                    upstream_trim=correction.get('five_prime_upstream_trim', 0),
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
            if correction.get('oc_terminal_base'):
                hc_modified |= extend_read_3prime_for_overcall_rescue(
                    read,
                    correction['strand'],
                    correction['oc_homopolymer_extension'],
                    correction['oc_overcall_count'],
                    correction['oc_terminal_base'],
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
            if correction.get('oc_terminal_base'):
                sc_modified |= extend_read_3prime_for_overcall_rescue(
                    read,
                    correction['strand'],
                    correction['oc_homopolymer_extension'],
                    correction['oc_overcall_count'],
                    correction['oc_terminal_base'],
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
