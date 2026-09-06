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
import copy
import logging
import pysam

from ...utils.genome import get_chrom_sequence
from ..splice.overhang_informativeness import is_canonical_junction
from .read_edits import (
    clip_read_to_corrected_3prime,
    softclip_read_to_corrected_3prime,
    extend_read_5prime_for_junction_rescue,
    extend_read_3prime_for_overcall_rescue,
    extend_read_3prime_for_softclip_rescue,
    softclip_intronic_tail_5prime,
    reroute_intronic_tail_5prime_via_junction,
    projected_5prime_rescue_intron_edge,
    realign_exon_blocks,
    _apply_reanchor_from_clip_len,
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
            i_5p_e2    = hdr.index('five_prime_exon2_prefix')    if 'five_prime_exon2_prefix'     in hdr else -1
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
                five_prime_e2      = int(parts[i_5p_e2])     if i_5p_e2 >= 0 and len(parts) > i_5p_e2 and parts[i_5p_e2] else 0
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
                    'five_prime_exon2_prefix':    five_prime_e2,
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


def _n_op_intervals(read: pysam.AlignedSegment) -> Tuple[Tuple[int, int], ...]:
    """(start, end) genomic interval of every N-op in the live CIGAR.

    0-based half-open, matching the pool/annotation convention.
    """
    if not read.cigartuples or read.reference_start is None:
        return ()
    out = []
    pos = read.reference_start
    for op, length in read.cigartuples:
        if op == 3:
            out.append((pos, pos + length))
            pos += length
        elif op in (0, 2, 7, 8):
            pos += length
    return tuple(out)


# A 5' rescue whose N-op overshoots the acceptor by a base or two is the read
# MISSING those bases at the exon-2 start, not an invented junction: `extend`
# ends its N-op at the read's own alignment start, so every exon-2 base the
# basecaller dropped is swallowed into the intron. Pulling the acceptor back onto
# the canonical position and representing the dropped bases as D recovers the
# real junction. Bounded deliberately: beyond a few bases the disagreement is
# about placement, not about lost terminal bases, and the junction should be
# refused rather than massaged.
_MAX_ACCEPTOR_REPAIR_BP = 3


def _repair_acceptor_overshoot(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    n_index: int,
    start: int,
    end: int,
) -> bool:
    """Pull a non-canonical rescued N-op back onto the nearest canonical acceptor.

    ``n_index`` is the position of the N op in ``read.cigartuples``. The recovered
    reference bases become a D op on the BODY side of the N — which is the side
    the acceptor abuts on either strand: for a plus-strand rescue the CIGAR is
    ``[exon..., N, body...]`` and the acceptor is the N's HIGH edge, for a minus-
    strand rescue it is ``[body..., N, exon...]`` and the acceptor is the LOW
    edge. Returns True when a repair was applied.

    The smallest qualifying shift wins, so a nearer canonical position is never
    passed over for a farther one (on the bundled upf1d cat3_plus_2 read both
    delta=2 — the annotated GT-AG — and delta=4 are canonical).
    """
    cigar = list(read.cigartuples or [])
    if not (0 <= n_index < len(cigar)) or cigar[n_index][0] != 3:
        return False
    length = cigar[n_index][1]
    # Which edge abuts the read body? The exon ops the rescue wrote sit on the
    # other side of the N.
    body_is_after = any(op in (0, 7, 8) for op, _l in cigar[n_index + 1:])
    for delta in range(1, _MAX_ACCEPTOR_REPAIR_BP + 1):
        if length - delta <= 0:
            break
        cand = (start, end - delta) if body_is_after else (start + delta, end)
        if not is_canonical_junction(chrom_seq, cand[0], cand[1], atac=True):
            continue
        new_cigar = list(cigar)
        new_cigar[n_index] = (3, length - delta)
        if body_is_after:
            new_cigar.insert(n_index + 1, (2, delta))
        else:
            new_cigar.insert(n_index, (2, delta))
            # reference_start is unchanged: the D consumes reference that the N
            # used to cover, on the body side.
        read.cigartuples = new_cigar
        logger.debug(
            "5' rescue acceptor repaired for %s: N %d-%d -> %d-%d (+%dD)",
            read.query_name, start, end, cand[0], cand[1], delta,
        )
        return True
    return False


def _revert_selfinflicted_noncanonical_n(
    read: pysam.AlignedSegment,
    genome: Optional[Dict[str, str]],
    pre_nops: Tuple[Tuple[int, int], ...],
    pre_cigar,
    pre_start: Optional[int],
) -> bool:
    """Undo a 5' rescue surgery that INVENTED a non-canonical intron.

    Returns True when the read was reverted.

    This is a guard against a **self-inflicted junction**, not a motif filter on
    aligner evidence.  Only N-ops that were absent from ``pre_nops`` — i.e. that
    this writer's own 5' surgery drew — are judged; an N-op the aligner already
    reported is never touched, whatever its motif, because RECTIFY's job is to
    place the read's evidence, not to censor the aligner's.  The rescue that
    requested the surgery scored a junction from the annotation pool; the writer
    can end up drawing the N-op at *different* coordinates (ISSUE-002: `extend`
    derives its intron length from the live alignment edge and never reads
    ``five_prime_intron_clip_pos``; ISSUE-007: the reroute merge can fuse a
    partially trimmed N with the new one).  Nothing scored the motif at the
    coordinates that actually reach the BAM — this does.

    "Canonical" is the human/yeast splice grammar GT-AG, GC-AG and AT-AC
    (:func:`overhang_informativeness.is_canonical_junction` with ``atac=True``,
    which encodes all six forward-genome dinucleotide pairs, both strands).
    Motivation: on the Sumner human RNA004 chr5 slice the 5' rescue added
    ~24.9k N-ops of which ~74% were non-canonical, and the tester's clip-length
    table shows a length floor alone is not enough — even 15-29 nt rescues are
    40-58% non-canonical — so the floor is paired with this destination guard.

    A read whose genome sequence cannot be resolved is left alone (the guard
    cannot judge what it cannot read); that is logged at debug level so a
    silent no-op is not mistaken for "nothing to revert".
    """
    if genome is None or pre_start is None:
        return False
    new_nops = [n for n in _n_op_intervals(read) if n not in pre_nops]
    if not new_nops:
        return False
    chrom_seq, _key = get_chrom_sequence(genome, read.reference_name)
    if not chrom_seq:
        logger.debug(
            "5' rescue canonical guard skipped for %s: no genome sequence for %s",
            read.query_name, read.reference_name,
        )
        return False
    for start, end in new_nops:
        if is_canonical_junction(chrom_seq, start, end, atac=True):
            continue
        # A near miss is the read MISSING a base or two at the exon-2 start, not
        # an invented junction. Try to pull the acceptor back onto the canonical
        # position first; only a junction that cannot be recovered that way is
        # reverted. (Bundled upf1d cat3_plus_2: extend drew GT-TC 148194-148284
        # because the read lacks the 2 exon-2 bases TC, and shrinking by 2
        # recovers the ANNOTATED GT-AG 148194-148282. Reverting instead cost that
        # read its rescue and, through merge_corrected_tsvs' HP-edit-distance
        # scoring of the corrected reads, handed the merged row to a different
        # aligner with the wrong 5' position.)
        _idx = None
        _pos = read.reference_start
        for _i, (_op, _l) in enumerate(read.cigartuples or []):
            if _op == 3 and _pos == start:
                _idx = _i
                break
            if _op in (0, 2, 3, 7, 8):
                _pos += _l
        if _idx is not None and _repair_acceptor_overshoot(
                read, chrom_seq, _idx, start, end):
            continue
        logger.debug(
            "5' rescue reverted for %s: writer-created N-op %d-%d is "
            "non-canonical (%s..%s) and not recoverable within %d bp",
            read.query_name, start, end,
            chrom_seq[start:start + 2].upper(), chrom_seq[end - 2:end].upper(),
            _MAX_ACCEPTOR_REPAIR_BP,
        )
        read.cigartuples = pre_cigar
        read.reference_start = pre_start
        return True
    return False


# Refusal reason tokens emitted into the corrected-TSV `five_prime_rescue_refused`
# column. '' means the rescue was drawn into the BAM. Keep these stable — they are
# a user-visible audit trail, not log text.
REFUSAL_EXTEND = 'extend_refused'          # extend declined (no 5' clip, contig edge)
REFUSAL_REROUTE = 'reroute_refused'        # icp gate sent it to reroute/softclip; both declined
REFUSAL_NONCANONICAL = 'noncanonical_destination'   # the writer's own N-op was not GT-AG/GC-AG/AT-AC
# The intronic bases WERE hidden as a soft clip at the true acceptor, but no
# junction was drawn. The 5' position correction stands; only the junction claim
# is retracted. Unlike the three above this is a PARTIAL downgrade — see
# bam_processor, which must keep `five_prime_rescued` set or the writer would
# skip the soft-clip surgery this token is reporting.
REFUSAL_SOFTCLIP_ONLY = 'softclipped_no_junction'
# The surgery would have placed the read off the contig (ISSUE-015).
REFUSAL_OFF_CONTIG = 'placement_off_contig'


def _placement_within_contig(read: pysam.AlignedSegment,
                             genome: Optional[Dict[str, str]]) -> bool:
    """Does *read*'s current alignment fit inside its contig?

    `extend_read_5prime_for_junction_rescue` has carried its own off-edge
    refusal since 17d1c35 (planning/719), but that check lives INSIDE extend and
    the ISSUE-015 read never went through it: its rescue published
    `five_prime_intron_clip_pos` 5270, so the icp gate routed it to
    `reroute_intronic_tail_5prime_via_junction`, which computes its own
    `new_ref_start = five_prime_position - exon_ref_span + 1` (2990 - 8834 = -5843)
    with no bounds test. Rather than copy the guard into every helper, assert the
    invariant ONCE on the result of whichever helper ran.

    A single negative-POS record makes the whole BAM unindexable and the failure
    is silent until `samtools index` — on the 145k run that cost 1.5 h and a
    non-zero exit with an unindexed corrected BAM.
    """
    if read.reference_start is not None and read.reference_start < 0:
        return False
    end = read.reference_end
    if end is None:
        return True
    length = None
    try:
        if read.reference_name is not None:
            length = read.header.get_reference_length(read.reference_name)
    except Exception:
        length = None
    if length is None and genome is not None:
        seq, _k = get_chrom_sequence(genome, read.reference_name)
        length = len(seq) if seq else None
    return length is None or end <= length


def apply_5prime_rescue_surgery(
    read: pysam.AlignedSegment,
    correction: Dict,
    genome: Optional[Dict[str, str]] = None,
) -> Tuple[bool, str]:
    """Draw the 5' junction rescue into *read*. Returns ``(modified, refusal)``.

    THE single implementation of the 5' surgery. It used to be copy-pasted into
    ``apply_corrected_edits_to_read``, ``write_softclipped_bam`` and
    ``write_dual_bam``; the ISSUE-002 fix initially landed in only the first, so
    two writers kept fabricating junctions. One function, three call sites.

    Exactly ONE of the three helpers may run — they are alternative geometries
    for the same rescue, not stages of one.

    ISSUE-002: ``extend_read_5prime_for_junction_rescue`` derives its intron
    length from the LIVE alignment edge and never reads
    ``five_prime_intron_clip_pos``. That is correct only when the edge it will
    use IS the icp (or when the rescue published no icp at all). Otherwise it
    fabricates an N-op running to the read's OLD 5' edge: 22 of the 26 F1+F2 rows
    on the Sumner human panel, 15 of them replacing the annotated GT-AG the TSV
    names with a novel non-canonical junction.

    The gate is NOT ``icp < 0``: two panel rows (cea5e842, 7d297145) have
    icp >= 0 and reach the correct annotated junction through extend plus
    ``five_prime_upstream_trim``, a geometry reroute cannot express. Ask extend's
    own projection where its N-op would land instead.

    ``refusal`` is '' when a surgery stuck (or when the row claims no rescue),
    else one of the ``REFUSAL_*`` tokens. ``bam_processor`` consults it through
    :func:`predict_5prime_rescue_refusal` so the corrected TSV never advertises a
    rescue the corrected BAM does not carry.
    """
    # Snapshot BEFORE the surgery so _revert_selfinflicted_noncanonical_n can
    # tell an N-op this writer invented from one the aligner already had.
    _pre_nops = _n_op_intervals(read)
    _pre_cigar = read.cigartuples
    _pre_start = read.reference_start

    _icp = correction.get('five_prime_intron_clip_pos', -1)
    _exon_cig = correction.get('five_prime_exon_cigar', '')
    _rescued_flag = bool(correction.get('five_prime_rescued'))
    _rescued = _rescued_flag and correction.get('five_prime_position') is not None

    # ISSUE-026 invariant D: junction-side clip bases over exon-2 positions,
    # drawn as M after the N-op so it ends at the reported acceptor.
    _e2 = int(correction.get('five_prime_exon2_prefix', 0) or 0)
    _extend_ok = True
    if _rescued and _icp >= 0:
        _edge = projected_5prime_rescue_intron_edge(
            read,
            correction['five_prime_soft_clip'],
            correction['strand'],
            correction.get('five_prime_upstream_trim', 0),
            exon2_prefix=_e2,
        )
        _extend_ok = (_edge is not None and _edge == _icp)

    modified = False
    refusal = ''
    if _rescued and _extend_ok:
        modified = extend_read_5prime_for_junction_rescue(
            read,
            correction['five_prime_position'],
            correction['five_prime_soft_clip'],
            correction['strand'],
            exon_cigar_str=_exon_cig,
            upstream_trim=correction.get('five_prime_upstream_trim', 0),
            exon2_prefix=_e2,
        )
        if not modified:
            refusal = REFUSAL_EXTEND
    elif _rescued:
        # The N-op has to end at _icp and extend cannot express that. Reroute
        # draws it there; if reroute refuses, hide the intronic bases as a soft
        # clip at the true acceptor rather than fabricating a junction.
        if _exon_cig:
            modified = reroute_intronic_tail_5prime_via_junction(
                read,
                clip_boundary=_icp,
                five_prime_position=correction['five_prime_position'],
                exon_cigar_str=_exon_cig,
                strand=correction['strand'],
            )
        if not modified and _icp >= 0:
            modified = softclip_intronic_tail_5prime(
                read,
                clip_boundary=_icp,
                strand=correction['strand'],
            )
        if not modified:
            refusal = REFUSAL_REROUTE
    elif _rescued_flag and _icp >= 0 and not _exon_cig:
        modified = softclip_intronic_tail_5prime(
            read,
            clip_boundary=_icp,
            strand=correction['strand'],
        )

    # Hard invariant, checked on whichever helper ran (ISSUE-015): a rescue may
    # not place the read off its contig. Nothing downstream can recover from a
    # negative POS — it makes the entire BAM unindexable.
    if modified and not _placement_within_contig(read, genome):
        logger.warning(
            "5' rescue refused for %s: surgery placed the read off %s "
            "(start=%s end=%s)",
            read.query_name, read.reference_name,
            read.reference_start, read.reference_end,
        )
        read.cigartuples = _pre_cigar
        read.reference_start = _pre_start
        return False, REFUSAL_OFF_CONTIG

    if _revert_selfinflicted_noncanonical_n(
            read, genome, _pre_nops, _pre_cigar, _pre_start):
        modified = False
        refusal = REFUSAL_NONCANONICAL

    # The rescue's claim is a JUNCTION, so "did the surgery stick" is not the
    # right question — `softclip_intronic_tail_5prime` succeeds by hiding the
    # intronic bases without drawing one. Report that separately rather than
    # letting the row keep advertising a junction the BAM does not contain.
    if _rescued and not refusal and modified:
        if not any(n not in _pre_nops for n in _n_op_intervals(read)):
            refusal = REFUSAL_SOFTCLIP_ONLY

    return modified, refusal


def predict_5prime_rescue_refusal(
    read: pysam.AlignedSegment,
    correction: Dict,
    genome: Optional[Dict[str, str]] = None,
) -> str:
    """Refusal token :func:`apply_5prime_rescue_surgery` will produce for *read*,
    computed on a throwaway copy so the caller's read is untouched.

    ``bam_processor`` calls this while assembling the corrected-TSV row, because
    the TSV is written a whole pipeline stage BEFORE any BAM writer runs. Without
    it a refused surgery leaves the row still claiming `five_prime_rescued=1`
    with an exon CIGAR and a junction — and a live TSV consumer (the browser)
    draws a junction the BAM does not contain.

    This is a prediction only in the scheduling sense: it runs the SAME function
    on the SAME read after the SAME two pre-passes, so it cannot drift from what
    the writer does. It is self-consistent in the other direction too — a row
    whose fields this causes to be cleared makes the writer skip the surgery
    entirely, which is the outcome the refusal already produced.
    """
    if correction is None or not correction.get('five_prime_rescued'):
        return ''
    probe = copy.deepcopy(read)
    # Mirror apply_corrected_edits_to_read's pre-passes exactly (order matters:
    # both mutate the CIGAR the surgery then measures).
    if genome is not None:
        _decode_eq_seq_inplace(probe, genome)
    _rcl = correction.get('reanchor_clip_len', 0)
    if _rcl > 0:
        _apply_reanchor_from_clip_len(probe, _rcl)
    if genome is not None:
        realign_exon_blocks(probe, genome)
    _, refusal = apply_5prime_rescue_surgery(probe, correction, genome)
    return refusal


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
    _rcl = correction.get('reanchor_clip_len', 0)
    if _rcl > 0:
        modified |= _apply_reanchor_from_clip_len(read, _rcl)

    # Homopolymer CIGAR surgery: re-align exon blocks with X ops at homopolymer
    # positions so under-called DRS homopolymers are represented as indels.
    if genome is not None:
        modified |= realign_exon_blocks(read, genome)

    _5p_modified, _ = apply_5prime_rescue_surgery(read, correction, genome)
    modified |= _5p_modified

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
            _rcl = correction.get('reanchor_clip_len', 0)
            if _rcl > 0:
                modified |= _apply_reanchor_from_clip_len(read, _rcl)

            # Homopolymer CIGAR surgery: re-align exon blocks.
            if genome is not None:
                modified |= realign_exon_blocks(read, genome)

            # 5' junction rescue (Cat3 / Cases 1/2/2b/4) — the SAME routing and
            # canonical-destination guard the hard-clip writer uses. This block
            # used to be a copy of it and silently missed the ISSUE-002 fix.
            _5p_mod, _ = apply_5prime_rescue_surgery(read, correction, genome)
            modified |= _5p_mod

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
            _rcl = correction.get('reanchor_clip_len', 0)
            if _rcl > 0:
                shared_modified |= _apply_reanchor_from_clip_len(read, _rcl)

            # Homopolymer CIGAR surgery: re-align exon blocks.
            if genome is not None:
                shared_modified |= realign_exon_blocks(read, genome)

            # 5' rescue (Cat3 / Cases 1/2/2b/4) — identical for both BAMs, and
            # the SAME routing + canonical-destination guard the hard-clip writer
            # uses. This block used to be a copy of it and silently missed the
            # ISSUE-002 fix.
            _5p_mod, _ = apply_5prime_rescue_surgery(read, correction, genome)
            shared_modified |= _5p_mod

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
