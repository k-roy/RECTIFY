"""Protocol-agnostic 3'-end walkback using read-vs-reference comparison.

Background
----------
Existing Module 2E ("A-tract walk-back") consults *only* the reference: it
flags an ambiguity window wherever the genome has a runs of A/T near the 3'
alignment end, but it doesn't validate that the *read* actually contains a
poly(A) extension. The result is correct ambiguity detection but no
correction — so for many protocols (QuantSeq REV in particular) ~92% of reads
emit ``correction_applied="none"`` while their alignments are silently
drifting 1–N bp into a genomic A-stretch.

The cDNA pipeline (``rectify/core/cdna_correct_command.py::walk_back_anchor_
and_tail``) already implements the correct algorithm — compare the read base
to the reference base at each aligned position, walk past stop-base (A/T)
matches, and anchor at the first non-stop match. The logic is protocol-
agnostic once a few wiring decisions are made (which side of the alignment is
the 3' end; what the "stop base" is on the read; how the gene strand maps to
``read.is_reverse``).

This module is the consolidated home for that algorithm. Each protocol-
specific entry point in ``rectify/core/correct/protocols/`` calls
:func:`walkback_3prime` with the appropriate wiring.

Author: Kevin R. Roy
Issue: ``.github/ISSUE_walkback_integration.md`` (NEW-075)
"""
from __future__ import annotations

from typing import Optional, Tuple

import pysam


THREE_PRIME_SIDE_RIGHT = "right"
THREE_PRIME_SIDE_LEFT = "left"

# Outcome string written to ``correction_applied`` column of corrected_3ends.tsv
APPLIED_WALKBACK = "polya_walkback_readgenome"
APPLIED_NONE = "none"


def walkback_3prime(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    three_prime_side: str,
    stop_base: str = "A",
) -> Tuple[int, int, str]:
    """Walk back from the 3' alignment end through a stop-base tract.

    The walk-back scans aligned (query, reference) pairs starting at the side
    of the alignment that holds the 3' end (``three_prime_side``) and moves
    inward toward the gene body. At each position the read base is compared
    to the reference base. The first position where they agree *and* the
    base is not the stop base is the canonical cleavage anchor.

    One terminal gate is applied upfront to skip the walk-back in the
    unambiguously correct case:

    1. Terminal is a non-stop base AND matches the reference → already
       anchored; no walk-back.

    All other cases (terminal is a stop-base, or terminal mismatches the
    reference) proceed to the walk-back scan. Critically, a terminal stop-base
    (A) that happens to match a genomic A is NOT skipped — non-genomically
    encoded poly-A tails frequently align over genomic A-tracts (the exact
    scenario RECTIFY was designed to correct), so the scan must walk inward to
    find the true cleavage anchor regardless of genomic A-content.

    Parameters
    ----------
    read
        ``pysam.AlignedSegment`` to correct. Must be mapped (no upstream
        check is performed).
    chrom_seq
        Reference sequence for ``read.reference_name``. 0-based, on the
        ``+`` genomic strand. The caller is responsible for selecting the
        right chromosome from the cache.
    three_prime_side
        Either :data:`THREE_PRIME_SIDE_RIGHT` (poly-A at the rightmost ref
        position, i.e. ``reference_end - 1``) or
        :data:`THREE_PRIME_SIDE_LEFT` (poly-A at ``reference_start``). The
        wrapping protocol module decides this from ``read.is_reverse`` and
        protocol chemistry.
    stop_base
        Nucleotide to walk past, on the **read** side. Default ``"A"``.
        Protocol wrappers select this from chemistry + BAM strand. The
        ``is_reverse`` → ``stop_base`` mapping is symmetric across all
        three protocols (DRS, cDNA, QuantSeq REV):

        - ``is_reverse=True`` → side=RIGHT, ``stop_base="A"`` (BAM SEQ is
          pysam RC of the basecall, so basecalled polyA T's become A's on
          the right).
        - ``is_reverse=False`` → side=LEFT, ``stop_base="T"`` (BAM SEQ is
          the basecall in alignment orientation; basecalled polyA T's stay
          at the left).

        Protocols differ only in the BAM→gene-strand mapping: DRS and cDNA
        are sense (BAM strand passes through); QuantSeq REV is antisense
        (BAM strand inverts).

    Returns
    -------
    original_3prime
        Reference position of the original 3' end before any correction.
    corrected_3prime
        Reference position of the canonical cleavage anchor. Equal to
        ``original_3prime`` when no correction was applied.
    applied
        One of :data:`APPLIED_WALKBACK` or :data:`APPLIED_NONE`.

    Notes
    -----
    *Strand handling.* This function does not know or care about gene
    strand. It returns positions in reference (``+`` genomic) coordinates
    only. Mapping ``read.is_reverse`` ↔ gene strand belongs to the
    protocol wrapper, which then passes ``three_prime_side`` accordingly.

    *Intron handling.* ``get_aligned_pairs(matches_only=False)`` is used
    with a manual ``qp is not None and rp is not None`` filter — this
    drops insertions, deletions, and ``N``-cigar (splice) gaps so the scan
    naturally crosses introns without special-case code.

    *No walk-back cap.* The aligned-pair list is bounded by read length
    (~150 bp for QuantSeq, ~longer for ONT), which serves as the natural
    cap.

    *V-primer tip artifact.* QuantSeq libraries sometimes terminate with a
    spurious non-A base (e.g. ``...AAAAAAAG``). The terminal-mismatch gate
    falls into case 3 and the scan continues through the A-run until a
    valid genomic anchor is found.
    """
    if three_prime_side not in (THREE_PRIME_SIDE_RIGHT, THREE_PRIME_SIDE_LEFT):
        raise ValueError(
            f"three_prime_side must be 'right' or 'left', got {three_prime_side!r}"
        )
    stop_base = stop_base.upper()
    if len(stop_base) != 1 or stop_base not in "ACGT":
        raise ValueError(f"stop_base must be one of A/C/G/T, got {stop_base!r}")

    if three_prime_side == THREE_PRIME_SIDE_RIGHT:
        original_3prime = read.reference_end - 1
    else:
        original_3prime = read.reference_start

    pairs = [
        (qp, rp)
        for qp, rp in read.get_aligned_pairs(matches_only=False)
        if qp is not None and rp is not None
    ]
    if not pairs:
        return original_3prime, original_3prime, APPLIED_NONE

    if three_prime_side == THREE_PRIME_SIDE_RIGHT:
        scan_pairs = list(reversed(pairs))
    else:
        scan_pairs = pairs

    qs = read.query_sequence
    if qs is None:
        return original_3prime, original_3prime, APPLIED_NONE

    # --- Terminal-position gate: skip walkback in the obvious-correct cases. ---
    first_qp, first_rp = scan_pairs[0]
    first_rb = qs[first_qp].upper()
    first_gb = chrom_seq[first_rp].upper()

    # Gate: non-stop read base matching genome → already anchored.
    if first_rb != stop_base and first_rb == first_gb:
        return original_3prime, original_3prime, APPLIED_NONE

    # Walk back until a non-stop read-vs-ref match is found.
    # This covers terminal stop-base (A) regardless of genomic A-content:
    # poly-A tails routinely align over genomic A-tracts and must be walked
    # back to the true cleavage anchor.
    corrected = original_3prime
    changed = False
    for qp, rp in scan_pairs:
        read_base = qs[qp].upper()
        ref_base = chrom_seq[rp].upper()
        if read_base == ref_base and read_base != stop_base:
            corrected = rp
            changed = corrected != original_3prime
            break

    return original_3prime, corrected, (APPLIED_WALKBACK if changed else APPLIED_NONE)


def walkback_3prime_with_qpos(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    three_prime_side: str,
    stop_base: str = "A",
) -> Tuple[int, int, str, int]:
    """Same scan as :func:`walkback_3prime`, but also returns the anchor's
    query position. Callers that need to inspect read sequence at/around the
    anchor (e.g. the cDNA pipeline counting poly-A tail length) use this
    variant; everything else should use :func:`walkback_3prime`.

    Returns ``(original_3prime, corrected_3prime, applied, anchor_qpos)``.
    ``anchor_qpos`` is the read coordinate paired with ``corrected_3prime``,
    or ``-1`` if no anchor was found (e.g. unmapped, empty pairs, missing
    sequence).
    """
    if three_prime_side not in (THREE_PRIME_SIDE_RIGHT, THREE_PRIME_SIDE_LEFT):
        raise ValueError(
            f"three_prime_side must be 'right' or 'left', got {three_prime_side!r}"
        )
    stop_base = stop_base.upper()
    if len(stop_base) != 1 or stop_base not in "ACGT":
        raise ValueError(f"stop_base must be one of A/C/G/T, got {stop_base!r}")

    if three_prime_side == THREE_PRIME_SIDE_RIGHT:
        original_3prime = read.reference_end - 1
    else:
        original_3prime = read.reference_start

    pairs = [
        (qp, rp)
        for qp, rp in read.get_aligned_pairs(matches_only=False)
        if qp is not None and rp is not None
    ]
    if not pairs:
        return original_3prime, original_3prime, APPLIED_NONE, -1

    if three_prime_side == THREE_PRIME_SIDE_RIGHT:
        scan_pairs = list(reversed(pairs))
    else:
        scan_pairs = pairs

    qs = read.query_sequence
    if qs is None:
        return original_3prime, original_3prime, APPLIED_NONE, -1

    first_qp, first_rp = scan_pairs[0]
    first_rb = qs[first_qp].upper()
    first_gb = chrom_seq[first_rp].upper()
    if first_rb != stop_base and first_rb == first_gb:
        # Terminal gate fires — already anchored at first_qp.
        return original_3prime, original_3prime, APPLIED_NONE, first_qp

    corrected = original_3prime
    anchor_qp = -1
    changed = False
    for qp, rp in scan_pairs:
        read_base = qs[qp].upper()
        ref_base = chrom_seq[rp].upper()
        if read_base == ref_base and read_base != stop_base:
            corrected = rp
            anchor_qp = qp
            changed = corrected != original_3prime
            break

    return (
        original_3prime,
        corrected,
        (APPLIED_WALKBACK if changed else APPLIED_NONE),
        anchor_qp,
    )


# =============================================================================
# Guarded walkback — DRS-grade implementation with three artifact guards
# =============================================================================
#
# This is the full-fidelity port of ``indel_corrector.find_polya_boundary``,
# parameterized by ``(side, stop_base)`` so it serves DRS (both strands), the
# cDNA pipeline, and QuantSeq REV from one implementation. The biology is the
# same in all three protocols — a stretch of read bases equal to ``stop_base``
# aligns over genomic bases of the same letter and the true 3' anchor lies
# inward at the first read=ref agreement on a non-stop base. The artifact
# guards defend against three poly-A-region noise patterns:
#
#   1. Large-deletion pre-scan (v2.9.3) — when the aligner introduces a ≥
#      ``large_del_min_bp`` deletion within ``poly_noise_window`` of the 3'
#      end to bridge a poly-A tail back to a genomic match, the backward scan
#      is restricted to before/after that deletion so the walkback can't stop
#      at a coincidental non-stop match inside the noise zone.
#
#   2. Intron-boundary guard (v2.9.4) — for ``side='left'`` walkbacks (the
#      walkback runs left-to-right toward the gene body), the scan cannot
#      cross the first N-op (intron) boundary; the terminal exon is the
#      leftmost block. Without this, all-poly-T terminal exons let the scan
#      hit a spurious match in the next exon.
#
#      For ``side='right'`` walkbacks the analogous behavior is gated by
#      the ``right_side_bridging_guard`` flag (default off). When enabled
#      it fires only if the rightmost N-op end lies within
#      ``poly_noise_window`` of the 3' end (aligner-inserted intron
#      bridging poly-A noise back to a distant exon) — in that case the
#      scan is clipped to the rightmost (terminal) exon. DRS keeps this
#      off so Cat4-style false N-ops near the 3' end are still absorbed
#      (walkback can anchor a few bp upstream in real gene-body sequence).
#      QSrev turns it on because dT-primer bridging artifacts produce a
#      synthetic tiny pre-N exon that the legacy scan would falsely
#      anchor on.
#
#   3. Tail-context false-stop (v3.0.3) — at every candidate anchor, look
#      ``tail_context_k`` positions further inward; if all of those positions
#      are read-stop-base AND at least one shows a read-vs-ref mismatch, the
#      candidate is a false stop inside the poly-A tail (a trailing T at the
#      end of a basecalled poly-A run that happens to match a genomic T).
#      Skip and continue scanning.
#
# The scan uses pre-allocated module-level scratch arrays to avoid per-call
# heap allocation (a real bottleneck on the rDNA locus where 78k long reads
# all end at the poly-A IGS). Safe under multiprocessing because each worker
# imports its own copy of this module.

import array as _array_mod

_GUARD_SCRATCH_N = 1200
_g_rp   = _array_mod.array('i', [0] * _GUARD_SCRATCH_N)  # read_pos; -1 = deletion
_g_refp = _array_mod.array('i', [0] * _GUARD_SCRATCH_N)
_g_rb   = _array_mod.array('B', [0] * _GUARD_SCRATCH_N)  # ASCII ord of read base
_g_gb   = _array_mod.array('B', [0] * _GUARD_SCRATCH_N)  # ASCII ord of genome base


def walkback_3prime_guarded(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    three_prime_side: str,
    stop_base: str = "A",
    *,
    artifact_n_ref_starts: Optional[set] = None,
    large_del_min_bp: int = 5,
    poly_noise_window: int = 50,
    tail_context_k: int = 4,
    max_scan_depth: int = 1000,
    early_exit_homopolymer_check: bool = True,
    early_exit_min_homopolymer_len: int = 4,
) -> Optional[dict]:
    """Walk back through a stop-base tract with three artifact guards.

    Same algorithm as :func:`walkback_3prime` (read-vs-reference comparison
    anchoring at the first non-stop match), but with three guards that defend
    against poly-A-region noise:

      * **Large-deletion pre-scan** — restrict the scan range so the loop
        does not stop at a coincidental match across an aligner-introduced
        ≥``large_del_min_bp``-bp deletion within ``poly_noise_window`` of
        the 3' end.
      * **Intron-boundary guard (symmetric)** — N-ops (introns) classified
        as **real** clip the scan at the boundary so it cannot anchor on
        the far side of the intron. N-ops classified as **artifact**
        (aligner-inserted bridges over poly-A noise) are crossed silently
        — the scratch array already skips them transparently. The caller
        passes the classification via ``artifact_n_ref_starts``; biology,
        not protocol identity, decides whether to cross.
      * **Tail-context false-stop** — lookahead ``tail_context_k`` positions
        inward of each candidate anchor; if all are read-stop-base AND any
        of them mismatch the genome, treat the candidate as a false stop
        and continue scanning.

    Returns a dict matching :func:`rectify.core.correct.indel_corrector.find_polya_boundary`:

        ``{'corrected_pos': int, 'original_pos': int,
           'polya_aligned_bp': int, 'correction_bp': int}``

    or ``None`` when no correction is needed (terminal already at a
    non-stop match, scan completed without an inward anchor, missing
    sequence / CIGAR, or genome out-of-range).

    Parameters
    ----------
    read, chrom_seq, three_prime_side, stop_base
        As in :func:`walkback_3prime`.
    artifact_n_ref_starts
        Set of ``ref_pos`` start positions of N-ops the caller has
        classified as poly-A-artifact bridges (typically via
        :func:`rectify.core.splice.false_junction_filter.filter_polya_artifact_junctions`).
        Artifact N-ops are crossed silently; any N-op NOT in this set is
        treated as a **real** intron and clips the scan at the boundary.

        ``None`` means "no classification was performed" — for back-compat
        with callers that haven't been migrated, every N-op is treated as
        real on ``side='left'`` and every N-op within
        ``poly_noise_window`` of the 3' end is treated as real on
        ``side='right'``. New production callers (DRS, QSrev, the
        ``find_polya_boundary`` shim) pass an explicit set.
    large_del_min_bp
        Minimum deletion length (bp) considered an aligner artifact in the
        pre-scan. Set to 0 to disable.
    poly_noise_window
        Reference-coordinate window from the 3' end. Used by the large-
        deletion pre-scan, the right-side real-N-op clipping window, and
        the back-compat ``artifact_n_ref_starts=None`` semantics on the
        right side. Default 50 bp.
    tail_context_k
        Lookahead length for the tail-context false-stop check. Set to 0 to
        disable.
    max_scan_depth
        Reference-coordinate window from the 3' end inside which positions
        are recorded into the scratch arrays. Default 1000 bp; positions
        outside this window are skipped to bound memory.

    Notes
    -----
    This is the DRS production walkback. It supersedes
    :func:`rectify.core.correct.indel_corrector.find_polya_boundary`; the two
    functions return byte-identical results on the bundled Cat1–9
    validation reads (verified in ``tests/test_walkback_readvsref.py``).
    """
    if three_prime_side not in (THREE_PRIME_SIDE_RIGHT, THREE_PRIME_SIDE_LEFT):
        raise ValueError(
            f"three_prime_side must be 'right' or 'left', got {three_prime_side!r}"
        )
    stop_base = stop_base.upper()
    if len(stop_base) != 1 or stop_base not in "ACGT":
        raise ValueError(f"stop_base must be one of A/C/G/T, got {stop_base!r}")

    cigar = read.cigartuples
    seq = read.query_sequence
    if not cigar or not seq:
        return None
    if not chrom_seq:
        return None

    stop_ord = ord(stop_base)

    # Early-exit homopolymer check (right-side / plus-strand-style only).
    #
    # For side='right' walkback, the poly-A tail aligns over a genomic
    # A-tract at the 3' end. If the genome immediately around the 3' end
    # does not contain a ``stop_base × N`` homopolymer, there's no
    # poly-A-over-A-tract artifact to correct — skip the expensive scan
    # for performance AND for conservatism (don't correct reads where the
    # artifact is not biologically plausible).
    #
    # For side='left' walkback we deliberately do NOT apply this check.
    # On the genome ``+`` strand at a minus-strand-gene CPA site, the
    # bases preceding the alignment start can be A's (mRNA strand T's) or
    # T's depending on local sequence context, so a genomic stop_base
    # homopolymer is not a reliable proxy for whether a correction is
    # needed there. (Historical note: this matches the legacy
    # ``find_polya_boundary`` behavior for minus-strand reads.)
    if (
        early_exit_homopolymer_check
        and three_prime_side == THREE_PRIME_SIDE_RIGHT
        and early_exit_min_homopolymer_len > 0
    ):
        _raw_3p = read.reference_end - 1
        _w_start = max(0, _raw_3p - 1)
        _w_end = min(
            len(chrom_seq), _raw_3p + early_exit_min_homopolymer_len + 1
        )
        _hp_window = chrom_seq[_w_start:_w_end].upper()
        if stop_base * early_exit_min_homopolymer_len not in _hp_window:
            return None

    # -------------------------------------------------------------------
    # Build the scratch arrays for the local-to-3'-end region.
    # CIGAR ops: M=0, I=1, D=2, N=3, S=4, H=5; equivalents 7 (=) and 8 (X).
    # -------------------------------------------------------------------
    ref_pos = read.reference_start
    read_pos = 0
    if cigar[0][0] == 4:  # leading soft-clip
        read_pos = cigar[0][1]

    if three_prime_side == THREE_PRIME_SIDE_RIGHT:
        original_3prime = read.reference_end - 1
        scan_start_ref = max(0, original_3prime - max_scan_depth)
        scan_end_ref = len(chrom_seq)
    else:
        original_3prime = read.reference_start
        scan_start_ref = 0
        scan_end_ref = read.reference_start + max_scan_depth

    # real_nops: list of (n_ref_start, n_ref_end, scratch_idx_just_after)
    # for every N-op the caller has NOT classified as a poly-A-artifact bridge.
    # Artifact N-ops contribute nothing — they're skipped during scratch
    # building (already true: the N-op branch never increments n) and ignored
    # during the pre-scan. LEFT clips scan_hi at the smallest
    # scratch_idx_just_after; RIGHT clips scan_lo at the largest
    # scratch_idx_just_after within poly_noise_window of the 3' end.
    real_nops: list = []
    n = 0

    for op, length in cigar:
        if op == 4:  # soft-clip
            continue
        elif op in (0, 7, 8):  # M, =, X
            for _ in range(length):
                if ref_pos < len(chrom_seq):
                    if scan_start_ref <= ref_pos <= scan_end_ref and n < _GUARD_SCRATCH_N:
                        _g_rp[n] = read_pos
                        _g_refp[n] = ref_pos
                        _g_rb[n] = ord(seq[read_pos].upper()) if read_pos < len(seq) else 78  # 'N'
                        _g_gb[n] = ord(chrom_seq[ref_pos].upper())
                        n += 1
                read_pos += 1
                ref_pos += 1
        elif op == 1:  # insertion — consumes query only
            read_pos += length
        elif op == 2:  # deletion — consumes reference only
            for _ in range(length):
                if ref_pos < len(chrom_seq):
                    if scan_start_ref <= ref_pos <= scan_end_ref and n < _GUARD_SCRATCH_N:
                        _g_rp[n] = -1
                        _g_refp[n] = ref_pos
                        _g_rb[n] = 0
                        _g_gb[n] = ord(chrom_seq[ref_pos].upper())
                        n += 1
                ref_pos += 1
        elif op == 3:  # N (intron) — consumes reference only
            n_start = ref_pos
            n_end = ref_pos + length
            is_artifact = (
                artifact_n_ref_starts is not None
                and n_start in artifact_n_ref_starts
            )
            if not is_artifact:
                # scratch_idx_just_after = current n (next M base lands here)
                real_nops.append((n_start, n_end, n))
            ref_pos += length

    if n == 0:
        return None

    # -------------------------------------------------------------------
    # Scan range determination per side, with the three guards.
    # -------------------------------------------------------------------
    if three_prime_side == THREE_PRIME_SIDE_RIGHT:
        scan_lo = 0
        scan_hi = n  # exclusive upper bound; main loop iterates [scan_lo, scan_hi) in reverse

        # Large-deletion pre-scan: identify any large deletion within
        # poly_noise_window of the 3' end. The check fires only when the
        # alignment ends in a stop_base context — otherwise the pre-scan
        # would shrink the window in cases where there is no poly-A artifact.
        last_gb_chr: Optional[str] = None
        for _k in range(n - 1, -1, -1):
            if _g_rp[_k] != -1:
                last_gb_chr = chr(_g_gb[_k])
                break
        if large_del_min_bp > 0 and last_gb_chr == stop_base:
            _i = n - 1
            while _i >= 0:
                refp_i = _g_refp[_i]
                rp_i = _g_rp[_i]
                if original_3prime - refp_i > poly_noise_window:
                    break
                if rp_i == -1:  # deletion
                    _j = _i
                    while _j > 0 and _g_rp[_j - 1] == -1:
                        _j -= 1
                    if _i - _j + 1 >= large_del_min_bp:
                        scan_hi = _j
                    _i = _j - 1
                else:
                    _i -= 1

        # Real-N-op clipping: every N-op the caller classified as real is
        # an exon boundary the scan cannot cross. Use the LARGEST
        # scratch_idx_just_after among real N-ops whose end lies within
        # poly_noise_window of the 3' end — that is the last real intron
        # exit on the rightmost exon's left side. Anything beyond that is
        # an upstream exon and not eligible to anchor a 3' end.
        for _n_start, _n_end, _post_idx in real_nops:
            if original_3prime - _n_end <= poly_noise_window and _post_idx > scan_lo:
                scan_lo = _post_idx

        # Terminal gate: check the rightmost real position before any walk.
        # find_polya_boundary doesn't have this gate (its loop simply stops
        # at the first non-stop match, which == the same outcome), but we
        # mirror it here for symmetry with walkback_3prime.
        # (Implicitly handled by the loop below — no special case needed.)

        # Main scan: walk right-to-left through scratch arrays in [scan_lo, scan_hi).
        true_cpa: Optional[int] = None
        for i in range(scan_hi - 1, scan_lo - 1, -1):
            if _g_rp[i] == -1:
                continue
            refp = _g_refp[i]
            rb = _g_rb[i]
            gb = _g_gb[i]
            if rb == gb and gb != stop_ord:
                # Tail-context false-stop guard.
                if tail_context_k > 0:
                    ctx_all_stop = True
                    ctx_has_mismatch = False
                    ctx_n = 0
                    for _j in range(i - 1, scan_lo - 1, -1):
                        if _g_rp[_j] == -1:
                            continue
                        jrb = _g_rb[_j]
                        jgb = _g_gb[_j]
                        if jrb != stop_ord:
                            ctx_all_stop = False
                            break
                        if jgb != jrb:
                            ctx_has_mismatch = True
                        ctx_n += 1
                        if ctx_n >= tail_context_k:
                            break
                    if ctx_n >= tail_context_k and ctx_all_stop and ctx_has_mismatch:
                        continue  # false stop — keep scanning
                true_cpa = refp
                break

        if true_cpa is None:
            return None
        polya_len = original_3prime - true_cpa
        if polya_len <= 0:
            return None
        return {
            'corrected_pos': true_cpa,
            'original_pos': original_3prime,
            'polya_aligned_bp': polya_len,
            'correction_bp': polya_len,
        }

    # ---- side == 'left' (minus-strand DRS, cDNA orient=rev, QSrev plus-gene) ----
    scan_lo = 0
    scan_hi = n  # main loop iterates [scan_lo, scan_hi) forward

    # Real-N-op clipping: clip scan_hi to the smallest scratch_idx_just_after
    # among real N-ops — i.e., the first real-intron exit on the leftmost
    # (terminal-3') exon's right side. The scan cannot cross the boundary.
    # first_real_n is also retained for the fallback snap below.
    first_real_n: Optional[tuple] = None
    if real_nops:
        first_real_n = real_nops[0]  # CIGAR order = ref-left-to-right
        for _n_start, _n_end, _post_idx in real_nops:
            if _post_idx < scan_hi:
                scan_hi = _post_idx

    # Find the first real (non-deletion) genome base in the scan window.
    first_gb_chr: Optional[str] = None
    for _k in range(scan_lo, scan_hi):
        if _g_rp[_k] != -1:
            first_gb_chr = chr(_g_gb[_k])
            break

    # For left-side walkback, find_polya_boundary returns None when the
    # alignment doesn't start in a stop_base context — there's no poly-A
    # over-alignment to correct. Preserve that behavior.
    if first_gb_chr != stop_base:
        return None

    # Large-deletion pre-scan (forward variant).
    if large_del_min_bp > 0:
        _i = scan_lo
        while _i < scan_hi:
            rp_i = _g_rp[_i]
            refp_i = _g_refp[_i]
            if refp_i - original_3prime > poly_noise_window:
                break
            if rp_i == -1:
                _j = _i
                while _j < scan_hi - 1 and _g_rp[_j + 1] == -1:
                    _j += 1
                if _j - _i + 1 >= large_del_min_bp:
                    scan_lo = _j + 1
                _i = _j + 1
            else:
                _i += 1

    # Main scan: walk left-to-right through scratch arrays.
    true_cpa: Optional[int] = None
    for i in range(scan_lo, scan_hi):
        if _g_rp[i] == -1:
            continue
        refp = _g_refp[i]
        rb = _g_rb[i]
        gb = _g_gb[i]
        if rb == gb and gb != stop_ord:
            if tail_context_k > 0:
                ctx_all_stop = True
                ctx_has_mismatch = False
                ctx_n = 0
                for _j in range(i + 1, scan_hi):
                    if _g_rp[_j] == -1:
                        continue
                    jrb = _g_rb[_j]
                    jgb = _g_gb[_j]
                    if jrb != stop_ord:
                        ctx_all_stop = False
                        break
                    if jgb != jrb:
                        ctx_has_mismatch = True
                    ctx_n += 1
                    if ctx_n >= tail_context_k:
                        break
                if ctx_n >= tail_context_k and ctx_all_stop and ctx_has_mismatch:
                    continue
            true_cpa = refp
            break

    # Fallback: scan found nothing AND scan_hi was clipped at a real N-op AND
    # the last pre-N base is still in the stop-base zone → use the N-op
    # start as the CPA anchor (the intron is the natural exon end).
    if (
        true_cpa is None
        and first_real_n is not None
        and scan_hi > 0
    ):
        last_pre_n_gb: Optional[str] = None
        for _k in range(scan_hi - 1, -1, -1):
            if _g_rp[_k] != -1:
                last_pre_n_gb = chr(_g_gb[_k])
                break
        if last_pre_n_gb == stop_base:
            true_cpa = first_real_n[0]  # n_ref_start

    if true_cpa is None:
        return None
    polya_len = true_cpa - original_3prime
    if polya_len <= 0:
        return None
    return {
        'corrected_pos': true_cpa,
        'original_pos': original_3prime,
        'polya_aligned_bp': polya_len,
        'correction_bp': polya_len,
    }


def is_minus_strand_dRNA(read: pysam.AlignedSegment) -> bool:
    """Return ``True`` if a direct-RNA read is on the genome ``-`` strand.

    For ONT direct-RNA (dRNA-seq) the BAM strand passes through unchanged:
    ``is_reverse=False`` is a ``+`` strand RNA, ``is_reverse=True`` is a
    ``-`` strand RNA.
    """
    return read.is_reverse


def walkback_drs(
    read: pysam.AlignedSegment, chrom_seq: str
) -> Tuple[int, int, str, str]:
    """Walk-back wrapper for ONT direct-RNA sequencing (dRNA-seq).

    Strand convention: BAM strand == gene strand.

    Pysam returns ``query_sequence`` in alignment orientation, so an
    ``is_reverse=True`` (minus-strand) read has its sequence
    reverse-complemented relative to the basecall. The RNA poly(A) tail
    (A's at the basecall 3' end) therefore appears as **T's at the LEFT
    of query_sequence** in minus-strand reads. Symmetrically, plus-strand
    reads keep A's at the right.

    Stop base:
      - Plus strand  (``is_reverse=False``): walk past read A's on the right
      - Minus strand (``is_reverse=True``):  walk past read T's on the left

    Returns ``(original_3prime, corrected_3prime, applied, gene_strand)``.

    This uses the simple :func:`walkback_3prime` core (no artifact guards) —
    suitable for unit tests on synthetic reads. The DRS production path in
    :mod:`rectify.core.bam.bam_processor` calls :func:`walkback_drs_full`, which
    delegates to :func:`walkback_3prime_guarded` with the same three
    artifact guards that the legacy ``indel_corrector.find_polya_boundary``
    used to apply (now a thin alias for :func:`walkback_3prime_guarded`).
    """
    if is_minus_strand_dRNA(read):
        gene_strand = "-"
        side = THREE_PRIME_SIDE_LEFT
        stop_base = "T"
    else:
        gene_strand = "+"
        side = THREE_PRIME_SIDE_RIGHT
        stop_base = "A"
    orig, corr, applied = walkback_3prime(read, chrom_seq, side, stop_base=stop_base)
    return orig, corr, applied, gene_strand


def _classify_artifact_nops(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    gene_strand: str,
    artifact_analyses=None,
) -> set:
    """Return the set of N-op ``ref_pos`` starts the caller treats as
    poly-A artifact bridges (to be crossed silently during walkback).

    If ``artifact_analyses`` is provided (caller already ran the filter
    upstream — e.g., ``bam_processor.correct_read_3prime``), extract
    junction starts from it. Otherwise run
    :func:`rectify.core.splice.false_junction_filter.filter_polya_artifact_junctions`
    here.
    """
    if artifact_analyses is not None:
        return {a.junction_start for a in artifact_analyses}
    from rectify.core.splice.false_junction_filter import (
        filter_polya_artifact_junctions,
        _GenomeDictReference,
    )
    _, analyses = filter_polya_artifact_junctions(
        read, _GenomeDictReference({read.reference_name: chrom_seq}), gene_strand
    )
    return {a.junction_start for a in analyses}


def walkback_drs_full(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    *,
    artifact_analyses=None,
    large_del_min_bp: int = 5,
    poly_noise_window: int = 50,
    tail_context_k: int = 4,
    max_scan_depth: int = 1000,
) -> Optional[dict]:
    """DRS production walkback — :func:`walkback_3prime_guarded` with the
    DRS strand→side/stop_base mapping pre-wired and poly-A-artifact N-op
    classification routed in.

    Returns the same dict shape as
    :func:`rectify.core.correct.indel_corrector.find_polya_boundary`:

        ``{'corrected_pos', 'original_pos', 'polya_aligned_bp', 'correction_bp'}``

    or ``None`` when no correction is needed. This is the drop-in
    replacement for the legacy ``find_polya_boundary`` in
    :mod:`rectify.core.bam.bam_processor`.

    Pass ``artifact_analyses`` when the caller has already run
    :func:`rectify.core.splice.false_junction_filter.filter_polya_artifact_junctions`
    (e.g. ``bam_processor.correct_read_3prime`` runs it for the 3'SS-rescue
    path); otherwise this function runs it itself.
    """
    if is_minus_strand_dRNA(read):
        side = THREE_PRIME_SIDE_LEFT
        stop_base = "T"
        gene_strand = "-"
    else:
        side = THREE_PRIME_SIDE_RIGHT
        stop_base = "A"
        gene_strand = "+"
    artifact_starts = _classify_artifact_nops(
        read, chrom_seq, gene_strand, artifact_analyses
    )
    return walkback_3prime_guarded(
        read,
        chrom_seq,
        side,
        stop_base=stop_base,
        artifact_n_ref_starts=artifact_starts,
        large_del_min_bp=large_del_min_bp,
        poly_noise_window=poly_noise_window,
        tail_context_k=tail_context_k,
        max_scan_depth=max_scan_depth,
    )
