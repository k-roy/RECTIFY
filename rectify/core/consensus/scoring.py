"""
Per-alignment scoring primitives + 5' soft-clip rescue.

Operates on a single ``AlignmentInfo`` (from ``extract.py``) plus the
backing pysam read / genome / candidate-junction context. The selection
logic that compares scored AlignmentInfos across aligners lives in
``select.py``; the streaming orchestrator lives in ``consensus.py``.

Note: ``_rescue_5prime_softclip`` mutates the passed ``AlignmentInfo``
only through its return value (no in-place writes to ``alignment``).
The ``score_alignment`` entry point writes the computed score back to
``alignment.junction_score`` as a side effect.
"""

import logging
from typing import Dict, List, Optional, Set, Tuple

import pysam

from .extract import AlignmentInfo, extract_junctions_from_cigar, get_softclip_lengths

logger = logging.getLogger(__name__)


def detect_false_3prime_junction(
    read: pysam.AlignedSegment,
    junctions: List[Tuple[int, int]],
    genome: Dict[str, str],
    min_polya_for_false: int = 3,
) -> bool:
    """
    Detect if the 3'-most junction is a false junction from poly(A) artifacts.

    False junction pattern:
    - Junction is at the 3' end of the read
    - Region before junction ends with A's
    - Region after junction starts with A's
    - Read has no substantial sequence after the junction (just poly(A))

    Returns True if the 3'-most junction appears to be a poly(A) artifact.
    """
    if not junctions:
        return False

    chrom = read.reference_name
    if chrom not in genome:
        return False

    seq = genome[chrom]
    strand = '-' if read.is_reverse else '+'

    # Get 3'-most junction (last for + strand, first for - strand)
    if strand == '+':
        junc_start, junc_end = junctions[-1]
    else:
        junc_start, junc_end = junctions[0]

    # Check if this junction is near the 3' end of the read
    read_3prime = read.reference_end if strand == '+' else read.reference_start

    if strand == '+':
        # For plus strand, junction should be near reference_end
        dist_to_3prime = read.reference_end - junc_end
    else:
        # For minus strand, junction should be near reference_start
        dist_to_3prime = junc_start - read.reference_start

    # If junction is not near 3' end, it's likely real
    if dist_to_3prime > 50:
        return False

    # Check for A-tract pattern
    # Before junction (exon side)
    if strand == '+':
        before_seq = seq[max(0, junc_start - 10):junc_start]
    else:
        before_seq = seq[junc_end:min(len(seq), junc_end + 10)]

    # Count trailing A's before junction
    trailing_a = 0
    for base in reversed(before_seq):
        if base.upper() == 'A':
            trailing_a += 1
        else:
            break

    # After junction
    if strand == '+':
        after_seq = seq[junc_end:min(len(seq), junc_end + 10)]
    else:
        after_seq = seq[max(0, junc_start - 10):junc_start]

    # Count leading A's after junction
    leading_a = 0
    for base in after_seq:
        if base.upper() == 'A':
            leading_a += 1
        else:
            break

    # Pattern for false junction: A's on both sides
    if trailing_a >= min_polya_for_false and leading_a >= min_polya_for_false:
        logger.debug(f"Detected false 3' junction at {chrom}:{junc_start}-{junc_end}")
        return True

    return False


def _rescue_5prime_softclip(
    alignment: AlignmentInfo,
    genome: Dict[str, str],
    candidate_junctions: Set[Tuple[str, int, int, str]],
    max_edit_frac: float = 0.2,
    search_window_bp: int = 300,
    rescue_seq_override: str = "",
) -> bool:
    """Check whether a 5' soft-clip (or MPB terminal error region) is explained
    by a missed upstream intron.

    Algorithm (2-pass rescue):
      Pass 1 — candidate junction pool is built by the caller from annotated
               junctions UNION junctions detected by any aligner for this read.
      Pass 2 — For each candidate junction that is UPSTREAM of the alignment's
               5' end (i.e., intron_end ≤ align_5prime for + strand, or
               intron_start ≥ reference_end for − strand), within search_window_bp:
                 1. Fetch the upstream exon-end sequence (last N bases of the exon
                    before the intron donor, where N = rescue sequence length).
                 2. Compute edit distance between rescue query bases and that
                    reference sequence.
                 3. If distance / clip_len ≤ max_edit_frac → rescue.

    Rescue sequence priority:
      rescue_seq_override (explicit caller-supplied)
        > five_prime_softclip_seq (explicit soft-clip, e.g. minimap2/gapmm2)
        > effective_five_prime_clip_seq (MPB forced-mismatch/indel terminal region)

    Coordinate conventions:
      + strand: alignment 5' end = reference_start (leftmost mapped base).
                Upstream intron has intron_end ≤ reference_start.
                Exon upstream of intron: genome[intron_start - N : intron_start].
      − strand: alignment 5' end = reference_end − 1 (rightmost mapped base).
                Upstream intron (in transcript) has intron_start ≥ reference_end.
                Exon upstream of intron (in transcript) = genome[intron_end : intron_end + N].
                BAM stores reverse-strand query_sequence in reference orientation,
                so no reverse-complement is needed for the comparison.

    Args:
        alignment: AlignmentInfo with five_prime_softclip_seq (and optionally
            effective_five_prime_clip_seq) populated
        genome: Dict mapping chrom to sequence
        candidate_junctions: Set of (chrom, intron_start, intron_end, strand)
        max_edit_frac: Max edit distance / clip_len to declare a rescue
        search_window_bp: Max distance from 5' alignment boundary to intron edge
        rescue_seq_override: Explicit sequence to use instead of alignment fields.
            When supplied, skips the field-selection logic entirely.

    Returns:
        True if the clip/error region is explained by a missed intron (no penalty),
        False if unexplained (apply normal penalty).
    """
    # Priority: explicit override > soft-clip seq > MPB terminal error seq
    clip_seq = (
        rescue_seq_override
        or alignment.five_prime_softclip_seq
        or alignment.effective_five_prime_clip_seq
    )
    if not clip_seq:
        return False

    clip_len = len(clip_seq)
    chrom = alignment.chrom

    if chrom not in genome:
        return False

    from ..spikein_filter import edit_distance
    genome_seq = genome[chrom]
    clip_seq_upper = clip_seq.upper()  # Hoist out of per-junction loop

    # Cap clip length to avoid O(clip_len²) edit_distance for very long soft-clips
    # (e.g. uLTRA reads with 200bp+ clips: 200²=40k ops/junction → 50²=2.5k ops, 16×).
    # Take the portion nearest the alignment 5' boundary (most informative for junction matching).
    _MAX_RESCUE_SEQ = 50
    if clip_len > _MAX_RESCUE_SEQ:
        if alignment.strand == '+':
            clip_seq_upper = clip_seq_upper[-_MAX_RESCUE_SEQ:]
        else:
            clip_seq_upper = clip_seq_upper[:_MAX_RESCUE_SEQ]
        clip_len = _MAX_RESCUE_SEQ

    if alignment.strand == '+':
        # 5' alignment start = leftmost mapped base
        align_5prime = alignment.reference_start
        for (j_chrom, intron_start, intron_end, j_strand) in candidate_junctions:
            if j_chrom != chrom:
                continue
            # Upstream intron: its 3'SS (intron_end) must be at or before align_5prime.
            # Directional: intron_end ≤ align_5prime (not abs, to exclude internal junctions).
            dist = align_5prime - intron_end
            if dist < 0 or dist > search_window_bp:
                continue
            # Exon1-end sequence: last clip_len bases of exon before the intron donor
            exon_end_start = intron_start - clip_len
            if exon_end_start < 0:
                continue
            exon_end_seq = genome_seq[exon_end_start:intron_start].upper()
            ed = edit_distance(clip_seq_upper, exon_end_seq)
            if ed / clip_len <= max_edit_frac:
                logger.debug(
                    f"5' rescue (+): {chrom}:{align_5prime} clip={clip_len}bp "
                    f"matches exon1 end before intron {intron_start}-{intron_end} "
                    f"(edit={ed}/{clip_len})"
                )
                return True
    else:
        # Minus strand: 5' end is at reference_end - 1 (rightmost mapped base).
        # In transcript orientation, "upstream" means higher reference coordinates.
        align_5prime = alignment.reference_end - 1
        reference_end = alignment.reference_end
        for (j_chrom, intron_start, intron_end, j_strand) in candidate_junctions:
            if j_chrom != chrom:
                continue
            # Upstream intron (in transcript): intron_start must be ≥ reference_end.
            # Directional: intron_start > align_5prime (excludes internal junctions).
            dist = intron_start - align_5prime
            if dist == 0 or dist > search_window_bp:
                continue
            # Exon upstream of intron (in transcript) = exon to the right of intron_end
            # in reference. BAM reverse-strand query_sequence is in reference orientation,
            # so compare directly without reverse-complementing.
            exon_end_seq = genome_seq[intron_end:intron_end + clip_len].upper()
            if len(exon_end_seq) < clip_len:
                continue
            ed = edit_distance(clip_seq_upper, exon_end_seq)
            if ed / clip_len <= max_edit_frac:
                logger.debug(
                    f"5' rescue (−): {chrom}:{align_5prime} clip={clip_len}bp "
                    f"matches exon upstream of intron {intron_start}-{intron_end} "
                    f"(edit={ed}/{clip_len})"
                )
                return True

    return False


def _cigar_terminal_errors(
    cigartuples,
    query_seq: str,
    ref_seq: str,
    ref_start: int,
    clip_len: int,
    scan_from_right: bool,
    scan_bp: int,
    polya_base: str = '',
) -> List[int]:
    """Build error array for the first scan_bp aligned bases near one end.

    Replaces ``read.get_aligned_pairs()`` for terminal-clip scanning.  Walking
    the CIGAR directly avoids allocating a full O(read_length) list of tuples
    and supports early termination after ``scan_bp`` positions.

    Args:
        cigartuples:     read.cigartuples
        query_seq:       read.query_sequence
        ref_seq:         genome[chrom]
        ref_start:       read.reference_start
        clip_len:        soft-clip length to skip (five_clip or three_clip)
        scan_from_right: True  → scan from high genomic coords (5' for − strand,
                                  3' for + strand).
                         False → scan from low genomic coords (5' for + strand,
                                  3' for − strand).
        scan_bp:         maximum error positions to collect before stopping.
        polya_base:      if non-empty ('A' or 'T'), insertions/mismatches of this
                         base are counted as 0 errors (poly-A filter for 3' end).
    """
    errors: List[int] = []
    if not cigartuples or not query_seq or not ref_seq:
        return errors

    ref_len = len(ref_seq)
    query_len = len(query_seq)

    if not scan_from_right:
        # Forward walk: scan left→right, skipping the first clip_len query bases.
        ref_pos = ref_start
        query_pos = 0
        for op, length in cigartuples:
            if len(errors) >= scan_bp:
                break
            if op == 5:                              # H: no bases consumed
                pass
            elif op == 4:                            # S: query only, skip
                query_pos += length
            elif op in (0, 7, 8):                    # M, =, X: both consumed
                skip_n = max(0, clip_len - query_pos)  # bases still in left clip
                for i in range(skip_n, length):
                    if len(errors) >= scan_bp:
                        break
                    rp = ref_pos + i
                    qp = query_pos + i
                    if rp < ref_len:
                        ref_b = ref_seq[rp].upper()
                        read_b = query_seq[qp].upper()
                        if polya_base and read_b == polya_base:
                            errors.append(0)
                        else:
                            errors.append(1 if (ref_b != 'N' and read_b != ref_b) else 0)
                    else:
                        errors.append(0)
                query_pos += length
                ref_pos += length
            elif op == 1:                            # I: query only (insertion)
                for i in range(length):
                    qp = query_pos + i
                    if qp >= clip_len and len(errors) < scan_bp:
                        read_b = query_seq[qp].upper()
                        if polya_base and read_b == polya_base:
                            errors.append(0)
                        else:
                            errors.append(1)
                query_pos += length
            elif op in (2, 3):                       # D, N: ref only
                ref_pos += length
    else:
        # Reversed walk: scan right→left, skipping the last clip_len query bases.
        # Pre-compute reference_end by a quick forward pass over CIGAR.
        ref_pos = ref_start
        for op, length in cigartuples:
            if op in (0, 2, 3, 7, 8):               # ops that consume reference
                ref_pos += length
        # ref_pos is now reference_end (exclusive)

        query_pos = query_len
        cutoff_qp = query_len - 1 - clip_len        # last aligned pos before right clip

        for op, length in reversed(cigartuples):
            if len(errors) >= scan_bp:
                break
            if op == 5:                              # H
                pass
            elif op == 4:                            # S: right-side soft-clip
                query_pos -= length
            elif op in (0, 7, 8):                    # M, =, X
                skip_n = max(0, query_pos - 1 - cutoff_qp)  # bases still in right clip
                for i in range(skip_n, length):
                    if len(errors) >= scan_bp:
                        break
                    qp = query_pos - 1 - i
                    rp = ref_pos - 1 - i
                    if 0 <= rp < ref_len:
                        ref_b = ref_seq[rp].upper()
                        read_b = query_seq[qp].upper()
                        if polya_base and read_b == polya_base:
                            errors.append(0)
                        else:
                            errors.append(1 if (ref_b != 'N' and read_b != ref_b) else 0)
                    else:
                        errors.append(0)
                query_pos -= length
                ref_pos -= length
            elif op == 1:                            # I: insertion
                for i in range(length):
                    qp = query_pos - 1 - i
                    if qp <= cutoff_qp and len(errors) < scan_bp:
                        read_b = query_seq[qp].upper()
                        if polya_base and read_b == polya_base:
                            errors.append(0)
                        else:
                            errors.append(1)
                query_pos -= length
            elif op in (2, 3):                       # D, N
                ref_pos -= length

    return errors


def _get_effective_5prime_clip(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    scan_bp: int = 25,
    window: int = 8,
    error_threshold: float = 0.40,
    min_errors: int = 2,
) -> int:
    """
    Compute effective 5' clip length including both explicit soft-clips and
    terminal mismatch/indel regions.

    Some aligners (notably mapPacBio) substitute forced mismatches/indels for
    soft-clips at splice junction boundaries — the read cannot be aligned there
    due to a missed or partially-resolved intron, but instead of clipping, the
    aligner records mismatches. This is structurally the same situation as a
    soft-clip and should receive the same scoring penalty.

    Algorithm:
      1. Collect per-position error flags (mismatch or insertion) for the first
         `scan_bp` aligned bases from the 5' end (after explicit soft-clips).
      2. Scan with a sliding window of size `window`. The terminal error region
         extends as long as the leading window(s) have density >= error_threshold.
         Stop at the first clean window (greedy from 5' end).
      3. Return max(explicit_soft_clip, explicit_soft_clip + terminal_error_length),
         provided at least `min_errors` errors were found in the terminal region.

    The explicit soft-clip sequence in `five_prime_softclip_seq` is NOT updated —
    the sequence-based rescue uses that field and operates on the clip sequence only.

    Args:
        read: pysam aligned read
        genome: chromosome → sequence dict (for mismatch detection without MD tags)
        scan_bp: aligned positions to scan from 5' end (after explicit soft-clips)
        window: sliding window size for density estimation
        error_threshold: fraction of errors per window to qualify as terminal region
        min_errors: minimum total errors to apply terminal clip extension

    Returns:
        Effective 5' clip length (>= explicit soft-clip).
    """
    five_clip, _ = get_softclip_lengths(read)

    if not read.query_sequence or not read.cigartuples:
        return five_clip

    chrom = read.reference_name
    if not chrom or chrom not in genome:
        return five_clip

    ref_seq = genome[chrom]

    # Build error array using CIGAR walk (avoids O(read_length) get_aligned_pairs call).
    # scan_from_right=True for minus strand (5' is at high genomic coords).
    errors = _cigar_terminal_errors(
        read.cigartuples,
        read.query_sequence,
        ref_seq,
        read.reference_start,
        clip_len=five_clip,
        scan_from_right=read.is_reverse,
        scan_bp=scan_bp,
    )

    if len(errors) < window:
        return five_clip

    # Greedy scan from 5' end: extend terminal boundary while leading windows are
    # high-error. Stop at the first clean window — terminal errors are contiguous.
    terminal_end = 0
    for i in range(len(errors) - window + 1):
        w = errors[i:i + window]
        n_err = sum(w)
        if n_err / window >= error_threshold:
            terminal_end = i + window
        else:
            break  # first clean window; stop extending

    if terminal_end == 0:
        return five_clip

    total_errors = sum(errors[:terminal_end])
    if total_errors < min_errors:
        return five_clip

    return five_clip + terminal_end


def _get_effective_3prime_clip(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    scan_bp: int = 25,
    window: int = 8,
    error_threshold: float = 0.40,
    min_errors: int = 2,
) -> int:
    """
    Compute effective 3' clip length for non-poly(A) terminal errors.

    Aligners should only soft-clip the poly(A) tail at the 3' end. Clipping or
    force-mismatching real exon sequence (non-A on + strand, non-T on - strand in
    read coordinates) means the aligner stopped before the true 3' end and should
    be penalized.

    Key distinction from 5' end: poly(A) base errors/clips are EXPECTED at the 3'
    end and are already handled by the A-tract depth penalty. Only non-A terminal
    errors are counted here.

    Algorithm:
      1. Collect per-position error flags for the first `scan_bp` aligned bases
         from the 3' end, scanning inward. Mark a position as an error only if it
         is a mismatch/insertion AND the read base is non-A (+ strand) or non-T
         (- strand in BAM coords, where poly(A) appears as poly(T)).
      2. Apply the same greedy sliding-window scan as _get_effective_5prime_clip.
      3. Return max(explicit_3prime_clip, explicit_3prime_clip + terminal_end).

    The explicit 3' soft-clip is also filtered: only non-A/T bases within the
    clip sequence count toward the penalty.

    Args:
        read: pysam aligned read
        genome: chromosome → sequence dict
        scan_bp: aligned positions to scan inward from 3' end
        window: sliding window size
        error_threshold: fraction of errors per window (non-A only)
        min_errors: minimum non-A errors to trigger penalty

    Returns:
        Effective 3' clip length (>= explicit soft-clip).
    """
    _, three_clip = get_softclip_lengths(read)

    if not read.query_sequence or not read.cigartuples:
        return three_clip

    chrom = read.reference_name
    if not chrom or chrom not in genome:
        return three_clip

    ref_seq = genome[chrom]

    # Poly(A) base in read coordinates:
    # - Plus strand: A in query = expected at 3' end (poly(A) tail)
    # - Minus strand: poly(A) tail is stored as poly(T) in the reverse-complemented BAM query
    polya_base = 'T' if read.is_reverse else 'A'

    # Build error array using CIGAR walk (avoids O(read_length) get_aligned_pairs call).
    # scan_from_right=True for plus strand (3' is at high genomic coords).
    errors = _cigar_terminal_errors(
        read.cigartuples,
        read.query_sequence,
        ref_seq,
        read.reference_start,
        clip_len=three_clip,
        scan_from_right=not read.is_reverse,
        scan_bp=scan_bp,
        polya_base=polya_base,
    )

    if len(errors) < window:
        return three_clip

    # Same greedy scan as _get_effective_5prime_clip
    terminal_end = 0
    for i in range(len(errors) - window + 1):
        w = errors[i:i + window]
        n_err = sum(w)
        if n_err / window >= error_threshold:
            terminal_end = i + window
        else:
            break

    if terminal_end == 0:
        return three_clip

    total_errors = sum(errors[:terminal_end])
    if total_errors < min_errors:
        return three_clip

    return max(three_clip, three_clip + terminal_end)


def _is_homopolymer_position(ref_seq: str, rp: int, min_run: int = 3) -> bool:
    """Return True if reference position rp is within a homopolymer run >= min_run."""
    if rp < 0 or rp >= len(ref_seq):
        return False
    base = ref_seq[rp].upper()
    if base == 'N':
        return False
    left = rp
    while left > 0 and ref_seq[left - 1].upper() == base:
        left -= 1
    right = rp + 1
    while right < len(ref_seq) and ref_seq[right].upper() == base:
        right += 1
    return (right - left) >= min_run


def _count_junction_proximity_errors(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    junction_window_bp: int = 5,
) -> float:
    """
    Count mismatches and indels within junction_window_bp of any splice junction.

    For each N operation (intron), inspect the `junction_window_bp` aligned
    bases on the exon side of both junction boundaries.  Errors (mismatches,
    insertions, deletions) in those windows indicate an aligner that placed
    indels or mismatches right at the junction rather than resolving the
    splice cleanly.

    Insertions are attributed to the preceding aligned reference position;
    deletions are attributed to their reference position directly.

    Errors at homopolymer reference positions are weighted 0.5 instead of 1.0.
    Nanopore DRS commonly undercalls homopolymer lengths; such errors are
    sequencing artifacts rather than evidence of misalignment.

    Returns:
        Total weighted error count summed across all junction-proximal windows.
    """
    if not read.cigartuples or not read.query_sequence:
        return 0

    chrom = read.reference_name
    if not chrom or chrom not in genome:
        return 0

    junctions = extract_junctions_from_cigar(read)
    if not junctions:
        return 0

    ref_seq = genome[chrom]
    query_seq = read.query_sequence

    # Build set of reference positions within junction_window_bp of any boundary.
    prox: set = set()
    for junc_start, junc_end in junctions:
        for rp in range(max(0, junc_start - junction_window_bp), junc_start):
            prox.add(rp)
        for rp in range(junc_end, min(len(ref_seq), junc_end + junction_window_bp)):
            prox.add(rp)

    ref_len = len(ref_seq)
    errors: float = 0.0
    ref_pos = read.reference_start
    query_pos = 0
    prev_rp: Optional[int] = None

    for op, length in read.cigartuples:
        if op == 5:    # H: hard-clip — no bases consumed
            pass
        elif op == 4:  # S: soft-clip — query only
            query_pos += length
        elif op in (0, 7, 8):  # M, =, X: both consumed
            for i in range(length):
                rp = ref_pos + i
                qp = query_pos + i
                if rp in prox and rp < ref_len:
                    ref_b = ref_seq[rp].upper()
                    read_b = query_seq[qp].upper()
                    if ref_b != 'N' and read_b != ref_b:
                        weight = 0.5 if _is_homopolymer_position(ref_seq, rp) else 1.0
                        errors += weight
                prev_rp = rp
            query_pos += length
            ref_pos += length
        elif op == 1:  # I: insertion — attribute to preceding ref position
            if prev_rp is not None and prev_rp in prox:
                weight = 0.5 if _is_homopolymer_position(ref_seq, prev_rp) else 1.0
                errors += weight
            query_pos += length
        elif op == 2:  # D: deletion
            for i in range(length):
                rp = ref_pos + i
                if rp in prox:
                    weight = 0.5 if _is_homopolymer_position(ref_seq, rp) else 1.0
                    errors += weight
                prev_rp = rp
            ref_pos += length
        elif op == 3:  # N: intron skip
            prev_rp = ref_pos + length - 1
            ref_pos += length

    return errors


def score_alignment(
    alignment: AlignmentInfo,
    genome: Dict[str, str],
    candidate_junctions: Optional[Set[Tuple[str, int, int, str]]] = None,
) -> float:
    """
    Score an alignment based on junction quality.

    Scoring factors:
    - 5' soft-clip penalty: -2 per unexplained clipped base; 0 if sequence-rescue
      confirms the aligner found an intron but could not align the upstream exon end
    - 3' non-poly(A) terminal error penalty: -2 per non-A/T base in trailing
      soft-clip or terminal mismatch region (indicates missed exon coverage)
    - Junction proximity errors: -1 per mismatch/indel within 5 bp of a junction

    Neither canonical splice site motifs (GT/AG) nor annotated junction
    matches are scored here, to avoid biasing against novel junctions.
    Both are used only as tiebreakers in select_best_alignment().

    Note: 3' raw endpoint proximity to the CPA is deliberately NOT scored.
    Aligners maximise their own score without knowing where the poly-A tail
    begins; find_polya_boundary() assigns the true CPA regardless of where
    the aligner stopped.  Penalising deeper 3' extension would disfavour
    aligners (e.g. mapPacBio) that correctly resolve HP insertions at the
    exon/poly-A boundary.

    Args:
        alignment: AlignmentInfo for this aligner
        genome: Dict mapping chrom to sequence
        candidate_junctions: Per-read junction pool (annotated + all aligners'
            observed junctions). Used for 5' soft-clip rescue. When provided,
            a soft-clip that matches an upstream exon end is NOT penalized.
    """
    score = 0.0

    # 5' terminal penalty: penalize both explicit soft-clips and terminal mismatch/
    # indel regions equivalently. `effective_five_prime_clip` = max(explicit_clip,
    # terminal_error_length) so aligners that force mismatches instead of clipping
    # (mapPacBio) are scored on equal footing with aligners that soft-clip (minimap2).
    #
    # Rescue uses whichever sequence is available: explicit soft-clip bytes
    # (minimap2/gapmm2/uLTRA) or MPB's forced-mismatch/indel terminal region
    # (effective_five_prime_clip_seq). Both represent the same structural
    # situation — bases that should align to the upstream exon across the intron.
    effective_clip = alignment.effective_five_prime_clip
    if effective_clip > 0:
        rescue_seq = (
            alignment.five_prime_softclip_seq
            or alignment.effective_five_prime_clip_seq
        )
        clip_rescued = (
            bool(rescue_seq)
            and candidate_junctions is not None
            and genome
            and _rescue_5prime_softclip(alignment, genome, candidate_junctions)
        )
        if not clip_rescued:
            score -= effective_clip * 2

    # NOTE: three_prime_atract_depth is intentionally NOT used for scoring.
    # Aligners do not know or care whether a 3' A-run is genomically encoded or
    # poly(A) tail — they simply maximise alignment score.  An aligner that extends
    # further into the A-rich 3'UTR (e.g. mapPacBio resolving an HP insertion to
    # reach the last non-A exon base) is MORE informative, not worse.
    # find_polya_boundary() correctly assigns the CPA regardless of where the raw
    # alignment ends — it walks back from the raw 3' end to the last non-A exon base.
    # Penalising raw 3' endpoint proximity would disfavour the most accurate aligner
    # and is the wrong criterion for consensus selection.

    # 3' non-poly(A) terminal error penalty: penalizes aligners that stop before the
    # true 3' end. Aligners should only soft-clip the poly(A) tail; non-A/T clipping or
    # forced mismatches at the 3' end indicate missed exon coverage.
    # Same -2/base scale as the 5' penalty; capped at 10 to avoid overwhelming junction scores.
    if alignment.effective_three_prime_clip > 0:
        score -= min(alignment.effective_three_prime_clip * 2, 10)

    # Junction-proximity mismatch/indel penalty: penalize aligners that place
    # errors right at splice junction boundaries.  -1 per proximal error, capped
    # at 10 to avoid overwhelming the junction/clip scores.  This favors aligners
    # (e.g. mapPacBio) that produce clean junctions with no forced errors in the
    # flanking exon sequence, relative to those (e.g. minimap2, gapmm2) that
    # accumulate mismatches/indels within a few bp of each junction.
    if alignment.junction_proximity_errors > 0:
        score -= min(alignment.junction_proximity_errors, 10)

    # NOTE: Canonical junction motifs (GT-AG) and annotated junction matches are
    # deliberately not scored here. They are used only as tiebreakers in
    # select_best_alignment() to avoid biasing against novel non-canonical junctions.

    alignment.junction_score = score
    return score
