#!/usr/bin/env python3
"""
Restore poly(A) + adapter bases as soft-clips in the rectified soft-clip BAM.

After DRS pre-trimming (drs_trim_command.py) and rectify correct, the
rectified_pA_softclip.bam contains soft-clipped bases only from the
correction pipeline (e.g., residual poly-A from indel correction). The
original poly(A) tail + adapter stub that was trimmed before alignment is
absent from the BAM.

This module adds those trimmed bases BACK as soft-clip ops so the full
poly(A) + adapter is visible in IGV ("Show soft-clipped bases") and available
for downstream analysis.

Strandedness:
  Plus strand (is_reverse=False):
    RNA 5'->3' = genomic left->right. Poly(A) is at the 3' end = RIGHT of
    query_sequence. Append trimmed_3prime_seq to the right; extend or add a
    trailing S op in the CIGAR.

  Minus strand (is_reverse=True):
    BAM stores query_sequence as reverse-complement of RNA. The RNA 3' end
    (poly(A)) is at the LEFT of query_sequence (stored as poly(T)). Prepend
    reverse_complement(trimmed_3prime_seq) to the left; extend or add a
    leading S op. reference_start is NOT changed (left S ops do not consume
    reference coordinates).

Existing soft-clips from the correction pipeline are preserved:
  - Plus strand: existing trailing S is extended with the new bases.
  - Minus strand: existing leading S is extended with the new bases.
  - The opposite-end S ops (e.g., 5' soft-clip on plus strand) are untouched.

Author: Kevin R. Roy
Date: 2026-04-17
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from ..utils.genome import reverse_complement
from .drs_trim_command import load_trim_metadata

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------

_CIGAR_S = 4   # soft-clip op code
_CIGAR_H = 5   # hard-clip op code


def _extend_or_add_trailing_s(cigar: List[Tuple[int, int]], extra: int) -> List[Tuple[int, int]]:
    """Add `extra` bases to a trailing S op, or append a new one.

    Hard-clips at the far right are preserved after the S op (BAM spec allows
    H after S at the 3' end).

    Args:
        cigar:  List of (op_code, length) tuples.
        extra:  Number of additional soft-clipped bases to add.

    Returns:
        Updated CIGAR list.
    """
    if extra <= 0:
        return cigar
    cigar = list(cigar)

    # Strip trailing H (preserve for re-attachment after S)
    trailing_h = 0
    while cigar and cigar[-1][0] == _CIGAR_H:
        trailing_h += cigar.pop()[1]

    if cigar and cigar[-1][0] == _CIGAR_S:
        op, length = cigar.pop()
        cigar.append((_CIGAR_S, length + extra))
    else:
        cigar.append((_CIGAR_S, extra))

    if trailing_h:
        cigar.append((_CIGAR_H, trailing_h))

    return cigar


def _extend_or_add_leading_s(cigar: List[Tuple[int, int]], extra: int) -> List[Tuple[int, int]]:
    """Add `extra` bases to a leading S op, or prepend a new one.

    Hard-clips at the far left are preserved before the S op.

    Args:
        cigar:  List of (op_code, length) tuples.
        extra:  Number of additional soft-clipped bases to add.

    Returns:
        Updated CIGAR list.
    """
    if extra <= 0:
        return cigar
    cigar = list(cigar)

    # Strip leading H (preserve for re-attachment before S)
    leading_h = 0
    while cigar and cigar[0][0] == _CIGAR_H:
        leading_h += cigar.pop(0)[1]

    if cigar and cigar[0][0] == _CIGAR_S:
        op, length = cigar.pop(0)
        cigar.insert(0, (_CIGAR_S, length + extra))
    else:
        cigar.insert(0, (_CIGAR_S, extra))

    if leading_h:
        cigar.insert(0, (_CIGAR_H, leading_h))

    return cigar


# ---------------------------------------------------------------------------
# Per-read restoration
# ---------------------------------------------------------------------------


def _restore_read_polya(
    read: pysam.AlignedSegment,
    trimmed_3prime_seq: str,
    trimmed_3prime_quals: List[int],
) -> bool:
    """Add trimmed poly(A)+adapter back as soft-clip ops on a single read.

    Modifies `read` in-place.

    Args:
        read:                 pysam.AlignedSegment from rectified_pA_softclip.bam.
        trimmed_3prime_seq:   Poly(A)+adapter sequence in RNA 5'->3' orientation.
        trimmed_3prime_quals: Quality scores matching trimmed_3prime_seq (RNA order).

    Returns:
        True if the read was modified; False if skipped (no trimmed seq, or unmapped).
    """
    if not trimmed_3prime_seq:
        return False
    if read.is_unmapped or read.query_sequence is None:
        return False
    cigar = list(read.cigartuples or [])
    if not cigar:
        return False

    n = len(trimmed_3prime_seq)
    seq = read.query_sequence
    quals = list(read.query_qualities or [])
    has_quals = len(quals) == len(seq)

    if not read.is_reverse:
        # Plus strand: append trimmed_3prime_seq to the right
        new_seq = seq + trimmed_3prime_seq
        new_cigar = _extend_or_add_trailing_s(cigar, n)

        if has_quals:
            append_quals = trimmed_3prime_quals if len(trimmed_3prime_quals) == n else [0] * n
            new_quals = quals + append_quals
        else:
            new_quals = None

    else:
        # Minus strand: prepend RC(trimmed_3prime_seq) to the left
        # The BAM stores the read as RC of RNA; the RNA 3'-end poly(A) therefore
        # appears as poly(T) at the leftmost positions of query_sequence.
        prepend_seq = reverse_complement(trimmed_3prime_seq)
        new_seq = prepend_seq + seq
        new_cigar = _extend_or_add_leading_s(cigar, n)
        # reference_start is NOT changed — leading S ops do not consume reference

        if has_quals:
            # Quality scores for the prepend region: reverse the RNA-ordered quals
            # so they match the RC orientation stored in the BAM
            prepend_quals = list(reversed(trimmed_3prime_quals)) if len(trimmed_3prime_quals) == n else [0] * n
            new_quals = prepend_quals + quals
        else:
            new_quals = None

    read.query_sequence = new_seq
    read.cigartuples = new_cigar
    if new_quals is not None:
        read.query_qualities = pysam.qualitystring_to_array(
            ''.join(chr(min(q + 33, 126)) for q in new_quals)
        )

    return True


# ---------------------------------------------------------------------------
# Main restore function
# ---------------------------------------------------------------------------


def restore_polya_softclips(
    softclip_bam_path: str,
    metadata_path: str,
    output_bam_path: str,
    threads: int = 1,
) -> Dict:
    """Add trimmed poly(A)+adapter bases back as soft-clip ops.

    Reads rectified_pA_softclip.bam, looks up each read's trimmed poly(A)
    metadata, and appends (plus strand) or prepends (minus strand) the
    original poly(A)+adapter bases as soft-clip ops. Writes to a new BAM.

    Reads without metadata (e.g., supplementary, reads where nothing was
    trimmed) are written unchanged.

    Args:
        softclip_bam_path:  rectified_pA_softclip.bam from `rectify correct`.
        metadata_path:      Parquet or TSV from `rectify trim-polya`.
        output_bam_path:    Destination BAM (rectified_pA_softclip_full.bam).
        threads:            pysam thread count.

    Returns:
        Stats dict: total, restored, unchanged, skipped_no_meta, skipped_no_trim.
    """
    stats: Dict = {
        'total': 0,
        'restored': 0,
        'unchanged': 0,
        'skipped_no_meta': 0,
        'skipped_no_trim': 0,
    }

    print(f"Loading trim metadata from: {metadata_path}")
    metadata = load_trim_metadata(Path(metadata_path))
    print(f"  Loaded {len(metadata):,} read records")

    with pysam.AlignmentFile(softclip_bam_path, 'rb', threads=threads) as bam_in, \
         pysam.AlignmentFile(output_bam_path, 'wb', header=bam_in.header, threads=threads) as bam_out:

        for read in bam_in:
            stats['total'] += 1

            meta = metadata.get(read.query_name)
            if meta is None:
                bam_out.write(read)
                stats['skipped_no_meta'] += 1
                continue

            trimmed_seq = str(meta.get('trimmed_3prime_seq', '') or '')
            if not trimmed_seq:
                bam_out.write(read)
                stats['skipped_no_trim'] += 1
                continue

            # Parse quality scores from stored comma-separated string
            raw_quals = meta.get('trimmed_3prime_quals', '')
            if isinstance(raw_quals, list):
                trimmed_quals = [int(q) for q in raw_quals]
            elif isinstance(raw_quals, str) and raw_quals:
                try:
                    trimmed_quals = [int(q) for q in raw_quals.split(',')]
                except ValueError:
                    trimmed_quals = []
            else:
                trimmed_quals = []

            modified = _restore_read_polya(read, trimmed_seq, trimmed_quals)
            if modified:
                stats['restored'] += 1
            else:
                stats['unchanged'] += 1
            bam_out.write(read)

    logger.info(
        "restore_polya_softclips: total=%d restored=%d unchanged=%d "
        "skipped_no_meta=%d skipped_no_trim=%d",
        stats['total'], stats['restored'], stats['unchanged'],
        stats['skipped_no_meta'], stats['skipped_no_trim'],
    )
    return stats


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def run(args) -> int:
    """Entry point called from rectify/cli.py for `rectify restore-softclip`."""
    import sys

    softclip_bam = Path(args.softclip_bam)
    if not softclip_bam.exists():
        print(f"ERROR: Soft-clip BAM not found: {softclip_bam}", file=sys.stderr)
        return 1

    metadata_path = Path(args.trim_metadata)
    if not metadata_path.exists():
        # Try sibling extension
        alt = metadata_path.with_suffix('.tsv') if metadata_path.suffix == '.parquet' \
              else metadata_path.with_suffix('.parquet')
        if alt.exists():
            metadata_path = alt
        else:
            print(f"ERROR: Trim metadata not found: {metadata_path}", file=sys.stderr)
            return 1

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Soft-clip BAM: {softclip_bam}")
    print(f"Trim metadata: {metadata_path}")
    print(f"Output BAM:    {output_path}")

    stats = restore_polya_softclips(
        softclip_bam_path=str(softclip_bam),
        metadata_path=str(metadata_path),
        output_bam_path=str(output_path),
    )

    print(f"\nRestore summary:")
    print(f"  Total reads:        {stats['total']:,}")
    print(f"  Restored (polyA):   {stats['restored']:,}")
    print(f"  Unchanged:          {stats['unchanged']:,}")
    print(f"  No metadata:        {stats['skipped_no_meta']:,}")
    print(f"  No trim (polya=0):  {stats['skipped_no_trim']:,}")

    return 0
