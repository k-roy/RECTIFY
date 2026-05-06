#!/usr/bin/env python3
"""
Restore poly(A) + adapter bases as soft-clips using the winning aligner's raw BAM.

After DRS pre-trimming (drs_trim_command.py), alignment (all aligners), per-aligner
correction, and consensus selection, this module produces corrected_polya.bam for
IGV validation. For each read:

  1. Look up the winning aligner from corrected_3ends.tsv (winning_aligner column).
  2. Fetch that read's record from the winning aligner's RAW (pre-correction) BAM.
  3. Append the original Dorado-called poly(A)+adapter sequence from the parquet
     trim metadata as a 3' soft clip.

The result shows the full read exactly as Dorado sequenced it (raw alignment +
poly-A tail), enabling direct IGV comparison against corrected.bam to see
precisely what Rectify changed.

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
        read:                 pysam.AlignedSegment from rectified_pA_tail_trimmed.bam.
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
    corrected_tsv_path: str,
    aligner_bam_paths: Dict[str, str],
    metadata_path: str,
    output_bam_path: str,
    threads: int = 1,
) -> Dict:
    """Reconstruct full Dorado reads by pulling winning aligner's raw BAM record
    and restoring the original poly(A) tail from parquet trim metadata.

    For each read in corrected_3ends.tsv:
      1. Read the winning_aligner column to find which raw BAM to pull from.
      2. Fetch that read's record from the winning aligner's raw (pre-correction) BAM.
      3. Append the Dorado-called poly(A)+adapter sequence from parquet as a 3' soft clip.

    The output (corrected_polya.bam) shows the full read as Dorado sequenced it.
    Compare against corrected.bam in IGV to see exactly what Rectify changed.

    Args:
        corrected_tsv_path:  corrected_3ends.tsv with winning_aligner column.
        aligner_bam_paths:   dict mapping aligner name → path to raw aligner BAM.
        metadata_path:       Parquet from `rectify trim-polya`.
        output_bam_path:     Destination BAM (corrected_polya.bam).
        threads:             pysam thread count.

    Returns:
        Stats dict: total, restored, unchanged, skipped_no_meta, skipped_no_trim,
                    skipped_no_aligner_bam, skipped_no_read.
    """
    import pandas as pd

    stats: Dict = {
        'total': 0,
        'restored': 0,
        'unchanged': 0,
        'skipped_no_meta': 0,
        'skipped_no_trim': 0,
        'skipped_no_aligner_bam': 0,
        'skipped_no_read': 0,
    }

    print(f"Loading trim metadata from: {metadata_path}")
    metadata = load_trim_metadata(Path(metadata_path))
    print(f"  Loaded {len(metadata):,} read records")

    print(f"Loading corrected TSV from: {corrected_tsv_path}")
    tsv = pd.read_csv(corrected_tsv_path, sep='\t', usecols=lambda c: c in ('read_id', 'winning_aligner'))
    if 'winning_aligner' not in tsv.columns:
        raise ValueError(
            f"corrected_3ends.tsv at {corrected_tsv_path} has no winning_aligner column. "
            "Re-run rectify consensus / run-all to regenerate with the updated corrected_consensus.py."
        )
    read_to_aligner: Dict[str, str] = dict(zip(tsv['read_id'], tsv['winning_aligner']))
    print(f"  {len(read_to_aligner):,} reads, {len(set(read_to_aligner.values()))} aligners")

    # Build index: aligner → set of read_ids that aligner won
    from collections import defaultdict
    aligner_to_reads: Dict[str, set] = defaultdict(set)
    for read_id, aligner in read_to_aligner.items():
        aligner_to_reads[aligner].add(read_id)

    # Use the first available raw BAM to get the BAM header for output
    header_bam_path = next(iter(aligner_bam_paths.values()))
    with pysam.AlignmentFile(str(header_bam_path), 'rb', threads=threads) as _hdr:
        header = _hdr.header.to_dict()

    with pysam.AlignmentFile(output_bam_path, 'wb', header=header, threads=threads) as bam_out:
        for aligner_name, winning_read_ids in aligner_to_reads.items():
            raw_bam_path = aligner_bam_paths.get(aligner_name)
            if not raw_bam_path or not Path(str(raw_bam_path)).exists():
                logger.warning("Raw BAM not found for aligner %s — skipping %d reads",
                               aligner_name, len(winning_read_ids))
                stats['skipped_no_aligner_bam'] += len(winning_read_ids)
                continue

            with pysam.AlignmentFile(str(raw_bam_path), 'rb', threads=threads) as bam_in:
                for read in bam_in:
                    if read.query_name not in winning_read_ids:
                        continue
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue

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
        "skipped_no_meta=%d skipped_no_trim=%d skipped_no_aligner_bam=%d",
        stats['total'], stats['restored'], stats['unchanged'],
        stats['skipped_no_meta'], stats['skipped_no_trim'], stats['skipped_no_aligner_bam'],
    )
    return stats


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def run(args) -> int:
    """Entry point called from rectify/cli.py for `rectify restore-softclip`."""
    import sys

    corrected_tsv = Path(args.corrected_tsv)
    if not corrected_tsv.exists():
        print(f"ERROR: Corrected TSV not found: {corrected_tsv}", file=sys.stderr)
        return 1

    # Parse aligner:path pairs from --aligner-bams
    aligner_bam_paths: Dict[str, str] = {}
    for token in getattr(args, 'aligner_bams', []):
        if ':' not in token:
            print(
                f"ERROR: --aligner-bams must be in 'aligner:path' format, got: {token!r}",
                file=sys.stderr,
            )
            return 1
        aligner, bam_path = token.split(':', 1)
        aligner_bam_paths[aligner] = bam_path

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

    print(f"Corrected TSV: {corrected_tsv}")
    print(f"Aligner BAMs:  {aligner_bam_paths}")
    print(f"Trim metadata: {metadata_path}")
    print(f"Output BAM:    {output_path}")

    stats = restore_polya_softclips(
        corrected_tsv_path=str(corrected_tsv),
        aligner_bam_paths=aligner_bam_paths,
        metadata_path=str(metadata_path),
        output_bam_path=str(output_path),
        threads=getattr(args, 'threads', 1),
    )

    print(f"\nRestore summary:")
    print(f"  Total reads:        {stats['total']:,}")
    print(f"  Restored (polyA):   {stats['restored']:,}")
    print(f"  Unchanged:          {stats['unchanged']:,}")
    print(f"  No metadata:        {stats['skipped_no_meta']:,}")
    print(f"  No trim (polya=0):  {stats['skipped_no_trim']:,}")

    return 0
