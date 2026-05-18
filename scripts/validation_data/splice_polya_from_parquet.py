#!/usr/bin/env python3
"""Splice the poly(A) tail back onto a rectified soft-clipped BAM.

The validation bundle's ``rectified_pA_tail_soft_clipped.bam`` historically was
produced by :func:`rectify.core.commands.restore_polya_command.restore_polya_softclips`,
which appends the trimmer-determined poly(A) bases (from the
``rectify trim-polya`` parquet/TSV metadata) as a 3' soft-clip on each read.

A recent rebuild path that writes the soft-clipped BAM directly (without going
through ``restore_polya_softclips``) preserved the CIGAR rescue surgery but
dropped the poly(A) bases. This script splices them back in.

The parquet/TSV is the **authoritative** source of the poly(A) bases — NOT the
3' soft-clip of ``validation_reads_dorado_source.bam`` (which can include
adapter or non-poly(A) flanking bases that the trimmer separated out).

Strand convention (matches ``restore_polya_softclips``):
  * Plus  (``read.is_reverse == False``): append ``trimmed_3prime_seq`` to the
    END of ``query_sequence``; extend or add a TRAILING ``{N}S`` op.
  * Minus (``read.is_reverse == True``):  prepend ``reverse_complement(
    trimmed_3prime_seq)`` to the START of ``query_sequence``; extend or add a
    LEADING ``{N}S`` op. ``reference_start`` is NOT changed (leading S ops do
    not consume reference coordinates).

Usage::

    splice_polya_from_parquet.py \\
        --in-bam  rectify/data/validation/rectified/rectified_pA_tail_soft_clipped.bam \\
        --parquet path/to/polya_trim_metadata.parquet \\
        --out-bam rectify/data/validation/rectified/rectified_pA_tail_soft_clipped.with_polya.bam

The output is coordinate-sorted and indexed.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pysam

# Allow running as a standalone script *or* via -m
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from rectify.core.commands.drs_trim_command import load_trim_metadata  # noqa: E402
from rectify.utils.genome import reverse_complement  # noqa: E402

logger = logging.getLogger(__name__)

_CIGAR_S = 4  # soft-clip
_CIGAR_H = 5  # hard-clip


# ---------------------------------------------------------------------------
# CIGAR helpers (mirror restore_polya_command._extend_or_add_*_s)
# ---------------------------------------------------------------------------


def _extend_or_add_trailing_s(
    cigar: List[Tuple[int, int]], extra: int
) -> List[Tuple[int, int]]:
    if extra <= 0:
        return cigar
    cigar = list(cigar)
    trailing_h = 0
    while cigar and cigar[-1][0] == _CIGAR_H:
        trailing_h += cigar.pop()[1]
    if cigar and cigar[-1][0] == _CIGAR_S:
        _, length = cigar.pop()
        cigar.append((_CIGAR_S, length + extra))
    else:
        cigar.append((_CIGAR_S, extra))
    if trailing_h:
        cigar.append((_CIGAR_H, trailing_h))
    return cigar


def _extend_or_add_leading_s(
    cigar: List[Tuple[int, int]], extra: int
) -> List[Tuple[int, int]]:
    if extra <= 0:
        return cigar
    cigar = list(cigar)
    leading_h = 0
    while cigar and cigar[0][0] == _CIGAR_H:
        leading_h += cigar.pop(0)[1]
    if cigar and cigar[0][0] == _CIGAR_S:
        _, length = cigar.pop(0)
        cigar.insert(0, (_CIGAR_S, length + extra))
    else:
        cigar.insert(0, (_CIGAR_S, extra))
    if leading_h:
        cigar.insert(0, (_CIGAR_H, leading_h))
    return cigar


# ---------------------------------------------------------------------------
# Per-read splice
# ---------------------------------------------------------------------------


def splice_polya(
    read: pysam.AlignedSegment,
    trimmed_3prime_seq: str,
    trimmed_3prime_quals: List[int],
) -> bool:
    """Splice trimmed poly(A) bases back onto ``read`` as a 3' soft-clip.

    Mirrors ``rectify.core.commands.restore_polya_command._restore_read_polya``:
    no idempotency guard. Pre-existing soft-clips on the 3' end (from rescue
    surgery — adapter, recovered body bases, etc.) are EXTENDED by ``polya_len``
    rather than treated as "already poly-A". Invoking this script twice on the
    same BAM will double-append poly(A); it is meant for one-shot regeneration.

    Returns:
        True if the read was modified.
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
        # Plus strand: 3' end of RNA = right of query_sequence
        new_seq = seq + trimmed_3prime_seq
        new_cigar = _extend_or_add_trailing_s(cigar, n)
        if has_quals:
            append_quals = (
                trimmed_3prime_quals if len(trimmed_3prime_quals) == n else [0] * n
            )
            new_quals = quals + append_quals
        else:
            new_quals = None
    else:
        # Minus strand: BAM stores RC(RNA); RNA 3' end -> left of query_sequence
        prepend_seq = reverse_complement(trimmed_3prime_seq)
        new_seq = prepend_seq + seq
        new_cigar = _extend_or_add_leading_s(cigar, n)
        if has_quals:
            prepend_quals = (
                list(reversed(trimmed_3prime_quals))
                if len(trimmed_3prime_quals) == n
                else [0] * n
            )
            new_quals = prepend_quals + quals
        else:
            new_quals = None

    read.query_sequence = new_seq
    read.cigartuples = new_cigar
    if new_quals is not None:
        read.query_qualities = pysam.qualitystring_to_array(
            "".join(chr(min(q + 33, 126)) for q in new_quals)
        )
    return True


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def splice_polya_from_metadata(
    in_bam: Path,
    parquet: Path,
    out_bam: Path,
    threads: int = 1,
) -> Dict[str, int]:
    """Splice poly(A) tails from ``parquet`` onto reads in ``in_bam``.

    Writes a coordinate-sorted, indexed BAM at ``out_bam``.
    """
    logger.info("Loading trim metadata: %s", parquet)
    metadata = load_trim_metadata(parquet)
    logger.info("  %d trim records", len(metadata))

    stats = {
        "total": 0,
        "spliced": 0,
        "skipped_no_meta": 0,
        "skipped_no_polya": 0,
        "skipped_no_change": 0,
        "skipped_secondary_or_unmapped": 0,
    }

    unsorted_out = out_bam.with_suffix(out_bam.suffix + ".unsorted")
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    with pysam.AlignmentFile(str(in_bam), "rb", threads=threads) as bam_in, \
            pysam.AlignmentFile(
                str(unsorted_out), "wb", template=bam_in, threads=threads
            ) as bam_out:
        for read in bam_in:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                bam_out.write(read)
                stats["skipped_secondary_or_unmapped"] += 1
                continue

            stats["total"] += 1
            meta = metadata.get(read.query_name)
            if meta is None:
                bam_out.write(read)
                stats["skipped_no_meta"] += 1
                continue

            trimmed_seq = str(meta.get("trimmed_3prime_seq", "") or "")
            if not trimmed_seq:
                bam_out.write(read)
                stats["skipped_no_polya"] += 1
                continue

            raw_quals = meta.get("trimmed_3prime_quals", [])
            if isinstance(raw_quals, list):
                trimmed_quals = [int(q) for q in raw_quals]
            elif isinstance(raw_quals, str) and raw_quals:
                try:
                    trimmed_quals = [int(q) for q in raw_quals.split(",")]
                except ValueError:
                    trimmed_quals = []
            else:
                trimmed_quals = []

            modified = splice_polya(read, trimmed_seq, trimmed_quals)
            if modified:
                stats["spliced"] += 1
            else:
                stats["skipped_no_change"] += 1
            bam_out.write(read)

    # Sort + index
    pysam.sort("-@", str(threads), "-o", str(out_bam), str(unsorted_out))
    unsorted_out.unlink()
    pysam.index(str(out_bam))

    logger.info("splice summary: %s", stats)
    return stats


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--in-bam", required=True, type=Path,
                   help="Input BAM (e.g. rectified_pA_tail_soft_clipped.bam)")
    p.add_argument("--parquet", required=True, type=Path,
                   help="Parquet or TSV metadata from `rectify trim-polya`")
    p.add_argument("--out-bam", required=True, type=Path,
                   help="Output BAM (coordinate-sorted, indexed)")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if not args.in_bam.exists():
        print(f"ERROR: in-bam not found: {args.in_bam}", file=sys.stderr)
        return 1
    if not args.parquet.exists():
        # load_trim_metadata also tries sibling .tsv/.parquet
        alt = (
            args.parquet.with_suffix(".tsv")
            if args.parquet.suffix == ".parquet"
            else args.parquet.with_suffix(".parquet")
        )
        if not alt.exists():
            print(
                f"ERROR: trim metadata not found: {args.parquet} (also tried {alt})",
                file=sys.stderr,
            )
            return 1

    stats = splice_polya_from_metadata(
        in_bam=args.in_bam,
        parquet=args.parquet,
        out_bam=args.out_bam,
        threads=args.threads,
    )
    print("Splice summary:")
    for k, v in stats.items():
        print(f"  {k:32s} {v:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
