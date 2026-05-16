#!/usr/bin/env python3
"""
Bedgraph writers driven by ``corrected_reads.tsv``.

Both functions accumulate fractional counts per (chrom, strand, corrected_3prime)
and emit a pair of strand-specific bedGraph files (``{prefix}.plus.bedgraph`` and
``{prefix}.minus.bedgraph``) via the shared netseq writer:

- ``write_netseq_assigned_bedgraph`` — Cat6 NET-seq-assigned positions only
  (rows whose ``correction_applied`` contains ``'netseq_refinement'``).
- ``write_corrected_reads_bedgraph`` — all Cat1–6 corrected 3' ends.
"""

from typing import Dict, Tuple
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def write_netseq_assigned_bedgraph(
    corrected_tsv_path: str,
    output_prefix: str,
) -> Dict[str, int]:
    """
    Write strand-specific bedGraph files for NET-seq-assigned 3' ends (Cat6).

    Reads all rows from *corrected_tsv_path* where ``correction_applied``
    contains ``'netseq_refinement'``, accumulates fractional counts per
    ``(chrom, strand, corrected_3prime)`` position, and writes two bedGraph
    files (plus and minus strand) named::

        {output_prefix}.plus.bedgraph
        {output_prefix}.minus.bedgraph

    Counts are **not** RPM-normalised (``normalize_rpm=False``) because the
    fractions already represent proportional signal and a single-sample RPM
    factor is meaningless for the subset of NET-seq-assigned reads.

    Args:
        corrected_tsv_path: Path to ``corrected_reads.tsv`` produced by
            ``rectify correct``.
        output_prefix: Path prefix for output bedGraph files (no extension).

    Returns:
        Dict with counts per strand: ``{'plus': n_positions, 'minus': n_positions}``.
    """
    from ..netseq.netseq_output import write_bedgraph as _write_bedgraph

    counts: Dict[Tuple[str, str, int], float] = {}
    n_netseq_rows = 0

    try:
        with open(corrected_tsv_path) as _f:
            hdr = _f.readline().strip().split('\t')
            try:
                i_chrom  = hdr.index('chrom')
                i_strand = hdr.index('strand')
                i_pos    = hdr.index('corrected_3prime')
                i_corr   = hdr.index('correction_applied')
            except ValueError as exc:
                raise ValueError(
                    f"Required column missing in {corrected_tsv_path}: {exc}"
                ) from exc

            i_frac = hdr.index('fraction') if 'fraction' in hdr else -1

            for line in _f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(i_chrom, i_strand, i_pos, i_corr):
                    continue
                corrections_str = parts[i_corr]
                if 'netseq_refinement' not in corrections_str:
                    continue
                try:
                    chrom  = parts[i_chrom]
                    strand = parts[i_strand]
                    pos    = int(parts[i_pos])
                    frac   = float(parts[i_frac]) if i_frac >= 0 else 1.0
                except (ValueError, IndexError):
                    continue
                key = (chrom, strand, pos)
                counts[key] = counts.get(key, 0.0) + frac
                n_netseq_rows += 1
    except OSError as exc:
        raise OSError(
            f"Cannot read corrected TSV {corrected_tsv_path}: {exc}"
        ) from exc

    if n_netseq_rows == 0:
        logger.info(
            "write_netseq_assigned_bedgraph: no netseq_refinement rows in %s; "
            "bedgraph files not written.",
            corrected_tsv_path,
        )
        return {'plus': 0, 'minus': 0}

    result_counts: Dict[str, int] = {}
    for strand_char, strand_name in [('+', 'plus'), ('-', 'minus')]:
        bg_path = Path(f"{output_prefix}.{strand_name}.bedgraph")
        _write_bedgraph(
            counts,
            bg_path,
            strand_char,
            total_reads=1,           # RPM normalisation disabled
            normalize_rpm=False,
            track_name=bg_path.stem,
        )
        n_pos = sum(1 for (_, s, _) in counts if s == strand_char)
        logger.info(
            "  Wrote %s (%d positions, %.1f fractional reads)",
            bg_path,
            n_pos,
            sum(v for (_, s, _), v in counts.items() if s == strand_char),
        )
        result_counts[strand_name] = n_pos

    return result_counts


def write_corrected_reads_bedgraph(
    corrected_tsv_path: str,
    output_prefix: str,
) -> Dict[str, int]:
    """
    Write strand-specific bedGraph files for all corrected 3' ends (Cat1–6).

    Reads every row from *corrected_tsv_path*, accumulates fractional counts
    per ``(chrom, strand, corrected_3prime)`` position using the ``fraction``
    column (defaults to 1.0 when absent), and writes two bedGraph files::

        {output_prefix}.plus.bedgraph
        {output_prefix}.minus.bedgraph

    Cat6 multi-peak rows contribute their fractional values (< 1.0) so the
    bedgraph shows proportional signal at each peak.  All other reads
    contribute 1.0.

    Counts are **not** RPM-normalised.

    Returns:
        Dict with counts per strand: ``{'plus': n_positions, 'minus': n_positions}``.
    """
    from ..netseq.netseq_output import write_bedgraph as _write_bedgraph

    counts: Dict[Tuple[str, str, int], float] = {}
    n_rows = 0

    try:
        with open(corrected_tsv_path) as _f:
            hdr = _f.readline().strip().split('\t')
            try:
                i_chrom  = hdr.index('chrom')
                i_strand = hdr.index('strand')
                i_pos    = hdr.index('corrected_3prime')
            except ValueError as exc:
                raise ValueError(
                    f"Required column missing in {corrected_tsv_path}: {exc}"
                ) from exc

            i_frac = hdr.index('fraction') if 'fraction' in hdr else -1

            for line in _f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(i_chrom, i_strand, i_pos):
                    continue
                try:
                    chrom  = parts[i_chrom]
                    strand = parts[i_strand]
                    pos    = int(parts[i_pos])
                    frac   = float(parts[i_frac]) if i_frac >= 0 else 1.0
                except (ValueError, IndexError):
                    continue
                key = (chrom, strand, pos)
                counts[key] = counts.get(key, 0.0) + frac
                n_rows += 1
    except OSError as exc:
        raise OSError(
            f"Cannot read corrected TSV {corrected_tsv_path}: {exc}"
        ) from exc

    if n_rows == 0:
        logger.info(
            "write_corrected_reads_bedgraph: no rows in %s; "
            "bedgraph files not written.",
            corrected_tsv_path,
        )
        return {'plus': 0, 'minus': 0}

    result_counts: Dict[str, int] = {}
    for strand_char, strand_name in [('+', 'plus'), ('-', 'minus')]:
        bg_path = Path(f"{output_prefix}.{strand_name}.bedgraph")
        _write_bedgraph(
            counts,
            bg_path,
            strand_char,
            total_reads=1,
            normalize_rpm=False,
            track_name=bg_path.stem,
        )
        n_pos = sum(1 for (_, s, _) in counts if s == strand_char)
        logger.info(
            "  Wrote %s (%d positions, %.1f fractional reads)",
            bg_path,
            n_pos,
            sum(v for (_, s, _), v in counts.items() if s == strand_char),
        )
        result_counts[strand_name] = n_pos

    return result_counts
