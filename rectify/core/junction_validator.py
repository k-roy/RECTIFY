#!/usr/bin/env python3
"""
Cross-sample junction validation for RECTIFY.

Three-pass COMPASS-inspired architecture:
  Pass 1 (per-sample): extract_sample_junctions()
  Pass 2 (cross-sample): filter_cross_sample_junctions()
  Pass 3 (apply back): apply_junction_filter()

Background
----------
Per-read, per-sample junction confidence scoring is insufficient.
Multi-aligner agreement on the same wrong alignment still gets tagged
"high confidence". This module implements three-pass cross-sample
validation inspired by the COMPASS pipeline:

  1. Per-sample: extract junctions with quality metrics from each BAM.
  2. Cross-sample: aggregate across all samples, apply minimum-support
     and splice-motif filters to build a validated junction set.
  3. Back to per-sample: stream through each BAM and downgrade the XC
     confidence tag for any read whose junctions are not validated.

CLI entry point (wiring is separate)
-------------------------------------
  rectify extract-junctions  --bam <bam> --sample-id <id> [--genome <fa>]
                              -o <out.tsv>

  rectify validate-junctions  --junction-tsvs <tsv1> [<tsv2> ...]
                               [--min-samples 2] [--min-total-reads 3]
                               [--max-intron 10000] [--no-canonical]
                               -o <validated.tsv>

  rectify apply-junction-filter  --input-bam <bam>
                                  --validated-junctions <validated.tsv>
                                  -o <output.bam>

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-04-01
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Dict, FrozenSet, List, Optional, Tuple

import pandas as pd

from ..utils.genome import standardize_chrom_name

try:
    import pysam
    HAS_PYSAM = True
except ImportError:  # pragma: no cover
    HAS_PYSAM = False

logger = logging.getLogger(__name__)

# CIGAR operation codes (pysam constants)
_BAM_CMATCH = 0   # M
_BAM_CREF_SKIP = 3  # N  — intron

# Canonical splice motifs
_CANONICAL_MOTIFS = frozenset({'GT-AG', 'AT-AC', 'GC-AG'})

# Confidence downgrade ladder
_DOWNGRADE = {'high': 'medium', 'medium': 'low', 'low': 'low'}


# =============================================================================
# Pass 1 — per-sample junction extraction
# =============================================================================

def extract_sample_junctions(
    bam_path: str,
    sample_id: str,
    genome: Optional[object] = None,
) -> pd.DataFrame:
    """
    Extract splice junctions from a BAM file.

    Iterates all primary (non-supplementary, non-secondary) alignments.
    For each N (intron skip) CIGAR operation the intron coordinates,
    strand, and optional splice motif are recorded.

    Parameters
    ----------
    bam_path:
        Path to coordinate-sorted, indexed BAM file.
    sample_id:
        Identifier for this sample, stored in the output DataFrame.
    genome:
        Optional open pysam.FastaFile (or compatible object with a
        ``fetch(chrom, start, end)`` method).  When provided, the 2-bp
        splice motifs at both ends of each intron are looked up and
        classified as canonical (GT-AG / AT-AC) or non-canonical.

    Returns
    -------
    pd.DataFrame with columns:
        chrom, intron_start, intron_end, strand, sample_id,
        read_count, motif, is_canonical
    """
    if not HAS_PYSAM:
        raise ImportError("pysam is required for extract_sample_junctions()")

    # Junction accumulator: (chrom, start, end, strand) -> read count
    junction_counts: Dict[Tuple[str, int, int, str], int] = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            # Skip secondary / supplementary / unmapped / no CIGAR
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_unmapped:
                continue
            if not read.cigartuples:
                continue

            strand = '-' if read.is_reverse else '+'
            chrom = standardize_chrom_name(read.reference_name) or read.reference_name
            ref_pos = read.reference_start

            for op, length in read.cigartuples:
                if op == _BAM_CREF_SKIP:
                    intron_start = ref_pos
                    intron_end = ref_pos + length
                    if length > 0:  # Skip malformed zero-length N-ops
                        junction_counts[(chrom, intron_start, intron_end, strand)] += 1
                # Advance reference position for all reference-consuming ops
                if op in (0, 2, 3, 7, 8):  # M, D, N, =, X
                    ref_pos += length

    if not junction_counts:
        return pd.DataFrame(columns=[
            'chrom', 'intron_start', 'intron_end', 'strand',
            'sample_id', 'read_count', 'motif', 'is_canonical',
        ])

    # Build parallel lists then construct DataFrame once — avoids creating a
    # list-of-dicts (one Python dict per row) which uses 5–10× more memory.
    out: Dict[str, list] = {
        'chrom': [], 'intron_start': [], 'intron_end': [], 'strand': [],
        'sample_id': [], 'read_count': [], 'motif': [], 'is_canonical': [],
    }
    for (chrom, intron_start, intron_end, strand), count in junction_counts.items():
        motif, is_canonical = _get_splice_motif(
            genome, chrom, intron_start, intron_end, strand
        )
        out['chrom'].append(chrom)
        out['intron_start'].append(intron_start)
        out['intron_end'].append(intron_end)
        out['strand'].append(strand)
        out['sample_id'].append(sample_id)
        out['read_count'].append(count)
        out['motif'].append(motif)
        out['is_canonical'].append(is_canonical)

    return pd.DataFrame(out)


def _get_splice_motif(
    genome: Optional[object],
    chrom: str,
    intron_start: int,
    intron_end: int,
    strand: str,
) -> Tuple[str, bool]:
    """
    Return (motif_string, is_canonical) for a given intron.

    The motif is 'donor2nt-acceptor2nt', e.g. 'GT-AG'.
    If no genome is provided returns ('unknown', False).

    Coordinates are 0-based half-open (pysam convention):
      donor 2-bp  : genome[intron_start : intron_start + 2]
      acceptor 2-bp: genome[intron_end - 2 : intron_end]
    """
    if genome is None:
        return ('unknown', False)

    try:
        donor2 = genome.fetch(chrom, intron_start, intron_start + 2).upper()
        acceptor2 = genome.fetch(chrom, intron_end - 2, intron_end).upper()
    except (ValueError, KeyError):
        return ('unknown', False)

    if strand == '-':
        # On the minus strand the donor is at intron_end (genomic right)
        # and acceptor is at intron_start (genomic left).
        # Reverse-complement both to get motif in transcript orientation.
        donor2 = _revcomp(acceptor2)
        acceptor2 = _revcomp(
            genome.fetch(chrom, intron_start, intron_start + 2).upper()
        )

    motif = f"{donor2}-{acceptor2}"
    is_canonical = motif in _CANONICAL_MOTIFS
    return (motif, is_canonical)


def _revcomp(seq: str) -> str:
    """Reverse-complement a DNA sequence (upper case)."""
    table = str.maketrans('ACGT', 'TGCA')
    return seq.translate(table)[::-1]


# =============================================================================
# Pass 2 — cross-sample aggregation and filtering
# =============================================================================

def filter_cross_sample_junctions(
    per_sample_dfs: List[pd.DataFrame],
    min_samples: int = 2,
    min_total_reads: int = 3,
    max_intron: int = 10_000,
    require_canonical: bool = True,
) -> FrozenSet[Tuple[str, int, int]]:
    """
    Aggregate per-sample junction DataFrames and return a validated set.

    The COMPASS approach: a junction is genuine only when it is
    independently observed in multiple samples with sufficient read
    support.  This eliminates junctions supported by a single noisy
    sample (including cases where multiple aligners agree on the same
    wrong alignment within one sample).

    Parameters
    ----------
    per_sample_dfs:
        List of DataFrames returned by ``extract_sample_junctions()``.
    min_samples:
        Minimum number of samples in which the junction must appear.
    min_total_reads:
        Minimum total read count summed across all samples.
    max_intron:
        Maximum allowed intron length in bp.  Junctions spanning more
        than this are discarded (likely mis-alignments or structural
        variants).
    require_canonical:
        When True, keep only junctions with a canonical splice motif
        (GT-AG or AT-AC).  When motif is 'unknown' (no genome provided)
        the junction is **kept** — absence of information is not evidence
        of non-canonicality.

    Returns
    -------
    frozenset of (chrom, intron_start, intron_end) tuples that passed
    all filters.  Strand is intentionally excluded so that the downstream
    filter step can apply a strand-agnostic lookup when needed.
    """
    if not per_sample_dfs:
        return frozenset()

    # Streaming aggregation — process one DataFrame at a time to avoid
    # loading all samples into memory simultaneously (OOM risk with many samples).
    # acc: (chrom, intron_start, intron_end, strand) ->
    #       {'samples': set[str], 'total_reads': int,
    #        'motif_counts': dict[str, int], 'is_canonical': bool}
    acc: Dict[Tuple[str, int, int, str], dict] = {}

    for df in per_sample_dfs:
        if df.empty:
            continue
        # Use numpy arrays instead of itertuples() — avoids materialising every
        # row as a Python namedtuple object (5–10× cheaper for large DataFrames).
        chroms = df['chrom'].to_numpy()
        starts = df['intron_start'].to_numpy()
        ends = df['intron_end'].to_numpy()
        strands = df['strand'].to_numpy()
        sids = df['sample_id'].to_numpy()
        rcounts = df['read_count'].to_numpy()
        motifs = df['motif'].to_numpy()
        canonicals = df['is_canonical'].to_numpy()

        for i in range(len(chroms)):
            key = (str(chroms[i]), int(starts[i]), int(ends[i]), str(strands[i]))
            if key not in acc:
                acc[key] = {
                    'samples': set(),
                    'total_reads': 0,
                    'motif_counts': defaultdict(int),
                    'is_canonical': False,
                }
            entry = acc[key]
            entry['samples'].add(str(sids[i]))
            entry['total_reads'] += int(rcounts[i])
            entry['motif_counts'][str(motifs[i])] += 1
            entry['is_canonical'] = entry['is_canonical'] or bool(canonicals[i])

    if not acc:
        return frozenset()

    validated_list = []
    for (chrom, intron_start, intron_end, strand), entry in acc.items():
        # Size filter
        if (intron_end - intron_start) > max_intron:
            continue
        # Multi-sample and read-depth filters
        if len(entry['samples']) < min_samples:
            continue
        if entry['total_reads'] < min_total_reads:
            continue
        # Pick most-common motif
        motif = max(entry['motif_counts'], key=entry['motif_counts'].__getitem__)
        is_canonical = entry['is_canonical']
        # Canonical motif filter — only exclude if motif is known non-canonical
        if require_canonical and not is_canonical and motif != 'unknown':
            continue
        validated_list.append((chrom, intron_start, intron_end))

    validated = frozenset(validated_list)

    logger.info(
        "Junction validation: %d junctions passed cross-sample filters "
        "(min_samples=%d, min_total_reads=%d, max_intron=%d, "
        "require_canonical=%s)",
        len(validated), min_samples, min_total_reads, max_intron,
        require_canonical,
    )
    return validated


# =============================================================================
# Pass 3 — apply filter back to per-sample BAMs
# =============================================================================

def apply_junction_filter(
    input_bam: str,
    output_bam: str,
    validated_junctions: FrozenSet[Tuple[str, int, int]],
    downgrade_tag: str = 'XC',
) -> Dict[str, int]:
    """
    Stream through a BAM and downgrade XC confidence for reads whose
    junctions are not in the validated set.

    Reads are **never discarded** — only the confidence tag is lowered.
    This preserves data for downstream inspection while signalling that
    the junction evidence is weaker than originally assigned.

    Downgrade ladder:
      'high'   → 'medium'
      'medium' → 'low'
      'low'    → 'low'  (no further downgrade)

    Parameters
    ----------
    input_bam:
        Path to input BAM (must be readable by pysam).
    output_bam:
        Path for the output BAM (will be created/overwritten).
    validated_junctions:
        Set of (chrom, intron_start, intron_end) tuples from
        ``filter_cross_sample_junctions()``.
    downgrade_tag:
        BAM tag that stores confidence level (default 'XC').

    Returns
    -------
    dict with keys: total, downgraded, unchanged
    """
    if not HAS_PYSAM:
        raise ImportError("pysam is required for apply_junction_filter()")

    counts = {'total': 0, 'downgraded': 0, 'unchanged': 0}

    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        header = bam_in.header
        with pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:
            for read in bam_in.fetch(until_eof=True):
                counts['total'] += 1

                if _read_has_invalid_junction(read, validated_junctions):
                    changed = _downgrade_read(read, downgrade_tag)
                    if changed:
                        counts['downgraded'] += 1
                    else:
                        counts['unchanged'] += 1
                else:
                    counts['unchanged'] += 1

                bam_out.write(read)

    logger.info(
        "apply_junction_filter: total=%d downgraded=%d unchanged=%d",
        counts['total'], counts['downgraded'], counts['unchanged'],
    )
    return counts


def _read_has_invalid_junction(
    read: "pysam.AlignedSegment",
    validated_junctions: FrozenSet[Tuple[str, int, int]],
) -> bool:
    """
    Return True if the read contains at least one N operation whose
    (chrom, start, end) is NOT in validated_junctions.

    Reads without any N operations are considered fully valid.
    Unmapped or cigar-less reads are also considered valid (no junctions
    to check).
    """
    if read.is_unmapped or not read.cigartuples:
        return False

    chrom = standardize_chrom_name(read.reference_name) or read.reference_name
    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op == _BAM_CREF_SKIP:
            intron_start = ref_pos
            intron_end = ref_pos + length
            key = (chrom, intron_start, intron_end)
            if key not in validated_junctions:
                return True
        if op in (0, 2, 3, 7, 8):  # M, D, N, =, X — reference-consuming
            ref_pos += length

    return False


def _downgrade_read(read: "pysam.AlignedSegment", tag: str) -> bool:
    """
    Lower the confidence tag value by one level in-place.

    If the tag is absent the read is left unchanged (no tag is set,
    to avoid introducing spurious annotations on reads that were never
    assigned a confidence level).

    Returns
    -------
    bool
        True if the tag value was actually changed, False otherwise
        (tag absent or already at lowest level).
    """
    try:
        current = read.get_tag(tag)
    except KeyError:
        return False  # Tag not present — nothing to downgrade

    new_value = _DOWNGRADE.get(current, current)
    if new_value != current:
        read.set_tag(tag, new_value)
        return True
    return False
