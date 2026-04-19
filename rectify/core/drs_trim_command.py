#!/usr/bin/env python3
"""
DRS poly(A) and adapter pre-trimming for RECTIFY.

Trims the poly(A) tail + adapter stub from Dorado-aligned DRS BAM files
BEFORE re-alignment. This gives cleaner splice-junction and 3'-end alignments
by removing the non-genomic poly(A) sequence from reads prior to minimap2.

Algorithm from docs/DRS_POLYA_ADAPTER_ANALYSIS.md:
  - Pass 0: Dorado already removed adapter; poly(A) runs to end of read
  - Pass 1: regex T[CT]{0,10}$ detects residual adapter stub
  - Pass 2: iterative peel rescues stubs with A-basecalling errors

Reads are oriented to RNA 5'->3' for detection, then stored in the unaligned
output BAM in that same orientation (for minimap2 -uf re-alignment).

Metadata written alongside the unaligned BAM records, for each read:
  - polya_len, adapter_seq, adapter_pass, pt_tag
  - trimmed_3prime_seq and trimmed_3prime_quals (RNA-oriented, for soft-clip restore)

Author: Kevin R. Roy
Date: 2026-04-17
"""

import re
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from ..utils.genome import reverse_complement

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Algorithm constants (from DRS_POLYA_ADAPTER_ANALYSIS.md)
# ---------------------------------------------------------------------------

_ADAPTER_RE = re.compile(r'T[CT]{0,10}$')
_PASS2_MAX_STUB = 15       # max stub length to try peeling in pass 2
_MIN_POLYA_PASS2 = 5       # min poly(A) length required to accept a pass-2 call
_ADAPTER_WINDOW = 150      # last N bases of RNA-oriented sequence to analyse

# ---------------------------------------------------------------------------
# Core detection functions
# ---------------------------------------------------------------------------


def _scan_polya(seq: str, max_error_rate: float = 0.0, max_consecutive_non_a: int = 1) -> int:
    """Right-to-left poly(A) scan anchored at the actual 3' end.

    Walks from seq[-1] leftward and stops when either:
      (a) the cumulative non-A rate exceeds max_error_rate, OR
      (b) consecutive non-A bases without an intervening A exceed
          max_consecutive_non_a.

    The consecutive guard is the critical correctness condition: it prevents
    consuming pure non-A flanking sequence (e.g. CTT gene-body bases upstream
    of the poly-A tail) even when max_error_rate > 0.  Without it, a cumulative
    rate of 3/20 = 15% would greedily absorb "CTT" + 17 A's because the errors
    are diluted by the large A-run to the right.

    max_error_rate alone is insufficient: it only limits the *density* of errors
    over the full scanned window, not whether any A's remain to the left of a
    non-A run.  A non-A base should only be bridged when there are genuinely
    more A's on the other side (i.e. it looks like an isolated sequencing error
    *within* the poly-A tail).  Two consecutive non-A bases are a reliable signal
    that we have left the poly-A zone and entered flanking sequence.

    The default max_error_rate=0.0 with max_consecutive_non_a=1:
      - Stops at the first non-A (rate 1/n > 0 for any n).
      - The consecutive guard is irrelevant in strict mode.
    With max_error_rate > 0 and max_consecutive_non_a=1 (recommended):
      - Bridges an isolated single-base sequencing error (e.g. AAATAAAA).
      - Stops at the first run of ≥2 consecutive non-A bases.

    Args:
        seq: RNA-oriented sequence (A-rich at 3' end for polyadenylated reads).
        max_error_rate: Tolerated fraction of non-A bases (default 0.0 = strict).
        max_consecutive_non_a: Stop when consecutive non-A count exceeds this
            value (default 1 = stop at 2+ consecutive non-A).

    Returns:
        Number of bases belonging to the poly(A) tail (0 if none detected).
    """
    n = len(seq)
    errors = 0
    total = 0
    polya_start = n
    consecutive_non_a = 0
    for i in range(n - 1, -1, -1):
        total += 1
        if seq[i] != 'A':
            errors += 1
            consecutive_non_a += 1
        else:
            consecutive_non_a = 0
        if errors / total > max_error_rate:
            break
        if consecutive_non_a > max_consecutive_non_a:
            break
        polya_start = i

    # The consecutive guard may leave polya_start pointing at a non-A base:
    # when the (max_consecutive_non_a+1)-th consecutive non-A triggers the break,
    # polya_start was already advanced to the first non-A in that run on the
    # previous iteration.  Backtrack to the first A to avoid including any
    # leading non-A bases that were speculatively accepted before the break.
    while polya_start < n and seq[polya_start] != 'A':
        polya_start += 1

    return n - polya_start


def find_polya_and_adapter(
    seq: str,
    max_error_rate: float = 0.0,
    max_consecutive_non_a: int = 1,
    min_polya: int = 1,
    adapter_window: int = _ADAPTER_WINDOW,
) -> Tuple[int, str, str, int]:
    """Detect poly(A) tail and adapter stub from a 3'-oriented read window.

    Three-pass algorithm:
      Pass 0 — no regex needed: Dorado already removed adapter; poly(A) runs
               to end of read (seq ends in A).
      Pass 1 — regex T[CT]{0,10}$ detects residual stub after poly(A).
      Pass 2 — iterative peel rescues stubs with A-basecalling errors.

    Args:
        seq: Full RNA-oriented read sequence (5'->3').
        max_error_rate: Tolerated non-A fraction for poly(A) scan.
        min_polya: Minimum poly(A) length to report (default 1).
        adapter_window: Only inspect the last N bases (default 150).

    Returns:
        (polya_len, adapter_seq, last_base, adapter_pass) where:
          polya_len    - detected poly(A) tail length in bp (0 if none)
          adapter_seq  - detected adapter stub (empty string for pass 0)
          last_base    - last base of poly(A)-trimmed sequence (quality check)
          adapter_pass - 0/1/2 (see above)
    """
    # Restrict analysis to the last adapter_window bases for efficiency
    if len(seq) > adapter_window:
        window = seq[-adapter_window:]
    else:
        window = seq

    n = len(window)
    if n == 0:
        return 0, '', '', 0

    # --- Pass 1: regex adapter detection ---
    m = _ADAPTER_RE.search(window)
    if m:
        adapter_seq = window[m.start():]
        window = window[:m.start()]
        adapter_pass = 1
    else:
        adapter_seq = ''
        adapter_pass = 0

    n = len(window)
    if n == 0:
        return 0, adapter_seq, '', adapter_pass

    last_base = window[-1]
    polya_len = _scan_polya(window, max_error_rate, max_consecutive_non_a)

    if polya_len >= min_polya:
        return polya_len, adapter_seq, last_base, adapter_pass

    # --- Pass 2: iterative peel for stubs with A-basecalling errors ---
    # Only activates when last_base is not 'A'
    if last_base != 'A':
        for k in range(1, _PASS2_MAX_STUB + 1):
            if len(window) <= k:
                break
            candidate = window[:-k]
            if not candidate or candidate[-1] != 'A':
                continue
            stub = window[-k:]
            if 'T' not in stub:
                # No adapter boundary T -> likely genuine gene body termination
                continue
            polya_len2 = _scan_polya(candidate, max_error_rate, max_consecutive_non_a)
            if polya_len2 >= _MIN_POLYA_PASS2:
                combined_adapter = (stub + adapter_seq)[:adapter_window]
                return polya_len2, combined_adapter, candidate[-1], 2

    return 0, adapter_seq, last_base, adapter_pass


# ---------------------------------------------------------------------------
# BAM record helpers
# ---------------------------------------------------------------------------


def _make_unaligned_record(
    source_read: pysam.AlignedSegment,
    trimmed_seq: str,
    trimmed_quals: Optional[List[int]],
    header: pysam.AlignmentHeader,
) -> pysam.AlignedSegment:
    """Create an unaligned BAM record from a trimmed RNA-oriented sequence.

    The returned record has:
      - flag = 4 (unmapped), is_reverse = False
      - No alignment info (reference_id = -1, etc.)
      - query_sequence in RNA 5'->3' orientation (for minimap2 -uf)
      - All original tags preserved except mapping-related ones

    Args:
        source_read:  Original aligned read (used for name and tags).
        trimmed_seq:  Trimmed sequence in RNA 5'->3' orientation.
        trimmed_quals: Trimmed base qualities matching trimmed_seq, or None.
        header:       AlignmentHeader for the new unaligned BAM.

    Returns:
        New pysam.AlignedSegment (unaligned).
    """
    new_read = pysam.AlignedSegment(header)
    new_read.query_name = source_read.query_name
    new_read.flag = 4  # unmapped
    new_read.reference_id = -1
    new_read.reference_start = 0
    new_read.mapping_quality = 0
    new_read.cigar = None
    new_read.next_reference_id = -1
    new_read.next_reference_start = 0
    new_read.template_length = 0
    new_read.query_sequence = trimmed_seq
    if trimmed_quals is not None and len(trimmed_quals) == len(trimmed_seq):
        new_read.query_qualities = pysam.qualitystring_to_array(
            ''.join(chr(min(q + 33, 126)) for q in trimmed_quals)
        )
    # Preserve the pt tag if present
    try:
        pt = source_read.get_tag('pt')
        new_read.set_tag('pt', pt)
    except KeyError:
        pass
    return new_read


def _make_unaligned_header(original_header: pysam.AlignmentHeader) -> pysam.AlignmentHeader:
    """Create a minimal header for the unaligned output BAM."""
    return pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6', 'SO': 'unsorted'}})


# ---------------------------------------------------------------------------
# Metadata I/O
# ---------------------------------------------------------------------------


def _write_metadata(rows: List[Dict], path: Path) -> None:
    """Write per-read trim metadata to parquet (preferred) or TSV (fallback)."""
    import pandas as pd

    df = pd.DataFrame(rows)

    path = Path(path)
    if path.suffix == '.parquet':
        try:
            df.to_parquet(path, index=False, engine='pyarrow')
            logger.info("Trim metadata written to %s (%d reads)", path, len(df))
            return
        except Exception as e:
            logger.warning("Cannot write parquet (%s); falling back to TSV", e)
            tsv_path = path.with_suffix('.tsv')
            df.to_csv(tsv_path, sep='\t', index=False)
            logger.info("Trim metadata written to %s (%d reads)", tsv_path, len(df))
            return

    df.to_csv(path, sep='\t', index=False)
    logger.info("Trim metadata written to %s (%d reads)", path, len(df))


def load_trim_metadata(path: Path) -> Dict[str, Dict]:
    """Load trim metadata into a read_id -> row dict.

    Accepts .parquet or .tsv path (auto-detects based on suffix; also tries
    the sibling .tsv if the .parquet file does not exist).

    Returns:
        Dict mapping read_id -> dict with trim fields.
    """
    import pandas as pd

    path = Path(path)
    if not path.exists():
        # Try sibling extension
        alt = path.with_suffix('.tsv') if path.suffix == '.parquet' else path.with_suffix('.parquet')
        if alt.exists():
            path = alt
        else:
            raise FileNotFoundError(f"Trim metadata not found: {path} (also tried {alt})")

    if path.suffix == '.parquet':
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path, sep='\t')

    # Convert trimmed_3prime_quals from string representation to list if needed
    if 'trimmed_3prime_quals' in df.columns:
        def _parse_quals(v):
            if isinstance(v, list):
                return v
            if isinstance(v, str) and v:
                try:
                    return [int(x) for x in v.split(',')]
                except ValueError:
                    return []
            return []
        df['trimmed_3prime_quals'] = df['trimmed_3prime_quals'].apply(_parse_quals)

    return {row['read_id']: row.to_dict() for _, row in df.iterrows()}


# ---------------------------------------------------------------------------
# Main trim function
# ---------------------------------------------------------------------------


def trim_drs_bam_polya(
    input_bam_path: str,
    output_bam_path: str,
    metadata_path: str,
    max_error_rate: float = 0.0,
    max_consecutive_non_a: int = 1,
    adapter_window: int = _ADAPTER_WINDOW,
    threads: int = 1,
) -> Dict:
    """Trim poly(A) tails and adapter stubs from a Dorado-aligned DRS BAM.

    Reads are oriented to RNA 5'->3' for detection. The unaligned output BAM
    stores all reads in RNA 5'->3' orientation (compatible with minimap2 -uf).

    The default max_error_rate=0.0 (strict mode) trims only unambiguously pure
    poly(A) bases. Any ambiguous bases — internal T's, runs bridging genomic
    A-tracts — are intentionally left in the trimmed read so that
    find_polya_boundary() can resolve them post-alignment using genomic context.
    Using a non-zero error rate risks over-trimming into genomic sequence, which
    would require a walk-forward correction step to undo; strict mode avoids
    this entirely by leaving all ambiguous cases for the walk-back step.

    Args:
        input_bam_path:  Dorado-aligned BAM (with optional pt:i: tags).
        output_bam_path: Destination for unaligned trimmed BAM.
        metadata_path:   Destination for per-read metadata (.parquet or .tsv).
        max_error_rate:        Tolerated non-A fraction in poly(A) scan (default 0.0).
                               0.0 = strict pure-A trim (recommended).
        max_consecutive_non_a: Stop scan after this many consecutive non-A bases
                               (default 1 = stop at ≥2 consecutive). Guards against
                               consuming upstream gene-body sequence when
                               max_error_rate > 0. Has no effect in strict mode.
        adapter_window:        Bases from 3' end to search for adapter (default 150).
        threads:         pysam thread count for BAM I/O.

    Returns:
        Stats dict: total, trimmed, untrimmed, unmapped, pass_counts {0,1,2}.
    """
    stats: Dict = {
        'total': 0,
        'trimmed': 0,
        'untrimmed': 0,
        'unmapped': 0,
        'pass_counts': {0: 0, 1: 0, 2: 0},
    }
    metadata_rows: List[Dict] = []

    input_bam_path = str(input_bam_path)
    output_bam_path = str(output_bam_path)
    metadata_path = Path(metadata_path)

    with pysam.AlignmentFile(input_bam_path, 'rb', threads=threads) as bam_in:
        unaligned_header = _make_unaligned_header(bam_in.header)

        with pysam.AlignmentFile(
            output_bam_path, 'wb', header=unaligned_header, threads=threads
        ) as bam_out:

            for read in bam_in:
                stats['total'] += 1

                seq = read.query_sequence
                if read.is_unmapped or seq is None:
                    stats['unmapped'] += 1
                    # Write a minimal unaligned record preserving the name and seq
                    if seq is not None:
                        rec = _make_unaligned_record(read, seq, list(read.query_qualities or []), unaligned_header)
                    else:
                        rec = _make_unaligned_record(read, '', None, unaligned_header)
                    bam_out.write(rec)
                    continue

                strand = '-' if read.is_reverse else '+'
                quals = list(read.query_qualities or [])

                # Orient to RNA 5'->3' for poly(A) detection
                if strand == '-':
                    rna_seq = reverse_complement(seq)
                    rna_quals = quals[::-1]
                else:
                    rna_seq = seq
                    rna_quals = quals

                # Get Dorado signal-based poly(A) estimate
                try:
                    pt_tag = int(read.get_tag('pt'))
                except (KeyError, ValueError):
                    pt_tag = -1

                # Detect poly(A) + adapter
                polya_len, adapter_seq, _last_base, adapter_pass = find_polya_and_adapter(
                    rna_seq,
                    max_error_rate=max_error_rate,
                    max_consecutive_non_a=max_consecutive_non_a,
                    adapter_window=adapter_window,
                )
                stats['pass_counts'][adapter_pass] += 1

                total_trim = polya_len + len(adapter_seq)

                # Only trim when poly(A) is actually detected (polya_len >= 1).
                # Reads with polya_len=0 are genuine non-polyadenylated transcripts
                # (~2% of DRS reads) — do not trim even if adapter stub was detected.
                if polya_len >= 1 and total_trim < len(rna_seq):
                    # Record the exact bases removed (RNA-oriented)
                    trimmed_3prime_seq = rna_seq[-total_trim:]
                    trimmed_3prime_quals = rna_quals[-total_trim:] if rna_quals else []

                    trimmed_rna_seq = rna_seq[:-total_trim]
                    trimmed_rna_quals = rna_quals[:-total_trim] if rna_quals else []
                    stats['trimmed'] += 1
                else:
                    trimmed_3prime_seq = ''
                    trimmed_3prime_quals = []
                    trimmed_rna_seq = rna_seq
                    trimmed_rna_quals = rna_quals
                    stats['untrimmed'] += 1

                # Write unaligned record in RNA 5'->3' orientation (for minimap2 -uf)
                rec = _make_unaligned_record(
                    read,
                    trimmed_rna_seq,
                    trimmed_rna_quals,
                    unaligned_header,
                )
                bam_out.write(rec)

                # Store quals as comma-separated string for portability across formats
                quals_str = ','.join(str(q) for q in trimmed_3prime_quals) if trimmed_3prime_quals else ''

                metadata_rows.append({
                    'read_id': read.query_name,
                    'strand': strand,
                    'polya_len': polya_len,
                    'adapter_seq': adapter_seq[:50],
                    'adapter_pass': adapter_pass,
                    'pt_tag': pt_tag,
                    'trimmed_3prime_seq': trimmed_3prime_seq,
                    'trimmed_3prime_quals': quals_str,
                    'original_seq_len': len(seq),
                    'trimmed_seq_len': len(trimmed_rna_seq),
                })

    _write_metadata(metadata_rows, metadata_path)

    logger.info(
        "trim_drs_bam_polya: total=%d trimmed=%d untrimmed=%d unmapped=%d "
        "pass0=%d pass1=%d pass2=%d",
        stats['total'], stats['trimmed'], stats['untrimmed'], stats['unmapped'],
        stats['pass_counts'][0], stats['pass_counts'][1], stats['pass_counts'][2],
    )
    return stats


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def run(args) -> int:
    """Entry point called from rectify/cli.py for `rectify trim-polya`."""
    import sys

    input_bam = Path(args.input_bam)
    if not input_bam.exists():
        print(f"ERROR: Input BAM not found: {input_bam}", file=sys.stderr)
        return 1

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_stem = getattr(args, 'prefix', None) or input_bam.stem
    output_bam = output_dir / f"{sample_stem}.bam"
    use_tsv = getattr(args, 'tsv', False)
    meta_ext = '.tsv' if use_tsv else '.parquet'
    metadata_path = output_dir / f"{sample_stem}_polya_trim_metadata{meta_ext}"

    print(f"Input BAM:   {input_bam}")
    print(f"Output BAM:  {output_bam}")
    print(f"Metadata:    {metadata_path}")
    print(f"max_error_rate:          {args.max_error_rate}")
    print(f"max_consecutive_non_a:   {args.max_consecutive_non_a}")
    print(f"adapter_window:          {args.adapter_window}")

    stats = trim_drs_bam_polya(
        input_bam_path=str(input_bam),
        output_bam_path=str(output_bam),
        metadata_path=str(metadata_path),
        max_error_rate=args.max_error_rate,
        max_consecutive_non_a=args.max_consecutive_non_a,
        adapter_window=args.adapter_window,
    )

    print(f"\nTrim summary:")
    print(f"  Total reads:      {stats['total']:,}")
    print(f"  Trimmed:          {stats['trimmed']:,}")
    print(f"  Untrimmed:        {stats['untrimmed']:,}")
    print(f"  Unmapped (skip):  {stats['unmapped']:,}")
    print(f"  Pass 0 (no stub): {stats['pass_counts'][0]:,}")
    print(f"  Pass 1 (regex):   {stats['pass_counts'][1]:,}")
    print(f"  Pass 2 (peel):    {stats['pass_counts'][2]:,}")

    return 0
