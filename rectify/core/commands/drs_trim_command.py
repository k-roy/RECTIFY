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
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from ...utils.genome import reverse_complement
from ..splice.repeat_expansion import (
    dominant_repeat_period,
    DEFAULT_K_RANGE as _REPEAT_K_RANGE,
    _MIN_REPEATS,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Algorithm constants (from DRS_POLYA_ADAPTER_ANALYSIS.md)
# ---------------------------------------------------------------------------

# DNA adapter stub pattern (Pass 1) — also handles RNA notation (U instead of T)
_ADAPTER_RE = re.compile(r'[TU][CTU]{0,10}$')
_PASS2_MAX_STUB = 15       # max stub length to try peeling in pass 2
_MIN_POLYA_PASS2 = 5       # min poly(A) length required to accept a pass-2 call
_ADAPTER_WINDOW = 150      # last N bases of RNA-oriented sequence to analyse

# Terminal (AAG/GAA)n triplet-repeat basecaller-artifact strip (Dorado RNA004).
# The poly-A homopolymer is mis-basecalled as a low-period multi-base repeat that
# the poly-A scan cannot see (no terminal A-run). We peel that block before the
# poly-A scan so the artifact is not carried into alignment. Strict no-op on clean
# poly-A: the block must contain periodic non-A bases (see _find_terminal_repeat_block).
# min_len=15 chosen from the CNTL_21.8 block-length distribution (2026-06-14): real
# artifacts are long (median 59 bp, p5=17), only ~2.9% of fired reads fall in 9-14 bp.
# 15 keeps ~97% of true artifacts while sparing genuine short GAA/AAG genomic tracts —
# important in the SMA (repeat-disease) patient context where short expansions are biology.
_REPEAT_MIN_LEN = 15       # >= 5 tandem copies of a triplet
_REPEAT_MIN_FRAC = 0.8     # canonical-motif purity required to peel

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


def _find_terminal_repeat_block(
    window: str,
    min_len: int = _REPEAT_MIN_LEN,
    min_frac: float = _REPEAT_MIN_FRAC,
    k_range: Tuple[int, ...] = _REPEAT_K_RANGE,
) -> Tuple[int, str]:
    """Locate a terminal low-period MULTI-BASE repeat block (the (AAG/GAA)n Dorado
    RNA004 basecaller artifact) at the 3' end of an RNA-oriented ``window``.

    Returns ``(block_len, canonical_motif)`` — the number of 3'-terminal bases to
    peel and the dominant repeat motif — or ``(0, '')`` when no block qualifies.

    The block's 5' boundary is anchored on the leftmost *on-phase non-A* base, so a
    pure poly-A run immediately 5' of the artifact is **preserved** for the
    downstream poly-A scan (the 4.8% "blocked poly-A" recovery case). This anchoring
    is the key correctness condition: poly-A satisfies *every* period trivially, so a
    naive period walk would merge the genuine poly-A run into the artifact block and
    destroy poly-A recovery.

    Strict no-op on clean poly-A/poly-T: the walk is anchored on a non-A base and the
    final gate rejects single-base motifs (``len(set(motif)) <= 1`` on
    :func:`dominant_repeat_period`'s result — the same multi-base rule as
    :func:`is_repeat_expansion`), so a pure homopolymer tail never qualifies. Off-phase
    basecall errors are tolerated (break only at 2 consecutive, mirroring
    ``_scan_polya``); the purity gate then bounds how far the block may absorb
    non-periodic flanking sequence.
    """
    n = len(window)
    if n < min_len:
        return 0, ''

    # Anchor phase on the rightmost non-A base. The artifact carries the periodic
    # non-A signal; trailing A's of the final repeat unit (e.g. the "AA" of a
    # terminal ...GAA) sit 3' of r and are included via window[lo:].
    r = -1
    for i in range(n - 1, -1, -1):
        if window[i] != 'A':
            r = i
            break
    if r < 0:
        return 0, ''  # pure-A window: no multi-base repeat at any period

    def _walk(k: int) -> int:
        """Leftmost (5') boundary of the terminal period-k repeat block: the
        index of the furthest on-phase non-A base reachable from r without
        crossing a full period of A's (poly-A/body) or 2 consecutive off-phase
        non-A bases (left the repeat; a lone off-phase base is a tolerated
        basecall error)."""
        phase = r % k
        lo = r
        consec_a = 0
        off_phase = 0
        i = r - 1
        while i >= 0:
            b = window[i]
            if b == 'A':
                consec_a += 1
                off_phase = 0
                if consec_a >= k:
                    break
            else:
                consec_a = 0
                if i % k == phase:
                    lo = i
                    off_phase = 0
                else:
                    off_phase += 1
                    if off_phase > 1:
                        break
            i -= 1
        return lo

    ks = [k for k in k_range if n >= k * _MIN_REPEATS]
    if not ks:
        return 0, ''

    # Stage 1 — coarse extent: the longest terminal repeat over candidate periods
    # (smallest lo). A wrong period yields a shorter or accidentally-purer
    # sub-block; using the longest extent guarantees the true repeat is fully
    # captured so dominant_repeat_period can identify its real period below.
    coarse_lo = min(_walk(k) for k in ks)
    if n - coarse_lo < min_len:
        return 0, ''
    coarse = dominant_repeat_period(window[coarse_lo:], min_frac, k_range)
    if coarse is None or len(set(coarse[1])) <= 1:
        return 0, ''  # not a multi-base repeat (e.g. pure poly-A) -> strict no-op

    # Stage 2 — precise boundary at the true period, then gate.
    final_lo = _walk(coarse[0])
    block = window[final_lo:]
    if len(block) < min_len:
        return 0, ''
    res = dominant_repeat_period(block, min_frac, k_range)
    if res is None or len(set(res[1])) <= 1 or res[2] < min_frac:
        return 0, ''
    return len(block), res[1]


def find_polya_and_adapter(
    seq: str,
    max_error_rate: float = 0.0,
    max_consecutive_non_a: int = 1,
    min_polya: int = 1,
    adapter_window: int = _ADAPTER_WINDOW,
    strip_repeat: bool = True,
    repeat_min_len: int = _REPEAT_MIN_LEN,
    repeat_min_frac: float = _REPEAT_MIN_FRAC,
) -> Tuple[int, str, str, int, int, str]:
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
        (polya_len, adapter_seq, last_base, adapter_pass, repeat_len, repeat_motif) where:
          polya_len    - detected poly(A) tail length in bp (0 if none); measured on
                         the window *after* the repeat block (if any) is peeled
          adapter_seq  - detected adapter stub (empty string for pass 0)
          last_base    - last base of poly(A)-trimmed sequence (quality check)
          adapter_pass - 0/1/2 (see above)
          repeat_len   - length of the terminal (AAG/GAA)n artifact block peeled (0
                         if none / strip disabled). Sits 3' of polya_len; callers fold
                         it into total_trim = repeat_len + polya_len + len(adapter_seq).
          repeat_motif - canonical motif of the peeled block ('' if none)
    """
    # Restrict analysis to the last adapter_window bases for efficiency
    if len(seq) > adapter_window:
        window = seq[-adapter_window:]
    else:
        window = seq

    n = len(window)
    if n == 0:
        return 0, '', '', 0, 0, ''

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
        return 0, adapter_seq, '', adapter_pass, 0, ''

    # --- Terminal (AAG/GAA)n triplet-repeat artifact strip ---
    # Peel the basecaller artifact block (if any) BEFORE the poly-A scan: the
    # mis-basecalled tail has no terminal A-run, so _scan_polya cannot see it and the
    # artifact would otherwise be force-aligned. Peeling exposes any genuine poly-A
    # run 5' of the artifact for the scan below.
    repeat_len = 0
    repeat_motif = ''
    if strip_repeat:
        repeat_len, repeat_motif = _find_terminal_repeat_block(
            window, min_len=repeat_min_len, min_frac=repeat_min_frac,
        )
        if repeat_len:
            window = window[:-repeat_len]

    n = len(window)
    if n == 0:
        return 0, adapter_seq, '', adapter_pass, repeat_len, repeat_motif

    last_base = window[-1]
    polya_len = _scan_polya(window, max_error_rate, max_consecutive_non_a)

    if polya_len >= min_polya:
        return polya_len, adapter_seq, last_base, adapter_pass, repeat_len, repeat_motif

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
                return polya_len2, combined_adapter, candidate[-1], 2, repeat_len, repeat_motif

    return 0, adapter_seq, last_base, adapter_pass, repeat_len, repeat_motif


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
    strip_repeat: bool = True,
    repeat_min_len: int = _REPEAT_MIN_LEN,
    repeat_min_frac: float = _REPEAT_MIN_FRAC,
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
                # 668: secondary/supplementary records share the primary's QNAME; the trim
                # output is rebuilt un-flagged, so passing them through creates duplicate
                # QNAMEs the aligner panel's RN-keyed merge would silently collapse (and
                # SEQ-less secondaries were being written as empty "unmapped" records here).
                # One record per molecule: the primary. Same fix as _trim_region_task.
                if read.is_secondary or read.is_supplementary:
                    continue
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

                # Detect poly(A) + adapter (+ terminal (AAG/GAA)n artifact)
                polya_len, adapter_seq, _last_base, adapter_pass, repeat_len, repeat_motif = find_polya_and_adapter(
                    rna_seq,
                    max_error_rate=max_error_rate,
                    max_consecutive_non_a=max_consecutive_non_a,
                    adapter_window=adapter_window,
                    strip_repeat=strip_repeat,
                    repeat_min_len=repeat_min_len,
                    repeat_min_frac=repeat_min_frac,
                )
                stats['pass_counts'][adapter_pass] += 1

                total_trim = polya_len + repeat_len + len(adapter_seq)

                # Trim when a poly(A) tail (polya_len >= 1) OR a terminal repeat-
                # expansion artifact (repeat_len > 0) is detected. The artifact case
                # often has polya_len=0 (the (GAA)n replaced the tail outright) yet
                # must still be peeled so it is not force-aligned. Reads with neither
                # are genuine non-polyadenylated transcripts — left untrimmed.
                if (polya_len >= 1 or repeat_len > 0) and total_trim < len(rna_seq):
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
                    'repeat_len': repeat_len,
                    'repeat_motif': repeat_motif,
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
# FASTQ trimmer (for DRS reads without a Dorado-aligned BAM)
# ---------------------------------------------------------------------------


def trim_drs_fastq_polya(
    input_fastq_path: str,
    output_fastq_path: str,
    metadata_path: str,
    max_error_rate: float = 0.0,
    max_consecutive_non_a: int = 1,
    adapter_window: int = _ADAPTER_WINDOW,
    strip_repeat: bool = True,
    repeat_min_len: int = _REPEAT_MIN_LEN,
    repeat_min_frac: float = _REPEAT_MIN_FRAC,
) -> dict:
    """Trim poly(A) tails from a DRS FASTQ where reads are already 5'→3' RNA.

    Dorado writes DRS reads in 5'→3' orientation; poly(A) is always at the
    3' end of the stored sequence.  No alignment or strand flip is needed.
    The pt:i: tag is unavailable for FASTQ input; pt_tag is recorded as -1.

    Args:
        input_fastq_path:  Input FASTQ (.fastq or .fastq.gz / BGZF).
        output_fastq_path: Output trimmed FASTQ (.fastq.gz).
        metadata_path:     Destination for per-read metadata (.parquet or .tsv).

    Returns:
        Stats dict: total, trimmed, untrimmed, pass_counts {0,1,2}.
    """
    import gzip as _gzip

    stats: dict = {
        'total': 0,
        'trimmed': 0,
        'untrimmed': 0,
        'unmapped': 0,
        'pass_counts': {0: 0, 1: 0, 2: 0},
    }
    metadata_rows: list = []

    def _open_fq(path):
        s = str(path)
        if s.endswith('.gz') or s.endswith('.bgz'):
            return _gzip.open(s, 'rt')
        return open(s, 'r')

    def _readline(fh):
        try:
            return fh.readline()
        except EOFError:
            return ''

    with _open_fq(input_fastq_path) as fin, \
         _gzip.open(output_fastq_path, 'wt') as fout:
        while True:
            header = _readline(fin)
            if not header:
                break
            header = header.rstrip('\n')
            try:
                seq   = _readline(fin).rstrip('\n')
                plus  = _readline(fin).rstrip('\n')
                qual  = _readline(fin).rstrip('\n')
            except EOFError:
                break
            if not seq:
                break

            stats['total'] += 1
            read_id = header[1:].split()[0] if header.startswith('@') else header[1:]

            # DRS reads are already 5'→3'; poly-A is at the 3' end
            polya_len, adapter_seq, _last_base, adapter_pass, repeat_len, repeat_motif = find_polya_and_adapter(
                seq,
                max_error_rate=max_error_rate,
                max_consecutive_non_a=max_consecutive_non_a,
                adapter_window=adapter_window,
                strip_repeat=strip_repeat,
                repeat_min_len=repeat_min_len,
                repeat_min_frac=repeat_min_frac,
            )
            stats['pass_counts'][adapter_pass] += 1
            total_trim = polya_len + repeat_len + len(adapter_seq)

            if (polya_len >= 1 or repeat_len > 0) and total_trim < len(seq):
                trimmed_3prime_seq  = seq[-total_trim:]
                trimmed_3prime_qual = qual[-total_trim:] if qual else ''
                trimmed_seq  = seq[:-total_trim]
                trimmed_qual = qual[:-total_trim] if qual else ''
                stats['trimmed'] += 1
            else:
                trimmed_3prime_seq  = ''
                trimmed_3prime_qual = ''
                trimmed_seq  = seq
                trimmed_qual = qual
                stats['untrimmed'] += 1

            fout.write(f"{header}\n{trimmed_seq}\n{plus}\n{trimmed_qual}\n")

            quals_str = ','.join(str(ord(c) - 33) for c in trimmed_3prime_qual) if trimmed_3prime_qual else ''
            metadata_rows.append({
                'read_id':             read_id,
                'strand':              '+',
                'polya_len':           polya_len,
                'adapter_seq':         adapter_seq[:50],
                'adapter_pass':        adapter_pass,
                'repeat_len':          repeat_len,
                'repeat_motif':        repeat_motif,
                'pt_tag':              -1,
                'trimmed_3prime_seq':  trimmed_3prime_seq,
                'trimmed_3prime_quals': quals_str,
                'original_seq_len':    len(seq),
                'trimmed_seq_len':     len(trimmed_seq),
            })

    _write_metadata(metadata_rows, Path(metadata_path))
    logger.info(
        "trim_drs_fastq_polya: total=%d trimmed=%d untrimmed=%d "
        "pass0=%d pass1=%d pass2=%d",
        stats['total'], stats['trimmed'], stats['untrimmed'],
        stats['pass_counts'][0], stats['pass_counts'][1], stats['pass_counts'][2],
    )
    return stats


# ---------------------------------------------------------------------------
# Per-read helper shared by sequential and parallel paths
# ---------------------------------------------------------------------------


def _process_mapped_read(
    read: pysam.AlignedSegment,
    unaligned_header: pysam.AlignmentHeader,
    max_error_rate: float,
    max_consecutive_non_a: int,
    adapter_window: int,
    strip_repeat: bool = True,
    repeat_min_len: int = _REPEAT_MIN_LEN,
    repeat_min_frac: float = _REPEAT_MIN_FRAC,
) -> Tuple[pysam.AlignedSegment, Dict, int, bool]:
    """Trim poly-A from one mapped read; return (rec, row, adapter_pass, was_trimmed)."""
    seq = read.query_sequence or ""
    strand = "-" if read.is_reverse else "+"
    quals = list(read.query_qualities or [])
    rna_seq = reverse_complement(seq) if strand == "-" else seq
    rna_quals = quals[::-1] if strand == "-" else quals
    try:
        pt_tag = int(read.get_tag("pt"))
    except (KeyError, ValueError):
        pt_tag = -1
    polya_len, adapter_seq, _last, adapter_pass, repeat_len, repeat_motif = find_polya_and_adapter(
        rna_seq, max_error_rate=max_error_rate,
        max_consecutive_non_a=max_consecutive_non_a, adapter_window=adapter_window,
        strip_repeat=strip_repeat, repeat_min_len=repeat_min_len,
        repeat_min_frac=repeat_min_frac,
    )
    total_trim = polya_len + repeat_len + len(adapter_seq)
    if (polya_len >= 1 or repeat_len > 0) and total_trim < len(rna_seq):
        trimmed_rna_seq = rna_seq[:-total_trim]
        trimmed_rna_quals = rna_quals[:-total_trim] if rna_quals else []
        trimmed_3p_seq = rna_seq[-total_trim:]
        trimmed_3p_quals = rna_quals[-total_trim:] if rna_quals else []
        was_trimmed = True
    else:
        trimmed_rna_seq = rna_seq
        trimmed_rna_quals = rna_quals
        trimmed_3p_seq = ""
        trimmed_3p_quals = []
        was_trimmed = False
    rec = _make_unaligned_record(read, trimmed_rna_seq, trimmed_rna_quals, unaligned_header)
    quals_str = ",".join(str(q) for q in trimmed_3p_quals) if trimmed_3p_quals else ""
    row: Dict = {
        "read_id": read.query_name,
        "strand": strand,
        "polya_len": polya_len,
        "adapter_seq": adapter_seq[:50],
        "adapter_pass": adapter_pass,
        "repeat_len": repeat_len,
        "repeat_motif": repeat_motif,
        "pt_tag": pt_tag,
        "trimmed_3prime_seq": trimmed_3p_seq,
        "trimmed_3prime_quals": quals_str,
        "original_seq_len": len(seq),
        "trimmed_seq_len": len(trimmed_rna_seq),
    }
    return rec, row, adapter_pass, was_trimmed


# ---------------------------------------------------------------------------
# Parallel BAM trimming
# ---------------------------------------------------------------------------


def _trim_region_task(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    out_bam_path: str,
    max_error_rate: float,
    max_consecutive_non_a: int,
    adapter_window: int,
    strip_repeat: bool = True,
    repeat_min_len: int = _REPEAT_MIN_LEN,
    repeat_min_frac: float = _REPEAT_MIN_FRAC,
) -> Tuple[List[Dict], int, int]:
    """Trim one BAM region in a worker process.

    Returns (metadata_rows, n_trimmed, n_untrimmed).
    Unmapped reads are skipped here; handled by _trim_unmapped_task.
    """
    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unsorted"}})
    rows: List[Dict] = []
    n_trimmed = n_untrimmed = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam_in, \
         pysam.AlignmentFile(out_bam_path, "wb", header=hdr) as bam_out:
        for read in bam_in.fetch(chrom, start, end):
            if read.is_unmapped or not read.query_sequence:
                continue
            # 668: two duplicate-QNAME sources, both real (caught by align's RN-collision
            # guard: 130,311 dup QNAMEs / 3.2% on wt_rep1):
            #  (a) fetch() returns reads OVERLAPPING the region, so a boundary-spanning read
            #      is emitted by BOTH neighbouring workers — the same defect class fixed in
            #      cdna_correct (624). Partition on reference_start: regions tile disjointly,
            #      so each read starts in exactly one region.
            #  (b) secondary/supplementary records share the primary's QNAME, and the trim
            #      output is rebuilt UN-flagged, so they can't be filtered downstream. Trim
            #      wants one record per molecule: the primary.
            if read.is_secondary or read.is_supplementary:
                continue
            if not (start <= read.reference_start < end):
                continue
            rec, row, _, was_trimmed = _process_mapped_read(
                read, hdr, max_error_rate, max_consecutive_non_a, adapter_window,
                strip_repeat=strip_repeat, repeat_min_len=repeat_min_len,
                repeat_min_frac=repeat_min_frac,
            )
            bam_out.write(rec)
            rows.append(row)
            if was_trimmed:
                n_trimmed += 1
            else:
                n_untrimmed += 1
    return rows, n_trimmed, n_untrimmed


def _trim_unmapped_task(
    bam_path: str,
    out_bam_path: str,
) -> List[Dict]:
    """Pass unmapped reads through to output BAM unchanged. Returns metadata rows."""
    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unsorted"}})
    rows: List[Dict] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam_in, \
         pysam.AlignmentFile(out_bam_path, "wb", header=hdr) as bam_out:
        try:
            it = bam_in.fetch(contig="*")
        except (ValueError, KeyError):
            return rows
        for read in it:
            seq = read.query_sequence
            rec = _make_unaligned_record(
                read, seq or "", list(read.query_qualities or []) if seq else None, hdr
            )
            bam_out.write(rec)
            if seq is not None:
                try:
                    pt_tag = int(read.get_tag("pt"))
                except (KeyError, ValueError):
                    pt_tag = -1
                rows.append({
                    "read_id": read.query_name,
                    "strand": "+",
                    "polya_len": 0,
                    "adapter_seq": "",
                    "adapter_pass": 0,
                    "repeat_len": 0,
                    "repeat_motif": "",
                    "pt_tag": pt_tag,
                    "trimmed_3prime_seq": "",
                    "trimmed_3prime_quals": "",
                    "original_seq_len": len(seq),
                    "trimmed_seq_len": len(seq),
                })
    return rows


def _trim_drs_bam_polya_parallel(
    input_bam_path: str,
    output_bam_path: str,
    metadata_path: Path,
    max_error_rate: float = 0.0,
    max_consecutive_non_a: int = 1,
    adapter_window: int = _ADAPTER_WINDOW,
    workers: int = 4,
    strip_repeat: bool = True,
    repeat_min_len: int = _REPEAT_MIN_LEN,
    repeat_min_frac: float = _REPEAT_MIN_FRAC,
) -> Dict:
    """Region-parallel poly-A trimmer.

    Splits the input BAM by coordinate region, processes each region in a
    worker process, handles unmapped reads in-process, then concatenates
    outputs via ``samtools cat``.

    Axis-2 safety: ``plan_regions`` opens and closes pysam before the pool is
    created. Workers receive a string path and open their own fresh handles.
    """
    from rectify.core.bam.regions import plan_regions

    input_bam_path = str(input_bam_path)
    bai = input_bam_path + ".bai"
    alt_bai = input_bam_path.replace(".bam", ".bai")
    if not Path(bai).exists() and not Path(alt_bai).exists():
        logger.info("No BAM index found — indexing %s", input_bam_path)
        subprocess.run(["samtools", "index", input_bam_path], check=True)

    with tempfile.TemporaryDirectory(prefix="drs_trim_par_") as tmp_str:
        tmp = Path(tmp_str)
        plans = plan_regions(input_bam_path, tmp)

        futures: Dict = {}
        per_region_bams: List[str] = []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for plan in plans:
                out_bam = str(tmp / f"{plan.region_id}.bam")
                per_region_bams.append(out_bam)
                fut = pool.submit(
                    _trim_region_task,
                    input_bam_path, plan.chrom, plan.start, plan.end,
                    out_bam, max_error_rate, max_consecutive_non_a, adapter_window,
                    strip_repeat, repeat_min_len, repeat_min_frac,
                )
                futures[fut] = plan.region_id

            all_rows: List[Dict] = []
            n_trimmed = n_untrimmed = 0
            pass_counts: Dict[int, int] = {0: 0, 1: 0, 2: 0}
            for fut in as_completed(futures):
                rows, nt, nu = fut.result()
                all_rows.extend(rows)
                n_trimmed += nt
                n_untrimmed += nu
                for row in rows:
                    p = row["adapter_pass"]
                    pass_counts[p] = pass_counts.get(p, 0) + 1

        unmapped_bam = str(tmp / "unmapped.bam")
        unmap_rows = _trim_unmapped_task(input_bam_path, unmapped_bam)
        all_rows.extend(unmap_rows)

        existing_bams = [b for b in per_region_bams + [unmapped_bam] if Path(b).exists()]
        if existing_bams:
            subprocess.run(
                ["samtools", "cat", "-o", output_bam_path] + existing_bams,
                check=True,
            )
        else:
            with pysam.AlignmentFile(input_bam_path, "rb") as _src:
                _hdr = _make_unaligned_header(_src.header)
            with pysam.AlignmentFile(output_bam_path, "wb", header=_hdr):
                pass

        _write_metadata(all_rows, metadata_path)

        stats: Dict = {
            "total": n_trimmed + n_untrimmed + len(unmap_rows),
            "trimmed": n_trimmed,
            "untrimmed": n_untrimmed,
            "unmapped": len(unmap_rows),
            "pass_counts": pass_counts,
        }
        logger.info(
            "trim_drs_bam_polya (parallel %d workers): total=%d trimmed=%d untrimmed=%d unmapped=%d",
            workers, stats["total"], stats["trimmed"], stats["untrimmed"], stats["unmapped"],
        )
        return stats


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

_FASTQ_SUFFIXES = {'.fastq', '.fq', '.fastq.gz', '.fq.gz'}


def _is_fastq(path: Path) -> bool:
    return path.suffix in {'.fastq', '.fq'} or str(path).endswith(('.fastq.gz', '.fq.gz'))


def run(args) -> int:
    """Entry point called from rectify/cli.py for `rectify trim-polya`."""
    import sys
    from datetime import datetime as _dt, timezone as _tz
    from time import perf_counter as _perf

    input_path = Path(args.input_bam)  # accepts BAM or FASTQ
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
        return 1

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_stem = getattr(args, 'prefix', None) or input_path.stem
    # strip double extension for .fastq.gz
    if sample_stem.endswith('.fastq') or sample_stem.endswith('.fq'):
        sample_stem = Path(sample_stem).stem

    # --- Resume skip-check ---
    _input_key = 'input_fastq' if _is_fastq(input_path) else 'input_bam'
    _current_inputs = {_input_key: str(input_path)}
    from rectify.core.provenance import can_skip_stage
    from rectify.utils.version import get_rectify_git_sha as _get_sha
    _drs_trim_sha = _get_sha()
    _skip = can_skip_stage(
        stage='drs_trim',
        sample_output=output_dir,
        sample_id=sample_stem,
        current_inputs=_current_inputs,
        current_argv=sys.argv,
        rectify_git_sha=_drs_trim_sha,
        force=getattr(args, 'force_all', False),
        force_stages=set(getattr(args, 'force_stage', '').split(','))
                     if getattr(args, 'force_stage', None) else None,
        accept_prior_provenance=getattr(args, 'accept_prior_provenance', False),
    )
    logger.info("[RESUME] sample=%s stage=drs_trim decision=%s reason=%s",
                sample_stem, 'SKIP' if _skip.skip else 'RUN', _skip.reason)
    if _skip.skip:
        return 0
    if getattr(args, 'dry_run_resume', False):
        print(f"[dry-run-resume] stage=drs_trim decision=RUN reason={_skip.reason}")
        return 0

    _stage_started_at = _dt.now(_tz.utc).isoformat()
    _t_start = _perf()

    use_tsv = getattr(args, 'tsv', False)
    meta_ext = '.tsv' if use_tsv else '.parquet'
    metadata_path = output_dir / f"{sample_stem}_polya_trim_metadata{meta_ext}"

    print(f"Input:       {input_path}")
    print(f"Output dir:  {output_dir}")
    print(f"Metadata:    {metadata_path}")
    print(f"max_error_rate:          {args.max_error_rate}")
    print(f"max_consecutive_non_a:   {args.max_consecutive_non_a}")
    print(f"adapter_window:          {args.adapter_window}")

    _strip_repeat = getattr(args, 'strip_repeat_expansion', True)
    _repeat_min_len = getattr(args, 'repeat_min_len', _REPEAT_MIN_LEN)
    _repeat_min_frac = getattr(args, 'repeat_min_frac', _REPEAT_MIN_FRAC)
    print(f"strip_repeat_expansion:  {_strip_repeat}")
    if _strip_repeat:
        print(f"repeat_min_len:          {_repeat_min_len}")
        print(f"repeat_min_frac:         {_repeat_min_frac}")

    _workers = max(1, getattr(args, 'workers', 1) or 1)

    if _is_fastq(input_path):
        # FASTQ path: reads are 5'→3' RNA, poly-A at 3' end — no alignment needed
        output_fastq = output_dir / f"{sample_stem}_trimmed.fastq.gz"
        print(f"Mode:        FASTQ (DRS 5'→3', pt_tag=-1)")
        print(f"Output FASTQ: {output_fastq}")
        stats = trim_drs_fastq_polya(
            input_fastq_path=str(input_path),
            output_fastq_path=str(output_fastq),
            metadata_path=str(metadata_path),
            max_error_rate=args.max_error_rate,
            max_consecutive_non_a=args.max_consecutive_non_a,
            adapter_window=args.adapter_window,
            strip_repeat=_strip_repeat,
            repeat_min_len=_repeat_min_len,
            repeat_min_frac=_repeat_min_frac,
        )
        _outputs = {'output_fastq': str(output_fastq), 'metadata': str(metadata_path)}
    else:
        # BAM path: standard Dorado-aligned DRS BAM
        output_bam = output_dir / f"{sample_stem}.bam"
        print(f"Mode:        BAM")
        print(f"Output BAM:  {output_bam}")
        if _workers > 1:
            print(f"Workers:     {_workers} (parallel regions)")
            stats = _trim_drs_bam_polya_parallel(
                input_bam_path=str(input_path),
                output_bam_path=str(output_bam),
                metadata_path=metadata_path,
                max_error_rate=args.max_error_rate,
                max_consecutive_non_a=args.max_consecutive_non_a,
                adapter_window=args.adapter_window,
                workers=_workers,
                strip_repeat=_strip_repeat,
                repeat_min_len=_repeat_min_len,
                repeat_min_frac=_repeat_min_frac,
            )
        else:
            stats = trim_drs_bam_polya(
                input_bam_path=str(input_path),
                output_bam_path=str(output_bam),
                metadata_path=str(metadata_path),
                max_error_rate=args.max_error_rate,
                max_consecutive_non_a=args.max_consecutive_non_a,
                adapter_window=args.adapter_window,
                strip_repeat=_strip_repeat,
                repeat_min_len=_repeat_min_len,
                repeat_min_frac=_repeat_min_frac,
            )
        _outputs = {'output_bam': str(output_bam), 'metadata': str(metadata_path)}

    print(f"\nTrim summary:")
    print(f"  Total reads:      {stats['total']:,}")
    print(f"  Trimmed:          {stats['trimmed']:,}")
    print(f"  Untrimmed:        {stats['untrimmed']:,}")
    if stats.get('unmapped', 0):
        print(f"  Unmapped (skip):  {stats['unmapped']:,}")
    print(f"  Pass 0 (no stub): {stats['pass_counts'][0]:,}")
    print(f"  Pass 1 (regex):   {stats['pass_counts'][1]:,}")
    print(f"  Pass 2 (peel):    {stats['pass_counts'][2]:,}")

    # Write stage sidecar (non-fatal on failure)
    try:
        from rectify.core.provenance import ProvenanceRecord, write_stage_sidecar
        _sc_record = ProvenanceRecord.from_components(
            stage='drs_trim',
            sample_id=sample_stem,
            sample_output_dir=output_dir,
            started_at=_stage_started_at,
            completed_at=_dt.now(_tz.utc).isoformat(),
            exit_status=0,
            inputs=_current_inputs,
            outputs=_outputs,
            stats={
                'wall_seconds': _perf() - _t_start,
                'total': stats['total'],
                'trimmed': stats['trimmed'],
                'untrimmed': stats['untrimmed'],
                'unmapped': stats.get('unmapped', 0),
            },
            argv=sys.argv,
            rectify_git_sha=_drs_trim_sha,
        )
        write_stage_sidecar(_sc_record, sample_output=output_dir)
    except Exception as _sc_exc:
        logger.warning("Failed to write drs_trim sidecar: %s", _sc_exc)

    return 0


def create_trim_polya_parser(subparsers):
    """Wire the `trim-polya` subcommand into the given subparsers group."""
    import argparse
    # =========================================================================
    # trim-polya command (DRS poly(A)+adapter pre-trimming)
    # =========================================================================
    trim_polya_parser = subparsers.add_parser(
        'trim-polya',
        help='Trim poly(A) tail + adapter from DRS BAM or FASTQ before alignment',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    trim_polya_parser.add_argument(
        'input_bam',
        help=(
            'Input file: Dorado-aligned BAM (with optional pt:i: tags), '
            'OR a raw DRS FASTQ/FASTQ.GZ where reads are already 5\'→3\' RNA orientation. '
            'FASTQ mode: pt_tag stored as -1; poly-A detected from sequence directly.'
        ),
    )
    trim_polya_parser.add_argument(
        '-o', '--output-dir',
        dest='output_dir',
        required=True,
        help='Output directory for trimmed BAM and metadata parquet',
    )
    trim_polya_parser.add_argument(
        '--adapter-window',
        type=int,
        default=150,
        help='Bases from 3\' end to search for poly(A)+adapter',
    )
    trim_polya_parser.add_argument(
        '--max-error-rate',
        type=float,
        default=0.0,
        help='Max cumulative non-A fraction for poly(A) scan (default 0.0 = strict pure-A)',
    )
    trim_polya_parser.add_argument(
        '--max-consecutive-non-a',
        type=int,
        default=1,
        dest='max_consecutive_non_a',
        help=(
            'Stop poly(A) scan after this many consecutive non-A bases (default 1 = '
            'stop at 2+ consecutive). Prevents consuming upstream gene-body sequence '
            'even when --max-error-rate > 0.'
        ),
    )
    trim_polya_parser.add_argument(
        '--strip-repeat-expansion',
        dest='strip_repeat_expansion',
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            'Peel a terminal (AAG/GAA)n triplet-repeat basecaller artifact (Dorado '
            'RNA004) before the poly(A) scan, exposing any genuine poly(A) behind it. '
            'Strict no-op on clean poly(A) tails. Use --no-strip-repeat-expansion to '
            'disable.'
        ),
    )
    trim_polya_parser.add_argument(
        '--repeat-min-len',
        dest='repeat_min_len',
        type=int,
        default=_REPEAT_MIN_LEN,
        help='Minimum terminal repeat-block length (bp) to peel (default: %(default)s).',
    )
    trim_polya_parser.add_argument(
        '--repeat-min-frac',
        dest='repeat_min_frac',
        type=float,
        default=_REPEAT_MIN_FRAC,
        help='Minimum canonical-motif purity to peel a repeat block (default: %(default)s).',
    )
    trim_polya_parser.add_argument(
        '--prefix',
        default=None,
        help='Output file prefix (default: derived from input BAM stem)',
    )
    trim_polya_parser.add_argument(
        '--tsv',
        action='store_true',
        default=False,
        help='Write metadata as TSV instead of parquet',
    )
    trim_polya_parser.add_argument(
        '--workers',
        type=int,
        default=1,
        help='Number of parallel worker processes for BAM trimming (default: 1 = sequential). '
             'Requires a sorted, indexed input BAM. Auto-indexes if .bai is absent.',
    )
    from rectify.core.provenance.cli import add_resume_args
    add_resume_args(trim_polya_parser)
