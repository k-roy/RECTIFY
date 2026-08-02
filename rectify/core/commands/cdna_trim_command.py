#!/usr/bin/env python3
"""
cDNA poly(A) pre-trimming for RECTIFY.

Trims poly(A) tails + adapter stubs from Nanopore cDNA FASTQ files BEFORE
alignment.  Also handles optional 5' poly(T) trimming for antisense reads.

Unlike ``rectify trim-polya`` (which requires Dorado-aligned BAMs with ``pt:i:``
tags), this module operates directly on raw FASTQ files — suitable for cDNA data
where Dorado unaligned BAMs with poly-A signal estimates are not available.

Algorithm
---------
3' end  (RNA-sense reads, poly-A at 3'):
  Pass 0 — pure poly-A tail (no adapter stub) — trim only the poly-A.
  Pass 1 — adapter anchor detected after poly-A: trim poly-A + stub.
  Pass 2 — iterative peel for A-basecalling errors near the adapter boundary.

5' end  (antisense reads, poly-T at 5'):
  Optional (``--trim-5p-polyt``): trim poly-T run from the 5' end.  Mirrors the
  3' algorithm but scans from the left using the complement base (T instead of A).

Adapter defaults
----------------
3' adapter (after poly-A on RNA-sense reads):
  CRTA  = CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG  (45 bp)
  Seed  = first 12 bp = CTTGCGGGCGGC  (highly specific, tolerates 1 mismatch)

5' adapter (before poly-T on antisense reads):
  CRTA_RC = CTTACCTGCGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG  (47 bp)
  Seed    = LAST 12 bp = GCCGCCCGCAAG  (the boundary that abuts the poly-T)

  The 5' seed is taken from the adapter's 3' end, mirroring the 3' path (which
  seeds on the adapter's 5' end because that is the boundary abutting the
  poly-A).  Measured on PBM17777 WT_BY4742_rep1 (planning/541): the last-12
  seed matches 92.5% of antisense reads, the first-12 seed 0% — the observed
  PCB114 primer diverges from CRTA_RC over its 5' ~12 bp but is exact over the
  35 bp abutting the poly-T.

Pass why these are "longer adapters": the DRS trim-polya module uses a 1-11 bp
stub regex (``T[CT]{0,10}$``).  The cDNA CRTA adapter is 45 bp — substantially
longer — and requires a k-mer anchor strategy rather than a short regex.

Author: Kevin R. Roy
Date: 2026-05-06
"""

import gzip
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .drs_trim_command import _scan_polya, _write_metadata

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default adapter sequences (CRTA template-switching adapter and its RC)
# ---------------------------------------------------------------------------

_CRTA         = "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG"
_CRTA_RC      = "CTTACCTGCGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG"

# Default seed lengths and mismatch tolerance for adapter detection
_SEED_LEN          = 12
_SEED_MAX_MISMATCHES = 1
_ADAPTER_WINDOW    = 200   # bp from 3' (or 5') end to search in
_PASS2_MAX_STUB    = 20    # max stub length to try peeling in pass 2
_MIN_POLYA_PASS2   = 5     # min poly-A length to accept a pass-2 call

# Minimum tail length required before a tail is trusted as an ORIENTATION label
# (as opposed to merely being trimmed).  Short A/T runs are common as genomic
# sequence: on PBM17777 WT_BY4742_rep1, 9.4% of reads carry a 1-4 nt 3' "poly-A"
# with no poly-T, which is far too weak to call a read sense.  Trimming keeps its
# own (looser) min_polya / min_polyt thresholds; only the label uses this one.
_MIN_ORIENTATION_TAIL = 10

# Minimum poly-T run accepted when NO 5' adapter anchor was found.  With an
# anchor the run position is known and any run length is meaningful; without
# one we are scanning naked sequence and a short T-run is more likely genomic
# than a reverse-transcription priming site, so require a real homopolymer.
_MIN_POLYT_UNANCHORED = 10


# ---------------------------------------------------------------------------
# Adapter anchor detection helpers
# ---------------------------------------------------------------------------

def _hamming(a: str, b: str) -> int:
    """Count mismatches between two equal-length strings."""
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def _find_adapter_anchor(seq: str, seed: str, max_mismatches: int) -> int:
    """Find the leftmost position of *seed* in *seq* with ≤ max_mismatches.

    Returns the start index of the match, or -1 if not found.

    The search is restricted to ``len(seed)``-long windows from the left so
    that a full seed comparison is always possible.
    """
    k = len(seed)
    n = len(seq)
    for i in range(n - k + 1):
        if _hamming(seq[i:i + k], seed) <= max_mismatches:
            return i
    return -1


def _find_adapter_anchor_rightmost(seq: str, seed: str, max_mismatches: int) -> int:
    """Rightmost position of *seed* in *seq* with ≤ max_mismatches (-1 if none).

    Used for the 5' end: the poly-T abuts the LAST adapter copy, so on a
    chimeric read carrying two adapter copies the rightmost one is the correct
    boundary.  (The 3' path wants the leftmost for the mirror-image reason.)
    """
    k = len(seed)
    for i in range(len(seq) - k, -1, -1):
        if _hamming(seq[i:i + k], seed) <= max_mismatches:
            return i
    return -1


def _scan_polyt(window: str, start: int, max_error_rate: float,
                max_consecutive_non_t: int) -> int:
    """Length of the poly-T run beginning at ``window[start]``.

    Same tolerance semantics as :func:`_scan_polya` (which scans a poly-A from
    the right), but scans rightward from an arbitrary offset.
    """
    n = len(window)
    if start >= n:
        return 0
    errors = 0
    total = 0
    polyt_end = start
    consecutive_non_t = 0
    for i in range(start, n):
        total += 1
        if window[i] != 'T':
            errors += 1
            consecutive_non_t += 1
        else:
            consecutive_non_t = 0
        if errors / total > max_error_rate:
            break
        if consecutive_non_t > max_consecutive_non_t:
            break
        polyt_end = i + 1
    # Back off to the last T if the error budget over-consumed
    while polyt_end > start and window[polyt_end - 1] != 'T':
        polyt_end -= 1
    return polyt_end - start


# ---------------------------------------------------------------------------
# Core detection: 3' poly-A + adapter
# ---------------------------------------------------------------------------

def find_cdna_3p_polya_and_adapter(
    seq: str,
    adapter_seed: str = _CRTA[:_SEED_LEN],
    adapter_full: str = _CRTA,
    seed_max_mismatches: int = _SEED_MAX_MISMATCHES,
    max_error_rate: float = 0.0,
    max_consecutive_non_a: int = 1,
    adapter_window: int = _ADAPTER_WINDOW,
    min_polya: int = 1,
) -> Tuple[int, str, str, int]:
    """Detect poly-A tail and 3' adapter stub (e.g. CRTA) from a read.

    Three-pass algorithm identical in structure to ``find_polya_and_adapter``
    in ``drs_trim_command``, but uses a configurable k-mer seed (tolerating
    up to *seed_max_mismatches*) instead of the short T[CT] DRS stub regex.

    Args:
        seq: Full read sequence (as-sequenced; RNA 5'→3' for RNA-sense reads,
             RC orientation for antisense reads that start with poly-T).
        adapter_seed: First ``_SEED_LEN`` bases of the expected 3' adapter.
        adapter_full: Full adapter sequence (for trimming the detected stub).
        seed_max_mismatches: Mismatches allowed in seed match (default 1).
        max_error_rate: Tolerated non-A fraction for poly-A scan (default 0.0).
        max_consecutive_non_a: Stop poly-A scan at ≥(this+1) consecutive non-A.
        adapter_window: Bases from 3' end to inspect (default 200).
        min_polya: Minimum poly-A length to report (default 1).

    Returns:
        (polya_len, adapter_seq, last_base, adapter_pass) where:
          polya_len    – number of poly-A bases trimmed (0 if none)
          adapter_seq  – detected adapter stub (empty string if none found)
          last_base    – last base of poly-A-trimmed sequence
          adapter_pass – 0 (no stub), 1 (seed match), or 2 (iterative peel)
    """
    n = len(seq)
    if n == 0:
        return 0, '', '', 0

    # Restrict analysis window to the last adapter_window bases
    offset = max(0, n - adapter_window)
    window = seq[offset:]
    wn = len(window)

    # --- Pass 1: adapter seed detection ---
    anchor_pos = _find_adapter_anchor(window, adapter_seed, seed_max_mismatches)
    if anchor_pos >= 0:
        adapter_seq   = window[anchor_pos:]
        window        = window[:anchor_pos]
        adapter_pass  = 1
    else:
        adapter_seq   = ''
        adapter_pass  = 0

    if not window:
        return 0, adapter_seq, '', adapter_pass

    last_base = window[-1]
    polya_len = _scan_polya(window, max_error_rate, max_consecutive_non_a)

    if polya_len >= min_polya:
        return polya_len, adapter_seq, last_base, adapter_pass

    # --- Pass 2: iterative peel for A-basecalling errors near adapter ---
    # Only activates when the last base is not 'A' (otherwise polya_len≥1 above)
    if last_base != 'A':
        for k in range(1, _PASS2_MAX_STUB + 1):
            if len(window) <= k:
                break
            candidate  = window[:-k]
            if not candidate or candidate[-1] != 'A':
                continue
            stub = window[-k:]
            # Require at least one non-A base that could be an adapter boundary
            if all(b == 'A' for b in stub):
                continue
            polya_len2 = _scan_polya(candidate, max_error_rate, max_consecutive_non_a)
            if polya_len2 >= _MIN_POLYA_PASS2:
                combined_adapter = stub + adapter_seq
                return polya_len2, combined_adapter, candidate[-1], 2

    return 0, adapter_seq, last_base, adapter_pass


# ---------------------------------------------------------------------------
# 5' poly-T detection (antisense reads)
# ---------------------------------------------------------------------------

def find_cdna_5p_polyt(
    seq: str,
    adapter_seed: str = _CRTA_RC[-_SEED_LEN:],
    seed_max_mismatches: int = _SEED_MAX_MISMATCHES,
    max_error_rate: float = 0.0,
    max_consecutive_non_t: int = 1,
    adapter_window: int = _ADAPTER_WINDOW,
    min_polyt: int = 1,
    min_polyt_unanchored: int = _MIN_POLYT_UNANCHORED,
) -> Tuple[int, int]:
    """Detect the poly-T run at the 5' end of an antisense cDNA read.

    Returns ``(polyt_len, polyt_start)`` — the run length and its 0-based start
    index in *seq*.  ``(0, -1)`` when no run is detected.

    Two passes, mirroring :func:`find_cdna_3p_polya_and_adapter`:

    Pass 1 (**anchored**) — locate the 5' adapter seed (the LAST ``_SEED_LEN``
    bases of the 5' adapter, i.e. the boundary abutting the poly-T) anywhere in
    the first *adapter_window* bases, then scan the poly-T rightward from the
    end of that match.  This is the pass that fires in practice.

    Pass 2 (**unanchored**) — no adapter seed found: scan the window for the
    first T-run of at least *min_polyt_unanchored* bases.  Covers reads whose
    adapter was trimmed upstream, and reads whose adapter basecall was too poor
    to seed.  The stricter length floor keeps genomic T-tracts out.

    .. note::
       Before 2026-08-01 this function walked from ``seq[0]`` and stopped at the
       first non-T, so it could only fire when the read *literally began* with
       T's.  On SQK-PCB114 the poly-T sits behind ~120-140 nt of
       adapter/barcode/UMI, so it fired on **0.11%** of reads instead of ~37%,
       and the ONT-cDNA orientation label was effectively absent.  See
       ``planning/541_ont_cdna_strand_fix.md`` section 4.

    Args:
        seq: Read sequence (antisense reads carry poly-T behind the 5' adapter).
        adapter_seed: Seed for the 5' adapter boundary (default: last 12 bases
            of CRTA_RC).  Pass ``''`` to force the unanchored pass.
        seed_max_mismatches: Mismatches allowed in the seed match (default 1).
        max_error_rate: Tolerated non-T fraction (default 0.0 = strict).
        max_consecutive_non_t: Stop at ≥(this+1) consecutive non-T (default 1).
        adapter_window: Bases from 5' end to inspect (default 200).
        min_polyt: Minimum poly-T length to report when anchored (default 1).
        min_polyt_unanchored: Minimum length to report when NOT anchored.

    Returns:
        ``(polyt_len, polyt_start)``.
    """
    window = seq[:adapter_window]
    if not window:
        return 0, -1

    # --- Pass 1: adapter-anchored ---
    if adapter_seed:
        anchor = _find_adapter_anchor_rightmost(window, adapter_seed, seed_max_mismatches)
        if anchor >= 0:
            start = anchor + len(adapter_seed)
            run = _scan_polyt(window, start, max_error_rate, max_consecutive_non_t)
            if run >= min_polyt:
                return run, start
            # Adapter present but no poly-T behind it — this is a sense read
            # whose 5' end carries the adapter, not an antisense read.  Do NOT
            # fall through to the unanchored scan: that would pick up a genomic
            # T-tract and mislabel the read's orientation.
            return 0, -1

    # --- Pass 2: unanchored ---
    n = len(window)
    i = 0
    while i < n:
        if window[i] != 'T':
            i += 1
            continue
        run = _scan_polyt(window, i, max_error_rate, max_consecutive_non_t)
        if run >= max(min_polyt, min_polyt_unanchored):
            return run, i
        i += max(run, 1)
    return 0, -1


# ---------------------------------------------------------------------------
# FASTQ I/O helpers
# ---------------------------------------------------------------------------

def _open_fastq(path: Path):
    """Return an open text-mode file handle for a FASTQ (plain or gzip)."""
    s = str(path)
    if s.endswith('.gz') or s.endswith('.gzip'):
        return gzip.open(s, 'rt')
    return open(s, 'r')


def _iter_fastq(path: Path):
    """Yield (name, seq, quals) tuples from a FASTQ file.

    *name* is the full header line without the leading '@'.
    """
    with _open_fastq(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq   = fh.readline().rstrip('\n')
            plus  = fh.readline()    # '+' line
            quals = fh.readline().rstrip('\n')
            if not (header and seq and plus and quals):
                break
            yield header[1:].rstrip('\n'), seq, quals


# ---------------------------------------------------------------------------
# Main trim function
# ---------------------------------------------------------------------------

def trim_cdna_fastq_polya(
    input_fastq_path: str,
    output_fastq_path: str,
    metadata_path: str,
    adapter_3p: str = _CRTA,
    adapter_5p: Optional[str] = _CRTA_RC,
    seed_len: int = _SEED_LEN,
    seed_max_mismatches: int = _SEED_MAX_MISMATCHES,
    max_error_rate: float = 0.0,
    max_consecutive_non_a: int = 1,
    adapter_window: int = _ADAPTER_WINDOW,
    trim_5p_polyt: bool = False,
    min_polya: int = 1,
    min_polyt: int = 1,
) -> Dict:
    """Trim poly-A tails and adapter stubs from a cDNA Nanopore FASTQ file.

    Reads are processed in FASTQ orientation (no RC is applied — the read is
    assumed to be already in its as-sequenced form).  Poly-A detection scans
    the 3' end; if ``trim_5p_polyt=True``, poly-T detection also scans the 5'
    end for antisense reads.

    The read name written to the output FASTQ is the bare UUID (everything
    before the first space in the FASTQ header), consistent with the bare-UUID
    policy throughout the pipeline.

    Args:
        input_fastq_path:  Source FASTQ (plain or gzip).
        output_fastq_path: Destination FASTQ (plain or gzip, auto-detected
                           from extension).
        metadata_path:     Per-read trim metadata (.parquet or .tsv).
        adapter_3p:        Full 3' adapter sequence (CRTA by default).
        adapter_5p:        Full 5' adapter sequence (CRTA_RC by default, used
                           only for seed derivation; set to None to disable 5'
                           adapter seed detection).
        seed_len:          k-mer seed length for adapter detection (default 12).
        seed_max_mismatches: Mismatches allowed in k-mer seed (default 1).
        max_error_rate:    Non-A fraction tolerance for poly-A scan (0.0=strict).
        max_consecutive_non_a: Stop poly-A scan at ≥(this+1) consecutive non-A.
        adapter_window:    Bases from each end to inspect (default 200).
        trim_5p_polyt:     Also trim poly-T from the 5' end (default False).
        min_polya:         Minimum poly-A length to trim (default 1).
        min_polyt:         Minimum poly-T length to trim (default 1; requires
                           ``trim_5p_polyt=True``).

    Returns:
        Stats dict: total, trimmed_3p, trimmed_5p, untrimmed, pass_counts {0,1,2}.
    """
    stats: Dict = {
        'total':       0,
        'trimmed_3p':  0,
        'trimmed_5p':  0,
        'untrimmed':   0,
        'pass_counts': {0: 0, 1: 0, 2: 0},
        # Read-orientation calls from tail evidence: S=sense (3' poly-A),
        # A=antisense (5' poly-T), B=both.  Unresolved reads are not counted.
        'orientation': {'S': 0, 'A': 0, 'B': 0},
    }
    metadata_rows: List[Dict] = []

    # Derive seeds from the full adapter sequences.  The 3' seed is the
    # adapter's FIRST seed_len bases and the 5' seed its LAST seed_len bases:
    # in both cases that is the end of the adapter abutting the homopolymer,
    # which is the boundary we need to locate.
    adapter_3p_seed = adapter_3p[:seed_len] if adapter_3p else ''
    adapter_5p_seed = adapter_5p[-seed_len:] if adapter_5p else ''

    # Open output FASTQ
    out_path = str(output_fastq_path)
    if out_path.endswith('.gz') or out_path.endswith('.gzip'):
        out_fh = gzip.open(out_path, 'wt')
    else:
        out_fh = open(out_path, 'w')

    try:
        for header, seq, quals in _iter_fastq(Path(input_fastq_path)):
            stats['total'] += 1

            # Bare UUID: everything before the first whitespace in the header
            bare_id = header.split()[0]

            orig_len = len(seq)
            trim_3p_polya  = 0
            trim_3p_pass   = 0

            # ── 3' end: poly-A + adapter detection ───────────────────────────
            polya_len, adapter_stub, _last_base, adapter_pass = \
                find_cdna_3p_polya_and_adapter(
                    seq,
                    adapter_seed          = adapter_3p_seed,
                    adapter_full          = adapter_3p,
                    seed_max_mismatches   = seed_max_mismatches,
                    max_error_rate        = max_error_rate,
                    max_consecutive_non_a = max_consecutive_non_a,
                    adapter_window        = adapter_window,
                    min_polya             = min_polya,
                )

            stats['pass_counts'][adapter_pass] += 1
            total_trim_3p = polya_len + len(adapter_stub)

            # Captured sequences/quals for parquet (empty strings if nothing trimmed)
            saved_polya_seq    = ''
            saved_polya_quals  = ''
            saved_adapt_seq    = ''
            saved_adapt_quals  = ''
            saved_5p_seq       = ''
            saved_5p_quals     = ''
            saved_5p_prefix_seq   = ''
            saved_5p_prefix_quals = ''

            if polya_len >= 1 and total_trim_3p < len(seq):
                trim_3p_polya  = polya_len
                trim_3p_pass   = adapter_pass
                # Split the removed tail: [poly-A portion][adapter portion]
                # The poly-A portion is the rightmost polya_len bases of the
                # pre-trimmed sequence; the adapter portion is to the right of that.
                pre_trim_end   = len(seq)
                adapt_start    = pre_trim_end - len(adapter_stub)
                polya_start    = adapt_start - polya_len
                saved_polya_seq   = seq[polya_start:adapt_start]
                saved_polya_quals = quals[polya_start:adapt_start]
                saved_adapt_seq   = seq[adapt_start:]
                saved_adapt_quals = quals[adapt_start:]
                seq   = seq[:polya_start]
                quals = quals[:polya_start]
                stats['trimmed_3p'] += 1

            # ── 5' end: poly-T detection (antisense reads) ───────────────────
            # The poly-T sits BEHIND the 5' adapter/barcode/UMI (~120-140 nt on
            # SQK-PCB114), so everything from the read start through the end of
            # the poly-T run is removed.  The prefix 5' of the run is kept in
            # the sidecar (it carries the UMI) rather than discarded.
            n_polyt = 0
            polyt_start = -1
            if trim_5p_polyt and seq:
                n_polyt, polyt_start = find_cdna_5p_polyt(
                    seq,
                    adapter_seed          = adapter_5p_seed,
                    seed_max_mismatches   = seed_max_mismatches,
                    max_error_rate        = max_error_rate,
                    max_consecutive_non_t = max_consecutive_non_a,  # symmetric
                    adapter_window        = adapter_window,
                    min_polyt             = min_polyt,
                )
                cut = polyt_start + n_polyt if n_polyt >= 1 else 0
                if n_polyt >= 1 and cut < len(seq):
                    saved_5p_prefix_seq   = seq[:polyt_start]
                    saved_5p_prefix_quals = quals[:polyt_start]
                    saved_5p_seq   = seq[polyt_start:cut]
                    saved_5p_quals = quals[polyt_start:cut]
                    seq   = seq[cut:]
                    quals = quals[cut:]
                    stats['trimmed_5p'] += 1
                else:
                    n_polyt, polyt_start = 0, -1

            if trim_3p_polya == 0 and (not trim_5p_polyt or not saved_5p_seq):
                stats['untrimmed'] += 1

            # Read orientation, from the tail evidence.  A 3' poly-A means the
            # read is the mRNA sense strand; a 5' poly-T means it is the reverse
            # complement.  Emitted as a SAM-style FASTQ comment so that
            # `minimap2 -y` carries it into the BAM as `ro:A:S|A|U`, where
            # `rectify correct --ONT-cDNA` consumes it to resolve RNA strand
            # per read instead of applying the (DRS) fixed rule.
            _has_a = trim_3p_polya  >= _MIN_ORIENTATION_TAIL
            _has_t = len(saved_5p_seq) >= _MIN_ORIENTATION_TAIL
            if _has_a and _has_t:
                orientation = 'B'          # both tails — resolved as sense downstream
            elif _has_t:
                orientation = 'A'          # antisense
            elif _has_a:
                orientation = 'S'          # sense
            else:
                orientation = 'U'          # unresolved
            if orientation != 'U':
                stats['orientation'][orientation] += 1

            # Write trimmed read (bare UUID as read name + the orientation tag)
            out_fh.write(f'@{bare_id}\tro:A:{orientation}\n{seq}\n+\n{quals}\n')

            metadata_rows.append({
                'read_id':            bare_id,
                'original_len':       orig_len,
                'trimmed_len':        len(seq),
                # 3' poly-A tail
                'polya_3p_len':       trim_3p_polya,
                'polya_3p_seq':       saved_polya_seq,
                'polya_3p_quals':     saved_polya_quals,
                # 3' adapter stub (after poly-A)
                'adapter_3p_len':     len(saved_adapt_seq),
                'adapter_3p_seq':     saved_adapt_seq,
                'adapter_3p_quals':   saved_adapt_quals,
                'adapter_3p_pass':    trim_3p_pass,
                # 5' poly-T (antisense reads)
                'polyt_5p_len':       len(saved_5p_seq),
                'polyt_5p_seq':       saved_5p_seq,
                'polyt_5p_quals':     saved_5p_quals,
                # 5' prefix removed ahead of the poly-T (adapter/barcode/UMI)
                'polyt_5p_offset':    polyt_start,
                'prefix_5p_len':      len(saved_5p_prefix_seq),
                'prefix_5p_seq':      saved_5p_prefix_seq,
                'prefix_5p_quals':    saved_5p_prefix_quals,
                # Read orientation label: S=sense, A=antisense, B=both, U=unresolved
                'orientation':        orientation,
            })

    finally:
        out_fh.close()

    _write_metadata(metadata_rows, Path(metadata_path))

    logger.info(
        "trim_cdna_fastq_polya: total=%d trimmed_3p=%d trimmed_5p=%d untrimmed=%d "
        "pass0=%d pass1=%d pass2=%d | orientation sense=%d antisense=%d both=%d "
        "unresolved=%d",
        stats['total'], stats['trimmed_3p'], stats['trimmed_5p'], stats['untrimmed'],
        stats['pass_counts'][0], stats['pass_counts'][1], stats['pass_counts'][2],
        stats['orientation']['S'], stats['orientation']['A'], stats['orientation']['B'],
        stats['total'] - sum(stats['orientation'].values()),
    )
    return stats


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def run(args) -> int:
    """Entry point called from rectify/cli.py for ``rectify trim-cdna-polya``."""
    from datetime import datetime as _dt, timezone as _tz
    from time import perf_counter as _perf
    _stage_started_at = _dt.now(_tz.utc).isoformat()
    _t_start = _perf()
    input_fastq = Path(args.input_fastq)
    if not input_fastq.exists():
        print(f"ERROR: Input FASTQ not found: {input_fastq}", file=sys.stderr)
        return 1

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_stem   = getattr(args, 'prefix', None) or input_fastq.name
    # Strip common FASTQ extensions to get a clean stem
    for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
        if sample_stem.endswith(ext):
            sample_stem = sample_stem[:-len(ext)]
            break

    # Output FASTQ preserves the input compression format
    if str(input_fastq).endswith('.gz'):
        out_fastq = output_dir / f"{sample_stem}.polya_trimmed.fastq.gz"
    else:
        out_fastq = output_dir / f"{sample_stem}.polya_trimmed.fastq"

    use_tsv     = getattr(args, 'tsv', False)
    meta_ext    = '.tsv' if use_tsv else '.parquet'
    meta_path   = output_dir / f"{sample_stem}_cdna_polya_trim_metadata{meta_ext}"

    adapter_3p  = getattr(args, 'adapter_3p',  _CRTA)
    adapter_5p  = getattr(args, 'adapter_5p',  _CRTA_RC)
    seed_len    = getattr(args, 'seed_len',    _SEED_LEN)
    seed_mm     = getattr(args, 'seed_max_mismatches', _SEED_MAX_MISMATCHES)

    print(f"Input FASTQ:     {input_fastq}")
    print(f"Output FASTQ:    {out_fastq}")
    print(f"Metadata:        {meta_path}")
    print(f"3' adapter seed: {adapter_3p[:seed_len]!r} ({seed_len} bp, ≤{seed_mm} mm)")
    if adapter_5p:
        print(f"5' adapter seed: {adapter_5p[:seed_len]!r} ({seed_len} bp)")
    print(f"max_error_rate:  {args.max_error_rate}")
    print(f"trim_5p_polyt:   {args.trim_5p_polyt}")

    stats = trim_cdna_fastq_polya(
        input_fastq_path      = str(input_fastq),
        output_fastq_path     = str(out_fastq),
        metadata_path         = str(meta_path),
        adapter_3p            = adapter_3p,
        adapter_5p            = adapter_5p if adapter_5p else None,
        seed_len              = seed_len,
        seed_max_mismatches   = seed_mm,
        max_error_rate        = args.max_error_rate,
        max_consecutive_non_a = args.max_consecutive_non_a,
        adapter_window        = args.adapter_window,
        trim_5p_polyt         = args.trim_5p_polyt,
        min_polya             = args.min_polya,
        min_polyt             = getattr(args, 'min_polyt', 1),
    )

    total   = stats['total']
    pct_3p  = 100 * stats['trimmed_3p'] / total if total else 0.0
    pct_5p  = 100 * stats['trimmed_5p'] / total if total else 0.0

    print(f"\nTrim summary:")
    print(f"  Total reads:              {total:,}")
    print(f"  3' poly-A trimmed:        {stats['trimmed_3p']:,}  ({pct_3p:.1f}%)")
    print(f"  5' poly-T trimmed:        {stats['trimmed_5p']:,}  ({pct_5p:.1f}%)")
    print(f"  Untrimmed:                {stats['untrimmed']:,}")
    print(f"  Pass 0 (no adapter stub): {stats['pass_counts'][0]:,}")
    print(f"  Pass 1 (seed match):      {stats['pass_counts'][1]:,}")
    print(f"  Pass 2 (peel):            {stats['pass_counts'][2]:,}")
    print(f"\nRestoration columns in metadata parquet:")
    print(f"  polya_3p_seq / polya_3p_quals   — poly-A tail only")
    print(f"  adapter_3p_seq / adapter_3p_quals — adapter stub only")
    print(f"  polyt_5p_seq / polyt_5p_quals   — 5' poly-T (antisense reads)")

    # ---- provenance sidecar -------------------------------------------------
    # trim-cdna-polya REWRITES the reads, so its code version has to be pinned:
    # every downstream coordinate depends on what this stage removed. Both bugs
    # fixed on 2026-08-01 (the non-adapter-anchored 5' poly-T detector, and the
    # missing orientation label) lived HERE, and this stage wrote no provenance
    # at all — so the affected FASTQs carry no record of which code produced
    # them. Non-fatal on failure, per the other stages' convention.
    try:
        from rectify.core.provenance import ProvenanceRecord, write_stage_sidecar
        from rectify.utils.version import get_rectify_git_sha as _get_sha_trim
        _sc = ProvenanceRecord.from_components(
            stage='trim_cdna_polya',
            stage_subtype='ont_cdna',
            sample_id=sample_stem,
            sample_output_dir=output_dir,
            started_at=_stage_started_at,
            completed_at=_dt.now(_tz.utc).isoformat(),
            exit_status=0,
            inputs={'fastq': input_fastq},
            outputs={'trimmed_fastq': out_fastq, 'trim_metadata': meta_path},
            stats={
                'wall_seconds': _perf() - _t_start,
                'total': stats['total'],
                'trimmed_3p': stats['trimmed_3p'],
                'trimmed_5p': stats['trimmed_5p'],
                'untrimmed': stats['untrimmed'],
                'pass0': stats['pass_counts'][0],
                'pass1': stats['pass_counts'][1],
                'pass2': stats['pass_counts'][2],
                # Orientation calls -- the signal `correct --ONT-cDNA` consumes.
                'orientation_sense': stats.get('orientation', {}).get('S', 0),
                'orientation_antisense': stats.get('orientation', {}).get('A', 0),
                'orientation_both': stats.get('orientation', {}).get('B', 0),
            },
            argv=sys.argv,
            rectify_git_sha=_get_sha_trim(),
        )
        write_stage_sidecar(_sc, sample_output=output_dir)
    except Exception as _sc_exc:
        logger.warning("Failed to write trim_cdna_polya sidecar: %s", _sc_exc)

    return 0


def create_trim_cdna_polya_parser(subparsers):
    """Wire the `trim-cdna-polya` subcommand into the given subparsers group."""
    import argparse
    # =========================================================================
    # trim-cdna-polya command (cDNA FASTQ poly(A)+adapter pre-trimming)
    # =========================================================================
    trim_cdna_parser = subparsers.add_parser(
        'trim-cdna-polya',
        help='Trim poly(A) tails + adapter stubs from cDNA Nanopore FASTQ files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            'Pre-trim poly(A) tails and CRTA adapter stubs from raw cDNA Nanopore '
            'FASTQ files before alignment.  Unlike rectify trim-polya (which requires '
            'Dorado-aligned BAMs with pt:i: tags), this module operates on plain FASTQ '
            'files — suitable when only FASTQ data is available.\n\n'
            'Default 3\' adapter: CRTA  (CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG)\n'
            'Default 5\' adapter: CRTA_RC (CTTACCTGCGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG)\n\n'
            'Adapter detection uses a configurable k-mer seed (default: first 12 bp of the\n'
            'adapter, ≤1 mismatch allowed) — more robust than the short stub regex used for\n'
            'DRS data.'
        ),
    )
    trim_cdna_parser.add_argument(
        'input_fastq',
        help='Input FASTQ file (plain or gzip; .fastq, .fastq.gz, .fq, .fq.gz)',
    )
    trim_cdna_parser.add_argument(
        '-o', '--output-dir',
        dest='output_dir',
        required=True,
        help='Output directory for trimmed FASTQ and trim metadata',
    )
    trim_cdna_parser.add_argument(
        '--adapter-3p',
        dest='adapter_3p',
        default='CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG',
        metavar='SEQ',
        help='Full 3\' adapter sequence expected after the poly(A) tail on RNA-sense reads '
             '(default: CRTA template-switching adapter)',
    )
    trim_cdna_parser.add_argument(
        '--adapter-5p',
        dest='adapter_5p',
        default='CTTACCTGCGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG',
        metavar='SEQ',
        help='Full 5\' adapter sequence expected before the poly(T) run on antisense reads '
             '(default: CRTA_RC; used only for seed derivation). Pass "" to disable.',
    )
    trim_cdna_parser.add_argument(
        '--seed-len',
        type=int,
        default=12,
        dest='seed_len',
        metavar='INT',
        help='k-mer seed length for adapter detection (default: 12 bp from adapter 5\' end)',
    )
    trim_cdna_parser.add_argument(
        '--seed-max-mismatches',
        type=int,
        default=1,
        dest='seed_max_mismatches',
        metavar='INT',
        help='Mismatches allowed in the adapter seed match (default: 1)',
    )
    trim_cdna_parser.add_argument(
        '--adapter-window',
        type=int,
        default=200,
        dest='adapter_window',
        metavar='INT',
        help='Bases from each end to inspect for poly-A/T + adapter (default: 200)',
    )
    trim_cdna_parser.add_argument(
        '--max-error-rate',
        type=float,
        default=0.0,
        dest='max_error_rate',
        help='Max cumulative non-A fraction for poly-A scan (0.0 = strict pure-A)',
    )
    trim_cdna_parser.add_argument(
        '--max-consecutive-non-a',
        type=int,
        default=1,
        dest='max_consecutive_non_a',
        help='Stop poly-A scan after this many consecutive non-A bases (default 1 = stop at ≥2)',
    )
    trim_cdna_parser.add_argument(
        '--trim-5p-polyt',
        action='store_true',
        default=False,
        dest='trim_5p_polyt',
        help='Also trim poly-T runs from the 5\' end (antisense reads starting with poly-T)',
    )
    trim_cdna_parser.add_argument(
        '--min-polya',
        type=int,
        default=1,
        dest='min_polya',
        metavar='INT',
        help='Minimum 3\' poly-A run length to trim (default: 1)',
    )
    trim_cdna_parser.add_argument(
        '--min-polyt',
        type=int,
        default=1,
        dest='min_polyt',
        metavar='INT',
        help='Minimum 5\' poly-T run length to trim (requires --trim-5p-polyt; default: 1)',
    )
    trim_cdna_parser.add_argument(
        '--prefix',
        default=None,
        help='Output file prefix (default: derived from input FASTQ filename)',
    )
    trim_cdna_parser.add_argument(
        '--tsv',
        action='store_true',
        default=False,
        help='Write trim metadata as TSV instead of parquet',
    )
