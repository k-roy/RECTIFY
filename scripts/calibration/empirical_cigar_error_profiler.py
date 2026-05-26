#!/usr/bin/env python3
"""
Empirical CIGAR Error Profiler for Nanopore Reads

Derives empirical per-base error rates (deletion, insertion, substitution) from
regions where multiple independent aligners agree on the mismatch/indel profile.

Strategy
--------
For each read present in ALL N aligners:
  1. Extract exonic (non-N-op) reference intervals from each aligner's CIGAR.
  2. Intersect those intervals → positions that every aligner judges as exonic.
  3. At each agreed position, check whether all aligners call the same CIGAR op
     (match, mismatch/substitution, deletion) and the same reference base.
  4. Look up the homopolymer run length at that reference position.
  5. Accumulate counts: {(op_type, ref_base, hp_run_length): count}.
     When --strand-split is given, also accumulate a parallel strand_counts dict
     with keys (op_type, ref_base, hp_run_length, strand) where strand is '+'/'-'.

After processing all reads/chunks:
  - Compute per-op error rates as a function of HP run length.
  - Normalize to sub-rate at HP=1 to derive penalty scores comparable to the
    existing heuristic table in junction_refiner._hp_edit_distance().
  - Write TSVs and a summary plot.

Usage (quick test — chunk 0, first 5 000 reads)
-------
    python empirical_cigar_error_profiler.py \\
        --run-dir /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/dev_runs/wt_by4742_rep1_chunked_20260412 \\
        --reference /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa \\
        --output-dir ./error_profile_wt_by4742_rep1 \\
        --chunks 0 \\
        --max-reads 5000

Usage (full run, 8 parallel workers — DRS)
-------------------------------------------
    python empirical_cigar_error_profiler.py \\
        --run-dir /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/dev_runs/wt_by4742_rep1_chunked_20260412 \\
        --reference /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa \\
        --output-dir ./error_profile_wt_by4742_rep1 \\
        --isolation-flank 10 --union --str-repeat --workers 8

Protocol-aware runs (cDNA / QSrev)
-----------------------------------
    # ONT PCR-cDNA: 5-aligner long-read panel, mask 50 bp of both read ends
    python empirical_cigar_error_profiler.py \\
        --protocol cdna \\
        --run-dir /path/to/wt_cdna_chunked_run \\
        --reference /path/to/genome.fasta \\
        --output-dir ./error_profile_cdna \\
        --isolation-flank 10 --union --str-repeat --workers 8

    # QuantSeq REV (Illumina short-read): bbmap + bwa panel, mask 30 bp 5' + 50 bp 3'
    python empirical_cigar_error_profiler.py \\
        --protocol qsrev \\
        --run-dir /path/to/wt_qsrev_chunked_run \\
        --reference /path/to/genome.fasta \\
        --output-dir ./error_profile_qsrev \\
        --isolation-flank 10 --workers 8

Or with explicit per-aligner BAMs (any dataset):
    python empirical_cigar_error_profiler.py \\
        --aligner-bams minimap2:/path/to/mm2.bam mapPacBio:/path/to/mp.bam \\
        --reference /path/to/genome.fasta \\
        --output-dir ./error_profile

Author: Kevin R. Roy
Email:  kevinrjroy@gmail.com
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pysam
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-7s  %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DEFAULT_ALIGNERS = ['minimap2', 'mapPacBio', 'deSALT', 'uLTRA', 'gapmm2']

# Per-protocol aligner panels (gated by --protocol).
# DRS / cDNA share the 5-aligner long-read panel; QSrev uses short-read aligners.
_DEFAULT_ALIGNERS_BY_PROTOCOL: Dict[str, List[str]] = {
    'drs':   _DEFAULT_ALIGNERS,
    'cdna':  _DEFAULT_ALIGNERS,
    'qsrev': ['bbmap', 'bwa'],
}

# Per-protocol read-coordinate exclusion defaults (5'_bp, 3'_bp), used when the
# user does not pass --exclude-5prime-bp / --exclude-3prime-bp explicitly.
# Note: DRS defaults to (0, 0) — not (0, 50) as suggested by an earlier draft of
# the cDNA/QSrev handoff — to preserve byte-identical reproducibility of the
# bundled DRS tables (calibrated 2026-04-22). Pass --exclude-3prime-bp 50 to
# regenerate DRS tables with the same poly-A-boundary mask used for cDNA.
_PROTOCOL_EXCLUSION_DEFAULTS: Dict[str, Tuple[int, int]] = {
    'drs':   (0, 0),
    'cdna':  (50, 50),
    # QSrev: 76-bp Han 2023 reads forced a smaller budget than the handoff's
    # (30, 50). (10, 16) keeps ~50 bp interior while still masking the dT-primer
    # adapter remnant (5') and the residual A-tract walkback near the 3' end.
    # Approved by Kevin 2026-05-17.
    'qsrev': (10, 16),
}

# CIGAR op codes (pysam integers)
_OP_M   = 0   # alignment match (M) — may be match OR mismatch
_OP_I   = 1   # insertion to reference
_OP_D   = 2   # deletion from reference
_OP_N   = 3   # intron skip
_OP_S   = 4   # soft clip
_OP_H   = 5   # hard clip
_OP_EQ  = 7   # sequence match (=)
_OP_X   = 8   # sequence mismatch (X)

# Sets of ops by what they consume
_REF_CONSUMING   = frozenset([_OP_M, _OP_D, _OP_N, _OP_EQ, _OP_X])
_QUERY_CONSUMING = frozenset([_OP_M, _OP_I, _OP_S, _OP_EQ, _OP_X])

# Ops that move through reference without producing an error event.
# (N moves ref_pos but we record it as an intron boundary, not a CIGAR op.)
_EXON_REF_CONSUMING = frozenset([_OP_M, _OP_D, _OP_EQ, _OP_X])

# Label strings for output tables / plots
_OP_LABELS = {'M': 'match', 'D': 'deletion', 'X': 'substitution', 'I': 'insertion'}
_OP_COLORS = {'deletion': '#e6194b', 'substitution': '#4363d8', 'insertion': '#f58231'}

# Op-code encoding for the vectorised ops_matrix in profile_chunk
# uint8 values: M=0, X=1, D=2, absent=255
_OP_CODE = {'M': np.uint8(0), 'X': np.uint8(1), 'D': np.uint8(2)}
_CODE_OP = {0: 'M', 1: 'X', 2: 'D'}


# ---------------------------------------------------------------------------
# UMI cluster-size binning
# ---------------------------------------------------------------------------

# A bin is a tuple (label, kind, n) where:
#   kind == 'eq' → matches reads with UMI cluster size == n
#   kind == 'ge' → matches reads with UMI cluster size >= n  (e.g. '3+')
# Plain data representation (no lambdas) so the spec pickles cleanly across
# multiprocessing fork/spawn boundaries.
UmiBinSpec = Tuple[str, str, int]


def _sanitize_bin_label(label: str) -> str:
    """Filename-safe form of a bin label: '3+' → '3plus', '1' → '1'."""
    return label.replace('+', 'plus')


def parse_umi_bins(spec: str) -> List[UmiBinSpec]:
    """Parse a CLI bin spec like '1,2,3+' into a list of UmiBinSpec tuples.

    Each token is either an integer (interpreted as eq, exact match) or an
    integer followed by '+' (interpreted as ge, ≥). Whitespace is tolerated.

    Returns a list of (label, kind, n) tuples; labels preserve the original
    token (e.g. '3+' not '3plus') — sanitize with _sanitize_bin_label() at
    write time. Raises ValueError on malformed input.
    """
    bins: List[UmiBinSpec] = []
    if not spec or not spec.strip():
        raise ValueError('--umi-bins must contain at least one bin token')
    for raw in spec.split(','):
        tok = raw.strip()
        if not tok:
            raise ValueError(f'empty bin token in --umi-bins {spec!r}')
        if tok.endswith('+'):
            try:
                n = int(tok[:-1])
            except ValueError as exc:
                raise ValueError(f'malformed bin token {tok!r}: expected <int>+') from exc
            bins.append((tok, 'ge', n))
        else:
            try:
                n = int(tok)
            except ValueError as exc:
                raise ValueError(f'malformed bin token {tok!r}: expected <int> or <int>+') from exc
            bins.append((tok, 'eq', n))
    return bins


def assign_bin(umi_size: Optional[int], bins: List[UmiBinSpec]) -> Optional[str]:
    """Return the label of the first matching bin, or None for no match.

    None inputs (missing tag) always return None — they are pooled-only.
    Bins are checked in order; the first matching bin wins. The caller is
    responsible for ensuring bins are non-overlapping if exact ownership
    matters (1, 2, 3+ is non-overlapping by construction).
    """
    if umi_size is None:
        return None
    for label, kind, n in bins:
        if kind == 'eq' and umi_size == n:
            return label
        if kind == 'ge' and umi_size >= n:
            return label
    return None


# ---------------------------------------------------------------------------
# Homopolymer utilities  (ported from junction_refiner._hp_run_length)
# ---------------------------------------------------------------------------

def hp_run_length(seq: str, pos: int) -> int:
    """Return the length of the homopolymer run that contains *pos*."""
    if pos < 0 or pos >= len(seq):
        return 1
    c = seq[pos]
    left = pos
    while left > 0 and seq[left - 1] == c:
        left -= 1
    right = pos
    while right < len(seq) - 1 and seq[right + 1] == c:
        right += 1
    return right - left + 1


def str_repeat_info(
    ref_seq: str,
    pos: int,
    max_unit: int = 3,
    min_copies: int = 3,
) -> Tuple[Optional[str], int]:
    """Return (canonical_unit, n_copies) if *pos* is within a short tandem repeat.

    Only considers repeat units of length 2–max_unit with ≥ 2 distinct
    characters (mononucleotide repeats are handled by hp_run_length and
    excluded here).  The canonical unit is the lexicographically smallest
    rotation of the raw unit (e.g. TATATA → canonical 'AT', not 'TA').

    Returns (None, 0) if the position is not in a qualifying STR context.

    Algorithm:
        For each unit length L in [2, max_unit], try all L phases that could
        place *pos* within a repeat unit.  Walk left to find the leftmost
        occurrence, count copies rightward, verify *pos* falls within the run.
        Keep the best (most copies) result across all unit lengths and phases.
    """
    n = len(ref_seq)
    best_unit: Optional[str] = None
    best_copies = 0

    for unit_len in range(2, max_unit + 1):
        for phase in range(unit_len):
            unit_start = pos - phase
            if unit_start < 0 or unit_start + unit_len > n:
                continue
            unit = ref_seq[unit_start:unit_start + unit_len]
            if len(set(unit)) == 1:          # mononucleotide — skip
                continue
            # Walk left to find leftmost copy of this exact unit
            left = unit_start
            while left >= unit_len and ref_seq[left - unit_len:left] == unit:
                left -= unit_len
            # Count all copies from left
            copies, p = 0, left
            while p + unit_len <= n and ref_seq[p:p + unit_len] == unit:
                copies += 1
                p += unit_len
            # Verify pos is inside this repeat run
            if not (left <= pos < left + copies * unit_len):
                continue
            if copies >= min_copies and copies > best_copies:
                best_copies = copies
                best_unit = min(unit[i:] + unit[:i] for i in range(unit_len))

    return (best_unit, best_copies) if best_unit is not None else (None, 0)


# ---------------------------------------------------------------------------
# Precomputed lookup tables for HP and STR context (replaces per-call scanning)
# ---------------------------------------------------------------------------

def _precompute_hp_array(seq: str) -> np.ndarray:
    """Build uint16 hp_array[i] = HP run length containing position i.

    Single O(n) numpy pass — replaces repeated hp_run_length() calls during
    profile_chunk inner loops.  For a 12 Mb yeast genome this takes <10 ms.
    """
    n = len(seq)
    if n == 0:
        return np.zeros(0, dtype=np.uint16)
    seq_arr = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    # Mark positions where the base changes (= start of a new run)
    change = np.empty(n, dtype=bool)
    change[0] = True
    change[1:] = seq_arr[1:] != seq_arr[:-1]
    # run_id[i] = 1-based index of the run containing position i
    run_id = np.cumsum(change)
    # Count positions per run, then map back to each position
    _, run_lengths = np.unique(run_id, return_counts=True)
    return run_lengths[run_id - 1].astype(np.uint16)


def _precompute_str_map(
    ref_seq: str,
    max_unit: int = 3,
    min_copies: int = 3,
) -> Dict[int, Tuple[str, int]]:
    """Build {pos: (canonical_unit, n_copies)} for all STR-context positions.

    Single O(n × max_unit) sweep — replaces per-position str_repeat_info() calls.
    Mononucleotide runs are excluded (those positions are handled by hp_run_length /
    _precompute_hp_array).  When multiple overlapping STRs claim the same position,
    the entry with the highest copy count wins (same tie-breaking as str_repeat_info).
    """
    n = len(ref_seq)
    result: Dict[int, Tuple[str, int]] = {}
    for unit_len in range(2, max_unit + 1):
        pos = 0
        while pos <= n - unit_len * min_copies:
            unit = ref_seq[pos:pos + unit_len]
            if len(set(unit)) == 1:          # mononucleotide repeat — skip
                pos += 1
                continue
            # Count consecutive copies of this exact unit starting at pos
            copies = 0
            end = pos
            while end + unit_len <= n and ref_seq[end:end + unit_len] == unit:
                copies += 1
                end += unit_len
            if copies >= min_copies:
                canon = min(unit[i:] + unit[:i] for i in range(unit_len))
                for p in range(pos, end):
                    if p not in result or copies > result[p][1]:
                        result[p] = (canon, copies)
                pos = end     # jump past this repeat; sub-repeats already covered
            else:
                pos += 1
    return result


# ---------------------------------------------------------------------------
# Exonic interval extraction & intersection
# ---------------------------------------------------------------------------

def get_exonic_ranges(read: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    """Return sorted list of (ref_start, ref_end) for non-N exonic regions."""
    if not read.cigartuples:
        return []

    ranges: List[Tuple[int, int]] = []
    ref_pos = read.reference_start
    seg_start: Optional[int] = None

    for op, length in read.cigartuples:
        if op == _OP_N:
            if seg_start is not None:
                ranges.append((seg_start, ref_pos))
                seg_start = None
            ref_pos += length
        elif op in _REF_CONSUMING:
            if seg_start is None:
                seg_start = ref_pos
            ref_pos += length
        # S, H, I: don't consume reference → no effect on exonic range

    if seg_start is not None:
        ranges.append((seg_start, ref_pos))

    return ranges


def get_exonic_ranges_inner(
    read: pysam.AlignedSegment,
    exclude_5p_bp: int,
    exclude_3p_bp: int,
) -> List[Tuple[int, int]]:
    """Like :func:`get_exonic_ranges` but trims to read-coord interior.

    Returns ref intervals where the corresponding read-coordinate (in BAM
    orientation, *after* soft-clip stripping) falls inside the inner window
    ``[qa_start + exclude_5p_bp, qa_end - exclude_3p_bp)``.

    Used to mask out protocol-specific end-of-read noise (e.g. template-switch
    adapter bleed for cDNA, CSP A-tract walkback for QSrev) before consensus
    error scoring. Deletion ops anchor to the preceding query position; insert
    and clip ops don't contribute ref intervals.

    Falls back to :func:`get_exonic_ranges` when both exclusion budgets are 0.
    """
    if not read.cigartuples:
        return []
    if exclude_5p_bp == 0 and exclude_3p_bp == 0:
        return get_exonic_ranges(read)

    qa_start = read.query_alignment_start
    qa_end = read.query_alignment_end  # half-open
    inner_q_lo = qa_start + exclude_5p_bp
    inner_q_hi = qa_end - exclude_3p_bp
    if inner_q_hi <= inner_q_lo:
        return []

    ranges: List[Tuple[int, int]] = []
    ref_pos = read.reference_start
    qry_pos = 0  # position in query_sequence (includes soft clips)
    seg_start: Optional[int] = None

    for op, length in read.cigartuples:
        if op == _OP_S:
            qry_pos += length
        elif op == _OP_H:
            pass  # hard clip: not in query_sequence
        elif op == _OP_I:
            qry_pos += length  # consumes query only
        elif op == _OP_N:
            if seg_start is not None:
                ranges.append((seg_start, ref_pos))
                seg_start = None
            ref_pos += length
        elif op == _OP_D:
            # Anchor deletion to the last seen query position
            anchor_q = qry_pos - 1 if qry_pos > 0 else 0
            if inner_q_lo <= anchor_q < inner_q_hi:
                if seg_start is None:
                    seg_start = ref_pos
                ref_pos += length
            else:
                if seg_start is not None:
                    ranges.append((seg_start, ref_pos))
                    seg_start = None
                ref_pos += length
        elif op == _OP_M or op == _OP_EQ or op == _OP_X:
            # Sub-range of this op that falls inside the inner zone
            q_op_lo, q_op_hi = qry_pos, qry_pos + length
            sub_lo = max(q_op_lo, inner_q_lo)
            sub_hi = min(q_op_hi, inner_q_hi)
            if sub_lo >= sub_hi:
                # Entire op outside inner zone
                if seg_start is not None:
                    ranges.append((seg_start, ref_pos))
                    seg_start = None
            else:
                live_ref_lo = ref_pos + (sub_lo - q_op_lo)
                live_ref_hi = ref_pos + (sub_hi - q_op_lo)
                if sub_lo > q_op_lo and seg_start is not None:
                    # Leading dead chunk; flush prior segment
                    ranges.append((seg_start, ref_pos))
                    seg_start = None
                if seg_start is None:
                    seg_start = live_ref_lo
                if sub_hi < q_op_hi:
                    # Trailing dead chunk; close segment here
                    ranges.append((seg_start, live_ref_hi))
                    seg_start = None
            ref_pos += length
            qry_pos += length

    if seg_start is not None:
        ranges.append((seg_start, ref_pos))

    return ranges


def _pairwise_intersect(
    a: List[Tuple[int, int]],
    b: List[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """Intersect two sorted, non-overlapping interval lists."""
    result: List[Tuple[int, int]] = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo < hi:
            result.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return result


def intersect_all(
    lists: List[List[Tuple[int, int]]],
) -> List[Tuple[int, int]]:
    """Intersect an arbitrary number of sorted interval lists."""
    if not lists:
        return []
    result = lists[0]
    for other in lists[1:]:
        result = _pairwise_intersect(result, other)
        if not result:
            return []
    return result


# ---------------------------------------------------------------------------
# Per-position CIGAR op map
# ---------------------------------------------------------------------------

def build_pos_op_map(
    read: pysam.AlignedSegment,
    ref_seq: str,
) -> Tuple[Dict[int, Tuple[str, str]], Dict[int, int]]:
    """Build per-reference-position op maps from a single alignment.

    Handles two CIGAR styles:
      - Sequence-bearing reads (minimap2, mapPacBio, deSALT, uLTRA):
          Uses get_aligned_pairs() + reference FASTA to distinguish M vs X.
      - No-sequence reads (gapmm2): CIGAR must use explicit = / X ops.
          Falls back to CIGAR-walk using reference only; ambiguous M ops skipped.

    Returns
    -------
    pos_ops   : {ref_pos: (op_char, ref_base)}
                op_char is 'M' (match), 'X' (mismatch/sub), or 'D' (deletion).
    ins_before: {ref_pos: inserted_base_count}
                Number of inserted query bases immediately before ref_pos.
    """
    if read.query_sequence:
        return _build_pos_op_map_with_seq(read, ref_seq)
    else:
        return _build_pos_op_map_from_cigar(read, ref_seq)


def _build_pos_op_map_with_seq(
    read: pysam.AlignedSegment,
    ref_seq: str,
) -> Tuple[Dict[int, Tuple[str, str]], Dict[int, int]]:
    """Op map using get_aligned_pairs() — requires query_sequence."""
    pos_ops: Dict[int, Tuple[str, str]] = {}
    ins_before: Dict[int, int] = {}
    pending_ins = 0
    query_seq = read.query_sequence

    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos is None:
            pending_ins += 1
        elif qpos is None:
            # Deletion
            if pending_ins:
                ins_before[rpos] = ins_before.get(rpos, 0) + pending_ins
                pending_ins = 0
            if 0 <= rpos < len(ref_seq):
                pos_ops[rpos] = ('D', ref_seq[rpos])
        else:
            # Match or mismatch
            if pending_ins:
                ins_before[rpos] = ins_before.get(rpos, 0) + pending_ins
                pending_ins = 0
            if 0 <= rpos < len(ref_seq) and qpos < len(query_seq):
                ref_base = ref_seq[rpos]
                read_base = query_seq[qpos].upper()
                # SAM spec allows '=' in SEQ to mean "matches reference at this
                # position" — samtools/minimap2 sometimes emit this as a size
                # compression. pysam returns it verbatim. Without this branch,
                # every '=' position is mis-classified as a mismatch ('X') and
                # `m_mask` ends up empty for the entire chunk — which silently
                # kills all M counts (and therefore all counts, since the
                # denominator uses M + D + X). Verified on cDNA stage1_consensus
                # BAMs 2026-05-17.
                if read_base == '=':
                    op = 'M'
                else:
                    op = 'M' if read_base == ref_base else 'X'
                pos_ops[rpos] = (op, ref_base)

    return pos_ops, ins_before


def _build_pos_op_map_from_cigar(
    read: pysam.AlignedSegment,
    ref_seq: str,
) -> Tuple[Dict[int, Tuple[str, str]], Dict[int, int]]:
    """Op map from CIGAR tuples alone — for reads without stored query sequence.

    Requires the aligner to use explicit = (EQ) and X ops; ambiguous M ops
    are skipped because without the query sequence we cannot resolve them.
    gapmm2 always writes = / X, so this is safe for gapmm2 BAMs.
    """
    pos_ops: Dict[int, Tuple[str, str]] = {}
    ins_before: Dict[int, int] = {}
    ref_pos = read.reference_start
    pending_ins = 0

    for op, length in read.cigartuples:
        if op == _OP_N:
            ref_pos += length
        elif op == _OP_D:
            for _ in range(length):
                if pending_ins:
                    ins_before[ref_pos] = ins_before.get(ref_pos, 0) + pending_ins
                    pending_ins = 0
                if 0 <= ref_pos < len(ref_seq):
                    pos_ops[ref_pos] = ('D', ref_seq[ref_pos])
                ref_pos += 1
        elif op == _OP_EQ:
            if pending_ins:
                ins_before[ref_pos] = ins_before.get(ref_pos, 0) + pending_ins
                pending_ins = 0
            for _ in range(length):
                if 0 <= ref_pos < len(ref_seq):
                    pos_ops[ref_pos] = ('M', ref_seq[ref_pos])
                ref_pos += 1
        elif op == _OP_X:
            if pending_ins:
                ins_before[ref_pos] = ins_before.get(ref_pos, 0) + pending_ins
                pending_ins = 0
            for _ in range(length):
                if 0 <= ref_pos < len(ref_seq):
                    pos_ops[ref_pos] = ('X', ref_seq[ref_pos])
                ref_pos += 1
        elif op == _OP_M:
            # Ambiguous: can't resolve match vs mismatch without query sequence.
            # Advance ref_pos but record nothing (conservative — avoids false M counts).
            ref_pos += length
        elif op == _OP_I:
            pending_ins += length
        # S, H: no effect on ref_pos

    return pos_ops, ins_before


# ---------------------------------------------------------------------------
# Isolation filter helpers
# ---------------------------------------------------------------------------

def _precompute_flank_lengths(
    pos_ops: Dict[int, Tuple[str, str]],
    ranges: List[Tuple[int, int]],
) -> Tuple[Dict[int, int], Dict[int, int]]:
    """Compute consecutive-M run lengths to the left and right of each position.

    left_m[p]  = consecutive M's ending immediately at p-1 within the same range.
    right_m[p] = consecutive M's starting immediately at p+1 within the same range.

    M-run counts do not cross range boundaries (exon-intron junctions reset the run).
    Used for the --isolation-flank filter: an error run [s, e] is "isolated" iff
    left_m[s] >= N and right_m[e] >= N in ALL seq-bearing aligners.

    Implemented with numpy cumsum to avoid per-position Python dict lookups.
    """
    left_m: Dict[int, int] = {}
    right_m: Dict[int, int] = {}
    for rng_start, rng_end in ranges:
        length = rng_end - rng_start
        if length == 0:
            continue
        positions = range(rng_start, rng_end)
        # Boolean is_M array: True where pos_ops[p] exists and op == 'M'
        is_m = np.array(
            [pos_ops.get(p, ('',))[0] == 'M' for p in positions],
            dtype=np.int32,
        )
        # Left run lengths: consecutive M's to the LEFT of each position.
        # left_m[i] = number of consecutive M's immediately before position i.
        # Reset to 0 whenever the preceding position is not M.
        # Equivalent to: cumsum of is_m, minus the cumsum at the last non-M reset.
        # Use a group-reset approach via np.maximum.accumulate on negated positions.
        lm = np.zeros(length, dtype=np.int32)
        if length > 1:
            # lm[0] = 0 always; lm[i] = lm[i-1]+1 if is_m[i-1] else 0
            lm[1:] = np.where(is_m[:-1], 0, -1)   # -1 = reset marker
            # Replace non-reset entries (0) with cumulative left M count
            # Use cumsum + reset trick: assign cumsum of is_m, then subtract
            # the running max at each reset point.
            cum = np.cumsum(is_m[:-1])              # cumulative M count (0-based, up to i-1)
            reset_vals = np.where(is_m[:-1], 0, cum)
            reset_max = np.maximum.accumulate(reset_vals)
            lm[1:] = cum - reset_max

        # Right run lengths: consecutive M's STARTING AT p+1 for each position p.
        # Computed by applying the same left-run formula to the reversed is_m array,
        # then reversing the result back.  Verified correct by symmetry with left_m.
        rm = np.zeros(length, dtype=np.int32)
        if length > 1:
            rev_is_m = is_m[::-1]
            cum_rev = np.cumsum(rev_is_m[:-1])
            reset_rev = np.where(rev_is_m[:-1], 0, cum_rev)
            reset_max_rev = np.maximum.accumulate(reset_rev)
            left_rev = np.zeros(length, dtype=np.int32)
            left_rev[1:] = cum_rev - reset_max_rev
            rm = left_rev[::-1].copy()

        # Write results to output dicts
        lm_list = lm.tolist()
        rm_list = rm.tolist()
        for i, p in enumerate(positions):
            left_m[p] = lm_list[i]
            right_m[p] = rm_list[i]
    return left_m, right_m


def _group_pos_into_runs(positions: Set[int]) -> List[Tuple[int, int]]:
    """Group a set of positions into sorted contiguous runs [(start, end), ...]."""
    if not positions:
        return []
    sorted_pos = sorted(positions)
    runs: List[Tuple[int, int]] = []
    s = prev = sorted_pos[0]
    for p in sorted_pos[1:]:
        if p == prev + 1:
            prev = p
        else:
            runs.append((s, prev))
            s = prev = p
    runs.append((s, prev))
    return runs


# ---------------------------------------------------------------------------
# BAM file discovery
# ---------------------------------------------------------------------------

def find_chunk_bams(
    run_dir: Path,
    chunk_idx: int,
    n_chunks: int,
    aligners: List[str],
) -> Dict[str, Path]:
    """Return {aligner: bam_path} for one chunk index.

    Handles both .bam and .sam outputs (deSALT occasionally writes SAM).
    """
    chunk_str = f'chunk_{chunk_idx:03d}_of_{n_chunks:03d}'
    # Derive sample name from what's in the minimap2 directory
    align_dir = run_dir / 'align'
    bams: Dict[str, Path] = {}

    for aligner in aligners:
        aligner_dir = align_dir / aligner
        if not aligner_dir.exists():
            log.warning('  Aligner directory not found: %s', aligner_dir)
            continue

        # Search for .bam first, then .sam
        for ext in ('.bam', '.sam'):
            candidates = sorted(aligner_dir.glob(f'*{chunk_str}*{aligner}{ext}'))
            # Exclude rectified BAMs and empty files
            candidates = [p for p in candidates
                          if 'rectified' not in p.name and p.stat().st_size > 0]
            if candidates:
                bams[aligner] = candidates[0]
                break
        else:
            log.warning('  No BAM/SAM for aligner=%s chunk=%d', aligner, chunk_idx)

    return bams


def discover_chunk_count(run_dir: Path) -> int:
    """Infer number of chunks from BAM filenames in any aligner directory."""
    align_dir = run_dir / 'align'
    for aligner_dir in align_dir.iterdir():
        if not aligner_dir.is_dir():
            continue
        for p in aligner_dir.iterdir():
            if '_of_' in p.name and (p.suffix in ('.bam', '.sam')):
                # Extract N from '_of_NNN'
                try:
                    n_str = p.name.split('_of_')[1].split('.')[0]
                    return int(n_str)
                except (IndexError, ValueError):
                    pass
    raise RuntimeError(f'Cannot determine chunk count from {align_dir}')


# ---------------------------------------------------------------------------
# Core profiling: one chunk
# ---------------------------------------------------------------------------

def open_bam(path: Path, retries: int = 5, base_delay: float = 2.0) -> pysam.AlignmentFile:
    """Open BAM or SAM file with retry logic for transient Oak NFS errors.

    Oak NFS occasionally returns a corrupt/empty view of a file for a fraction
    of a second ('file does not contain alignment data').  Retrying with
    exponential backoff is sufficient to recover in practice.
    """
    import time
    mode = 'rb' if path.suffix == '.bam' else 'r'
    last_exc: Exception = RuntimeError("unreachable")
    for attempt in range(retries):
        try:
            return pysam.AlignmentFile(str(path), mode)
        except (ValueError, OSError) as exc:
            last_exc = exc
            delay = base_delay * (2 ** attempt)
            log.warning("open_bam attempt %d/%d failed for %s (%s); retrying in %.1fs",
                        attempt + 1, retries, path.name, exc, delay)
            time.sleep(delay)
    raise RuntimeError(f"open_bam failed after {retries} attempts for {path}: {last_exc}") from last_exc


_MAPPACBIO_SUFFIX_RE = re.compile(r'_pt:i:\d+$')


def _normalise_read_name(name: str) -> str:
    """Strip aligner-appended suffixes so names match across aligners.

    Two normalisations are applied:
      * Whitespace truncation: bbmap (and some other tools) write the entire
        FASTQ header into QNAME (e.g. ``SRR22434624.1029737 1029737 length=76``)
        while bwa / minimap2 / etc. strip at the first whitespace. SAM spec
        forbids whitespace in QNAME, so the strip is the correct normalisation.
        Without this, bbmap+bwa BAMs have an empty intersection at scale
        (verified 2026-05-17 on Han 2023 QSrev — 783k bbmap, 834k bwa, 0 common
        before fix; 779k common after).
      * mapPacBio's `_pt:i:N` poly-A-length suffix.
    """
    name = name.split(None, 1)[0]
    return _MAPPACBIO_SUFFIX_RE.sub('', name)


def load_reads(bam_path: Path, retries: int = 5, base_delay: float = 2.0) -> Dict[str, pysam.AlignedSegment]:
    """Load all mapped, cigar-having reads into {read_name: read}.

    query_sequence is NOT required — gapmm2 BAMs store no sequence but still
    carry valid = / X CIGAR ops that _build_pos_op_map_from_cigar() can use.
    Read names are normalised to strip aligner-specific suffixes (e.g. the
    mapPacBio '_pt:i:N' poly-A tail annotation).

    Retries the entire load on OSError (Oak NFS can return truncated BAM blocks
    transiently; the error surfaces during iteration, not just file opening).
    """
    import time
    last_exc: Exception = RuntimeError("unreachable")
    for attempt in range(retries):
        try:
            reads: Dict[str, pysam.AlignedSegment] = {}
            with open_bam(bam_path) as bam:
                for read in bam:
                    if (
                        not read.is_unmapped
                        and read.cigartuples
                        and not read.is_secondary
                        and not read.is_supplementary
                    ):
                        reads[_normalise_read_name(read.query_name)] = read
            return reads
        except OSError as exc:
            last_exc = exc
            delay = base_delay * (2 ** attempt)
            log.warning("load_reads attempt %d/%d failed for %s (%s); retrying in %.1fs",
                        attempt + 1, retries, bam_path.name, exc, delay)
            time.sleep(delay)
    raise RuntimeError(f"load_reads failed after {retries} attempts for {bam_path}: {last_exc}") from last_exc


def profile_chunk(
    chunk_bams: Dict[str, Path],
    ref_seqs: Dict[str, str],
    min_aligners: int,
    max_reads: Optional[int],
    counts: Dict[Tuple[str, str, int], int],
    *,
    isolation_flank: int = 0,
    exclude_5p_bp: int = 0,
    exclude_3p_bp: int = 0,
    union_counts: Optional[Dict[Tuple[str, str, int], int]] = None,
    str_counts: Optional[Dict[Tuple[str, str, int], int]] = None,
    strand_counts: Optional[Dict[Tuple[str, str, int, str], int]] = None,
    umi_tag: Optional[str] = None,
    umi_bins: Optional[List[UmiBinSpec]] = None,
    counts_by_bin: Optional[Dict[str, Dict[Tuple[str, str, int], int]]] = None,
    union_counts_by_bin: Optional[Dict[str, Dict[Tuple[str, str, int], int]]] = None,
    str_counts_by_bin: Optional[Dict[str, Dict[Tuple[str, str, int], int]]] = None,
    reads_per_bin: Optional[Dict[str, int]] = None,
    umi_lookup: Optional[Dict[str, int]] = None,
    usable_mask: Optional[Dict[str, "np.ndarray"]] = None,
) -> Tuple[int, int]:
    """Profile one chunk, updating *counts* (and optionally *union_counts*) in place.

    usable_mask (optional, for diploid/outbred genomes):
        dict[chrom] -> bool ndarray (length = chrom length). A reference
        position is counted only where usable_mask[chrom][pos] is True. Built
        from --exclude-vcf (mask known variant sites, indels expanded) and/or
        --high-conf-bed (restrict to high-confidence intervals). When a chrom
        is absent from the mask it is treated as fully excluded if a mask was
        supplied at all. This prevents heterozygous alleles from being counted
        as substitution "errors" -- critical for the HP=1 substitution baseline
        that normalizes the whole penalty table. Unmodified/isogenic data
        (e.g. yeast, IVT) can omit it.

    Two classes of aligner are handled:

    *Sequence-bearing* (minimap2, mapPacBio, deSALT, uLTRA):
        Stores query_sequence; every reference position is classified as
        match (M), substitution (X), or deletion (D).

    *Cigar-only* (gapmm2):
        Does NOT store query_sequence.  Uses explicit = / X CIGAR ops for
        mismatches but writes ambiguous M for most matched positions — so
        matched positions are absent from its pos_op map.

    Agreement rules:
      - D / X / I positions: require at least *min_aligners* in total.
      - M (match) positions: require all *seq-bearing* aligners to agree;
        cigar-only aligners are allowed to abstain (they hold no opinion on
        M vs X at those positions).  This gives us the denominator (total
        covered positions) without discarding the majority of the data.

    isolation_flank (default 0 = off):
        When > 0, only count an error event if every seq-bearing aligner has
        ≥ isolation_flank consecutive exact matches on BOTH sides of the
        contiguous error run (e.g. 10=2D10= requires isolation_flank=10).
        The check is applied at the level of the full error run, so a 2D
        deletion is accepted as long as its outer flanks pass — interior
        positions of the run are not penalised for being adjacent to other
        deleted positions.

    union_counts (default None):
        When provided, populated with "any aligner" error events: a D/X is
        counted if ANY (not all) aligner calls it at that position.  The
        isolation_flank filter is applied to these union events as well,
        using the outer flanks of the contiguous union error run.  The M
        denominator is always taken from the strict-consensus counts.

    str_counts (default None):
        When provided, positions that are in a short tandem repeat (STR)
        context of unit length 2–3 with ≥ 3 copies are counted here instead
        of in *counts*.  The key is (op_type, canonical_unit, n_copies) rather
        than (op_type, ref_base, hp_len).  This separates dinucleotide /
        trinucleotide repeat slippage errors from true HP errors, preventing
        STR events from inflating the HP=1 denominator and baseline.
        Only positions where hp_run_length == 1 are eligible for STR
        classification; positions in genuine HP runs (hp ≥ 2) are always
        counted in *counts*.

    Returns (n_processed, n_common_reads).
    """
    if len(chunk_bams) < min_aligners:
        log.warning('  Only %d aligners available, need %d — skipping chunk',
                    len(chunk_bams), min_aligners)
        return 0, 0

    # Load all reads from each aligner; skip any that fail persistently
    aligner_reads: Dict[str, Dict[str, pysam.AlignedSegment]] = {}
    for aligner, bam_path in chunk_bams.items():
        log.info('    Loading %s ...', bam_path.name)
        try:
            aligner_reads[aligner] = load_reads(bam_path)
            log.info('      %d reads loaded', len(aligner_reads[aligner]))
        except (RuntimeError, OSError) as exc:
            log.warning('    Skipping %s — failed to load after retries: %s', bam_path.name, exc)

    if len(aligner_reads) < min_aligners:
        log.warning('  Only %d aligners loaded, need %d — skipping chunk',
                    len(aligner_reads), min_aligners)
        return 0, 0

    aligners = list(aligner_reads.keys())

    # Detect which aligners store query_sequence
    seq_bearing: set = set()
    for aligner, rdict in aligner_reads.items():
        sample = list(rdict.values())[:200]
        if any(r.query_sequence is not None for r in sample):
            seq_bearing.add(aligner)
    cigar_only = set(aligners) - seq_bearing

    n_seq = len(seq_bearing)
    log.info('  Seq-bearing aligners (%d): %s', n_seq, sorted(seq_bearing))
    if cigar_only:
        log.info('  Cigar-only aligners  (%d): %s  (abstain for M positions)',
                 len(cigar_only), sorted(cigar_only))

    # Index into aligners list
    seq_idxs   = [i for i, a in enumerate(aligners) if a in seq_bearing]
    # all idxs = range(len(aligners)), but useful to have separate for D/X/I

    # Find read names present in ALL aligner dicts
    name_sets = [set(d.keys()) for d in aligner_reads.values()]
    common_names = set.intersection(*name_sets)
    log.info('  %d reads common to all %d aligners', len(common_names), len(aligners))

    n_processed = 0
    n_no_agreed_range = 0
    n_chrom_mismatch = 0
    _lookup_hits = 0
    _lookup_misses = 0

    # Precompute HP run-length arrays for every reference chromosome (fast numpy).
    # STR maps are computed lazily per chromosome (slower Python sweep, ~2–5 s each).
    hp_arrays: Dict[str, np.ndarray] = {
        chrom: _precompute_hp_array(seq) for chrom, seq in ref_seqs.items()
    }
    str_maps: Dict[str, Dict[int, Tuple[str, int]]] = {}

    for read_name in common_names:
        if max_reads is not None and n_processed >= max_reads:
            break

        reads = {a: aligner_reads[a][read_name] for a in aligners}

        # Determine read strand (used when strand_counts is active).
        # All seq-bearing aligners must agree; skip read if they disagree.
        if strand_counts is not None:
            strand_flags = {aligner_reads[a][read_name].is_reverse
                            for a in (seq_bearing if seq_bearing else aligners)}
            if len(strand_flags) != 1:
                continue  # aligners disagree on strand — skip
            read_strand = '-' if next(iter(strand_flags)) else '+'
        else:
            read_strand = None

        # All aligners must agree on chromosome
        chroms = {r.reference_name for r in reads.values()}
        if len(chroms) != 1:
            n_chrom_mismatch += 1
            continue
        chrom = next(iter(chroms))

        ref_seq = ref_seqs.get(chrom)
        if ref_seq is None:
            continue

        # Per-chrom usable-position mask (variant / high-conf filtering).
        # None => no masking for this chrom. If a mask was supplied globally
        # but this chrom is absent, skip the whole read group (conservative:
        # we have no confidence calls there).
        if usable_mask is not None:
            usable_chrom = usable_mask.get(chrom)
            if usable_chrom is None:
                continue
        else:
            usable_chrom = None

        # --- UMI bin lookup (once per read, before per-range loops) ---
        # The bin label drives parallel per-bin counters; reads with no UMI
        # tag (or no matching bin) get bin_label=None — those still feed the
        # pooled `counts` (and union/str), but NOT any per-bin output.
        bin_label: Optional[str] = None
        if umi_bins and counts_by_bin is not None and umi_tag:
            umi_val: Optional[int] = None
            if umi_lookup is not None:
                # Sidecar BAM path: use prebuilt name→tag dict (read_name is
                # already normalised, matching umi_lookup keys).
                raw = umi_lookup.get(read_name)
                if raw is not None:
                    try:
                        umi_val = int(raw)
                        _lookup_hits += 1
                    except (TypeError, ValueError):
                        pass
                else:
                    _lookup_misses += 1
            else:
                # Inline tag path: pick any aligner that carries the tag —
                # they should all match but only one BAM tends to propagate
                # aux tags from the input (e.g. minimap2 with -y).
                for _a in aligners:
                    _r = reads[_a]
                    try:
                        umi_val = int(_r.get_tag(umi_tag))
                        break
                    except KeyError:
                        continue
                    except (TypeError, ValueError):
                        continue
            bin_label = assign_bin(umi_val, umi_bins)
            if bin_label is not None and reads_per_bin is not None:
                reads_per_bin[bin_label] = reads_per_bin.get(bin_label, 0) + 1

        # Retrieve precomputed HP array; build STR map lazily (only when needed)
        hp_arr = hp_arrays[chrom]
        if str_counts is not None and chrom not in str_maps:
            str_maps[chrom] = _precompute_str_map(ref_seq)
        str_map = str_maps.get(chrom, {})

        # Exonic ranges per aligner, then intersect.
        # When exclude_*_bp > 0, ranges are pre-trimmed to the read-coord interior
        # (cuts protocol-specific end-of-read noise: template-switch bleed for cDNA,
        # CSP A-tract walkback for QSrev). DRS defaults to no exclusion.
        if exclude_5p_bp > 0 or exclude_3p_bp > 0:
            exonic_per_aligner = [
                get_exonic_ranges_inner(reads[a], exclude_5p_bp, exclude_3p_bp)
                for a in aligners
            ]
        else:
            exonic_per_aligner = [get_exonic_ranges(reads[a]) for a in aligners]
        agreed_ranges = intersect_all(exonic_per_aligner)
        if not agreed_ranges:
            n_no_agreed_range += 1
            continue

        # Build pos-op maps for each aligner: [(pos_ops_dict, ins_before_dict), ...]
        all_pos_ops = [build_pos_op_map(reads[a], ref_seq) for a in aligners]

        # Precompute per-aligner M-run flank lengths (used by isolation filter)
        if isolation_flank > 0:
            flank_data = [
                _precompute_flank_lengths(all_pos_ops[i][0], agreed_ranges)
                for i in seq_idxs
            ]
        else:
            flank_data = None

        # -----------------------------------------------------------------------
        # Process each agreed exonic range
        # -----------------------------------------------------------------------
        for rng_start, rng_end in agreed_ranges:

            # --- Phase 1: classify every position in this range (vectorised) ---
            # Build a dense uint8 ops_matrix of shape (n_aligners, rng_len):
            #   0=M, 1=X, 2=D, 255=absent (position not in aligner's pos_ops).
            # Vectorised numpy consensus then replaces the original per-position
            # Python loop (list allocs + set builds per position eliminated).
            rng_len = rng_end - rng_start
            ops_matrix = np.full((len(aligners), rng_len), 255, dtype=np.uint8)
            for _i, (_pod, _) in enumerate(all_pos_ops):
                for _p in range(rng_start, rng_end):
                    _entry = _pod.get(_p)
                    if _entry is not None:
                        _code = _OP_CODE.get(_entry[0])
                        if _code is not None:
                            ops_matrix[_i, _p - rng_start] = _code

            # Presence / consensus masks
            present      = ops_matrix != 255                            # (n_aln, rng_len)
            all_present  = np.all(present, axis=0)                     # (rng_len,)
            first_op     = ops_matrix[0]                               # (rng_len,)
            all_agree    = all_present & np.all(
                ops_matrix == first_op[np.newaxis, :], axis=0
            )

            # Seq-bearing subset present at each position
            if seq_idxs:
                seq_present_mask = np.all(present[seq_idxs], axis=0)
            else:
                seq_present_mask = all_present

            # Strict M: all agree, op=M, all seq-bearing present
            m_mask = all_agree & (first_op == 0) & seq_present_mask

            # Strict D/X: all agree, at least min_aligners present
            min_aln_mask = np.sum(present, axis=0) >= min_aligners
            d_mask = all_agree & (first_op == 2) & min_aln_mask
            x_mask = all_agree & (first_op == 1) & min_aln_mask

            # Union any-aligner masks (computed only when needed)
            if union_counts is not None:
                any_d_mask = np.any(ops_matrix == 2, axis=0)
                any_x_mask = (
                    np.any(ops_matrix[seq_idxs] == 1, axis=0)
                    if seq_idxs else np.zeros(rng_len, dtype=bool)
                )
                any_del_pos: Set[int] = set(
                    (np.where(any_d_mask)[0] + rng_start).tolist()
                )
                any_sub_pos: Set[int] = set(
                    (np.where(any_x_mask)[0] + rng_start).tolist()
                )
                if usable_chrom is not None:
                    any_del_pos = {p for p in any_del_pos if usable_chrom[p]}
                    any_sub_pos = {p for p in any_sub_pos if usable_chrom[p]}
            else:
                any_del_pos = set()
                any_sub_pos = set()

            # Count strict M positions (iterate only matched positions — small subset)
            for _rel in np.where(m_mask)[0]:
                pos = int(_rel) + rng_start
                if usable_chrom is not None and not usable_chrom[pos]:
                    continue
                ref_base = ref_seq[pos]
                hp_len = int(hp_arr[pos])
                if hp_len == 1 and str_counts is not None:
                    su, sc = str_map.get(pos, (None, 0))
                    if su is not None:
                        str_counts[('M', su, sc)] += 1
                        if bin_label is not None and str_counts_by_bin is not None:
                            str_counts_by_bin[bin_label][('M', su, sc)] += 1
                        if union_counts is not None:
                            union_counts[('M', ref_base, hp_len)] += 1
                            if bin_label is not None and union_counts_by_bin is not None:
                                union_counts_by_bin[bin_label][('M', ref_base, hp_len)] += 1
                        continue   # exclude from HP counts (cleaner baseline)
                counts[('M', ref_base, hp_len)] += 1
                if bin_label is not None and counts_by_bin is not None:
                    counts_by_bin[bin_label][('M', ref_base, hp_len)] += 1
                if strand_counts is not None and read_strand is not None:
                    strand_counts[('M', ref_base, hp_len, read_strand)] += 1
                if union_counts is not None:
                    union_counts[('M', ref_base, hp_len)] += 1
                    if bin_label is not None and union_counts_by_bin is not None:
                        union_counts_by_bin[bin_label][('M', ref_base, hp_len)] += 1

            # Collect strict D/X error events (iterate only error positions)
            strict_error_events: List[Tuple[int, str, str]] = []
            for _rel in np.where(d_mask | x_mask)[0]:
                pos = int(_rel) + rng_start
                if usable_chrom is not None and not usable_chrom[pos]:
                    continue
                ref_base = ref_seq[pos]
                op_type = 'D' if d_mask[_rel] else 'X'
                strict_error_events.append((pos, op_type, ref_base))

            # --- Phase 2: group strict error events into contiguous runs ---
            # This lets us apply the isolation filter to the OUTER flanks of
            # a multi-position deletion (e.g. 10=2D10=) rather than per-position,
            # which would incorrectly reject interior positions of the run.
            strict_runs: List[List[Tuple[int, str, str]]] = []
            if strict_error_events:
                cur: List[Tuple[int, str, str]] = [strict_error_events[0]]
                for item in strict_error_events[1:]:
                    pos, op_t, rb = item
                    ppos, pop, _ = cur[-1]
                    if op_t == pop and pos == ppos + 1:
                        cur.append(item)
                    else:
                        strict_runs.append(cur)
                        cur = [item]
                strict_runs.append(cur)

            # --- Phase 3: apply isolation filter and count strict events ---
            for run in strict_runs:
                run_start_p = run[0][0]
                run_end_p   = run[-1][0]

                # Isolation check: all seq-bearing aligners must have >= N M's
                # immediately outside the contiguous error run.
                isolated = True
                if isolation_flank > 0 and flank_data:
                    for lf, rf in flank_data:
                        if (lf.get(run_start_p, 0) < isolation_flank or
                                rf.get(run_end_p, 0) < isolation_flank):
                            isolated = False
                            break

                if isolated:
                    # Check STR context using the first position in the run.
                    # A run is classified as STR only when every position in the
                    # run has hp_len=1 (genuine HP runs always go to counts).
                    run_su: Optional[str] = None
                    run_sc: int = 0
                    if str_counts is not None:
                        first_pos = run[0][0]
                        if all(hp_arr[p] == 1 for p, _, _ in run):
                            run_su, run_sc = str_map.get(first_pos, (None, 0))
                    for pos, op_type, ref_base in run:
                        if run_su is not None:
                            str_counts[(op_type, run_su, run_sc)] += 1
                            if bin_label is not None and str_counts_by_bin is not None:
                                str_counts_by_bin[bin_label][(op_type, run_su, run_sc)] += 1
                        else:
                            hp_len = int(hp_arr[pos])
                            counts[(op_type, ref_base, hp_len)] += 1
                            if bin_label is not None and counts_by_bin is not None:
                                counts_by_bin[bin_label][(op_type, ref_base, hp_len)] += 1
                            if strand_counts is not None and read_strand is not None:
                                strand_counts[(op_type, ref_base, hp_len, read_strand)] += 1

            # --- Phase 4: union error events (any aligner) ---
            if union_counts is not None:
                for pos_set, op_char in [(any_del_pos, 'D'), (any_sub_pos, 'X')]:
                    for run_start_p, run_end_p in _group_pos_into_runs(pos_set):
                        isolated = True
                        if isolation_flank > 0 and flank_data:
                            for lf, rf in flank_data:
                                if (lf.get(run_start_p, 0) < isolation_flank or
                                        rf.get(run_end_p, 0) < isolation_flank):
                                    isolated = False
                                    break
                        if isolated:
                            for pos in range(run_start_p, run_end_p + 1):
                                if 0 <= pos < len(ref_seq):
                                    ref_base = ref_seq[pos]
                                    hp_len = int(hp_arr[pos])
                                    union_counts[(op_char, ref_base, hp_len)] += 1
                                    if bin_label is not None and union_counts_by_bin is not None:
                                        union_counts_by_bin[bin_label][(op_char, ref_base, hp_len)] += 1

            # --- Phase 5: insertions (strict and union) ---
            strict_ins_before: Set[int] = set()
            any_ins_before: Set[int] = set()

            for pos in range(rng_start + 1, rng_end + 1):
                ins_vals = [ins.get(pos, 0) for _, ins in all_pos_ops]
                n_with   = sum(1 for v in ins_vals if v > 0)
                if n_with >= min_aligners:
                    strict_ins_before.add(pos)
                if union_counts is not None and n_with >= 1:
                    any_ins_before.add(pos)

            for pos_set, target_counts, target_by_bin in [
                (strict_ins_before, counts, counts_by_bin),
                (any_ins_before, union_counts, union_counts_by_bin),
            ]:
                if target_counts is None:
                    continue
                for pos in pos_set:
                    # Isolation check for insertion before `pos`:
                    #   left  flank = M's ending at pos-1  → lf[pos]
                    #   right flank = M's starting at pos  → rf[pos-1]
                    isolated = True
                    if isolation_flank > 0 and flank_data:
                        for lf, rf in flank_data:
                            if (lf.get(pos, 0) < isolation_flank or
                                    rf.get(pos - 1, 0) < isolation_flank):
                                isolated = False
                                break
                    if isolated:
                        hp_left  = int(hp_arr[pos - 1]) if pos > 0 else 1
                        hp_right = int(hp_arr[pos]) if pos < len(ref_seq) else 1
                        hp_len   = max(hp_left, hp_right)
                        hp_base  = (ref_seq[pos - 1] if hp_left >= hp_right
                                    else ref_seq[pos])
                        target_counts[('I', hp_base, hp_len)] += 1
                        if bin_label is not None and target_by_bin is not None:
                            target_by_bin[bin_label][('I', hp_base, hp_len)] += 1

        n_processed += 1

    log.info('  Processed %d reads  (chrom_mismatch=%d  no_agreed_range=%d)',
             n_processed, n_chrom_mismatch, n_no_agreed_range)
    if umi_lookup is not None and umi_bins:
        total_looked = _lookup_hits + _lookup_misses
        pct = 100.0 * _lookup_hits / total_looked if total_looked else 0.0
        log.info('  Sidecar UMI lookup: %d hits / %d total (%.1f%% hit rate)',
                 _lookup_hits, total_looked, pct)
    return n_processed, len(common_names)


# ---------------------------------------------------------------------------
# Rate computation & penalty derivation
# ---------------------------------------------------------------------------

def compute_rates(counts: Dict[Tuple[str, str, int], int]) -> pd.DataFrame:
    """Compute error rates per (op_type, hp_length) from raw counts.

    Denominator for deletions and substitutions:
        total_positions(hp) = sum of (M + D + X) counts at hp_length=hp
        This is the number of agreed reference positions at that HP context.

    Denominator for insertions:
        Same total_positions(hp) — insertions are counted as "extra events"
        relative to the positions surveyed.

    Returns a DataFrame with columns:
        op_type, hp_length, base, count, total_positions, rate
    """
    rows = []
    for (op_type, ref_base, hp_len), cnt in counts.items():
        rows.append({
            'op_type': op_type,
            'ref_base': ref_base,
            'hp_length': hp_len,
            'count': cnt,
        })

    if not rows:
        return pd.DataFrame(columns=['op_type', 'ref_base', 'hp_length',
                                     'count', 'total_positions', 'rate'])

    df = pd.DataFrame(rows)

    # Compute denominator = M + D + X counts at each (ref_base, hp_length)
    # (positions that were actually covered and agreed upon)
    denom_ops = {'M', 'D', 'X'}
    denom_df = (
        df[df['op_type'].isin(denom_ops)]
        .groupby(['ref_base', 'hp_length'], as_index=False)['count']
        .sum()
        .rename(columns={'count': 'total_positions'})
    )

    df = df.merge(denom_df, on=['ref_base', 'hp_length'], how='left')
    df['total_positions'] = df['total_positions'].fillna(0).astype(int)
    df['rate'] = df['count'] / df['total_positions'].clip(lower=1)

    return df.sort_values(['op_type', 'hp_length', 'ref_base'])


def _isotonic_rates(rates: np.ndarray) -> np.ndarray:
    """Enforce monotonically non-decreasing rates via pool-adjacent-violators.

    Nanopore deletion rates should increase (or stay flat) with HP run length.
    When empirical estimates violate this at high HP lengths (sparse data), we
    pool adjacent blocks to enforce the expected relationship.

    Args:
        rates: 1-D array of rates indexed 0 = HP 1, 1 = HP 2, etc.

    Returns:
        Smoothed rates array with the same shape, values non-decreasing.
    """
    n = len(rates)
    # Represent each element as a block: (mean_value, weight=1)
    blocks = [[float(rates[i]), 1] for i in range(n)]
    i = 0
    while i < len(blocks) - 1:
        if blocks[i][0] > blocks[i + 1][0]:
            # Merge block i and i+1 into a single pooled block
            combined_sum = blocks[i][0] * blocks[i][1] + blocks[i + 1][0] * blocks[i + 1][1]
            combined_weight = blocks[i][1] + blocks[i + 1][1]
            blocks[i] = [combined_sum / combined_weight, combined_weight]
            blocks.pop(i + 1)
            # Step back to check the newly merged block against its predecessor
            if i > 0:
                i -= 1
        else:
            i += 1
    # Expand blocks back to original length
    result = []
    for val, weight in blocks:
        result.extend([val] * int(weight))
    return np.array(result)


def derive_penalties(
    df: pd.DataFrame,
    default_ins: float = 1.25,
    min_count_flag: int = 50,
) -> pd.DataFrame:
    """Derive penalty scores from empirical rates, split by base class (AT vs CG).

    **Base-class splitting:**

    A/T and C/G homopolymers have dramatically different error rates (up to 40×
    at HP=8), so penalties are derived separately for each base class.  The output
    has a ``base_class`` column ('AT' or 'CG') so the junction refiner can look up
    the correct cost per reference position.

    **Normalization (per base class):**

    - Substitutions (X): normalized so that ``sub(HP=1, base_class) = 1.0``.
    - Deletions (D): normalized to the same per-class sub(HP=1) reference.
    - Insertions (I): normalized independently so that ``ins(HP=1) = default_ins``
      using the per-class insertion baseline.
    - Matches (M): included for diagnostics only.

    **Monotonicity:**  Pool-adjacent-violators isotonic regression is applied
    independently to deletion rates within each base class.

    **Low-count flagging:**  Rows with ``count_total < min_count_flag`` are
    marked ``low_count = True``.

    Returns a DataFrame with columns:
        op_type, base_class, hp_length, rate_mean, count_total, penalty_score, low_count
    """
    df = df.copy()
    df['base_class'] = df['ref_base'].map(lambda b: 'AT' if b in ('A', 'T') else 'CG')

    # Aggregate rates across individual bases within each class
    agg = (
        df.groupby(['op_type', 'base_class', 'hp_length'])
        .apply(lambda g: pd.Series({
            'rate_mean':   np.average(g['rate'], weights=g['count'].clip(lower=1)),
            'count_total': g['count'].sum(),
        }))
        .reset_index()
    )

    # --- Per-class normalization reference rates ---
    ref_sub: Dict[str, float] = {}
    ref_ins: Dict[str, Optional[float]] = {}
    for bc in ('AT', 'CG'):
        sub_hp1 = agg.loc[
            (agg['op_type'] == 'X') & (agg['base_class'] == bc) & (agg['hp_length'] == 1),
            'rate_mean',
        ]
        if sub_hp1.empty or float(sub_hp1.iloc[0]) == 0:
            log.warning('No sub data at HP=1 for base_class=%s — using 0.02 fallback', bc)
            ref_sub[bc] = 0.02
        else:
            ref_sub[bc] = float(sub_hp1.iloc[0])

        ins_hp1 = agg.loc[
            (agg['op_type'] == 'I') & (agg['base_class'] == bc) & (agg['hp_length'] == 1),
            'rate_mean',
        ]
        ref_ins[bc] = float(ins_hp1.iloc[0]) if (not ins_hp1.empty and float(ins_hp1.iloc[0]) > 0) else None

    log.info('Reference rates: sub(HP=1) AT=%.4f CG=%.4f  ins(HP=1) AT=%s CG=%s',
             ref_sub.get('AT', 0), ref_sub.get('CG', 0),
             f'{ref_ins["AT"]:.4f}' if ref_ins.get('AT') else 'n/a',
             f'{ref_ins["CG"]:.4f}' if ref_ins.get('CG') else 'n/a')

    # --- Apply isotonic regression per base class for deletions ---
    result_parts = []
    non_del = agg[agg['op_type'] != 'D']
    result_parts.append(non_del)
    for bc in ('AT', 'CG'):
        del_mask = (agg['op_type'] == 'D') & (agg['base_class'] == bc)
        if not del_mask.any():
            continue
        del_agg = agg[del_mask].sort_values('hp_length').copy()
        smoothed = _isotonic_rates(del_agg['rate_mean'].values)
        n_adj = int((smoothed != del_agg['rate_mean'].values).sum())
        if n_adj > 0:
            log.info('Isotonic regression adjusted %d deletion rate(s) for base_class=%s', n_adj, bc)
        del_agg['rate_mean'] = smoothed
        result_parts.append(del_agg)
    agg = pd.concat(result_parts, ignore_index=True)

    # --- Compute penalty scores ---
    def _penalty(row: pd.Series) -> float:
        bc  = row['base_class']
        op  = row['op_type']
        rate = row['rate_mean']
        if rate == 0:
            return 10.0
        if op == 'I':
            ri = ref_ins.get(bc)
            if ri is None:
                return default_ins
            return min(default_ins * ri / rate, 10.0)
        else:
            rs = ref_sub.get(bc, 0.02)
            return min(rs / rate, 10.0)

    agg['penalty_score'] = agg.apply(_penalty, axis=1)
    agg['low_count'] = agg['count_total'] < min_count_flag

    return agg.sort_values(['op_type', 'base_class', 'hp_length']).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

_OP_DISPLAY = {'D': 'deletion', 'X': 'substitution', 'I': 'insertion'}


def plot_error_rates(penalties: pd.DataFrame, output_dir: Path) -> None:
    """Line plot: error rate vs HP run length, one line per op type."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # --- Left: rate vs HP length ---
    ax = axes[0]
    for op_char, label in _OP_DISPLAY.items():
        sub = penalties[penalties['op_type'] == op_char]
        if sub.empty:
            continue
        sub = sub.sort_values('hp_length')
        # Only plot HP lengths with reasonable counts
        sub = sub[sub['count_total'] >= 10]
        ax.plot(sub['hp_length'], sub['rate_mean'],
                marker='o', label=label, color=_OP_COLORS[label], linewidth=2)

    ax.set_xlabel('Homopolymer run length')
    ax.set_ylabel('Empirical error rate')
    ax.set_title('Nanopore error rate vs HP context\n(multi-aligner agreed positions)')
    ax.legend()
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)

    # --- Right: penalty score vs HP length ---
    ax = axes[1]
    for op_char, label in _OP_DISPLAY.items():
        sub = penalties[penalties['op_type'] == op_char]
        if sub.empty:
            continue
        sub = sub.sort_values('hp_length')
        sub = sub[sub['count_total'] >= 10]
        ax.plot(sub['hp_length'], sub['penalty_score'],
                marker='o', label=label, color=_OP_COLORS[label], linewidth=2)

    # Reference lines for current heuristic
    ax.axhline(1.00, color='gray', linestyle='--', linewidth=1,
               label='current del/sub (HP=1) = 1.0')
    ax.axhline(0.70, color='orange', linestyle='--', linewidth=1,
               label='current ins (HP=1) = 0.70')
    ax.axhline(0.50, color='red', linestyle=':', linewidth=1,
               label='current del (HP≥4) = 0.50')

    ax.set_xlabel('Homopolymer run length')
    ax.set_ylabel('Penalty score  (sub, HP=1 → 1.0)')
    ax.set_title('Empirical penalty scores vs HP context\n(dashed = current heuristics)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    out_path = output_dir / 'error_rate_vs_hp_length.png'
    plt.savefig(str(out_path), dpi=150, bbox_inches='tight')
    plt.close()
    log.info('Plot saved: %s', out_path)


def plot_hp_length_distribution(
    df: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Bar plot: count of agreed error events by HP run length."""
    fig, ax = plt.subplots(figsize=(10, 4))

    error_ops = {'D': 'deletion', 'X': 'substitution', 'I': 'insertion'}
    bottoms = defaultdict(float)
    max_hp = int(df['hp_length'].max()) if not df.empty else 20
    x = np.arange(1, max_hp + 1)

    for op_char, label in error_ops.items():
        sub = (
            df[df['op_type'] == op_char]
            .groupby('hp_length')['count']
            .sum()
            .reindex(x, fill_value=0)
        )
        ax.bar(x, sub.values, bottom=[bottoms[xi] for xi in x],
               label=label, color=_OP_COLORS[label], alpha=0.8, width=0.6)
        for xi, v in zip(x, sub.values):
            bottoms[xi] += v

    ax.set_xlabel('Homopolymer run length')
    ax.set_ylabel('Agreed error events (count)')
    ax.set_title('Distribution of agreed error events by HP context')
    ax.legend()
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, axis='y')

    out_path = output_dir / 'hp_length_distribution.png'
    plt.savefig(str(out_path), dpi=150, bbox_inches='tight')
    plt.close()
    log.info('Plot saved: %s', out_path)


def print_penalty_table(penalties: pd.DataFrame, max_hp: int = 12) -> None:
    """Print a summary penalty table to stdout, showing AT and CG base classes."""
    print('\n' + '=' * 88)
    print('EMPIRICAL PENALTY TABLE  (sub at HP=1 → 1.00 per base class)')
    print('=' * 88)
    hdr = (f"{'HP':>4}  "
           f"{'del_AT':>8}  {'del_CG':>8}  "
           f"{'sub_AT':>8}  {'sub_CG':>8}  "
           f"{'ins_AT':>8}  {'ins_CG':>8}")
    print(hdr)
    print('-' * 88)

    has_bc = 'base_class' in penalties.columns

    def _get(op, bc, hp):
        if has_bc:
            match = penalties[(penalties['op_type'] == op) &
                              (penalties['base_class'] == bc) &
                              (penalties['hp_length'] == hp)]
        else:
            match = penalties[(penalties['op_type'] == op) & (penalties['hp_length'] == hp)]
        return f'{match["penalty_score"].iloc[0]:.3f}' if not match.empty else '   —  '

    for hp in range(1, max_hp + 1):
        print(f"{hp:>4}  "
              f"{_get('D','AT',hp):>8}  {_get('D','CG',hp):>8}  "
              f"{_get('X','AT',hp):>8}  {_get('X','CG',hp):>8}  "
              f"{_get('I','AT',hp):>8}  {_get('I','CG',hp):>8}")

    print('=' * 88)
    print()
    print('Current heuristic (junction_refiner._hp_edit_distance defaults):')
    print('  HP < 4 : del=1.00, ins=0.70, sub=1.00')
    print('  HP = 4 : del=0.50, ins=0.35, sub=1.00')
    print('  HP = 8 : del=0.25, ins=0.18, sub=1.00')
    print()


def _write_str_outputs(
    str_counts: Dict[Tuple[str, str, int], int],
    output_dir: Path,
    max_hp: int = 20,
) -> None:
    """Compute and write STR error rate TSV + print a console summary.

    str_counts keys are (op_type, canonical_unit, n_copies) where op_type is
    one of 'M' (denominator), 'D', 'X', 'I'.  Rates are computed per
    (unit, n_copies) group as error_count / total_positions (M+D+X).
    """
    if not str_counts:
        return

    rows = []
    for (op_type, unit, n_copies), cnt in str_counts.items():
        rows.append({'op_type': op_type, 'unit': unit, 'n_copies': n_copies, 'count': cnt})

    df = pd.DataFrame(rows)
    if df.empty:
        return

    # Denominator: M + D + X per (unit, n_copies)
    denom_df = (
        df[df['op_type'].isin({'M', 'D', 'X'})]
        .groupby(['unit', 'n_copies'], as_index=False)['count']
        .sum()
        .rename(columns={'count': 'total_positions'})
    )
    df = df.merge(denom_df, on=['unit', 'n_copies'], how='left')
    df['total_positions'] = df['total_positions'].fillna(0).astype(int)
    df['rate'] = df['count'] / df['total_positions'].clip(lower=1)

    # --- STR penalty scores (normalize to overall HP=1 baseline from main counts) ---
    # Use a fallback reference rate of 0.02 if we can't derive one from STR data.
    sub_str = df[(df['op_type'] == 'X')]
    ref_rate_sub_str = (
        float(sub_str.loc[sub_str['total_positions'] > 50, 'rate'].max())
        if not sub_str.empty else 0.02
    )
    ref_rate_sub_str = max(ref_rate_sub_str, 0.001)  # floor to avoid division issues

    ins_str = df[(df['op_type'] == 'I') & (df['n_copies'] == df['n_copies'].min())]
    ref_rate_ins_str = float(ins_str['rate'].mean()) if not ins_str.empty else None

    def _str_penalty(row: pd.Series) -> float:
        rate = row['rate']
        if rate == 0:
            return 10.0
        if row['op_type'] == 'I':
            if ref_rate_ins_str is None:
                return 1.25
            return min(1.25 * ref_rate_ins_str / rate, 10.0)
        else:
            return min(ref_rate_sub_str / rate, 10.0)

    df['penalty_score'] = df.apply(_str_penalty, axis=1)
    df['low_count'] = df['count'] < 50

    pen_path = output_dir / 'str_penalty_scores.tsv'
    (df[df['op_type'].isin({'D', 'X', 'I'})]
     .sort_values(['op_type', 'unit', 'n_copies'])
     .to_csv(pen_path, sep='\t', index=False, float_format='%.4f'))
    log.info('Wrote: %s', pen_path)

    out_path = output_dir / 'str_error_rates.tsv'
    df.sort_values(['op_type', 'unit', 'n_copies']).to_csv(
        out_path, sep='\t', index=False, float_format='%.6f')
    log.info('Wrote: %s', out_path)

    # Console summary: top STR contexts by deletion rate
    del_df = (
        df[df['op_type'] == 'D']
        .sort_values('rate', ascending=False)
        .head(20)
    )
    if not del_df.empty:
        print('\n' + '=' * 72)
        print('STR DELETION RATES  (top contexts by rate)')
        print('Positions in HP=1 dinucleotide/trinucleotide repeats ≥3 copies')
        print('=' * 72)
        print(f"{'unit':>6}  {'copies':>6}  {'count_D':>8}  {'total_pos':>10}  {'del_rate':>10}")
        print('-' * 72)
        for _, row in del_df.iterrows():
            m_row = df[(df['op_type'] == 'M') & (df['unit'] == row['unit'])
                       & (df['n_copies'] == row['n_copies'])]
            m_cnt = int(m_row['count'].sum()) if not m_row.empty else 0
            print(f"{row['unit']:>6}  {int(row['n_copies']):>6}  "
                  f"{int(row['count']):>8}  {int(row['total_positions']):>10}  "
                  f"{row['rate']:>10.4f}")
        print('=' * 72)
        print()


def _print_comparison_table(
    strict_df: pd.DataFrame,
    union_df: pd.DataFrame,
    max_hp: int = 12,
) -> None:
    """Print strict vs union penalty side-by-side for quick visual validation.

    If strict ≈ union at each HP length, positional ambiguity and aligner
    artifacts are not significantly biasing the strict consensus profiler.
    """
    print('\n' + '=' * 88)
    print('STRICT vs UNION COMPARISON  (isolation filter applied to both)')
    print('If values are similar, positional ambiguity is not a major concern.')
    print('=' * 88)
    hdr = (f"{'HP':>4}  "
           f"{'del strict':>10}  {'del union':>10}  {'del ratio':>10}  "
           f"{'sub strict':>10}  {'sub union':>10}  {'sub ratio':>10}")
    print(hdr)
    print('-' * 88)

    for hp in range(1, max_hp + 1):
        row = {}
        for label, df in [('strict', strict_df), ('union', union_df)]:
            for op_char, name in [('D', 'del'), ('X', 'sub')]:
                match = df[(df['op_type'] == op_char) & (df['hp_length'] == hp)]
                key = f'{name}_{label}'
                row[key] = match['penalty_score'].iloc[0] if not match.empty else None

        def fmt(v):
            return f'{v:.3f}' if v is not None else '   —  '

        def ratio(s, u):
            if s is None or u is None or u == 0:
                return '   —  '
            return f'{s / u:.3f}'

        print(f"{hp:>4}  "
              f"{fmt(row.get('del_strict')):>10}  "
              f"{fmt(row.get('del_union')):>10}  "
              f"{ratio(row.get('del_strict'), row.get('del_union')):>10}  "
              f"{fmt(row.get('sub_strict')):>10}  "
              f"{fmt(row.get('sub_union')):>10}  "
              f"{ratio(row.get('sub_strict'), row.get('sub_union')):>10}")

    print('=' * 88)
    print()
    print('Ratio < 1: union detects more events than strict (some positional ambiguity).')
    print('Ratio ≈ 1: strict consensus captures the same events as union — results reliable.')
    print('Ratio > 1: strict overcounts relative to union (unexpected; check data).')
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument(
        '--run-dir', type=Path,
        help='Rectify dev-run directory containing align/<aligner>/ subdirs. '
             'Chunk BAMs are auto-discovered.',
    )
    src.add_argument(
        '--aligner-bams', nargs='+', metavar='ALIGNER:PATH',
        help='Explicit list of aligner:path pairs (bypasses --run-dir). '
             'E.g. "minimap2:/path/to/mm2.bam mapPacBio:/path/to/mp.bam"',
    )

    p.add_argument('--reference', type=Path, required=True,
                   help='Reference genome FASTA (must be indexed with samtools faidx).')
    p.add_argument('--output-dir', type=Path, required=True,
                   help='Output directory for TSVs and plots.')
    p.add_argument('--exclude-vcf', type=Path, default=None, metavar='VCF',
                   help='Mask reference positions overlapping records in this VCF '
                        '(bgzipped+tabix; indels expanded by REF length). For diploid/'
                        'outbred samples this excludes heterozygous/variant alleles that '
                        'would otherwise be counted as substitution "errors", biasing the '
                        'HP=1 substitution baseline. Not needed for isogenic/unmodified '
                        'data (yeast, IVT).')
    p.add_argument('--high-conf-bed', type=Path, default=None, metavar='BED',
                   help='Restrict counting to reference positions inside these intervals '
                        '(e.g. a GIAB high-confidence-regions BED). Positions outside, and '
                        'whole chromosomes absent from the BED, are excluded.')
    p.add_argument('--protocol', choices=['drs', 'cdna', 'qsrev'], default='drs',
                   help='Sequencing protocol. Controls (a) default --aligners panel '
                        '(qsrev=bbmap+bwa; drs/cdna=5-aligner long-read panel) and '
                        '(b) default read-coordinate exclusion: drs=(5p=0,3p=0) — '
                        'matches bundled DRS tables; cdna=(50,50); qsrev=(30,50). '
                        'Explicit --aligners / --exclude-Nprime-bp override.')
    p.add_argument('--aligners', nargs='+', default=None,
                   metavar='ALIGNER',
                   help='Aligners to use. Default depends on --protocol: drs/cdna='
                        f'{" ".join(_DEFAULT_ALIGNERS_BY_PROTOCOL["drs"])}; qsrev='
                        f'{" ".join(_DEFAULT_ALIGNERS_BY_PROTOCOL["qsrev"])}.')
    p.add_argument('--exclude-5prime-bp', type=int, default=None, metavar='N',
                   dest='exclude_5p_bp',
                   help='Exclude positions whose read-coordinate (BAM orientation, '
                        'after soft-clip stripping) falls within N bp of the read 5-prime '
                        'end. Default driven by --protocol. Set 0 to disable.')
    p.add_argument('--exclude-3prime-bp', type=int, default=None, metavar='N',
                   dest='exclude_3p_bp',
                   help='Exclude positions whose read-coordinate (BAM orientation, '
                        'after soft-clip stripping) falls within N bp of the read 3-prime '
                        'end. Default driven by --protocol. Set 0 to disable.')
    p.add_argument('--min-aligners', type=int, default=None,
                   help='Minimum aligners that must agree (default: all provided).')
    p.add_argument('--chunks', nargs='+', type=int, default=None,
                   metavar='N',
                   help='Which chunk indices to process (default: all). '
                        'E.g. --chunks 0 1 2')
    p.add_argument('--max-reads', type=int, default=None,
                   help='Max reads per chunk (useful for quick tests).')
    p.add_argument('--max-hp', type=int, default=20,
                   help='Maximum HP run length to include in outputs (default: 20).')
    p.add_argument('--default-ins', type=float, default=1.25,
                   help='Insertion penalty at HP=1 (default: 1.25, matching the junction_refiner '
                        'heuristic baseline). Insertion penalties at longer HP lengths are scaled '
                        'relative to this value. Prevents all insertion rows clipping to 10.0 '
                        'when insertions are far rarer than substitutions.')
    p.add_argument('--min-count-flag', type=int, default=50,
                   help='Mark rows with fewer than this many agreed events as low_count=True '
                        'in penalty_scores.tsv (default: 50). HpPenaltyTable.from_tsv() uses '
                        'this threshold to skip unreliable high-HP entries.')
    p.add_argument('--isolation-flank', type=int, default=0, metavar='N',
                   help='Only count error events flanked by ≥N consecutive exact matches '
                        'on both sides in all seq-bearing aligners (default: 0 = off). '
                        'E.g. --isolation-flank 10 requires patterns like 10=2D10= or '
                        '10=1X10=.  Isolation is checked at the outer edges of contiguous '
                        'multi-position error runs, so 2D is accepted with 10= on each '
                        'outer side even though interior positions adjoin each other.')
    p.add_argument('--union', action='store_true', default=False,
                   help='In addition to the strict all-aligner consensus counts, also '
                        'accumulate "union" counts where ANY aligner sees D/X/I at a '
                        'position.  The same --isolation-flank filter is applied to both. '
                        'Outputs a second set of TSVs and penalty tables prefixed with '
                        '"union_" for comparison.  If strict ≈ union rates, positional '
                        'ambiguity within HP runs is not inflating/deflating the profiler '
                        'results.')
    p.add_argument('--workers', type=int, default=1, metavar='N',
                   help='Number of parallel worker processes for chunk profiling '
                        '(default: 1 = sequential).  On Linux, workers share the '
                        'reference genome via copy-on-write fork at no memory cost. '
                        'Each worker loads ~6–7 GB of BAM data, so limit workers '
                        'to available_RAM_GB / 7.  E.g. --workers 8 on a 192 GB node.')
    p.add_argument('--str-repeat', action='store_true', default=False,
                   help='Separately profile dinucleotide and trinucleotide repeat (STR) '
                        'errors (e.g. TATATA→TATA).  Positions with hp_run_length=1 that '
                        'fall within an STR of unit length 2–3 with ≥3 copies are counted '
                        'in a separate str_error_counts.tsv keyed by (unit, n_copies) '
                        'instead of (ref_base, hp_len).  This prevents STR slippage errors '
                        'from inflating the HP=1 baseline used to normalize the penalty '
                        'table.  The HP penalty table from the main output will reflect '
                        'true non-STR HP=1 error rates after this correction.')
    p.add_argument('--umi-tag', type=str, default='XC', metavar='TAG',
                   help='BAM aux tag carrying the integer UMI cluster size '
                        '(default: XC). When the tag is present, reads are '
                        'routed into per-bin error-count tables (see '
                        '--umi-bins) in addition to the pooled output. Reads '
                        'missing this tag still contribute to the pooled '
                        'counts but to no per-bin output.')
    p.add_argument('--umi-bins', type=str, default='1,2,3+', metavar='SPEC',
                   help='Comma-separated UMI cluster-size bins for stratified '
                        "calibration (default: '1,2,3+'). Each token is "
                        "either '<N>' (exact match) or '<N>+' (≥ N). For "
                        'each bin the profiler emits '
                        "error_counts_umi<label>.tsv alongside the pooled "
                        "error_counts.tsv. Pass '1' to disable stratification "
                        '(single bin matching only singletons).')
    p.add_argument('--strand-split', action='store_true', default=False,
                   help='In addition to the pooled counts, also accumulate per-strand '
                        'counts keyed by (op_type, ref_base, hp_length, strand).  Outputs '
                        'strand_error_rates_by_base.tsv with an extra "strand" (+/-) column.  '
                        'Useful for verifying that + and - strand reads show symmetric error '
                        'profiles (a significant asymmetry would suggest strand-specific '
                        'sequencing artifacts or systematic misalignment).')
    p.add_argument('--umi-source-bam', type=Path, default=None, metavar='PATH',
                   help='Sidecar BAM carrying the UMI tag (e.g. stage1_consensus.bam) '
                        'when the per-aligner BAMs lack it. At startup, read names and '
                        'their --umi-tag values are loaded from this BAM into an in-memory '
                        'lookup and used instead of querying the aligned BAMs. Falls back '
                        'to the aligned BAM tag for reads not found in the sidecar. '
                        'Read names are normalised the same way as per-aligner BAMs.')

    return p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _worker(args_tuple):
    """Multiprocessing worker: profile one chunk and return local count dicts.

    Returns a tuple of:
        (counts, union_counts, str_counts, strand_counts,
         counts_by_bin, union_counts_by_bin, str_counts_by_bin, reads_per_bin)

    Each *_counts is a plain dict (or None when the corresponding flag is off);
    each *_by_bin is a dict of plain dicts keyed by bin label.

    Uses process-local defaultdicts so no shared state is needed. On Linux
    (fork start method) ref_seqs is shared copy-on-write at no cost. The
    `umi_bins` spec is a list of plain tuples — no lambdas — so it pickles
    cleanly under both fork and spawn start methods.
    """
    (chunk_idx, run_dir, n_chunks, aligners, ref_seqs,
     min_aligners, max_reads, isolation_flank, exclude_5p_bp, exclude_3p_bp,
     do_union, do_str, do_strand, umi_tag, umi_bins, umi_lookup,
     usable_mask) = args_tuple

    local_counts: Dict[Tuple[str, str, int], int] = defaultdict(int)
    local_union: Optional[Dict[Tuple[str, str, int], int]] = (
        defaultdict(int) if do_union else None
    )
    local_str: Optional[Dict[Tuple[str, str, int], int]] = (
        defaultdict(int) if do_str else None
    )
    local_strand: Optional[Dict[Tuple[str, str, int, str], int]] = (
        defaultdict(int) if do_strand else None
    )

    if umi_bins:
        local_counts_by_bin: Optional[Dict[str, Dict[Tuple[str, str, int], int]]] = {
            label: defaultdict(int) for label, _, _ in umi_bins
        }
        local_union_by_bin: Optional[Dict[str, Dict[Tuple[str, str, int], int]]] = (
            {label: defaultdict(int) for label, _, _ in umi_bins} if do_union else None
        )
        local_str_by_bin: Optional[Dict[str, Dict[Tuple[str, str, int], int]]] = (
            {label: defaultdict(int) for label, _, _ in umi_bins} if do_str else None
        )
        local_reads_per_bin: Optional[Dict[str, int]] = {label: 0 for label, _, _ in umi_bins}
    else:
        local_counts_by_bin = None
        local_union_by_bin = None
        local_str_by_bin = None
        local_reads_per_bin = None

    chunk_bams = find_chunk_bams(run_dir, chunk_idx, n_chunks, aligners)
    if not chunk_bams:
        log.warning('Worker: no BAMs for chunk %d — skipping', chunk_idx)
        return dict(local_counts), None, None, None, None, None, None, None

    log.info('[worker pid=%d] chunk %d: %s', os.getpid(), chunk_idx,
             list(chunk_bams.keys()))
    profile_chunk(chunk_bams, ref_seqs, min_aligners, max_reads, local_counts,
                  isolation_flank=isolation_flank,
                  exclude_5p_bp=exclude_5p_bp,
                  exclude_3p_bp=exclude_3p_bp,
                  union_counts=local_union,
                  str_counts=local_str,
                  strand_counts=local_strand,
                  umi_tag=umi_tag,
                  umi_bins=umi_bins,
                  counts_by_bin=local_counts_by_bin,
                  union_counts_by_bin=local_union_by_bin,
                  str_counts_by_bin=local_str_by_bin,
                  reads_per_bin=local_reads_per_bin,
                  umi_lookup=umi_lookup,
                  usable_mask=usable_mask)

    def _freeze(d):
        return {k: dict(v) for k, v in d.items()} if d is not None else None

    return (dict(local_counts),
            dict(local_union)  if local_union  is not None else None,
            dict(local_str)    if local_str    is not None else None,
            dict(local_strand) if local_strand is not None else None,
            _freeze(local_counts_by_bin),
            _freeze(local_union_by_bin),
            _freeze(local_str_by_bin),
            dict(local_reads_per_bin) if local_reads_per_bin is not None else None)


def _merge_counts(
    target: Dict[Tuple[str, str, int], int],
    source: Optional[Dict[Tuple[str, str, int], int]],
) -> None:
    """Add all values from *source* into *target* in place."""
    if source:
        for k, v in source.items():
            target[k] += v


def _merge_counts_by_bin(
    target: Optional[Dict[str, Dict[Tuple[str, str, int], int]]],
    source: Optional[Dict[str, Dict[Tuple[str, str, int], int]]],
) -> None:
    """Add per-bin source dicts into per-bin target in place."""
    if not target or not source:
        return
    for bin_label, sub in source.items():
        if bin_label not in target:
            target[bin_label] = defaultdict(int)
        for k, v in sub.items():
            target[bin_label][k] += v


def load_usable_mask(
    ref_seqs: Dict[str, str],
    vcf_path: Optional[Path],
    bed_path: Optional[Path],
) -> Optional[Dict[str, np.ndarray]]:
    """Build per-chrom boolean 'usable' arrays from an exclude-VCF and/or a
    high-conf BED. Returns None when neither is given (no masking).

    Memory: one bool array per chrom (1 byte/base, ~3.1 GB for GRCh38). Built
    once in the parent; on Linux (fork) it is shared copy-on-write with --run-dir
    workers, but is also pickled per chunk — so for large genomes prefer the
    in-process --aligner-bams path (no pickling)."""
    if vcf_path is None and bed_path is None:
        return None
    default_val = bed_path is None  # BED given => start excluded; else start usable
    mask: Dict[str, np.ndarray] = {
        chrom: np.full(len(seq), default_val, dtype=bool)
        for chrom, seq in ref_seqs.items()
    }
    if bed_path is not None:
        n_iv = 0
        with open(bed_path) as fh:
            for line in fh:
                if not line or line[0] == '#':
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 3:
                    continue
                arr = mask.get(f[0])
                if arr is not None:
                    arr[int(f[1]):int(f[2])] = True
                    n_iv += 1
        log.info('  high-conf BED: %d intervals loaded', n_iv)
    if vcf_path is not None:
        vf = pysam.VariantFile(str(vcf_path))
        n_var = 0
        for rec in vf.fetch():
            arr = mask.get(rec.chrom)
            if arr is not None:
                reflen = len(rec.ref) if rec.ref else 1
                arr[rec.start:rec.start + reflen] = False
                n_var += 1
        log.info('  exclude-VCF: %d records masked', n_var)
    usable_bp = sum(int(a.sum()) for a in mask.values())
    log.info('  usable positions after masking: %.1f Mb', usable_bp / 1e6)
    return mask


def main(argv: Optional[List[str]] = None) -> None:
    args = build_parser().parse_args(argv)

    # --- Resolve protocol-driven defaults ---
    if args.aligners is None:
        args.aligners = list(_DEFAULT_ALIGNERS_BY_PROTOCOL[args.protocol])
    default_5p, default_3p = _PROTOCOL_EXCLUSION_DEFAULTS[args.protocol]
    if args.exclude_5p_bp is None:
        args.exclude_5p_bp = default_5p
    if args.exclude_3p_bp is None:
        args.exclude_3p_bp = default_3p
    log.info('Protocol: %s  |  aligners: %s  |  read-coord exclusion: 5\'=%d bp, 3\'=%d bp',
             args.protocol, ' '.join(args.aligners),
             args.exclude_5p_bp, args.exclude_3p_bp)

    # --- Parse UMI bins (plain-data tuples; pickles cleanly for workers) ---
    try:
        args.umi_bin_specs: List[UmiBinSpec] = parse_umi_bins(args.umi_bins)
    except ValueError as exc:
        sys.exit(f'ERROR: invalid --umi-bins {args.umi_bins!r}: {exc}')
    log.info('UMI tag: %s  |  bins: %s',
             args.umi_tag,
             ', '.join(f'{lbl}({kind} {n})' for lbl, kind, n in args.umi_bin_specs))

    # --- Build sidecar UMI lookup (optional) ---
    umi_lookup: Optional[Dict[str, int]] = None
    if args.umi_source_bam is not None:
        if not args.umi_source_bam.exists():
            sys.exit(f'ERROR: --umi-source-bam not found: {args.umi_source_bam}')
        log.info('Building UMI lookup from sidecar BAM: %s', args.umi_source_bam)
        umi_lookup = {}
        with pysam.AlignmentFile(str(args.umi_source_bam), 'rb', check_sq=False) as _bam:
            for _r in _bam:
                if _r.is_unmapped or _r.is_secondary or _r.is_supplementary:
                    continue
                try:
                    _val = int(_r.get_tag(args.umi_tag))
                except (KeyError, TypeError, ValueError):
                    continue
                umi_lookup[_normalise_read_name(_r.query_name)] = _val
        log.info('Sidecar UMI lookup built: %d entries', len(umi_lookup))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = args.output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)

    # --- Reference genome ---
    log.info('Loading reference genome: %s', args.reference)
    fasta = pysam.FastaFile(str(args.reference))
    ref_seqs: Dict[str, str] = {}
    for chrom in fasta.references:
        ref_seqs[chrom] = fasta.fetch(chrom).upper()
    log.info('Loaded %d chromosomes (%.1f Mb total)',
             len(ref_seqs),
             sum(len(s) for s in ref_seqs.values()) / 1e6)

    # --- Variant / high-conf position mask (optional; diploid/outbred genomes) ---
    usable_mask = load_usable_mask(ref_seqs, args.exclude_vcf, args.high_conf_bed)
    if usable_mask is not None:
        log.info('Variant/high-conf masking ACTIVE (exclude-vcf=%s, high-conf-bed=%s)',
                 args.exclude_vcf, args.high_conf_bed)

    # --- Accumulators ---
    counts: Dict[Tuple[str, str, int], int] = defaultdict(int)
    union_counts: Optional[Dict[Tuple[str, str, int], int]] = (
        defaultdict(int) if args.union else None
    )
    str_counts: Optional[Dict[Tuple[str, str, int], int]] = (
        defaultdict(int) if args.str_repeat else None
    )
    strand_counts: Optional[Dict[Tuple[str, str, int, str], int]] = (
        defaultdict(int) if args.strand_split else None
    )

    # Per-bin accumulators (parallel to pooled counts/union/str above).
    # strand_counts is intentionally pooled-only — strand symmetry analysis
    # doesn't need UMI stratification at this stage.
    counts_by_bin: Dict[str, Dict[Tuple[str, str, int], int]] = {
        label: defaultdict(int) for label, _, _ in args.umi_bin_specs
    }
    union_counts_by_bin: Optional[Dict[str, Dict[Tuple[str, str, int], int]]] = (
        {label: defaultdict(int) for label, _, _ in args.umi_bin_specs}
        if args.union else None
    )
    str_counts_by_bin: Optional[Dict[str, Dict[Tuple[str, str, int], int]]] = (
        {label: defaultdict(int) for label, _, _ in args.umi_bin_specs}
        if args.str_repeat else None
    )
    reads_per_bin: Dict[str, int] = {label: 0 for label, _, _ in args.umi_bin_specs}

    if args.isolation_flank > 0:
        log.info('Isolation filter: --isolation-flank %d  (require ≥%d M-flanking bases)',
                 args.isolation_flank, args.isolation_flank)
    if args.union:
        log.info('Union mode enabled: accumulating any-aligner error events alongside strict')
    if args.str_repeat:
        log.info('STR profiling enabled: dinucleotide/trinucleotide repeats counted separately')

    # --- Route: explicit --aligner-bams vs --run-dir ---
    profile_kwargs = dict(
        isolation_flank=args.isolation_flank,
        exclude_5p_bp=args.exclude_5p_bp,
        exclude_3p_bp=args.exclude_3p_bp,
        union_counts=union_counts,
        str_counts=str_counts,
        strand_counts=strand_counts,
        umi_tag=args.umi_tag,
        umi_bins=args.umi_bin_specs,
        counts_by_bin=counts_by_bin,
        union_counts_by_bin=union_counts_by_bin,
        str_counts_by_bin=str_counts_by_bin,
        reads_per_bin=reads_per_bin,
        umi_lookup=umi_lookup,
        usable_mask=usable_mask,
    )

    if args.aligner_bams:
        # Parse "aligner:path" pairs
        chunk_bams: Dict[str, Path] = {}
        for item in args.aligner_bams:
            if ':' not in item:
                sys.exit(f'ERROR: --aligner-bams entries must be "aligner:path", got: {item!r}')
            aligner, path_str = item.split(':', 1)
            chunk_bams[aligner] = Path(path_str)

        min_aligners = args.min_aligners or len(chunk_bams)
        log.info('Processing single BAM set with %d aligners (min_aligners=%d)',
                 len(chunk_bams), min_aligners)
        profile_chunk(chunk_bams, ref_seqs, min_aligners, args.max_reads, counts,
                      **profile_kwargs)

    else:
        # Auto-discover chunks from --run-dir
        run_dir = args.run_dir
        n_chunks = discover_chunk_count(run_dir)
        chunk_indices = args.chunks if args.chunks is not None else list(range(n_chunks))

        log.info('Run dir: %s', run_dir)
        log.info('Total chunks: %d  |  Processing: %s', n_chunks, chunk_indices)

        min_aligners = args.min_aligners or len(args.aligners)

        worker_args = [
            (chunk_idx, run_dir, n_chunks, args.aligners, ref_seqs,
             min_aligners, args.max_reads, args.isolation_flank,
             args.exclude_5p_bp, args.exclude_3p_bp,
             args.union, args.str_repeat, args.strand_split,
             args.umi_tag, args.umi_bin_specs, umi_lookup, usable_mask)
            for chunk_idx in chunk_indices
        ]

        n_workers = min(args.workers, len(chunk_indices))
        if n_workers > 1:
            log.info('Parallel mode: %d workers for %d chunks', n_workers, len(chunk_indices))
            with multiprocessing.Pool(processes=n_workers) as pool:
                results = pool.map(_worker, worker_args)
        else:
            log.info('Sequential mode (--workers 1)')
            results = [_worker(a) for a in worker_args]

        for (chunk_counts, chunk_union, chunk_str, chunk_strand,
             chunk_counts_by_bin, chunk_union_by_bin, chunk_str_by_bin,
             chunk_reads_per_bin) in results:
            _merge_counts(counts, chunk_counts)
            if union_counts is not None:
                _merge_counts(union_counts, chunk_union)
            if str_counts is not None:
                _merge_counts(str_counts, chunk_str)
            if strand_counts is not None:
                _merge_counts(strand_counts, chunk_strand)
            _merge_counts_by_bin(counts_by_bin, chunk_counts_by_bin)
            if union_counts_by_bin is not None:
                _merge_counts_by_bin(union_counts_by_bin, chunk_union_by_bin)
            if str_counts_by_bin is not None:
                _merge_counts_by_bin(str_counts_by_bin, chunk_str_by_bin)
            if chunk_reads_per_bin:
                for k, v in chunk_reads_per_bin.items():
                    reads_per_bin[k] = reads_per_bin.get(k, 0) + v

    if not counts:
        log.error('No counts accumulated — check inputs.')
        sys.exit(1)

    log.info('Total strict agreed events: %d', sum(counts.values()))
    if union_counts:
        log.info('Total union events:         %d', sum(union_counts.values()))
    if str_counts:
        log.info('Total STR events:           %d', sum(str_counts.values()))

    # --- Helper to write one full set of outputs ---
    def _write_outputs(
        cnt: Dict[Tuple[str, str, int], int],
        prefix: str,
        label: str,
    ) -> None:
        log.info('--- %s outputs ---', label)
        rate_df = compute_rates(cnt)
        rate_df = rate_df[rate_df['hp_length'] <= args.max_hp]

        pen_df = derive_penalties(
            rate_df,
            default_ins=args.default_ins,
            min_count_flag=args.min_count_flag,
        )

        # Raw counts
        cnt_rows = [
            {'op_type': op, 'ref_base': base, 'hp_length': hp, 'count': c}
            for (op, base, hp), c in cnt.items()
            if hp <= args.max_hp
        ]
        cnt_path = args.output_dir / f'{prefix}error_counts.tsv'
        pd.DataFrame(cnt_rows).sort_values(
            ['op_type', 'hp_length', 'ref_base']
        ).to_csv(cnt_path, sep='\t', index=False)
        log.info('Wrote: %s', cnt_path)

        rate_path = args.output_dir / f'{prefix}error_rates.tsv'
        rate_df.to_csv(rate_path, sep='\t', index=False, float_format='%.6f')
        log.info('Wrote: %s', rate_path)

        base_path = args.output_dir / f'{prefix}error_rates_by_base.tsv'
        rate_df[rate_df['op_type'].isin({'D', 'X', 'I'})].to_csv(
            base_path, sep='\t', index=False, float_format='%.6f')
        log.info('Wrote: %s', base_path)

        pen_path = args.output_dir / f'{prefix}penalty_scores.tsv'
        pen_df.to_csv(pen_path, sep='\t', index=False, float_format='%.4f')
        log.info('Wrote: %s', pen_path)

        return rate_df, pen_df

    # --- Rates & penalties: strict ---
    log.info('Computing rates ...')
    rate_df_strict, penalty_df_strict = _write_outputs(counts, '', 'Strict')

    # --- Per-bin error counts (Phase A: counts only, no derived penalties) ---
    # bundle_protocol_tables.py (Phase B) consumes these to derive per-bin
    # penalty tables. We deliberately do NOT call _write_outputs() per bin
    # because that path also re-derives penalties / rates / plots — the
    # design decision is that bin-aware penalty derivation belongs in the
    # bundler, where the min_count fallback logic across bins lives.
    log.info('--- Per-UMI-bin error counts ---')
    print('\nUMI-bin summary  (--umi-tag %s):' % args.umi_tag)
    for bin_label, _kind, _n in args.umi_bin_specs:
        bin_counts = counts_by_bin.get(bin_label, {})
        n_events = sum(bin_counts.values())
        n_reads = reads_per_bin.get(bin_label, 0)
        suffix = _sanitize_bin_label(bin_label)
        out_path = args.output_dir / f'error_counts_umi{suffix}.tsv'
        cnt_rows = [
            {'op_type': op, 'ref_base': base, 'hp_length': hp, 'count': c}
            for (op, base, hp), c in bin_counts.items()
            if hp <= args.max_hp
        ]
        df_out = pd.DataFrame(
            cnt_rows,
            columns=['op_type', 'ref_base', 'hp_length', 'count'],
        )
        if not df_out.empty:
            df_out = df_out.sort_values(['op_type', 'hp_length', 'ref_base'])
        df_out.to_csv(out_path, sep='\t', index=False)
        log.info('Wrote: %s  (n_reads=%d, n_events=%d)', out_path, n_reads, n_events)
        print(f'  umi{suffix:<7}  n_reads={n_reads:>8d}  n_events={n_events:>10d}'
              f'  -> {out_path.name}')
    print()

    # --- Plots (strict) ---
    log.info('Generating plots ...')
    plot_error_rates(penalty_df_strict, plots_dir)
    plot_hp_length_distribution(rate_df_strict, plots_dir)

    # --- Console summary: strict ---
    print_penalty_table(penalty_df_strict, max_hp=min(args.max_hp, 15))

    # --- Rates & penalties: union (if requested) ---
    if union_counts:
        rate_df_union, penalty_df_union = _write_outputs(
            union_counts, 'union_', 'Union (any-aligner)')

        # --- Comparison table ---
        _print_comparison_table(penalty_df_strict, penalty_df_union,
                                max_hp=min(args.max_hp, 12))

    # --- STR error rates (if requested) ---
    if str_counts:
        _write_str_outputs(str_counts, args.output_dir, args.max_hp)

    # --- Strand-split error rates (if requested) ---
    if strand_counts:
        log.info('--- Strand-split outputs ---')
        strand_rows = []
        for (op, base, hp, strand), cnt in strand_counts.items():
            if hp <= args.max_hp:
                strand_rows.append({
                    'op_type': op, 'ref_base': base,
                    'hp_length': hp, 'strand': strand, 'count': cnt,
                })
        sc_df = pd.DataFrame(strand_rows)
        # Denominator: M + D + X per (ref_base, hp_length, strand)
        denom_ops = {'M', 'D', 'X'}
        denom_df = (
            sc_df[sc_df['op_type'].isin(denom_ops)]
            .groupby(['ref_base', 'hp_length', 'strand'], as_index=False)['count']
            .sum()
            .rename(columns={'count': 'total_positions'})
        )
        sc_df = sc_df.merge(denom_df, on=['ref_base', 'hp_length', 'strand'], how='left')
        sc_df['total_positions'] = sc_df['total_positions'].fillna(0).astype(int)
        sc_df['rate'] = sc_df['count'] / sc_df['total_positions'].clip(lower=1)
        sc_df = sc_df.sort_values(['op_type', 'strand', 'hp_length', 'ref_base'])

        strand_path = args.output_dir / 'strand_error_rates_by_base.tsv'
        sc_df[sc_df['op_type'].isin({'D', 'X', 'I'})].to_csv(
            strand_path, sep='\t', index=False, float_format='%.6f')
        log.info('Wrote: %s', strand_path)

    log.info('Done. Results in %s', args.output_dir)


if __name__ == '__main__':
    main()
