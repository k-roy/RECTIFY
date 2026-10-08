#!/usr/bin/env python3
"""
Literal read bases for the multialigned BAM.

Every aligner BAM in the panel goes through ``samtools calmd -e`` (the overhang
resolver counts mismatches as non-``=`` bases, so that encoding stays), and the
consensus record keeps its arm's SEQ. Decoding ``=`` from the final record is
not exact, because the record is not always the one ``calmd`` saw:

- uLTRA can write a record's SEQ in the orientation opposite to its own flag;
  ``calmd`` then marks chance matches and the record's bases, and its NM, are
  wrong (e.g. NM 759 on a 944-nt read; 87 with the read oriented by its flag).
- later CIGAR edits can leave ``=`` outside the aligned blocks, e.g. inside a
  soft clip, where no reference base defines it.

The input reads are the ground truth: every aligner was given exactly those
bases. ``restore_literal_sequences`` rewrites each record's SEQ/QUAL from the
input FASTQ, oriented by the record's own flag (hard clips respected), after
checking it against every base the record spells out. Alignments, flags and
tags are untouched, so consensus selection and the resolver see nothing new.

Author: Kevin R. Roy
Date: 2026-10-08
"""

import gzip
import logging
import os
import subprocess
import zlib
from collections import Counter
from typing import Dict, Optional, Tuple

import numpy as np
import pysam

logger = logging.getLogger(__name__)
_normalize = None   # consensus._normalize_bam_read_name, bound on first use (avoids an import cycle at module load)
_EQ, _N = ord('='), ord('N')

_RC = str.maketrans('ACGTNacgtnRYKMSWBDHVrykmswbdhv', 'TGCANtgcanYRMKSWVHDByrmkswvhdb')
_CIGAR_H = 5
_QUERY_OPS = (0, 1, 4, 7, 8)          # M I S = X consume query bases stored in SEQ


def _fastq_name(header: str) -> str:
    """First whitespace token of a FASTQ header, without the '@' (what aligners put in QNAME)."""
    return header[1:].split(None, 1)[0][:254] if len(header) > 1 else ''


def _bam_name(name: str) -> str:
    """QNAME as the consensus writes it (aligner comment suffixes stripped)."""
    global _normalize
    if _normalize is None:
        from ..consensus.consensus import _normalize_bam_read_name
        _normalize = _normalize_bam_read_name
    return _normalize(name or '')


def _spelled_agree(stored: str, candidate: str) -> bool:
    """True when every base ``stored`` spells out (not ``=``, not ``N``) equals ``candidate`` at that position."""
    if stored == candidate:
        return True
    a = np.frombuffer(stored.encode(), dtype=np.uint8)
    b = np.frombuffer(candidate.encode(), dtype=np.uint8)
    m = (a != _EQ) & (a != _N)
    return bool(np.array_equal(a[m], b[m]))


def _bucket(name: str, n_buckets: int) -> int:
    return zlib.crc32(name.encode()) % n_buckets if n_buckets > 1 else 0


def _iter_fastq(path: str):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip('\n')
            fh.readline()
            qual = fh.readline().rstrip('\n')
            if header.startswith('@'):
                yield _fastq_name(header.rstrip('\n')), seq, qual


def rebuild_record_sequence(read: pysam.AlignedSegment, seq: str, qual: Optional[str]) -> str:
    """Set ``read``'s SEQ/QUAL from the input read, oriented by its flag. Returns a status word.

    ``seq``/``qual`` are the read as given to the aligners. Statuses: ``rebuilt`` (spelled bases agree),
    ``rebuilt_flipped`` (they agree only with the opposite orientation, i.e. the stored SEQ was flipped),
    ``unchanged_inconsistent`` (neither orientation agrees: left as is), ``unchanged_length`` (CIGAR, with
    hard clips, does not span the read: left as is).
    """
    stored = read.query_sequence or ''
    rev = read.is_reverse
    oriented = seq.translate(_RC)[::-1] if rev else seq
    oq = (qual[::-1] if rev else qual) if qual else None
    if read.is_unmapped or not read.cigartuples:
        lo, hi = 0, len(oriented)
    else:
        ct = read.cigartuples
        h_left = ct[0][1] if ct[0][0] == _CIGAR_H else 0
        h_right = ct[-1][1] if len(ct) > 1 and ct[-1][0] == _CIGAR_H else 0
        qlen = sum(n for op, n in ct if op in _QUERY_OPS)
        if h_left + qlen + h_right != len(oriented):
            return 'unchanged_length'
        lo, hi = h_left, h_left + qlen
    new = oriented[lo:hi]
    if stored and len(stored) == len(new):
        if _spelled_agree(stored, new):
            status = 'rebuilt'
        elif _spelled_agree(stored, (seq if rev else seq.translate(_RC)[::-1])[lo:hi]):
            status = 'rebuilt_flipped'
        else:
            return 'unchanged_inconsistent'
    elif stored and len(stored) != len(new):
        return 'unchanged_length'
    else:
        status = 'rebuilt'
    read.query_sequence = new
    if oq is not None and len(oq) == len(oriented):
        read.query_qualities = pysam.qualitystring_to_array(oq[lo:hi])
    return status


def restore_literal_sequences(bam_path: str, reads_path: str, threads: int = 1,
                              max_reads_in_memory: int = 1_000_000,
                              max_inconsistent_fraction: float = 0.001) -> Dict[str, int]:
    """Rewrite every record of ``bam_path`` (in place, order kept, re-indexed) with literal bases from ``reads_path``.

    Reads are held in memory in name-hash buckets of at most ``max_reads_in_memory`` reads; with more than one
    bucket the BAM is read once per bucket and the output is coordinate-sorted again. Raises if more than
    ``max_inconsistent_fraction`` of records cannot be rebuilt or lack their read (the BAM is then left as it was).
    """
    n_reads = sum(1 for _ in _iter_fastq(reads_path))
    n_buckets = max(1, -(-n_reads // max(1, max_reads_in_memory)))
    tmp = f'{bam_path}.literal_tmp.bam'
    st: Counter = Counter()
    with pysam.AlignmentFile(bam_path, 'rb') as probe:
        header = probe.header.to_dict()
        coord_sorted = header.get('HD', {}).get('SO') == 'coordinate'
    header.setdefault('PG', []).append({'ID': f'rectify-literal-seq.{len(header.get("PG", []))}', 'PN': 'rectify-literal-seq',
                                        'CL': f'SEQ/QUAL rewritten from {os.path.basename(str(reads_path))}, oriented by each record flag'})
    out_path = tmp if n_buckets == 1 else f'{bam_path}.literal_unsorted.bam'
    with pysam.AlignmentFile(out_path, 'wb', header=header, threads=max(1, threads)) as out:
        for b in range(n_buckets):
            reads: Dict[str, Tuple[str, str]] = {}
            for name, seq, qual in _iter_fastq(reads_path):
                if _bucket(name, n_buckets) == b:
                    reads[name] = (seq, qual)
            with pysam.AlignmentFile(bam_path, 'rb', threads=max(1, threads)) as inp:
                for r in inp:
                    name = _bam_name(r.query_name)
                    if _bucket(name, n_buckets) != b:
                        continue
                    st['records'] += 1
                    src = reads.get(name)
                    if src is None:
                        st['no_read'] += 1
                    else:
                        st[rebuild_record_sequence(r, src[0], src[1])] += 1
                    out.write(r)
            reads = {}
    bad = st['no_read'] + st['unchanged_inconsistent']
    if st['records'] and bad / st['records'] > max_inconsistent_fraction:
        for p in (tmp, f'{bam_path}.literal_unsorted.bam'):
            if os.path.exists(p):
                os.remove(p)
        raise RuntimeError(f'literal-sequence restore: {bad} of {st["records"]} records lack their read or disagree with it ({dict(st)})')
    if n_buckets > 1:
        if coord_sorted:
            subprocess.run(['samtools', 'sort', '-@', str(max(1, threads)), '-o', tmp, out_path], check=True)
            os.remove(out_path)
        else:
            os.replace(out_path, tmp)
    os.replace(tmp, bam_path)
    if coord_sorted:
        pysam.index(bam_path)
    st['fastq_reads'] = n_reads
    st['buckets'] = n_buckets
    return dict(st)
