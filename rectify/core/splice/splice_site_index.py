"""Precomputed splice-site index: candidate LOOKUP, not a scan.

planning/641 SPEC §2.3: resolving an overhang becomes a binary-search range
query over candidate donor/acceptor positions within the information-bounded
window W, instead of a base-by-base DP sweep. Complexity drops from
``O(W x DP)`` to ``O(candidates_in_W x DP)``.

Coordinate conventions (0-based, half-open introns ``[start, end)``):

===============  =============================  ===========================
array            forward-genome dinucleotide    stored coordinate
===============  =============================  ===========================
``don_gt_plus``  ``GT`` at [p, p+2)             p  = intron START
``don_gc_plus``  ``GC`` at [p, p+2)             p  = intron START
``acc_plus``     ``AG`` at [e-2, e)             e  = intron END (exclusive)
``don_minus``    ``AC``|``GC`` at [e-2, e)      e  = intron END (exclusive)
``acc_minus``    ``CT`` at [s, s+2)             s  = intron START
===============  =============================  ===========================

(minus-strand introns read GT..AG on the transcript, i.e. CT..AC / CT..GC on
the forward genome.)

The index is cached beside the genome as ``<genome>.splice_sites.npz`` with a
provenance fingerprint (realpath + size + mtime_ns + format version); a
mismatch triggers a silent rebuild. Yeast builds in <1 s; large genomes load
memory-mapped.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np

logger = logging.getLogger(__name__)

_FORMAT_VERSION = 1

_KINDS = ('don_gt_plus', 'don_gc_plus', 'acc_plus', 'don_minus', 'acc_minus')


def _dinuc_positions(seq_u8: np.ndarray, a: str, b: str) -> np.ndarray:
    """Sorted positions p where seq[p] == a and seq[p+1] == b (uppercase)."""
    if seq_u8.size < 2:
        return np.empty(0, dtype=np.uint32)
    hits = (seq_u8[:-1] == ord(a)) & (seq_u8[1:] == ord(b))
    return np.flatnonzero(hits).astype(np.uint32)


class SpliceSiteIndex:
    """Sorted per-chromosome splice-site position arrays with range queries."""

    def __init__(self, arrays: Dict[str, np.ndarray], fingerprint: Optional[dict] = None):
        # keys are f"{chrom}|{kind}"
        self._arrays = arrays
        self.fingerprint = fingerprint or {}

    # -- construction -----------------------------------------------------

    @classmethod
    def build(cls, genome: Dict[str, str], fingerprint: Optional[dict] = None) -> 'SpliceSiteIndex':
        arrays: Dict[str, np.ndarray] = {}
        for chrom, seq in genome.items():
            u8 = np.frombuffer(seq.upper().encode('ascii'), dtype=np.uint8)
            gt = _dinuc_positions(u8, 'G', 'T')
            gc = _dinuc_positions(u8, 'G', 'C')
            ag = _dinuc_positions(u8, 'A', 'G')
            ac = _dinuc_positions(u8, 'A', 'C')
            ct = _dinuc_positions(u8, 'C', 'T')
            arrays[f'{chrom}|don_gt_plus'] = gt
            arrays[f'{chrom}|don_gc_plus'] = gc
            # acceptor stored as intron END (exclusive) => dinuc position + 2
            arrays[f'{chrom}|acc_plus'] = (ag + 2).astype(np.uint32)
            # minus donor at intron END: AC or GC => merge, store pos + 2
            don_m = np.sort(np.concatenate([ac, gc])).astype(np.uint32)
            arrays[f'{chrom}|don_minus'] = (don_m + 2).astype(np.uint32)
            arrays[f'{chrom}|acc_minus'] = ct
        return cls(arrays, fingerprint=fingerprint)

    @staticmethod
    def _fingerprint_for(genome_path: str) -> dict:
        st = os.stat(genome_path)
        return {
            'format_version': _FORMAT_VERSION,
            'genome_realpath': os.path.realpath(genome_path),
            'size': st.st_size,
            'mtime_ns': st.st_mtime_ns,
        }

    @classmethod
    def load_or_build(
        cls,
        genome_path: str,
        genome: Dict[str, str],
        cache_path: Optional[str] = None,
    ) -> 'SpliceSiteIndex':
        """Load the cached index if its fingerprint matches ``genome_path``;
        otherwise build from ``genome`` and write the cache (best-effort —
        an unwritable cache dir degrades to build-only, never an error)."""
        cache = Path(cache_path) if cache_path else Path(str(genome_path) + '.splice_sites.npz')
        want = cls._fingerprint_for(genome_path)
        if cache.exists():
            try:
                npz = np.load(cache, mmap_mode='r', allow_pickle=False)
                got = json.loads(str(npz['__fingerprint__'][0]))
                if got == want:
                    arrays = {k: npz[k] for k in npz.files if k != '__fingerprint__'}
                    return cls(arrays, fingerprint=got)
                logger.info('splice-site index cache stale (fingerprint mismatch); rebuilding: %s', cache)
            except Exception as exc:  # corrupt cache => rebuild
                logger.warning('splice-site index cache unreadable (%s); rebuilding: %s', exc, cache)
        idx = cls.build(genome, fingerprint=want)
        try:
            tmp = cache.with_name(cache.name + '.tmp')
            # File object, not a path: np.savez appends '.npz' to bare paths,
            # which would break the atomic-rename dance.
            with open(tmp, 'wb') as fh:
                np.savez_compressed(
                    fh, __fingerprint__=np.array([json.dumps(want)]), **idx._arrays,
                )
            os.replace(tmp, cache)
        except OSError as exc:
            logger.warning('splice-site index cache not written (%s): %s', exc, cache)
        return idx

    # -- queries ----------------------------------------------------------

    def sites_in(self, chrom: str, kind: str, lo: int, hi: int) -> np.ndarray:
        """Sorted positions of ``kind`` in ``[lo, hi)`` on ``chrom``.

        ``kind`` is one of don_gt_plus / don_gc_plus / acc_plus / don_minus /
        acc_minus, or the convenience union ``don_plus`` (GT+GC).
        """
        if kind == 'don_plus':
            a = self.sites_in(chrom, 'don_gt_plus', lo, hi)
            b = self.sites_in(chrom, 'don_gc_plus', lo, hi)
            return np.sort(np.concatenate([a, b]))
        arr = self._arrays.get(f'{chrom}|{kind}')
        if arr is None or arr.size == 0:
            return np.empty(0, dtype=np.uint32)
        i = np.searchsorted(arr, lo, side='left')
        j = np.searchsorted(arr, hi, side='left')
        return np.asarray(arr[i:j])

    def n_sites(self, chrom: str, kind: str) -> int:
        arr = self._arrays.get(f'{chrom}|{kind}')
        return 0 if arr is None else int(arr.size)

    def chroms(self) -> Iterable[str]:
        return sorted({k.split('|', 1)[0] for k in self._arrays})
