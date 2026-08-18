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

**Extended (Prp18-class) acceptors** — Roy et al. 2023 NAR (PMID 37956322,
gkad968) measured genome-wide activation of NON-YAG 3'SS in prp18 mutants,
and the same junction classes accumulate in plain upf1Δ (NMD stabilizes the
non-productive isoforms). The alternative acceptor terminal dinucleotides are
the **BG class (TG / CG / GG)** plus the **non-G HAU class (AT)**; donors
stay GT/GC (4 of 1,833 published alt-3'SS junctions had a non-canonical
donor — the donor index never needs extension). The index therefore ALWAYS
carries two extra arrays (built unconditionally; querying is opt-in via the
``acc_plus_all`` / ``acc_minus_all`` union kinds):

===================  ==========================================  ============
``acc_plus_ext``     TG|CG|GG|AT at [e-2, e)                     e = intron
                                                                 END (excl.)
``acc_minus_ext``    CA|CG|CC|AT at [s, s+2)  (revcomp reading)  s = intron
                                                                 START
===================  ==========================================  ============

Measured price (planning/722b): acceptor density 5.83 %/bp -> 22.3 %/bp
(x4.82); measured benefit on the published junction set: enumerable
utilized alt-3'SS 47 % -> 88 % in prp18, 62 % -> 85 % in upf1Δ-only.

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

_FORMAT_VERSION = 2  # v2: + acc_plus_ext / acc_minus_ext (Prp18 classes)

_KINDS = ('don_gt_plus', 'don_gc_plus', 'acc_plus', 'don_minus', 'acc_minus',
          'acc_plus_ext', 'acc_minus_ext')


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
            # Prp18-class extended acceptors (module docstring). Plus strand:
            # transcript acceptor dinuc in {TG, CG, GG, AT} at the intron end.
            tg = _dinuc_positions(u8, 'T', 'G')
            cgg = _dinuc_positions(u8, 'G', 'G')
            at = _dinuc_positions(u8, 'A', 'T')
            cg = _dinuc_positions(u8, 'C', 'G')
            acc_ext_p = np.sort(np.concatenate([tg, cg, cgg, at])).astype(np.uint32)
            arrays[f'{chrom}|acc_plus_ext'] = (acc_ext_p + 2).astype(np.uint32)
            # Minus strand: transcript acceptor is at the intron START read in
            # revcomp, so the forward-genome dinucleotides are the revcomps
            # {CA, CG, CC, AT}, stored at the intron start like acc_minus.
            ca = _dinuc_positions(u8, 'C', 'A')
            cc = _dinuc_positions(u8, 'C', 'C')
            acc_ext_m = np.sort(np.concatenate([ca, cg, cc, at])).astype(np.uint32)
            arrays[f'{chrom}|acc_minus_ext'] = acc_ext_m
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
        acc_minus / acc_plus_ext / acc_minus_ext, or a convenience union:
        ``don_plus`` (GT+GC), ``acc_plus_all`` / ``acc_minus_all``
        (canonical AG-class + the Prp18 extended classes).
        """
        if kind == 'don_plus':
            a = self.sites_in(chrom, 'don_gt_plus', lo, hi)
            b = self.sites_in(chrom, 'don_gc_plus', lo, hi)
            return np.sort(np.concatenate([a, b]))
        if kind in ('acc_plus_all', 'acc_minus_all'):
            # (A v1 cache without the ext arrays can never be loaded here —
            # the format-version bump makes its fingerprint mismatch, which
            # forces a rebuild in load_or_build.)
            base = kind.replace('_all', '')
            a = self.sites_in(chrom, base, lo, hi)
            b = self.sites_in(chrom, base + '_ext', lo, hi)
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
