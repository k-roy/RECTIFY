"""Genome-wide clip mappability — clause C of the resolver v3 refusal set.

Why this module exists (planning/769d, planning/771 A3)
-------------------------------------------------------
The overhang resolver's clause B (a within-window null) refuses a placement
when some position in the resolver's own +/-``max_intron`` search window scores
better. That catches the **P02** artifact — a 161 nt Ty1-delta-LTR clip forced
into a ``180M1852N327M1S`` alignment at chrIII:149657, with ~48 positions in the
resolver's own window scoring better than the one it picked.

It structurally **cannot** catch **P07** (``chrIII:148616-151517``). P07 scores
``gain = 0.0`` under clause B: its placement really is rank-1 within +/-5 kb.
But the clip places at hp-ED **2.0** at chrIII:84298 / 84814 / 90439 — **60-67 kb
away**, far outside any window clause B can see. Both clips are Ty1-region
repeat sequence and both best-place 17-67 kb from where the resolver put them.

So clause C asks a different question: *is this clip's accepted placement its
best placement ANYWHERE in the genome?*

🔴 Uniqueness is computed on the CLIP, never on the read's MAPQ
---------------------------------------------------------------
All four artifact reads carry **MAPQ 60** on the minimap2 record — their MAPQ
reflects the clean multi-hundred-base anchor, not the clip. The MAPQ 0 quoted in
planning/769b was the clip **re-mapped on its own**. A MAPQ pre-filter on the
read was tried and withdrawn (planning/769d correction 1); do not reintroduce it.

Algorithm — pigeonhole seed-and-extend, NOT variant enumeration
----------------------------------------------------------------
Enumerating ed<=2 variants of a 20-mer is ~1,770 substitutions alone and scales
``O(X^e * 3^e)`` — it gets *worse* as X grows, and X is the knob we want to
raise. Instead the genome's K-mers are indexed once (bucket-sorted by 2-bit
code), seeds are probed from the query, and surviving candidate offsets are
verified with the resolver's own hp-edit distance. Lookups are O(1) and constant
in query length.

🔴 Dust/entropy-filter every seed and cap its genome frequency. A poly-A 10-mer
hits everywhere; without this the Ty1 trap comes straight back in a new costume.

The same index is what arm 1's spec calls for, and clip mappability comes free
with it (planning/770).
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

_FORMAT_VERSION = 1

#: K-mer length. 4**10 = 1,048,576 buckets — ~11 positions per bucket on a 12 Mb
#: yeast genome, and an 8 MB offset table. K=12 would be 16.7M buckets / 134 MB
#: of offsets for no sensitivity gain at this genome size.
DEFAULT_K = 10

#: A seed occurring more than this many times genome-wide is dropped as
#: uninformative (the poly-A / Ty1 trap). Scaled from the genome length at build
#: time; this is the floor.
FREQ_CAP = 200

#: Candidate offsets to verify with the DP, ranked by seed-vote count.
MAX_CANDIDATES = 60

_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b'ACGT'):
    _LUT[_b] = _i
    _LUT[_b + 32] = _i          # lower case


def _seed_is_informative(seed: str) -> bool:
    """Reject low-complexity seeds before they are ever looked up.

    Two cheap tests that between them kill the failure mode: a homopolymer-
    dominated seed (poly-A tails, A-rich yeast intergenic) and a seed built from
    too few distinct bases (dinucleotide repeats, (CA)n / (GT)n microsatellite).
    """
    n = len(seed)
    if n == 0:
        return False
    counts: Dict[str, int] = {}
    for c in seed:
        counts[c] = counts.get(c, 0) + 1
    if len(counts) < 3:
        return False
    if max(counts.values()) / n >= 0.7:
        return False
    return True


class KmerIndex:
    """Bucket-sorted genome K-mer index over the concatenated contigs.

    ``starts[c]`` .. ``starts[c+1]`` indexes ``order`` for K-mer code ``c``;
    ``order`` holds GLOBAL positions in the concatenated sequence, which
    :meth:`locate` maps back to ``(chrom, pos)``.
    """

    def __init__(self, k: int, order: np.ndarray, starts: np.ndarray,
                 names: List[str], offsets: np.ndarray, lengths: np.ndarray,
                 fingerprint: Optional[dict] = None):
        self.k = int(k)
        self.order = order
        self.starts = starts
        self.names = list(names)
        self.offsets = offsets          # global start of each contig
        self.lengths = lengths
        self.fingerprint = fingerprint
        self._name_to_i = {n: i for i, n in enumerate(self.names)}

    # ------------------------------------------------------------------ build
    @classmethod
    def build(cls, genome: Dict[str, str], k: int = DEFAULT_K,
              fingerprint: Optional[dict] = None) -> 'KmerIndex':
        names = sorted(genome)
        lengths = np.array([len(genome[n]) for n in names], dtype=np.int64)
        offsets = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)
        total = int(offsets[-1])
        # A separator run of invalid bases keeps K-mers from spanning contigs.
        buf = np.empty(total, dtype=np.uint8)
        for i, n in enumerate(names):
            seq = genome[n].encode('ascii', 'replace')
            buf[offsets[i]:offsets[i + 1]] = np.frombuffer(seq, dtype=np.uint8)
        val = _LUT[buf]
        n_pos = total - k + 1
        if n_pos <= 0:
            raise ValueError('genome shorter than K')

        code = np.zeros(n_pos, dtype=np.int64)
        for i in range(k):
            code = code * 4 + val[i:i + n_pos].astype(np.int64)
        # validity: no non-ACGT inside the window, and no contig boundary
        bad = (val >= 4).astype(np.int32)
        cbad = np.concatenate([[0], np.cumsum(bad)])
        valid = (cbad[k:] - cbad[:-k]) == 0
        # a K-mer must lie wholly inside one contig
        starts_in = np.searchsorted(offsets, np.arange(n_pos), side='right') - 1
        ends_in = np.searchsorted(offsets, np.arange(n_pos) + k - 1, side='right') - 1
        valid &= (starts_in == ends_in)

        pos = np.nonzero(valid)[0].astype(np.int64)
        codes = code[pos]
        srt = np.argsort(codes, kind='stable')
        order = pos[srt].astype(np.int64)
        counts = np.bincount(codes[srt], minlength=4 ** k)
        starts = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
        logger.info('clip-mappability K-mer index: k=%d, %d valid positions, '
                    '%d contigs', k, order.size, len(names))
        return cls(k, order, starts, names, offsets[:-1], lengths, fingerprint)

    # ------------------------------------------------------------------ cache
    @staticmethod
    def _fingerprint_for(genome_path: str, k: int) -> dict:
        st = os.stat(genome_path)
        return {
            'format_version': _FORMAT_VERSION,
            'k': int(k),
            'genome_realpath': os.path.realpath(genome_path),
            'size': st.st_size,
            'mtime_ns': st.st_mtime_ns,
        }

    @classmethod
    def load_or_build(cls, genome_path: str, genome: Dict[str, str],
                      k: int = DEFAULT_K,
                      cache_path: Optional[str] = None) -> 'KmerIndex':
        """Same contract as :meth:`SpliceSiteIndex.load_or_build` — an
        unwritable cache degrades to build-only, never an error."""
        cache = Path(cache_path) if cache_path else Path(f'{genome_path}.kmer{k}.npz')
        want = cls._fingerprint_for(str(genome_path), k)
        if cache.exists():
            try:
                npz = np.load(cache, allow_pickle=False)
                got = json.loads(str(npz['__fingerprint__'][0]))
                if got == want:
                    return cls(k, npz['order'], npz['starts'],
                               [str(x) for x in npz['names']],
                               npz['offsets'], npz['lengths'], want)
                logger.info('clip-mappability cache stale; rebuilding: %s', cache)
            except Exception as exc:
                logger.warning('clip-mappability cache unreadable (%s); '
                               'rebuilding: %s', exc, cache)
        idx = cls.build(genome, k=k, fingerprint=want)
        try:
            # PID-suffixed: concurrent array tasks on a cold cache would
            # otherwise race one another on a single fixed .tmp path.
            tmp = cache.with_name(f'{cache.name}.{os.getpid()}.tmp')
            with open(tmp, 'wb') as fh:
                np.savez(fh, order=idx.order, starts=idx.starts,
                         names=np.array(idx.names),
                         offsets=idx.offsets, lengths=idx.lengths,
                         __fingerprint__=np.array([json.dumps(want)]))
            os.replace(tmp, cache)
        except Exception as exc:
            logger.warning('clip-mappability cache not written (%s): %s', exc, cache)
        return idx

    # ------------------------------------------------------------------ query
    def _encode(self, seed: str) -> int:
        c = 0
        for ch in seed:
            v = _LUT[ord(ch) & 0xFF]
            if v >= 4:
                return -1
            c = c * 4 + int(v)
        return c

    def positions(self, seed: str) -> np.ndarray:
        c = self._encode(seed)
        if c < 0:
            return np.empty(0, dtype=np.int64)
        return self.order[self.starts[c]:self.starts[c + 1]]

    def locate(self, gpos: int) -> Tuple[str, int]:
        i = int(np.searchsorted(self.offsets, gpos, side='right') - 1)
        return self.names[i], int(gpos - self.offsets[i])

    def global_pos(self, chrom: str, pos: int) -> int:
        return int(self.offsets[self._name_to_i[chrom]]) + int(pos)


def best_out_of_window(
    index: KmerIndex,
    genome: Dict[str, str],
    query: str,
    chrom: str,
    exclude_lo: int,
    exclude_hi: int,
    accepted_ed: float,
    margin: float,
    scorer,
    freq_cap: int = FREQ_CAP,
    max_candidates: int = MAX_CANDIDATES,
) -> Tuple[float, Optional[str], int]:
    """Best placement of ``query`` OUTSIDE ``chrom[exclude_lo:exclude_hi]``.

    Returns ``(gain, chrom, pos)`` where ``gain = accepted_ed - best_ed``, or
    ``(0.0, None, -1)`` when nothing beats ``accepted_ed - margin``. ``scorer``
    is the edit-distance function the gate itself uses (``hp_edit_distance_
    bounded``), so the comparison is apples to apples with the accepted ED.
    """
    cutoff = accepted_ed - margin
    if cutoff < 0 or len(query) < index.k:
        return 0.0, None, -1
    k = index.k
    votes: Dict[int, int] = {}
    n = len(query)
    # Stride by k for non-overlapping seeds; a clip at hp-ED 13-35 has no
    # pigeonhole guarantee, so sensitivity comes from probing many seeds rather
    # than from a proof about three of them.
    for off in range(0, n - k + 1, k):
        seed = query[off:off + k]
        if not _seed_is_informative(seed):
            continue
        hits = index.positions(seed)
        if hits.size == 0 or hits.size > freq_cap:
            continue          # 🔴 the frequency cap IS the Ty1 guard
        for gp in hits.tolist():
            start = gp - off
            if start < 0:
                continue
            votes[start] = votes.get(start, 0) + 1
    if not votes:
        return 0.0, None, -1

    lo_g = index.global_pos(chrom, max(0, exclude_lo))
    hi_g = index.global_pos(chrom, exclude_hi)
    ranked = sorted(votes.items(), key=lambda kv: -kv[1])
    best_ed: Optional[float] = None
    best_gp = -1
    seen = 0
    for start, _v in ranked:
        if lo_g <= start <= hi_g:
            continue                      # inside the window: clause B's job
        cn, cp = index.locate(start)
        seq = genome.get(cn)
        if seq is None or cp + n > len(seq):
            continue
        ed = scorer(query, seq[cp:cp + n], cutoff)
        seen += 1
        if ed <= cutoff and (best_ed is None or ed < best_ed):
            best_ed, best_gp = ed, start
            if best_ed <= 0.0:
                break
        if seen >= max_candidates:
            break
    if best_ed is None:
        return 0.0, None, -1
    cn, cp = index.locate(best_gp)
    return accepted_ed - best_ed, cn, cp
