"""Read-side tri/di-nucleotide repeat-expansion artifact detection.

Nanopore DRS (pronounced on RNA004 chemistry) emits long low-period repeat
expansions — (AAG/CTT)n and similar di/tri-nucleotide runs, tens to hundreds of bp
— in read soft-clips and insertions. They arise from the sequencing motor stalling
at RNA secondary structure and slipping on a short repeat tract, NOT from a missed
upstream intron. Attempting 5'SS-rescue on such a clip is wasted work: the
low-complexity sequence is compared against every candidate junction in the pool
(~10k+ on human) and can match many spuriously. This module detects the expansion
so the rescue path can short-circuit and the read can be flagged.

Cross-organism + chemistry-specific (2026-05-26): yeast and human RNA004 DRS show
the same k=3 (AAG·CTT) dominance; R9.4 (SG-NEx A549) does not. Method matches the
canonical ``multimer_scan.py`` ``dominant_kmer_rna`` (canonical k-mer family covers
>= ``min_frac`` of windows, segment spans >= 3 full repeats).

Motif-agnostic by construction: detects repeat STRUCTURE at any period, never gates
on a specific motif (no GT-AG-style gating). Homopolymers (poly-A/poly-T, period 1)
are deliberately NOT flagged here — those are the 3' poly-A tail handled by the
poly-A walkback, so the dominant repeat motif is required to be multi-base.
"""

from __future__ import annotations

from collections import Counter
from typing import Optional, Tuple

# Defaults match the validated yeast/human profile (dev/repeat_expansion_profile.py).
DEFAULT_MIN_LEN = 30
DEFAULT_MIN_FRAC = 0.75
DEFAULT_K_RANGE: Tuple[int, ...] = (2, 3, 4, 5)  # multi-base periods; period-1 homopolymer excluded
_MIN_REPEATS = 3  # require at least 3 full repeat units in the segment


def _canonical(kmer: str) -> str:
    """Lexicographically smallest rotation (AAG, AGA, GAA -> AAG)."""
    return min(kmer[i:] + kmer[:i] for i in range(len(kmer)))


def dominant_repeat_period(
    seq: str,
    min_frac: float = DEFAULT_MIN_FRAC,
    k_range: Tuple[int, ...] = DEFAULT_K_RANGE,
    min_repeats: int = _MIN_REPEATS,
) -> Optional[Tuple[int, str, float]]:
    """Smallest period ``k`` in ``k_range`` whose dominant canonical k-mer family
    covers >= ``min_frac`` of overlapping windows, with the segment spanning >=
    ``min_repeats`` full repeats.

    Returns ``(k, canonical_motif, frac)`` or ``None``. The motif may be a
    homopolymer (e.g. ``"AA"``) when ``k_range`` includes such periods over a
    poly-A run — callers wanting multi-base expansions only should use
    :func:`is_repeat_expansion`, which rejects single-base motifs.
    """
    s = seq.upper()
    n = len(s)
    for k in k_range:
        if n < k * min_repeats:
            continue
        n_kmers = n - k + 1
        counts: Counter = Counter(s[i:i + k] for i in range(n_kmers))
        by_canon: Counter = Counter()
        for kmer, c in counts.items():
            if "N" in kmer:
                continue
            by_canon[_canonical(kmer)] += c
        if not by_canon:
            continue
        motif, cnt = by_canon.most_common(1)[0]
        frac = cnt / n_kmers
        if frac >= min_frac:
            return k, motif, frac
    return None


def is_repeat_expansion(
    seq: str,
    min_len: int = DEFAULT_MIN_LEN,
    min_frac: float = DEFAULT_MIN_FRAC,
    k_range: Tuple[int, ...] = DEFAULT_K_RANGE,
) -> bool:
    """True if ``seq`` (a soft-clip or insertion) is a multi-base low-period repeat
    expansion of at least ``min_len`` bp.

    A pure poly-A/poly-T run is rejected even though it satisfies a period-2/3
    canonical family (``"AA"``/``"AAA"``): the dominant motif must contain at least
    two distinct bases. This keeps genuine (AAG)n/(CTT)n expansions — including
    A-rich ones like ``AAAAAGAAGAAG`` — while excluding homopolymer tails and
    poly-A runs with sparse basecall errors (whose dominant family is ``"AAA"``).
    """
    if not seq or len(seq) < min_len:
        return False
    res = dominant_repeat_period(seq, min_frac, k_range)
    if res is None:
        return False
    _k, motif, _frac = res
    return len(set(motif)) > 1
