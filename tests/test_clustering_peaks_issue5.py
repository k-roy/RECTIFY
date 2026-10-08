"""Issue #5: `_identify_peaks` / `_find_valleys_between_peaks` must give the same output as the original
O(positions x peaks) versions, and finish in seconds where those took hours."""
import time

import numpy as np
import pytest

from rectify.core.analyze.clustering import (
    DEFAULT_MIN_PEAK_SEPARATION,
    _find_valleys_between_peaks,
    _identify_peaks,
)


def _identify_peaks_reference(positions, counts, min_separation):
    """The pre-#5 implementation, verbatim."""
    if len(positions) == 0:
        return []
    sorted_indices = np.argsort(counts)[::-1]
    peaks = []
    for idx in sorted_indices:
        pos = positions[idx]
        if all(abs(pos - p) >= min_separation for p in peaks):
            peaks.append(int(pos))
    peaks.sort()
    return peaks


def _find_valleys_reference(positions, counts, peaks):
    """The pre-#5 implementation, verbatim."""
    valleys = []
    for i in range(len(peaks) - 1):
        left_peak, right_peak = peaks[i], peaks[i + 1]
        mask = (positions > left_peak) & (positions < right_peak)
        between_positions, between_counts = positions[mask], counts[mask]
        if len(between_positions) == 0:
            valleys.append((left_peak + right_peak) // 2)
        else:
            valleys.append(int(between_positions[np.argmin(between_counts)]))
    return valleys


@pytest.mark.parametrize("n,span,sep,max_count,seed", [
    (1, 10, 5, 3, 0),
    (2, 3, 5, 3, 1),
    (500, 1_500_000, 5, 50, 2),
    (3000, 1_500_000, 5, 50, 3),
    (3000, 20_000, 5, 30, 4),     # dense: many separations collide
    (5000, 30_000, 12, 30, 5),
    (4000, 8_000, 5, 3, 6),       # heavy count ties
    (2000, 2_000, 1, 2, 7),       # every position adjacent, sep 1
])
def test_same_peaks_and_valleys_as_reference(n, span, sep, max_count, seed):
    rng = np.random.default_rng(seed)
    pos = np.sort(rng.choice(span, size=n, replace=False))
    cnt = rng.integers(1, max_count + 1, size=n)
    peaks = _identify_peaks(pos, cnt, sep)
    assert peaks == _identify_peaks_reference(pos, cnt, sep)
    assert _find_valleys_between_peaks(pos, cnt, peaks) == _find_valleys_reference(pos, cnt, peaks)


def test_empty_and_single():
    empty = np.array([], dtype=np.int64)
    assert _identify_peaks(empty, empty, 5) == []
    assert _find_valleys_between_peaks(empty, empty, []) == []
    one = np.array([100]); c = np.array([7])
    assert _identify_peaks(one, c, 5) == [100]
    assert _find_valleys_between_peaks(one, c, [100]) == []


def test_valleys_accept_unsorted_positions():
    rng = np.random.default_rng(11)
    pos = rng.choice(100_000, size=2000, replace=False)          # unsorted, unique counts -> no tie ambiguity
    cnt = rng.permutation(2000) + 1
    peaks = _identify_peaks(pos, cnt, 5)
    assert _find_valleys_between_peaks(pos, cnt, peaks) == _find_valleys_reference(pos, cnt, peaks)


def test_scales_to_deep_libraries():
    # 60,000 distinct positions on one strand at the default separation: the pre-#5 code took ~3.5 min here.
    rng = np.random.default_rng(0)
    pos = np.sort(rng.choice(1_500_000, size=60_000, replace=False))
    cnt = rng.integers(1, 50, size=60_000)
    t = time.perf_counter()
    peaks = _identify_peaks(pos, cnt, DEFAULT_MIN_PEAK_SEPARATION)
    _find_valleys_between_peaks(pos, cnt, peaks)
    assert time.perf_counter() - t < 10.0
    assert len(peaks) > 50_000
