"""UMI / position clustering for the cdna pipeline.

The graph-clustering algorithms themselves are NOT cdna-specific (they operate on
plain ``List[str]`` + Levenshtein), so they now live in
:mod:`rectify.core.umi.clustering` and are shared with the short-read UMI path.
They are re-exported here unchanged, so this module's public API is identical to
before the move:

  * :func:`umi_components` -- connected-components (legacy, pre-v1.13).
  * :func:`umi_components_directional` -- umi_tools-style directional (v1.13+).
  * :func:`position_components_directional` -- v1.17 isoform grouping on positions.
  * :func:`canonical_umi` -- pick a cluster's canonical UMI.

What stays here is what genuinely depends on the cdna :class:`ReadInfo` dataclass:
:func:`position_components` and :func:`umi_canon_placeholder`.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Tuple

# Re-exported for backward compatibility -- these used to be defined in this file.
from rectify.core.umi.clustering import (  # noqa: F401
    canonical_umi,
    position_components_directional,
    umi_components,
    umi_components_directional,
)

from .read_info import ReadInfo

__all__ = [
    "umi_components",
    "umi_components_directional",
    "position_components_directional",
    "position_components",
    "canonical_umi",
    "umi_canon_placeholder",
]


def position_components(reads: List[ReadInfo]) -> List[List[int]]:
    """v1.15: exact (aln_start, aln_end) grouping for Type-2 (SSP-less) reads.

    No UMI is available for Type-2 reads, so we collapse exact-position duplicates
    only. This is conservative — molecules sharing the SAME aln_start AND aln_end
    are treated as PCR/sequencing duplicates of one decay/truncation event.
    Wobble of even ±1 bp produces a new cluster; that's fine because real distinct
    decay-intermediate species typically have well-defined termini.
    """
    pos_to_indices: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i, r in enumerate(reads):
        pos_to_indices[(r.aln_start, r.aln_end)].append(i)
    return list(pos_to_indices.values())


def umi_canon_placeholder(cluster: List[ReadInfo]) -> str:
    """Return the canonical UMI of a Type-1 cluster (most-frequent UMI string).
    For Type-2 clusters this returns empty string."""
    counts: defaultdict[str, int] = defaultdict(int)
    for r in cluster:
        counts[r.umi] += 1
    if not counts: return ""
    return min(counts, key=lambda u: (-counts[u], u))
