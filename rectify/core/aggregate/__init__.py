"""
Aggregate modules for RECTIFY.

These modules aggregate per-read data into datasets for downstream analysis:
- three_prime: 3' end (CPA) clustering with 5' gene attribution
- five_prime: 5' end (TSS) clustering with 3' gene attribution
- junctions: Junction aggregation with multi-aligner consensus

Author: Kevin R. Roy
"""

from .three_prime import (
    cluster_3prime_ends,
    aggregate_3prime_clusters,
)

from .five_prime import (
    cluster_5prime_ends,
    aggregate_5prime_clusters,
)

from .junctions import (
    aggregate_junctions,
    merge_with_partial_evidence,
)

__all__ = [
    'cluster_3prime_ends',
    'aggregate_3prime_clusters',
    'cluster_5prime_ends',
    'aggregate_5prime_clusters',
    'aggregate_junctions',
    'merge_with_partial_evidence',
]
