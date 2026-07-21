"""Protocol-agnostic UMI support for RECTIFY.

Standard fixed-length random N-mer UMIs (Illumina short-read chemistries such as
Lexogen CORALL V2), as distinct from the ONT PCB114 structured-motif UMI handled
in ``core/cdna/``.

Closes the gap documented in ``docs/comparison.md``: *"RECTIFY does not consume
Lexogen UMIs (use umi_tools dedup before rectify correct)."*

Layout:
  ``clustering`` -- shared directional/components graph clustering (also re-exported
                    by ``core/cdna/umi.py``, which is where it used to live).
  ``extract``    -- pull the UMI off a FASTQ record (read-id / r1-start / r2-start).
  ``grouping``   -- paired-end fragment-span keys (soft-clip corrected, pair-status aware).
"""
from __future__ import annotations

from .clustering import (
    canonical_umi,
    position_components_directional,
    umi_components,
    umi_components_directional,
)
from .extract import (
    CORALL_UMI_LENGTH,
    DEFAULT_UMI_SEPARATOR,
    UMI_LOCATIONS,
    ExtractedUmi,
    extract_umi,
    extract_umi_from_read_id,
    extract_umi_from_sequence,
    is_valid_umi,
    read_name_of,
    strip_umi_from_read_id,
)
from .grouping import (
    FragmentKey,
    PairStatus,
    fragment_key,
    fragment_strand,
    five_prime_unclipped,
    is_spliced,
    leading_soft_clip,
    trailing_soft_clip,
)
from .cli import (
    UMI_CLUSTERING_CHOICES,
    UMI_COLLAPSE_CHOICES,
    add_umi_args,
    validate_umi_args,
)
from .dedup import (
    FragmentRecord,
    UmiDedupStats,
    dedup_bam,
    extract_fragments,
    select_molecules,
)

__all__ = [
    # clustering
    "canonical_umi",
    "position_components_directional",
    "umi_components",
    "umi_components_directional",
    # extract
    "CORALL_UMI_LENGTH",
    "DEFAULT_UMI_SEPARATOR",
    "UMI_LOCATIONS",
    "ExtractedUmi",
    "extract_umi",
    "extract_umi_from_read_id",
    "extract_umi_from_sequence",
    "is_valid_umi",
    "read_name_of",
    "strip_umi_from_read_id",
    # grouping
    "FragmentKey",
    "PairStatus",
    "fragment_key",
    "fragment_strand",
    "five_prime_unclipped",
    "is_spliced",
    "leading_soft_clip",
    "trailing_soft_clip",
    # cli
    "UMI_CLUSTERING_CHOICES",
    "UMI_COLLAPSE_CHOICES",
    "add_umi_args",
    "validate_umi_args",
    # dedup
    "FragmentRecord",
    "UmiDedupStats",
    "dedup_bam",
    "extract_fragments",
    "select_molecules",
]
