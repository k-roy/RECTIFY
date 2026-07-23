"""CMA (compressed-multialign) — pre-correct multi-aligner footprint compression.

RECTIFY runs a panel of aligners per read and retains the pre-correct per-aligner
BAMs (the only lossless store that supports add-an-aligner and resume-correct).
The per-aligner SEQ+QUAL is byte-identical across aligners and, measured on
Hoffman2 DRS data, is ~81% of per-aligner BAM bytes
(planning/256_rectify_M0_measurement.md). CMA stores that read payload **once**
plus one lightweight record per distinct cross-aligner placement.

Public surface:
    build_cma(read_stream, template_header, out_path, panel)   # cma_writer
    expand(cma_path[, genome])                                 # cma_reader
    validate_cma(cma_path)                                     # validate

Design: planning/254_rectify_multialign_compression_plan.md (§2 record model,
§2.4 SEQ/QUAL reconstruction frame, §3.1 expand()).
"""

from .cma_schema import (
    SCHEMA_VERSION,
    TAG_ALIGNERS,
    TAG_PAYLOAD,
    TAG_K,
    TAG_QSLICE,
    TAG_MAPQ,
    collapse_key,
    revcomp,
)
from .cma_writer import build_cma, load_aligner_records
from .cma_reader import expand
from .validate import validate_cma

__all__ = [
    "SCHEMA_VERSION",
    "TAG_ALIGNERS",
    "TAG_PAYLOAD",
    "TAG_K",
    "TAG_QSLICE",
    "TAG_MAPQ",
    "collapse_key",
    "revcomp",
    "build_cma",
    "load_aligner_records",
    "expand",
    "validate_cma",
]
