"""
RECTIFY: Unified RNA 3' End Correction Framework

A modular framework for correcting 3' end mapping artifacts in poly(A)-tailed RNA sequencing data.

Modules:
- A-tract ambiguity detection (universal)
- AG mispriming screening (oligo-dT methods)
- Poly(A) tail trimming and indel correction (direct RNA-seq)
- NET-seq refinement (optional)
- Visualization (optional, requires matplotlib)

Features (v0.9.0 — first public release):
- Multi-aligner consensus (minimap2, mapPacBio, gapmm2, uLTRA, deSALT) with per-read selection
- 3' end indel correction and soft-clip rescue at homopolymer boundaries
- 5' splice-junction rescue (Cat3) and intronic-tail rerouting
- False junction handling for poly(A) tails landing on genomic A-tracts
- Region-based parallel BAM processing with coverage-gap splitting
- Streaming output mode for large BAM files
- SLURM-aware CPU detection to prevent oversubscription
- NET-seq NNLS deconvolution for A-tract refinement
- Adaptive valley-based CPA clustering and DESeq2 at gene + cluster resolution
- Export to bedGraph/bigWig; metagene and genome-browser visualization

Author: Kevin R. Roy
License: MIT
"""

__version__ = "0.9.0"
__author__ = "Kevin R. Roy"
__email__ = "kevinrjroy@gmail.com"

from . import core, utils, slurm
from .slurm import get_available_cpus, set_thread_limits, is_slurm_job

# Conditionally import visualize module if matplotlib is available
try:
    from . import visualize
    _VISUALIZE_AVAILABLE = True
except ImportError:
    _VISUALIZE_AVAILABLE = False

__all__ = [
    "core",
    "utils",
    "slurm",
    "get_available_cpus",
    "set_thread_limits",
    "is_slurm_job",
    "__version__",
]

if _VISUALIZE_AVAILABLE:
    __all__.append("visualize")
