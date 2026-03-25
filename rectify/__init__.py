"""
RECTIFY: Unified RNA 3' End Correction Framework

A modular framework for correcting 3' end mapping artifacts in poly(A)-tailed RNA sequencing data.

Modules:
- A-tract ambiguity detection (universal)
- AG mispriming screening (oligo-dT methods)
- Poly(A) tail trimming and indel correction (direct RNA-seq)
- NET-seq refinement (optional)
- Visualization (optional, requires matplotlib)

Features (v2.6.0):
- Region-based parallel BAM processing with coverage gap splitting
- SLURM-aware CPU detection to prevent oversubscription
- Streaming output mode for large BAM files
- Export to bedGraph/bigWig format
- Metagene signal aggregation and multi-track visualization

Author: Kevin R. Roy
License: MIT
"""

__version__ = "2.6.0"
__author__ = "Kevin R. Roy"
__email__ = "kevinroy@stanford.edu"

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
