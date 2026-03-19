"""
RECTIFY: Unified RNA 3' End Correction Framework

A modular framework for correcting 3' end mapping artifacts in poly(A)-tailed RNA sequencing data.

Modules:
- A-tract ambiguity detection (universal)
- AG mispriming screening (oligo-dT methods)
- Poly(A) tail trimming and indel correction (direct RNA-seq)
- NET-seq refinement (optional)

Features (v2.3.0):
- Region-based parallel BAM processing with coverage gap splitting
- SLURM-aware CPU detection to prevent oversubscription
- Streaming output mode for large BAM files
- Export to bedGraph/bigWig format

Author: Kevin R. Roy
License: MIT
"""

__version__ = "2.3.0"
__author__ = "Kevin R. Roy"
__email__ = "kevinroy@stanford.edu"

from . import core, utils, slurm
from .slurm import get_available_cpus, set_thread_limits, is_slurm_job

__all__ = [
    "core",
    "utils",
    "slurm",
    "get_available_cpus",
    "set_thread_limits",
    "is_slurm_job",
    "__version__",
]
