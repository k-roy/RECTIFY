"""FASTQ chunking support utilities."""

from .sidecar import (
    ReadNumSidecar,
    ReadNumSidecarWriter,
    SidecarRow,
    VerificationResult,
)

__all__ = [
    'ReadNumSidecar',
    'ReadNumSidecarWriter',
    'SidecarRow',
    'VerificationResult',
]
