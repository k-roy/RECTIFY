"""
Helpers shared across the ``rectify run-all`` runners.

- ``_resolve_reference_paths``: organism/bundled-data lookup that fills
  ``args.genome`` / ``args.annotation`` / ``args.go_annotations`` in place.
- ``_rectified_bam_path``: canonical rectified-BAM path for a sample.
- ``_collect_per_aligner_bams``: enumerate per-aligner BAMs that exist on disk.
- ``_bam_has_md_tags``: sample the first mapped reads for an MD tag.
- ``_validate_bam_integrity``: existence + ``.bai`` + ``samtools quickcheck``.
"""

import sys
from pathlib import Path
from typing import Dict, Optional


# ---------------------------------------------------------------------------
# Reference path resolution
# ---------------------------------------------------------------------------

def _resolve_reference_paths(args) -> None:
    """
    Legacy in-place reference resolver for ``rectify run-all`` / ``batch``.

    Thin wrapper around :func:`rectify.data.resolve_reference_paths` preserved
    for direct internal callers (``run_command.run``, ``batch_command`` slurm
    generators, ``tests/test_run_command_wiring.py``) that need the strict
    "no genome → exit" behavior. New subcommands should rely on the global
    hook in :func:`rectify.cli.main`, which calls
    ``resolve_reference_paths(..., require_genome=False)``.
    """
    from rectify.data import resolve_reference_paths

    resolve_reference_paths(args, require_genome=True, verbose=True)


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def _rectified_bam_path(sample_id: str, sample_output_dir: Path) -> Path:
    """Canonical path for the rectified BAM for a sample."""
    return sample_output_dir / f"{sample_id}.rectified.bam"


_ALIGNER_NAMES = ['minimap2', 'mapPacBio', 'gapmm2', 'uLTRA', 'deSALT']


def _collect_per_aligner_bams(
    sample_id: str,
    sample_output_dir: Path,
) -> Dict[str, Path]:
    """Return per-aligner BAM paths that exist on disk (keyed by aligner name)."""
    bams: Dict[str, Path] = {}
    for aligner in _ALIGNER_NAMES:
        bam = sample_output_dir / f"{sample_id}.{aligner}.bam"
        if bam.exists():
            bams[aligner] = bam
    return bams


# ---------------------------------------------------------------------------
# BAM integrity checks
# ---------------------------------------------------------------------------

def _bam_has_md_tags(bam_path: Path) -> bool:
    """Check if a BAM file's reads have MD tags (sample first 10 mapped reads)."""
    import logging
    logger = logging.getLogger(__name__)
    try:
        import pysam
        with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
            checked = 0
            for read in bam:
                if read.is_unmapped or read.is_secondary:
                    continue
                if read.has_tag('MD'):
                    return True
                checked += 1
                if checked >= 10:
                    break
    except FileNotFoundError:
        raise
    except Exception as e:
        logger.warning(f"Could not check MD tags: {e}")
    return False


def _validate_bam_integrity(bam_path: Path) -> bool:
    """
    Validate BAM file integrity before reuse.

    Returns True only if:
    - BAM file exists
    - BAM index (.bai) exists
    - samtools quickcheck passes

    Prevents reuse of corrupt or truncated BAMs from prior crashed runs.
    """
    if not bam_path.exists():
        return False

    # Check for .bai index (required for pysam.fetch)
    bai_path = Path(str(bam_path) + '.bai')
    alt_bai = bam_path.with_suffix('.bam.bai')
    if not bai_path.exists() and not alt_bai.exists():
        return False

    # Run samtools quickcheck to detect truncation/corruption
    result = _subprocess.run(
        ['samtools', 'quickcheck', str(bam_path)],
        capture_output=True,
        timeout=30,
    )
    return result.returncode == 0
