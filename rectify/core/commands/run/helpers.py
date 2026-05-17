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
    Resolve genome/annotation/GO annotation paths from explicit args or bundled data.

    Reads args.genome, args.annotation, args.organism.
    Updates args.genome, args.annotation, and args.go_annotations in place.

    Raises SystemExit if a genome cannot be resolved (nothing to align/correct against).
    """
    from rectify.data import ensure_reference_data, get_bundled_go_annotations_path, normalize_organism

    organism = getattr(args, 'organism', None)
    genome_arg = getattr(args, 'genome', None)
    annotation_arg = getattr(args, 'annotation', None)

    if organism is None and genome_arg is None:
        print(
            "ERROR: No reference genome provided.\n"
            "  Supply --genome /path/to/genome.fa,\n"
            "  or use --Scer (S. cerevisiae bundled data),\n"
            "  or use --organism <name> for another supported organism.",
            file=sys.stderr,
        )
        sys.exit(1)

    if organism is not None:
        genome_path, annotation_path, data_source = ensure_reference_data(
            organism=organism,
            custom_genome=genome_arg,
            custom_annotation=annotation_arg,
            verbose=True,
        )
    else:
        # Explicit paths only — validate they exist
        genome_path = Path(genome_arg) if genome_arg else None
        annotation_path = Path(annotation_arg) if annotation_arg else None
        data_source = 'custom' if genome_path else 'none'

    if genome_path is None:
        print(
            f"ERROR: No bundled genome available for organism '{organism}'. "
            "Use --genome to provide a custom reference.",
            file=sys.stderr,
        )
        sys.exit(1)

    args.genome = genome_path
    args.annotation = annotation_path

    # Resolve bundled GO annotations if not already set
    if not getattr(args, 'go_annotations', None) and organism:
        go_path = get_bundled_go_annotations_path(normalize_organism(organism))
        if go_path:
            args.go_annotations = go_path


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
