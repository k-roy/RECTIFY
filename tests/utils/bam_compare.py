"""QNAME-sort + record-by-record BAM equivalence helper.

Used by all parallel-vs-legacy equivalence tests across Commits A-F.

The primary entry point is ``assert_bams_equivalent(bam_a, bam_b)``.
It QNAME-sorts both BAMs and compares record-by-record on:
  - query_name
  - reference_name (chrom)
  - reference_start (POS)
  - flag
  - mapping_quality (MAPQ)
  - cigarstring
  - all BAM tags by name

Differences in @PG header lines are ignored by default (the parallel path
adds extra @PG lines from each region's merge).

Author: Kevin R. Roy
"""

from typing import List, Optional, Tuple
import pysam


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _sort_key(read: pysam.AlignedSegment) -> Tuple:
    """Stable sort key for QNAME ordering (reproducible across BAMs)."""
    return (
        read.query_name or '',
        read.reference_name or '',
        read.reference_start if not read.is_unmapped else -1,
        int(read.is_reverse),
        read.flag,
    )


def _load_reads(bam_path: str) -> List[pysam.AlignedSegment]:
    """Load all reads from a BAM into a list."""
    reads = []
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            reads.append(read)
    return reads


def _tags_dict(read: pysam.AlignedSegment) -> dict:
    """Return a dict of {tag_name: tag_value} for all tags on a read."""
    return dict(read.get_tags(with_value_type=False))


def _compare_one_read(
    read_a: pysam.AlignedSegment,
    read_b: pysam.AlignedSegment,
) -> List[str]:
    """Compare two reads field-by-field.

    Returns a list of difference descriptions (empty if reads are equivalent).
    """
    diffs: List[str] = []

    # Core alignment fields
    if read_a.query_name != read_b.query_name:
        diffs.append(f"query_name: {read_a.query_name!r} != {read_b.query_name!r}")
    if read_a.reference_name != read_b.reference_name:
        diffs.append(f"reference_name: {read_a.reference_name!r} != {read_b.reference_name!r}")
    if read_a.reference_start != read_b.reference_start:
        diffs.append(f"reference_start: {read_a.reference_start} != {read_b.reference_start}")
    if read_a.flag != read_b.flag:
        diffs.append(f"flag: {read_a.flag} != {read_b.flag}")
    if read_a.mapping_quality != read_b.mapping_quality:
        diffs.append(f"mapping_quality: {read_a.mapping_quality} != {read_b.mapping_quality}")
    if read_a.cigarstring != read_b.cigarstring:
        diffs.append(f"cigarstring: {read_a.cigarstring!r} != {read_b.cigarstring!r}")

    # Tags
    tags_a = _tags_dict(read_a)
    tags_b = _tags_dict(read_b)
    all_tag_names = sorted(set(tags_a) | set(tags_b))
    for tag in all_tag_names:
        if tag not in tags_a:
            diffs.append(f"tag {tag!r}: present in B only (value={tags_b[tag]!r})")
        elif tag not in tags_b:
            diffs.append(f"tag {tag!r}: present in A only (value={tags_a[tag]!r})")
        elif tags_a[tag] != tags_b[tag]:
            diffs.append(f"tag {tag!r}: {tags_a[tag]!r} != {tags_b[tag]!r}")

    return diffs


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def assert_bams_equivalent(
    bam_a: str,
    bam_b: str,
    *,
    ignore_pg: bool = True,
) -> None:
    """Assert two BAMs are equivalent modulo:
      - sort tie order (we QNAME-sort both before comparison)
      - the @PG chain (if ignore_pg=True; the parallel path adds extra @PG lines)

    Compares: CIGAR, FLAG, POS, MAPQ, all tags, query_name, reference_name.

    Raises AssertionError with a precise location on first divergence.

    Args:
        bam_a:     Path to the first BAM (e.g. legacy single-threaded output).
        bam_b:     Path to the second BAM (e.g. parallel output).
        ignore_pg: If True, skip @PG header comparison (default True).

    Raises:
        AssertionError: On any difference, with the offending QNAME and a
            summary of the first difference found.
    """
    reads_a = sorted(_load_reads(bam_a), key=_sort_key)
    reads_b = sorted(_load_reads(bam_b), key=_sort_key)

    # Read-count check first — asymmetric drops are the most common failure.
    if len(reads_a) != len(reads_b):
        raise AssertionError(
            f"BAM read count mismatch: {bam_a} has {len(reads_a)} reads, "
            f"{bam_b} has {len(reads_b)} reads."
        )

    # Record-by-record comparison (QNAME-sorted).
    for i, (ra, rb) in enumerate(zip(reads_a, reads_b)):
        diffs = _compare_one_read(ra, rb)
        if diffs:
            qname = ra.query_name or "(unknown)"
            diff_summary = "; ".join(diffs[:5])  # limit to first 5 diffs
            if len(diffs) > 5:
                diff_summary += f" ... (+{len(diffs) - 5} more)"
            raise AssertionError(
                f"BAM divergence at read index {i} (QNAME={qname!r}):\n  {diff_summary}"
            )
