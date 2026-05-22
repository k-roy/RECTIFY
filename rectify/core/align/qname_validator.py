"""Post-alignment FASTQ↔BAM QNAME validator + auto-sanitizer.

When an aligner mutates QNAMEs in a way that the K-way consensus merge
cannot recover from, this module catches it at alignment time rather than
at consensus time. By default it also **repairs** the mutation in place:
streams the BAM through pysam, normalizes every ``query_name`` via
``_normalize_bam_read_name`` (plus a whitespace strip), and writes back to
the same path. Only when the repair pass fails to make the BAM round-trip
to the source FASTQ does the validator raise.

This auto-sanitize behavior is on by default so that users without CLI
experience get a working pipeline. Set ``RECTIFY_NO_AUTO_QNAME_SANITIZE=1``
to disable repair (validator becomes detect-and-raise only). Set
``RECTIFY_SKIP_POST_ALIGN_VALIDATION=1`` to skip the validator entirely.

Three checks per BAM:
  1. SAM-spec compliance — no QNAME contains a space or tab character.
  2. Length cap — no QNAME exceeds 254 chars.
  3. Round-trip normalization — sampled BAM QNAMEs, after
     ``_normalize_bam_read_name``, must appear in the source FASTQ's
     normalized-QNAME set.

Sampling-based for runtime; full coverage of the K-way merge invariant
remains ``_check_read_name_compatibility`` at consensus time.
"""

from __future__ import annotations

import gzip
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


_SAMPLE_SIZE_DEFAULT = 200
_LOOKAHEAD_MULTIPLIER = 10
_MISMATCH_THRESHOLD = 0.05


# ─── helpers ──────────────────────────────────────────────────────────────────


def _strip_fastq_whitespace(raw: str) -> str:
    """Cut the FASTQ header body at the first space or tab."""
    cut = len(raw)
    for ws in (' ', '\t'):
        i = raw.find(ws)
        if 0 <= i < cut:
            cut = i
    return raw[:cut]


def _load_fastq_qname_set(fastq_path: str, max_records: Optional[int] = None) -> set:
    from rectify.core.consensus.consensus import _normalize_bam_read_name

    opener = gzip.open if str(fastq_path).endswith('.gz') else open
    out: set = set()
    with opener(fastq_path, 'rt') as fin:
        seen = 0
        while max_records is None or seen < max_records:
            header = fin.readline()
            if not header:
                break
            fin.readline()
            fin.readline()
            fin.readline()
            if not header.startswith('@'):
                continue
            raw = header[1:].rstrip('\n')
            out.add(_normalize_bam_read_name(_strip_fastq_whitespace(raw)))
            seen += 1
    return out


@dataclass
class _ScanResult:
    n_checked: int = 0
    n_missing: int = 0
    n_with_ws: int = 0
    n_overlength: int = 0
    n_needs_normalize: int = 0  # round-trip OK but QNAME would change under _normalize_bam_read_name
    example_missing: Optional[str] = None
    example_ws: Optional[str] = None
    example_overlength: Optional[str] = None
    example_needs_normalize: Optional[str] = None

    def is_clean(self) -> bool:
        """No detectable mutation in the sample.

        Includes the "cosmetic" case where round-trip succeeds but a QNAME
        would change under normalization (e.g., `uuid_pt:i:25` → `uuid`).
        The merge would still work without rewriting, but end users running
        `samtools view` on the per-aligner BAM would see the mutated form.
        """
        if self.n_with_ws or self.n_overlength or self.n_needs_normalize:
            return False
        if self.n_checked == 0:
            return True
        return (self.n_missing / self.n_checked) <= _MISMATCH_THRESHOLD

    def is_recoverable(self) -> bool:
        """All sampled mutations are repairable by streaming normalize.

        Round-trip failures (n_missing > threshold) indicate the normalizer
        regex can't recover the QNAME — repair will not fix the merge.
        Whitespace, overlength, and needs_normalize are all fixable.
        """
        if self.n_checked == 0:
            return True
        return (self.n_missing / self.n_checked) <= _MISMATCH_THRESHOLD

    def summary(self) -> str:
        parts = []
        if self.n_with_ws:
            parts.append(f"{self.n_with_ws} with whitespace")
        if self.n_overlength:
            parts.append(f"{self.n_overlength} >254 chars")
        if self.n_needs_normalize:
            parts.append(f"{self.n_needs_normalize} cosmetic (normalize-only)")
        if self.n_missing:
            parts.append(
                f"{self.n_missing}/{self.n_checked} not round-tripping"
            )
        return ', '.join(parts) or 'no issues'


def _scan_bam(bam_path: str, fastq_qnames: set, sample_size: int) -> _ScanResult:
    import pysam
    from rectify.core.consensus.consensus import _normalize_bam_read_name

    r = _ScanResult()
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            qn = read.query_name or ''
            normed = _normalize_bam_read_name(_strip_fastq_whitespace(qn))
            if ' ' in qn or '\t' in qn:
                r.n_with_ws += 1
                if r.example_ws is None:
                    r.example_ws = qn
            if len(qn) > 254:
                r.n_overlength += 1
                if r.example_overlength is None:
                    r.example_overlength = qn
            if normed not in fastq_qnames:
                r.n_missing += 1
                if r.example_missing is None:
                    r.example_missing = qn
            elif normed != qn:
                # Round-trip OK but QNAME contains a suffix the normalizer
                # would strip — rewrite so downstream `samtools view` sees
                # the bare form.
                r.n_needs_normalize += 1
                if r.example_needs_normalize is None:
                    r.example_needs_normalize = qn
            r.n_checked += 1
            if r.n_checked >= sample_size:
                break
    return r


def _rewrite_bam_normalize_qnames(bam_path: str) -> int:
    """Stream the BAM through pysam, normalizing every QNAME in place.

    Strips post-whitespace and applies ``_normalize_bam_read_name`` to each
    record. Writes to a sibling ``.qname_fixed.tmp.bam`` then atomically
    renames over the original. Returns the count of records whose QNAME
    actually changed.
    """
    import pysam
    from rectify.core.consensus.consensus import _normalize_bam_read_name

    src = Path(bam_path)
    tmp = src.with_suffix('.qname_fixed.tmp.bam')
    n_changed = 0
    with pysam.AlignmentFile(str(src), 'rb') as bam_in, \
         pysam.AlignmentFile(str(tmp), 'wb', header=bam_in.header) as bam_out:
        for read in bam_in:
            qn = read.query_name or ''
            new_qn = _normalize_bam_read_name(_strip_fastq_whitespace(qn))
            if new_qn != qn:
                read.query_name = new_qn
                n_changed += 1
            bam_out.write(read)
    tmp.replace(src)
    return n_changed


# ─── public entry point ──────────────────────────────────────────────────────


def validate_post_alignment_qnames(
    bam_path: str,
    fastq_path: str,
    aligner_name: str,
    sample_size: int = _SAMPLE_SIZE_DEFAULT,
    auto_sanitize: bool = True,
) -> None:
    """Validate (and optionally auto-repair) per-aligner BAM QNAMEs.

    By default, if the sample reveals QNAME mutations the validator first
    attempts a streaming in-place repair (normalize every record's QNAME),
    then re-validates. Only raises ``RuntimeError`` if the repair fails to
    make the BAM round-trip to the source FASTQ.

    Args:
        bam_path: Output BAM from the aligner wrapper.
        fastq_path: Source FASTQ that was passed to the wrapper.
        aligner_name: Identity string for log/error messages (e.g.
            ``"mapPacBio"``).
        sample_size: Number of primary BAM records to sample.
        auto_sanitize: When True (default), attempt in-place repair on
            detected mutations before raising. Disable via env var
            ``RECTIFY_NO_AUTO_QNAME_SANITIZE=1`` or by passing False.

    Env vars:
        ``RECTIFY_SKIP_POST_ALIGN_VALIDATION=1`` — skip the validator entirely.
        ``RECTIFY_NO_AUTO_QNAME_SANITIZE=1`` — disable auto-repair (raise on first detection).
    """
    if os.environ.get('RECTIFY_SKIP_POST_ALIGN_VALIDATION') == '1':
        logger.debug(f"[{aligner_name}] post-alignment QNAME validation skipped via env")
        return

    if os.environ.get('RECTIFY_NO_AUTO_QNAME_SANITIZE') == '1':
        auto_sanitize = False

    try:
        # The BAM is usually coordinate-sorted here, so its first sampled
        # alignments are not guaranteed to appear near the start of the FASTQ.
        fastq_qnames = _load_fastq_qname_set(fastq_path)
    except (OSError, IOError) as exc:
        logger.warning(
            f"[{aligner_name}] cannot read source FASTQ for validation "
            f"({fastq_path}): {exc}. Skipping post-alignment check."
        )
        return

    if not fastq_qnames:
        logger.warning(
            f"[{aligner_name}] source FASTQ produced no QNAMEs to validate against "
            f"({fastq_path}). Skipping post-alignment check."
        )
        return

    scan = _scan_bam(bam_path, fastq_qnames, sample_size)
    if scan.is_clean():
        logger.info(
            f"[{aligner_name}] QNAME post-alignment validation OK "
            f"(sampled {scan.n_checked} primary reads, 0 mutations detected)"
        )
        return

    if auto_sanitize and scan.is_recoverable():
        logger.warning(
            f"[{aligner_name}] QNAME mutations detected ({scan.summary()}) — "
            f"running auto-sanitizer..."
        )
        n_changed = _rewrite_bam_normalize_qnames(bam_path)
        logger.warning(
            f"[{aligner_name}] auto-sanitizer rewrote {n_changed} BAM record(s); "
            f"re-validating..."
        )
        scan = _scan_bam(bam_path, fastq_qnames, sample_size)
        if scan.is_clean():
            logger.info(
                f"[{aligner_name}] post-repair validation OK "
                f"(sampled {scan.n_checked} primary reads, mutations recovered)"
            )
            return

    _raise_for_scan(scan, aligner_name, repair_attempted=auto_sanitize)


def _raise_for_scan(
    scan: _ScanResult,
    aligner_name: str,
    repair_attempted: bool,
) -> None:
    suffix = (
        "\nAuto-sanitizer ran but could not repair the BAM. The aligner is "
        "mutating QNAMEs in a way the normalizer regex does not cover. "
        "Report this at https://github.com/anthropics/rectify/issues with the "
        "example QNAME above."
        if repair_attempted else
        "\nAuto-sanitizer was disabled "
        "(RECTIFY_NO_AUTO_QNAME_SANITIZE=1) — unset to attempt in-place repair."
    )

    if scan.n_with_ws:
        raise RuntimeError(
            f"[{aligner_name}] {scan.n_with_ws}/{scan.n_checked} BAM QNAMEs "
            f"contain whitespace (SAM spec violation). "
            f"Example: {scan.example_ws!r}.{suffix}"
        )
    if scan.n_overlength:
        raise RuntimeError(
            f"[{aligner_name}] {scan.n_overlength}/{scan.n_checked} BAM QNAMEs "
            f"exceed 254 chars (SAM spec limit). "
            f"Example: {scan.example_overlength!r}.{suffix}"
        )
    raise RuntimeError(
        f"[{aligner_name}] {scan.n_missing}/{scan.n_checked} "
        f"({100 * scan.n_missing / scan.n_checked:.1f}%) BAM QNAMEs do not "
        f"normalize back to any FASTQ QNAME. "
        f"Example: {scan.example_missing!r}.{suffix}"
    )
