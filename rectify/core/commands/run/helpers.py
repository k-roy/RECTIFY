"""
Helpers shared across the ``rectify run-all`` runners.

- ``_resolve_reference_paths``: organism/bundled-data lookup that fills
  ``args.genome`` / ``args.annotation`` / ``args.go_annotations`` in place.
- ``_multialigned_bam_path``: canonical pre-correction merged multi-aligner BAM
  path for a sample (formerly ``<sample>.rectified.bam``; renamed 2026-06-25 so
  ``.rectified.bam`` can name the actually-rectified final output).
- ``_final_rectified_bam_path``: canonical final corrected BAM path
  (``<sample>.rectified.bam``; formerly ``corrected_consensus.bam``).
- ``_emit_legacy_consensus_symlink``: drop the back-compat
  ``corrected_consensus.bam`` symlink next to the final output.
- ``_collect_per_aligner_bams``: enumerate per-aligner BAMs that exist on disk.
- ``_bam_has_md_tags``: sample the first mapped reads for an MD tag.
- ``_validate_bam_integrity``: existence + ``.bai`` + ``samtools quickcheck``.
"""

import subprocess as _subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Optional


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

def _multialigned_bam_path(sample_id: str, sample_output_dir: Path) -> Path:
    """Canonical path for the pre-correction merged multi-aligner BAM for a sample.

    This is the alignment-stage artifact (per-placement aligner votes in the
    ``Xa`` tag) that the ``run`` reuse-gate keys on — NOT the final corrected
    output.  It was named ``<sample>.rectified.bam`` before 2026-06-25; the
    ``.rectified.bam`` token now names the actually-rectified final output
    (see ``_final_rectified_bam_path``).
    """
    return sample_output_dir / f"{sample_id}.multialigned.bam"


# Back-compat alias: legacy callers/imports referenced ``_rectified_bam_path``
# for the *pre-correction* artifact.  Keep the old name pointing at the
# multialigned path so any external importer keeps resolving the alignment
# artifact (its historical meaning), not the new final output.
_rectified_bam_path = _multialigned_bam_path


def _resolve_consensus_tag_bam(
    explicit: Optional[str],
    sample_id: str,
    search_dirs,
    fallback: Optional[Path] = None,
) -> Optional[Path]:
    """BAM to read ``Xa``/``Xc``/``Xn``/``Xt`` from for the merged TSV.

    Order: an explicit ``--consensus-bam``; else ``<sample>.multialigned.bam``
    (then the pre-2026-06-25 ``<sample>.consensus.bam``) in each of
    *search_dirs* — the align stage may write next to the per-aligner arms in
    scratch, in ``--bam-dir/<sample>``, or in the sample output dir; else
    *fallback*, the BAM the correct stage was actually handed.

    Returns None when nothing is found, and the merge then leaves the four
    ``consensus_*`` columns empty exactly as before.  Always prints which path
    was chosen: a silently-absent BAM and a genuinely untagged one produce the
    same empty columns downstream, and only the log distinguishes them.
    """
    if explicit:
        path = Path(explicit)
        if path.exists():
            print(f"    Consensus tags: {path} (--consensus-bam)")
            return path
        print(
            f"    WARNING: --consensus-bam {path} does not exist; "
            f"the consensus_* columns will be empty",
            file=sys.stderr,
        )
        return None

    seen = []
    for directory in search_dirs:
        if directory is None:
            continue
        directory = Path(directory)
        if directory in seen:
            continue
        seen.append(directory)
        for candidate, label in (
            (_multialigned_bam_path(sample_id, directory), 'auto'),
            (directory / f"{sample_id}.consensus.bam", 'auto, legacy name'),
        ):
            if candidate.exists():
                print(f"    Consensus tags: {candidate} ({label})")
                return candidate

    if fallback is not None and Path(fallback).exists():
        print(f"    Consensus tags: {fallback} (the corrected input BAM)")
        return Path(fallback)

    print(
        "    Consensus tags: no multialigned BAM found next to the arms "
        f"({', '.join(str(d) for d in seen) or 'no search dirs'}) — "
        "the consensus_* columns will be empty"
    )
    return None


def _final_rectified_bam_path(sample_id: str, sample_output_dir: Path) -> Path:
    """Canonical path for the FINAL corrected (rectified) BAM for a sample.

    This is the winner-take-all corrected output, formerly written as the bare
    ``corrected_consensus.bam`` (no sample prefix).  As of 2026-06-25 it carries
    the sample prefix like the per-aligner BAMs and is named
    ``<sample>.rectified.bam``.
    """
    return sample_output_dir / f"{sample_id}.rectified.bam"


def _emit_legacy_consensus_symlink(final_bam: Path) -> None:
    """Create back-compat ``corrected_consensus.bam`` symlink -> ``final_bam``.

    Downstream consumers (the workshop ``stage_igv.sh``, ``deseq_from_corrected``,
    ``splice_efficiency``, IGV, the shared ``processed/alignments/`` store) and
    existing on-disk runs reference the old bare name ``corrected_consensus.bam``.
    To keep old globs working for one release we drop a *relative* symlink
    (``corrected_consensus.bam`` and its ``.bai``) pointing at the new
    ``<sample>.rectified.bam`` in the same directory.

    The symlink target is a bare basename so the link stays valid when the
    per-sample directory is moved/synced (e.g. into ``processed/alignments/``).
    Idempotent: an existing link/file at the legacy path is replaced.

    NOTE: we deliberately do NOT emit a ``*.rectified.bam`` alias for the
    pre-correction merged file — that would put two files under the
    ``*.rectified.bam`` glob and reintroduce the swap-in-place hazard the rename
    was meant to remove.  Exactly one file matches ``*.rectified.bam`` per run.
    """
    final_bam = Path(final_bam)
    if not final_bam.exists():
        return
    legacy = final_bam.with_name("corrected_consensus.bam")
    for link_path, target_name in (
        (legacy, final_bam.name),
        (Path(str(legacy) + ".bai"), final_bam.name + ".bai"),
    ):
        # Only create the .bai link if the real index exists.
        if target_name.endswith(".bai") and not (final_bam.parent / target_name).exists():
            continue
        try:
            if link_path.is_symlink() or link_path.exists():
                link_path.unlink()
            link_path.symlink_to(target_name)  # relative target (same dir)
        except OSError:
            # Symlinks unsupported (rare FS); fall back to a copy so old globs
            # still resolve.  Non-fatal — the new name is the source of truth.
            try:
                import shutil as _shutil
                _shutil.copy2(final_bam.parent / target_name, link_path)
            except OSError:
                pass


#: Every aligner name ``align_command`` can emit a ``<sample>.<aligner>.bam`` for.
#:
#: 🔴 This list is what ``_collect_per_aligner_bams`` globs, so an arm missing here is INVISIBLE to
#: correction, to the Module-2H junction pool and to the merged ``corrected_reads.tsv`` — with no
#: error. The COMPASS arms (STAR x2, HISAT2 x2, magicblast, gsnap) were all missing, so
#: ``run-all --short-read`` (TruSeq SE and PE alike) collected the **bbmap arm only** while the
#: 5- or 7-way consensus sat unused in ``multialigned.bam`` (planning 833 C-2). Keep it a superset
#: of ``align_command``'s ``--aligners`` choices; ``COMPASS_{SE,PE}_ALIGNERS`` is the panel, this is
#: the emitter inventory.
_ALIGNER_NAMES = [
    # Long-read arms
    'minimap2', 'mapPacBio', 'gapmm2', 'uLTRA', 'deSALT',
    'winnowmap2', 'minisplice_mm2',
    # Short-read 3'-end arms (--dT-primed-cDNA panel)
    'bbmap', 'bwa',
    # COMPASS short-read panel (833 C-2)
    'STAR_default', 'STAR_noncanonical',
    'HISAT2_default', 'HISAT2_noncanonical',
    'magicblast', 'gsnap',
]


def _collect_per_aligner_bams(
    sample_id: str,
    sample_output_dir: Path,
    run_provenance: Optional[Dict[str, Any]] = None,
    trust_existing_bams: bool = False,
) -> Dict[str, Path]:
    """Return per-aligner BAM paths that exist on disk (keyed by aligner name).

    When ``run_provenance`` is provided and ``trust_existing_bams`` is False,
    each BAM is validated against its sidecar. BAMs with no sidecar or a
    mismatched sidecar are excluded (logged at WARNING level). Callers receive
    a smaller pool; downstream stages proceed with whatever matches.
    """
    bams: Dict[str, Path] = {}
    for aligner in _ALIGNER_NAMES:
        bam = sample_output_dir / f"{sample_id}.{aligner}.bam"
        if not bam.exists():
            continue
        if run_provenance is not None and not trust_existing_bams:
            try:
                from rectify.utils.bam_provenance import (
                    expected_provenance_for_aligner,
                    matches_strict,
                    read_sidecar,
                )
                import logging as _logging
                _log = _logging.getLogger(__name__)
                expected = expected_provenance_for_aligner(run_provenance, aligner)
                stored = read_sidecar(bam)
                ok, reason = matches_strict(stored, expected)
                if not ok:
                    _log.warning(
                        "Excluding %s per-aligner BAM — provenance mismatch (%s): %s. "
                        "Pass --trust-existing-bams to reuse anyway.",
                        aligner, reason, bam,
                    )
                    continue
            except Exception as _exc:
                import logging as _logging
                _logging.getLogger(__name__).warning(
                    "Could not check provenance for %s BAM (%s); including it anyway.",
                    aligner, _exc,
                )
        bams[aligner] = bam

    # Overhang-resolver substitution (Station A): when the resolver ran,
    # its output SUBSTITUTES the minimap2 arm for every downstream stage —
    # align_command.run_align performs the same swap in-process for
    # consensus, but this collector re-derives arms from DISK, and without
    # the swap per-aligner correction would silently run on the RAW
    # minimap2 BAM, undoing the resolution for the correct stage. Same
    # freshness gate as run_align's reuse check (newer than its input arm,
    # non-trivial size). The raw minimap2 BAM stays on disk for the XB
    # delta census; it is just not handed downstream.
    if 'minimap2' in bams:
        _resolver_bam = sample_output_dir / f"{sample_id}.overhang_resolver.bam"
        if (_resolver_bam.exists() and _resolver_bam.stat().st_size > 2000
                and _resolver_bam.stat().st_mtime
                >= bams['minimap2'].stat().st_mtime):
            import logging as _logging
            _logging.getLogger(__name__).info(
                "per-aligner collection: substituting the minimap2 arm with "
                "the overhang-resolver output (%s).", _resolver_bam)
            bams['minimap2'] = _resolver_bam
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
