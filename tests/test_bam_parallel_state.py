"""Byte-equivalence safety net for the ``bam_parallel`` worker-state refactor.

The Commit C-II.0 refactor changes ``_init_region_worker_state`` to hold
**only** the fields shared across aligners (``genome``, ``polya_model``,
``annotated_junctions``, ``pool_chrom_index``, ``gene_interval_trees``,
config flags, ``netseq_dir``) and moves the per-aligner fields
(``bam_path``, ``variant_aware_rescue``) into per-task tuples processed
by ``_process_region_worker_from_state``.

This test exercises ``process_bam_file_parallel`` on a tiny validation
BAM at ``n_threads=2`` (forcing the multi-process path), hashes a
sorted-by-read_id JSON of the result list, and asserts the hash matches
a recorded golden value. Pre-refactor and post-refactor hashes must be
byte-identical.

The golden hash is captured by running this test once with
``RECTIFY_RECORD_GOLDEN=1`` in the environment — it prints the observed
hash and fails the assertion so the developer knows to paste it back.

Why this matters: subtle bugs in worker state isolation (the same class
as the AGENT_FIXES "heap corruption at ≥4M-read scale" entry) often
don't surface until cluster scale. A small deterministic comparison
here catches state-leak regressions before they escape to H2.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest


# Recorded golden hash for the validation minimap2 BAM at n_threads=2.
# Capture/refresh via:
#     RECTIFY_RECORD_GOLDEN=1 pytest tests/test_bam_parallel_state.py -s
GOLDEN_HASH_VALIDATION_MINIMAP2_NT2 = (
    "a41ec734e189dcd0a2e2b96ceb092e490cbac5a6d148bba43261c69bc2ba03c6"
)
# Re-recorded 2026-06-29 (drs-validation-rebuild): the walkback homopolymer-undercall
# guard (walkback.py large-deletion pre-scan now preserves a deletion flanked 3' by a
# run of genuine read=ref matches — cat2_plus_1's 9D+39= over a 24-A genomic run) changed
# per-read minimap2 correction output. This is an intentional correction-semantics change;
# prior golden 4f326833… (2026-06-19).
# Re-recorded 2026-06-19 (drs-validation-rebuild). The prior golden
# (b267b1da…, recorded 2026-05-22 at 55089f7) went stale when later 3'SS-rescue
# commits legitimately changed process_bam_file_parallel output —
# bd20f9e (narrow 3'SS-rescue candidates), 961c844 (interval-tree pool index),
# cf5ebb9 (skip 3'SS rescue on rDNA/Pol III), 0c1773b (DRS 3'-end fix).
# Output is deterministic (stable observed hash across runs) and independent of
# the short-read consensus locus-index change (verified: baseline select.py
# produces the identical hash). 36 reads from validation_reads.minimap2.bam.


def _stable_hash_results(results):
    """Hash a list[dict] from process_bam_file_parallel deterministically.

    Sort by ``(read_id, chrom, original_3prime)`` so multi-region returns
    don't depend on region-completion order. Serialize each dict with
    sorted keys so dict iteration order can't perturb the hash.
    """
    sortable = sorted(
        results,
        key=lambda r: (
            r.get('read_id', ''),
            r.get('chrom', ''),
            r.get('original_3prime', 0),
        ),
    )
    # Lists inside dicts (correction_applied, qc_flags) get sorted for
    # idempotence too — order within a row shouldn't change output meaning.
    canonical = []
    for row in sortable:
        normalized = {}
        for k in sorted(row.keys()):
            v = row[k]
            if isinstance(v, list):
                v = sorted(v)
            normalized[k] = v
        canonical.append(normalized)
    payload = json.dumps(canonical, sort_keys=True, default=str).encode()
    return hashlib.sha256(payload).hexdigest()


@pytest.fixture(scope='module')
def validation_bam(tmp_path_factory):
    """Locate the validation minimap2 BAM relative to the rectify repo root."""
    here = Path(__file__).resolve()
    repo_root = here.parent.parent
    bam = repo_root / 'rectify' / 'data' / 'validation' / 'aligners' / 'validation_reads.minimap2.bam'
    if not bam.exists():
        pytest.skip(f"Validation BAM not found at {bam}")
    return str(bam)


@pytest.fixture(scope='module')
def validation_genome():
    here = Path(__file__).resolve()
    repo_root = here.parent.parent
    genome = (
        repo_root / 'rectify' / 'data' / 'genomes' /
        'saccharomyces_cerevisiae' / 'S288C_reference_sequence_R64-5-1_20240529.fsa.gz'
    )
    if not genome.exists():
        pytest.skip(f"Validation genome not found at {genome}")
    return str(genome)


def test_process_bam_file_parallel_deterministic(validation_bam, validation_genome):
    """process_bam_file_parallel must produce a deterministic result set.

    This test pins the byte-equivalence contract that Commit C-II.0
    must preserve. If you intentionally change the output schema or
    correction semantics, re-record the golden hash via
    ``RECTIFY_RECORD_GOLDEN=1``.
    """
    from rectify.core.bam.parallel import process_bam_file_parallel

    results = process_bam_file_parallel(
        bam_path=validation_bam,
        genome_path=validation_genome,
        n_threads=2,  # forces the mp.Pool path
        apply_atract=True,
        apply_ag_mispriming=False,
        ag_threshold=0.65,
        apply_polya_trim=False,
        apply_indel_correction=False,
        output_path=None,
        show_progress=False,
        return_stats=False,
        variant_aware=False,
        apply_3ss_rescue=True,
        polya_model_path=None,
        dt_primed_cDNA=False,
        min_mapq=0,
        min_aligned_length=0,
    )

    assert isinstance(results, list)
    assert len(results) > 0, "validation BAM should produce at least one result"

    observed = _stable_hash_results(results)

    if os.environ.get('RECTIFY_RECORD_GOLDEN') == '1':
        print(f"\n[RECORD_GOLDEN] hash = {observed}\n"
              f"[RECORD_GOLDEN] n_results = {len(results)}\n"
              f"  → paste into GOLDEN_HASH_VALIDATION_MINIMAP2_NT2")
        pytest.fail("RECORD_GOLDEN mode — paste the printed hash into the test constant.")

    assert observed == GOLDEN_HASH_VALIDATION_MINIMAP2_NT2, (
        f"process_bam_file_parallel output hash changed.\n"
        f"  observed: {observed}\n"
        f"  golden:   {GOLDEN_HASH_VALIDATION_MINIMAP2_NT2}\n"
        f"  n_results: {len(results)}\n"
        f"If this is an intentional semantic change, re-record with "
        f"RECTIFY_RECORD_GOLDEN=1 pytest {__file__} -s"
    )


def test_process_bam_file_parallel_single_threaded_matches_parallel(
    validation_bam, validation_genome
):
    """n_threads=1 (no Pool) must match n_threads=2 (Pool) byte-for-byte.

    This pins the contract that the per-task tuple shape doesn't perturb
    the single-threaded fallback. The fallback path at parallel.py:683
    calls ``_process_region_worker`` directly without going through
    ``_init_region_worker_state`` — both code paths must converge.
    """
    from rectify.core.bam.parallel import process_bam_file_parallel

    common_kwargs = dict(
        bam_path=validation_bam,
        genome_path=validation_genome,
        apply_atract=True,
        apply_ag_mispriming=False,
        ag_threshold=0.65,
        apply_polya_trim=False,
        apply_indel_correction=False,
        output_path=None,
        show_progress=False,
        return_stats=False,
        variant_aware=False,
        apply_3ss_rescue=True,
        polya_model_path=None,
        dt_primed_cDNA=False,
        min_mapq=0,
        min_aligned_length=0,
    )

    serial = process_bam_file_parallel(n_threads=1, **common_kwargs)
    parallel = process_bam_file_parallel(n_threads=2, **common_kwargs)

    assert _stable_hash_results(serial) == _stable_hash_results(parallel), (
        "n_threads=1 and n_threads=2 produced different hashes — worker "
        "state model has a serial/parallel divergence bug."
    )
