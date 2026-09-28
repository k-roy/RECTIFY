"""The correction TSV's paired alignment state (ISSUE-079).

``rectify correct`` decides every correction row against ``bam_to_process`` — the input BAM
AFTER Module 2H has refined its junctions. The deferred consumers (lazy HP scoring in
``merge_corrected_tsvs``, the final ``write_corrected_consensus_bam``) replay those rows onto a
BAM record, so they must be handed the SAME alignment state. Before this module they were handed
the original arm BAM: every 2H placement was silently undone in the final rectified BAM while the
merged TSV kept describing it (four of eight witness molecules, 16 query bases each, 16/16 exact
matches -> 6/16).

The TSV cannot carry the state itself — a junction string does not say where the query bases
went — so ``correct --retain-writer-input`` keeps the exact refined BAM beside its TSV and
publishes a receipt LAST. Consumers resolve the BAM through :func:`resolve_writer_input`, which
refuses anything the receipt does not vouch for; a legacy TSV with no receipt is rebuilt, never
replayed onto raw geometry.

Receipt states::

    ready                 the named BAM (or the arm's own BAM, kind='source') pairs with the TSV
    evicted_after_final   the retained BAM was deleted once the final outputs were complete
"""

import hashlib
import json
import logging
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional

logger = logging.getLogger(__name__)

SCHEMA = 1
ROLE = 'post_preprocessing_pre_correction'
STATE_READY = 'ready'
STATE_EVICTED = 'evicted_after_final'
#: The writer input IS the arm BAM the caller already holds (2H did not run, or failed open).
KIND_SOURCE = 'source'
#: The writer input is the retained post-2H BAM named by the receipt.
KIND_RETAINED = 'retained'

_CANONICAL_SUFFIXES = ('.manifest.tsv', '.region_000.tsv', '.tsv')


class WriterInputError(RuntimeError):
    """The receipt is absent, stale, evicted or disagrees with the bytes on disk."""


def _stem(corrected_tsv: Path) -> str:
    name = Path(corrected_tsv).name
    for suffix in _CANONICAL_SUFFIXES:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return Path(corrected_tsv).stem


def receipt_path(corrected_tsv: Path) -> Path:
    """Receipt beside the TSV; the flat TSV, its manifest and its region TSV share one."""
    corrected_tsv = Path(corrected_tsv)
    return corrected_tsv.parent / f'{_stem(corrected_tsv)}.writer_input.json'


def retained_bam_path(corrected_tsv: Path) -> Path:
    corrected_tsv = Path(corrected_tsv)
    return corrected_tsv.parent / f'{_stem(corrected_tsv)}.writer_input.bam'


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _identity(path) -> Dict[str, object]:
    path = Path(str(path))
    return {'name': path.name, 'size': path.stat().st_size if path.exists() else None}


def _pool_identity(pool_bams: Optional[Iterable]) -> List[Dict[str, object]]:
    return sorted((_identity(p) for p in pool_bams or []), key=lambda d: str(d['name']))


def _table_files(corrected_tsv: Path) -> List[Path]:
    """The canonical entry point plus every region TSV a manifest names."""
    corrected_tsv = Path(corrected_tsv)
    files = [corrected_tsv]
    if corrected_tsv.name.endswith('.manifest.tsv'):
        from .tsv_partition import load_manifest
        files.extend(Path(str(e['tsv_path'])) for e in load_manifest(corrected_tsv))
    return files


def invalidate(corrected_tsv: Path) -> None:
    """Withdraw the receipt BEFORE a correction run touches its TSV.

    A run that later fails leaves no receipt, so the catch-and-continue callers cannot mistake
    the previous run's TSV for this run's success.
    """
    receipt = receipt_path(corrected_tsv)
    retained = retained_bam_path(corrected_tsv)
    for path in (receipt, retained, Path(str(retained) + '.bai')):
        if path.exists():
            path.unlink()


def publish(
    corrected_tsv: Path,
    writer_input_bam: Path,
    *,
    owned: bool,
    origin_bam: Path,
    pool_bams: Optional[Iterable] = None,
    refinement: Optional[Dict[str, object]] = None,
    producer: Optional[str] = None,
) -> Path:
    """Bind ``writer_input_bam`` to ``corrected_tsv`` and publish the receipt last.

    ``owned`` says the BAM is this run's own 2H intermediate: it is MOVED into the retained name
    (same directory, no second copy). A BAM the run does not own is never moved, linked or
    copied — the receipt records ``kind='source'`` and consumers use the arm BAM they hold, whose
    records it shares (a coordinate-sorted or spike-in-filtered view selects the same records by
    read id).
    """
    corrected_tsv = Path(corrected_tsv)
    writer_input_bam = Path(writer_input_bam)
    receipt = receipt_path(corrected_tsv)

    body: Dict[str, object] = {
        'schema': SCHEMA,
        'role': ROLE,
        'state': STATE_READY,
        'producer': producer,
        'origin': _identity(origin_bam),
        'pool': _pool_identity(pool_bams),
        'refinement': refinement or {},
    }
    if owned:
        retained = retained_bam_path(corrected_tsv)
        os.replace(str(writer_input_bam), str(retained))
        index = Path(str(writer_input_bam) + '.bai')
        if index.exists():
            os.replace(str(index), str(retained) + '.bai')
        body.update(kind=KIND_RETAINED, bam=retained.name,
                    bam_size=retained.stat().st_size, bam_sha256=_sha256(retained))
    else:
        body.update(kind=KIND_SOURCE, bam=None, bam_size=None, bam_sha256=None)

    body['tables'] = [
        {'name': os.path.relpath(str(p), str(receipt.parent)),
         'size': p.stat().st_size, 'sha256': _sha256(p)}
        for p in _table_files(corrected_tsv)
    ]

    tmp = receipt.with_name(receipt.name + '.tmp')
    with open(tmp, 'w') as fh:
        json.dump(body, fh, indent=1, sort_keys=True)
        fh.write('\n')
    os.replace(str(tmp), str(receipt))
    return receipt


def load(
    corrected_tsv: Path,
    *,
    origin_bam: Optional[Path] = None,
    pool_bams: Optional[Iterable] = None,
    allow_evicted: bool = False,
    verify_bytes: bool = True,
) -> Dict[str, object]:
    """Return the validated receipt, or raise :class:`WriterInputError` naming the reason.

    ``verify_bytes=False`` checks presence and size only — for eviction, which is about to delete
    the BAM and has no use for a second full read of it.
    """
    corrected_tsv = Path(corrected_tsv)
    receipt = receipt_path(corrected_tsv)
    if not receipt.exists():
        raise WriterInputError(f'no writer-input receipt beside {corrected_tsv}')
    try:
        with open(receipt) as fh:
            body = json.load(fh)
    except (OSError, ValueError) as exc:
        raise WriterInputError(f'unreadable writer-input receipt {receipt}: {exc}') from exc
    if body.get('schema') != SCHEMA or body.get('role') != ROLE:
        raise WriterInputError(f'{receipt}: unknown schema/role')

    state = body.get('state')
    if state == STATE_EVICTED and not allow_evicted:
        raise WriterInputError(f'{receipt}: paired input evicted after final outputs')
    if state not in (STATE_READY, STATE_EVICTED):
        raise WriterInputError(f'{receipt}: state {state!r}')

    for entry in body.get('tables') or []:
        path = receipt.parent / str(entry['name'])
        if not path.exists():
            raise WriterInputError(f'{receipt}: correction table {path.name} is missing')
        if path.stat().st_size != entry['size'] or (
                verify_bytes and _sha256(path) != entry['sha256']):
            raise WriterInputError(f'{receipt}: correction table {path.name} changed')
    if not body.get('tables'):
        raise WriterInputError(f'{receipt}: binds no correction table')

    if origin_bam is not None and _identity(origin_bam)['size'] != body['origin']['size']:
        raise WriterInputError(f'{receipt}: arm BAM {Path(str(origin_bam)).name} changed')
    if pool_bams is not None and _pool_identity(pool_bams) != body.get('pool'):
        raise WriterInputError(f'{receipt}: candidate-pool BAMs changed')

    if state == STATE_READY and body.get('kind') == KIND_RETAINED:
        retained = receipt.parent / str(body['bam'])
        if not retained.exists():
            raise WriterInputError(f'{receipt}: retained BAM {retained.name} is missing')
        if retained.stat().st_size != body['bam_size'] or (
                verify_bytes and _sha256(retained) != body['bam_sha256']):
            raise WriterInputError(f'{receipt}: retained BAM {retained.name} changed')
    elif state == STATE_READY and body.get('kind') != KIND_SOURCE:
        raise WriterInputError(f'{receipt}: kind {body.get("kind")!r}')
    return body


def resolve_writer_input(
    corrected_tsv: Path,
    arm_bam: Path,
    *,
    pool_bams: Optional[Iterable] = None,
) -> Path:
    """The BAM the deferred consumers must replay ``corrected_tsv`` onto.

    ``arm_bam`` is the arm's own alignment BAM; it is returned only when the receipt says the
    correction was computed on its records (kind='source'). Raises :class:`WriterInputError`
    otherwise — never falls back to raw geometry.
    """
    body = load(corrected_tsv, origin_bam=arm_bam, pool_bams=pool_bams)
    if body['kind'] == KIND_RETAINED:
        return receipt_path(corrected_tsv).parent / str(body['bam'])
    if not Path(str(arm_bam)).exists():
        raise WriterInputError(f'arm BAM {arm_bam} is missing')
    return Path(str(arm_bam))


def resolve_writer_inputs(
    per_aligner_tsvs: Dict[str, Path],
    per_aligner_bams: Dict[str, Path],
) -> Dict[str, str]:
    """:func:`resolve_writer_input` for every corrected arm, against the complete pool."""
    pool = list(per_aligner_bams.values())
    return {
        aligner: str(resolve_writer_input(tsv, per_aligner_bams[aligner], pool_bams=pool))
        for aligner, tsv in per_aligner_tsvs.items()
    }


def evict(corrected_tsv: Path) -> bool:
    """Delete a retained BAM once the final outputs exist; mark the receipt first.

    State changes BEFORE the delete: a crash in between leaves extra bytes, never a ``ready``
    receipt pointing at nothing. Returns True when a BAM was removed.
    """
    corrected_tsv = Path(corrected_tsv)
    receipt = receipt_path(corrected_tsv)
    try:
        body = load(corrected_tsv, verify_bytes=False)
    except WriterInputError:
        return False
    if body['kind'] != KIND_RETAINED:
        return False
    body['state'] = STATE_EVICTED
    tmp = receipt.with_name(receipt.name + '.tmp')
    with open(tmp, 'w') as fh:
        json.dump(body, fh, indent=1, sort_keys=True)
        fh.write('\n')
    os.replace(str(tmp), str(receipt))
    retained = receipt.parent / str(body['bam'])
    for path in (retained, Path(str(retained) + '.bai')):
        if path.exists():
            path.unlink()
    return True


def evict_all(per_aligner_tsvs: Dict[str, Path]) -> int:
    """:func:`evict` every arm; returns how many retained BAMs were released."""
    return sum(1 for tsv in per_aligner_tsvs.values() if evict(tsv))
