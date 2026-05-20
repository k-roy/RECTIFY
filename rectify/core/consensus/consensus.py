"""
Multi-Aligner Consensus Module for RECTIFY.

Streams through name-sorted BAMs from multiple aligners and selects the
best alignment per read. The per-read scoring + selection layer lives in
peer modules:

- ``extract.py``  — ``AlignmentInfo`` / ``ConsensusResult`` dataclasses,
                    ``extract_alignment_info``, ``check_canonical_splice_sites``,
                    plus the small CIGAR-extract helpers.
- ``scoring.py``  — ``score_alignment`` + the per-alignment scoring
                    primitives, including the 5' soft-clip rescue.
- ``select.py``   — ``select_best_alignment`` (per-read winner picking).

This module holds the streaming/orchestration layer:

- ``run_consensus_selection`` — top-level streaming entry point
- ``_process_and_write_batch`` — per-batch dispatch + output writing
- name-grouped K-way merge across name-sorted BAMs
  (``_iter_name_grouped_bams`` and friends)
- ``_restore_sequence_from_aligner_reads`` (kept here so the patched
  ``"rectify.core.consensus.consensus"`` logger in
  ``tests/test_gapmm2_seq_restore.py`` continues to attach to its logs)
- ``load_annotated_junctions`` (GFF/GTF intron loader)
- ``merge_slurm_array_bams`` (post-array-job consolidation)

Re-exports the names that external callers and tests import from
``rectify.core.consensus.consensus``: ``AlignmentInfo``, ``ConsensusResult``,
``extract_alignment_info``, ``extract_junctions_from_cigar``,
``get_softclip_lengths``, ``check_canonical_splice_sites``,
``score_alignment``, ``_rescue_5prime_softclip``, ``select_best_alignment``.

Scoring priorities:
1. Prefer alignments that splice through junctions vs soft-clipping (5' rescue):
   - If a 5' soft-clip can be explained by a missed intron (the clipped bases
     match the upstream exon end, checked via edit distance), the soft-clip is
     "rescued" and carries no penalty. This avoids penalizing aligners that
     correctly identify a junction but soft-clip a few upstream exon bases,
     relative to aligners that simply start mapping AFTER the junction.
   - Per-read junction pool = annotated junctions UNION all aligners' observed
     junctions for this read.
2. Prefer alignments whose 3' end lands outside a downstream A-tract (3' end quality)
3. Prefer junctions supported by multiple aligners
4. Tiebreaker: prefer aligner whose corrected 3' position agrees with majority
5. Tiebreaker: canonical splice site motifs (GT/AG) and annotated junctions

Note: A-tract 3' correction is applied to each aligner pre-scoring using genome
sequence only. Full indel correction (MD-tag dependent) is applied post-consensus
as a refinement step.

Author: Kevin R. Roy
"""

import logging
import os
import hashlib
import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Set

import pysam

# Re-exports for backwards compatibility with external callers / tests that
# import from rectify.core.consensus.consensus.  Listed in __all__ below.
from .extract import (  # noqa: F401  (re-exported)
    AlignmentInfo,
    ConsensusResult,
    extract_alignment_info,
    extract_junctions_from_cigar,
    get_softclip_lengths,
    check_canonical_splice_sites,
)
from .scoring import (  # noqa: F401  (re-exported)
    score_alignment,
    _rescue_5prime_softclip,
)
from .select import select_best_alignment  # noqa: F401  (re-exported)

logger = logging.getLogger(__name__)

__all__ = [
    # Re-exported from extract/scoring/select for backwards compatibility
    'AlignmentInfo',
    'ConsensusResult',
    'extract_alignment_info',
    'extract_junctions_from_cigar',
    'get_softclip_lengths',
    'check_canonical_splice_sites',
    'score_alignment',
    '_rescue_5prime_softclip',
    'select_best_alignment',
    # Native to this module
    '_process_and_write_batch',
    '_restore_sequence_from_aligner_reads',
    'run_consensus_selection',
    'merge_slurm_array_bams',
    'load_annotated_junctions',
]


def _ensure_name_sorted(bam_path: str) -> str:
    """
    Ensure a BAM file is name-sorted. If not, create a name-sorted copy.

    Returns path to name-sorted BAM (may be same as input if already sorted).
    """
    bam = pysam.AlignmentFile(bam_path, 'rb')
    header = bam.header.to_dict()
    bam.close()

    sort_order = header.get('HD', {}).get('SO', 'unknown')
    if sort_order == 'queryname':
        logger.debug(f"BAM already name-sorted: {bam_path}")
        return bam_path

    sorted_path = bam_path.replace('.bam', '.namesorted.bam')
    if os.path.exists(sorted_path):
        if os.path.getmtime(sorted_path) > os.path.getmtime(bam_path):
            logger.info(f"Using existing name-sorted BAM: {sorted_path}")
            return sorted_path

    logger.info(f"Name-sorting BAM: {bam_path} -> {sorted_path}")
    # Cap samtools sort memory: without -m, samtools uses 768MB × all threads.
    # With 5 BAMs sorted sequentially, Python's allocator retains each peak
    # in RSS, compounding to ~60GB on a 16-core node. 1G per sort is ample
    # for typical per-sample BAM sizes.
    pysam.sort('-n', '-m', '1G', '-o', sorted_path, bam_path)
    return sorted_path


def _read_id_hash(read_id: str, n_buckets: int) -> int:
    """Deterministic hash of read_id for SLURM array splitting."""
    h = hashlib.md5(read_id.encode()).hexdigest()
    return int(h, 16) % n_buckets


def _filtered_read_iterator(bam: pysam.AlignmentFile):
    """Yield only primary, mapped reads from a BAM file."""
    for read in bam:
        if not (read.is_unmapped or read.is_secondary or read.is_supplementary):
            yield read


def _natural_sort_key(s: str) -> list:
    """Key for natural (version) sort matching samtools queryname:natural order.

    Samtools natural sort compares runs of digits as integers rather than
    lexicographically.  Example: ``98297e97`` sorts before ``0633141e``
    because 98297 < 633141, even though ``'9' > '0'`` lexicographically.
    The K-way merge must use the same ordering as the BAM iterators so that
    reads present in only a subset of aligners do not desynchronise the merge.
    """
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', s)]


# Underscore-encoded comment suffixes that samtools sort may produce by
# converting in-QNAME spaces to underscores. Two patterns from our pipeline:
#   mapPacBio:                 '<uuid>_pt:i:<N>'
#   BBmap + retained comment:  '<accession>_<record_num>_length=<read_len>'
# We match specific forms (not bare '_<anything>') to avoid mangling qnames
# that legitimately contain underscores (e.g. Illumina-style flow-cell ids).
_UNDERSCORE_COMMENT_RE = re.compile(
    r'(?:_pt:i:\d+|_\d+_length=\d+).*$'
)


def _normalize_bam_read_name(name: str) -> str:
    """Strip aligner-specific comment suffixes from BAM read names.

    Two failure modes this guards against in the K-way consensus merge:

    1. **BBmap-retained FASTQ comments** (the v0.9.1 wt_R1 bug). BBmap by
       default emits the full FASTQ header line into QNAME, including the
       trailing comment: ``'SRR22434624.1654499 1654499 length=76'``.
       BWA truncates to the first whitespace token: ``'SRR22434624.1654499'``.
       Unnormalized, the two never join in consensus → no deduplication.

    2. **mapPacBio pt:i:N suffix** (pre-existing). mapPacBio embeds the
       poly-A tail length from the FASTQ auxiliary tag into the BAM read
       name. The separator depends on processing stage:
       - Space-separated ('UUID pt:i:25'): direct mapPacBio output, pre-sort.
       - Underscore-separated ('UUID_pt:i:25'): after samtools sort (BAM spec
         forbids spaces in QNAME, so samtools converts space → underscore).

    Strategy: strip everything after the first whitespace (handles BBmap's
    retained comment and ONT runid/sampleid comments while qnames are still
    in raw post-aligner form), then strip specific underscore-encoded
    suffixes produced by samtools sort.
    """
    # Strip after first whitespace (space or tab).
    for sep in (' ', '\t'):
        idx = name.find(sep)
        if idx != -1:
            name = name[:idx]
    # Strip recognized underscore-encoded suffixes left by samtools sort.
    return _UNDERSCORE_COMMENT_RE.sub('', name)


def _check_read_name_compatibility(
    bam_paths: Dict[str, str], n_sample: int = 100
) -> None:
    """Sample the first N reads from each input BAM and verify normalized
    names overlap meaningfully across aligners.

    Catches the v0.9.1 wt_R1 failure mode where BBmap retained FASTQ
    comments while BWA truncated them — the K-way merge then silently
    emitted both rows per read instead of one consensus winner.

    Raises with a clear error if normalized-name overlap is <50% across
    any pair of aligners, suggesting a read_id format mismatch that the
    consensus join would never recover from.
    """
    if len(bam_paths) < 2:
        return
    samples: Dict[str, set] = {}
    for aligner, path in bam_paths.items():
        names = set()
        with pysam.AlignmentFile(path, 'rb') as bam:
            for i, read in enumerate(bam):
                if i >= n_sample:
                    break
                names.add(_normalize_bam_read_name(read.query_name or ''))
        samples[aligner] = names
    aligners = list(samples.keys())
    for i, a in enumerate(aligners):
        for b in aligners[i + 1:]:
            inter = samples[a] & samples[b]
            union = samples[a] | samples[b]
            if not union:
                continue
            overlap = len(inter) / len(union)
            if overlap < 0.5:
                raise RuntimeError(
                    f"Read-id format mismatch between aligners {a!r} and {b!r}: "
                    f"only {len(inter)}/{len(union)} ({100*overlap:.1f}%) of the "
                    f"first {n_sample} normalized read names overlap. The K-way "
                    f"consensus merge will not match reads across these BAMs. "
                    f"Inspect QNAME formats in the input BAMs; if BBmap retained "
                    f"a FASTQ comment, re-run BBmap with 'trd=t' or extend "
                    f"_normalize_bam_read_name() to strip the offending suffix.\n"
                    f"  {a} sample: {sorted(samples[a])[:3]}\n"
                    f"  {b} sample: {sorted(samples[b])[:3]}"
                )


def _iter_name_grouped_bams(bam_paths: Dict[str, str]):
    """
    K-way merge across name-sorted BAMs, yielding all alignments per read.

    Memory: O(n_aligners) per read instead of O(total_reads * n_aligners).
    """
    _check_read_name_compatibility(bam_paths)
    bams = {}
    iterators = {}
    for aligner, path in bam_paths.items():
        bam = pysam.AlignmentFile(path, 'rb')
        bams[aligner] = bam
        iterators[aligner] = _filtered_read_iterator(bam)

    current_reads = {}
    for aligner, it in iterators.items():
        try:
            current_reads[aligner] = next(it)
        except StopIteration:
            current_reads[aligner] = None

    try:
        while any(r is not None for r in current_reads.values()):
            # Use normalized names as the merge key so mapPacBio's UUID_pt:i:N
            # reads group with bare UUID reads from other aligners.
            min_read_id = min(
                (_normalize_bam_read_name(r.query_name)
                 for r in current_reads.values() if r is not None),
                key=_natural_sort_key,
            )
            group = {}
            for aligner in list(current_reads.keys()):
                read = current_reads[aligner]
                if read is not None and _normalize_bam_read_name(read.query_name) == min_read_id:
                    # Drain ALL same-name records from this aligner (handles the
                    # rare case where mapPacBio emits two primary alignments for
                    # the same read, e.g. a truncated exon-1-only alignment plus
                    # the full spliced alignment).  Among duplicates, prefer the
                    # alignment with the most N-ops (most complete splice chain),
                    # then highest MAPQ.
                    candidates = [read]
                    try:
                        nxt = next(iterators[aligner])
                        while (nxt is not None
                               and _normalize_bam_read_name(nxt.query_name) == min_read_id):
                            candidates.append(nxt)
                            nxt = next(iterators[aligner])
                        current_reads[aligner] = nxt
                    except StopIteration:
                        current_reads[aligner] = None
                    if len(candidates) == 1:
                        group[aligner] = candidates[0]
                    else:
                        def _n_ops(r: pysam.AlignedSegment) -> int:
                            return sum(1 for op, _ in (r.cigartuples or []) if op == 3)
                        group[aligner] = max(
                            candidates,
                            key=lambda r: _n_ops(r),
                        )
            yield min_read_id, group
    finally:
        for bam in bams.values():
            bam.close()


def _cigar_query_length(read: pysam.AlignedSegment) -> int:
    """Return the total number of query-consuming bases implied by the CIGAR."""
    if not read.cigartuples:
        return 0
    # ops that consume query: M=0, I=1, S=4, =7, X=8
    query_ops = {0, 1, 4, 7, 8}
    return sum(length for op, length in read.cigartuples if op in query_ops)


def _restore_sequence_from_aligner_reads(
    best_read: pysam.AlignedSegment,
    aligner_reads: Dict[str, pysam.AlignedSegment],
) -> None:
    """Copy query sequence and quality scores to best_read from another aligner.

    gapmm2 outputs PAF which carries no read sequence, so _paf_to_bam() leaves
    query_sequence=None on every gapmm2 BAM record.  When gapmm2 wins consensus
    selection the output BAM would contain SEQ=* records that break all
    downstream steps (indel correction, poly-A trimming, etc.).

    This function looks through the other aligners' reads for the same read_id
    and copies the first non-None sequence whose length matches the CIGAR's
    expected query length.  Donors with mismatched lengths (e.g. hard-clipped
    records from deSALT) are skipped to prevent samtools from rejecting the
    BAM with "CIGAR and query sequence lengths differ".

    If no donor with the correct length is found, best_read is left as SEQ=*
    and a warning is logged.

    Args:
        best_read: The winning pysam.AlignedSegment (modified in place).
        aligner_reads: Dict mapping aligner name to pysam.AlignedSegment for
                       the same read_id.
    """
    expected_len = _cigar_query_length(best_read)
    for donor_read in aligner_reads.values():
        seq = donor_read.query_sequence
        if seq is None:
            continue
        if expected_len > 0 and len(seq) != expected_len:
            logger.debug(
                f"Skipping donor for '{best_read.query_name}': "
                f"sequence length {len(seq)} != CIGAR query length {expected_len}"
            )
            continue
        best_read.query_sequence = seq
        best_read.query_qualities = donor_read.query_qualities
        return
    logger.warning(
        f"No aligner has query_sequence for read '{best_read.query_name}'; "
        "writing SEQ=* record"
    )


def _process_and_write_batch(read_batch, raw_read_batch, genome, annotated_junctions, out_bam, stats, use_chimeric=False):
    """Process a batch of reads and write best alignments to output BAM."""
    if use_chimeric:
        from .chimeric_consensus import select_best_chimeric, build_chimeric_read

    for i, (read_id, alignments) in enumerate(read_batch):
        _, aligner_reads = raw_read_batch[i]

        if use_chimeric:
            chimeric_result = select_best_chimeric(aligner_reads, genome, annotated_junctions)

            # Pick a template read with a valid sequence (gapmm2 yields None).
            # Pass 1: prefer a read whose sequence length matches the chimeric CIGAR
            # query length (prevents "CIGAR and query sequence lengths differ" crash).
            # Pass 2 (fallback): accept any read with a sequence; the CIGAR may be
            # slightly wrong in length, but losing the read entirely is worse.
            query_ops = {0, 1, 4, 7, 8}  # M, I, S, =, X
            expected_len = sum(
                length for op, length in chimeric_result.chimeric_cigar if op in query_ops
            ) if chimeric_result.chimeric_cigar else 0
            template = None
            for r in aligner_reads.values():
                seq = r.query_sequence
                if seq is not None and (expected_len == 0 or len(seq) == expected_len):
                    template = r
                    break
            if template is None:
                # Fallback: accept any read with a sequence even if length mismatches.
                for r in aligner_reads.values():
                    if r.query_sequence is not None:
                        template = r
                        break
            if template is None:
                logger.warning(
                    f"No valid template read for chimeric assembly of '{read_id}'; skipping"
                )
                continue

            out_read = build_chimeric_read(
                template_read=template,
                ref_start=chimeric_result.chimeric_ref_start,
                cigar_tuples=chimeric_result.chimeric_cigar,
                chimeric_result=chimeric_result,
                header=out_bam.header,
            )
            out_read.flag &= ~0x900  # enforce primary

            # Validate CIGAR/sequence consistency.  The Pass-2 fallback above can
            # produce a mismatched template whose sequence length differs from the
            # chimeric CIGAR query span; pysam writes it silently but samtools sort
            # then crashes with "truncated file".  Drop the read rather than corrupt
            # the output BAM.
            if out_read.query_sequence is not None and out_read.cigartuples is not None:
                _q_ops = {0, 1, 4, 7, 8}  # M I S = X
                _cig_span = sum(l for op, l in out_read.cigartuples if op in _q_ops)
                if _cig_span != len(out_read.query_sequence):
                    logger.warning(
                        f"Chimeric read '{read_id}' CIGAR/sequence mismatch "
                        f"(cigar_query_span={_cig_span}, "
                        f"seq_len={len(out_read.query_sequence)}); dropping read"
                    )
                    continue

            if chimeric_result.confidence == 'high':
                stats['consensus_high'] += 1
            elif chimeric_result.confidence == 'medium':
                stats['consensus_medium'] += 1
            else:
                stats['consensus_low'] += 1
            if chimeric_result.is_chimeric:
                stats['chimeric_reads'] += 1
            for _pos, winner, _qs, _qe in (chimeric_result.segment_winners or []):
                stats['by_aligner'][winner] += 1
            unique_winners = frozenset(
                w[1] for w in (chimeric_result.segment_winners or [])
            )
            stats['by_aligner_combo'][unique_winners] += 1

            out_bam.write(out_read)

        else:
            result = select_best_alignment(alignments, genome, annotated_junctions)
            if result.confidence == 'high':
                stats['consensus_high'] += 1
            elif result.confidence == 'medium':
                stats['consensus_medium'] += 1
            else:
                stats['consensus_low'] += 1
            if result.was_5prime_rescued:
                stats['5prime_rescued'] += 1
            if result.n_tied_score > 1:
                stats['tied_score'] += 1
            stats['by_aligner'][result.best_aligner] += 1
            stats['by_aligner_combo'][frozenset(result.aligners_compared)] += 1

            if result.best_aligner in aligner_reads:
                best_read = aligner_reads[result.best_aligner]

                # Enforce exactly one primary per read: clear secondary (0x100) and
                # supplementary (0x800) bits so the winning record is always primary.
                best_read.flag &= ~0x900

                # Normalize mapPacBio's 'UUID pt:i:N' / 'UUID_pt:i:N' read name to bare
                # UUID so all downstream tools see a consistent read name regardless of
                # which aligner won the consensus selection.
                qn = best_read.query_name or ''
                if ' pt:i:' in qn or '_pt:i:' in qn:
                    best_read.query_name = _normalize_bam_read_name(best_read.query_name)

                # gapmm2 PAF→BAM conversion does not preserve read sequences;
                # restore SEQ from another aligner's record for the same read.
                if best_read.query_sequence is None:
                    _restore_sequence_from_aligner_reads(best_read, aligner_reads)

                # Aligner-selection metadata tags — lowercase second-letter to
                # avoid colliding with X[upper] tags emitted by upstream tools
                # (e.g. the cDNA pipeline writes XU=UMI, XC=cluster_size,
                # XA=tail_len, XR=read_ids, XF=full-length tier on FASTQ
                # comments propagated via minimap2 `-y`).
                best_read.set_tag('Xa', result.best_aligner)
                best_read.set_tag('Xc', result.confidence)
                best_read.set_tag('Xn', result.n_aligners_agree)
                if result.tied_aligners:
                    best_read.set_tag('Xt', ','.join(sorted(result.tied_aligners)))
                if result.was_5prime_rescued:
                    best_read.set_tag('Xj', 1)
                if result.false_junction_removed:
                    best_read.set_tag('Xv', 1)
                out_bam.write(best_read)


def run_consensus_selection(
    bam_paths: Dict[str, str],
    genome: Dict[str, str],
    output_bam: str,
    annotated_junctions: Optional[Set[Tuple[str, int, int, str]]] = None,
    write_all_to_tag: bool = True,
    n_workers: int = 0,
    batch_size: int = 10000,
    slurm_array_task: Optional[int] = None,
    slurm_array_total: Optional[int] = None,
    use_chimeric: bool = False,
    checkpoint_dir: Optional[str] = None,
) -> Dict[str, int]:
    """
    Run consensus selection across multiple BAM files.

    Streams through name-sorted BAMs to avoid loading all reads into memory.
    Supports SLURM array job splitting for cluster-scale parallelism.

    Memory usage: O(batch_size * n_aligners) instead of O(total_reads * n_aligners).

    Args:
        bam_paths: Dict mapping aligner name to BAM path
        genome: Dict mapping chrom to sequence
        output_bam: Output path for rectified BAM
        annotated_junctions: Optional set of annotated junctions
        write_all_to_tag: If True, write all aligner info to BAM tags
        n_workers: Number of worker processes (0 = auto-detect, 1 = single-threaded)
        batch_size: Number of read groups to accumulate before processing
        slurm_array_task: Current SLURM array task ID (0-indexed).
                          When set, only reads where
                          hash(read_id) % slurm_array_total == slurm_array_task
                          are processed.
        slurm_array_total: Total number of SLURM array tasks.

    Returns:
        Summary statistics dict
    """
    from ...slurm import get_available_cpus, get_slurm_info

    # Auto-detect SLURM array settings from environment.
    # Only activates when RECTIFY_CONSENSUS_ARRAY_MODE=1 is explicitly set.
    # In run-all mode, SLURM array indices are for sample parallelism, not
    # read-level partitioning — do not auto-activate there.
    if slurm_array_task is None and slurm_array_total is None:
        slurm_info = get_slurm_info()
        if (slurm_info.get('array_task_id') is not None
                and os.environ.get('RECTIFY_CONSENSUS_ARRAY_MODE') == '1'):
            try:
                slurm_array_task = int(slurm_info['array_task_id'])
                slurm_array_total = int(os.environ.get(
                    'SLURM_ARRAY_TASK_COUNT',
                    os.environ.get('SLURM_ARRAY_TASK_MAX', '0')
                ))
                if slurm_array_total > 0:
                    task_min = int(os.environ.get('SLURM_ARRAY_TASK_MIN', '0'))
                    task_step = int(os.environ.get('SLURM_ARRAY_TASK_STEP', '1'))
                    if 'SLURM_ARRAY_TASK_COUNT' not in os.environ:
                        slurm_array_total = (slurm_array_total - task_min) // task_step + 1
                    logger.info(
                        f"SLURM array detected: task {slurm_array_task} of {slurm_array_total}"
                    )
                else:
                    slurm_array_task = None
                    slurm_array_total = None
            except (ValueError, TypeError):
                slurm_array_task = None
                slurm_array_total = None

    use_slurm_filter = (
        slurm_array_task is not None and
        slurm_array_total is not None and
        slurm_array_total > 1
    )

    # Auto-detect workers
    if n_workers <= 0:
        n_workers = get_available_cpus()

    # Ensure BAMs are name-sorted
    import time as _time
    _t_total = _time.perf_counter()
    logger.info("Ensuring BAMs are name-sorted...")
    _t_ns = _time.perf_counter()
    sorted_bam_paths = {}
    for aligner, path in bam_paths.items():
        sorted_bam_paths[aligner] = _ensure_name_sorted(path)
    logger.info(f"[TIMING] Name-sort: {_time.perf_counter() - _t_ns:.1f}s")

    # Get header from first BAM
    first_bam_path = list(sorted_bam_paths.values())[0]
    first_bam = pysam.AlignmentFile(first_bam_path, 'rb')
    header = first_bam.header.to_dict()
    first_bam.close()

    # Add program group for RECTIFY consensus
    if 'PG' not in header:
        header['PG'] = []
    header['PG'].append({
        'ID': 'RECTIFY',
        'PN': 'RECTIFY',
        'VN': '2.0',
        'CL': f'consensus selection from {",".join(bam_paths.keys())}',
    })

    # Modify output path for SLURM array tasks
    if use_slurm_filter:
        base, ext = os.path.splitext(output_bam)
        output_bam = f"{base}.task{slurm_array_task}{ext}"
        logger.info(f"SLURM array task {slurm_array_task}: writing to {output_bam}")

    # Initialize stats
    stats = {
        'total_reads': 0,
        'reads_skipped_slurm_filter': 0,
        'consensus_high': 0,
        'consensus_medium': 0,
        'consensus_low': 0,
        '5prime_rescued': 0,
        'tied_score': 0,
        'chimeric_reads': 0,
        'by_aligner': defaultdict(int),
        'by_aligner_combo': defaultdict(int),  # frozenset of available aligners → count
    }

    # ── Checkpoint setup ─────────────────────────────────────────────────────
    # When checkpoint_dir is set:
    #   - Per-batch BAMs are written to checkpoint_dir/consensus_batch_NNNNNN.bam
    #   - A sentinel .done file marks each completed batch
    #   - consensus_checkpoint.json stores stats for resume
    #   - Final coordinate-sort + index run on checkpoint_dir (scratch), then
    #     the sorted BAM is copied to output_bam (Oak) to avoid NFS write hangs
    # When checkpoint_dir is None: original single-file behaviour (may hang on Oak NFS)
    _bam_header = pysam.AlignmentHeader.from_dict(header)  # reused for per-batch opens
    _batch_bam_paths: List[str] = []
    n_reads_to_skip = 0
    _ckpt_batch_num = 0

    def _write_ckpt_json() -> None:
        import json as _json
        ckpt = {
            'n_batches': _ckpt_batch_num,
            'total_reads_processed': stats['total_reads'],
            'stats': {
                k: (dict(v) if isinstance(v, defaultdict) else v)
                for k, v in stats.items()
                if k != 'by_aligner_combo'  # frozenset keys not JSON-serialisable
            },
        }
        _p = os.path.join(checkpoint_dir, 'consensus_checkpoint.json')
        with open(_p, 'w') as _f:
            _json.dump(ckpt, _f, indent=2)

    if checkpoint_dir:
        import glob as _glob, json as _json_ckpt
        os.makedirs(checkpoint_dir, exist_ok=True)

        # Collect completed batches from a previous (interrupted) run
        _done_files = sorted(_glob.glob(
            os.path.join(checkpoint_dir, 'consensus_batch_*.done')
        ))
        for _df in _done_files:
            _bam_f = _df.replace('.done', '.bam')
            if os.path.exists(_bam_f):
                _batch_bam_paths.append(_bam_f)
        _ckpt_batch_num = len(_batch_bam_paths)

        # Restore stats from checkpoint JSON if present
        _ckpt_json_path = os.path.join(checkpoint_dir, 'consensus_checkpoint.json')
        if _ckpt_batch_num > 0 and os.path.exists(_ckpt_json_path):
            with open(_ckpt_json_path) as _f:
                _ckpt_data = _json_ckpt.load(_f)
            # Use stored total_reads_processed as authoritative skip count
            n_reads_to_skip = _ckpt_data.get('total_reads_processed',
                                              _ckpt_batch_num * batch_size)
            _saved = _ckpt_data.get('stats', {})
            for _k, _v in _saved.items():
                if _k in stats:
                    if isinstance(stats[_k], defaultdict) and isinstance(_v, dict):
                        stats[_k].update(_v)
                    else:
                        stats[_k] = _v
            logger.info(
                f"Resuming consensus from checkpoint: {_ckpt_batch_num} batches done, "
                f"{n_reads_to_skip} reads to skip"
            )
        elif _ckpt_batch_num > 0:
            # No JSON but .done files exist → all completed batches are full
            n_reads_to_skip = _ckpt_batch_num * batch_size
            logger.info(
                f"Resuming consensus (no checkpoint JSON): "
                f"{_ckpt_batch_num} batches done, skipping {n_reads_to_skip} reads"
            )

        out_bam = None  # opened lazily at the start of each new batch
    else:
        # Original behaviour: single output file opened for the whole run
        out_bam = pysam.AlignmentFile(output_bam, 'wb', header=_bam_header)

    try:
        # Stream through name-sorted BAMs
        _t_stream = _time.perf_counter()
        logger.info(f"Streaming consensus selection (batch_size={batch_size})...")
        if use_slurm_filter:
            logger.info(
                f"  SLURM array filter: task {slurm_array_task}/{slurm_array_total}"
            )
        if checkpoint_dir and n_reads_to_skip > 0:
            logger.info(f"  Skipping {n_reads_to_skip} already-processed reads (resume)")

        # Accumulate batches for processing
        read_batch = []
        raw_read_batch = []
        n_batches = 0
        _n_skipped = 0  # reads skipped for checkpoint resume

        for read_id, aligner_reads in _iter_name_grouped_bams(sorted_bam_paths):
            # SLURM array filtering
            if use_slurm_filter:
                if _read_id_hash(read_id, slurm_array_total) != slurm_array_task:
                    stats['reads_skipped_slurm_filter'] += 1
                    continue

            # Skip already-checkpointed reads on resume
            if _n_skipped < n_reads_to_skip:
                _n_skipped += 1
                continue

            stats['total_reads'] += 1

            # Lazily open the current batch BAM when the first read of a new batch arrives
            if checkpoint_dir and out_bam is None:
                _cur_batch_path = os.path.join(
                    checkpoint_dir, f'consensus_batch_{_ckpt_batch_num:06d}.bam'
                )
                out_bam = pysam.AlignmentFile(_cur_batch_path, 'wb', header=_bam_header)

            # Extract alignment info for scoring (non-chimeric path only).
            # select_best_chimeric operates directly on raw pysam reads; it never
            # consults AlignmentInfo objects.  Calling extract_alignment_info in
            # chimeric mode costs 1M expensive pysam operations for nothing.
            if not use_chimeric:
                alignments = {
                    aligner: extract_alignment_info(read, aligner, genome)
                    for aligner, read in aligner_reads.items()
                }
            else:
                alignments = {}  # unused by _process_and_write_batch in chimeric mode

            read_batch.append((read_id, alignments))
            raw_read_batch.append((read_id, aligner_reads))

            # Process batch when full
            if len(read_batch) >= batch_size:
                _process_and_write_batch(
                    read_batch, raw_read_batch, genome,
                    annotated_junctions, out_bam, stats,
                    use_chimeric=use_chimeric,
                )
                read_batch = []
                raw_read_batch = []
                n_batches += 1

                if checkpoint_dir:
                    out_bam.close()
                    _sentinel = os.path.join(
                        checkpoint_dir, f'consensus_batch_{_ckpt_batch_num:06d}.done'
                    )
                    open(_sentinel, 'w').close()
                    _batch_bam_paths.append(_cur_batch_path)
                    _ckpt_batch_num += 1
                    out_bam = None  # will be lazily opened for next batch
                    _write_ckpt_json()
                    logger.info(
                        f"  Checkpoint: batch {_ckpt_batch_num} written "
                        f"({stats['total_reads']:,} reads total)"
                    )

                if stats['total_reads'] % 100000 == 0:
                    logger.info(f"  Processed {stats['total_reads']:,} reads...")

        # Process remaining reads
        if read_batch:
            if checkpoint_dir and out_bam is None:
                _cur_batch_path = os.path.join(
                    checkpoint_dir, f'consensus_batch_{_ckpt_batch_num:06d}.bam'
                )
                out_bam = pysam.AlignmentFile(_cur_batch_path, 'wb', header=_bam_header)
            _process_and_write_batch(
                read_batch, raw_read_batch, genome,
                annotated_junctions, out_bam, stats,
                use_chimeric=use_chimeric,
            )
            n_batches += 1
            if checkpoint_dir:
                out_bam.close()
                _sentinel = os.path.join(
                    checkpoint_dir, f'consensus_batch_{_ckpt_batch_num:06d}.done'
                )
                open(_sentinel, 'w').close()
                _batch_bam_paths.append(_cur_batch_path)
                _ckpt_batch_num += 1
                out_bam = None
                _write_ckpt_json()

    except Exception:
        if out_bam is not None:
            try:
                out_bam.close()
            except Exception:
                pass
        if not checkpoint_dir:
            # Remove partial output BAM so callers don't see an incomplete file
            try:
                os.unlink(output_bam)
            except OSError:
                pass
        # With checkpoint_dir: leave batch files intact for resume
        raise

    # Close output (non-checkpoint path)
    if out_bam is not None:
        out_bam.close()
    logger.info(f"[TIMING] Streaming ({stats['total_reads']:,} reads, {n_batches} batches): {_time.perf_counter() - _t_stream:.1f}s")
    if stats['total_reads'] > 0:
        _reads_per_sec = stats['total_reads'] / max(_time.perf_counter() - _t_stream, 0.001)
        logger.info(f"[TIMING] Throughput: {_reads_per_sec:,.0f} reads/sec")

    # Sort output by coordinate and index
    _t_sort = _time.perf_counter()
    logger.info("Coordinate-sorting output BAM...")
    if checkpoint_dir:
        # ── Scratch-based sort: merge all batch BAMs on scratch, sort, index,
        #    then copy the final sorted BAM to output_bam (Oak).
        #    This avoids pysam.sort writing directly to Oak NFS which can hang
        #    under concurrent array-job I/O load.
        import shutil as _shutil_sort
        _merged_path = os.path.join(checkpoint_dir, 'consensus_merged.bam')
        _sorted_path = os.path.join(checkpoint_dir, 'consensus_sorted.bam')

        if not _batch_bam_paths:
            # No reads processed at all — write an empty coordinate-sorted BAM
            _empty = pysam.AlignmentFile(_sorted_path, 'wb', header=_bam_header)
            _empty.close()
        elif len(_batch_bam_paths) == 1:
            pysam.sort('-m', '2G', '-o', _sorted_path, _batch_bam_paths[0])
        else:
            # Concatenate (samtools cat) then coordinate-sort
            pysam.cat('-o', _merged_path, *_batch_bam_paths)
            pysam.sort('-m', '2G', '-o', _sorted_path, _merged_path)
            try:
                os.unlink(_merged_path)
            except OSError:
                pass

        pysam.index(_sorted_path)

        # Copy sorted BAM + index to Oak output path
        _shutil_sort.copy2(_sorted_path, output_bam)
        _bai_src = _sorted_path + '.bai'
        if os.path.exists(_bai_src):
            _shutil_sort.copy2(_bai_src, output_bam + '.bai')

        # Clean up scratch sort files
        for _p in (_sorted_path, _bai_src):
            try:
                os.unlink(_p)
            except OSError:
                pass
        # Clean up per-batch BAMs (already merged+sorted)
        for _b in _batch_bam_paths:
            try:
                os.unlink(_b)
            except OSError:
                pass
    else:
        # Original behaviour: sort in-place at output_bam location
        # Pre-sort validation: sample first 1000 records to detect CIGAR/sequence mismatches
        # (avoid expensive sort failure on corrupt BAM written by batch writer)
        _mismatch_count = 0
        try:
            with pysam.AlignmentFile(output_bam, 'rb') as _check_bam:
                for _i, _r in enumerate(_check_bam.fetch(until_eof=True)):
                    if _i >= 1000:
                        break
                    if _r.query_sequence is not None and _r.cigartuples is not None:
                        _q_ops = {0, 1, 4, 7, 8}  # M I S = X
                        _cig_span = sum(l for op, l in _r.cigartuples if op in _q_ops)
                        if _cig_span != len(_r.query_sequence):
                            _mismatch_count += 1
            if _mismatch_count > 0:
                raise RuntimeError(
                    f"Pre-sort validation failed: {_mismatch_count}/1000 sampled reads have "
                    f"CIGAR/sequence length mismatches in {output_bam} — BAM is corrupt/truncated"
                )
        except RuntimeError:
            raise
        except Exception as e:
            logger.warning(f"Pre-sort CIGAR validation skipped: {e}")

        sorted_output = output_bam.replace('.bam', '.sorted.bam')
        pysam.sort('-m', '1G', '-o', sorted_output, output_bam)
        os.replace(sorted_output, output_bam)
        pysam.index(output_bam)
    logger.info(f"[TIMING] Coordinate-sort + index: {_time.perf_counter() - _t_sort:.1f}s")

    # Drop a provenance sidecar so a future run-all reuse gate can verify the
    # rectified BAM was produced by the same rectify version before reusing it.
    try:
        import sys as _sys_prov
        from rectify.utils.bam_provenance import compute_run_provenance, write_sidecar
        _prov = compute_run_provenance(command=_sys_prov.argv)
        write_sidecar(output_bam, _prov, aligner_name="consensus")
    except Exception as _prov_err:
        logger.warning("Failed to write provenance sidecar for rectified BAM (%s); "
                       "BAM emitted without sidecar", _prov_err)

    # Log summary
    logger.info(f"\nConsensus selection complete:")
    logger.info(f"  Total reads processed: {stats['total_reads']}")
    if use_slurm_filter:
        logger.info(f"  Reads skipped (other SLURM tasks): {stats['reads_skipped_slurm_filter']}")
    logger.info(f"  High confidence: {stats['consensus_high']}")
    logger.info(f"  Medium confidence: {stats['consensus_medium']}")
    logger.info(f"  Low confidence: {stats['consensus_low']}")
    logger.info(f"  5' rescued: {stats['5prime_rescued']}")
    logger.info(f"  Tied score (tiebreaker used): {stats['tied_score']}")
    if use_chimeric:
        logger.info(f"  Chimeric reads (multi-aligner segments): {stats['chimeric_reads']}")
    logger.info(f"  By aligner: {dict(stats['by_aligner'])}")
    logger.info(f"  By aligner combo: { {'+'.join(sorted(k)): v for k, v in stats['by_aligner_combo'].items()} }")
    logger.info(f"  Batches processed: {n_batches}")
    logger.info(f"[TIMING] run_consensus_selection total: {_time.perf_counter() - _t_total:.1f}s")

    return stats


def merge_slurm_array_bams(
    output_bam_pattern: str,
    n_tasks: int,
    merged_output: str,
):
    """
    Merge BAM files from SLURM array tasks into a single output.

    Call this after all array tasks have completed.

    Args:
        output_bam_pattern: Pattern with {task} placeholder
        n_tasks: Number of array tasks
        merged_output: Path for merged output BAM
    """
    task_bams = []
    for task_id in range(n_tasks):
        bam_path = output_bam_pattern.format(task=task_id)
        if os.path.exists(bam_path):
            task_bams.append(bam_path)
        else:
            logger.warning(f"Missing SLURM array task BAM: {bam_path}")

    if not task_bams:
        raise FileNotFoundError("No SLURM array task BAMs found")

    logger.info(f"Merging {len(task_bams)} SLURM array task BAMs...")
    pysam.merge('-f', merged_output, *task_bams)

    sorted_output = merged_output.replace('.bam', '.sorted.bam')
    pysam.sort('-m', '1G', '-o', sorted_output, merged_output)
    os.replace(sorted_output, merged_output)
    pysam.index(merged_output)

    logger.info(f"Merged output: {merged_output}")

    for bam_path in task_bams:
        idx_path = bam_path + '.bai'
        if os.path.exists(idx_path):
            os.remove(idx_path)
        os.remove(bam_path)

    logger.info("SLURM array merge complete")


def load_annotated_junctions(annotation_path: str) -> Set[Tuple[str, int, int, str]]:
    """
    Load annotated junctions from GFF/GTF file.

    Returns set of (chrom, intron_start, intron_end, strand) tuples where chrom is in
    standardized canonical format (chrI, chrII, etc.) so that junction lookups
    match the standardized chrom names used during correction.
    """
    from ...utils.genome import standardize_chrom_name

    junctions = set()

    import gzip as _gzip
    _open = _gzip.open if str(annotation_path).endswith('.gz') else open
    with _open(annotation_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2].lower()

            # Look for intron features — match any subtype (intron,
            # five_prime_UTR_intron, three_prime_UTR_intron, etc.)
            if 'intron' in feature_type:
                chrom = standardize_chrom_name(parts[0])
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])  # Already exclusive in GFF end
                strand = parts[6] if parts[6] in ('+', '-') else '+'
                junctions.add((chrom, start, end, strand))

    if len(junctions) == 0:
        logger.warning(
            "load_annotated_junctions: 0 junctions loaded from %s. "
            "Check that the file exists, is readable, and contains 'intron' "
            "feature records (column 3). Junction-guided scoring will be disabled.",
            annotation_path,
        )
    else:
        logger.info(f"Loaded {len(junctions)} annotated junctions from {annotation_path}")

    return junctions
