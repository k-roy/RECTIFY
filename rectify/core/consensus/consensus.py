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
2. Prefer junctions supported by multiple aligners
3. Tiebreaker: prefer aligner whose corrected 3' position agrees with majority
4. Tiebreaker: canonical splice site motifs (GT/AG) and annotated junctions

Note: A-tract 3' correction is applied to each aligner pre-scoring using genome
sequence only. Full indel correction (MD-tag dependent) is applied post-consensus
as a refinement step.

Author: Kevin R. Roy
"""

import logging
import os
import hashlib
import re
import shutil
from collections import defaultdict
from .intron_sanity import truncate_impossible_introns, max_reportable_intron_from_env
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


def _fsync_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path)) or '.'
    try:
        fd = os.open(parent, os.O_RDONLY)
    except OSError as exc:
        logger.debug("Could not open parent directory for fsync: %s (%s)", parent, exc)
        return
    try:
        os.fsync(fd)
    except OSError as exc:
        logger.debug("Could not fsync parent directory: %s (%s)", parent, exc)
    finally:
        os.close(fd)


def _fsync_file(path: str) -> None:
    with open(path, 'rb') as fh:
        os.fsync(fh.fileno())


def _atomic_write_text(path: str, text: str = '') -> None:
    tmp_path = f"{path}.tmp"
    with open(tmp_path, 'w') as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp_path, path)
    _fsync_parent_dir(path)


def _copy2_and_fsync(src: str, dst: str) -> None:
    shutil.copy2(src, dst)
    _fsync_file(dst)
    _fsync_parent_dir(dst)


def _validate_bam_sample(path: str, *, max_records: int = 1000) -> None:
    try:
        pysam.quickcheck(path)
    except Exception as exc:
        raise RuntimeError(f"BAM quickcheck failed for checkpoint batch {path}: {exc}") from exc

    mismatch_count = 0
    past_contig = []
    try:
        with pysam.AlignmentFile(path, 'rb') as bam:
            contig_lens = dict(zip(bam.references, bam.lengths))
            for i, read in enumerate(bam.fetch(until_eof=True)):
                if i >= max_records:
                    break
                if read.query_sequence is None or read.cigartuples is None:
                    continue
                query_ops = {0, 1, 4, 7, 8}  # M I S = X
                cigar_query_span = sum(length for op, length in read.cigartuples if op in query_ops)
                if cigar_query_span != len(read.query_sequence):
                    mismatch_count += 1
                # Invariant: an alignment may not run off EITHER edge of its
                # contig. `intron_sanity` truncates past-END records before
                # they are written, so this should be unreachable — which is
                # the point. Truncation fixes the alignments we know about;
                # this is what stops the class returning silently years from
                # now (arg. from the 668-drs-arm session). See planning/684c:
                # deSALT produced 36 past-contig alignments per 400k reads on
                # DRS before the fix. The check is TWO-SIDED (planning/719):
                # the 5' rescue once walked a POS-0 read 25 bp off the START
                # of chrXV, producing a negative POS that killed `samtools
                # index` a stage later — the original one-sided guard never
                # fired on it.
                if not read.is_unmapped and read.reference_end is not None:
                    clen = contig_lens.get(read.reference_name)
                    _off_start = read.reference_start < 0
                    _off_end = clen is not None and read.reference_end > clen
                    if _off_start or _off_end:
                        if len(past_contig) < 5:
                            past_contig.append(
                                f"{read.query_name} {read.reference_name}:"
                                f"{read.reference_start}-{read.reference_end}"
                                f" outside [0, {clen}]"
                            )
    except Exception as exc:
        raise RuntimeError(f"Could not read checkpoint batch BAM {path}: {exc}") from exc

    if mismatch_count:
        raise RuntimeError(
            f"Checkpoint batch validation failed for {path}: "
            f"{mismatch_count}/{max_records} sampled reads have CIGAR/sequence length mismatches"
        )
    if past_contig:
        raise RuntimeError(
            f"Checkpoint batch validation failed for {path}: alignment(s) whose "
            f"reference span runs OFF A CONTIG EDGE (past the end, or a negative "
            f"POS off the start) — malformed by definition; intron_sanity / the "
            f"5' rescue edge guard should have prevented them before write. Examples: "
            + '; '.join(past_contig)
        )


def _commit_checkpoint_batch(bam_path: str, done_path: str) -> None:
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"Checkpoint batch BAM does not exist: {bam_path}")
    _fsync_file(bam_path)
    _fsync_parent_dir(bam_path)
    _validate_bam_sample(bam_path)
    _atomic_write_text(done_path, 'ok\n')


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


def _ensure_rn_sorted(bam_path: str) -> str:
    """
    Ensure a BAM is sorted by its integer ``RN:i`` tag (the K-way merge key).

    The RN-keyed merge requires every input stream to yield RN in
    non-decreasing order. Name-sort does NOT guarantee that: RN is assigned
    in FASTQ record order (``_build_qname_rn_map``), and natural QNAME order
    is unrelated to FASTQ order for plain-uuid reads (ONT DRS). With
    name-sorted inputs the merge desynchronises wherever a read is missing
    from a subset of aligners and emits split groups (~2.2 primary records
    per read on the 668 DRS panel). Sorting by the merge key itself is the
    only ordering correct for every QNAME shape.
    """
    sorted_path = bam_path.replace('.bam', '.rnsorted.bam')
    if os.path.exists(sorted_path):
        if os.path.getmtime(sorted_path) > os.path.getmtime(bam_path):
            logger.info(f"Using existing RN-sorted BAM: {sorted_path}")
            return sorted_path

    logger.info(f"RN-sorting BAM: {bam_path} -> {sorted_path}")
    # -t RN sorts by the integer tag value; -n makes ties break by natural
    # QNAME (deterministic). Same -m rationale as _ensure_name_sorted.
    pysam.sort('-t', 'RN', '-n', '-m', '1G', '-o', sorted_path, bam_path)
    return sorted_path


def _read_id_hash(read_id: str, n_buckets: int) -> int:
    """Deterministic hash of read_id for SLURM array splitting."""
    h = hashlib.md5(str(read_id).encode()).hexdigest()
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
# converting in-QNAME whitespace (spaces or tabs) to underscores. Patterns
# observed in our pipeline:
#   mapPacBio:                  '<uuid>_pt:i:<N>'
#   BBmap + retained comment:   '<accession>_<record_num>_length=<read_len>'
#   tab-aux leak (minimap2 -y): '<uuid>_XA:Z:foo', '<uuid>_XC:i:42', etc.
#   Dorado metadata leak:       '<uuid>_runid=<v>_ch=<v>_start_time=<v>...'
#
# Matched specifically (not bare '_<anything>') to avoid mangling qnames that
# legitimately contain underscores (e.g. Illumina-style flow-cell ids).
#
# The generic SAM-aux pattern `_<2char>:[AZifHB]:` covers any underscore-encoded
# aux tag (XA:Z:, XC:i:, ts:A:, etc.) — defends against future aligner outputs
# without needing per-tag regex updates. The Dorado metadata keys are
# enumerated because their `=`-separated form is too easy to overmatch.
_UNDERSCORE_COMMENT_RE = re.compile(
    r'(?:'
    r'_pt:i:\d+'                                # mapPacBio poly-A annotation
    r'|_\d+_length=\d+'                         # BBmap + retained SRA-style comment
    r'|_[A-Za-z][A-Za-z0-9]:[AZifHB]:'          # generic SAM aux tag (XA:Z:, etc.)
    r'|_(?:runid|ch|start_time|sampleid'
    r'|flow_cell_id|read|parent_read_id'
    r'|model_version_id|barcode)='              # Dorado metadata keys
    r').*$'
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




def _first_primary_has_rn(bam_path: str) -> bool:
    """Return True if the first primary mapped record carries an RN tag."""
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            return read.has_tag('RN')
    return False


def _all_inputs_have_rn(bam_paths: Dict[str, str]) -> bool:
    if not bam_paths:
        return False
    return all(_first_primary_has_rn(path) for path in bam_paths.values())


def _merge_key(read: pysam.AlignedSegment, use_rn_key: bool):
    if use_rn_key:
        try:
            return int(read.get_tag('RN'))
        except KeyError:
            return None
    return _normalize_bam_read_name(read.query_name or '')


def _merge_key_sort_value(key):
    if isinstance(key, int):
        return key
    return _natural_sort_key(str(key))


def _parse_fastq_comment_tags(comment: str) -> List[Tuple[str, str, object]]:
    """Parse SAM-format FASTQ comment tags into ``(tag, type, value)`` tuples."""
    tags = []
    for tok in (comment or '').split():
        parts = tok.split(':', 2)
        if len(parts) != 3:
            continue
        tag, typ, raw = parts
        if len(tag) != 2 or tag == 'RN':
            continue
        try:
            if typ in ('i', 'I'):
                value = int(raw)
            elif typ == 'f':
                value = float(raw)
            elif typ in ('Z', 'H'):
                value = raw
            elif typ == 'A' and len(raw) == 1:
                value = raw
            else:
                continue
        except ValueError:
            continue
        tags.append((tag, typ, value))
    return tags


_CDNA_COMMENT_TAGS = ('XU', 'XO', 'XT', 'XY', 'XC', 'XF', 'XM', 'XB', 'XR', 'XA')


def _restore_comment_tags_from_siblings(best_read, aligner_reads):
    """Ensure the consensus winner carries the authoritative cDNA FASTQ-comment
    aux tags, copying them from a sibling aligner's record if the winner lacks
    them — so a read another aligner mapped is NEVER dropped downstream just
    because the *winning* aligner did not propagate the FASTQ comment.

    The cDNA Stage-1 writer puts per-molecule metadata (XU=UMI, XO=orient,
    XT/XY=read type, XC=cluster_size, XF=full-length tier, +XM/XB/XR/XA) on the
    FASTQ comment; these reach the BAM ONLY via an aligner run with FASTQ-comment
    pass-through (minimap2 / minimap2-family presets with ``-y``). A winner from a
    junction aligner that does not pass FASTQ comments (uLTRA) either LACKS these
    tags or carries its OWN colliding tags — uLTRA emits XA/XC with splice
    semantics (e.g. ``XC:Z:NO_SPLICE``), which then crash cdna-analyze's
    ``int(XC)`` and mis-key downstream joins.

    Robust for ANY aligner combination (not hardcoded to minimap2): the UMI tag
    ``XU`` is emitted only by the Stage-1 comment (uLTRA never produces it), so a
    record carrying ``XU`` is an authoritative comment-propagating source. If the
    winner already carries ``XU`` it is authoritative and left untouched; else we
    copy the full comment-tag set from a sibling that carries ``XU`` (preferring
    ``minimap2`` for determinism), OVERWRITING any colliding tags on the winner.
    If no sibling carries ``XU`` (a rare molecule mapped only by a non-propagating
    aligner) the winner is left as-is and cdna-analyze drops it via its defensive
    tag guard.
    """
    if best_read.has_tag('XU'):
        return  # winner is itself from a comment-propagating aligner
    # Pick an authoritative source sibling: any record with XU, minimap2 first.
    source = None
    ordered = sorted(
        (a for a in aligner_reads if aligner_reads[a] is not None),
        key=lambda a: (a != 'minimap2', a),
    )
    for a in ordered:
        rec = aligner_reads[a]
        if rec is not best_read and rec.has_tag('XU'):
            source = rec
            break
    if source is None:
        return
    for t in _CDNA_COMMENT_TAGS:
        if source.has_tag(t):
            val, vt = source.get_tag(t, with_value_type=True)
            best_read.set_tag(t, val, value_type=vt)


def _restore_sidecar_tags(read: pysam.AlignedSegment, read_num_sidecar) -> None:
    """Restore FASTQ-comment aux tags from a read-num sidecar, if available."""
    if read_num_sidecar is None:
        return
    try:
        if read.has_tag('RN'):
            row = read_num_sidecar.lookup(int(read.get_tag('RN')))
        else:
            row = read_num_sidecar.lookup_by_qname(
                _normalize_bam_read_name(read.query_name or '')
            )
    except KeyError:
        logger.warning(
            "Read-num sidecar has no row for consensus read %r; writing without restored tags",
            read.query_name,
        )
        return
    for tag, typ, value in _parse_fastq_comment_tags(row.fastq_comment):
        # 659: the FASTQ comment is the AUTHORITATIVE source for the Stage-1 cDNA tags, so
        # for those the sidecar must OVERWRITE a colliding aligner-emitted tag — uLTRA writes
        # its own XC:Z:NO_SPLICE / XA:Z: with splice semantics, and the old don't-clobber
        # guard let them survive; cdna-analyze's int(XC) then dropped every uLTRA-won record
        # as "missing required tags" (7.0% of rna15aa_rep1, BIASED toward spliced molecules).
        # Same rule as _restore_comment_tags_from_siblings. Non-cDNA tags keep the guard.
        if read.has_tag(tag) and tag not in _CDNA_COMMENT_TAGS:
            continue
        read.set_tag(tag, value, value_type=typ)

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


def _detect_duplicate_molecule_names(
    bam_paths: Dict[str, str], use_rn_key: bool, sample: int = 200000
) -> None:
    """Fail loud on the duplicate-molecule-name pathology that silently
    collapses the K-way merge.

    Both the integer-RN key and the normalized-QNAME key are derived 1:1 from
    the read name, so a Stage-1 FASTQ whose ``cluster_<cid>`` names collide
    across regions (e.g. per-chromosome ``correct-cdna --region`` output
    concatenated without a region prefix) makes DISTINCT molecules share one
    merge key. The merge then collapses ~N distinct molecules to one output
    record and drops the rest — historically ~87% of reads (see planning/251,
    /250a-c). The fix is globally-unique Stage-1 names; this guard ensures a
    regression fails immediately instead of silently.

    Detection is cheap and specific: scan up to ``sample`` primary reads from
    one name-sorted input BAM (colliding names are ADJACENT under name-sort) and
    flag names that appear on ≥2 UNPAIRED primary records at DIFFERENT
    (reference, position). Multimappers are excluded (only primary reads are
    iterated; the cDNA presets run ``--secondary=no``); paired reads are exempt
    (R1/R2 legitimately share a name and are mate-split downstream). Raises when
    >1% of the sampled distinct names collide to distinct loci.
    """
    if not bam_paths:
        return
    probe_path = next(iter(bam_paths.values()))
    name_loci: Dict[object, set] = defaultdict(set)
    n_seen = 0
    with pysam.AlignmentFile(probe_path, 'rb') as bam:
        for read in _filtered_read_iterator(bam):
            if read.is_paired:
                return  # paired data: name-sharing is expected; guard N/A
            key = _merge_key(read, use_rn_key)
            if key is None:
                continue
            name_loci[key].add((read.reference_id, read.reference_start))
            n_seen += 1
            if n_seen >= sample:
                break
    n_names = len(name_loci)
    if n_names == 0:
        return
    n_collide = sum(1 for loci in name_loci.values() if len(loci) > 1)
    frac = n_collide / n_names
    if frac > 0.01:
        raise RuntimeError(
            "Duplicate molecule names detected in the consensus input "
            f"({probe_path}): {n_collide}/{n_names} ({100 * frac:.1f}%) of the "
            f"first {n_seen} primary read names map to ≥2 distinct genomic loci. "
            "The K-way consensus merge keys on the read name/RN and would "
            "silently COLLAPSE these distinct molecules into one record "
            "(dropping the rest — the ~87% cDNA align2 loss). Root cause is "
            "non-globally-unique Stage-1 consensus names — e.g. per-region "
            "`rectify correct-cdna --region <chr>` output concatenated without a "
            "region prefix. Regenerate Stage-1 with globally-unique names "
            "(fixed: the --region path now prefixes cluster names with the "
            "region; see planning/251) and re-run align."
        )


def _iter_name_grouped_bams(bam_paths: Dict[str, str], use_rn_key: Optional[bool] = None):
    """
    K-way merge across merge-key-sorted BAMs, yielding all alignments per read.

    Memory: O(n_aligners) per read instead of O(total_reads * n_aligners).
    When every input BAM carries RN tags, the merge key is the integer RN
    and inputs MUST be RN-sorted (``_ensure_rn_sorted``); otherwise the
    normalized-QNAME key path is used and inputs must be name-sorted.

    The merge is only correct when each input stream yields keys in
    non-decreasing order — a stream sorted on a DIFFERENT key desynchronises
    wherever a read is missing from a subset of aligners, silently splitting
    (QNAME, RN) groups into multiple output records (the 668 DRS defect).
    A per-stream monotonicity guard therefore fails loud on out-of-order keys.
    """
    if use_rn_key is None:
        use_rn_key = _all_inputs_have_rn(bam_paths)
    if use_rn_key:
        logger.info("Consensus K-way merge: using RN:i keys")
    else:
        _check_read_name_compatibility(bam_paths)
        logger.info("Consensus K-way merge: using normalized QNAME keys")
    _detect_duplicate_molecule_names(bam_paths, use_rn_key)

    def _monotone_checked(aligner: str, path: str, it):
        """Yield reads, raising if the merge key ever decreases (fail loud
        instead of silently splitting groups)."""
        prev_sort_value = None
        prev_key = None
        for read in it:
            key = _merge_key(read, use_rn_key)
            if key is not None:
                sort_value = _merge_key_sort_value(key)
                if prev_sort_value is not None and sort_value < prev_sort_value:
                    raise RuntimeError(
                        f"Consensus input BAM for aligner {aligner!r} ({path}) is "
                        f"not sorted by the K-way merge key "
                        f"({'RN:i' if use_rn_key else 'normalized QNAME'}): key "
                        f"{key!r} follows {prev_key!r}. The merge would silently "
                        f"split read groups and emit multiple primary records per "
                        f"read. Re-sort the input with "
                        f"{'samtools sort -t RN -n' if use_rn_key else 'samtools sort -n'} "
                        f"(rectify does this automatically in run_consensus_selection)."
                    )
                prev_sort_value = sort_value
                prev_key = key
            yield read

    bams = {}
    iterators = {}
    for aligner, path in bam_paths.items():
        bam = pysam.AlignmentFile(path, 'rb')
        bams[aligner] = bam
        iterators[aligner] = _monotone_checked(
            aligner, path, _filtered_read_iterator(bam))

    current_reads = {}
    for aligner, it in iterators.items():
        try:
            current_reads[aligner] = next(it)
        except StopIteration:
            current_reads[aligner] = None

    def _n_ops(r: pysam.AlignedSegment) -> int:
        return sum(1 for op, _ in (r.cigartuples or []) if op == 3)

    def _mate_key(r: pysam.AlignedSegment):
        # Paired short-read: R1 and R2 share one RN/QNAME but are DISTINCT reads
        # and must be adjudicated separately. Without this split the per-aligner
        # collapse below would keep only the more-spliced mate via max(_n_ops) and
        # silently drop the other — ~half the paired junction evidence, biased
        # against the less-spliced mate. Unpaired data → None (single bucket →
        # byte-identical to the original single-record-per-aligner behaviour).
        return r.is_read2 if r.is_paired else None

    try:
        while any(r is not None for r in current_reads.values()):
            min_read_id = min(
                (_merge_key(r, use_rn_key) for r in current_reads.values() if r is not None),
                key=_merge_key_sort_value,
            )
            # Collect ALL records with this key per aligner. Both mates of a pair
            # share the key and are adjacent under name-sort, so they are gathered
            # here together (this does NOT rely on cross-aligner mate ordering).
            per_aligner: Dict[str, List[pysam.AlignedSegment]] = {}
            for aligner in list(current_reads.keys()):
                read = current_reads[aligner]
                if read is not None and _merge_key(read, use_rn_key) == min_read_id:
                    candidates = [read]
                    try:
                        nxt = next(iterators[aligner])
                        while (nxt is not None
                               and _merge_key(nxt, use_rn_key) == min_read_id):
                            candidates.append(nxt)
                            nxt = next(iterators[aligner])
                        current_reads[aligner] = nxt
                    except StopIteration:
                        current_reads[aligner] = None
                    per_aligner[aligner] = candidates

            # Partition each aligner's records by mate, then emit ONE group per
            # mate (one record per aligner — the existing downstream contract).
            # Same-mate duplicates (e.g. an aligner's multimappers) still collapse
            # to the most-spliced record. Unpaired → a single None bucket.
            mate_buckets: Dict[object, Dict[str, pysam.AlignedSegment]] = defaultdict(dict)
            for aligner, cands in per_aligner.items():
                by_mate: Dict[object, List[pysam.AlignedSegment]] = defaultdict(list)
                for c in cands:
                    by_mate[_mate_key(c)].append(c)
                for mk, lst in by_mate.items():
                    mate_buckets[mk][aligner] = lst[0] if len(lst) == 1 else max(lst, key=_n_ops)

            # Deterministic mate order: None (unpaired) first, then R1 (False),
            # then R2 (True). Avoids None-vs-bool comparison errors.
            # Unpaired data yields the bare key (byte-identical to the original
            # contract); paired data yields (read_id, mate) so R1/R2 stay distinct.
            for mate_key in sorted(mate_buckets, key=lambda m: (m is not None, m)):
                out_key = min_read_id if mate_key is None else (min_read_id, mate_key)
                yield out_key, mate_buckets[mate_key]
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
        # pysam reverse-complements `query_sequence` for `is_reverse=True`
        # records, so a donor mapped on the opposite strand would inject
        # the RC of the correct sequence. Skip mismatched-strand donors.
        if donor_read.is_reverse != best_read.is_reverse:
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


def _credit_tied_aligners(stats, result) -> None:
    """Credit every aligner tied for the top junction score in the win tally.

    On a score tie, select.py picks a single ``best_aligner`` via an arbitrary
    ASCII-name fallback, but per Kevin's directive the tally must keep ALL tied
    top aligners — so ``aligner_wins_*`` reflects "was among the top-scoring
    aligners", not one arbitrary winner. ``tied_aligners == [best_aligner]`` on a
    unique win, so single-winner reads (and the empty-alignments degenerate case,
    ``tied_aligners == []`` / ``best_aligner == 'none'``, handled by the fallback)
    are counted exactly as before; an N-way tie credits all N. This can make
    ``by_aligner`` sum to >100% of reads — ``stats['tied_score']`` separately
    counts the multi-credited reads.

    Tally-only: this does NOT change which record is written to the consensus BAM
    (``best_aligner`` alone still selects the written alignment).
    """
    for aligner in (result.tied_aligners or [result.best_aligner]):
        stats['by_aligner'][aligner] += 1


def _enforce_intron_sanity(read, out_bam, max_intron_bp, stats=None):
    """Soft-clip a winning alignment at its first physically impossible intron.

    uLTRA/deSALT emit N-ops of hundreds of kb and the scorer SELECTS them —
    "spans more query" outranks "is physically possible" (planning/684c: 268 of
    400,001 reads with an N-op > 10 kb, max 261,350 bp, 3 running off the contig
    end, in a genome whose longest annotated intron is ~1 kb; the minimap2 arm,
    constrained by -G, produced zero). This is the last point before the record
    reaches `multialigned.bam`, which is what `cdna-analyze` consumes.

    Tags the read `Xi:i:<offending bp>` so downstream viewers (rbrowse) can
    flag it, and counts it. Truncating rather than dropping keeps the
    well-supported 5' portion — see intron_sanity for why.
    """
    if max_intron_bp <= 0 or read.is_unmapped or not read.cigartuples:
        return
    try:
        contig_len = out_bam.get_reference_length(read.reference_name)
    except (KeyError, ValueError):
        contig_len = None
    new_cigar, offending_bp = truncate_impossible_introns(
        read.cigartuples, max_intron_bp, read.reference_start or 0, contig_len,
    )
    if not offending_bp:
        return
    # Xi, NOT Xn — `Xn` is already `n_aligners_agree`, set ~150 lines below at
    # the aligner-selection block, and _enforce_intron_sanity runs AFTER it, so
    # writing Xn here silently overwrote the consensus concordance count on
    # every read the guard fired for. Caught by the rbrowse session, which
    # consumes Xn. The comment on that block spells out the convention this
    # violated: lowercase second letter, and do not collide.
    read.set_tag('Xi', int(offending_bp), value_type='i')
    if stats is not None:
        stats['impossible_intron_truncated'] += 1
    if new_cigar is not None:
        read.cigartuples = new_cigar


def _process_and_write_batch(read_batch, raw_read_batch, genome, annotated_junctions, out_bam, stats, use_chimeric=False, read_num_sidecar=None, chimeric_stats=None, tiebreak='rectify', pool_min_anchor_bp=0, pool_max_intron_len=0):
    """Process a batch of reads and write best alignments to output BAM."""
    _max_reportable_intron = max_reportable_intron_from_env()
    if use_chimeric:
        from .chimeric_consensus import select_best_chimeric, build_chimeric_read

    for i, (read_id, alignments) in enumerate(read_batch):
        _, aligner_reads = raw_read_batch[i]

        if use_chimeric:
            chimeric_result = select_best_chimeric(aligner_reads, genome, annotated_junctions)

            # Measure the SELECTOR (segment wins / fallback rate / loss-reason)
            # independent of whether the BAM write below succeeds — a dropped
            # read still reflects what select_best_chimeric decided.
            if chimeric_stats is not None:
                chimeric_stats.update(chimeric_result)

            # 🔴 The template supplies the output record's CONTIG
            # (`chimeric_consensus.build_chimeric_read` sets
            # ``out.reference_id = template_read.reference_id``) while ``ref_start`` and the
            # CIGAR come from the WINNING aligner. On the true-chimeric path that is safe:
            # ``select_best_chimeric`` refuses to assemble unless every arm is on ONE contig.
            # But its cross-contig FALLBACK (``_fallback_simple_selection``) returns the
            # winner's ``chimeric_ref_start`` with the arms still on DIFFERENT contigs — so a
            # template picked by dict order emits <template's chrom> x <winner's position>, a
            # coordinate that never existed.
            #
            # Measured on a QuantSeq smoke, rectify master @ fd2e2d2, 2026-09-03
            # (planning/862): read D00689:118:C890GANXX:8:2204:16881:55011 was placed by bwa at
            # chrIX:300228 (2S48M) and by bbmap at chrXIII:480619 (1=1X47=1X); the consensus
            # wrote chrIX:480619, which is past the end of chrIX (439,888 nt). It surfaced only
            # because chrIX is SHORTER than the borrowed position — with a longer template
            # contig the fabricated locus is written SILENTLY, on the wrong chromosome.
            #
            # Restrict the template pool to reads on the winning aligner's contig. When every
            # arm is already on one contig (the normal case, and every true-chimeric case) the
            # pool is unchanged and this is a no-op.
            _winner_name = (
                chimeric_result.segment_winners[0][1]
                if chimeric_result.segment_winners else None
            )
            _winner_read = aligner_reads.get(_winner_name) if _winner_name else None
            _winner_chrom = _winner_read.reference_name if _winner_read is not None else None
            if _winner_chrom is not None:
                _template_pool = [
                    r for r in aligner_reads.values()
                    if r.reference_name == _winner_chrom
                ] or list(aligner_reads.values())
            else:
                _template_pool = list(aligner_reads.values())

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
            for r in _template_pool:
                seq = r.query_sequence
                if seq is not None and (expected_len == 0 or len(seq) == expected_len):
                    template = r
                    break
            if template is None:
                # Fallback: accept any read with a sequence even if length mismatches.
                for r in _template_pool:
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

            # Same normalization rule as the non-chimeric path at line 489 —
            # the chimeric template can be any aligner's record, so without
            # this, mapPacBio/uLTRA-templated chimeric reads carry mutated
            # QNAMEs into every downstream QNAME-keyed join.
            out_read.query_name = _normalize_bam_read_name(out_read.query_name or '')
            _restore_sidecar_tags(out_read, read_num_sidecar)
            _enforce_intron_sanity(out_read, out_bam, _max_reportable_intron, stats)
            out_bam.write(out_read)

        else:
            result = select_best_alignment(
                alignments, genome, annotated_junctions, tiebreak=tiebreak,
                pool_min_anchor_bp=pool_min_anchor_bp,
                pool_max_intron_len=pool_max_intron_len,
            )
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
            _credit_tied_aligners(stats, result)
            stats['by_aligner_combo'][frozenset(result.aligners_compared)] += 1

            if result.best_aligner in aligner_reads:
                best_read = aligner_reads[result.best_aligner]

                # Enforce exactly one primary per read: clear secondary (0x100) and
                # supplementary (0x800) bits so the winning record is always primary.
                best_read.flag &= ~0x900

                # Normalize the winning read's QNAME so all downstream tools see a
                # consistent bare form regardless of which aligner won. Without this,
                # the consensus.bam contains a mix of QNAME shapes (uLTRA's
                # `uuid_runid=abc_ch=42`, bbmap's retained comment, mapPacBio's pt:i,
                # etc.) depending on which aligner happened to win each read.
                # Symmetric with the merge key in _iter_name_grouped_bams.
                best_read.query_name = _normalize_bam_read_name(best_read.query_name or '')

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
                _restore_sidecar_tags(best_read, read_num_sidecar)
                _restore_comment_tags_from_siblings(best_read, aligner_reads)
                _enforce_intron_sanity(best_read, out_bam, _max_reportable_intron, stats)
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
    keep_checkpoints: bool = False,
    read_num_sidecar: Optional[str] = None,
    tiebreak: str = 'rectify',
    pool_min_anchor_bp: int = 0,
    pool_max_intron_len: int = 0,
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
        checkpoint_dir: If set, write resumable per-batch checkpoints here
                        (consensus_batch_*.bam + .done sentinels +
                        consensus_checkpoint.json) and finalise via scratch.
        keep_checkpoints: If False (default), delete the dead checkpoint resume
                          state (per-batch BAMs, .done sentinels, JSON) after the
                          final BAM is written successfully, leaving ~0 files.
                          If True, retain the full resumable state for debugging.
                          Only affects the success path; failures always keep
                          checkpoints for resume.

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

    # Ensure BAMs are sorted by the K-way merge key. RN-keyed merges need
    # RN-sorted inputs; name-sort is only correct for the QNAME-keyed path
    # (RN is assigned in FASTQ order, which natural QNAME order does not
    # follow for plain-uuid reads — the 668 DRS group-split defect).
    import time as _time
    _t_total = _time.perf_counter()
    use_rn_key = _all_inputs_have_rn(bam_paths)
    _sorter = _ensure_rn_sorted if use_rn_key else _ensure_name_sorted
    logger.info(
        f"Ensuring BAMs are {'RN' if use_rn_key else 'name'}-sorted...")
    _t_ns = _time.perf_counter()
    sorted_bam_paths = {}
    for aligner, path in bam_paths.items():
        sorted_bam_paths[aligner] = _sorter(path)
    logger.info(f"[TIMING] Merge-key sort: {_time.perf_counter() - _t_ns:.1f}s")

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
        # Consensus selected an alignment carrying a physically impossible
        # intron (planning/684c); truncated at that junction before writing.
        'impossible_intron_truncated': 0,
        'by_aligner': defaultdict(int),
        'by_aligner_combo': defaultdict(int),  # frozenset of available aligners → count
    }

    # Per-region / loss-reason instrumentation for the chimeric path.  Only
    # instantiated when use_chimeric is on; a JSON sidecar + log summary are
    # emitted at the end so a GMAP-segment-win measurement is reproducible.
    chimeric_stats = None
    if use_chimeric:
        from .chimeric_consensus import ChimericStats
        chimeric_stats = ChimericStats()

    sidecar_reader = None
    if read_num_sidecar:
        from rectify.core.chunking.sidecar import ReadNumSidecar
        sidecar_reader = ReadNumSidecar.open(read_num_sidecar)
        logger.info("Read-num sidecar loaded: %s (%d rows)", read_num_sidecar, len(sidecar_reader))

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
        import json as _json_ckpt
        os.makedirs(checkpoint_dir, exist_ok=True)

        # Collect only the contiguous valid prefix from a previous run. Skipping
        # input reads by count is only safe when batch N implies all batches < N.
        _resume_batch_idx = 0
        while True:
            _bam_f = os.path.join(
                checkpoint_dir, f'consensus_batch_{_resume_batch_idx:06d}.bam'
            )
            _done_f = os.path.join(
                checkpoint_dir, f'consensus_batch_{_resume_batch_idx:06d}.done'
            )
            if not os.path.exists(_done_f):
                break
            if not os.path.exists(_bam_f):
                logger.warning(
                    "Ignoring consensus checkpoint from batch %d onward: done marker exists "
                    "without BAM (%s)",
                    _resume_batch_idx, _done_f,
                )
                break
            try:
                _validate_bam_sample(_bam_f)
            except Exception as exc:
                logger.warning(
                    "Ignoring consensus checkpoint from batch %d onward: invalid batch %s (%s)",
                    _resume_batch_idx, _bam_f, exc,
                )
                break
            _batch_bam_paths.append(_bam_f)
            _resume_batch_idx += 1
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

        for read_id, aligner_reads in _iter_name_grouped_bams(
                sorted_bam_paths, use_rn_key=use_rn_key):
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
                    read_num_sidecar=sidecar_reader,
                    chimeric_stats=chimeric_stats,
                    tiebreak=tiebreak,
                    pool_min_anchor_bp=pool_min_anchor_bp,
                    pool_max_intron_len=pool_max_intron_len,
                )
                read_batch = []
                raw_read_batch = []
                n_batches += 1

                if checkpoint_dir:
                    out_bam.close()
                    _sentinel = os.path.join(
                        checkpoint_dir, f'consensus_batch_{_ckpt_batch_num:06d}.done'
                    )
                    _commit_checkpoint_batch(_cur_batch_path, _sentinel)
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
                read_num_sidecar=sidecar_reader,
                chimeric_stats=chimeric_stats,
                tiebreak=tiebreak,
                pool_min_anchor_bp=pool_min_anchor_bp,
                pool_max_intron_len=pool_max_intron_len,
            )
            n_batches += 1
            if checkpoint_dir:
                out_bam.close()
                _sentinel = os.path.join(
                    checkpoint_dir, f'consensus_batch_{_ckpt_batch_num:06d}.done'
                )
                _commit_checkpoint_batch(_cur_batch_path, _sentinel)
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
        _copy2_and_fsync(_sorted_path, output_bam)
        _bai_src = _sorted_path + '.bai'
        if os.path.exists(_bai_src):
            _copy2_and_fsync(_bai_src, output_bam + '.bai')

        # Clean up scratch sort files
        for _p in (_sorted_path, _bai_src):
            try:
                os.unlink(_p)
            except OSError:
                pass
        # ── Success-path checkpoint cleanup ─────────────────────────────────
        # The final coordinate-sorted BAM (+ .bai) has now been written to
        # output_bam above. The remaining checkpoint artefacts — per-batch BAMs
        # (already merged+sorted), their .done sentinels, and
        # consensus_checkpoint.json — are now DEAD resume state: a chunk that
        # has finalised its output never resumes. On a large cluster panel
        # these dominate the file/inode count, so remove them by default.
        # keep_checkpoints=True (--keep-checkpoints) retains the full, coherent
        # resume state (batch BAMs + sentinels + JSON) for debugging.
        #
        # This runs ONLY on the success path. The except-handler above closes
        # out_bam and re-raises WITHOUT touching any checkpoint file, so an
        # interrupted/failed run leaves the batch BAMs + sentinels + JSON intact
        # and a requeued task resumes from the last committed batch.
        if not keep_checkpoints:
            _ckpt_cleanup: List[str] = []
            for _b in _batch_bam_paths:
                _ckpt_cleanup.append(_b)                            # consensus_batch_NNNNNN.bam
                _ckpt_cleanup.append(_b[:-len('.bam')] + '.done')   # matching .done sentinel
            _ckpt_cleanup.append(
                os.path.join(checkpoint_dir, 'consensus_checkpoint.json')
            )
            for _p in _ckpt_cleanup:
                try:
                    os.unlink(_p)
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

    # ── Chimeric instrumentation: human-readable summary + JSON sidecar ──────
    # The sidecar (`<output_bam>.chimeric_stats.json`) carries the per-region
    # segment wins, fallback rate, and contested-segment loss-reason breakdown
    # so a GMAP-segment measurement is reproducible and the annotation-off
    # confound is auditable after the fact.
    if chimeric_stats is not None:
        chimeric_stats.log_summary()
        try:
            import json as _json_cs
            _cs_path = output_bam + '.chimeric_stats.json'
            with open(_cs_path, 'w') as _csf:
                _json_cs.dump(chimeric_stats.summary_dict(), _csf, indent=2)
            logger.info(f"  Chimeric stats sidecar: {_cs_path}")
        except Exception as _cs_err:
            logger.warning("Failed to write chimeric stats sidecar (%s)", _cs_err)

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
    # Buffer exon coordinates per transcript so we can DERIVE introns when the
    # annotation has no explicit `intron` records (the common case: GENCODE /
    # Ensembl exon-GTFs).  Without this, an exon-only GTF loads 0 junctions and
    # silently disables junction-guided scoring (the annotated +8 reward in the
    # chimeric score_segment path goes dead).  The exon path is a strict
    # FALLBACK: if any `intron` record is present (e.g. the yeast SGD GTF/GFF),
    # the intron set is returned unchanged and exon derivation never runs — so
    # intron-bearing annotations stay byte-identical.
    exons_by_tx: Dict[Tuple[str, str, str], List[Tuple[int, int]]] = defaultdict(list)

    def _parse_tx_id(attr: str) -> str:
        # GTF: transcript_id "X";   GFF3: transcript_id=X / Parent=X
        for key in ('transcript_id', 'Parent'):
            idx = attr.find(key)
            if idx == -1:
                continue
            rest = attr[idx + len(key):].lstrip(' =')
            if rest.startswith('"'):
                end = rest.find('"', 1)
                return rest[1:end] if end != -1 else rest[1:]
            # unquoted: up to ';' or whitespace
            for sep in (';', '\t', ' '):
                j = rest.find(sep)
                if j != -1:
                    rest = rest[:j]
            return rest.strip()
        return ''

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
            elif feature_type == 'exon':
                tid = _parse_tx_id(parts[8])
                if tid:
                    chrom = standardize_chrom_name(parts[0])
                    strand = parts[6] if parts[6] in ('+', '-') else '+'
                    exons_by_tx[(chrom, tid, strand)].append(
                        (int(parts[3]), int(parts[4]))
                    )

    # FALLBACK: no explicit intron records → derive from adjacent exons.
    if len(junctions) == 0 and exons_by_tx:
        for (chrom, _tid, strand), exons in exons_by_tx.items():
            exons.sort()
            for i in range(len(exons) - 1):
                # intron spans [prev_exon_end, next_exon_start-1] in 1-based
                # inclusive GTF coords → (0-based start, exclusive end):
                #   start_0 = prev_exon_end   (= prev_exon_end_1based, since
                #             0-based intron start is the base after the exon)
                #   end     = next_exon_start - 1
                intron_start = exons[i][1]          # 0-based: first intronic base
                intron_end = exons[i + 1][0] - 1    # exclusive end
                if intron_end > intron_start:
                    junctions.add((chrom, intron_start, intron_end, strand))
        if junctions:
            logger.info(
                "load_annotated_junctions: 0 explicit intron records; derived "
                "%d introns from %d exon-bearing transcripts in %s",
                len(junctions), len(exons_by_tx), annotation_path,
            )

    if len(junctions) == 0:
        logger.warning(
            "load_annotated_junctions: 0 junctions loaded from %s. "
            "Check that the file exists, is readable, and contains 'intron' "
            "or 'exon' feature records (column 3). Junction-guided scoring "
            "will be disabled.",
            annotation_path,
        )
    else:
        logger.info(f"Loaded {len(junctions)} annotated junctions from {annotation_path}")

    return junctions
