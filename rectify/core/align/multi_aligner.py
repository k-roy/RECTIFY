"""
Multi-aligner alignment pipeline for RECTIFY.

Runs multiple aligners on the same reads and provides utilities for consensus
analysis.  Aligners are grouped into three tiers by read type:

Tier 1 — Long-read (default, nanopore direct RNA-seq / dT-primed cDNA):
  - minimap2:   Fast seed-and-chain baseline with junction annotation support
                (`--junc-bed` improves splice-junction accuracy)
  - mapPacBio:  BBTools long-read aligner with splice-aware mode
  - gapmm2:     minimap2 wrapper with terminal-exon homopolymer refinement

Tier 2 — Long-read, opt-in (`--tier2-aligners deSALT uLTRA`):
  - deSALT:     High-sensitivity splice aligner; vendored binary bundled for
                Linux x86_64 (`rectify/data/bin/linux_x86_64/deSALT`)
  - uLTRA:      Annotation-guided aligner; requires `--annotation` GFF/GTF

Short-read mode (`--short-read`, Illumina/Aviti ≤150 bp):
  - bbmap:      BBTools splice-aware aligner (`intronlen=40` for short introns)
  - bwa:        BWA-MEM aligner (not splice-aware; use bbmap for spliced data)
  When `--short-read` is active the Tier 1 long-read panel is replaced by
  bbmap + bwa, and poly(A)-tail modules are disabled.

Junction annotations are used to IMPROVE alignment quality but scoring
remains BLIND to annotations (novel junctions can still be detected).

Author: Kevin R. Roy
"""

import gzip
import re
import subprocess
import logging
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Mapping
from dataclasses import dataclass, field

from rectify.core.align.qname_validator import validate_post_alignment_qnames

# Pre-compiled regexes for gapmm2 PAF cs-tag parsing (avoid per-call recompilation)
_CS_TOK_RE = re.compile(r':[0-9]+|[*][a-z][a-z]|[+][a-z]+|[-][a-z]+|~[a-z]{2}[0-9]+[a-z]{2}')
_INTRON_LEN_RE = re.compile(r'[0-9]+')

# Reverse-complement table for sequence injection into gapmm2 BAM records.
_RC_TABLE = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def _clean_fastq(reads_path: Path, output_dir: Path) -> Path:
    """Write a deduplicated, tag-stripped FASTQ suitable for deSALT and gapmm2.

    Strips DRS auxiliary tags from read names, skips reads with empty sequences
    (Dorado placeholders), and skips duplicate UUIDs.  Returns the path to the
    temporary FASTQ file; the caller is responsible for deleting it.

    Handles both standard gzip and BGZF-compressed FASTQ.  BGZF (used by
    samtools/dorado output) appends a 28-byte EOF terminator block that
    Python's gzip module misidentifies as a truncation error; we catch that
    EOFError and treat it as a clean end-of-file.
    """
    opener = gzip.open if str(reads_path).endswith('.gz') else open
    tmp_fq = _tempfile.NamedTemporaryFile(suffix='.fastq', dir=output_dir, delete=False)
    tmp_path = Path(tmp_fq.name)
    seen_uuids: set = set()
    n_skipped_empty = 0
    n_skipped_dup = 0
    with opener(reads_path, 'rt') as src, open(tmp_path, 'w') as dst:
        while True:
            try:
                header = src.readline()
            except EOFError:
                break
            if not header:
                break
            try:
                seq  = src.readline()
                plus = src.readline()
                qual = src.readline()
            except EOFError:
                break
            clean_name = header[1:].split()[0]
            if not seq.rstrip():
                n_skipped_empty += 1
                continue
            if clean_name in seen_uuids:
                n_skipped_dup += 1
                continue
            seen_uuids.add(clean_name)
            dst.write(f'@{clean_name}\n{seq}{plus}{qual}')
    logger.info(
        f"  Cleaned FASTQ: {len(seen_uuids)} reads written, "
        f"{n_skipped_empty} empty-seq skipped, {n_skipped_dup} duplicate-UUID skipped"
    )
    return tmp_path


def _load_fastq_sequences(reads_path: str) -> Dict[str, Tuple[str, str]]:
    """Load read sequences and Phred quality strings from a FASTQ file.

    Handles plain and gzip-compressed FASTQ.  Strips DRS auxiliary tags
    (e.g. ``pt:i:25``) that samtools fastq embeds in the read name after a
    tab or space, so the returned keys are bare UUIDs that match the query
    names minimap2 writes into PAF (minimap2 truncates the name at the first
    whitespace).

    Returns:
        {read_name: (sequence, qual_string)}  — forward (5'→3') orientation.
    """
    opener = gzip.open if str(reads_path).endswith('.gz') else open
    seqs: Dict[str, Tuple[str, str]] = {}
    with opener(reads_path, 'rt') as fh:
        while True:
            try:
                header = fh.readline()
            except EOFError:
                break
            if not header:
                break
            try:
                seq = fh.readline().rstrip()
                fh.readline()            # discard '+' line
                qual = fh.readline().rstrip()
            except EOFError:
                break
            # Strip leading '@', then take only the part before any whitespace
            # to match minimap2's PAF query-name truncation behaviour.
            read_name = header[1:].split()[0]
            seqs[read_name] = (seq, qual)
    return seqs



def _load_fastq_rn_map(reads_path: str) -> Dict[str, int]:
    """Return QNAME -> RN map for an input FASTQ.

    When FASTQ headers carry ``RN:i:`` tags (rectify split output), those are
    used verbatim.  When no ``RN:i:`` tags are found, sequential integers
    (0, 1, 2, …) are assigned by read order.  The sequential fallback ensures
    all aligner BAMs receive consistent ``RN:i`` tags even for non-split
    FASTQs (e.g. DRS reads from ``samtools fastq``), which enables the
    RN-keyed K-way consensus merge and bypasses the QNAME-set compatibility
    check that fails when aligners include different read subsets.
    """
    from rectify.core.chunking.sidecar import parse_rn_from_fastq_header

    opener = gzip.open if str(reads_path).endswith('.gz') else open
    tag_map: Dict[str, int] = {}   # built from RN:i: header tags (if present)
    seq_map: Dict[str, int] = {}   # sequential fallback
    n_dup = 0
    seq_counter = 0
    n_records = 0                  # total reads seen (for duplicate-QNAME guard)

    with opener(reads_path, 'rt') as fh:
        while True:
            try:
                header = fh.readline()
            except EOFError:
                break
            if not header:
                break
            try:
                fh.readline()  # seq
                fh.readline()  # +
                fh.readline()  # qual
            except EOFError:
                break
            if not header.startswith('@'):
                continue
            qname, rn = parse_rn_from_fastq_header(header)
            qname = qname[:254]
            n_records += 1
            if rn is not None:
                if qname in tag_map and tag_map[qname] != rn:
                    n_dup += 1
                else:
                    tag_map[qname] = rn
            if qname not in seq_map:
                seq_map[qname] = seq_counter
                seq_counter += 1

    if tag_map:
        if n_dup:
            logger.warning(
                "RN map for %s saw %d duplicate QNAME(s) with conflicting RN tags; "
                "kept the first mapping", reads_path, n_dup,
            )
        return tag_map

    # No RN:i: tags found — assign sequential integers. But if the FASTQ has
    # pervasive DUPLICATE QNAMEs, the per-QNAME sequential map assigns the SAME
    # RN to distinct molecules, and the RN-keyed K-way consensus merge then
    # collapses them to one record — silently dropping the rest (the ~87% cDNA
    # align2 loss; see planning/251, /250a-c). Fail LOUD instead: distinct-name
    # << record-count is never legitimate for a consensus/aligner input.
    n_distinct = len(seq_map)
    if n_records and (n_records - n_distinct) / n_records > 0.01:
        raise RuntimeError(
            f"FASTQ {reads_path}: {n_records} reads but only {n_distinct} distinct "
            f"QNAMEs ({n_records - n_distinct} duplicates, "
            f"{100 * (n_records - n_distinct) / n_records:.1f}%). Assigning "
            "sequential RN per QNAME would give distinct molecules the same RN and "
            "the RN-keyed consensus merge would silently COLLAPSE them (dropping "
            "reads). Root cause is non-globally-unique read names — e.g. per-region "
            "`rectify correct-cdna --region <chr>` output concatenated without a "
            "region prefix. Regenerate with globally-unique names (the --region "
            "path now prefixes cluster names with the region; see planning/251)."
        )
    logger.info(
        "FASTQ %s: no RN:i tags found; assigned sequential RN for %d reads",
        reads_path, n_distinct,
    )
    return seq_map


def _load_fastq_umi_map(reads_path: str) -> Dict[str, str]:
    """Return QNAME -> UMI map from an input FASTQ carrying ``RX:Z:`` header tags.

    This is the short-read UMI counterpart of :func:`_load_fastq_rn_map`, and it
    exists for the same reason: STAR/HISAT2/magicblast/gsnap do NOT copy FASTQ
    comments into the BAM (only ``minimap2 -y`` does), so a UMI written into the
    comment by ``rectify split --umi`` would be silently lost for those aligners.
    We re-read it here, keyed on QNAME, and re-attach it as ``RX:Z`` after
    alignment via :func:`_inject_rn_into_bam`.

    Returns ``{}`` when no ``RX:Z:`` tags are present (non-UMI runs), which makes
    downstream injection a no-op -- so this is safe to call unconditionally in
    every short-read runner. Deliberately additive: it does not touch the RN map,
    whose duplicate-QNAME guards and loud-fail must stay exactly as they are.
    """
    from rectify.core.chunking.sidecar import parse_rx_from_fastq_header

    opener = gzip.open if str(reads_path).endswith('.gz') else open
    umi_map: Dict[str, str] = {}
    n_conflict = 0
    with opener(reads_path, 'rt') as fh:
        while True:
            try:
                header = fh.readline()
            except EOFError:
                break
            if not header:
                break
            try:
                fh.readline()  # seq
                fh.readline()  # +
                fh.readline()  # qual
            except EOFError:
                break
            if not header.startswith('@'):
                continue
            qname, umi = parse_rx_from_fastq_header(header)
            if umi is None:
                continue
            qname = qname[:254]
            existing = umi_map.get(qname)
            if existing is None:
                umi_map[qname] = umi
            elif existing != umi:
                # Both mates should carry the same fragment UMI; a conflict means
                # a malformed header. Keep the first and count it.
                n_conflict += 1
    if n_conflict:
        logger.warning(
            "UMI map for %s saw %d QNAME(s) with conflicting RX tags; kept the first",
            reads_path, n_conflict,
        )
    return umi_map


def _load_fastq_comment_tag_map(reads_path: str) -> Dict[str, list]:
    """Return QNAME -> [(tag, type, value), ...] for FASTQ-comment aux tags.

    🔴 WHY THIS EXISTS. Only minimap2 is invoked with ``-y``; mapPacBio, gapmm2,
    uLTRA and deSALT do not propagate FASTQ comments at all, and several actively
    strip them (``_sanitize_mpb_fastq`` / ``_clean_fastq``). Measured on a real
    ONT-cDNA library (planning/561, 140,000-record consensus):

        winner      share     pl/ro present
        mapPacBio   48.4%     0.00%
        minimap2    40.9%     100.00%
        gapmm2      10.6%     0.00%
        ALL                   40.94%

    So on the DEFAULT 3-aligner panel **59% of records lost both ``pl`` (poly-A
    tail length) and ``ro`` (read orientation)** -- and losing ``ro`` means
    ``correct --ONT-cDNA`` cannot resolve per-read RNA strand, i.e. it silently
    corrupts strand assignment, not merely tail length.

    Restoring here -- inside ``_inject_rn_into_bam``, which already runs on EVERY
    aligner's BAM -- fixes every aligner and every combination at one choke point,
    rather than chasing per-aligner comment flags that may not exist.
    """
    from rectify.core.chunking.sidecar import split_fastq_header

    opener = gzip.open if str(reads_path).endswith('.gz') else open
    out: Dict[str, list] = {}
    with opener(reads_path, 'rt') as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            fh.readline(); fh.readline(); fh.readline()
            if not header.startswith('@'):
                continue
            qname, comment = split_fastq_header(header)
            if not comment:
                continue
            tags = []
            for tok in comment.replace(' ', '\t').split('\t'):
                bits = tok.split(':')
                # SAM aux syntax TAG:TYPE:VALUE. Skip RN -- it is injected
                # separately and authoritatively below.
                if len(bits) < 3 or len(bits[0]) != 2 or bits[0] == 'RN':
                    continue
                tag, typ, val = bits[0], bits[1], ':'.join(bits[2:])
                if typ == 'i':
                    try:
                        val = int(val)
                    except ValueError:
                        continue
                elif typ == 'f':
                    try:
                        val = float(val)
                    except ValueError:
                        continue
                elif typ not in ('Z', 'A'):
                    continue
                tags.append((tag, typ, val))
            if tags:
                out[qname[:254]] = tags
    return out


def _inject_rn_into_bam(bam_path: str, qname_to_rn: Mapping[str, int],
                        qname_to_umi: Optional[Mapping[str, str]] = None,
                        qname_to_tags: Optional[Mapping[str, list]] = None,
                        reads_path: Optional[str] = None) -> int:
    """Stream-rewrite a BAM, adding ``RN:i`` (and optionally ``RX:Z`` UMI) per QNAME.

    ``qname_to_umi``, when supplied, attaches the SAM-standard raw-UMI tag
    ``RX:Z`` to every record whose QNAME is in the map -- so BOTH mates of a pair
    receive the fragment's UMI (the UMI identifies the fragment, not one read).
    Passing an empty/None UMI map preserves the original RN-only behaviour exactly.

    ``qname_to_tags`` restores FASTQ-comment aux tags (``ro``, ``pl``, ``XU``,
    ``XO``, …) on aligners that do not propagate comments. Stage-1 cDNA tags
    (``_CDNA_COMMENT_TAGS``) OVERWRITE same-named aligner-emitted tags — the
    comment is authoritative for them, and uLTRA emits colliding ``XC:Z``/``XA:Z``
    (659). All other existing tags are never overwritten; on minimap2 ``-y``
    records the cDNA-set overwrite is value-identical, so that path is unchanged
    in effect.
    """
    # Build the comment-tag map here rather than at each call site: there are a
    # dozen callers, and requiring each to pass it is precisely the copy-paste
    # pattern that let the ONT-cDNA pre-trim exist in one run-all entry point and
    # not the other. One choke point, every aligner, no per-caller drift.
    if qname_to_tags is None and reads_path:
        try:
            qname_to_tags = _load_fastq_comment_tag_map(str(reads_path))
        except Exception as exc:                      # never fail an alignment over this
            logger.warning("[%s] comment-tag map unavailable (%s); "
                           "FASTQ-comment tags will not be restored", bam_path, exc)
            qname_to_tags = None

    has_umi = bool(qname_to_umi)
    has_tags = bool(qname_to_tags)
    if not qname_to_rn and not has_umi and not has_tags:
        return 0
    import pysam
    from rectify.core.consensus.consensus import (
        _CDNA_COMMENT_TAGS,
        _normalize_bam_read_name,
    )

    src = Path(bam_path)
    tmp = src.with_suffix('.rn_injected.tmp.bam')
    n_tagged = 0
    n_missing = 0
    n_umi_tagged = 0
    n_tag_restored = 0
    with pysam.AlignmentFile(str(src), 'rb') as bam_in, \
         pysam.AlignmentFile(str(tmp), 'wb', header=bam_in.header) as bam_out:
        for read in bam_in:
            qn = _normalize_bam_read_name(read.query_name or '')[:254]
            rn = qname_to_rn.get(qn) if qname_to_rn else None
            if rn is not None:
                read.set_tag('RN', int(rn), value_type='i')
                n_tagged += 1
            else:
                n_missing += 1
            if has_umi:
                umi = qname_to_umi.get(qn)
                if umi:
                    read.set_tag('RX', umi, value_type='Z')
                    n_umi_tagged += 1
            if has_tags:
                for _tag, _typ, _val in qname_to_tags.get(qn, ()):
                    # 659: for the Stage-1 cDNA tag set the FASTQ comment is AUTHORITATIVE and
                    # must OVERWRITE a colliding aligner-emitted tag — uLTRA writes its own
                    # XC:Z:NO_SPLICE / XA:Z: with splice semantics, and the old blanket
                    # don't-clobber guard let them shadow XC:i/XA:i; cdna-analyze's int(XC)
                    # then dropped every uLTRA-won molecule (7.0% of rna15aa_rep1, biased
                    # toward spliced reads). On minimap2 -y records the overwrite is a
                    # value-identical no-op, so that path stays authoritative in effect.
                    # Non-cDNA tags keep the original guard.
                    if read.has_tag(_tag) and _tag not in _CDNA_COMMENT_TAGS:
                        continue
                    try:
                        read.set_tag(_tag, _val, value_type=_typ)
                        n_tag_restored += 1
                    except (TypeError, ValueError):
                        continue
            bam_out.write(read)
    tmp.replace(src)
    if has_tags:
        logger.info("[%s] FASTQ-comment tag restore: %d tags added", bam_path, n_tag_restored)
    if n_missing and qname_to_rn:
        logger.warning(
            "[%s] RN injection: %d/%d records had no QNAME->RN mapping",
            bam_path, n_missing, n_tagged + n_missing,
        )
    if has_umi:
        logger.info("[%s] RX (UMI) injection: %d records tagged", bam_path, n_umi_tagged)
    return n_tagged

# Maximum wall-clock seconds to wait for any single aligner subprocess.
# mapPacBio on 9.7M nanopore reads needs ~3 h; uLTRA/deSALT ~30-60 min.
# Set to 6 h as a safe upper bound before treating a run as hung.
ALIGNER_TIMEOUT = 21600

#: Fallback max intron (bp) when no annotation is available to derive one
#: from. This is the historical S. cerevisiae constant (minimap2 -G 5000);
#: kept only as the no-annotation fallback — see derive_max_intron().
DEFAULT_MAX_INTRON = 5000

#: Bounds on the annotation-derived max intron. The floor keeps a sparse or
#: intron-poor annotation from strangling the aligners; the ceiling matches
#: the long-standing "use 500000 or larger for human" guidance and caps the
#: compute/false-positive bill an outlier annotated intron would otherwise buy.
_DERIVED_MAX_INTRON_BOUNDS = (1000, 500_000)

# uLTRA --reduce_read_ployA threshold high enough that poly-A reduction never
# fires (longest plausible long read « this), so uLTRA never truncates the
# emitted SEQ. See run_ultra for the no-trim rationale.
_ULTRA_DISABLE_POLYA_REDUCE = 10_000_000

from ...utils.junction_bed import get_minimap2_junc_args, generate_junction_bed

from rectify.core.align.mpb_split_reads import split_long_reads, stitch_split_bam, MAX_MPB_READ_LENGTH

logger = logging.getLogger(__name__)

import os as _os
import platform as _platform
import tempfile as _tempfile


def _get_vendored_binary(name: str) -> Optional[str]:
    """Return path to a vendored binary bundled at rectify/data/bin/<os>_<arch>/<name>.

    Returns None if no matching binary exists or it is not executable.
    """
    system = _platform.system().lower()
    machine = _platform.machine().lower()
    if machine == 'arm64':
        machine = 'aarch64'

    try:
        import importlib.resources as _res
        data_pkg = _res.files('rectify').joinpath('data')
        bin_path = data_pkg.joinpath('bin', f'{system}_{machine}', name)
        with _res.as_file(bin_path) as p:
            candidate = str(p)
    except (AttributeError, TypeError, FileNotFoundError):
        here = Path(__file__).parent.parent
        candidate = str(here / 'data' / 'bin' / f'{system}_{machine}' / name)

    if _os.path.isfile(candidate) and _os.access(candidate, _os.X_OK):
        logger.debug(f"Using vendored {name}: {candidate}")
        return candidate

    logger.debug(f"No vendored {name} for {system}/{machine} at {candidate}")
    return None


def _get_vendored_desalt() -> Optional[str]:
    """Return path to the vendored deSALT binary, or None if unavailable."""
    return _get_vendored_binary('deSALT')


@dataclass
class AlignerConfig:
    """Configuration for a single aligner."""
    name: str
    enabled: bool = True
    path: str = ""  # Path to executable (empty = use PATH)
    threads: int = 8
    extra_args: List[str] = field(default_factory=list)


@dataclass
class MultiAlignerConfig:
    """Configuration for multi-aligner pipeline."""
    minimap2: AlignerConfig = field(default_factory=lambda: AlignerConfig(
        name="minimap2",
        enabled=True,
        path="minimap2"
    ))
    map_pacbio: AlignerConfig = field(default_factory=lambda: AlignerConfig(
        name="mapPacBio",
        enabled=True,
        path="mapPacBio.sh"
    ))
    gapmm2: AlignerConfig = field(default_factory=lambda: AlignerConfig(
        name="gapmm2",
        enabled=True,
        path="gapmm2"
    ))
    bbmap: AlignerConfig = field(default_factory=lambda: AlignerConfig(
        name="bbmap",
        enabled=False,  # Opt-in: only appropriate for short reads (Illumina/Aviti)
        path="bbmap.sh"
    ))
    bwa_mem: AlignerConfig = field(default_factory=lambda: AlignerConfig(
        name="bwa",
        enabled=False,  # Opt-in: only appropriate for short reads (Illumina/Aviti)
        path="bwa"
    ))
    # Information-bounded splice-overhang resolver (planning/641): the
    # mapPacBio-role replacement. Native (no external binary); consumes the
    # minimap2 arm BAM, so minimap2 must run in the same invocation.
    overhang_resolver: AlignerConfig = field(default_factory=lambda: AlignerConfig(
        name="overhang_resolver",
        enabled=False,  # Opt-in until the planning/641 T3-T6 acceptance runs pass
        path=""
    ))

    # Junction annotation options
    use_junction_annotation: bool = True
    junc_bonus: int = 9

    # Output options
    keep_sam: bool = False


def check_aligner_available(aligner: str) -> bool:
    """Check if an aligner is available in PATH.

    Returns False for a None exec path. The COMPASS short-read aligners
    (STAR/HISAT2/magicblast/gsnap) resolve their binary inside the wrapper via
    ``_require_binary`` and pass ``exec_path=None`` to the dispatcher's
    availability pre-check; ``shutil.which(None)`` would raise TypeError, so we
    short-circuit (callers already guard, this is belt-and-suspenders)."""
    if aligner is None:
        return False
    return shutil.which(aligner) is not None


def _write_bare_qname_fastq(src: str, dst: str) -> None:
    """Copy a FASTQ to ``dst`` (decompressing gz), rewriting every header to its
    bare QNAME — the first whitespace token.

    Neutralises multi-token FASTQ headers that Magic-BLAST mangles under
    ``-no_query_id_trim``: rectify-split injects ``RN:i:N`` plus the original
    Casava comment, and magicblast spills those tokens into SAM columns 2-3,
    shifting the mandatory fields so samtools rejects the SAM. RN is re-applied
    afterwards from the qname->RN map in ``_finalize_short_read_bam``.
    """
    fin = gzip.open(src, 'rt') if _is_gz(src) else open(src, 'rt')
    with fin, open(dst, 'wt') as fout:
        for i, line in enumerate(fin):
            fout.write(line.split(None, 1)[0] + '\n' if i % 4 == 0 else line)


def run_minimap2(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    annotation_path: Optional[str] = None,
    junc_bonus: int = 9,
    cache_dir: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
    max_intron: int = 5000,
) -> str:
    """Run minimap2 alignment with optional junction annotation.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_bam: Path for output BAM file
        threads: Number of threads
        annotation_path: Optional GFF/GTF for junction annotation
        junc_bonus: Bonus score for annotated junctions
        cache_dir: Directory to cache junction BED
        extra_args: Additional minimap2 arguments

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # Base minimap2 command for direct RNA
    cmd = [
        'minimap2',
        '-ax', 'splice',
        '-uf',  # Stranded (forward)
        '-k14',  # Smaller k-mer for sensitivity
        '-G', str(max_intron),  # Max intron size
        '--splice-flank=no',  # Disable for compatibility
        '--secondary=no',
        '--MD',  # Include MD tag for indel correction and alignment identity
        '-y',   # Copy FASTQ comment fields (SAM-format tags) to BAM aux records
        '-t', str(threads),
    ]

    # Add junction annotation if available
    if annotation_path:
        junc_args = get_minimap2_junc_args(
            annotation_path=annotation_path,
            junc_bonus=junc_bonus,
            cache_dir=cache_dir
        )
        cmd.extend(junc_args)

    # Add extra arguments
    if extra_args:
        cmd.extend(extra_args)

    # Add reference and reads
    cmd.extend([genome_path, reads_path])

    logger.info(f"Running minimap2: {' '.join(cmd[:10])}...")

    # Run minimap2, pipe directly to name-sorted BAM.
    # Name-sort (-n) so consensus selection can stream across aligners without
    # a secondary sort step. Coordinate index is not created — not needed until
    # after consensus selects the best alignment.
    sam_output = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    sort_cmd = [
        'samtools', 'sort',
        '-n',           # name-sort
        '-@', str(threads),
        '-o', str(output_bam)
    ]

    sort_proc = subprocess.Popen(
        sort_cmd,
        stdin=sam_output.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    sam_output.stdout.close()

    # Drain minimap2 stderr in a background thread to prevent the OS pipe buffer
    # (~64 KB) from filling and deadlocking the pipeline: minimap2 blocks on
    # stderr write → cannot write more SAM → sort_proc.communicate() never returns.
    import threading
    mm2_stderr_chunks = []
    def _drain_stderr():
        mm2_stderr_chunks.append(sam_output.stderr.read())
    drain_thread = threading.Thread(target=_drain_stderr, daemon=True)
    drain_thread.start()

    try:
        _, stderr = sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        sam_output.kill()
        raise RuntimeError(f"minimap2/samtools timed out after {ALIGNER_TIMEOUT}s")
    finally:
        drain_thread.join()

    # Check minimap2 exit code (sort failing is often a symptom of minimap2 failing)
    sam_output.wait()
    if sam_output.returncode != 0:
        mm2_err = mm2_stderr_chunks[0].decode(errors='replace') if mm2_stderr_chunks else ''
        raise RuntimeError(f"minimap2 failed (exit {sam_output.returncode}): {mm2_err[-500:]}")

    if sort_proc.returncode != 0:
        raise RuntimeError(f"samtools sort failed: {stderr.decode()}")

    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    logger.info(f"minimap2 complete: {output_bam}")
    validate_post_alignment_qnames(str(output_bam), str(reads_path), 'minimap2')
    qname_to_rn = _load_fastq_rn_map(str(reads_path))
    _inject_rn_into_bam(str(output_bam), qname_to_rn, reads_path=str(reads_path))
    return str(output_bam)


def extract_fastq_chunk(
    input_fastq: str,
    output_fastq: str,
    chunk_idx: int,
    n_chunks: int,
) -> int:
    """Write reads where read_index % n_chunks == chunk_idx to output_fastq.

    Interleaved distribution ensures read-length is evenly spread across chunks
    without a pre-count pass.  Returns number of reads written.
    """
    import gzip as _gzip
    _open_in = _gzip.open if str(input_fastq).endswith('.gz') else open
    n_written = 0
    read_idx = 0
    with _open_in(input_fastq, 'rt') as fin, open(output_fastq, 'w') as fout:
        while True:
            try:
                header = fin.readline()
            except EOFError:
                break
            if not header:
                break
            try:
                seq  = fin.readline()
                plus = fin.readline()
                qual = fin.readline()
            except EOFError:
                break
            if read_idx % n_chunks == chunk_idx:
                fout.write(header if header.endswith('\n') else header + '\n')
                fout.write(seq   if seq.endswith('\n')    else seq   + '\n')
                fout.write('+\n')
                fout.write(qual  if qual.endswith('\n')   else qual  + '\n')
                n_written += 1
            read_idx += 1
    logger.info(
        "FASTQ chunk %d/%d extracted: %d reads → %s",
        chunk_idx, n_chunks, n_written, output_fastq,
    )
    return n_written


def _merge_bams(input_bams: List[str], output_bam: str, threads: int = 1) -> None:
    """Merge a list of name-sorted BAMs into one name-sorted BAM via samtools merge."""
    tmp = str(output_bam) + '.merging_tmp'
    result = subprocess.run(
        ['samtools', 'merge', '-n', '-f', '-@', str(threads), tmp] + list(input_bams),
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        Path(tmp).unlink(missing_ok=True)
        raise RuntimeError(f"samtools merge failed: {result.stderr}")
    Path(tmp).rename(output_bam)


def _apply_calmd_eq(bam_path: Path, genome_path: str, threads: int = 1) -> None:
    """Convert M CIGAR ops to =/X in-place using ``samtools calmd -e``.

    Processes the BAM in place: writes to a ``.calmd_tmp.bam`` sidecar, then
    atomically renames it over the original.  Works on both coordinate-sorted
    and name-sorted BAMs.  Running on a BAM that already has =/X is safe and
    idempotent.

    Args:
        bam_path: Path to the BAM to convert (modified in-place).
        genome_path: Reference FASTA (must have a companion ``.fai`` index).
        threads: samtools thread count (passed as ``-@``).
    """
    tmp = bam_path.with_suffix('.calmd_tmp.bam')
    try:
        with open(tmp, 'wb') as out_fh:
            result = subprocess.run(
                ['samtools', 'calmd', '-e', '-b', f'-@{threads}', str(bam_path), genome_path],
                stdout=out_fh,
                stderr=subprocess.PIPE,
            )
        if result.returncode != 0:
            raise RuntimeError(
                f"samtools calmd failed for {bam_path.name}: "
                f"{result.stderr.decode(errors='replace')[-500:]}"
            )
        tmp.rename(bam_path)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise


# deSALT v1.5.6 segfaults (SIGSEGV, exit 139) or is OOM-killed (exit 137)
# during second-pass "Loop-ProcessReads" when specific pseudo-exon structures
# are inferred from the input.  The crash is deterministic for a given input
# batch — retries never succeed.  When this happens we emit an empty
# name-sorted BAM so the merge step proceeds with the other 4 aligners.
# Upstream bug filed: github.com/ydLiu-HIT/deSALT/issues/49
# 139/137 = shell (128+signal); -11/-9 = Python subprocess (negative signal number)
_DESALT_CRASH_EXITS = frozenset({139, 137, -11, -9})


def _create_empty_name_sorted_bam(
    output_bam: Path,
    genome_path: Optional[str] = None,
) -> None:
    """Write a valid empty name-sorted BAM (header-only) to ``output_bam``.

    When ``genome_path`` is provided, ``@SQ`` lines are synthesised from
    ``<genome_path>.fai`` and prepended after the ``@HD`` line. This is
    required for downstream consumers that open the BAM with pysam's default
    ``check_sq=True`` (e.g. the consensus name-sort step in
    ``rectify align``) — a placeholder with only ``@HD`` raises
    ``ValueError: file has no sequences defined``. Callers that don't supply
    ``genome_path`` get the historical ``@HD``-only behaviour.
    """
    header_lines: List[bytes] = [b'@HD\tVN:1.6\tSO:queryname']
    if genome_path:
        fai_path = Path(str(genome_path) + '.fai')
        if fai_path.exists():
            with open(fai_path) as fh:
                for line in fh:
                    parts = line.rstrip('\n').split('\t', 2)
                    if len(parts) >= 2:
                        header_lines.append(
                            f'@SQ\tSN:{parts[0]}\tLN:{parts[1]}'.encode()
                        )
    payload = b'\n'.join(header_lines) + b'\n'
    result = subprocess.run(
        ['samtools', 'view', '-bS', '-o', str(output_bam)],
        input=payload,
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to create empty BAM at {output_bam}: "
            f"{result.stderr.decode(errors='replace')}"
        )


# Tolerance for BBMap's TranslateColorspaceRead.realign_new AssertionError.
# When a spliced sub-read's genomic extent (greflen) exceeds ALIGN_COLUMNS-80
# (7520), BBMap's padding-shrinkage loops drive extraPadLeft/Right negative
# and the JVM assertion fires.  This drops exactly 1 read per thread crash.
# We accept the partial SAM when the discrepancy is ≤ this many reads.
_MAX_MPB_LOST_READS = 3


def _parse_mpb_lost_reads(stderr: str) -> Optional[int]:
    """Return number of reads lost to BBMap's read-count assertion, or None.

    BBMap's AbstractMapper.printOutputStats fires when reads_out != reads_in.
    The message includes a line like ``316339 != 316340``; we extract the
    delta so callers can decide whether to accept the partial SAM.
    """
    if 'number of reads out does not add up' not in stderr:
        return None
    m = re.search(r'(\d[\d+]*) = (\d+) != (\d+)', stderr)
    if not m:
        return None
    reads_out, reads_in = int(m.group(2)), int(m.group(3))
    if reads_in > reads_out:
        return reads_in - reads_out
    return None


def _sanitize_mpb_fastq(src_path: str, dst_path: str) -> None:
    """Rewrite a FASTQ for mapPacBio with SAM-spec-compliant QNAMEs.

    Strips everything after the first whitespace (space OR tab) and truncates
    to 254 chars. Runs unconditionally per-record — a per-file gate misses
    mixed-comment files whose first read is bare-UUID while later reads carry
    comments, and missing the tab case leaks minimap2 `-y` aux tags into the
    BAM QNAME as `id_XA:Z:foo` after samtools sort.

    Accepts gzipped or plain FASTQ for ``src_path``.
    """
    opener = gzip.open if str(src_path).endswith('.gz') else open
    with opener(src_path, 'rt') as fin, open(dst_path, 'w') as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            sep = fin.readline()
            qual = fin.readline()
            qname = header[1:].rstrip('\n')
            cut = len(qname)
            for ws in (' ', '\t'):
                i = qname.find(ws)
                if 0 <= i < cut:
                    cut = i
            fout.write('@' + qname[:cut][:254] + '\n')
            fout.write(seq)
            fout.write(sep)
            fout.write(qual)


def _bbtools_xmx_gb(genome_path: str) -> int:
    """Java heap (GB) for BBTools (mapPacBio / bbmap), sized to the genome.

    BBTools' wrapper scripts set ``-Xms`` equal to ``-Xmx``, so an oversized
    ``-Xmx`` force-commits that much at JVM startup and OOM-crashes the JVM under
    a cgroup memory limit. This was observed crashing mapPacBio on yeast (the old
    hard-coded ``-Xmx32g`` while gmap/uLTRA were also resident on a 48 GB task).
    Yeast-scale references need only a few GB; mammalian BBMap indexes need
    ~24-32 GB. Size from the (possibly gzipped) FASTA on disk: >80 MB ~ mammalian.
    Override with the RECTIFY_BBTOOLS_XMX_GB env var if needed.
    """
    import os as _os
    _env = _os.environ.get('RECTIFY_BBTOOLS_XMX_GB')
    if _env and _env.isdigit():
        return int(_env)
    try:
        gsize = Path(genome_path).stat().st_size
    except OSError:
        gsize = 0
    return 32 if gsize > 80_000_000 else 8


def run_map_pacbio(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    extra_args: Optional[List[str]] = None,
    chunk_idx: Optional[int] = None,
    n_chunks: Optional[int] = None,
    max_intron: int = 5000,
) -> str:
    """Run BBTools mapPacBio alignment.

    mapPacBio is splice-aware for long reads. The mapPacBio.sh script default
    of fastareadlen=6000 causes AssertionErrors on reads >6019 bp; we patch
    the script default to 100000 at install time (see multi_aligner.py).
    intronlen controls the maximum intron length recognized as an N-op in CIGAR;
    use >=100k for mammalian RNA-seq (human introns up to ~2.5 Mbp).

    Chunked mode (chunk_idx / n_chunks):
        When chunk_idx is provided, only the reads where
        read_index % n_chunks == chunk_idx are processed.  The output is written
        to {output_bam_stem}.chunk_{chunk_idx}_of_{n_chunks}.bam instead of the
        final output_bam path.

        When chunk_idx is None but n_chunks > 1, the function looks for all N
        chunk BAMs and merges them into output_bam.  If any chunk BAM is missing
        it falls back to a full single-pass alignment.

    Args:
        reads_path: Path to FASTQ file (or FASTQ.GZ)
        genome_path: Path to genome FASTA
        output_bam: Path for output BAM file (or chunk BAM when chunk_idx set)
        threads: Number of threads
        extra_args: Additional mapPacBio arguments
        chunk_idx: 0-based index of this chunk (None = all reads)
        n_chunks: Total number of chunks (None or 1 = no chunking)

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # ── Chunk merge mode: all chunk BAMs exist → merge and return ──────────
    if chunk_idx is None and n_chunks and n_chunks > 1:
        chunk_bams = [
            output_bam.parent / f"{output_bam.stem}.chunk_{k}_of_{n_chunks}.bam"
            for k in range(n_chunks)
        ]
        if all(p.exists() for p in chunk_bams):
            logger.info(
                "All %d mapPacBio chunk BAMs found — merging into %s",
                n_chunks, output_bam,
            )
            _merge_bams([str(p) for p in chunk_bams], str(output_bam), threads=threads)
            _apply_calmd_eq(output_bam, genome_path, threads=threads)
            logger.info("mapPacBio merge complete: %s", output_bam)
            return str(output_bam)
        else:
            # HARD ERROR, not a silent fallback (planning/641 §3, planning/632):
            # "--mapPacBio-chunks N" without "--mapPacBio-chunk-idx" is the
            # MERGE verb. Falling back to a full single-pass alignment here
            # turned an intended parallel run into a ~41 h monolith once
            # already. Fail loudly with the two legitimate next steps.
            missing = [p.name for p in chunk_bams if not p.exists()]
            raise FileNotFoundError(
                f"mapPacBio chunk merge: {len(missing)}/{n_chunks} chunk BAMs "
                f"missing ({', '.join(missing)}). --mapPacBio-chunks N without "
                f"--mapPacBio-chunk-idx MERGES existing chunks — it does not "
                f"parallelise. Either run each chunk first "
                f"(--mapPacBio-chunks {n_chunks} --mapPacBio-chunk-idx K for "
                f"K=0..{n_chunks - 1}), or drop --mapPacBio-chunks for a full "
                f"single-pass alignment."
            )

    # ── Chunk extraction mode: write only reads for this chunk ─────────────
    actual_reads_path = reads_path
    chunk_tmp_fq = None
    if chunk_idx is not None and n_chunks and n_chunks > 1:
        # Redirect output to the chunk BAM path
        output_bam = output_bam.parent / f"{output_bam.stem}.chunk_{chunk_idx}_of_{n_chunks}.bam"
        chunk_tmp_fq = str(output_bam.with_suffix('.chunk_input.fastq'))
        extract_fastq_chunk(reads_path, chunk_tmp_fq, chunk_idx, n_chunks)
        actual_reads_path = chunk_tmp_fq

    qname_to_rn = _load_fastq_rn_map(str(actual_reads_path))

    # Check if mapPacBio is available
    map_pacbio_path = shutil.which('mapPacBio.sh')
    if not map_pacbio_path:
        if chunk_tmp_fq:
            Path(chunk_tmp_fq).unlink(missing_ok=True)
        raise FileNotFoundError("mapPacBio.sh not found in PATH")

    # mapPacBio outputs SAM to a temp file alongside the output BAM
    sam_path = output_bam.with_suffix('.sam')

    # Store BBTools index alongside the genome so all jobs share it
    mpb_index_dir = Path(genome_path).parent / 'bbmap_index'
    mpb_index_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        map_pacbio_path,
        f'ref={genome_path}',
        f'in={actual_reads_path}',
        f'out={sam_path}',
        f'threads={threads}',
        f'path={mpb_index_dir}',
        'fastareadlen=100000',       # Belt-and-suspenders: also patched in mapPacBio.sh default
        # BBMap `intronlen` is the MIN deletion length relabeled D->N (Cufflinks
        # convention), NOT a max — keep it small or real introns stay D-ops.
        'intronlen=40',
        # `maxindel` caps the largest gap BBMap will SEARCH for; the default 16000
        # cannot span mammalian introns (read gets soft-clipped instead). Decoupled
        # from `max_intron` (which feeds minimap2 -G / gapmm2): BBMap needs >=200k
        # for human RNA. See docs/aligners/mapPacBio.md.
        f'maxindel={max(200000, max_intron)}',
        'minratio=0.4',
        # Heap sized to the genome (BBTools sets -Xms==-Xmx; an oversized value
        # force-commits at JVM start and OOM-crashes under a cgroup). See helper.
        f'-Xmx{_bbtools_xmx_gb(genome_path)}g',
    ]

    if extra_args:
        cmd.extend(extra_args)

    # ── Split long reads (>6 kb) for mapPacBio ──
    mpb_split_fq = Path(sam_path).with_suffix('.split.fastq')
    mpb_chunk_map, mpb_n_split = split_long_reads(
        actual_reads_path, str(mpb_split_fq),
        max_length=MAX_MPB_READ_LENGTH,
    )
    if mpb_n_split > 0:
        cmd = [f'in={mpb_split_fq}' if c.startswith('in=') else c for c in cmd]

    # ── Sanitise QNAME for SAM spec compliance ──
    # ONT Dorado FASTQs embed run metadata in the read description:
    #   @<uuid> runid=... ch=... start_time=... flow_cell_id=...
    # minimap2 -y / cDNA-pipeline FASTQs use tab-separated SAM-format aux tags:
    #   @<uuid>\tXA:Z:foo\tXC:i:42
    # mapPacBio.sh copies the full header verbatim into the SAM QNAME field,
    # which violates the SAM spec 254-char limit and causes ``samtools view``
    # to exit 1 with "query name too long". The sanitiser runs per-record
    # (see ``_sanitize_mpb_fastq``) and strips both space and tab comments.
    _mpb_in = next(
        (c.split('=', 1)[1] for c in cmd if c.startswith('in=')),
        str(actual_reads_path),
    )
    _mpb_san_fq = Path(sam_path).with_suffix('.mpb_san.fastq')
    _sanitize_mpb_fastq(_mpb_in, str(_mpb_san_fq))
    cmd = [f'in={_mpb_san_fq}' if c.startswith('in=') else c for c in cmd]

    logger.info(
        "Running mapPacBio: %s",
        ' '.join(cmd[:5]) + (f' [chunk {chunk_idx}/{n_chunks}]' if chunk_idx is not None else ''),
    )

    # The sanitised FASTQ is only needed while mapPacBio runs. Unlink it in a
    # finally: it is uncompressed and was measured at 55.1% of a whole panel
    # output directory when leaked (planning/641 §3 / planning/586).
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)
    finally:
        _mpb_san_fq.unlink(missing_ok=True)
        if chunk_tmp_fq:
            Path(chunk_tmp_fq).unlink(missing_ok=True)
    if result.returncode != 0:
        lost = _parse_mpb_lost_reads(result.stderr)
        if (lost is not None and 0 < lost <= _MAX_MPB_LOST_READS
                and sam_path.exists() and sam_path.stat().st_size > 0):
            # BBMap's TranslateColorspaceRead.realign_new assertion: a spliced
            # sub-read with greflen >= ALIGN_COLUMNS-80 (7520) drives padding
            # negative and crashes the thread, losing exactly 1 read.
            # The partial SAM is valid for all other reads — accept it.
            logger.warning(
                "mapPacBio: accepting partial SAM — %d read(s) lost to "
                "realign_new AssertionError (BBMap spliced greflen overflow; "
                "upstream bug report: github.com/bbushnell/BBTools/issues/19)",
                lost,
            )
        else:
            mpb_split_fq.unlink(missing_ok=True)
            raise RuntimeError(f"mapPacBio failed: {result.stderr}")

    # Convert SAM to name-sorted BAM. Name-sort (-n) so consensus selection
    # can stream across aligners without a secondary sort step.
    #
    # stderr=subprocess.DEVNULL for view_proc: using PIPE without a reader
    # causes a deadlock when samtools view writes enough warnings to fill the
    # OS pipe buffer (~64 KB) — view_proc blocks on stderr write, stops
    # forwarding data to sort_proc, and sort_proc's communicate() hangs.
    view_proc = subprocess.Popen(
        ['samtools', 'view', '-bS', str(sam_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-n', '-@', str(threads), '-o', str(output_bam)],
        stdin=view_proc.stdout,
        stderr=subprocess.PIPE,
    )
    view_proc.stdout.close()
    try:
        _sort_out, sort_stderr = sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
        if sort_proc.returncode != 0:
            raise RuntimeError(
                f"samtools sort failed: {sort_stderr.decode(errors='replace') if sort_stderr else ''}"
            )
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        view_proc.kill()
        sort_proc.communicate()
        view_proc.communicate()
        raise RuntimeError("samtools sort (mapPacBio) timed out")
    sam_path.unlink(missing_ok=True)
    view_proc.wait()
    if view_proc.returncode not in (0, -9):  # -9 = SIGKILL on normal exit path
        raise RuntimeError(f"samtools view (mapPacBio) failed with exit code {view_proc.returncode}")

    # ── Strip mapPacBio _pt:i:N suffix from QNAMEs (Read ID Purity Policy) ──
    # mapPacBio appends its poly-A call to the QNAME (e.g. "UUID pt:i:25" in
    # native SAM, "UUID_pt:i:25" after samtools sort converts spaces to
    # underscores per SAM spec).  Strip it here so all downstream consumers
    # (bam_processor, corrected_reads.tsv, corrected_consensus.bam) see only
    # bare UUIDs — consistent with every other aligner and the parquet metadata.
    import pysam
    _pt_fixed = output_bam.with_suffix('.ptfix.bam')
    _n_stripped = 0
    with pysam.AlignmentFile(str(output_bam), 'rb') as _src, \
         pysam.AlignmentFile(str(_pt_fixed), 'wb', header=_src.header) as _dst:
        for _read in _src:
            _qn = _read.query_name or ''
            if '_pt:i:' in _qn:
                _read.query_name = _qn.split('_pt:i:')[0]
                _n_stripped += 1
            _dst.write(_read)
    output_bam.unlink()
    _pt_fixed.rename(output_bam)
    if _n_stripped:
        logger.info("mapPacBio: stripped _pt:i:N suffix from %d QNAMEs", _n_stripped)

    # ── Stitch split mapPacBio chunks back into full-read alignments ──
    if mpb_n_split > 0:
        mpb_pre_stitch = output_bam.with_suffix('.pre_stitch.bam')
        output_bam.rename(mpb_pre_stitch)
        stitch_split_bam(
            str(mpb_pre_stitch), str(output_bam),
            mpb_chunk_map, threads=threads,
        )
        mpb_pre_stitch.unlink(missing_ok=True)
    mpb_split_fq.unlink(missing_ok=True)

    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    logger.info("mapPacBio complete: %s", output_bam)
    validate_post_alignment_qnames(str(output_bam), str(actual_reads_path), 'mapPacBio')
    _inject_rn_into_bam(str(output_bam), qname_to_rn, reads_path=str(actual_reads_path))
    return str(output_bam)


def run_bbmap(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    extra_args: Optional[List[str]] = None,
    bbmap_path: Optional[str] = None,
    reads2_path: Optional[str] = None,
) -> str:
    """Run vanilla BBMap for short-read splice-aware alignment (Illumina/Aviti).

    Uses intronlen=40 so any reference gap ≥40 bp is encoded as an N-op
    (intron skip) rather than a D-op (deletion).  This is critical for
    rectify's junction_refiner, which scans for N-ops.

    Unlike mapPacBio, this invokes align2.BBMap (not align2.BBMapPacBio),
    which is tuned for short reads and does not apply the PacBio-specific
    long-read chunking logic.  The bbmap_index directory is shared with
    mapPacBio (same BBTools suite, same reference build).

    Args:
        reads_path: Path to FASTQ (or FASTQ.GZ)
        genome_path: Path to genome FASTA
        output_bam: Path for output BAM
        threads: Number of threads
        extra_args: Additional bbmap.sh arguments

    Returns:
        Path to output BAM (name-sorted)
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    qname_to_rn = _load_fastq_rn_map(str(reads_path))
    qname_to_umi = _load_fastq_umi_map(str(reads_path))  # RX:Z UMI, no-op if absent

    # Resolve bbmap.sh: caller-provided path wins (absolute or PATH-resolvable);
    # otherwise fall back to PATH lookup of bare "bbmap.sh".
    if bbmap_path and bbmap_path != 'bbmap.sh':
        bbmap_resolved = bbmap_path if Path(bbmap_path).is_file() else shutil.which(bbmap_path)
    else:
        bbmap_resolved = shutil.which('bbmap.sh')
    if not bbmap_resolved:
        raise FileNotFoundError(
            f"bbmap.sh not found (tried path={bbmap_path!r}, also not on PATH). "
            "Pass --bbmap-path /full/path/to/bbmap.sh, or activate a conda env "
            "with bbtools installed."
        )
    bbmap_path = bbmap_resolved

    sam_path = output_bam.with_suffix('.sam')
    bbmap_index_dir = Path(genome_path).parent / 'bbmap_index'
    bbmap_index_dir.mkdir(parents=True, exist_ok=True)

    # Paired-end (COMPASS short-read): in1=/in2= with mate-2 supplied. Single-end
    # keeps the legacy in= form.
    if reads2_path:
        _in_args = [f'in1={reads_path}', f'in2={reads2_path}']
    else:
        _in_args = [f'in={reads_path}']
    # Intron sizing: COMPASS's human panel uses maxindel=200000 pairlen=200000
    # (MAX_INTRON_LENGTH=200000). The single-end yeast path keeps the smaller
    # cap. pairlen (bbmap default ~32kb) must match maxindel for paired mode or
    # mates spanning long introns are not flagged proper-pair / are placed
    # differently than the rest of the COMPASS panel.
    if reads2_path:
        _maxindel = 200000
        _bbmap_intron_args = [f'maxindel={_maxindel}', f'pairlen={_maxindel}']
    else:
        _bbmap_intron_args = ['maxindel=100000']
    cmd = [
        bbmap_path,
        f'ref={genome_path}',
        *_in_args,
        f'out={sam_path}',
        f'threads={threads}',
        f'path={bbmap_index_dir}',
        'intronlen=40',      # Gaps ≥40 bp → N-op (intron skip) in CIGAR
        *_bbmap_intron_args,  # maxindel (+ pairlen when paired)
        'minratio=0.56',     # Default short-read sensitivity
        'ambiguous=best',
        'trd=t',             # Trim read description: keep only the first
                             # whitespace token of the FASTQ header in QNAME.
                             # Without this, BBmap retains the full header
                             # (e.g. 'SRR.123 123 length=76') while BWA truncates
                             # to the bare accession — consensus K-way merge
                             # then fails to join cross-aligner reads.
        f'-Xmx{_bbtools_xmx_gb(genome_path)}g',   # genome-sized heap (see helper)
    ]

    if extra_args:
        cmd.extend(extra_args)

    logger.info("Running bbmap: %s", ' '.join(cmd[:5]) + '...')
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)
    if result.returncode != 0:
        raise RuntimeError(f"bbmap failed: {result.stderr[-1000:]}")

    # SAM → name-sorted BAM (same pattern as mapPacBio)
    view_proc = subprocess.Popen(
        ['samtools', 'view', '-bS', str(sam_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-n', '-@', str(threads), '-o', str(output_bam)],
        stdin=view_proc.stdout,
        stderr=subprocess.PIPE,
    )
    view_proc.stdout.close()
    try:
        _, sort_stderr = sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
        if sort_proc.returncode != 0:
            raise RuntimeError(
                f"samtools sort (bbmap) failed: {sort_stderr.decode(errors='replace')}"
            )
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        view_proc.kill()
        sort_proc.communicate()
        view_proc.communicate()
        raise RuntimeError("samtools sort (bbmap) timed out")
    finally:
        sam_path.unlink(missing_ok=True)

    view_proc.wait()
    if view_proc.returncode not in (0, -9):
        raise RuntimeError(
            f"samtools view (bbmap) failed with exit code {view_proc.returncode}"
        )

    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    logger.info("bbmap complete: %s", output_bam)
    validate_post_alignment_qnames(str(output_bam), str(reads_path), 'bbmap')
    _inject_rn_into_bam(str(output_bam), qname_to_rn, qname_to_umi, reads_path=str(reads_path))
    return str(output_bam)


def run_bwa_mem(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    extra_args: Optional[List[str]] = None,
    bwa_path: Optional[str] = None,
) -> str:
    """Run BWA-MEM for short-read alignment (Illumina/Aviti).

    BWA-MEM is not splice-aware: intron-spanning reads are split into
    supplementary alignments rather than getting a single N-op CIGAR.
    For 3'-end sequencing data (QuantSeq, MACE-seq) where reads are
    predominantly in UTRs, this is acceptable — the multi-aligner
    consensus will prefer bbmap over bwa for the minority of reads that
    span introns.

    A BWA index is built automatically the first time, alongside the
    genome FASTA.

    Args:
        reads_path: Path to FASTQ (or FASTQ.GZ)
        genome_path: Path to genome FASTA
        output_bam: Path for output BAM
        threads: Number of threads
        extra_args: Additional bwa mem arguments

    Returns:
        Path to output BAM (name-sorted)
    """
    import threading as _threading

    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    qname_to_rn = _load_fastq_rn_map(str(reads_path))
    qname_to_umi = _load_fastq_umi_map(str(reads_path))  # RX:Z UMI, no-op if absent

    # Resolve bwa: caller-provided path wins; otherwise PATH lookup.
    if bwa_path and bwa_path != 'bwa':
        bwa_resolved = bwa_path if Path(bwa_path).is_file() else shutil.which(bwa_path)
    else:
        bwa_resolved = shutil.which('bwa')
    if not bwa_resolved:
        raise FileNotFoundError(
            f"bwa not found (tried path={bwa_path!r}, also not on PATH). "
            "Pass --bwa-path /full/path/to/bwa, or `module load bwa`."
        )
    bwa_path = bwa_resolved

    # Build BWA index if not present (index files sit next to the FASTA)
    bwt_path = Path(str(genome_path) + '.bwt')
    if not bwt_path.exists():
        logger.info("Building BWA index: %s (one-time, ~1 min for yeast)", genome_path)
        idx_result = subprocess.run(
            ['bwa', 'index', str(genome_path)],
            capture_output=True, text=True,
        )
        if idx_result.returncode != 0:
            raise RuntimeError(f"bwa index failed: {idx_result.stderr[-500:]}")
        logger.info("BWA index built: %s", genome_path)

    cmd = [
        bwa_path, 'mem',
        '-t', str(threads),
        '-M',   # Mark split-hit secondary alignments (samtools-compatible)
    ]
    if extra_args:
        cmd.extend(extra_args)
    cmd.extend([str(genome_path), str(reads_path)])

    logger.info("Running bwa mem: %s", ' '.join(cmd[:6]) + '...')

    bwa_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-n', '-@', str(threads), '-o', str(output_bam)],
        stdin=bwa_proc.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    bwa_proc.stdout.close()

    bwa_stderr_chunks: list = []
    def _drain_bwa_stderr():
        bwa_stderr_chunks.append(bwa_proc.stderr.read())
    drain_thread = _threading.Thread(target=_drain_bwa_stderr, daemon=True)
    drain_thread.start()

    try:
        _, sort_stderr = sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        bwa_proc.kill()
        sort_proc.communicate()
        raise RuntimeError("bwa mem / samtools sort timed out")
    finally:
        drain_thread.join()

    bwa_proc.wait()
    if bwa_proc.returncode != 0:
        bwa_err = bwa_stderr_chunks[0].decode(errors='replace') if bwa_stderr_chunks else ''
        raise RuntimeError(
            f"bwa mem failed (exit {bwa_proc.returncode}): {bwa_err[-500:]}"
        )
    if sort_proc.returncode != 0:
        raise RuntimeError(
            f"samtools sort (bwa mem) failed: {sort_stderr.decode(errors='replace')}"
        )

    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    logger.info("bwa mem complete: %s", output_bam)
    validate_post_alignment_qnames(str(output_bam), str(reads_path), 'bwa')
    _inject_rn_into_bam(str(output_bam), qname_to_rn, qname_to_umi, reads_path=str(reads_path))
    return str(output_bam)


# ══════════════════════════════════════════════════════════════════════════
# COMPASS short-read aligner panel (paired-end)
# ══════════════════════════════════════════════════════════════════════════
# STAR_default / STAR_noncanonical, HISAT2_default / HISAT2_noncanonical,
# Magic-BLAST, GSNAP. Invocations + flags are ported verbatim from the COMPASS
# process_reads_and_align.sh aligner block. Three deliberate differences:
#
#   1. gzipped input. COMPASS fed UNCOMPRESSED renumbered reads (@N_R1); RECTIFY
#      feeds gzipped chunk FASTQs (with RN:i comments), so gz-aware flags are
#      added per tool — STAR `--readFilesCommand zcat`, GSNAP `--gunzip`; HISAT2,
#      Magic-BLAST and bbmap read .gz natively.
#   2. No samfixcigar. COMPASS reformatted CIGARs to =/X for its OWN junction
#      comparison code. RECTIFY's consensus scoring reads the genome directly and
#      handles plain `M` ops, so that step is dropped; `_apply_calmd_eq`
#      normalizes MD/NM (and =/X) natively, same as the bbmap/bwa wrappers.
#   3. Read identity. COMPASS renumbered headers destructively; RECTIFY keeps the
#      original QNAME and re-derives the shared RN via `_inject_rn_into_bam`
#      (both mates share one bare QNAME + RN — see split_fastq_paired).
#
# Index layout mirrors COMPASS exactly so the prebuilt indices under
# genome_references/ are reused with no extra arguments:
#   REFERENCE_DIR   = <genome FASTA parent>
#   GENOME_VERSION  = basename(FASTA) minus the fasta extension
#   STAR_GENOME_DIR = REFERENCE_DIR/STAR_annotated_<read_length>_bp_SJDB_index
#   HISAT2_INDEX    = REFERENCE_DIR/HISAT2_annotated_index/<GENOME_VERSION>
#   SPLICE_SITES    = REFERENCE_DIR/HISAT2_annotated_index/<GV>_splice_sites.txt
#   BLAST_INDEX     = REFERENCE_DIR/BLAST/<GENOME_VERSION>
#   GSNAP_GENOME_DIR= REFERENCE_DIR/GSNAP  (-d <GENOME_VERSION>)

from collections import namedtuple as _namedtuple

_CompassIndexPaths = _namedtuple(
    '_CompassIndexPaths',
    'ref_dir genome_version star_dir hisat2_index splice_sites blast_index gsnap_dir',
)

_FASTA_EXTS = ('.fasta.gz', '.fa.gz', '.fna.gz', '.fasta', '.fa', '.fna')


def _genome_version(genome_path) -> str:
    """basename(FASTA) with the fasta extension removed (COMPASS GENOME_VERSION)."""
    name = Path(genome_path).name
    for ext in _FASTA_EXTS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return Path(genome_path).stem


def _compass_index_paths(genome_path, read_length: int = 150) -> _CompassIndexPaths:
    """Derive the COMPASS index paths from the genome FASTA + read length."""
    ref_dir = Path(genome_path).resolve().parent
    gv = _genome_version(genome_path)
    hisat2_dir = ref_dir / 'HISAT2_annotated_index'
    return _CompassIndexPaths(
        ref_dir=ref_dir,
        genome_version=gv,
        star_dir=ref_dir / f'STAR_annotated_{read_length}_bp_SJDB_index',
        hisat2_index=hisat2_dir / gv,
        splice_sites=hisat2_dir / f'{gv}_splice_sites.txt',
        blast_index=ref_dir / 'BLAST' / gv,
        gsnap_dir=ref_dir / 'GSNAP',
    )


def _is_gz(p) -> bool:
    return str(p).endswith('.gz')


# ── Pure command builders (unit-testable without the binaries) ───────────────

def _build_star_cmd(reads_r1, reads_r2, star_genome_dir, out_prefix, threads,
                    read_length: int = 150, noncanonical: bool = False,
                    max_intron: int = 0) -> List[str]:
    # reads_r2=None -> single-end (TruSeq-style SE COMPASS subset)
    reads_in = [str(reads_r1)] + ([str(reads_r2)] if reads_r2 else [])
    cmd = [
        'STAR',
        '--runThreadN', str(threads),
        '--genomeDir', str(star_genome_dir),
        '--sjdbOverhang', str(read_length - 1),
        '--readFilesIn', *reads_in,
        '--outFileNamePrefix', str(out_prefix),
        '--alignEndsType', 'EndToEnd',
        '--outSAMattributes', 'NH', 'HI', 'NM', 'MD', 'AS', 'nM', 'jM', 'jI', 'XS',
    ]
    # planning/833 G-8: without this STAR uses its own winBinNbits-derived
    # default (~589 kb), so it emits junctions rectify itself calls impossible
    # and they enter the consensus candidate pool. `--alignMatesGapMax` is the
    # paired counterpart -- STAR treats an uncapped mate gap the same way.
    if max_intron and max_intron > 0:
        cmd += ['--alignIntronMax', str(max_intron)]
        if reads_r2:
            cmd += ['--alignMatesGapMax', str(max_intron)]
    if _is_gz(reads_r1):
        cmd += ['--readFilesCommand', 'zcat']
    if noncanonical:
        cmd += ['--scoreGapNoncan', '0']
    return cmd


def _build_hisat2_cmd(reads_r1, reads_r2, hisat2_index, splice_sites, out_sam,
                     threads, min_intron, max_intron, novel_sj_out, summary_out,
                     noncanonical: bool = False) -> List[str]:
    cmd = [
        'hisat2',
        '--known-splicesite-infile', str(splice_sites),
        '--no-softclip',
        '--threads', str(threads),
        '--time', '--reorder',
        '-x', str(hisat2_index),
    ]
    # reads_r2=None -> single-end (-U); else paired (-1/-2)
    if reads_r2:
        cmd += ['-1', str(reads_r1), '-2', str(reads_r2)]
    else:
        cmd += ['-U', str(reads_r1)]
    cmd += [
        '-S', str(out_sam),
        '--min-intronlen', str(min_intron),
        '--max-intronlen', str(max_intron),
    ]
    if noncanonical:
        cmd += ['--pen-noncansplice', '0']
    cmd += [
        '--rna-strandness', 'RF',
        '--novel-splicesite-outfile', str(novel_sj_out),
        '--summary-file', str(summary_out),
        '--new-summary',
    ]
    return cmd


def _build_magicblast_cmd(reads_r1, reads_r2, blast_index, out_sam,
                         threads: int = 12, max_intron: int = 0) -> List[str]:
    cmd = [
        'magicblast',
        '-query', str(reads_r1),
        '-query_mate', str(reads_r2),
        '-db', str(blast_index),
        '-md_tag', '-fr', '-no_query_id_trim',
        '-infmt', 'fastq',
        '-num_threads', str(threads),
        '-max_db_word_count', '10',
    ]
    # planning/833 G-8: Magic-BLAST's own default is 500,000 bp.
    if max_intron and max_intron > 0:
        cmd += ['-max_intron_length', str(max_intron)]
    cmd += ['-out', str(out_sam)]
    return cmd


_GSNAP_AMBIG_SUPPORT = None


def _gsnap_supports_ambig_noclip() -> bool:
    """Whether the gsnap on PATH accepts ``--ambig-splice-noclip``.

    GSNAP 2024-11-20 removed this flag (present through the COMPASS-pinned
    2021-05-27). Probed once via ``gsnap --help`` and cached; on any probe error
    we assume supported (the legacy default) so behaviour is unchanged when the
    binary is the pinned one.
    """
    global _GSNAP_AMBIG_SUPPORT
    if _GSNAP_AMBIG_SUPPORT is None:
        try:
            r = subprocess.run(['gsnap', '--help'], capture_output=True,
                               text=True, timeout=30)
            _GSNAP_AMBIG_SUPPORT = '--ambig-splice-noclip' in (r.stdout + r.stderr)
        except Exception:
            _GSNAP_AMBIG_SUPPORT = True
    return _GSNAP_AMBIG_SUPPORT


def _build_gsnap_cmd(reads_r1, reads_r2, gsnap_dir, genome_version, out_sam,
                    threads, ambig_splice_noclip: bool = True,
                    max_intron: int = 0) -> List[str]:
    cmd = [
        'gsnap',
        '-D', str(gsnap_dir),
        '-d', str(genome_version),
        f'--use-splicing={genome_version}',
        str(reads_r1), str(reads_r2),
        '--output-file', str(out_sam),
        f'--nthreads={threads}',
    ]
    # GSNAP >=2024 dropped --ambig-splice-noclip; omit it when unsupported.
    if ambig_splice_noclip:
        cmd.append('--ambig-splice-noclip')
    cmd += [
        '--novelsplicing=1',
        '--add-paired-nomappers',
        '--sam-extended-cigar',
        '--format=sam',
    ]
    # planning/833 G-8: GSNAP's own --localsplicedist default is 200,000 bp.
    # (--pairmax-rna, the mate-distance cap, is deliberately left alone: it
    # bounds fragment length + intron, not the intron, so max_intron is not the
    # right value for it.)
    if max_intron and max_intron > 0:
        cmd.append(f'--localsplicedist={max_intron}')
    if _is_gz(reads_r1):
        cmd.append('--gunzip')
    return cmd


def _sam_to_name_sorted_bam(sam_path, output_bam, threads: int) -> None:
    """``samtools view -bS | samtools sort -n`` then delete the SAM."""
    sam_path = Path(sam_path)
    view_proc = subprocess.Popen(
        ['samtools', 'view', '-bS', str(sam_path)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-n', '-@', str(threads), '-o', str(output_bam)],
        stdin=view_proc.stdout, stderr=subprocess.PIPE,
    )
    view_proc.stdout.close()
    try:
        _, sort_stderr = sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
        if sort_proc.returncode != 0:
            raise RuntimeError(
                f"samtools sort failed for {sam_path.name}: "
                f"{sort_stderr.decode(errors='replace')[-500:]}"
            )
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        view_proc.kill()
        sort_proc.communicate()
        view_proc.communicate()
        raise RuntimeError(f"samtools sort timed out for {sam_path.name}")
    finally:
        sam_path.unlink(missing_ok=True)
    view_proc.wait()
    if view_proc.returncode not in (0, -9):
        raise RuntimeError(
            f"samtools view failed for {sam_path.name} (exit {view_proc.returncode})"
        )


def _require_binary(name: str) -> str:
    resolved = shutil.which(name)
    if not resolved:
        raise FileNotFoundError(
            f"{name} not found on PATH. Activate the conda env with the COMPASS "
            f"short-read aligner panel installed (STAR, hisat2, magicblast, gsnap)."
        )
    return resolved


def _require_compass_index(label: str, kind: str, path) -> None:
    """Refuse to launch a COMPASS arm whose prebuilt index is absent.

    🔴 Not cosmetic. Nothing in RECTIFY builds the COMPASS indices
    (``_compass_index_paths`` only *derives* their location), and launching an
    aligner against a missing one is not a clean failure:

    * ``STAR`` does fail fast and cleanly (``could not open genome file …``).
    * ``hisat2`` does **not**. Its Perl wrapper re-spawns itself; when
      ``hisat2-align-s`` dies on the missing index the direct child exits but the
      re-spawned wrapper is reparented to init and keeps the stdout/stderr pipe
      write ends open, so the ``subprocess.run(..., capture_output=True)`` here
      never sees EOF. With ``ALIGNER_TIMEOUT = 21600`` (6 h) the whole
      ``run-all`` sits at 0 % CPU until the scheduler kills it, with no
      diagnostic. Measured on Hoffman2, ``run-all --short-read --Scer``
      (planning 861 S3: two ``[perl] <defunct>`` children, two orphan
      ``perl …/hisat2`` at PPID 1, parent in ``poll_schedule_timeout`` holding
      four pipe read fds).

    Checking first turns that hang into the ``DROPPED-ALIGNER`` the caller
    already knows how to report. It is a MITIGATION, not a cure: the same
    deadlock is reachable for any fast-failing wrapper-script aligner (corrupt
    index, unreadable path), so the subprocess plumbing is the real fix
    (KNOWN_ISSUES).
    """
    p = Path(path)
    if kind == 'star':
        ok = (p / 'SAindex').exists() and (p / 'genomeParameters.txt').exists()
        what = f"STAR genome dir {p} (needs SAindex + genomeParameters.txt)"
    elif kind == 'hisat2':
        ok = Path(str(p) + '.1.ht2').exists() or Path(str(p) + '.1.ht2l').exists()
        what = f"HISAT2 index prefix {p} (needs <prefix>.1.ht2)"
    elif kind == 'blast':
        ok = any(Path(str(p) + e).exists() for e in ('.nin', '.nal', '.nsq'))
        what = f"BLAST database {p} (needs <db>.nin/.nsq)"
    elif kind == 'gsnap':
        ok = p.is_dir() and any(p.iterdir())
        what = f"GSNAP genome dir {p}"
    else:  # pragma: no cover - defensive
        ok, what = p.exists(), str(p)
    if not ok:
        raise FileNotFoundError(
            f"{label}: prebuilt COMPASS index missing — {what}. RECTIFY does not "
            f"build COMPASS indices; create it under the genome FASTA's parent "
            f"directory using the COMPASS layout, or drop this arm from "
            f"--base-aligners."
        )


def _finalize_short_read_bam(output_bam: Path, genome_path: str, reads_r1: str,
                            aligner_name: str, qname_to_rn, threads: int,
                            qname_to_umi=None) -> str:
    """Shared post-alignment steps: calmd → QNAME validate → RN (+ optional RX/UMI) inject.

    This is the single choke point for the COMPASS PE panel (STAR/HISAT2/
    magicblast/gsnap). ``qname_to_umi`` carries the standard short-read UMI to the
    BAM as ``RX:Z`` -- necessary because none of these aligners pass FASTQ comments
    through (unlike ``minimap2 -y``). None/empty ⇒ no UMI, original behaviour.
    """
    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    validate_post_alignment_qnames(str(output_bam), str(reads_r1), aligner_name)
    _inject_rn_into_bam(str(output_bam), qname_to_rn, qname_to_umi, reads_path=str(reads_r1))
    logger.info("%s complete: %s", aligner_name, output_bam)
    return str(output_bam)


def run_star(
    reads_path: str,
    reads2_path: Optional[str],
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    read_length: int = 150,
    noncanonical: bool = False,
    star_genome_dir: Optional[str] = None,
    max_intron: int = 0,
) -> str:
    """STAR splice-aware alignment (COMPASS default / non-canonical).

    Paired-end when ``reads2_path`` is given, single-end when it is None
    (TruSeq-style SE mode). noncanonical=True adds ``--scoreGapNoncan 0``
    (COMPASS STAR_noncanonical).
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    qname_to_rn = _load_fastq_rn_map(str(reads_path))
    qname_to_umi = _load_fastq_umi_map(str(reads_path))  # RX:Z UMI, no-op if absent
    _require_binary('STAR')

    paths = _compass_index_paths(genome_path, read_length=read_length)
    star_dir = Path(star_genome_dir) if star_genome_dir else paths.star_dir
    label_pre = 'STAR_noncanonical' if noncanonical else 'STAR_default'
    _require_compass_index(label_pre, 'star', star_dir)
    out_prefix = str(output_bam.with_suffix('')) + '.'  # STAR appends 'Aligned.out.sam'
    cmd = _build_star_cmd(
        reads_path, reads2_path, star_dir, out_prefix, threads,
        read_length=read_length, noncanonical=noncanonical,
        max_intron=max_intron,
    )
    label = 'STAR_noncanonical' if noncanonical else 'STAR_default'
    logger.info("Running %s: %s ...", label, ' '.join(cmd[:6]))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed: {result.stderr[-1000:]}")

    sam_path = Path(out_prefix + 'Aligned.out.sam')
    _sam_to_name_sorted_bam(sam_path, output_bam, threads)
    return _finalize_short_read_bam(
        output_bam, genome_path, reads_path, label, qname_to_rn, threads, qname_to_umi
    )


def run_hisat2(
    reads_path: str,
    reads2_path: Optional[str],
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    min_intron: int = 20,
    max_intron: int = DEFAULT_MAX_INTRON,
    noncanonical: bool = False,
    hisat2_index: Optional[str] = None,
    splice_sites: Optional[str] = None,
) -> str:
    """HISAT2 splice-aware alignment (COMPASS default / non-canonical).

    ``max_intron`` defaults to rectify's own cap, NOT hisat2's 200,000 --
    planning/833 G-8: both dispatchers used to omit it, so HISAT2 kept its own
    default and emitted junctions rectify itself calls impossible (a recurring
    5,584 nt N-op on yeast), which then entered the consensus candidate pool.

    Paired-end when ``reads2_path`` is given, single-end (``-U``) when it is
    None (TruSeq-style SE mode). noncanonical=True adds
    ``--pen-noncansplice 0`` (COMPASS HISAT2_noncanonical).
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    qname_to_rn = _load_fastq_rn_map(str(reads_path))
    qname_to_umi = _load_fastq_umi_map(str(reads_path))  # RX:Z UMI, no-op if absent
    _require_binary('hisat2')

    paths = _compass_index_paths(genome_path)
    idx = hisat2_index if hisat2_index else str(paths.hisat2_index)
    ss = splice_sites if splice_sites else str(paths.splice_sites)
    sam_path = output_bam.with_suffix('.sam')
    label = 'HISAT2_noncanonical' if noncanonical else 'HISAT2_default'
    _require_compass_index(label, 'hisat2', idx)
    # hisat2 accepts a non-existent --known-splicesite-infile SILENTLY (exit 0,
    # normal alignment rate) -- so a half-built index dir yields a green run with
    # zero annotated junctions. Warn loudly rather than fail: novel splicing is
    # still on, and the arm is usable, just not "annotated" (planning 861 S3).
    if not Path(ss).exists():
        logger.warning(
            "%s: known-splicesite file %s does not exist; hisat2 accepts this "
            "silently and the arm will run with NO annotated junctions.", label, ss)
    novel_sj = str(output_bam.with_suffix('.novel_splicesites.txt'))
    summary = str(output_bam.with_suffix('.summary.txt'))
    cmd = _build_hisat2_cmd(
        reads_path, reads2_path, idx, ss, sam_path, threads,
        min_intron, max_intron, novel_sj, summary, noncanonical=noncanonical,
    )
    logger.info("Running %s: %s ...", label, ' '.join(cmd[:6]))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed: {result.stderr[-1000:]}")

    _sam_to_name_sorted_bam(sam_path, output_bam, threads)
    return _finalize_short_read_bam(
        output_bam, genome_path, reads_path, label, qname_to_rn, threads, qname_to_umi
    )


def run_magicblast(
    reads_path: str,
    reads2_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 12,
    blast_index: Optional[str] = None,
    max_intron: int = 0,
) -> str:
    """Magic-BLAST paired-end alignment (splice-site agnostic, COMPASS panel).

    Threads default to 12 — COMPASS deliberately caps Magic-BLAST threads to
    avoid OOM from its high per-thread memory use.
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    qname_to_rn = _load_fastq_rn_map(str(reads_path))
    qname_to_umi = _load_fastq_umi_map(str(reads_path))  # RX:Z UMI, no-op if absent
    _require_binary('magicblast')

    paths = _compass_index_paths(genome_path)
    idx = blast_index if blast_index else str(paths.blast_index)
    _require_compass_index('magicblast', 'blast', idx)
    sam_path = output_bam.with_suffix('.sam')

    # Magic-BLAST gz support is version-dependent (older builds choke on a gzipped
    # -query). Decompress defensively to scratch beside the output so the wrapper
    # is version-independent. RN map / QNAME validation still use the originals.
    _tmp_inputs: List[Path] = []

    def _ensure_plain(p) -> str:
        # Magic-BLAST mangles multi-token FASTQ headers under -no_query_id_trim:
        # rectify-split injects `RN:i:N` plus the original Casava comment, and
        # magicblast spills those tokens into SAM columns 2-3, shifting the
        # mandatory fields so samtools rejects the SAM (FLAG becomes "RN:i:0").
        # Rewrite every header to its bare QNAME (first whitespace token); RN is
        # re-applied afterwards from the qname->RN map in
        # _finalize_short_read_bam. Also decompresses gz (older magicblast chokes
        # on a gzipped -query), so this runs for plain inputs too.
        stem = Path(str(p)[:-3] if _is_gz(p) else str(p)).name
        plain = output_bam.parent / f'.mb_{output_bam.stem}_{stem}.fastq'
        _write_bare_qname_fastq(str(p), str(plain))
        _tmp_inputs.append(plain)
        return str(plain)

    try:
        r1_plain = _ensure_plain(reads_path)
        r2_plain = _ensure_plain(reads2_path)
        cmd = _build_magicblast_cmd(r1_plain, r2_plain, idx, sam_path,
                                    threads=threads, max_intron=max_intron)
        logger.info("Running magicblast: %s ...", ' '.join(cmd[:6]))
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)
        if result.returncode != 0:
            raise RuntimeError(f"magicblast failed: {result.stderr[-1000:]}")
        _sam_to_name_sorted_bam(sam_path, output_bam, threads)
    finally:
        for t in _tmp_inputs:
            t.unlink(missing_ok=True)

    return _finalize_short_read_bam(
        output_bam, genome_path, reads_path, 'magicblast', qname_to_rn, threads, qname_to_umi
    )


def run_gsnap(
    reads_path: str,
    reads2_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    gsnap_genome_dir: Optional[str] = None,
    genome_version: Optional[str] = None,
    max_intron: int = 0,
) -> str:
    """GSNAP paired-end splice-aware alignment (COMPASS panel)."""
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    qname_to_rn = _load_fastq_rn_map(str(reads_path))
    qname_to_umi = _load_fastq_umi_map(str(reads_path))  # RX:Z UMI, no-op if absent
    _require_binary('gsnap')

    paths = _compass_index_paths(genome_path)
    gdir = gsnap_genome_dir if gsnap_genome_dir else str(paths.gsnap_dir)
    gv = genome_version if genome_version else paths.genome_version
    _require_compass_index('gsnap', 'gsnap', Path(gdir) / gv)
    sam_path = output_bam.with_suffix('.sam')
    cmd = _build_gsnap_cmd(reads_path, reads2_path, gdir, gv, sam_path, threads,
                           ambig_splice_noclip=_gsnap_supports_ambig_noclip(),
                           max_intron=max_intron)
    logger.info("Running gsnap: %s ...", ' '.join(cmd[:6]))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)
    if result.returncode != 0:
        raise RuntimeError(f"gsnap failed: {result.stderr[-1000:]}")

    _sam_to_name_sorted_bam(sam_path, output_bam, threads)
    return _finalize_short_read_bam(
        output_bam, genome_path, reads_path, 'gsnap', qname_to_rn, threads, qname_to_umi
    )


def _cs_long_to_cigar(cs: str, query_len: int, query_start: int, query_end: int,
                      strand: str) -> str:
    """Convert gapmm2 long-form cs tag to CIGAR string.

    gapmm2 outputs PAF with long-form cs tags (e.g. :4*at:70-t+g~nn500nn).
    The cs tag only covers the aligned block; query_start and
    query_len - query_end give the soft-clip lengths on each end.

    cs operations:
      :n   = n matches
      *xy  = substitution (1 bp)
      +seq = insertion (query only)
      -seq = deletion (ref only)
      ~nnNnn (or ~NNintronlenNN) = splice junction (N in CIGAR)
    """
    left_clip = query_start
    right_clip = query_len - query_end

    ops = []

    # Soft-clip at 5' end (left for +, right for -)
    if strand == '+':
        if left_clip:
            ops.append((4, left_clip))   # S
    else:
        if right_clip:
            ops.append((4, right_clip))  # S

    # Parse cs string (uses module-level pre-compiled regexes)
    for tok in _CS_TOK_RE.findall(cs):
        if tok[0] == ':':
            ops.append((0, int(tok[1:])))       # M
        elif tok[0] == '*':
            ops.append((8, 1))                  # X (mismatch)
        elif tok[0] == '+':
            ops.append((1, len(tok) - 1))       # I
        elif tok[0] == '-':
            ops.append((2, len(tok) - 1))       # D
        elif tok[0] == '~':
            # ~nnNNNnn where NNN is intron length digits
            intron_len = int(_INTRON_LEN_RE.search(tok).group())
            ops.append((3, intron_len))         # N (splice)

    # Soft-clip at 3' end
    if strand == '+':
        if right_clip:
            ops.append((4, right_clip))  # S
    else:
        if left_clip:
            ops.append((4, left_clip))   # S

    if not ops:
        return '*'

    return ''.join(f'{length}{chr(ord("MIDNSHP=X"[op]))}' for op, length in ops)


def _paf_to_bam(
    paf_path: Path,
    output_bam: Path,
    genome_path: str,
    reads_path: Optional[str] = None,
    threads: int = 8,
):
    """Convert gapmm2 PAF output (with cs tags) to a sorted, indexed BAM.

    gapmm2 only supports PAF/GFF3 output, not SAM. This function converts
    its PAF (which includes long-form cs tags) to BAM via pysam so the
    consensus selection module can work with it normally.

    Args:
        paf_path: Path to the gapmm2 PAF file.
        output_bam: Destination BAM path (name-sorted, indexed).
        genome_path: Genome FASTA (must have a .fai companion).
        reads_path: Optional path to the FASTQ that was aligned.  When
            provided, sequences and base qualities are injected into every
            BAM record (required for ``rectify correct`` downstream).
            DRS headers with ``pt:i:N`` suffixes are stripped automatically
            so name look-ups match the bare UUIDs in the PAF.
        threads: samtools sort thread count.
    """
    import pysam

    seq_dict: Dict[str, Tuple[str, str]] = {}
    if reads_path is not None:
        logger.info(f"Loading read sequences for gapmm2 BAM from {reads_path}")
        seq_dict = _load_fastq_sequences(reads_path)
        logger.info(f"  Loaded {len(seq_dict):,} sequences")

    # Load chromosome lengths from genome FASTA index (fast: reads .fai, not full FASTA)
    chrom_lengths = {}
    with pysam.FastaFile(str(genome_path)) as fa:
        for ref in fa.references:
            chrom_lengths[ref] = fa.get_reference_length(ref)

    # PAF is in input (FASTQ) order — write directly then name-sort.
    # Name-sort so consensus selection can stream across aligners without
    # a secondary sort step.
    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': chrom, 'LN': length} for chrom, length in chrom_lengths.items()],
        'PG': [{'ID': 'gapmm2', 'PN': 'gapmm2'}],
    })

    tmp_unsorted = output_bam.with_suffix('.unsorted.bam')

    n_skipped = 0
    with pysam.AlignmentFile(str(tmp_unsorted), 'wb', header=header) as out_bam:
        with open(paf_path) as fh:
            for line in fh:
                fields = line.rstrip().split('\t')
                if len(fields) < 12:
                    continue

                read_name = fields[0]
                query_len = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                strand = fields[4]
                target_name = fields[5]
                target_start = int(fields[7])
                mapq = int(fields[11])

                if target_name not in chrom_lengths:
                    continue

                # Extract cs tag (build tag dict for O(1) lookup)
                tags = {f[:4]: f[5:] for f in fields[12:] if len(f) > 5}
                cs_tag = tags.get('cs:Z')
                nm_raw = tags.get('NM:i')
                nm_tag = int(nm_raw) if nm_raw is not None else None

                if cs_tag is None:
                    continue

                cigar_str = _cs_long_to_cigar(cs_tag, query_len, query_start, query_end, strand)
                if cigar_str == '*':
                    continue

                flag = 0
                if strand == '-':
                    flag |= 16  # reverse complement

                # PAF tp:A: tag encodes alignment type: 'P'=primary, 'S'=secondary,
                # 'I'=inversion. Mark non-primary PAF records as secondary (FLAG 0x100)
                # so downstream _filtered_read_iterator skips them and only one primary
                # per read_id comes out of the gapmm2 BAM.
                tp_tag = tags.get('tp:A')
                if tp_tag is not None and tp_tag != 'P':
                    flag |= 0x100  # secondary alignment

                # Validate position before writing BAM record (gapmm2 can
                # produce artifacts with negative or out-of-bounds positions)
                if target_start < 0:
                    n_skipped += 1
                    continue
                ref_id = header.get_tid(target_name)
                if ref_id < 0:
                    n_skipped += 1
                    continue

                seg = pysam.AlignedSegment(header)
                seg.query_name = read_name
                seg.flag = flag
                seg.reference_id = ref_id
                seg.reference_start = target_start
                seg.mapping_quality = mapq
                seg.cigarstring = cigar_str
                if nm_tag is not None:
                    seg.set_tag('NM', nm_tag)
                seg.set_tag('cs', cs_tag)

                # Inject sequence and qualities from FASTQ when available.
                # pysam expects query_sequence in alignment orientation:
                # plus-strand → forward FASTQ sequence; minus-strand → rev-comp.
                if seq_dict:
                    entry = seq_dict.get(read_name)
                    if entry is not None:
                        fwd_seq, fwd_qual = entry
                        # Validate CIGAR query-length against the actual sequence.
                        # gapmm2 occasionally emits cs tags that over-consume 1–4
                        # query bases past query_end (gapmm2 bug). If the CIGAR
                        # query length != len(fwd_seq), pysam will reject the record.
                        # Skip these reads rather than attempt surgery on the CIGAR.
                        cigar_qlen = sum(
                            n for n, op in (
                                (int(c[:-1]), c[-1])
                                for c in
                                re.findall(r'\d+[MIDNSHP=X]', cigar_str)
                            )
                            if op in 'MISX='
                        )
                        if cigar_qlen != len(fwd_seq):
                            n_skipped += 1
                            continue
                        if strand == '-':
                            seg.query_sequence = _reverse_complement(fwd_seq)
                            seg.query_qualities = pysam.qualitystring_to_array(fwd_qual)[::-1]
                        else:
                            seg.query_sequence = fwd_seq
                            seg.query_qualities = pysam.qualitystring_to_array(fwd_qual)

                out_bam.write(seg)

    if n_skipped > 0:
        logger.warning(
            f"Skipped {n_skipped} gapmm2 records (invalid position or cs/seq length mismatch) "
            f"in {paf_path.name}"
        )

    # Name-sort the unsorted BAM
    pysam.sort('-n', '-@', str(threads), '-o', str(output_bam), str(tmp_unsorted))
    tmp_unsorted.unlink(missing_ok=True)


def run_gapmm2(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    extra_args: Optional[List[str]] = None,
    max_intron: int = 5000,
) -> str:
    """Run gapmm2 alignment with terminal exon refinement.

    gapmm2 is a minimap2 wrapper that uses edlib to refine alignments at
    terminal exons, improving 5' and 3' end accuracy.

    gapmm2 only outputs PAF (not SAM). The PAF is written to a .paf sidecar
    file and then converted to a sorted BAM via _paf_to_bam() so downstream
    consensus selection can treat it like any other aligner.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_bam: Path for output BAM file
        threads: Number of threads
        extra_args: Additional gapmm2 arguments

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # Check if gapmm2 is available
    gapmm2_path = shutil.which('gapmm2')
    if not gapmm2_path:
        raise FileNotFoundError("gapmm2 not found in PATH")

    # gapmm2 builds a mappy sequence index from the FASTQ and looks up each
    # aligned read's sequence by PAF query name via query_idx.seq(name).
    # Two issues require a cleaned FASTQ to be written before passing to gapmm2:
    #
    # 1. DRS auxiliary-tag headers: `samtools fastq -T pt` produces headers like
    #    "@UUID\tpt:i:25".  minimap2 truncates at the first whitespace → PAF has
    #    bare UUID.  mappy strips tabs similarly, so seq() works in practice.
    #    We strip tags anyway for robustness.
    #
    # 2. Duplicate UUIDs with empty sequences: DRS-trimmed FASTQs occasionally
    #    contain two entries for the same UUID — one with an empty sequence (the
    #    dorado placeholder) and one with the real sequence.  When mappy indexes
    #    a FASTQ with duplicate names, seq() returns None for BOTH entries,
    #    causing a TypeError in gapmm2's refinement loop.  Fix: skip reads with
    #    empty sequences and skip any UUID seen more than once.
    qname_to_rn = _load_fastq_rn_map(str(reads_path))
    logger.info("Preparing deduplicated FASTQ for gapmm2 compatibility")
    tmp_fastq = _clean_fastq(reads_path, output_bam.parent)
    reads_input = str(tmp_fastq)

    # gapmm2 outputs PAF natively — name it .paf (not .sam)
    paf_path = output_bam.with_suffix('.paf')

    cmd = [
        gapmm2_path,
        '-t', str(threads),
        # `-m 1` was here as "Mode 1: direct RNA" — that comment is wrong.
        # `-m`/`--min-mapq` is gapmm2's min-mapq filter (default 1). Some
        # 25.4.5 wheels ship without `type=int` on the argparse arg, so
        # passing `-m 1` made min_mapq the string "1" and TypeError'd on
        # `h.mapq < min_mapq`. Default (int 1) is what we want anyway.
        '-i', str(max_intron),  # Max intron size
        '-o', str(paf_path),
        genome_path,
        reads_input,
    ]

    if extra_args:
        cmd.extend(extra_args)

    logger.info(f"Running gapmm2: {' '.join(cmd[:5])}...")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)
    finally:
        tmp_fastq.unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(f"gapmm2 failed: {result.stderr}")

    if not paf_path.exists() or paf_path.stat().st_size == 0:
        raise RuntimeError(f"gapmm2 produced empty PAF: {paf_path}")

    # Convert PAF → sorted BAM (gapmm2 has no SAM output mode).
    # Pass the original reads_path (with pt:i tags intact) so _paf_to_bam can
    # look up sequences using stripped names (_load_fastq_sequences strips them).
    logger.info(f"Converting gapmm2 PAF to BAM: {output_bam}")
    _paf_to_bam(paf_path, output_bam, genome_path, reads_path=reads_path, threads=threads)

    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    logger.info(f"gapmm2 complete: {output_bam}")
    validate_post_alignment_qnames(str(output_bam), str(reads_path), 'gapmm2')
    _inject_rn_into_bam(str(output_bam), qname_to_rn, reads_path=str(reads_path))
    return str(output_bam)


def run_winnowmap2(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    repetitive_kmers: Optional[str] = None,
    cache_dir: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
    max_intron: int = 5000,
) -> str:
    """Run winnowmap2 alignment for long-read RNA-seq.

    Winnowmap2 uses weighted minimizers to suppress false alignments in repetitive
    regions (e.g., SMN1/SMN2 on human chr5).  It requires a list of high-frequency
    k-mers built from the reference genome via ``meryl``.

    The repetitive k-mer list is cached in cache_dir (or adjacent to the genome)
    and reused on subsequent runs.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_bam: Path for output name-sorted BAM
        threads: Number of threads
        repetitive_kmers: Pre-computed repetitive k-mers file (skips meryl build)
        cache_dir: Directory to cache meryl database and repetitive k-mers
        extra_args: Additional winnowmap arguments
        max_intron: Maximum intron size (-G)

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    genome_path = str(genome_path)

    qname_to_rn = _load_fastq_rn_map(str(reads_path))

    genome_stem = Path(genome_path).name.split('.')[0]
    cache = Path(cache_dir) if cache_dir else Path(genome_path).parent / '.winnowmap_cache'
    rep_kmers_path = Path(repetitive_kmers) if repetitive_kmers else (
        cache / f'{genome_stem}_repetitive_k15.txt'
    )

    if not rep_kmers_path.exists():
        meryl_bin = shutil.which('meryl')
        if not meryl_bin:
            raise RuntimeError(
                "meryl not found on PATH; required for winnowmap2 repetitive-kmer build. "
                "Install: conda install -c bioconda meryl"
            )
        cache.mkdir(parents=True, exist_ok=True)
        meryl_db = cache / f'{genome_stem}_merylDB'
        logger.info(f"winnowmap2: building meryl k15 database at {meryl_db}...")
        count_result = subprocess.run(
            [meryl_bin, 'count', 'k=15', f'output={meryl_db}', genome_path],
            capture_output=True, text=True, timeout=ALIGNER_TIMEOUT,
        )
        if count_result.returncode != 0:
            raise RuntimeError(
                f"meryl count failed (exit {count_result.returncode}): "
                f"{count_result.stderr[-500:]}"
            )
        logger.info("winnowmap2: extracting repetitive k-mers...")
        with open(rep_kmers_path, 'w') as fh:
            print_result = subprocess.run(
                [meryl_bin, 'print', 'greater-than', 'distinct=0.9998', str(meryl_db)],
                stdout=fh, stderr=subprocess.PIPE, text=True, timeout=ALIGNER_TIMEOUT,
            )
        if print_result.returncode != 0:
            rep_kmers_path.unlink(missing_ok=True)
            raise RuntimeError(
                f"meryl print failed (exit {print_result.returncode}): "
                f"{print_result.stderr[-500:]}"
            )
        logger.info(f"winnowmap2: repetitive k-mers → {rep_kmers_path}")

    winnowmap_bin = shutil.which('winnowmap') or shutil.which('winnowmap2')
    if not winnowmap_bin:
        raise RuntimeError(
            "winnowmap not found on PATH. Install: conda install -c bioconda winnowmap"
        )

    cmd = [
        winnowmap_bin,
        '-W', str(rep_kmers_path),
        '-ax', 'splice',
        '-uf',
        # k MUST match the meryl repetitive-kmer build (k=15 above / _repetitive_k15.txt);
        # winnowmap errors "input list of k-mers and winnowmap parameter k are inconsistent"
        # otherwise. 15 is also winnowmap's documented map-ont default.
        '-k15',
        '-G', str(max_intron),
        '--secondary=no',
        '--MD',
        '-y',
        '-t', str(threads),
    ]
    if extra_args:
        cmd.extend(extra_args)
    cmd.extend([genome_path, str(reads_path)])

    logger.info(f"Running winnowmap2: {' '.join(cmd[:8])}...")

    sam_output = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-n', '-@', str(threads), '-o', str(output_bam)],
        stdin=sam_output.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    sam_output.stdout.close()

    import threading
    _wm_err: list = []
    def _drain_wm():
        _wm_err.append(sam_output.stderr.read())
    _dt = threading.Thread(target=_drain_wm, daemon=True)
    _dt.start()

    try:
        _, sort_stderr = sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        sam_output.kill()
        raise RuntimeError(f"winnowmap2/samtools timed out after {ALIGNER_TIMEOUT}s")
    finally:
        _dt.join()

    sam_output.wait()
    if sam_output.returncode != 0:
        wm_err = _wm_err[0].decode(errors='replace') if _wm_err else ''
        raise RuntimeError(
            f"winnowmap2 failed (exit {sam_output.returncode}): {wm_err[-500:]}"
        )
    if sort_proc.returncode != 0:
        raise RuntimeError(f"samtools sort failed: {sort_stderr.decode()}")

    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    logger.info(f"winnowmap2 complete: {output_bam}")
    validate_post_alignment_qnames(str(output_bam), str(reads_path), 'winnowmap2')
    _inject_rn_into_bam(str(output_bam), qname_to_rn, reads_path=str(reads_path))
    return str(output_bam)


def run_minisplice_mm2(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    model_path: Optional[str] = None,
    model_cali_path: Optional[str] = None,
    splice_scores: Optional[str] = None,
    cache_dir: Optional[str] = None,
    annotation_path: Optional[str] = None,
    junc_bonus: int = 9,
    extra_args: Optional[List[str]] = None,
    max_intron: int = 5000,
) -> str:
    """Run minisplice + minimap2 with deep-learning splice-site scores.

    Minisplice scores candidate splice sites across the genome using a trained
    neural network model, then passes those scores to minimap2 via ``--spsc``
    to improve junction accuracy on ONT DRS data.

    The predict step (``minisplice predict``) is cached per genome in cache_dir
    (or adjacent to the genome) and reused on subsequent runs.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_bam: Path for output name-sorted BAM
        threads: Number of threads
        model_path: Path to minisplice model file (e.g. vi2-7k.kan)
        model_cali_path: Optional calibration file (-c flag)
        splice_scores: Pre-computed splice site scores TSV (skips predict step)
        cache_dir: Directory to cache predicted splice scores
        annotation_path: Optional GFF/GTF for minimap2 junction annotation
        junc_bonus: minimap2 junction bonus score
        extra_args: Additional minimap2 arguments
        max_intron: Maximum intron size (-G)

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    genome_path = str(genome_path)

    qname_to_rn = _load_fastq_rn_map(str(reads_path))

    genome_stem = Path(genome_path).name.split('.')[0]
    cache = Path(cache_dir) if cache_dir else Path(genome_path).parent / '.minisplice_cache'
    scores_path = (
        Path(splice_scores) if splice_scores
        else cache / f'{genome_stem}_splice_scores.tsv'
    )

    if not scores_path.exists():
        if not model_path:
            raise RuntimeError(
                "minisplice_mm2 requires either --minisplice-model (path to a "
                "minisplice model file, e.g. vi2-7k.kan) or --minisplice-scores "
                "(pre-computed splice site scores TSV)."
            )
        minisplice_bin = shutil.which('minisplice')
        if not minisplice_bin:
            raise RuntimeError(
                "minisplice not found on PATH. Install: conda install -c bioconda minisplice"
            )
        cache.mkdir(parents=True, exist_ok=True)
        predict_cmd = [minisplice_bin, 'predict', '-t', str(threads)]
        if model_cali_path:
            predict_cmd.extend(['-c', model_cali_path])
        predict_cmd.extend([model_path, genome_path])
        logger.info(f"minisplice_mm2: predicting splice sites → {scores_path}...")
        with open(scores_path, 'w') as fh:
            pred_result = subprocess.run(
                predict_cmd, stdout=fh, stderr=subprocess.PIPE,
                text=True, timeout=ALIGNER_TIMEOUT,
            )
        if pred_result.returncode != 0:
            scores_path.unlink(missing_ok=True)
            raise RuntimeError(
                f"minisplice predict failed (exit {pred_result.returncode}): "
                f"{pred_result.stderr[-500:]}"
            )
        logger.info(f"minisplice_mm2: splice scores → {scores_path}")

    cmd = [
        shutil.which('minimap2') or 'minimap2',
        f'--spsc={scores_path}',
        '-ax', 'splice',
        '-uf',
        '-k14',
        '-G', str(max_intron),
        '--splice-flank=no',
        '--secondary=no',
        '--MD',
        '-y',
        '-t', str(threads),
    ]
    if annotation_path:
        cmd.extend(get_minimap2_junc_args(
            annotation_path=annotation_path,
            junc_bonus=junc_bonus,
            cache_dir=str(cache),
        ))
    if extra_args:
        cmd.extend(extra_args)
    cmd.extend([genome_path, str(reads_path)])

    logger.info(f"Running minisplice_mm2: {' '.join(cmd[:8])}...")

    sam_output = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-n', '-@', str(threads), '-o', str(output_bam)],
        stdin=sam_output.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    sam_output.stdout.close()

    import threading
    _ms_err: list = []
    def _drain_ms():
        _ms_err.append(sam_output.stderr.read())
    _dt = threading.Thread(target=_drain_ms, daemon=True)
    _dt.start()

    try:
        _, sort_stderr = sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        sam_output.kill()
        raise RuntimeError(f"minisplice_mm2/samtools timed out after {ALIGNER_TIMEOUT}s")
    finally:
        _dt.join()

    sam_output.wait()
    if sam_output.returncode != 0:
        ms_err = _ms_err[0].decode(errors='replace') if _ms_err else ''
        raise RuntimeError(
            f"minisplice_mm2 (minimap2) failed (exit {sam_output.returncode}): {ms_err[-500:]}"
        )
    if sort_proc.returncode != 0:
        raise RuntimeError(f"samtools sort failed: {sort_stderr.decode()}")

    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    logger.info(f"minisplice_mm2 complete: {output_bam}")
    validate_post_alignment_qnames(str(output_bam), str(reads_path), 'minisplice_mm2')
    _inject_rn_into_bam(str(output_bam), qname_to_rn, reads_path=str(reads_path))
    return str(output_bam)


def _dedup_gtf_attrs(attr_str: str) -> str:
    """Remove duplicate attribute keys from a GTF attribute string.

    SGD GTFs repeat ``transcript_id`` on mRNA lines (once for the gene
    name, once for a RefSeq/SGD accession).  gffutils rejects these with
    a ``ValueError: more than one value`` when building its database.
    Keep the *first* occurrence of each key.
    """
    import re
    seen: set = set()
    parts: list = []
    for m in re.finditer(r'(\w+)\s+"([^"]+)"', attr_str):
        key = m.group(1)
        if key not in seen:
            seen.add(key)
            parts.append(f'{key} "{m.group(2)}"')
    return '; '.join(parts) + ';' if parts else attr_str


def _subtract_introns_gtf(tx_start: int, tx_end: int, intron_list: list) -> list:
    """Return exon intervals (1-based closed) by subtracting introns from a transcript span."""
    if not intron_list:
        return [(tx_start, tx_end)]
    exons = []
    pos = tx_start
    for i_start, i_end in sorted(intron_list):
        if i_start > pos:
            exons.append((pos, i_start - 1))
        pos = max(pos, i_end + 1)
    if pos <= tx_end:
        exons.append((pos, tx_end))
    return exons


def _normalize_gtf_for_ultra(gtf_path: str, out_path: str) -> None:
    """Convert an SGD-style or gffread-produced GTF to a uLTRA-compatible GTF.

    uLTRA requires ``gene``, ``transcript``, and ``exon`` features with a
    ``gene_id`` attribute on every exon line, and ``gene_id != transcript_id``
    on every exon.  Three common failure modes are handled here:

    1. SGD-style GTFs: ``mRNA`` feature instead of ``transcript``; exons
       absent and must be derived from transcript span minus ``intron`` intervals.
    2. gffread-produced GTFs: ``gene`` features omitted entirely; their absence
       prevents gffutils from building the gene→transcript→exon hierarchy,
       causing ``KeyError: 'gene_id'`` at ``create_augmented_gene.py:323``.
       Synthetic gene records are generated from transcript extents.
    3. gffread-produced GTFs: ``gene_id == transcript_id`` (e.g., both
       ``"YAL069W_mRNA"``).  The isoform suffix is stripped from gene_id so
       the two attributes are distinct.

    SGD naming convention (critical for intron–transcript matching):
    - ``mRNA`` lines carry ``transcript_id "YAL030W"`` (gene name) and
      ``Name "YAL030W_id001"`` (isoform name).
    - ``intron`` and ``CDS`` lines carry ``transcript_id "YAL030W_id001"``
      which equals the mRNA's ``Name`` field, not its ``transcript_id``.

    Changes applied:
    - ``mRNA`` feature → ``transcript``
    - For each transcript: use existing ``exon`` intervals if present;
      otherwise derive from transcript span minus any ``intron`` intervals.
    - Ensure every emitted exon line carries ``gene_id``, sourced from the
      parent transcript's attributes.

    The result is written to *out_path* and can be cached alongside the
    source GTF for repeated use.
    """
    import re
    import collections

    def _get_attr(attr_str: str, key: str) -> str:
        m = re.search(rf'{key} "([^"]+)"', attr_str)
        return m.group(1) if m else ''

    def _ensure_gene_id(attr_str: str, gene_id: str) -> str:
        if not gene_id:
            return attr_str
        if re.search(r'gene_id "', attr_str):
            return attr_str
        return f'gene_id "{gene_id}"; ' + attr_str

    # Suffixes that gffread appends to gene IDs when naming transcripts.
    # Used to recover the parent gene ID from transcript_id when Name is absent.
    _ISOFORM_SUFFIXES = ('_mRNA', '_id001', '_id002', '_id003', '_id004', '_id005')

    # First pass: collect transcript records, introns, and existing exons.
    # Keyed by isoform Name (SGD convention) or transcript_id (gffread GTFs).
    transcripts: dict = {}   # isoform_name → list of 9 GTF fields
    introns: dict = collections.defaultdict(list)   # isoform_name → [(start, end)]
    exons: dict = collections.defaultdict(list)     # isoform_name → [(start, end)]
    gene_lines: list = []

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                gene_lines.append(line)
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature == 'gene':
                gene_lines.append(line)
            elif feature in ('mRNA', 'transcript'):
                # Key by Name (isoform ID) so introns can be joined later.
                # Fall back to transcript_id if Name is absent.
                name = _get_attr(parts[8], 'Name') or _get_attr(parts[8], 'transcript_id')
                if name:
                    transcripts[name] = parts
            elif feature == 'intron':
                tid = _get_attr(parts[8], 'transcript_id')
                if tid:
                    introns[tid].append((int(parts[3]), int(parts[4])))
            elif feature == 'exon':
                tid = _get_attr(parts[8], 'transcript_id')
                if tid:
                    exons[tid].append((int(parts[3]), int(parts[4])))

    # gffread-produced GTFs omit 'gene' features entirely.  uLTRA's gffutils
    # parser requires them to build the gene→transcript→exon hierarchy; without
    # them, exon.attributes["gene_id"] raises KeyError at create_augmented_gene.py:323.
    # Synthesize gene records from transcript extents when none are present.
    _has_gene_features = any('\tgene\t' in ln for ln in gene_lines if not ln.startswith('#'))
    if not _has_gene_features and transcripts:
        _gene_extents: dict = {}  # bare_gene_id → [seqname, min_start, max_end, strand]
        for _n, _p in transcripts.items():
            _gid = _get_attr(_p[8], 'gene_id') or _n
            for _sfx in _ISOFORM_SUFFIXES:
                if _gid.endswith(_sfx):
                    _gid = _gid[:-len(_sfx)]
                    break
            _s, _e = int(_p[3]), int(_p[4])
            if _gid not in _gene_extents:
                _gene_extents[_gid] = [_p[0], _s, _e, _p[6]]
            else:
                _gene_extents[_gid][1] = min(_gene_extents[_gid][1], _s)
                _gene_extents[_gid][2] = max(_gene_extents[_gid][2], _e)
        for _gid, (_seq, _s, _e, _strand) in sorted(
                _gene_extents.items(), key=lambda x: (x[1][0], x[1][1])):
            gene_lines.append(
                f'{_seq}\t.\tgene\t{_s}\t{_e}\t.\t{_strand}\t.'
                f'\tgene_id "{_gid}"; transcript_id "{_gid}";\n'
            )

    with open(out_path, 'w') as fh:
        for line in gene_lines:
            fh.write(line)
        for name, parts in transcripts.items():
            tx_start = int(parts[3])
            tx_end = int(parts[4])
            clean_attrs = _dedup_gtf_attrs(parts[8])
            # uLTRA requires transcript_id != gene_id.  SGD GTFs use the gene
            # name for both, so replace transcript_id with the isoform Name.
            clean_attrs = re.sub(r'transcript_id "[^"]+"', f'transcript_id "{name}"',
                                 clean_attrs, count=1)
            gene_id = _get_attr(clean_attrs, 'gene_id')
            # For gffread-produced GTFs, gene_id carries an isoform suffix
            # (e.g., "YAL069W_mRNA") making it equal to transcript_id.  Strip
            # the suffix so uLTRA sees distinct gene_id and transcript_id.
            if gene_id and gene_id == name:
                for _sfx in _ISOFORM_SUFFIXES:
                    if gene_id.endswith(_sfx):
                        gene_id = gene_id[:-len(_sfx)]
                        clean_attrs = re.sub(r'gene_id "[^"]+"', f'gene_id "{gene_id}"',
                                             clean_attrs, count=1)
                        break
            # Write transcript line
            tx_parts = parts[:]
            tx_parts[2] = 'transcript'
            tx_parts[8] = clean_attrs
            fh.write('\t'.join(tx_parts) + '\n')
            # Use existing exon intervals if present; otherwise derive from
            # transcript span minus introns.  Always ensure gene_id is set.
            if name in exons:
                exon_intervals = sorted(exons[name])
            elif name in introns:
                exon_intervals = _subtract_introns_gtf(tx_start, tx_end, introns[name])
            else:
                exon_intervals = [(tx_start, tx_end)]
            for ex_start, ex_end in exon_intervals:
                ex_parts = parts[:]
                ex_parts[2] = 'exon'
                ex_parts[3] = str(ex_start)
                ex_parts[4] = str(ex_end)
                ex_parts[8] = _ensure_gene_id(clean_attrs, gene_id)
                fh.write('\t'.join(ex_parts) + '\n')




def derive_max_intron(
    annotation_path: Optional[str],
    fallback: int = DEFAULT_MAX_INTRON,
) -> int:
    """Derive a max-intron cap (bp) from an annotation, organism-agnostically.

    Rule: ``2 x longest annotated intron``, rounded UP to the nearest 100,
    clamped to ``_DERIVED_MAX_INTRON_BOUNDS``. The 2x margin admits real
    introns modestly longer than anything annotated while still amputating
    the parameter-cliff class (deSALT's own -I default of 200,000 bp on
    yeast produced a >10 kb junction population whose 99th percentile sat at
    196,914 bp — planning/694: a parameter echo, not biology).

    For the bundled S. cerevisiae R64 annotation the longest annotated
    intron is 2,483 bp (the chrmt Q0060 group-II introns), so the derived
    value is 2*2483 = 4966 -> **5,000** — bit-identical to the historical
    yeast constant, i.e. zero behavior change on existing cohorts while the
    rule itself generalizes to other organisms.

    Intron lengths are taken from explicit ``intron`` features when the
    annotation has them (the R64 GFF does); otherwise they are derived from
    per-transcript exon gaps (the GTF case). Returns ``fallback`` when no
    annotation is given or no intron can be found.
    """
    if not annotation_path:
        return fallback
    path = str(annotation_path)
    try:
        opener = gzip.open if path.endswith('.gz') else open
        max_len = 0
        exons_by_parent: Dict[str, List[Tuple[int, int]]] = {}
        with opener(path, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                if line.startswith('>'):
                    break  # trailing FASTA section (GFF3 ##FASTA convention)
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                ftype = parts[2]
                if ftype == 'intron':
                    max_len = max(max_len, int(parts[4]) - int(parts[3]) + 1)
                elif ftype == 'exon' and max_len == 0:
                    # Only pay the exon-gap bookkeeping while no explicit
                    # intron feature has been seen.
                    attrs = parts[8]
                    parent = None
                    if 'Parent=' in attrs:
                        parent = attrs.split('Parent=')[1].split(';')[0]
                    elif 'transcript_id "' in attrs:
                        parent = attrs.split('transcript_id "')[1].split('"')[0]
                    if parent:
                        exons_by_parent.setdefault(parent, []).append(
                            (int(parts[3]), int(parts[4])))
        if max_len == 0:
            for ivs in exons_by_parent.values():
                if len(ivs) < 2:
                    continue
                ivs.sort()
                for (s1, e1), (s2, e2) in zip(ivs, ivs[1:]):
                    gap = s2 - e1 - 1
                    if gap > max_len:
                        max_len = gap
        if max_len <= 0:
            return fallback
        derived = -(-2 * max_len // 100) * 100  # 2x, rounded up to 100
        lo, hi = _DERIVED_MAX_INTRON_BOUNDS
        return max(lo, min(hi, derived))
    except Exception as exc:  # unreadable/malformed annotation: fall back
        logger.warning(
            "derive_max_intron: could not derive from %s (%s); "
            "falling back to %d", path, exc, fallback)
        return fallback


def run_ultra(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    annotation_path: str,
    threads: int = 8,
    max_intron: int = 5000,
    extra_args: Optional[List[str]] = None
) -> str:
    """Run uLTRA annotation-guided alignment.

    uLTRA uses collinear chaining over a genome graph and excels at aligning
    reads that span small exons (11-20 nt) that seed-chain aligners miss.
    Requires a GTF or GFF annotation file.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_bam: Path for output BAM file
        annotation_path: Path to GFF/GTF annotation (required)
        threads: Number of threads
        extra_args: Additional uLTRA arguments

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    qname_to_rn = _load_fastq_rn_map(str(reads_path))

    ultra_path = shutil.which('uLTRA')
    if not ultra_path:
        raise FileNotFoundError("uLTRA not found in PATH")

    # Ensure namfinder is available — uLTRA calls it as a subprocess.
    # If not on PATH, prepend the vendored binary's directory to the env.
    _extra_env: Optional[dict] = None
    if not shutil.which('namfinder'):
        vendored_namfinder = _get_vendored_binary('namfinder')
        if vendored_namfinder:
            logger.info(f"namfinder not in PATH; using vendored binary: {vendored_namfinder}")
            _extra_env = _os.environ.copy()
            _extra_env['PATH'] = str(Path(vendored_namfinder).parent) + _os.pathsep + _extra_env.get('PATH', '')
        else:
            raise FileNotFoundError(
                "namfinder not found in PATH and no vendored binary available. "
                "Install with: conda install -c bioconda namfinder"
            )

    # uLTRA writes to an output directory; the SAM lives at <out_dir>/<prefix>.sam
    ultra_out_dir = output_bam.parent / f"{output_bam.stem}_ultra_tmp"
    ultra_out_dir.mkdir(parents=True, exist_ok=True)

    # If a prior run left a corrupted cache (refs_sequences.fa empty or missing),
    # remove the stale database so uLTRA re-indexes from scratch.  A valid cache
    # (non-empty refs_sequences.fa) is intentionally preserved — the index is
    # genome/GTF-derived and can be shared across chunks.
    _refs_fa = ultra_out_dir / 'refs_sequences.fa'
    _db = ultra_out_dir / 'database.db'
    if _db.exists() and (not _refs_fa.exists() or _refs_fa.stat().st_size == 0):
        logger.warning(f"uLTRA: stale/empty cache detected in {ultra_out_dir}; removing to force re-index")
        import shutil as _shutil_pre
        _shutil_pre.rmtree(ultra_out_dir)
        ultra_out_dir.mkdir(parents=True, exist_ok=True)
    prefix = "ultra"

    # uLTRA does not support gzipped inputs and requires GTF (not GFF) format.
    _tmp_dir = None
    ref_path = genome_path

    def _is_gzipped(p: str) -> bool:
        return str(p).endswith('.gz')

    def _is_gff(p: str) -> bool:
        return any(str(p).endswith(ext) for ext in ('.gff', '.gff3', '.gff.gz', '.gff3.gz'))

    # Resolve annotation: GFF(gz) → look for sibling .gtf; GTF(gz) → decompress.
    ann_p = Path(annotation_path)
    if _is_gff(annotation_path):
        stem = ann_p.stem if not ann_p.stem.endswith('.gz') else ann_p.stem[:-3]
        stem = stem.rsplit('.', 1)[0] if stem.endswith('.gff') or stem.endswith('.gff3') else stem
        candidate_gtf = ann_p.parent / (stem + '.gtf')
        if not candidate_gtf.exists():
            for ext in ('.gff.gz', '.gff3.gz', '.gff', '.gff3'):
                if str(ann_p).endswith(ext):
                    candidate_gtf = ann_p.parent / (str(ann_p.name)[:-len(ext)] + '.gtf')
                    break
        if not candidate_gtf.exists():
            raise FileNotFoundError(
                f"uLTRA requires a sibling GTF alongside the GFF but none was found: "
                f"{annotation_path}. Expected: {candidate_gtf}. "
                f"For RECTIFY built-in genomes this file is bundled; for custom genomes "
                f"provide a GTF from the same annotation source (NOT gffread -T output, "
                f"which omits gene features required by uLTRA)."
            )
        ann_path = str(candidate_gtf)
        logger.info(f"uLTRA: using GTF annotation at {ann_path}")
    else:
        ann_path = annotation_path

    # Decompress genome if gzipped (uLTRA cannot read gzip FASTA)
    if _is_gzipped(genome_path):
        _tmp_dir = _tempfile.mkdtemp(prefix='ultra_decomp_')
        ref_dest = Path(_tmp_dir) / Path(genome_path).stem
        with gzip.open(genome_path, 'rb') as f_in, open(ref_dest, 'wb') as f_out:
            f_out.write(f_in.read())
        ref_path = str(ref_dest)

    # Decompress annotation if gzipped GTF
    if _is_gzipped(ann_path):
        if _tmp_dir is None:
            _tmp_dir = _tempfile.mkdtemp(prefix='ultra_decomp_')
        ann_dest = Path(_tmp_dir) / Path(ann_path).stem
        with gzip.open(ann_path, 'rb') as f_in, open(ann_dest, 'wb') as f_out:
            f_out.write(f_in.read())
        ann_path = str(ann_dest)

    # Normalize GTF for uLTRA if:
    #   (a) the annotation came from GFF discovery (the sibling .gtf may have
    #       been generated by gffread, which produces exon features but omits
    #       gene_id on ~10% of them — crashing uLTRA at
    #       create_augmented_gene.py:323), or
    #   (b) the GTF contains no 'exon' features at all (SGD-style GTFs use
    #       'mRNA' instead of 'transcript' and derive exons from intron spans).
    # Cache the normalized GTF alongside the source so it is only built once.
    _needs_norm = _is_gff(annotation_path)  # always normalize GFF-sourced GTFs
    if not _needs_norm:
        with open(ann_path) as _fh:
            for _line in _fh:
                if _line.startswith('#'):
                    continue
                _parts = _line.split('\t')
                if len(_parts) >= 3 and _parts[2] == 'exon':
                    break
            else:
                _needs_norm = True

    if _needs_norm:
        _ann_p = Path(ann_path)
        _norm_gtf = _ann_p.parent / (_ann_p.stem + '.ultra_norm.gtf')
        # Rebuild cache if source GTF is newer than the cached normalized GTF
        if not _norm_gtf.exists() or _ann_p.stat().st_mtime > _norm_gtf.stat().st_mtime:
            logger.info(f"uLTRA: normalizing GTF (mRNA→transcript, deriving exons) → {_norm_gtf}")
            _normalize_gtf_for_ultra(ann_path, str(_norm_gtf))
        else:
            logger.info(f"uLTRA: using cached normalized GTF at {_norm_gtf}")
        ann_path = str(_norm_gtf)

    # uLTRA pipeline = index + align in one step
    # --disable_infer: skip gffutils gene-boundary inference (avoids crash on complex GFFs)
    cmd = [
        ultra_path, 'pipeline',
        '--ont',
        '--disable_infer',
        '--t', str(threads),
        '--prefix', prefix,
    ]

    # No-trim policy: aligners must not DROP query bases — unaligned bases belong
    # in soft-clips (counted against the aligner, then rectify rescues), never
    # deleted from SEQ. uLTRA's --reduce_read_ployA (default 8) reduces poly-A
    # runs >X bp to 5 bp "before MEM finding"; despite its help claiming this
    # "does not affect final read alignments", it actually TRUNCATES the emitted
    # SEQ (a 640 bp read came back 632 bp — 8 bases gone, not soft-clipped).
    # Disable it (value > any read length) so the full query is preserved; the
    # previously-dropped bases reappear as soft-clip. Verified: this only extends
    # the terminal soft-clip — the aligned core, reference_start, and all other
    # reads' alignments are byte-identical. Skip if the caller set it explicitly.
    if not (extra_args and any('--reduce_read_ployA' in str(a) for a in extra_args)):
        cmd += ['--reduce_read_ployA', str(_ULTRA_DISABLE_POLYA_REDUCE)]

    # Max intron. uLTRA and deSALT were the ONLY two wrappers in the panel that
    # never received rectify's ``max_intron`` — every other one does (minimap2
    # -G, BBMap maxindel, STAR --max-intronlen, gapmm2 -i, GMAP
    # --max-intronlength-*). They are also exactly the two arms that emit
    # physically impossible introns, and that is not a coincidence: uncapped,
    # deSALT runs at its own -I default of 200,000 bp on a genome whose longest
    # real intron is ~1 kb. Measured across 31 cDNA samples (planning/694), the
    # >10 kb junction population has its 99th percentile at 196,914 bp and only
    # 0.55 % above 200,000 — a parameter cliff, not biology. Capping here removes
    # the class at source, which is why the downstream 10 kb guard (d0e3a0f)
    # should not have to be the primary control.
    if not (extra_args and any('--max_intron' in str(a) for a in extra_args)):
        cmd += ['--max_intron', str(max_intron)]

    cmd += [
        ref_path,
        ann_path,
        reads_path,
        str(ultra_out_dir),
    ]

    if extra_args:
        cmd.extend(extra_args)

    logger.info(f"Running uLTRA: {' '.join(cmd[:6])}...")

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT,
                            env=_extra_env)  # _extra_env is None (inherit) or has namfinder prepended

    # Clean up decompressed temp files
    if _tmp_dir:
        import shutil as _shutil2
        _shutil2.rmtree(_tmp_dir, ignore_errors=True)

    if result.returncode != 0:
        raise RuntimeError(f"uLTRA failed: {result.stderr}")

    sam_path = ultra_out_dir / f"{prefix}.sam"
    if not sam_path.exists() or sam_path.stat().st_size == 0:
        raise RuntimeError(f"uLTRA produced no output at {sam_path}")

    # Convert SAM → coordinate-sorted BAM
    view_proc = subprocess.Popen(
        ['samtools', 'view', '-bS', str(sam_path)],
        stdout=subprocess.PIPE
    )
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-@', str(threads), '-o', str(output_bam)],
        stdin=view_proc.stdout
    )
    view_proc.stdout.close()
    try:
        sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        view_proc.kill()
        raise RuntimeError(f"samtools sort (uLTRA) timed out after {ALIGNER_TIMEOUT}s")

    view_proc.wait()
    if view_proc.returncode != 0:
        raise RuntimeError(f"samtools view (uLTRA) failed with exit code {view_proc.returncode}")
    if sort_proc.returncode != 0:
        raise RuntimeError(f"samtools sort (uLTRA) failed with exit code {sort_proc.returncode}")

    import shutil as _shutil
    _shutil.rmtree(ultra_out_dir, ignore_errors=True)

    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    logger.info(f"uLTRA complete: {output_bam}")
    validate_post_alignment_qnames(str(output_bam), str(reads_path), 'uLTRA')
    _inject_rn_into_bam(str(output_bam), qname_to_rn, reads_path=str(reads_path))
    return str(output_bam)


def _dedup_desalt_bam(input_bam: Path, output_bam: Path, threads: int = 1) -> None:
    """Remove duplicate alignments from deSALT output.

    deSALT has a known bug: it outputs each alignment N times where N is the
    number of secondary alignment slots (-N flag, default 4). This removes
    duplicates by keeping the first occurrence of each (read_name, flag,
    chrom, pos, cigar) combination.
    """
    import pysam

    seen: set = set()
    n_total = 0
    n_kept = 0

    with pysam.AlignmentFile(str(input_bam), 'rb') as bam_in, \
         pysam.AlignmentFile(str(output_bam), 'wb', header=bam_in.header) as bam_out:
        for read in bam_in:
            n_total += 1
            key = (
                read.query_name,
                read.flag,
                read.reference_name,
                read.reference_start,
                read.cigarstring,
            )
            if key not in seen:
                seen.add(key)
                bam_out.write(read)
                n_kept += 1

    n_removed = n_total - n_kept
    if n_removed > 0:
        logger.info(
            f"deSALT dedup: {n_total} → {n_kept} alignments "
            f"({n_removed} duplicates removed)"
        )


def run_desalt(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    annotation_path: Optional[str] = None,
    threads: int = 8,
    index_path: Optional[str] = None,
    max_intron: int = 5000,
    extra_args: Optional[List[str]] = None
) -> str:
    """Run deSALT De Bruijn graph splice aligner.

    deSALT requires a pre-built RdBG index (build with: deSALT index <ref.fa> <index_dir>).
    If index_path is not given, looks for a 'desalt_index' directory adjacent to genome_path.

    Note: deSALT has a known bug where it outputs each primary alignment
    N times. This function deduplicates automatically.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA (used only for default index discovery)
        output_bam: Path for output BAM file
        annotation_path: Optional GTF annotation (accepted for API compatibility; unused —
                         deSALT's -G flag causes SIGSEGV on yeast GTF)
        threads: Number of threads
        index_path: Path to pre-built deSALT RdBG index directory
        extra_args: Additional deSALT arguments

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    desalt_exec = shutil.which('deSALT')
    if not desalt_exec:
        # Fall back to vendored binary bundled with the package
        desalt_exec = _get_vendored_desalt()
    if not desalt_exec:
        raise FileNotFoundError(
            "deSALT not found in PATH and no compatible vendored binary available. "
            "Install with: rectify install-aligners --desalt"
        )

    # Resolve index directory. The index must contain a non-empty `ref.seq`
    # (deBGA's main reference output) — an empty `desalt_index/` placeholder
    # silently picks up here otherwise and deSALT fails late with a confusing
    # "index load" error.
    def _is_built_desalt_index(p: Path) -> bool:
        return p.is_dir() and (p / 'ref.seq').exists() and (p / 'ref.seq').stat().st_size > 0

    if index_path:
        index_dir = Path(index_path)
        if not _is_built_desalt_index(index_dir):
            raise FileNotFoundError(
                f"deSALT index at {index_dir} is missing or empty. "
                f"Build it with: deSALT index <genome.fa> {index_dir}"
            )
    else:
        genome_p = Path(genome_path)
        candidates = [
            genome_p.parent / 'desalt_index',
            genome_p.parent / f"{genome_p.name.split('.')[0]}_desalt_index",
        ]
        index_dir = next((c for c in candidates if _is_built_desalt_index(c)), None)
        if index_dir is None:
            tried = ', '.join(str(c) for c in candidates)
            raise FileNotFoundError(
                f"deSALT index not found (or empty) adjacent to {genome_path}. "
                f"Looked at: {tried}. "
                f"Build it with: deSALT index <genome.fa> {candidates[0]}"
            )

    qname_to_rn = _load_fastq_rn_map(str(reads_path))

    # deSALT v1.5.6 silently misparses gzipped FASTQ (reads binary gz header as
    # plain text, producing garbage "reads" and exit code 0) and stops at the
    # first empty-sequence record (Dorado placeholder reads in DRS-trimmed FASTQs).
    # Apply the same cleaning as run_gapmm2(): decompress if needed, strip DRS
    # auxiliary tags from read names, skip empty-sequence and duplicate UUID records.
    logger.info("Preparing deduplicated FASTQ for deSALT compatibility")
    tmp_cleaned_fastq = _clean_fastq(reads_path, output_bam.parent)
    reads_input = str(tmp_cleaned_fastq)

    sam_path = output_bam.with_suffix('.sam')
    # deSALT requires its -f tmp file on a local (non-NFS) filesystem to avoid
    # "double free or corruption" crashes caused by memory-mapping over NFS.
    tmp_file = Path(_tempfile.gettempdir()) / f"desalt_tmp_{_os.getpid()}_{output_bam.stem}.bin"

    cmd = [
        desalt_exec, 'aln',
        '-t', str(threads),
        '-f', str(tmp_file),
        '-o', str(sam_path),
    ]

    # NOTE: deSALT's -G annotation flag causes a SIGSEGV when loading yeast GTF.
    # Skip -G entirely; deSALT de novo splice detection is sufficient.
    _ = annotation_path

    # Max intron (-I). deSALT's own default is 200,000 bp; uncapped it is the
    # dominant source of the impossible-intron class — 11 of the 12 longest
    # junctions in a cDNA sample were deSALT-only, single-read, anchored by
    # 15-57 bp of exon against a ~200 kb gap (planning/694). See the matching
    # note in run_ultra() for why these two wrappers were the ones missing it.
    if not (extra_args and any(str(a) in ('-I', '--max-intron-len') for a in extra_args)):
        cmd += ['-I', str(max_intron)]

    if extra_args:
        cmd.extend(extra_args)

    cmd.extend([str(index_dir), reads_input])

    logger.info(f"Running deSALT: {' '.join(cmd[:6])}...")

    # Strip LD_LIBRARY_PATH so deSALT uses system glibc only (conda allocator
    # incompatibility causes "double free or corruption" in parallel launches).
    desalt_env = _os.environ.copy()
    desalt_env.pop('LD_LIBRARY_PATH', None)

    result = subprocess.run(cmd, capture_output=True, text=True, env=desalt_env, timeout=ALIGNER_TIMEOUT)
    tmp_file.unlink(missing_ok=True)
    for sidecar in tmp_file.parent.glob(f"{tmp_file.name}*"):
        sidecar.unlink(missing_ok=True)
    tmp_cleaned_fastq.unlink(missing_ok=True)
    if result.returncode != 0:
        sam_path.unlink(missing_ok=True)
        logger.debug("deSALT exit code: %d (in crash set: %s)",
                     result.returncode, result.returncode in _DESALT_CRASH_EXITS)
        if result.returncode in _DESALT_CRASH_EXITS:
            # Deterministic SIGSEGV/OOM in second-pass Loop-ProcessReads.
            # Retrying the same input never helps.  Emit an empty BAM so
            # merge_aligners proceeds with the other 4 aligners for this chunk.
            logger.warning(
                "deSALT crashed (exit %d, likely SIGSEGV in Loop-ProcessReads) — "
                "emitting empty BAM; chunk will use 4-aligner consensus. "
                "Upstream bug: github.com/ydLiu-HIT/deSALT",
                result.returncode,
            )
            _create_empty_name_sorted_bam(output_bam, genome_path)
            return str(output_bam)
        raise RuntimeError(f"deSALT failed (exit {result.returncode}): {result.stderr}")

    if not sam_path.exists() or sam_path.stat().st_size == 0:
        sam_path.unlink(missing_ok=True)
        raise RuntimeError(f"deSALT produced no output at {sam_path}")

    # Convert SAM → coordinate-sorted BAM, then deduplicate
    raw_bam = output_bam.with_suffix('.raw.bam')
    view_proc = subprocess.Popen(
        ['samtools', 'view', '-bS', str(sam_path)],
        stdout=subprocess.PIPE
    )
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-@', str(threads), '-o', str(raw_bam)],
        stdin=view_proc.stdout
    )
    view_proc.stdout.close()
    try:
        sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        view_proc.kill()
        raise RuntimeError(f"samtools sort (deSALT) timed out after {ALIGNER_TIMEOUT}s")

    view_proc.wait()
    if view_proc.returncode != 0 or sort_proc.returncode != 0:
        # deSALT occasionally writes malformed SAM (CIGAR/seq length mismatch)
        # that samtools rejects.  Same deterministic-input behaviour as the
        # SIGSEGV case — retrying never recovers.  Emit an empty BAM.
        raw_bam.unlink(missing_ok=True)
        sam_path.unlink(missing_ok=True)
        logger.warning(
            "deSALT wrote invalid SAM (samtools view exit %d / sort exit %d — "
            "likely CIGAR mismatch; upstream bug: github.com/ydLiu-HIT/deSALT/issues/49) — "
            "emitting empty BAM; chunk will use 4-aligner consensus.",
            view_proc.returncode, sort_proc.returncode,
        )
        _create_empty_name_sorted_bam(output_bam, genome_path)
        return str(output_bam)

    sam_path.unlink(missing_ok=True)

    _dedup_desalt_bam(raw_bam, output_bam, threads=threads)
    raw_bam.unlink(missing_ok=True)

    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    logger.info(f"deSALT complete: {output_bam}")
    validate_post_alignment_qnames(str(output_bam), str(reads_path), 'deSALT')
    _inject_rn_into_bam(str(output_bam), qname_to_rn, reads_path=str(reads_path))
    return str(output_bam)


def run_gmap(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    annotation_path: Optional[str] = None,
    threads: int = 8,
    gmap_db: Optional[str] = None,
    gmap_path: str = 'gmap',
    max_intron: int = 5000,
    extra_args: Optional[List[str]] = None,
) -> str:
    """Run GMAP, the splice-aware long-read / cDNA aligner (GSNAP's sibling).

    GMAP needs a pre-built database (build once with:
    ``gmap_build -D <dir> -d <name> genome.fa``). Pass the built db directory via
    ``gmap_db`` (``<dir>/<name>``); otherwise a ``gmap_db/<genome_stem>`` directory
    adjacent to ``genome_path`` is used. The db is NOT auto-built here — building is
    a one-time setup step (slow, and racy if many users trigger it at once).

    On any GMAP/samtools failure this emits an empty name-sorted BAM (like deSALT)
    so the consensus proceeds with the remaining aligners rather than aborting.

    Args:
        reads_path: Path to FASTQ (gz ok; cleaned/decompressed internally)
        genome_path: Genome FASTA (used for default db discovery + calmd)
        output_bam: Output BAM path
        annotation_path: Accepted for API compatibility; unused (GMAP is de-novo)
        threads: Threads
        gmap_db: Path to a built GMAP db directory (``<-D dir>/<-d name>``)
        gmap_path: gmap executable (default 'gmap')
        max_intron: Max intron length (yeast-safe default 5000)
        extra_args: Extra gmap arguments

    Returns:
        Path to output BAM
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    gmap_exec = shutil.which(gmap_path) or shutil.which('gmap')
    if not gmap_exec:
        raise FileNotFoundError(
            "gmap not found in PATH. Install with: rectify install-aligners --gmap "
            "(or conda install -c bioconda gmap)"
        )

    # Resolve the built GMAP db. A built db dir <dir>/<name>/ contains
    # <name>.genomecomp (deterministic marker of a complete gmap_build).
    def _is_built_gmap_db(p: Path) -> bool:
        return p.is_dir() and (
            (p / f"{p.name}.genomecomp").exists() or any(p.glob('*.genomecomp'))
        )

    if gmap_db:
        db_dir = Path(gmap_db)
        if not _is_built_gmap_db(db_dir):
            raise FileNotFoundError(
                f"GMAP db at {db_dir} is missing or incomplete. "
                f"Build it with: gmap_build -D {db_dir.parent} -d {db_dir.name} <genome.fa>"
            )
    else:
        genome_p = Path(genome_path)
        stem = genome_p.name.split('.')[0]
        candidates = [
            genome_p.parent / 'gmap_db' / stem,
            genome_p.parent / 'gmap_db' / genome_p.parent.name,
        ]
        db_dir = next((c for c in candidates if _is_built_gmap_db(c)), None)
        if db_dir is None:
            tried = ', '.join(str(c) for c in candidates)
            raise FileNotFoundError(
                f"GMAP db not found adjacent to {genome_path}. Looked at: {tried}. "
                f"Build it with: gmap_build -D <dir> -d <name> <genome.fa>"
            )

    qname_to_rn = _load_fastq_rn_map(str(reads_path))

    # GMAP reads plain FASTQ; reuse the shared cleaner (decompress, strip DRS
    # auxiliary read-name tags, drop empty-seq / duplicate-UUID records).
    logger.info("Preparing FASTQ for GMAP")
    tmp_cleaned_fastq = _clean_fastq(reads_path, output_bam.parent)
    reads_input = str(tmp_cleaned_fastq)

    sam_path = output_bam.with_suffix('.sam')
    cmd = [
        gmap_exec,
        '-D', str(db_dir.parent),
        '-d', db_dir.name,
        '-t', str(threads),
        '-f', 'samse',                 # SAM output, single-end
        '--npaths', '1',               # best path only (avoid multimapper flood)
        '--max-intronlength-middle', str(max_intron),
        '--max-intronlength-ends', str(max_intron),
    ]
    if extra_args:
        cmd.extend(extra_args)
    cmd.append(reads_input)

    logger.info(f"Running GMAP: {' '.join(cmd[:8])}...")

    try:
        with open(sam_path, 'w') as _sam_out:
            result = subprocess.run(
                cmd, stdout=_sam_out, stderr=subprocess.PIPE, text=True,
                timeout=ALIGNER_TIMEOUT,
            )
    finally:
        tmp_cleaned_fastq.unlink(missing_ok=True)

    if result.returncode != 0:
        sam_path.unlink(missing_ok=True)
        logger.warning(
            "GMAP failed (exit %d) — emitting empty BAM; chunk uses the remaining "
            "aligners. stderr: %s",
            result.returncode, (result.stderr or '')[-500:],
        )
        _create_empty_name_sorted_bam(output_bam, genome_path)
        return str(output_bam)

    if not sam_path.exists() or sam_path.stat().st_size == 0:
        sam_path.unlink(missing_ok=True)
        logger.warning("GMAP produced no output — emitting empty BAM.")
        _create_empty_name_sorted_bam(output_bam, genome_path)
        return str(output_bam)

    # SAM → coordinate-sorted BAM
    view_proc = subprocess.Popen(
        ['samtools', 'view', '-bS', str(sam_path)], stdout=subprocess.PIPE
    )
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-@', str(threads), '-o', str(output_bam)],
        stdin=view_proc.stdout,
    )
    view_proc.stdout.close()
    try:
        sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        view_proc.kill()
        raise RuntimeError(f"samtools sort (GMAP) timed out after {ALIGNER_TIMEOUT}s")
    view_proc.wait()
    if view_proc.returncode != 0 or sort_proc.returncode != 0:
        sam_path.unlink(missing_ok=True)
        logger.warning(
            "GMAP SAM→BAM failed (view %s / sort %s) — emitting empty BAM.",
            view_proc.returncode, sort_proc.returncode,
        )
        _create_empty_name_sorted_bam(output_bam, genome_path)
        return str(output_bam)

    sam_path.unlink(missing_ok=True)
    _apply_calmd_eq(output_bam, genome_path, threads=threads)
    logger.info(f"GMAP complete: {output_bam}")
    validate_post_alignment_qnames(str(output_bam), str(reads_path), 'gmap')
    _inject_rn_into_bam(str(output_bam), qname_to_rn, reads_path=str(reads_path))
    return str(output_bam)


def run_multi_aligner(
    reads_path: str,
    genome_path: str,
    output_dir: str,
    sample_name: str,
    annotation_path: Optional[str] = None,
    config: Optional[MultiAlignerConfig] = None,
    threads: int = 8,
    aligners: Optional[List[str]] = None,
    reads2_path: Optional[str] = None,
    read_length: int = 150,
    max_intron: int = DEFAULT_MAX_INTRON,
) -> Dict[str, str]:
    """Run multiple aligners on the same reads.

    Args:
        reads_path: Path to FASTQ file (R1 for paired-end)
        genome_path: Path to genome FASTA
        output_dir: Output directory for BAM files
        sample_name: Sample name for output files
        annotation_path: Optional GFF/GTF for junction annotation
        config: Optional MultiAlignerConfig
        threads: Number of threads per aligner
        aligners: List of aligners to run (default: all enabled)
        reads2_path: Mate-2 (R2) FASTQ for paired-end short-read aligners
            (bbmap, STAR_*, HISAT2_*, magicblast, gsnap). When None those
            aligners run single-end where supported.
        read_length: Read length for STAR sjdbOverhang / index selection.
        max_intron: Intron-length cap handed to every splice-aware arm
            (planning/833 G-8). Callers should pass the value derived from the
            annotation by ``derive_max_intron``; the default is the historical
            S. cerevisiae constant, never an aligner's own default.

    Returns:
        Dict mapping aligner name to output BAM path
    """
    if config is None:
        config = MultiAlignerConfig()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = output_dir / '.cache'

    results = {}

    # Determine which aligners to run
    if aligners is None:
        aligners = []
        if config.minimap2.enabled:
            aligners.append('minimap2')
        if config.map_pacbio.enabled:
            aligners.append('mapPacBio')
        if config.gapmm2.enabled:
            aligners.append('gapmm2')
        if config.bbmap.enabled:
            aligners.append('bbmap')
        if config.bwa_mem.enabled:
            aligners.append('bwa')
        if config.overhang_resolver.enabled:
            # Must come after minimap2 (consumes its BAM).
            aligners.append('overhang_resolver')

    # COMPASS short-read aligners are run paired-end; they need a mate-2 FASTQ.
    _paired_required = {
        'STAR_default', 'STAR_noncanonical', 'HISAT2_default',
        'HISAT2_noncanonical', 'magicblast', 'gsnap',
    }

    for aligner in aligners:
        output_bam = output_dir / f"{sample_name}.{aligner}.sorted.bam"

        if aligner in _paired_required and not reads2_path:
            logger.error(
                "Aligner %s requires a mate-2 (R2) FASTQ (reads2_path); skipping. "
                "Pass paired chunk FASTQs from `rectify split --read2`.", aligner,
            )
            continue

        try:
            if aligner == 'minimap2':
                results['minimap2'] = run_minimap2(
                    reads_path=reads_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    annotation_path=annotation_path if config.use_junction_annotation else None,
                    junc_bonus=config.junc_bonus,
                    cache_dir=str(cache_dir),
                    extra_args=config.minimap2.extra_args
                )
            elif aligner == 'mapPacBio':
                results['mapPacBio'] = run_map_pacbio(
                    reads_path=reads_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    extra_args=config.map_pacbio.extra_args
                )
            elif aligner == 'gapmm2':
                results['gapmm2'] = run_gapmm2(
                    reads_path=reads_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    extra_args=config.gapmm2.extra_args
                )
            elif aligner == 'overhang_resolver':
                # Native mapPacBio-role replacement (planning/641): re-places
                # terminal soft clips of the minimap2 arm across junctions
                # under an information bound, so it needs that BAM first.
                if 'minimap2' not in results:
                    logger.error(
                        "overhang_resolver requires the minimap2 arm from this "
                        "invocation — list minimap2 before overhang_resolver; "
                        "skipping."
                    )
                    continue
                from .overhang_resolver import run_overhang_resolver
                results['overhang_resolver'] = run_overhang_resolver(
                    base_bam=results['minimap2'],
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                )
            elif aligner == 'bbmap':
                results['bbmap'] = run_bbmap(
                    reads_path=reads_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    extra_args=config.bbmap.extra_args,
                    reads2_path=reads2_path,
                )
            elif aligner == 'bwa':
                results['bwa'] = run_bwa_mem(
                    reads_path=reads_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    extra_args=config.bwa_mem.extra_args
                )
            elif aligner in ('STAR_default', 'STAR_noncanonical'):
                results[aligner] = run_star(
                    reads_path=reads_path,
                    reads2_path=reads2_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    read_length=read_length,
                    noncanonical=(aligner == 'STAR_noncanonical'),
                    max_intron=max_intron,
                )
            elif aligner in ('HISAT2_default', 'HISAT2_noncanonical'):
                results[aligner] = run_hisat2(
                    reads_path=reads_path,
                    reads2_path=reads2_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    noncanonical=(aligner == 'HISAT2_noncanonical'),
                    max_intron=max_intron,
                )
            elif aligner == 'magicblast':
                results['magicblast'] = run_magicblast(
                    reads_path=reads_path,
                    reads2_path=reads2_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=min(threads, 12),
                    max_intron=max_intron,
                )
            elif aligner == 'gsnap':
                results['gsnap'] = run_gsnap(
                    reads_path=reads_path,
                    reads2_path=reads2_path,
                    genome_path=genome_path,
                    output_bam=str(output_bam),
                    threads=threads,
                    max_intron=max_intron,
                )
            else:
                logger.warning(f"Unknown aligner: {aligner}")

        except Exception as e:
            logger.error(f"Aligner {aligner} failed: {e}")

    if 'minimap2' in aligners and 'minimap2' not in results:
        succeeded = ', '.join(sorted(results)) or 'none'
        raise RuntimeError(
            "Required baseline aligner minimap2 failed; cannot proceed with "
            f"consensus selection (succeeded: {succeeded})"
        )
    if not results:
        raise RuntimeError("All requested aligners failed; no BAMs are available for consensus selection")

    logger.info(f"Multi-aligner pipeline complete: {len(results)} aligners succeeded")
    return results
