"""
Multi-aligner alignment pipeline for RECTIFY.

Runs multiple aligners (minimap2, mapPacBio, gapmm2) on the same reads
and provides utilities for consensus analysis.

Strategy:
- minimap2: Fast seed-and-chain baseline with junction annotation support
- mapPacBio: BBTools long-read aligner with splice-aware mode
- gapmm2: minimap2 wrapper with terminal exon refinement

Junction annotations are used to IMPROVE alignment quality but scoring
remains BLIND to annotations (novel junctions can still be detected).

Author: Kevin R. Roy
"""

import re
import subprocess
import logging
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field

# Pre-compiled regexes for gapmm2 PAF cs-tag parsing (avoid per-call recompilation)
_CS_TOK_RE = re.compile(r':[0-9]+|[*][a-z][a-z]|[+][a-z]+|[-][a-z]+|~[a-z]{2}[0-9]+[a-z]{2}')
_INTRON_LEN_RE = re.compile(r'[0-9]+')

# Maximum wall-clock seconds to wait for any single aligner subprocess.
# uLTRA and deSALT on a large genome can take 30–60 min; 7200 s (2 h) is a
# safe upper bound before treating a run as hung.
ALIGNER_TIMEOUT = 7200

from ..utils.junction_bed import get_minimap2_junc_args, generate_junction_bed

from rectify.core.mpb_split_reads import split_long_reads, stitch_split_bam, MAX_MPB_READ_LENGTH

logger = logging.getLogger(__name__)


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

    # Junction annotation options
    use_junction_annotation: bool = True
    junc_bonus: int = 9

    # Output options
    keep_sam: bool = False


def check_aligner_available(aligner: str) -> bool:
    """Check if an aligner is available in PATH."""
    return shutil.which(aligner) is not None


def run_minimap2(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    annotation_path: Optional[str] = None,
    junc_bonus: int = 9,
    cache_dir: Optional[str] = None,
    extra_args: Optional[List[str]] = None
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
        '-G', '5000',  # Max intron size
        '--splice-flank=no',  # Disable for compatibility
        '--secondary=no',
        '--MD',  # Include MD tag for indel correction and alignment identity
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
    try:
        _, stderr = sort_proc.communicate(timeout=ALIGNER_TIMEOUT)
    except subprocess.TimeoutExpired:
        sort_proc.kill()
        sam_output.kill()
        raise RuntimeError(f"minimap2/samtools timed out after {ALIGNER_TIMEOUT}s")

    if sort_proc.returncode != 0:
        raise RuntimeError(f"minimap2/samtools failed: {stderr.decode()}")

    logger.info(f"minimap2 complete: {output_bam}")
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
            header = fin.readline()
            if not header:
                break
            seq  = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
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


def run_map_pacbio(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    extra_args: Optional[List[str]] = None,
    chunk_idx: Optional[int] = None,
    n_chunks: Optional[int] = None,
) -> str:
    """Run BBTools mapPacBio alignment.

    mapPacBio is splice-aware for long reads. The mapPacBio.sh script default
    of fastareadlen=6000 causes AssertionErrors on reads >6019 bp; we patch
    the script default to 100000 at install time (see multi_aligner.py).
    intronlen=50 converts deletions >=50bp to N (intron) CIGAR operations.

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
            logger.info("mapPacBio merge complete: %s", output_bam)
            return str(output_bam)
        else:
            missing = [p.name for p in chunk_bams if not p.exists()]
            logger.warning(
                "%d/%d chunk BAMs missing (%s) — falling back to full alignment",
                len(missing), n_chunks, ', '.join(missing),
            )
            # fall through to full alignment below

    # ── Chunk extraction mode: write only reads for this chunk ─────────────
    actual_reads_path = reads_path
    chunk_tmp_fq = None
    if chunk_idx is not None and n_chunks and n_chunks > 1:
        # Redirect output to the chunk BAM path
        output_bam = output_bam.parent / f"{output_bam.stem}.chunk_{chunk_idx}_of_{n_chunks}.bam"
        chunk_tmp_fq = str(output_bam.with_suffix('.chunk_input.fastq'))
        extract_fastq_chunk(reads_path, chunk_tmp_fq, chunk_idx, n_chunks)
        actual_reads_path = chunk_tmp_fq

    # Check if mapPacBio is available
    map_pacbio_path = shutil.which('mapPacBio.sh')
    if not map_pacbio_path:
        # Try common anaconda location
        map_pacbio_path = '/home/groups/larsms/users/kevinroy/anaconda3/bin/mapPacBio.sh'
        if not Path(map_pacbio_path).exists():
            if chunk_tmp_fq:
                Path(chunk_tmp_fq).unlink(missing_ok=True)
            raise FileNotFoundError("mapPacBio.sh not found")

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
        'fastareadlen=100000',  # Belt-and-suspenders: also patched in mapPacBio.sh default
        'intronlen=50',         # Convert large deletions to introns
        'minratio=0.4',
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

    logger.info(
        "Running mapPacBio: %s",
        ' '.join(cmd[:5]) + (f' [chunk {chunk_idx}/{n_chunks}]' if chunk_idx is not None else ''),
    )

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)
    if chunk_tmp_fq:
        Path(chunk_tmp_fq).unlink(missing_ok=True)
    if result.returncode != 0:
        mpb_split_fq.unlink(missing_ok=True)
        raise RuntimeError(f"mapPacBio failed: {result.stderr}")

    # Convert SAM to name-sorted BAM. Name-sort (-n) so consensus selection
    # can stream across aligners without a secondary sort step.
    view_proc = subprocess.Popen(
        ['samtools', 'view', '-bS', str(sam_path)],
        stdout=subprocess.PIPE,
    )
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-n', '-@', str(threads), '-o', str(output_bam)],
        stdin=view_proc.stdout,
    )
    view_proc.stdout.close()
    sort_proc.communicate()
    sam_path.unlink(missing_ok=True)
    if sort_proc.returncode != 0:
        raise RuntimeError(
            f"samtools sort (mapPacBio) failed with code {sort_proc.returncode}"
        )

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

    logger.info("mapPacBio complete: %s", output_bam)
    return str(output_bam)


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


def _paf_to_bam(paf_path: Path, output_bam: Path, genome_path: str, threads: int = 8):
    """Convert gapmm2 PAF output (with cs tags) to a sorted, indexed BAM.

    gapmm2 only supports PAF/GFF3 output, not SAM. This function converts
    its PAF (which includes long-form cs tags) to BAM via pysam so the
    consensus selection module can work with it normally.
    """
    import pysam

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

                out_bam.write(seg)

    if n_skipped > 0:
        logger.warning(
            f"Skipped {n_skipped} gapmm2 records with invalid positions in {paf_path.name}"
        )

    # Name-sort the unsorted BAM
    pysam.sort('-n', '-@', str(threads), '-o', str(output_bam), str(tmp_unsorted))
    tmp_unsorted.unlink(missing_ok=True)


def run_gapmm2(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    extra_args: Optional[List[str]] = None
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
        gapmm2_path = '/home/groups/larsms/users/kevinroy/anaconda3/bin/gapmm2'
        if not Path(gapmm2_path).exists():
            raise FileNotFoundError("gapmm2 not found")

    # gapmm2 outputs PAF natively — name it .paf (not .sam)
    paf_path = output_bam.with_suffix('.paf')

    cmd = [
        gapmm2_path,
        '-t', str(threads),
        '-m', '1',   # Mode 1: direct RNA
        '-i', '5000',  # Max intron size
        '-o', str(paf_path),
        genome_path,
        reads_path,
    ]

    if extra_args:
        cmd.extend(extra_args)

    logger.info(f"Running gapmm2: {' '.join(cmd[:5])}...")

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)
    if result.returncode != 0:
        raise RuntimeError(f"gapmm2 failed: {result.stderr}")

    if not paf_path.exists() or paf_path.stat().st_size == 0:
        raise RuntimeError(f"gapmm2 produced empty PAF: {paf_path}")

    # Convert PAF → sorted BAM (gapmm2 has no SAM output mode)
    logger.info(f"Converting gapmm2 PAF to BAM: {output_bam}")
    _paf_to_bam(paf_path, output_bam, genome_path, threads=threads)

    logger.info(f"gapmm2 complete: {output_bam}")
    return str(output_bam)


def run_ultra(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    annotation_path: str,
    threads: int = 8,
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
    import gzip
    import tempfile

    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    ultra_path = shutil.which('uLTRA')
    if not ultra_path:
        ultra_path = '/home/users/kevinroy/.local/bin/uLTRA'
        if not Path(ultra_path).exists():
            raise FileNotFoundError("uLTRA not found in PATH or ~/.local/bin")

    # uLTRA writes to an output directory; the SAM lives at <out_dir>/<prefix>.sam
    ultra_out_dir = output_bam.parent / f"{output_bam.stem}_ultra_tmp"
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
                f"uLTRA requires GTF annotation but only GFF was found: {annotation_path}. "
                f"Generate a GTF with: rectify build-gtf --annotation {annotation_path} "
                f"-o {candidate_gtf}"
            )
        ann_path = str(candidate_gtf)
        logger.info(f"uLTRA: using GTF annotation at {ann_path}")
    else:
        ann_path = annotation_path

    # Decompress genome if gzipped (uLTRA cannot read gzip FASTA)
    if _is_gzipped(genome_path):
        _tmp_dir = tempfile.mkdtemp(prefix='ultra_decomp_')
        ref_dest = Path(_tmp_dir) / Path(genome_path).stem
        with gzip.open(genome_path, 'rb') as f_in, open(ref_dest, 'wb') as f_out:
            f_out.write(f_in.read())
        ref_path = str(ref_dest)

    # Decompress annotation if gzipped GTF
    if _is_gzipped(ann_path):
        if _tmp_dir is None:
            _tmp_dir = tempfile.mkdtemp(prefix='ultra_decomp_')
        ann_dest = Path(_tmp_dir) / Path(ann_path).stem
        with gzip.open(ann_path, 'rb') as f_in, open(ann_dest, 'wb') as f_out:
            f_out.write(f_in.read())
        ann_path = str(ann_dest)

    # uLTRA pipeline = index + align in one step
    # --disable_infer: skip gffutils gene-boundary inference (avoids crash on complex GFFs)
    cmd = [
        ultra_path, 'pipeline',
        '--ont',
        '--disable_infer',
        '--t', str(threads),
        '--prefix', prefix,
        ref_path,
        ann_path,
        reads_path,
        str(ultra_out_dir),
    ]

    if extra_args:
        cmd.extend(extra_args)

    logger.info(f"Running uLTRA: {' '.join(cmd[:6])}...")

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=ALIGNER_TIMEOUT)

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

    import shutil as _shutil
    _shutil.rmtree(ultra_out_dir, ignore_errors=True)

    logger.info(f"uLTRA complete: {output_bam}")
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
        desalt_exec = '/home/users/kevinroy/bin/deSALT'
        if not Path(desalt_exec).exists():
            raise FileNotFoundError("deSALT not found in PATH or ~/bin")

    # Resolve index directory
    if index_path:
        index_dir = Path(index_path)
    else:
        genome_p = Path(genome_path)
        candidates = [
            genome_p.parent / 'desalt_index',
            genome_p.parent / f"{genome_p.name.split('.')[0]}_desalt_index",
        ]
        index_dir = next((c for c in candidates if c.exists()), None)
        if index_dir is None:
            raise FileNotFoundError(
                f"deSALT index not found adjacent to {genome_path}. "
                f"Build it with: deSALT index <genome.fa> {candidates[0]}"
            )

    sam_path = output_bam.with_suffix('.sam')
    # deSALT requires its -f tmp file on a local (non-NFS) filesystem to avoid
    # "double free or corruption" crashes caused by memory-mapping over NFS.
    import tempfile as _tmpmod
    import os as _os
    tmp_file = Path(_tmpmod.gettempdir()) / f"desalt_tmp_{_os.getpid()}_{output_bam.stem}.bin"

    cmd = [
        desalt_exec, 'aln',
        '-t', str(threads),
        '-f', str(tmp_file),
        '-o', str(sam_path),
    ]

    # NOTE: deSALT's -G annotation flag causes a SIGSEGV when loading yeast GTF.
    # Skip -G entirely; deSALT de novo splice detection is sufficient.
    _ = annotation_path

    if extra_args:
        cmd.extend(extra_args)

    cmd.extend([str(index_dir), reads_path])

    logger.info(f"Running deSALT: {' '.join(cmd[:6])}...")

    # Strip LD_LIBRARY_PATH so deSALT uses system glibc only (conda allocator
    # incompatibility causes "double free or corruption" in parallel launches).
    desalt_env = _os.environ.copy()
    desalt_env.pop('LD_LIBRARY_PATH', None)

    result = subprocess.run(cmd, capture_output=True, text=True, env=desalt_env, timeout=ALIGNER_TIMEOUT)
    tmp_file.unlink(missing_ok=True)
    for sidecar in tmp_file.parent.glob(f"{tmp_file.name}*"):
        sidecar.unlink(missing_ok=True)
    if result.returncode != 0:
        raise RuntimeError(f"deSALT failed: {result.stderr}")

    if not sam_path.exists() or sam_path.stat().st_size == 0:
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
    sam_path.unlink(missing_ok=True)

    _dedup_desalt_bam(raw_bam, output_bam, threads=threads)
    raw_bam.unlink(missing_ok=True)

    logger.info(f"deSALT complete: {output_bam}")
    return str(output_bam)


def run_multi_aligner(
    reads_path: str,
    genome_path: str,
    output_dir: str,
    sample_name: str,
    annotation_path: Optional[str] = None,
    config: Optional[MultiAlignerConfig] = None,
    threads: int = 8,
    aligners: Optional[List[str]] = None
) -> Dict[str, str]:
    """Run multiple aligners on the same reads.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_dir: Output directory for BAM files
        sample_name: Sample name for output files
        annotation_path: Optional GFF/GTF for junction annotation
        config: Optional MultiAlignerConfig
        threads: Number of threads per aligner
        aligners: List of aligners to run (default: all enabled)

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

    for aligner in aligners:
        output_bam = output_dir / f"{sample_name}.{aligner}.sorted.bam"

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
            else:
                logger.warning(f"Unknown aligner: {aligner}")

        except Exception as e:
            logger.error(f"Aligner {aligner} failed: {e}")

    logger.info(f"Multi-aligner pipeline complete: {len(results)} aligners succeeded")
    return results
