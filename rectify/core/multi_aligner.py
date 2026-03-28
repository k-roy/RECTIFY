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

import subprocess
import logging
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field

from ..utils.junction_bed import get_minimap2_junc_args, generate_junction_bed

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

    # Run minimap2, pipe to samtools
    sam_output = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    # Sort and convert to BAM
    sort_cmd = [
        'samtools', 'sort',
        '-@', str(threads),
        '-o', str(output_bam)
    ]

    sort_proc = subprocess.Popen(
        sort_cmd,
        stdin=sam_output.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    # Wait for completion
    sam_output.stdout.close()
    _, stderr = sort_proc.communicate()

    if sort_proc.returncode != 0:
        raise RuntimeError(f"minimap2/samtools failed: {stderr.decode()}")

    # Index BAM
    subprocess.run(['samtools', 'index', str(output_bam)], check=True)

    logger.info(f"minimap2 complete: {output_bam}")
    return str(output_bam)


def run_map_pacbio(
    reads_path: str,
    genome_path: str,
    output_bam: str,
    threads: int = 8,
    extra_args: Optional[List[str]] = None
) -> str:
    """Run BBTools mapPacBio alignment.

    mapPacBio is splice-aware for long reads. fastareadlen=100000 overrides
    the default 6000 limit that causes AssertionErrors on reads >6019 bp.
    intronlen=50 converts deletions >=50bp to N (intron) CIGAR operations.

    Args:
        reads_path: Path to FASTQ file
        genome_path: Path to genome FASTA
        output_bam: Path for output BAM file
        threads: Number of threads
        extra_args: Additional mapPacBio arguments

    Returns:
        Path to output BAM file
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # Check if mapPacBio is available
    map_pacbio_path = shutil.which('mapPacBio.sh')
    if not map_pacbio_path:
        # Try common anaconda location
        map_pacbio_path = '/home/groups/larsms/users/kevinroy/anaconda3/bin/mapPacBio.sh'
        if not Path(map_pacbio_path).exists():
            raise FileNotFoundError("mapPacBio.sh not found")

    # mapPacBio outputs SAM to stdout
    sam_path = output_bam.with_suffix('.sam')

    cmd = [
        map_pacbio_path,
        f'ref={genome_path}',
        f'in={reads_path}',
        f'out={sam_path}',
        f'threads={threads}',
        'fastareadlen=100000',  # Override 6000 default; avoids AssertionError on reads >6019 bp
        'intronlen=50',  # Convert large deletions to introns
        'nodisk',
        'minratio=0.4',
    ]

    if extra_args:
        cmd.extend(extra_args)

    logger.info(f"Running mapPacBio: {' '.join(cmd[:5])}...")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"mapPacBio failed: {result.stderr}")

    # Convert to sorted BAM
    view_proc = subprocess.Popen(
        ['samtools', 'view', '-bS', str(sam_path)],
        stdout=subprocess.PIPE
    )
    sort_proc = subprocess.Popen(
        ['samtools', 'sort', '-@', str(threads), '-o', str(output_bam)],
        stdin=view_proc.stdout
    )
    view_proc.stdout.close()
    sort_proc.communicate()

    # Index
    subprocess.run(['samtools', 'index', str(output_bam)], check=True)

    # Clean up SAM
    sam_path.unlink(missing_ok=True)

    logger.info(f"mapPacBio complete: {output_bam}")
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
    import re

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

    # Parse cs string
    for tok in re.findall(r':[0-9]+|[*][a-z][a-z]|[+][a-z]+|[-][a-z]+|~[a-z]{2}[0-9]+[a-z]{2}', cs):
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
            intron_len = int(re.search(r'[0-9]+', tok).group())
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
    import gzip

    # Load chromosome lengths from genome FASTA for BAM header
    chrom_lengths = {}
    opener = gzip.open if str(genome_path).endswith('.gz') else open
    current_chrom = None
    current_len = 0
    with opener(genome_path, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_chrom:
                    chrom_lengths[current_chrom] = current_len
                current_chrom = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
        if current_chrom:
            chrom_lengths[current_chrom] = current_len

    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'unsorted'},
        'SQ': [{'SN': chrom, 'LN': length} for chrom, length in chrom_lengths.items()],
        'PG': [{'ID': 'gapmm2', 'PN': 'gapmm2'}],
    })

    tmp_unsorted = output_bam.with_suffix('.unsorted.bam')

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
                target_len = int(fields[6])
                target_start = int(fields[7])
                target_end = int(fields[8])
                mapq = int(fields[11])

                if target_name not in chrom_lengths:
                    continue

                # Extract cs tag
                cs_tag = None
                nm_tag = None
                for f in fields[12:]:
                    if f.startswith('cs:Z:'):
                        cs_tag = f[5:]
                    elif f.startswith('NM:i:'):
                        nm_tag = int(f[5:])

                if cs_tag is None:
                    continue

                cigar_str = _cs_long_to_cigar(cs_tag, query_len, query_start, query_end, strand)
                if cigar_str == '*':
                    continue

                # Build SAM flag
                flag = 0
                if strand == '-':
                    flag |= 16  # reverse complement

                seg = pysam.AlignedSegment(header)
                seg.query_name = read_name
                seg.flag = flag
                seg.reference_id = header.get_tid(target_name)
                seg.reference_start = target_start
                seg.mapping_quality = mapq
                seg.cigarstring = cigar_str
                if nm_tag is not None:
                    seg.set_tag('NM', nm_tag)
                seg.set_tag('cs', cs_tag)

                out_bam.write(seg)

    # Sort and index
    pysam.sort('-@', str(threads), '-o', str(output_bam), str(tmp_unsorted))
    pysam.index(str(output_bam))
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

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"gapmm2 failed: {result.stderr}")

    if not paf_path.exists() or paf_path.stat().st_size == 0:
        raise RuntimeError(f"gapmm2 produced empty PAF: {paf_path}")

    # Convert PAF → sorted BAM (gapmm2 has no SAM output mode)
    logger.info(f"Converting gapmm2 PAF to BAM: {output_bam}")
    _paf_to_bam(paf_path, output_bam, genome_path, threads=threads)

    logger.info(f"gapmm2 complete: {output_bam}")
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
