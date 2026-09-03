"""``rectify run-all --netseq`` — the NET-seq / mNET-seq FASTQ→tracks pipeline.

Three stages and *no* ``correct``/consensus stage::

    [1] trim    cutadapt: 3' linker, quality trim, min length.  THE RANDOMER STAYS IN THE READ.
    [2] align   STAR with the planning-829 §4 ABSOLUTE match floors, optionally against a
                combined genome + spike-in index; spike-in counted from idxstats.
    [3] netseq  `rectify netseq --rna3p-at read5p` on the primary-genome BAM, with the
                planning-836 defaults (donor-side junction rescue on, tracks on the CORRECTED
                3' end, evidence-gated walkback).

**Why no `correct` stage** (Kevin, planning 835 §2, 2026-09-02 17:40). The junction-overhang
resolver is a *discovery* stage wrapping minimap2 and belongs to the long-read protocols; the
1-bp-and-up soft-clip RE-PLACEMENT that NET-seq needs is `correct`'s logic, and planning 836 built
exactly that logic — donor-side junction-pool rescue against the exon-2 start, randomer-tolerant
remainder, evidence-gated poly(A) walkback — *inside* ``rectify netseq``, where it operates on the
NET-seq geometry directly. Running the COMPASS panel plus a `correct` arm on 35-nt NET-seq reads
would add no placement information and cost the `correct` arm, which is 87x the whole align stage
per BAM (planning 586 / 728).

Everything a test needs to check without STAR or cutadapt on PATH is a **pure command builder**
(``build_cutadapt_cmd``, ``build_star_index_cmd``, ``build_netseq_star_cmd``,
``parse_idxstats_spikein``, ``build_netseq_namespace``). The 829 §4 floors are the expensive part
of this file — they cost 42 % of one library's "unique" alignments to learn — so they live in one
named table, ``NETSEQ_STAR_FLOORS``, that a unit test asserts flag by flag.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = [
    'CHURCHMAN_3P_LINKER',
    'NETSEQ_STAR_FLOORS',
    'build_cutadapt_cmd',
    'build_star_index_cmd',
    'build_netseq_star_cmd',
    'build_combined_genome',
    'parse_idxstats_spikein',
    'build_netseq_namespace',
    'run_netseq_pipeline',
]


#: Churchman 2011 3' linker, measured in 81.4 % of PRJNA1521488 reads (planning 829 §2).
CHURCHMAN_3P_LINKER = 'ATCTCGTATGCCGTCTTCTGCTTG'


#: 🔴 The planning-829 §4 ABSOLUTE match floors, and the reason this mode exists at all.
#:
#: STAR's *relative* thresholds (``--outFilter{Score,Match}NminOverLread 0.3``) let a trimmed 35-nt
#: NET-seq read be placed on an **11-nt chance match**. Measured on PRJNA1521488: 24 % of
#: asf1_rep1's "unique" alignments were <= 20 nt, 9.4 % of its reads sat on ONE chrVI locus, and
#: 54 % of its "spliced" reads had a 9-nt block. The absolute floors below cost WT ~2 % of
#: alignments and removed the chrVI hotspot entirely (antisense-to-CDS 0.544 -> 0.866).
#:
#: The relative thresholds are KEPT at 0.3 alongside them -- STAR applies both, and dropping them
#: would change the validated command line. Ordered as a flag->args mapping so a test can assert
#: each floor by name rather than pattern-matching a command string.
NETSEQ_STAR_FLOORS: Dict[str, Tuple[str, ...]] = {
    # Absolute floors (829 §4) — the chance-match fix
    '--outFilterMatchNmin': ('20',),            # >= 20 matched nt, regardless of read length
    '--outFilterMismatchNoverLmax': ('0.06',),  # <= 6 % mismatches
    '--outFilterMismatchNmax': ('3',),          # ... and never more than 3
    '--alignSJoverhangMin': ('12',),            # >= 12 nt on the short side of a junction
    '--outFilterIntronMotifs': ('RemoveNoncanonical',),
    # Relative thresholds, kept for parity with the validated 829 command line
    '--outFilterScoreMinOverLread': ('0.3',),
    '--outFilterMatchNminOverLread': ('0.3',),
    # Multimappers: undepleted Churchman NET-seq is ~63 % chrXII rDNA multimappers (829 §3), and
    # `ambiguous=best`-style arbitrary placement manufactures signal. Drop them.
    '--outFilterMultimapNmax': ('1',),
    # Local, not EndToEnd: the randomer and any exon-2 overhang MUST be soft-clipped and left for
    # the planning-836 donor-side rescue to re-place. EndToEnd would force them onto the genome.
    '--alignEndsType': ('Local',),
    # Yeast introns; the mammalian defaults let a 35-nt read span a fabricated 500-kb "intron".
    '--alignIntronMin': ('20',),
    '--alignIntronMax': ('2000',),
}

#: SAM attributes STAR must write. MD is load-bearing: the terminal-oligo(A) X-mismatch correction
#: in ``netseq_bam_processor`` falls back to parsing MD for ``M`` ops.
_NETSEQ_STAR_SAM_ATTRS = ('NH', 'HI', 'AS', 'nM', 'NM', 'MD')

_FASTA_EXTS = ('.fasta.gz', '.fa.gz', '.fna.gz', '.fasta', '.fa', '.fna')


# ══════════════════════════════════════════════════════════════════════════════
# Pure command builders
# ══════════════════════════════════════════════════════════════════════════════

def build_cutadapt_cmd(
    reads: Path,
    output: Optional[Path],
    adapter: str = CHURCHMAN_3P_LINKER,
    min_length: int = 18,
    nextseq_trim: Optional[int] = 20,
    quality_cutoff: Optional[int] = None,
    threads: int = 4,
    cutadapt_path: str = 'cutadapt',
) -> List[str]:
    """cutadapt command for the NET-seq 3' linker trim.

    🔴 **There is no ``-u``.** The 5' randomer is LEFT IN THE READ. planning 829 §4 measured
    PRJNA1521488 as a MIXTURE — 58-60 % of reads carry no randomer at all, 21-22 % carry a 6-nt one
    — so an unconditional ``-u 6`` deletes six genuine nucleotides from the majority class, and the
    randomer-free class is exactly the one whose 3' end is already correct. The aligner soft-clips
    the randomer (``--alignEndsType Local``) and ``rectify netseq --umi-length N`` accounts for it
    downstream, per-read, as a non-templated remainder.

    ``output=None`` writes to stdout, which is how the aligner stage streams trimmed reads into
    STAR through a process substitution without a trimmed FASTQ ever hitting disk.
    """
    cmd = [cutadapt_path, '-j', str(threads), '-a', adapter]
    if quality_cutoff is not None:
        cmd += ['-q', str(quality_cutoff)]
    elif nextseq_trim is not None:
        # Two-colour-aware quality trim: on NextSeq/NovaSeq a no-signal base reads as a
        # high-quality G, so plain -q does not remove the poly-G tail. 829 used this on all
        # 52 libraries.
        cmd += [f'--nextseq-trim={nextseq_trim}']
    cmd += ['-m', str(min_length), '--report=minimal', '-o', str(output) if output else '-',
            str(reads)]
    return cmd


def build_star_index_cmd(
    genome_fasta: Path,
    index_dir: Path,
    threads: int = 4,
    sa_index_nbases: int = 11,
    star_path: str = 'STAR',
) -> List[str]:
    """``STAR --runMode genomeGenerate`` for a small (yeast-scale) genome.

    ``--genomeSAindexNbases 11`` is the STAR-manual value for a ~12 Mb genome
    (``min(14, log2(genomeLength)/2 - 1)``); the mammalian default 14 makes STAR allocate for a
    3-Gb suffix array and it will not build. No ``--sjdbGTFfile``: the 829 index is
    annotation-free and junctions are found de novo under ``--outFilterIntronMotifs
    RemoveNoncanonical``. Annotation enters later, in the rescue pool.
    """
    return [
        star_path,
        '--runMode', 'genomeGenerate',
        '--runThreadN', str(threads),
        '--genomeDir', str(index_dir),
        '--genomeFastaFiles', str(genome_fasta),
        '--genomeSAindexNbases', str(sa_index_nbases),
        '--outFileNamePrefix', str(index_dir) + '/',
    ]


def build_netseq_star_cmd(
    reads_in: str,
    index_dir: Path,
    out_prefix: Path,
    threads: int = 4,
    gzipped: bool = False,
    star_path: str = 'STAR',
    extra: Optional[Sequence[str]] = None,
) -> List[str]:
    """STAR alignment for NET-seq, carrying every ``NETSEQ_STAR_FLOORS`` entry.

    ``reads_in`` is a string, not a Path, so the caller can pass a process-substitution path
    (``/dev/fd/63``) and stream cutadapt straight into STAR.
    """
    cmd = [
        star_path,
        '--runThreadN', str(threads),
        '--genomeDir', str(index_dir),
        '--readFilesIn', reads_in,
        '--outSAMtype', 'BAM', 'Unsorted',
        '--outSAMunmapped', 'None',
        '--outSAMattributes', *_NETSEQ_STAR_SAM_ATTRS,
        '--outFileNamePrefix', str(out_prefix),
    ]
    if gzipped:
        cmd += ['--readFilesCommand', 'zcat']
    for flag, values in NETSEQ_STAR_FLOORS.items():
        cmd += [flag, *values]
    if extra:
        cmd += list(extra)
    return cmd


def _open_maybe_gz(path: Path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'r')


def build_combined_genome(
    genome_fasta: Path,
    spikein_fasta: Path,
    output_fasta: Path,
    spikein_prefix: str = 'Sp_',
) -> Dict[str, int]:
    """Write ``genome + prefixed(spikein)`` to ``output_fasta``; return contig counts.

    Ported from ``planning/829_scripts/829_star_index.sh`` (validated on 52 libraries): the
    spike-in headers are rewritten ``>I Schizosaccharomyces_pombe`` -> ``>Sp_I`` so that a single
    STAR index carries both organisms and ``samtools idxstats`` separates them by contig prefix.
    The primary genome is copied VERBATIM, so its contig names — and therefore every downstream
    coordinate — are unchanged.

    Fails loudly rather than writing an index that silently contains one organism.
    """
    n_primary = n_spikein = 0
    tmp = output_fasta.with_suffix(output_fasta.suffix + '.part')
    with open(tmp, 'w') as out:
        with _open_maybe_gz(genome_fasta) as fh:
            for line in fh:
                if line.startswith('>'):
                    name = line[1:].split()[0] if len(line) > 1 else ''
                    if name.startswith(spikein_prefix):
                        raise ValueError(
                            f"primary genome contig {name!r} already starts with the spike-in "
                            f"prefix {spikein_prefix!r} — the idxstats split would count it as "
                            "spike-in. Pass a different --spikein-prefix."
                        )
                    n_primary += 1
                out.write(line)
        with _open_maybe_gz(spikein_fasta) as fh:
            for line in fh:
                if line.startswith('>'):
                    n_spikein += 1
                    out.write('>' + spikein_prefix + line[1:].split()[0] + '\n')
                else:
                    out.write(line)
    if n_primary == 0:
        raise ValueError(f"no contigs read from the primary genome {genome_fasta}")
    if n_spikein == 0:
        raise ValueError(f"no contigs read from the spike-in genome {spikein_fasta}")
    tmp.replace(output_fasta)
    return {'primary_contigs': n_primary, 'spikein_contigs': n_spikein}


def parse_idxstats_spikein(idxstats_text: str, spikein_prefix: str = 'Sp_') -> Dict[str, object]:
    """Split a ``samtools idxstats`` table into primary / spike-in mapped-read counts.

    idxstats columns: ``name  length  mapped  unmapped``. The ``*`` row (unmapped, unplaced) is
    excluded from both. Returns the counts, the fraction, and the contig lists — the fraction is
    the spike-in normalisation factor the Bryll & Peterson headline depends on, which RPM erases
    by construction (planning 834 gap D).
    """
    primary = spikein = 0
    primary_contigs: List[str] = []
    spikein_contigs: List[str] = []
    for line in idxstats_text.splitlines():
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 3:
            continue
        name = fields[0]
        if name == '*':
            continue
        try:
            mapped = int(fields[2])
        except ValueError:
            continue
        if name.startswith(spikein_prefix):
            spikein += mapped
            spikein_contigs.append(name)
        else:
            primary += mapped
            primary_contigs.append(name)
    total = primary + spikein
    return {
        'primary_mapped': primary,
        'spikein_mapped': spikein,
        'spikein_fraction': (spikein / total) if total else 0.0,
        'primary_contigs': primary_contigs,
        'spikein_contigs': spikein_contigs,
    }


# ══════════════════════════════════════════════════════════════════════════════
# The netseq Namespace — built from the sub-parser's OWN defaults
# ══════════════════════════════════════════════════════════════════════════════

#: run-all flag -> `rectify netseq` dest. Named, so a new netseq knob cannot be silently dropped
#: the way `short_read` was from `correct`'s hand-built Namespace (planning 832 G-1 / 834 §6.4).
_NETSEQ_FORWARDED: Tuple[Tuple[str, str], ...] = (
    ('netseq_umi_length', 'umi_length'),
    ('netseq_dedup', 'dedup'),
    ('netseq_umi_source', 'umi_source'),
    ('netseq_min_mapq', 'min_mapq'),
    ('netseq_max_reads', 'max_reads'),
    ('netseq_output_format', 'output_format'),
    ('netseq_track_position', 'track_position'),
    ('netseq_junction_pool', 'junction_pool'),
    ('netseq_rescue_max_intronic', 'rescue_max_intronic'),
    ('netseq_rescue_min_k', 'rescue_min_k'),
    ('netseq_rescue_min_k_with_remainder', 'rescue_min_k_with_remainder'),
    ('netseq_pool_include_trna', 'pool_include_trna'),
    ('netseq_pool_include_organellar', 'pool_include_organellar'),
    ('netseq_no_deconvolution', 'no_deconvolution'),
    ('netseq_min_atract_length', 'min_atract_length'),
    ('netseq_walkback_unconditional', None),   # inverted below
    ('netseq_no_tail_detection', 'no_tail_detection'),
    ('netseq_exclude_mito', 'exclude_mito'),
)


def build_netseq_namespace(
    args,
    bam_path: Path,
    genome_path: Optional[Path],
    annotation_path: Optional[Path],
    output_dir: Path,
) -> argparse.Namespace:
    """Namespace for ``rectify netseq``, seeded from ``add_netseq_parser``'s own defaults.

    🔴 Deliberately NOT hand-built. ``stages._run_correction`` hand-builds `correct`'s Namespace
    and dropped ``short_read`` for months (planning 832 G-1 / 833 C-1), then would have dropped
    ``netseq`` (834 §6.4); ``rectify netseq --Scer`` silently had no annotation because the parser
    had no ``--annotation`` slot (836 CP0b). Three incidents, one shape. Here the sub-parser is
    instantiated and parsed, so EVERY netseq default arrives automatically and only the values
    this mode overrides are assigned.
    """
    from ..netseq_command import add_netseq_parser

    parser = argparse.ArgumentParser(prog='rectify')
    sub = parser.add_subparsers(dest='command')
    add_netseq_parser(sub)
    ns = parser.parse_args(['netseq', str(bam_path), '-o', str(output_dir)])

    ns.input = [Path(bam_path)]
    ns.output_dir = Path(output_dir)
    ns.genome = Path(genome_path) if genome_path else None
    ns.gff = Path(annotation_path) if annotation_path else None
    ns.annotation = ns.gff
    ns.organism = getattr(args, 'organism', None)

    # ── The mode's own defaults, each with a reason ────────────────────────────────────────
    # Churchman geometry: the RNA 3' end is the read 5' terminus and the gene strand is the
    # INVERSE of the BAM strand. Verified three ways (netseq_bam_processor.py:227-256).
    ns.rna3p_at = 'read5p'
    # Tracks on the CORRECTED end — the planning-836 rescue/tail/walkback output. `raw` was the
    # pre-2026-09 behaviour and discards every correction this mode exists to make.
    ns.track_position = getattr(args, 'netseq_track_position', 'corrected')
    # 🔴 Pol III INCLUDED by default. `--include-pol3=False` drops every `snRNA` GFF feature
    # (netseq_command.py:105), and the planning-829 §4 primary QC observable is the U5/SNR7 3' end
    # at chrVII:939,521 (-), 3.1 % of all reads — an snRNA. Excluding it deletes the one position
    # that proves the strand convention is right. rDNA is still excluded (--include-rdna off).
    ns.include_pol3 = not getattr(args, 'netseq_exclude_pol3', False)
    # Counts, not RPM: NNLS does not conserve mass so any counting must use raw counts
    # (829 §4: deconv 3.14 M vs 2.39 M reads), the `--netseq-dir` refiner is calibrated on raw
    # signal (834 gap H), and the validated 836 build wrote counts.
    ns.no_rpm_normalize = not getattr(args, 'netseq_rpm', False)

    for run_dest, netseq_dest in _NETSEQ_FORWARDED:
        if netseq_dest is None:
            continue
        if hasattr(args, run_dest):
            value = getattr(args, run_dest)
            if value is not None:
                setattr(ns, netseq_dest, value)
    # Inverted pair: `netseq`'s dest is walkback_requires_clip_a (default True).
    if getattr(args, 'netseq_walkback_unconditional', False):
        ns.walkback_requires_clip_a = False

    if ns.dedup and not ns.umi_length:
        raise ValueError(
            "--netseq-dedup requires --netseq-umi-length N. And verify the randomer is UNIVERSAL "
            "first (aligned 5'-clip histogram, not the kit): on a mixed library the overshoot "
            "correction shifts the randomer-free class by the full UMI length (planning 829 §4)."
        )
    return ns


# ══════════════════════════════════════════════════════════════════════════════
# Stage runners
# ══════════════════════════════════════════════════════════════════════════════

def _require(binary: str, what: str) -> str:
    path = shutil.which(binary)
    if not path:
        raise RuntimeError(
            f"{what} requires `{binary}` on PATH and it was not found. "
            f"run-all --netseq runs {binary} directly; install it or pass an explicit path."
        )
    return path


def _run(cmd: Sequence[str], log_path: Optional[Path] = None, **kwargs) -> None:
    printable = ' '.join(str(c) for c in cmd)
    print(f"    $ {printable[:400]}{'…' if len(printable) > 400 else ''}")
    if log_path is not None:
        with open(log_path, 'w') as log:
            rc = subprocess.run([str(c) for c in cmd], stdout=log, stderr=subprocess.STDOUT,
                                **kwargs).returncode
    else:
        rc = subprocess.run([str(c) for c in cmd], **kwargs).returncode
    if rc != 0:
        tail = ''
        if log_path is not None and log_path.exists():
            tail = '\n'.join(log_path.read_text(errors='replace').splitlines()[-25:])
        raise RuntimeError(f"command failed (rc={rc}): {printable}\n{tail}")


def netseq_trim(
    args,
    reads: Path,
    out_dir: Path,
    sample_id: str,
    threads: int,
) -> Path:
    """Stage 1: 3' linker trim. Returns the trimmed FASTQ path."""
    out_dir.mkdir(parents=True, exist_ok=True)
    trimmed = out_dir / f"{sample_id}.trimmed.fastq.gz"
    report = out_dir / f"{sample_id}.cutadapt.txt"
    if trimmed.exists() and trimmed.stat().st_size > 0:
        print(f"    reusing existing trimmed FASTQ: {trimmed.name}")
        return trimmed
    cutadapt = getattr(args, 'cutadapt_path', None) or _require('cutadapt', 'run-all --netseq')
    # 🔴 The temp name MUST keep the `.fastq.gz` suffix. cutadapt picks its output compression
    # from the file extension, so a `<name>.part` temp is written as PLAIN TEXT; renaming it to
    # `.gz` afterwards then makes STAR's `--readFilesCommand zcat` read **zero reads and exit 0**
    # (observed 2026-09-03, job 14648677: `Number of input reads | 0`, rc=0, empty BAM). That is a
    # green run over an empty input, which this codebase treats as worse than a crash.
    part = trimmed.with_name(trimmed.name.replace('.fastq.gz', '.part.fastq.gz'))
    cmd = build_cutadapt_cmd(
        reads=reads,
        output=part,
        adapter=getattr(args, 'netseq_adapter', CHURCHMAN_3P_LINKER),
        min_length=getattr(args, 'netseq_min_length', 18),
        nextseq_trim=getattr(args, 'netseq_nextseq_trim', 20),
        quality_cutoff=getattr(args, 'netseq_quality_cutoff', None),
        threads=threads,
        cutadapt_path=cutadapt,
    )
    _run(cmd, log_path=report)
    if not part.exists() or part.stat().st_size == 0:
        raise RuntimeError(f"cutadapt produced no output at {part}; see {report}")
    with open(part, 'rb') as fh:
        if fh.read(2) != b'\x1f\x8b':
            raise RuntimeError(
                f"{part} is named .gz but is not gzip — the downstream `zcat` would read zero "
                "reads and STAR would exit 0 with an empty BAM. Check the cutadapt version's "
                "extension handling."
            )
    n_out = _cutadapt_out_reads(report)
    if n_out is not None:
        print(f"    cutadapt kept {n_out:,} reads")
        if n_out == 0:
            raise RuntimeError(f"cutadapt kept 0 reads; see {report}")
    part.replace(trimmed)
    print(f"    trimmed -> {trimmed.name}   (report: {report.name})")
    return trimmed


def _cutadapt_out_reads(report: Path) -> Optional[int]:
    """``out_reads`` from a cutadapt ``--report=minimal`` table, or None if unparseable."""
    try:
        lines = [ln for ln in report.read_text().splitlines() if ln.strip()]
        header, values = lines[0].split('\t'), lines[1].split('\t')
        return int(values[header.index('out_reads')])
    except Exception:
        return None


def _star_input_reads(log_final: Path) -> Optional[int]:
    """``Number of input reads`` from STAR's Log.final.out, or None."""
    try:
        for line in log_final.read_text().splitlines():
            if 'Number of input reads' in line:
                return int(line.split('|')[1].strip())
    except Exception:
        pass
    return None


def ensure_star_index(
    args,
    genome_path: Path,
    work_dir: Path,
    threads: int,
) -> Tuple[Path, Optional[Path]]:
    """Stage 2a: resolve (or build) the STAR index. Returns ``(index_dir, combined_fasta)``.

    An explicit ``--star-index`` is used as-is and never rebuilt. Otherwise the index lives at
    ``--star-index-dir`` (default ``<work>/star_index``) and is built once and REUSED — the
    ``SAindex`` file is the completion sentinel, so a half-built index from a killed job is
    rebuilt rather than silently used.
    """
    explicit = getattr(args, 'star_index', None)
    if explicit:
        index_dir = Path(explicit)
        if not (index_dir / 'SAindex').exists():
            raise FileNotFoundError(f"--star-index {index_dir} has no SAindex")
        print(f"    using prebuilt STAR index: {index_dir}")
        return index_dir, None

    spikein = getattr(args, 'spikein_fasta', None)
    prefix = getattr(args, 'spikein_prefix', 'Sp_')
    index_dir = Path(getattr(args, 'star_index_dir', None) or (work_dir / 'star_index'))
    combined: Optional[Path] = None
    if spikein:
        combined = index_dir.parent / (
            f"{Path(genome_path).stem}__{Path(spikein).stem}.fa")
        combined.parent.mkdir(parents=True, exist_ok=True)
        if not combined.exists():
            counts = build_combined_genome(Path(genome_path), Path(spikein), combined, prefix)
            print(f"    combined genome: {counts['primary_contigs']} primary + "
                  f"{counts['spikein_contigs']} spike-in contigs (prefix {prefix!r}) -> "
                  f"{combined.name}")
        genome_for_index: Path = combined
    else:
        genome_for_index = Path(genome_path)

    if (index_dir / 'SAindex').exists():
        print(f"    reusing STAR index: {index_dir}")
        return index_dir, combined

    star = getattr(args, 'star_path', None) or _require('STAR', 'run-all --netseq')
    index_dir.mkdir(parents=True, exist_ok=True)
    print(f"    building STAR index -> {index_dir}")
    _run(build_star_index_cmd(
        genome_for_index, index_dir, threads=threads,
        sa_index_nbases=getattr(args, 'star_sa_index_nbases', 11), star_path=star,
    ), log_path=index_dir / 'genomeGenerate.log')
    if not (index_dir / 'SAindex').exists():
        raise RuntimeError(f"STAR genomeGenerate wrote no SAindex in {index_dir}")
    return index_dir, combined


def netseq_align(
    args,
    trimmed_fastq: Path,
    index_dir: Path,
    out_dir: Path,
    sample_id: str,
    threads: int,
) -> Tuple[Path, Dict[str, object]]:
    """Stage 2b: STAR + sort + index + idxstats + spike-in split.

    Returns ``(primary_genome_bam, spikein_stats)``. The BAM handed to stage 3 carries the
    PRIMARY genome contigs only: spike-in contigs are dropped so they cannot enter the tracks,
    the exclusion regions or the deconvolution, and they are accounted for separately in
    ``<sample>.spikein.tsv``.
    """
    import pysam

    out_dir.mkdir(parents=True, exist_ok=True)
    all_bam = out_dir / f"{sample_id}.all.sorted.bam"
    prefix = getattr(args, 'spikein_prefix', 'Sp_')
    star = getattr(args, 'star_path', None) or _require('STAR', 'run-all --netseq')

    if not (all_bam.exists() and all_bam.stat().st_size > 0):
        star_prefix = out_dir / f"{sample_id}."
        cmd = build_netseq_star_cmd(
            reads_in=str(trimmed_fastq), index_dir=index_dir, out_prefix=star_prefix,
            threads=threads, gzipped=str(trimmed_fastq).endswith('.gz'), star_path=star,
        )
        _run(cmd, log_path=out_dir / f"{sample_id}.star.log")
        # 🔴 STAR exits 0 on an unreadable/empty input, reporting `Number of input reads | 0`.
        # Check the INPUT count, not just the output BAM: an empty BAM from zero mapped reads and
        # an empty BAM from zero input reads are the same file but completely different bugs.
        n_in = _star_input_reads(out_dir / f"{sample_id}.Log.final.out")
        if n_in is not None:
            print(f"    STAR input reads: {n_in:,}")
            if n_in == 0:
                raise RuntimeError(
                    f"STAR read ZERO input reads from {trimmed_fastq} and still exited 0. The "
                    "usual cause is a file named .gz that is not gzip (--readFilesCommand zcat "
                    f"then yields nothing). See {out_dir / (sample_id + '.Log.final.out')}."
                )
        unsorted = out_dir / f"{sample_id}.Aligned.out.bam"
        if not unsorted.exists() or unsorted.stat().st_size == 0:
            raise RuntimeError(f"STAR wrote no alignments at {unsorted}")
        pysam.sort('-@', str(max(1, threads - 1)), '-o', str(all_bam), str(unsorted))
        unsorted.unlink()
    if not Path(str(all_bam) + '.bai').exists():
        pysam.index(str(all_bam))

    idxstats_text = pysam.idxstats(str(all_bam))
    (out_dir / f"{sample_id}.idxstats.tsv").write_text(idxstats_text)
    stats = parse_idxstats_spikein(idxstats_text, prefix)

    spikein_tsv = out_dir / f"{sample_id}.spikein.tsv"
    spikein_tsv.write_text(
        "sample\tprimary_mapped\tspikein_mapped\tspikein_fraction\tspikein_prefix\n"
        f"{sample_id}\t{stats['primary_mapped']}\t{stats['spikein_mapped']}\t"
        f"{stats['spikein_fraction']:.6f}\t{prefix}\n"
    )
    print(f"    mapped: primary {stats['primary_mapped']:,}   "
          f"spike-in {stats['spikein_mapped']:,}   "
          f"fraction {stats['spikein_fraction']:.4f}   -> {spikein_tsv.name}")

    if int(stats['primary_mapped']) == 0:
        raise RuntimeError(
            f"no reads mapped to the primary genome in {all_bam} — refusing to write a green "
            "sentinel over an empty run"
        )

    if not stats['spikein_contigs']:
        return all_bam, stats

    primary_bam = out_dir / f"{sample_id}.bam"
    if not (primary_bam.exists() and primary_bam.stat().st_size > 0):
        pysam.view('-@', str(max(1, threads - 1)), '-b', '-o', str(primary_bam), str(all_bam),
                   *stats['primary_contigs'], catch_stdout=False)
    if not Path(str(primary_bam) + '.bai').exists():
        pysam.index(str(primary_bam))
    print(f"    primary-genome BAM (spike-in contigs dropped): {primary_bam.name}")
    return primary_bam, stats


def run_netseq_pipeline(
    args,
    input_path: Path,
    input_type: str,
    output_dir: Path,
    work_dir: Path,
    sample_id: str,
    genome_path: Optional[Path],
    annotation_path: Optional[Path],
    threads: int = 4,
) -> int:
    """Run the whole ``run-all --netseq`` mode. Returns a process exit code."""
    from ..netseq_command import run_netseq

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(work_dir)
    summary: Dict[str, object] = {'sample': sample_id, 'mode': 'run-all --netseq'}

    print("\n" + "=" * 70)
    print("RECTIFY run-all --netseq  (NET-seq / mNET-seq nascent 3'-end pipeline)")
    print("=" * 70)
    print("  Stages: [1] cutadapt 3' linker trim  ->  [2] STAR (absolute match floors)  ->  "
          "[3] rectify netseq")
    print("  NOTE: there is deliberately NO `correct`/consensus stage. The 1-10 nt soft-clip "
          "re-placement\n        NET-seq needs is the donor-side junction rescue inside `rectify "
          "netseq` (planning 836);\n        the COMPASS panel adds no placement information on "
          "35-nt reads. See --help.")
    print(f"  NOTE: Pol III / snRNA loci are INCLUDED by default "
          f"({'ON' if not getattr(args, 'netseq_exclude_pol3', False) else 'OFF'}) — excluding "
          "them drops\n        every snRNA feature, and the U5/SNR7 3' end is the QC observable "
          "that proves the\n        strand convention. rDNA is still excluded. "
          "--netseq-exclude-pol3 to change.")

    if input_type in ('fastq', 'fastq.gz'):
        if getattr(args, 'netseq_skip_trim', False):
            print("\n[Step 1/3] Trimming SKIPPED (--netseq-skip-trim): input is already "
                  "linker-trimmed.")
            trimmed = Path(input_path)
        else:
            print("\n[Step 1/3] Trimming (cutadapt; randomer LEFT IN PLACE)...")
            print("-" * 50)
            trimmed = netseq_trim(args, input_path, work_dir / 'netseq_trim', sample_id, threads)

        print("\n[Step 2/3] Aligning (STAR, planning-829 §4 absolute match floors)...")
        print("-" * 50)
        if not genome_path:
            raise ValueError("run-all --netseq needs --genome (or --Scer/--organism) to align")
        index_dir, combined = ensure_star_index(args, Path(genome_path), work_dir, threads)
        bam_path, spikein = netseq_align(
            args, trimmed, index_dir, work_dir / 'netseq_align', sample_id, threads)
        summary['spikein'] = spikein
        summary['star_index'] = str(index_dir)
        if combined:
            summary['combined_genome'] = str(combined)
    else:
        print("\n[Steps 1-2/3] Skipped — input is already an aligned BAM.")
        bam_path = Path(input_path)

    print("\n[Step 3/3] rectify netseq (donor-side junction rescue, corrected-end tracks)...")
    print("-" * 50)
    netseq_args = build_netseq_namespace(
        args, bam_path, genome_path, annotation_path, output_dir)
    rc = run_netseq(netseq_args)
    if rc != 0:
        raise RuntimeError(f"rectify netseq failed (rc={rc})")

    summary['bam'] = str(bam_path)
    summary['netseq_flags'] = {
        'rna3p_at': netseq_args.rna3p_at,
        'track_position': netseq_args.track_position,
        'include_pol3': netseq_args.include_pol3,
        'no_rpm_normalize': netseq_args.no_rpm_normalize,
        'umi_length': netseq_args.umi_length,
        'dedup': netseq_args.dedup,
        'output_format': list(netseq_args.output_format),
        'walkback_requires_clip_a': netseq_args.walkback_requires_clip_a,
        'rescue_min_k': netseq_args.rescue_min_k,
        'rescue_min_k_with_remainder': netseq_args.rescue_min_k_with_remainder,
        'pool_include_trna': netseq_args.pool_include_trna,
        'pool_include_organellar': netseq_args.pool_include_organellar,
    }
    (output_dir / f"{sample_id}.runall_netseq.json").write_text(json.dumps(summary, indent=1))

    print("\n" + "=" * 70)
    print("run-all --netseq complete")
    print(f"  outputs: {output_dir}")
    print("=" * 70)
    return 0
