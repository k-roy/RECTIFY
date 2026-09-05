"""``rectify run-all --short-read`` — the Illumina adapter-trim stage.

🔴 **Why this stage exists.** Before it, ``run-all --short-read`` had no trimming at all, and the
adapter was *invisible in the output BAM*: HISAT2 (97 % of consensus winners on the planning-861
TruSeq smoke) runs ``--no-softclip`` and STAR runs ``EndToEnd`` on that path, so an adapter-bearing
read is turned into mismatches or an unaligned read rather than a soft clip. Measured on a real
library (planning 861 §S0c, Peterson RNA-seq wt_rep1, 200 k reads):

===========================================  ==========================================
reads carrying the R1 TruSeq adapter          **13.0 %** (25,911 / 200,000)
of which full-length dimers (whole read)      **11.5 %** (23,016 — trimmed length 75 → 0)
records showing ANY terminal clip in the BAM  **0.01 %**
===========================================  ==========================================

QuantSeq-class libraries (``--short-read --dT-primed-cDNA``) have the same exposure with a
different adapter, and there the trimmed 3' end *is* the measured RNA 3' end, so untrimmed adapter
moves called cleavage sites rather than merely losing reads.

Everything a unit test needs is a **pure command builder** (:func:`build_short_read_cutadapt_cmd`)
that never touches the filesystem or requires cutadapt on ``PATH``.

The defaults are the Illumina TruSeq / Nextera universal stems, which are also the prefix of the
NEBNext, Lexogen QuantSeq and SMARTer adapters — trimming the stem removes the whole 3' remainder:

* R1 ``AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`` (TruSeq read-1 / universal)
* R2 ``AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`` (TruSeq read-2)

⚠️ **This is NOT the NET-seq trim.** ``run-all --netseq`` has its own stage with the Churchman
linker and the "no ``-u``, the randomer stays in the read" rule
(:mod:`rectify.core.commands.run.netseq_pipeline`); the two must not be merged.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

__all__ = [
    'TRUSEQ_R1_ADAPTER',
    'TRUSEQ_R2_ADAPTER',
    'build_short_read_cutadapt_cmd',
    'parse_cutadapt_report',
    'run_short_read_trim',
]

#: Illumina TruSeq read-1 / universal adapter stem (also the prefix of NEBNext and Lexogen QuantSeq).
TRUSEQ_R1_ADAPTER = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
#: Illumina TruSeq read-2 adapter stem.
TRUSEQ_R2_ADAPTER = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'


def build_short_read_cutadapt_cmd(
    reads: Path,
    output: Path,
    *,
    read2: Optional[Path] = None,
    output2: Optional[Path] = None,
    adapter: str = TRUSEQ_R1_ADAPTER,
    adapter2: str = TRUSEQ_R2_ADAPTER,
    min_length: int = 20,
    quality_cutoff: Optional[int] = None,
    nextseq_trim: Optional[int] = None,
    threads: int = 4,
    cutadapt_path: str = 'cutadapt',
    report_path: Optional[Path] = None,
) -> List[str]:
    """Build the cutadapt command for an Illumina short-read library.

    Single-end uses ``-a``; paired-end adds ``-A`` and ``-p`` so both mates are trimmed in one pass
    and stay synchronised (trimming the mates separately desynchronises the FASTQs whenever one
    mate falls below ``-m`` — the failure mode that makes a PE aligner exit non-zero mid-run).

    ``min_length`` defaults to 20 rather than cutadapt's 0: a 0-length record is written as an empty
    read, which bbmap and magicblast reject with an error and STAR silently drops. On a library with
    11.5 % full-length dimers that is not a corner case.

    ``nextseq_trim`` selects cutadapt's two-colour-aware quality trim (a no-signal base reads as a
    high-quality G on NextSeq/NovaSeq, so plain ``-q`` does not remove the poly-G tail); it is
    mutually exclusive with ``quality_cutoff`` and the caller must not pass both.

    Returns the argv list; raises ``ValueError`` on an incoherent paired-end request.
    """
    if (read2 is None) != (output2 is None):
        raise ValueError('paired-end trimming needs both read2 and output2')
    if quality_cutoff is not None and nextseq_trim is not None:
        raise ValueError('pass either quality_cutoff or nextseq_trim, not both')
    cmd: List[str] = [cutadapt_path, '-j', str(threads), '-a', adapter]
    if read2 is not None:
        cmd += ['-A', adapter2]
    if quality_cutoff is not None:
        cmd += ['-q', str(quality_cutoff)]
    elif nextseq_trim is not None:
        cmd += [f'--nextseq-trim={nextseq_trim}']
    cmd += ['-m', str(min_length), '--report=minimal', '-o', str(output)]
    if read2 is not None:
        cmd += ['-p', str(output2)]
    cmd += [str(reads)]
    if read2 is not None:
        cmd += [str(read2)]
    return cmd


def parse_cutadapt_report(text: str) -> dict:
    """Parse ``cutadapt --report=minimal`` output into a dict of ints/floats.

    The minimal report is two tab-separated lines (header, values). Unknown or unparsable input
    returns ``{}`` rather than raising — a trim whose *report* could not be read still produced a
    valid FASTQ, and the caller logs the miss instead of failing the run.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        return {}
    keys = lines[0].split('\t')
    vals = lines[1].split('\t')
    if len(keys) != len(vals):
        return {}
    out = {}
    for k, v in zip(keys, vals):
        try:
            out[k] = int(v)
        except ValueError:
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = v
    return out


def run_short_read_trim(
    input_path: Path,
    work_dir: Path,
    sample_id: str,
    *,
    read2: Optional[Path] = None,
    adapter: str = TRUSEQ_R1_ADAPTER,
    adapter2: str = TRUSEQ_R2_ADAPTER,
    min_length: int = 20,
    quality_cutoff: Optional[int] = None,
    nextseq_trim: Optional[int] = None,
    threads: int = 4,
    cutadapt_path: Optional[str] = None,
    log=None,
) -> Tuple[Path, Optional[Path], dict]:
    """Trim an Illumina library and return ``(trimmed_r1, trimmed_r2_or_None, report_dict)``.

    Writes into ``work_dir/short_read_trim/``. Raises ``FileNotFoundError`` when cutadapt is not on
    ``PATH`` and no explicit path was given, and ``RuntimeError`` when cutadapt exits non-zero — a
    silently-skipped trim is exactly the failure this stage exists to prevent, so it is never
    swallowed here; the caller decides whether to abort or continue.
    """
    exe = cutadapt_path or shutil.which('cutadapt')
    if not exe:
        raise FileNotFoundError(
            'cutadapt not found on PATH. Install it, pass --cutadapt-path, or pass '
            '--skip-short-read-trim if the input FASTQ is already adapter-trimmed.'
        )
    out_dir = Path(work_dir) / 'short_read_trim'
    out_dir.mkdir(parents=True, exist_ok=True)
    r1_out = out_dir / f'{sample_id}_R1.trimmed.fastq.gz'
    r2_out = out_dir / f'{sample_id}_R2.trimmed.fastq.gz' if read2 is not None else None
    cmd = build_short_read_cutadapt_cmd(
        Path(input_path), r1_out, read2=Path(read2) if read2 is not None else None,
        output2=r2_out, adapter=adapter, adapter2=adapter2, min_length=min_length,
        quality_cutoff=quality_cutoff, nextseq_trim=nextseq_trim, threads=threads,
        cutadapt_path=exe,
    )
    if log is not None:
        log.write('short-read trim: ' + ' '.join(cmd) + '\n')
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f'cutadapt failed (exit {proc.returncode}): {proc.stderr.strip()[:500]}'
        )
    report = parse_cutadapt_report(proc.stdout)
    (out_dir / f'{sample_id}.cutadapt.txt').write_text(proc.stdout + '\n' + proc.stderr)
    if log is not None and report:
        log.write(f'short-read trim report: {report}\n')
    n_in = report.get('in_reads')
    n_out = report.get('out_reads')
    if n_in and n_out is not None and n_out < 0.5 * n_in:
        msg = (f'[{sample_id}] short-read trim kept only {n_out:,}/{n_in:,} reads '
               f'({100.0 * n_out / n_in:.1f} %) — check --short-read-adapter and '
               f'--short-read-min-length')
        print(f'WARNING: {msg}', file=sys.stderr)
        if log is not None:
            log.write(f'WARNING: {msg}\n')
    return r1_out, r2_out, report
