"""Regression tests for two `rectify split --generate-slurm` generator bugs.

BUG#2 (blank-line continuation): the per-chunk `rectify correct` array script
interpolated optional args (``--junction-penalty-table`` / ``--str-penalty-table``)
as their own template lines. When those tables are absent (the common case) the
slots render as *blank lines* immediately after ``--annotation "..." \\``. In
shell, a blank line after a ``\\``-continued line terminates the logical command,
so every argument below it — including the trailing ``-o``/``--output`` — is
orphaned, and every correct task fails with "error: -o/--output required".
Note ``bash -n`` does NOT catch this (the orphaned lines are syntactically valid
separate commands), so these tests assert the *runtime* invariant: after
reconstructing bash's line-continuation join, the ``rectify correct`` logical
command still carries ``-o`` and all its arguments.

BUG#1 (minimap2 ``-y`` aux leak): ``rectify align`` runs minimap2 ``-y``, which
copies the FASTQ comment verbatim into the SAM aux column. An ENA/uuid read
description (e.g. ``<uuid>/1``) or a Dorado ``runid=... ch=...`` comment is not a
valid SAM aux field, so htslib rejects it and samtools sort aborts (SIGPIPE 141).
``rectify split`` sanitizes the comment to well-formed aux tokens at chunk-write
time. These tests assert every split-written chunk header carries only valid SAM
aux tokens, for all input types that route through ``split_fastq`` /
``split_fastq_paired``.
"""
import argparse
import glob
import gzip
import os
from pathlib import Path

import pytest

from rectify.core.chunking.sidecar import (
    _SAM_AUX_TOKEN_RE,
    format_fastq_header_with_rn,
    sanitize_fastq_comment_for_aux,
    split_fastq_header,
)
from rectify.core.commands import split_command as sc
from rectify.core.commands.split_command import (
    _correct_array_body,
    _strip_blank_continuation_lines,
    split_fastq,
    split_fastq_paired,
)


# ── helpers ──────────────────────────────────────────────────────────────────

def _bash_logical_lines(script: str):
    """Reconstruct bash's line-continuation join: a trailing ``\\`` splices the
    next physical line onto the current logical line; a blank line's own
    (unescaped) newline ends the logical command."""
    logical, buf = [], ''
    for line in script.split('\n'):
        if line.rstrip().endswith('\\'):
            buf += line.rstrip()[:-1]  # drop the escaping backslash, keep going
        else:
            logical.append(buf + line)
            buf = ''
    if buf:
        logical.append(buf)
    return logical


def _correct_logical_command(script: str) -> str:
    return next(l for l in _bash_logical_lines(script) if 'rectify correct' in l)


def _count_blank_after_continuation(script: str) -> int:
    lines = script.split('\n')
    return sum(
        1
        for i in range(1, len(lines))
        if lines[i].strip() == '' and lines[i - 1].rstrip().endswith('\\')
    )


def _make_correct_body(jpt='', spt='', wcb=False, aligners=('gapmm2', 'uLTRA', 'deSALT')):
    return _correct_array_body(
        aligner=aligners[0], n_chunks=14, all_aligners=list(aligners),
        sample_prefix='S', output_dir=Path('/out'),
        genome='/g.fa', annot='/a.gff', python_path='/py', rectify_src='/src',
        junction_penalty_table=jpt, str_penalty_table=spt,
        write_per_aligner_corrected_bams=wcb,
    )


# ── BUG#2: blank-line continuation in the correct-array emitter ───────────────

@pytest.mark.parametrize('jpt', ['', '/pen/junc.tsv'])
@pytest.mark.parametrize('spt', ['', '/pen/str.tsv'])
@pytest.mark.parametrize('wcb', [False, True])
def test_correct_array_command_survives_line_continuation(jpt, spt, wcb):
    """The rendered ``rectify correct`` command keeps every arg after joining
    continuations — most importantly ``-o`` (the arg BUG#2 dropped)."""
    body = _make_correct_body(jpt=jpt, spt=spt, wcb=wcb)

    # No blank line may follow a backslash-continued line (would truncate).
    assert _count_blank_after_continuation(body) == 0

    cmd = _correct_logical_command(body)
    required = [
        '-o "$CHUNK_OUTDIR/corrected_reads.tsv"',
        '--junction-pool-cache "$POOL_PKL"',
        '--aligner-bams "gapmm2:',
        '--aligner-bams "uLTRA:',
        '--aligner-bams "deSALT:',
        '--threads "$CORRECT_CPUS"',
        '--streaming',
        '--checkpoint-dir "$CHECKPOINT_DIR"',
    ]
    if jpt:
        required.append('--junction-penalty-table "/pen/junc.tsv"')
    if spt:
        required.append('--str-penalty-table "/pen/str.tsv"')
    if wcb:
        required.append('--write-corrected-bam')
    for token in required:
        assert token in cmd, f'{token!r} missing from correct command (jpt={jpt!r} spt={spt!r} wcb={wcb})'


def test_correct_array_omits_penalty_flags_when_absent():
    """The optional penalty-table flags must NOT appear when no table is given
    (the case that previously rendered blank lines)."""
    cmd = _correct_logical_command(_make_correct_body(jpt='', spt=''))
    assert '--junction-penalty-table' not in cmd
    assert '--str-penalty-table' not in cmd


@pytest.mark.parametrize('scheduler', ['slurm', 'uge'])
def test_generate_scripts_emits_no_truncating_blanks(tmp_path, scheduler):
    """End-to-end: every script written by ``_generate_scripts`` (which is where
    the ``_write`` belt-and-suspenders scrubber runs) is free of truncating
    blanks, and the correct-array command still carries ``-o`` — for both
    schedulers."""
    args = argparse.Namespace(
        output_dir=tmp_path, genome=None, annotation=None,
        scheduler=scheduler, python_path='/py', rectify_src='/src',
        other_aligners=['gapmm2', 'uLTRA', 'deSALT'], skip_map_pacbio=True,
        slurm_partition=None, slurm_account=None,
        uge_queue=None, uge_pe='smp', pbs_queue='workq',
    )
    sc._generate_scripts(args, n_chunks=4, sample_prefix='S')

    scripts = sorted(glob.glob(str(tmp_path / '*.sh')))
    assert scripts, 'no scripts generated'
    for path in scripts:
        text = Path(path).read_text()
        assert _count_blank_after_continuation(text) == 0, (
            f'{os.path.basename(path)} has a blank line after a continuation'
        )

    for path in glob.glob(str(tmp_path / 'run_array_correct_*.sh')):
        cmd = _correct_logical_command(Path(path).read_text())
        assert '-o "$CHUNK_OUTDIR/corrected_reads.tsv"' in cmd
        assert '--junction-pool-cache "$POOL_PKL"' in cmd


# ── BUG#2 Layer 2: the standalone continuation scrubber ───────────────────────

def test_strip_blank_continuation_removes_truncating_blank():
    buggy = (
        'rectify correct \\\n'
        '    --genome g \\\n'
        '\n'            # <- truncating blank (an empty optional-arg slot)
        '    -o out\n'
    )
    fixed = _strip_blank_continuation_lines(buggy)
    assert _count_blank_after_continuation(fixed) == 0
    assert '-o out' in _bash_logical_lines(fixed)[0]


def test_strip_blank_continuation_handles_consecutive_blanks():
    buggy = (
        'cmd \\\n'
        '    --a x \\\n'
        '\n'            # both jpt and spt slots empty
        '\n'
        '    --b y \\\n'
        '    -o out\n'
    )
    fixed = _strip_blank_continuation_lines(buggy)
    joined = _bash_logical_lines(fixed)[0]
    for tok in ('--a x', '--b y', '-o out'):
        assert tok in joined


def test_strip_blank_continuation_preserves_legit_blanks():
    """A blank line that does NOT follow a continuation is a real block
    separator and must be kept."""
    clean = 'echo one\n\necho two \\\n    --flag\n'
    assert _strip_blank_continuation_lines(clean) == clean


# ── BUG#1: split-write sanitizes FASTQ comments to valid SAM aux ──────────────

def test_sanitize_fastq_comment_for_aux_drops_invalid_tokens():
    # ENA/uuid read description + Dorado run metadata are NOT valid aux → dropped;
    # genuine SAM aux tokens are kept (tab-joined).
    assert sanitize_fastq_comment_for_aux('<uuid>/1') == ''
    assert sanitize_fastq_comment_for_aux('runid=abc ch=5 start_time=x') == ''
    assert sanitize_fastq_comment_for_aux('XA:Z:foo NM:i:3') == 'XA:Z:foo\tNM:i:3'
    assert sanitize_fastq_comment_for_aux('<uuid>/1 XC:i:42') == 'XC:i:42'
    assert sanitize_fastq_comment_for_aux('') == ''


def test_format_fastq_header_with_rn_sanitizes_comment():
    header = format_fastq_header_with_rn('read_uuid', '<uuid>/1 XA:Z:tail', 7)
    # RN tag is always present; the invalid <uuid>/1 token is stripped.
    assert 'RN:i:7' in header
    assert '<uuid>/1' not in header
    assert 'XA:Z:tail' in header
    # Every whitespace-delimited comment token is a valid SAM aux field.
    _, comment = split_fastq_header(header)
    for tok in comment.split():
        assert _SAM_AUX_TOKEN_RE.match(tok), f'{tok!r} is not a valid SAM aux token'


def _write_fastq(path: Path, records):
    with gzip.open(path, 'wt') as fh:
        for name, comment, seq in records:
            header = f'@{name} {comment}' if comment else f'@{name}'
            fh.write(f'{header}\n{seq}\n+\n{"I" * len(seq)}\n')


def _assert_chunk_headers_all_valid_aux(chunk_paths):
    seen = 0
    for p in chunk_paths:
        with gzip.open(p, 'rt') as fh:
            for i, line in enumerate(fh):
                if i % 4 != 0:
                    continue
                seen += 1
                _, comment = split_fastq_header(line)
                for tok in comment.split():
                    assert _SAM_AUX_TOKEN_RE.match(tok), (
                        f'{tok!r} in {p} is not a valid SAM aux token'
                    )
    assert seen > 0


def test_split_fastq_single_end_sanitizes_ena_and_dorado_comments(tmp_path):
    """DRS + single-end short-read (QuantSeq) both route through ``split_fastq``;
    every chunk header must carry only valid SAM aux tokens."""
    src = tmp_path / 'reads.fastq.gz'
    _write_fastq(src, [
        ('ce7f9uuid-0001/1', '<uuid>/1', 'ACGTACGT'),
        ('read2', 'runid=abc ch=5 start_time=2024', 'TTTTGGGG'),
        ('read3', 'XA:Z:polyA NM:i:2', 'CCCCAAAA'),
        ('read4', '', 'GGGGCCCC'),
    ])
    chunks = split_fastq(src, tmp_path / 'chunks', n_chunks=2, prefix='s')
    _assert_chunk_headers_all_valid_aux(chunks)


def test_split_fastq_paired_sanitizes_comments(tmp_path):
    """Paired short-read (TruSeq) routes through ``split_fastq_paired``; both
    mates' chunk headers must carry only valid SAM aux tokens."""
    r1 = tmp_path / 'r1.fastq.gz'
    r2 = tmp_path / 'r2.fastq.gz'
    _write_fastq(r1, [
        ('pairuuid-1', '1:N:0:ACGT <uuid>/1', 'ACGTACGT'),
        ('pairuuid-2', 'runid=x ch=1', 'TTTTGGGG'),
    ])
    _write_fastq(r2, [
        ('pairuuid-1', '2:N:0:ACGT <uuid>/2', 'ACGTACGT'),
        ('pairuuid-2', 'runid=x ch=1', 'TTTTGGGG'),
    ])
    r1_chunks, r2_chunks = split_fastq_paired(r1, r2, tmp_path / 'chunks', n_chunks=2, prefix='s')
    _assert_chunk_headers_all_valid_aux(list(r1_chunks) + list(r2_chunks))
