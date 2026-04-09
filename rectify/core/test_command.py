#!/usr/bin/env python3
"""
rectify test — Installation smoke-test and integration test.

Quick mode (default): runs `rectify correct` on the bundled validation BAM
and verifies all six correction categories.  No external aligner tools needed.

Full mode (--full): starts from the bundled FASTQ, runs `rectify align` to
generate per-aligner BAMs, then runs correction.  Requires minimap2 (and
optionally mapPacBio / gapmm2) on PATH.

Author: Kevin R. Roy
"""

import argparse
import csv
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Terminal colours (degrade gracefully)
# ---------------------------------------------------------------------------

def _supports_color() -> bool:
    return hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()


def _green(s: str) -> str:
    return f'\033[32m{s}\033[0m' if _supports_color() else s


def _red(s: str) -> str:
    return f'\033[31m{s}\033[0m' if _supports_color() else s


def _yellow(s: str) -> str:
    return f'\033[33m{s}\033[0m' if _supports_color() else s


def _bold(s: str) -> str:
    return f'\033[1m{s}\033[0m' if _supports_color() else s


PASS = _green('PASS')
FAIL = _red('FAIL')
SKIP = _yellow('SKIP')

# ---------------------------------------------------------------------------
# Bundled data helpers
# ---------------------------------------------------------------------------

def _data_dir() -> Path:
    return Path(__file__).parent.parent / 'data'


def _bundled_bam() -> Path:
    return _data_dir() / 'validation_reads.bam'


def _bundled_fastq() -> Path:
    return _data_dir() / 'validation_reads.fastq.gz'


def _get_genome() -> Optional[Path]:
    try:
        from rectify.data import get_bundled_genome_path
        return get_bundled_genome_path('saccharomyces_cerevisiae')
    except Exception:
        return None


def _get_annotation() -> Optional[Path]:
    try:
        from rectify.data import get_bundled_annotation_path
        return get_bundled_annotation_path('saccharomyces_cerevisiae')
    except Exception:
        return None

# ---------------------------------------------------------------------------
# Subprocess helpers
# ---------------------------------------------------------------------------

def _run(cmd: List[str], label: str) -> subprocess.CompletedProcess:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'\n  [{FAIL}] {label}')
        print(f'  Command: {" ".join(cmd)}')
        print(f'  stderr:\n{result.stderr[-2000:]}')
    return result


def _rectify(*args) -> List[str]:
    return [sys.executable, '-m', 'rectify.cli'] + list(args)

# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------

def _load_tsv_first_row(tsv: Path) -> Dict[str, Dict]:
    """Return first output row per read_id."""
    rows: Dict[str, Dict] = {}
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rid = row.get('read_id', '')
            if rid not in rows:
                rows[rid] = row
    return rows


def _load_tsv_all_rows(tsv: Path) -> Dict[str, List[Dict]]:
    """Return all output rows per read_id (for Cat6 multi-peak)."""
    rows: Dict[str, List[Dict]] = {}
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rid = row.get('read_id', '')
            rows.setdefault(rid, []).append(row)
    return rows


def _read_names_by_label(bam: Path) -> Dict[str, str]:
    """Map XV label → read_id from a BAM."""
    try:
        import pysam
    except ImportError:
        return {}
    result = {}
    with pysam.AlignmentFile(str(bam), 'rb') as f:
        for r in f:
            if r.is_secondary or r.is_supplementary:
                continue
            try:
                result[r.get_tag('XV')] = r.query_name
            except KeyError:
                pass
    return result

# ---------------------------------------------------------------------------
# Check helpers
# ---------------------------------------------------------------------------

class _Results:
    def __init__(self):
        self.passed: int = 0
        self.failed: int = 0
        self.skipped: int = 0

    def check(self, label: str, ok: bool, msg: str = '', skip: bool = False,
              detail: str = '') -> bool:
        indent = '    '
        if skip:
            self.skipped += 1
            print(f'{indent}[{SKIP}] {label}')
            return False
        if ok:
            self.passed += 1
            suffix = f'  {detail}' if detail else ''
            print(f'{indent}[{PASS}] {label}{suffix}')
        else:
            self.failed += 1
            print(f'{indent}[{FAIL}] {label}  ← {msg}')
        return ok

# ---------------------------------------------------------------------------
# Per-category checks
# ---------------------------------------------------------------------------

def _check_cat1(res: _Results, label_map: Dict, corrected: Dict) -> None:
    print('\n  Category 1 — poly-A walkback (A-tract ambiguity):')
    for label, strand in [
        ('cat1_plus_1', '+'), ('cat1_plus_2', '+'),
        ('cat1_minus_1', '-'), ('cat1_minus_2', '-'),
    ]:
        rid = label_map.get(label)
        row = corrected.get(rid) if rid else None
        if row is None:
            res.check(label, False, skip=True, msg='')
            continue
        orig = int(row['original_3prime'])
        corr = int(row['corrected_3prime'])
        shifted = orig != corr
        direction_ok = (corr < orig) if strand == '+' else (corr > orig)
        shift = corr - orig
        ok = shifted and direction_ok
        res.check(label, ok,
                  msg=f'corrected={corr} original={orig} (no shift or wrong direction)',
                  detail=f'{shift:+d} bp')


def _check_cat2(res: _Results, label_map: Dict, corrected: Dict) -> None:
    print('\n  Category 2 — soft-clip rescue:')
    for label, strand in [
        ('cat2_plus_1', '+'), ('cat2_plus_2', '+'),
        ('cat2_minus_1', '-'), ('cat2_minus_2', '-'),
    ]:
        rid = label_map.get(label)
        row = corrected.get(rid) if rid else None
        if row is None:
            res.check(label, False, skip=True, msg='')
            continue
        orig = int(row['original_3prime'])
        corr = int(row['corrected_3prime'])
        shift = abs(corr - orig)
        ok = 1 <= shift <= 20
        res.check(label, ok,
                  msg=f'shift={shift} bp (expected 1–20)',
                  detail=f'{corr - orig:+d} bp')


def _check_cat3(res: _Results, label_map: Dict, corrected: Dict) -> None:
    print('\n  Category 3 — 5\u2019 junction rescue:')
    raw5 = {
        'cat3_plus_1': 125270,
        'cat3_plus_2': 168808,
        'cat3_minus_1': 900766,
        'cat3_minus_2': 60192,
    }
    for label, strand in [
        ('cat3_plus_1', '+'), ('cat3_plus_2', '+'),
        ('cat3_minus_1', '-'), ('cat3_minus_2', '-'),
    ]:
        rid = label_map.get(label)
        row = corrected.get(rid) if rid else None
        if row is None:
            res.check(label, False, skip=True, msg='')
            continue
        five = row.get('five_prime_position', '')
        if five in ('', 'NA', 'None', '-1', 'nan'):
            res.check(label, False, msg='five_prime_position not set (annotation missing?)')
            continue
        rescued = int(float(five))
        expected_raw = raw5[label]
        ok = rescued != expected_raw
        res.check(label, ok,
                  msg=f'five_prime_position={rescued} unchanged from raw {expected_raw}',
                  detail=f'rescued to {rescued}')


def _check_cat4(res: _Results, label_map: Dict, corrected: Dict) -> None:
    print('\n  Category 4 — false N op (FJF):')
    for label, strand in [
        ('cat4_plus_1', '+'), ('cat4_plus_2', '+'),
        ('cat4_minus_1', '-'), ('cat4_minus_2', '-'),
    ]:
        rid = label_map.get(label)
        row = corrected.get(rid) if rid else None
        if row is None:
            res.check(label, False, skip=True, msg='')
            continue
        orig = int(row['original_3prime'])
        corr = int(row['corrected_3prime'])
        shifted = orig != corr
        direction_ok = (corr < orig) if strand == '+' else (corr > orig)
        ok = shifted and direction_ok
        res.check(label, ok,
                  msg=f'corrected={corr} original={orig} (no shift or wrong direction)',
                  detail=f'{corr - orig:+d} bp')


def _check_cat5(res: _Results, label_map: Dict, bam: Path) -> None:
    print('\n  Category 5 — chimeric reconstruction (BAM tag check):')
    try:
        import pysam
        reads_by_name: Dict[str, object] = {}
        with pysam.AlignmentFile(str(bam), 'rb') as f:
            for r in f:
                if r.is_secondary or r.is_supplementary:
                    continue
                reads_by_name[r.query_name] = r
    except ImportError:
        for label in ['cat5_plus_3aligner', 'cat5_plus_2aligner',
                      'cat5_minus_long', 'cat5_minus_short']:
            res.check(label, False, skip=True, msg='')
        return

    for label in ['cat5_plus_3aligner', 'cat5_plus_2aligner',
                  'cat5_minus_long', 'cat5_minus_short']:
        rid = label_map.get(label)
        r = reads_by_name.get(rid) if rid else None
        if r is None:
            res.check(label, False, skip=True, msg='')
            continue
        try:
            xk = r.get_tag('XK')
            xa = r.get_tag('XA')
            xs = r.get_tag('XS') if r.has_tag('XS') else 0
        except KeyError as e:
            res.check(label, False, msg=f'missing tag {e}')
            continue
        ok = xk == 1 and ',' in xa and xs >= 3
        res.check(label, ok,
                  msg=f'XK={xk} XA={xa!r} XS={xs}',
                  detail=f'XA={xa} XS={xs}')


def _check_cat6(res: _Results, label_map: Dict, all_rows: Dict[str, List]) -> None:
    print('\n  Category 6 — NET-seq A-tract refinement:')
    singles = ['cat6_plus_single', 'cat6_minus_single']
    multis  = ['cat6_plus_multi',  'cat6_minus_multi']

    for label in singles:
        rid = label_map.get(label)
        rows = all_rows.get(rid, []) if rid else []
        if not rows:
            res.check(label, False, skip=True, msg='')
            continue
        frac = float(rows[0].get('fraction', 1.0))
        ok = len(rows) == 1 and abs(frac - 1.0) < 0.01
        res.check(label, ok,
                  msg=f'{len(rows)} rows, fraction={frac:.3f} (expected 1 row, fraction=1.0)',
                  detail=f'1 row, fraction={frac:.3f}')

    for label in multis:
        rid = label_map.get(label)
        rows = all_rows.get(rid, []) if rid else []
        if not rows:
            res.check(label, False, skip=True, msg='')
            continue
        total = sum(float(r.get('fraction', 1.0)) for r in rows)
        ok = len(rows) >= 2 and abs(total - 1.0) < 0.02
        res.check(label, ok,
                  msg=f'{len(rows)} rows, sum={total:.3f} (expected ≥2 rows summing to 1.0)',
                  detail=f'{len(rows)} rows, sum={total:.3f}')

# ---------------------------------------------------------------------------
# gene_id attribution check
# ---------------------------------------------------------------------------

# Reads confirmed to overlap an annotated gene (from VALIDATION_READS.md coords).
# cat1_plus_2, cat2_minus_1, cat2_minus_2 are genuinely intergenic — excluded.
_EXPECTED_GENIC = {
    'cat1_plus_1', 'cat1_minus_1', 'cat1_minus_2',
    'cat2_plus_1', 'cat2_plus_2',
    'cat3_plus_1', 'cat3_plus_2', 'cat3_minus_1', 'cat3_minus_2',
    'cat4_plus_1', 'cat4_plus_2', 'cat4_minus_1', 'cat4_minus_2',
    # Cat5: only cat5_plus_3aligner overlaps an annotated gene (YAL066W, chrI:9806-10460+);
    # the other three Cat5 reads land in intergenic or wrong-strand regions.
    'cat5_plus_3aligner',
    'cat6_plus_single', 'cat6_plus_multi', 'cat6_minus_single', 'cat6_minus_multi',
}


def _check_gene_id(res: _Results, label_map: Dict, corrected: Dict) -> None:
    print('\n  Gene attribution (gene_id per read body):')
    for label in sorted(_EXPECTED_GENIC):
        rid = label_map.get(label)
        row = corrected.get(rid) if rid else None
        if row is None:
            res.check(label, False, skip=True, msg='')
            continue
        gid = row.get('gene_id', '')
        ok = bool(gid and gid not in ('nan', 'None', ''))
        res.check(label, ok,
                  msg='gene_id is empty (annotation not loaded or interval tree mismatch)',
                  detail=gid)


# ---------------------------------------------------------------------------
# Alignment step (--full mode)
# ---------------------------------------------------------------------------

def _run_align(fastq: Path, genome: Path, annotation: Optional[Path],
               out_dir: Path, threads: int) -> Optional[Path]:
    """Run rectify align and return path to the output rectified BAM."""
    align_dir = out_dir / 'aligners'
    align_dir.mkdir(parents=True, exist_ok=True)

    cmd = _rectify(
        'align',
        str(fastq),
        '--genome', str(genome),
        '--output', str(out_dir),
        '--threads', str(threads),
    )
    if annotation:
        cmd += ['--annotation', str(annotation)]

    print(f'  Running rectify align (threads={threads})...')
    t0 = time.time()
    r = _run(cmd, 'rectify align')
    elapsed = time.time() - t0
    if r.returncode != 0:
        return None
    print(f'  Alignment finished in {elapsed:.0f}s')

    # rectify align writes <output>/rectified.bam
    out_bam = out_dir / 'rectified.bam'
    if not out_bam.exists():
        print(f'  [{FAIL}] Expected rectified.bam at {out_bam}')
        return None
    return out_bam

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run(args: argparse.Namespace) -> int:
    from rectify import __version__
    t_start = time.time()

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    width = 62
    print('=' * width)
    print(_bold(f'  RECTIFY installation test') + f'  (v{__version__})')
    print('=' * width)

    # --- Locate bundled data ------------------------------------------------
    genome = Path(args.genome) if args.genome else _get_genome()
    annotation = Path(args.annotation) if args.annotation else _get_annotation()

    if genome is None or not genome.exists():
        print(f'[{FAIL}] Bundled genome not found. '
              f'Pass --genome or reinstall RECTIFY.')
        return 1

    ann_label = str(annotation) if annotation else 'not available'
    print(f'\nBundled data:')
    print(f'  genome     : {genome}')
    print(f'  annotation : {ann_label}')
    print(f'  output dir : {out_dir}')

    # --- Choose input BAM ---------------------------------------------------
    if args.full:
        fastq = _bundled_fastq()
        if not fastq.exists():
            print(f'[{FAIL}] Bundled FASTQ not found: {fastq}')
            return 1

        bam = _run_align(fastq, genome, annotation, out_dir, args.threads)
        if bam is None:
            return 1
    else:
        bam = _bundled_bam()
        if not bam.exists():
            print(f'[{FAIL}] Bundled validation BAM not found: {bam}')
            return 1
        print(f'\nUsing pre-aligned validation BAM (quick mode; use --full to '
              f'test alignment)')

    # --- Map XV labels → read_ids -------------------------------------------
    label_map = _read_names_by_label(_bundled_bam())
    if not label_map:
        print(f'[{FAIL}] Could not read XV labels from validation BAM '
              f'(pysam not installed?)')
        return 1

    # --- Run rectify correct ------------------------------------------------
    tsv_path = out_dir / 'corrected_3ends.tsv'
    cmd = _rectify(
        'correct',
        str(bam),
        '--genome', str(genome),
        '-o', str(tsv_path),
        '--threads', str(args.threads),
    )
    if annotation:
        cmd += ['--annotation', str(annotation)]

    print(f'\nRunning rectify correct...')
    t0 = time.time()
    r = _run(cmd, 'rectify correct')
    elapsed = time.time() - t0

    if r.returncode != 0:
        return 1
    print(f'  Correction finished in {elapsed:.0f}s  → {tsv_path.name}')

    corrected = _load_tsv_first_row(tsv_path)
    all_rows  = _load_tsv_all_rows(tsv_path)

    # --- Run checks ---------------------------------------------------------
    print(f'\nChecking correction categories:')
    res = _Results()
    _check_cat1(res, label_map, corrected)
    _check_cat2(res, label_map, corrected)
    _check_cat3(res, label_map, corrected)
    _check_cat4(res, label_map, corrected)
    _check_cat5(res, label_map, bam)
    _check_cat6(res, label_map, all_rows)
    _check_gene_id(res, label_map, corrected)

    # --- Summary ------------------------------------------------------------
    total = time.time() - t_start
    n_total = res.passed + res.failed + res.skipped
    status_str = _green(f'{res.passed}/{n_total} passed') if res.failed == 0 \
        else _red(f'{res.failed}/{n_total} FAILED')

    print(f'\n{"=" * width}')
    print(f'  Result  : {status_str}', end='')
    if res.skipped:
        print(f'  ({res.skipped} skipped)', end='')
    print()
    print(f'  Time    : {total:.1f}s')
    print(f'  Output  : {out_dir}')
    print('=' * width)

    if not args.keep and not args.full:
        shutil.rmtree(out_dir, ignore_errors=True)
    elif not args.keep and args.full:
        # Keep aligners/ for inspection but remove TSV
        tsv_path.unlink(missing_ok=True)

    return 0 if res.failed == 0 else 1


# ---------------------------------------------------------------------------
# CLI parser
# ---------------------------------------------------------------------------

def create_test_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    p = subparsers.add_parser(
        'test',
        help='Run installation smoke-test on bundled validation reads',
        description=(
            'Verify RECTIFY is correctly installed by running correction on '
            '24 bundled validation reads (all six correction categories) and '
            'checking that expected corrections are applied.\n\n'
            'Quick mode (default): uses the pre-aligned validation BAM — '
            'no external aligner tools required.\n\n'
            'Full mode (--full): starts from the bundled FASTQ and runs the '
            'complete align → correct pipeline. Requires minimap2 on PATH.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Examples:\n'
               '  rectify test                   # quick smoke-test\n'
               '  rectify test --full            # full pipeline test\n'
               '  rectify test --output-dir /tmp/rectify_test --keep\n',
    )
    p.add_argument(
        '--output-dir', default='rectify_test',
        metavar='DIR',
        help='Directory for test outputs (default: rectify_test/)',
    )
    p.add_argument(
        '--genome', default=None, metavar='FASTA',
        help='Genome FASTA (default: bundled S. cerevisiae genome)',
    )
    p.add_argument(
        '--annotation', default=None, metavar='GFF',
        help='Annotation GFF (default: bundled S. cerevisiae annotation)',
    )
    p.add_argument(
        '--full', action='store_true',
        help='Run full pipeline from FASTQ including alignment step',
    )
    p.add_argument(
        '--keep', action='store_true',
        help='Retain output directory after test completes',
    )
    p.add_argument(
        '--threads', type=int, default=4, metavar='N',
        help='Threads for correction/alignment (default: 4)',
    )
    return p
