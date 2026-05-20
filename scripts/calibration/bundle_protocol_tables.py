#!/usr/bin/env python3
"""Merge per-sample profiler outputs and bundle protocol-specific penalty
tables into rectify/data/genomes/<organism>/penalty_tables/.

Inputs:
    --inputs DIR [DIR ...]
        One or more profiler output directories (each containing
        error_counts.tsv and optionally str_error_counts.tsv). Counts are
        summed across all inputs before penalty derivation, so independent
        replicates can be pooled at the raw-event level (denominators stay
        consistent — see compute_rates docstring).

    --protocol {drs, cdna, qsrev}
        Suffix applied to the output filenames and the protocol key under
        which they're surfaced via BUNDLED_GENOMES[...]['protocols'].

    --overhang PATH (optional)
        Path to an overhang TSV to bundle as junction_overhang_table_<protocol>.tsv.
        Omit for QSrev (per handoff §6 — short-read overhang falls back to DRS).

Writes:
    rectify/data/genomes/<organism>/penalty_tables/penalty_scores_<protocol>.tsv
    rectify/data/genomes/<organism>/penalty_tables/str_penalty_scores_<protocol>.tsv
    [optional] junction_overhang_table_<protocol>.tsv

The BUNDLED_GENOMES['saccharomyces_cerevisiae']['protocols'][protocol] map
already points at these filenames — once they land on disk, the
get_bundled_*(organism, protocol=...) resolvers will pick them up.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Tuple

import pandas as pd

# Re-use the profiler's penalty derivation
sys.path.insert(0, str(Path(__file__).parent))
from empirical_cigar_error_profiler import (  # noqa: E402
    compute_rates, derive_penalties, _write_str_outputs,
)


ORGANISM = 'saccharomyces_cerevisiae'  # only supported organism for now
REPO_ROOT = Path(__file__).resolve().parents[2]
BUNDLE_DIR = REPO_ROOT / 'rectify' / 'data' / 'genomes' / ORGANISM / 'penalty_tables'


def _read_counts(
    input_dirs: Iterable[Path],
    hp_filename: str = 'error_counts.tsv',
) -> Tuple[Dict, Dict]:
    """Pool ``hp_filename`` (HP) and ``str_error_counts.tsv`` (STR) across runs.

    The HP filename is parameterized so per-UMI-bin variants (e.g.
    ``error_counts_umi1.tsv``) can be pooled without touching the STR side —
    Phase A's profiler keeps STR pooled-only, so there is no per-bin STR file.
    """
    hp_counts: Dict[Tuple[str, str, int], int] = defaultdict(int)
    str_counts: Dict[Tuple[str, str, int], int] = defaultdict(int)

    for d in input_dirs:
        hp_tsv = d / hp_filename
        if not hp_tsv.exists():
            raise FileNotFoundError(f'{hp_tsv} missing — was the profiler run on {d}?')
        df = pd.read_csv(hp_tsv, sep='\t')
        for _, row in df.iterrows():
            hp_counts[(row['op_type'], row['ref_base'], int(row['hp_length']))] += int(row['count'])

        str_tsv = d / 'str_error_counts.tsv'
        if str_tsv.exists():
            sdf = pd.read_csv(str_tsv, sep='\t')
            for _, row in sdf.iterrows():
                # STR rows are keyed by (op_type, unit, n_copies)
                str_counts[(row['op_type'], row['unit'], int(row['n_copies']))] += int(row['count'])
    return hp_counts, str_counts


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--protocol', choices=['drs', 'cdna', 'qsrev'], required=True)
    p.add_argument('--inputs', type=Path, nargs='+', required=True,
                   help='Per-sample profiler output dirs')
    p.add_argument('--umi-bin', dest='umi_bin', default=None,
                   help='Optional per-UMI-bin label (typically umi1/umi2/umi3plus). '
                        'When set, reads error_counts_<LABEL>.tsv from each --inputs '
                        'dir and writes penalty_scores_<protocol>_<LABEL>.tsv. STR '
                        'output is skipped (Phase A profiler keeps STR pooled-only).')
    p.add_argument('--overhang', type=Path, default=None,
                   help='Path to junction_overhang_table.tsv to bundle (omit for qsrev '
                        'or for per-UMI-bin runs)')
    p.add_argument('--max-hp', type=int, default=20)
    p.add_argument('--default-ins', type=float, default=1.25)
    p.add_argument('--min-count-flag', type=int, default=50)
    p.add_argument('--dry-run', action='store_true',
                   help='Print actions without writing files')
    args = p.parse_args(argv)

    if not BUNDLE_DIR.exists():
        sys.exit(f'ERROR: {BUNDLE_DIR} does not exist. Run from rectify repo root.')

    # --- Pool counts ---
    hp_filename = (f'error_counts_{args.umi_bin}.tsv'
                   if args.umi_bin else 'error_counts.tsv')
    print(f'Pooling counts across {len(args.inputs)} input(s) '
          f'(HP source: {hp_filename}) ...')
    hp_counts, str_counts = _read_counts(args.inputs, hp_filename=hp_filename)
    # Per-bin runs: STR is pooled-only at the profiler stage, so suppress STR
    # output regardless of what _read_counts found.
    if args.umi_bin:
        str_counts = {}
    print(f'  HP counts:  {sum(hp_counts.values())} events across '
          f'{len(hp_counts)} (op, base, hp) keys')
    print(f'  STR counts: {sum(str_counts.values())} events across '
          f'{len(str_counts)} (op, unit, n_copies) keys')

    if not hp_counts:
        sys.exit('ERROR: no HP counts pooled — check --inputs.')

    # --- Derive penalties (re-use profiler logic) ---
    print('Computing rates + deriving penalties ...')
    rate_df = compute_rates(hp_counts)
    rate_df = rate_df[rate_df['hp_length'] <= args.max_hp]
    pen_df = derive_penalties(rate_df,
                              default_ins=args.default_ins,
                              min_count_flag=args.min_count_flag)

    # --- Write outputs ---
    suffix = f'_{args.protocol}'
    if args.umi_bin:
        suffix = f'{suffix}_{args.umi_bin}'
    out_pen  = BUNDLE_DIR / f'penalty_scores{suffix}.tsv'
    out_str  = BUNDLE_DIR / f'str_penalty_scores{suffix}.tsv'
    out_ovh  = BUNDLE_DIR / f'junction_overhang_table{suffix}.tsv'

    if args.dry_run:
        print(f'[dry-run] would write {out_pen}')
        print(f'[dry-run] would write {out_str}')
        if args.overhang:
            print(f'[dry-run] would copy {args.overhang} -> {out_ovh}')
        return

    pen_df.to_csv(out_pen, sep='\t', index=False, float_format='%.4f')
    print(f'Wrote: {out_pen} ({len(pen_df)} rows)')

    if str_counts:
        # Reuse the profiler's STR-output writer for format parity
        _write_str_outputs(str_counts, BUNDLE_DIR, args.max_hp)
        # The writer emits str_penalty_scores.tsv (no suffix); rename it
        emitted = BUNDLE_DIR / 'str_penalty_scores.tsv'
        if emitted.exists() and emitted != out_str:
            shutil.move(emitted, out_str)
            print(f'Wrote: {out_str}')
        # Also tidy up str_error_counts.tsv / str_error_rates.tsv side-effects
        for stray in ('str_error_counts.tsv', 'str_error_rates.tsv'):
            stray_path = BUNDLE_DIR / stray
            if stray_path.exists():
                stray_path.unlink()
    else:
        print(f'  (no STR counts — skipping {out_str.name})')

    if args.overhang:
        if not args.overhang.exists():
            sys.exit(f'ERROR: --overhang {args.overhang} does not exist')
        shutil.copy2(args.overhang, out_ovh)
        print(f'Copied: {args.overhang} -> {out_ovh}')

    print('\nDone. The protocols entry in rectify/data/__init__.py:BUNDLED_GENOMES')
    print('already points at these filenames, so --Scer + the matching protocol')
    print('flags (--dT-primed-cDNA / --short-read) will now resolve to them.')


if __name__ == '__main__':
    main()
