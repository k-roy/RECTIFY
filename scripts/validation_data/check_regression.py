#!/usr/bin/env python3
"""
Check for regressions across all 36 validation reads.

Compares two corrected_reads.tsv files (current vs baseline) and reports any
per-read differences in correction-relevant columns. Joins XV (per-read label)
and XG (category) tags from the input BAM so each diff row is identified by
its category label, not just UUID.

Use cases:
  - Quick A/B check after a code change: regenerate corrected_reads.tsv into a
    sidecar path and diff it against the committed baseline.
  - Spot regressions before committing: any unintended shift on a non-target
    read indicates the change wasn't fully scoped.

Usage:
  python scripts/validation_data/check_regression.py \\
      [--current rectify/data/validation/corrected_reads.tsv] \\
      [--baseline <path-or-git-ref>] \\
      [--bam rectify/data/validation/validation_reads.bam]

  --baseline accepts either a filesystem path, or a git-ref
  (HEAD, HEAD~3, <commit-sha>) for the bundled validation TSV.

Output:
  Per-column diff table; exit 0 if no differences, 1 otherwise.
"""
from __future__ import annotations
import argparse
import io
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pysam


# Columns whose changes are functionally meaningful and worth flagging.
# Reorder reflects priority for surfacing the most important shifts first.
DIFF_COLS = [
    'corrected_3prime',
    'five_prime_position',
    'five_prime_rescued',
    'correction_applied',
    'confidence',
    'n_junctions',
    'junctions',
    'sc_rescued_seq',
    'sc_homopolymer_extension',
    'oc_terminal_base',
    'oc_overcall_count',
    'oc_homopolymer_extension',
    'original_3prime',
]


def _load_tsv(spec: str, tsv_relpath: str = 'rectify/data/validation/corrected_reads.tsv') -> pd.DataFrame:
    """Load a TSV either from filesystem or from `git show <ref>:<relpath>`."""
    p = Path(spec)
    if p.exists():
        return pd.read_csv(p, sep='\t')
    # Treat as git ref
    try:
        out = subprocess.check_output(['git', 'show', f'{spec}:{tsv_relpath}'])
    except subprocess.CalledProcessError as exc:
        raise SystemExit(
            f'check_regression: --baseline {spec!r} is neither a file '
            f'nor a resolvable git ref ({tsv_relpath} at {spec}): {exc}'
        ) from exc
    return pd.read_csv(io.BytesIO(out), sep='\t')


def _load_xv_xg(bam_path: Path) -> dict[str, tuple[str, str]]:
    """Map read_id → (XV label, XG category) from input BAM."""
    out: dict[str, tuple[str, str]] = {}
    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        for r in bam.fetch():
            try:
                out[r.query_name] = (r.get_tag('XV'), r.get_tag('XG'))
            except KeyError:
                pass
    return out


def diff_tsvs(current: pd.DataFrame, baseline: pd.DataFrame,
              xv: dict[str, tuple[str, str]]) -> list[dict]:
    cols = [c for c in DIFF_COLS if c in current.columns and c in baseline.columns]
    cur = current.set_index('read_id')
    base = baseline.set_index('read_id')
    shared = cur.index.intersection(base.index)
    only_cur = cur.index.difference(base.index)
    only_base = base.index.difference(cur.index)

    diffs: list[dict] = []
    for rid in shared:
        for col in cols:
            a = str(base.loc[rid, col])
            b = str(cur.loc[rid, col])
            if a != b:
                diffs.append({
                    'read_id': rid,
                    'label':   xv.get(rid, ('?', '?'))[0],
                    'category': xv.get(rid, ('?', '?'))[1],
                    'column':  col,
                    'baseline': a,
                    'current':  b,
                })
    return diffs, list(only_cur), list(only_base)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--current',
                    default='rectify/data/validation/corrected_reads.tsv',
                    help='Path to current corrected_reads.tsv')
    ap.add_argument('--baseline',
                    default='HEAD',
                    help='Baseline: filesystem path OR git ref (default: HEAD)')
    ap.add_argument('--bam',
                    default='rectify/data/validation/validation_reads.bam',
                    help='Input BAM (for XV/XG labels)')
    ap.add_argument('--category',
                    help='Restrict diff output to this XG category (repeatable)',
                    action='append')
    args = ap.parse_args()

    cur = _load_tsv(args.current)
    base = _load_tsv(args.baseline)
    xv = _load_xv_xg(Path(args.bam))

    diffs, only_cur, only_base = diff_tsvs(cur, base, xv)
    if args.category:
        filt = set(args.category)
        diffs = [d for d in diffs if d['category'] in filt]

    print(f'current:  {len(cur)} reads  ({args.current})')
    print(f'baseline: {len(base)} reads  ({args.baseline})')
    if only_cur:
        print(f'  only in current  ({len(only_cur)}): '
              + ', '.join(rid[:8] for rid in only_cur[:10]))
    if only_base:
        print(f'  only in baseline ({len(only_base)}): '
              + ', '.join(rid[:8] for rid in only_base[:10]))

    if not diffs:
        print('\nNo differences.')
        return 0

    # Group by category for readability
    print(f'\n{len(diffs)} differences across '
          f'{len({d["read_id"] for d in diffs})} reads:\n')

    diffs.sort(key=lambda d: (d['category'], d['label'], DIFF_COLS.index(d['column'])
                              if d['column'] in DIFF_COLS else 99))
    last_label = None
    for d in diffs:
        if d['label'] != last_label:
            print(f"  {d['read_id'][:8]}  {d['label']:22s}  ({d['category']})")
            last_label = d['label']
        print(f"      {d['column']:30s}  {d['baseline']!r}  →  {d['current']!r}")
    return 1


if __name__ == '__main__':
    sys.exit(main())
