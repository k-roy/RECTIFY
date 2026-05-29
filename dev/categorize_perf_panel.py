#!/usr/bin/env python3
"""Task #13: categorize the genome-wide expensive-read panel into cost archetypes.

Reads the per-aligner panel TSVs emitted by profile_correct_reads.py --cost-threshold-ms
(one row per read whose correct_read_3prime cost >= threshold) plus the per-aligner
whole-set summary TSVs, assigns each panel read to ONE disjoint archetype by precedence,
and reports per-archetype read counts AND share of total rescue (panel) time, per aligner
and pooled.

Archetypes (precedence high->low; a read is assigned to the FIRST that matches):
  excluded_BUG           : is_excluded == 1        (BUG SIGNAL: rRNA/PolIII reads
                                                    should SKIP rescue at bam_processor
                                                    L423-426 and never reach the >500ms
                                                    panel. Any time here = exclusion-gate
                                                    leak; flagged first to surface it.)
  bigN_artifact_intron   : n_N > 2000              (spurious mega-intron / mis-relabel)
  large_5prime_clip      : five_clip >= 100        (long rescue_seq -> 200x200 DP ceiling)
  over_cap_candidates    : n_pool_nearby > 25      (exceeds the K=25 proximity cap)
  high_cigar_op          : n_cigar_ops >= 50       (garbled / high-mismatch, many N-ops)
  dead_end_with_work     : is_noop==1 AND (n_pool_nearby>0 OR n_op_intervals>0)
                                                    (ran the DP, emitted nothing)
  productive_candidate   : (five_rescued OR three_shifted) AND n_pool_nearby>0
  other                  : everything else in the panel

Usage:
  python dev/categorize_perf_panel.py --panel-dir dev/perf_panel_t13_results \
      [--threshold-ms 500]
"""
import argparse
import glob
import os
from collections import defaultdict


ARCHETYPE_ORDER = [
    'excluded_BUG',
    'bigN_artifact_intron',
    'large_5prime_clip',
    'over_cap_candidates',
    'high_cigar_op',
    'dead_end_with_work',
    'productive_candidate',
    'other',
]


def classify(r):
    n_N = int(r['n_N'])
    five_clip = int(r['five_clip'])
    n_pool = int(r['n_pool_nearby'])
    n_ops = int(r['n_cigar_ops'])
    n_iv = int(r['n_op_intervals'])
    is_excl = int(r['is_excluded'])
    is_noop = int(r['is_noop'])
    five_rescued = int(r['five_rescued'])
    three_shifted = int(r['three_shifted'])
    if is_excl:
        return 'excluded_BUG'   # should never appear: exclusion gate skips rescue
    if n_N > 2000:
        return 'bigN_artifact_intron'
    if five_clip >= 100:
        return 'large_5prime_clip'
    if n_pool > 25:
        return 'over_cap_candidates'
    if n_ops >= 50:
        return 'high_cigar_op'
    if is_noop and (n_pool > 0 or n_iv > 0):
        return 'dead_end_with_work'
    if (five_rescued or three_shifted) and n_pool > 0:
        return 'productive_candidate'
    return 'other'


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            vals = line.rstrip('\n').split('\t')
            if len(vals) != len(header):
                continue
            yield dict(zip(header, vals))


def read_summary(path):
    """Return dict of the single data row."""
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        line = fh.readline().rstrip('\n')
        if not line:
            return None
        return dict(zip(header, line.split('\t')))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel-dir', required=True)
    ap.add_argument('--threshold-ms', type=float, default=None)
    ap.add_argument('--md-out', default='',
                    help='if set, also emit a markdown table block to this file.')
    args = ap.parse_args()

    panel_files = sorted(glob.glob(os.path.join(args.panel_dir, 'profile_*.panel.tsv')))
    summary_files = {os.path.basename(p).replace('summary_', '').replace('.tsv', ''): p
                     for p in glob.glob(os.path.join(args.panel_dir, 'summary_*.tsv'))}

    # per-aligner: archetype -> [count, ms]
    per_aln = {}
    summaries = {}
    pooled = defaultdict(lambda: [0, 0.0])
    pooled_total_scanned = 0
    pooled_total_time = 0.0
    pooled_panel_ms = 0.0

    md_lines = []

    for pf in panel_files:
        aln = os.path.basename(pf).replace('profile_', '').replace('.panel.tsv', '')
        arch = defaultdict(lambda: [0, 0.0])
        panel_ms_aln = 0.0
        panel_n_aln = 0
        for r in read_tsv(pf):
            cat = classify(r)
            ms = float(r['elapsed_ms'])
            arch[cat][0] += 1
            arch[cat][1] += ms
            panel_ms_aln += ms
            panel_n_aln += 1
            pooled[cat][0] += 1
            pooled[cat][1] += ms
        per_aln[aln] = (arch, panel_n_aln, panel_ms_aln)
        pooled_panel_ms += panel_ms_aln

        s = read_summary(summary_files[aln]) if aln in summary_files else None
        summaries[aln] = s
        if s:
            pooled_total_scanned += int(float(s['reads_scanned']))
            pooled_total_time += float(s['total_time_s'])

    # ---- report ----
    def emit(line=''):
        print(line)
        md_lines.append(line)

    emit(f"# Task #13 expensive-read panel categorization")
    emit(f"panel-dir: {args.panel_dir}")
    emit('')
    emit("## Whole-set + panel sizing (per aligner)")
    emit('')
    emit("| Aligner | reads scanned | stride | total correct time (s) | "
         "threshold (ms) | panel reads | panel time (s) | panel share of total |")
    emit("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for aln in sorted(per_aln):
        s = summaries.get(aln)
        arch, pn, pms = per_aln[aln]
        if s:
            emit(f"| {aln} | {int(float(s['reads_scanned'])):,} | {s['stride']} | "
                 f"{float(s['total_time_s']):.1f} | {float(s['threshold_ms']):.0f} | "
                 f"{pn:,} | {pms/1000:.1f} | {float(s['panel_share_pct']):.1f}% |")
        else:
            emit(f"| {aln} | ? | ? | ? | ? | {pn:,} | {pms/1000:.1f} | ? |")
    emit(f"| **POOLED** | {pooled_total_scanned:,} | - | {pooled_total_time:.1f} | - | "
         f"{sum(v[0] for v in pooled.values()):,} | {pooled_panel_ms/1000:.1f} | "
         f"{100*pooled_panel_ms/pooled_total_time if pooled_total_time else 0:.1f}% |")
    emit('')

    # per-aligner archetype tables (count + % of that aligner's panel time)
    emit("## Archetype breakdown per aligner (count | % of panel time)")
    emit('')
    header = "| Archetype | " + " | ".join(sorted(per_aln)) + " | POOLED |"
    emit(header)
    emit("| --- | " + " | ".join(["---:"] * (len(per_aln) + 1)) + " |")
    for cat in ARCHETYPE_ORDER:
        cells = []
        for aln in sorted(per_aln):
            arch, pn, pms = per_aln[aln]
            c, ms = arch.get(cat, [0, 0.0])
            pct = 100 * ms / pms if pms else 0
            cells.append(f"{c} ({pct:.1f}%)")
        pc, pmsv = pooled.get(cat, [0, 0.0])
        ppct = 100 * pmsv / pooled_panel_ms if pooled_panel_ms else 0
        cells.append(f"{pc} ({ppct:.1f}%)")
        emit(f"| `{cat}` | " + " | ".join(cells) + " |")
    emit('')

    # pooled time-share summary (the headline table)
    emit("## POOLED archetype time-share (the headline)")
    emit('')
    emit("| Archetype | panel reads | panel time (s) | % of panel time |")
    emit("| --- | ---: | ---: | ---: |")
    for cat in ARCHETYPE_ORDER:
        c, ms = pooled.get(cat, [0, 0.0])
        pct = 100 * ms / pooled_panel_ms if pooled_panel_ms else 0
        emit(f"| `{cat}` | {c:,} | {ms/1000:.1f} | {pct:.1f}% |")
    emit(f"| **TOTAL** | {sum(v[0] for v in pooled.values()):,} | "
         f"{pooled_panel_ms/1000:.1f} | 100.0% |")
    emit('')

    # feature medians per archetype (pooled) for write-up grounding
    emit("## Pooled feature medians per archetype")
    emit('')
    feats = ['elapsed_ms', 'five_clip', 'n_pool_nearby', 'n_cigar_ops', 'n_N',
             'n_op_intervals', 'mapq']
    by_cat_feat = defaultdict(lambda: defaultdict(list))
    for pf in panel_files:
        for r in read_tsv(pf):
            cat = classify(r)
            for f in feats:
                try:
                    by_cat_feat[cat][f].append(float(r[f]))
                except (ValueError, KeyError):
                    pass
    emit("| Archetype | n | " + " | ".join(f"med_{f}" for f in feats) + " |")
    emit("| --- | ---: | " + " | ".join(["---:"] * len(feats)) + " |")
    for cat in ARCHETYPE_ORDER:
        vals = by_cat_feat.get(cat)
        if not vals:
            continue
        n = len(vals[feats[0]]) if feats[0] in vals else 0
        meds = []
        for f in feats:
            v = sorted(vals.get(f, []))
            meds.append(f"{v[len(v)//2]:.1f}" if v else "-")
        emit(f"| `{cat}` | {n} | " + " | ".join(meds) + " |")
    emit('')

    if args.md_out:
        with open(args.md_out, 'w') as fh:
            fh.write('\n'.join(md_lines) + '\n')
        print(f"\nwrote markdown: {args.md_out}")


if __name__ == '__main__':
    main()
