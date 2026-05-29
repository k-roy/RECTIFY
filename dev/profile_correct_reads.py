#!/usr/bin/env python3
"""Per-read timing profiler for `rectify correct`'s Module 2F 3'SS rescue.

Task #12 (perf): pull out the top ~1% most computationally intensive reads in
`correct_read_3prime` POST all the 2026-05-26 fixes (rRNA exclusion `cf5ebb9`,
align_clip_to_exon 500bp cap `da18965`, terminal-peel K=20 `fe1bb69`, baseline
body K=25 `25d7a30`, offset-loop ED=0 early-exit `5ceb243`) and deep-dive what
STILL dominates + which reads are dead ends.

This is a STANDALONE diagnostic — NOT a core hot-path change. It replicates the
real `correct` setup from `correct_command.py`:
  - load genome + annotation
  - build the cross-aligner junction pool from the 5 per-aligner BAMs
  - build pool_chrom_index via bam_processor._build_pool_chrom_index
  - build the ExclusionRegionDetector (rDNA/Pol III) the same way correct does
  - build gene_interval_trees
then iterate a representative BAM, timing correct_read_3prime() per read with
full per-read feature extraction, sorts, and writes the top-1% to a TSV.

Usage:
  python dev/profile_correct_reads.py \
      --bam        <representative.bam> \
      --genome     <genome.fsa[.gz]> \
      --annotation <annotation.gff[.gz]> \
      --aligner-bams b1.bam b2.bam ... \
      --max-reads  4000 \
      --out-prefix /tmp/profile_correct
"""
import argparse
import time
import sys
from pathlib import Path

import pysam

from rectify.core.bam import bam_processor
from rectify.core.bam.bam_processor import (
    correct_read_3prime,
    get_read_5prime_position,
    get_read_3prime_position,
    _build_pool_chrom_index,
)
from rectify.core.splice.splice_aware_5prime import (
    _get_5prime_softclip_len,
    _get_n_op_intervals,
)
from rectify.utils.genome import load_genome, standardize_chrom_name, register_genome_contigs_from_fasta


def cigar_complexity(read):
    """Return (n_ops, n_S, n_I, n_D, n_N, n_X, ref_len, query_len, max_clip)."""
    if not read.cigartuples:
        return (0, 0, 0, 0, 0, 0, 0, 0, 0)
    n_ops = len(read.cigartuples)
    counts = {1: 0, 2: 0, 3: 0, 4: 0, 8: 0}  # I D N S X
    max_clip = 0
    for op, ln in read.cigartuples:
        if op in counts:
            counts[op] += ln
        if op == 4 and ln > max_clip:
            max_clip = ln
    ref_len = read.reference_length or 0
    qlen = read.query_length or 0
    return (n_ops, counts[4], counts[1], counts[2], counts[3], counts[8],
            ref_len, qlen, max_clip)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True)
    ap.add_argument('--genome', required=True)
    ap.add_argument('--annotation', required=True)
    ap.add_argument('--aligner-bams', nargs='*', default=[])
    ap.add_argument('--max-reads', type=int, default=4000)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--junction-max-size', type=int, default=10000)
    ap.add_argument('--include-rdna', action='store_true',
                    help='disable rRNA/Pol III exclusion (default: excluded, like correct)')
    ap.add_argument('--aligner-label', default='',
                    help='label recorded in the aligner column (e.g. minimap2)')
    ap.add_argument('--sample-stride', type=int, default=1,
                    help='take 1 of every N reads (spread the sample across all '
                         'chromosomes instead of just chrI). 1 = first max-reads.')
    ap.add_argument('--cost-threshold-ms', type=float, default=None,
                    help='task #13 panel mode: stream EVERY read whose '
                         'correct_read_3prime elapsed_ms >= this threshold to the '
                         'panel TSV, while accumulating whole-set aggregates '
                         '(total reads scanned, total time, percentile histogram) '
                         'instead of buffering all rows. Use with --sample-stride to '
                         'sweep the whole BAM genome-wide. When set, --max-reads is '
                         'interpreted as a cap on reads PROFILED (post-stride), not a '
                         'hard stop, and the top1pct/all TSVs are not written.')
    ap.add_argument('--pyspy-read', default='',
                    help='if set, run correct_read_3prime on ONLY this read_id, '
                         'repeated --pyspy-repeat times (for py-spy targeting). '
                         'Skips the normal profile.')
    ap.add_argument('--pyspy-repeat', type=int, default=200)
    ap.add_argument('--summary-out', default='',
                    help='task #13 panel mode: write a one-row TSV of whole-set '
                         'aggregates (reads scanned, total time, panel size + share).')
    ap.add_argument('--apply-indel-correction', action='store_true',
                    help='match production correct (2C indel correction ON). Off by '
                         'default so the profile isolates the 3SS-rescue cost.')
    args = ap.parse_args()

    t0 = time.perf_counter()
    print(f"[setup] loading genome {args.genome}", file=sys.stderr)
    genome = load_genome(args.genome)
    register_genome_contigs_from_fasta(args.genome)

    # --- annotated junctions (Module 2F) ---
    from rectify.core.consensus.consensus import load_annotated_junctions
    annotated_junctions = load_annotated_junctions(args.annotation)
    print(f"[setup] {len(annotated_junctions):,} annotated junctions", file=sys.stderr)

    # --- cross-aligner junction pool + pool_chrom_index (mirrors correct_command) ---
    pool_chrom_index = None
    if args.aligner_bams:
        from rectify.core.splice.junction_refiner import build_junction_pool
        pool, annot_set = build_junction_pool(
            args.aligner_bams, annotated_junctions,
            max_junction_size=args.junction_max_size,
        )
        print(f"[setup] junction pool: {len(pool):,} ({len(annot_set):,} annotated)",
              file=sys.stderr)
        jmax = args.junction_max_size
        if jmax is not None:
            pool = {
                j for j in pool
                if (j in annot_set or (len(j) >= 3 and int(j[2]) - int(j[1]) <= int(jmax)))
            }
        pool_chrom_index = _build_pool_chrom_index(pool)

    # --- gene interval trees (per-read attribution) ---
    gene_interval_trees = None
    try:
        from rectify.core.analyze.gene_attribution import build_cds_interval_tree
        from rectify.core.analyze.loaders import load_annotation as _load_annotation_for_trees
        _ann_df = _load_annotation_for_trees(args.annotation, normalize_chroms=False)
        gene_interval_trees = build_cds_interval_tree(_ann_df)
    except Exception as e:
        print(f"[setup] gene trees unavailable: {e}", file=sys.stderr)

    # --- exclusion detector (rDNA / Pol III), default-on like correct ---
    exclusion_detector = None
    if not args.include_rdna:
        try:
            from rectify.core.exclusion_regions import ExclusionRegionDetector
            exclusion_detector = ExclusionRegionDetector(flanking_bp=100)
            n_excl = exclusion_detector.load_from_gff(
                Path(args.annotation),
                exclude_tRNA=True, exclude_snRNA=True,
                exclude_rDNA=True, exclude_mito=False,
            )
            print(f"[setup] exclusion regions: {n_excl}", file=sys.stderr)
            if n_excl == 0:
                exclusion_detector = None
        except Exception as e:
            print(f"[setup] exclusion detector unavailable: {e}", file=sys.stderr)
            exclusion_detector = None

    print(f"[setup] done in {time.perf_counter()-t0:.1f}s; profiling reads...",
          file=sys.stderr)

    bam = pysam.AlignmentFile(args.bam, 'rb')

    # --- py-spy targeting mode: run ONLY the named read, repeated, then exit ---
    if args.pyspy_read:
        target = None
        for read in bam.fetch(until_eof=True):
            if read.query_name == args.pyspy_read and not read.is_secondary \
                    and not read.is_supplementary and not read.is_unmapped:
                target = read
                break
        if target is None:
            print(f"[pyspy] read {args.pyspy_read} not found", file=sys.stderr)
            sys.exit(1)
        print(f"[pyspy] found {args.pyspy_read}; repeating {args.pyspy_repeat}x "
              f"(attach py-spy NOW to this pid)", file=sys.stderr)
        t = time.perf_counter()
        for _ in range(args.pyspy_repeat):
            # fresh copy each iter — correct_read_3prime mutates CIGAR in-place then restores
            correct_read_3prime(
                target, genome,
                annotated_junctions=annotated_junctions,
                pool_chrom_index=pool_chrom_index,
                apply_3ss_rescue=True,
                apply_indel_correction=args.apply_indel_correction,
                gene_interval_trees=gene_interval_trees,
                exclusion_detector=exclusion_detector,
            )
        dt = time.perf_counter() - t
        print(f"[pyspy] {args.pyspy_repeat}x in {dt:.2f}s "
              f"= {1000*dt/args.pyspy_repeat:.1f}ms/call", file=sys.stderr)
        return

    cols = ['read_id', 'elapsed_ms', 'chrom', 'pos5', 'strand', 'mapq', 'aligner',
            'ref_len', 'query_len', 'n_cigar_ops', 'five_clip', 'max_clip',
            'n_S', 'n_I', 'n_D', 'n_N', 'n_X', 'n_op_intervals', 'n_pool_nearby',
            'is_excluded', 'five_rescued', 'three_shifted', 'corr_category', 'is_noop']

    # ---- task #13 panel mode: stream reads above a cost threshold ----
    panel_mode = args.cost_threshold_ms is not None
    panel_fh = None
    panel_n = 0          # reads written to panel
    panel_ms = 0.0       # total time of panel reads
    if panel_mode:
        out_panel = f"{args.out_prefix}.panel.tsv"
        panel_fh = open(out_panel, 'w')
        panel_fh.write('\t'.join(cols) + '\n')
        # log10-ish histogram buckets (ms) for whole-set time concentration
        hist_edges = [0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0,
                      250.0, 500.0, 1000.0, 5000.0, 1e9]
        hist_cnt = [0] * (len(hist_edges) - 1)
        hist_ms = [0.0] * (len(hist_edges) - 1)

    rows = []
    n = 0
    seen = 0
    n_total_time = 0.0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        seen += 1
        if args.sample_stride > 1 and (seen % args.sample_stride) != 0:
            continue
        n += 1
        if n > args.max_reads:
            break

        chrom_std = standardize_chrom_name(read.reference_name)
        strand = '-' if read.is_reverse else '+'
        five_p = get_read_5prime_position(read, strand)
        # features (cheap, computed before timing the rescue)
        (n_ops, n_S, n_I, n_D, n_N, n_X, ref_len, qlen, max_clip) = cigar_complexity(read)
        five_clip = _get_5prime_softclip_len(read)
        n_n_intervals = len(_get_n_op_intervals(read))
        is_excl = bool(exclusion_detector is not None and five_p is not None
                       and exclusion_detector.is_excluded(chrom_std, five_p))
        # n nearby pool junctions (same window correct uses)
        n_pool_nearby = 0
        if pool_chrom_index is not None and five_p is not None:
            tree = pool_chrom_index.get(chrom_std)
            if tree:
                w = bam_processor._POOL_FETCH_HALF_WINDOW + five_clip
                n_pool_nearby = len(tree.overlap(five_p - w, five_p + w + 1))

        t_read = time.perf_counter()
        try:
            result = correct_read_3prime(
                read, genome,
                annotated_junctions=annotated_junctions,
                pool_chrom_index=pool_chrom_index,
                apply_3ss_rescue=True,
                apply_indel_correction=args.apply_indel_correction,
                gene_interval_trees=gene_interval_trees,
                exclusion_detector=exclusion_detector,
            )
        except Exception as e:
            result = []
            print(f"[warn] read {read.query_name} raised {type(e).__name__}: {e}",
                  file=sys.stderr)
        elapsed_ms = (time.perf_counter() - t_read) * 1000.0
        n_total_time += elapsed_ms

        # outcome classification (from emitted dict)
        five_rescued = False
        corr_cat = None
        orig3 = None
        corr3 = None
        if result:
            r0 = result[0]
            five_rescued = bool(r0.get('five_prime_rescued', False))
            corr_cat = r0.get('correction_category')
            orig3 = r0.get('original_3prime')
            corr3 = r0.get('corrected_3prime')
        three_shifted = (orig3 is not None and corr3 is not None and orig3 != corr3)
        is_noop = (not five_rescued) and (not three_shifted) and (corr_cat is None)

        row = {
            'read_id': read.query_name,
            'elapsed_ms': elapsed_ms,
            'chrom': chrom_std,
            'pos5': five_p if five_p is not None else -1,
            'strand': strand,
            'mapq': read.mapping_quality,
            'aligner': args.aligner_label,
            'ref_len': ref_len,
            'query_len': qlen,
            'n_cigar_ops': n_ops,
            'five_clip': five_clip,
            'max_clip': max_clip,
            'n_S': n_S, 'n_I': n_I, 'n_D': n_D, 'n_N': n_N, 'n_X': n_X,
            'n_op_intervals': n_n_intervals,
            'n_pool_nearby': n_pool_nearby,
            'is_excluded': int(is_excl),
            'five_rescued': int(five_rescued),
            'three_shifted': int(three_shifted),
            'corr_category': corr_cat if corr_cat is not None else '',
            'is_noop': int(is_noop),
        }

        if panel_mode:
            # whole-set aggregates (don't buffer every row).
            # n_total_time already accumulated above (line ~264).
            # histogram bucket
            for _b in range(len(hist_edges) - 1):
                if elapsed_ms < hist_edges[_b + 1]:
                    hist_cnt[_b] += 1
                    hist_ms[_b] += elapsed_ms
                    break
            if elapsed_ms >= args.cost_threshold_ms:
                panel_fh.write('\t'.join(str(row[c]) for c in cols) + '\n')
                panel_n += 1
                panel_ms += elapsed_ms
                if (panel_n % 200) == 0:
                    panel_fh.flush()
        else:
            rows.append(row)
    bam.close()

    # ---- panel-mode summary + early return ----
    if panel_mode:
        panel_fh.close()
        print(f"\n===== PANEL SUMMARY (aligner={args.aligner_label or '?'}) =====")
        print(f"reads scanned (post-stride): {n}   stride={args.sample_stride}")
        print(f"total correct time: {n_total_time/1000:.2f}s")
        print(f"threshold: >= {args.cost_threshold_ms:.1f} ms")
        print(f"panel reads (>= threshold): {panel_n}   "
              f"panel time: {panel_ms/1000:.2f}s "
              f"({100*panel_ms/n_total_time if n_total_time else 0:.1f}% of total)")
        print("time concentration histogram (ms bucket: n_reads  sum_s  %total_time):")
        for _b in range(len(hist_edges) - 1):
            lab = f"[{hist_edges[_b]:.2f},{hist_edges[_b+1]:.2f})"
            pct = 100 * hist_ms[_b] / n_total_time if n_total_time else 0
            print(f"  {lab:>22}: {hist_cnt[_b]:>8}  {hist_ms[_b]/1000:>8.2f}s  {pct:5.1f}%")
        if args.summary_out:
            with open(args.summary_out, 'w') as sfh:
                sfh.write("aligner\treads_scanned\tstride\ttotal_time_s\tthreshold_ms\t"
                          "panel_reads\tpanel_time_s\tpanel_share_pct\n")
                sfh.write(f"{args.aligner_label}\t{n}\t{args.sample_stride}\t"
                          f"{n_total_time/1000:.3f}\t{args.cost_threshold_ms}\t"
                          f"{panel_n}\t{panel_ms/1000:.3f}\t"
                          f"{100*panel_ms/n_total_time if n_total_time else 0:.2f}\n")
        print(f"wrote: {args.out_prefix}.panel.tsv")
        if args.summary_out:
            print(f"wrote: {args.summary_out}")
        return

    rows.sort(key=lambda r: r['elapsed_ms'], reverse=True)
    total_ms = sum(r['elapsed_ms'] for r in rows)
    n_reads = len(rows)
    top1pct_n = max(1, n_reads // 100)
    top = rows[:top1pct_n]
    top_ms = sum(r['elapsed_ms'] for r in top)

    # write top-1% TSV  (cols defined above, before the read loop)
    out_top = f"{args.out_prefix}.top1pct.tsv"
    with open(out_top, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in top:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # also write ALL rows for downstream analysis
    out_all = f"{args.out_prefix}.all.tsv"
    with open(out_all, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # summary
    def pctl(vals, p):
        if not vals:
            return 0.0
        s = sorted(vals)
        i = min(len(s) - 1, int(p / 100.0 * len(s)))
        return s[i]
    all_ms = [r['elapsed_ms'] for r in rows]
    print(f"\n===== PROFILE SUMMARY (aligner={args.aligner_label or '?'}) =====")
    print(f"reads profiled: {n_reads}   total time: {total_ms/1000:.2f}s")
    print(f"p50={pctl(all_ms,50):.3f}ms  p95={pctl(all_ms,95):.3f}ms  "
          f"p99={pctl(all_ms,99):.3f}ms  max={max(all_ms):.1f}ms")
    print(f"top-1% ({top1pct_n} reads): {top_ms/1000:.2f}s = {100*top_ms/total_ms:.1f}% of all time")
    # breakdown of top-1% by outcome
    n_top_noop = sum(r['is_noop'] for r in top)
    ms_top_noop = sum(r['elapsed_ms'] for r in top if r['is_noop'])
    n_top_excl = sum(r['is_excluded'] for r in top)
    n_top_5resc = sum(r['five_rescued'] for r in top)
    n_top_3shift = sum(r['three_shifted'] for r in top)
    print(f"  top-1% no-ops: {n_top_noop}/{top1pct_n} reads = "
          f"{ms_top_noop/1000:.2f}s ({100*ms_top_noop/top_ms:.1f}% of top-1% time)")
    print(f"  top-1% excluded(rRNA/PolIII): {n_top_excl}")
    print(f"  top-1% five_rescued: {n_top_5resc}   three_shifted: {n_top_3shift}")
    # feature medians in top-1% vs overall
    def med(rs, key):
        v = sorted(r[key] for r in rs)
        return v[len(v)//2] if v else 0
    print(f"  top-1% median n_pool_nearby={med(top,'n_pool_nearby')} "
          f"(overall={med(rows,'n_pool_nearby')})")
    print(f"  top-1% median five_clip={med(top,'five_clip')} "
          f"(overall={med(rows,'five_clip')})")
    print(f"  top-1% median max_clip={med(top,'max_clip')} "
          f"(overall={med(rows,'max_clip')})")
    print(f"  top-1% median n_op_intervals={med(top,'n_op_intervals')} "
          f"(overall={med(rows,'n_op_intervals')})")
    print(f"  top-1% median ref_len={med(top,'ref_len')} "
          f"(overall={med(rows,'ref_len')})")
    # whole-set no-op burden
    n_noop = sum(r['is_noop'] for r in rows)
    ms_noop = sum(r['elapsed_ms'] for r in rows if r['is_noop'])
    print(f"all no-op reads: {n_noop}/{n_reads} = {ms_noop/1000:.2f}s "
          f"({100*ms_noop/total_ms:.1f}% of all time)")
    print(f"wrote: {out_top}")
    print(f"wrote: {out_all}")


if __name__ == '__main__':
    main()
