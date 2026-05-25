#!/usr/bin/env python3
"""
aligner_contribution_analysis.py — per-category aligner contribution benchmark

Usage:
    # Single sample (validation bundle):
    python scripts/benchmark/aligner_contribution_analysis.py \
        --samples rectify/data/validation \
        --output /tmp/aligner_contrib_validation

    # Multiple production samples (set2 on Sherlock):
    python scripts/benchmark/aligner_contribution_analysis.py \
        --samples /path/to/wt_rep1 /path/to/ysh1_rep1 /path/to/rna15_rep1 \
        --output /tmp/aligner_contrib_set2 \
        --annotation /path/to/annotation.gtf

Outputs (all TSVs in --output):
    win_rate_by_category.tsv      — winning_aligner × correction_applied cross-tab
    win_rate_by_sample.tsv        — winning_aligner distribution per sample
    leave_one_out_summary.tsv     — per dropped aligner: reads changed, corrections lost
    leave_one_out_detail.tsv      — per dropped aligner × new_winner × new correction_applied
    novel_junction_summary.tsv    — per aligner: unique/annotated/novel junctions
                                    (only when --annotation provided)

correction_applied → cat1-9 approximate mapping:
    indel_correction               → cat1_indel
    softclip_rescue                → cat2_softclip
    five_prime_rescued             → cat3_junction
    (false junc / junction refine) → cat4/cat9 — not directly in correction_applied
    atract_ambiguity               → cat8_netseq_refine
    none                           → cat5/cat6/cat7 or clean read (no correction needed)
"""

import argparse
import csv
import glob
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Sample discovery
# ---------------------------------------------------------------------------

def find_aligner_summary_files(sample_dir: str) -> List[str]:
    """Return all aligner_summary.tsv paths for a sample directory.

    Validation bundle: <sample_dir>/rectified/per_aligner_summary.tsv
    Production:        <sample_dir>/chunks/bypass_2aln_*/chunks/*/aligner_summary.tsv
    """
    val_path = os.path.join(sample_dir, "rectified", "per_aligner_summary.tsv")
    if os.path.exists(val_path):
        return [val_path]

    pattern = os.path.join(
        sample_dir, "chunks", "bypass_2aln_*", "chunks", "*", "aligner_summary.tsv"
    )
    return sorted(glob.glob(pattern))


def find_corrected_reads(sample_dir: str) -> Optional[str]:
    for candidate in [
        os.path.join(sample_dir, "rectified", "corrected_reads.tsv"),
        os.path.join(sample_dir, "corrected_reads.tsv"),
    ]:
        if os.path.exists(candidate):
            return candidate
    return None


# ---------------------------------------------------------------------------
# Win-rate analysis — streaming, O(unique combinations) memory
# ---------------------------------------------------------------------------

def accumulate_win_rates(
    corrected_reads_path: str,
    sample_name: str,
    cross: Dict[Tuple, int],
    by_sample: Dict[Tuple, int],
    totals: Dict[str, int],
    junc_per_aligner: Optional[Dict[str, Dict[Tuple, int]]],
    annotated_junctions: Optional[set],
) -> None:
    """Stream corrected_reads.tsv and accumulate counts in-place."""
    with open(corrected_reads_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            winner = row.get("winning_aligner", "").strip() or "unknown"
            ca = row.get("correction_applied", "none").strip() or "none"
            cross[(winner, ca)] += 1
            by_sample[(sample_name, winner)] += 1
            totals[winner] += 1

            # Junction analysis (optional)
            if junc_per_aligner is not None:
                junc_str = row.get("junctions", "").strip()
                if junc_str:
                    chrom = row.get("chrom", "")
                    strand = row.get("strand", "+")
                    per_win = junc_per_aligner.setdefault(winner, defaultdict(int))
                    for j in junc_str.split(","):
                        j = j.strip()
                        if "-" not in j:
                            continue
                        parts = j.split("-")
                        if len(parts) != 2:
                            continue
                        try:
                            start, end = int(parts[0]), int(parts[1])
                            per_win[(chrom, start, end, strand)] += 1
                        except ValueError:
                            continue


# ---------------------------------------------------------------------------
# Leave-one-out simulation — streaming chunk by chunk
# ---------------------------------------------------------------------------

# Compact tuple for LOO: (aligner, is_winner, rank, correction_applied)
# Avoids the ~300-byte overhead of a full dict per row
_LOO_ALIGNER = 0
_LOO_IS_WINNER = 1
_LOO_RANK = 2
_LOO_CA = 3


def _row_to_tuple(row: dict) -> tuple:
    aligner = row.get("_aligner", "")
    is_winner = row.get("_is_winner", "").lower() == "true"
    # Try _conf_rank (production), fall back to hp_edit_distance (validation bundle)
    try:
        rank = int(row.get("_conf_rank", 999) or 999)
    except (ValueError, TypeError):
        try:
            rank = int(float(row.get("hp_edit_distance", 9999.0) or 9999.0))
        except (ValueError, TypeError):
            rank = 9999
    ca = row.get("correction_applied", "none") or "none"
    return (aligner, is_winner, rank, ca)


def _best_remaining_tuple(tuples: List[tuple], dropped: str) -> Optional[tuple]:
    """Return best remaining candidate tuple after dropping one aligner."""
    candidates = [t for t in tuples if t[_LOO_ALIGNER] != dropped]
    if not candidates:
        return None
    candidates.sort(key=lambda t: (t[_LOO_RANK], t[_LOO_ALIGNER]))
    return candidates[0]


def run_leave_one_out(
    aligner_summary_paths: List[str],
    chunk_size: int = 50_000,
) -> Tuple[Dict, Dict]:
    """
    Simulate leaving out each aligner and count how many reads change winner.

    Streams aligner_summary files in read-grouped chunks to stay memory-bounded.
    Stores compact tuples (not full dicts) to minimize per-read memory.

    Returns:
        loo_summary: {dropped_aligner: {"n_total", "n_changed", "n_lost_all", "n_lost_correction"}}
        loo_detail:  {(dropped, new_winner, new_correction_applied): count}
    """
    loo_summary: Dict[str, Dict[str, int]] = {}
    loo_detail: Dict[Tuple, int] = defaultdict(int)

    def _process_chunk(chunk: Dict[str, List[tuple]]) -> None:
        all_aligners: set = set()
        for tuples in chunk.values():
            for t in tuples:
                all_aligners.add(t[_LOO_ALIGNER])

        for dropped in all_aligners:
            if dropped not in loo_summary:
                loo_summary[dropped] = {
                    "n_total": 0, "n_changed": 0,
                    "n_lost_all": 0, "n_lost_correction": 0,
                }
            stats = loo_summary[dropped]

            for read_id, tuples in chunk.items():
                stats["n_total"] += 1
                winners = [t for t in tuples if t[_LOO_IS_WINNER]]
                if not winners:
                    continue
                orig = winners[0]
                orig_aligner = orig[_LOO_ALIGNER]
                orig_ca = orig[_LOO_CA]

                if orig_aligner != dropped:
                    loo_detail[(dropped, orig_aligner, orig_ca)] += 1
                    continue

                new_winner = _best_remaining_tuple(tuples, dropped)
                if new_winner is None:
                    stats["n_lost_all"] += 1
                    stats["n_changed"] += 1
                    loo_detail[(dropped, "__no_aligner__", orig_ca)] += 1
                    continue

                new_aligner = new_winner[_LOO_ALIGNER]
                new_ca = new_winner[_LOO_CA]
                stats["n_changed"] += 1
                if orig_ca != "none" and new_ca == "none":
                    stats["n_lost_correction"] += 1
                loo_detail[(dropped, new_aligner, new_ca)] += 1

    # Stream all paths, buffering up to chunk_size unique reads at a time.
    pending: Dict[str, List[tuple]] = {}
    n_reads_pending = 0

    for path in aligner_summary_paths:
        print(f"[INFO]   LOO: {os.path.basename(os.path.dirname(os.path.dirname(path)))}/{os.path.basename(path)}", file=sys.stderr)
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rid = row["read_id"]
                if rid not in pending:
                    n_reads_pending += 1
                    pending[rid] = []
                pending[rid].append(_row_to_tuple(row))

                if n_reads_pending >= chunk_size:
                    _process_chunk(pending)
                    pending = {}
                    n_reads_pending = 0

    if pending:
        _process_chunk(pending)

    return loo_summary, dict(loo_detail)


# ---------------------------------------------------------------------------
# Novel junction analysis
# ---------------------------------------------------------------------------

def load_annotated_junctions(annotation_path: str) -> set:
    """Load annotated intron coordinates from a GTF file."""
    junctions: set = set()
    with open(annotation_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "intron":
                continue
            chrom, start, end, strand = (
                parts[0], int(parts[3]) - 1, int(parts[4]), parts[6]
            )
            junctions.add((chrom, start, end, strand))
    return junctions


def summarize_junctions(
    junc_per_aligner: Dict[str, Dict[Tuple, int]],
    annotated_junctions: set,
) -> List[dict]:
    results = []
    for aligner in sorted(junc_per_aligner.keys()):
        junc_counts = junc_per_aligner[aligner]
        unique = len(junc_counts)
        total_calls = sum(junc_counts.values())
        annotated = sum(1 for j in junc_counts if j in annotated_junctions)
        novel = unique - annotated
        results.append({
            "aligner": aligner,
            "n_junction_calls": total_calls,
            "n_unique_junctions": unique,
            "n_annotated": annotated,
            "n_novel": novel,
            "pct_annotated": round(100 * annotated / unique, 1) if unique else 0,
            "pct_novel": round(100 * novel / unique, 1) if unique else 0,
        })
    return results


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_tsv(path: str, rows: List[dict], fieldnames: Optional[List[str]] = None) -> None:
    if not rows:
        with open(path, "w") as fh:
            fh.write("(no data)\n")
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_cross_tab(
    path: str,
    cross: Dict[Tuple, int],
    totals: Dict[str, int],
    all_corrections: List[str],
    all_aligners: List[str],
) -> None:
    rows = []
    for aligner in sorted(all_aligners):
        row = {"winning_aligner": aligner, "TOTAL": totals.get(aligner, 0)}
        for ca in all_corrections:
            row[ca] = cross.get((aligner, ca), 0)
        rows.append(row)
    total_row = {"winning_aligner": "TOTAL", "TOTAL": sum(totals.values())}
    for ca in all_corrections:
        total_row[ca] = sum(cross.get((a, ca), 0) for a in all_aligners)
    rows.append(total_row)

    fieldnames = ["winning_aligner", "TOTAL"] + all_corrections
    write_tsv(path, rows, fieldnames=fieldnames)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--samples", nargs="+", required=True, metavar="DIR",
        help="Sample directories (each containing corrected_reads.tsv and chunk aligner_summary.tsv files)",
    )
    p.add_argument(
        "--output", "-o", required=True, metavar="DIR",
        help="Output directory for TSV reports",
    )
    p.add_argument(
        "--annotation", metavar="GTF",
        help="Annotation GTF for novel junction analysis (optional)",
    )
    p.add_argument(
        "--chunk-size", type=int, default=50_000, metavar="N",
        help="Reads per LOO simulation chunk (default: 50000)",
    )
    p.add_argument(
        "--skip-loo", action="store_true",
        help="Skip leave-one-out simulation (faster if only win-rate needed)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    annotated_junctions: Optional[set] = None
    if args.annotation:
        print(f"[INFO] Loading annotation from {args.annotation}...", file=sys.stderr)
        annotated_junctions = load_annotated_junctions(args.annotation)
        print(
            f"[INFO]   {len(annotated_junctions)} annotated introns loaded.", file=sys.stderr
        )

    # Streaming accumulators — O(unique combinations), not O(n reads)
    cross: Dict[Tuple, int] = defaultdict(int)
    by_sample: Dict[Tuple, int] = defaultdict(int)
    totals: Dict[str, int] = defaultdict(int)
    junc_per_aligner: Optional[Dict] = {} if annotated_junctions is not None else None

    all_aligner_summary_paths: List[str] = []
    all_samples: List[str] = []

    for sample_dir in args.samples:
        sample_name = os.path.basename(sample_dir.rstrip("/"))
        cr_path = find_corrected_reads(sample_dir)
        if cr_path is None:
            print(f"[WARN] No corrected_reads.tsv in {sample_dir}, skipping", file=sys.stderr)
            continue

        print(f"[INFO] {sample_name}: {cr_path}", file=sys.stderr)
        accumulate_win_rates(
            cr_path, sample_name, cross, by_sample, totals,
            junc_per_aligner, annotated_junctions,
        )
        all_samples.append(sample_name)

        aligner_paths = find_aligner_summary_files(sample_dir)
        if not aligner_paths:
            print(f"[WARN] No aligner_summary.tsv for {sample_name}", file=sys.stderr)
        else:
            print(f"[INFO]   {len(aligner_paths)} aligner_summary chunk(s)", file=sys.stderr)
        all_aligner_summary_paths.extend(aligner_paths)

    if not totals:
        print("[ERROR] No data loaded.", file=sys.stderr)
        sys.exit(1)

    # --- Win-rate cross-tab ---
    print("[INFO] Writing win-rate tables...", file=sys.stderr)
    all_aligners = sorted(totals.keys())
    all_corrections = sorted(set(ca for (_, ca) in cross.keys()))

    cross_path = os.path.join(args.output, "win_rate_by_category.tsv")
    write_cross_tab(cross_path, dict(cross), dict(totals), all_corrections, all_aligners)
    print(f"[INFO] -> {cross_path}", file=sys.stderr)

    sample_rows = []
    for sample in all_samples:
        total_s = sum(by_sample.get((sample, a), 0) for a in all_aligners)
        for aligner in all_aligners:
            count = by_sample.get((sample, aligner), 0)
            sample_rows.append({
                "sample": sample,
                "winning_aligner": aligner,
                "n_reads": count,
                "pct": round(100 * count / total_s, 2) if total_s else 0,
            })
    sample_path = os.path.join(args.output, "win_rate_by_sample.tsv")
    write_tsv(sample_path, sample_rows)
    print(f"[INFO] -> {sample_path}", file=sys.stderr)

    # Console summary
    n_total_reads = sum(totals.values())
    print(f"\n=== Win Rate by Correction Category ({n_total_reads:,} reads) ===")
    cols = all_corrections[:6]  # cap at 6 cols for readability
    hdr = f"{'aligner':<15} {'TOTAL':>9}"
    for ca in cols:
        hdr += f"  {ca[:18]:>18}"
    if len(all_corrections) > 6:
        hdr += f"  {'(+others)':>10}"
    print(hdr)
    print("-" * len(hdr))
    for aligner in all_aligners:
        t = totals.get(aligner, 0)
        pct_t = 100 * t / n_total_reads if n_total_reads else 0
        row_str = f"{aligner:<15} {t:>7} ({pct_t:4.1f}%)"
        for ca in cols:
            cnt = cross.get((aligner, ca), 0)
            pct = 100 * cnt / t if t else 0
            row_str += f"  {cnt:>10} ({pct:4.1f}%)"
        print(row_str)

    # --- Leave-one-out ---
    if not args.skip_loo:
        if not all_aligner_summary_paths:
            print("[WARN] No aligner_summary files; skipping LOO.", file=sys.stderr)
        else:
            print(
                f"\n[INFO] LOO simulation: {len(all_aligner_summary_paths)} chunk file(s)...",
                file=sys.stderr,
            )
            loo_summary, loo_detail = run_leave_one_out(
                all_aligner_summary_paths, chunk_size=args.chunk_size
            )

            loo_rows = []
            for dropped, stats in sorted(loo_summary.items()):
                n_total = stats["n_total"]
                n_changed = stats["n_changed"]
                loo_rows.append({
                    "dropped_aligner": dropped,
                    "n_total_reads": n_total,
                    "n_reads_changed_winner": n_changed,
                    "pct_changed": round(100 * n_changed / n_total, 1) if n_total else 0,
                    "n_lost_all_aligners": stats["n_lost_all"],
                    "n_lost_correction": stats["n_lost_correction"],
                })
            write_tsv(os.path.join(args.output, "leave_one_out_summary.tsv"), loo_rows)
            print(f"[INFO] -> {os.path.join(args.output, 'leave_one_out_summary.tsv')}", file=sys.stderr)

            detail_rows = []
            for (dropped, new_winner, new_ca), count in sorted(loo_detail.items()):
                detail_rows.append({
                    "dropped_aligner": dropped,
                    "new_winning_aligner": new_winner,
                    "new_correction_applied": new_ca,
                    "n_reads": count,
                })
            write_tsv(os.path.join(args.output, "leave_one_out_detail.tsv"), detail_rows)
            print(f"[INFO] -> {os.path.join(args.output, 'leave_one_out_detail.tsv')}", file=sys.stderr)

            print("\n=== Leave-One-Out Simulation ===")
            print(f"{'dropped':<15} {'n_total':>10} {'n_changed':>10} {'%changed':>9} "
                  f"{'n_lost_correction':>18}")
            print("-" * 68)
            for r in loo_rows:
                print(
                    f"{r['dropped_aligner']:<15} {r['n_total_reads']:>10,} "
                    f"{r['n_reads_changed_winner']:>10,} {r['pct_changed']:>9} "
                    f"{r['n_lost_correction']:>18,}"
                )

    # --- Novel junctions ---
    if junc_per_aligner is not None and annotated_junctions is not None:
        junc_rows = summarize_junctions(junc_per_aligner, annotated_junctions)
        write_tsv(os.path.join(args.output, "novel_junction_summary.tsv"), junc_rows)
        print(f"[INFO] -> {os.path.join(args.output, 'novel_junction_summary.tsv')}", file=sys.stderr)

        print("\n=== Novel Junction Contribution ===")
        print(f"{'aligner':<15} {'n_calls':>8} {'n_unique':>9} {'annotated':>10} "
              f"{'novel':>7} {'%novel':>7}")
        print("-" * 62)
        for r in junc_rows:
            print(
                f"{r['aligner']:<15} {r['n_junction_calls']:>8,} "
                f"{r['n_unique_junctions']:>9,} {r['n_annotated']:>10,} "
                f"{r['n_novel']:>7,} {r['pct_novel']:>6}%"
            )

    print(f"\n[INFO] Done. Outputs in {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
