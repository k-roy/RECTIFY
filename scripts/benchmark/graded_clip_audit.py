#!/usr/bin/env python3
"""Graded 3'-clip penalty audit — Option-B EER-ED winner-reshuffle on real DRS data.

Background (see rectify/data/validation/CASE_STUDIES.md, cat2_plus_1): the
HP-aware winner edit distance charges a FLAT 1.0/base for a 3' soft/hard-clip.
That penalises an aligner for HONESTLY soft-clipping a poly-A tail, while an
aligner that force-aligns the tail (insertions+mismatches) pays less and wins.
Proposed fix: discount the clip penalty by the estimated poly-A-tail length
inside the clip.

The naive "whole-clip A-fraction" discount over-fires on long non-tail clips
that happen to be moderately A-rich (validation read 00a1e01e: a 76-bp clip,
36% A overall, 25% A in its terminal window — a discarded 2nd-exon segment, NOT
a tail). The discriminator is the **terminal** stop-base run, not the whole-clip
average. This script estimates the tail length from the 3' terminus inward
(expanding while the cumulative stop-base fraction stays >= TAIL_FRAC) and only
forgives that many bases; the rest of the clip keeps full penalty.

It recomputes, per read, the winner under (a) the current flat-clip ED and
(b) the graded ED, then reports the reshuffle, classifies each flip
(tail-flip = the newly-winning aligner's clip is a genuine tail; suspicious =
otherwise), and runs a REGRESSION CHECK: any read whose flip is driven by
discounting a LOW-terminal-richness clip is flagged (the 00a1e01e failure mode).

Usage:
    python graded_clip_audit.py --per-aligner-dir DIR --genome GENOME.fsa \
        [--aligners minimap2 gapmm2 mapPacBio deSALT uLTRA] \
        [--variant softclipped] [--tail-frac 0.6] [--term-window 12] \
        [--out report.tsv]

DIR must contain <aligner>.<variant>.bam (corrected per-aligner BAMs, e.g. from
`rectify correct --write-softclipped-bam`). The softclipped variant is REQUIRED
(the trimmed/hard-clipped variant drops the clipped bases, so tail composition
cannot be assessed).
"""
import argparse
import sys
from pathlib import Path

import pysam

# Reuse the production HP penalty table + HP-run helper so del/ins costs match
# the real winner-selection scorer exactly. Import from the canonical hp_penalty
# module (junction_refiner re-exports them, but older deployed rectify trees may
# not) with a fallback.
try:
    from rectify.core.splice.hp_penalty import HpPenaltyTable, _hp_run_length
except ImportError:  # pragma: no cover - older trees re-export via junction_refiner
    from rectify.core.splice.junction_refiner import HpPenaltyTable, _hp_run_length


def estimate_tail_len(clip: str, stop: str, tail_frac: float, term_window: int) -> int:
    """Estimate the poly-A-tail length within a 3' soft-clip.

    Expand a window from the 3' terminus (end of ``clip``) inward while the
    cumulative stop-base fraction stays >= ``tail_frac``. Returns the length of
    the longest such terminal run (capped at ``term_window``). A degenerate
    AT-rich tail (e.g. cat2_plus_1 ``AAATAAAAT``, 78% A) clears tail_frac=0.6;
    a non-tail clip whose terminus is genomic (00a1e01e ``...TCCATACATCCT``,
    25% A) does not, so almost nothing is forgiven.
    """
    if not clip:
        return 0
    best = 0
    stop_n = 0
    rev = clip[::-1]  # iterate from the 3' terminus inward
    for i, base in enumerate(rev[:term_window], start=1):
        if base.upper() == stop:
            stop_n += 1
        if stop_n / i >= tail_frac:
            best = i
    return best


def cigar_ed(read, genome_seq, pt, *, graded, tail_frac, term_window):
    """HP-aware edit distance for one corrected read. ``graded=False`` is the
    current production model (flat 1.0/clip-base); ``graded=True`` forgives the
    estimated poly-A-tail portion of the trailing 3' soft-clip.

    Returns (ed, tail_len, clip_len) — tail_len/clip_len describe the trailing
    3' soft-clip so the caller can classify flips.
    """
    ct = read.cigartuples
    if ct is None:
        return None
    seq = read.query_sequence
    if seq is None or read.is_secondary or read.is_supplementary:
        return None
    stop = "A" if not read.is_reverse else "T"  # plus->A tail, minus->T (BAM SEQ is RC)
    rp = read.reference_start
    qp = 0
    total = 0.0
    n = len(ct)
    tail_len = 0
    clip_len = 0
    for ci, (op, length) in enumerate(ct):
        if op == 7:                       # =
            rp += length; qp += length
        elif op == 8:                     # X
            total += length; rp += length; qp += length
        elif op == 0:                     # M (compare bases)
            total += sum(1 for k in range(length)
                         if genome_seq[rp + k] != seq[qp + k].upper())
            rp += length; qp += length
        elif op == 2:                     # D — HP-aware
            total += sum(pt.del_cost(_hp_run_length(genome_seq, rp + i),
                                     genome_seq[rp + i]) for i in range(length))
            rp += length
        elif op == 1:                     # I — HP-aware
            total += length * pt.ins_cost(_hp_run_length(genome_seq, rp),
                                          genome_seq[rp])
            qp += length
        elif op == 4:                     # S (soft-clip)
            is_trailing_3p = (ci == n - 1)
            if is_trailing_3p:
                clip_len = length
                clip = seq[qp:qp + length]
                if graded:
                    tail_len = estimate_tail_len(clip, stop, tail_frac, term_window)
                    total += (length - tail_len) * 1.0   # forgive the tail portion
                else:
                    total += length * 1.0
            else:
                total += length * 1.0     # 5' clip: always full penalty
            qp += length
        elif op == 5:                     # H (hard-clip)
            total += length * 1.0
    return total, tail_len, clip_len


def winner(eds, idx):
    cand = {a: v[idx][0] for a, v in eds.items() if v[idx] is not None}
    return min(cand, key=cand.get) if cand else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-aligner-dir", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--aligners", nargs="+",
                    default=["minimap2", "gapmm2", "mapPacBio", "deSALT", "uLTRA"])
    ap.add_argument("--variant", default="softclipped",
                    help="corrected-BAM variant (must retain soft-clip bases)")
    ap.add_argument("--tail-frac", type=float, default=0.6)
    ap.add_argument("--term-window", type=int, default=12)
    ap.add_argument("--penalty-tsv", default=None,
                    help="penalty_scores.tsv (default: bundled with --genome's tree)")
    ap.add_argument("--str-penalty-tsv", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    pdir = Path(args.per_aligner_dir)
    # Resolve penalty tables (default to the genome tree's penalty_tables/).
    if args.penalty_tsv:
        pen, strpen = args.penalty_tsv, args.str_penalty_tsv
    else:
        pt_dir = Path(args.genome).resolve().parent / "penalty_tables"
        pen = str(pt_dir / "penalty_scores.tsv")
        strpen = str(pt_dir / "str_penalty_scores.tsv")
    pt = HpPenaltyTable.from_tsv(pen, strpen if strpen and Path(strpen).exists() else None)

    fa = pysam.FastaFile(args.genome)
    chrom_seq = {c: fa.fetch(c).upper() for c in fa.references}

    # eds[read_id][aligner] = ((ed_flat, tail, clip), (ed_graded, tail, clip))
    eds = {}
    for al in args.aligners:
        bam_path = pdir / f"{al}.{args.variant}.bam"
        if not bam_path.exists():
            print(f"  ! missing {bam_path}; skipping", file=sys.stderr)
            continue
        with pysam.AlignmentFile(str(bam_path)) as bam:
            for r in bam.fetch():
                cs = chrom_seq.get(r.reference_name)
                if cs is None:
                    continue
                flat = cigar_ed(r, cs, pt, graded=False,
                                tail_frac=args.tail_frac, term_window=args.term_window)
                grad = cigar_ed(r, cs, pt, graded=True,
                                tail_frac=args.tail_frac, term_window=args.term_window)
                if flat is None or grad is None:
                    continue
                eds.setdefault(r.query_name, {})[al] = (flat, grad)

    rows = []
    n_flip = n_tailflip = n_suspicious = 0
    for rid, d in eds.items():
        w0, w1 = winner(d, 0), winner(d, 1)
        if w0 == w1 or w1 is None:
            continue
        n_flip += 1
        # Classify: did the new winner's trailing clip look like a genuine tail?
        _, tail_w1, clip_w1 = d[w1][1]
        tail_share = (tail_w1 / clip_w1) if clip_w1 else 0.0
        kind = "tail-flip" if tail_share >= 0.5 else "SUSPICIOUS"
        if kind == "tail-flip":
            n_tailflip += 1
        else:
            n_suspicious += 1
        rows.append((rid, w0, round(d[w0][0][0], 2), w1, round(d[w1][1][0], 2),
                     clip_w1, tail_w1, round(tail_share, 2), kind))

    rows.sort(key=lambda r: (r[8] != "SUSPICIOUS", r[0]))  # suspicious first
    hdr = ("read_id", "old_winner", "old_flatED", "new_winner", "new_gradedED",
           "newwin_clip", "newwin_tail", "tail_share", "kind")
    out = sys.stdout if not args.out else open(args.out, "w")
    print("\t".join(hdr), file=out)
    for r in rows:
        print("\t".join(str(x) for x in r), file=out)
    if args.out:
        out.close()

    print(f"\n=== {n_flip}/{len(eds)} reads change winner under graded clip penalty "
          f"({n_tailflip} tail-flips, {n_suspicious} SUSPICIOUS) ===", file=sys.stderr)
    print(f"params: tail_frac={args.tail_frac} term_window={args.term_window} "
          f"variant={args.variant}", file=sys.stderr)
    if n_suspicious:
        print("REGRESSION WARNING: SUSPICIOUS flips discount a low-terminal-richness "
              "clip (00a1e01e failure mode) — inspect before adopting.", file=sys.stderr)


if __name__ == "__main__":
    main()
