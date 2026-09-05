"""
rectify pool-gate — Station C v0: the pool-level junction admission REPORT.

Censuses junctions from a consensus/corrected BAM and writes a per-junction
admission table (admit_candidate / review / demote_orthogonal_evidence) from
the measured phase-2 scorer: short-exon-side overhang quality x repeat-context
flags x within-sample recurrence, on the two-track (canonical-in-class vs
non-canonical) discipline. Report-only — annotates, never deletes.

Design + measurements: dev/REALIGNER_LANDSCAPE_AND_PATH_20260809.md §3,
dev/PHASE2_OVERHANG_644H_20260811.md, dev/STATIONC_REPEAT_FLAG_644I_20260811.md.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def run_pool_gate(args: argparse.Namespace) -> int:
    from rectify.data import resolve_reference_paths
    from rectify.utils.genome import load_genome
    from rectify.core.consensus.station_c import (
        PoolGateConfig, derive_pool_gate_max_intron,
        find_bundled_background_sv_bed,
        find_bundled_selfhom_bed, pool_gate,
        write_pool_gate_outputs,
    )

    resolve_reference_paths(args, require_genome=True)
    if not getattr(args, 'annotation', None):
        print("ERROR: rectify pool-gate needs --annotation (or an --organism/"
              "--Scer bundle providing one).", file=sys.stderr)
        return 1

    genome = load_genome(Path(args.genome))

    selfhom = Path(args.selfhom_bed) if args.selfhom_bed else \
        find_bundled_selfhom_bed(Path(args.genome))
    background_sv = Path(args.background_sv_bed) if args.background_sv_bed \
        else find_bundled_background_sv_bed(Path(args.genome))

    # Length pre-gate bound: explicit flag wins; otherwise a high QUANTILE of
    # the annotated intron-length distribution (ISSUE-013). The aligner's rule
    # (2x the single longest intron) saturates the 500,000 clamp on human, so
    # the pre-gate never fired; yeast is unchanged at 5,000 either way.
    _max_intron = getattr(args, 'max_intron', None)
    _derived_by = 'explicit --max-intron'
    if _max_intron is None:
        _max_intron = derive_pool_gate_max_intron(
            Path(args.annotation),
            quantile=args.max_intron_quantile,
            multiplier=args.max_intron_multiplier,
        )
        _derived_by = (f'{args.max_intron_multiplier:g}x p'
                       f'{args.max_intron_quantile * 100:g} of the annotated '
                       f'intron lengths')

    cfg = PoolGateConfig(
        q_canon=args.q_canon,
        q_noncanon=args.q_noncanon,
        min_support=args.min_support,
        max_intron=_max_intron,
        adj_indel_max_ops=args.adj_indel_max_ops,
        adj_indel_max_bp=args.adj_indel_max_bp,
    )
    rows, summary = pool_gate(
        args.input, genome, Path(args.annotation), cfg=cfg, selfhom_bed=selfhom,
        background_sv_bed=background_sv,
    )
    tsv, js = write_pool_gate_outputs(rows, summary, Path(args.output))

    v = summary['verdicts']
    print(f"Station C (pool-gate): {summary['n_junctions_censused']} junctions "
          f"censused ({summary['n_annotated']} annotated); "
          f"admit_candidate={v.get('admit_candidate', 0)} "
          f"review={v.get('review', 0)} "
          f"demote_orthogonal_evidence={v.get('demote_orthogonal_evidence', 0)}")
    if selfhom:
        print(f"  self-homology track: {selfhom}")
    else:
        print("  self-homology track: none (annotation flag only — see "
              "dev/STATIONC_REPEAT_FLAG_644I_20260811.md to build one)")
    if background_sv:
        print(f"  background-SV track: {background_sv}")

    # A dead demotion term must never be mistaken for a clean result. Three of
    # the four are unavailable on human input (no bundled self-homology /
    # background-SV BED, no REPEAT_FEATURE_TYPES equivalent in a GENCODE GTF),
    # and the table used to look identical to a fully-gated yeast run.
    unavailable = summary.get('tracks_unavailable') or []
    if unavailable:
        print(f"  WARNING: {len(unavailable)} of 3 flag tracks UNAVAILABLE for "
              f"this genome ({', '.join(unavailable)}) — those columns read "
              f"'track_unavailable', NOT a clean result. Junctions here are "
              f"gated on canonical_in_class x q_max x support (+ the length "
              f"pre-gate) only.", file=sys.stderr)
    if summary.get('n_annotated_introns_parsed') == 0:
        print(f"  WARNING: 0 annotated introns parsed from {args.annotation} — "
              f"every annotated junction is being reported as a discovery "
              f"candidate. Check the annotation has 'intron' features or exons "
              f"with Parent=/transcript_id attributes.", file=sys.stderr)
    print(f"  length pre-gate: {_max_intron:,} bp ({_derived_by})")
    n_adj = summary.get('n_junctions_with_adjacent_indel')
    if n_adj:
        print(f"  {n_adj} censused junction(s) carry a boundary-adjacent indel "
              f"(adj_indel_l/adj_indel_r columns) — the CIGAR signature of a "
              f"refiner-moved boundary")
    print(f"  table: {tsv}\n  summary: {js}")
    return 0


def create_pool_gate_parser(subparsers) -> argparse.ArgumentParser:
    p = subparsers.add_parser(
        'pool-gate',
        help='Station C (report-only): per-junction admission table from '
             'overhang quality x repeat flags x recurrence, two-track',
    )
    p.add_argument('input', help='Consensus / corrected BAM to census')
    p.add_argument('-o', '--output', required=True,
                   help='Output prefix (writes <prefix>.pool_gate.tsv/.json)')
    p.add_argument('--genome', help='Reference FASTA (or use an organism bundle)')
    p.add_argument('--annotation', help='GFF/GTF annotation (or organism bundle)')
    p.add_argument('--selfhom-bed', default=None,
                   help='Genome self-homology BED (default: auto-detect beside '
                        'the genome / bundled yeast track)')
    p.add_argument('--background-sv-bed', default=None,
                   help='Known background-SV regions of the reference, BED '
                        'col4=name (default: auto-detect *background_sv.bed '
                        'beside the genome / bundled R64 track). Junctions '
                        'overlapping one demote on BOTH tracks — the '
                        'reference is known wrong there (e.g. R64 chrIII '
                        'SRD1 flank-A, yKR888 T2T)')
    p.add_argument('--q-canon', type=float, default=40.0,
                   help='Canonical-track overhang-quality admit threshold, bits '
                        '(default 40 — measured 644h)')
    p.add_argument('--q-noncanon', type=float, default=80.0,
                   help='Non-canonical-track threshold, bits (default 80)')
    p.add_argument('--min-support', type=int, default=2,
                   help='Within-sample read-support gate (default 2)')
    p.add_argument('--adj-indel-max-ops', type=int, default=2, metavar='N',
                   help='Census anchor walk: intervening I/D ops tolerated '
                        'between an N-op and the aligned run that anchors it, '
                        'per side (default 2). Module 2H plants a compensating '
                        'indel beside every junction it moves, so a walk of 0 '
                        'censuses only 16/121 of the junctions RECTIFY created '
                        'on the Sumner panel. On that panel THIS is the binding '
                        'budget: 2 -> 80/121, 3 -> 85, 4 -> 90 (at 30 bp). '
                        'Raising it admits 2F rescued-exon shapes where no '
                        'contiguous anchor exists — see _anchor_run.')
    p.add_argument('--adj-indel-max-bp', type=int, default=30, metavar='BP',
                   help='Census anchor walk: summed length of those '
                        'stepped-over indel ops, per side (default 30). At the '
                        'default op limit, 20 and 30 bp are identical on the '
                        'Sumner panel (80/121); 50 bp gives 87/121.')
    p.add_argument('--max-intron', type=int, default=None, metavar='BP',
                   help='Length pre-gate: junctions longer than this demote '
                        'before the verdict (planning/684c). Default: derived '
                        'from --annotation as 2x the p99.5 of the annotated '
                        'intron lengths (yeast derives 5,000 — identical to '
                        'the historical value; GENCODE v48 chr5 derives '
                        '310,100, where 2x the longest intron would saturate '
                        'the 500,000 ceiling and never fire)')
    p.add_argument('--max-intron-quantile', type=float, default=0.995,
                   metavar='Q',
                   help='Quantile of the annotated intron-length distribution '
                        'the length pre-gate is derived from (default 0.995). '
                        'p99.9 was measured and rejected: on GENCODE chr5 it '
                        'lands at 524,000, back above the ceiling and inert.')
    p.add_argument('--max-intron-multiplier', type=int, default=2, metavar='N',
                   help='Multiplier on that quantile (default 2 — the same 2x '
                        'margin multi_aligner.derive_max_intron applies to the '
                        'longest annotated intron)')

    from rectify.data import add_organism_args
    add_organism_args(p)
    return p
