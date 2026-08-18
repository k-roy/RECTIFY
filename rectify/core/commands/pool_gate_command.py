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
        PoolGateConfig, find_bundled_selfhom_bed, pool_gate,
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

    # Length pre-gate bound: explicit flag wins; otherwise derived from the
    # annotation (2x the longest annotated intron — yeast derives 5000).
    _max_intron = getattr(args, 'max_intron', None)
    if _max_intron is None:
        from rectify.core.align.multi_aligner import derive_max_intron
        _max_intron = derive_max_intron(str(args.annotation))

    cfg = PoolGateConfig(
        q_canon=args.q_canon,
        q_noncanon=args.q_noncanon,
        min_support=args.min_support,
        max_intron=_max_intron,
    )
    rows, summary = pool_gate(
        args.input, genome, Path(args.annotation), cfg=cfg, selfhom_bed=selfhom,
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
    p.add_argument('--q-canon', type=float, default=40.0,
                   help='Canonical-track overhang-quality admit threshold, bits '
                        '(default 40 — measured 644h)')
    p.add_argument('--q-noncanon', type=float, default=80.0,
                   help='Non-canonical-track threshold, bits (default 80)')
    p.add_argument('--min-support', type=int, default=2,
                   help='Within-sample read-support gate (default 2)')
    p.add_argument('--max-intron', type=int, default=None, metavar='BP',
                   help='Length pre-gate: junctions longer than this demote '
                        'before the verdict (planning/684c). Default: derived '
                        'from --annotation as 2x the longest annotated intron '
                        '(yeast derives 5000)')

    from rectify.data import add_organism_args
    add_organism_args(p)
    return p
