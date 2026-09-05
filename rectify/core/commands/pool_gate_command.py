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
        PoolGateConfig, find_bundled_background_sv_bed,
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
    p.add_argument('--max-intron', type=int, default=None, metavar='BP',
                   help='Length pre-gate: junctions longer than this demote '
                        'before the verdict (planning/684c). Default: derived '
                        'from --annotation as 2x the longest annotated intron '
                        '(yeast derives 5000)')

    from rectify.data import add_organism_args
    add_organism_args(p)
    return p
