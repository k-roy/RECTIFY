"""
rectify triage — classify consensus alignments; re-align the triaged minority.

The CLI face of the consensus-triage layer (``rectify.core.consensus.triage``;
design: dev/DENOVO_ALIGNER_REASSESSMENT_20260809.md §3). Classification is
read-evidence-only; the junction re-align leg is motif-blind refinement with
strict hp_ed re-entry (a candidate never auto-replaces the aligner's placement).

Outputs, under ``-o``:
  triage.tsv           per-read label + reasons + signals + realign outcome
  <input>.triaged.bam  the input with accepted re-alignments swapped in
                       (only with --realign; sorted + indexed)
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def run_triage(args: argparse.Namespace) -> int:
    from rectify.data import resolve_reference_paths
    from rectify.utils.genome import load_genome
    from rectify.core.consensus.consensus import load_annotated_junctions
    from rectify.core.consensus.triage import TriagePolicy, triage_realign_bam

    resolve_reference_paths(args, require_genome=True)
    if not getattr(args, 'annotation', None):
        print("ERROR: rectify triage needs --annotation (or an --organism/--Scer "
              "bundle providing one) for junction annotation status.", file=sys.stderr)
        return 1

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    genome = load_genome(Path(args.genome))
    annotated = load_annotated_junctions(str(args.annotation))

    policy = TriagePolicy(
        max_junction_proximal_errors=args.max_junction_proximal_errors,
        max_clip_5p=args.max_clip_5p,
        max_clip_3p=args.max_clip_3p,
        triage_unannotated_junctions=not args.no_triage_unannotated,
        clip_legs_enable=getattr(args, 'clip_legs', False),
    )

    penalty_table = None
    if getattr(args, 'penalty_table', None):
        from rectify.core.splice.hp_penalty import HpPenaltyTable
        penalty_table = HpPenaltyTable.from_tsv(args.penalty_table)

    in_bam = Path(args.input)
    out_bam = out_dir / (in_bam.stem + '.triaged.bam')

    rows, stats = triage_realign_bam(
        str(in_bam), str(out_bam),
        genome=genome,
        annotated_junctions=annotated,
        policy=policy,
        pool_bams=list(args.pool_bams) if args.pool_bams else None,
        original_bams=(list(args.original_bams)
                       if getattr(args, 'original_bams', None) else None),
        penalty_table=penalty_table,
        realign=args.realign,
    )

    tsv_path = out_dir / 'triage.tsv'
    fieldnames = ['read_id', 'label', 'reasons', 'junction_proximal_errors',
                  'clip_5p', 'clip_3p', 'n_junctions', 'n_unannotated',
                  'realigned', 'accepted', 'reverted_to_original']
    with open(tsv_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        for row in rows:
            w.writerow(row)

    print(f"Triage: {stats['classified']} classified — "
          f"{stats['high_confidence']} high-confidence (bypass), "
          f"{stats['triaged']} triaged "
          f"({stats['junction_leg']} junction-leg)")
    if args.realign:
        print(f"Re-align: {stats['realigned']} moved by the motif-blind refiner; "
              f"{stats['accepted']} accepted by hp_ed re-entry -> {out_bam}")
        if stats['orig_leg']:
            print(f"Pre-correction candidates: {stats['orig_leg']} offered, "
                  f"{stats['orig_proposed']} differed from the incumbent, "
                  f"{stats['orig_accepted']} reverted to the original "
                  f"alignment by hp_ed re-entry")
            if stats['orig_skipped_unknown_chrom']:
                print(f"  WARNING: {stats['orig_skipped_unknown_chrom']} "
                      f"pre-correction record(s) skipped — reference name not "
                      f"in the input BAM header (chromosome naming mismatch)",
                      file=sys.stderr)
    else:
        print("Re-align: skipped (--no-realign)")
    print(f"Per-read table: {tsv_path}")
    return 0


def create_triage_parser(subparsers) -> argparse.ArgumentParser:
    p = subparsers.add_parser(
        'triage',
        help='Classify consensus alignments (high-confidence bypass vs triage); '
             'motif-blind re-align the triaged junction reads with hp_ed re-entry',
    )
    p.add_argument('input', help='Consensus / corrected BAM to triage')
    p.add_argument('-o', '--output', required=True, help='Output directory')
    p.add_argument('--genome', help='Reference FASTA (or use an organism bundle)')
    p.add_argument('--annotation', help='GFF/GTF annotation (or organism bundle)')
    p.add_argument('--pool-bams', nargs='+', default=None,
                   help='BAMs whose junction evidence seeds the re-align pool '
                        '(default: the input BAM itself). Pass the per-aligner '
                        'panel BAMs when available — the pool must reflect FULL '
                        'evidence, never the triaged subset.')
    p.add_argument('--original-bams', nargs='+', default=None,
                   help='PRE-CORRECTION BAM(s) — the aligner output before '
                        '`rectify correct`. Each triaged read\'s original '
                        'alignment is offered as a CANDIDATE against the same '
                        'strict hp_ed re-entry arbiter, so Station B can undo '
                        '2F/2H damage when the original strictly wins. Without '
                        'this the only proposer is Module 2H itself, which '
                        're-derives its own fixed point and can never offer '
                        'back the placement it moved away from.')
    p.add_argument('--penalty-table', default=None,
                   help='Empirical HP penalty table TSV for the hp_ed re-entry '
                        'metric (default: flat costs)')
    p.add_argument('--realign', dest='realign', action='store_true', default=True,
                   help='Run the junction re-align leg (default)')
    p.add_argument('--no-realign', dest='realign', action='store_false',
                   help='Classify only; write triage.tsv and skip re-alignment')
    p.add_argument('--max-junction-proximal-errors', type=float, default=1.0,
                   help='HP-weighted junction-proximal error budget for the '
                        'high-confidence bypass (default 1.0)')
    p.add_argument('--max-clip-5p', type=int, default=30,
                   help='5\' soft-clip bases tolerated before triage (default 30)')
    p.add_argument('--max-clip-3p', type=int, default=30,
                   help='3\' soft-clip bases tolerated before triage (default 30)')
    p.add_argument('--no-triage-unannotated', action='store_true', default=False,
                   help='Do not triage reads solely for carrying an unannotated junction')
    p.add_argument('--clip-legs', dest='clip_legs', action='store_true',
                   default=False,
                   help='Enable the clip legs: terminal clips to the overhang '
                        'resolver, 5\' clips to Cat3 rescue, one refusal '
                        'discipline, hp_ed re-entry as the arbiter (default off)')

    from rectify.data import add_organism_args
    add_organism_args(p)
    return p
