#!/usr/bin/env python3
"""Non-canonical preservation control for the resolver's v2 junction re-arbitration.

Question (the discovery-FDR gap in planning/644b's verdict criteria): does the
rearb pass — in particular the splicing-grammar tiebreak, which snaps a
non-canonical-class incumbent to a canonical-class alternative at
match-or-beat ED — FLATTEN genuine non-canonical junctions when a canonical
decoy sits a few bp away? The noncanon_sim mixture panel is by-construction
ground truth for exactly this (WT canonical + cryptic non-canonical reads per
locus, canonical decoy near the cryptic acceptor).

Arms (same minimap2 BAM in):
  mm2    — raw minimap2 (baseline)
  noarb  — resolver, arb_enable=False   (clip resolution only)
  v2     — resolver, arb_enable=True    (clips + re-arbitration + grammar)

Flattening signal = cryptic-cell recovery drop mm2 -> v2 (noarb isolates the
clip stage); canonical/INTRONFREE cells are the do-no-harm controls.
"""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

TREE = Path('/private/tmp/claude-501/-Users-kevinroy-work-rectify/'
            '12349562-4dab-426f-bb0c-0f5ca6d90c91/scratchpad/rearb_control')
WORK = TREE.parent / 'rearb_control_out'
PY = sys.executable
SIM = TREE / 'scripts' / 'benchmark' / 'noncanon_sim'

sys.path.insert(0, str(TREE))


def sh(cmd, **kw):
    print('+', ' '.join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def import_from(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    WORK.mkdir(exist_ok=True)

    # 1+2. panel + reads (fallback injector; WT+cryptic mixture, decoy present)
    if not (WORK / 'sim_ref.fa').exists():
        sh([PY, SIM / 'build_panel.py', '--panel', 'smoke', '--out-dir', WORK])
    if not (WORK / 'reads.fastq').exists():
        sh([PY, SIM / 'gen_reads.py', '--out-dir', WORK, '--error-regime', 'mid'])

    # 3. minimap2 (harness helper: annotation-blind, -uf, G=5000)
    run_arms = import_from(SIM / 'run_arms.py', 'run_arms_mod')
    mm2_bam = WORK / 'arm_mm2.bam'
    if not mm2_bam.exists():
        run_arms.align_reads(WORK / 'reads.fastq', WORK / 'sim_ref.fa', mm2_bam)

    # 4. resolver arms
    from rectify.core.align.overhang_resolver import (ResolverConfig,
                                                      run_overhang_resolver)
    for name, arb in (('noarb', False), ('v2', True)):
        out = WORK / f'arm_{name}.bam'
        if out.exists():
            continue
        cfg = ResolverConfig(arb_enable=arb)
        run_overhang_resolver(str(mm2_bam), str(WORK / 'sim_ref.fa'), str(out),
                              config=cfg)
        stats = run_overhang_resolver.last_stats.as_dict()
        (WORK / f'stats_{name}.json').write_text(json.dumps(stats, indent=1))
        print(f'[{name}] stats: ' + json.dumps(stats), flush=True)
        subprocess.run(['samtools', 'index', str(out)], check=False)

    # 5. score all three arms against the by-construction truth
    sh([PY, SIM / 'score_trade.py', '--work-dir', WORK,
        '--arm', f'mm2={mm2_bam}',
        '--arm', f'noarb={WORK / "arm_noarb.bam"}',
        '--arm', f'v2={WORK / "arm_v2.bam"}',
        '--out', WORK / 'trade_curve.json'])

    print('\nDONE ->', WORK / 'trade_curve.json')


if __name__ == '__main__':
    main()
