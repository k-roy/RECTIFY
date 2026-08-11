"""Arm-Ff: the FAITHFUL C1 one-time gap-OPEN discount (run-latch _score_hp_affine_del,
staircase fixed) wired into the re-placer.  Emits BOTH:
  arm_Ff.bam   — standalone (motif_blind, NO guard)
  arm_Ffg.bam  — guard-combined (motif_blind + hp_drift_margin=3.0)
Reuses the existing aligned.sorted.bam + one shared junction pool (single arm at a
time, n_workers=4).  lambda in the re-placer's own del-cost units (del_hp=0.5).

Usage: _make_arm_ff.py <work_dir> <lam>   (e.g. mix_fair_out 0.1)
"""
import sys; sys.path.insert(0, '.'); sys.path.insert(0, '../../..')
import multiprocessing as _mp
try:
    _mp.set_start_method("fork")
except RuntimeError:
    pass
from pathlib import Path
from run_arms import load_sim_genome, DEFAULT_GFF, DEFAULT_PENALTY_TABLE, build_junction_pool
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.core.splice.junction_refiner import refine_bam_junctions

wd = Path(sys.argv[1])
lam = float(sys.argv[2])
guard_margin = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
aligned = wd / 'aligned.sorted.bam'
genome = load_sim_genome(wd / 'sim_ref.fa')
annot = load_annotated_junctions(str(DEFAULT_GFF))
pool, annot_set = build_junction_pool([str(aligned)], annot)   # once, shared

common = dict(
    input_bam=str(aligned), aligner_bams=[str(aligned)], annotated_junctions=annot,
    genome=genome, penalty_table_path=None,           # flat legacy del/ins (arm-B base)
    prebuilt_junction_pool=pool, prebuilt_annotated_set=annot_set,
    sort_and_index=True, n_workers=4, motif_blind=True,
    del_open_penalty_table_path=str(DEFAULT_PENALTY_TABLE),  # rate_mean carrier
    del_open_lam=lam,
)

# standalone (no guard)
s1 = refine_bam_junctions(output_bam=str(wd / 'arm_Ff.bam'), hp_drift_margin=0.0, **common)
print(f"arm_Ff  (lam={lam}, no guard):     refined {s1.get('refined')} of {s1.get('n_op_reads')}",
      flush=True)
# guard-combined
s2 = refine_bam_junctions(output_bam=str(wd / 'arm_Ffg.bam'), hp_drift_margin=guard_margin, **common)
print(f"arm_Ffg (lam={lam}, guard={guard_margin}): refined {s2.get('refined')} of {s2.get('n_op_reads')}",
      flush=True)
print("ARM-Ff DONE")
