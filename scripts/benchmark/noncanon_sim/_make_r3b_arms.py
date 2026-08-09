"""r3b panel: emit margin-matched arm-E (guard-only) + faithful arm-Ff/arm-Ffg so the
guard comparison isolates the del_open term.  Reuses aligned.sorted.bam + one pool.
  arm_E_m3p0.bam — motif_blind + hp_drift_margin=3.0 (no del_open)   [guard baseline]
  arm_Ff.bam     — faithful del_open (lam), NO guard
  arm_Ffg.bam    — faithful del_open (lam) + hp_drift_margin=3.0
Usage: _make_r3b_arms.py <work_dir> <lam>
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

wd = Path(sys.argv[1]); lam = float(sys.argv[2]); M = 3.0
aligned = wd / 'aligned.sorted.bam'
genome = load_sim_genome(wd / 'sim_ref.fa')
annot = load_annotated_junctions(str(DEFAULT_GFF))
pool, annot_set = build_junction_pool([str(aligned)], annot)
base = dict(input_bam=str(aligned), aligner_bams=[str(aligned)], annotated_junctions=annot,
            genome=genome, penalty_table_path=None, prebuilt_junction_pool=pool,
            prebuilt_annotated_set=annot_set, sort_and_index=True, n_workers=4, motif_blind=True)
PT = str(DEFAULT_PENALTY_TABLE)

s = refine_bam_junctions(output_bam=str(wd/'arm_E_m3p0.bam'), hp_drift_margin=M, **base)
print(f"arm_E_m3p0 (guard={M}): refined {s.get('refined')} of {s.get('n_op_reads')}", flush=True)
s = refine_bam_junctions(output_bam=str(wd/'arm_Ff.bam'), hp_drift_margin=0.0,
                         del_open_penalty_table_path=PT, del_open_lam=lam, **base)
print(f"arm_Ff (lam={lam}, no guard): refined {s.get('refined')} of {s.get('n_op_reads')}", flush=True)
s = refine_bam_junctions(output_bam=str(wd/'arm_Ffg.bam'), hp_drift_margin=M,
                         del_open_penalty_table_path=PT, del_open_lam=lam, **base)
print(f"arm_Ffg (lam={lam}, guard={M}): refined {s.get('refined')} of {s.get('n_op_reads')}", flush=True)
print("R3B ARMS DONE")
