import sys; sys.path.insert(0,'.'); sys.path.insert(0,'../../..')
import multiprocessing as _mp
try: _mp.set_start_method("fork")
except RuntimeError: pass
from pathlib import Path
from run_arms import load_sim_genome, DEFAULT_GFF, build_junction_pool
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.core.splice.junction_refiner import refine_bam_junctions
wd = Path(sys.argv[1]); hpd = float(sys.argv[2])
aligned = wd/'aligned.sorted.bam'
genome = load_sim_genome(wd/'sim_ref.fa')
annot = load_annotated_junctions(str(DEFAULT_GFF))
pool, annot_set = build_junction_pool([str(aligned)], annot)
stats = refine_bam_junctions(
    input_bam=str(aligned), output_bam=str(wd/'arm_E.bam'),
    aligner_bams=[str(aligned)], annotated_junctions=annot, genome=genome,
    penalty_table_path=None, prebuilt_junction_pool=pool, prebuilt_annotated_set=annot_set,
    sort_and_index=True, n_workers=4, motif_blind=True, hp_drift_margin=hpd)
print('arm_E refined:', stats.get('refined'), 'of', stats.get('n_op_reads'))
