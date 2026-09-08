# Read-level review — how to inspect one read

The CLAUDE.md rule ("vet alignments at the INDIVIDUAL-READ level at EVERY pipeline stage") mechanized for a review bundle
(`stock.bam` + one `<sha>.bam` per arm + `manifest.tsv`, e.g. `175_rectified/rbrowse_bundles/review_2f_<sha>/`).

## How to inspect one read

1. **The numbers** — `python3 dev/todo_run_20260905/replay_scripts/inspect_read.py <bundle_dir> <read8> [arm ...]`
   prints, per arm, the CIGAR, soft clips, every exon block with `=`/X/I/D counts and identity, every intron with its motif and
   GENCODE status, the 5′ block base by base, and the 5′ clip's best ungapped placement against annotated exon ends within 5 kb.
2. **The picture** — `python -m rectify.visualize.read_junction --bundle <bundle_dir> --read <read8> --genome <FASTA> --gtf <GTF>
   --out <read8>.png [--all-junctions] [--window 32]`: one panel per junction the arms disagree about (else the 5′-most; `--all-junctions`
   for every junction), each with the reference letters, the annotated exon structures, one letter-level row per arm and a
   numbers-only strip per arm. A `.md` sidecar with the inspector's numbers is written beside the PNG.
3. **The five questions to answer per read** (no thresholds here — they are being modeled):
   - **Identity of the 5′ block** — matched/aligned of the RNA-5′-most block (`5′ 9/32 28 %`): is the rescued exon real sequence or a forced placement?
   - **Junction-proximal clean run** — consecutive matches walking away from the junction on its 5′ and 3′ sides (`run 12|30`): does the read actually anchor the splice site?
   - **Leading I/D = a clip in disguise** — an insertion or deletion at the very start of the 5′ block (`I8 D0` on a 6-nt block) is a soft clip the writer had to pad.
   - **Clip fit to the annotated exon ends within 5 kb** — `clip 63 fit 18/40 @12,975,731 (−1)`: where the clipped letters would have gone, and how well.
   - **MAPQ** — the aligner's own confidence in the body the rescue hangs from (`MQ 1` vs `MQ 60`).
