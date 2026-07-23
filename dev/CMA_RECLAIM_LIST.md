# Reclaim delete list — `aligner_chunks/` (AWAITING KEVIN'S APPROVAL)

**Generated:** 2026-07-23 · **Status:** PROPOSED — nothing deleted. Approve the list before any removal.
**Verification:** H2 `~/m0_scratch/reclaim_verify_all_20260723_1240.txt` (read-only). Method: per sample, for
every aligner, `Σ(raw per-chunk *.<aligner>.bam primary reads) == merged_bams/<sample>.<aligner>.bam primary
reads`, AND `consensus/corrected_consensus.bam` present. SAFE only if ALL aligners match AND consensus exists.

Run: `/u/scratch/k/kevinroy/drs_decay_mutants_rectify_20260702/`

## SAFE to delete — explicit named directories (≈ 101.6 GB)

| # | Directory (delete the whole `aligner_chunks/` tree) | Size |
|--|---|--:|
| 1 | `…/Dcp2_AID_25C_repA/aligner_chunks` | 11 G |
| 2 | `…/Dcp2_AID_25C_repB/aligner_chunks` | 6.5 G |
| 3 | `…/WT_30_1_repA/aligner_chunks` | 8.0 G |
| 4 | `…/WT_W303_AID_control_25C_repA/aligner_chunks` | 4.9 G |
| 5 | `…/WT_W303_AID_control_25C_repB/aligner_chunks` | 15 G |
| 6 | `…/Xrn1_AID_25C_repA/aligner_chunks` | 2.2 G |
| 7 | `…/Xrn1_AID_25C_repB/aligner_chunks` | 28 G |
| 8 | `…/total_mRNA_WT_W303_digested_xrn1_repA/aligner_chunks` | 13 G |
| 9 | `…/total_mRNA_WT_W303_undigested_repA/aligner_chunks` | 13 G |
|   | **TOTAL** | **≈ 101.6 G** |

(`…` = `/u/scratch/k/kevinroy/drs_decay_mutants_rectify_20260702`)

Each `aligner_chunks/` holds raw per-chunk per-aligner BAMs (→ complete in `merged_bams/`), per-chunk
`.junction_refined_*.bam` (→ pool in `junction_pool.pkl`), and per-chunk `.rectified.bam` (→ final
`consensus/corrected_consensus.bam`). All three durable products verified present.

## EXCLUDED — do NOT delete
- `…/Xrn1_AID_25C_repB_4aln/aligner_chunks` (3.6 G) — **UNSAFE**: no `consensus/corrected_consensus.bam`,
  and deSALT/gapmm2/uLTRA have **no merged_bams** (incomplete/aborted 4-aligner run). Its chunks are the only
  copy — deleting them loses data.
- `408d_plants` (127 G) — **no `aligner_chunks`/`merged_bams`** (different pipeline layout: `bam/ fastq/ out/
  raw/ trim/`). Not reclaimable via this lever; needs a separate inventory if space is wanted there.

## Residual caveats (my recommendation: run item 1 before deleting the 28 G Xrn1_repB)
1. Read-count equality proves no reads were **dropped**; it does not prove merged records are **field-identical**
   to the chunk records. Cheap spot-check on one sample: compare a merged BAM vs its concatenated chunks on the
   normalized field view before deleting.
2. "`junction_pool.pkl` makes `.junction_refined_*.bam` redundant" is asserted (the pkl exists) but not proven
   to be their sole consumer.

## Recommended deletion method (rm-safety compliant — NO `rm -rf`, NO globs)
Per approved directory, delete files then empty dirs, scoped to the one named path:
```
D="…/Dcp2_AID_25C_repA/aligner_chunks"      # one explicit named dir at a time
find "$D" -type f -delete && find "$D" -depth -type d -delete
```
Safer alternative (reversible until final purge): `mv "$D" "$D.TRASH"` first, confirm the run still resumes,
then delete `.TRASH`. Kevin chooses; I execute only the approved list.
