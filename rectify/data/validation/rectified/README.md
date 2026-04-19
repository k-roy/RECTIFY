# rectified

**Description:** RECTIFY corrected 3' end positions

**Created:** 2026-04-19T11:29:44.605584
**Last updated:** 2026-04-19T11:29:46.731729
**Total runs:** 1

## Files (2)

| File | Size | Modified |
|------|------|----------|
| corrected_3ends.tsv | 6.3 KB | 2026-04-19 |
| corrected_3ends_stats.tsv | 1.1 KB | 2026-04-19 |

## Recent Runs

### Run 1 (2026-04-19)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 4 -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam
```

Outputs: 2 files

---
*Full provenance in `PROVENANCE.json`*
