# validation

**Description:** RECTIFY corrected 3' end positions

**Created:** 2026-04-11T17:36:00.417595
**Last updated:** 2026-04-21T12:42:07.209433
**Total runs:** 20

## Files (2)

| File | Size | Modified |
|------|------|----------|
| corrected_3ends.tsv | 7.1 KB | 2026-04-21 |
| corrected_3ends_stats.tsv | 1.1 KB | 2026-04-21 |

## Recent Runs

### Run 16 (2026-04-14)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/output/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff --write-corrected-bam rectify/data/output/rectified/rectified_pA_hardclip.bam --write-softclipped-bam rectify/data/output/rectified/rectified_pA_softclip.bam --write-bedgraph rectify/data/output/rectified/rectified_pA_hardclip -o rectify/data/corrected_3ends.tsv
```

Outputs: 2 files

### Run 17 (2026-04-14)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/output/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --netseq-dir rectify/data -o rectify/data/corrected_3ends.tsv --report rectify/data/corrected_3ends_stats.tsv --threads 4
```

Outputs: 2 files

### Run 18 (2026-04-16)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/output/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff --Scer -o rectify/data/corrected_3ends.tsv --write-corrected-bam rectify/data/output/aligners/validation_reads.rectified.bam --write-softclipped-bam rectify/data/output/aligners/validation_reads.rectified_softclip.bam --write-bedgraph rectify/data/output/aligners/validation_reads.rectified -j 8
```

Outputs: 2 files

### Run 19 (2026-04-16)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/output/validation_reads.bam --Scer --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz -o rectify/data/corrected_3ends.tsv --write-corrected-bam rectify/data/output/aligners/validation_reads.rectified.bam --write-softclipped-bam rectify/data/output/aligners/validation_reads.rectified_softclip.bam --write-bedgraph rectify/data/output/aligners/validation_reads.rectified
```

Outputs: 2 files

### Run 20 (2026-04-21)

```
rectify correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 4 -o rectify/data/validation/corrected_3ends.tsv
```

Outputs: 2 files

---
*Full provenance in `PROVENANCE.json`*
