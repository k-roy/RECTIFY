# rectified

**Description:** RECTIFY corrected 3' end positions

**Created:** 2026-04-19T11:29:44.605584
**Last updated:** 2026-04-22T16:45:57.376172
**Total runs:** 10

## Files (2)

| File | Size | Modified |
|------|------|----------|
| corrected_3ends.tsv | 7.1 KB | 2026-04-22 |
| corrected_3ends_stats.tsv | 1.1 KB | 2026-04-22 |

## Recent Runs

### Run 6 (2026-04-21)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/cli.py correct /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/validation_reads.bam --genome /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 4 --aligner-bams minimap2:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams deSALT:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/aligners/validation_reads.deSALT.bam --aligner-bams uLTRA:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/aligners/validation_reads.uLTRA.bam -o /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/rectified/rectified_corrected_3end
```

Outputs: 2 files

### Run 7 (2026-04-22)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 4 -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam
```

Outputs: 2 files

### Run 8 (2026-04-22)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff -o rectify/data/validation/rectified/corrected_3ends.tsv --threads 8
```

Outputs: 2 files

### Run 9 (2026-04-22)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --threads 8
```

Outputs: 2 files

### Run 10 (2026-04-22)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --threads 8
```

Outputs: 2 files

---
*Full provenance in `PROVENANCE.json`*
