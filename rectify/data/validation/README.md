# validation

**Description:** RECTIFY corrected 3' end positions

**Created:** 2026-04-11T17:36:00.417595
**Last updated:** 2026-04-22T17:02:38.415788
**Total runs:** 29

## Files (2)

| File | Size | Modified |
|------|------|----------|
| corrected_3ends.tsv | 7.1 KB | 2026-04-22 |
| corrected_3ends_stats.tsv | 1.1 KB | 2026-04-22 |

## Recent Runs

### Run 25 (2026-04-22)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 8 -o rectify/data/validation/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/rectified_corrected_3end --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

### Run 26 (2026-04-22)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 8 -o rectify/data/validation/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/rectified_corrected_3end --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

### Run 27 (2026-04-22)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/validation_reads.bam --genome /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 8 -o /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/corrected_3ends.tsv --write-corrected-bam /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/rectified/rectified_corrected_3end --aligner-bams minimap2:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

### Run 28 (2026-04-22)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 8 -o rectify/data/validation/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/rectified_corrected_3end --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

### Run 29 (2026-04-22)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 8 -o rectify/data/validation/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/rectified_corrected_3end --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

---
*Full provenance in `PROVENANCE.json`*
