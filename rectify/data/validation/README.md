# validation

**Description:** RECTIFY corrected 3' end positions

**Created:** 2026-04-11T17:36:00.417595
**Last updated:** 2026-04-27T18:22:21.049312
**Total runs:** 37

## Files (2)

| File | Size | Modified |
|------|------|----------|
| corrected_3ends.tsv | 6.9 KB | 2026-04-27 |
| corrected_3ends_stats.tsv | 1.1 KB | 2026-04-27 |

## Recent Runs

### Run 33 (2026-04-24)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff --threads 8 -o rectify/data/validation/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/rectified_corrected_3end --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

### Run 34 (2026-04-24)

```
replace_cat1_reads.py + replace_cat2_reads.py + rebuild_from_committed.py (Cat1/Cat2 DRS replacement) + rectify correct --aligner-bams (regenerate corrected_3ends.tsv)
```

Outputs: 2 files

### Run 35 (2026-04-27)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff --threads 8 -o rectify/data/validation/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

### Run 36 (2026-04-27)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff --threads 8 -o rectify/data/validation/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/corrected_3ends --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

### Run 37 (2026-04-27)

```
/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff --threads 8 -o rectify/data/validation/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/corrected_3ends --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

---
*Full provenance in `PROVENANCE.json`*
