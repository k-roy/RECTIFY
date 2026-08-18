# rectified

**Description:** RECTIFY corrected 3' end positions

**Created:** 2026-04-19T11:29:44.605584
**Last updated:** 2026-04-27T18:07:30.932546
**Total runs:** 16

## Files (2)

| File | Size | Modified |
|------|------|----------|
| corrected_3ends.tsv | 6.9 KB | 2026-04-27 |
| corrected_3ends_stats.tsv | 1.1 KB | 2026-04-27 |

## Recent Runs

### Run 12 (2026-04-24)

```
rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --threads 8
```

Outputs: 2 files

### Run 13 (2026-04-24)

```
rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff --threads 8 -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam
```

Outputs: 2 files

### Run 14 (2026-04-24)

```
rectify/cli.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 4 --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam
```

Outputs: 2 files

### Run 15 (2026-04-27)

```
rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff --threads 8 -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/corrected_3ends --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

### Run 16 (2026-04-27)

```
rectify/__main__.py correct rectify/data/validation/validation_reads.bam --genome rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa --annotation rectify/data/saccharomyces_cerevisiae_R64-5-1_20240529.gff --threads 8 -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/corrected_3ends --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam
```

Outputs: 2 files

---
*Full provenance in `PROVENANCE.json`*
