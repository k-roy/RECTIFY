# rectified

**Description:** RECTIFY corrected 3' end positions

**Created:** 2026-04-19T11:29:44.605584
**Last updated:** 2026-04-21T16:47:01.518775
**Total runs:** 4

## Files (2)

| File | Size | Modified |
|------|------|----------|
| corrected_3ends.tsv | 7.1 KB | 2026-04-21 |
| corrected_3ends_stats.tsv | 1.1 KB | 2026-04-21 |

## Recent Runs

### Run 1 (2026-04-19)

```
rectify correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 4 -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam
```

Outputs: 2 files

### Run 2 (2026-04-20)

```
rectify correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 4 -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam
```

Outputs: 2 files

### Run 3 (2026-04-21)

```
rectify correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 4 -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/rectified_corrected_3end
```

Outputs: 2 files

### Run 4 (2026-04-21)

```
rectify correct rectify/data/validation/validation_reads.bam --genome rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz --annotation rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz --threads 4 --aligner-bams minimap2:rectify/data/validation/aligners/validation_reads.minimap2.bam --aligner-bams gapmm2:rectify/data/validation/aligners/validation_reads.gapmm2.bam --aligner-bams mapPacBio:rectify/data/validation/aligners/validation_reads.mapPacBio.bam --aligner-bams deSALT:rectify/data/validation/aligners/validation_reads.deSALT.bam --aligner-bams uLTRA:rectify/data/validation/aligners/validation_reads.uLTRA.bam -o rectify/data/validation/rectified/corrected_3ends.tsv --write-corrected-bam rectify/data/validation/rectified/rectified_corrected_3end.bam --write-softclipped-bam rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam --write-bedgraph rectify/data/validation/rectified/rectified_corrected_3end
```

Outputs: 2 files

---
*Full provenance in `PROVENANCE.json`*
