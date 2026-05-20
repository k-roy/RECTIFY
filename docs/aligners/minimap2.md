# Minimap2 alignment: junction annotation

The `rectify align` command generates a junction BED from the GFF
annotation and passes it to minimap2 via `--junc-bed`. This improves
splice junction accuracy for long-read RNA-seq.

The junction BED is cached as `annotation.junc.bed` in the sample output
directory. Exact minimap2 command:

```
minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD
    -t <threads>
    --junc-bed <sample_dir>/annotation.junc.bed
    --junc-bonus 9
    <genome.fsa.gz> <reads.fastq.gz>
```

Key flags:
- `-uf`: forward-strand-only (correct for direct RNA / cDNA sense reads).
- `-k14`: smaller k-mer for sensitivity on noisy nanopore reads.
- `-G 5000`: max intron size (yeast introns < 1 kb).
- `--splice-flank=no`: disables GT-AG bonus (important for 3' end
  accuracy).
- `--MD`: required for indel artifact correction downstream.
