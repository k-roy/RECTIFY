# Magic-BLAST — candidate aligner (read-indexed BLAST lineage)

**Status (2026-06-14): EXPLORATORY, smoke-test gated.** Algorithmically the most
unusual seeding in the survey, but **no published ONT direct-RNA track record** —
the smoke-test must first confirm it produces usable spliced alignments on noisy
ONT reads at all (Sherlock job 29546795).

## Why it's interesting (algorithmic orthogonality)

Magic-BLAST inverts the indexing direction of every other panel member:

- **Indexes the READS, not the genome** — builds a lookup table of word locations
  across the read batch, then streams genome sequence past it.
- Seeds are extended by **BLAST local gapped extension**, then stitched by
  **optimizing a spliced-alignment score** (its own splice model — distinct from
  minimap2's kernel, so it adds splice-model diversity).

It is the only **read-indexed** aligner considered; orthogonal to every
genome-indexed member (minimap2/uLTRA/deSALT/BBMap/GMAP/Graphmap2). Functionally
it is seed-and-extend, somewhat closer to BBMap's k-mer+DP family than GMAP is, so
claim "orthogonal lineage" but not "clean orthogonality."

## Splice / output

- Splice-aware: emits SAM with **N-op introns**, MAPQ, MD corrected around
  introns; also a BTOP/extended-BTOP tabular splice-signal format.
- Paper claims best-in-class intron discovery over wide conditions and best
  mapping for reads >250 bp from any platform; robust to high mismatch / extreme
  base composition.

## ONT suitability (UNKNOWN — the key risk)

- Aligns nanopore reads in principle, but the paper's benchmarks emphasize
  PacBio/Illumina/454. **No dedicated published ONT-DRS accuracy numbers.**
- This is exactly the BBMap lesson: "aligns long reads" ≠ "handles ONT direct-RNA
  error model." Treat ONT-DRS suitability as **unproven** until the smoke-test +
  per-aligner sole-win quality measurement on our data.

## Invocation

```bash
makeblastdb -in GRCh38_chr5.fa -dbtype nucl -out <db>          # build BLAST db
magicblast -query reads.fastq -db <db> -infmt fastq \
           -out out.sam -num_threads <n> -no_unaligned
```

## Install / license / maintenance

- bioconda: `conda install -c bioconda magicblast` (v1.7.2, 2023).
- NCBI public domain; NCBI-maintained.

## Sources

- Boratyn et al. 2019, *BMC Bioinformatics* ([PMC6659269](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6659269/)).
- Output format: [ncbi.github.io/magicblast](https://ncbi.github.io/magicblast/doc/output.html).
- 2026-06-14 orthogonality panel: memory `project-aligner-orthogonality-panel`.
