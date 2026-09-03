# rectify netseq

Optional refinement step. Process NET-seq BAM files: extract 3' end positions, apply oligo(A)-spreading deconvolution, and write strand-specific BigWig files for use with `--netseq-dir` in `rectify correct`.

---

## Usage

```bash
rectify netseq <input.bam> [options] -o <output_dir>
```

## Examples

```bash
# Basic NET-seq processing
rectify netseq netseq.bam \
    --genome genome.fa.gz \
    --gff genes.gff.gz \
    -o netseq_output/

# With exclusion regions
rectify netseq netseq.bam \
    --genome genome.fa.gz \
    --gff genes.gff.gz \
    --exclude-mito \
    -o netseq_output/

# Include rDNA and Pol III genes
rectify netseq netseq.bam \
    --Scer \
    --include-rdna \
    --include-pol3 \
    -o netseq_output/
```

---

## Arguments

### Required

| Argument | Description |
|----------|-------------|
| `input` | Input NET-seq BAM file |
| `-o, --output-dir` | Output directory |

### Reference data

| Argument | Description |
|----------|-------------|
| `--genome` | Reference genome FASTA |
| `--gff` | Gene annotation GFF file |
| `--Scer` | Bundled *S. cerevisiae* data |

### Exclusion regions

| Argument | Default | Description |
|----------|---------|-------------|
| `--exclude-mito` | off | Exclude mitochondrial reads |
| `--include-rdna` | off | Include rDNA loci (excluded by default) |
| `--include-pol3` | off | Include Pol III transcripts (excluded by default) |
| `--pol3-flanking` | 100 | Flanking bp around Pol III genes to exclude |

### Donor-side junction rescue and 3'-tail calling

| Argument | Default | Description |
|----------|---------|-------------|
| `--junction-rescue` / `--no-junction-rescue` | ON when a junction source exists | Re-place a 1-10 nt exon-2 overhang across the 5' splice site |
| `--junction-pool FILE` | — | External junction TSV (`chrom`/`donor`/`acceptor`/`strand`) merged with the annotated introns |
| `--pool-include-trna` | off | Keep tRNA introns in the pool (they are dropped by default — excised by the tRNA endonuclease, not the spliceosome, at Pol III loci) |
| `--rescue-max-intronic` | 10 | How far past the donor the aligned RNA 3' end may sit |
| `--rescue-min-k` | 1 | Minimum recovered exon-2 length when there is NO non-templated remainder (`r == 0`) |
| `--rescue-min-k-with-remainder` | 4 | The floor when a randomer remainder is invoked (`r > 0`) — that is where the chance matches live |
| `--walkback-requires-clip-a` | **ON** | Only walk back when the clip carries a non-templated A (nascent 3' ends have no tail by default) |
| `--walkback-unconditional` | off | Restore invariant-7 walkback, for a poly(A)-selected input |
| `--no-tail-detection` | off | Disable the tail call entirely |

See [NET-seq Refinement](../../algorithms/netseq_refinement.md#donor-side-junction-rescue-rectify-netseq)
for the geometry, the acceptance rule and the chance-match null.

### Deconvolution

| Argument | Default | Description |
|----------|---------|-------------|
| `--no-deconvolution` | off | Disable NNLS deconvolution (output raw 3' ends only) |
| `--min-atract-length` | 3 | Minimum downstream A's to trigger deconvolution |

### Output options

| Argument | Default | Description |
|----------|---------|-------------|
| `--output-format` | parquet bedgraph | One or more of `parquet`, `bedgraph`, `bigwig`, `tsv`. `parquet` falls back to `tsv` with a warning when pyarrow is not installed |
| `--track-position` | corrected | Which 3' end drives the primary track (deconvolution input): `corrected` or `raw` |
| `--no-rpm-normalize` | off | Disable RPM normalization for bedgraph/bigwig output |

### Processing

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-mapq` | 0 | Minimum mapping quality |
| `--max-reads` | — | Cap total reads processed (for testing) |
| `-v, --verbose` | off | Verbose output |

---

## Output files

Defaults: per-read records as parquet plus per-position signal as bedgraph (one strand-specific file each, raw and deconvolved). Use `--output-format bigwig` to produce BigWig files in addition (requires `pyBigWig`); `--output-format tsv` writes a plain-text equivalent of the parquet.

| File | Description |
|------|-------------|
| `{sample}.unified_reads.parquet` | Per-read 3' end records (default) |
| `{sample}.raw.plus.bedgraph` | Plus strand 3' end signal (default; RPM-normalized unless `--no-rpm-normalize`) |
| `{sample}.raw.minus.bedgraph` | Minus strand 3' end signal (default) |
| `{sample}.deconv.plus.bedgraph` | NNLS-deconvolved plus strand signal (default) |
| `{sample}.deconv.minus.bedgraph` | NNLS-deconvolved minus strand signal (default) |
| `{sample}.corrected.plus.bedgraph` | Plus strand signal on the CORRECTED 3' end (tail walkback + junction rescue) |
| `{sample}.corrected.minus.bedgraph` | Minus strand corrected signal |
| `{sample}.netseq_summary.json` | Per-sample correction counters: rescue classes, rescued-by-k (split into the clean and randomer channels), the decoy-acceptor null, tail-length histogram, and the fraction of ends that moved |
| `{sample}.{raw,corrected,deconv}.{plus,minus}.bw` | BigWig variants of the above (only with `--output-format bigwig`) |

`--track-position raw|corrected` (default `corrected`) selects which track the NNLS deconvolution
consumes; **both** the `raw` and `corrected` bedgraphs are always written, so before/after is
diffable without a second pass.

---

## Oligo(A) deconvolution

NET-seq 3' ends are oligo-adenylated, creating downstream signal spreading. RECTIFY applies NNLS deconvolution using a point-spread function (PSF) derived from 5,000+ zero-A control sites to recover true CPA positions.

The deconvolved BigWigs are used by `rectify correct --netseq-dir` for A-tract ambiguity resolution in direct RNA experiments.

See [NET-seq Refinement](../../algorithms/netseq_refinement.md) for details on the deconvolution algorithm.

---

## Notes

- rDNA and Pol III genes are excluded by default because their extremely high signal density distorts the deconvolution PSF
- The PSF is derived from positions with zero downstream A's — sites where the true CPA is unambiguous
- Use `--include-rdna` and `--include-pol3` only if you specifically want to analyze those loci
