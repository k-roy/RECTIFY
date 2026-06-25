# rectify correct

Correct 3' (and 5') end positions in a BAM or FASTQ file.

This is the core RECTIFY command. It applies the walk-back indel correction algorithm, A-tract ambiguity detection, AG-mispriming screening, and optionally NET-seq refinement.

---

## Usage

```bash
rectify correct <input> [options] -o <output>
```

## Examples

```bash
# Bundled yeast data — no external files needed
rectify correct reads.bam --Scer -o corrected.tsv

# Custom genome + annotation
rectify correct reads.bam \
    --genome genome.fa.gz \
    --annotation genes.gff.gz \
    -o corrected.tsv

# With NET-seq refinement
rectify correct reads.bam --Scer --netseq-dir my_netseq/ -o corrected.tsv

# Streaming mode for large BAMs (>2 GB)
rectify correct large.bam --Scer --streaming -o corrected.tsv

# Remove spike-in reads
rectify correct reads.bam --Scer --filter-spikein ENO2 -o corrected.tsv
```

---

## Arguments

### Required

| Argument | Description |
|----------|-------------|
| `input` | Input BAM or FASTQ/FASTQ.GZ file |
| `-o, --output` | Output TSV file path |

### Reference data

| Argument | Description |
|----------|-------------|
| `--genome` | Reference genome FASTA (optionally gzipped) |
| `--annotation` | Gene annotation file (GTF or GFF, optionally gzipped) |
| `--Scer` | Use bundled *S. cerevisiae* S288C data (genome + annotation + NET-seq) |
| `--organism` | Organism name for bundled data (e.g. `yeast`) |

### NET-seq refinement

| Argument | Default | Description |
|----------|---------|-------------|
| `--netseq-dir` | — | Directory of NET-seq BigWig files for A-tract refinement |
| `--netseq-samples` | all | Specific NET-seq samples to use |

### Technology / protocol

| Argument | Default | Description |
|----------|---------|-------------|
| `--dT-primed-cDNA` | off | Input was generated with oligo-dT priming (QuantSeq, dT-primed cDNA). Poly(A) tail is NOT in the read. Enables indel artifact correction and poly(A) trimming modules. Do NOT use for ONT direct RNA. |
| `--ONT-cDNA` | off | Input is Oxford Nanopore PCR-cDNA (e.g. SQK-PCB114). Poly(A) tail IS present as a 3' soft-clip. Enables poly(A) trimming + indel correction; disables AG-mispriming. Do NOT combine with `--dT-primed-cDNA`. |
| `--short-read` | off | Input is short-read data (Illumina/Aviti ≤150 bp). Disables poly(A)-tail trimming, A-tract correction, and indel modules. Pair with `rectify align --short-read` (bbmap + bwa). |

*(deprecated alias retained for backwards compatibility: `--polya-sequenced` (maps to `--dT-primed-cDNA`). Use the current flags above for new pipelines. Note: `--no-polya-sequenced` is a deprecated alias for `rectify run-all`, not for `rectify correct`.)*

### Poly(A) handling

| Argument | Default | Description |
|----------|---------|-------------|
| `--polya-model` | built-in | Pre-trained poly(A) tail model (JSON from `rectify train-polya`) |
| `--min-polya-score` | — | Minimum poly(A) model confidence score (0–1) to mark `polya_pass=1`; reads below are flagged but still written. Requires `--polya-model` or `--dT-primed-cDNA`. |
| `--use-dorado-polya` | off | Treat Dorado's `pt:i` tail-length estimate as the authoritative `polya_length` (otherwise measured from alignment/soft-clip; the `pt:i` value is always recorded in the `pt_tag` column). |

### Module selection

| Argument | Description |
|----------|-------------|
| `--skip-atract-check` | Skip A-tract ambiguity detection |
| `--skip-ag-check` | Skip AG-mispriming screening |
| `--skip-polya-trim` | Skip poly(A) tail trimming |
| `--skip-indel-correction` | Skip indel artifact correction |
| `--skip-3ss-rescue` | Skip Cat3 5' splice-site truncation rescue |
| `--skip-variant-aware` | Skip variant-aware homopolymer rescue (two-pass) |
| `--min-mapq` | Minimum MAPQ to include a read (default 0; use 1 to drop multi-mappers) |
| `--min-aligned-length` | Minimum reference bases consumed by the alignment (default 0; 30 recommended for dT-primed short-read) |

### Filtering

| Argument | Description |
|----------|-------------|
| `--filter-spikein GENE [GENE ...]` | Remove reads from spike-in genes by name |

### Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--ag-threshold` | 17.0 | Weighted AG-richness composite score threshold for mispriming flag. Range 0.0–34.5; default 17.0 is calibrated for yeast DRS data. This is NOT a simple fraction — do not use values in the 0–1 range. |

### Output options

| Argument | Description |
|----------|-------------|
| `--output-bam` | Write a poly(A)-trimmed BAM alongside the corrected TSV (requires `--dT-primed-cDNA`) |
| `--write-corrected-bam` | Write a corrected BAM hard-clipped at every read's corrected 3' end (sorted + indexed) |
| `--write-softclipped-bam` | Like `--write-corrected-bam` but poly(A) bases are soft-clipped (visible in IGV) |
| `--write-bedgraph PREFIX` | Write strand-specific bedGraph files (`PREFIX.plus.bedgraph` / `PREFIX.minus.bedgraph`) for NET-seq fractional pileups |
| `--report` | Write QC report to this path (`.html` or `.pdf`) |

### Junction refinement (Module 2H)

Provide per-aligner BAMs to activate post-consensus N-op refinement: every intron N-op in
the winning BAM is re-scored against all candidate junctions observed across all aligner BAMs
within a search radius, and replaced with the best sequence-supported junction.

| Argument | Default | Description |
|----------|---------|-------------|
| `--aligner-bams BAM` | (none) | Per-aligner BAM for junction candidate pool. Repeat for each aligner: `--aligner-bams minimap2.bam --aligner-bams gapmm2.bam`. Accepts plain paths or `aligner:path` pairs. When supplied, enables Module 2H. |
| `--junction-hp-pen FLOAT` | 0.25 | Homopolymer indel penalty for junction scoring (0 < value ≤ 1.0). Lower values tolerate poly-T/A undercalling errors near splice sites. 0.25 is recommended for Nanopore DRS. |
| `--junction-search-radius BP` | 5000 | Search radius (bp) around each N-op for candidate junctions. Covers the longest known yeast intron (~1 kb) with margin. |
| `--junction-window BP` | 15 | Edit-distance window half-width for split-alignment scoring (bp on each side of the junction). |
| `--junction-max-slide BP` | 10 | Maximum query split displacement for split-alignment scoring (±bp). |
| `--junction-max-boundary-shift BP` | 50 | Maximum allowed shift of either intron boundary when applying a junction replacement. Prevents false-positive matches from junctions in neighbouring genes. |
| `--junction-penalty-table PATH` | (heuristic) | Path to empirical HP-context penalty table (`penalty_scores.tsv`) produced by `empirical_cigar_error_profiler.py`. Overrides heuristic del/ins costs with per-HP-length values derived from multi-aligner agreement on this dataset. |
| `--str-penalty-table PATH` | (none) | Path to STR penalty table (`str_penalty_scores.tsv`). Scores dinucleotide/trinucleotide repeat contexts with empirical STR-specific penalties alongside `--junction-penalty-table`. |

See [Multi-Aligner Consensus](../../algorithms/multi_aligner_consensus.md) and
[docs/EMPIRICAL_HP_PENALTY_SCORING.md](../../EMPIRICAL_HP_PENALTY_SCORING.md) for scoring details.

### Checkpoint / resume / sidecar

| Argument | Default | Description |
|----------|---------|-------------|
| `--checkpoint-dir DIR` | (none) | Directory for per-region checkpoint sentinels and the variant-scan pickle. Enables resume: completed regions are skipped on re-run, the variant scan pickle is reloaded, and the partial output TSV is appended to. Has no effect without `--streaming`. On H2/Sherlock, use `$L_SCRATCH/rectify_ckpt`. |
| `--tmp-dir DIR` | `$TMPDIR` | Scratch directory for intermediate per-region BAM files used by the parallel BAM writer. On H2/Sherlock, pass `$L_SCRATCH/rectify_regions` for fast local disk. |
| `--variant-scan-cache PKL` | (none) | Pre-computed variant scan pickle from `rectify prescan`. Skips the Pass-1 all-reads scan (Module 2D). Use when running per-chunk correction with a shared scan built from the merged BAM. |
| `--junction-pool-cache PKL` | (none) | Pre-computed junction pool pickle from `rectify prescan`. Skips the aligner-BAM junction collection step (Module 2H). |
| `--emit-merged-tsv` | off | Also emit a legacy concatenated `corrected_reads.tsv` alongside the `corrected_reads.manifest.tsv`. For downstream scripts not yet migrated to accept manifests. |
| `--force-all` | off | Ignore all sidecars; rerun every stage unconditionally. |
| `--force-stage NAME[,NAME...]` | (none) | Force-rerun specific stages (e.g. `analyze,correct`). Stages downstream of a forced stage are also forced. |
| `--accept-prior-provenance` | off | Treat git SHA mismatch between prior and current run as non-blocking. Use only for cosmetic commits (docs, comments) that don't affect output bytes. |
| `--dry-run-resume` | off | Print the SKIP/RUN decision for each stage (with diff of which input or argv changed) and exit. No work done. |
| `--legacy-single-threaded` | off | Force the pre-parallel single-threaded BAM writer. For debugging correctness regressions only. **DEPRECATED** — will be removed in a future release. |

### Performance

| Argument | Default | Description |
|----------|---------|-------------|
| `-j, --threads` | auto | Number of processing threads (0 = auto-detect) |
| `--streaming` | off | Streaming output mode — keeps peak RAM ~4–5 GB for any BAM size |
| `--chunk-size` | 10000 | Reads per output chunk (streaming mode only) |

---

## Output files

| File | Description |
|------|-------------|
| `<output>.tsv` | Per-read corrected positions (see [Output Formats](../output_formats.md)) |
| `<output>.manifest.tsv` | Region manifest pointing at the per-region TSV(s); written by default in the streaming/parallel path. Pass `--emit-merged-tsv` to also write the concatenated `<output>.tsv`, or run `rectify export-merged-tsv <manifest>` later |
| `<output>_stats.tsv` | Processing QC summary |
| `<output>_potential_variants.tsv` | Candidate genomic variants flagged during variant-aware homopolymer rescue (when enabled) |

---

## Notes

- For FASTQ input, pre-align first with `rectify align`, then pass the BAM to `rectify correct`
- To run alignment + correction in one step, use `rectify run-all`
- `rectify correct` runs post-consensus on the winning aligner's BAM; use `--aligner-bams` to supply per-aligner BAMs as a junction candidate pool for Module 2H (post-consensus N-op refinement)
- `--streaming` is recommended for BAMs larger than 2 GB; it is the default in the bundled SLURM profiles
- In the streaming/parallel path, `correct` writes a `corrected_reads.manifest.tsv` by default and defers the concatenated TSV; pass `--emit-merged-tsv` or call `rectify export-merged-tsv` when a single flat TSV is needed
