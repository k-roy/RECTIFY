# uLTRA alignment

uLTRA uses collinear chaining over a genome-graph built from the GTF annotation.
It excels at reads that span small exons (11–20 nt) that seed-chain aligners
(minimap2, mapPacBio, gapmm2) miss because no seed lands in the exon.

---

## GTF annotation requirements

uLTRA's index builder (`prep_splicing`) requires a GTF where:

- Column 3 contains `gene`, `transcript`, and `exon` features.
- **Every `exon` line carries a `gene_id` attribute.**  A single exon line
  without `gene_id` crashes uLTRA at `create_augmented_gene.py:323` with
  `KeyError: 'gene_id'`, which RECTIFY surfaces as:

  ```
  RuntimeError: uLTRA failed: ...
    File "create_augmented_gene.py", line 323, in create_graph_from_exon_parts
      exon_gene_ids = exon.attributes["gene_id"]
  KeyError: 'gene_id'
  ```

### Why gffread-produced GTFs fail

`gffread -T` produces GTFs with `exon` features, but for tRNAs, ncRNAs, and
some multi-exon mRNAs it omits `gene_id`, retaining only `transcript_id` and
`gene_name`.  A yeast R64 GTF generated this way has ~11 % of exon lines
missing `gene_id` (578 out of 5090 in the test case).

### How RECTIFY handles this automatically

When `annotation_path` is a GFF (`.gff`, `.gff3`, `.gff.gz`, `.gff3.gz`),
`run_ultra()` locates the sibling `.gtf` in the same directory and then
**always runs `_normalize_gtf_for_ultra()`** regardless of whether the GTF
already contains `exon` features.  The normalizer:

1. Collects existing `exon` intervals from the GTF (preserves real coordinates).
2. For transcripts with no exon lines, derives them from transcript span minus
   any `intron` intervals (SGD-style GTF path).
3. Ensures every emitted `exon` line carries `gene_id`, copied from the parent
   `transcript` line if absent.
4. Writes a cached `<stem>.ultra_norm.gtf` alongside the source; rebuilds
   only when the source is newer.

You never need to pre-process the GTF manually — RECTIFY does it.

---

## Sibling GTF requirement

uLTRA does not accept GFF input.  When the RECTIFY annotation is a GFF file,
`run_ultra()` looks for `<annotation_stem>.gtf` in the same directory.

```
saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz   ← annotation_path
saccharomyces_cerevisiae_R64-5-1_20240529.gtf       ← required sibling
```

If the sibling is missing, RECTIFY raises `FileNotFoundError` with the exact
`gffread` command to create it:

```bash
gffread saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz \
    -T -o saccharomyces_cerevisiae_R64-5-1_20240529.gtf
```

The bundled yeast genome in `rectify/data/genomes/saccharomyces_cerevisiae/`
does **not** ship with a sibling GTF.  The `scripts/setup/setup_rectify_env.sh`
bootstrap script does not generate one either — it must be created once per
checkout or shared env.  See [Environment setup](../installation.md) for the
full bootstrap checklist.

---

## namfinder dependency

uLTRA calls `namfinder` as a subprocess.  It is a separate binary not installed
as a transitive dependency of `ultra-bioinformatics`.  RECTIFY's `run_ultra()`
checks for `namfinder` on `PATH` and falls back to the vendored binary in
`rectify/data/bin/<platform>/namfinder` if needed.

The vendored binary is built against a newer GLIBCXX than some HPC nodes
(Hoffman2's older RHEL).  If you see:

```
namfinder: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.30' not found
```

install the conda-built binary instead:

```bash
conda install -c bioconda namfinder
```

This is performed automatically by `scripts/setup/setup_rectify_env.sh`.

---

## SGD yeast GTF assertion patch

uLTRA 0.1's `modules/create_augmented_gene.py` raises an `AssertionError`
on SGD yeast GTF (overlapping per-gene exons trigger it near 89–99 % of
iteration).  `setup_rectify_env.sh` comments the assertion out in-place
with an idempotent `awk` patch.  The algorithm completes correctly past it
for yeast use cases.

Without the patch, `rectify align` will silently drop all uLTRA results
(the exception is caught by `_run_one_aligner()` and logged, but the
pipeline continues with the remaining aligners).

---

## Verifying uLTRA works

```bash
# Quick smoke: align 50 reads, confirm BAM has mapped records
rectify align \
    --reads <reads.fastq.gz> \
    --genome <genome.fsa.gz> \
    --annotation <annotation.gff.gz> \
    --aligners uLTRA \
    --output /tmp/ultra_smoke \
    --threads 4

samtools flagstat /tmp/ultra_smoke/*.uLTRA.bam
```

Expected: nonzero mapped count.  If mapped = 0 or the BAM is absent, check
`rectify align` logs for `uLTRA failed:` lines.

---

## Failure modes quick-reference

| Symptom | Cause | Fix |
|---------|-------|-----|
| `database.db` created; no `.uLTRA.bam` emitted | `KeyError: 'gene_id'` in `create_augmented_gene.py` | RECTIFY ≥ post-2026-05-19 patch auto-fixes; for older installs: check sibling GTF for exons missing `gene_id` |
| `namfinder not found` | namfinder not on PATH | `conda install -c bioconda namfinder` |
| `GLIBCXX_3.4.30 not found` | Vendored namfinder built against newer GLIBCXX | `conda install -c bioconda namfinder` |
| `AssertionError` near 89–99 % of annotation | SGD overlapping exon triggers uLTRA internal assert | Run `scripts/setup/setup_rectify_env.sh` to patch |
| `FileNotFoundError: uLTRA requires a sibling GTF` | No `.gtf` next to `.gff.gz` | `gffread <ann.gff.gz> -T -o <ann.gtf>` |

---

## Primary-alignment & duplicate handling

`rectify align` passes **no secondary-suppression flag** to uLTRA — its command is
`pipeline --ont --disable_infer -t <n> --prefix ...` (`multi_aligner.py:2066-2076`),
with no uLTRA secondary control and no post-hoc dedup. Whether uLTRA emits secondary
records by default is not determinable from rectify source. This is safe for 3′-end
counting because `rectify correct` skips `is_secondary`/`is_supplementary` records
regardless of which aligner produced them.

What `rectify correct` does **not** guard against is duplicate **primary** records
in an external BAM (e.g. a doubled input FASTQ) → 2× double-counted 3′ ends. See the
canonical writeup and cross-aligner table in
[minimap2.md](minimap2.md#-duplicate-primary-alignments--2-double-counted-3-ends-external-bam-hazard).
