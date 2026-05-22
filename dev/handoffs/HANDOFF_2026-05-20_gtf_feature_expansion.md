# GTF/GFF Feature Expansion For Genomic Distribution

Date: 2026-05-20

Status: implemented in Codex working tree on 2026-05-20. Kept as a design note
for review and future agents.

## Problem

`rectify.core.analyze.loaders._parse_gtf()` currently keeps only `gene` rows.
That loses feature-level annotations needed by genomic-distribution code.

Current relevant files:

- `rectify/core/analyze/loaders.py`
- `rectify/core/analyze/genomic_distribution.py`
- `rectify/core/analyze/gene_attribution.py`
- `rectify/core/commands/analyze_command.py`
- `rectify/core/analyze/manifest.py`

Observed impact:

- Parsed GTF/GFF falls back to broad gene-body heuristics.
- Feature categories such as `UTR3`, `UTR5`, `CDS`, `CUT`, `SUT`, `XUT`,
  `snoRNA`, and other ncRNA/rRNA feature types can be absent even when present
  in the annotation file.
- Genomic distribution plots and summaries can under-report UTR/CDS/ncRNA/rDNA
  classes and over-report generic gene body/intergenic classes.

## Desired Behavior

`load_annotation()` should return an annotation DataFrame that preserves both:

- gene-level rows for gene attribution and nearest-gene distance, and
- feature-level rows for genomic-distribution classification.

Expected columns should remain compatible with existing callers:

- `chrom`
- `start` 0-based inclusive
- `end` 0-based exclusive
- `strand`
- `gene_id`
- `gene_name`
- `feature_type`

Additional optional columns are fine if useful:

- `transcript_id`
- `parent_id`
- `source`
- `raw_feature`

## Coordinate Rules

GTF/GFF coordinates are 1-based inclusive. Convert to 0-based half-open:

```text
start = int(fields[3]) - 1
end = int(fields[4])
```

Do not apply `end - 1` in annotation tables. Position queries should use
`start <= pos < end`.

## Feature Mapping Guidance

Start conservative. Preserve raw feature names, then map common names into the
categories expected by `genomic_distribution.py`.

Likely mappings:

- `gene` -> `gene`
- `CDS` -> `CDS`
- `three_prime_UTR`, `3UTR`, `UTR3`, `three_prime_utr` -> `UTR3`
- `five_prime_UTR`, `5UTR`, `UTR5`, `five_prime_utr` -> `UTR5`
- `exon` -> `exon`
- `intron` -> `intron`
- `snoRNA`, `snRNA`, `tRNA`, `rRNA`, `ncRNA` -> preserve or map to `ncRNA`
  depending on downstream expectations.
- Yeast-specific `CUT`, `SUT`, `XUT` should be preserved as feature types.

For GFF3 attributes:

- `ID`
- `Name`
- `gene`
- `Parent`
- `gene_id`
- `gene_name`

For GTF attributes:

- `gene_id`
- `gene_name`
- `transcript_id`

## Implementation Plan

1. Split `_parse_gtf()` into smaller helpers:
   - attribute parser,
   - feature-type normalizer,
   - row builder.
2. Parse all rows with usable coordinates and strand, not just `gene`.
3. Ensure every feature row has `gene_id` and `gene_name` when possible:
   - direct attributes first,
   - `Parent` fallback,
   - feature `ID` fallback,
   - coordinate-derived fallback only as last resort.
4. Keep gene rows, because `annotate_clusters_with_genes()` expects gene-like
   intervals.
5. Update genomic-distribution tests to verify feature-level classification.
6. Check whether callers that expect one row per gene need filtering to
   `feature_type == 'gene'` or broader accepted gene types.

## Test Fixtures To Add

Create a tiny synthetic GFF/GTF with:

- one plus-strand gene,
- one minus-strand gene,
- `CDS`,
- `UTR3`,
- `UTR5`,
- one ncRNA feature such as `snoRNA`,
- one yeast-specific transcript class such as `CUT`.

Assertions:

- Coordinates are converted to half-open exactly.
- Feature rows are present in `load_annotation()` output.
- Gene rows remain present.
- `calculate_genomic_distribution()` classifies positions inside UTR/CDS/ncRNA
  features as those features, not generic gene body.
- Existing gene attribution tests still pass.

## Validation Commands

Targeted:

```bash
pytest tests/test_analyze.py tests/test_metagene_loaders.py tests/test_visualize_gene_track.py -q
```

Broader analysis pass:

```bash
pytest tests/test_analyze.py tests/test_cdna_analyze.py tests/test_metagene.py -q
```

## Notes From 2026-05-20 Audit-Fix Tranche

Related coordinate fixes already landed:

- plus-strand gene 3' distance should use `gene['end'] - 1`,
- rDNA/exclusion filters should use half-open `pos < end`,
- bedgraph emits corrected positions as `[pos, pos + 1)`.

Keep this same convention when expanding feature parsing.
