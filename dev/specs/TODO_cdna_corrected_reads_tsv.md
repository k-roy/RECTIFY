# TODO — cDNA pipeline should emit a `corrected_reads.tsv` (DRS-equivalent)

**Filed:** 2026-05-30 (Kevin's request). **Effort:** small — the data already exists.

## Problem
DRS (`rectify correct`) emits a per-read **`corrected_reads.tsv`** with columns incl.
`chrom, strand, original_3prime, corrected_3prime, gene_id, …`. Every downstream
3′-end analysis (metagenes, pileups, the cross-modality TRT/ART framework) consumes
that one TSV and never re-parses the BAM.

The **cDNA pipeline** (`rectify cdna-analyze`) emits no such file. Its 3′ ends are
only recoverable from `clusters.tsv` via the `anchor`/`orient` columns:
- `anchor` = gene-strand 3′ end (CPA site): `orient=fwd` → `anchor==aln_end`;
  `orient=rev` → `anchor==aln_start`.
- `orient` ∈ {fwd, rev} → strand {+, −}.

So consumers must hand-map `anchor`/`orient` → `corrected_3prime`/`strand` per analysis
(done ad hoc for the 4-WT batch comparison, 2026-05-30). Brittle and non-uniform.

## Ask
Have the cDNA pipeline emit a **`corrected_reads.tsv`** with the same schema as DRS
(at minimum `chrom, strand, corrected_3prime`; ideally also `original_3prime`,
`gene_id`, a read/molecule id, `tail_len`), **one row per UMI-consensus molecule**
(= one row per `clusters.tsv` cluster). Then DRS and cDNA 3′-end analyses share one loader.

## Where it wires in
- `rectify/core/commands/cdna_analyze_command.py` — alongside the existing
  `clusters.tsv` / `isoforms.tsv` writes.
- `rectify/core/cdna/cluster.py` — clusters already carry `chrom, anchor, orient,
  n_reads, tail_len, gene`; emitting the TSV is a projection:
  `corrected_3prime=anchor`, `strand={fwd:+, rev:-}[orient]`.

## Notes
- Unit: per-**molecule** (UMI-collapsed), not per-read — that's the cDNA equivalent of
  DRS per-read and the correct granularity for molecule-level 3′-end pileups.
- Keep column names identical to DRS so `load_index(usecols=['chrom','strand','corrected_3prime'])`
  works unchanged across modalities.
