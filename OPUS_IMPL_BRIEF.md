# Implementation Brief — RECTIFY cDNA Pipeline Refactor
**For:** Opus 4.7 implementing agent  
**Repo:** `/Users/kevinroy/work/rectify`  
**Spec:** `PROPOSED_README_UPDATE.md` in the same repo (read this first — it is the authoritative reference for all biology, naming, and architectural decisions)  
**Date:** 2026-05-14

---

## What RECTIFY does (brief)

RECTIFY is a read-correction pipeline for 3'-end sequencing data. It corrects
poly-A tail mis-positioning via a read-vs-reference walkback, then calls
cleavage-and-polyadenylation (CPA) sites. It supports three protocol tracks:
DRS (ONT direct RNA-seq), ONT cDNA (PCB114.24 SSP+UMI), and QuantSeq REV
(short-read antisense). This brief concerns the cDNA track only.

---

## Current cDNA pipeline (what exists today)

`rectify correct-cdna` is implemented in `rectify/core/cdna_correct_command.py`
(1962 lines). Per chromosome it does:

1. Parse raw input BAM → classify reads as Type-1 (SSP+UMI captured) or Type-2
   (SSP-truncated, no UMI) via `parse_read_info()`
2. Bucket reads by 3' anchor position → UMI-directional cluster within each
   bucket → one cluster per molecule
3. For each cluster: build abPOA consensus (strand-split-then-merge if both
   is_reverse=True and is_reverse=False reads are present)
4. Re-align consensus with **minimap2 only** via `realign_consensus()` (uses
   mappy Python bindings)
5. Assign gene (XG tag), sense/antisense (XS tag)
6. Group Stage-1 clusters into isoforms: Type-1 uses 5'+3' position clustering;
   Type-2 uses 3' position only
7. Link Type-1 ↔ Type-2 clusters from the same molecule ("same-molecule strand
   pairing") via `reconcile_t1_t2_pairs()`
8. Write output BAM + `clusters.tsv` + `isoforms.tsv` + `t1t2_pairs.tsv`

The consensus sequence passed to `realign_consensus()` is **untrimmed** — it
still contains SSP/UMI/GGG at the 5' end and poly-A/adapter at the 3' end.
The aligner soft-clips these; post-hoc correction (`walk_forward_tss`,
`walk_back_anchor_and_tail`) then adjusts the terminal positions.

---

## Target architecture (what to build)

The goal is to align the cDNA downstream workflow with DRS so both tracks share
the same multi-aligner → chimeric consensus → walkback pipeline after their
respective upstream steps.

```
[cDNA-specific upstream]
  Raw input BAM (pre-aligned by user with minimap2)
    ↓
  parse_read_info() → Type-1a / Type-1b / Type-2 classification
    ↓
  UMI directional clustering
    ↓
  abPOA consensus (strand-split-then-merge)
    ↓
[Shared downstream — matches DRS from here]
  Pre-trim consensus:
    A. Strip SSP/UMI/GGG from 5' (or SSP_RC/UMI_RC/CCC from 5' for orient=rev)
    B. Strip poly-A + adapter from 3' (record trim length → XA tag)
    ↓
  Multi-aligner (all three on trimmed consensus):
    minimap2  (already present)
    mapPacBio (add)
    gapmm2    (add)
    ↓
  chimeric_consensus.select_best(minimap2, mpb, gapmm2)
    ↓
  Restore XA tag from pre-trim length
  walkback (3' CPA correction) + walk_forward_tss (5' TSS correction)
  Assign XG, XS, XY tags
    ↓
  Write corrected BAM (per-read tags: XU XC XR XF XA XG XS XY XL)

[Downstream analysis — separate step, shared with DRS]
  Read corrected BAM
    ↓
  5'/3' isoform clustering (Type-1a/1b: TSS+CPA; Type-2: CPA only)
  Same-molecule strand pairing (T1↔T2)
    ↓
  Write clusters.tsv, isoforms.tsv, t1t2_pairs.tsv
```

---

## Implementation tasks (in priority order)

### Task 1 — Add `read_subtype` to `ReadInfo` and BAM output

**Why:** The README uses Type-1a/1b/2 as primary terminology. The code must
match so downstream agents and users can read the code and docs together.

**What to change:**

In `ReadInfo` dataclass (search for `@dataclass` near line 258):
```python
read_subtype: str   # "1a", "1b", or "2"
```

In `parse_read_info()` (line ~558), set it after the orient/read_type
determination:
```python
if read_type == 1:
    read_subtype = "1a" if orient == "fwd" else "1b"
else:
    read_subtype = "2"
```

In `make_consensus_record()` and the Type-2 record-writing path, emit:
```python
rec.set_tag("XY", ri.read_subtype, value_type="Z")
```

Add `read_subtype` column to `clusters.tsv` (from the first read in each
cluster) and `isoforms.tsv`.

**Do NOT** add a new `read_type` integer value. `read_type` stays as 1 or 2.
`read_subtype` is the human-readable label. Both 1a and 1b route identically
through UMI clustering.

---

### Task 2 — Add consensus pre-trimming before re-alignment

**Why:** Without pre-trim, all three aligners see SSP/UMI/GGG at the 5' and
poly-A/adapter at the 3'. This degrades junction calls near the TSS and CPA —
exactly where accuracy matters most. The DRS track pre-trims poly-A before
alignment; the cDNA track needs both a 5' and 3' trim.

**5' trim — SSP/UMI/GGG strip:**

The GGG/CCC boundary is already determined during `parse_read_info()` as the
`pos5_corrected` field, but that is a reference coordinate. The query-level
trim boundary is the `tss_qpos` computed inside `walk_forward_tss()` (or the
equivalent query position for orient=rev: the position right after SSP_RC in
the read). You need to:

1. During `parse_read_info()`, also record the query-level SSP/UMI/GGG trim
   position (call it `q_trim_5`) — this is the query index of the first mRNA
   base after the GGG bridge. For orient=fwd it's `p + len(SSP_FWD) + UMI_LEN
   + 3` (the +3 is the GGG bridge). For orient=rev it's the analogous position
   before SSP_RC.
2. Before `realign_consensus()`, slice the consensus string:
   `trimmed = consensus_seq[q_trim_5:]` for fwd,
   `trimmed = consensus_seq[:q_trim_5_rev]` for rev.
3. Pass `trimmed` to the aligners instead of the full `consensus_seq`.
4. When building the output BAM record, add back the soft-clip for the trimmed
   prefix/suffix so the CIGAR is valid.

**3' trim — poly-A + adapter strip:**

Re-use the same logic already used in the DRS track. Look at how `run_command.py`
calls the poly-A pre-trim before alignment and restores `XA` afterward.
Specifically:
- Find the poly-A start position in the trimmed consensus using
  `detect_full_length_tier()` / `_find_adapter_anchor_pos()` (already in this
  file).
- Strip from that position to the end; record the count of A's stripped as
  `pretrim_pa_len`.
- After multi-aligner + walkback, set `XA:i` = walkback-corrected tail length
  (if walkback ran) or `pretrim_pa_len` (fallback).

**Test:** add a fixture in `tests/test_cdna_correct_command.py` where the raw
consensus is `"TTTCTGTTGGTGCTGATATTGCT" + UMI + "GGG" + "ATGCATGC" + "AAAAAAAA" + adapter`.
Assert that the sequence passed to the aligner is exactly `"ATGCATGC"` (8 nt mRNA
body only).

---

### Task 3 — Add mapPacBio + gapmm2 to `correct-cdna`

**Why:** DRS uses three aligners (minimap2 + mapPacBio + gapmm2). cDNA should
match after pre-trimming lands clean mRNA sequence on the aligners.

**Reference implementations:**
- `rectify/core/multi_aligner.py` — `run_map_pacbio()` and `run_gapmm2()` both
  exist and are used by the DRS track. Read these before implementing.
- `rectify/core/chimeric_consensus.py` — `select_best()` — already
  protocol-agnostic.
- `rectify/core/run_command.py` lines ~140–220 — the DRS triple-aligner call
  pattern. Mirror this for cDNA.

**What to change in `cdna_correct_command.py`:**

Replace the single `realign_consensus(poa_seq, mp_aligner)` call (~line 1229)
with:
```python
mm2_hit = realign_consensus(trimmed_seq, mp_aligner)
mpb_hit = run_map_pacbio_on_consensus(trimmed_seq, ref_path, chrom, ...)
gmm_hit = run_gapmm2_on_consensus(trimmed_seq, ref_path, chrom, ...)
best_hit = chimeric_consensus.select_best(mm2_hit, mpb_hit, gmm_hit)
```

mapPacBio and gapmm2 write to temp BAM files; wrap them to return the same
`(seq, cigar, ref_start, qual)` tuple format as `realign_consensus()`, or
adapt `select_best()` to consume their native output directly (check how DRS
does it).

Re-alignment applies to **all** consensus reads — singletons (XR=1) and
multi-read clusters (XR≥2) alike.

**Test:** add a synthetic consensus where minimap2 miscalls a splice junction
but mapPacBio rescues it; verify `select_best` returns the mapPacBio alignment.

---

### Task 4 — Move isoform analysis to a downstream step

**Why:** Currently `assign_isoforms()` and `reconcile_t1_t2_pairs()` run on
minimap2-only alignments. Moving them downstream means they run on the
best-available alignment (multi-aligner consensus winner), giving more accurate
TSS and CPA positions for isoform grouping. It also makes the cDNA analysis
pipeline structurally parallel to DRS.

**Scope of this task:**
- Remove the `assign_isoforms()` and `reconcile_t1_t2_pairs()` calls from the
  per-chromosome processing loop in `cdna_correct_command.py`.
- The `correct-cdna` output BAM retains all per-read tags (XU, XC, XR, XF,
  XA, XG, XS, XY, XL) **except** isoform ID (XI) and t1t2 pair info — those
  move downstream.
- Create a new subcommand `rectify cdna-analyze` (or extend `rectify analyze`
  with a `--cdna` mode) that reads the corrected BAM and produces:
  - `clusters.tsv` — per-cluster manifest
  - `isoforms.tsv` — isoform-level aggregation
  - `t1t2_pairs.tsv` — same-molecule strand pairing results
  The functions `assign_isoforms()` and `reconcile_t1_t2_pairs()` move into
  this new subcommand unchanged; only their call site changes.
- Update `run_v4_h2.py` on H2 (at `/u/scratch/k/kevinroy/ont_cdna_v4/`) to
  call the new downstream step after merging per-chrom BAMs.

**Important:** `clusters.tsv` currently includes per-chrom sub-files that are
concatenated. The merge/concat step in `run_v4_h2.py` already handles this.
If clusters.tsv moves to the downstream step, the per-chrom step no longer
writes it — update the concat logic accordingly.

---

### Task 5 — Fix stale Type-2 docstring

In `cdna_correct_command.py`, the v1.15 changelog block (~line 27) says:
> "position-based deduplication (no UMI available)"

Replace with:
> "independent observations carried forward without deduplication (no UMI
> available; reads are grouped by CPA position for isoform counting only)"

---

## Key guardrails — do NOT do these

1. **Do not change `read_type` integer values.** `read_type=1` and
   `read_type=2` remain. Clustering logic checks `read_type == 1` in ~20
   places; changing the integer breaks them all silently.

2. **Do not introduce `Type1a` or `Type1b` as symbol names.** Use
   `read_subtype="1a"` / `"1b"` / `"2"` (strings) and the `XY:Z` BAM tag.

3. **Do not rename `reconcile_t1_t2_pairs` to something else internally.**
   The function name is fine; only its call site moves. The user-facing name
   in TSV headers and README is "same-molecule strand pairing."

4. **Do not remove backward-compatible columns from TSV outputs.** Add
   `read_subtype` as a new column; do not rename or drop existing columns.

5. **Do not skip the pre-trim work before adding multi-aligner.** The
   aligners must receive clean mRNA sequence (no SSP/UMI/GGG, no poly-A).
   Adding mapPacBio + gapmm2 on untrimmed sequence will degrade rather than
   improve junction calls.

6. **Do not implement Task 4 without first having Tasks 1–3 passing tests.**
   The isoform move is the riskiest change; the pre-trim + multi-aligner work
   is a prerequisite.

---

## Files to read before writing any code

| File | Why |
|------|-----|
| `rectify/core/cdna_correct_command.py` | Main target — read in full |
| `rectify/core/multi_aligner.py` | `run_map_pacbio()`, `run_gapmm2()` implementations |
| `rectify/core/chimeric_consensus.py` | `select_best()` signature and contract |
| `rectify/core/run_command.py` | DRS triple-aligner pattern to mirror (lines 140–220) |
| `tests/test_cdna_correct_command.py` | Existing test patterns to follow |
| `PROPOSED_README_UPDATE.md` | Full biology spec, naming decisions, figure specs |

---

## Key biology facts (do not re-derive — these were settled in the design session)

- **Type-1a** (`read_type=1, orient="fwd"`): SSP-first read. SSP+UMI+GGG at
  basecalled-5'. Full-length molecule, UMI identified at read start.
- **Type-1b** (`read_type=1, orient="rev"`): pA-first read that traveled far
  enough to capture SSP_RC+UMI_RC at basecalled-3'. Also UMI-anchored.
- **Type-2** (`read_type=2`): pA-first read truncated before reaching the UMI.
  3' poly-A is the only reliable anchor. **Not deduplicated** — each read is
  an independent molecule observation.
- **Same-molecule strand pairing** (not "reconciliation"): links a Type-2
  cluster to its Type-1 partner cluster from the same physical dsDNA molecule,
  using both 5' TSS and 3' CPA proximity (|Δ5'| ≤ 5 bp AND |Δ3'| ≤ 5 bp).
- **UMI merging asymmetry:** all reads in a UMI cluster (Type-1a + Type-1b)
  merge regardless of 3' truncation. Type-2 reads at the same CPA site are
  NOT merged — without a UMI they cannot be confirmed as the same molecule.
- **gapmm2 is kept** in the DRS aligner set and should also be in cDNA.
  It adds terminal-exon homopolymer refinement on top of minimap2/mapPacBio.

---

## Deliverables

When complete:
1. All existing tests pass (`pytest tests/`)
2. New tests added for pre-trim (Task 2) and multi-aligner rescue (Task 3)
3. `rectify correct-cdna --help` still works
4. A dry-run on chrI of any cDNA sample (e.g. the bundled test data) completes
   without error and the output BAM has `XY:Z` tags
5. `PROPOSED_README_UPDATE.md` implementation gaps 1–3 can be crossed off
   (Task 4 — isoform move — may be a separate PR if scope is too large)
