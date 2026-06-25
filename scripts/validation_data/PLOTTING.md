# Validation-read plotting reference

Per-read alignment visualization for the RECTIFY validation bundle. Pure
pysam + matplotlib — no IGV, no screencapture. ~0.5–2 s per read.

This doc is the operating reference for the **plotter session**: a Claude
instance whose sole job is rendering. A separate **debugger session**
owns production-code fixes. Read this end-to-end before touching the
rendering code — most layout choices have already been visual-iterated.

**Role boundary.** The plotter does not edit `rectify/core/**`. When you
spot something that looks wrong (a CIGAR that doesn't match the TSV, a
missing N-op, a soft-clip pattern that contradicts the trim parquet),
log it in
`<google_drive>/Chanfreau Lab/validation_read_review/cat{N}_*_findings.md`
for the debugger to pick up. The plotter regenerates PNGs whenever the
bundle changes.

---

## Files

| File | Role |
| --- | --- |
| `scripts/validation_data/render_read_alignment.py` | Single-read renderer. Library + CLI. |
| `scripts/validation_data/generate_review_report.py` | Wrapper: loops the renderer over a category and writes a per-category `.md` report. |
| `scripts/validation_data/audit_polya_trim.py` | Audit: runs `find_polya_and_adapter` on every read in the dorado-source BAM and emits a TSV + invariant flag summary. |
| `scripts/validation_data/restore_polya_from_parquet.py` | Re-attaches the poly(A) tail to `rectified_pA_tail_soft_clipped.bam` (parent of the pA-rest track). NOT a renderer — pipeline plumbing. Accepts both new and legacy "splice" function names. |
| `scripts/validation_data/regen_pa_rest_bundle.py` | End-to-end bundle regen via `rectify correct` (Module 2H) + filter + parquet pA restore. Use when the renderer's input BAMs need refreshing. |

Output lives at
`<google_drive>/Chanfreau Lab/validation_read_review/`:

- `cat{N}_*_review.md` (one per category)
- `cat{N}_*_review_pngs/<xv>.png` (one PNG per read, embedded by the `.md`)
- `STATUS.md` (folder landing page; points back to this doc)

The renderer **no longer emits HTML reports**. The `.md` + embedded PNG
flow renders cleanly in Markdown viewers and Google Drive's preview;
HTML had been a 2026-05 detour and was retired when the rectify-correct
`.md` output stabilized.

---

## What the renderer produces

A single PNG per read, vertical stack:

1. **Title** + figure-level **color legend** (single horizontal strip:
   `= match · X mismatch · D deletion · I insertion · S soft-clip · N intron`).
   M is intentionally absent — `_split_m_into_eqx` expands every M op
   into runs of `=`/`X` before painting. If you ever see loud-yellow
   (`#FFEB3B`) segments, the bundle is stale.
2. **Overview painted-CIGAR rows** — one row per aligner GROUP. Aligners
   with identical post-M-expansion-and-coalesce CIGARs collapse into a
   single row whose label comma-separates the contributing aligners
   (`winner: deSALT, gapmm2, minimap2, uLTRA`). The "winner: …" row is
   bolded; no separate accent shape.
3. **EER ED + raw ED columns** anchored on the display-right of the
   overview, one value per group row. ("EER ED" = empirical error-rate
   ED, the metric formerly called "HP-ED" — renamed because the metric
   covers more than homopolymers.)
4. **Cross-chrom mini-panel(s)** (only when an aligner placed the read
   on a different chromosome — most commonly deSALT picking a paralog).
   Each cc aligner gets its own painted-CIGAR strip + tick row in ITS
   chromosome's coords, anchored on the body start to line up visually
   with the main overview. The cc row is also listed in a red-italic
   banner at the bottom of the figure.
5. **Junction zoom panel(s)** — one per N-op in the corrected BAM. See
   "Junction zoom" below.
6. **3'-end pileup bar chart** (from `rectified/corrected_3ends.{plus,minus}.bedgraph`)
   with ↓ markers above the ref row for `orig 3p`, `corr 3p`,
   `samtools 3p` (where a naive `samtools view | awk` would place the
   read's 3' end — sourced from the unrectified minimap2 BAM), plus
   `orig 5p` and `corr 5p` when they fall in the rendered window.
7. **Reference row** — colored genomic bases, complemented and
   x-axis-flipped for minus-strand reads so the panel reads in RNA 5'→3'.
8. **Per-aligner alignment rows** (per-base character display). One row
   per aligner GROUP plus a top **"minimap2 (unrectified)"** baseline
   row (see below). Each row renders:
   - matches: light grey background, base char colored by identity
   - mismatches: pink background
   - deletions: orange `-`
   - insertions: purple-bordered pill above the row with a black ▼
     pointing to the exact insertion site; nearby insertions stagger
     vertically. Same pill style in the junction zoom.
   - soft-clip: **two-tone** — grey (`SOFTCLIP_BG #F0F0F0`) for the
     aligner's soft-clip of the trimmed body, green (`SOFTCLIP_PA_BG
     #C8E6C9`) for the bases that `restore_polya_from_parquet.py`
     re-attached (= polya_len + adapter_seq from the trim parquet).
   - N-op (intron): grey `|` separator.
9. **Position-tick row** at the bottom.

### "minimap2 (unrectified)" row

Top-of-the-stack baseline row that shows minimap2's alignment of the
**untrimmed** read (still carrying the pA tail). The body typically
extends 1–3 bp into the genomic A-tract before soft-clipping the rest of
the tail — the visual motivator for RECTIFY. Sourced per-read:

1. Preferred: `rectify/data/validation/validation_reads_dorado_source.bam`
   (Dorado uses minimap2 internally, so this IS minimap2's natural
   behavior on the full read).
2. Fallback: `rectify/data/validation/rectified/per_aligner/minimap2_unrectified.bam`
   (the pA-restored variant of the trimmed alignment — used when the
   dorado-source BAM doesn't contain this read; body won't extend into
   the A-tract here, so the boundary-extension nuance is invisible for
   those reads).

The fallback is per-read (file-existence is not enough — the dorado BAM
has historically gone stale before the bundle did; the renderer probes
for the read and falls through if absent). Two reads currently use the
fallback because `build_dorado_source.py` Phase 2 pulled them from a
post-trim source — `cat1_plus_1` and `cat1_minus_1`.

The **pA-tail green coloring is suppressed for this row** — minimap2 was
aligning the full read with the pA still attached, so its soft-clip is
not the result of a rectify trim + restore. The whole soft-clip renders
grey here.

---

## Color palette

Module-level constants at top of `render_read_alignment.py`:

```python
MISMATCH_BG       = "#FFD0D0"
MATCH_BG          = "#FAFAFA"
SOFTCLIP_BG       = "#F0F0F0"   # aligner soft-clip of trimmed body
SOFTCLIP_FG       = "#999999"
SOFTCLIP_PA_BG    = "#C8E6C9"   # pA tail (restored from trim parquet)
SOFTCLIP_PA_FG    = "#2E7D32"
DELETION_BG       = "#FFE0B2"
INTRON_BG         = "#E0E0E0"
```

Overview painted-CIGAR palette (`CIGAR_OV_COLOR` inside `render_overview`):
- `=` neutral grey `#cfd8dc`
- `X` pink `#e91e63`
- `D` orange `#ff9800`
- `I` purple tick `#9c27b0`
- `S` faded grey `#bbbbbb`
- `N` thin grey gap connector `#aaaaaa`
- `M` loud yellow `#FFEB3B` — should never appear; visible as alarm if a
  stale BAM slips through.

---

## Module-level toggles

```python
SHOW_DIVERGENCE_SHADING = False
```

When `True`, the overview track overlays an amber background patch +
darker amber top-edge ribbon over every ref column where ≥2 aligner
groups disagree on the op (or pending-insertion length). Off by default
because the visual gets noisy on high-ED cat1/cat2 reads. Flip to `True`
for ad-hoc reviews where divergence extent is the question.

---

## Strand-aware axis layout

Minus-strand reads use `ax.invert_xaxis()` so the RNA 5'→3' direction
reads left-to-right (higher genomic coord on the left, lower on the
right). To keep the ED-column padding (and the divergence-shading
ribbon) on the display-right after that flip, `render_overview` and
`render_axis_ticks` accept an `is_reverse` parameter and extend xlim to
the LEFT of the bars (`xlim=(-n*0.14, n)`) instead of the right
(`xlim=(0, n*1.14)`).

`suppress_last_tick` flips polarity correspondingly — on plus strand it
drops `tick_xs[-1]` (the rightmost tick, which would bleed under the ED
columns); on minus strand it drops `tick_xs[0]` (which is the position
that lands closest to the ED columns after invert).

Cross-chrom mini-panels are excluded from `invert_xaxis()` — their data
is on a different chromosome and main-read strand convention doesn't
apply. They always use the plus-strand layout.

---

## Junction zoom panel

For every N-op in the corrected BAM's CIGAR, render one zoom panel above
the 3'-end pileup. Visual layout (158 columns wide):

```
 0 .. 49     upstream exon end          (ref [donor-50, donor))
50 .. 74     intron stub donor side     (ref [donor, donor+25))
75 .. 82     gap region                 (intron length label, no bases)
83 .. 107    intron stub acceptor side  (ref [acceptor-25, acceptor))
108 .. 157   downstream exon start      (ref [acceptor, acceptor+50))
```

Constants:

```python
ZOOM_EXON_FLANK = 50   # bp of exon shown on each side
ZOOM_INTRON_STUB = 25  # bp of intron shown on each side
ZOOM_GAP_VISUAL = 8    # visual columns for the "≈ NNN bp" label
ZOOM_TOTAL_COLS = 158  # derived
```

- Exon bases use the same look as the per-base ref row.
- Intron-stub bases get a `#EEEEF6` background + italic glyph.
- Donor + acceptor markers (`↓` triangles) sit at the exon-intron
  boundary with the splice-motif dinucleotide label
  (`donor GT`, `acceptor AG`).
- Intron-length label `≈NNN bp` centered in the gap with a rounded
  yellow background.
- Insertion pills carry over from the per-base row.

**Strand correctness:** for plus-strand reads, RNA donor = genomic
intron_start (lower coord), RNA acceptor = genomic intron_end-1 (higher
coord). For minus-strand reads the roles flip; donor and acceptor
markers + the genomic→visual mapping are re-keyed to place the RNA
donor on the (RNA-)5' side. See `render_zoom_ref_row()` for the exact
convention.

**Multi-junction reads:** one stacked panel per N-op. Cat5/Cat9 reads
with 3–4 junctions get 3–4 panels.

---

## orig/corr label overlap handling

Most reads have `orig_3p == corr_3p` (zero-shift correction) or
`|orig - corr| ≤ 3`. Naive rendering puts both text labels at the same x.

| Condition | Behavior |
| --- | --- |
| `orig_3p == corr_3p` | Single olive `orig=corr` marker. |
| `0 < \|orig - corr\| ≤ 3` | Stacked labels: `orig` at `txt_y_lo=1.20`, `corr` at `txt_y_hi=1.34`. Both triangles still at `tri_y=0.95`. |
| `\|orig - corr\| > 3` | Render both at `txt_y_lo`. |
| `samtools 3p` within ±5 cols of orig/corr | Stagger to `txt_y_st=1.50` (third tier). |

Ref-row `ylim` is `1.75` to fit tier 3.

---

## Bundle paths

```
rectify/data/genomes/saccharomyces_cerevisiae/
    S288C_reference_sequence_R64-5-1_20240529.fsa.gz       <- genome
rectify/data/validation/
    validation_reads.bam                  <- bundle (XV/XG-tagged 36 reads)
    validation_reads_dorado_source.bam    <- untrimmed minimap2 baseline
    aligners/validation_reads.<aligner>.bam               <- per-aligner raw
    rectified/
        rectified_pA_tail_trimmed.bam     <- merged-winner corrected (no pA)
        rectified_pA_tail_soft_clipped.bam <- merged-winner + pA restored
        corrected_3ends.{plus,minus}.bedgraph
        per_aligner_summary.tsv           <- HP-ED, eff_group, pick info
        per_aligner/
            <aligner>.trimmed.bam         <- corrected (hard-clipped)
            <aligner>.softclipped.bam     <- corrected (soft-clipped)
            <aligner>.pA_rest.bam         <- corrected + pA restored
            minimap2_unrectified.bam      <- pre-RECTIFY baseline (pA-spliced)
scripts/validation_data/rebuild_2026_05/trimmed/
    validation_reads_polya_trim_metadata.tsv   <- pA + adapter per read
```

---

## CLI

### Per-category report (the common case)

```bash
.venv/bin/python scripts/validation_data/generate_review_report.py \
    --arm drs --per-category --format md
```

Writes:

- `<google_drive>/Chanfreau Lab/validation_read_review/cat{N}_*_review.md`
  (one file per category)
- `<google_drive>/.../cat{N}_*_review_pngs/<xv>.png`

Common flags:

- `--category cat3_junction` (repeatable) — limit to specific categories.
- `--format md` — current default; HTML output is deprecated.

### Single-read rendering

```bash
.venv/bin/python scripts/validation_data/render_read_alignment.py \
    --qname 79f61403-cf63-4522-b555-569590dc4304 \
    --chrom chrI --start 143330 --end 143430 \
    --orig-3p 143380 --corr-3p 143380 \
    --mapped-bam     rectify/data/validation/aligners/validation_reads.deSALT.bam \
    --minimap2-bam   rectify/data/validation/aligners/validation_reads.minimap2.bam \
    --winner-label   "deSALT" \
    --corrected-bam  rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam \
    --parestore-bam  rectify/data/validation/rectified/rectified_pA_tail_soft_clipped.bam \
    --genome-fa      rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz \
    --bg-plus        rectify/data/validation/rectified/corrected_3ends.plus.bedgraph \
    --bg-minus       rectify/data/validation/rectified/corrected_3ends.minus.bedgraph \
    --aligner-bam-dir rectify/data/validation/aligners \
    --summary-tsv    rectify/data/validation/rectified/per_aligner_summary.tsv \
    --title "cat3_plus_2 — diagnostic view" \
    --out /path/to/output.png
```

Window guidance: pick a 100 bp window centered on `corr_3p`. The
overview panel auto-renders for N-op-spanning alignments or windows
> 200 bp. Junction-zoom panels are independent of window choice.

Suppress junction zoom: `--no-junction-zoom`.

### Audit

```bash
.venv/bin/python scripts/validation_data/audit_polya_trim.py \
    rectify/data/validation/validation_reads_dorado_source.bam \
    --output "<google_drive>/Chanfreau Lab/validation_read_review/polya_trim_audit.tsv"
```

Emits a 36-row TSV (`XV`, `qname8`, `strand`, `last_50bp_RNA`,
`polya_len`, `adapter_seq`, `pass`, `trim_applied`, `raw_last_base`,
`post_trim_boundary`, `source`) plus a summary of four invariant flags
(A: adapter detected + trim skipped; B: trim_applied + post_trim
boundary == 'A' — the "trim to first non-A" violation; C: Pass 2
triggered; D: no stub but raw_last_base non-A). Expected: 36/36
trim_applied, flag B = 0.

---

## Bundle regen (when input BAMs are stale)

```bash
.venv/bin/python scripts/validation_data/regen_pa_rest_bundle.py
# Skip the parquet restore (e.g. when the parquet is missing):
.venv/bin/python scripts/validation_data/regen_pa_rest_bundle.py --skip-polya-restore
```

`regen_pa_rest_bundle.py` Step 5 also rebuilds
`validation_reads_dorado_source.bam` via `build_dorado_source.py` so the
"minimap2 (unrectified)" track always reflects the current 36-read set.
The script handles backups automatically (`*.pre_regen`).

After regen, re-run `generate_review_report.py --per-category` to
refresh the per-category `.md` files and PNGs.

---

## Known visualization quirks (NOT bugs in the renderer)

These are real artifacts of the underlying BAMs / pipeline. The renderer
is faithfully displaying what's in the corrected BAM. Don't "fix" them
in the plotter — surface them to the debugger via the findings docs.

1. **Acceptor label reads "TA" instead of "AG" for some deSALT-winning
   Cat3 reads**: the rescue picks the canonical annotated junction but
   `extend_read_5prime_for_junction_rescue` truncates `intron_end` to
   `read.reference_start` (deSALT's pre-rescue alignment start).
   Acceptor canonicality is lost. Fix is in `read_edits.py` — honor the
   rescue's `rescued_junction[2]` rather than pinning to ref_start.

2. **Refined upstream exon CIGAR is mangled** into noise (e.g.
   `22D1M21I1M…`) instead of the clean `14=1D9=…` that the refiner
   emits internally. Some downstream pass (`indel_corrector` or
   `realign_exon_blocks` is the likely culprit) re-aligns and produces
   garbage. Debugger task.

3. **"donor GG" / "acceptor TA" for non-canonical reads is HONEST**:
   when the BAM-called N-op coords aren't canonical, the renderer
   correctly shows the actual motif. This is signal, not noise — it
   flags reads where the rescue + refiner haven't converged on
   canonical GT/AG.

4. **Loud-yellow segments in painted CIGAR** = a residual `M` op
   survived `_split_m_into_eqx`. Should never happen post-rectify (the
   `=`-SEQ decoder in `bam_writer.py` ensures the chrom_seq + query_seq
   are populated). If you see yellow: bundle is stale or a fresh aligner
   BAM hasn't been processed.

---

## Renderer source organization

`render_read_alignment.py` reading order (current line numbers — refresh
after refactors):

| Symbol | Line | Purpose |
| --- | --- | --- |
| Color constants (`MISMATCH_BG`, `SOFTCLIP_PA_BG`, …) | ~60–85 | Palette |
| `SHOW_DIVERGENCE_SHADING` | ~65 | Module-level toggle |
| `POLYA_TRIM_TSV` | ~70 | Path to parquet metadata |
| `_load_polya_lengths` | 146 | Cache read_id → len(trimmed_3prime_seq) |
| `_mark_polya_softclip` | 178 | Relabel the restored-chunk portion of the soft-clip |
| `walk_cigar` | 224 | Per-window read base + state extraction. Foundation. |
| `_decode_eq_chars`, `compute_mismatches` | 345, 383 | =-SEQ decode + mismatch detection |
| `load_bedgraph_window` | 397 | 3'-end pileup data |
| `_color_for_base` | 419 | State → (bg, fg) — handles `softclip_pa` |
| `render_alignment_row` | 434 | Per-base row (matches/mismatches/clips/inserts/dels/introns + pA-tail two-tone) |
| `render_ref_row` | 534 | Ref bases + orig/corr/samtools/5p markers |
| `render_bedgraph` | 650 | 3'-end pileup bar chart |
| `render_axis_ticks` | 680 | Bottom tick row (strand-aware xlim) |
| Junction-zoom helpers | 747–786 | `_zoom_visual_x`, `_zoom_ref_positions`, `_collect_corrected_junctions` |
| `_coalesce_adjacent_same_op`, `_split_m_into_eqx` | 814, 829 | CIGAR normalization for the overview |
| `_cigar_raw_edit_distance` | 917 | raw ED for the ED columns |
| `load_correction_panel_data` | 1018 | Per-aligner ED + winner info from summary_tsv |
| `render_correction_panel` | 1215 | Legacy text panel — currently unused (`panel_data` is always None at render call site) |
| `render_zoom_ref_row` | 1419 | Junction-zoom ref row (RNA-strand-correct donor/acceptor) |
| `render_zoom_alignment_row` | 1500 | Junction-zoom alignment row (one per track per junction) |
| `render_zoom_ticks` | 1618 | Junction-zoom tick row |
| `render_overview` | 1642 | Schematic painted-CIGAR overview (strand-aware xlim; ED columns; divergence shading gated by toggle) |
| `needs_overview`, `overview_window` | 2081, 2092 | Auto-decide overview window |
| `render` | 2101 | Top-level orchestrator, builds gridspec, dispatches per-panel |
| `main` | 2869 | CLI |

---

## Open plotting items (for the next plotter session)

1. **Cat5 per-segment provenance overlay** — DEFERRED. Needs data-side
   restoration of the chimeric segment-provenance tags
   (`Xa`/`Xg`/`Xz` per the validation-BAM tag schema). See
   `validation_read_review/cat5_chimeric_findings.md` for the unblock
   path. Once the debugger emits per-segment provenance (either as BAM
   tags or as a sidecar TSV), the plotter overlay is a moderate
   addition.

2. **Splice motif tier visualization** — color the donor/acceptor
   labels by canonical_tier (green = GT/AG, amber = GC/AG, red =
   non-canonical). Currently all the same purple regardless of
   canonicality.

3. **5' markers in zoom panels** — the ref-row 5' orig/corr ↓ markers
   don't echo into junction-zoom ref rows. Easy add.

4. **Insertion-pill stagger inside zoom panels** — `render_zoom_alignment_row`
   has the pills but the per-base anti-collision stagger logic isn't
   reused. Rare in practice but worth noting if multiple insertions
   cluster within a zoom window.

5. **`--divergence-shading` CLI flag** — currently a module constant
   (`SHOW_DIVERGENCE_SHADING`). If you want per-render control without
   editing source, plumb a flag through `generate_review_report.py` →
   `render` → `render_overview`.
