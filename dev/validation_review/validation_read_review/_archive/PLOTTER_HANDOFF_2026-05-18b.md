# Plotter Session Handoff — 2026-05-18 (evening, second session)

For the **next plotter agent**. Read this FIRST, then `PLOTTER_HANDOFF_2026-05-18.md`
(the earlier handoff this same day — describes Phase 1 redesign), then
`scripts/validation_data/PLOTTING.md` (legacy reference; partial sections
still useful).

**Role boundary unchanged**: you are the plotter. The debugger session owns
`rectify/core/**` fixes — don't edit there. Log findings to
`validation_read_review/cat{N}_*_findings.md`. The debugger has continued
to ship significant fixes in parallel during this session (commits between
`14d9add` and `0653172` cover walkback / 5'-rescue / chimera flag logic
across cat1–3) — re-pulling and re-regen-ing the bundle is the standard
opening move.

---

## TL;DR for someone walking in cold

The renderer is mature. Most surface-level layout choices are
visual-iterated, documented in section 6–9 below. **Three things still
worth doing**:

1. **cat3_plus_1 / cat3_minus_1 mapPacBio rescue** is queued for the
   debugger (task #12). Plotter-side it's surfaced via the painted
   CIGAR + the divergence band; nothing more to do until the debugger
   ships.
2. **Cat5 segment-provenance overlay** (task #5) is blocked on the
   debugger restoring chimeric tags to the bundle. Detailed unblock
   path in `cat5_chimeric_findings.md`.
3. **deSALT cross-chrom alignment perception**: the math is correct
   (aln_start anchor, same window length, same xlim) but the user
   reported the visual alignment looking slightly off in multiple
   iterations. Anchor choice was iterated (aln_start ↔ full_aln_start)
   and current setting is aln_start. May warrant another pass if it
   re-surfaces. See "Known issues" below.

Everything else is shipped + style-reviewed. Run
`generate_review_report.py --category catN_…` and iterate.

---

## What landed in THIS session (after the morning Phase 1 handoff)

### 1. Cross-chrom outlier handling (deSALT-on-paralog case)

cat1_plus_1's deSALT places the read on `chrVI:8526-8704` (paralog) while
the other 4 aligners land on `chrXIV:10435-10611`. Previously the
renderer painted deSALT's row against chrXIV's reference (wrong chrom)
producing a wall of fake mismatches. New behavior:

- **Per-base view filtering**: when an aligner's `reference_name`
  differs from the rendered chrom, that aligner is excluded from the
  per-base aligned_set (avoids cross-chrom base comparison). Tracked
  in `cross_chrom_rows` + `cross_chrom_bams`.
- **Top-of-figure italic banner** lists each cross-chrom outlier with
  its actual coords + ED values, in red italic.
- **Dedicated mini-panel** at the top of the figure (between the
  overview block and the bedgraph): painted-CIGAR row with its OWN
  ref-coord x-axis + tick row, against its actual chromosome. Re-uses
  `render_overview` machinery so all the overview features (M→=/X
  expansion, divergence highlighting, EER/raw ED label columns, etc.)
  work seamlessly. The aligner's CIGAR appears clean against the
  correct genome.
- **Query-position anchoring**: the cross-chrom panel's `cc_win_start`
  is set so the read's `aln_start` (body start) lands at the same axis
  x as the MAIN winner's `aln_start`. Window length matches
  `ov_end - ov_start`, so both axes share the same scale. Anchor
  choice (aln_start vs. full_aln_start) was iterated — aln_start was
  chosen because the BODY is the dominant visual element. Aligners
  with leading soft-clip have their soft-clip extend left of axis x
  =20; bodies still align at x=20 across panels.
- The cross-chrom mini-panel is excluded from the minus-strand
  x-axis-invert pass (its data is on a different chrom — main-read
  strand convention doesn't apply).

### 2. Inline EER ED / raw ED columns in the overview's left margin

Replaced the standalone ED-table panel with two right-aligned text
columns inside the overview's left margin, alongside the row labels.
Column headers ("EER ED", "raw ED") sit at the same y as the
"overview" tag at the top.

- **Renamed** `HP-ED` → `EER ED` (empirical error-rate ED — user's
  preferred term; "HP" was misleading since the metric covers more than
  homopolymers).
- **Single value per row** (NOT comma-separated). Aligners in a group
  share their post-normalization CIGAR by definition, so their EER/raw
  ED are equal — `_eds_for_label` picks the first aligner in each
  group's name list. This was a late-session refinement; previously
  showed `30.3 · 30.3 · 30.3` which pushed the left margin out
  awkwardly.
- **Font size** in ED columns is 6.5pt; row label is 7.5pt. ED columns
  use monospace; row labels use the default sans.
- **Final layout (post-iteration)**: ED columns moved to the RIGHT of
  the painted CIGAR (not the left margin). The overview's xlim is
  extended 14% past the data range (`xlim = (0, n * 1.14)`) so the ED
  values have their own real estate without competing with the bars.
  - Row label at transAxes `x = -0.02` (just left of axis — matches
    where matplotlib y-tick labels sit by default; vertically aligns
    with the zoom-panel labels below).
  - EER ED at data `x = n * 1.07` (right-aligned, fontsize 8).
  - raw ED at data `x = n * 1.13` (right-aligned, fontsize 8).
  - Headers at `y = len(aligned_set) + 0.55` (just above top row).
- **Left margin** back to `0.12` (close to default 0.11) since labels
  no longer claim extra space on the left.
- **Headers suppressed on cross-chrom panels** via
  `render_overview(show_column_headers=False)` — the main overview
  already shows them.
- **Cross-chrom rows** also display ED. The label parser strips a
  trailing `[chrom]` tag before looking up the aligner name in
  `aligner_ed_map`.

### 3. ED-loader bug fix — raw ED on cross-chrom reads

`load_correction_panel_data` previously passed a single `chrom_seq`
(the winning aligner's chrom) into `_cigar_raw_edit_distance` for ALL
aligners. Cross-chrom aligners (deSALT on chrVI) were scored against
the wrong chromosome → inflated raw ED (was reading 158, should be 39
for this read).

Fix: `load_correction_panel_data` now accepts `genome_fa: Path` and
re-opens it per aligner to fetch each read's actual chrom sequence.
Falls back to the passed `chrom_seq` when `genome_fa` is unavailable.

### 4. Polish on alignment-row internals

- **Insertion pills**: dropped the dark purple border (was redundant
  with the ▼ tick + colored letters inside). Light purple fill only.
- **Letter size**: insertion-bases match normal base size (`char_width`
  1.0, was 0.55), with generous padding (`pad` 0.20 each side, was 0.10).
  Pill height bumped 0.55 → 0.70 for vertical buffer.
- **Zoom-panel insertions** now use the same pill style as the per-base
  alignment row (matching colored letters + ▼ tick + stagger). Was the
  user's biggest specific ask on this iteration. ylim of
  `render_zoom_alignment_row` bumped to `(0, 2.60)` and the zoom-row
  height_ratio bumped from 2.0 to 3.2 to fit the pills.

### 5. Divergence-band visibility

Earlier `#FFD54F` @ alpha 0.50 amber background bands were drowned by
the painted bars (bars are zorder=1 with alpha=0.90). Replaced with
a two-layer system:

- **Wide background patch** (low-alpha amber, zorder=0, full row
  height) for tracing extent.
- **Thin top-edge ribbon** (`#F57C00` @ alpha 0.95, zorder=3, ~0.10
  row units tall) ABOVE the bars — guarantees even 1-bp divergences
  stay perceptible.
- 1-bp ranges expand to a minimum visual width (0.4% of panel) so
  they're spottable at wide overview scales.

### 6. Background banding removed; winner = green left-edge accent + bold

Per style-review recommendation (spawned 2026-05-18 evening): full-row
beige + green bands were dropped in favor of a cleaner white
background. Winner now indicated by:
- A **3-px-wide solid green (`#2E7D32`) left-edge accent bar** flush
  with the row label (drawn in transAxes coords, x=-0.185 to -0.173).
- **Bold typeface on the winner's row label** (other rows use normal
  weight).

Same approach in the zoom panels — light-green accent behind the
y-tick label area.

### 7. Painted-CIGAR op-length labels removed + color legend added

User asked to try a version without the op-length numbers (`29`, `13`,
`+2I`, etc.) inside the painted CIGAR — replaced with a single
top-of-figure color legend showing what each segment color means.
Insertion `+N` labels are kept (length info there can't be derived
from the per-base view).

The legend (added via `fig.legend()` at `loc="upper center"`,
`bbox_to_anchor=(0.5, 0.978)`) lists the seven op codes:
`= match`, `M ambig.`, `X mismatch`, `D deletion`, `I insertion`,
`S soft-clip`, `N intron`. One Patch handle per op with its color
swatch + label.

### 8. Font consolidation — Helvetica Neue + Menlo + 3 sizes

- `matplotlib.rcParams['font.sans-serif'] = ['Helvetica Neue',
  'Helvetica', 'Arial', 'DejaVu Sans']`
- `matplotlib.rcParams['font.monospace'] = ['Menlo', 'Monaco',
  'Consolas', 'DejaVu Sans Mono']`
- Module-level constants: `SIZE_TITLE = 11`, `SIZE_TEXT = 9`,
  `SIZE_SMALL = 6.5`. Most labels / headers / ED values / italic
  banners use `SIZE_TEXT`; CIGAR-op-length labels (now removed) and
  the intron-length labels in the overview's gap connectors use
  `SIZE_SMALL`. Per-base character rendering (where each base sits in
  a 1-data-unit column) keeps its explicit `fontsize=8/8.5` to fit
  the column width.

### 7. Tick-row spacing

- Overview ticks: height_ratio 0.30 → 0.50.
- Main ticks: 0.35 → 0.55.
- Tick mark y=[0.78, 0.98] (top of panel), label at y=0.00 (bottom) so
  coord numbers have visible breathing room.

### 8. 5'-end ↓ markers on ref row

Two new markers alongside the existing 3'-end ones:
- `orig 5p` (faded blue `#42a5f5`)
- `corr 5p` (deep blue `#1565c0`)

Sourced from the winner's raw aligner BAM + the rectified BAM. Only
appear when the rendered window happens to include the 5' position
(rare for 3'-end-centered window; useful in overview spans + zoom).

### 9. Cleanup

- Italic subtitle `"per-aligner rows: corrected CIGAR + soft-clipped
  poly(A) tail spliced back from the trim parquet"` lowered from
  y=0.978 → y=0.962 so it no longer overlaps the title.
- Red cross-chrom annotation overlay above the deSALT mini-panel
  removed (was duplicating the top banner + colliding with overview
  ticks).

---

## Files modified in THIS session

```
scripts/validation_data/render_read_alignment.py    — substantial
                                                       (insertions, cross-chrom panel,
                                                        inline ED columns, divergence
                                                        layers, 5' markers, …)
```

No changes to `regen_pa_rest_bundle.py`, `bam_writer.py`, or any
`rectify/core/**` file in this session.

The bundle was regenerated multiple times by the debugger (most recent
mtime 2026-05-18 13:45). The plotter doesn't need to regen — just
re-render whenever the bundle updates.

---

## Findings docs updated this session

All in `validation_read_review/`:

```
cat1_walkback_findings.md       — unchanged
cat2_softclip_findings.md       — added long-deletion-extension policy
                                   (task #14) re cat2_minus_2
cat3_junction_findings.md       — added mpb 5'-anchor reanchor finding
                                   (task #12) re cat3_minus_1 + cat3_plus_1;
                                   chimera-exemption interaction (task #11)
cat5_chimeric_findings.md       — Cat5 segment provenance blocker
                                   (task #5)
PLOTTER_HANDOFF_2026-05-18.md   — morning handoff (Phase 1 redesign)
PLOTTER_HANDOFF_2026-05-18b.md  — THIS doc
```

The debugger has been consuming these and shipping. Specifically
addressed since the morning handoff:
- cat1 walkback fixes (`09e4627`, `a1728eb`, `20aa1c8`)
- cat2 expected-value updates (`c725819`)
- cat3 2D-Nop CIGAR fix (`6d2cf59`, `0653172`)
- cat3_plus_1 silent-False (`62241ea`)
- Bedgraph regen step + per-aligner effective_utility column (`75b0338`)

Still on the debugger's queue (per cat findings docs):
- **#11**: chimera exemption should cover canonical annotated junctions
  (not just 5'-rescued reads). Resolves the cat3 dispute where mpb has
  cleanest CIGAR but loses winner-selection.
- **#12**: pre-rescue 5'-anchor reanchor for mapPacBio's mismatch/
  indel-heavy 5' clusters (cat3_minus_1, cat3_plus_1).
- **#14**: long-deletion-extension policy (cat2_minus_2 + related;
  shares principle with #12).

---

## Open tasks at session-end

### Plotter-owned, still pending

- **#5 — Cat5 per-segment provenance overlay** — BLOCKED on data
  restoration. Cat5 segment-provenance tags (XA/XS/Xz) are missing or
  wrong in the regenerated bundle BAMs (per the tag-rename refactor in
  `9e1b9a1`). Detailed unblock path + rendering plan in
  `cat5_chimeric_findings.md`. Next plotter session can pick this up
  once the debugger restores the tags OR emits a sidecar TSV.

### Debugger-owned (NOT plotter)

- **#11**, **#12**, **#14** (above).

### Plotter-owned, completed in this session

- **#4** — 5' orig/corr markers on render_ref_row.
- **#8** — Per-aligner effective utility tracking (read-level eff
  column + sample-wide rollup shipped by debugger in `75b0338`).
- **#9** — Raw ED column (winner CIGAR footer replaced by painted
  CIGAR overview in Phase 1).
- **#10 Phase 1 + Phase 2** — Per-aligner painted CIGAR in overview
  + zoom + per-base view; identical-CIGAR grouping; divergence
  highlights; winner box; inline EER/raw ED columns.
- **#13** — Insertion pills in junction-zoom panels.
- **#15** — Separate-x-axis cross-chrom mini-panel.

---

## How to test / iterate

```bash
cd /Users/kevinroy/work/rectify

# When the bundle is fresh — usually it already is (debugger reruns it).
# Heavy and only needed when input aligner BAMs or rectify correct
# logic changes. ~3 min on M1.
.venv/bin/python scripts/validation_data/regen_pa_rest_bundle.py

# Re-render PNGs + MD reports for one or more categories. Fast (~10 sec
# per category of 4 reads).
.venv/bin/python scripts/validation_data/generate_review_report.py \
    --arm drs --per-category \
    --category cat1_indel \
    --format md

# Outputs:
#   /Users/kevinroy/Library/CloudStorage/.../validation_read_review/
#     cat1_indel_review.md
#     cat1_indel_review_pngs/cat1_*.png
```

Specific reads useful as test fixtures:
- **cat1_plus_1** (`0cb5a111`): exercises the cross-chrom panel
  (deSALT on chrVI vs others on chrXIV). 5 aligners, 1 cross-chrom.
- **cat3_plus_1** (`0a28167d`): 1 intron, multiple aligner CIGARs
  diverge; exercises divergence bands, insertion pills in zoom
  panels, and winner-box highlight.
- **cat2_minus_2** (`b313b50d`): no introns; clean test of inline
  ED column layout for a 5-aligner-grouped row.

---

(Earlier mid-session layout map removed — see the "Final layout map
(post-session)" section near the end of this doc for the current state.)

---

## Known issues / wishlist (still open)

1. **Label collision in dense CIGAR regions** — *now mostly moot* since
   CIGAR op-length labels were removed (color legend at top replaces
   them). The remaining label-collision risk is on `+N` insertion
   labels when multiple insertions cluster within ~3 ref bp. Lateral
   redistribution attempt was reverted as too fragile.

2. **Zoom-panel painted CIGAR** — zoom panels still render per-base
   characters (which is useful for base-level inspection). The overview
   already uses painted CIGAR. Could add a third "painted-zoom" view
   for visual parity if requested.

3. **5' markers in zoom panels** — the ref-row 5' markers don't echo
   into junction-zoom ref rows. Easy add if useful.

4. **Cat5 provenance** — see `cat5_chimeric_findings.md`. Data side
   needs restoration first.

5. **Sample-wide effective-utility surfacing** — the debugger now
   writes `effective_group` + `effectively_matched_winner` to the
   summary TSV and logs a rollup; the plotter doesn't surface that
   anywhere in the rendered output. Could be a `sample_summary.md`
   artifact next to the per-category review docs.

6. **deSALT cross-chrom alignment perception** — the math is verified
   correct (both panels use the same n and xlim; deSALT's aln_start
   lands at axis x = `main_anchor_offset` = 20 on its panel, same x
   as the main winner's aln_start). User has reported visually-
   perceived misalignment multiple times across iterations. Anchor
   was iterated between `aln_start` (current) and `full_aln_start`
   (= aln_start − leading 5'S). Differences in leading/trailing
   soft-clip lengths between aligners DO cause visible offsets in
   the leftmost/rightmost rendered content (which is biologically
   honest), but the bodies align. If this re-surfaces, options to
   try:
   - Render cross-chrom panels WITHOUT the leading/trailing
     soft-clip patches (body-only).
   - Rescale the cc panel's x-axis non-uniformly so the cc full
     alignment span exactly matches the main winner's span (loses
     per-bp scale equivalence).
   - Show the cc bar in a different visual style (outlined / dashed)
     to make clear it's on a different chromosome and not directly
     comparable bp-for-bp.

7. **Painted CIGAR overview "M" segments** — when the M→=/X expansion
   can't run (no genome OR no query_seq), M blocks render in soft
   blue (`#90caf9`) and are noted in the legend as "ambig.". With
   the bam_writer `=`-SEQ decoder fix, post-rectify BAMs should
   always allow expansion, so blue should be rare; if you start
   seeing lots of blue, suspect a fresh aligner BAM the bundle
   regen hasn't processed yet.

---

## Final layout map (post-session)

Top → bottom panel order in a typical PNG, after all this session's
changes:

```
title  (Helvetica Neue 11pt)
  color legend:  = M X D I S N    (top center, just below title)

— overview ————————————————————————————————————————————————————
overview                                          EER ED   raw ED
| winner: <name>      [painted CIGAR row]         XX.X      NN
| <other group(s)>    [painted CIGAR row]         XX.X      NN
[3-px green ▌ left-edge accent on winner row only; bold name on winner]
[divergence bands (amber background + top ribbon) span ≥2-aligner-
 disagreement ref ranges]
       10,449  10,489 …                          (data coord ticks,
                                                  rightmost suppressed)

— cross-chrom mini-panel(s) (each cc aligner gets its own row) ——
deSALT [chrVI]      [painted CIGAR on chrVI]       XX.X      NN
       8,546  8,586 …                              (chrVI coord ticks)

— 3' end pileup ————————————————————————————————————————————————

— ref row + ↓ markers ——————————————————————————————————————————
ref     [colored genomic bases, RNA-strand-oriented]
        orig=corr ↓  samtools ↓  orig 5' ↓ corr 5' ↓
                       (markers stagger when close in x)

— junction zoom panels (one per N-op) ——————————————————————————
ref (junc N)  [donor ↓ ▼ ≈intron-len bp ▼ ↓ acceptor]
<per-aligner zoom rows with per-base chars + purple insertion pills>
[zoom ticks]

— per-aligner alignment rows (per-base char display) ——————————

— main ticks ———————————————————————————————————————————————————

[fig.text italic banner: "per-aligner rows: corrected CIGAR + ..."]
[fig.text red bold italic: "cross-chrom (excluded from per-base view): ..."]
```

---

## One-line orientation for a fresh agent

You're the plotter. Read this doc + the morning handoff
(`PLOTTER_HANDOFF_2026-05-18.md`) + the cat findings docs. Don't touch
`rectify/core/**`. Render PNGs via `generate_review_report.py`. Log
findings to `validation_read_review/cat{N}_*_findings.md`. The
debugger session sees your findings via that folder; they'll ship
pipeline fixes and re-regen the bundle on their own cadence.
