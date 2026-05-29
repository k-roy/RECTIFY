# Plotter Session Handoff — 2026-05-18 evening

For the **next plotter agent**. Read this AND `scripts/validation_data/PLOTTING.md`
before touching `render_read_alignment.py`. PLOTTING.md is the legacy reference
(some sections now stale, see "Layout changes" below); this doc captures
everything that's changed in the 2026-05-18 session.

**Role boundary**: you are the plotter. A separate debugger session owns
production-code fixes (walkback / indel_corrector / scoring / junction
refiner / chimera filter / bam_writer). Don't edit `rectify/core/**` —
log findings to `validation_read_review/cat{N}_*_findings.md` for the
debugger to pick up. The exception was the `=`-SEQ decoder in
`bam_writer.py` (commits in this session), which had to land at the
producer for the plotter to render correct HP-ED. Don't repeat that
pattern lightly — discuss with the user first.

---

## Layout changes shipped in this session (PLOTTING.md is stale on these)

### 1. End-correction summary table — REMOVED (Phase 1 redesign)

The big monospace text table that previously sat above the alignment
rows has been dropped. The metadata it carried (HP-ED, raw ED, eff
cluster, pick badges, per-aligner 5'/3' coords) is now visible directly
in the alignment rows themselves via:

- Per-aligner-grouped painted CIGAR in the overview track.
- Per-aligner 5'-orig/corr/samtools and 3'-orig/corr ↓ markers on the ref
  row.
- Effective-utility info written to the bundle summary TSV (sample-wide
  rollup logged by the debugger's `merge_corrected_tsvs`).

If a future session needs the metadata back, it's still computed by
`load_correction_panel_data()` — just not rendered.

### 2. Per-base view: 5 per-aligner rows (was 3-row mapped/corrected/pA-rest)

`track_sources` now contains one entry per aligner that has a
per-aligner pA-rest BAM in `<bundle>/per_aligner/<aligner>.pA_rest.bam`,
plus an "unrectified minimap2" row sourced from
`validation_reads_dorado_source.bam` (true raw minimap2 output on the
untrimmed read). The "winner: <name>" row renders first, others
alphabetical.

Aligners with identical post-M-expansion-and-coalesce CIGARs are
collapsed into a single row whose label comma-separates the contributing
aligners ("winner: deSALT, gapmm2, minimap2, uLTRA"). The `+1 1 1` vs.
`2X` style encoding differences across aligners are normalized via the
expand + coalesce pipeline before grouping (see
`_split_m_into_eqx` + `_coalesce_adjacent_same_op`).

### 3. Overview track: painted CIGAR with divergence highlighting

Each row's painted CIGAR is per-op-colored:
- `=` light grey, `M` soft blue (when M-expansion couldn't run),
  `X` pink, `D` orange, `I` purple tick, `S/H` faded grey,
  `N` thin grey gap connector with intron-length label.

Labels:
- Length number only (op letter dropped — color carries the type).
- Insertion labels prefixed with `+` to signal "no ref consumed".
- `D` and `X` labels always render (event ops); `M`/`=` only when wide.
- Single-line placement (no anti-collision stagger — tried and reverted
  in earlier session because data-coord width estimates weren't stable).
- Labels match the segment color for `D` and `X` so you can re-associate
  small ticks with their labels.

Divergence highlights — when ≥2 aligners disagree on the op at a ref
position (or pending-insertion-length at that position), the column is
flagged. Adjacent divergent columns coalesce into ranges. Each range
renders as:
- A wide low-opacity amber background patch (`#FFD54F` @ alpha 0.40)
  behind all bars for context.
- A dark amber ribbon (`#F57C00` @ alpha 0.95, ~0.10 row units tall) at
  the top of the row group, full opacity, on TOP of bars — guarantees
  even 1-bp divergences stay perceptible.
- 1-bp ranges expand to a minimum visual width (0.4% of panel) so
  they're spottable in wide overview windows.

Winner row also gets a faint teal box (`#00897B` @ linewidth 1.4) around
its painted CIGAR for instant identification.

### 4. New ↓ markers on the ref row

Beyond the original orig_3p (red) and corr_3p (green):
- **samtools 3p** (brown / `#795548`): where a naive
  `samtools view | awk` would place the read's 3' end. Sourced from the
  unrectified minimap2 BAM (or the Dorado source BAM when available).
- **orig 5p** (faded blue / `#42a5f5`) + **corr 5p** (deep blue /
  `#1565c0`): the original and corrected 5' positions, computed from
  the winner's raw aligner BAM + the rectified BAM.

3-tier label stagger: tier 1 (txt_y_lo=1.20), tier 2 (txt_y_hi=1.34),
tier 3 (txt_y_st=1.50). samtools auto-stagger when within ±5 cols of
orig/corr. 5' markers share the same stagger logic.

ylim bumped 1.55 → 1.75 to fit tier 3.

### 5. Junction zoom panels

- Per-aligner-group rows (same `track_sources` as the per-base view).
- Per-row painted-style is NOT applied in zoom — kept the per-base
  character rendering since zoom IS the base-level inspection panel.
- Insertions now render in the zoom: purple vertical tick at the
  insertion ref position + small `+N` label hanging above the row.
  (Was missing; added in this session.)
- Winner-row teal-box highlight: same as overview, spans the zoom width.

### 6. Tick-row spacing

Both overview-ticks (height_ratio 0.30→0.50) and main-ticks
(0.35→0.55) were bumped so coord numbers have breathing room below the
alignment baseline. Tick mark itself moved up (y=[0.78, 0.98]) and
label moved down (y=0.00).

---

## Files modified in this session

```
scripts/validation_data/render_read_alignment.py   — major refactor (Phase 1)
scripts/validation_data/regen_pa_rest_bundle.py    — persist per-aligner pA_rest BAMs,
                                                      splice pA into raw minimap2 BAM
scripts/validation_data/generate_review_report.py  — pass SUMMARY_TSV + ALIGNER_DIR
rectify/core/bam/bam_writer.py                     — =-SEQ decoder (producer-side fix)
```

The bundle was regenerated; debugger has continued to regenerate it
multiple times since. Bundle artifacts referenced by the renderer:

```
rectify/data/validation/
  validation_reads.bam                             — bundled test set (XV/XG tags)
  validation_reads_dorado_source.bam               — true raw minimap2 alignment
  aligners/validation_reads.<aligner>.bam          — uncorrected per-aligner BAMs
  rectified/rectified_pA_tail_trimmed.bam          — merged-winner corrected BAM
  rectified/rectified_pA_tail_soft_clipped.bam     — merged-winner pA-restored BAM
  rectified/per_aligner_summary.tsv                — comparison summary (HP-ED, eff_group, etc.)
  rectified/per_aligner/<aligner>.trimmed.bam      — per-aligner corrected (hard-clipped)
  rectified/per_aligner/<aligner>.softclipped.bam  — per-aligner corrected (soft-clipped)
  rectified/per_aligner/<aligner>.pA_rest.bam      — per-aligner pA-restored
  rectified/per_aligner/minimap2_unrectified.bam   — pre-RECTIFY baseline (pA-spliced)
```

---

## Findings docs the debugger is consuming

These live in `validation_read_review/` and are the canonical record of
plotter→debugger handoffs. The debugger has already shipped fixes
addressing most of cat1's findings; cat2 and cat3 still have open
items (search the docs for "Update — 2026-05-18, formal cat1–3 review").

```
cat1_walkback_findings.md       — cat1_minus_1 walkback under-shoot,
                                   cat1_plus_1 over-extension policy,
                                   cat1_plus_2 HP-aware ins
cat2_softclip_findings.md       — cat2_plus_1 HP-ED parsimony mis-ranking,
                                   cat2_minus_2 long-deletion-extension policy,
                                   cat2_plus_2 effective-utility feature spec
cat3_junction_findings.md       — cat3 2D-Nop CIGAR (FIXED by 6d2cf59 + 0653172),
                                   stale 3' pileup (FIXED by 75b0338),
                                   mpb 5'-anchor reanchor (cat3_minus_1 + cat3_plus_1),
                                   chimera exemption interaction
cat5_chimeric_findings.md       — BLOCKED: segment provenance tags missing
                                   from regenerated bundle BAMs
```

---

## Open tasks at session-end

### Plotter-owned, open

- **#5 (Cat5 per-segment provenance overlay)** — BLOCKED on tag
  restoration in the bundle (see `cat5_chimeric_findings.md` for the
  unblock path).

### Debugger-owned, open (do NOT touch these from plotter)

- **#11** — Chimera exemption should cover canonical annotated junctions,
  not just 5'-rescued reads. Resolves the cat3 dispute where mpb has the
  cleanest CIGAR + lowest HP-ED but loses winner-selection.
- **#12** — Pre-rescue 5'-anchor reanchor for mapPacBio's
  mismatch/indel-heavy 5' clusters (cat3_minus_1, cat3_plus_1).
- **#14** — Long-deletion-extension policy (cat2_minus_2 + related).
  Shares principle with #12.

### Plotter-owned, completed in this session

- #4: 5' orig/corr markers on render_ref_row.
- #8: Per-aligner effective utility (read-level eff column shipped;
  sample-wide rollup also shipped by debugger commit 75b0338).
- #9: Raw ED column + winner CIGAR footer (raw ED column survived;
  CIGAR footer replaced by painted-CIGAR overview in Phase 1).
- #10: Phase 1 redesign (drop summary table, per-aligner painted CIGAR
  in overview, per-aligner rows in per-base + zoom, identical-CIGAR
  grouping, winner box).
- #13: Insertions rendered in base-level zoom panels.

---

## How to test / iterate

Standard regen + render workflow:

```bash
cd /Users/kevinroy/work/rectify

# Regenerate the bundle (per-aligner BAMs, summary TSV, etc.). Takes
# ~3 minutes. Only needed when input aligner BAMs change or rectify
# correct logic changes.
.venv/bin/python scripts/validation_data/regen_pa_rest_bundle.py

# Re-render PNGs + MD reports. Choose category subset for fast iteration.
.venv/bin/python scripts/validation_data/generate_review_report.py \
    --arm drs --per-category \
    --category cat3_junction \
    --format md

# Output lands in:
#   /Users/kevinroy/Library/CloudStorage/.../validation_read_review/
#     cat3_junction_review.md
#     cat3_junction_review_pngs/cat3_*.png
```

For one-off single-read rendering see the CLI section of `PLOTTING.md`
(still accurate for `render_read_alignment.py` CLI usage; the underlying
behavior just looks different now).

---

## Known issues / wishlist

1. **Label collision in dense CIGAR regions** — labels still pile up
   when multiple short ops cluster within ~3-5 ref bases. Tried
   anti-collision stagger; reverted because data-coord width estimates
   weren't reliable enough. Workable in practice because the eye can
   still follow the bar colors; could be revisited if matplotlib
   exposes a better renderer-time text bbox API.

2. **Junction zoom — painted CIGAR** would be a natural extension if
   anyone wants per-aligner-group painted bars in the zoom too. The
   per-base view is more informative at zoom scale, but if the user
   prefers visual parity with the overview it's a moderate change.

3. **5' markers + zoom rendering** — the 5'-orig/corr ↓ markers only
   appear in the per-base ref row when they fall in the rendered
   window. Zoom panels don't currently echo them. Easy add if useful.

4. **Cat5 provenance** — see `cat5_chimeric_findings.md`. Needs data-
   side restoration first.

5. **Sample-wide artifacts** — `merge_corrected_tsvs` now emits
   `effective_group` + `effectively_matched_winner` columns and a
   per-aligner rollup at `logger.info`. Plotter doesn't surface these
   yet (no sample-wide summary view in the renderer). Could be a new
   `validation_read_review/sample_summary.md` artifact for future
   reviews.

---

## One-line orientation for a fresh agent

You're the plotter. Read this doc + `scripts/validation_data/PLOTTING.md`
+ the cat findings docs. Don't touch `rectify/core/**`. Render PNGs via
`generate_review_report.py`. Log findings to
`validation_read_review/cat{N}_*_findings.md`. The debugger session
sees your findings via that folder; they'll ship pipeline fixes and
re-regen the bundle on their own cadence.
