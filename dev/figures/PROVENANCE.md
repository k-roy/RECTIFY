# RECTIFY Figure Provenance

**Last updated:** 2026-04-30
**Author:** Kevin R. Roy / Claude (Cowork Desktop App)

All figures in this directory are the canonical source of truth for RECTIFY's
README and documentation. The SVGs are version-controlled; PNGs are rendered
derivatives.

---

## Design System

All hand-crafted figures follow a shared visual language:

- **Standard width:** 760px (viewBox-based; height varies per figure)
- **Font stack:** `Helvetica Neue, Helvetica, Arial, sans-serif`
- **Nucleotide color palette:**
  - A = `#166534` (green) on `#dcfce7` background
  - T = `#9a3412` (orange) on `#ffedd5` background
  - G = `#1e40af` (blue) on `#dbeafe` background
  - C = `#6b21a8` (purple) on `#f3e8ff` background
- **Box dimensions:** 26x24px per nucleotide, 3px border-radius
- **Error highlighting:** `#fee2e2` background, `#dc2626` stroke
- **Correction highlighting:** `#ccfbf1` background, `#0d9488` stroke
- **Section headers:** colored bold text with letter-spacing
- **Separator lines:** `#e2e8f0`, 1px

---

## Figure Inventory

### 1. indel_correction.svg (760x370, 10,250 bytes)

**Title:** 3' End Walk-Back: Recovering the True CPA Site
**Type:** Hand-crafted SVG
**Shows:** Three rows (Genome / Read / Corrected). Read has a deletion in an
A-tract causing the 3' end to overshoot past the true CPA site. Walk-back
arrow shows skipping A's, deletion, and T error to recover the CPA.

**Key elements:**
- Genome row: 5 non-A bases + 9-base A-tract (green borders) with one G error (red border) + 3 non-A bases (amber border)
- Read row: matching alignment with deletion shown as "-" at the G position, T error at end
- Walk-back arrow (blue): labeled "skip A's, deletion, T error -> stop at C"
- Corrected row: trimmed to CPA site with green marker

---

### 2. polya_pretrim.svg (760x420, 12,954 bytes)

**Title:** DRS Poly(A) Pre-Trimming
**Type:** Hand-crafted SVG
**Shows:** Three-pass poly(A) + adapter detection pipeline (RNA 5'->3' orientation).

**Key elements:**
- Pass 0 (green): Pure A scan from 3' end, stops at first non-A (T in this case)
- Pass 1 (amber): Regex `T[CT]{0,10}$` strips adapter stub, then scans A's
- Pass 2 (red): Iterative peel rescue for error stubs that break the regex (G in stub)
- Output: Trimmed read ready for alignment, ambiguous T noted
- Yellow-highlighted T at boundary with "seq error or genomic?" annotation
- Note: "poly(A) length, adapter sequence, and pass number saved to metadata"

---

### 3. softclip_rescue.svg (760x350, 15,693 bytes)

**Title:** Soft-Clip Rescue at Homopolymer Boundaries
**Type:** Hand-crafted SVG
**Shows:** Before/after correction of a read whose 3' end is misplaced due to
nanopore homopolymer errors in a T-tract.

**Key elements:**
- Reference: 11-bp T-tract + C G A
- Read (before): Only 8 T's called, then C A A... soft-clipped — apparent 3' end at position 8
- Corrected: 8 aligned T's + 3 rescued T's (teal border) + C (CPA base) + poly(A) tail (blue border)
- "+4 bp rescued (3 T's + CPA base)" annotation

---

### 4. false_junction_walkback.svg (760x360, 8,403 bytes)

**Title:** False Junction Walk-Back
**Type:** Hand-crafted SVG
**Shows:** A read aligned with a false N-op (skip) in an A-tract region,
and how walk-back discards the false junction.

**Key elements:**
- Genome: non-A bases + A-tract + gap + A-tract + non-A bases
- Read: same but with red N-op box bridging the A-tract gap, T error at end
- Walk-back arrow: "skip A's + T error, discard false N"
- Corrected: trimmed to true CPA at the start of the first A-tract

---

### 5. 5prime_junction_rescue.svg (760x560, 6,792 bytes)

**Title:** Unified 5' End and Junction Correction
**Type:** Hand-crafted SVG (2-panel figure)

**Part 1 — 5' Soft-Clip Junction Rescue (Cat3):**
- Genome: Exon 1 + intron (GT...AG) + Exon 2
- Read: soft-clipped region (matches Exon 1) + aligned to Exon 2
- Semi-global NW alignment rescues bases into Exon 1
- Corrected: rescued Exon 1 block + N (splice) + Exon 2

**Part 2 — N-op Junction Refinement (Module 2H):**
- Read with junction-proximal mismatches near annotated N-op
- HP-aware re-scoring: "sequence match > stability > GT-AG > annotation > shift"
- Corrected: junction shifted to sequence-optimal position with canonical GT-AG
- Note: "only N-ops with a mismatch or indel within 10 bp of the junction boundary are re-scored; clean alignments are skipped"

---

### 6. multi_aligner_consensus.svg (760x680, 12,937 bytes)

**Title:** Multi-Aligner Rectification Pipeline
**Type:** Hand-crafted SVG (3-stage figure)

**Stage 1 — Per-Aligner Rectification:**
- Three parallel columns: minimap2 (blue), mapPacBio (amber), gapmm2 (green)
- Each: BAM -> rectify correct -> corrected BAM + 3' ends TSV
- Left annotation: "each aligner independently: 5' junction rescue, 3' walk-back, soft-clip rescue"

**Stage 2 — Consensus Selection:**
- Three corrected BAMs converge into `rectify consensus`
- Priority chain: 5' rescued -> confidence -> 3' agreement -> span -> n_junctions
- Output: consensus BAM + corrected 3' ends TSV (tags: XA, XC, XN)

**Stage 3 — Chimeric Reconstruction (optional):**
- Genome with Exon A / intron 1 / Exon B / intron 2 / Exon C
- minimap2 gets junction 1 right but misses junction 2
- mapPacBio gets junction 2 right but misses junction 1
- Chimeric stitch recovers both junctions

---

### 7. adaptive_clustering.svg (760x420, 11,253 bytes)

**Title:** Adaptive Valley Clustering
**Type:** Hand-crafted SVG (histogram/chart figure)

**Plot area:** x=70..730, y=40..340 (660w x 300h)
**Axes:** Y = 0..50 (3' End Read Counts), X = 0..100 bp (Genomic Position)
**Scale:** x = 70 + pos*6.6, y = 340 - count*6

**Three clusters with histogram bars:**
- Cluster 1 (blue #6c9fd8): positions 15-27, peak at pos 22 (38 counts), 192 reads
- Cluster 2 (orange #e8a848): positions 42-53, peak at pos 49 (26 counts), 170 reads
- Cluster 3 (green #8bc46a): positions 76-87, peak at pos 80 (16 counts), 80 reads

**Markers:**
- Peak markers (red inverted triangles) at pos 22, 49, 80
- Valley markers (cyan upward triangles) at pos 38, 68
- Cluster boundaries (dashed lines) at x=149, 321, 519, 677

**Legend:** Horizontal bar under title with Peak/Valley/Boundary items
**Algorithm note:** "find peaks -> find valleys -> boundary = midpoint(peak, valley) +/- 10 bp cap"

---

### 8. oligo_a_spreading.svg (422.7x242.6pt, 56,512 bytes)

**Title:** (from NET-seq oligo-A analysis)
**Type:** Matplotlib-generated SVG export
**Note:** This is NOT hand-crafted. It was exported from the matplotlib script
`plot_oligo_a_spreading.py` on Sherlock. Uses `pt` units and SVG 1.1 DOCTYPE.

---

### 9. oligo_a_deconvolution.svg (710.4x241.1pt, 84,743 bytes)

**Title:** (from NET-seq oligo-A NNLS deconvolution)
**Type:** Matplotlib-generated SVG export
**Note:** This is NOT hand-crafted. Same provenance as oligo_a_spreading.

---

## PNG Rendering

PNGs are rendered at 2x resolution (retina) using `render_all_pngs.py`:

```bash
pip install cairosvg
python render_all_pngs.py              # all 9 figures
python render_all_pngs.py -f <stem>    # single figure
python render_all_pngs.py -s 3         # 3x scale
```

The script extracts width from each SVG's viewBox and renders at `width * scale`.
For the 7 hand-crafted SVGs (760px width), output is 1520px at 2x.

## Deployment to Sherlock

```bash
rclone copy 'stanford_google_drive:ClaudeSherlockCode/rectify_figures/' \
    /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/docs/figures/ \
    --include '*.png' --include '*.svg' -v
```
