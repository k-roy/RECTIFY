# RECTIFY Figure Style Guide

**Purpose:** Concrete rules for generating README schematic figures with a consistent, professional aesthetic. This guide documents the design system used in `generate_figures_v2.py` and explains *why* each choice was made, so any agent can reproduce the style.

---

## Side-by-Side: What Went Wrong

Comparing the other agent's `indel_correction.png` (left column) with the correct version (right column) illustrates every problem this guide addresses:

| Issue | Other agent | Correct version |
|-------|------------|-----------------|
| **Corrected row** | Shows 17 identical A boxes stretching across full width — implies bases were modified | Shows only `CAGTC` with a vertical CPA marker — the output is a *position*, not modified sequence |
| **Green success banner** | Full-width rounded rectangle at bottom: "True 3' end recovered . Indel artifacts removed . Poly(A) tail length preserved" | No banner. The figure speaks for itself |
| **Yellow warning box** | "For minus strand: poly(A) appears as poly(T)..." in a yellow rounded rectangle | Not present — strand details belong in documentation, not a schematic |
| **Section framing** | "The Problem:" / "RECTIFY's Solution:" — prescriptive, marketing-tone headers | "ALIGNED READ" / "AFTER WALK-BACK" — descriptive, factual headers |
| **Vertical density** | ~730px tall with large empty regions around the banners | ~400px tall, no wasted space |
| **Nucleotide count** | Read row has 20 boxes (including 7 soft-clip A's), corrected row has 17 all-A boxes | Read row has 13 boxes (pre-trimmed), corrected row has 4 boxes. Minimal example that shows the concept |
| **Annotation collision** | "poly(A) tail start" label overlaps the section header text | Labels positioned with explicit y-offsets, no overlaps |

The `multi_aligner_consensus.png` had a different but related set of problems: it used pastel-colored blocks with no internal labels (you had to read the legend), repeated the genome track three times across two redundant panels, and added a scoring criteria box with a green "Best Alignment Selected" result banner. The correct version uses a flowchart with labeled nodes and one clean chimeric reconstruction example.

---

## Core Principles

1. **Show the minimum example that communicates the concept.** If 5 nucleotide boxes demonstrate the A-tract problem, don't use 12. If the corrected output is a position, show just the position marker — don't render the entire modified sequence.

2. **No success banners, warning boxes, or marketing framing.** These are technical figures for a README, not slides for a sales pitch. "The Problem / The Solution" framing with green checkmark banners looks like a SaaS landing page. Use factual, descriptive section headers instead ("ALIGNED READ", "AFTER WALK-BACK").

3. **Restrained color — 2-3 accent colors per figure.** Every color must earn its place. If a nucleotide box is green-filled, that's because its identity (A) matters to the story. If something is red, it's an error or artifact. Don't use color for decoration.

4. **White space is structure.** Separate sections with a thin divider line and vertical padding — not colored background panels. Let the eye rest between the problem and the solution.

5. **The figure should be readable at 680px width on GitHub.** All text >= 8.5pt. No text that requires zooming. Every label must be legible at the rendered size.

---

## Design Tokens (from `generate_figures_v2.py`)

Every figure must use these exact values. No ad hoc colors or sizes.

### Dimensions

```python
FIG_W = 760          # all figures same width (renders at 680px on GitHub)
BW, BH, BR = 26, 24, 3  # nucleotide box: width, height, border-radius
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"
```

### Color Palette

```python
PAL = dict(
    # --- Text hierarchy (4 levels) ---
    title   = "#1e293b",   # figure title only
    heading = "#475569",   # section headers
    label   = "#64748b",   # row labels (Genome, Read, Corrected)
    muted   = "#94a3b8",   # secondary annotations, intron labels

    # --- Structural ---
    border  = "#cbd5e1",   # default nucleotide box border
    divider = "#e2e8f0",   # horizontal section dividers
    bg      = "#ffffff",   # figure background

    # --- Accents (use sparingly) ---
    blue    = "#2563eb",   # primary accent (arrows, CPA markers)
    blue_l  = "#dbeafe",   # light blue fill
    green   = "#059669",   # correct / rescued elements
    green_l = "#d1fae5",
    red     = "#dc2626",   # errors, artifacts, deletions
    red_l   = "#fee2e2",
    teal    = "#0d9488",   # rescued / corrected regions
    teal_l  = "#ccfbf1",
    orange  = "#d97706",   # warnings, non-A boundaries
    orange_l= "#fef3c7",

    # --- Exon blocks ---
    exon    = "#6366f1",   # indigo fill (distinctive, professional)
    exon_t  = "#ffffff",   # white text on exon
    exon_l  = "#e0e7ff",   # light indigo (aligned regions)
)
```

### Nucleotide Colors

Each base has a muted fill and a dark text color. These are soft pastels, not saturated primaries.

```python
# Fill / Text pairs
A:  "#dcfce7" / "#166534"   # soft green
T:  "#ffedd5" / "#9a3412"   # soft orange
G:  "#dbeafe" / "#1e40af"   # soft blue
C:  "#f3e8ff" / "#6b21a8"   # soft purple
-:  "#fee2e2" / "#dc2626"   # red (deletion)
N:  "#fee2e2" / "#dc2626"   # red (false junction)
```

**Why muted pastels?** Saturated base colors (bright red T, bright green A) create visual noise when you have 10+ boxes in a row. Muted fills let the eye scan across the sequence while still distinguishing bases.

---

## Typography Hierarchy

| Element | Size | Weight | Color | Anchor |
|---------|------|--------|-------|--------|
| Figure title | 14px | bold | `title` (#1e293b) | center |
| Section header | 10px | bold, letter-spacing 0.8 | `heading` (#475569) | left, x=60 |
| Row label | 10px | normal | `label` (#64748b) | left, x=14 |
| Annotation | 8.5-9px | normal or 600 | context-dependent | varies |
| Nucleotide | 11px | 600 | base-specific | center in box |

**Rules:**
- One font family throughout (`Helvetica Neue` stack)
- Never use font sizes below 8.5px
- Section headers are UPPERCASE with letter-spacing for visual distinction without needing a background panel
- Row labels ("Genome", "Read", "Corrected") always at x=14, vertically centered with the row

---

## Layout Rules

### Vertical Spacing

```
Title:           y = 22
First section:   y = title + 30
  Section head:  y (relative to section top)
  First row:     section head + 28
  Second row:    first row + 50 (genome-to-read gap)
  Annotation:    row + BH + 14
Divider:         last annotation + 20
Next section:    divider + 25
```

Minimum 20px between the last element in a section and the divider line. This prevents the cramped feeling.

### Horizontal Layout

```
Row label:       x = 14
First base box:  x = 80-100 (enough clearance for "Corrected" label)
Box spacing:     BW + 2 (= 28px center-to-center)
Right margin:    last box should not exceed x = 700
```

### Section Dividers

Use a single 1px line in `divider` color (#e2e8f0), spanning x=50 to x=710:

```python
def hdivider(y, x1=50, x2=710):
    return f'<line stroke="{PAL["divider"]}" stroke-width="1" x1="{x1}" x2="{x2}" y1="{y}" y2="{y}"/>'
```

**Never use:** colored background panels for sections, rounded rectangles with fills as section containers, or thick separator lines.

---

## Element Catalog

### Nucleotide Boxes

Default: muted fill, 0.7px border in `border` color. Highlighted: accent fill + 1.2px accent border.

```python
def base(x, y, b, highlight_fill=None, highlight_border=None):
    fill = highlight_fill or BASE_FILL.get(b, "#f1f5f9")
    bdr  = highlight_border or PAL["border"]
    sw   = "1.2" if highlight_border else "0.7"
    # rect + centered text
```

**When to highlight:** Only to draw attention to a specific base that matters to the story (e.g., a T sequencing error, a G interruption in an A-tract). If a base is just context, use the default muted styling.

### Braces (Above / Below)

Thin brackets with 0.6 opacity, 8.5px label. Use to annotate a region of bases.

```python
brace_above(x1, x2, y, color, label)   # bracket opens upward
brace_below(x1, x2, y, color, label)   # bracket opens downward
```

**When to use:** To label a contiguous region (e.g., "genomic A-tract", "adapter stub"). Don't stack more than 2 braces on the same row — it gets cluttered.

### Arrows

Horizontal arrows with optional centered label above. 1.5px stroke, solid triangular arrowhead.

```python
h_arrow(x1, x2, y, color, label=None)
```

**When to use:** To show direction of a scan or walk-back operation. One arrow per figure is usually enough.

### Vertical Markers

2px vertical lines marking a position (e.g., CPA site, true 3' end). Label on right or left side.

```python
vert_marker(x, y1, y2, color, label=None, label_side="right")
```

**When to use:** To mark a key genomic position that spans multiple rows (genome + read alignment).

### Exon Blocks

Solid fill with white text, 4px border-radius. Indigo by default.

```python
exon(x, y, w, h, label, color=None)
```

### Intron Lines

Dashed line in `muted` color with centered "intron" label above.

```python
intron(x1, x2, y, label="intron")
```

### Soft-Clip Blocks

Dashed red border, light red fill, red label text. Visually distinct from aligned regions.

```python
softclip_block(x, y, w, h, label="soft-clipped")
```

### Aligned Region Blocks

Semi-transparent (opacity 0.18) fill with 1px solid border. Used for schematic aligned regions when you don't need individual bases.

```python
aligned_block(x, y, w, h, label, color=None)
```

### Rescued Region Blocks

Teal fill/border. Shows regions that RECTIFY recovered.

```python
rescued_block(x, y, w, h, label)
```

---

## Anti-Patterns (Things to Avoid)

### 1. Success / Result Banners

**Bad:**
```
┌─────────────────────────────────────────────────────┐
│  True 3' end recovered  ·  Indel artifacts removed  │
└─────────────────────────────────────────────────────┘
```
A full-width green rounded rectangle announcing the result. This is marketing copy, not a technical figure. The reader can see the result from the figure itself.

**Good:** End the figure after the corrected row with a brief italic annotation in `muted` color (e.g., "poly(A) was removed during pre-trimming; CPA position recorded in corrected TSV").

### 2. "The Problem" / "RECTIFY's Solution" Framing

**Bad:** Section headers that frame the algorithm as a product pitch.

**Good:** Descriptive headers: "ALIGNED READ (pre-trimmed, poly(A) already removed)" / "AFTER WALK-BACK". Tell the reader what they're looking at, not how to feel about it.

### 3. Warning / Note Boxes

**Bad:** Yellow rounded rectangles with strand-handling notes embedded in the figure.

**Good:** Strand details and edge cases belong in the README text, not in a schematic. Keep the figure focused on one concept.

### 4. Redundant Tracks

**Bad:** Repeating the genome reference track in every panel of a multi-panel figure.

**Good:** Show the genome once at the top. If subsequent panels need it for alignment context, show it, but only if the genomic context is different.

### 5. Unlabeled Color Blocks

**Bad:** Pastel rectangles where you need to match aligner name (in the row label) to block color (elsewhere in the figure) to understand what you're looking at.

**Good:** Put labels inside or immediately adjacent to every block. The reader should never need to build a mental color legend.

### 6. Overly Long Sequences

**Bad:** 17 identical A boxes in the corrected row to show that "the poly(A) tail length is preserved."

**Good:** If the concept is about a *position* (the CPA site), show the position. If you need to show the poly(A) tail, 5-6 A's followed by "..." or a labeled block is sufficient.

---

## Rendering Pipeline

```python
import cairosvg

def save(name, svg_str, h):
    svg_path = os.path.join(OUTDIR, f"{name}.svg")
    png_path = os.path.join(OUTDIR, f"{name}.png")
    with open(svg_path, "w") as f:
        f.write(svg_str)
    cairosvg.svg2png(
        url=svg_path,
        write_to=png_path,
        output_width=FIG_W * 2,    # 2x for retina
        output_height=h * 2,
    )
```

- Always render PNG at 2x resolution (1520px wide from 760px SVG)
- Save both SVG and PNG; the README references the PNG
- The SVG is the source of truth for edits

---

## Checklist for New Figures

Before committing a figure, verify:

- [ ] Width is exactly 760px
- [ ] Uses only colors from `PAL` dict — no ad hoc hex values
- [ ] Font family is the `FONT` stack throughout
- [ ] No text below 8.5px
- [ ] Title is centered at y=22, 14px bold
- [ ] Section headers are uppercase with letter-spacing
- [ ] Row labels are left-aligned at consistent x position
- [ ] No success banners, warning boxes, or marketing framing
- [ ] No overlapping text (check all annotation y-offsets)
- [ ] Nucleotide sequences use minimum bases needed to show the concept
- [ ] All colored blocks have labels (inside or adjacent)
- [ ] Renders legibly at 680px width on GitHub
- [ ] 2x PNG output for retina displays
