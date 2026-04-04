# Coordinate System

RECTIFY uses **0-based, half-open coordinates** throughout — the same convention as pysam, BED files, and the Python `range()` function.

---

## Definitions

| Term | Definition |
|------|------------|
| `reference_start` | 0-based leftmost aligned position (inclusive) |
| `reference_end` | 0-based rightmost position + 1 (exclusive) |
| `5' end` | Transcription start site (TSS) end of the RNA fragment |
| `3' end` | Cleavage/polyadenylation site (CPA) end of the RNA fragment |

---

## Strand-specific coordinate mapping

| Strand | 5' end (TSS) | 3' end (CPA) |
|--------|--------------|--------------|
| **+** (forward) | `reference_start` | `reference_end - 1` |
| **-** (reverse) | `reference_end - 1` | `reference_start` |

---

## Visual diagrams

### Plus strand (+)

RNA is transcribed left → right.

```
5' ────────────────────────────────────────────── 3'  (RNA direction)

Genomic:  100       110       120       130       140
          |---------|---------|---------|---------|
Read:            [===================]
          start=105                  end=135 (exclusive)

5' end = reference_start = 105   (leftmost = TSS for plus strand)
3' end = reference_end - 1 = 134 (rightmost = CPA for plus strand)
```

### Minus strand (-)

RNA is transcribed right → left.

```
3' ────────────────────────────────────────────── 5'  (RNA direction)

Genomic:  100       110       120       130       140
          |---------|---------|---------|---------|
Read:            [===================]
          start=105                  end=135 (exclusive)

5' end = reference_end - 1 = 134  (rightmost = TSS for minus strand)
3' end = reference_start = 105    (leftmost = CPA for minus strand)
```

**Key insight:** For minus strand genes, 5' end is at **higher** genomic coordinates and 3' end is at **lower** coordinates. This is opposite to plus strand.

---

## Python implementation (pysam)

```python
def get_read_5prime_position(read: pysam.AlignedSegment, strand: str = None) -> int:
    """Get 5' end (TSS) position."""
    if strand is None:
        strand = '-' if read.is_reverse else '+'
    if strand == '+':
        return read.reference_start      # Leftmost = 5' for plus
    else:
        return read.reference_end - 1    # Rightmost = 5' for minus

def get_read_3prime_position(read: pysam.AlignedSegment, strand: str = None) -> int:
    """Get 3' end (CPA) position."""
    if strand is None:
        strand = '-' if read.is_reverse else '+'
    if strand == '+':
        return read.reference_end - 1    # Rightmost = 3' for plus
    else:
        return read.reference_start      # Leftmost = 3' for minus
```

---

## Validation examples (*S. cerevisiae*)

| Gene | Strand | Genomic coords (0-based) | Expected 3' ends near | Expected 5' ends near |
|------|--------|--------------------------|----------------------|----------------------|
| CDC19 (YAL038W) | + | 71,787–73,289 | 73,289 (right side) | 71,787 (left side) |
| TFC3 (YAL001C) | - | 151,097–151,166 | 151,097 (left side) | 151,166 (right side) |

---

## GFF/GTF coordinate conversion

GFF and GTF files use **1-based, closed** coordinates (both start and end are inclusive).

```python
# Convert GFF 1-based to 0-based half-open:
pos_0based = gff_pos - 1        # For a single position
start_0based = gff_start - 1    # For interval start
end_0based = gff_end            # For interval end (no change — GFF end is inclusive = 0-based exclusive)
```

RECTIFY converts GFF coordinates to 0-based on load. All internal coordinates are 0-based.

---

## Output files

All positions in RECTIFY output files are **0-based**:

- `corrected_3ends.tsv`: `original_position`, `corrected_position` — 0-based
- `cpa_clusters.tsv`: `start`, `end` — 0-based half-open
- BigWig/BedGraph: 0-based half-open (standard for these formats)

**Exception:** When visualizing in IGV, coordinates are displayed as 1-based. IGV performs this conversion automatically.
