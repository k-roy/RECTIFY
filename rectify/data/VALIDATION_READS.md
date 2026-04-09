# RECTIFY Validation Reads — Reference

`validation_reads.bam` contains 24 real reads from `wt_by4742_rep1.bam`
(S. cerevisiae wild-type, direct RNA nanopore) covering all six correction
categories. Each read carries:

- **XV**: label (e.g. `cat1_plus_1`) — used as primary key in tests
- **XG**: category name (e.g. `cat1_indel`)
- **XK**: chimeric flag (`1` = chimeric, Cat5 only)
- **XA**: comma-separated aligner list (Cat5 only)
- **XS**: segment count (Cat5 only)

---

## Category 1 — polya_walkback (`cat1_indel`)

Reads end in genomic A-tracts; the aligner's 3' boundary is shifted inward
(plus: left; minus: right) by polya_walkback. All four reads confirmed to
produce `corrected_3prime != original_3prime` after `rectify correct`.

| Label | Coords (0-based half-open) | Strand | Shift |
|---|---|---|---|
| `cat1_plus_1` | chrVIII:508599–508794 | + | −3 bp |
| `cat1_plus_2` | chrIX:22544–22850 | + | −3 bp |
| `cat1_minus_1` | chrII:9855–10533 | − | +6 bp |
| `cat1_minus_2` | chrII:9813–10539 | − | +4 bp |

---

## Category 2 — soft-clip rescue (`cat2_softclip`)

Reads with A-rich 3' soft-clips that match the reference. The soft-clip is
rescued and the corrected position shifts by 1–20 bp relative to the raw
alignment end.

| Label | Coords | Strand | Shift |
|---|---|---|---|
| `cat2_plus_1` | chrI:31119–31567 | + | −6 bp |
| `cat2_plus_2` | chrI:11324–14730 | + | −1 bp |
| `cat2_minus_1` | chrI:32932–34333 | − | +1 bp |
| `cat2_minus_2` | chrI:31417–32958 | − | +1 bp |

---

## Category 3 — 5' junction rescue (`cat3_junction`)

Reads with a 5' soft-clip that lands within `junction_proximity_bp` (default
10 bp) of an annotated splice-site 3'SS. The 5' end is rescued to the intron
boundary. Requires `--annotation` to be passed to `rectify correct`.

| Label | Coords | Strand | Gene / Intron | Rescued 5' end |
|---|---|---|---|---|
| `cat3_plus_1` | chrII:125270–126182 | + | intron 125154–125270 | 125153 |
| `cat3_plus_2` | chrII:168808–169478 | + | intron 168424–168808 | 168423 |
| `cat3_minus_1` | chrXV:900071–900767 | − | RPL20B intron 900767–901193 | 901193 |
| `cat3_minus_2` | chrII:59715–60193 | − | intron 60193–60697 | 60697 |

For plus-strand reads the raw 5' end equals `reference_start`; for minus-strand
reads it equals `reference_end − 1`. The test asserts `five_prime_position ≠
raw_5prime`.

---

## Category 4 — false N op near 3' end (`cat4_false_junc`)

Reads where the aligner inserts a spurious N (intron) CIGAR op in the
poly-T / poly-A region close to the 3' end. The N is absorbed and the
corrected position walks back past it (plus: corrected < original;
minus: corrected > original).

| Label | Coords (0-based half-open) | Strand | Spurious N op | dist from 3' end | Shift |
|---|---|---|---|---|---|
| `cat4_plus_1` | chrXI:19592–22073 | + | 20527–22047 (1520 bp) | 26 bp | −45 bp |
| `cat4_plus_2` | chrX:392246–393837 | + | 393725–393825 (100 bp) | 12 bp | −72 bp |
| `cat4_minus_1` | chrI:128094–129063 | − | 128521–129021 (500 bp) | 427 bp | +3 bp |
| `cat4_minus_2` | chrIX:76016–77313 | − | 76027–76250 (223 bp) | 11 bp | +11 bp |

---

## Category 5 — chimeric reconstruction (`cat5_chimeric`)

Chimeric reads whose consensus was reconstructed from multiple aligner
outputs. No re-correction is applied; tests verify BAM tags only.

| Label | Coords (0-based half-open) | Strand | Aligners (XA) | Segments (XS) | Gene |
|---|---|---|---|---|---|
| `cat5_plus_3aligner` | chrI:9411–9965 | + | minimap2,mapPacBio,gapmm2 | 9 | YAL066W (partial, 159 bp overlap) |
| `cat5_plus_2aligner` | chrI:7109–9030 | + | minimap2,mapPacBio | 32 | intergenic (YAL067C is −strand) |
| `cat5_minus_long` | chrI:4923–6316 | − | minimap2,mapPacBio | 28 | intergenic |
| `cat5_minus_short` | chrI:5044–5683 | − | minimap2,mapPacBio | 18 | intergenic |

All four carry `XK=1`.

---

## Category 6 — NET-seq A-tract refinement (`cat6_netseq_refine`)

Reads ending in genomic A-tracts at positions with NET-seq signal are
assigned fractional weights across nearby A-tract peaks. Single-peak reads
produce exactly 1 output row with `fraction=1.0`; multi-peak reads produce
≥2 rows with fractions in (0, 1) summing to 1.0.

| Label | Coords | Strand | Expected output |
|---|---|---|---|
| `cat6_plus_single` | chrII:297998–300059 | + | 1 row, fraction=1.0 |
| `cat6_plus_multi` | chrIV:232523–234067 | + | ≥2 rows, fractions sum=1.0 |
| `cat6_minus_single` | chrIV:1169625–1172005 | − | 1 row, fraction=1.0 |
| `cat6_minus_multi` | chrVIII:100505–101004 | − | ≥2 rows, fractions sum=1.0 |

---

## Regenerating the BAM

All reads are sourced from `wt_by4742_rep1.bam`. If the validation BAM needs
to be rebuilt, use the read names stored in the XV tags as query names to
`samtools view` the source BAM, then re-apply the XV/XG/XK/XA/XS tags with
a pysam script. The selection criteria for each category are documented in
the session transcript.
