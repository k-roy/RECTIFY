# RECTIFY Validation Reads — Reference

`validation_reads.bam` contains 32 real reads from `wt_by4742_rep1.bam`
(S. cerevisiae wild-type, direct RNA nanopore) covering eight correction
categories with 4 reads each (2 plus-strand, 2 minus-strand). Each read carries:

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
| `cat1_plus_1` | chrVIII:508599–508794 | + | −4 bp |
| `cat1_plus_2` | chrIX:22544–22850 | + | −3 bp |
| `cat1_minus_1` | chrII:9855–10533 | − | +6 bp |
| `cat1_minus_2` | chrII:9813–10539 | − | +4 bp |

---

## Category 2 — soft-clip rescue (`cat2_softclip`)

Reads where the aligner ended inside a genomic homopolymer (under-called
by the basecaller) and left matching downstream bases as a 3' soft-clip.
`rescue_softclip_at_homopolymer` detects the homopolymer at the alignment
boundary, skips over the uncalled homopolymer extension in the reference,
and matches the soft-clipped bases to the downstream reference sequence.
The corrected position shifts **outward** (corrected > original for +,
corrected < original for −).

| Label | Coords (0-based half-open) | Strand | Shift | Homopolymer base | D op | M op |
|---|---|---|---|---|---|---|
| `cat2_plus_1`  | chrI:68707–69833   | + | +12 bp | T | 10 | GC (2 bp) |
| `cat2_plus_2`  | chrI:198764–199921 | + | +9 bp  | T |  8 | G (1 bp)    |
| `cat2_minus_1` | chrI:65622–66228   | − | −12 bp | T |  8 | ATAA (4 bp) |
| `cat2_minus_2` | chrI:128113–129052 | − | −11 bp | T | 10 | A (1 bp)    |

Note: "Shift" is outward (away from the gene body). For + strand reads
`corrected_3prime = original_3prime + shift`; for − strand reads
`corrected_3prime = original_3prime − |shift|`.

**D op** = number of under-called reference homopolymer bases (reference deletion).
**M op** = non-poly-A/T bases from the soft-clip matched to reference.

CIGAR surgery (v0.9.1): `extend_read_3prime_for_softclip_rescue()` converts the
3' soft-clip to `{D}D {M}M {poly-A}S|H`, making the true RNA end visible in IGV.
The bundled corrected BAMs are `rectified_pA_hardclip.bam` and
`rectified_pA_softclip.bam` (renamed from `rectified.bam` / `rectified_softclip.bam`
in v0.9.1 to convey poly-A tail handling).

Rescue guard (v0.9.2): The soft-clip rescue guard now uses a consecutive-run
threshold (≥ `min_homopolymer_len` = 3) rather than a single-base stop.  Isolated
A's (plus strand) or T's (minus strand) in the soft-clip that match the reference
are rescued; only a run of ≥3 consecutive A/T triggers rollback and stop.  This
fixed cat2_minus_1 (rescues 'ATAA' not 'AA').

Terminal-A strip (v0.9.3): After rescue, any outermost rescued base that is
genomic A (plus strand) or genomic T (minus strand) is stripped, since the
corrected 3′ end must not land on a base that appears as A in the RNA.  This
reduced cat2_plus_2 from +10 ('GA') to +9 ('G'), and cat2_minus_2 from −12
('TA') to −11 ('A').

---

## Category 3 — 5' junction rescue (`cat3_junction`)

Reads where 4/5 aligners soft-clip the 5' exon near an annotated 3'SS.
mapPacBio may span the intron or map upstream without soft-clipping; the other
four aligners (minimap2, gapmm2, deSALT, uLTRA) produce a 5' soft-clip.
The 5' end is rescued to the intron boundary using the affine-gap semi-global
aligner (`local_aligner.py`). Source alignment: minimap2 (soft-clip present).
Requires `--annotation`.

| Label | Coords (0-based half-open) | Strand | 5' soft-clip | Raw 5' end |
|---|---|---|---|---|
| `cat3_plus_1`  | chrII:168808–169462 | + | 10S | 168808 (reference_start) |
| `cat3_plus_2`  | chrI:142618–143383  | + | 22S | 142618 (reference_start) |
| `cat3_minus_1` | chrXV:900071–900767 | − | —   | 900766 (reference_end−1) |
| `cat3_minus_2` | chrII:365845–366503 | − | 25S | 366502 (reference_end−1) |

For plus-strand reads the raw 5' end equals `reference_start`; for minus-strand
reads it equals `reference_end − 1`. The test asserts `five_prime_position ≠ raw_5prime`.

---

## Category 4 — false N op near 3' end (`cat4_false_junc`)

Reads where the aligner inserts a spurious N (intron) CIGAR op in the
poly-T / poly-A region close to the 3' end. The N is absorbed and the
corrected position walks back past it (plus: corrected < original;
minus: corrected > original).

| Label | Coords (0-based half-open) | Strand | Spurious N op | dist from 3' end | Corrected 3' end | Shift |
|---|---|---|---|---|---|---|
| `cat4_plus_1` | chrXI:19592–22073 | + | 20527–22047 (1520 bp) | 26 bp | 22070 | −3 bp |
| `cat4_plus_2` | chrX:392246–393837 | + | 393725–393825 (100 bp) | 12 bp | 393721 | −115 bp |
| `cat4_minus_1` | chrI:128094–129063 | − | 128521–129021 (500 bp) | 427 bp | 128098 | +4 bp |
| `cat4_minus_2` | chrIX:76016–77313 | − | 76027–76250 (223 bp) | 11 bp | 76027 | +11 bp |

**Notes on Cat4 behavior:**
- `cat4_plus_1`: The post-N exon (22047–22073) is `AAAAAAAAAAAAAAAAAAAAAAAATC`. polya_walkback + NET-seq places the corrected 3' end at 22070 (inside the poly-A run). The N-snap does not fire because the walkback does not enter the N region. `n_junctions=1` (N retained as a reported junction).
- `cat4_plus_2`: polya_walkback lands at 393721 (before N at 393725). Ambiguity window [393721,393836] is clipped to [393721,393724] so NET-seq cannot place signal inside the N. corrected=393721. `n_junctions=1`.
- `cat4_minus_1`: N is far from the 3' end (outside FJF window). Treated as a real junction. Normal polya_walkback + NET-seq → 128098. `n_junctions=1`.
- `cat4_minus_2`: polya_walkback lands inside the N (76027–76250). Snap fires → corrected=76027 (junction_start). All leading poly-T M bases consumed. N remains at the left CIGAR boundary. `n_junctions=1`.

---

## Category 5 — multi-aligner chimeric reconstruction (`cat5_chimeric`)

Reads where **no single aligner** reconstructs the full read correctly.
The read spans a two-intron gene; each aligner only finds one of the two GT-AG
introns. Chimeric consensus combines the per-aligner corrected outputs to
produce an alignment with both introns.

All four carry `XK=1`. `XA` lists the two aligners that each contribute a
unique intron. The source alignment (stored in the BAM) comes from one of
those aligners and contains exactly one N op; the other aligner's intron is
visible only after chimeric merging.

| Label | Coords (0-based half-open) | Strand | Source aligner | Source intron (bp) | Complementary aligner | Complementary intron (bp) | Read name prefix |
|---|---|---|---|---|---|---|---|
| `cat5_plus_1`  | chrV:423577–424442   | + | mapPacBio | 423590–423951 (361) | minimap2  | 424421–425030 (609) | `040195ff` |
| `cat5_plus_2`  | chrII:332822–334047  | + | minimap2  | 332875–333386 (511) | mapPacBio | 334050–334122  (72) | `4d1e5c19` |
| `cat5_minus_1` | chrVII:437400–438403 | − | gapmm2    | 437941–438397 (456) | minimap2  | 436480–437396 (916) | `8f86cb34` |
| `cat5_minus_2` | chrIII:177359–178238 | − | mapPacBio | 177906–178213 (307) | gapmm2    | 176709–177362 (653) | `02165816` |

All introns are GT-AG. Verified from `S288C_reference_sequence_R64-5-1_20240529.fsa`.

**Per-read structure:**
- `cat5_plus_1`: The read covers three exons + two introns on chrV (plus strand).
  mapPacBio spans the first exon + intron 1 (423590–423951) but soft-clips at the
  second exon. minimap2 starts at the second exon and catches intron 2
  (424421–425030). Source BAM: mapPacBio chunk_0.
- `cat5_plus_2`: Three exons + two introns on chrII (plus strand). minimap2 spans
  exon1 + intron 1 (332875–333386). mapPacBio starts at exon2 and catches intron 2
  (334050–334122, 72 bp). Source BAM: minimap2 chunk_0.
- `cat5_minus_1`: Three exons + two introns on chrVII (minus strand). gapmm2 spans
  the 5'-proximal intron (437941–438397) near reference_end. minimap2 spans the
  3'-proximal intron (436480–437396) near reference_start. Source BAM: gapmm2 chunk_0.
- `cat5_minus_2`: Three exons + two introns on chrIII (minus strand). mapPacBio
  catches the 5'-proximal intron (177906–178213). gapmm2 catches the 3'-proximal
  intron (176709–177362). minimap2 finds neither. Source BAM: mapPacBio chunk_0.

**Candidate selection method:** `identify_cat5_candidates()` in
`corrected_consensus.py` was run on per-aligner corrected TSVs from chunk_0 of
`wt_by4742_rep1`. Candidates with unique junctions >50 bp apart were inspected;
GT-AG was verified from genome sequence; source reads were confirmed to contain
the expected N op in their CIGAR.

---

## Category 6 — simple chimeric / single-aligner intron rescue (`cat6_chimeric`)

Reads where **one aligner** (mapPacBio) correctly spans the 5' intron while
minimap2 and gapmm2 soft-clip the same region. The reads are stored in
`validation_reads.bam` using the mapPacBio alignment (`XU=1`), making the
junction directly visible. These contrast with Cat3 (all aligners miss the
intron) and Cat5 (two or more aligners each contribute different information).

For plus-strand reads, the 5' intron is near `reference_start` (left side of
CIGAR); for minus-strand reads, it is near `reference_end` (right side of CIGAR,
corresponding to the RNA 5' end).

| Label | Coords | Strand | Intron (mapPacBio) | mm2 soft-clip |
|---|---|---|---|---|
| `cat6_plus_1`  | chrII:125140–126182  | + | 125154–125270 (116 bp) | 14S |
| `cat6_plus_2`  | chrII:168401–169493  | + | 168424–168808 (384 bp) | 20S |
| `cat6_minus_1` | chrII:59701–60697    | − | 60193–60697 (504 bp)   | 9S  |
| `cat6_minus_2` | chrIV:306811–307796  | − | 307333–307765 (432 bp) | 31S |

**Source:** mapPacBio BAM from `dev_runs/wt_by4742_rep1_chunked_20260412/`.
Each read carries `XG=cat6_chimeric`, `XU=1`.

---

## Category 7 — Non-canonical, unannotated splice junctions (`cat7_alt_splice`)

Reads containing a single, biologically plausible splice junction that is
(1) non-canonical (not GT-AG, AT-AC, or GC-AG), (2) not within 50 bp of any
annotated intron boundary, and (3) aligned by mapPacBio (annotation-agnostic).
These test that RECTIFY preserves genuine alternative splice junctions during
3' end correction without attempting 5' rescue.

Motif notation: `{5'SS dinucleotide}-{3'SS dinucleotide}` on the RNA strand.
All reads carry `XG=cat7_alt_splice`, `XU=1` (single aligner: mapPacBio).

| Label | Coords | Strand | Gene region | Junction (intron) | Motif | Read ID prefix |
|---|---|---|---|---|---|---|
| `cat7_plus_1` | chrIII:137716–139164 | + | YCR012W | 138856–138946 (90 bp) | **GT-AT** | `4e43165e` |
| `cat7_plus_2` | chrXII:593396–599029 | + | — | 595736–595852 (116 bp) | **GC-AT** | `0f021462` |
| `cat7_minus_1` | chrII:442949–444728 | − | YBR101C | 443720–443833 (113 bp) | **GT-CG** | `c79f1fb9` |
| `cat7_minus_2` | chrVII:882287–883838 | − | YGR192C | 882352–882702 (350 bp) | **GC-AT** | `5c59f0bc` |

**Biological rationale for motif choices:**
- **GT-AT**: canonical GT donor (standard U2 5'SS), non-canonical AT acceptor
- **GC-AT**: canonical GC donor (known U2 alternative), non-canonical AT acceptor
- **GT-CG**: canonical GT donor, CG acceptor (documented in rare yeast transcripts)

**Source BAM:** `wt_by4742_rep1.mapPacBio.bam` from
`dev_runs/wt_by4742_rep1_chunked_20260412/merged/`. mapPacBio was chosen
over minimap2 because minimap2 applies GT-AG splice bonuses that bias
junction detection toward canonical motifs.

**Expected correction output:**
- `n_junctions = 1` for all four reads
- `five_prime_rescued = 0` (cat7 reads are not 5' truncated at a splice site)
- `correction_applied` includes `atract_ambiguity` and/or `netseq_refinement`
  (normal 3' end processing), but not `five_prime_rescued`

---

## Category 8 — NET-seq A-tract refinement (`cat8_netseq_refine`)

Formerly Cat6. Reads ending in genomic A-tracts at positions with NET-seq signal
are assigned fractional weights across nearby A-tract peaks. Single-peak reads
produce exactly 1 output row with `fraction=1.0`; multi-peak reads produce ≥2
rows with fractions in (0, 1) summing to 1.0.

| Label | Coords | Strand | Expected output |
|---|---|---|---|
| `cat8_plus_single` | chrII:297998–300059 | + | 1 row, fraction=1.0 |
| `cat8_plus_multi` | chrIV:232523–234067 | + | primary corrected_3prime at non-A base (chrIV:234059=G); fractions sum=1.0 |
| `cat8_minus_single` | chrIV:1169625–1172005 | − | 1 row, fraction=1.0 |
| `cat8_minus_multi` | chrVIII:100505–101004 | − | primary corrected_3prime at non-T base (chrVIII:100520=A ref); fractions sum=1.0 |

---

## Regenerating the BAM

All reads are sourced from `wt_by4742_rep1.bam`. If the validation BAM needs
to be rebuilt, use the read names stored in the XV tags as query names to
`samtools view` the source BAM, then re-apply the XV/XG/XK/XA/XS tags with
a pysam script. The selection criteria for each category are documented in
the session transcript.
