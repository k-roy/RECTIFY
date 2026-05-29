# RECTIFY Validation Reads — upf1Δ parallel set

`validation_reads_upf1d.bam` contains **32 real reads from upf1Δ rep1**
(S. cerevisiae *upf1Δ*, direct-RNA nanopore) covering **8 correction categories
with 4 reads each (2 plus-strand, 2 minus-strand)**. It is a **parallel,
purely-additive** companion to the wild-type `validation_reads.bam` (mirrors
`validation_reads_cdna.bam` / `_quantseq_rev.bam`); the wt 36-read set and
`tests/test_validation_reads.py` are untouched. Tests: `tests/test_validation_reads_upf1d.py`
(77 passing). Each read carries:

- **XV**: label (e.g. `cat1_plus_1`) — primary key in tests
- **XG**: category name (e.g. `cat1_indel`)
- **Xg**: segment count (`#N-ops + 1`; cat5/cat6 ≥ 2) — rectify-internal namespace
- **Xm**: aligner-consensus count (≥ 1) — rectify-internal namespace

---

## Provenance

```
drs_trim_upf1d/upf1d_rep1.bam   (Dorado, raw, 742,308 reads — immutable source)
  └─ rectify trim-polya ─────►  upf1d_rep1_trimmed.fastq.gz (+ trim metadata)
       └─ subsample ─────────►  upf1d_rep1_sub100k.fastq.gz (98,848 reads)
            └─ rectify align ─► 5-aligner panel  (minimap2, gapmm2, mapPacBio,
                                deSALT, uLTRA)  — $SCRATCH/upf1d_fullpanel_20260527/align/
```

`validation_reads_upf1d_dorado_source.bam` is the immutable raw-Dorado archive
(reads exactly as basecalled, XV/XG tags added). The aligner BAMs are built from
**poly-A-trimmed** reads. Ground truth was generated on M1 (authoritative code:
rescue DP cap + cat3 early-exit fix) via the correct-first pipeline
(`rectify correct` per aligner → `merge_corrected_tsvs`).

Rebuild: `python scripts/validation_data/upf1d_2026_05/build_upf1d_validation.py`
(manifest at `scripts/validation_data/upf1d_2026_05/manifest.tsv`).

### Category 4 deliberately omitted
The poly-A-trimmed panel produces **zero** `false_junction` (cat4) artifacts
across all 96 k reads / 5 aligners: the cat4 artifact (aligner maps the poly-A
**tail** to a distant genomic A-tract via an N-op) requires the untrimmed tail,
which the conservative trimmer removes. cat4 is therefore **not a scenario the
real trimmed pipeline encounters**; it remains covered by the wt set and will be
added here only if genuine leak-through reads appear in another dataset.

### Note on mapPacBio multi-primary
mapPacBio emits multiple records all flagged **PRIMARY** for multi-mapping reads
(see `docs/aligners/mapPacBio.md`). The builder disambiguates by the manifest
`chrom` column (preferring the N-op-bearing record). Multi-mapping reads were
avoided for the structural categories.

---

## Category 1 — A-tract indel / polya_walkback (`cat1_indel`)
Reads with a small indel near the 3′ end over a genomic A-tract; corrected via
`indel_correction` (± `polya_walkback`). Test: `indel_correction` applied +
exact `corrected_3prime`.

| Label | Coords | Strand | orig→corr 3′ |
|---|---|---|---|
| cat1_plus_1 | chrXIV | + | 708337→708336 |
| cat1_plus_2 | chrIV | + | 1239967→1239967 |
| cat1_minus_1 | chrIV | − | 1162817→1162861 |
| cat1_minus_2 | chrVIII | − | 381917→381917 |

## Category 2 — soft-clip rescue at an under-called homopolymer (`cat2_softclip`)
Aligner ended inside a genomic homopolymer leaving matching downstream bases as a
3′ soft-clip; `softclip_rescue` extends to the true end. Test: `softclip_rescue`
applied + exact `corrected_3prime`.

| Label | Coords | Strand | corr 3′ |
|---|---|---|---|
| cat2_plus_1 | chrII | + | 723682 |
| cat2_plus_2 | chrIV | + | 66060 |
| cat2_minus_1 | chrIV | − | 233917 |
| cat2_minus_2 | chrIV | − | 236091 |

## Category 3 — 5′ junction rescue (`cat3_junction`)
4/5 aligners 5′-soft-clip near an annotated 3′SS; the 5′ end is rescued across
the intron (exercises the rescue DP cap + early-exit fix). Test: `five_prime_rescued`
applied + `five_prime_position` moved from raw + exact value.

| Label | Coords | Strand | rescued 5′ | intron (GT-AG, annotated) |
|---|---|---|---|---|
| cat3_plus_1 | chrXVI | + | 282139 | 282140–282665 |
| cat3_plus_2 | chrV | + | 148193 | 148194–148282 |
| cat3_minus_1 | chrIV | − | 307765 | 307333–307765 |
| cat3_minus_2 | chrXII | − | 713155 | 712725–713155 |

## Category 5 — multi-intron reads (`cat5_multi_intron`)
Reads spanning a 2-intron transcript (≥2 annotated introns, Xg ≥ 3). The strict
"chimeric, no single aligner reconstructs" variant is **rare** in upf1Δ (1 found
genome-wide), so this set uses genuine 2-intron transcripts as the multi-segment
exemplars. Test: structural (present, Xg ≥ 2, stored alignment has an N-op).

| Label | Coords | Strand | introns |
|---|---|---|---|
| cat5_plus_1 | chrII | + | 462209–462289, 462429–462499 |
| cat5_plus_2 | chrIII | + | 107033–107110, 107191–107287 |
| cat5_minus_1 | chrIV | − | 1319617–1319697, 1319720–1319816 |
| cat5_minus_2 | chrVII | − | 364964–365432, 365526–365985 |

## Category 6 — single-aligner intron span (`cat6_chimeric`)
mapPacBio spans the intron with no 5′ soft-clip where minimap2/gapmm2 5′-clip.
Test: Xm ≥ 1, stored mapPacBio has the N-op + no 5′ clip, corrected output
contains the annotated intron, n_junctions = 1.

| Label | Coords | Strand | intron (GT-AG) |
|---|---|---|---|
| cat6_plus_1 | chrII | + | 45644–45977 |
| cat6_plus_2 | chrII | + | 45644–45977 |
| cat6_minus_1 | chrII | − | 89132–89440 |
| cat6_minus_2 | chrII | − | 592416–592768 |

## Category 7 — non-canonical, unannotated splice junctions (`cat7_alt_splice`)
Single biologically-plausible non-canonical junction (not GT-AG/AT-AC/GC-AG),
> 50 bp from any annotated boundary. Test: n_junctions ≥ 1, junction coords,
`five_prime_rescued` NOT applied.

| Label | Coords | Strand | junction | motif (RNA) |
|---|---|---|---|---|
| cat7_plus_1 | chrII | + | 246695–246734 | TT-TC |
| cat7_plus_2 | chrII | + | 248644–248687 | TA-CA |
| cat7_minus_1 | chrIV | − | 1489311–1489390 | GT-TA |
| cat7_minus_2 | chrXII | − | 476768–476896 | TC-GC |

## Category 8 — A-tract 3′ end (`cat8_atract`)
Reads ending in a genomic A-tract; corrected 3′ lands on the first non-A base
(`atract_ambiguity`). NET-seq refinement is opt-in, so each read yields a single
output row (fraction = 1.0). Test: 1 row, fraction = 1.0, `atract_ambiguity` applied.

| Label | Coords | Strand | corr 3′ (non-A base) |
|---|---|---|---|
| cat8_plus_1 | chrII | + | 290447 |
| cat8_plus_2 | chrII | + | 24929 |
| cat8_minus_1 | chrII | − | 601356 |
| cat8_minus_2 | chrI | − | 128117 |

## Category 9 — Module-2H N-op boundary refinement (`cat9_junction_refine`)
With `--aligner-bams` + `--annotation`, Module 2H re-scores an imprecise N-op
boundary and snaps it to the annotated junction. Test: corrected (Module-2H)
junction == annotated coords.

| Label | Coords | Strand | refined intron (annotated) |
|---|---|---|---|
| cat9_plus_1 | chrVII | + | 439093–439323 |
| cat9_plus_2 | chrII | + | 110879–110948 |
| cat9_minus_1 | chrII | − | 414753–415259 |
| cat9_minus_2 | chrII | − | 407027–407122 |
