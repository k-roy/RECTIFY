# The `intronlen=40` D→N threshold manufactures novel junctions — and raising it is the WRONG fix

**Found 2026-08-12 overnight, jointly by `668-drs-arm` (measurement) and `rectify-realigner`
(mechanism + annotation analysis). Confirmed decisively. NOT FIXED — the lever is a design
decision for Kevin.**

## The observation

Station C junction counts vary **46×** across DRS samples at matched ~1M depth (741 → 34,122)
while *annotated* junction counts stay flat at 222–270. The entire spread is in NOVEL junctions,
and in the high-count samples **96.2 % of them are singletons** — ~70 % of junction-carrying reads
each produce a junction seen exactly once.

Those singletons are **not** low-quality: median `q_max` 93.3 (vs 28.5 in low-count samples),
12 % below q<50 (vs 62 %), and *less* repeat context (12.6 % vs 42.6 %). They are also **81.5 %
NON-canonical** — and real yeast splicing is essentially all GT-AG.

## The mechanism: they are DELETIONS RELABELLED AS INTRONS

The panel converts a reference gap into an N-op purely on LENGTH. It is explicit in source:

```
multi_aligner.py:1016   "BBMap `intronlen` is the MIN deletion length relabeled D->N"
multi_aligner.py:1018   'intronlen=40'
multi_aligner.py:1177   "Uses intronlen=40 so any reference gap >= 40 bp is encoded as an N-op"
overhang_resolver.py:85 min_intron: int = 40   # matches BBMap intronlen=40 D->N semantics
                        (--min-intronlen default 20 on the STAR/HISAT2 path)
```

An N-op consumes REFERENCE without consuming query, so it cannot absorb errors *in* the read — it
encodes a stretch where the read skips reference. DRS's dominant error mode is deletion relative to
reference. **A read carrying a ~40–100 bp net reference gap gets an N-op by the length rule alone,
however the gap arose.**

### The histogram that confirmed it (`wt_rep1`, 33,030 singleton novel junctions)

```
  len < 20:    1  (0.003 %)
  20-24:       0
  25-29:       0
  30-34:      19
  35-39:     156   ###
  40-44:     852   ####################          <-- intronlen=40
  45-49:    1998   ################################################
  50-54:    2471   ############################################################  PEAK
  55-59:    2314 · 60-64: 2254 · 65-69: 2130 · 70-74: 1970 · ... smooth decay
```

**99.5 % sit at or above 40 bp, with essentially nothing between 20 and 39 and a step change at
exactly the threshold.** A splicing-length distribution has no reason to know about 40.
`rna15_rep1` reproduces the shape (0.45 % below 40, step at 40–44, peak at 50–54).

## 🔴 Why raising `intronlen` is the WRONG fix

The obvious lever is to raise the threshold. **The bundled R64-5-1 annotation says do not.**
378 annotated introns:

| threshold | annotated introns lost | of which spliceosomal (non-tRNA) |
|---|---:|---:|
| < 40 bp | 60 (15.9 %) | **2** |
| < 50 bp | 62 (16.4 %) | ~4 |
| < 60 bp | 67 (17.7 %) | ~9 |
| < 100 bp | **167 (44.2 %)** | **~107** |

**58 of the 60 sub-40 bp introns are tRNA introns** (spliced by the tRNA endonuclease, not the
spliceosome, and largely irrelevant to mRNA DRS/cDNA). So `intronlen=40` is well chosen *as a
threshold* — it costs only 2 of 318 spliceosomal introns.

But the false population **peaks at 50–54 bp with a long tail**, so removing most of it would
require raising the threshold near or above 100 bp — which would suppress **~107 real
spliceosomal introns, about a third of them**. The threshold cannot separate the two populations
because they overlap almost completely in length.

## ✅ The right lever: require a canonical motif for SHORT gaps

⚠️ **READ THE QUALIFICATION BELOW FIRST — an earlier draft of this section overstated the
evidence.** The motif is the right LEVER, but the 81.5 %-non-canonical figure is **not** diagnostic
of the D→N class; it is the baseline for *novel* junctions generally in this data. The cliff at 40
is the evidence for the mechanism. See "What the motif figure does and does not show".

Real spliceosomal introns are essentially all GT-AG. So a motif requirement is a usable separator
where length is not:

> A reference gap below some length may be called an intron (N) only if it carries a canonical
> donor/acceptor; otherwise it stays a deletion (D).

Removes ~81.5 % of the false population at ~zero cost to real introns. Two placements, in
increasing invasiveness:

1. **Reporting-side (cheapest, reversible):** filter in the junction table — a short NOVEL
   junction must be canonical to be censused. Station C already computes `canonical_in_class`, so
   this is nearly free and changes no alignment.
2. **Alignment-side (correct at source, higher blast radius):** apply the motif test at the D/N
   decision. Changes CIGARs for every dataset — needs its own validation pass.

**Neither is implemented.** Nothing landed on 2026-08-12 touches this.

## Consequences to carry regardless

- **`review` counts are NOT comparable across samples in either direction.** The same column holds
  short D→N reclassifications in high-count samples and long spurious junctions in low-count ones
  (`ski7d_rep1` core junctions: median **24,313 bp**, 72.8 % > 1 kb — the separate impossible-
  junction class addressed by `intron_sanity.py`).
- **`support >= 2` is working as intended**, correctly refusing ~33,000 one-off reference gaps.
- **This is library-dependent, not genotype-dependent** — `wt_rep1`, the WT-AA *control*, has the
  most. Whatever drives per-read deletion rate (basecaller version, chemistry, read quality) drives
  this, which is why "deposit" fitted better than "strain background" throughout.
- **The sub-40 residue** (19 in 30–34, 156 in 35–39, zero in 20–29) is consistent with a second
  path at a different threshold — STAR/HISAT2 default to `--min-intronlen 20`. A **per-arm**
  histogram, each arm showing a cliff at its own threshold, would confirm. Not yet run.

## 🔴 THE TWO PHENOMENA HAVE OPPOSITE ARM SIGNATURES — minimap2-only is NOT a refuge

Per-arm, `wt_rep1`, 300k reads/arm (measured by `668-drs-arm`):

| arm | short D→N class (~40–120 bp) | long impossible class (>10 kb) | N-ops < 20 bp |
|---|---:|---:|---:|
| minimap2 | **~34,000** | **0** | 0 |
| uLTRA | ~35,000 | 394 | **132** |
| deSALT | ~35,000 | 2,067 | 0 |
| overhang_resolver | ~34,000 (byte-identical to minimap2) | **0** | 0 |

**The long class is uLTRA/deSALT-specific and minimap2 is clean. The short class is PANEL-WIDE and
minimap2 participates equally.** So "minimap2 is the clean arm" is true only of the long class.

⚠️ **The natural mitigation — switch to minimap2-only to avoid impossible junctions — carries the
short D→N population entirely untouched.** State this wherever the long class is discussed.

`overhang_resolver` is byte-identical to minimap2 here, confirming again that Station A inherits
minimap2 placements and contributes nothing of its own to either class.

### Correction: the 40 boundary is NOT a shared parameter we set

A per-arm histogram was predicted to show each arm cliffing at ITS OWN threshold (40 BBMap-family,
20 STAR/HISAT2). **It does not — every arm cliffs at ~40, including minimap2**, which is not on the
BBMap path. And rectify passes minimap2 **no min-intron parameter at all**:

```
run_minimap2: -ax splice -uf -k14 -G <max_intron> --splice-flank=no --secondary=no --MD -y
```

So this is not a rule we configure on that arm, and not (as first supposed) a shared
post-alignment D→N normalisation — no such step exists in the tree.

**HYPOTHESIS (untested): it is a SCORING CROSSOVER intrinsic to splice-aware alignment.** A
deletion of length L costs `gap_open + L * gap_extend` (affine, grows with L) while an intron costs
a roughly fixed splice penalty. Above a crossover length the intron is simply cheaper, so gaps
flip from D to N there — no threshold required. That crossover lands in the 30–40 bp region under
default long-read parameters, which is also where BBMap's explicit `intronlen=40` rule sits, so the
panel looks like one boundary while arriving at it by two different routes. It also explains why
the transition is a steep RISE rather than a hard zero: a cost crossover is soft, a rule is not.

**TEST: vary minimap2's gap costs (`-O`/`-E`) or the splice penalty and re-histogram.** If the
crossover MOVES, it is scoring; if it stays at ~40 regardless, there is a hard rule somewhere still
to be found. Cheap — one small realignment.

**If the crossover explanation holds, it strengthens the motif fix decisively:** you cannot
parameterise your way out uniformly, because each arm reaches the boundary by different means. A
motif requirement is arm-independent.

The one arm that differs is **uLTRA**, the only one with a sub-20 population (132 N-ops < 20 bp;
16 at 10–14, 8 at 15–19) — consistent with a lower or absent floor on that path.

## ✅ CLOSED: the driver is per-library basecall quality, and it is quantitatively linked

**Strain background REFUTED outright.** Same anchor-away background, different deposit:

| sample | background | deposit | junctions | annotated |
|---|---|---|---:|---:|
| `wtaa_rep1/2/3` | **HHY168-AA** | `wtaa_drs_2025` | **1,266 / 1,323 / 1,366** | 262 |
| `wt_rep1/2/3` | **HHY168-AA** | `PRJNA1229592` | **34,122 / 34,680 / 19,095** | 264–269 |

Identical background, **22× apart** on the means. Genotype is refuted too — `wt_rep1` (the WT-AA
*control*) sits at 34k while `wtaa`, also WT-AA, sits at 1.3k.

**The extreme samples are the `pod5_skip` batch.** `ysh1_rep1` (178,469) and `ysh1_rep2` (196,643)
are the corpus's only two `pod5_skip` libraries — two barcodes of one flow cell, already flagged
INDEPENDENTLY for untailed-read excess, and split the same way by three other QC metrics and by the
pipeline's own PCA/correlation clustering. Their own third replicate sits at 23,915.

### The two QC signals scale together ~1:1 — strong evidence of a single upstream cause

| sample | junctions | × sibling | untailed reads | × sibling | ratio of folds |
|---|---:|---:|---:|---:|---:|
| `ysh1_rep1` | 178,469 | **7.46×** | 14.7 % | **7.74×** | **0.96** |
| `ysh1_rep2` | 196,643 | **8.22×** | 13.9 % | **7.32×** | **1.12** |
| `ysh1_rep3` | 23,915 | 1.00× | 1.9 % | 1.00× | — |

Two independently-derived QC signals — one from poly(A) tailing, one from junction counting —
track each other within ~12 % across an 8× range. That is what a common upstream driver
(basecall quality → per-read deletion rate) looks like, and it is hard to explain as coincidence.

**`annotated` stays 262–271 across all nine samples, spanning a 155× range of totals** — the
internal control that real signal is detected consistently while only the false population moves.

### 🔴 ACTIONABLE: quarantine `ysh1_rep1` and `ysh1_rep2` from junction-level analysis

178k–197k censused with 99k–108k in `review`, against 24k for their own third replicate. **Any
cross-genotype junction comparison including them compares library quality, not biology.** The
corpus now has two independent reasons to exclude them.

### A useful inversion

If junction count tracks basecall quality this tightly, **junction count is itself a cheap
library-quality QC metric** — screenable before a junction table is trusted. The artefact doubles
as a diagnostic.

**Still open:** whether `pod5_skip` per se or a correlated property (chemistry, basecaller version,
flow-cell health) is proximate. Directly testable as per-read deletion rate for these libraries vs
`wtaa`, which would convert "library quality" into a measurable covariate.

## What the motif figure does and does NOT show (correction, `668-drs-arm`)

**The ~80 % non-canonical fraction is NOT specific to the short D→N class.** Measured on the same
data, novel junctions **> 150 bp** are `wt_rep1` **79.5 %** and `ysh1_rep1` **79.3 %**
non-canonical — essentially the same as the short class's 81.5 %.

⇒ **Non-canonical distinguishes NOVEL from ANNOTATED, not the D→N class from anything else.**
81.5 % is therefore *not elevated* relative to other novel junctions, and an earlier framing here
(and in my first relay of it) implied that it was. Corrected:

- **Evidence for the diagnosis = the cliff at 40**, plus the arm-independence and the per-library
  scaling. Not the motif fraction.
- **The motif is a LEVER, not a diagnostic.** It happens to remove ~82 % of the short class, which
  is useful, but it would remove ~80 % of long novel junctions too.
- ⇒ **Restricting the filter to short gaps is a POLICY CHOICE, not something the data forces.** If
  the target is the D→N artefact specifically, keep the length restriction. If the target is
  spurious junctions generally, the same rule hits the long class — **which overlaps
  `intron_sanity.py`'s truncation guard, and the interaction of the two needs thinking about
  before both are enabled.**

## Reporting-side placement has a STRUCTURALLY ZERO cost to known introns

Verified in source: Station C **excludes annotated junctions from its table before writing**

```python
station_c.py:371   if (chrom, s, e) in annotated:
                       n_annotated += 1
                       continue        # counted, then SKIPPED — never reaches the rows
```

so **every row is novel by construction, and a motif filter applied there cannot remove a single
annotated intron.** Not "removes few" — *cannot*. That is a stronger safety property than the
annotation-side analysis above had to assume, and it separates the two placements sharply:

| placement | risk to annotated introns | reversible | changes CIGARs |
|---|---|---|---|
| **1. reporting-side (Station C table)** | **structurally zero** | yes | no |
| 2. alignment-side (D/N decision) | real — acts before annotation is consulted; the ~2-of-318 sub-40 figure applies | no | yes, every dataset |

Also verified: `canonical_in_class` is a genuine motif test — *"True iff ANY member of the
junction's ambiguity class is canonical"* via `is_canonical_junction` over the ambiguity window.
It is **conservative in the safe direction**: a junction counts as canonical if any shifted
position in its window is GT-AG, so the non-canonical fractions quoted here are lower bounds on
how non-canonical the population really is.

### Measured benefit on DRS (novel junctions ≤ 150 bp a canonical requirement would remove)

| sample | novel ≤150 bp | non-canonical | removed |
|---|---:|---:|---:|
| `wt_rep1` | 28,007 | 22,933 | **81.9 %** |
| `ysh1_rep1` | 119,393 | 98,717 | **82.7 %** |
| `wtaa_rep1` (clean library) | 200 | 128 | **64.0 %** |

The filter is **least** effective on the clean library and most on the degraded ones — consistent
with it removing the false population rather than a fixed share of everything.

## Method note

Four hypotheses were proposed and tested overnight; three were refuted by measurement
(coordinate jitter, intron length as originally framed, and low-quality noise — the last
*inverted*: the excess is higher-quality than the baseline). The one that survived was grounded in
a constant read out of the source rather than in a plausible story, which is why it was
falsifiable in a single histogram.
