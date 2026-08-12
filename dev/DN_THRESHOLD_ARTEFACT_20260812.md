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

Real spliceosomal introns are essentially all GT-AG; the false population is **81.5 %
non-canonical**. So a motif requirement separates them where length cannot:

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

## Method note

Four hypotheses were proposed and tested overnight; three were refuted by measurement
(coordinate jitter, intron length as originally framed, and low-quality noise — the last
*inverted*: the excess is higher-quality than the baseline). The one that survived was grounded in
a constant read out of the source rather than in a plausible story, which is why it was
falsifiable in a single histogram.
