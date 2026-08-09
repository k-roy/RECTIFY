# C6 — variant-aware junctions: gate plan + VERDICT (2026-06-30)

The discovery-FDR-GUARD facet. A standing germline/somatic variant near a splice
site can FABRICATE a spurious "novel" intron: the haploid-reference aligner
re-expresses a GT..AG-flanked **deletion variant** as an **N-op (intron)** instead
of a **D-op (deletion)**, inflating the de-novo (esp. non-canonical) discovery FDR.
C6's premise: a **variant-aware emission** would score the locus as a deletion and
recover the TRUE placement.

Mirrors the C1–C4 gate discipline: establish by MEASUREMENT against TRUTH that the
incumbent is BELOW CEILING on *recoverable* reads AND the gap is **addressable from
READ EVIDENCE** (not the zero-evidence trap) BEFORE proposing any decoder, with a
pre-committed null, a zero-evidence guard, and a SPECIFICITY fence. Fitness = truth,
NEVER an internal score.

## The incumbent + the stratum (both pinned firsthand)

- Incumbent = `minimap2 -ax splice -uf --eqx` (the shipped scorer/smoke (E) locus
  metric is minimap2-based; the faithful native-junction caller).
- Stratum = `controlled.gen_variant_stratum` (REUSED, untouched):
  - **SPLICE_MIMIC_DEL drivers** — single-exon contig with a `GT..AG`-flanked block
    of length {40,60,100} bp (≥40 = driver), random flanks/interior, het/hom/VAF
    varied. Read = `left+right` (block deleted). **Truth = a DELETION; no intron.**
  - **Controls** — `SNP_NEAR_JUNC` (real intron + SNP 3 bp into exon1) and
    `SNP_DISTANT`; plus short-deletion (<40 bp) blocks that stay a correct D.
- Added (probe-local, append-only, default-OFF re shipped code): a **length-matched
  REAL-intron arm** at the SAME lengths {40,60,100} bp (GT..AG, random interior,
  spliced reads, **truth = the intron**) — the specificity arm. No shipped code
  touched; smoke default byte-identical (re-run GREEN, (E) unchanged).

## The decisive measurements (DECODER-FREE; `scripts/benchmark/c6_headroom.py`)

1. **(a) below-ceiling** — fraction of drivers where mm2 fabricates a variant-adjacent
   FP intron instead of the truth D. ceiling = the D is in-principle placeable (the
   read genuinely lacks the block → a D-op reconstructs it 0-mismatch by construction).
2. **(b) zero-evidence NM test** (the C5/(G)-equivalent) — for each driver FP intron,
   is mm2's N-op alignment **NM==0**? If yes, the intron is an EQUALLY-GOOD (0-mismatch)
   alignment of the SAME query as the truth D → no read-level evidence separates them →
   a motif-blind score cannot prefer truth. (Largely construction-predetermined; a
   confirmation, **not** the load-bearing measure.)
3. **(c) THE TENSION** (a DISCLOSED NULL, *not* an at-risk falsifiable test — panel
   correction) — the ONLY read-level rule that recovers a driver is "a long GT..AG gap
   is a DELETION, not an intron." Applied to the length-matched REAL introns it converts
   `k` real introns to FN, ~1:1. **Be honest about what this is:** the matched-intron arm
   is drawn from the *identical* generator as the drivers (same `L=200` random flanks,
   same `GT+_rand_unique+AG` block at the same lengths), so the two arms are IID and
   `mm2_FP_N == mm2_TP_N` is FORCED by distributional symmetry — the REFUTE branch was
   the only achievable outcome, and the PROCEED branch is structurally unreachable on
   this corpus. That symmetry *is* the zero-evidence proof (two constructs identical
   under the read-generating process cannot be separated by ANY function of the read),
   but it is zero-evidence **BY CONSTRUCTION**, a disclosed null — NOT a measurement that
   could have gone the other way. The at-risk, information-carrying measures are **(a)**
   (the 100% FP rate) and **(d)** (the ~1.0 real-junction recall).
4. **(d) specificity (positive)** — mm2 RETAINS real junctions (recall ~1.0 on
   SNP_NEAR_JUNC and on the matched real introns) → mm2 is NOT blunt; a blunt
   "abstain/relabel near any variant" rule would wrongly suppress real junctions.
   **Scope (panel correction):** (d) refutes a *blunt/read-evidence* rule. It does NOT
   refute a **catalog-gated** rule ("relabel N→D only where a catalogued deletion allele
   overlaps"): the matched real introns carry NO `VariantTruth`, so on this very corpus a
   catalog-gated rule spares 180/180 real introns → specificity 1.0. The catalog-gated
   rule's *recall*, however, is circular here (the catalog == truth on the SIM) — which
   is why it is a DEFER-to-real-data integration, not a corpus PROCEED.

**The flat-arm trap (advisor, integrated).** The C1/C3 flat-affine DP arm has NO
N-op, so it emits D for EVERY long gap — it "recovers" every driver *by being
intron-blind*, and that SAME arm converts every real intron to a giant deletion
(catastrophic specificity failure). **Flat-arm recovery is NOT a C6-addressability
signal; it IS the tension.** Hence minimap2-only, and the recovery rule is scored
against real introns, never read as a member win.

## RESULT — (2026-06-30, reps=60, seed=7; `c6_headroom_result.txt`)

| measure | value |
| --- | --- |
| (a) drivers n / placed | 180 / 180 |
| (a) ceiling (D placeable) | 1.0000 |
| (a) incumbent (mm2 emits truth D) | **0.0000** (0/180) |
| (a) **HEADROOM** (mm2 fabricates variant-adjacent FP intron) | **1.0000** (180/180) |
| (a) scorer `fp_variant_adjacent` / indel concordance | 180 / 0.0000 |
| (b) driver FP introns with **NM==0** | **180/180** |
| (c) tension, per length {40,60,100} | driver FP-N = real-intron TP-N = 60 each (1:1) |
| (d) SNP_NEAR_JUNC real-junction recall | 0.9833 |
| (d) matched real-intron recall | 0.9889 (N on 180/180; NM==0 on 178) |

minimap2 fabricates a variant-adjacent FP intron on **100%** of the ≥40 bp drivers
(the smoke (E) effect at scale — the incumbent IS below ceiling). **But** all 180 FP
introns are **NM==0** — byte-identical 0-mismatch alignments to the truth D — and the
sole read-level recovery rule converts the 180 length-matched real introns to FN
**1:1**. Drivers and real introns are indistinguishable at the read level (same
length, GT..AG, random interior). Specificity (d) confirms mm2 already retains real
junctions, so a blunt near-variant abstain/relabel rule would *suppress* them.

## VERDICT — **DEFER** (active matched-VCF integration milestone); the read-evidence-only member is REFUTED

Headline = **DEFER**, not REFUTE. The DEFER vs REFUTE test: a named real-data trigger +
a well-defined fix + a real (non-phantom) FP population ⇒ DEFER. C6 has all three (180/180
variant-adjacent FP, a catalog-gated N→D emission, a matched-germline-VCF trigger) — the
gap is REAL, merely unvalidatable on this SIM. (Contrast C3, a true REFUTE: the incumbent
was AT CEILING — no gap at all.) The **REFUTE** is the precision, scoped to the
*read-evidence-only* member: no per-read evidence separates a GT..AG deletion from a
GT..AG cryptic intron, so a member that consumes only the read cannot close the gap.

- The below-ceiling gap is the **ZERO-EVIDENCE trap**: a GT..AG deletion-block and a
  GT..AG cryptic intron produce identical 0-mismatch alignments of the same read (mm2's
  N even sits at the EXACT block coords `[L, L+dlen)` — only the op LABEL is wrong, so
  there is not even a positional signal); NO read-evidence-weighing method (incl. any C6
  member) can separate them per-read. This is the pre-committed null, **measured** (NM==0
  on 180/180) and zero-evidence **by construction** (the IID matched-intron symmetry of
  (c)), not asserted.
- The only deciding input the premise actually needs is an **EXTERNAL variant catalog
  (a matched germline VCF)**. On THIS corpus the VCF == the truth label → **circular**;
  the SIM cannot validate a PROCEED for a VCF-consulting emission.
- **Do NOT build** a variant-aware emission *as a read-evidence member* against this
  corpus. **Scope correction (panel):** the "specificity forbids it" argument applies to
  a *blunt* abstain/relabel-near-any-variant rule (which (d) shows suppresses real
  junctions, recall ~1.0). It does NOT forbid a **catalog-gated** N→D relabel — that rule
  has specificity 1.0 even on this corpus (real introns are absent from the catalog); its
  only defect here is circular recall (catalog == truth on the SIM). The catalog-gated
  emission is therefore a **DEFER-to-real-data integration milestone** (residue #1 below),
  not a corpus PROCEED and not a "never build." C6 joins C2/C3/C4 as a gate-refuted (as a
  read-evidence member) / deferred (as a VCF integration) facet.

### Named real-data residue (DEFERRED to a measured trigger — mirrors C4)
- The SIM is zero-evidence **by construction**: random block interior, single read, no
  splice-site strength beyond GT..AG, no across-read zygosity/retention, no matched VCF.
  Real C6 has signals the SIM omits. C6 becomes live again ONLY when a **MEASURED**
  real-data trigger appears AND the deciding input is genuinely INDEPENDENT of truth:
  1. **Matched germline VCF** (DNA-seq of the same sample) — a *catalog-gated* emission
     (emit D, not N, ONLY at a catalogued deletion allele) is an *integration*, not a
     read-discovery gap. This is an **ACTIVE milestone, not dormant**: the FDR inflation
     is real and total (180/180, fp_variant_adjacent=180), the fix is well-defined, and
     matched VCFs are routine in the stated targets (SMA SMN1/SMN2 patient data,
     known-variant cell lines). Gate the BUILD decision on a measured real-data FDR
     reduction (where the VCF is INDEPENDENT evidence, not the truth label), not on this
     SIM. **Required specificity check (do not skip):** the catalog-gated rule's
     specificity is 1.0 *only on this corpus* — corpus-limited, because the matched real
     introns are variant-free BY CONSTRUCTION, so the SIM structurally cannot exhibit the
     rule's true failure mode: **a real (novel) junction that COINCIDES with a catalogued
     deletion allele → N→D relabel converts it to an FN.** Deletion polymorphisms do occur
     in/near introns, so this case is real. Any catalog-gated build MUST be specificity-
     tested against "real junctions overlapping catalogued variants" before shipping.
  2. **Across-read DNA-level evidence** (het allele balance / retention split across
     many reads at the locus) — a multi-read signal the per-read SIM **structurally
     cannot model** (`coverage=1`, controlled.py:656). Half of C6's designed evidence is
     thus untested here, which is itself a reason the corpus REFUTE is scoped, not global.
  3. **Splice-site strength** beyond the GT..AG dinucleotide (branch point, PPT,
     U2/U12 PWM) — if a real intron scores strong and a deletion-block scores weak on a
     calibrated motif model, that IS read-level separation. The SIM's random interior
     makes both equally (un)splice-like, so it cannot test this. A dedicated
     splice-strength stratum is the precise, named revisit locus.
- Until such a stratum/VCF is built AND a variant-aware model is shown to reduce real-data
  FDR WITHOUT FNing real junctions: **do not build a read-evidence member** (the catalog-
  gated integration is the active milestone above, gated on real-data, not this SIM).

## Adversarial panel (3 distinct angles + advisor) — ran, CONVERGED on REFUTE/DEFER
- **Probe-faithfulness:** FIDELITY CONFIRMED. minimap2 flags byte-identical to the smoke
  (`-ax splice -uf --eqx -k 14`); `gen_variant_stratum` reused unmodified; driver
  partition identical to smoke `_is_driver`; matched-intron arm faithful to
  `gen_junction_discovery_stratum`; HEADROOM 180 == scorer `fp_variant_adjacent` 180. One
  conservative non-fatal divergence (no canonicity-retry; can only understate real-intron
  recall, verdict uses the has-N count).
- **Measurement-validity:** conclusion CORRECT, honestly scoped; one framing fix
  integrated — (c) is a *disclosed null* (IID arms force the 1:1), not an at-risk test;
  the at-risk measures are (a) and (d). No answer leak; ceiling=1.0 legitimate; the
  label-vs-placement note (N at exact block coords) STRENGTHENS zero-evidence.
- **Steelman-PROCEED:** does NOT overturn the read-evidence REFUTE (a VCF is not read
  evidence; catalog-gating = PROCEED-as-integration = DEFER residue #1). It DID defeat the
  secondary "specificity forbids building anything" overreach — under catalog-gating
  specificity is 1.0 on this corpus. Integrated as the scope correction + active milestone.

## Files (this gate)
- `scripts/benchmark/c6_headroom.py` — the decoder-free headroom + tension + specificity probe.
- `scripts/benchmark/c6_headroom_result.txt` — the result table (reps=60, seed=7).
- M1-light, decoder-free, fitness=truth. NO shipped/product code added (reused
  `gen_variant_stratum`; matched-intron arm is probe-local). Smoke gate GREEN, (E) unchanged.
