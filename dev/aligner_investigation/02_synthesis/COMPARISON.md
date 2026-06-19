# Long-Read Splice Aligner Comparison — Similarities & Differences

> **Build note:** code-level claims here were verified against `master`; see `../CORRECTIONS_vs_DRS_BUILD.md` for re-verification vs `origin/drs-validation-rebuild`. The aligner-flag lists below are CONFIRMED on the build (deSALT no `-x`/`-G`; minimap2 set + `-y`; gapmm2 `-i`; uLTRA `--ont --disable_infer`). **Corrected inline:** mapPacBio `intronlen=50` → `intronlen=10` and explicit `maxindel=max(200000, max_intron)` (`multi_aligner.py:749,754`); penalty tables are bundled and `--Scer`-auto-resolved (not absent); winner-selection runs HP-edit-distance (Path A) in production.

**Synthesis scope.** Five aligners in RECTIFY's correct-first ensemble for ONT
DRS/cDNA *S. cerevisiae* RNA: **minimap2, deSALT, uLTRA, gapmm2, mapPacBio (BBMap)**.
Cross-referenced against the PacBio (pbmm2/Iso-Seq) and ONT-DRS ecosystem dossiers.

**Sources.** `01_investigation/{minimap2,deSALT,uLTRA,gapmm2,mapPacBio_bbmap,pacbio_ecosystem,ont_drs_ecosystem}.md`;
`CLAUDE.md` (correct-first win rates: **deSALT 78.9%, mapPacBio 18.2%, uLTRA 2%, gapmm2 0.8%, minimap2 0.1%**);
ground-truth invocations in `rectify/core/align/multi_aligner.py` (`run_minimap2` L217, `run_map_pacbio` L418,
`run_bbmap` L616, `run_gapmm2` L1052, `run_ultra` L1265, `run_desalt` L1506 — all verified against the dossiers).

**Convention.** *FACT* = traceable to a source dossier/paper/RECTIFY source. *INFER* = synthesis reasoning.
Disagreements and shaky claims are flagged **⚠ FOR REVIEW** for the adversarial pass.

---

## 1. Master Comparison Table

| Axis | **minimap2** | **deSALT** | **uLTRA** | **gapmm2** | **mapPacBio (BBMap)** | **GMAP** *(contrast — not in panel)* |
|---|---|---|---|---|---|---|
| **Seeding strategy** | (w,k)-**minimizers**; splice preset w=5, RECTIFY k=14 | **MEMs** = "Match Blocks" recovered from RdBG unipath hits (l=15, s=5, a=8) | **NAMs via strobemers** (namfinder; current default) — legacy MEMs (slaMEM) optional | inherits **minimap2 minimizers** (delegates body alignment) | **all-k-mer** index (every reference k-mer), k=13 default; **multi-k-mer seed-and-extend** voting | **sampled oligomer (k-mer) hash**, k=15, **sampling every 3 bp** (fixed stride, positionally deterministic — not window-minimum) |
| **Index data structure** | open hash table: minimizer-hash → (pos,strand) list; bucket_bits=14 | **Reduced de Bruijn graph (RdBG)** of unipaths (deBGA, k=22) + hash over unipaths; 3 lookups | **annotation segmentation**: GTF→"parts"→"segments" tiling + splice-graph junction tables; pickled `database.db` + augmented `refs_sequences.fa` | minimap2 `.mmi`/mappy in-memory index | hash: k-mer key → list of genomic sites (all k-mers); disk-cached, shareable (≈6 B/ref base) | **gamma-pointer / offset hash** (basesize=12 prefix + 3 bp gamma layer) + positions array + **IIT** (Interval Index Tree) splice/SNP files; ~4 GB build @ k=15 (mammalian; trivial on yeast) |
| **Chaining / colinear algorithm** | **SDP** (sparse DP) `f(i)=max[f(j)+α−β]`, h=50 cap; 2021 RMQ long-gap path | **SDP** over Match Blocks → coarse "alignment skeleton" (pass 1) | **two-pass colinear chaining** (Mäkinen Alg 15.1); RMaxQST segment trees, O(n log n); **MAM chaining with overlaps+gaps** (pass 2) | minimap2 SDP for the body; **no own chainer** | **multi-k-mer hit voting** (tally densest sites) — not formal SDP; seed clusters nucleate DP | **oligomer chaining + diagonalization** (group seeds by genomic−query diagonal → candidate regions; chain same-diagonal HSPs allowing intron/indel jumps) |
| **Base-level DP engine** | **KSW2** (Suzuki–Kasahara diagonal, SIMD/SSE); `ksw_exts2_sse` splice DP | ksw2-style **SIMD banded SW** (own build), Z-drop=400, bw=500 | **parasail** (`sg_trace_scan_16`, semi-global SIMD) + **edlib** (MAM scoring; HW mode for refs >20 kb) | **edlib** (`mode=HW`, k=0, `task=path`) for **terminal-exon refinement only**; body DP is minimap2/KSW2 | **full affine-transform banded global DP** (custom matrix: match +100 / mismatch −127 / 2nd consec −51) | **sandwich DP** (two DPs from both anchored exon ends toward the gap, joint base+splice optimum) + **nucleotide-level DP** (`dynprog*.c`); single-threaded per alignment, no SIMD in classic path |
| **Splice / intron model** | dedicated intron DP state `x2` (open q2=32, **extend e2=0** → free intron); two-piece affine | two-tier gap `-O 2,32 / -E 1,0` (cheap indel + near-flat long-gap = intron); Z-drop decides exon end | **exon-segment tiling**: gap between MAMs classified intron (annotated or >10 bp) vs deletion (≤10 bp) in `find_exons()` | synthesizes `~gtNNNNag` / `~ctNNNNac` intron tokens around recovered terminal exons | **no splice scorer**: introns discovered as **scored long deletions**, then `D→N` reclassified at `intronlen≥10` (build; was `50` on master) | **sequence-based splice-site model** evaluated on the genome: sandwich DP slides the exon/intron boundary to the joint (base-alignment + splice-strength) optimum; `--min-intronlength` (9) = N-vs-D boundary; `-K` max intron (man 1000000 / src 200000 ⚠) |
| **Junction signal (GT-AG / annot / cross-read)** | **GT-AG donor/acceptor scoring** (d(i)/a(i) tables); `--junc-bed`+`--junc-bonus 9` (soft); `--splice-flank=no` in RECTIFY | **GT-AG via `-R 9`** non-canonical penalty; **cross-read exon pooling** (all skeletons integrated → consensus boundaries) — *unique* | **exact annotated junctions** snapped to GTF `valid_introns_sites`; novel = gap>10 bp; FSM reads privileged | **canonical-dinucleotide bias** in edlib slide-scan (`find_all_splice_CT`, GT/AG analogues) | **none** (no GT-AG model); RECTIFY's Module 2H refines its `N`-ops downstream | **canonical-as-tiebreaker reward** (GT-AG > GC-AG/AT-AC > non-canonical, never gated); adaptive `--canonical-mode 0/1/2` (conditional down-weights bonus on clean reads); MaxEnt novel-splice scoring in modern builds; **per-read only, no cross-read consensus** |
| **3'-end / poly-A behavior** | `--end-bonus 0` → **soft-clips poly-A**; 3'=`reference_end`, drifts in A-runs | **no poly-A awareness**; aligns A-tail as ordinary seq → can slip 3' into genomic A-run | **strips poly-A** (`remove_read_polyA_ends`, ≤100 bp window) before align → 3' = trimmed end, no CPA model | inherits minimap2 soft-clip; edlib only fixes the orphan terminal exon | **column-optimal DP boundary** → best raw 3'-end precision of the panel (its edge) | **none** — no poly-A awareness; aligns the A-run as ordinary sequence → can slip 3' into a genomic A-run; splice strength helps internal junctions, **nothing** for the unspliced terminal exon / CPA |
| **Error-model assumptions** | low A=1/B=2 (noisy-read scoring); ONT `splice` preset, not `splice:hq` | `-x` presets (ccs/clr/ont1d/ont2d/null); **RECTIFY runs default `null` ~13%** ⚠ (no ONT preset passed) | `--ont`: min_acc=0.6, s=9, mm2_ksize=14 (lower acc threshold than IsoSeq) | minimap2 `splice` preset error model | long-read tuned, **indel-dominated** affine weights ("recommended for Nanopore") | strong **SNP/indel + divergence tolerance** (DP absorbs mismatches/small indels; SNP-tolerant IIT mode; `--cross-species`); designed/validated on **mRNA/EST + PacBio**, *not* ONT homopolymer profile |
| **Annotation requirement** | optional (`--junc-bed`); novel-junction capable | **GTF-free** de novo exon inference (RECTIFY deliberately omits `-G`: yeast GTF→SIGSEGV) | **REQUIRED** (GTF with gene/transcript/exon; `--disable_infer`); mm2 fallback for >10% out-of-annotation reads | none (uses minimap2 splice preset) | none (D→N threshold only) | **not required** — fundamentally annotation-free (exon/intron structure from genome + splice model); optional known-splicesite IIT can assist |
| **Output format** | SAM (→ `samtools sort -n`) | SAM (→ sort → **dedup** → calmd) | `ultra.sam` (→ sort → calmd) | **PAF + `cs` only** (no SEQ/QUAL) → RECTIFY rebuilds BAM, injects sequences | SAM with **`=`/`X` CIGAR ops** (→ sort -n) | **SAM/GFF/PSL/etc.** (`-f samse/sampe/gff3_gene/…`); Iso-Seq recipe uses `-f samse -n 0` |
| **Time / memory profile** | **fastest**; near-linear; yeast index trivial | slower (2-pass + SIMD DP); RdBG ~35 GB on human, **tiny on yeast** | **heaviest/brittle**: pickled DB + augmented ref + per-gene MAM sets; namfinder+mm2 deps | cheap marginal (only ~7% of termini refined) | **slow** (JVM + exhaustive DP); RAM-hungry (≈6 B/base + per-thread DP), `-Xmx32g` cap; 6 kb read ceiling | **slowest of all — disqualifying**: per-read nucleotide DP, ≥**30× slower than minimap2** (Li 2018); read-parallel only, no within-alignment SIMD; ~4 GB index build @ k=15 (mammalian) |
| **Determinism** | deterministic (single primary; `--secondary=no`) | needs **post-hoc dedup** (emits each primary N×) | mapq from competing alns; deterministic per index | name-sort; **~0.02% reads dropped** (cs-overrun); duplicate-UUID crash without dedup | `ambiguous=best` → deterministic single record; emits **unmapped duplicates** (need guard) | **deterministic per read** (fixed sampling stride, deterministic diagonalization + DP tie-breaks; thread count changes throughput only, never alignments) |

---

## 2. What They SHARE

**(a) The seed–chain–align paradigm.** All five are instances of the same three-stage skeleton:
find short exact/approximate anchors (seeds), order collinear anchors into a chain, then fill gaps
with base-level DP. minimap2, deSALT, uLTRA, and (transitively) gapmm2 implement it explicitly;
mapPacBio is the loosest fit — its "chaining" is multi-k-mer hit-density voting rather than a formal
SDP recurrence, but the seed→localize→DP structure is the same *(INFER)*. **GMAP** (the external
contrast aligner, **not in RECTIFY's panel**) is also a clean instance of this paradigm: sampled-oligomer
hash seeds → **diagonalization + oligomer chaining** → sandwich/nucleotide DP. So the seed→chain→align
skeleton is the shared backbone of *every* long-read splice aligner surveyed here, panel or not.

**(b) Sparse-DP chaining lineage.** minimap2 and deSALT share a near-identical SDP recurrence —
`f(i)=max[f(j)+α(j,i)−β(j,i)]` with a gap term that makes long reference gaps cheap (the intron-friendly
`log2`/two-tier shape). uLTRA's two-pass colinear chaining descends from the same colinear-chaining
theory (Mäkinen *Genome-Scale Algorithmic Design* Alg 15.1) but **generalizes it to allow overlaps and
gaps** — its stated novel contribution over forbid-overlap SDP. So three of the five sit on one
algorithmic family tree.

**(c) Two-piece / two-tier affine gap cost = the shared intron trick.** minimap2 (`-O 2,32 -E 1,0`),
deSALT (`-O 2,32 -E 1,0`), and pbmm2's ISOSEQ preset all use the *identical* two-piece gap shape: a cheap
linear piece for ordinary indels and a near-flat second piece (open ~32, extend 0) so a single intron is
penalized by a constant open cost, not per base. This is what lets one banded DP emit a long `N` without
the score collapsing. mapPacBio achieves the same end (long scored gap → reclassify to `N`) with a
different, single affine matrix + an `intronlen` post-filter.

**(d) GT-AG canonical prior.** Four of five bias toward canonical splice signals: minimap2 (donor/acceptor
penalty tables), deSALT (`-R 9` non-canonical penalty), uLTRA (snaps to annotated GT-AG coordinates),
gapmm2 (synthesizes `gt..ag`/`ct..ac` tokens). **mapPacBio is the sole exception** — it has no GT-AG model
at all. Notably, none treat GT-AG as a hard constraint; all permit non-canonical junctions at a cost,
preserving novel-junction discovery — the same design principle RECTIFY's Module 2H enforces ("annotation/
canonical are tie-breakers, never gates"). **GMAP is the strongest external precedent for exactly this
philosophy:** its sandwich DP adds a *reward* (never a gate) when the chosen boundary yields a canonical
dinucleotide, ranks `GT-AG > GC-AG/AT-AC > non-canonical`, and its `--canonical-mode 2` (conditional)
*down-weights* the canonical bonus on clean reads so it won't bend a good alignment to manufacture a GT-AG
site — the same "canonical-as-tie-breaker, sequence-first" stance RECTIFY codifies in Module 2H's
`_CANONICAL_HP_PRIOR = 0.5`. GMAP demonstrates that this prior can live *inside* the aligner, per read.

**(e) Modern CIGAR + calmd.** deSALT, uLTRA, gapmm2, and mapPacBio outputs all pass through RECTIFY's
`_apply_calmd_eq` to normalize MD / `=` / `X` ops; minimap2 emits `--MD` directly. The pipeline
standardizes every aligner onto explicit match/mismatch CIGARs so the downstream CIGAR-surgery /
correction code has uniform ops to operate on.

---

## 3. What Genuinely DIFFERS

**(a) De Bruijn graph index vs minimizer/k-mer hash.** deSALT is the only aligner whose index is a
**graph** (RdBG of collapsed unipaths). The decisive property: one hash lookup returns *every* genomic
copy of a matched unipath, surfacing all candidate loci (repeats, paralogs, short exons) cheaply, where
minimizer sampling (minimap2) can simply miss the true exon's seed. This breadth-of-candidates is the
root of deSALT's sensitivity *(deSALT dossier §1, INFER→78.9%)*. mapPacBio's all-k-mer index is the
opposite extreme — exhaustive (no sampling) but flat (no graph structure, no splice awareness).

**(b) Annotation-guided vs de-novo.** uLTRA is **annotation-REQUIRED** — it tiles the GTF into exon
segments and chains the read across them, *injecting all annotated exons of the candidate gene into Pass-2
even if no seed landed in them* (the headline 8-nt micro-exon recovery). This gives exact canonical
junctions on FSM reads but caps it structurally: reads diverging from the GTF fall through to a minimap2
safety net, so uLTRA *cannot beat minimap2* on novel structure. deSALT does GTF-free de-novo exon
inference; minimap2/gapmm2/mapPacBio are de-novo with at most a soft annotation hint. On compact,
well-annotated yeast this caps uLTRA at ~2%.

**(c) Full-SW vs heuristic gap-fill.** mapPacBio runs a **full affine-transform DP of every candidate
window** and reports the globally column-optimal boundary. minimap2/deSALT/gapmm2 use *heuristic*,
chain-anchored, banded extension that stops at the first good-enough chain — which can place an indel a few
bases off when a homopolymer admits several near-equal placements. This is exactly why mapPacBio recovers
**precise 3' ends / CPA boundaries** (its 18.2% niche) despite having no splice model.

**(d) Two-pass cross-read junction refinement — deSALT's exclusive mechanism.** deSALT integrates the
alignment skeletons of *every read in the batch* (B=655350) before re-aligning each read: pooled MB
projections define consensus exon boundaries that get applied back to individual noisy reads, producing
**homogeneous** junctions (same CIGAR placement across reads at a locus). No other aligner in the panel
looks across reads — minimap2/gapmm2/mapPacBio commit per-read in isolation; uLTRA pools the *annotation*,
not the reads. This is the single biggest structural differentiator and the mechanistic core of the win-rate
spread *(deSALT dossier §3, INFER)*.

**(e) Wrapper vs native, and output completeness.** gapmm2 is a **thin minimap2 wrapper** that only adds an
edlib terminal-exon refinement pass and emits **PAF (no SEQ/QUAL)** — forcing RECTIFY to reconstruct BAMs
and inject sequences. mapPacBio is a wholly independent Java engine. This wrapper status is why gapmm2's
body alignment *is* minimap2 — so wherever minimap2 is right, gapmm2 ≈ minimap2, and a dedicated aligner
usually scores ≥ as well (the structural reason for gapmm2's 0.8%).

**(f) Splice-placement *philosophy* — GMAP as the "splice-model-first" pole (external contrast).** Across the
panel, junction placement is driven by four distinct primary criteria: minimap2 is **gap-cost-first** (an
intron is emitted when the two-piece affine break-even is crossed; GT-AG is only a per-position bonus); uLTRA
is **annotation-first** (exon coordinates supplied by the GTF); deSALT is **cross-read-consensus-first** (pooled
exon projections fix the boundary before re-alignment); mapPacBio is **score-first / model-free** (long gap by
DP score, then `D→N`). **GMAP** (the external contrast, **not in the panel**, so **no win-rate**) occupies the
remaining, purest pole: **splice-model-first** — the sandwich DP brackets the intron between two anchored exon
ends and slides the boundary to the position that *jointly* maximizes base-alignment score + sequence
splice-strength, even when basecall errors sit at the junction. This is the cleanest external example of a
sequence-first, canonical-aware splice solver, and a strong precedent for RECTIFY's own Module 2H — except that
GMAP does it per read, inside the aligner, at alignment time, whereas RECTIFY decouples it as a post-hoc
refinement (and so avoids GMAP's speed cost). **GMAP's sandwich DP — a *joint* base+splice optimization that
slides the boundary onto the best-supported site — is the closest external analogue to RECTIFY's sequence-first
junction refiner (`junction_refiner.py`); every panel aligner instead treats the splice signal as a bolt-on
penalty/bonus on top of a base-alignment that is decided separately.** *(INFER; GMAP is investigated only as an
algorithmic contrast / candidate aligner, not a panel member — see `01_investigation/gmap.md` §8.)*

---

## 4. Algorithm → Empirical Win-Rate Mapping

| Aligner | Win rate | Mechanistic explanation |
|---|---:|---|
| **deSALT** | **78.9%** | RdBG breadth-of-candidates (rarely misses true/short/repeat exons) **+** two-pass coarse-then-precise (commits junctions only after exon structure is known) **+** *cross-read exon pooling → homogeneous, consensus junctions*. RECTIFY selects per-read winners on **corrected 3'-end position** and Module 2H rewards sequence-supported, consistent junctions — precisely the trait deSALT optimizes. deSALT does a mini cross-read junction-correction *internally* before RECTIFY's correction runs → compounding head start. *(INFER, deSALT §2-5)* |
| **mapPacBio** | **18.2%** | **Full affine-scored global DP** of every candidate window → **column-optimal indel & terminal-exon boundary** = best raw 3'-end (CPA) precision in the panel. Indel-dominated long-read error weights match ONT's spectrum; multi-k-mer voting localizes noisy reads robustly. Loses to deSALT on *spliced* reads (no GT-AG/graph model — finds introns only as scored gaps), wins on **unspliced terminal exons & exact 3' ends**. *(FACT DP + INFER, mapPacBio §"Why second-best")* |
| **uLTRA** | **2%** | Wins only reads on **well-annotated multi-exon isoforms** where snapping to GTF coordinates beats noisier de-novo aligners (zero junction wobble on FSM reads). Structurally capped: novel/intergenic reads → minimap2 fallback (can't beat mm2 there); yeast is compact + well-annotated + few isoforms → narrow niche. Strips poly-A, offers no CPA precision — contributes a *junction backbone*, not a terminus. *(FACT+INFER, uLTRA §"Weaknesses")* |
| **gapmm2** | **0.8%** | Body alignment **is** minimap2, so for ~93% of reads it offers nothing distinctive; it only changes ~7% of termini (edlib recovery, 5'-biased). RECTIFY's junction-proximity penalty docks minimap2/gapmm2 for near-junction mismatches, so even on refined termini a clean-junction aligner wins. Value = the rare read where its edlib terminal recovery is **uniquely** correct (independent failure mode). *(FACT+INFER, gapmm2 §6-7)* |
| **minimap2** | **0.1%** | `--end-bonus 0` + noisy reads → **soft-clips poly-A**, `reference_end` drifts; HP indel noise distributes indels arbitrarily in A-runs (low A=1/B=2 doesn't pin a unique boundary); intron-boundary jitter ±1-few bp. Fast, accurate *junction baseline* but **not optimized for single-base 3'-end precision** — the exact gap RECTIFY fills. It loses to every dedicated aligner on the corrected-3'-end metric RECTIFY scores. *(FACT+INFER, minimap2 §6)* |

**The unifying logic.** RECTIFY selects on **corrected 3'-end position**, not raw alignment score. Two
properties dominate that metric: (1) *junction sharpness/consistency* (so the corrected CPA sits at a clean
exon edge, not in an intron) and (2) *terminal-boundary precision*. deSALT owns (1) via cross-read pooling;
mapPacBio owns (2) via full SW. Everything else is a baseline (minimap2), a wrapper of the baseline
(gapmm2), or a narrow annotation specialist (uLTRA). The 79/18/2/0.8/0.1 spread is the direct numerical
shadow of these complementary mechanisms — and the spread is *why an ensemble + correct-first selection
works*: no single aligner is uniformly best at the noisy 3' terminal exon, and the failure modes are
complementary (minimap2 over-extends/soft-clips; deSALT/uLTRA recover small/terminal exons; mapPacBio nails
the boundary).

---

## 5. Cross-Reference: PacBio vs ONT (why this panel is ONT-tuned)

- **Poly-A timing is the architectural fork.** PacBio Iso-Seq trims poly-A **before** alignment
  (`isoseq refine --require-polya`, `--min-polya-length 20`); ONT DRS sequences the tail **in-read**.
  RECTIFY's `--drs` Step 0/Step 4 *manufactures* the PacBio "tail-already-trimmed" representation for ONT
  while preserving tail-length metadata PacBio discards. The five aligners here all see the in-read A-tail
  (or, under `--dT-primed-cDNA`, an off-read A-context), which is *why* every one of them needs RECTIFY's
  downstream 3'-end correction. *(pacbio §3, ont_drs §2)*
- **Tuning direction is opposite to `splice:hq`.** RECTIFY runs the base `splice` preset (B=2, O=2,32, C=9)
  + `-k14` + `-G 5000` + `--splice-flank=no` — the *noisy-ONT* direction. PacBio's `splice:hq`/ISOSEQ raises
  mismatch/gap penalties and lowers non-canonical cost (trust the accurate read). The empirical penalty
  tables and `_CANONICAL_HP_PRIOR=0.5` are "R10.4.1-specific" and must **not** transfer to HiFi.
- **HiFi has no systematic homopolymer-deletion bias** — the exact artifact this panel + correction stack
  exists to repair. On HiFi most RECTIFY corrections would be near-no-ops *(pacbio §4, INFER)*.

---

## 6. Disagreements & Shaky Claims (flagged for adversarial review)

1. **⚠ deSALT error preset.** Dossier confirms (and RECTIFY source L1560 verifies) deSALT runs `aln` with
   **no `-x`**, i.e. the `null` default (~13% error model) — *not* an ONT preset (`ont1d`/`ont2d`), despite
   processing ONT DRS data. The dossier flags this as possibly sub-optimal; it is the single most concrete
   "is this intentional?" item. (Counter-point: deSALT still wins 78.9%, so the null model is evidently
   good enough — but no one has A/B-tested `-x ont2d`.)

2. **⚠ AT-vs-CG deletion-rate direction.** `CLAUDE.md`'s penalty table says **AT runs have HIGHER deletion
   rates than CG** (D/AT HP=1=0.37 lower penalty = more expected deletions); the broader ONT literature
   reports A/T homopolymers are *basecalled more accurately* overall. The ont_drs dossier argues these are
   **not contradictory** (deletion-specific rate vs aggregate accuracy measure different operations). This
   reconciliation is *inference*, not a cited measurement — worth an adversarial double-check.

3. **⚠ mapPacBio as "seed-chain-align."** Calling BBMap's multi-k-mer voting a "chain" stage is a
   convenient taxonomy but BBMap implements no formal SDP recurrence. The seed→DP structure holds; the
   "chaining" label is loose *(INFER)*.

4. **⚠ gapmm2 `-i` semantics.** The gapmm2 dossier argues `-i 5000` sizes gapmm2's *terminal-refinement
   window*, and is *not necessarily* passed to minimap2 as `-G`. The minimap2 sub-pass may therefore run at
   its **splice-preset default max-intron (200k)**, not 5 kb — a subtle divergence from the other aligners'
   explicit `-G 5000`. Reconstructed from README + RECTIFY code, not a line-by-line read of the wheel.

5. **⚠ Win rates are single-dataset.** The 78.9/18.2/2/0.8/0.1 figures trace to one validated lineage
   (mini-1k sidecar / chunk003); they are presented as biology but are one experiment. The pre-v3.3.0
   `index_col` merge bug meant minimap2/gapmm2/uLTRA read_ids were corrupted and *never properly competed*
   — so any win rate measured before that fix is suspect. ⚠ Confirm the cited rates are post-fix.

6. **⚠ "GTF-free deSALT" + yeast GTF SIGSEGV.** deSALT *can* take `-G`, but RECTIFY disables it because the
   yeast GTF crashes it. The "GTF-free is a strength" framing is partly a workaround rationalization.
