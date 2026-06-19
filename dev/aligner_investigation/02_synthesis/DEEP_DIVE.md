# Long-Read Splice Aligner Deep Dive — Granular Operation, Per-Aligner

**Scope.** Function/equation-level operational walkthroughs for the five aligners in RECTIFY's
correct-first ensemble, distilled from the `01_investigation/` dossiers and verified against
`rectify/core/align/multi_aligner.py`. Then a cross-cutting advantage/disadvantage matrix and the
crux section: **failure modes at the 3' end / CPA / splice boundary**, and why the ensemble +
correct-first selection recovers truth.

**Convention.** *FACT* = source/paper/RECTIFY-source. *INFER* = synthesis. ⚠ = flagged for review.
Correct-first win rates (`CLAUDE.md`): **deSALT 78.9% · mapPacBio 18.2% · uLTRA 2% · gapmm2 0.8% · minimap2 0.1%**.

---

## 1. minimap2 — fast minimizer seed-chain-align baseline

**RECTIFY invocation (verified L217):**
`minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD -y -t N [--junc-bed … --junc-bonus 9] genome reads`

### Seed → Chain → Align
1. **Seed (minimizers).** A (w,k)-minimizer = smallest-hash k-mer in each window of `w` consecutive
   k-mers. Splice preset sets w=5 (≈1 minimizer per 3 bp), RECTIFY forces k=14 (denser, more sensitive
   on noisy DRS, helps short terminal exons). Index = open hash table, minimizer-hash → (pos,strand) list.
   Repetitive minimizers (top `-f 2e-4`) filtered before chaining.
2. **Chain (SDP).** Anchors sorted by reference pos; chain score
   `f(i)=max{ max_{j<i}[f(j)+α(j,i)−β(j,i)], w_i }` where `α`=matched bases (capped at seed length),
   `β=γ_c(Δref−Δquery)` with gap cost `γ_c(l)=0.01·w̄·|l|+0.5·log2|l|` (the `log2` term makes one long
   intron-gap cheap). Inner loop bounded: start at `i−1`, stop after **h=50** non-improving predecessors →
   near-linear. 2021 RMQ path handles longer gaps.
3. **Align (KSW2 splice DP).** Between adjacent anchors, `ksw_exts2_sse` fills gaps with a banded diagonal
   DP using a **dedicated intron state `x2`** (open cost q2=32, **extend e2=0** → free intron extension).
   Per cell `z=max(z,a,b,a2a)` chooses normal-deletion vs open-intron. **GT-AG signal:** donor `d(i)`=0 for
   `GTA/GTG`, p/2 for `GTC/GTT`, p else; acceptor symmetric (`CAG/TAG`=0); p=`-C`=9. `--junc-bed`+
   `--junc-bonus 9` add a *soft* bonus at annotated sites (tie-breaker, novel-junction-capable).
   **Z-drop** (`zdrop=200, zdrop_inv=100`) breaks the alignment where score drops too fast but tolerates a
   single long gap (introns survive). Splice gaps → `N` ops.

### Key RECTIFY flag rationale
- `-uf`: DRS is stranded (3'→5'), forward-only halves splice DP work; matches ONT's own DRS recipe.
- `--splice-flank=no`: disables the GT**A/G**..**C/T**AG extra-base bonus that would nudge boundaries ±1 bp
  toward the statistically-favored flank — RECTIFY wants junction placement handed to its own Module 2H,
  not nudged by an uncontrolled aligner prior. (CLAUDE.md: "important for 3' end accuracy.")
- `--end-bonus 0` (default): **no reward for extending to read ends** → no incentive to push the 3'
  alignment into a poly-A run → soft-clips it. *(Root of its 3'-end drift.)*

### Advantages / Disadvantages
- **+** Fastest; principled splice DP; annotation-assisted yet novel-capable; de-facto ONT standard.
- **−** `--end-bonus 0` soft-clips poly-A → `reference_end` drifts; HP indel noise placed arbitrarily
  (low A=1/B=2 doesn't pin a boundary); intron-boundary jitter ±1-few bp; `--secondary=no` discards
  alternatives. **Wins 0.1%** — a junction baseline, not a CPA-precision tool.

---

## 2. deSALT — RdBG index + two-pass cross-read refinement (the winner, 78.9%)

**RECTIFY invocation (verified L1506):** `deSALT aln -t N -f <tmpfile> -o <sam> <index_dir> <reads>` —
**no `-x` preset** (runs `null` ~13% model ⚠), **no `-G`** (yeast GTF→SIGSEGV).

### Index (RdBG, from deBGA)
A **Reduced de Bruijn Graph** over reference (k+1)-mers (k=22): non-branching k-mer chains (unipaths) are
collapsed into single unitig nodes; multi-copy unitigs record all genomic start positions. A hash over
unipaths exposes 3 lookups: (i) k-mer→unipath+offset, (ii) unipath→**all** genomic positions, (iii)
unipath→sequence. **One lookup returns every genomic copy** of a matched unipath → all candidate loci
(repeats, paralogs, short exons) surfaced cheaply — the root of deSALT's sensitivity *(INFER)*.

### Pass 1 — alignment skeleton
1. **Match Blocks (MBs):** maximal exact matches recovered as unipath hits (l=15, s=5, a=8, n=50).
2. **SDP chaining** composes MBs into a coarse, gapped **alignment skeleton** — DP only *between*
   non-overlapping seeds (cheap edit estimate, not full base-DP). Large gaps = candidate introns.
   (c=30 min score, g=2000 read gap, I=200000 max intron.)

### Pass 2 — refined alignment (+ the decisive cross-read step)
1. **Exon inference by integrating ALL skeletons** — MBs from *every read in the batch* (B=655350) are
   projected onto the genome; pooled projections define exon regions. **GTF-free, de novo.**
2. **Re-align each read** against the inferred exons; a **local hash query with k=8** catches small/
   micro-exons and weak-seed exons Pass 1 missed (`-e 5` extension). Output = full-length spliced alignment.

### Base DP & splice scoring
ksw2-style **SIMD banded SW**, two-tier gap `-O 2,32 / -E 1,0` (cheap indel + near-flat long-gap = intron),
Z-drop=400, bw=500. **GT-AG via `-R 9`** non-canonical penalty (canonical = free).

### Why it wins (the cross-read mechanism)
A single noisy read may have a frayed exon end; when dozens of reads at a locus contribute MBs, their
projections **agree** on the true boundary, and that consensus is applied back to the noisy read →
**homogeneous junctions** (same CIGAR placement across reads). RECTIFY selects on **corrected 3'-end** and
Module 2H rewards sequence-supported, consistent junctions — exactly deSALT's output trait. deSALT performs
a mini cross-read junction-correction *internally before* RECTIFY's correction runs → compounding head
start, hence ~4/5 of all reads *(INFER)*.

### Advantages / Disadvantages
- **+** Best sensitivity + junction homogeneity; GTF-free; micro-exon rescue (k=8); GT-AG-biased DP; tiny/
  fast on yeast.
- **−** RdBG ~35 GB on human (irrelevant for yeast); slower than minimap2; parameter-sensitive;
  **engineering-fragile in RECTIFY** (emits each primary N× → dedup; misparses gzip FASTQ → decompress
  first; `-G` SIGSEGV; `-f` temp must be local non-NFS; conda allocator double-free → strip
  `LD_LIBRARY_PATH`). **No poly-A awareness** → can slip 3' into a genomic A-run (RECTIFY's 2B/2E/2G fix the
  tail; deSALT nails the structure). Complementary by design.

---

## 3. mapPacBio (BBMap) — full affine-SW seed-and-extend (the precision aligner, 18.2%)

**NOT minimap2's `map-pb` preset** — this is Bushnell's Java `align2.BBMapPacBio`, a distinct engine.

**RECTIFY invocation (verified L500-516):** `mapPacBio.sh ref=… in=… out=… path=<shared bbmap_index>
fastareadlen=100000 intronlen=50 minratio=0.4 -Xmx32g` (+ pre-split reads >6 kb, stitch after).

### Seed → vote → DP
1. **Index:** **all-k-mer** hash (every reference k-mer, k=13), key→genomic-site list; ≈6 B/ref base;
   disk-cached + shared across parallel jobs.
2. **Multi-k-mer seed-and-extend** ("growing polycrystalline silicon"): collect many k-mer hits from a
   read, tally which reference sites accumulate the most (`minhits` floor), nucleate outward from dense
   clusters. Robust to errors — a single seed can be destroyed but dozens survive. `minratio=0.40`: accept
   alignments ≥40% of perfect score (long-read tolerant).
3. **Full affine-transform banded global DP** of the read vs each candidate window, custom empirical matrix:
   **match +100, single mismatch −127, 2nd consecutive mismatch additional −51** (affine, indel-class-
   specific). Finds the **globally column-optimal** placement of every indel — vs minimap2's heuristic
   chain-anchored extension that stops at first-good-enough.

### Splice handling (no GT-AG model)
The affine DP opens long reference gaps when score favors it (subject to `maxindel`); **`intronlen=50`
reclassifies any `D` ≥50 bp into an `N`** op afterward. Introns are discovered as *scored gaps*, then
relabeled — no canonical-site scorer. Emits **`=`/`X` CIGARs** (directly consumable by RECTIFY surgery).

### Why precise at the 3' end
A full DP scores the whole window and snaps the terminal-exon boundary to the single best-scoring column;
for the CPA — the exact reference base where the terminal exon ends — a globally-optimal boundary lands
closer to the true site than a heuristically-extended one. Long-read indel-dominated error weights match
ONT's spectrum. **This is its 18.2% niche: unspliced terminal exons & exact 3' ends.**

### Advantages / Disadvantages
- **+** Column-optimal indel/end placement (best raw CPA precision); tunable D→N; `=`/`X` CIGARs; shareable
  index; deterministic (`ambiguous=best`).
- **−** **6 kb read ceiling** (hard AssertionError ~6019 bp → RECTIFY pre-splits/stitches); RAM-hungry +
  JVM overhead → slow, can OOM neighbors (`-Xmx32g` cap); **no splice-graph/GT-AG** → loses to deSALT on
  heavily-spliced reads; imports DRS quirks (`pt:i:N` into QNAME → strip; **unmapped duplicate records** →
  `is_unmapped` guard). HP boundaries still admit near-equal placements (no aligner nails every HP).

---

## 4. uLTRA — annotation-segmentation colinear chaining (the specialist, 2%)

**RECTIFY invocation (verified L1405):** `uLTRA pipeline --ont --disable_infer --t N --prefix P
<ref.fa> <annotation.gtf> <reads.fastq> <outdir>` (after `_normalize_gtf_for_ultra`; namfinder vendored).

### Index = annotation segmentation (not the genome as a flat string)
GTF parsed via gffutils (`--disable_infer` requires pre-existing gene/transcript features → RECTIFY
normalizes SGD GTFs first). Exons clustered into **"parts"** (new part if `exon.start−1 > active_stop+20`),
then `get_canonical_segments()` cuts each part at every annotated boundary into **"segments"** (packed
`array("L",[chr,p1,p2]).tobytes()`). A gene = linear tiling of disjoint segments; an isoform = subset path.
Splice-graph junction tables (`all_splice_pairs_annotations` = `valid_introns_sites`) + intergenic `flanks`.
Whole index pickled to `database.db`; segment sequences → augmented `refs_sequences.fa`.

### Seed → two-pass colinear chain → align
1. **Seed:** **NAMs via strobemers** (namfinder; current default; legacy MEMs optional). `--ont`: min_acc=0.6,
   s=9, mm2_ksize=14.
2. **Pass 1 (read coverage):** `n_logn_read_coverage` — colinear chaining of seeds maximizing read coverage
   via two **RMaxQST segment trees** (O(n log n)); selects the candidate gene(s). DP:
   `C_a`=best disjoint predecessor, `C_b`=best overlapping (only non-overlapping tail credited),
   `C[j+1]=max(C_a,C_b)`.
3. **Pass 2 (MAM chaining with overlaps+gaps):** chains **MAMs** (Maximal *Approximate* Matches, edlib-
   scored) allowing fractional overlap penalty `0.1·gap`. **The decisive trick:** inject **all annotated
   exons of the candidate gene — including seedless ones** — into the chaining DP, so the chain can tile a
   micro-exon (headline **8-nt exon**) that minimap2/deSALT drop.
4. **Align:** `find_exons()` classifies each MAM-gap: in `valid_introns_sites` → annotated `N` at exact GTF
   coords; `>10 bp` → novel `N`; else `D`. Per-exon CIGAR via **parasail** (`sg_trace_scan_16`, match+2/
   mismatch−2/open3/ext1) or **edlib** (HW for refs >20 kb). FSM reads privileged (accuracy filter waived).
   **Poly-A stripped** (`remove_read_polyA_ends`, ≤100 bp window) before align.

### Why only 2%
Wins reads on **well-annotated multi-exon isoforms** where GTF-snapped junctions beat noisier de-novo
aligners (zero wobble on FSM). But: novel/intergenic reads (>10% out-of-annotation, `--genomic_frac`) fall
through to a **minimap2 safety net** → can't beat mm2 there. Yeast is compact + well-annotated + few
isoforms → narrow niche. Strips poly-A → **no CPA model**: contributes a junction backbone, not a terminus.

### Advantages / Disadvantages
- **+** Micro-exon recovery (signature); exact canonical junctions on FSM; principled overlap/gap O(n log n)
  chaining; graceful novel-read fallback (no silent mis-forcing).
- **−** **Annotation-bound** (quality capped by GTF; small upside on yeast); novel-isoform blind spot;
  poly-A stripped (no 3'-end precision); **operationally heavy/brittle** (GTF-not-GFF, no gzip, namfinder +
  mm2 deps, large pickled DB); memory/runtime inflated.

---

## 5. gapmm2 — minimap2 wrapper + edlib terminal-exon refinement (the rare-rescue, 0.8%)

**RECTIFY invocation (verified L1110):** `gapmm2 -t N -i 5000 -o out.paf genome cleaned_reads` (omits `-m`
deliberately — argparse `--min-mapq` type bug on some wheels). PAF→BAM with sequence injection.

### Operation
gapmm2 implements **no seed-chain-align core** — it runs minimap2 `splice` (mappy or, since v25.4.13, the
`minimap2` binary directly after mappy segfaults) for the body, then for each primary hit `h`:
- `if h.q_st > 0:` → unaligned 5'/left query bases → **left refinement**.
- `if len(seq) > h.q_en:` → unaligned 3'/right bases → **right refinement**.

**Refinement (edlib):** carve a reference window near the alignment edge bounded by `maxlen=-i` (RECTIFY
5000, sizing the terminal search ⚠ *not necessarily* minimap2's `-G`), run
`edlib.align(orphan_query, ref, mode="HW", k=0, task="path", additionalEqualities=degenNuc)` — **HW = infix/
semi-global** (query fully consumed, free reference end-gaps = "where does this orphan exon sit upstream of
the intron"). A `slide` scan biases the synthesized intron toward canonical `gt..ag`/`ct..ac`, stitched
into `cs` as `~gtNNNNag` / `~ctNNNNac`. README benchmark: 409 left / 63 right refined per 6926 alns (≈7%,
5'-biased — short 5' exons are the dropped ones).

### RECTIFY reconstruction
PAF carries **no SEQ/QUAL** → `_paf_to_bam` parses `cs` (`:n`→M, `*xy`→X, `+`→I, `-`→D, `~nnNNNnn`→N),
recomputes soft-clips from coords, **injects sequences** from the original FASTQ (RC for minus strand), and
**skips records where `cigar_qlen ≠ len(seq)`** (gapmm2 cs-overrun ~0.02%). Mandatory `_clean_fastq` dedup
(empty-seq + duplicate-UUID → mappy `seq()` returns `None` for both copies → `len(None)` `TypeError` crash).

### Why only 0.8%
Body **is** minimap2 → ≈minimap2 wherever minimap2 is right (93% of reads it changes nothing); only ~7% of
termini refined; RECTIFY's junction-proximity penalty docks minimap2/gapmm2 for near-junction mismatches.
Value = the rare read where its edlib terminal recovery is **uniquely** correct (independent failure mode).

### Advantages / Disadvantages
- **+** Recovers short terminal exons minimap2 drops (13 bp 5' exon demo); canonical-junction-aware terminal
  placement; HW edlib + degenerate equalities (tolerant on noisy ends); cheap marginal cost; independent
  failure mode → occasional unique correct 3' end.
- **−** PAF-only (reconstruction fragility); duplicate-UUID crash; cs-overrun drops reads; mappy segfaults;
  argparse type bug; out-of-bounds coords; rarely wins (body=minimap2); default `-i 500` too small for yeast
  (RECTIFY overrides to 5000).

---

## 5b. GMAP — sandwich-DP splice-model-first aligner (external contrast — NOT in the panel, no win-rate)

**Role.** GMAP is investigated as an **algorithmic contrast / candidate**, *not* a member of RECTIFY's
5-aligner ensemble — it has **no win-rate** here. It is the historically dominant annotation-free Iso-Seq
mapper (the de-facto PacBio mapper before minimap2) and the cleanest external example of **splice-model-first**
junction placement. *(See `01_investigation/gmap.md`.)*

**Iso-Seq invocation (ecosystem recipe, not a RECTIFY wrapper):**
`gmap -d <db> -D <dbdir> -f samse -n 0 --no-chimeras -t N input.fasta > aln.sam`

### Index = sampled oligomer hash + IIT
A *sampled* genomic k-mer lookup table (not a suffix tree): default **k=15**, two-level scheme with a
**basesize=12 offset prefix** resolved by a compressed **gamma-pointer** layer over the remaining 3 positions,
oligomers **sampled every 3 bp** (fixed stride). Build RAM scales with k (k=15 → **~4 GB** mammalian; trivial
on yeast). **IIT** (Interval Index Tree) files hold optional splice-site / SNP auxiliary data for tolerant modes.

### Stage 1 — approximate mapping (oligomer hash → diagonalization)
Each sampled query 15-mer is looked up in the genomic hash → candidate genomic positions. GMAP
**diagonalizes** these hits — groups them by `(genomic − query)` offset, so all seeds of one ungapped exon
segment fall on (nearly) the same diagonal. Well-supported diagonals define **candidate genomic regions**.
Up to `-n/--npaths = 5` paths survive (`-n 0` ⇒ 1, or 2 if a chimera is found); best = primary.

### Stage 2 — oligomer chaining + **sandwich DP** (the distinctive stage)
Within a region, scattered same-diagonal seeds are chained into ordered HSPs (ungapped exon blocks),
allowing diagonal *jumps*: a jump where genomic offset grows far faster than query offset is a **candidate
intron**; a small balanced jump is an indel. When two anchored exon segments flank an unresolved gap, GMAP
runs **two dynamic programs simultaneously from both anchored ends toward the middle** (the donor-side DP
extending right, the acceptor-side DP extending left) and solves for the single intron position that
**maximizes (alignment score + splice-site score)**. This is the structural contrast with the panel: every
panel aligner emits the intron as a long gap op inside *one* left-to-right banded DP (cost = two-piece affine
gap + a GT-AG bonus); GMAP **brackets the gap between two anchors and searches the boundary jointly from both
sides** — letting it slide the exon/intron boundary onto the position that is simultaneously the best base
alignment *and* the strongest splice signal, even when basecall errors sit at the junction.
`--min-intronlength` (9) is GMAP's own **N-vs-D** decision boundary; `-K` bounds intron span (man 1000000 /
src 200000 ⚠).

### Stage 3 — nucleotide DP + microexon + chimera
Nucleotide-level DP (`dynprog*.c`) resolves remaining mismatches/indels and the exact intron boundaries
within the Stage-2 skeleton, producing the base-perfect alignment and SAM CIGAR (`N` for introns). The DP
absorbs substantial SNPs/small indels without breaking the exon chain (a SNP-tolerant IIT mode and
`--cross-species` extend this). **De novo microexon detection:** a tiny internal exon is inserted only if a
flanking **splice-site probability** exceeds `--microexon-spliceprob` (src 0.95 / man 0.90 ⚠) — uLTRA-class
small-exon competence, but annotation-free, from the splice model. **Chimera / read-through detection:** if a
query margin (`--chimera-margin`, src 30 / man 40 ⚠) can't align to the same locus, GMAP searches for the
remainder elsewhere and reports a two-part alignment (`-n 0` forces 2 paths).

### Splice model = the algorithmic signature
Boundary placement is **sequence-based** and evaluated on the genome, with the canonical signal as a
**reward, never a gate**: `GT-AG > GC-AG/AT-AC > non-canonical`, magnitude set by `--canonical-mode 0/1/2`
where **2 (conditional)** *down-weights* the canonical bonus on high-identity reads (won't bend a clean
alignment to fake GT-AG) and *up-weights* it on noisy reads (pull an ambiguous boundary onto real signal).
Modern builds add **MaxEnt** novel-splice scoring (GSNAP→GMAP), superseding the 2005 "no probabilistic model"
core. **No cross-read consensus** — each read is refined independently against genome + model.

### Why it is not a panel fit
- **Speed is disqualifying.** Per-read nucleotide/sandwich DP, single-threaded per alignment; independent
  benchmarks put **minimap2 ≥ ~30× faster** (Li 2018). Throughput is read-parallel only (no within-alignment
  SIMD in the classic path) → wall-time-prohibitive at modern Iso-Seq/ONT depth — fatal for RECTIFY's
  chunked-alignment HPC mandate.
- **No 3'-end benefit.** GMAP has **no poly-A awareness** — it aligns the in-read A-run as ordinary sequence,
  so on ONT DRS the DP can slip the 3' exon end into a genomic A-run, leaving the CPA imprecise — *the exact
  artifact RECTIFY corrects*. Its splice model helps internal junctions but does **nothing** for the unspliced
  3' terminal exon / CPA where RECTIFY's value lies. (Same limitation as deSALT/minimap2; GMAP is no better.)
- **Per-read isolation.** No cross-read exon consensus (unlike deSALT) → junctions can be placed a few bp
  apart between reads of the same isoform — the inhomogeneity correct-first selection penalizes.

### Relevance to RECTIFY (INFER)
GMAP's sandwich DP is the **closest external analogue to RECTIFY's sequence-first junction refiner** (Module 2H,
`junction_refiner.py`): both score junctions sequence-first, treat canonical as a tie-breaker/prior never a
gate, and slide the boundary onto the best-supported site. GMAP proves the policy can live inside the aligner;
RECTIFY decouples it as a post-hoc refinement and avoids GMAP's speed cost. As a candidate aligner GMAP would
likely earn only a **low win-rate** specialist slot (its de novo microexon + divergence/SNP tolerance, on
difficult loci, competing with uLTRA's niche but annotation-free), if any — speed + no 3'-end benefit argue
against inclusion in the production panel.

---

## 6. Cross-Cutting Advantages & Disadvantages Matrix

| Aligner | Signature advantage | Signature disadvantage | Independent failure mode (why kept) |
|---|---|---|---|
| **minimap2** | Fast, principled splice DP, ONT-standard baseline | Soft-clips poly-A; 3'-end drifts; HP/junction jitter | Junction baseline; soft-clip-rescue source for consensus |
| **deSALT** | Cross-read homogeneous junctions; RdBG sensitivity; micro-exon rescue | No poly-A awareness; fragile integration; slow | Best junction structure — wins most reads |
| **mapPacBio** | Column-optimal full-SW boundary = best CPA precision | 6 kb ceiling; RAM/JVM-slow; no GT-AG model | Exact unspliced 3' ends deSALT can't beat |
| **uLTRA** | Exact GTF junctions; 8-nt micro-exon recovery | Annotation-bound; poly-A stripped; brittle | Annotated multi-exon isoform junctions |
| **gapmm2** | edlib terminal-exon recovery (short 5' exons) | Body=minimap2; PAF-only; rarely distinctive | The unique edlib-recovered terminus no one else found |
| **GMAP** *(contrast — NOT in panel, no win-rate)* | Accurate **annotation-free** de novo splice sites (sandwich DP joint base+splice optimum); strong **SNP/indel + divergence tolerance**; de novo microexon; **sandwich-DP precedent** for RECTIFY's sequence-first refiner | **Disqualifying speed** (≥30× slower than minimap2; **single-threaded per-alignment DP**, read-parallel only); **no poly-A / ONT-homopolymer / 3'-end awareness**; **per-read isolation** (no cross-read consensus) | n/a — not in panel; investigated as algorithmic contrast / candidate (low expected win-rate if ever added) |

---

## 7. Failure Modes at the 3' End / CPA / Splice Boundary (the crux)

The ONT error model concentrates error at **exactly** the 3' boundary: homopolymer deletion is the dominant
residual error (even R10.4.1), and it lands in the A-run immediately 5' of the CPA. Each aligner errs
differently there — the failure modes are **complementary**, which is the entire justification for the
ensemble + correct-first selection.

### Where each aligner errs (and why)
- **minimap2 — over-extends / soft-clips poly-A.** `--end-bonus 0` gives no incentive to extend cleanly
  through an A-run; it soft-clips trailing homopolymer bases and `reference_end` drifts. HP indel noise is
  distributed arbitrarily within the run (low A=1/B=2 + affine gaps don't pin a unique column). Intron `N`
  boundaries slide ±1-few bp when exon-end+intron-start resembles intron-end+exon2-start (common at HP-rich
  junctions). *(minimap2 §6)*
- **deSALT — slips the 3' into a genomic A-run.** No poly-A awareness: it aligns the A-tail as ordinary
  sequence, so the banded DP can extend/slip the reported CPA into a downstream genomic A-run. *But its
  exon/junction structure is correct* — so RECTIFY's tail correction (2B/2E/2G) on top of deSALT's clean
  backbone is the dominant winning combination. *(deSALT §7)*
- **mapPacBio — column-optimal but HP-ambiguous & splice-blind.** Its full DP gives the best raw boundary,
  but long-read HP indel noise still admits multiple near-equal DP placements; and with no GT-AG model it
  can mis-place a junction it found only as a scored gap. *(mapPacBio §"Weaknesses")*
- **uLTRA — no 3' end at all.** It *strips* poly-A before aligning; the 3' end is just the trimmed
  alignment end. It offers a clean junction backbone but **zero CPA precision**. *(uLTRA §"Weaknesses")*
- **gapmm2 — inherits minimap2's 3'-end drift** for the body, fixing only orphan terminal exons via edlib;
  the cs-overrun bug occasionally over-runs the query span at the terminus. *(gapmm2 §5.2, §7)*
- **GMAP — splice-strong but 3'-end-blind (external contrast, not in panel).** Despite the strongest
  *internal*-junction placement of any aligner surveyed (sandwich DP, joint base+splice optimum), GMAP offers
  **NO 3'-end / CPA benefit to RECTIFY**: it has no poly-A model and aligns the in-read A-run as ordinary
  sequence, so on ONT DRS the DP slips the unspliced terminal-exon end into a downstream genomic A-run exactly
  as deSALT/minimap2 do. Its splice strength does nothing for the terminal exon where RECTIFY's value lies —
  it would be just another upstream aligner whose 3' end RECTIFY's correct-first logic must clean. *(gmap.md §6, §8.2)*

### Why ensemble + correct-first recovers truth
1. **Aligner choice is read-dependent.** No single aligner is uniformly best at the noisy 3' terminal exon.
   The failure modes are *orthogonal*: minimap2 over-extends/soft-clips; deSALT/uLTRA recover small/terminal
   exons but slip/strip the tail; mapPacBio nails the boundary but misses splices. Running all five surfaces
   *a* correct 3' end for almost every read somewhere in the panel. *(ont_drs §4.3, INFER)*
2. **Correct-FIRST selection.** RECTIFY corrects each aligner's 3' end *before* choosing a winner — using HP-
   aware modules: **2E** A-tract walk-back (uses the **reference** genome, length-independent → robust to HP
   miscounting), **2C** indel correction (needs MD tags), **2G** soft-clip rescue at homopolymer boundaries
   (extends 3' outward, takes priority over 2E to avoid cancelling corrections), **2H** junction refinement
   (rescores every `N`-op with HP-anchored semi-global DP + 0.5 canonical-HP prior). Winners are chosen on
   **corrected** positions, so the noisy raw alignment never decides — the aligner that produced the best
   *correctable* structure wins. This is why deSALT's homogeneous junctions (→ corrected CPA at a clean exon
   edge) and mapPacBio's column-optimal boundary dominate the corrected-3' metric. *(CLAUDE.md, deSALT §3)*
3. **Module 2H feeds on the ensemble.** `--aligner-bams` passes all 5 per-aligner BAMs to every correction
   run; the candidate junction pool is the union across aligners. A junction one aligner placed precisely
   rescues the others — cross-aligner junction recovery that no single aligner could do.
4. **Orthogonal validation closes the loop.** Corrected CPAs are refined/validated against NET-seq tables
   (the dedicated-3'-assay gold standard) — the field's best practice for 3'-end accuracy. *(ont_drs §5.1)*

**Net.** The win-rate spread (79/18/2/0.8/0.1) is the numerical shadow of complementary failure modes:
deSALT supplies the best *junction structure* (→ best corrected CPA after tail correction), mapPacBio the
best *terminal boundary*, and the rest cover residual independent failures. Correct-first selection on HP-
corrected positions is the mechanism that turns five individually-imperfect 3'-end estimates into one
bp-accurate per-read CPA.

---

## 8. Flagged for Adversarial Review

1. **⚠ deSALT runs no `-x` preset** (verified L1560: `null` ~13% model on ONT DRS data). Possibly
   sub-optimal vs `-x ont2d`; never A/B-tested. Most concrete "is this intentional?" item.
2. **⚠ gapmm2 `-i 5000`** likely sizes only gapmm2's edlib terminal window; the minimap2 sub-pass may use
   the **splice default 200k max-intron**, diverging from the other aligners' explicit `-G 5000`.
   Reconstructed from README + RECTIFY code, not a wheel read.
3. **⚠ AT>CG deletion rate** (CLAUDE.md table) vs "A/T basecalled more accurately overall" (literature) —
   reconciled by *inference* (deletion-specific vs aggregate), not a cited measurement.
4. **⚠ mapPacBio "chaining"** is multi-k-mer voting, not formal SDP — the seed-chain-align taxonomy is
   loose for BBMap.
5. **⚠ Win rates are single-dataset** and a pre-v3.3.0 `index_col` merge bug corrupted minimap2/gapmm2/
   uLTRA read_ids (they "never properly competed"). Confirm the cited rates are **post-fix**.
6. **⚠ deSALT "duplicate primary alignments" + mapPacBio "unmapped duplicates"** are both observed RECTIFY
   behaviors patched post-hoc; worth confirming they are aligner bugs vs invocation artifacts.

---

*All FACT claims trace to the seven `01_investigation/` dossiers and `multi_aligner.py` (invocations
re-verified). INFER claims — especially the win-rate causal chains and the 3'-end failure-mode
complementarity — are synthesis, not assertions by the tool authors.*
