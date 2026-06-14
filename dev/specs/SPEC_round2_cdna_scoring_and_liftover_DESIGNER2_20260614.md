# SPEC — Round-2 cDNA Alignment: Scoring, cDNA→Genome Lift-over, Integration, Validation

**Author:** Designer 2 (2-designer + 1-adversary scoping panel)
**Date:** 2026-06-14
**Scope of THIS doc:** the *scoring scheme* (genome-vs-cDNA, fairly), the *cDNA→genome
lift-over* (N-op insertion), *integration* with the existing consensus + anchor gate,
and *acceptance/validation/falsification*. The isoform-library construction and the
read→cDNA alignment step are Designer 1's slice; this doc consumes Designer 1's
per-cDNA **transcript→genome block map** and the per-read cDNA alignment record.

**Status:** design only — no production code. Concrete enough to build from.

---

## 0. The one-sentence thesis

> **Don't invent a new "fair" genome-vs-cDNA metric. RECTIFY's HP-aware edit distance
> (HP-ED) is *already* gap-neutral (N-ops cost 0). The fair comparison is HP-ED computed
> on the cDNA alignment's *lifted genome CIGAR* — after lift-over, introns are free on
> both candidates, so a cDNA can only win by *strictly lowering* exonic error, and its
> transferred N-ops are policed by the very same perfect-match anchor gate that fixed the
> mapPacBio gaming.**

Everything below is the consequence of that thesis.

---

## 1. THE SCORING SCHEME — genome alignment vs cDNA alignment, fairly

### 1.1 Why a naive "fewer gaps wins" metric is the trap (and why we don't build it)

A cDNA alignment is contiguous by construction (introns pre-removed). Any metric that
credits "fewer N-ops / fewer gaps / higher contiguous identity" makes the cDNA win every
read trivially — the *exact* gaming the anchor gate just fixed (`_add_chimera_flag`,
mapPacBio relabeling 3' soft-clips as free introns), run in reverse. We refuse to build
such a metric.

### 1.2 The actual metric: HP-ED on the *lifted* CIGAR

RECTIFY already ranks aligners per read by:

```
sort key (use_hp_ed): ['read_id', '_effective_chimera_ok', 'hp_edit_distance', '_span']
ascending            : [   —     ,         True            ,       True       ,  False ]
```
(`corrected_consensus.py` lines ~1378-1381). HP-ED costs (`_cigar_hp_edit_distance`,
lines 59-145):

| CIGAR op | cost |
| --- | --- |
| `=` match | 0 |
| `X` / `M`-mismatch | 1.0 |
| `D` deletion | HP-aware `del_cost(hp, base)` (fallback 1.0) |
| `I` insertion | HP-aware `ins_cost(hp, base)` (fallback 1.25) |
| **`N` intron** | **0 — free pass** |
| `S` soft-clip | 1.0/base |
| `H` hard-clip | 1.0/base |

**The cDNA candidate is treated as just another aligner.** We:

1. Take the read's alignment to the padded cDNA (Designer 1).
2. **Lift it to a genome CIGAR** (§2) — this *inserts* the cDNA's exon-junction N-ops.
3. Compute `hp_edit_distance`, `aligned_bases` (`_cigar_aligned_bases`),
   `min_junction_anchor` (`_cigar_min_junction_anchor`) on that lifted record —
   **the identical functions** the genome aligners go through in `_read_hp_edit_distances`.
4. Let it compete in `merge_corrected_tsvs` under the existing sort key.

After lift-over, **both** candidates have free N-ops. The *only* things that can differ are:

- (a) exonic bases the genome aligner soft-clipped / mismatched at a junction it failed to
  span, which the cDNA recovers as `=` (HP-ED drops); and
- (b) the junction anchors — already policed by `_chimera_ok`.

That is apples-to-apples by construction, with **zero new scoring code** — only the
lift-over (§2) is new.

### 1.3 The trivial-win defeat (REQUIRED trace — put in the doc/tests)

**Case: a read the genome aligner *already splices correctly*.**

- Genome winner: `120= 87N 95=`, all exonic bases match → HP-ED = 0, one N (free),
  left anchor 120, right anchor 95.
- The same read maps to the cDNA spanning the same two exons; lift-over reinserts the
  **same** 87N at the **same** boundary with the **same** flanking `=` runs → lifted
  HP-ED = 0, same anchor, same span.
- Sort key ties on `_chimera_ok` (both 0), `hp_edit_distance` (both 0), `_span` (equal).
  **The cDNA does NOT displace the genome winner.** No tie-break grants it a win.

A cDNA candidate can win **only when it strictly lowers HP-ED** — i.e. it genuinely
resolved exonic sequence the genome round could not. This trace *is* the apples-to-apples
argument; ship it as a unit test.

### 1.4 The explicit decision rule (the "win-guard") — do NOT rely on sort stability alone

HP-ED is a *sum*, so a shorter alignment can have a lower absolute HP-ED without being
better. Sort-key ties are necessary but not sufficient protection. State the eligibility
rule explicitly:

> A cDNA-derived candidate may replace the genome winner for a read **IFF**, after lift-over:
> 1. it passes the anchor gate: `_effective_chimera_ok == 0` (its lifted N-ops cleared
>    `_add_chimera_flag` exactly as a genome aligner's would), **AND**
> 2. **strict** improvement: `hp_ed_cdna ≤ hp_ed_genome_winner − ε` (ε = 1.0 — must beat
>    by at least one full mismatch-equivalent; ties keep genome), **AND**
> 3. **no-shrink guard:** `aligned_bases_cdna ≥ aligned_bases_genome_winner` (it did not
>    "win" by aligning *fewer* exonic bases).

Tie, or fewer aligned bases, or failed gate → **keep genome.** Condition (3) is the guard
the §1.3 reframe alone does not cover; it is mandatory.

Implementation note: (1) is automatic once the cDNA candidate is a row in `rep_df` (the
sort puts `_effective_chimera_ok` ahead of HP-ED). (2)+(3) are enforced by a small
**post-sort veto** restricted to reads whose winner is the cDNA candidate: re-fetch the
best *genome* candidate for that read and revert the winner to it unless (2)∧(3) hold.
Equivalently, encode (2) as an `ε`-quantized HP-ED column and (3) as a sort tier — but the
explicit veto is clearer and auditable. Prefer the veto.

### 1.5 Worked examples — chosen so existing rescue does NOT already cover them

The genome round already has `_rescue_5prime_softclip` (a plain 5' exon-1 soft-clip is
*already* converted to matches → HP-ED ≈ 0 without any cDNA round). So the win examples
below are cases the genome path **cannot** fix on its own.

**WIN #1 — micro-exon the splice-aware genome aligner could not seed.**
Read crosses exon A — 9 bp micro-exon — exon B. Genome aligners (seed-and-extend) cannot
seed a 9 bp exon; they emit `…M (mismatches over the micro-exon region) S` or skip it with
a single fused intron that swallows the micro-exon → those 9 bp + flank become `X`/`S`
(HP-ED ≈ 12). The cDNA contains exon A+micro+B contiguously; the read maps clean. Lift-over
emits `…= 5xxN 9= 3yyN =…` (two N-ops). Both N-ops have flanking `=` anchors ≥ K → gate
passes. Lifted HP-ED ≈ 0. (2) `0 ≤ 12−1` ✓, (3) aligned_bases higher (recovered 9 bp) ✓
→ **cDNA wins, correctly.** Attribution: "placed micro-exon."

**WIN #2 — multi-junction read whose genome seed fragmented.**
A 4-exon read; one genome aligner spans 3 junctions but the third junction's exon-4 anchor
is 4 bp and noisy → either `_chimera_ok=1` (flagged) or carries `X`s near the junction
(`_count_junction_proximity_errors` penalty in `score_alignment`). The cDNA (the exact
4-exon isoform) maps the whole read clean; lift-over emits 3 N-ops each with ≥K anchors.
Gate passes, HP-ED strictly lower, aligned_bases ≥ → **cDNA wins.** Attribution: "spanned
multi-junction."

**LOSE #1 — intron-retention isoform (THE protective example).**
Read retains intron 2 (a real biological IR transcript, or a mutant). The library's spliced
cDNA does **not** contain that intron. Forcing the read onto the spliced cDNA soft-clips /
mismatches the retained-intron bases (they have nowhere to go in the cDNA) → lifted cDNA
HP-ED **higher** than the genome alignment, which represents the retained intron in-line as
`M`/`=` across the genomic interval. (2) fails (cDNA not strictly better) → **genome wins.**
This doubles as proof the metric protects novelty: a transcript the library lacks is left to
the genome round, not force-fit.

**LOSE #2 — fully-spliced, both correct (the §1.3 trace).** Tie → genome kept. No win.

**LOSE #3 — cDNA wins by aligning *fewer* bases.** A read maps to only the 3' half of a
cDNA with HP-ED 2, vs a genome alignment over the full read with HP-ED 3. Absolute HP-ED
favors cDNA, but `aligned_bases_cdna < aligned_bases_genome` → (3) vetoes → **genome kept.**

### 1.6 Shaky assumptions (flagged for the adversary)

- **A1.** "Lifting then HP-ED is sufficient" assumes the lift-over is *faithful* — a buggy
  lift that drops boundary indels or mis-attributes I-bases could fabricate a lower HP-ED.
  Mitigation: §2 edge-case table + a round-trip invariant test (lifted CIGAR's
  query-consuming length == read length; ref-consuming span == sum of block lengths
  touched).
- **A2.** ε = 1.0 and the no-shrink guard are heuristics; the right ε may be HP-ED-scale
  dependent. Falsification (§4) measures win attribution to tune, not guess.
- **A3.** HP-ED treats `N` as free on BOTH sides — correct here, but it means the cDNA
  round cannot be *penalized* for inserting an extra (real) junction the genome missed;
  the *only* thing stopping a bad inserted junction is the anchor gate. Hence §3's hard
  precondition.

---

## 2. cDNA→GENOME LIFT-OVER (N-op insertion)

### 2.1 Inputs (from Designer 1)

- **Read→cDNA alignment** (a pysam record vs the padded cDNA reference): cDNA-local CIGAR,
  `cdna_ref_start` (0-based offset into the padded cDNA), read query sequence, strand of
  the *mapping* to the cDNA.
- **Block map for that cDNA** — ordered list of blocks, each:
  ```
  {transcript_coord:int,   # 0-based start in the PADDED cDNA sequence
   genome_chrom:str, genome_coord:int,  # 0-based genome start of the block
   length:int,
   strand:'+'|'-',         # gene/RNA strand
   kind:'pad5'|'exon'|'pad3'}
  ```
  Blocks tile the padded cDNA contiguously in transcript coordinates. **Adjacent `exon`
  blocks whose genome coords are non-contiguous define the junctions** (the N-ops to
  insert). `pad5`/`pad3` blocks are genomic flank — contiguous with the first/last exon in
  genome space (no N between a pad and its abutting exon).

### 2.2 Output

A genome BAM record: `reference_name = genome_chrom`, `reference_start`, a genome CIGAR
(`=`/`X`/`M`/`I`/`D` within blocks, `N` at spanned exon→exon boundaries, `S` for read ends
beyond the padded reference), and the `is_reverse` flag matching gene strand. Then it goes
through §1.2 step 3 (HP-ED/anchor/span) like any aligner.

### 2.3 Algorithm (step by step)

Work in **transcript coordinate space**, then translate each maximal "same-block" run to
genome space, emitting N at every block boundary the read *spans*.

```
1. Normalize strand. ASSUMPTION (confirm with Designer 1): the padded cDNA is stored
   in mRNA-SENSE — i.e. for a '-' gene the cDNA reference is revcomp(genome). This is
   what the 5'/3' pad language already implies (pad5 = upstream of TSS in transcript).
   A read maps FORWARD to this mRNA-sense cDNA, i.e. in TRANSCRIPT orientation — NOT
   reference orientation. (Caution: do NOT reuse the `_rescue_5prime_softclip`
   "no revcomp needed" comment here — that comment is about a read ALREADY in a
   genome-aligned BAM, hence already reference-forward. A read aligned to an mRNA-sense
   cDNA is a different category and IS transcript-oriented.)

   To emit a valid genome BAM record for a '-' gene transcript you MUST:
     (i)   reverse-complement the read query_sequence (so it is reference-forward),
     (ii)  emit CIGAR ops in REVERSE order (transcript pos 0 maps to the HIGHEST
           genome coord; the genome CIGAR must read low→high coordinate),
     (iii) set is_reverse = True,
     (iv)  iterate blocks genome-ASCENDING when building the CIGAR.
   For a '+' gene: query as-is, CIGAR in transcript order, is_reverse = False,
   blocks already genome-ascending. Steps 2–4 below operate on the (possibly
   revcomp'd, possibly reversed) read+CIGAR so they are strand-agnostic thereafter.

2. Walk the cDNA-local CIGAR op by op, maintaining (t_pos = current transcript
   coord, q_pos = read query pos). For each op, intersect its transcript span
   [t_pos, t_pos+len) with the block tiling:

   a. Leading/trailing S in the cDNA CIGAR (read overhangs the padded cDNA ends):
      emit S directly (read ran off the pad — true novel TSS/CPA beyond the flank).

   b. For M/=/X/D/I (consume transcript and/or query):
      - Find which block t_pos falls in. Map t_pos → genome:
          genome_pos = block.genome_coord + (t_pos - block.transcript_coord)   ['+' gene]
          (mirror for '-' gene, decreasing)
      - Emit the op against the genome, but SPLIT the op at block boundaries:
        if a single M/=/X run crosses from block_i (exon) into block_{i+1} (exon)
        and their genome coords are non-contiguous, emit
            <left part as =/X/M>  N(gap = genome gap between blocks)  <right part>
        The matched bases on each side of the N ARE the perfect-match anchor that
        _cigar_min_junction_anchor will read. This split is the CORE step.

   c. pad→exon and exon→pad boundaries are genome-CONTIGUOUS → NO N; the op simply
      continues across the boundary as ordinary M/=/X/D/I.

3. reference_start = genome coord of the first reference-consuming op (after leading S).

4. Merge adjacent identical ops (= = -> =) but NEVER merge across an N.
```

### 2.4 Edge-case table

| # | Case | Handling |
| --- | --- | --- |
| E1 | Read spans 0 junctions (within one exon, or entirely in a pad) | No N. Pure `=`/`X`/`M`(+I/D). Contiguous — the trivial genome-equivalent; usually loses on tie. |
| E2 | Read spans 1 junction | One N at the spanned exon→exon boundary; `=` runs on both flanks = anchors. |
| E3 | Read spans many junctions | One N per spanned boundary; each gets independently anchor-gated (weakest flank of weakest N sets `min_junction_anchor`). |
| E4 | **M-op split exactly at a boundary** (common) | Split the run: `M_left, N, M_right`. The two halves carry the matched bases the gate measures. This is E2/E3's mechanics. |
| E5 | Deletion (`D`) abutting a boundary | **Convention: D belongs to the EXON side, never merged into the N.** A D immediately 5' of a junction stays as `…=, D, =, N, …`; a D immediately 3' stays `…N, =, D, …`. Never `…D N…` collapsed into one larger N (that would fabricate a junction shift and corrupt the anchor count). |
| E6 | Insertion (`I`) abutting a boundary | I consumes query only. **Convention: I attaches to the 5' (upstream) exon side** of the boundary: emit `…=, I, N, =…` (the inserted read bases are placed at the donor exon end, never inside the intron). Deterministic; matches how genome aligners place boundary insertions. |
| E7 | Read 5' end falls in `pad5` (novel/variant TSS upstream of annotated start) | The pad is genomic → aligns as `=`/`X` contiguously into exon 1. **This is the feature**: TSS variation is captured as ordinary matched bases, NOT soft-clipped. Soft-clip only if the read extends *beyond* the pad5 block. |
| E8 | Read 3' end falls in `pad3` (CPA variation / poly-A) | Same: pad3 aligns contiguously; the 3' end lands at the true position. Downstream `find_polya_boundary`/walkback then assigns the CPA exactly as for any genome alignment — **lift-over does NOT do CPA correction**, it only produces the raw genome record. |
| E9 | Read extends beyond pad (longer than padded cDNA) | Emit `S` for the overhang (E2.a). If this is frequent, Designer 1's pad is too short — surface as a metric, not a silent clip. |
| E10 | **Novel junction**: read maps to an exon combination present in the library cDNA, but whose *specific* junction the genome round never called | Lift-over emits the N **unconditionally** (it is implied by the block map). It then **must clear the same anchor gate** — and as a single-read novel junction it gets **no cross-read support relaxation** (`_compute_junction_stats` needs ≥`_MIN_JUNCTION_SUPPORT` observations), so it must satisfy the raw `≥ min_junction_anchor_bp` perfect-match anchor. **This is correct, not a bug.** A genuinely novel junction with strong flanking matches passes; a smuggled bad one is rejected exactly like a bad genome N. |
| E11 | `-` gene strand | **revcomp the read query AND reverse the CIGAR op order** (transcript pos 0 → highest genome coord), iterate blocks genome-ascending, set `is_reverse=True` (§2.3 step 1). Do NOT skip the revcomp — the source read is in transcript orientation (vs an mRNA-sense cDNA), not reference orientation. |
| E12 | cDNA-local alignment has its own N (shouldn't, but defensively) | A no-splice cDNA aligner should never emit N; if one appears, treat as a hard error / drop the candidate (it violates the contiguity premise). |

### 2.5 Lift-over invariants (assert in code / tests)

- **Query-length invariant:** sum of query-consuming ops (`M/=/X/I/S`) in the lifted CIGAR
  == read length. (Catches dropped/duplicated bases — defends A1.)
- **Reference-span invariant:** `reference_end − reference_start` == sum of ref-consuming
  spans over touched blocks + sum of inserted N gaps. (Catches mis-placed N gaps.)
- **No-N-merge invariant:** no two `=`/`M` runs were merged across an N.
- **Anchor sanity:** `_cigar_min_junction_anchor(lifted) ≥ 0` and finite when N present.
- **Orientation invariant (CRITICAL — catches the §2.3-step-1 strand bug):** the FIRST
  matched read base of the lifted record must equal `genome[reference_start]`, and the
  LAST matched read base must equal `genome[reference_end−1]` (allowing for genuine
  mismatch positions — test on a known-clean read). The query-length and ref-span
  invariants above are **revcomp-invariant** and CANNOT catch a mis-oriented minus-strand
  record; this base-identity check is the one that does. A mis-oriented record would feed a
  garbage/deflated HP-ED into selection and pass every other invariant.

---

## 3. INTEGRATION with the existing consensus + anchor gate

### 3.1 The choice: a candidate in `merge_corrected_tsvs`, NOT a post-consensus rescue

**Decision: lift each cDNA win to a per-aligner corrected BAM under a synthetic aligner
name (e.g. `cdna_round2`) and feed it into the *existing* `merge_corrected_tsvs`
winner-selection as just another candidate.** Argument:

1. **Reuse = the entire safety story.** As a `rep_df` row, the cDNA candidate
   automatically gets `_chimera_ok` (the anchor gate runs on its lifted N-ops via
   `_add_chimera_flag` → `_cigar_min_junction_anchor`), `hp_edit_distance`, `aligned_bases`,
   `_span`, and the minimal-anchor hard-filter. A post-consensus rescue would re-implement
   each of these in a parallel path that *will* drift out of sync — and the gate is the one
   thing standing between us and reverse-gaming. We refuse to duplicate it.
2. **One winner-selection invariant.** The §1.4 win-guard (strict-improve + no-shrink) is a
   thin post-sort veto on top of the existing sort, not a second selection engine.
3. **Provenance/audit for free.** The cDNA candidate appears in the same summary TSV
   (`_is_winner`, `_chimera_ok`, `hp_edit_distance`) so per-read win attribution (§4) is a
   straight read of existing columns.

Post-consensus rescue is rejected: it would either bypass the gate (unsafe) or re-implement
it (drift risk), and it loses the unified audit trail.

### 3.2 Anchor gate MUST be ON whenever Round 2 is active — HARD PRECONDITION

`_MIN_JUNCTION_ANCHOR = 0` by default (gate OFF; the yeast mapPacBio panel relies on it
being off). **The cDNA round's entire safety depends on the gate being ON.** Therefore:

> **Precondition (enforced at config-load):** if the cDNA Round-2 feature is enabled, the
> pipeline MUST set `min_junction_anchor_bp > 0` (recommended 10, matching the human-DRS
> setting). If Round 2 is enabled with the gate at 0, **fail fast** with a clear error —
> never run Round 2 with the gate off (on yeast that would let a cDNA transfer smuggle a
> bad junction past selection).

This is the single most important integration constraint. Surface it loudly.

### 3.3 Where Round 2 sits in the pipeline

```
Round 1: per-aligner align → per-aligner rectify correct (5' rescue, walkback) → corrected BAMs
   │
   ├─(NEW) Round-2 entry gate (per read; §3.4): is the Round-1 result weak?
   │        if NO  → skip Round 2 for this read (clean winner can't be beaten by a tie)
   │        if YES → align read to cDNA library (Designer 1) → lift-over (§2)
   │                 → emit cdna_round2 corrected record
   │
   └─ merge_corrected_tsvs / corrected_consensus  ← genome aligners + cdna_round2 compete
           under existing sort key + anchor gate + §1.4 win-guard veto
   │
   └─ Stage 2/3 aggregate/analyze (CPA, bedgraphs) — unchanged; consumes the consensus BAM
```

Round 2 runs **after Round-1 per-aligner correction**, **before** consensus selection
(so its output is a peer candidate), and is **gated on Round-1 weakness** for efficiency.

### 3.4 Entry-gating on Round-1 weakness (perf) — and its one limitation

Entry-gating on the **clip + anchor** conditions is correctness-neutral: a Round-1 winner
with no effective clips AND clean anchors that *also* has HP-ED 0 can at best *tie* a cDNA
candidate, and §1.3/§1.4 keep the genome on ties. So the safe (recall-preserving) entry gate
skips Round 2 only for reads whose Round-1 consensus winner has **all** of:
- no 5'/3' effective soft-clip (`effective_five_prime_clip == 0` and `…three… == 0`),
- `min_junction_anchor ≥ min_junction_anchor_bp` (or no junctions),
- `hp_edit_distance == 0` (τ = 0).

**Caveat on τ > 0:** raising the HP-ED entry threshold to τ > 0 (e.g. skip reads with HP-ED
≤ 2) is NOT correctness-neutral — a Round-1 winner at HP-ED 2 (two mismatches) *can* be
strictly beaten by a cDNA at HP-ED 0, and τ=2 would silently forgo that legitimate win. So
τ > 0 trades **recall for compute**, it does not just save compute. **Default τ = 0**
(gate only on clips + anchors, which genuinely cannot lose on a tie); expose τ as an opt-in
compute-vs-recall knob, documented as such — do not set τ=2 by default.

Reads failing any condition → eligible for Round 2.

**Limitation (hand to adversary):** a read whose Round-1 winner *looks* clean but is on the
**wrong paralog** (correct-looking CIGAR against the wrong locus) never gets the cDNA
chance under entry-gating. Two options: (a) accept it (paralog mis-assignment is a separate
problem the cross-chrom `_n_agree` tier already partly addresses); (b) offer a
`--round2-all-reads` mode that disables entry-gating for paralog-sensitive runs at a compute
cost. Recommend (a) default, (b) opt-in.

---

## 4. ACCEPTANCE / VALIDATION / FALSIFICATION

The adversary's #1 attack: **the cDNA library is annotation-derived → annotation bias,
which conflicts head-on with RECTIFY's motif-agnostic, novel-junction-for-mutants
philosophy.** This is load-bearing, not a footnote. Two non-negotiable rules:

- **R1 — annotated-junction fraction is a VALIDATION LABEL ONLY, never a selection input.**
  It must NOT enter HP-ED, the gate, or the win-guard (mirrors the existing rule that
  `score_alignment` does not score annotated-junction matches — see its docstring,
  lines ~750-752). Using it in selection would reintroduce exactly the bias the genome
  round avoids.
- **R2 — Round 2 must not suppress genuine novel junctions the genome round found.** A
  mutant read with a novel exon combination absent from the library must be *unharmed*:
  the cDNA loses (LOSE #1 mechanics), the genome wins, the novel junction survives into the
  consensus. This is a positive test, not an aspiration.

### 4.1 Metrics (net-win test)

| Metric | Tool / source | Net-win direction |
| --- | --- | --- |
| **Per-read win attribution** | summary TSV (`_is_winner`, `_chimera_ok`, `hp_edit_distance`, aligned_bases) + lift-over reason tag | EVERY `cdna_round2` win must carry an attributable cause: *recovered-clip*, *placed-micro-exon*, or *spanned-multi-junction*. Compute the fraction of cDNA wins with **no** attributable cause. |
| **CPA proximity** | PolyASite atlas (human) / lab CPA clusters (yeast); distance of corrected 3' end to nearest annotated CPA | cDNA-won reads should be **no worse** (ideally better) than the genome winner they replaced. |
| **Annotated-junction fraction** (LABEL ONLY, R1) | reference annotation (GENCODE / SGD) | Among cDNA wins, junction calls should be ≥ genome baseline — but this is a *report*, never a gate. |
| **Novel-junction preservation** (R2) | spike-in / mutant reads with library-absent exon combos | 100% of such reads keep the genome winner; 0% force-fit onto a library isoform. |
| **Net HP-ED reduction** | consensus HP-ED with vs without Round 2 | Total corrected-read HP-ED should drop; if it doesn't, Round 2 adds nothing. |

### 4.2 Falsification plan — when do we conclude the cDNA round does NOT help?

Round 2 is **rejected** (turned off by default) if ANY of:

1. **Trivial-win leak.** A meaningful fraction (> ~5%) of cDNA wins have **no attributable
   cause** (§4.1 row 1). That means the metric is letting contiguity-by-construction win —
   the reverse-gaming we set out to prevent. (Primary falsifier.)
2. **CPA regression.** cDNA-won reads land *further* from annotated CPAs than the genome
   winners they displaced (median distance worse, signif.). Means lift-over/3'-handling is
   degrading the very 3'-end accuracy RECTIFY exists for.
3. **Novelty suppression (R2 fails).** Any library-absent novel-junction read gets force-fit
   onto a library isoform (cDNA wrongly wins). Even a small rate is disqualifying — it
   violates the core philosophy.
4. **No net HP-ED reduction.** Consensus HP-ED with Round 2 ≈ without. Then Round 2 is pure
   compute cost with no benefit → off.
5. **Gate dependence proven harmful.** If, with the gate correctly ON, cDNA wins still carry
   bad junctions (annotated-junction fraction of cDNA wins *below* genome baseline), the
   anchor gate is not sufficient protection for transferred N-ops → halt and redesign.

If (1)–(5) all pass on a held-out sample (human DRS A549 chr5 + a yeast mutant with known
novel splicing), Round 2 ships **on by default for the protocols tested**, gate forced ON.

### 4.3 Provenance

Per the lab standard, every Round-2 corrected BAM/TSV gets a `PROVENANCE.json` sidecar:
library version (Designer 1), block-map build, `min_junction_anchor_bp`, ε / no-shrink
guard values, entry-gate τ, and the commit. The cDNA candidate's rows in the summary TSV
carry an `aligner = cdna_round2` tag so any downstream audit can isolate its contributions.

---

## 5. Open questions / handoff to the BUILD agent and adversary

- **For Designer 1:** confirm the block-map field names (§2.1) and whether pad blocks are
  emitted explicitly or implied by exon-block extents; confirm `-` gene block ordering.
- **For the adversary:** the three sharpest assumptions are A1 (lift faithfulness — mitigated
  by §2.5 invariants), the annotation-bias philosophy conflict (§4 R1/R2 — the central
  defense), and the gate-OFF-on-yeast hazard (§3.2 hard precondition). The entry-gating
  paralog blind spot (§3.4) is a known, accepted limitation with an opt-in escape hatch.
- **ε and τ** are heuristics (A2); §4.1 win-attribution data tunes them rather than guessing.
```
