# Red-Team Audit — De-Novo / Tier-1 Aligner Dossiers

**Auditor role:** skeptical adversarial fact-checker.
**Targets:** `minimap2.md`, `deSALT.md`, `gapmm2.md`, `mapPacBio_bbmap.md`
(in `dev/aligner_investigation/01_investigation/`).
**Date:** 2026-06-19.

**Method.** Each algorithmically-substantive claim was classified FACT /
PLAUSIBLE INFERENCE / SUSPECT, then high-risk claims were verified against
primary sources: minimap2 `options.c` + README (lh3/minimap2), the deSALT Genome
Biology 2019 paper + README (ydLiu-HIT/deSALT), gapmm2 README (nextgenusfs/gapmm2),
BBMap `mapPacBio.sh` + `BBMapGuide.txt` (BioInfoTools/BBMap mirror), and the RECTIFY
wrapper source `rectify/core/align/multi_aligner.py` (read directly).

**Bottom line up front.** All four dossiers are **unusually well-sourced and largely
accurate**. The external/primary-source facts I spot-checked verified almost
verbatim. The weak points are *not* fabricated facts but (a) a small number of
overstated "default" / "splice-aware" characterizations, (b) one internal
source-comment vs. dossier-rationale mismatch, and (c) the recurring causal leap from
a real mechanism to the **78.9% / 18.2% / … win-rate** numbers, which every dossier
labels INFERENCE but states with more confidence than the evidence supports. The win
rates themselves are **un-sourced beyond CLAUDE.md** — no dossier shows the underlying
measurement.

---

## 1. Audited-Claims Table

Verdict key: **FACT** (verified against primary source) · **FACT✓src** (verified in
RECTIFY source) · **PI** (plausible inference, reasonable) · **OVERSTATED** (true
kernel, too-strong wording) · **SUSPECT** (questionable / unverifiable) · **CONTRA**
(internal contradiction).

### minimap2.md

| # | Claim | Verdict | Evidence / source |
|---|-------|---------|-------------------|
| M1 | splice preset = `k=15, w=5, A=1, B=2, O=2,32, E=1,0, C=9, z=200,100, junc_bonus=9, junc_pen=5` | **FACT** | `options.c::mm_set_opt`: verified k=15,w=5,a=1,b=2,q=2,e=1,q2=32,e2=0,noncan=9,zdrop=200,zdrop_inv=100,junc_bonus=9,junc_pen=5. Exact match. |
| M2 | splice sets `MM_F_SPLICE_FOR | MM_F_SPLICE_REV` → default `-u b` (both strands) | **FACT** | `options.c`: flags = `MM_F_SPLICE | MM_F_SPLICE_FOR | MM_F_SPLICE_REV | MM_F_SPLICE_FLANK`. Confirms `-uf` is a RECTIFY override of the both-strand default. |
| M3 | splice:hq differs by `B=4, O=6,24, C=5` | **FACT** | `options.c` splice:hq: `noncan=5, b=4, q=6, q2=24`. Exact. |
| M4 | `-k14` is the documented official ONT-DRS value, not RECTIFY-invented | **FACT** | minimap2 README: `minimap2 -ax splice -uf -k14 ref.fa reads.fa # noisy Nanopore Direct RNA-seq`; "smaller k-mer … for first/last exons." Dossier correctly says k=15 is the preset default and -k14 is the recommended DRS tweak. |
| M5 | map-pb sets HPC (`MM_I_HPC`) + k=19; HPC hurts ONT | **FACT** (params) / **PI** (hurts-ONT) | `options.c`: map-pb `flag |= MM_I_HPC, k=19`. "HPC hurts Nanopore" is the documented cookbook rationale (widely stated by lh3) — accept as FACT-adjacent. |
| M6 | bw=200000; dossier preset string lists `bw=200k, max_gap=2000, max_gap_ref=200000` | **MOSTLY FACT / minor** | `options.c` shows `bw=200000, bw_long=200000, max_gap=2000`. The dossier's `max_gap_ref=200000` label most likely = `bw_long`; the *number* is right, the *name* `max_gap_ref` is not a literal splice-preset field. Cosmetic. |
| M7 | `--splice-flank=no` "hands intron-boundary/3'-end placement to RECTIFY" for 3'-end accuracy | **OVERSTATED / CONTRA(minor)** | The **source comment** at `multi_aligner.py:252` says `--splice-flank=no  # Disable for compatibility`, NOT "for 3' end accuracy." CLAUDE.md's minimap2 section *does* say "important for 3' end accuracy." So the rationale is asserted by CLAUDE.md but the live code comment says "compatibility." Dossier §2.3 builds an elaborate mechanistic story on the CLAUDE.md framing; flag that the code author's stated reason differs. Splice-flank semantics themselves (GTA/GTG→0 etc.) are FACT. |
| M8 | "minimap2 wins only 0.1%" drives the Weaknesses framing | **UNSOURCED** | Only source is CLAUDE.md. No measurement shown. The mechanistic weakness claims (poly-A soft-clip, HP indel jitter) are independently plausible, but the 0.1% number must be treated as a RECTIFY-internal assertion, not an established fact. |
| M9 | `--end-bonus 0` default → no incentive to extend 3' into poly-A → soft-clips tail | **PI (good)** | Labeled INFERENCE in dossier. end-bonus default is 0; the causal chain is reasonable and properly hedged. |
| M10 | Chaining recurrence, KSW2 two-piece affine, Z-drop, intron `x2` state, GT-AG donor/acceptor tables | **FACT** | Matches paper + `ksw2_exts2_sse.c` description; internally consistent. Not independently re-derived line-by-line but consistent with published algorithm. |

### deSALT.md

| # | Claim | Verdict | Evidence / source |
|---|-------|---------|-------------------|
| D1 | All default params (`-k22[20-28], -l15, -s5, -a8, -n50, -c30, -g2000, -I200000, -i20, -d10, -e5, -m1, -M2, -O2,32, -E1,0, -z400, -w500, -R9, -B655350`) | **FACT** | deSALT README options table — every value verified verbatim. |
| D2 | deSALT wraps deBGA; RdBG-index; `START_POS_REF=0` patch; ≤500 bp line width; index memory hg38≈35 GB / mouse≈31 GB / fly≈3.5 GB | **FACT** | README verified verbatim (incl. memory table and the START_POS_REF patch). |
| D3 | Read-type presets `ccs/clr/ont1d/ont2d/null`; RECTIFY runs no `-x` → `null` (~13% err) | **FACT** | README preset list verified; `run_desalt()` builds `['aln','-t','-f','-o']` with no `-x` (source confirms). |
| **D4** | **Cross-read pooling: exon inference "integrates ALL skeletons of ALL reads"** | **FACT-PAPER (verified — NOT an overstatement)** | Genome Biology 2019: *"All the alignment skeletons of all the reads are then integrated to comprehensively detect the exon regions."* The dossier's central mechanism claim is **accurate**, not inflated. |
| D5 | Paper markets "homogeneous" alignments | **FACT** | Paper: "accurate and homogeneous full-length alignments." Verified. |
| D6 | GTF-free de-novo exon inference (`-G` optional, RECTIFY omits it) | **FACT** | Paper + README; `run_desalt()` passes no `-G`. |
| **D7** | **"This is why deSALT wins 78.9% … one aligner takes ~4/5 of reads" causal chain** | **PI but OVERSTATED** | Labeled INFERENCE, but stated with near-certainty ("the decisive factor," "structural head start"). Real mechanism (D4/D5) ✔, but the leap to the *specific* 78.9% and the claim that homogeneity is THE cause **conflates a plausible mechanism with an un-decomposed outcome**. Win rate is CLAUDE.md-only; no per-mechanism attribution was measured. |
| D8 | Engineering fragilities (N-times duplicate primaries → dedup; gz misparse; `-G` yeast SIGSEGV; `-f` must be local disk; LD_LIBRARY_PATH strip) | **FACT✓src** | All present in `run_desalt()` / `_dedup_desalt_bam()` / docstring (source read). |
| D9 | "Dual affine gap = same two-cost model as minimap2/ksw2; intron = near-flat long gap" | **PI (good)** | Param shape (O=2,32/E=1,0) is FACT; "same engine as ksw2" is labeled INFERENCE. deSALT does use an SSE banded SW; reasonable. |
| D10 | RdBG "breadth of candidates … root of sensitivity" | **PI** | Labeled INFERENCE; the three index lookups are FACT-PAPER; the causal framing is reasonable but not asserted by authors. |

### gapmm2.md

| # | Claim | Verdict | Evidence / source |
|---|-------|---------|-------------------|
| G1 | Wraps minimap2/mappy + edlib terminal refinement; deps mappy/edlib/natsort | **FACT** | README: "wrapper for minimap2 … uses edlib to find the terminal exon positions." |
| G2 | Defaults `-t 3, -m/--min-mapq 1, -i/--max-intron 500` | **FACT** | README options verified verbatim. |
| G3 | 13 bp 5' exon worked example at 409609–409621; benchmark 6926 aln → 409 left / 63 right | **FACT** | README API output `{'n':6926,'low-mapq':0,'refine-left':409,'refine-right':63}` and coords verified verbatim. |
| G4 | v25.4.13: runs `minimap2` directly (not mappy) due to mappy segfaults; cmd `minimap2 -x splice --cs=short` | **FACT (cmd partially)** | README confirms the segfault note and `minimap2 -x splice --cs`. The exact `--cs=short` (vs `--cs`) wasn't visible in the README text I fetched — dossier says `--cs=short`; downgrade to **near-FACT**, verify the `short` qualifier against the installed wheel. |
| G5 | RECTIFY `-i 5000`, omits `-m` (argparse `type=int` bug on 25.4.5 wheels) | **FACT✓src** | `run_gapmm2` lines 1108–1120: `-i 5000`, `-m` commented out with exact rationale. Verified. |
| G6 | mode="HW", k=0, task="path", degenNuc; `maxlen`=`-i`; cs stitch `~gtNNNNag`/`~ctNNNNac` | **FACT (README/src excerpts)** | README + RECTIFY `_cs_long_to_cigar`. The dossier itself flags (header) that gapmm2 source is "sparse … reconstructed from README excerpts" — honest caveat. Accept, but these are README-level not line-by-line wheel reads. |
| G7 | Duplicate-UUID → `seq()` returns None for BOTH → `TypeError` at `align.py:883`, ~22k of 695k; cs-overrun ~0.02% | **FACT✓src (RECTIFY characterization)** | `_clean_fastq` + `_paf_to_bam` `cigar_qlen != len(fwd_seq)` skip verified in source and CLAUDE.md. The specific gapmm2 line number (`align.py:883`) and "~22k of 695k" come from RECTIFY's own diagnosis, not upstream gapmm2 docs — correctly attributed to RECTIFY. |
| G8 | "Win rate ≈ 0.8%"; junction-proximity penalty docks gapmm2/minimap2 | **UNSOURCED (rate) / FACT✓src (penalty exists)** | 0.8% = CLAUDE.md only. `scoring.py:745` penalty referenced — plausibly real (not independently re-read here). |
| G9 | "Body alignment IS minimap2 so where minimap2 is right gapmm2 ≈ minimap2" → low marginal value | **PI (good)** | Sound reasoning, labeled INFER. |

### mapPacBio_bbmap.md

| # | Claim | Verdict | Evidence / source |
|---|-------|---------|-------------------|
| B1 | mapPacBio = `align2.BBMapPacBio` (BBMap), NOT minimap2 map-pb | **FACT** | `mapPacBio.sh` invokes `align2.BBMapPacBio`. The dossier's strongly-worded disambiguation is correct and valuable. |
| B2 | mapPacBio defaults: `build=1 overwrite=true minratio=0.40 fastareadlen=6000 ambiguous=best minscaf=100 startpad=10000 stoppad=10000 midpad=6000` | **FACT** | `mapPacBio.sh` literal default-args line verified verbatim. |
| B3 | Affine weights +100 / −127 / −51; 6 bytes/base; usemodulo ≈3 bytes; "recommended for Nanopore"; "dominated by short indels" | **FACT** | BBMapGuide.txt verified verbatim for all five. |
| **B4** | **"Default k = 13"** | **OVERSTATED** | BBMapGuide.txt mentions `k=13` only as a **memory-estimation example** (`stats.sh … k=13`), not as the index default. BBMap's actual default index k *is* 13, but the cited guide passage does not assert it as the default — the `[FACT]` label leans on a passage that doesn't say what's claimed. Soften to "k=13 is the standard/typical build; guide uses it as the canonical example." |
| **B5** | **mapPacBio sets `intronlen=10` and an explicit `maxindel=max(200000, max_intron)`** | **FACT✓src — RESOLVED** | Source confirms: `run_map_pacbio` sets `intronlen=10` (`multi_aligner.py:749`) and `maxindel=max(200000, max_intron)` (`multi_aligner.py:754`), alongside `minratio=0.4, fastareadlen=100000, -Xmx32g`. The earlier "no `maxindel` cap on the long-read path" concern is fully resolved: mapPacBio does not rely on BBMap's soft ~16000 default — its `maxindel` is ≥200 kb and scales with `max_intron`, so yeast introns (<1 kb) and even long mammalian introns are safely within the searched range. The D→N threshold is `intronlen=10`, not 50. |
| B6 | "BBMap is splice-aware"; splices found as scored deletions then D→N reclassified via intronlen; no GT-AG model | **FACT/PI** | BBMap is marketed "splice-aware" (Bushnell title) and the D→N intronlen mechanism is documented. "No GT-AG model in the core scorer" is INFERENCE but consistent with BBMap docs (no splice-signal scoring described). Reasonable. |
| B7 | 6 kb ceiling / AssertionError ~6019 bp; RECTIFY patches fastareadlen=100000 + pre-splits at MAX_MPB_READ_LENGTH=6000 + stitches | **FACT✓src** | `run_map_pacbio` lines 513, 522–609; `MAX_MPB_READ_LENGTH=6000`. Verified. The "~6019 bp" exact threshold is RECTIFY-docstring-sourced, fine. |
| B8 | `pt:i:N` is Dorado's poly-A tag, not BBMap output; mapPacBio copies header→QNAME; RECTIFY strips `_pt:i:` | **FACT✓src** | pt strip at lines 578–598 verified; Dorado `--estimate-poly-a` pt tag is correct attribution. |
| B9 | Unmapped-duplicate records pass `not is_secondary and not is_supplementary`; need `is_unmapped` guard | **FACT (CLAUDE.md v3.1.8)** | RECTIFY-side, correctly attributed. The "[INFERENCE]" that this stems from BBMap writing all reads to one stream is reasonable. |
| B10 | minratio=0.40 mapPacBio default "vs ~0.56 short-read default" | **MINOR SUSPECT** | mapPacBio 0.40 is FACT. The "~0.56" short-read default is uncertain — web sources cite BBMap default minratio as **0.56**, **0.76**, and skimmer 0.40 inconsistently across versions. Don't over-anchor the comparison number; 0.40 for mapPacBio is the load-bearing fact and it's solid. |
| B11 | "Wins 18.2%, second-best … globally optimal DP boundary beats heuristic extension" | **PI / UNSOURCED rate** | 18.2% = CLAUDE.md only. Full-DP-vs-heuristic mechanism is plausible and hedged as INFERENCE, but again not measured per-mechanism. |

---

## 2. Corrections Required (concrete)

1. **minimap2 §2.3 / table (M7) — splice-flank rationale mismatch.** The live code
   comment is `--splice-flank=no  # Disable for compatibility` (`multi_aligner.py:252`),
   not "for 3' end accuracy." Either (a) note that the *code author's* stated reason is
   "compatibility" while CLAUDE.md asserts "3' end accuracy," or (b) reconcile them.
   Do not present the elaborate 3'-end mechanistic story as the established reason
   without flagging that the source comment disagrees.

2. **mapPacBio B4 — "Default k=13".** Reword. The BBMap guide uses k=13 only as a
   memory-calc example; change `[FACT — default k=13]` to
   `[FACT — k=13 is the standard build / canonical example; guide does not state it as a hard default]`.

3. **mapPacBio B5 — RESOLVED.** mapPacBio sets `intronlen=10` (`multi_aligner.py:749`) and
   an explicit `maxindel=max(200000, max_intron)` (`multi_aligner.py:754`). The cap is ≥200 kb
   and scales with `max_intron`, so it does not rely on BBMap's soft ~16000 default; yeast
   introns (<1 kb) are comfortably within the searched range. The dossier now states this
   directly. No further action.

4. **All four — win-rate numbers (M8, D7, G8, B11).** Add one sentence per dossier:
   "The 78.9% / 18.2% / 2% / 0.8% / 0.1% win rates are RECTIFY-internal measurements
   (CLAUDE.md); this dossier does not reproduce the underlying benchmark, and the
   per-mechanism causal attributions to those exact numbers are INFERENCE, not a
   decomposition of the measured result." This is the single most important calibration
   fix — the mechanisms are real; the *attribution of a specific percentage to a specific
   mechanism* is not demonstrated.

5. **deSALT D7 — soften the causal certainty.** Phrases like "the decisive factor" and
   "structural head start" overstate. The cross-read pooling mechanism (D4) is genuinely
   FACT-PAPER and well-verified; the inference that it *causes* 78.9% should read as "is
   consistent with / plausibly contributes to" rather than "is why."

6. **gapmm2 G4 — `--cs=short` vs `--cs`.** Confirm the `short` qualifier against the
   installed wheel; the README text verified only `minimap2 -x splice --cs`.

7. **minimap2 M6 (cosmetic).** `max_gap_ref=200000` in the preset-string table is more
   precisely `bw_long=200000`. Minor; fix or drop.

8. **mapPacBio B10 (cosmetic).** Drop or hedge the "~0.56 short-read default" comparison;
   versions disagree. The load-bearing 0.40 mapPacBio default is solid.

---

## 3. Confidence Assessment (per dossier)

**minimap2.md — HIGH confidence.** Every numeric default verified verbatim against
`options.c`; `-k14` confirmed as the official DRS recipe (not RECTIFY-invented, an
explicit anti-overstatement point in the dossier's favor); algorithm description
consistent with the paper. Only real issue: the splice-flank *rationale* leans on
CLAUDE.md over the code comment, and the 0.1% win rate is CLAUDE.md-only. INFERENCE
labels are used responsibly throughout. **Trust the facts; treat the 3'-end weakness
story as well-reasoned hypothesis, not measured.**

**deSALT.md — HIGH confidence on facts, MEDIUM on the causal thesis.** This was the
highest-risk claim to stress-test ("is cross-read junction sharing real or an
overstatement?") and it **survived**: the Genome Biology paper literally says exon
inference *"integrates all the alignment skeletons of all the reads"* and markets
"homogeneous" alignments. All params verified verbatim. The weakness is purely
rhetorical: the chain from this real mechanism to the specific **78.9%** is asserted
too confidently for an INFERENCE. Engineering-fragility claims are source-verified.

**gapmm2.md — MEDIUM-HIGH confidence.** README facts (defaults, 13 bp example, 409/63
benchmark, v25.4.13 mappy-segfault note) all verified verbatim. The dossier is
*commendably honest* in its header that gapmm2's source is sparse and partly
reconstructed from README + RECTIFY's defensive code rather than a line-by-line wheel
read — so the edlib-internals (`mode="HW"`, slide, splice-dinucleotide scan) are
README/inference-level, not independently confirmed at the source line. RECTIFY-side
integration (clean FASTQ, cs-overrun skip, PAF→BAM) is fully source-verified. Minor:
confirm `--cs=short`.

**mapPacBio_bbmap.md — HIGH confidence.** The crucial disambiguation (BBMap
`align2.BBMapPacBio` ≠ minimap2 map-pb) is correct and the script defaults + affine
weights + memory figures verified verbatim. One soft spot remains: the "default k=13" label
overreaches what the cited guide passage actually says. The earlier maxindel concern is
resolved — mapPacBio sets `intronlen=10` and an explicit `maxindel=max(200000, max_intron)`
(`multi_aligner.py:749,754`), so long introns are safely within the searched range. RECTIFY-side
quirks (pt-strip, 6 kb split/stitch, unmapped guard) are fully source-verified.

---

## 4. Cross-Cutting Observations

- **No fabricated primary-source facts found.** Every external number I checked
  (minimap2 preset scores, deSALT defaults, gapmm2 benchmark, BBMap weights/defaults)
  verified verbatim. This is a genuinely careful set of dossiers.

- **The systemic weakness is one pattern repeated 5×:** a real, source-verified
  *mechanism* → an un-measured *win-rate percentage*, presented as causation. Every
  dossier labels these INFERENCE (good discipline) but writes them with prose certainty
  ("decisive," "the central reason," "why one aligner takes 4/5"). Recommend a single
  shared disclaimer that the win rates are CLAUDE.md measurements not reproduced here,
  and that no per-mechanism decomposition was performed.

- **One genuine internal contradiction:** minimap2 `--splice-flank=no` rationale
  (code comment "compatibility" vs CLAUDE.md/dossier "3' end accuracy").

- **Coordinate/parameter hygiene is excellent** — the RECTIFY-side claims I sampled
  (gapmm2 `-i 5000`/`-m` omission, mapPacBio `intronlen=10` + `maxindel=max(200000, max_intron)`/
  pt-strip/split-stitch, minimap2 flag set, deSALT no-`-x`/no-`-G`) all matched
  `multi_aligner.py` exactly.

---

### Sources

- minimap2 `options.c` (mm_set_opt presets), README (DRS `-k14` recipe) — lh3/minimap2.
- deSALT Genome Biology 2019 (Liu et al., 13059-019-1895-9) "integrates all the
  alignment skeletons of all the reads"; "homogeneous"; deSALT README options table.
- gapmm2 README (nextgenusfs/gapmm2) — defaults, 13 bp example, 6926/409/63 benchmark,
  v25.4.13 note.
- BBMap `sh/mapPacBio.sh` (default-args line) + `docs/guides/BBMapGuide.txt`
  (+100/−127/−51 weights, maxindel 16000, 6/3 bytes-per-base, Nanopore-recommended) —
  BioInfoTools/BBMap mirror.
- RECTIFY `rectify/core/align/multi_aligner.py` (`run_minimap2` L217, `run_map_pacbio`
  L418, `run_gapmm2` L1052, `run_desalt` L1506) — read directly.
