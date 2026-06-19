# Red-Team Adversarial Audit — Aligner-Ecosystem Dossiers

**Auditor role:** Skeptical adversarial fact-checker / red-team reviewer (bioinformatics).
**Date:** 2026-06-19.
**Targets (in `dev/aligner_investigation/01_investigation/`):**
- `uLTRA.md`
- `pacbio_ecosystem.md`
- `ont_drs_ecosystem.md`

**Method:** Each substantive claim classified **FACT / PLAUSIBLE INFERENCE / SUSPECT-WRONG**.
Shaky/numeric claims verified against primary sources (paper PDFs, repo source, man pages,
tool docs) via web search + fetch. Source quotes embedded. Calibrated skepticism: well-supported
claims are credited explicitly; only genuinely shaky items are escalated.

**Headline verdict:** The three dossiers are **unusually high-quality and well-sourced**. The
hard technical cores (uLTRA's algorithm, minimap2/pbmm2 preset numerics, the Iso-Seq pipeline
order, the ONT error model, dorado tag semantics) are **verified correct**. The defects are
concentrated in (1) **two specific unverifiable numbers** (uLTRA "8-nt exon", and the
98.67% RNA004 identity), and (2) **a handful of INFERENCE claims dressed with `[FACT]`-adjacent
confidence** — most importantly the ont_drs "AT-vs-CG reconciliation" and the "no mainstream
tool corrects per-read 3' ends" universal negative. None of the corrections overturn a dossier's
central thesis.

---

## 1. Audited-Claims Table

Legend: ✅ FACT (verified) · 🟡 PLAUSIBLE INFERENCE (reasoned, not directly sourced; label is
fair) · 🟠 OVER-LABELED (claim is sound but tagged `[FACT]` when it is really inference) ·
🔴 SUSPECT-WRONG / UNVERIFIABLE.

### 1A. uLTRA.md

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| U1 | `colinear_solver.py` cites "Algorithm 15.1 in *Genome-Scale Algorithmic Design*, Mäkinen et al." | ✅ FACT | Verified in source: docstring literally reads `"Algorithm 15.1 in Genome scale algorithmic design, Makinen et al."` (raw GitHub `modules/colinear_solver.py`). |
| U2 | Four functions: `read_coverage`, `n_logn_read_coverage`, `read_coverage_mam_score`, `n_logn_read_coverage_mams`; segment trees `RMaxQST` | ✅ FACT | All four present in source; file imports `range_query_max_search_tree as RMaxQST` and calls `construct_tree/range_query/update`. |
| U3 | MAM chaining allows overlaps with fractional penalty `0.1*(gap)` | ✅ FACT | Source: `T_values = [(j_prime, c_val - 0.1*(v.c - mams[j_prime].d - 1))...]` in `read_coverage_mam_score`; `ovl_penalty = ... + 0.0001` in n_logn variant. The `0.1*gap` form is exact. |
| U4 | `--ont` sets `min_acc=0.6, s=9, mm2_ksize=14`; `--isoseq` sets `min_acc=0.8, s=10` | ✅ FACT | Verified against `uLTRA` argparse source: `--ont` → min_acc 0.6 / s 9 / mm2_ksize 14; `--isoseq` → min_acc 0.8 / s 10 (mm2_ksize stays default 15). Exactly as written. |
| U5 | `--genomic_frac` default 0.10; reads >10% outside indexed regions fall to minimap2; `--disable_mm2` disables both genomic detection and QC | ✅ FACT | Source help text: "threshold for fraction of aligned portion outside uLTRA indexed regions to be considered genomic (default 0.1)"; `--disable_mm2`: "Disables utilizing minimap2 to detect genomic primary alignments and to quality check uLTRAs primary alignments." |
| U6 | `--disable_infer` is an indexing-only gffutils speedup requiring pre-existing gene/transcript features | ✅ FACT | uLTRA README: speeds indexing "but it only works if you have the gene feature and transcript feature in your GTF file." Indexing-only scope is correct. |
| U7 | namfinder (NAMs via strobemers) supersedes slaMEM MEMs; namfinder is current default | ✅ FACT | namfinder repo: "Finds Non-overlapping Approximate Matches (NAMs) ... using strobemers"; uLTRA deps list namfinder; PyPI confirms NAMs "constructed from overlapping strobemer seeds." slaMEM was the paper-era MEM finder. |
| U8 | parasail (sg_trace_scan_16, match+2/mismatch−2/gap_open3/gap_ext1) + edlib base alignment; edlib for refs >20 kb (HW mode) | 🟡 PLAUSIBLE (mostly FACT) | parasail + edlib dependency confirmed (PyPI/Bioconda). The exact parasail numerics (2/−2/3/1) and the 20 kb edlib switch are cited to `help_functions.py`; I did not independently re-open that file, but the dossier flags them as source-grounded and they are internally consistent. **Spot-check recommended** but no red flag. |
| U9 | **Paper's "headline 8-nt exon example"** missed by minimap2 & deSALT | 🔴 **UNVERIFIABLE / LIKELY WRONG NUMBER** | The paper's small-exon *result* is binned: ~60% accuracy for exons **≤10 nt**, ~90% for 11–20 nt (verified, biorxiv v2 + Oxford Academic). The widely-cited tiny-exon failure case in the minimap2 tracker is **6 nt** (lh3/minimap2 issue #253: "Alignment ... to a tiny exon (6 nucleotides)"). I found **no source for a specific "8-nt" headline example.** The qualitative claim (uLTRA recovers micro-exons others miss) is solid; the **specific "8 nt" figure is unsupported** and should be removed or corrected to "≤10 nt (paper bin)". |
| U10 | uLTRA strips poly-A via `remove_read_polyA_ends()`; offers no CPA/3'-end precision, only junction backbone | 🟡 PLAUSIBLE (well-reasoned) | Poly-A trimming behaviour is consistent with the tool's design; the "no 3'-end precision → junction backbone only" is the dossier's own inference, **correctly labeled [INFERENCE]**. Reasonable. |
| U11 | ~2% win rate is because uLTRA can only beat de-novo aligners on well-annotated FSM isoforms; on compact yeast its niche is narrow | 🟡 PLAUSIBLE INFERENCE | This is **correlation-as-mechanism but labeled honestly** as INFERENCE. The causal chain (annotation-bound → small upside on simple yeast genes → 2%) is plausible and internally consistent, but it is a *story*, not a measurement. The dossier does not overclaim it. Acceptable. |
| U12 | FSM reads exempt from the `alignment_score < 2*threshold*len` reject filter; novel-junction `>10 bp` gap rule | 🟡 PLAUSIBLE | Cited to `align.py`/`classify_alignment2.py`. Consistent with uLTRA's SQANTI-style design. Not independently opened; no contradiction with paper. Spot-check optional. |

### 1B. pacbio_ecosystem.md

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| P1 | minimap2 `splice:hq` = `splice` then override `-C5 -O6,24 -B4` | ✅ FACT | minimap2 man page (Debian): `splice:hq` = "-xsplice -C5 -O6,24 -B4". Exact match. |
| P2 | `splice` preset: k15/w5/-B2/-O2,32/-C9/-G200000/--junc-bonus 9/--splice-flank=yes | ✅ FACT | Consistent with minimap2 man page / options.c; widely corroborated. |
| P3 | pbmm2 links minimap2 C API directly (no shell-out); native pbbam I/O; `-Y` softclip, `--eqx` | ✅ FACT | pbmm2 README ("a minimap2 SMRT wrapper ... native PacBio data in ⇨ native PacBio BAM out"). |
| P4 | HPC minimizers ON for SUBREAD/UNROLLED, OFF for CCS/ISOSEQ (`-u`) | ✅ FACT | pbmm2 README: "-H ... is always on for SUBREAD and UNROLLED presets and can be disabled with -u." CCS/ISOSEQ use -u (off). |
| P5 | pbmm2 ≥v1.13.0 → minimap2 v2.26; <v1.13.0 → v2.15 | ✅ FACT | pbmm2 README verbatim. |
| P6 | ISOSEQ preset numeric string `-k15 -w5 -u -o2 -O32 -e1 -E0 -A1 -B2 -z200 -Z100 -r200000 -L0.5 -g2000 -C5 -G200000` | 🟡 PLAUSIBLE, version-flagged | pbmm2 README explicitly says preset numbers "vary by version" and are in `--help`; v26.1.0 changed ISOSEQ numbers. The dossier **correctly flags this as version-dependent and tells the reader to re-verify**. ISOSEQ exposes -G/-C/--no-splice-flank (confirmed). The k15/w5 and -u are confirmed-consistent; the full numeric string is **not independently verified to the digit** but is appropriately caveated. |
| P7 | Iso-Seq order: `ccs → lima → refine → cluster2 → pbmm2 align (ISOSEQ) → collapse → pigeon/SQANTI3`; alignment is a distinct step between clustering and collapse | ✅ FACT | IsoSeq docs (isoseq.how clustering CLI). collapse requires the sorted aligned BAM. Correct. |
| P8 | `isoseq refine --require-polya` REQUIRES ≥`--min-polya-length` (default 20) poly-A AND **strips the tail**; output = FLNC | ✅ FACT | IsoSeq clustering docs: refine with `--require-polya` keeps FL reads with ≥20 bp poly(A) "and removes the identified tail." The dossier's load-bearing "trimmed BEFORE alignment" claim is **correct**. |
| P9 | refine's exact poly-A detector unpublished; historical `trim_isoseq_polyA` used an HMM | 🟡 PLAUSIBLE, honestly hedged | Labeled [INFERENCE] for refine's internal algorithm. Fair. |
| P10 | HiFi reads have **no systematic homopolymer-deletion bias**; CCS averages stochastic per-pass errors | ✅ FACT (mechanistically standard) | Well-established; CCS = circular consensus over passes → random errors cancel. Correct contrast with ONT systematic deletion. |
| P11 | RNA004 indel proportion 7.19% (RNA002) → 0.88% (~8×); RNA004 median identity 98.67% | 🟠 MOSTLY FACT, one number soft | 7.19%→0.88% **verified exactly** (biorxiv 2025.05.01.651717, "Evaluation of Nanopore DRS updates"). The **98.67% median identity is NOT independently confirmed in the snippets** — the paper exists and discusses identity, but the exact decimal traces to a search snippet. Dossier already lists it under "re-verify against source PDFs." Acceptable with that caveat. |
| P12 | Homopolymers ≥5 bp correctly called ≈ HiFi 99% / raw ONT ~91% / CLR ~88% | 🟡 PLAUSIBLE, snippet-sourced | Directionally standard (HiFi ≫ ONT ≫ CLR on HP). Exact percentages are snippet-level, not traced to a single canonical table. Treat as approximate; dossier does not over-anchor them. |
| P13 | GMAP = minimal-sampling oligomer lookup → oligomer chaining → sandwich DP → microexon ID; Wu & Watanabe Bioinformatics 2005 21:1859–1875 | ✅ FACT | Verified: GMAP method = "minimal sampling ... oligomer chaining ... sandwich DP for splice site detection ... microexon identification"; citation & pages exact. |
| P14 | deSALT `-x` accepts `ccs/clr/ont1d/ont2d`; `-T` forward-strand splice | ✅ FACT | deSALT README: `-x ccs -O6,24 -M4` for CCS; `-x ont1d` for ONT. Preset names correct. |
| P15 | dorado aligner wraps minimap2 internally | ✅ FACT | dorado docs: "alignment uses minimap2 ... by default uses the lr:hq preset"; internal minimap2 v2.28. (Minor: dossier could note default `lr:hq`, not load-bearing.) |
| P16 | SQANTI3 intra-priming flag = ≥60% downstream A's → likely false 3' end; flags/filters whole isoforms post hoc vs RECTIFY per-read | 🟡 PLAUSIBLE | SQANTI3 intra-priming concept confirmed (Pardo-Palacios Nat Methods 2024). The "≥60%" threshold is a known SQANTI default. The per-read-vs-per-isoform contrast is sound inference. |
| P17 | §6 porting plan (flip to `splice:hq`/`ccs`/`isoseq`, disable HP/poly-A modules, don't reuse ONT penalty tables) | 🟡 PLAUSIBLE INFERENCE | Honestly tagged [INFERENCE]; analytically sound. The one caveat-worthy step: "drop `-k14`, HiFi wants k15" — correct, but `splice:hq` already sets k15, so the instruction is slightly redundant (flipping preset already does it). Minor. |

### 1C. ont_drs_ecosystem.md

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| O1 | DRS sequenced **3'→5**', poly-A in-read, no priming step (so no internal mispriming) | ✅ FACT | Standard ONT DRS biology; matches ONT docs + literature. The "no priming → AG module off for DRS" logic is correct and matches CLAUDE.md. |
| O2 | Motor falls off **~10–12 nt (≈11 nt)** before the 5' terminus; final ~6 nt unresolvable | 🟡 PLAUSIBLE (could not directly re-quote) | Cited to PMC8906549 (poly(A) isoform scaffolds). The ~11 nt motor fall-off is a **frequently-reported DRS figure** and is consistent with the literature, but I hit 403s re-fetching the primary and could not re-quote the exact "≈11 nt / ~6 nt" numbers. **Plausible and load-light** (used only to justify "5' = rescue, not trust"), but the specific numbers rest on a single citation — flag as not-independently-reconfirmed. |
| O3 | dorado `--estimate-poly-a` (off by default); `pt:i:N`(N>0)=length, `pt:i:0`=anchor found/length not estimable, `pt:i:-1`=anchor not found; supports cDNA(PCS/PCB)+RNA | ✅ FACT | dorado docs verbatim: "-1 ... if the primer anchor for the tail was not found, or ... 0 if the primer anchor was found, but the length could not be estimated." Exact. `--estimate-poly-a` off by default, cDNA+RNA support: confirmed. |
| O4 | Aligner wrapping: FLAIR→minimap2→BED12; Bambu→minimap2 `-ax splice`; NanoCount→minimap2-to-transcriptome+EM; Sicelore→minimap2; dorado aligner→minimap2 | ✅ FACT | FLAIR docs (minimap2→BED12, common TSS/TES); NanoCount (NAR 2022, minimap2 transcriptome EM); Sicelore README (minimap2 mapping); dorado (minimap2). All correct. |
| O5 | **"No mainstream tool corrects per-read 3' ends"** (the RECTIFY niche) | 🟠 OVER-STRONG UNIVERSAL NEGATIVE (but defensible) | This is the dossier's central differentiation claim and is **labeled [INFERENCE]**, which is appropriate. It is **probably true** for the surveyed tools (FLAIR/LAPA cluster; NanoCount/Bambu quantify; signal tools do tail-length). BUT a universal "no tool does X" is unfalsifiable from a finite survey and ignores adjacent work (e.g., signal-level 3'-end refinement, Nano3P-seq's protocol-level 3'-precision, `2passtools` junction filtering). **Soften** to "no surveyed mainstream pipeline performs per-read 3'-end *coordinate* correction against an HP error model." The dossier mostly does this in §5.2 ("No surveyed tool combines (a)+(b)+(c)") — that triad-uniqueness framing is the **defensible** version and should be the headline, not the blunter §6.4 wording. |
| O6 | **AT-vs-CG reconciliation:** literature says A/T basecalled *more accurately* than C/G; CLAUDE.md table says AT *deletes more* than CG; dossier argues "not contradictory" because accuracy ≠ deletion-operation | 🔴 **WEAKEST CLAIM — correlation/mechanism handwave** | Two halves: (a) "A/T HP basecalled more accurately than C/G" — **VERIFIED** (PMC10070092 benchmarking; multiple R10.4 studies). (b) The *reconciliation* — that AT can simultaneously have higher overall accuracy AND higher deletion-specific rate — is **logically possible but presented with more confidence than the evidence supports.** The proposed mechanism ("A/T-rich k-mers pass faster → under-segmented → more deletions") is the dossier's own **[INFERENCE consistent with...]**, not a sourced finding. **No primary source is cited that directly measures deletion-rate(AT) > deletion-rate(CG) at equal HP length** independent of RECTIFY's own table. So the dossier is using RECTIFY's table to validate RECTIFY's table. This is **circular**. Verdict: the two statements are *plausibly* non-contradictory, but the dossier should state plainly that the AT>CG **deletion** direction is supported **only by RECTIFY's internal empirical table**, not by external literature, and that the "ratcheting speed" mechanism is unproven. |
| O7 | minimap2 `-ax splice -uf -k14`: 10M human DRS reads <1 wall-hr/16 cores, 94.2% junctions annotation-consistent | 🟡 PLAUSIBLE, single-source | Specific numbers (94.2%, <1 hr) are cited to the minimap2 paper/arXiv. Directionally credible; exact figures not re-verified here. Used only to support "fast & good on junctions ≠ 3'-precise," which is sound. |
| O8 | deSALT 78.9 / mapPacBio 18.2 / uLTRA 2 / gapmm2 0.8 / minimap2 0.1 win rates | ✅ FACT (internal, RECTIFY-measured) | These are RECTIFY's own measured correct-first win rates (match CLAUDE.md). Internally sourced, not external — fine, and consistently labeled as RECTIFY's measurement. |
| O9 | deSALT "homogeneous full-length alignments," small-exon recovery; two-pass de Bruijn skeleton | ✅ FACT | deSALT paper (Genome Biology 2019) language matches ("small exons, serious sequencing errors, consensus spliced alignment"). |
| O10 | scPAISO: "~½ of PASs span <20 bp"; LAPA/APALORD cluster 3' ends | 🟡 PLAUSIBLE, snippet-sourced | Tool roles (cluster-level) are correct and well-known. The "<20 bp ~½" specific figure is snippet-level; not load-bearing for the argument. |
| O11 | DRS chimera artifacts (two molecules basecalled as one) per a 2024 genomic-LM study | 🟡 PLAUSIBLE | Cited to PMC12923543. Chimera artifacts in DRS are a real, reported phenomenon; specific paper not re-opened. Low risk. |

---

## 2. Corrections Required

**Must-fix (factual error / unsupported number):**

1. **uLTRA.md U9 — the "8-nt exon" headline.** No source supports a specific "8 nt" example.
   The paper reports binned small-exon accuracy (**~60% at ≤10 nt**, ~90% at 11–20 nt) and the
   canonical tiny-exon failure case in the minimap2 tracker is **6 nt** (issue #253). **Action:**
   replace "8-nt exon" (appears in §"Pass 2", Strength #1) with "micro-exons (paper: ~60%
   accuracy at exon length ≤10 nt; minimap2 documented to miss exons as small as 6 nt)" or drop
   the specific integer. This appears **twice** and should be fixed in both.

2. **ont_drs_ecosystem.md O6 — AT-vs-CG deletion reconciliation is circular.** The claim that
   "AT deletes more than CG at equal HP length" is sourced **only to RECTIFY's own penalty
   table**, then used to argue consistency with external "A/T more accurate" literature. **Action:**
   state explicitly that (i) external literature confirms A/T HP *overall accuracy* > C/G, (ii)
   the AT>CG *deletion-rate* direction is **RECTIFY-internal/empirical, not externally
   corroborated**, and (iii) the "ratcheting-speed" mechanism is an unproven hypothesis. Do not
   present the reconciliation as settled.

**Should-fix (over-labeling / over-strong wording):**

3. **ont_drs_ecosystem.md O5 — "No mainstream tool corrects per-read 3' ends."** Soften the blunt
   universal negative (§6.4 takeaway #4) to the triad-uniqueness framing already used in §5.2
   ("no *surveyed* tool combines ensemble + per-read HP-error correction + correct-first
   selection"). A finite survey cannot prove a universal negative.

4. **pacbio_ecosystem.md P11 — 98.67% RNA004 identity.** Keep, but the figure is snippet-only
   (the 7.19%→0.88% indel numbers ARE verified; the identity decimal is not). Already caveated;
   leave the caveat in and do not promote to unqualified [FACT].

5. **pacbio_ecosystem.md P17 — "-k14 → HiFi wants k15."** Minor: `splice:hq` already sets k15, so
   "drop the -k14 override" is automatically satisfied by flipping the preset; reword to avoid
   implying a separate manual step.

**Nice-to-have (precision):**

6. **ont_drs O2 (≈11 nt motor fall-off)** and **O7 (94.2% junction figure)** rest on single
   citations I could not independently re-quote (403s). Not errors — flag as
   "single-source, not independently reconfirmed in this audit."

7. **pacbio P15 / ont O4 (dorado aligner):** add that dorado's *default* alignment preset is
   `lr:hq` (minimap2 v2.28) for completeness — currently just "wraps minimap2."

---

## 3. Per-Dossier Confidence Assessment

### uLTRA.md — **Confidence: HIGH (9/10)**
The algorithmic core is **independently verified at the source level**: Algorithm 15.1 citation,
all four chaining functions, RMaxQST segment trees, the `0.1*gap` overlap penalty, the `--ont`
(0.6/9/14) and `--isoseq` (0.8/10) parameter sets, `--genomic_frac=0.10`, `--disable_mm2`, and
`--disable_infer` semantics — **every one checked out against the live repo/paper.** The
FACT/INFERENCE discipline is exemplary; the ~2% win-rate causal story is honestly fenced as
inference. **Only material defect:** the fabricated-looking "8-nt exon" number (U9). Fix that and
this dossier is essentially bulletproof. The parasail/edlib numerics (U8) and `align.py` filter
details (U12) are source-cited but not re-opened here — low risk, optional spot-check.

### pacbio_ecosystem.md — **Confidence: HIGH (8.5/10)**
The load-bearing technical claims are **verified**: `splice:hq = -C5 -O6,24 -B4` (exact),
pbmm2 HPC on/off rule, minimap2 version bundling, the Iso-Seq pipeline order, and the
**central thesis** that PacBio trims poly-A *before* alignment at `isoseq refine --require-polya`
(`--min-polya-length 20`, tail removed). GMAP/STARlong/deSALT/SQANTI3/TAMA/Cupcake roles and
preset names are real and current. The dossier is **scrupulous about version-sensitivity**
(ISOSEQ numeric string flagged for re-verification; HiFi accuracy decimals flagged as
snippet-extracted). Deductions: the 98.67% identity decimal (P11) and HP-accuracy percentages
(P12) are snippet-level not source-traced; the §6 porting plan is sound but is analysis, not
fact (correctly labeled). Nothing here is wrong — the soft spots are all already flagged by the
author.

### ont_drs_ecosystem.md — **Confidence: MEDIUM-HIGH (7.5/10)**
Strong where it counts: DRS 3'→5' directionality, dorado `pt:i` tag semantics (**exact**),
the aligner-wrapping table (FLAIR/Bambu/NanoCount/Sicelore/dorado all → minimap2, verified),
RNA002→RNA004 indel reduction, and the DRS-vs-cDNA mispriming logic. **Two weaknesses pull the
score down:** (1) **O6, the AT-vs-CG reconciliation, is circular** — it leans on RECTIFY's own
table to argue consistency with external literature, and the proposed ratcheting mechanism is
unsourced; this is the single most "correlation-as-mechanism" passage across all three dossiers.
(2) **O5, the universal "no tool corrects per-read 3' ends,"** is an unfalsifiable absolute that
the §5.2 triad-uniqueness framing already states more defensibly. The motor-fall-off (~11 nt)
and 94.2%-junction numbers are single-source and were not re-confirmable in this audit. The
central argument (the ONT error model justifies per-read 3'-correction) survives intact; the
dossier just over-reaches in its two flagship rhetorical claims.

---

## 4. Cross-Cutting Observations

- **FACT/INFERENCE labeling is mostly excellent** and is the main reason these dossiers earn high
  marks. The failures are where an **inference is dressed in fact-adjacent confidence** (O6's
  reconciliation, O5's universal negative) or where **a specific number lacks a traceable source**
  (U9's "8 nt", P11's 98.67%).
- **One genuine circularity** (O6): RECTIFY's empirical penalty table is cited as evidence *for*
  the biological claim that motivates that same table. This should be disclosed, not smoothed over.
- **No fabricated tool names, preset names, or pipeline stages were found.** Every aligner preset
  (`splice:hq`, `-x ccs`, `--ont`, ISOSEQ), every Iso-Seq stage (ccs/lima/refine/cluster2/
  collapse), and every tag semantic (`pt:i:-1/0/N`) that I checked is real and current.
- **Numbers most worth re-verifying before any external publication:** uLTRA "8 nt" (remove),
  RNA004 98.67% identity (source PDF), ~11 nt motor fall-off (source PDF), 94.2% junction
  consistency (minimap2 paper), HP ≥5 bp accuracy percentages (single canonical table).

---

## 5. Sources Consulted

- minimap2 man page (Debian) — `splice` / `splice:hq` presets: https://manpages.debian.org/testing/minimap2/minimap2.1.en.html
- minimap2 issue #769 (pbmm2↔minimap2 preset mapping): https://github.com/lh3/minimap2/issues/769
- minimap2 issue #253 (6-nt exon alignment): https://github.com/lh3/minimap2/issues/253
- pbmm2 README & changelog: https://github.com/PacificBiosciences/pbmm2
- IsoSeq clustering / CLI workflow docs: https://github.com/PacificBiosciences/IsoSeq , https://isoseq.how/clustering/cli-workflow.html
- uLTRA repo (`uLTRA` main, `modules/colinear_solver.py`, README): https://github.com/ksahlin/ultra
- namfinder repo: https://github.com/ksahlin/namfinder ; ultra-bioinformatics PyPI: https://pypi.org/project/ultra-bioinformatics/
- Sahlin & Mäkinen, *Bioinformatics* 2021 / bioRxiv 2020.09.02.279208 (small-exon bins): https://academic.oup.com/bioinformatics/article/37/24/4643/6327681 , https://www.biorxiv.org/content/10.1101/2020.09.02.279208v2.full
- deSALT repo & paper (Genome Biology 2019): https://github.com/ydLiu-HIT/deSALT , https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6913027/
- GMAP (Wu & Watanabe, *Bioinformatics* 2005;21:1859–1875): https://academic.oup.com/bioinformatics/article/21/9/1859/409207
- dorado Poly(A) Estimation docs (pt:i semantics): https://software-docs.nanoporetech.com/dorado/latest/basecaller/polya_estimation/
- dorado alignment docs (wraps minimap2, lr:hq): https://software-docs.nanoporetech.com/dorado/latest/basecaller/alignment/
- ONT R10.4/R9.4.1 benchmarking (A/T vs C/G HP accuracy): https://pmc.ncbi.nlm.nih.gov/articles/PMC10070092/
- RNA004 evaluation (indel 7.19%→0.88%): https://www.biorxiv.org/content/10.1101/2025.05.01.651717v1.full
- FLAIR docs: https://flair.readthedocs.io/en/latest/modules.html ; NanoCount (NAR 2022): https://academic.oup.com/nar/article/50/4/e19/6439677 ; Sicelore: https://github.com/ucagenomix/sicelore
</content>
</invoke>
