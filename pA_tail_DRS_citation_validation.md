# Citation Validation — pA tail DRS adapter executive summary

Verification pass over every external citation I recommended for the HTML executive summary. Each entry records (a) what I originally claimed, (b) what the paper actually says, and (c) whether to keep, revise, or drop the citation. Verbatim quotes from primary sources are ≤15 words each.

Today: 2026-04-18.

---

## 1. Tudek et al. 2021 — yeast poly(A) tail metabolism

**Verdict: VERIFIED — keep; correct the assay description.**

Correct citation:
Tudek A, Krawczyk PS, Mroczek S, Tomecki R, Turtola M, Matylla-Kulińska K, Jensen TH, Dziembowski A. *Global view on the metabolism of RNA poly(A) tails in yeast Saccharomyces cerevisiae*. Nat Commun 12, 4951 (2021). DOI [10.1038/s41467-021-25251-w](https://doi.org/10.1038/s41467-021-25251-w). PMC: [PMC8367983](https://pmc.ncbi.nlm.nih.gov/articles/PMC8367983/).

**Assay actually used:** nanopore direct RNA sequencing (DRS) with poly(A)+ oligo‑dT25 selection, poly(A) lengths estimated via Nanopolish. **Not PAL-seq and not TAIL-seq** — I got the assay name wrong in my earlier draft. This is important because it means the Tudek reference numbers come from the same DRS + Nanopolish stack Kevin is running, so they are directly comparable.

Reported yeast poly(A) lengths (verbatim, ≤15 words each):

| Context | Number | Quote |
|---|---|---|
| WT steady state, total poly(A)+ | mean 40 nt, median 37 nt | "mean and median poly(A) tails of 40 and 37 As, respectively" |
| Highly expressed mRNAs | mean 30 nt, median 26 nt | "Highly expressed mRNAs had relatively short tails" |
| Newly transcribed (Pap1-dependent) | mean 50 nt, max ~200 nt | "newly transcribed mRNAs are 50 adenosine long on average" |
| ncRNA range | 20–60 nt | "ncRNA poly(A) tails are 20-60 adenosines long" |
| Post-deadenylation steady state | ~40 nt | "trim the tails to 40 adenosines on average" |

**"25.6 nt" does NOT appear anywhere in the paper.** My earlier attribution was wrong on two axes: the number and the assay name. Withdraw both.

**Per-read vs per-gene caveat (important).** Tudek reports "mean and median poly(A) tails of 40 and 37 As" and later "Highly expressed mRNAs had relatively short tails with a mean of 30 As and a median of 26 As." Both sentences are textually ambiguous about whether the statistics are (a) the distribution over individual tail lengths (per-read) or (b) the distribution over per-gene mean/median values. Standard DRS reporting convention and the phrasing "mean and median poly(A) tails" lean toward per-read, but the paper does not explicitly say "per read." Fig. 3C — the scatter panel Kevin referenced — is unambiguously **per-gene** (each dot = one gene, y = mean tail per gene, r = -0.5357). So when comparing to Kevin's pt:i: distribution, the relevant question is which version of the "40/37" and "30/26" statistic the paper is quoting. Do not assume either without checking Tudek's source data (MOESM series on PMC) or re-reading the Fig. 1c legend.

---

## 2. Begik et al. Nano3P-seq — yeast median 23.5 nt

**Verdict: VERIFIED — keep.**

Correct citation:
Begik O, Diensthuber G, Liu H, Delgado‑Tejedor A, Kontur C, Niazi AM, et al. *Nano3P-seq: transcriptome-wide analysis of gene expression and tail dynamics using end-capture nanopore cDNA sequencing*. Nat Methods 20, 75–85 (2023). DOI [10.1038/s41592-022-01714-w](https://doi.org/10.1038/s41592-022-01714-w). PMC: [PMC9834059](https://pmc.ncbi.nlm.nih.gov/articles/PMC9834059/).

Primary-source quote (verbatim, ≤15 words):
"median polyA tail length of 23.5 nucleotides, in agreement with previous studies" (yeast, Fig. 3d).
Paired comparator in same panel: mouse brain 90 nt median; yeast is the shortest of four species.

**Caveats that matter for the Section 4 comparison:**
- Nano3P-seq is **cDNA with template-switch end-capture**, not DRS. It avoids the adapter-trim boundary problem entirely, which makes its numbers a useful independent ground truth.
- Note the "in agreement with previous studies" hedge — Nano3P-seq and Tudek 2021 report slightly different yeast medians (23.5 vs 37 nt), reflecting differences in selection (end-capture all RNA vs poly(A)+ enrichment) and tail-call algorithm.
- **Per-read vs per-gene caveat:** The Fig. 3d legend explicitly states the panel displays "polyA tail length distribution of yeast, zebrafish, and mouse mRNAs represented as single-transcript values (left) and per-gene medians (right)." The body-text quote reporting 23.5 nt does not specify panel. The 23.5 nt could therefore be from the per-read (single-transcript) distribution or the per-gene median distribution. Pin this down from the figure image itself before using the number as a direct comparator to a pt:i: statistic whose weighting is known.
- Any claim that Kevin's pt:i: medians "match" 23.5 nt is currently weight-indeterminate on both sides of the comparison; don't present it as a confirmation.

---

## 3. Chang et al. 2025 — Sequins benchmark of poly(A) tools

**Verdict: VERIFIED — keep; slightly tighten the phrasing.**

Correct citation:
Chang JJ, Yang X, Teng H, Reames B, Corbin V, Coin L. *Using synthetic RNA to benchmark poly(A) length inference from direct RNA sequencing*. GigaScience 14 (2025). DOI [10.1093/gigascience/giaf098](https://doi.org/10.1093/gigascience/giaf098). Preprint: [bioRxiv 2024.10.25.620206](https://www.biorxiv.org/content/10.1101/2024.10.25.620206v1.full).

Primary-source quote (verbatim, 15 words):
"All four tools generate mean tail-length estimates which lie within 12% of the correct value."

**Keep but phrase carefully.** The paper says "within 12% of the correct value" as an aggregate across four tools — it does NOT single out Dorado for ±12%. My earlier "±12%" reframing is slightly stronger than what the abstract supports; replace with "within 12%". The paper does recommend Dorado specifically for speed and integration, not because it had the smallest error among the four tools.

Spike-ins are **Sequins** (synthetic RNA standards with known tail lengths), as I cited. That's correct.

---

## 4. Dorado GitHub Issue #1131 — HAC vs SUP poly-A artifacts

**Verdict: VERIFIED — keep; update the description.**

Correct reference:
[nanoporetech/dorado Issue #1131](https://github.com/nanoporetech/dorado/issues/1131) — "Bases are emitted from strided blocks within polyA region (SQK-RNA004)". Filed 2024-11-13 by `magmir71`. Still **open** as of this validation.

What the issue actually documents (verbatim excerpts, ≤15 words each):
- "in default HAC model, many samples from presumable polyA tail region emit bases"
- "With the SUP model, there was one strided block... right at the edge of the... adapter and polyA-tail region"
- Models referenced explicitly: `rna004_130bps_hac@v5.1.0` and `rna004_130bps_sup@v5.1.0`.
- Uses SQK-RNA004 mouse direct RNA data.

**Caveat:** the issue describes artifactual *indel emission* inside the poly-A region in HAC, not specifically the pt:i: / pa:f: tag estimates. Rephrase any claim that frames the issue as "pt:i: tag estimation is affected" — the direct evidence is for basecalled-sequence artifacts at the tail/adapter boundary, which could propagate into downstream poly-A length inference but the issue itself doesn't pin down that tag-level effect.

---

## 5. Bauer et al. 2025 MRA — adapter sequence disclosure

**Verdict: CONDITIONALLY OK — keep only if you want a citable source for the 11-bp stub sequence; see caveat below.**

Correct citation:
Bauer M, et al. *Duplexed direct RNA sequencing protocol using polyadenylation and polyuridylation*. Microbiol Resour Announc 14(3):e01041-24 (2025). DOI [10.1128/mra.01041-24](https://doi.org/10.1128/mra.01041-24). PMC: [PMC11895459](https://pmc.ncbi.nlm.nih.gov/articles/PMC11895459/). Note: published 2025, not 2024.

Primary-source quote of the relevant sequence (verbatim, 14 words):
"OligoB (5'-/5Phos/GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC-3')... was used in place of the RTA reagent"

**Caveat:** this paper uses the sequence `GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC` in a **custom poly(U) duplex protocol** replacing the standard RTA reagent — it's not a reference for what the stock SQK-RNA002/SQK-RNA004 RTA oligo is. That said, the 11-bp prefix `GGCTTCTTCTT` exactly matches the first 11 bases of the sequence Kevin is observing as a pre-trim stub. The Bauer paper is useful as one of the few openly-published sources that discloses this sequence verbatim. If a stronger citation is needed for the standard-kit RTA sequence, look for an ONT protocol document or an earlier method paper; do not lean on Bauer as the sole authority.

---

## 6. TAILcaller — Bioinformatics Advances 2025

**Verdict: VERIFIED — keep.**

Correct citation:
Maździarz M. *TAILcaller: an R package for analyzing differences in poly(A) tail length for Oxford Nanopore RNA sequencing*. Bioinform Adv 5(1):vbaf235 (2025). DOI [10.1093/bioadv/vbaf235](https://doi.org/10.1093/bioadv/vbaf235). PMC: [PMC12527348](https://pmc.ncbi.nlm.nih.gov/articles/PMC12527348/). Published 2025-09-26.

Verbatim supporting quote (≤15 words):
"TAILcaller supports direct RNA sequencing data generated by dorado (--estimate-poly-a)"

- Input: BAM files with pt:i: / pa:f: from `dorado --estimate-poly-a`.
- Statistical tests: Wilcoxon, Student's t, Welch's t, Kruskal-Wallis, ANOVA, Welch ANOVA, Games-Howell post-hoc.
- Does NOT benchmark Dorado's estimate accuracy itself, and does NOT compare HAC vs SUP. Cite it as tooling for downstream comparison only, not as an accuracy assessment.

---

## 7. Ninetails — non-adenosine detection in poly(A) tails

**Verdict: VERIFIED — keep.**

Primary citation:
Gumińska N, et al. *Direct profiling of non-adenosines in poly(A) tails of endogenous and therapeutic mRNAs with Ninetails*. Nat Commun 16 (2025). PMC: [PMC11920217](https://pmc.ncbi.nlm.nih.gov/articles/PMC11920217/). Repo: [github.com/LRB-IIMCB/ninetails](https://github.com/LRB-IIMCB/ninetails).

Only load-bearing claim: Ninetails is the tool for resolving non-A residues (U, C, G) inside poly(A) tails at the DRS squiggle level. Relevant to the executive summary as context for why "poly(A)" tails are not pure A-runs, but not directly relevant to the adapter-boundary issue unless the report is making that broader point.

---

## 8. BMC Genomics PMC11134706 — DRS systematic errors

**Verdict: DROP the citation as I originally framed it. Paper is real but does not support what I claimed.**

Correct citation:
Liu-Wei W, et al. *Sequencing accuracy and systematic errors of nanopore direct RNA sequencing*. BMC Genomics (2024). DOI [10.1186/s12864-024-10440-w](https://doi.org/10.1186/s12864-024-10440-w). PMC: [PMC11134706](https://pmc.ncbi.nlm.nih.gov/articles/PMC11134706/).

What it actually says (verbatim ≤15 words each):
- "deletions significantly outnumbered mismatches and insertions"
- "C bases are more likely to be basecalled as U"
- "G bases are more often confused with A"
- Adapter note: "Guppy began basecalling at the beginning of the adapter" — BUT "we did not observe that [this] negatively impacts the error rate of the subsequent RNA basecalls"

My earlier wording — "elevated C→T and G→A substitution rates when RNA-trained models process non-RNA signal" — is wrong on three points: (1) the substitutions reported are C→U and G→A, not C→T; (2) the paper explicitly states adapter‑region basecalls do **not** degrade downstream RNA accuracy; (3) there's no "RNA-trained model on non-RNA signal" framing in the paper. **Drop this citation** from any section that frames it as evidence of adapter-contamination-induced substitution error. If you want an evidence-based citation for DRS systematic error patterns (deletion-dominated, C/G homopolymer penalty), this paper works fine — just re-phrase to match what it actually says.

---

## Summary table — action per citation

| # | Citation | Action | Notes |
|---|---|---|---|
| 1 | Tudek 2021 | KEEP, correct assay; flag weighting | DRS + Nanopolish, not PAL-seq. Numbers 40/37, 50, 20–60 nt. Per-read vs per-gene ambiguous in text. |
| 2 | Begik 2023 Nano3P-seq | KEEP, flag panel | 23.5 nt yeast median. Fig. 3d has BOTH per-read and per-gene panels; which one the 23.5 refers to is unspecified in body text. |
| 3 | Chang 2025 GigaScience | KEEP, tighten phrasing | "within 12%" not "±12%". Four-tool aggregate, not Dorado-specific. |
| 4 | Dorado Issue #1131 | KEEP, re-describe | Artifact is base emission in tail region, not pt:i: directly. Still open. |
| 5 | Bauer 2025 MRA | KEEP with caveat | Shows 11-bp stub sequence but in a non-stock duplex protocol. |
| 6 | TAILcaller 2025 | KEEP | Downstream tooling; not an accuracy benchmark. |
| 7 | Ninetails 2025 | KEEP | Non-A residue detection — contextual, not load-bearing. |
| 8 | BMC PMC11134706 | DROP as framed | Paper is real but I mischaracterized the finding. |

---

## Cross-cutting caveat: per-read vs per-gene statistics

This ambiguity recurs across citations 1, 2, and Kevin's own pt:i: numbers, so it deserves its own section.

A poly(A) tail-length "median" or "mean" can be computed in two non-equivalent ways:

1. **Per-read** (aggregate over all reads in the library). The distribution weights are proportional to transcript abundance × sequencing depth. A per-read median is pulled toward whichever transcripts supply the most reads — typically a small number of highly abundant mRNAs.
2. **Per-gene** (compute a per-gene summary first, then summarize the per-gene distribution). Each gene contributes one data point regardless of how many reads it had. A per-gene median reflects the typical gene, not the typical read.

When the underlying per-gene tail-length distribution is anti-correlated with abundance (as Tudek's Fig. 3C shows at r = -0.54 for yeast), the two statistics **diverge**: per-read statistics run below per-gene statistics. Apples-to-apples requires the same weighting on both sides of the comparison.

What's known and unknown for the numbers in this sheet:

| Source | Number | Weighting |
|---|---|---|
| Tudek 2021, "40 / 37 nt" WT steady-state | 40 mean / 37 median | **Ambiguous in text.** Phrasing "mean and median poly(A) tails" leans per-read; paper does not explicitly say. Check MOESM source data before citing. |
| Tudek 2021, "30 / 26 nt" highly expressed | 30 mean / 26 median | **Ambiguous in text.** Same paragraph as above. |
| Tudek 2021, Fig. 3C scatter | each dot = one gene | **Per-gene.** Figure legend and axis label confirm. Correlation r = -0.5357 between per-gene mean tail length and log2 abundance. |
| Begik 2023 Nano3P-seq, "23.5 nt" yeast median | 23.5 nt | **Ambiguous in text.** Fig. 3d legend explicitly shows BOTH per-read (left) and per-gene medians (right) for yeast/zebrafish/mouse. Body text cites "Fig. 3d" without panel. Inspect the figure image to determine which panel. |
| Kevin's pt:i: 20–26 nt | 20–26 nt | **Unknown — please confirm.** Was this a median of pt:i: values across all reads (per-read), a median of per-gene medians (per-gene), or something else? |

**Operational implication:** Before the Section 4 comparison can be tightened, we need (a) the Tudek Fig. 3C x-y coordinates (MOESM source data) to produce a per-gene reference distribution, (b) confirmation of the Nano3P-seq 23.5 nt panel, and (c) the weighting of Kevin's pt:i: statistic.

---

## Proposed replacement text for Section 4 of the HTML report

Because the published yeast poly(A) reference numbers (Tudek 2021 "40/37", "30/26"; Begik 2023 "23.5") are textually ambiguous between per-read and per-gene weighting, and because the weighting of the pt:i: median in this report has not been confirmed, a direct "our 20–26 nt matches/does not match published yeast values" comparison is not yet supported. What can be said with the primary sources available:

1. Per-gene, Tudek et al. 2021 (DRS + Nanopolish, *S. cerevisiae*) shows an r ≈ -0.54 anti-correlation between mean tail length and abundance (their Fig. 3C), with the most abundant mRNAs sitting at per-gene means of ~25–35 nt and the bulk of moderately expressed mRNAs at ~40 nt.
2. Nano3P-seq (cDNA with template-switch end-capture, Begik et al. 2023) reports the shortest yeast tail distribution among four species — 23.5 nt median — from either a per-read or a per-gene summary of Fig. 3d (panel to be confirmed).
3. The adapter-boundary effect seen in nanoporetech/dorado issue #1131 (HAC vs SUP basecaller model, still open) is an orthogonal source of pt:i: variability that is not accounted for by mass bias alone.

The cleaner Section 4 deliverable is a per-gene scatter of the pt:i: data on this library plotted against log2 TPM, overlaid on (or directly compared to) Tudek Fig. 3C coordinates for the same yeast genes. If the red "highly expressed" cluster reproduces at ~25–35 nt and the overall per-gene mean sits near 40 nt with r in the -0.4 to -0.6 range, the Dorado pt:i: estimates are accurate at gene level; the low global median is a mass-bias artifact and not a basecaller undercount. If the per-gene scatter is systematically shifted below Tudek's for matched genes, that quantifies the undercount. Same pod5 → HAC vs SUP → per-gene pt:i: delta is the within-library control for the issue #1131 effect.

Sources:
- [Tudek et al. 2021, Nat Commun — PMC8367983](https://pmc.ncbi.nlm.nih.gov/articles/PMC8367983/)
- [Begik et al. 2023 Nano3P-seq, Nat Methods — PMC9834059](https://pmc.ncbi.nlm.nih.gov/articles/PMC9834059/)
- [Chang et al. 2025 GigaScience](https://doi.org/10.1093/gigascience/giaf098)
- [Dorado Issue #1131](https://github.com/nanoporetech/dorado/issues/1131)
- [Bauer et al. 2025 MRA — PMC11895459](https://pmc.ncbi.nlm.nih.gov/articles/PMC11895459/)
- [TAILcaller 2025 — PMC12527348](https://pmc.ncbi.nlm.nih.gov/articles/PMC12527348/)
- [Ninetails 2025 — PMC11920217](https://pmc.ncbi.nlm.nih.gov/articles/PMC11920217/)
- [Liu-Wei et al. 2024 BMC Genomics — PMC11134706](https://pmc.ncbi.nlm.nih.gov/articles/PMC11134706/)
