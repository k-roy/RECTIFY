# Sumner SMA — RNA modification analysis scoping (2026-07-09)

**Confirmed:** the Dorado WGS BAMs carry per-read modification calls. Model
`rna004_130bps_sup@v5.2.0_inosine_m6A_2OmeA@v1` → THREE mods, all on adenosine (MM/ML tags, all 15 samples):
- `A+a` = **m6A** (N6-methyladenosine) — splicing, stability, translation, export
- `A+17596` = **inosine** (A-to-I editing, ADAR) — can CREATE/DESTROY splice sites, edit codons, alter miRNA sites
- `A+69426` = **2′-O-methyl-A (Am)** — ribose methylation, RNA structure/stability (snoRNA-guided)
No prior mod analysis (only an unrelated `fingerprint/pileup.txt`). modkit NOT installed (need it — ONT's
standard mod-pileup/differential tool; static Rust binary, easy).

## Why RNA mods are relevant to SMA (compelling)
SMN's core function is snRNP assembly; snRNAs are heavily modified. m6A/editing/Am regulate the SAME processes
(splicing, stability, translation) dysregulated in SMA → mod changes could be an unexamined SMA axis AND could
mechanistically EXPLAIN the splicing/APA leads (sweeps B/C).

## Analysis plan (SMA vs WT, good-RIN subset 4 SMA / 2 WT for robustness)
1. TOOL: install modkit → `modkit pileup` per sample → bedMethyl (per-site mod fraction, per mod code), with
   `--filter-threshold` on ML (Dorado probability) to call a base modified. Then differential SMA vs WT.
2. PER-SITE DIFFERENTIAL (per mod): sites with adequate coverage in ≥2 SMA + ≥1 WT → mod-fraction shift SMA vs WT.
3. GENE/REGION aggregation: mod burden per gene/transcript (m6A, inosine, Am separately); genes with differential mod.
4. ★ FUNCTIONAL-SITE FOCUS (highest-value, cross-modal):
   - **Inosine at/near SPLICE SITES** — A-to-I editing that creates/destroys a splice site could EXPLAIN a sweep-B
     non-canonical junction. The single most compelling cross-modal link (editing → splicing).
   - **m6A in 3′UTRs** — stability/APA crosstalk → ties to sweep-C APA shifts.
   - **m6A/Am/editing in SMN-pathway + lead genes**: SMN2, UBA1, SNRPN, GTF2H2, and the sweep-A/B hits.
5. SMN-PATHWAY / SMA-MODIFIER cross-ref (PLS3, NCALD, ...) for differential mods.

## Caveats
- Small n + variable depth → per-site power limited; lean on high-coverage sites + gene-level aggregation.
- RIN: mods are per-BASE (less RIN-sensitive than APA) but degradation biases coverage → use good-RIN subset.
- ML threshold matters (calibrate; modkit --filter-threshold). m6A false-positives → require DRACH motif; inosine
  vs m6A both on A → the model separates them, but validate context (inosine splice-adjacent, m6A DRACH).
- Editing hotspots (Alu) → high inosine background; focus on genic/functional sites, not repeat-driven.

## Deliverable: modkit pileup per sample → per-mod differential → ranked mod-lead table (site · gene · mod ·
fraction SMA/WT · Δ · functional context · caveat), merged into the SMA lead set. → cluster sma_leads/.
