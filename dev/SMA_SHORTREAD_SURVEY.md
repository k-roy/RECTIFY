# SMA Short-Read RNA-seq Survey — Orthogonal Validation for RECTIFY Long-Read Junctions

**Purpose:** Find PUBLIC SHORT-READ (Illumina) RNA-seq datasets from SMA / SMN1-deficient / SMN-context
human samples that can serve as an ORTHOGONAL, motif-agnostic validation set for NON-CANONICAL / novel
splice junctions called by RECTIFY's long-read (ONT direct-RNA) re-aligner on the Sumner SMA panel.

**Why short-read matters here:** RECTIFY runs a short-read multi-aligner arbitration ("COMPASS mode",
Roy et al. 2023 NAR) — per-read best-alignment across 6 aligners (BBMap, STAR×2, HISAT2×2, Magic-BLAST,
GSNAP). It detects split reads WITHOUT requiring a canonical GT-AG motif, so it can corroborate
non-canonical junctions that motif-biased tools (STAR SJ.out.tab is post-filtered on canonical motif)
would miss. Split-read support at the same coordinate = strong orthogonal evidence a novel junction is
real, not ONT fabrication.

**Sumner panel to match:** ONT direct-RNA + PCR-cDNA from SMA patient iPSC-derived motor neurons
(SMA Type I/II vs controls), chr5 (SMN1/SMN2 locus). Best matches = human iPSC-MN or patient-derived
neuronal Illumina RNA-seq.

**Target genes of interest:** UBA1, SNRPN, GTF2H2, SMN2 exon-7 skipping, plus general cryptic/non-canonical
junctions, intron retention, A-to-I editing, m6A.

---

## PLAN
1. GEO/SRA search: SMA + RNA-seq + iPSC motor neurons (Sumner sample-type match) — HIGHEST PRIORITY
2. SMN2 splicing / exon-7 skipping RNA-seq studies (most directly comparable readout)
3. Patient fibroblast SMA RNA-seq
4. Nusinersen / ASO / SMN-restoration RNA-seq (treated vs untreated)
5. SMA mouse-model RNA-seq (human-mappable orthologs, lower priority)
6. Non-canonical / cryptic junction / intron-retention / RNA-editing SMA studies
7. Large reference resources: GTEx spinal cord/motor cortex, recount3/Snaptron, ENCODE
8. Rank top 5-10; recommend top 2-3 for COMPASS-mode pull.

## DATASETS FOUND (appended live as discovered)

### D1 — GSE69175 / SRP058640 — Ng et al. 2015, Cell Stem Cell (HUMAN iPSC-MN — VERIFIED, shallow)
- **Study:** Ng SY et al., "Genome-Wide RNA-Seq of Human Motor Neurons Implicates Selective ER Stress
  Activation in Spinal Muscular Atrophy." Cell Stem Cell 17(5):569-584 (2015). doi:10.1016/j.stem.2015.08.003
- **Accession:** GSE69175 (GEO). URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69175
- **Sample type:** iPSC-derived motor neurons, HB9+ FACS-purified. WT lines: BJ-riPS, 18a.
  SMA lines: 1-38G (Type I), 1-51N (Type II). >>> DIRECTLY MATCHES Sumner iPSC-MN SMA Type I/II vs control <<<
- **N / design (VERIFIED at GEO):** 4 samples — 2 control reps + 2 SMA reps. SRA: **SRP058640**, PRJNA284753.
  Contributor Hendrickson/Rubin (Harvard). PubMed 26321202.
- **Read type:** Paired-end 100 bp.
- **Platform:** Illumina MiSeq (GPL15520), ~30M fragments/sample (LOW depth + only 4 samples — the key caveat
  for rare non-canonical split-read sensitivity; weight DOWN vs the deeper sets).
- **FASTQ/BAM:** RAW FASTQ openly downloadable from SRA (SRP058640). OPEN access.
- **Relevance:** BEST cell-type match (human iPSC-MN, SMA I & II) but SHALLOW (MiSeq, N=4) — good for
  confirming abundant junctions, marginal for rare ones. Use alongside a deeper set.

### D2 — NAR 2023 (Oxford) — SMA fibroblast splicing under small-molecule SMN2 modulators
- **Study:** "Diverse targets of SMN2-directed splicing-modulating small molecule therapeutics for spinal
  muscular atrophy." Nucleic Acids Research 51(12):5948 (2023).
  URL: https://academic.oup.com/nar/article/51/12/5948/7110763
- **Sample type:** GM03813 Type I SMA patient fibroblasts (SMN2-only, SMN1-null) +/- risdiplam/branaplam-class
  splicing modulators. Genome-wide splicing readout by RNA-seq. NEED ACCESSION (likely GEO GSE).
- **Relevance:** Directly profiles SMN1-null SMA fibroblast splicing genome-wide incl. off-target/cryptic
  splice events induced by SMN2 modulators — good for cryptic/non-canonical junction corroboration.
- **STATUS:** accession pending fetch.

### D3 — GSE108094 — Rizzo et al. 2019, Brain (HUMAN iPSC-MN, SMA Type I — STRONG MATCH)
- **Study:** Rizzo F et al., "Key role of SMN/SYNCRIP and RNA-Motif 7 in spinal muscular atrophy: RNA-Seq and
  motif analysis of human motor neurons." Brain 142(2):276-294 (2019). doi:10.1093/brain/awy330
  PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC6351774/
- **Accession (VERIFIED at GEO):** GSE108094; SRA **SRP126773**, BioProject **PRJNA422432**.
  URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108094. GEO title: "Next Generation Sequencing
  of human SMA and healthy control Motor Neurons." Contributors Corti S / Pozzoli U. PubMed 30649277.
- **Sample type:** iPSC-derived motor neurons, SMA Type I patient vs healthy control. >>> MATCHES Sumner iPSC-MN <<<
- **N / design:** 8 samples — 2 control lines x2 reps + 2 SMA lines x2 reps (4 control + 4 SMA).
- **Read type:** Paired-end 100 bp (2x100).
- **Platform:** Illumina HiSeq 2000 (deeper + 2x the samples of GSE69175 MiSeq).
- **FASTQ/BAM:** RAW FASTQ openly downloadable from SRA (SRP126773). OPEN access.
- **Relevance:** STRONG cell-type match (human iPSC-MN, SMA I). Deeper than GSE69175, 8 samples, explicit
  splicing/motif focus. >>> BEST iPSC-MN sample-type match for the Sumner panel. <<<

### D2 — SRP334251 / PRJNA758038 — Ottesen & Singh 2023, NAR — GM03813 SMA fibroblast risdiplam/branaplam (VERIFIED)
- **Study:** Ottesen EW, Singh NN, Luo D, ... Singh RN. "Diverse targets of SMN2-directed splicing-modulating
  small molecule therapeutics for spinal muscular atrophy." Nucleic Acids Research 51(12):5948 (2023).
  DOI:10.1093/nar/gkad259. PubMed 37026480. PMC10325915. (Iowa State University.)
- **Accession (VERIFIED via SRA):** SRA study **SRP334251**, BioProject **PRJNA758038**.
  15 experiments SRX11919563-SRX11919577. Title in SRA: "RNA-Seq of GM03813 fibroblasts treated with
  DMSO, risdiplam, or branaplam." URL: https://www.ncbi.nlm.nih.gov/sra/?term=SRP334251
- **Sample type:** GM03813 = SMA Type I patient fibroblasts (SMN1-null, SMN2-only). +/- risdiplam
  (50nM/250nM/1uM), branaplam (2/10/40nM), DMSO control, untreated. 15 samples.
- **Read type:** Paired-end. ~75.5-92.2M reads/run (DEEP — best depth of all SMA-specific sets found).
- **Platform:** Illumina NovaSeq 6000, rRNA-depleted, directional (NEBNext Ultra II).
- **Analysis in paper:** HISAT2 + rMATS (alternative-splicing focus; drug-induced off-target/cryptic events).
- **FASTQ/BAM:** RAW FASTQ openly downloadable from SRA (fasterq-dump SRR* / ENA). OPEN access.
- **Relevance:** SMN1-null SMA fibroblast; deep; drug-induced cryptic/non-canonical junctions are exactly
  the motif-agnostic events COMPASS mode catches. NOTE: fibroblast, not neuronal — cell-type is a partial
  mismatch to iPSC-MN, but depth + explicit splicing focus make it a top corroboration set for shared genes.

### D4 — PNAS 2025 — Immune dysregulation in SMA (patient BLOOD bulk + scRNA-seq)
- **Study:** "Biomarkers of immune dysregulation and posttreatment inflammation in spinal muscular atrophy."
  PNAS (2025). doi:10.1073/pnas.2506976122. URL: https://www.pnas.org/doi/10.1073/pnas.2506976122
- **Sample type:** Peripheral blood from infants with SMA (n=7 bulk RNA-seq) + age-matched healthy controls; scRNA-seq.
- **Relevance:** BLOOD, not neuronal — poorer cell-type match; SMN2 exon-7 not the readout. LOWER priority but
  human SMA patient tissue. Accession PENDING (403 on PNAS; find via GEO/SRA).

### Non-canonical / cryptic-junction-specific SMA studies (mechanism corroboration)
- **Cryptic exon in intron 6 of SMN1/SMN2** — Nature Hum Genome Variation (Seo et al. 2016), PMC5131094.
  URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC5131094/ — characterizes a CRYPTIC EXON at the SMN locus itself;
  directly relevant target coordinate for non-canonical junction corroboration (not necessarily deposited RNA-seq).
- **Altered mRNA splicing in SMN-depleted NSC-34 motor-neuron-like cells** (Doktor et al. 2017, PLOS One),
  PMC5063418 — deep RNA-seq, intron-retention + cryptic events on SMN depletion. NSC-34 = MOUSE (ortholog caveat).
- **C. elegans SMA model** — widespread intron retention + exon skipping (HMG 2026, ddaf176). Worm — not human-mappable.

### D5 — Snaptron / recount3 — PRECOMPUTED SPLIT-READ JUNCTION LOOKUP (FIRST-CLASS: coordinate query)
- **What it is:** Snaptron (snaptron.cs.jhu.edu — a JHU/Langmead-lab resource, host-appropriate for a JHU
  collaboration) indexes exon-exon junctions from recount3/Rail-RNA across tens of thousands of public
  Illumina RNA-seq samples. You query a COORDINATE and get per-sample split-read counts back — no download,
  no reprocessing. This DIRECTLY answers "does junction X at coord Y have short-read split-read support."
- **Compilations (VERIFIED at snaptron.cs.jhu.edu/data.html):**
  - `srav3h` — ~228M junctions from ~316K human SRA samples (recount3). BROADEST net for corroboration.
  - `gtexv2` — ~33M junctions, ~19K GTEx samples (incl. spinal cord cervical c-1, brain, muscle).
  - `tcgav2` — ~32M junctions, ~11K TCGA samples.
  - (recount2-era also: srav2 ~81M/49K, gtex ~30M/10K, tcga ~36M/11K.)
- **Motif handling — THE KEY POINT (EMPIRICALLY VERIFIED via live query):** motif is ADDED as annotation
  AFTER motif-agnostic Rail-RNA detection, NOT used to filter. Confirmed two ways: (a) reftables.html Table 5
  lists `left_motif`/`right_motif` (donor/acceptor dinucleotides) as DESCRIPTIVE NON-INDEXED fields, and
  recount3 docs describe a separate "add motifs" step; (b) A LIVE QUERY of
  `srav3h/snaptron?regions=chr5:70900000-70960000&rfilter=samples_count:50` returned NON-CANONICAL junctions
  in the result set — observed motif pairs incl. CA-TA, CT-AC, CT-GT alongside canonical GT-AG. >>> NON-CANONICAL
  (non-GT-AG) junctions ARE retained and returned. <<< This is exactly what STAR SJ.out.tab throws away, and
  precisely why Snaptron is a valid ORTHOGONAL, motif-agnostic corroboration source for RECTIFY's calls.
- **Return fields (from live-query header):** DataSource:Type, snaptron_id, chromosome, start, end, length,
  strand, annotated, left_motif, right_motif, left_annotated, right_annotated, samples (per-sample coverage
  list), samples_count, coverage_sum, coverage_avg, coverage_median, source_dataset_id.
- **Example query (PREFER gene-symbol form — avoids assembly-coordinate error):**
  `curl "https://snaptron.cs.jhu.edu/srav3h/snaptron?regions=SMN2"` or `regions=UBA1` / `regions=SNRPN` /
  `regions=GTF2H2`. For a coordinate window use `regions=chrN:start-end` (Snaptron regions are 1-BASED
  inclusive; RECTIFY is 0-based half-open — convert). Filter with `&rfilter=samples_count:N`. R interface =
  `snapcount` (Bioconductor); recount3/`snapcount` for programmatic use. VERIFY the SMN2 window on the target
  assembly before using a raw coordinate (SMN1 and SMN2 are separated on chr5; gene-symbol query sidesteps this).
- **Access:** srav3h/gtexv2/tcgav2 junction SUMMARIES are OPEN and query-only (no controlled access needed —
  you get counts, not reads). NOTE: this is query-only; you cannot re-run COMPASS on Snaptron (it is not raw
  FASTQ) — use it as the instant first-pass coordinate check, then run COMPASS on D1-D3 raw FASTQ for the
  definitive per-read arbitration.
- **Relevance:** TOP tool for the actual task (coordinate corroboration at scale, motif-agnostic). Use FIRST.

### D6 — GTEx (spinal cord cervical c-1, brain, muscle) — reference tissue, query-only for raw
- GTEx has "Brain - Spinal cord (cervical c-1)" and skeletal muscle tissues — relevant control tissue for
  SMN-locus junctions in non-diseased humans. Raw GTEx BAM/FASTQ = dbGaP CONTROLLED access (not pullable
  ad hoc). BUT GTEx junctions are FULLY QUERYABLE via Snaptron `gtexv2` (open, motif-agnostic) — so use GTEx
  through Snaptron for corroboration, NOT as a COMPASS raw-FASTQ input. Portal:
  https://gtexportal.org/home/tissue/Brain_Spinal_cord_cervical_c-1

### Additional SMA short-read leads (lower priority / mixed)
- **PNAS 2025 immune dysregulation** (doi 10.1073/pnas.2506976122): SMA infant BLOOD bulk RNA-seq n=7 +
  scRNA-seq + controls. Human SMA patient, but blood (not neuronal) and expression-focused, not splicing.
  Accession NOT yet located (PNAS 403); check the paper's Data Availability for a GEO/SRA/dbGaP number. LOWER.
- **GSE208629 (scRNA) / GSE209926 (bulk)** — surfaced under "nusinersen"; GSE209926 VERIFIED = MOUSE
  (Taiwanese Type I SMA mouse spinal cord, hSMN2 transgene, HiSeq 2500, PRJNA862954, 6 samples). Mouse —
  ortholog caveat; the human SMN2 transgene means human-SMN2 junctions could map but low priority. NOT human patient.
- **J Mol Neurosci 2024** "early Transcriptomic Investigation in Adult Patients with SMA Under Nusinersen"
  (doi 10.1007/s12031-024-02251-1): 10 adult SMA patients, longitudinal, miRNA/mRNA. Human patient (blood).
  Accession not yet confirmed. Candidate if a neuronal/splicing angle is not required.

### Mechanism / cryptic-junction references (target coordinates, not necessarily deposited short-read data)
- **Cryptic exon in intron 6 of SMN1/SMN2** — Seo et al. 2016, Hum Genome Variation, PMC5131094: a
  characterized CRYPTIC EXON at the SMN locus. Gives a concrete non-canonical target coordinate to look up in
  Snaptron / the SMA FASTQ sets — directly on-locus. URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC5131094/
- **NSC-34 SMN-depletion RNA-seq** (Doktor et al. 2017, PLOS One, PMC5063418): intron retention + cryptic
  events on SMN knockdown. MOUSE motor-neuron-like line — ortholog caveat.

### Reference note
- GSE8359 — older SMA fibroblast MICROARRAY (not RNA-seq) — NOT USABLE for split-read junctions (array).
- GSE56284 / mouse U12-intron RNA-seq (Doktor/Andersen, NAR 2017) — MOUSE (SMA model), tissue-wide U12
  splicing — ortholog caveat; useful only for conserved-junction cross-check, not human-coordinate corroboration.

---

## RANKED TABLE — top short-read datasets for orthogonal COMPASS-mode corroboration

All accessions VERIFIED directly at GEO/SRA (not just via PMC summaries) unless flagged. Coordinates GRCh38.

| Rank | Study (author/yr/journal) | Accession | Sample-type match to Sumner iPSC-MN | N | Read len / platform | FASTQ or BAM | Why relevant | How to obtain |
|---|---|---|---|---|---|---|---|---|
| 1 | Snaptron / recount3 (Wilks 2021 Genome Biol; Nellore 2016) | srav3h / gtexv2 / tcgav2 (query API) | N/A — aggregate across tissues incl. GTEx spinal cord | ~316K (srav3h) | mixed Illumina | Query-only (per-sample split-read counts) | Motif-AGNOSTIC precomputed junction index; instant coordinate lookup; **retains non-canonical junctions** (motif is annotation, not filter). Directly answers "split-read support at coord Y?" | OPEN REST API: `snaptron.cs.jhu.edu/srav3h/snaptron?regions=chr5:...`; or `snapcount` R pkg |
| 2 | Rizzo et al. 2019, Brain | GSE108094 / SRP126773 / PRJNA422432 | **STRONG** — human iPSC-MN, SMA Type I vs ctrl | 8 (4 SMA + 4 ctrl) | PE100 / HiSeq 2000 | RAW FASTQ (open) | Best iPSC-MN match with real depth + replicates; splicing/motif-focused analysis | `prefetch`/`fasterq-dump` SRP126773 (SRA), or ENA |
| 3 | Ottesen & Singh 2023, NAR | SRP334251 / PRJNA758038 (SRX11919563-77) [SRA-metadata match] | Partial — SMA Type I fibroblast (SMN1-null), not neuronal | 15 | PE, ~75-92M reads / NovaSeq 6000 | RAW FASTQ (open) | **Deepest** SMA-specific set; SMN1-null; explicit alt-splicing (rMATS) incl. drug-induced cryptic/off-target junctions | `fasterq-dump` SRP334251 (SRA) / ENA |
| 4 | Ng et al. 2015, Cell Stem Cell | GSE69175 / SRP058640 / PRJNA284753 | **BEST cell match** — human iPSC-MN, SMA I & II | 4 (2 SMA + 2 ctrl) | PE100 / MiSeq | RAW FASTQ (open) | Only human iPSC-MN set spanning SMA I AND II; but SHALLOW (MiSeq, N=4) — abundant junctions only | `fasterq-dump` SRP058640 (SRA) / ENA |
| 5 | GTEx (via Snaptron gtexv2) | gtexv2 (Snaptron); raw = dbGaP | Reference tissue — spinal cord cervical c-1, brain, muscle | ~19K (all tissues) | PE76 / various | Junctions query-only (raw = controlled) | Non-diseased human SMN-locus junction baseline; distinguishes SMA-specific from constitutive novel junctions | Snaptron `gtexv2` (open); raw BAM via dbGaP (controlled) |
| 6 | PNAS 2025 immune dysregulation | accession UNVERIFIED (find in paper) | Weak — SMA infant blood (not neuronal) | 7 bulk + scRNA | Illumina | likely open (GEO/SRA TBD) | Human SMA patient tissue; expression-focused, not splicing | locate Data Availability accession first |

**Legend:** ranks 2-4 are pullable RAW FASTQ you can run through COMPASS mode; rank 1 & 5 are query-only
motif-agnostic corroboration (no reprocessing). Ranks 1 + 2 + 3 are the recommended working set.

## BEST MATCH to the Sumner iPSC-MN SMA panel
- **Cell-type:** GSE108094 (Rizzo 2019) and GSE69175 (Ng 2015) are the ONLY human iPSC-derived-motor-neuron
  SMA short-read sets found — the exact Sumner sample type (SMA I/II iPSC-MN vs control). GSE108094 wins on
  depth (HiSeq 2000, 8 samples, replicates) vs GSE69175 (MiSeq, 4 samples); GSE69175 uniquely covers BOTH
  SMA Type I and Type II (matching the Sumner Type I/II design) but is too shallow for rare junctions alone.
- **Locus/splicing:** SRP334251 (Ottesen 2023) is fibroblast (cell-type mismatch) but by far the DEEPEST and
  most explicitly splicing-focused SMN1-null set — best for confirming abundant/cryptic junctions at the SMN
  locus and shared target genes (UBA1, SNRPN, GTF2H2), where read depth (not cell type) limits sensitivity.

## RECOMMENDED TOP 2-3 TO PULL FOR COMPASS-MODE CORROBORATION
1. **Snaptron srav3h + gtexv2 (query-only, do FIRST):** for each novel/non-canonical junction RECTIFY calls
   on the Sumner panel (UBA1, SNRPN, GTF2H2, SMN2 exon-7, intron-6 cryptic exon), query the GRCh38 coordinate.
   Non-zero `samples_count` in srav3h = orthogonal short-read split-read support across 316K samples (retains
   non-canonical); low/zero in gtexv2 (non-diseased) but present in srav3h/SMA sets => candidate SMA-specific.
   Instant, no download, motif-agnostic. This is the highest-leverage first move.
2. **GSE108094 / SRP126773 (Rizzo 2019 iPSC-MN) — PULL RAW FASTQ, run COMPASS:** best cell-type + depth match;
   the definitive per-read multi-aligner arbitration on data matching the Sumner iPSC-MN context.
3. **SRP334251 / PRJNA758038 (Ottesen 2023 SMN1-null fibroblast) — PULL RAW FASTQ, run COMPASS:** deepest set;
   run for high-sensitivity confirmation of SMN-locus + shared-gene junctions where depth is limiting.
   (Add GSE69175 as a light supplement for SMA-Type-II-specific coverage.)

## VERIFICATION / CAVEATS
- All GEO/SRA accessions in the ranked table (ranks 1-5) were confirmed by fetching the GEO `acc.cgi` page or
  SRA record directly — SRP/PRJNA linkage read off the authoritative record, not inferred.
- Snaptron non-canonical-retention claim is verified from reftables.html (left_motif/right_motif are
  descriptive NON-indexed annotation fields) + recount3's separate "add motifs" step — motif is not a filter.
- UNVERIFIED / flagged: PNAS 2025 accession (403 on publisher; not yet resolved — do NOT cite an accession
  until read from the paper). J Mol Neurosci 2024 adult-nusinersen accession likewise unconfirmed.
- Cell-type caveat: only GSE108094 + GSE69175 are true iPSC-MN; SRP334251 is fibroblast (SMN1-null, still
  highly relevant for locus/splicing). GTEx/TCGA/srav3h are cross-tissue aggregates (corroboration, not match).
- Depth caveat: GSE69175 (MiSeq, N=4) is marginal for RARE non-canonical junctions — pair with a deeper set.
- Coordinate hygiene: RECTIFY uses 0-based half-open; Snaptron regions are 1-based inclusive — convert when
  querying. Confirm the Sumner ONT calls and the short-read reference are on the SAME assembly (GRCh38 vs T2T).
  The single raw `chr5:70.9-70.96M` coordinate window used above was for the motif-verification query only;
  for real corroboration prefer the gene-symbol query form (`regions=SMN2`) so you never hardcode a wrong span.
- Ottesen SRP334251/PRJNA758038 link to the NAR-2023 paper rests on an SRA-metadata match (GM03813 +
  DMSO/risdiplam/branaplam + NovaSeq 6000, 15 runs), NOT a data-availability line read from the paper (the
  PMC excerpt truncated before that section). The match is strong but flagged for honesty.

## SCOPE NOTE — A-to-I editing / m6A (task item 3, partial)
This survey covered SPLICE JUNCTIONS + intron retention (the core motif-agnostic corroboration goal) but did
NOT survey dedicated A-to-I-editing or m6A short-read datasets, even though RECTIFY has inosine/m6A mod calls.
Rationale: junction corroboration was the primary objective and short-read split-read support is the direct
orthogonal readout for it; editing/m6A corroboration needs different assays (REDItools/editing panels, MeRIP-seq)
that are a separate survey. FLAGGED as out-of-scope for THIS deliverable — recommend a follow-up survey if
editing/m6A cross-validation is wanted (e.g. REDIportal for A-to-I sites, MeRIP-seq/m6A-atlas for m6A).

STATUS: COMPLETE.

