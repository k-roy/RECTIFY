# Non-canonical junction benchmark survey (2026-07-07, web survey agent)

**The PI's half-remembered study = LRGASP** (Pardo-Palacios et al., Nat Methods 2024): synthetic (NanoSim,
per-read origin truth) + real (Illumina, PacBio, ONT-cDNA, ONT-**dRNA**) + SIRV spike-in, scored on the
SQANTI3 FSM/ISM/**NIC/NNC** novel-junction taxonomy. The closest match to "synthetic+real, short+long, with a
non-canonical subset." Portal: cgl.gi.ucsc.edu/data/LRGASP.

## Ranked datasets (value axis = non-canonical ground truth in RECTIFY's ONT-dRNA modality)
1. **SIRV-Set 4** (Lexogen) — physical spike-in already run on ONT dRNA inside SG-NEx & LongBench. 7 genes/69
   isoforms; junctions **96.9% GT-AG, GC-AG 1.7%, AT-AC 0.6%, CT-AG 0.4%, CT-AC 0.4%** (verified, Lexogen
   design page). Absolute truth (GTF), immune to ONT error. BORROW-FIRST, but non-canonical set is thin + annotated.
2. **SG-NEx** (`s3://sg-nex-data/`, ENA PRJEB44348): real ONT dRNA RNA002 + cDNA + PacBio + short-read, with a
   combined GRCh38+SIRV+Sequin+ERCC annotation → SIRV non-canonical in a real human dRNA background.
3. **LongBench** (`s3://longbench-data/`, bioRxiv 2025.09.11.675724): the ONLY SIRV-in-real-dRNA on **RNA004**
   (RECTIFY's target chemistry). Highest-priority real-modality test. POD5 present; spike-in ref fetched separately.
4. **Sequins** (sequin.xyz): synthetic in-silico-chr spike-in; predominantly CANONICAL → weak non-canonical source.
5. **intropolis / Snaptron** (github.com/nellore/intropolis): 43-81M real Illumina junctions across SRA/GTEx/TCGA
   with per-junction sample counts → cross-sample corroboration ORACLE (a novel dRNA junction gains credibility
   if recurrent in intropolis). Not dRNA, not per-read truth.
6. **COMPASS yeast** (Roy et al. 2023 NAR gkad968, PI's own; GEO GSE131797 — VERIFY): real yeast non-canonical
   (non-YAG) 3'SS, RT-PCR-validated → reusable real-yeast non-canonical truth in RECTIFY's native organism.
7. **IAOD / U12DB** (human minor-intron AT-AC reference) — human arm only (yeast has no minor spliceosome).

## SPIKE-IN prior art (directly informs task #20)
- **Physical spike-in** (SIRV/Sequins): already gives non-canonical junctions in a real ONT-dRNA background
  (SG-NEx, LongBench) — but ANNOTATED, not unannotated-novel.
- **SQANTI-SIM** (Mestre-Tomás et al. 2023, github.com/ConesaLab/SQANTI-SIM): THE published precedent for the
  PI's idea — creates ground-truth NOVEL junctions by REMOVING real transcripts from the reference (making them
  NIC/NNC), then simulates reads (wraps NanoSim cDNA+dRNA, PBSIM3, Polyester). "Inject known novel junctions to
  get ground truth," validated + published. RECTIFY's own NIC/NNC injector is a hand-built cousin → ADOPT/align to it.
- **★ GENUINELY NOVEL DESIGN SPACE (the defensible contribution):** NO published work spikes UNANNOTATED
  NON-CANONICAL reads into a real ONT-dRNA background. SIRV = real-background but annotated; SQANTI-SIM =
  unannotated-novel but simulated background. **Combining them — SQANTI-SIM-style novel non-canonical junctions
  layered onto a SIRV-in-real-dRNA (LongBench RNA004 / SG-NEx) background — is the PI's spike-in idea and is new.**

## Metric convention to adopt (from LRGASP + the field)
- SQANTI3 **FSM/ISM/NIC/NNC** stratification + **motif-stratified** precision/recall/FDR (separate non-canonical
  from canonical). Position-exact AND windowed (±10 bp, minimap2 convention). Ambiguity-aware normalization
  (collapse to canonical rep within the repeat window before scoring — RECTIFY's chimeric_consensus mirrors this).
- **uLTRA's warning = the central pitfall:** junction-penalty aligners MANUFACTURE spurious GT-AG by snapping to
  the motif (the exact failure the empirical penalty tables must avoid; also why short-read Q2 is motif-circular).
- Simulator shape-deficiency trap (pbsim over-dispersion ~13x too low vs real SIRV) — the "0.09→1.07" trap; sim
  wins are necessary not sufficient → must transfer to SIRV/real.

## Recommended adopt/adapt (accessions in the agent's full report)
1. Score **LongBench RNA004** SIRV spike-in reads (real-modality non-canonical numbers) vs pbsim3 Tier-2.
2. Adopt **LRGASP + SQANTI3 NIC/NNC + motif-stratified** reporting as RECTIFY's metric convention.
3. Adopt/align to **SQANTI-SIM** for controlled-novelty NIC/NNC truth (don't hand-roll).
4. Reuse **COMPASS/GSE131797** as the real-yeast non-canonical truth (RT-PCR-anchored).
5. (human) **intropolis/Snaptron + IAOD** as corroboration oracles for real non-canonical / AT-AC.

## Flagged UNVERIFIED (confirm before citing): GSE131797 accession; LRGASP Synapse syn25683367/syn25683630;
Sequins GSE; LongBench spike-in identity (SIRV-Set4/Sequins per PI notes, not on the GitHub page); ENA PRJEB44348.
Correction: PMID 37956322 = the NAR Prp18p paper that INTRODUCED COMPASS (not a standalone COMPASS methods paper).
