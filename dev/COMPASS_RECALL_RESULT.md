# COMPASS RECALL BENCHMARK — the RECALL half of the make-or-break precision/recall question

**Agent:** COMPASS-RECALL agent · branch `worktree-agent-a25a2c1e784ad37dc` · started 2026-07-11
**Cluster:** Sherlock (`ssh sherlock`), env `compass` (STAR/hisat2/bbmap) or `rectify`.

## THE QUESTION
The spike-in (dev/SPIKEIN_RESULT.md) measured PRECISION (arm-B fabricates ~1.2% pool-dependent
non-HP k-bp drift; guard misses it) but STRUCTURALLY COULD NOT test RECALL, because pbsim
under-disperses ~13× so minimap2 never flattens in-sim (raw-mm2 recall flat ~0.83 across all
decoy distances). The RECALL test needs REAL data where minimap2 DOES flatten, with a
MOTIF-AGNOSTIC short-read truth (Snaptron/recount3 is CIRCULAR — STAR-built + motif-filtered).

**Q: does the long-read native re-placer (arm-B, motif-blind) RECOVER real non-canonical
junctions that raw minimap2 flattens, as confirmed by MOTIF-AGNOSTIC short-read data?**

## DATA (Sherlock, per brief)
- SMA short-read FASTQ `/scratch/users/kevinroy/sma_shortread/` — 2 SMA (SRR6376960, SRR6376962)
  + 2 CTRL (SRR6376956, SRR6376958), iPSC-MN, PE100, ~13-19GB each (GSE108094).
- Genome+GTF `/scratch/users/kevinroy/compass_a549/COMPASS/genome_references/GRCh38_gencode_v44.fasta(+.gtf)`.
- Long-read re-placer output `/scratch/users/kevinroy/sumner_gw/panel_deep/` (25% run) + sweep-B leads:
  UBA1 chrX:47211252-47212397, SNRPN chr15:24977922-24978198, CACNA2D3 chr3:54837229-54838566,
  PCBP2 chr12:53468808-53471613, L1CAM, CLU, HSPA5.

## PLAN (written before compute)
1. SUBSAMPLE each FASTQ to ~30M read pairs (seqtk sample, fixed seed) for feasibility.
2. Get MOTIF-AGNOSTIC short-read junctions: run BBMap (canonical human params: maxindel=200000
   pairlen=200000 intronlen=20 ambig=best) AND STAR 2-pass, extracting junctions from the RAW BAM
   N-cigar split reads (NOT SJ.out.tab which motif-filters). Junction = (chrom, donor, acceptor)
   with split-read count, per sample.
3. Build the short-read NON-CANONICAL junction TRUTH: junctions in expressed loci with >=K split
   reads in >=2 samples, classify canonical(GT-AG/CT-AC) vs non-canonical.
4. RECALL: for each short-read-confirmed NON-CANONICAL junction, check the long-read data — does raw
   minimap2 FLATTEN it (place at a nearby canonical decoy) and does arm-B (panel_deep) RECOVER it?
   Recall = recovered / (short-read-confirmed non-canonical that minimap2 flattens).
5. REVERSE / non-circular support-check: of the long-read re-placer's non-canonical calls (sweep-B
   leads), how many have motif-agnostic short-read support (resolves the Snaptron circularity)?
6. Calibrated controls: annotated junctions in expressed loci (must have high short-read support);
   SMN locus + GTF2H2/SMN2-exon7.

## ★ REVISED PLAN (advisor 2026-07-11) — key corrections
1. **Recall = raw-mm2 (sub.bam) AND arm-B (refined.bam) vs the SAME short-read truth set, side by side.**
   Measure raw-mm2's recall FIRST (one sample). The mission's "recall over the flattened subset" = arm-B's
   recall among junctions raw-mm2 missed — no separate flattening classifier needed. This directly extends
   the spike-in's raw-mm2-vs-arm-B table to REAL data.
   - raw-mm2 recall ≈ arm-B recall (both high) → real data behaves like sim, mm2 doesn't flatten (honest null).
   - raw-mm2 recall LOW, arm-B recall HIGH → the real-flattening regime the sim couldn't produce = the WIN.
2. **Cohort mismatch (iPSC-MN short-read vs Sumner long-read) is a hard confound on forward recall.** Restrict
   the recall DENOMINATOR to junctions with actual long-read coverage at the locus. Count the intersection
   (short-read non-canonical truth ∩ Sumner covered loci) as an EARLY feasibility gate. If small → forward
   recall underpowered → weight shifts to the reverse per-locus check.
3. **Reverse per-locus check = the most defensible deliverable (lead with it).** For each arm-B non-canonical
   call (revealed_noncanon + the sweep-B leads): does motif-agnostic short-read support arm-B's EXACT coord,
   or a CANONICAL neighbor within ±30bp drift window, or NEITHER? exact→recovery (STAR/Snaptron FN, the win);
   canonical-neighbor→drift (the ~1.2% fabrication, now non-circular); neither→inconclusive (low power).
4. **BBMap = PRIMARY motif-agnostic oracle** (canonical human params). STAR's --scoreGapNoncan (−8) penalizes
   non-canonical gaps in the ALIGNMENT itself → its raw BAM still inherits motif bias → STAR = SECONDARY
   corroboration only (or relax --scoreGapNoncan). Confidence gate: ≥2-aligner + ≥K-reads + ≥2-samples.
5. refined.bam is arm-B WITH guard (hp_drift_margin=3.0) = arm-E. For non-HP junctions arm-E==arm-B (spike-in),
   so a valid arm-B proxy — but NO arm-A / unguarded-arm-B on real data from panel_deep; don't claim arm-A-vs-B.

## CHECKPOINT 0 (2026-07-11) — data + assembly confirmed
- All 8 short-read FASTQ present, .dl_rc=0 (SRR6376956/58 CTRL, 60/62 SMA; 13-19GB each). ✓
- panel_deep: 15 .sub.bam (raw-mm2) + 13 .refined.bam (arm-B/E) + revealed_noncanon.tsv per sample. ✓
  revealed_noncanon.tsv = (chrom, donor, acceptor, sample) = ref_nn − raw_nn (arm-B non-canonical NOT in raw).
  = EXACTLY "what arm-B places non-canonically that raw-mm2 flattened." SMA_191 revealed_novel_noncanon=186,966.
- **★ ASSEMBLY IDENTICAL:** short-read genome (compass GRCh38_gencode_v44.fasta) and long-read genome
  (GRCh38.primary_assembly.genome.fa) are the SAME assembly — 194 contigs each, identical name+length md5
  (f1225b925ccff466400c4106aeeeb071), empty diff, chr5 len 181538259 both. Coordinates line up EXACTLY.
  Removes the coordinate-harmonization confound. Both chr-named + unplaced (GL000195.1 etc.) present.
- Full-genome STAR index + bbmap index BOTH pre-built at compass_a549/COMPASS/genome_references/. ✓
- pipeline: sub.bam=raw-mm2 downsampled(frac 0.25); refined.bam=refine_bam_junctions motif_blind=True
  hp_drift_margin=3.0 (=arm-E). Discovery genome-wide per sample (sumner_gw_discover.py).

## CHECKPOINT 1 (2026-07-11) — BBMap recall alignment SUBMITTED
- sbatch job **33531335** (4-task array, larsms, 16c/64G/12h): SRR6376960(SMA) SRR6376962(SMA)
  SRR6376956(CTRL) SRR6376958(CTRL). Per task: reformat.sh subsample→30M PAIRS (sampleseed=42) →
  BBMap (maxindel=200000 pairlen=200000 intronlen=20, ambig=best default, canonical human COMPASS
  params) → sort/index → extract_sr_junctions.py (motif-agnostic N-cigar split junctions, min-anchor 8,
  mapq 10, classify canonical GT-AG/CT-AC).
- Scripts: /scratch/users/kevinroy/sma_recall/{extract_sr_junctions.py, sma_recall_bbmap.sbatch, samples.txt}
- Sentinel: /scratch/users/kevinroy/sma_recall/.bbmap_done (one line per finished SRR).
- RESUME: ssh sherlock 'cat /scratch/users/kevinroy/sma_recall/.bbmap_done; ls -la /scratch/users/kevinroy/
  sma_recall/junc/*.junc.tsv'. 4 lines in .bbmap_done + 4 junc TSVs = done. If FAILED, tail logs/bbmap-*.log.
- Read length PE101, STAR idx sjdbOverhang 149 (fine). BBMap primary (motif-agnostic); STAR-raw = optional
  secondary (its --scoreGapNoncan biases even raw BAM; deferred unless BBMap alone is underpowered).

## CHECKPOINT 1b (2026-07-11) — set -u trap fixed, RESUBMITTED
- First submit (33531335/33531405) FAILED in 12-16s: conda env openjdk activate.d script references
  unbound `target_platform` under `set -u` (documented Sherlock trap). FIX: conda activate BEFORE `set -u`
  + `export target_platform=linux-64`. RESUBMITTED: BBMap **33531828**, LR-extract **33531830**.
- RESUME: ssh sherlock 'cat /scratch/users/kevinroy/sma_recall/.bbmap_done .lrjunc_done; sacct -j
  33531828,33531830 -X -o JobID,State,Elapsed,ExitCode'.

## CHECKPOINT 2 (2026-07-11) — LR junction extraction DONE; mechanism confirmed on REAL data
Job 33531830 COMPLETED (15 samples, sub + refined). My extractor (min-anchor 5, mapq 0, distinct junction
keys) on SMA_191:
  raw-mm2 (sub.bam)     : 200,070 canonical / 8,836 non-canonical (total 222,993 distinct junctions)
  arm-B   (refined.bam) : 195,312 canonical / 13,614 non-canonical (total 223,013)
→ arm-B converts ~4,758 canonical → ~4,778 non-canonical junctions = the "revealed" set = raw-mm2-flattens
  → arm-B-reveals mechanism, REPRODUCED on real data. (Counts differ slightly from discover script's per-read
  revealed_novel_noncanon=186,966 because I count DISTINCT junction keys + min-anchor 5, not per-read intron
  instances — direction identical.) The recall analysis now asks: of these revealed non-canonical junctions,
  which have motif-agnostic short-read support (recovery) vs a canonical neighbor (drift/fabrication)?
BBMap short-read (33531828): tasks 1-2 RUNNING, 3-4 pending (%4). Waiter birs04xok watching both terminal.

## CHECKPOINT 2b (2026-07-11) — min-anchor fix (advisor); re-run LR at anchor 0
UBA1 lead (chrX:47211252-47212397) is in 5 revealed_noncanon.tsv samples but was MISSING from my anchor-5
refined.junc extraction. CIGAR inspection (samtools view SMA_191.refined.bam chrX:47211250-47211255) showed
arm-B places it as `...72M3I119N3D190M...` — the N-op is flanked by 3I (left) / 3D (right) from RECTIFY's
CIGAR surgery. My extractor's CONTIGUOUS-anchor logic ({M,=,X} only) counts anchor=0 when an I/D abuts the
N, so anchor-5 systematically DROPPED arm-B's surgically-placed junctions (biasing forward recall AGAINST
arm-B). This is an EXTRACTOR ARTIFACT, not a real short-anchor finding — do NOT record "arm-B short anchor".
FIX (advisor): re-run LR extraction at --min-anchor 0 --mapq 0 (matches the discovery pipeline read_junctions,
which used no anchor filter to build revealed_noncanon). Re-run job **33532684**. KEEP --min-anchor 8 on the
SR/BBMap side (a 2bp short-read overhang over "non-canonical" is misalignment — the asymmetry is correct).
Also fixed classify_call: "drift" now requires a WELL-SUPPORTED canonical neighbor (>=min_split, >=min_samples),
and "inconclusive/neither" is explicitly UNMEASURABLE (cohort mismatch), NOT fabrication.
Reverse check (4b + leads) reads revealed_noncanon.tsv DIRECTLY (anchor-independent) = the primary deliverable.
VERIFIED: after re-run (33532684, anchor 0), SMA_191.refined.junc.tsv NOW contains UBA1 chrX:47211252-47212397
(count=1) — the fix recovers arm-B's surgically-placed junctions. LR redo COMPLETED 6/13, rest running.

## CHECKPOINT 2c (2026-07-11) — LR redo DONE; BBMap in flight; analysis auto-armed
LR extraction redo (anchor 0) COMPLETED all 15 tasks (.lrjunc_done=15). BBMap short-read (33531828):
tasks 1-2 RUNNING (~1h+, full-genome align of 30M pairs = slow), 3-4 pending. Background waiter b26ql4dwp
armed to AUTO-RUN recall_analyze.py --min-split 3 --min-samples 2 the moment BBMap goes terminal, writing
/scratch/users/kevinroy/sma_recall/RECALL_ANALYSIS.txt. On notification: paste its 4 sections here + verdict.

DRY-RUN validated the analysis end-to-end (rc=0, SR dir empty → all inconclusive, no crashes). STRUCTURAL
numbers already available (LR-only, 13 paired sub+refined samples):
  - LR distinct junction keys: raw-mm2 (sub) = 841,661 ; arm-B (refined) = 1,498,361 (arm-B places ~78% MORE
    distinct junction positions — the re-placer spreads reads across more coords, incl. non-canonical).
  - arm-B REVEALED non-canonical junctions (the reverse-check denominator): 142,498 in >=2 LR samples,
    56,891 in >=3, 17,728 in >=5.
  - 4 sweep-B leads (UBA1/SNRPN/CACNA2D3/PCBP2) all render, currently INCONCLUSIVE (awaiting SR data).
  When BBMap lands, the SR truth set populates recovery/drift/inconclusive.

## CHECKPOINTS (append-only; every number persisted the moment it lands)
- [x] data confirmed present on Sherlock + assembly identity verified (CHECKPOINT 0)
- [x] LR junction extraction (raw-mm2 + arm-B), anchor 0, all 15 samples — mechanism confirmed (CHECKPOINT 2/2c)
- [ ] subsampled ~30M pairs (1 sample first: SRR6376960 SMA)
- [ ] BBMap aligned (PRIMARY motif-agnostic) — split-read junctions extracted
- [ ] STAR-raw aligned (SECONDARY) — N-cigar junctions from raw BAM
- [ ] short-read non-canonical truth set (>=K reads, >=2 aligners; multi-sample after)
- [ ] FEASIBILITY GATE: |short-read non-canonical truth ∩ Sumner long-read covered loci|
- [ ] RECALL: raw-mm2 (sub.bam) vs arm-B (refined.bam) recall on the SAME truth set (side by side)
- [ ] REVERSE: sweep-B lead support (UBA1/SNRPN/CACNA2D3/PCBP2) — exact / canonical-neighbor / neither
- [ ] VERDICT + caveats

## ★★★ CHECKPOINT 3 (2026-07-11) — RECALL LANDED (§1–3 durable; §4b OOM, re-running w/ mem)
BBMap 33531828 COMPLETED (4/4 `*.bbmap.junc.tsv`); LR redo done (13 paired). Ran
`recall_analyze.py --min-split 3 --min-samples 2` (guard code `PYTHONPATH=/scratch/.../rectify_guard`).
Login-node memory cgroup SIGKILLed at §4b (exit 137) AFTER §1–3 printed. stdout→`recall_result.txt`.

**§1 SHORT-READ TRUTH (BBMap, motif-agnostic):** 165,140 confirmed (≥3 split, ≥2/4 samp) =
53,320 canonical + **111,820 NON-CANONICAL**. ⚠ 68% non-canonical → truth denominator is
noise-inflated (motif-agnostic SR admits split-align noise) → absolute recall UNDERCOUNTS; the
arm-B-vs-raw COMPARISON is the robust signal, not the absolute level.

**§2 FEASIBILITY GATE:** SR non-canon 111,820; with LR coverage 108,642 (coverage NOT the limiter).
LR junction keys: raw 841,661 ; armB 1,498,361.

**★ §3 RECALL (raw-mm2 sub.bam vs arm-B refined.bam, SAME truth set):**
| truth | raw-mm2 | arm-B | note |
| --- | --- | --- | --- |
| ALL SR non-canon (n=111,820) | 603 = **0.54%** | 19,522 = **17.46%** | arm-B recovers 18,941 raw misses |
| LR-covered non-canon (n=108,642) | 603 = 0.56% | 19,522 = **17.97%** | |
| SR CANONICAL control (n=53,320) | 50,459 = 94.63% | 50,339 = **94.41%** | arm-B holds canonical (−0.22%) |

**HEADLINE (non-circular WIN):** the motif-blind re-placer lifts non-canonical recall against an
INDEPENDENT motif-agnostic short-read truth **~32× (0.54%→17.46%)** while holding canonical flat
(94.63%→94.41%). Raw mm2 flattens 99.5% of the real non-canonical junctions an independent SR method
sees; the re-placer recovers ~18.9k. This is the real-flattening regime the sim structurally could NOT
produce (spike-in note) — the mission's WIN condition (raw LOW, arm-B HIGH), now on real data.
Caveats: absolute 17% capped by (a) iPSC-MN-SR vs Sumner-LR cohort mismatch, (b) noise-inflated §1
denominator. arm-B here is WITH the hp_drift_margin=3.0 guard (=arm-E); microhom guard NOT yet applied
(its removed drift is largely absent from SR truth anyway, so guard-on recall on this set ≈ preserved).

**§4b PENDING (the guard-threshold confirmation):** the AGGREGATE reverse check — classify arm-B's
recurrent revealed non-canonical junctions as recovery/drift/inconclusive — is what measures the
real-data FABRICATION rate non-circularly and confirms the microhom 0.5/8.0 operating point. OOM'd
before printing. RE-RUNNING on a compute node w/ memory (see RESUME §mem below).

## ⚠ §4b UPDATE (2026-07-11): job 33632262 TIMEOUT at 30min (mem was fine, TIME was not)
Root cause = ALGORITHMIC, not memory: `classify_call` linear-scans `sr_by_chrom.get(chrom, [])` for
EACH of ~142k recurrent revealed junctions (minrec=2) → O(142k × per-chrom-list-len) ≈ billions of ops.
§1–3 completed and are durable (above); §4/§4b never printed. **§4b is NO LONGER BLOCKING** — the audit
verdict is HOLD regardless (guard default stays OFF on the read-blind fault + incomplete-audit grounds,
`dev/MICROHOM_AUDIT_SYNTHESIS.md`), so the §4b real-data fabrication estimate is now context, not a gate.
**FIX before re-running §4b:** in `recall_analyze.py` replace the linear neighbor scan with a position-
binned index — build `sr_canon_by_bin[(chrom, pos//30)]` once (exact match → set lookup; ±30bp neighbor
→ check the 3 adjacent bins), then classify_call is O(1) per junction. Then `sbatch run_4b.sbatch`
(bump `--time=01:00:00` as belt-and-suspenders). Only pursue if we resume guard remediation.

## ★★★★ §4b RESULT (2026-07-14, job 33933333, binned-index fix `recall_analyze_fast.py`) — THE GATE NUMBER
The non-circular real-data fabrication rate on arm-E (`refined.bam`, HP-guard on / microhom-guard OFF)
revealed non-canonical junctions, classified vs independent motif-agnostic short-read (BBMap):
```
  recurrence   total     recovery(win)   DRIFT(fabrication)   inconclusive(unmeasurable)
  >=2 samp    142,498    0.5% (669)      27.3% (38,855)       72.3% (102,974)
  >=3 samp     56,891    0.6% (339)      26.6% (15,137)       72.8%
  >=5 samp     17,728    0.6% (102)      23.6% (4,182)        75.8%
```
**INTERPRETATION (careful):** among the ~25% of arm-E non-canonical calls that short-read CAN adjudicate
(recovery+drift), **~98% are DRIFT (fabrication)** vs ~2% recovery — drift:recovery ≈ 58:1. So on real SMA
data, guard-OFF non-canonical calls are overwhelmingly fabrication *where measurable*. **Fabrication is REAL
and substantial** (≫ the spike-in ~1.2%) → the guard track is JUSTIFIED, NOT solving a non-problem.
**CAVEATS:** (1) 72–76% are INCONCLUSIVE — no short-read signal either way (iPSC-MN-SR vs Sumner-LR cohort
mismatch); the 27% drift is only over the measurable minority, may not generalize. (2) "drift" = well-supported
canonical SR neighbor ≤±30bp + no exact SR support ⇒ parsimonious fabrication, not proof. (3) This is the
guard-OFF arm; it does NOT yet show the guard REMOVES this drift without hurting recall.
**LEADS:** CACNA2D3 = DRIFT (canonical 2bp away, in 2 samp — clear fabrication, in-window). UBA1/SNRPN/PCBP2 =
INCONCLUSIVE, and tellingly their nearest well-supported neighbors sit at **53 / 16 / 48 bp** — i.e. drift
distances that EXCEED the positional close's W=28 horizon (V5 finding). So real fabrication spans past W →
the current close would NOT catch much of it → the SCORER-LEVEL W-fix is the needed next work IF we finish the guard.
**⇒ DECISION-RELEVANT:** fabrication real → guard justified; current close incomplete at real drift distances
(>W=28) → scorer-level fix needed; DEFINITIVE confirmation = a close-ENABLED (guard-ON) COMPASS run measuring
fab-removed AND recall-held (the follow-up).

## ★ RESUME (mem) — finish §4/§4b on a compute node
**SUBMITTED: sbatch job 33632262** (32G, larsms, writes `recall_result.txt` + sentinel
`/scratch/users/kevinroy/sma_recall/.recall4b_rc`) → TIMED OUT (see §4b UPDATE above; needs the
binned-index fix, not just more mem/time). Poll:
`ssh sherlock 'cat /scratch/users/kevinroy/sma_recall/.recall4b_rc 2>/dev/null; sacct -j 33632262 -X -n -o State'`
— if `.recall4b_rc` == 0 → `cat /scratch/users/kevinroy/sma_recall/recall_result.txt` and read §4b's
recovery/drift/inconclusive split. If FAILED/OOM again → bump `--mem=64G` in `run_4b.sbatch`, resubmit.
Manual fallback (if the job route is unavailable):
`ssh sherlock` then:
`srun --account=larsms --partition=batch --mem=32G --time=30 --pty bash -c 'source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh; conda activate rectify; export PYTHONPATH=/scratch/users/kevinroy/rectify_guard; python -u /scratch/users/kevinroy/sma_recall/recall_analyze.py --min-split 3 --min-samples 2 > /scratch/users/kevinroy/sma_recall/recall_result.txt 2>&1'`
then `cat /scratch/users/kevinroy/sma_recall/recall_result.txt`. Read §4b's recovery/drift/inconclusive
split for the recurrent revealed set (≥2/≥3/≥5 samples) — drift% = the non-circular real-data
fabrication estimate to compare against the spike-in 1.31% and confirm the guard threshold.

## ★ RESUME (self-contained — run when BBMap 33531828 + LR-redo 33532684 both COMPLETED)
1. Check: `ssh sherlock 'cat /scratch/users/kevinroy/sma_recall/.bbmap_done; wc -l < /scratch/users/kevinroy/sma_recall/.lrjunc_done; ls /scratch/users/kevinroy/sma_recall/junc/'`
   Need 4 `*.bbmap.junc.tsv` in junc/ + 13 lines in .lrjunc_done. If a BBMap task FAILED: `tail
   /scratch/users/kevinroy/sma_recall/logs/bbmap-33531828_<task>.log` (reformat OOM / bbmap path); resubmit
   single task `sbatch --array=<N> sma_recall_bbmap.sbatch`.
2. Run the analysis (captures ALL numbers to the .md):
   `ssh sherlock 'source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh; conda activate compass;
    python /scratch/users/kevinroy/sma_recall/recall_analyze.py --min-split 3 --min-samples 2'`
   Also sweep --min-split 2 and 5 for robustness.
3. Paste the 4 sections (truth set / feasibility gate / recall raw-mm2-vs-arm-B / reverse leads+aggregate) into
   the checkpoints below, then write the VERDICT. LEAD with the reverse per-lead + aggregate (non-circular,
   cohort-robust); forward recall is the cohort-confounded secondary.
