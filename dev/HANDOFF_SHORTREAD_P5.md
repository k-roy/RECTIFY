════════════════════════════════════════════════════════════════════════════════
# ★★ CURRENT HANDOFF (2026-06-26 PM) — 2-corroborated CROSS-PLATFORM follow-up + de-novo-aligner Q
(Newest. Supersedes the DELIVERED block below for current state. The 111 deliverable itself is still
DELIVERED & locked on Oak — unchanged.)

## DONE (this session)
- **Characterized the 2 corroborated junctions** (forward NEXT item #1) into a verdict-ready packet,
  reusing the live rectify lib on-cluster. Tool: `dev/compass_corroborated_junctions_characterize.py`
  (deployed `$W/rectify_src/dev/`; out `$W/corroborated_2_characterization.json`). SELF-CHECK PASSED
  (recomputed short-read depths = 12 and 3 = locked json). Cross-platform table doc:
  `dev/COMPASS_2corroborated_CROSSPLATFORM.md`. Long-read probe: `dev/longread_probe.py`
  (out `$W/longread_probe_out.txt`).
  - **J1 chr5:140564954-140565547 (593bp, SLC35A4 +) → REAL novel junction (high confidence).** Canonical
    GT..AG, no annotated junction in-locus. Corroborated by minimap2(101)+deSALT(102)+uLTRA(81)+GMAP+
    COMPASS short(168) → **NOT actually gmap-only**: the "111 gmap-only" label was a 3bp ambiguity-
    normalization gap (others place the SAME 593bp intron 3bp shifted at 140564957-550, which strict
    base-equality normalize() didn't merge; adjudication depth 12 UNDERCOUNTS — true short-read ≈168).
  - **J2 chr5:179823051-179823857 (806bp, SQSTM1 +) → INCONCLUSIVE (leans not-confirmed).** A +2bp
    alt-acceptor off the ANNOTATED 804bp SQSTM1 intron (same donor 179823053; acceptor 859 vs 857).
    Short reads: 857=2258 sharp, 859=3 (0.13%, noise-level). Long reads SMEAR acceptor 855–862 (no clean
    859 peak) → platforms do NOT agree on a specific ~0.13% alt-acceptor freq (Kevin's cross-platform
    frequency test, the RIGHT test, is NOT satisfied). Canonical but not provably real nor a seq error.
- **Answered Kevin's de-novo-aligner "where were we" Q:** the idea (realign within the panel-mapped locus,
  scored by empirical ED tables) = **facet C1 of the RECTIFY native aligner member**. Status: DESIGNED, not
  built; build gated on a simulation benchmark that doesn't exist yet. Docs: `dev/ALIGNER_PROGRAM_SCOPING.md`
  + `dev/ALIGNER_MEMBER_DESIGN.md` (both 2026-06-18). Tables+realign DP both exist in prod but NOT wired
  (`realign_exon_blocks` read_edits.py:791 → `align_exon_block_global` local_aligner.py:522 still uses FLAT
  costs; TODO at local_aligner.py:37-60). Next concrete step = C1 Phase 0 equivalence harness, AFTER the
  benchmark gate.

## ⚠ RECOUNT RESULT (DONE) — the "109/111 artifacts" headline is NOT robust as stated
`dev/recount_111.py` (out `$W/recount_111_out.txt`): EXACT-normalized∩COMPASS = 2 (reproduces adjudication);
**RELAXED (same intron length, |Δstart|≤5)∩COMPASS = 4.** Two of the "109 artifacts" actually have STRONG
independent short-read support at a same-length ≤5bp-shifted coord the exact match missed:
  • `chr5:179824400-179832205` (7805bp, gmap 12) → **2959** COMPASS reads
  • `chr5:177592500-177593474` (974bp, gmap 23) → **323** COMPASS reads
→ corroboration is ≥4/111 (≥107 artifacts), AND the adjudication's exact-coordinate match SYSTEMATICALLY
undercounts (the 3bp normalization gap, same one that split J1). **Methodological fix needed**: length-aware /
windowed merge of same-length placements BEFORE the gmap-only test + motif-check AT the supported coordinate.
- **J1 REINTERPRETED (downgrade — earlier "clean real novel junction" was WRONG):** the bulk J1 support is at
  the NON-canonical placement 140564957-550 (donor AG, acceptor CA — NOT GT-AG, verified from genome), 156
  short + minimap2 101 + deSALT 102 + uLTRA 81; only the MINOR placement 140564954-547 is canonical GT-AG (12
  short + 8 GMAP). The two are different spliced products 3bp apart (GTG≠TCA, not an ambiguity slide). 4
  aligners + Illumina uniformly placing a 593bp gap at a non-canonical breakpoint ⇒ more consistent with a
  recurrent ~593bp GENOMIC DELETION/SV in A549 (CNV-rich line) or shared alignment artifact than splicing.
  J1 is NOT a confirmed canonical novel junction. (J2 verdict unchanged: inconclusive.)

## DE-NOVO ALIGNER (parallel track) — benchmark GATE built by dedicated agent ✅
User chose "simulation benchmark (the gate)" as first deliverable. A worktree-isolated agent BUILT +
validated + committed Deliverable A on branch **`worktree-agent-a25a2c1e784ad37dc`** (off
`drs-validation-rebuild`; NOT committed to the shared branch; all under benchmark-only paths
`rectify/core/benchmark/`, `scripts/benchmark/`, `dev/`). Scoped brief + the agent's own RESUME:
`dev/HANDOFF_ALIGNER_BENCHMARK.md`. Smoke green (1bp-shift→TP-not-FP; indel concordance 0.997; the live
flat-affine DP scores 0.980 on HP_HARD-noisy with 24/24 failures = the indel-vs-sub misplacement C1 targets
→ the gate CAN discriminate the member). **OPEN:** pbsim3 live run not executed (bioconda solve was slow on
Sherlock; wrapper code-complete, MAF→genome validated locally). **NEXT cycle (now unblocked, gated on this):**
prove length-law arm > flat arm on HP_HARD-noisy (harness ready); THEN C1 member code. For Kevin: review/merge
the worktree branch when ready.

## #2 circRNA check — DONE (2026-06-28): NOT circRNA (mechanism = non-canonical LINEAR splicing)
`dev/circ_check.py`/`circ_check2.py` on the long-read per-aligner BAMs. Across minimap2/deSALT/uLTRA/GMAP at
all 3 loci: **rolling-circle = 0, head-to-tail back-splice = 0.** Direct read inspection confirms FORWARD
N-op orientation (deSALT read N-op 179824391→179833031, donor<acceptor = linear). → forward linear RNA gaps
over present DNA = **non-canonical LINEAR splicing**, not back-splicing. Wrinkle: minimap2 renders the 7805
locus as 6447 CHIMERIC (SA) reads split to OTHER loci (not back into locus) → the 7805 reads may be part of
complex/chimeric transcripts; deSALT/uLTRA/GMAP call a clean forward intron. (forward-gap recount in
circ_check2 had a tolerance bug, but collect_junction_counts already proved the N-ops + direct dump confirms.)

## #2b FUSION/CHIMERA check — DONE (2026-06-28): NOT a distant fusion; all 3 INTRAGENIC
`dev/chimera_invest.py`: the 7805 locus's 6446 minimap2-chimeric reads map their SA segments back into SQSTM1
itself (chr5:179.8-179.9Mb) — local multi-exon read splitting, NO distant partner. All 3 are intragenic,
sense strand: 7805∈SQSTM1(+), 974∈TMED9(+), J1-593∈SLC35A4(+) (SQSTM1 hosts TWO incl. J2). Each dominant
NON-canonical coord has a canonical GT-AG/GC-AG 1-4bp away (ambiguity 0 = distinct placement). Accurate
SHORT reads place at the NON-canonical coord → leans genuinely-non-canonical (vs canonical-misplaced-by-
breakpoint-errors; resolve via motif-aware realign = the C1 feature / RT-PCR). → 3 = REAL intragenic novel
RNA junctions, DNA-confirmed, not deletions/circRNA/fusion.

## #1 v3 json — STAGED (generator ready, not yet locked): `$W/make_v3.py`
Run AFTER #3: `python $W/make_v3.py "<het-result-string>"` → writes `$W/adjudication_111_v3.json` (full
method chain + per-junction investigation + corrected summary: 107 artifacts / 3 real intragenic RNA
junctions / 1 inconclusive / 0 deletions/circRNA/fusion). Then copy to Oak + chmod 444 (chmod u+w $OAK first,
re-lock dir after) — same recipe as v2. Supersedes v1 ("109/2") AND v2 (STRUCTURAL_DELETION).

## ✅ GMAP TRACK FINALIZED (2026-07-02) — verdict corrected + consistent
Reval done (job 32220563): corrected (ambiguity-tolerant) gmap-only test → **201 candidates / 5 SUPPORTED /
~1 well-supported** across 5 chroms (vs inflated exact-match 349/13). VERIFIED KEY POINT: the 3 "real
intragenic junctions" (SLC35A4/TMED9/SQSTM1) are EXCLUDED from the corrected gmap-only set → they're real but
MULTI-ALIGNER, NOT gmap-unique. So GMAP's genuinely-unique well-supported yield ≈ 1/5-chrom (chr1:19219782).
**Final verdict: GMAP appropriate but LOW-unique-value; keep fenced; well-supported DROP candidate once C1
lands.** Fully consolidated in `dev/GMAP_PANEL_APPROPRIATENESS.md` (§F + rewritten VERDICT/RECOMMENDATION;
§D/§E numbers marked superseded). Files: `$W/reval_chr{5,1,11,17,19}.json`, `$W/agg_reval.py`, `$W/confirm_13.json`.
- OPTIONAL follow-ups (low priority, none blocking): (a) genome-wide confirm — harness ready (`--chrom all`)
  but needs a genome-wide GMAP+indep LR panel; 5-chrom is already decision-grade → low marginal value.
  (b) note in `adjudication_111_v3.json` that the 3 real junctions are multi-aligner not gmap-unique (reality
  unchanged, only gmap-credit). (c) C1 realigner on the non-canonical junctions (inbox note already sent).

## (DONE, see above) HARNESS FIX + RE-VALIDATION — job 32220563 (2026-06-30)
Fixed `dev/gmap_validate_harness.py` STEP-1: "others" (gmap-only) test is now AMBIGUITY-TOLERANT (exact OR
same-intron-length within ±win), not exact-normalized — closes the over-count the confirm exposed. Re-running
5 chroms (chr5 full panel + chr1/11/17/19 morechrom): `$W/reval5.sbatch` → `$W/reval_chr{5,1,11,17,19}.json`,
sentinel `$W/.gmap_reval5_rc`. EXPECT materially FEWER gmap-only candidates than the old 349 (the exact-match
inflated it; the truly-gmap-only 13→~3 pattern suggests a big drop).
- **RESUME:** `ssh sherlock 'cat /scratch/users/kevinroy/compass_a549/.gmap_reval5_rc; for c in chr5 chr1 chr11 chr17 chr19; do python3 -c "import json;d=json.load(open(\"/scratch/users/kevinroy/compass_a549/reval_$c.json\"));print(\"$c\",d[\"TOTAL\"] if \"TOTAL\" in d else d[\"per_chrom_summary\"])"; done'`
  - sentinel `0` → aggregate corrected candidates→SUPPORTED across 5 chroms; update `GMAP_PANEL_APPROPRIATENESS.md`
    §D/§E with the corrected (lower) gmap-only count. Job gone & no sentinel → `sacct -j 32220563 -X -o State`;
    re-`sbatch $W/reval5.sbatch`.

## ✅ 13-SUPPORTED CONFIRM — DONE (2026-06-30): none are SV; harness "gmap-only" is inflated
Ran WGS + reliable cross-aligner confirm on all 13 SUPPORTED (`$W/confirm_13.json`, `dev/confirm13*.py`).
- **WGS: 13/13 present normal-copy (ratio 0.72-1.33) → NONE a deletion/SV** → all RNA-level (splicing).
- **Cross-aligner (collect-based, the reliable method; manual-walk attempt was buggy): 10/13 are actually
  MULTI-aligner** (≥2 of minimap2/deSALT/uLTRA place them, a few bp off GMAP's coord) → the harness STEP-1
  "gmap-only" exact-match is INFLATED by the same normalization gap as the 111. **Genuinely gmap-only = 3/13**;
  of those only chr1:19219782 (sr99, non-canon) is well-supported (GMAP uniquely finding a real short-read-
  corroborated junction) — chr11:308314 (sr1), chr19:48965691 (sr2) are thin. GMAP's genuinely-unique
  well-corroborated yield across 5 chroms ≈ **1-2 junctions**. Folded into `GMAP_PANEL_APPROPRIATENESS.md` §E.
- **TWO FIXES for the harness before genome-wide run:** (1) apply ambiguity-tolerant (length-matched) matching
  to STEP-1's "others" set (currently exact-normalized → over-counts gmap-only); (2) the confirm's
  collect-cache OOMs if it holds many full-chrom dicts — free per-chrom (see `dev/confirm_chr19.py` lean pattern).

## ✅ GMAP-VALIDATION 5-CHROM — DONE (job 31903808 COMPLETED, 2026-06-30)
Harness ran chr5+chr1/11/17/19 (GMAP-2021 added to DRS morechrom panel, 999,228 primary aligns).
**RESULT: 349 gmap-only-recurrent-canonical-novel candidates → 13 short-read-SUPPORTED (~7 with real depth
sr≥30); well-supported calls mostly NON-canonical/complex-locus; ~0-1 well-supported canonical unique novel
per chromosome — chr5's "3-4" was the high end, pattern generalizes.** Folded into
`dev/GMAP_PANEL_APPROPRIATENESS.md` §D + verdict (GMAP = LOW-unique-value, kept because fully fenced; the
"needs genome-wide count" caveat is now substantially answered). Per-chrom jsons `$W/gmap_validate_chr{1,11,17,19}.json`;
aggregator `$W/agg5.py`; GMAP bam `$W/morechrom_trimmed.GMAP.bam`.
- **OPEN (cheap next step, NOT done):** the 13 SUPPORTED have short-read corroboration but not the full
  cross-aligner+WGS motif confirm (the chr5-grade adjudication). Close by feeding the 13 to `dev/lr_probe_4loci.py`
  + `dev/dna_split.py` + the WGS BAM. Also: extend to remaining autosomes for a true genome-wide count (needs
  GMAP+indep aligner on those chroms — morechrom only had chr1/11/17/19).
- NOTE: the job's `echo 0 > $SENT` sentinel line apparently didn't land (log shows `[done]` + all 4 jsons
  present + State COMPLETED) — result is complete regardless; verified via the jsons, not the sentinel.

## (superseded — see above) GMAP-VALIDATION EXTENDED to chr1/11/17/19 — job 31903808
Reuses the DRS agent's `morechrom` panel (READ-ONLY): minimap2/deSALT/uLTRA aligned across chr1/11/17/19
(coords == GRCh38, confirmed). GMAP was the only missing aligner → `$W/gmap_morechrom.sbatch` builds a
gmap-2021 index for morechrom_ref.fa + aligns GMAP → `$W/morechrom_trimmed.GMAP.bam` (MY workspace, NOT the
DRS dir) → runs `dev/gmap_validate_harness.py` per chrom → `$W/gmap_validate_chr{1,11,17,19}.json`. Sentinel
`$W/.gmap_morechrom_rc`. (DRS coordination note dropped in inbox.) Extends validation chr5 → 5 chromosomes.
- **⚠ SSH (2026-06-29):** the **M1's** Sherlock ControlMaster expired (Kerberos, ~daily) — agent ssh gets
  `Permission denied (gssapi-with-mic)`. The job runs ON Sherlock regardless. To re-check: USER runs
  `! ssh sherlock` **in the Claude session on the M1** + Duo (the iMac's live session does NOT refresh the
  M1 master). Then resume below.
- **RESUME:** `ssh sherlock 'cat /scratch/users/kevinroy/compass_a549/.gmap_morechrom_rc; for c in chr1 chr11 chr17 chr19; do echo $c; python3 -c "import json;d=json.load(open(\"/scratch/users/kevinroy/compass_a549/gmap_validate_$c.json\"));print(d[\"TOTAL\"])"; done'`
  - sentinel `0` → read the 4 per-chrom jsons: TOTAL candidates → SUPPORTED. Aggregate with chr5 (111→4) for
    a 5-chrom GMAP unique-real-novel estimate → update `dev/GMAP_PANEL_APPROPRIATENESS.md` (genome-wide-ward count).
  - `11`=gmap_build failed / `12`=align failed → `$W/gmap_morechrom.<JID>.out`, `$W/morechrom_gmap/gmap_*.log`.
  - sentinel absent & job gone (`sacct -j 31903808 -X -o State`) → re-`sbatch $W/gmap_morechrom.sbatch` (idempotent: skips built index / existing GMAP bam).
- NOTE: morechrom covers only 4 chroms (not full genome); true genome-wide still needs the remaining chroms
  aligned (GMAP + an independent aligner). This is the tractable next increment, not the whole genome.

## ✅ GOAL-2 SYNTHESIS + HARNESS + CASE STUDIES (2026-06-29)
- GMAP-appropriateness verdict: `dev/GMAP_PANEL_APPROPRIATENESS.md` — GMAP appropriate (FENCED); ~97% of its
  unique-junction output non-canonical; of 111 chr5 best unique-novel candidates only ~3-4 real. KEEP fenced;
  genome-wide validated count needed before any DROP (advisor corrected an earlier "marginal/30x-overstated"
  overreach — chr5 3-4 scales to ~tens genome-wide; fence-test "~100 genuine" is an unvalidated chr5 candidate count).
- **VALIDATION HARNESS built + self-tested:** `dev/gmap_validate_harness.py` (chrom-agnostic, `--chrom all`):
  candidate-gen → short-read ambiguity-tolerant validation → motif → (optional WGS). chr5 self-test reproduces
  the hand analysis EXACTLY (111 → 107 ARTIFACT / 4 SUPPORTED from raw BAMs). Out: `$W/gmap_validate_chr5.json`.
  **Genome-wide run BLOCKED on data:** need genome-wide A549 long-read alignments (GMAP + ≥1 independent
  aligner); current BAMs are chr5-trimmed. Short-read consensus (Oak) + WGS are genome-wide.
- Case studies: added "human non-canonical junction discovery (A549)" section + acceptance criteria to
  `rectify/data/validation/CASE_STUDIES.md` (points to the cross-platform doc; yeast footer preserved).
- C1 inbox note: `.claude/inbox/20260629T065724Z__from-compass-111-rna__c1-realign-testcase.md`.

## ✅ #3 + #1 DONE (2026-06-28) — het-SV EXCLUDED; v3 LOCKED on Oak
#3 deep WGS (job 31821023 COMPLETED, ~9-10x, `$W/wgs/a549_wgs_deep.bam`): all 3 interiors normal-copy
(interior/flank 0.99/0.78/1.0, 100% cov ~9-10x) → het-del EXCLUDED. `dev/dna_split.py`: DNA reads cross all
breakpoints continuously (0 soft-clip, 0 split/SA) → NO balanced rearrangement. **Every genomic alternative
ruled out → the 3 are real intragenic RNA splice junctions** (SQSTM1/TMED9/SLC35A4). #1: `$OAK/adjudication_111_v3.json`
(chmod 444, md5 439ec68c…) LOCKED — supersedes v1 ("109/2") + v2 (STRUCTURAL_DELETION). Generator `$W/make_v3.py`.
FINAL: 107 artifacts / 3 real intragenic RNA junctions / 1 inconclusive (J2) / 0 del/circRNA/fusion. ONLY open
science Q: genuinely-non-canonical vs canonical-placed-a-few-bp-off → motif-aware realign (the C1 feature) / RT-PCR.

## (DONE, see above) #3 DEEP WGS (het-SV exclusion) — job 31821023
Existing 1.3x already DISFAVORS het (interiors show no coverage halving: 7805 interior 1.63 ≈ genome 1.52).
Firming up: `$W/wgs_deep.sbatch` PARALLEL byte-range downloads first ~15GB of ENCODE R1 (~10x; single-stream
was 570KB/s bottleneck, ENCODE S3 honors HTTP Range = 206, full R1=31GB) → minimap2 → `$W/wgs/a549_wgs_deep.bam`
→ depth-ratio report `$W/wgs/wgs_deep_coverage.txt`. Sentinel `$W/.wgs_deep_rc`.
- **RESUME:** `ssh sherlock 'cat /scratch/users/kevinroy/compass_a549/.wgs_deep_rc; cat /scratch/users/kevinroy/compass_a549/wgs/wgs_deep_coverage.txt'`
  - sentinel `0` → read report: INTERIOR/flank meandepth ratio ~1.0 = normal copy (no del, het EXCLUDED);
    ~0.5 = het deletion; ~0 = homo (already excluded). Then proceed to #1 (re-lock v3 json).
  - `11`=download failed / `12`=align failed → `$W/wgs/wgs_deep.<JID>.out`, `$W/wgs/mm2_deep.log`.
  - sentinel absent & job gone (`sacct -j 31821023 -X -o State`) → re-`sbatch $W/wgs_deep.sbatch` (idempotent).
- **#1 (LAST, after #3):** re-lock `$OAK/adjudication_111_v3.json` — corrected verdict: 3 SUPPORTED junctions =
  real non-canonical LINEAR splice junctions (DNA-confirmed genomic, not deletions, not circRNA, not artifacts);
  J2 inconclusive alt-acceptor; ~107 artifacts. Supersede v2's STRUCTURAL_DELETION calls.

## ✅ A549 WGS DELETION-CHECK (homozygous) — DONE (2026-06-28): deletion REFUTED → real RNA splicing
RESULT: at ~1.3x A549 WGS (ENCODE ENCSR521ELB, `$W/wgs/a549_wgs_partial.bam`), all 3 candidate interiors
(7805/974/J1-593bp) have NORMAL DNA coverage (%cov 73/68/62 = genome-avg ~71) — NO dropout → **NOT deletions;
the genomic sequence is PRESENT.** Sequence present in DNA + removed in RNA = **splicing**, not deletion, not
artifact. So the 3 SUPPORTED_NONCANONICAL are **real non-canonical RNA junctions** (DNA-confirmed genomic,
multi-platform). REVERSES the "structural deletion" verdict. Full table + caveats in
`dev/COMPASS_2corroborated_CROSSPLATFORM.md` (top). Coverage report `$W/wgs/wgs_partial_coverage.txt`.
- **TODO (open, not started):** (a) `$OAK/adjudication_111_v2.json` STRUCTURAL_DELETION verdicts now WRONG →
  re-lock a v3 (DNA-confirmed-genomic / RNA-spliced / mechanism-open). (b) MECHANISM: non-canonical cis-splice
  vs **circRNA back-splice** (check head-to-tail/out-of-order junction orientation) vs trans-splice. (c) deeper
  WGS via PARALLEL-stream download (single-stream ENCODE/S3 ~570 KB/s was the bottleneck; no aria2c) to exclude
  het-SV; ONT A549 WGS (ENA PRJEB90580/ERR15968645) for DNA-level SV confirmation.
- Pipeline note: `wgs_del_check.sbatch` used `set -o pipefail` which tripped on the truncated-gz `zcat` AFTER
  minimap2+sort succeeded → false rc12; finalize the `.tmp` BAM manually or drop pipefail for that pipe.

## (cancelled) A549 WGS DELETION-CHECK — job 31810428 (too slow)
Testing whether the 3 SUPPORTED_NONCANONICAL "novels" (7805bp, 974bp, J1-593bp) are STRUCTURAL DELETIONS in
A549 vs real splicing — the RNA span-test was non-discriminating (real introns also look empty in mRNA-seq;
positive-control proved it). DNA settles it: a deletion = ~0 DNA coverage in the interval interior; real
genomic sequence (intron) = full DNA coverage.
- **Data:** ENCODE A549 WGS `ENCSR521ELB` (Illumina HiSeq X Ten 151bp; files ENCFF122NPY/ENCFF846WHK), raw
  FASTQ. **Sherlock login node has NO working curl** (binary 7.29 vs libcurl 8.9 mismatch → rc43) but wget &
  python urllib WORK → job downloads first 60M reads/file (~5.8x, first-N of Illumina FASTQ ≈ uniform genome
  sample) via python urllib stream, minimap2 -ax sr → BAM, samtools coverage at the 3 candidate interiors +
  flanks + a real-SQSTM1-intron control + intergenic control. Script `$W/wgs_del_check.sbatch` (idempotent:
  skips existing fastq/BAM on requeue).
- **⚠ DOWNLOAD TOO SLOW:** ENCODE/S3 single-stream from Sherlock = ~570 KB/s, so 60M reads/file ≈ 5h. CANCELLED
  job 31810428; salvaged the partial download `$W/wgs/sub1_partial.fastq.gz` (~11M reads ≈ 0.5x, R1 only) and
  submitted a quick align job **31813174** (`$W/wgs_align_partial.sbatch`) → `$W/wgs/a549_wgs_partial.bam` +
  coverage `$W/wgs/wgs_partial_coverage.txt`, sentinel `$W/.wgs_partial_rc`. 0.5x is enough for the BIG 7805bp
  locus (~26 reads if present vs 0 if deleted) but THIN for 974bp/J1-593bp. For deeper depth: no aria2c → use
  PARALLEL streams (download R1+R2 concurrently, or python byte-range threads) — single-stream curl/wget is the
  bottleneck. **RESUME (partial):** `ssh sherlock 'cat $W/.wgs_partial_rc; cat $W/wgs/wgs_partial_coverage.txt'`.
- (orig full job, cancelled) Sentinel `$W/.wgs_del_rc`; report `$W/wgs/wgs_coverage_report.txt`.
- **RESUME:** `ssh sherlock 'cat /scratch/users/kevinroy/compass_a549/.wgs_del_rc 2>/dev/null; \
  cat /scratch/users/kevinroy/compass_a549/wgs/wgs_coverage_report.txt 2>/dev/null'`
  - sentinel `0` → read report. **Candidate INTERIOR meandepth ≈0 while flanks+real-intron-control are full
    → homozygous DELETION confirmed** (→ these are NOT novel junctions; finalize v2 as structural). Interior
    ≈50% of flank → het deletion. Interior ≈ flank ≈ control → genomic sequence present → reconsider (could be
    real intron after all).
  - sentinel `11`=download failed / `12`=align failed → read `$W/wgs/wgs_del.<JID>.out` + `$W/wgs/minimap2.log`.
  - sentinel absent & job gone (`sacct -j 31810428 -X -o State`) → died; re-`sbatch $W/wgs_del_check.sbatch`
    (idempotent, resumes from downloaded fastqs).
- **Then:** if deletions confirmed, the locked `$OAK/adjudication_111_v2.json` STRUCTURAL_DELETION verdicts are
  DNA-confirmed (drop the "pending WGS" hedge); ONT A549 WGS (ENA PRJEB90580/ERR15968645) is the optional
  long-read SV confirmation.

## ✅ 111 RE-DERIVATION COMPLETE (2026-06-26 PM) — locked v2 on Oak
Done via `dev/rederive_111.py` (ambiguity-tolerant same-len ±5 + motif-at-supported-coord) +
`dev/lr_probe_4loci.py` (cross-aligner long-read placement = the decisive discriminator). Authoritative
writeup: `dev/COMPASS_2corroborated_CROSSPLATFORM.md` (top section). **CORRECTED VERDICT:** ~107/111 clean
GMAP artifacts; the 4 with independent short-read support are NOT validated novel splice junctions — **3 are
recurrent structural-DELETION signatures in A549** (multi-aligner+Illumina uniformly at a non-canonical
breakpoint; splice-aware aligners don't snap to the nearby GT-AG/GC-AG; standout = ~7.8kb deletion INSIDE
SQSTM1, 2959 short + ~400 long reads), **1 inconclusive +2 alt-acceptor (J2)**. **NONE rescued as a clean
novel splice junction.** The original "2 corroborated → likely real" is wrong in interpretation.
- **Corrected json LOCKED on Oak:** `$OAK/adjudication_111_v2.json` (chmod 444, md5 5deba55…, generated by
  `$W/make_v2.py` from `$W/rederive_111.json`). SUPERSEDES `$OAK/adjudication_111.json` (left untouched for
  provenance; its "109/2" headline is corrected by v2). Dir re-locked dr-xr-s (read-only).
- **NEXT (confirmation, for Kevin):** IGV the 7805bp (chr5:179824400-179832205) + 974bp (chr5:177592500-177593474)
  loci — deletion-breakpoint signatures (soft-clips/microhomology) vs clean exon/exon — and overlay an A549
  WGS/CNV track to confirm deletion-vs-splicing. Fix the 2 auto-classifier caveats before reusing
  `rederive_111.py` on other chroms (974 canonical-1bp-off mislabel; J2 alt-SS XOR-on-normalized-coord).
- (carried) Forward NEXT items #2 rep2/rep3, #3 split_command.py rm-rf/wildcard hardening (coordinate w/ DRS
  agent), #4 DRS/cDNA backlog (owned by DRS agent) — unchanged from block below.

## RESUME (concrete)
SSH `sherlock` (never tear down ControlMaster; retry transient sshd serially). `$W=/scratch/users/kevinroy/compass_a549`.
1. **Read the recount result:**
   `ssh sherlock 'cat $W/recount_111_out.txt'`
   - **empty?** still running (or died) — check `pgrep -fl recount_111` on sherlock; if gone & empty, re-run:
     `cd $W/rectify_src && PYTHONPATH=$W/rectify_src python dev/recount_111.py > $W/recount_111_out.txt 2>&1`
     (rectify env; ~2min collect over the 6.9G Oak BAM).
   - **"RELAXED … intersection: 2"** (== EXACT 2) → the **109/111 artifacts headline HOLDS**; J1 was the rare
     normalization-split case. Say so plainly to Kevin; ship the characterization as-is.
   - **RELAXED materially > 2** → "109 artifacts" is an OVERCOUNT (more junctions gain COMPASS support under
     the same-length ±5bp match). Report the new number + the gained-junction list to Kevin BEFORE he trusts 109.
2. **Then** (regardless) the J1/J2 packet + de-novo-aligner answer are already delivered to Kevin in-chat;
   the durable artifacts are the two `dev/COMPASS_*` / `dev/compass_*` files + the on-cluster JSON/txt.

## FILES (this session, all UNCOMMITTED, working tree)
- `dev/compass_corroborated_junctions_characterize.py`, `dev/longread_probe.py`, `dev/recount_111.py`
  (all deployed to `$W/rectify_src/dev/`), `dev/COMPASS_2corroborated_CROSSPLATFORM.md`.
- On cluster: `$W/corroborated_2_characterization.json`, `$W/longread_probe_out.txt`, `$W/recount_111_out.txt`.
- Long-read panel BAMs: `…/rectify_human_validation/sgnex_a549/alignments/a549_chr5_trimmed.{minimap2,deSALT,mapPacBio,GMAP,uLTRA,rectified}.bam`.
- NOTE: shared branch `drs-validation-rebuild` has a concurrent DRS agent — do NOT `git add -A`; stage only
  these explicit `dev/` paths if committing.
════════════════════════════════════════════════════════════════════════════════
# ★ HANDOFF (2026-06-26) — COMPASS short-read 111-adjudication: DELIVERED ✅
(Authoritative current state. Detailed chronological history is below this block.)

## DONE
- **COMPASS short-read pipeline (A549, paired Illumina) built, debugged, validated end-to-end at 7/7
  aligners.** 4 real bugs found+fixed (all in HEAD): env mismatch (run rectify in `rectify` py3.9 env, not
  py3.7 `compass`); `check_aligner_available(None)` TypeError; **consensus >2h hang** → memoized locus-index
  (`daea0fa`); gsnap worker symlinks + magicblast bare-QNAME header strip. Adjudication tool built+validated:
  `dev/compass_shortread_adjudicate_111.py`.
- **/goal all tests pass — MET:** full suite at HEAD `b08f1ac` = **1608 passed, 0 failed, 40 skipped, 1 xfailed**.
- **Aligner-version experiment — DONE:** version plays a NEGLIGIBLE role (6v6 STAR/HISAT2/magicblast/bbmap:
  jaccard 0.998, novel-Jaccard 1.0 @≥2 reads, 111∩=0 both). gsnap-2024 is not a drop-in (CLI change + index
  rebuild + SIGSEGV). `dev/compass_version_compare.py`.
- **chr5 single-aligner prefilter — REFUTED:** permissive STAR misses 43–55% of well-supported novel chr5
  junctions (consensus is magicblast-dominated). All-reads whole-genome alignment is REQUIRED, not optional.
- **DATA-LOSS RCA + FIX:** first full run finished but the generated merge did `rm -rf "$OUTDIR"` after copying
  the BAM to `$OAK_OUT` which DEFAULTED to `$OUTDIR/final` (a subdir) — self-deleting the deliverable. NOT a
  purge/agent/offload-bot (confirmed: 0 deletions in 1036 transcripts; imac-hub RCA agrees). Fixed two ways
  (`4ff87a2` guard + always pass real `--oak-output-dir`) and re-launched to durable Oak.
- **Human DRS/cDNA backlog surveyed + de-duped** (Explore survey was stale; concurrent agent already did
  anchor-gate CLI / gapmm2 pin / cluster_com / Dorado polyA / GMAP). Coordination inbox note dropped.

## VERIFIED
- 7/7 chunk run: all aligners SUCCESS, both mates survive per RN, consensus 187s, N-op junctions present.
- 1608/0 full suite at HEAD in isolated `rectify_src_dev` (no regressions from anyone).
- Fixed merge script generates `Oak out: /oak/...` + the subpath guard (tested before relaunch).

## ★ VERDICT — DELIVERED (2026-06-26, full-depth A549_rep1, durable read-only on Oak)
`$OAK/adjudication_111.json` (chmod 444): **compass chr5 junctions 21,532**; **POSITIVE control PASSED** —
7,941/18,450 annotated chr5 junctions recovered (43.0%); **NEGATIVE control clean** — 0/101 decoys;
**111 ∩ COMPASS = 2/111.** → **109/111 (98.2%) confirmed ARTIFACTS** (no independent short-read support);
**2/111 CORROBORATED** = likely real novel junctions the long-read panel missed:
  • `chr5:140564954-140565547` depth 12 (gmap 8) — solid; • `chr5:179823051-179823857` depth 3 (gmap 14) — modest.
(Refines Deliverable-B's single-pass STAR 0/111: the multi-aligner COMPASS panel at full depth rescues 2.)
FOLLOW-UP for Kevin: IGV the 2 corroborated (annotation/canonical/coverage) to accept-or-reject as real.
Merged BAM `$OAK/A549_rep1.consensus.bam` (6.9G) + json + sentinel LOCKED read-only (lock job 31412492).

## STATUS: COMPASS 111-deliverable COMPLETE & DELIVERED. SSH restored (168h ControlPersist; Kerberos tickets
expire ~daily → if `ssh sherlock` fails auth again, USER runs `! ssh sherlock` + Duo to re-create the master).

## NEXT (forward — pick up here)
1. **IGV the 2 corroborated junctions** (the live science follow-up): `chr5:140564954-140565547` (depth 12,
   strong) and `chr5:179823051-179823857` (depth 3, modest). Check annotation status, canonical GT-AG motif,
   and read coverage in the locked Oak BAM `$OAK/A549_rep1.consensus.bam` → accept/reject as real novel
   junctions. These are the 2/111 the multi-aligner panel rescued that single-pass STAR (Deliverable B) missed.
2. **(Optional) extend to A549 rep2/rep3** — this run was rep1 only. To run another rep/sample: reuse the
   pattern in `$W/cmp_sr_full_split_v2.sbatch` (re-split → array → merge → adjudicate → lock), but **ALWAYS
   pass `--oak-output-dir` on a REAL `/oak/...` path** (a fresh per-sample dir) and confirm the generated merge
   prints `Oak out: /oak/...` before launching. Pull the rep's FASTQ to `$W/COMPASS/fastq/` first.
3. **Rule-compliance hardening (small):** the rectify-generated SR final-merge still contains `rm -rf "$OUTDIR"`
   + a `rm -f …*_R1.fastq.gz` wildcard (`split_command.py` SR merge template, ~L2333). Safe today (scratch-only,
   Oak-verified) but violates the new global no-`rm -rf`/no-wildcard rule — rewrite to an explicit named-path
   removal. **Coordinate first** (shared `drs-validation-rebuild`, active DRS-arm agent) — don't clobber.
4. **Human DRS/cDNA dev backlog** (Round-2 cDNA realign incl. Phase-0 kill-gate + BLOCKER-1; cDNA
   `corrected_reads.tsv`; native-RNA004 validation of the IVT penalty tables; `--max-intron 500000` human
   default) — **owned by the concurrently-active DRS-arm agent.** Do NOT auto-build on the shared branch from
   here; coordinate via `.claude/inbox/`. Use `rectify_src_dev` for any isolated test runs.

## RE-READ THE DELIVERED RESULT (anytime)
`ssh sherlock 'OAK=/oak/stanford/groups/larsms/Users/kevinroy/compass_a549_out; cat $OAK/adjudication_111.json'`
(output is chmod-444 / dir non-writable — read-only, protected from accidental rm). All chain jobs COMPLETED
(`31165983` split, `31202917` array 500/500, `31202976` merge, `31202977` adjudicate, `31412492` lock).

## FILES / IDs
- Re-run chain (live): split `31165983` ✅ → array `31202917` (A549_rep1_sr, ~373/500) → merge `31202976`
  → adjudicate `31202977` → **lock `31412492`** (afterok; `chmod a-w` the Oak BAM/json/sentinel/inputs +
  dir, so the output is read-only & non-deletable). Merge verified Oak-safe (guard → rm -rf hits scratch only).
  FOLLOW-UP (rule compliance): the generated merge still uses `rm -rf "$OUTDIR"` + a `rm -f …*_R1.fastq.gz`
  wildcard — safe (scratch-only) but should be hardened to explicit-path removal per the new no-rm-rf rule.
- (orig) split+chain job `31165983` → array `A549_rep1_sr` → merge → adjudicate (afterok).
  W=`/scratch/users/kevinroy/compass_a549`; OAK=`/oak/stanford/groups/larsms/Users/kevinroy/compass_a549_out`
  (durable: `A549_rep1.consensus.bam`, `adjudication_111.json`, `.adjudication_111_rc`, `inputs/` FASTQ copies).
- Sbatch: `$W/cmp_sr_full_split_v2.sbatch` (re-split+chain, passes `--oak-output-dir`), `$W/adjudicate_oak.sbatch`.
- Tools: `dev/compass_shortread_adjudicate_111.py`, `dev/compass_version_compare.py`. 111 list:
  `/scratch/users/kevinroy/deliverable_b/rectify_src/dev/gmap_only_recurrent_novels_chr5.tsv`.
- Env: `rectify` (py3.9, pysam 0.23); aligner symlinks `~/.rectify/bin`; latest-aligner env `compass_latest`
  (micromamba) for the version experiment; genome `…/COMPASS/genome_references/GRCh38_gencode_v44.fasta`.
- Key commits (in HEAD `b08f1ac`+): `daea0fa` locus-index, `700709f` gsnap flag, `4ff87a2` merge-cleanup guard,
  `3d761d0`/`7a36c77` docs. Coordination: `.claude/inbox/20260625T0822*__from-compass-shortread__*`.
════════════════════════════════════════════════════════════════════════════════

# HANDOFF — RECTIFY short-read P5 RUN on Sherlock (2026-06-19)

Picks up `dev/HANDOFF_SHORTREAD.md` (P1–P4 code done, committed, locally tested; pre-flight green).
This doc tracks the **P5 cluster run** of the paired short-read COMPASS pipeline on A549 to adjudicate
the **111 GMAP-only recurrent novels** (`dev/gmap_only_recurrent_novels_chr5.tsv`).

## Key correction to the original handoff's assumption
The original handoff said run rectify **in the compass env**. That is WRONG: the `compass` conda env is
**Python 3.7.12 / pysam 0.15.3**, and rectify requires **py≥3.8** (it dies at import — `config.py:33` uses
py3.9 generic subscription). Resolution:
- **Run rectify in the `rectify` conda env** (`/home/groups/larsms/users/kevinroy/anaconda3/envs/rectify`,
  **py3.9.23 / pysam 0.23.3**) via `python -m rectify` (NOT a pip console-script — we do NOT `pip install`,
  to avoid mutating the shared env; `cd $RECTIFY_SRC` puts our branch on sys.path).
- The 4 specialized COMPASS aligners (**STAR, hisat2, magicblast, gsnap**) live ONLY in the `compass` env.
  They are exposed to the rectify env via **symlinks in `$HOME/.rectify/bin/`** (already on the generated
  script's PATH). Verified all 4 load via symlink with compass NOT on PATH (RPATH/$ORIGIN resolves). bbmap.sh,
  samtools, java come from the rectify env (fine).

## Done
- Deployed branch `drs-validation-rebuild` (local HEAD `b4d4b48`, incl. all P1–P5 commits through `2e0388c`)
  to Sherlock: `rsync --max-size=2M` (drops 2GB bundled yeast data) → `/scratch/users/kevinroy/compass_a549/rectify_src`.
- Bumped two generator constants in `rectify/core/commands/split_command.py` (local + redeployed):
  `SR_ALIGN_CORRECT_MEM_GB 32→64` (STAR loads a **29G** human index; sequential panel → peak ≈ one STAR;
  32G OOMs) and `SR_ALIGN_CORRECT_TIME '1:30:00'→'6:00:00'` (per-chunk 7-aligner index loads dominate;
  larsms is non-preempt). **NOT yet committed in git.**
- Generated + READ the smoke scripts (caught nothing wrong — the ANNOT path that looked off was just the
  COMPASS gtf symlink `.resolve()`-d to its real target `rectify_human_validation/.../gencode.v44.annotation.gtf`;
  identical file, 3,424,194 lines).

## Verified
- `rectify split --help` (deployed code) shows `-2/--read2 --read-length -n --short-read --generate-slurm
  --python-path --rectify-src --slurm-partition --slurm-account`.
- All 4 aligner symlinks load: STAR 2.7.10a, hisat2 2.2.1, magicblast 1.5.0, gsnap (runs; SSE4.2 note only).
- `for aligner in aligners:` in `multi_aligner.py:2802` is **sequential** → mem sizing = one STAR.
- Subsampled smoke input = **100k pairs** each (`A549_rep1_subsampled_{1,2}.fastq.gz`).

## CONSENSUS HANG — root-caused and FIXED (the big one; uncommitted)
Smoke attempt 2 ran the panel but then **hung ~2h in consensus** on 25k pairs. faulthandler traceback:
`select_best_alignment → score_alignment → _rescue_5prime_softclip → edit_distance`. Root cause
(`select.py:71-75`): the per-read 5' soft-clip rescue pool was built by scanning **every annotated
junction on the read's whole chromosome** (~tens of thousands) — O(reads × junctions_on_chrom), which is
fine for chr5-only long-read (deliverable_b) but hangs genome-wide at short-read scale (~50k reads). The
rescue only ever inspects junctions within `search_window_bp=300` of the read 5' boundary, so the chrom
scan was pure waste.
**Fix** (`select.py`): a memoized per-chrom **position index** (sorted start/end arrays + bisect) over
`annotated_junctions`; the candidate pool is now a locus window `[read_lo-350, read_hi+350]` query — a
**window-superset of the rescue's ≤300bp reach ⇒ byte-identical scoring**, O(log n) per read. The
`annotated_junctions` SET is kept intact for the `_n_annotated` tiebreak membership test. Helps the
long-read path identically.
**Verified in tight loop** (no SLURM burn): `repro2_consensus.py` on the 5 smoke BAMs went from >2h-hang →
**DONE 50.6s**, output BAM 64,586 records (= merge groups), **28,009 spliced N-op reads present**.

## gsnap + magicblast FAIL (panel currently 5/7) — env/tool issues, NOT code
- **gsnap**: compass env has GMAP workers (`gmap.sse42`,`gmap.nosimd`) but **NO gsnap worker binaries**
  (`gsnap.sse42`/`.nosimd` absent) — the `gsnap` dispatcher prints "does not exist" for every ISA and exits
  in 0.1s. The GSNAP half of the gmap/gsnap package was never installed. Fix = install complete gmap/gsnap
  (e.g. `conda install -n <env> -c bioconda gmap`, then re-point the `~/.rectify/bin/gsnap` symlink).
- **magicblast**: `.bam` has a valid `@HD` header but **0 records**; `samtools view failed ...magicblast.sam
  (exit 1)`. magicblast 1.5.0 (2019) emitted a SAM the rectify-env samtools rejected. Fix = inspect raw SAM /
  samtools-version compat / params in `run_magicblast`.
- 5 working = bbmap + STAR_default + STAR_noncanonical + HISAT2_default + HISAT2_noncanonical (the sensitive
  splice-aware workhorses). **USER CHOSE: fix env for full 7/7 before the run.**

### gsnap FIX (done, validated)
gmap/gsnap **2021-05-27** workers (`gsnap.sse42`,`gsnap.nosimd`,`gsnapl.*`) DO exist in compass/bin — the
dispatcher just looks for `gsnap.<isa>` **next to argv[0]** (= the `~/.rectify/bin` symlink dir), where only
`gsnap` itself was symlinked. Fix = symlinked the 4 workers into `~/.rectify/bin` too. `gsnap --version` via
that PATH now resolves `gsnap.sse42` cleanly (only a harmless `avx2 does not exist` note). **No env mutation.**
The GSNAP index was built with this same gmap 2021-05-27 → compatible.

### magicblast FIX (done, validated on subset; code change)
magicblast 1.5.0 under `-no_query_id_trim` keeps the whole multi-token FASTQ header (rectify-split injects
`RN:i:N` + the Casava comment) and spills those tokens into SAM cols 2-3, shifting the mandatory fields →
FLAG="RN:i:0" → samtools view exit 1 → 0 records. Dropping `-no_query_id_trim` makes magicblast fail
differently (rc=1). **Fix** (`multi_aligner.py::run_magicblast._ensure_plain`): rewrite every input-FASTQ
header to its bare QNAME (first whitespace token) before magicblast; RN is re-applied afterward from the
qname→RN map in `_finalize_short_read_bam` (keyed on bare qname). Verified on a 1k-pair subset: SAM columns
correct, `samtools view` rc=0, 2273 records. **Uncommitted.**

## ALIGNER-VERSION SENSITIVITY EXPERIMENT (user-requested 2026-06-20)
Goal: confidence aligner version plays a minimal role. Head-to-head on the SAME reads (100k-pair
A549_rep1 subsample), identical rectify pipeline, differing ONLY in aligner binaries/indices.
- Installed (COMPASS-pinned): STAR 2.7.10a, HISAT2 2.2.1, magicblast 1.5.0, gsnap 2021-05-27, bbmap 38.18.
  Latest: STAR 2.7.11b, HISAT2 2.2.1 (==, project dormant), magicblast 1.7.2, gmap/gsnap 2025.07.31, bbmap 39.79.
- Env `compass_latest` BUILT (via **micromamba** — classic conda solver OOMs on the login node; conda has
  no libmamba solver plugin and Sherlock curl is flaky, so micromamba binary was fetched on the M1 and
  rsync'd to `~/bin/micromamba`). Installed: STAR 2.7.11b, HISAT2 2.2.2, magicblast 1.6.0 (bioconda's
  newest is 1.7.0; solver picked 1.6.0), gsnap 2024-11-20, bbmap 39.26, samtools 1.21. Path:
  `/home/groups/larsms/users/kevinroy/anaconda3/envs/compass_latest`.
- **Invocation**: env has its own python 3.12, so run rectify via the ABSOLUTE rectify-env python with
  `PATH=compass_latest/bin:$PATH` (aligners + their wrapper-siblings — gsnap workers, hisat2-align-s,
  bbmap jar — all resolve from that one dir; no symlinks needed). Testing whether STAR 2.7.11b / gsnap
  2024 read the existing indices; rebuild only incompatible ones into `genome_references_latest/`.
- Baseline (COMPASS-pinned panel on the 100k subsample) → `$W/ver_cmp/compass/cmpver.rectified.bam`.
- Latest panel → `$W/ver_cmp/latest/latestver.rectified.bam`.
- Compare: `dev/compass_version_compare.py` (chr5 Jaccard, depth Spearman, novel-junction agreement,
  111∩each). High concordance ⇒ version minimal role. (Separate from the main full run job 30432422.)
- **Index compatibility tested**: STAR 2.7.11b READS the existing 2.7.4a index (reuse); HISAT2 2.2.2 reads
  2.2.1 index (reuse); BLAST db reused; bbmap builds on-the-fly. **gsnap 2024 CANNOT read the 2021 gmap
  index (rc=9)** → rebuilt with latest `gmap_build` (COMPASS recipe: gmap_build + `gtf_splicesites|awk
  '$4>9'|iit_store`). `genome_references_latest/` = symlinks to reusable + fresh GSNAP/.
- **Chained jobs**: `gsnap_latest_build` **30443075** (~1-2h, sentinel `$W/.gsnap_latest_build_rc`) →
  `latest_panel_cmp` **30443095** (`afterok` dep; runs latest panel on the 100k subsample then the
  comparison, writes `$W/ver_cmp/version_compare.json` + sentinel `$W/.latest_cmp_rc`). COMPASS-pinned
  baseline already done: `$W/ver_cmp/compass/cmpver.rectified.bam` (all 7 aligners, 311,614 reads).
- **RESUME (experiment)**: `cat $W/.latest_cmp_rc` → if `0`, read `$W/ver_cmp/version_compare.json`
  (jaccard_all, novel_jaccard, depth_spearman, 111_in_compass_pinned vs 111_in_latest, verdict). If a job
  failed: `sacct --name=gsnap_latest_build,latest_panel_cmp`; gsnap build logs `$W/logs/gsnap_latest_build_*`.
  Versions: STAR 2.7.10a→2.7.11b, HISAT2 2.2.1→2.2.2, magicblast 1.5.0→1.6.0, gsnap 2021-05-27→2024-11-20,
  bbmap 38.18→39.26.

## chr5 PREFILTER feasibility (user question 2026-06-20) — single aligner is NOT sufficient
Q: can a fast single-aligner pass prefilter reads to chr5 (avoid running the 7-panel on all reads, since
chr5 is only **4.8%** of mapped reads — measured)? Validated on the 100k subsample vs the all-reads truth
(`cmpver`); metric = chr5 junctions the prefilter MISSES.
- **STAR EndToEnd (COMPASS config)**: captured 2362/~7400 consensus chr5 pairs; **missed 543/1655 junctions
  (33%), 258/421 novel (61%)**. jaccard 0.67.
- **STAR permissive** (Local + `--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3
  --outFilterMultimapNmax 200 --seedSearchStartLmax 30`): better but still **missed 149/1655 (9%)**, and the
  misses are NOT low-depth noise — **novel junctions missed: 43% at ≥2 reads, 53% at ≥3, 55% at ≥10**. jaccard
  0.90, novel-Jaccard 0.74.
- **Mechanism / CONCLUSION**: the consensus chr5 read set is the UNION of 7 heterogeneous aligners and is
  **magicblast-dominated (magicblast wins 44% of consensus reads)**; magicblast (BLAST-based) maps spliced/
  divergent reads no STAR config reproduces. So **a single aligner cannot prefilter at chrom level accurately**
  for this multi-aligner consensus — it drops ~half the well-supported NOVEL junctions (the category the 111
  live in). The whole-genome ALL-READS alignment is therefore not merely conservative, it is **required** for
  faithful novel-junction recovery. A faithful prefilter would need a multi-aligner union INCLUDING magicblast
  (the slow, dominant one) → ~no net speedup. **Decision: keep the all-reads run; the prefilter shortcut is
  refuted by data.** (111∩=0 in both prefilter and all-reads here, but that is coincidental to the 111 being
  absent everywhere; for novel-junction discovery generally the single-aligner prefilter would mislead.)
- Killed the strict full-data prefilter job (30466878); no single-aligner fast-track relaunched.

## Test suite (user /goal: all tests pass) — GREEN after fixes
Full `pytest tests/` in the `rectify` env (py3.9), bundled data deployed (54M, indexes excluded):
- My 4 changes cause **zero regressions** — the 9 initial "failures" were all missing-bundled-data
  artifacts of the `--max-size=2M` deploy (yeast genome/gff); they pass once `rectify/data` is present.
- **Added regression tests**: `tests/test_consensus_locus_index.py` (window index correctness +
  far-junction scoring invariance), `tests/test_magicblast_header_strip.py` (bare-qname strip + gz +
  `check_aligner_available(None)` no-raise). Also hardened `check_aligner_available(None)→False` and
  refactored the magicblast strip into module-level `_write_bare_qname_fastq` (testable).
- **One pre-existing failure resolved, not mine**: `test_bam_parallel_state::…_deterministic` is a golden-
  hash test whose golden (recorded 2026-05-22 @55089f7) went stale when later 3'SS-rescue commits
  (bd20f9e/961c844/cf5ebb9/0c1773b) legitimately changed `process_bam_file_parallel` output. Proven not
  mine: baseline `select.py` yields the IDENTICAL observed hash; output is deterministic (stable across
  runs). Re-recorded the golden to `4f326833…` with documented rationale. (It was masked earlier because
  it's skipped when bundled data is absent.)

## ALIGNER-VERSION EXPERIMENT — RESULT (2026-06-20): version plays a NEGLIGIBLE role
Head-to-head on the same 100k-pair A549 subsample, identical rectify pipeline.
- **gsnap caveat**: gsnap 2024-11-20 removed `--ambig-splice-noclip` (rectify's shared wrapper passes it for
  gsnap 2021), so latest gsnap won't drop into the COMPASS flags unchanged → gsnap excluded from BOTH sides
  for a clean matched comparison (its version effect is UNtested, not shown-minimal). The other 4 distinct
  tools all upgraded cleanly: STAR 2.7.10a→2.7.11b, HISAT2 2.2.1→2.2.2, magicblast 1.5.0→1.6.0, bbmap 38.18→39.26.
- **CLEAN 6v6** (`$W/ver_cmp/version_compare_6v6.json`, bbmap+STAR×2+HISAT2×2+magicblast, pinned vs latest):
  - chr5 junctions 1661 vs 1662; **jaccard_all 0.9982**, **depth Spearman 0.9998**.
  - novel-junction Jaccard **0.9954 at ≥1 read**, **1.0 at ≥2/≥3/≥5/≥10 reads** (the single ≥1-read
    discordant junction is one 1-read call — detection-boundary noise, not algorithmic).
  - **111∩panel = 0 for BOTH** versions → the artifact conclusion is version-invariant.
- (As-is 7-pinned vs 6-latest, gsnap-confounded: jaccard 0.969, depth 0.971, novel 0.924 — conservative floor.)
- **Conclusion**: for STAR/HISAT2/magicblast/bbmap, latest vs COMPASS-pinned produce IDENTICAL well-supported
  junctions; aligner version does not change the 111 result.
- **IMPORTANT setup nuance**: `genome_references_latest/` fasta was a SYMLINK; rectify's
  `_compass_index_paths` does `Path(genome).resolve().parent`, and `.resolve()` followed the symlink back to
  the OLD dir — so the 6v6 "latest" panel actually used the OLD indices with the LATEST binaries (a clean
  binary-version isolation, still valid). Fixed by HARDLINKING the fasta+fai into `genome_references_latest/`
  so ref_dir resolves there (needed for gsnap's rebuilt index to be used).
- **gsnap 2024 — NOT a drop-in** (committed flag fix `_gsnap_supports_ambig_noclip`): (1) removed
  `--ambig-splice-noclip` (worked around); (2) requires a freshly-built index (rebuilt); (3) **~100x slower**
  — its new `localdb` mode allocs 12.4G and runs **13.6 queries/sec** (vs gsnap 2021 fast), so the full-100k
  rectify run hit a wall and exited nonzero at 512s. Standalone on 1k pairs it returns rc=0 / correct output,
  so the version is fine, just impractically slow as configured. Direct gsnap-2021-vs-2024 junction
  comparison on a 12k-pair subset: ATTEMPTED but **gsnap 2024 SIGSEGV**s on read
  `K00151:...:1112:3752:43128` (`Access_emergency_cleanup`), truncating output (12594 vs gsnap2021's 25828
  records; 55 vs 238 chr5 junctions) → comparison guardrail correctly refused. So gsnap 2024-11-20 is
  **broken on this A549 data** (crash bug), making its version effect untestable — which is itself the
  strongest reproducibility argument for KEEPING the COMPASS-pinned gsnap 2021-05-27.
- **FINAL ANSWER (user's question)**: aligner version plays a NEGLIGIBLE role for STAR/HISAT2/magicblast/bbmap
  (6v6: novel-Jaccard 1.0 at ≥2 reads, 111∩=0 on both). gsnap-latest is not a drop-in and crashes here, so
  pinning gsnap (as COMPASS does) is the right call; its junctions are corroborated by the other 6 aligners
  in consensus anyway. Net: the 111-artifact conclusion is robust to aligner version.

## Tool maintenance check (user-requested 2026-06-19) — both ACTIVELY maintained; ours are OLD
- **GMAP/GSNAP**: latest bioconda **2025.07.31** (Sep 2025); releases throughout 2024–2025; maintained by
  Thomas Wu / Genentech. **We run 2021-05-27** (~4 yr old). Functional (index built with same version).
- **Magic-BLAST**: GitHub commits into **Apr 2025**; latest tagged release **1.7.2 (Apr 2023)**. **We run
  1.5.0 (Aug 2019)** — ancient; its SAM-column mangling (now worked around) is plausibly fixed in 1.7.x.
- **Recommendation**: for the production/publishable run, consider upgrading both (esp. magicblast → 1.7.2)
  in a fresh env, then drop the magicblast header-strip workaround if 1.7.x handles multi-token headers.
  Not blocking — current fixes make the 2019/2021 versions produce correct 7/7 output.

## 7/7 END-TO-END VALIDATED (chunk-0, 2026-06-19 ~22:25)
Full `rectify align --short-read --read2 --aligners all` on smoke chunk 0 (25k pairs), 7m31s total:
- **All 7 aligners contribute**: By-aligner winners `{STAR_noncanonical 16822, STAR_default 1241,
  HISAT2_noncanonical 6067, bbmap 5206, gsnap 6583, magicblast 33459, HISAT2_default 5}` — gsnap+magicblast
  fixed and winning.
- **Consensus 187s for 69,383 (RN,mate) groups, no hang** (locus-index fix holds at full panel). 5' rescue
  intact (3,456 rescued, edit=0 matches in debug log).
- **Both mates survive**: records-per-RN dominated by 2 (10,223 RN); only 110 RN (0.5%) single-record →
  adversarial mate-drop bug absent. Rectified BAM 69,383 records (multimapper/secondary records inflate
  >2/RN; fine for junction extraction).
Output: `/scratch/users/kevinroy/compass_a549/chunk0_7aligner_test/chunk0test.rectified.bam`.

## Bug fixed before smoke could pass (commit-worthy, uncommitted)
First smoke (`30367486`) fast-failed: `align_command._run_one_aligner` sets `exec_path=None` for the
COMPASS panel (binary resolved inside the wrapper) but still called `check_aligner_available(None)` →
`shutil.which(None)` → TypeError. **Fix** (`align_command.py:475`): guard the check with
`exec_path is not None and not check_aligner_available(exec_path)`. This dispatch path had never executed
(handoff said so). Redeployed; COMPASS dispatch branches (583–629) otherwise present and correct.

## STATUS 2026-06-19 ~22:40 — pipeline 7/7 validated, tests green, NOT yet launched full run
Both smoke arrays (30367486, 30367768) were cancelled — the first hit the consensus hang; the second was
superseded by the direct chunk-0 7/7 validation (see below) once all bugs were fixed. **No cluster job is
currently running.** The pipeline is validated end-to-end at 7/7 on a real chunk; the FULL ~500-chunk A549
run has NOT been launched yet (awaiting go-ahead). The deployed code (`…/compass_a549/rectify_src`) carries
all fixes; `~/.rectify/bin` has all aligner symlinks incl. gsnap workers.

## NIGHT 2026-06-25 (autonomous) — conservative plan; shared branch has a CONCURRENT agent
The human-DRS/cDNA-arm agent is actively committing `drs-validation-rebuild` (last commit b08f1ac @01:13,
~minutes ago). So tonight I am NOT committing dev work to the shared branch (no Round-2/cDNA-TSV/BLOCKER-1) —
that would clobber their in-flight tree. My activity is confined to: (1) the COMPASS deliverable in
`/scratch/users/kevinroy/compass_a549` (isolated), (2) read-only test runs in `rectify_src_dev`, (3) this
handoff + the inbox note `.claude/inbox/20260625T0822*__from-compass-shortread__*`.
- **De-duped backlog (Explore survey was STALE — these are ALREADY DONE at HEAD by the other agent):**
  anchor-gate CLI (`--min-junction-anchor-bp`/`resolve_min_junction_anchor_bp`), gapmm2 pin (5c59f99),
  cluster_com (507c4ee), splice relabel (116bb28), Dorado pt:i polyA (879aa96), GMAP opt-in aligner (7b32fa0).
- **Still open (per docs, for Kevin/DRS-agent — DO NOT auto-build on shared branch):** Round-2 cDNA realign
  (zero-code; Phase-0 kill-gate FIRST; BLOCKER-1 = Cat3 exemption vs anchor gate — verify it still exists,
  `_effective_chimera_ok` not found at HEAD so may be stale/renamed), cDNA `corrected_reads.tsv` export,
  native-RNA004 validation of IVT penalty tables, `--max-intron 500000` human default, human 9-cat validation
  vetting. Recommend building Round-2 only when Kevin is awake to point it (advisor concurs — a run-and-trust
  GO/NO-GO harness is worse than none).
- **Re-run DELIVERABLE (to durable Oak): split+chain job `31165983`** (re-split → array → merge → adjudicate).
  Split RUNNING (~2h in of ~4h). Merge writes BAM to `/oak/.../compass_a549_out/A549_rep1.consensus.bam`
  (durable) then safe cleanup; adjudicate (afterok) → `$OAK/adjudication_111.json` + `$OAK/.adjudication_111_rc`.
  Source FASTQs copied to `$OAK/inputs/` as scratch-purge insurance.
  **RESUME**: `ssh sherlock 'cat /oak/stanford/groups/larsms/Users/kevinroy/compass_a549_out/.adjudication_111_rc'`
  — if `0`, read `$OAK/adjudication_111.json` (positive control / decoys / 111∩COMPASS together).
- **Test goal — MET**: full suite at HEAD b08f1ac (in isolated `rectify_src_dev`, my fixes + the concurrent
  agent's work) = **1608 passed, 0 failed, 40 skipped, 1 xfailed** (348s). No regressions from anyone; green.

## DATA-LOSS ROOT CAUSE + FIX (2026-06-25) — self-inflicted merge bug, NOT an agent
The first full run (30443246) finished all 500 chunks + merge (30443247 COMPLETED 06-21T17:48), but the
adjudication failed rc=2 (no BAM) and the whole `rectify_sr_full/` (incl. 500 chunk consensus BAMs + merged
BAM) was GONE.
- **Root cause**: the generated SR final-merge ends with `rm -rf "$OUTDIR"` to "clean up scratch" AFTER
  copying the merged BAM to `$OAK_OUT`. With **no `--oak-output-dir` given, `OAK_OUT` defaulted to
  `$OUTDIR/final`** (a SUBDIR of the working dir), so the cleanup deleted the only copy of the deliverable.
  Self-inflicted by the pipeline + my omission of a real Oak path.
- **Forensics (definitive)**: NO `rm`/`-delete` of `compass_a549`/`rectify_sr_full` in ANY of 1036 M1 session
  transcripts; the hard-drive-rescue/offload agents touch only the local 5TB drive / Drive, never Sherlock —
  **exonerated**. `compass_a549` mtime = 06-21T17:48:50 (13s after merge end) = the merge's own cleanup. The
  smoke dir survived only because its arrays were cancelled before merge ran. Not the scratch 90-day purge
  (the overdue 8TB old COMPASS_alignments_archive is still present).
- **FIX (committed)**: `split_command.py` merge cleanup now guards — only `rm -rf "$OUTDIR"` when `$OAK_OUT`
  is OUTSIDE it; else keep the merged BAM, reclaim only intermediates. Verified the generated merge script.
- **RE-LAUNCHED to durable Oak (job 31165983, 2026-06-25)**: re-split → array → merge → adjudicate, with
  `--oak-output-dir /oak/stanford/groups/larsms/Users/kevinroy/compass_a549_out`. Merged BAM lands at
  `$OAK/A549_rep1.consensus.bam` (durable); adjudication writes `$OAK/adjudication_111.json` + sentinel
  `$OAK/.adjudication_111_rc`. **RESUME**: `cat $OAK/.adjudication_111_rc`; if 0, read `$OAK/adjudication_111.json`.

## (FIRST ATTEMPT, data lost — see above) FULL RUN LAUNCHED 2026-06-20 — split+chain job `30432422`
(First attempt 30431657 fast-failed on the `set -u` + conda-activate trap — conda's java_home.sh
references unbound JAVA_HOME; fixed sbatch to `set -o pipefail` only. Resubmitted as 30432422.)
Submitted `/scratch/users/kevinroy/compass_a549/cmp_sr_full_split.sbatch`. It:
1. `rectify split` the full 42M-pair A549_rep1 (R1/R2) into **500 chunks** → `$W/rectify_sr_full/`
   (~4h at observed throughput; writes sentinel `$W/.sr_full_split_rc` with the rc).
2. On rc==0, chains `bash $W/rectify_sr_full/submit_pipeline.sh` → submits the **500-task array**
   (`A549_rep1_sr`, 64G/6h/task, idempotent `.consensus.bam` skip) + `afterok` final merge.

### RESUME — concrete branch logic for the full run
SSH `sherlock` open; never tear down ControlMaster; retry transient sshd serially. Env: `rectify` conda
env + `export PATH=$PATH:$HOME/.rectify/bin`. `$W=/scratch/users/kevinroy/compass_a549`.
```
ssh sherlock "squeue -u kevinroy -o '%.14i %.16j %.8T %.10M %R'; echo ---; \
  cat $W/.sr_full_split_rc 2>/dev/null; echo '--- chunks done ---'; \
  ls $W/rectify_sr_full/chunk_outputs/*.consensus.bam 2>/dev/null | wc -l; \
  ls $W/rectify_sr_full/final/ 2>/dev/null"
```
- **split job 30432422 RUNNING** → wait (~4h).
- **`.sr_full_split_rc` absent & job gone** → split died; check `$W/logs/cmp_sr_split_*.{out,err}`.
- **`.sr_full_split_rc` == 0** → array was chain-submitted; find it: `squeue`/`sacct --name=A549_rep1_sr`.
  - array tasks idempotent (`.consensus.bam` skip + atomic copy) → safe to requeue / re-`sbatch
    rectify_sr_full/run_array_short_read.sh` if some failed.
  - **all 500 `.consensus.bam` present + merge COMPLETED** → final merged BAM in `rectify_sr_full/final/`.
    Proceed to adjudication (below).
- **`.sr_full_split_rc` != 0** → split failed; do NOT expect an array. Read split err.

### ADJUDICATION (P5 step 4) — CLUSTER-DURABLE, chained on the merge
Tool: `dev/compass_shortread_adjudicate_111.py`. Submitted as SLURM job **30455867**
(`afterok:30443247` the merge), so it runs on the cluster independent of the M1 watcher →
writes `$W/rectify_sr_full/adjudication_111.json` + sentinel `$W/.adjudication_111_rc`.
Reports positive control (annotated chr5 junctions HIGH), negative (~0), and `111 ∩ COMPASS`.
Near-zero intersection ⇒ the 111 are artifacts — ONLY valid if the positive control passed.
**RESUME**: `cat $W/.adjudication_111_rc`; if `0`, read `$W/rectify_sr_full/adjudication_111.json`
(`POSITIVE_control_...`, `NEGATIVE_control_decoys_in_compass`, `INTERSECTION_111_in_compass`, `verdict`).
Chain: array `30443246` (500 tasks, ~5-wide × ~15-18min ≈ ~1 day) → merge `30443247` → adjudicate `30455867`.

## (superseded) Resume — concrete branch logic
SSH ControlMaster `sherlock` is open; never tear it down; retry transient sshd errors serially.
Check the smoke:
```
ssh sherlock "sacct -j 30367486 -X -o JobID,State,Elapsed,MaxRSS,Start,End; \
  ls -la /scratch/users/kevinroy/compass_a549/rectify_sr_smoke/chunk_outputs/ 2>&1"
```
- **If all 4 tasks COMPLETED** → run the acceptance checks (per-chunk log shows all 7 aligners succeeding;
  merged BAM has **two records per RN** — the adversarial-review mate-drop fix; N-op junctions present):
  ```
  L=/scratch/users/kevinroy/compass_a549/rectify_sr_smoke/logs
  grep -iE "aligner|STAR|hisat|magicblast|gsnap|bbmap|reads|error|0 reads" $L/30367486_0.out | head -60
  # both-mates check on a chunk consensus BAM:
  source .../conda.sh && conda activate rectify
  samtools view .../chunk_outputs/A549_rep1_subsampled_chunk_000_of_004.consensus.bam | \
    grep -oE "RN:i:[0-9]+" | sort | uniq -c | awk '{print $1}' | sort | uniq -c
  ```
  If both-mates and 7-aligners pass → proceed to FULL run.
- **If any task FAILED** → read `$L/30367486_<task>.err`; most likely causes: an aligner not found on PATH
  (symlink/PATH), STAR `versionGenome … INCOMPATIBLE` (regenerate index, ~30min), or OOM (bump mem again).

### FULL run (after smoke passes)
The split of the FULL 42M-pair fastq (7.5GB gz) is heavy — run it AS A JOB, not on the login node.
Generate with the same flags but input `$W/COMPASS/fastq/A549_rep1_R{1,2}.fastq.gz` and **no `-n`**
(auto-size ~500 chunks via `--target-reads-per-chunk`). Then `bash submit_pipeline.sh`. Drop a sentinel,
refresh THIS handoff.

### Deliverable (P5 step 4) — TOOLING NOT YET WRITTEN
Junction extraction from the merged consensus BAM + the 111-adjudication (positive control: annotated chr5
junctions HIGH; negative ≈0; `111 ∩ COMPASS`) has **no code yet**. WRITE IT DURING the full run.
Report the three numbers together — near-zero `111 ∩ COMPASS` only means "artifact" if the positive control passed.

## Files
- Deployed code: `/scratch/users/kevinroy/compass_a549/rectify_src/` (rsync target; re-sync after local edits).
- rectify env python: `/home/groups/larsms/users/kevinroy/anaconda3/envs/rectify/bin/python`.
- Aligner symlinks: `$HOME/.rectify/bin/{STAR,hisat2,magicblast,gsnap}` → compass env.
- Smoke: `/scratch/users/kevinroy/compass_a549/rectify_sr_smoke/`.
- Local generator edit (uncommitted): `rectify/core/commands/split_command.py` lines 108–110.
- Source of truth for aligner params: Sherlock `…/compass_a549/COMPASS/process_reads_and_align.sh`.
