# Deliverable B — real-data corroboration: C1 motif-aware realign of 3 A549 SUPPORTED_NONCANONICAL junctions

Agent: Deliverable-B (solo re-run). Worktree branch `worktree-agent-adf2053e3caf83776`
(copy of `worktree-agent-a25a2c1e784ad37dc`). Status: IN PROGRESS. REPORT-BACK, no commit.

## Question
COMPASS adjudication (rederive_111.json, class=SUPPORTED_NONCANONICAL) found 3 real human
A549 chr5 splice junctions whose **dominant supported placement is NON-canonical**, but a
**canonical GT-AG/GC-AG motif sits 1-4 bp away** (= where GMAP, a minority of reads, places it).
Per junction: does a motif-aware C1 empirical-del-law REALIGN of the supporting long reads
**SNAP** the junction to the nearby canonical motif (→ canonical, aligner-misplaced) or
**HOLD** the non-canonical placement (→ genuinely non-canonical)?

## The 3 junctions (from $W/rederive_111.json, $W/lr_probe_4loci_out.txt)
W = /scratch/users/kevinroy/compass_a549 (COMPASS agent's locked data, READ-ONLY)

| gene | canonical placement (gmap) | motif | noncanon dominant placement | motif | offset | shortrd | longrd-aligner depths |
| --- | --- | --- | --- | --- | --- | --- | --- |
| SQSTM1 | chr5:179824400-179832205 | GT..AG | 179824404-179832209 (depth 2419) | GT..GA | +4/+4 | 2959 | deSALT137 uLTRA107 GMAP149 mm2/mpb 0 |
| TMED9 | chr5:177592500-177593474 | GC..AG | 177592499-177593473 (depth 287) | CG..CA | -1/-1 | 323 | mm2 628 deSALT626 uLTRA485 GMAP491 |
| SLC35A4 | chr5:140564954-140565547 | GT..AG | 140564957-140565550 (depth 156) | AG..CA | +3/+3 | 168 | mm2 102 deSALT103 uLTRA81 |

Note SQSTM1: the donor is GT in BOTH placements (GT..AG vs GT..GA); only the **acceptor** differs
(canonical AG vs noncanonical GA, +4 bp). gmap_coord supported by only 12 GMAP reads; the bulk
of reads/aligners (2419) place +4 at the non-canonical GT..GA. The other 2 are full +/-N shifts.

## PRIOR / key caveat
- **Short reads (accurate, ~0.1% err) ALSO support the non-canonical placement** at all 3 loci
  (2959/323/168 reads). Accurate short reads do NOT have the homopolymer/deletion errors that
  would cause a snap → strong independent evidence these are **genuinely non-canonical**. So the
  honest prior is HOLD. A C1 realign that HOLDS *corroborates* the short reads (C1 agrees with
  truth); a realign that SNAPS would expose a C1 failure mode (yeast del-law over-pulling to
  canonical when truth is non-canonical).
- **YEAST penalty table on HUMAN reads** — `penalty_scores.tsv` is S. cerevisiae R10.4.1. The
  del-law SHAPE (del cheaper in/near homopolymers, monotone in HP length) is the hypothesis under
  test; absolute calibration is yeast. FLAGGED, not hidden. In-silico is SUGGESTIVE; RT-PCR/Sanger
  is the wet-lab gold standard.

## Method (C1 realign = the shipped del-law scorer)
The task's quoted signature `align_exon_block_global(..., penalty_table=LAW, lam=1.0)` does NOT
match the real API (that fn uses fixed match/mismatch/gap scoring, no penalty_table arg). The
ACTUAL del-law mechanism is `rectify.core.splice.hp_penalty._hp_edit_distance(s1, s2,
penalty_table=LAW)` (the HP-context-aware DP that `junction_refiner._score_junction` uses). I use
that.

Per junction, per supporting read (harvested from $W/morechrom_trimmed.GMAP.bam — read SEQ is
aligner-independent platform truth; GMAP BAM is just the read source + approx junction locator):
1. Extract the read's contiguous query string crossing the junction: bases aligned within W bp
   upstream (exon1 tail) ++ bases within W bp downstream (exon2 head).
2. For each hypothesis H ∈ {canonical, noncanon}, build the spliced reference
   `genome[donorH-W : donorH] ++ genome[acceptorH : acceptorH+W]`.
3. Score = `_hp_edit_distance(read_query, spliced_refH, penalty_table=LAW)`. Lower cost wins.
4. SNAP = canonical wins; HOLD = noncanon wins; TIE = equal (report margin distribution).
Aggregate snap/hold/tie counts per junction. Report WITH and WITHOUT the built-in
`_CANONICAL_HP_PRIOR=0.5` canonical discount (the honest raw-evidence vs shipped-policy split).
Sanity-check vs minimap2/deSALT/uLTRA/GMAP placements (lr_probe_4loci_out.txt).

## Cluster paths (mine, isolated)
- Outputs: /scratch/users/kevinroy/c1_realdata_dB/   (mine)
- Code:    /home/groups/larsms/users/kevinroy/c1_realdata_dB_code/  (mine, NOT shared)
- Genome:  $W/genome_references_latest/GRCh38_gencode_v44.fasta

## DATA-PROVENANCE BLOCKER discovered (2026-06-30) — reference mismatch
The note's read pointers don't resolve cleanly:
- `$W/morechrom_trimmed.GMAP.bam` contains only chr1/11/17/19 (the OTHER 111-candidate
  clusters), NOT chr5. The Jun-29 morechrom run REUSED that filename; the chr5 morechrom BAMs
  lr_probe used (Jun 26) were overwritten (scratch, gone).
- The per-aligner chr5 long-read BAMs that DO exist
  (`/scratch/users/kevinroy/rectify_human_validation/sgnex_a549/alignments/a549_chr5_trimmed.*.bam`,
  the LR source `span_test.py` points to) were aligned (prov.json) to
  `sumner_lab/references/GRCh38_chr5.fa` — a DIFFERENT chr5 build than rederive_111/lr_probe used
  (`genome_references_latest/GRCh38_gencode_v44.fasta`).
- The two chr5 builds DRIFT by local indels: the dominant long-read junction sits at a CONSTANT
  per-locus offset from my rederive coords — **+826 bp at SQSTM1** (LR junction 179825226-179833031,
  also GT..AG, also 7805 bp), **+176 bp at TMED9** (LR dominant 177592675-177593649 vs lr_probe
  177592499-177593473). Not a global constant → cumulative assembly-patch indel drift.
- CONSEQUENCE: my candidate canon/nonc coords (gencode_v44 space) do NOT map onto these
  sumner-aligned BAMs; harvesting at my coords yields 0 supporting reads. The reads EXIST (same
  A549 chr5 molecules, fastq at sumner_lab/.../trim/a549_chr5_trimmed.fastq.gz) — they're just in
  a different coordinate frame.

### Fix: re-align the chr5 reads to the SAME reference rederive used (gencode_v44), then realign.
Read SEQ is reference-independent; aligning the same molecules to gencode_v44 puts their junctions
at my candidate coords. Plan: minimap2 -ax splice (fast) + GMAP/deSALT (for SQSTM1, which minimap2
under-splices) of a549_chr5_trimmed.fastq.gz -> gencode_v44, into MY scratch
/scratch/users/kevinroy/c1_realdata_dB/align/, then run c1_realign_3junctions.py on those BAMs.

## RESULTS
(pending — awaiting re-alignment job 32180508)

### Resume (if this session drops)
IN FLIGHT: sbatch **32180508** (`c1_realign_align_chr5.sbatch`) re-aligns
a549_chr5_trimmed.fastq.gz -> gencode_v44 chr5 (minimap2 + deSALT, de-novo), outputs in
`/scratch/users/kevinroy/c1_realdata_dB/align/`. Sentinel `.align_rc` = `0` on success.
```
ssh sherlock 'sacct -j 32180508 -X -o State%12; cat /scratch/users/kevinroy/c1_realdata_dB/align/.align_rc 2>/dev/null'
```
If rc==0, run the realign on the NEW BAMs then it writes c1_realign_W15.json:
```
ssh sherlock 'source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh; conda activate rectify
 A=/scratch/users/kevinroy/c1_realdata_dB/align
 cd /home/groups/larsms/users/kevinroy/c1_realdata_dB_code
 python3 c1_realign_3junctions.py \
   --genome /scratch/users/kevinroy/compass_a549/genome_references_latest/GRCh38_gencode_v44.fasta \
   --bams $A/mm2.bam,$A/desalt.bam \
   --penalty-tsv /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv \
   --window 15 --out /scratch/users/kevinroy/c1_realdata_dB/c1_realign_W15.json'
```
FIRST sanity-check the new BAMs actually splice at the candidate coords (scan N-ops near each
junction) before trusting the realign n.

## Status / resume
(pending — see HANDOFF)

## RESULTS (2026-06-30, re-align job 32180508; minimap2 arm only — deSALT arm FAILED)
Job 32180508 FAILED at the deSALT step (`FAIL desalt_aln`), but minimap2 completed (mm2.bam, 440MB,
all 3 loci covered: SQSTM1 6798 / TMED9 1254 / SLC35A4 269 reads; ~6797/1248/266 SPLICED).
Ran the recovered scorer (c1_realign_3junctions.py, TOL=8) → **n_pooled=0 at ALL 3 candidate coords.**

ROOT CAUSE (decisive, NOT an under-splicing artifact): the reads splice at STRONG, CONSISTENT coords
that are OFFSET from the rederive_111 candidate coords — the frame mismatch PERSISTS even after re-aligning
to align/chr5.fa:
| junction | rederive candidate (canon/nonc) | reads' ACTUAL dominant junction (n) | offset | motif AT read coord |
| --- | --- | --- | --- | --- |
| TMED9 | 177592500 / 177593474 | **177592675-177593649 (1200 reads)** | +175/+175 | **GT..AG (CANONICAL, offset 0)** |
| SLC35A4 | 140564954 / 140565547 | **140564858-140565872 (159 reads)** | −96/+325 | **GT..AG (CANONICAL, offset 0)** |
| SQSTM1 | 179824400 / 179832205 | (none within 400bp; needs the failed deSALT arm) | large | — |

The TMED9 +175 offset is essentially IDENTICAL to the +176 the first attempt measured on the SUMNER build —
i.e. TWO independent references (sumner GRCh38_chr5.fa AND align/chr5.fa) AGREE with each other and both sit
+175 from the rederive_111 candidate. So the rederive_111 coordinate frame is the OUTLIER; align/chr5.fa is
NOT the gencode_v44 frame the candidates are in (the sbatch's ref fetch did not produce the rederive frame).

⚠️ CONSEQUENTIAL, UNRESOLVED INTERPRETATION (do NOT commit a verdict yet — it would contradict the COMPASS
"3 genuinely non-canonical junctions" claim on possibly-mis-framed evidence):
  (H1) The COMPASS non-canonical call was a COORDINATE-FRAME ARTIFACT — the reads actually splice at a clean
       canonical GT..AG, and the "non-canonical dominant" arose from scoring at mis-framed coords. → would
       REFUTE the 3 junctions as novel-non-canonical.
  (H2) align/chr5.fa is a DIFFERENT frame/locus, so 177592675-177593649 is a DIFFERENT (canonical) intron than
       the one COMPASS flagged — I am comparing the wrong junction and cannot conclude.
DECISIVE next step (frame reconciliation, deferred): obtain the EXACT reference rederive_111/lr_probe used
(the gencode_v44 build path in $W/genome_references_latest/), lift the candidate coords through the SAME
build, OR re-derive canon/nonc candidates IN align/chr5.fa's frame from the read-supported junction + its
±4bp motif neighbourhood, THEN realign canon-vs-nonc there. Also fix + rerun the deSALT arm for SQSTM1.
Result file: /scratch/users/kevinroy/c1_realdata_dB/c1_realign_mm2only.json (all n_pooled=0, documented).
