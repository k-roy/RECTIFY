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

## ⚠️ CORRECTION (advisor, 2026-06-30) — the above motif read is SNAPPING-CONFOUNDED; NO valid verdict yet
The "reads splice at canonical GT..AG" observation was measured on minimap2's N-op PLACEMENTS — but
minimap2 SNAPPING a non-canonical junction onto the nearest canonical GT-AG is the EXACT confound
Deliverable B exists to see past (OVERVIEW; smoke (G) non-canonical recall 0.000; THIS session's
c3_junction_headroom: mm2 snaps 46.7%). So "GT..AG in the mm2 BAM" is exactly what a GENUINELY
non-canonical junction looks like AFTER minimap2 snaps it (call it H3, MORE likely than H1). The
motif-at-mm2-placement evidence therefore CANNOT address whether the junction is truly non-canonical,
and H1 (frame artifact → junctions are canonical) is UNSUPPORTED. Retract any implication that the
COMPASS 3-junction non-canonical claim might be a frame artifact.
THE ACTUAL TEST NEVER RAN: n_pooled=0 → the C1 realign (which re-scores each read's RAW SEQUENCE
canon-vs-nonc, seeing past where minimap2 placed the N-op) scored ZERO reads. This is a NOT-YET-MEASURED
result, not a nuanced deferral. There is NO verdict to adjudicate until the realign runs in-frame.
Two facts that break the unified "frame-shift" narrative: SLC35A4 candidate span=593bp vs reads' span=1014bp
= a DIFFERENT intron (not a shifted one); SQSTM1 = deSALT FAILED + mm2 under-splices the 7805bp intron =
cannot be concluded at all. COMPASS's core evidence was accurate SHORT reads (untouched here).
DECISIVE PATH (the only thing that answers it): align the reads to the EXACT reference COMPASS used
($W/genome_references_latest/GRCh38_gencode_v44.fasta), NOT align/chr5.fa → the frame matches by
construction → run c1_realign_3junctions.py → read the per-read canon-vs-nonc LIKELIHOOD. Plus fix the
deSALT/splice-aware arm for SQSTM1.

## FRAME CHECK RESOLVED (2026-06-30) — reference IDENTICAL; MY re-alignment is the +176 OUTLIER, not COMPASS
Proved align/chr5.fa == the exact COMPASS reference (genome_references_latest/GRCh38_gencode_v44.fasta):
chr5 length identical (181,538,259) AND byte-identical sequence at both the candidate (177592498) and the
reads-actual (177592673) loci. So the +176 is NOT a reference/build shift.
DECISIVE ground truth from COMPASS's OWN in-frame data (lr_probe_4loci_out.txt, same reference):
  TMED9 974bp intron dominant = **177592499-177593473** (604 reads; minor neighbours 177592497-177593471:17,
  177592498-177593473:4 — the 1-4bp non-canonical/canonical cluster the COMPASS note describes), corroborated
  by accurate short reads (323). rederive_111.json gmap_coord = chr5:177592500-177593474 (≈ COMPASS dominant,
  the agent's candidate was RIGHT).
BUT my minimap2 re-align of sumner a549_chr5_trimmed.fastq put 1200 reads at 177592675-177593649 (+176, same
974bp length) → the OUTLIER. On an IDENTICAL reference the same reads cannot place 176bp apart, so my
re-aligned read set is NOT the reads COMPASS used (the sumner chr5-trimmed fastq ≠ COMPASS's read set / or a
pipeline shift), and 176bp is FAR too large to be minimap2 canonical-snapping (a few-bp effect). ⇒ MY minimap2
placement is unreliable; COMPASS's in-frame 177592499 placement is the ground truth, and the C1 realign must be
run on the reads that ACTUALLY splice there (COMPASS's long-read set), NOT my re-alignment. The n_pooled=0 was
caused by my reads splicing +176 off the (correct) candidate.
STATUS: still NO valid C1-realign verdict; but the reference-artifact / snapping hypotheses are BOTH ruled out
— the blocker is READ-SET PROVENANCE (get COMPASS's actual TMED9/SLC35A4/SQSTM1 supporting reads, e.g. from the
lr_probe source BAMs / the read IDs in lr_probe_4loci_out.txt, and realign-score THOSE). Advisor consult next.

## STOP + HONEST STATUS (2026-06-30, advisor-scoped) — "wrong read set" REFUTED; +176 is a same-reads aligner artifact
Decisive same-read test: QNAME 945df2f3-3e2e-4824-8709-89ea41bd98e3 exists in BOTH my re-align (mm2.bam) AND
COMPASS's source BAM (a549_chr5_trimmed.minimap2.bam, EXISTS at rectify_human_validation/sgnex_a549/alignments/)
→ SAME READS. So my prior "wrong read set" commit is REFUTED (same sample+reference ⇒ same splice coord
regardless of read set, as the advisor noted). The +176 is a SAME-READS ALIGNER/PIPELINE placement difference
(my single minimap2 run vs COMPASS's multi-aligner pipeline / a possible sumner-vs-gencode chr5 local indel —
the COMPASS-source read even shows a "176N" op near TMED9). It is NOT canonicity, NOT read set, NOT reference
build. Three of my earlier durable claims (canonical-frame-artifact, reference-drift, wrong-read-set) were each
overturned by the next datum — logged here so the record is honest.
RULED OUT: reference build (identical), read set (same QNAMEs), minimap2 canonical-snap (176bp ≫ few-bp snap).
REMAINING (an aligner/pipeline placement discrepancy) is a DISTRACTION from the deliverable per advisor.
NO valid C1-realign verdict was produced. IMPORTANT: this does NOT touch the COMPASS finding — COMPASS already
CROSS-PLATFORM-corroborated non-canonical (accurate SHORT reads + 4 long-read aligners agree); the C1 realign is
only CONFIRMATORY and its marginal value is low.
EXACT DECISIVE RESUME (single valid path; run ONCE, else leave open):
  Use COMPASS's OWN long-read BAM `a549_chr5_trimmed.minimap2.bam` (+ .deSALT/.uLTRA/.GMAP, same dir) which
  place TMED9 at COMPASS's coord. (1) confirm which reference it used (prov.json → sumner GRCh38_chr5.fa) and the
  TMED9/SLC35A4/SQSTM1 junction coords IN THAT frame (from lr_probe_4loci_out.txt: TMED9 dominant 177592499-177593473).
  (2) Region-select reads overlapping each locus; (3) score each read's RAW SEQUENCE canon (the GC-AG 1-4bp away)
  vs nonc (COMPASS dominant) via _hp_edit_distance + the C1 LAW table, report WITH and WITHOUT sequence-ambiguous
  (tie) reads. Honest prior = HOLD (short reads support non-canonical); a "canonical" realign result = a flag to
  re-examine the pipeline, NOT a refutation of COMPASS. Do NOT re-drill the +176.
