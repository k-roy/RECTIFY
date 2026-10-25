# cDNA `=` bug: cdna-analyze skipped the poly(A) walkback for five months, and #6 broke cDNA tail lengths

Session 35d7b9b5, 2026-10-10 (Kevin, 2026-10-09 23:4x: "That cDNA bug seems like it could have a pervasive impact. Please
dive deeper and sort this out."). Evidence on H2 under `/u/scratch/k/kevinroy/prod_20260928/route_smoke/pairs/` and
`.../cdnafix2/`; scripts in this folder.

## 1. Mechanism (verified per read)

- `rectify align` writes every aligner arm through `samtools calmd -e` (`multi_aligner._apply_calmd_eq`, present by
  f2b8020, 2026-05-06): every base that matches the reference is stored as `=`. Until #6 (`cb0737f`, merged as `83af3d1`,
  2026-10-08) the multialigned BAM kept the winning arm's SEQ, so it was `=`-encoded too (78,889 of 78,898 cDNA records on the
  smoke). Arm BAMs are still `=`-encoded on master.
- The shared walkback (`correct/walkback.py::walkback_3prime_with_qpos`) anchors the 3′ end at the first read base that equals
  the reference and is not the stop base. A stored `=` never equals a genome letter, so the scan never anchors and the caller
  falls back to the raw alignment end (tail count 0).
- `rectify correct` decodes `=` at intake (`bam_processor.correct_read_3prime`, `a1728eb`, 2026-05-18), so DRS 3′ ends are fine
  (#6 alone moves 11 of 42,627 DRS corrected ends: noise). `correct-cdna` stage 1 decodes since `47e3b39` (2026-08-24), and its
  input is normally literal anyway. **`cdna-analyze` (wired `7e408b1`, 2026-05-15) never decoded**: its input is by design the
  `=`-encoded multialigned BAM.
- Chanfreau planning/578 (2026-08-24) found exactly this root cause, fixed only stage 1, and deferred the shared-walkback fix.
  On 2026-08-21 `dfbda14` read "100 % `tail_len == 0`" as "pretrimmed records have nothing to walk" and fell back to the
  carried `XA`, which hid the anchor failure.
- Proof (`cdna_eq_verify.py`, H2 15119638): cdna-analyze's own `_read_info_from_bam_record` on each c287748 cDNA record, as
  stored vs decoded (same alignment): the 3′ anchor moves UPSTREAM for **20,160 of 78,898 molecules (25.6 %)**: 1 nt 7,928;
  2 nt 7,581; 3 nt 1,342; 4 nt 674; 5–10 nt 2,458; 11–50 nt 177; none downstream. Example `r000_3` (minus strand): genome
  `…GTCCAAAA`, stored SEQ all `=`, decoded the walkback steps back over `AAAA` to the C (3509 → 3513). `pos5_corrected` moves
  for 46 (0.06 %: the bridge GGG is almost always soft-clipped, so literal). The decoded anchor equals 83af3d1's (literal SEQ)
  for 78,852 of 78,884 same-placement records.

## 2. Second bug, live on master: #6 broke cDNA tail lengths

- `dfbda14`'s rule is "walkback tail count if > 0, else carried `XA`". With `=` SEQ the count was always 0, so the correct
  stage-1 `XA` won. With literal SEQ (#6) the walk counts the few A's the stage-1 pretrim left over genomic A's and replaces the
  carried tail: on the smoke `polya_length` falls for **20,080 of 77,138 molecules, median 33 → 5** (all molecules: median
  22 → 15). So master at `83af3d1`/`f15b523` writes wrong cDNA tail lengths. Kevin's #6 smoke was 3,000 DRS reads only.
- Stage 1 records the trim: `XK:i` = bases removed from the mRNA 3′ end (`cdna/consensus.py:273`), present on every record.

## 3. What the bug changes downstream (c287748 vs 83af3d1, same stage-1 molecules; `cdna_downstream_diff.py`, H2 15119678)

| cdna-analyze output | effect of the `=` bug |
| --- | --- |
| `corrected_reads.tsv` `corrected_3prime` | 19,655 of 77,138 move upstream: 1–2 nt 15,094; 3–4 nt 1,972; 5–10 nt 2,390; > 10 nt 199 |
| `clusters.tsv` `anchor` | 17,862 of 68,534 move |
| `polya_length` / `tail_len` / tagged `XA` | published values (c287748) are the carried stage-1 `XA`: correct by accident |
| gene assignment, gene counts | 8 molecules; **0 of 131 genes change count** |
| `isoforms.tsv` (3′ groups, `pos3_modal`) | unchanged: isoforms use the RAW alignment end by design, never the anchor |
| `t1t2_pairs.tsv` | ~2 % of pairs differ (102 / 103 of ~4,936) |
| `pos5_corrected` | 0.06 % |

Not affected: DRS `correct` outputs, gene-level quantification (DE, `cdna-length-correct`), isoform-level counts. Affected:
anything that uses cDNA per-molecule 3′ positions at nucleotide resolution (CPA site maps and pileups from cDNA, DRS-vs-cDNA
3′-end comparisons, cleavage-position windows, intronic PAS positions taken from cDNA ends), and T1/T2 pairing (~2 %).

## 4. Data blast radius (H2 inventory, `cdna_inventory_job.sh`, H2 15119741)

276 `consensus_tagged.bam` outputs under `/u/project/guillom/{kevinroy,shared}` (2026-06-25 … 2026-10-06), **all `=`-encoded**,
so every one carries unwalked cDNA 3′ ends: 1082 (60), processed (36), 933_overnight (32), 1080 (24), 936_rpl20b (19),
1031_d11 (18), 1088 (18), 587_ms2_cdna_consensus (18), 736_ms2_cdna_postfix_consensus (18), 1076 (11), 1090o (6),
268_aligner_benchmark (6), others ≤ 2. Chanfreau findings docs that consume cdna-analyze outputs include the 1090 manuscript
drafts (IPA; 2026-10-09), 1075/1076 (IPA QC, polysome IPA), 907c/d/e (intronic and internal PAS), 1021 (DRS vs cDNA
quantification), 1010/1010a (MS2), 936/925/942b (reporters), 1095b/1077 (length correction: gene-level, unaffected).

## 5. Other `=` consumers checked (code sweep, master `f15b523`)

- Safe: DRS `correct` (intake decode), stage 1 (decode, literal input), the overhang resolver (decodes), `rectify cma`
  (decodes), the landing's consensus scoring (decodes; at c287748 it mis-scored `=` records, which is P1's 973 pick changes),
  rbrowse (treats `=` as a match; slices from raw Dorado BAMs).
- NET-seq endpoint calls depended on SEQ encoding before the landing (its CHANGELOG); not quantified here.
- Fragile but loud: the rectify review renderer (`visualize/read_junction.py`) counts `=` as mismatches (a `=` bundle would show
  every base as a mismatch; the cards Kevin judged were literal).
- `extract --include-sequence` writes `=` into the sequence column when given an arm or pre-#6 multialigned BAM.
- `--filter-spikein` reads SEQ before `correct`'s intake decode; its known-gene path keys on the synthetic UTR (soft-clipped,
  literal), so it likely still works (unverified; opt-in, no current chain uses it).
- External: 135 Chanfreau planning scripts read SEQ, 40 of them handle `=`; any that read an arm, a pre-#6 multialigned BAM or a
  `consensus_tagged.bam` without decoding sees `=`. Not audited one by one.

## 5b. Third problem, introduced by the landing: master's `correct` aborts on old BAMs

The landing made the intake decode strict (`alignment_tags.decode_original_sequence` refuses a `=` inside a clip or
insertion) and `bam_processor.correct_read_3prime` calls it unguarded. Verified 2026-10-10 01:47 (H2 15119859 task 5):
master `f15b523`'s `rectify correct` on the c287748 DRS route-smoke BAM aborts with "Cannot decode original SEQ '=' for
read f2ee18b7…: equals in insertion or soft clip has no reference placement" → "Processing failed". So re-running `correct`
at master on any pre-#6 multialigned BAM fails outright (the smoke BAMs carry 1 such DRS record and 19 cDNA records). A
full `align` + `correct` at master is unaffected (its multialigned SEQ is literal).

## 6. Fix (branch `fix/cdna-analyze-eq-decode`, `.worktrees/fix-cdna-eq`; NOT pushed)

- `a72721d`: cdna-analyze decodes `=` at intake (as `correct` does); `walkback_3prime`, `walkback_3prime_with_qpos` and
  `walkback_3prime_guarded` resolve a stored `=` to the reference base, so no caller can skip the walkback silently again.
- `dd9f681`: a record with `XK > 0` and a carried `XA > 0` reports the carried `XA` (fixes #6's tail regression); records
  without `XK` keep `dfbda14`'s rule, as its test pins.
- `31db07f`: cdna-analyze no longer aborts on a record whose `=` cannot be decoded.
- `fa7f315`: `decode_original_sequence(..., unresolved='N')` decodes the aligned cells and writes N for a `=` inside a clip or
  insertion (an unknown base, never evidence); `correct`'s intake and cdna-analyze use it only after the strict decode refuses,
  with a WARNING; the strict default is unchanged for every other caller.
- `tests/test_cdna_analyze_eq_seq.py`: literal vs `=`-encoded molecules on both strands, the tagged record decoded, the shared
  core, the pretrimmed tail, `=` in a clip (cdna-analyze, the decoder's two modes, `correct`'s intake).
- Combined push candidate with the CI fix: `land/fixes-20261010` = `6e0b9ed` (`.worktrees/land-fixes-20261010`).

Verification (H2 `route_smoke/cdnafix{2,3,4}/`): the new tests all FAIL at `f15b523` (control) and pass at the fix (73 + the
decoder/writer suites); cdna-analyze at the fix on the OLD `=`-encoded c287748 cDNA BAM reproduces 83af3d1's anchors
(`corrected_3prime` differs for 95 of 77,138) and keeps the published tails (`polya_length` differs for 226: untrimmed
records now measured in SEQ), 19 undecodable records handled; DRS `correct` at the fix on the literal f15b523 BAM is
identical to f15b523 in every column (42,627 rows).

## 7. Remediation options (Kevin's call)

1. Land the fix on master (with the suite and the route smoke), so master writes walked-back cDNA ends AND correct tails.
2. Re-derive cDNA outputs without re-aligning: re-run only `cdna-analyze` (≈ 1.5 min on the 79k-molecule smoke) on the
   existing multialigned BAMs. On a c287748 + fix pin this keeps every consensus pick (D1's "stay on c287748 until ISSUE-087"
   holds) and changes only the anchors; tails stay the published carried `XA`.
3. Re-examine only the analyses that use cDNA per-molecule 3′ positions at nucleotide resolution (section 3); gene-level and
   isoform-level results stand.

## 8. Final verification (2026-10-10 02:41 PDT) — all green; nothing pushed

| check | result |
| --- | --- |
| control: the 7 new tests at unfixed `f15b523` | all 7 FAIL (each catches a real defect) |
| focused suites at `fa7f315` (cDNA, walkback, decoder, writer, 2F clip origin, both validation-read sets) | 214 passed, 0 failed |
| cdna-analyze at the fix on the OLD `=` c287748 cDNA BAM | anchors = 83af3d1's (95 of 77,138 differ); tails = published (226 differ); 19 records N-masked, no abort |
| DRS `correct` at the fix on the literal f15b523 BAM | identical to f15b523 in every column (42,627 rows) |
| DRS `correct` at the fix on the OLD c287748 BAM (master aborts there) | completes; 1 record N-masked; corrected_3prime identical to c287748 published for all 42,627; 52 junction / 51 clip-origin rows differ = the landing's own correct changes on old alignments |
| combined candidate `6e0b9ed` (CI fix + cDNA fix), full suite, Sherlock old stack | 4,306 passed; failure set identical to f15b523 (2 env tests) |
| `6e0b9ed` full suite under CI's stack (H2 venv: py3.12.9, matplotlib 3.11.2, pydeseq2 0.5.4) | 4,076 passed / 161 skipped; 2 H2-only artifacts (SGE NSLOTS in the loky test; no Liberation fonts) |

## 9. rbrowse exposure (checked 2026-10-10 02:50): none
rbrowse draws cDNA from `analyze/consensus_tagged.bam` records (raw alignments; `slice_locus.parse_alignment` treats `=` as a
match) and never joins cDNA `corrected_reads.tsv` (BUILD_MANUAL §7: "Do not fall back to corrected_reads.tsv for cDNA"; the
`aa-cdna` slug, 2026-09-24: "no corrected-TSV join"). Its `(strand, 3′-anchor)` molecule buckets are built from the alignments.
So no rbrowse view showed cdna-analyze's unwalked `corrected_3prime`, and the reviews Kevin did there are unaffected.

## 10. Production blast radius (2026-10-10 03:10; `cdna_eq_production.py`, H2 15120227; 10 % QNAME-hash samples, pre-fix code)

The smoke (one contig) UNDERSTATED it. cdna-analyze's own function, stored vs decoded, on three production libraries, split by
whether the walked-back end lies in an annotated intron on the molecule's strand:

| library (consensus_tagged.bam) | sampled molecules | other 3′ ends unwalked | intronic 3′ ends unwalked (n) |
| --- | ---: | ---: | ---: |
| 1076 tun_poly_bc03 (polysome) | 882,227 | 40.3 % | 56.6 % (555) |
| 1088 rrp6aa_rep2 (depletion cohort) | 1,347,203 | 28.7 % | 23.4 % (4,262) |
| 1082 wt_30min_rep1 (cold shock) | 725,186 | 27.8 % | 32.8 % (250) |

Shift sizes are alike in all three: about three quarters of the moved ends move 1–2 nt, ~16 % 3–4 nt, ~8 % 5–10 nt, ~0.5–1 %
more than 10 nt; never downstream. 99.8–99.9 % of the records are `=`-encoded.
