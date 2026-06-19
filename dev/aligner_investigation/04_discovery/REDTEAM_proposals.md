# RED-TEAM: Adversarial Review of the 6 Algorithmic-Improvement Proposals

**Reviewer role.** Hard-nosed adversarial reviewer. Goal: stress-test the six discovery proposals
for over-optimism, hidden costs, redundancy, conflicts, and whether they would *actually work* —
not to cheerlead.

**Inputs reviewed.** All six proposals in `04_discovery/`
(`improve_{seeding_indexing,chaining_dp,splice_junctions,3prime_cpa,ensemble_selection,ml_learned}.md`);
verified constraints from `02_synthesis/{COMPARISON.md,DEEP_DIVE.md}` and
`03_adversarial/redteam_winrates_selection.md`; and **the actual RECTIFY source** (every code claim
below was checked, not taken on the proposals' word).

> **Build note:** the code claims below were checked against the **`master`-derived working tree**;
> see `../CORRECTIONS_vs_DRS_BUILD.md` for re-verification vs `origin/drs-validation-rebuild`. **Two
> environment/code claims are SUPERSEDED on the build** and are corrected inline: (1) "the empirical
> penalty/STR tables are absent" was a **master-checkout artifact** — the tables are **bundled** under
> `rectify/data/genomes/*/penalty_tables/` and `--Scer`-auto-resolved (C13 below); (2) Path A
> (HP-edit-distance) is **wired and runs in production** (both merge call sites pass raw BAMs + genome),
> so the "wire the dead path on" framing (C1/C6, §1, roadmap E0/E1) is largely DONE. The Numba-import
> sub-claim (C11) is environment-dependent and unchanged.

**Environment note that colours several verdicts.** In this checkout: `edlib` imports; **`mappy`,
`parasail`, `lightgbm`, `scikit-learn`, and `numba` do NOT import**; `minimap2` is not on `PATH`;
the empirical `penalty_scores.tsv` / `str_penalty_scores.tsv` and the entire
`common/scripts/nanopore/` tree are **absent** *(in this **master** checkout)*. This matters: multiple
proposals assert "the Numba DP kernel already exists, only the call site is new" or "LightGBM is a
pip-install with no dependencies of consequence" — both understate real plumbing cost here.

> **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** the "penalty tables absent" half of this note
> is a **master-checkout artifact only.** On the build the tables ARE present and `--Scer`-resolved:
> `rectify/data/genomes/{saccharomyces_cerevisiae,homo_sapiens}/penalty_tables/{penalty_scores.tsv,
> penalty_scores_cdna*.tsv, penalty_scores_qsrev.tsv, str_penalty_scores.tsv, junction_overhang_table.tsv}`,
> auto-filled by `rectify/data/__init__.py:1188-1208` via the CLI hook `cli.py:198-199`. The
> **Numba-not-importable** sub-claim is genuinely environment-specific and is left standing.

---

## 1. Verified-vs-false code claims (the load-bearing ones)

| # | Claim (proposal) | Verdict | Evidence checked |
|---|---|---|---|
| C1 | Production selection runs the **legacy 5-level sort** (Path B); HP-edit-distance (Path A) is gated by `use_hp_ed = bool(per_aligner_corrected_bams)` and never wired | **TRUE** | `corrected_consensus.py:662` `use_hp_ed = bool(...)`; Path B key at `:749-752` (`_five_rescued,_chimera_ok,_conf_rank,_n_agree,_span,_n_junc`); Path A key at `:744-746` |
| C2 | `_n_agree` is a **popularity vote** (count of rows sharing `(read_id, corrected_3prime)`) | **TRUE** | `:693-698` `groupby(...).size()→_n_agree`; it is the dominant non-flag Path-B discriminator |
| C3 | **N-ops cost 0** ("free pass") in HP-edit-distance | **TRUE** | `_cigar_hp_edit_distance` `:68-69` `op == 3: ref_pos += length  # free pass`; docstring `:12` "N (intron): 0 — free pass" |
| C4 | `single_sample.py:495` passes **no BAMs / no genome / no penalty table** to the merge | **TRUE** | `:495-500` passes only `per_aligner_tsvs,output_tsv,summary_tsv,overhang_table` |
| C5 | `split_command.py:985` builds `corrected_bams` but does **not** feed them to the merge | **TRUE** | builds dict at `:970-978`, merge at `:985` omits it; uses BAMs only for the consensus BAM at `:1000-1014` |
| C6 | Per-aligner corrected BAM (`{stem}.rectified_corrected_3end.bam`) **already produced on disk** at the single-sample call site, "simply never collected" | **TRUE (with caveat)** | `stages.py:194,205` passes `corrected_bam=corrected_bam_path` into `correct_command.run`; but `_run_correction_per_aligner` (`:240-334`) **returns only TSV paths** — the BAM is written but not globbed back. The ensemble-selection §1 "wiring not computation" framing is accurate. Caveat: BAM-write is conditional on correction succeeding and on `correct_command` honouring `corrected_bam` for every aligner. |
| C7 | `--splice-flank=no` is a "compatibility" flag, **not** a proven 3'-accuracy mechanism | **TRUE** | `multi_aligner.py:252` literal comment `# Disable for compatibility` (contradicts CLAUDE.md's "important for 3' end accuracy"). chaining_dp's caveat is correct. |
| C8 | deSALT runs with **no `-x` ONT preset** (null ~13% model) | **TRUE** | `run_desalt` (`:1506+`) builds `aln`+index+`extra_args`; no `-x` hardcoded |
| C9 | mapPacBio sets `intronlen=50` but **no `maxindel`** | **TRUE** | `run_map_pacbio` `:513-515` sets `fastareadlen,intronlen=50,minratio=0.4`, no `maxindel`. (`run_bbmap` short-read at `:674` *does* set `maxindel=100000` — do not confuse the two.) |
| C10 | dorado `pt:i` tail length **already parsed/stored** | **TRUE** | `drs_trim_command.py:240,411,462` (`get_tag('pt')`, `pt_tag`, parquet column) |
| C11 | RECTIFY already owns `_score_hp_dp_numba` (Numba NW DP), `HpPenaltyTable.del_cost/ins_cost`, `local_aligner.py` (Gotoh affine) | **TRUE in code; PARTIALLY MISLEADING at runtime** | All three exist (`hp_penalty.py:44,261,270`; `local_aligner.py:21-27`). **But** `numba` is not importable here, so `_score_hp_dp_numba is None` → Module 2H runs the **pure-Python fallback**. "Only the call site is new" hides a real perf dependency. *(Build note: the Numba-availability concern is environment-specific and still stands; the separate "penalty tables absent" claim it was bundled with — see C13 — does NOT hold on the build.)* |
| C12 | `extra_args` passthrough exists on the minimap2 wrapper (so flag sweeps are zero-code) | **TRUE** | `multi_aligner.py:269-270` `if extra_args: cmd.extend(extra_args)` |
| C13 | Empirical penalty table exists on disk to consume | **FALSE on master · SUPERSEDED (TRUE) on build** | On the master checkout `penalty_scores.tsv` / `str_penalty_scores.tsv` and `common/scripts/nanopore/` are **absent**. **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** the tables ARE bundled — `rectify/data/genomes/{saccharomyces_cerevisiae,homo_sapiens}/penalty_tables/{penalty_scores.tsv, penalty_scores_cdna*.tsv, penalty_scores_qsrev.tsv, str_penalty_scores.tsv, junction_overhang_table.tsv}` — and `--Scer` auto-resolves them (`rectify/data/__init__.py:1188-1208`; `cli.py:198-199`); the generator `scripts/calibration/empirical_cigar_error_profiler.py` is also present. So on the build, proposals "reusing the empirical table" consume the **bundled** table; the open work is **validate/version**, not regenerate. |
| C14 | `--aligner-bams` wired into `rectify correct` enables Module 2H | **TRUE** | `correct_command.py:304,451,500` (`_run_2h` gated on `aligner_bams` + annotation) |
| C15 | NET-seq tables bundled (`netseq_{wt,pan}`, `atract`) | **TRUE** | `rectify/data/saccharomyces_cerevisiae_netseq_{wt,pan}.tsv.gz`, `..._atract_netseq.tsv.gz` |

**Net:** the proposals' factual base was unusually solid *against master* — the win-rate/selection
critique (C1-C6) was code-verified there, and the "RECTIFY already has the machinery" claims (C11, C14)
are real. **On the build, two of these are SUPERSEDED:** C1/C4/C5/C6 (Path A not wired) — both merge
call sites now pass raw BAMs + genome → `use_hp_ed=True` → **Path A runs in production** (see
`CORRECTIONS_vs_DRS_BUILD.md` §A); and C13 (tables absent) — the tables are **bundled and
`--Scer`-resolved**. The remaining soft spot is **C11 runtime** (numba absent → Python fallback), which
is environment-specific and unchanged. C2 (`_n_agree` popularity) and C3 (N=0) stay TRUE on the build —
but C2 describes the **fallback** path, not the production default, and C3's "only defense is binary" is
yeast-specific (human adds a junction-anchor gate).

---

## 2. Per-proposal verdict table

Legend: **KEEP** (sound, build it) · **KEEP-WITH-CAVEATS** · **DESCOPE** (good kernel, the
proposal over-reaches) · **REJECT**.

### File: improve_ensemble_selection.md

| Sub-idea | Verdict | One-line reason |
|---|---|---|
| §5 Win-rate harness + provenance + Path A/B ablation | **KEEP** | Cheapest, unblocks every accuracy claim in all six files. The single correct "do-first." |
| §1 Wire Path A (collect existing corrected BAMs, flip `use_hp_ed`) | **MOSTLY DONE on build** | **⚠️ BUILD CORRECTION:** this wiring already exists on the build via the lazy raw-BAM path — both call sites pass `per_aligner_raw_bams` + `genome` and `use_hp_ed = bool(per_aligner_corrected_bams or per_aligner_raw_bams)` (`corrected_consensus.py:1262`; `single_sample.py:238-250`; `split_command.py:1085-1094`). Remaining open work is **only the N-op calibrated cost** (C3 / splice P4), so Path A doesn't trade `_n_agree` herd bias for the N=0 free-intron exploit. Measure on §3 truth. |
| §2 Correlation-aware / lineage-weighted consensus (de-herd `_n_agree`) | **KEEP-WITH-CAVEATS** | The lineage insight (3 of 5 share minimap2 lineage) is real and the strongest novel idea in this file. But hand-set weights `{minimap2:0.5,...}` are themselves an unvalidated prior; only ship empirical-ρ weights from §3, and only as a tiebreaker term. |
| §3 Ground-truth-free per-aligner weighting (NET-seq + replicate + Dawid–Skene EM) | **KEEP-WITH-CAVEATS** | NET-seq concordance is the defensible signal; **but the EM annotator-reliability arm is circular** (it can converge to "trust the herd" — the exact bias being removed). Keep NET-seq + replicate; treat EM as advisory only, seeded by NET-seq, self excluded. Replicate concordance needs replicates that may not exist. |
| §4 Calibrated per-read confidence + abstain | **KEEP-WITH-CAVEATS** | High downstream value (APA), but calibration target is NET-seq → inherits §3's circularity; needs `sklearn.isotonic` which **does not import here**. Real, not free. |
| §6 Learned GBDT selector | **DESCOPE → advisory-only** | Duplicates ml_learned §1. See cross-file dedup. LightGBM not importable here; 36-read overfit risk; keep offline/advisory until it beats a *well-wired* Path A. |

### File: improve_ml_learned.md

| Sub-idea | Verdict | One-line reason |
|---|---|---|
| §1 Learned winner SELECTOR (GBDT, NET-seq labels) | **DESCOPE** | Same idea as ensemble §6; honest about overfit and "measurement first." But: (a) it must beat a wired Path A (§1) *which does not exist yet* — so it's premature; (b) LightGBM/sklearn absent here; (c) NET-seq labels are also the validation oracle → leakage unless held out by chromosome AND dataset. Build the heuristic Path A first; treat GBDT as a ceiling-probe, not a deliverable. |
| §2a `pt:i` as Bayesian prior in 2E walk-back | **KEEP** | Genuinely cheap, data already ingested (C10), orthogonal information. Strongest item in the file. |
| §2b Remora/dwell signal fusion | **REJECT (now) / defer to funded R&D** | Requires POD5 retention RECTIFY does not keep + signal↔ref anchoring + per-read signal I/O against the chunked-SLURM budget. Author concedes this. Not a discovery deliverable. |
| §3 SpliceAI/Pangolin junction scorer | **DESCOPE → metazoan-only, deferred** | Author correctly says NIL for yeast (~300 canonical introns). Keep only as the interface stub (shared with splice_junctions P5). |
| §4 full seq2seq/transformer correction | **REJECT** | Author agrees: 36-read training set, opaque, replaces auditable DP. |
| §4 narrow "learn the penalty table" | **KEEP-WITH-CAVEATS** | Defensible — but it is just *recalibrating the empirical table* (C13), which must be regenerated anyway. Fold into the table-regeneration work, not a separate ML track. |
| §5 foundation-model embeddings | **REJECT** | Author agrees; GPU cost violates SLURM design; one speculative chimera aside. |

### File: improve_3prime_cpa.md

| Sub-idea | Verdict | One-line reason |
|---|---|---|
| P1a Unified terminal re-scoring pass (per-aligner, pre-correction CPA solver) | **KEEP** | Highest-value, mostly consolidation of existing 2B/2E/2G into one boundary decision. **But see CONFLICT-1**: it claims to *collapse* 2E/2G whose priority ordering CLAUDE.md fought hard to get right (v2.9.x–v3.0.x). "Decide the boundary once" is the right idea but is a rewrite of a heavily-patched, validated code path, not a cheap add. Risk understated. |
| P1b true DP terminal state (KSW2/BBMap fork) | **REJECT (now)** | Author defers it; correct. mappy not importable, forks 5 aligners. |
| P2 Phase A `pt:i` fusion | **KEEP** | Same as ml_learned §2a — **dedup**. Orthogonal anchor, cheap. |
| P2 Phase B dwell fusion | **REJECT (now)** | Same as ml_learned §2b — **dedup**. POD5 plumbing. |
| P3 NET-seq soft prior in selection | **KEEP-WITH-CAVEATS** | Same as ensemble §3 selection-prior — **dedup**. Circularity is the whole risk; held-out `_wt`/`_pan` split is mandatory and is itself only a partial fix (see DEPENDENCY note). |
| P4 3'-tight / 5'-loose DP asymmetry | **KEEP-WITH-CAVEATS** | Mostly defensive consolidation of Cat3; low impact on the *3'* metric by the author's own admission. Worth it only as hygiene, ranked last. |
| P5 mispriming terminal veto (cDNA) | **KEEP** | Reuses the existing genomic-downstream-A detector; correctly DRS-exempt; bounded penalty not gate. Genuinely useful for the QuantSeq path and does not touch the contested win-rate machinery. |

### File: improve_chaining_dp.md

| Sub-idea | Verdict | One-line reason |
|---|---|---|
| P0 ground-truth harness + regenerate penalty table | **KEEP** | Same harness as everyone else (**dedup**); plus the table regeneration (C13) which is a real, non-trivial prerequisite the other files gloss. |
| P3 poly-A-aware two-state terminal DP (CPA = G→T transition column) | **KEEP-WITH-CAVEATS** | The most *principled* idea across all six files — replaces the 2B/2E/2G heuristic stack with one optimal boundary. **But it is the SAME target as 3prime_cpa P1a** (one is a 2-state DP, the other a re-scoring pass — they CONFLICT, see below) and it depends on the numba kernel that is **not compiled here** (C11). Pick ONE of {P3 two-state DP, P1a re-scoring}; do not build both. |
| P2 HP-aware full-SW realign of 3' window | **KEEP-WITH-CAVEATS** | `edlib` is importable (true); window is cheap. But it overlaps P1a/P3 heavily — it is a third spelling of "re-place the terminal indel with a better DP." Consolidate. |
| P1b empirical gap costs into realign DP | **DESCOPE** | Depends on the absent table (C13) and the absent numba kernel; chemistry-specific. Fold into P3 once the table exists. |
| P1a KSW2 fork | **REJECT** | Author defers; correct. |
| P4a harvest sub-optimal terminal chains (`--secondary=yes -N`) | **KEEP-WITH-CAVEATS** | Cheap flag change (C12), genuinely addresses "post-correction can't recover an unseeded exon." But more candidates → more chances for a wrong terminus to win Path B's popularity sort — **must** land after §1/§5, never before. |
| P4b 3'-end-anchored chaining bonus | **DESCOPE** | No minimap2 flag exists; needs a patch or RECTIFY re-chain. High risk (can pull 3' into a genomic A-run), fights P3. Defer. |
| P4c RMQ-vs-range for short reads | **KEEP (hygiene)** | Already implicitly true (separate short-read aligners); document + verify band. |
| P5 explicit cross-read CPA consensus (replace `_n_agree`) | **KEEP-WITH-CAVEATS** | Overlaps splice_junctions P1 + ensemble §2. **Dominant risk (author-flagged): over-homogenization erases real APA heterogeneity** — the entire point of a 3' study. Snap only within HP-noise radius. Cross-chunk peak pooling under SLURM arrays is non-trivial (not "plumbing"). |
| P6/P7/P8 deSALT `-x`, mapPacBio `maxindel`, minimap2 `--end-bonus` sweep | **KEEP** | All one-line config A/Bs (C8/C9); P8 explicitly diagnostic-only (do not ship a blunt end-bonus — it fights P3). |
| P9 Z-drop terminal model | **REJECT (subsumed)** | Author says "no separate work." |

### File: improve_splice_junctions.md

| Sub-idea | Verdict | One-line reason |
|---|---|---|
| P1 Cross-read junction consensus (first-class step) | **KEEP-WITH-CAVEATS** | Generalizes deSALT's mechanism to all 5 aligners — real. **Same over-homogenization risk as chaining_dp P5** (collapses genuine alt-5'/3'SS isoforms <50 bp apart, which exist in yeast RP genes). The per-read sequence veto + shared-boundary-anchoring mitigations are necessary, not optional. |
| P2 Unified calibrated junction scorer (MaxEnt + HP + cross-read + annotation prior, logistic fit) | **KEEP-WITH-CAVEATS** | The right direction (replace the brittle lexicographic tuple with a calibrated score). **But:** (a) MaxEnt is metazoan-trained → must retrain a yeast PWM; (b) `sklearn` for the logistic fit **does not import here**; (c) the calibration set (≥4-aligner concordance) is itself herd-defined → partial circularity; (d) wiring it into selection collides with the C1 problem — it only helps "if wired into selection," same dependency as ensemble §1. |
| P3 De-novo micro-exon recovery (≤10 nt) | **DESCOPE → metazoan** | Author honest: low yeast value (few micro-exons), no validation read exercises it, NET-seq can't validate it (3'-assay). Needs P2 first. Defer. |
| P4 Calibrated novel-junction handling (charge N-ops their improbability) | **KEEP** | Directly fixes C3 (N=0 free-intron exploit) the *right* way (calibrated cost, not a fixed gate). Mostly free given P2. Highest correctness-per-effort in this file. |
| P5 SpliceAI/Pangolin CNN | **DESCOPE → metazoan, interface-only** | Same as ml_learned §3 — **dedup**. Build the swappable interface, defer the CNN. |

### File: improve_seeding_indexing.md

| Sub-idea | Verdict | One-line reason |
|---|---|---|
| P0 accuracy oracle + seed-recall harness | **KEEP** | Same harness (**dedup**). Adds `seed_exon_recall` via minimap2 `--print-seeds` (binary not present here — a real setup cost). |
| 6a deSALT `-x ont2d` A/B | **KEEP** | = chaining_dp P6 (**dedup**). One-line, top-aligner, closes an open question. |
| 1a minimap2 `-k`/`-w` sweep | **KEEP** | Zero-code via `extra_args` (C12). Honest "syncmer-approx," good cheap win. |
| 5a CPA-anchor BED (3'-end-aware corrector seeding) | **KEEP-WITH-CAVEATS** | Overlaps 3prime_cpa P3 selection-prior (**partial dedup**). Circularity is the dominant risk; the **novel-CPA validation arm is mandatory** and the author says so. |
| 2 Scoped HPC re-seeding of poly-A-proximal soft-clips | **KEEP-WITH-CAVEATS** | Genuinely novel; but HPC anchors are ambiguous in A-rich yeast regions → the local-window constraint is load-bearing and the back-mapping (run-length↔coord) is fiddly. Medium risk, medium reward. |
| 3 De-novo strobemer aligner (namfinder, poly-A-preserving) | **DESCOPE** | A from-scratch chainer = the same correctness tail the dossiers catalogue for gapmm2. Yeast micro-exon upside is small (F2: real floor ~6 nt, rare in yeast). Research bet, not a near-term win. |
| 1b standalone open-syncmer seeder | **DESCOPE** | From-scratch aligner; gate behind 1a proving insufficient. |
| 4 Unified de-novo + annotation/CPA-prior graph index | **REJECT (as a graph build) / KEEP the layered approximation only** | Author already concedes the monolithic graph is a multi-month deSALT-fragility bet; the "layered ensemble + shared prior file" is just ensemble glue (= 5a + P5 consensus), not a new index. Do not build the graph. |
| 6b R10.4.1-learned spaced seeds | **REJECT** | Speculative, overfit-prone, needs a custom seeder that doesn't exist. Author ranks it last; agree. |

---

## 3. Cross-proposal analysis

### 3.1 REDUNDANCY (consolidate — these are the SAME idea wearing different hats)

1. **The ground-truth / NET-seq concordance harness appears in ALL SIX files** (ensemble §5, ml §"precondition", 3prime_cpa "common harness", chaining P0, splice "cross-cutting prerequisite", seeding P0). **Build it ONCE.** It is the universal "do-first." But see the circularity warning in 3.3.

2. **Poly-A-aware terminal DP / CPA boundary solver appears THREE times:**
   `improve_chaining_dp.md` **P3** (two-state G→T DP),
   `improve_3prime_cpa.md` **P1a** (terminal re-scoring pass),
   `improve_chaining_dp.md` **P2** (HP-aware window realign).
   These are three spellings of "re-place the terminal indel/boundary with a better, poly-A-aware DP."
   **Consolidate into ONE terminal CPA solver.** (They also mildly conflict — see 3.3.)

3. **`pt:i` signal fusion appears twice:** 3prime_cpa **P2 Phase A** ≡ ml_learned **§2a**. One item.

4. **Cross-read junction/CPA consensus appears THREE times:** splice **P1**, chaining **P5**,
   ensemble **§2** (the de-herding is the same operation viewed from the selection side). One module.

5. **Learned GBDT selector appears twice:** ensemble **§6** ≡ ml_learned **§1**. One item.

6. **NET-seq-as-selection-prior appears three times:** 3prime_cpa **P3**, ensemble **§3**, seeding **5a**.
   One soft-prior mechanism, applied at one site (the merge tiebreaker / terminal solver).

7. **SpliceAI/Pangolin interface appears twice:** splice **P5** ≡ ml_learned **§3**. One deferred stub.

8. **deSALT `-x ont2d` A/B appears twice:** chaining **P6** ≡ seeding **6a**. One config experiment.

**After dedup, the six files reduce to ~9 distinct work items** (see roadmap in 3.4).

### 3.2 CONFLICTS (proposals that fight each other or validated policy)

- **CONFLICT-1 — terminal-DP rewrite vs the validated 2E/2G priority dance.** 3prime_cpa P1a and
  chaining P3 both propose to *collapse* Modules 2B/2E/2G into one boundary call. CLAUDE.md documents
  an extensive lineage of bug-fixes that got the 2G-takes-priority-over-2E ordering, terminal-D/N
  stripping (v3.0.2/3.0.3), and `find_polya_boundary` N-op guards (v2.9.3/2.9.4) exactly right on real
  data. "Decide the boundary once" is architecturally cleaner but is a **rewrite of the most heavily
  patched code path in the project**, validated only against the 36-read set + cat1/cat2 shifts. The
  proposals frame this as "mostly consolidation / low effort." It is **high regression risk**, not low.

- **CONFLICT-2 — P3 two-state DP vs P2 window-realign vs P1a re-scoring** all claim the #1/#2 slot and
  all rewrite the terminal boundary differently. Building two of them double-corrects (the same class
  of mutual-cancellation bug they claim to remove). **Pick one.**

- **CONFLICT-3 — cross-read homogenization (splice P1 / chaining P5) vs biological APA signal.** This
  is the single most dangerous conflict in the set: the deliverable is *alternative polyadenylation*,
  and over-snapping reads to a consensus CPA/junction **destroys the exact heterogeneity being
  measured.** It is also self-reinforcing with the `_n_agree` herd bias (snap → more agreement →
  higher popularity → more snapping). The HP-noise-radius bound is not a nice-to-have; without it the
  feature is actively harmful.

- **CONFLICT-4 — any annotation/NET-seq/CPA prior vs the PERMANENT "canonical-as-tiebreaker-never-gate"
  policy.** ensemble §2 weights, 3prime_cpa P3/P5, splice P2 `w_ann`, seeding 5a/4 priors all risk
  re-introducing a gate. CLAUDE.md marks the no-gate policy PERMANENT (v3.1.7). Every one of these
  must stay a bounded tiebreaker within the HP-noise floor. The proposals all *say* this; the risk is
  in implementation drift, and it needs an explicit test (a novel-site read that a gate would
  mis-snap) for each prior.

- **CONFLICT-5 — chaining P8 blunt `--end-bonus` vs P3 poly-A DP.** A raw end-bonus extends into
  genomic A-runs (the opposite failure). Author flags it diagnostic-only; ensure it is never shipped
  as a fix alongside the poly-A solver.

- **No conflict with chunked-SLURM I/O policy** in any proposal — but chaining P5 / splice P1's
  cross-chunk peak pooling is described as "plumbing" when it actually requires a global merge barrier
  before pass 2 (like deSALT's batch), which complicates the SLURM-array DAG. Cost understated, not
  conflicting.

### 3.3 DEPENDENCY ORDER — and why "ground-truth harness first" is necessary but NOT sufficient

Every file says "build the NET-seq concordance harness first." That is correct as far as it goes, but
**the near-unanimous framing hides a second circularity the redteam already flagged and the proposals
under-weight:**

- **NET-seq is used as BOTH the training/tuning signal AND the validation oracle** in: ml §1 (GBDT
  labels), ensemble §3/§4 (weights + confidence calibration), 3prime_cpa P3 / seeding 5a (selection
  prior), splice P2 (logistic fit features). If the same NET-seq table tunes the model *and* scores it,
  concordance is guaranteed to "improve" by construction.
- The only partial fix offered is held-out `_wt` vs `_pan`. But `_wt` and `_pan` are **not independent
  experiments** — `_pan` is a pan/aggregate table over conditions; they share the same assay biology
  and likely the same peaks. Holding one out does **not** break the circularity the way a truly
  orthogonal assay (Quant-seq / 3'-seq / direct cDNA 3'-seq) would.
- **Worse: there is no Quant-seq/3'-seq truth set committed in the repo** — only NET-seq. So the
  "orthogonal CPA truth set" the entire program is predicated on **does not currently exist** beyond
  NET-seq, and NET-seq is also what RECTIFY *already refines against* (`netseq_refiner.py`), i.e. it is
  partly baked into the corrected positions being scored.

**Therefore the truly correct "do-first" is two-pronged and bigger than any file states:**
  1. **(a)** the win-rate provenance + Path A/B ablation (ensemble §5 / redteam exp 1-2) — cheap, pure
     code, no oracle needed. This alone settles whether the win rates are a metric artifact and whether
     to wire Path A. **This is the genuine cheapest unblocker.**
  2. **(b)** a **non-NET-seq** CPA truth set (Quant-seq / 3'-seq, or at minimum a held-out NET-seq
     sample from a *different run*, plus annotated-intron junction truth for spliced reads). Until (b)
     exists, NET-seq concordance is a sanity check, not an oracle, and **no ML/weight/prior proposal is
     creditable.** The proposals treat (b) as "curate the truth set" (low effort); it is the real
     gating cost of the whole program.

### 3.4 Corrected, dependency-ordered roadmap (de-duplicated)

```
TIER 0 — measurement, no oracle needed (DO FIRST, days)
  [E0] Win-rate provenance + Path A/B ablation harness        (ensemble §5; redteam exp 1-2,5)
       → answers "is the spread a metric artifact?" with ZERO accuracy assumptions.
  [E1] Wire Path A (collect existing corrected BAMs, C6) BUT only behind a flag,
       SHIPPED TOGETHER WITH the N-op calibrated cost (C3 fix = splice P4)         (ensemble §1 + splice P4)

TIER 1 — build the REAL oracle (gating cost, weeks; the program's true critical path)
  [O0] Non-NET-seq CPA truth set (Quant-seq/3'-seq or cross-run held-out) +
       annotated-intron junction truth + seed-exon-recall instrumentation.
       Regenerate the absent penalty table (C13) as part of this.                  (all P0s, merged)
  *** Until [O0] exists, everything below is a hypothesis. ***

TIER 2 — cheap, orthogonal, low-risk wins (config + already-ingested data)
  [C1] deSALT -x ont2d A/B; mapPacBio maxindel; minimap2 k/w sweep                 (chaining P6/P7 + seeding 1a)
  [C2] pt:i as bounded prior in 2E walk-back                                       (ml §2a / 3prime P2A)
  [C3] mispriming terminal veto for cDNA (reuses existing detector)               (3prime P5)

TIER 3 — the ONE terminal CPA solver (pick a single design; HIGH regression risk)
  [T1] Unified poly-A-aware terminal boundary solver — choose EITHER the
       two-state G→T DP OR the window re-scoring pass, not both. Gate on numba
       being available (C11). Validate cat1/cat2 shifts EXACTLY before merge.      (chaining P3 ≡ 3prime P1a ≡ chaining P2)

TIER 4 — calibrated selection/junction quality (depends on [O0] + [E1])
  [S1] Calibrated junction scorer (yeast PWM/MaxEnt, NOT sklearn-dependent fit)    (splice P2)
  [S2] Cross-read junction/CPA consensus — WITH hard APA-preservation radius       (splice P1 ≡ chaining P5 ≡ ensemble §2)
  [S3] Empirical lineage/reliability weights (NET-seq + replicate; NO EM herd)     (ensemble §3, de-scoped)

TIER 5 — ceiling probes / deferred (advisory only or metazoan)
  [X1] GBDT selector — advisory until it beats wired Path A on [O0]                (ml §1 ≡ ensemble §6)
  [X2] Calibrated per-read confidence + abstain                                    (ensemble §4)
  [X3] SpliceAI interface stub (metazoan); scoped HPC re-seeding (research)        (splice P5/ml §3; seeding P2)

REJECT outright: KSW2 fork, full transformer correction, FM embeddings, monolithic
  graph index, standalone syncmer/strobemer aligners, Remora/POD5 dwell fusion (now),
  R10.4.1 spaced seeds, blunt --end-bonus as a fix.
```

### 3.5 Genuinely novel/high-value vs repackaged existing behaviour

- **Genuinely novel & high-value:** the **N-op calibrated cost** (splice P4 — closes a real,
  code-verified exploit C3); the **poly-A-aware two-state terminal DP** (chaining P3 — the only
  *principled* replacement for the heuristic 2B/2E/2G stack); **lineage-weighted de-herding** of
  `_n_agree` (ensemble §2 — the lineage observation is correct and the fix is sound); the **mispriming
  terminal veto** (3prime P5 — useful, self-contained, low-risk).
- **High-value but mostly WIRING of existing behaviour (not new science):** wiring Path A
  (ensemble §1 — the BAMs already exist, C6); `pt:i` prior (data already parsed, C10); the harness
  (standard benchmarking).
- **Repackaged / over-claimed:** the "unified graph index" (seeding P4 — author concedes the shippable
  form is "ensemble glue, not a new index"); cross-read consensus is **deSALT's existing mechanism**
  re-implemented post-hoc (real value, but not novel); the GBDT/SpliceAI/foundation-model items are
  textbook methods whose authors mostly already mark them overkill for yeast.

---

## 4. THE SINGLE BIGGEST RISK

**The entire program is validated against an oracle (NET-seq) that is simultaneously (a) the tuning
signal for nearly every proposed model/weight/prior, and (b) already partially baked into the
corrected positions being scored (`netseq_refiner.py`), and (c) the *only* CPA truth set in the repo
— there is no committed orthogonal (Quant-seq/3'-seq) ground truth.** Holding out `_wt` vs `_pan` does
not break this because they are not independent experiments.

Consequence: **every accuracy claim — including the "deSALT 78.9%" the whole investigation is built
on — can be made to "improve" by construction, and the proposals' near-unanimous "ground-truth harness
first" gives false comfort because the ground truth is circular.** The correct sequencing is: do the
**zero-oracle** Path A/B + provenance ablation FIRST (it needs no truth set and settles the metric
question), and treat **procuring a genuinely orthogonal CPA truth set as the real, weeks-long gating
investment** — not the cheap "curate the truth set" line every file assumes. Any proposal that trains
on NET-seq and reports NET-seq concordance as its win is, until then, unfalsifiable.

Secondary risk, called out because it is mis-rated as "low effort": **CONFLICT-1/CONFLICT-3** — the
terminal-DP unification rewrites the most heavily-patched, real-data-validated code path in RECTIFY,
and the cross-read consensus can erase the APA signal that is the product's reason to exist. Both are
framed as "consolidation"; both are high-regression-risk and must gate on cat1/cat2 + a novel-site
APA-preservation test, not just the 36-read green suite.
