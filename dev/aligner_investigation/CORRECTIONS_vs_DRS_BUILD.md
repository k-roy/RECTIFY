# CORRECTIONS — Aligner Investigation re-verified against the REAL build (`origin/drs-validation-rebuild`)

**Why this file exists.** The investigation in `dev/aligner_investigation/` verified its
load-bearing CODE claims against the `master`-derived working tree. The build the team
actually runs is the git branch **`origin/drs-validation-rebuild`** (224 commits ahead of
the investigation's reference; hundreds of lines different in the key files). Several of the
most prominently-stated findings — including the single most-emphasized one ("HP-edit-distance
selection never runs in production") — are **now wrong on the real build** and must be
retracted or softened in the deliverables.

**Method.** All build code was read with `git show origin/drs-validation-rebuild:<path>`
(working tree untouched). Every verdict below cites build `file:line` + snippet.

**Verdict legend.** CONFIRMED = build matches the investigation's claim · CHANGED = claim was
true on master but the build differs materially · OUTDATED-FALSE = claim is now factually
wrong on the build and must be retracted. **CONFIRMED/INFERENCE** distinction is preserved
from the source docs.

---

## Claims table

| # | Claim (as stated in deliverables) | master verdict (per source) | BUILD verdict | Build evidence (`origin/drs-validation-rebuild`) |
|---|---|---|---|---|
| 1 | Empirical penalty / STR / overhang tables are **ABSENT from the checkout**; only the loader + docs exist → "heuristic fallback is what actually runs" (SPLICE_JUNCTION_PLACEMENT §2 L164-168; improve_splice_junctions J2) | TRUE on master | **OUTDATED-FALSE** | Tables ARE bundled: `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/{penalty_scores.tsv, penalty_scores_cdna.tsv, penalty_scores_cdna_umi1/2/3plus.tsv, penalty_scores_qsrev.tsv, str_penalty_scores.tsv, junction_overhang_table.tsv}` and a parallel `homo_sapiens/penalty_tables/` set. Generator present: `scripts/calibration/empirical_cigar_error_profiler.py`. |
| 1b | (implicit) tables are not auto-resolved by `--Scer` | n/a | **CHANGED — auto-resolved by default** | `rectify/data/__init__.py:1188-1208` `resolve_reference_paths()` fills `args.junction_penalty_table`/`str_penalty_table`/`junction_overhang_table` from bundled data when the slot is empty, **protocol-routed** (`_protocol_for_args` → drs/cdna/qsrev) via `get_bundled_junction_penalty_table` (`__init__.py:835`). The global CLI hook calls it: `rectify/cli.py:198-199`. So `--Scer` loads the DRS penalty+overhang tables by default. |
| 2 | **Both production call sites of `merge_corrected_tsvs` pass NO BAMs → `use_hp_ed=False` → legacy popularity-vote (Path B) runs. HP-edit-distance is "NOT-IN-PRODUCTION."** (redteam_winrates §a L56-75, claims 2/3/6; SPLICE §2 L189-214; "single most important experiment") | TRUE on master | **OUTDATED-FALSE — the headline finding is now wrong** | (a) Gate widened: `corrected_consensus.py:1262` `use_hp_ed = bool(per_aligner_corrected_bams or per_aligner_raw_bams)`. (b) NEW lazy HP path: passing `per_aligner_raw_bams` computes HP-ED in memory from raw BAM + TSV without materializing corrected BAMs (`corrected_consensus.py:1264-1364`). (c) **Single-sample call site passes raw BAMs + genome**: `rectify/core/commands/run/single_sample.py:238-250` — `merge_corrected_tsvs(..., per_aligner_raw_bams=_staged_s1 if not per_aligner_corrected_bams else None, genome=_merge_genome, ...)`. (d) **Chunked/split production call site passes raw BAMs + genome**: generated script body `rectify/core/commands/split_command.py:1085-1094` — `per_aligner_corrected_bams=...if corrected_bams else None, per_aligner_raw_bams=_staged if (raw_bams and not corrected_bams) else None, genome=genome`. → **Path A (HP-edit-distance) runs in production on the build.** |
| 3 | N-op is **cost 0** in `_cigar_hp_edit_distance` → "free false-intron exploit"; the overhang filter is the only (binary) defense (redteam_winrates §a note L77-81; SPLICE §2 item 2 L205-209; improve_splice_junctions J1) | TRUE on master | **CHANGED — N=0 confirmed, but a NEW junction-anchor gate backstops it** | N still free: `corrected_consensus.py:142-143` `elif op == 3: # N (intron skip) ... ref_pos += length # free pass`. NEW backstop: `_cigar_min_junction_anchor` (`:166-220`) measures the junction-local perfect-match anchor on both flanks of every N-op; `_add_chimera_flag` gates spurious-N winners at `:868-882` (`if anchor < min_junction_anchor_bp: ... row_is_chimeric = True` unless every junction is cross-read-supported). **BUT the gate is organism-gated**: `default_min_junction_anchor_bp` (`rectify/data/__init__.py:1053-1062`) = human 10, **yeast 0 (gate OFF)**, unknown 0; `_MIN_JUNCTION_ANCHOR = 0` (`corrected_consensus.py:650`). So on yeast the overhang filter is still the only defense; on human the anchor gate adds a graded structural defense. The "free intron / only binary defense" criticism holds for **yeast** but is **softened/false for human** and no longer describes the code as a whole. |
| 4a | deSALT runs with **NO `-x` ONT preset** (`null` ~13% model) and **no `-G`** (yeast GTF→SIGSEGV) (deSALT dossier; COMPARISON §6.1; redteam claim 10) | TRUE | **CONFIRMED** | `multi_aligner.py:2257-2266` cmd = `[desalt_exec, 'aln', '-t', ..., '-f', ..., '-o', ...]` — no `-x`. `-G` deliberately skipped: `:2264-2266` "deSALT's -G annotation flag causes a SIGSEGV ... Skip -G entirely". Crash handling + empty-BAM fallback `:2290-2331`; `_dedup_desalt_bam` present `:2131`. |
| 4b | minimap2 flags `-ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD --junc-bed --junc-bonus 9` | TRUE | **CONFIRMED** (one addition) | `multi_aligner.py:361-371` exactly these flags (`-G str(max_intron)`, `--splice-flank=no`, `--secondary=no`, `--MD`), plus `-y` (copy FASTQ comment tags). `--junc-bed`/`--junc-bonus` via `get_minimap2_junc_args` `:383-390`. |
| 4c | gapmm2 `-i 5000` (max intron) | TRUE | **CONFIRMED** | `multi_aligner.py:1403` `'-i', str(max_intron)`. (Note `-m 1` removed; min-mapq comment corrected `:1397-1402`.) |
| 4d | uLTRA `--ont --disable_infer` | TRUE | **CONFIRMED** | `multi_aligner.py:2066-2070` `[ultra_path, 'pipeline', '--ont', '--disable_infer', '--t', ...]`. |
| 4e | mapPacBio `intronlen=50`, **no `maxindel`** (dossier flagged the missing cap as a possibly load-bearing gap, redteam_denovo B5) | TRUE on master | **CHANGED** | `multi_aligner.py:749` now `intronlen=10` (NOT 50), and `:754` `f'maxindel={max(200000, max_intron)}'` — **the previously-missing maxindel cap is now set explicitly** (≥200 kb for human RNA). The "intronlen=50 / no maxindel" dossier statements are stale. |
| 5 | Module 2H scoring is sequence-first, bilateral `t1+t2`, `_CANONICAL_HP_PRIOR=0.5`, `tier_beats_alt` tuple flip, **no candidate guard** (PERMANENT v3.1.7), `search_radius=5000`, `max_boundary_shift=50` (SPLICE §2 L133-159) | TRUE | **CONFIRMED — materially unchanged** | `junction_scoring.py:293` `_CANONICAL_HP_PRIOR = 0.5`; `_score_junction` bilateral `score(k)=t1(k)+t2(k)` `:807-928` with `for k in range(L)` `:901`. `junction_refiner.py:462-663`: `search_radius=5000`, `max_boundary_shift=50` `:471-472`; `tier_beats_alt = current_tier >= 4` `:616`; `canonical_discount = _CANONICAL_HP_PRIOR if tier < 4 else 0.0` `:646-647`; tuple flip `:658-663`; no pre-scoring candidate guard. |
| 6 | `--splice-flank=no` code comment is **"Disable for compatibility"** (contradiction vs CLAUDE.md "important for 3' end accuracy") (SPLICE §1 L42-44) | TRUE | **CONFIRMED** | `multi_aligner.py:366` `'--splice-flank=no',  # Disable for compatibility`. Comment unchanged; the internal contradiction the dossier flagged still stands. |

---

## Findings that are NOW WRONG and must be retracted / softened in the deliverables

These are the load-bearing conclusions the investigation built on the master tree that **do not
hold on `origin/drs-validation-rebuild`**:

### A. (Claim 2 — the headline) "HP-edit-distance selection never runs in production; production runs the legacy popularity vote (Path B)." → RETRACT.
On the build, **both** production merge call sites pass per-aligner **raw** BAMs + genome, and
`merge_corrected_tsvs` gained a lazy HP-edit-distance path that activates on
`per_aligner_raw_bams` (no corrected-BAM materialization needed). `use_hp_ed` is therefore
**True** in production. Specifically the following statements are now false:

- redteam_winrates §a "Which path runs in production?" (L56-75) and its conclusion "the
  reported win rates were almost certainly produced by Path B (legacy sort)."
- redteam_winrates claim 2 ("LIKELY-BIASED / FACTUALLY WRONG: production sort never sorts on
  3'-end / HP-ED"), claim 3 ("HP-aware edit distance ... NOT-IN-PRODUCTION"), claim 6 (mechanism
  attribution premised on Path B).
- SPLICE_JUNCTION_PLACEMENT §2 "CRITICAL FINDING — junction quality does NOT drive production
  selection" (L189-214) and the table footnote L122-124 "measured under the legacy selection
  path, NOT a junction-quality metric."
- improve_splice_junctions / REDTEAM_proposals J1 step (b) "collect the already-produced
  per-aligner corrected BAMs ... and pass them to `merge_corrected_tsvs` so Path A actually
  runs" — **this wiring already exists on the build** (via the lazy raw-BAM path). J1(b) is
  largely DONE; only J1(a) (charging N-ops a calibrated cost) remains open.

> NOTE — *what survives*: the *legacy Path B code still exists* (`corrected_consensus.py:1387-1402`,
> `_n_agree` popularity term) and is still the fallback when no BAMs/genome are supplied (e.g.
> manual `merge_corrected_tsvs` calls, or if BAM staging fails and `ed_rows` is empty →
> `use_hp_ed=False`, `:1377`). So "Path B exists and is a popularity vote" is still TRUE as a
> *description of the fallback*; what is FALSE is "Path B is what production runs."

### B. (Claim 1) "Empirical penalty/STR/overhang tables are absent → only the heuristic fallback runs." → RETRACT.
The tables are bundled for both *S. cerevisiae* and *H. sapiens* with protocol variants, and
`--Scer` auto-resolves them by default (protocol-routed). SPLICE §2 L164-168 ("CRITICAL ...
ABSENT from this checkout ... heuristic fallback is what actually runs") and improve_splice
J2's premise ("currently absent") are wrong on the build. J2 should be re-scoped from
"regenerate + commit the missing tables" to "validate / version the already-bundled tables."

### C. (Claim 4e) mapPacBio "`intronlen=50`, no `maxindel`." → CORRECT THE NUMBERS.
Build uses `intronlen=10` and an explicit `maxindel=max(200000, max_intron)`. The
redteam_denovo B5 concern ("no maxindel — possibly load-bearing gap for long introns") is
RESOLVED on the build. Any sentence citing `intronlen=50` or "relies on BBMap's soft default
~16000 maxindel" is stale.

### D. (Claim 3, partial) "N-op cost 0 with only a binary overhang flag as defense." → SOFTEN, organism-qualify.
Still true for **yeast** (anchor gate defaults OFF). For **human** a graded junction-local
perfect-match **anchor gate** (`_cigar_min_junction_anchor` + `_add_chimera_flag` anchor branch,
default 10 bp) now backstops the free-N exploit. Statements that the overhang filter is "the
*only* defense" must be qualified "on yeast (anchor gate off by default); human enables a
structural anchor gate."

---

## Findings that STILL HOLD on the build (no change needed)

- **N-op is cost 0** in `_cigar_hp_edit_distance` (`corrected_consensus.py:142-143`). The
  *mechanism* the critique names is real; only its *consequence* ("only defense is binary") is
  now organism-qualified (see D).
- **The legacy 5-level sort (Path B) is a popularity vote** with `_n_agree` =
  `groupby(read_id, chrom, corrected_3prime).size()` (`corrected_consensus.py:1395-1407`,
  `:1451-1454`). Accurate as a description of the *fallback*. (Note: it now groups on
  `(read_id, chrom, corrected_3prime)`, not `(read_id, corrected_3prime)` — a paralog fix; minor.)
- **deSALT: no `-x`, no `-G`** (Claim 4a) — CONFIRMED.
- **minimap2 flag set** (Claim 4b) — CONFIRMED (+`-y`).
- **gapmm2 `-i max_intron`** (4c), **uLTRA `--ont --disable_infer`** (4d) — CONFIRMED.
- **Module 2H policy** (Claim 5): sequence-first, `_CANONICAL_HP_PRIOR=0.5`, `tier_beats_alt`,
  no candidate guard, `search_radius=5000`/`max_boundary_shift=50` — CONFIRMED unchanged.
- **`--splice-flank=no  # Disable for compatibility`** comment + the contradiction vs CLAUDE.md
  (Claim 6) — CONFIRMED.
- **Win-rate numbers (78.9/18.2/2/0.8/0.1) are single-dataset and un-committed** — still a valid
  caution (no `aligner_summary.tsv` artifact verified on the build for these exact figures).
  *However*, the specific argument "they came from Path B / legacy sort" is no longer
  supportable for a fresh build run (see A); the provenance caution should be reframed as
  "single dataset, un-committed artifact, and the selection metric on the current build is
  HP-ED Path A — so the legacy-sort lineage explanation no longer applies to new runs."

---

## Distinguishing CONFIRMED vs INFERENCE on the build

- **CONFIRMED (read directly in build source):** Claims 1, 1b, 2, 3 (both the N=0 cost and the
  anchor-gate code + its organism defaults), 4a–4e, 5, 6.
- **INFERENCE (not executed, reasoned from code):** that a *fresh* production run on the build
  *does* invoke Path A end-to-end depends on raw BAMs being present on disk at merge time
  (`split_command.py:1062-1067` globs `*.{aligner}.bam`; `single_sample.py:_sample_per_aligner_bams`).
  If those globs miss (aligner produced no BAM), that aligner falls out of HP scoring and, if
  *all* miss, `ed_rows` is empty → `use_hp_ed=False` → Path B fallback (`corrected_consensus.py:1377`).
  This is an operational edge, not the default. The win-rate *biology* claims remain UNTESTED
  (no ground-truth concordance committed) — that part of the redteam critique is unaffected by
  the build and still stands.

---

## Deliverable files to edit and exactly what to change

1. **`03_adversarial/redteam_winrates_selection.md`** — the most affected file.
   - Rewrite §(a) "Which path runs in production?" (L56-75): both call sites now pass
     `per_aligner_raw_bams` + `genome`; `use_hp_ed = bool(per_aligner_corrected_bams or
     per_aligner_raw_bams)` (`corrected_consensus.py:1262`); a lazy raw-BAM HP path exists. → On
     the build **Path A (HP-edit-distance) runs in production.** Keep Path B description but label
     it "fallback only."
   - Claims 2 and 3 in the table: change ratings from "FACTUALLY WRONG / NOT-IN-PRODUCTION" to
     "**SUPERSEDED on build — HP-ED now runs**." Claim 6's "penalty doesn't run on the
     correct-first path" stays true *for the raw-BAM `consensus.py` junction-proximity penalty*,
     but note HP-ED + overhang + anchor gate now DO run.
   - §(d) experiment 2 ("Path A vs Path B ablation ... single most important experiment"): note
     it is now an *ablation of two live paths*, not "wire a dead path on."
   - §(a) note on N=0 (L77-81): add the anchor-gate backstop + organism gating.

2. **`SPLICE_JUNCTION_PLACEMENT.md`**
   - §2 "Empirical `--junction-penalty-table` (FACT — and a CRITICAL gap)" L160-168: strike
     "ABSENT from this checkout / heuristic fallback is what actually runs." Tables are bundled
     and auto-resolved by `--Scer`.
   - §2 "CRITICAL FINDING — junction quality does NOT drive production selection" L189-214:
     retract; both call sites pass raw BAMs → `use_hp_ed=True`. Reframe item 2 (N-op cost 0) with
     the anchor-gate backstop (human=10, yeast=0).
   - §1 table footnote L122-124 ("measured under the legacy selection path"): retract for the build.
   - §5 J1/J2: J1(b) wiring is DONE on the build; J2's "currently absent" premise is wrong.

3. **`04_discovery/improve_splice_junctions.md`** and **`04_discovery/REDTEAM_proposals.md`**
   - J1/E0/C6 ("collect the already-on-disk per-aligner corrected BAMs and pass them to the
     merge so Path A runs"): mark as **already implemented** on the build (via lazy raw-BAM path).
     Remaining J1 work = the calibrated N-op open cost only.
   - J2 / C13 ("regenerate the absent penalty tables"): re-scope to "validate the bundled tables."
   - C11 ("Numba not importable / penalty tables absent here"): environment-specific to the
     master checkout; does not describe the build (tables present; Numba availability still
     environment-dependent — leave that sub-claim).

4. **Dossier `01_investigation/mapPacBio_bbmap.md`** (and any COMPARISON.md mention): correct
   `intronlen=50` → `intronlen=10`; add `maxindel=max(200000, max_intron)`; mark redteam_denovo
   B5 ("no maxindel") RESOLVED.

5. **`01_investigation/deSALT.md` / `uLTRA.md` / `gapmm2.md` / `minimap2.md`**: no change to the
   flag lists (CONFIRMED), but where they assert the penalty tables are absent or HP-ED is
   "not wired," add a build-correction pointer to this file.
