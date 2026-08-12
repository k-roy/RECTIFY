# HANDOFF — the Re-aligner line. Brief for a fresh agent.

**Written 2026-08-12 ~11:30 PDT by session `rectify-realigner`, closing an overnight run.**
You inherit the Re-aligner/Station line and the coordination role with three live sessions.
Detailed session state is in `HANDOFF.md` (five sections, current). This is the orientation.

---

## 1. What you own, and what you must not touch

**Yours:** the Re-aligner stations (A = overhang resolver, B = consensus triage, C = pool-gate),
the consensus module, junction-level correctness, and coordination with the DRS / cDNA / rbrowse
sessions.

**NOT yours — live sessions own these, each with its own handoff:**
- **DRS run-all** (`668-drs-arm`) — 48-sample DRS cohort, `682_drs1m`, Stations B/C sweeps
- **cDNA run-all** (`cdna-trim-fix`) — 31-sample cDNA cohort, `684_p1cdna_1M`
- **rbrowse** — the viewer; consumes our tags and schemas LIVE

Do not edit their trees, resubmit their jobs, or change a schema they consume without telling them.

🔴 **The 676 uBAM deletion is NOT done and must not be done without Kevin approving the file list.**
Data is safe on project (`/u/project/guillom/kevinroy/676_ubam_retire`, 313 G, 12 files), 320 GB of
scratch already reclaimed, `pt` tags confirmed preserved in every deposit. The count proof is
specified in `HANDOFF.md` §Resume — **md5 is the WRONG test** (27-byte `@PG` delta on 9 files; the
3 wtaa are genuine shard merges).

---

## 2. State of the tree

**`origin/master` = `0128840`**, full suite **2,312 passed / 0 failed**. Landed overnight:

| commit | what |
|---|---|
| `b3a8c35` | cDNA consensus trim fix (frame mismatch) — `XQ==0` 52.5→13.5 %, clip 65.8→3.8 nt |
| `4533de5` | `--require-aligners` + `DROPPED-ALIGNER` summary; resolver `threads` no-op now warns |
| `24c7805` | per-clip candidate ceiling + `refused_candidate_blowup` |
| `d0e3a0f` | impossible-intron guard: 10 kb bound, soft-clip at the junction, `Xi:i:<bp>` |
| `ae69e79` | write-time invariant: no alignment may run past its contig end |
| `d14cd22` | pool-gate `-o` silently overwrote versioned outputs (`with_suffix` REPLACES) |
| `cdf4bd9` | `analyze --manifest` no longer requires the unused positional |
| `0128840` | **`Xn` collision fix** — see §5 |

---

## 3. 🔴 THREE DECISIONS WAIT ON KEVIN. Do not implement without him.

All three change Station C's **scientific output**, and rbrowse consumes the schema live.

**(a) SUS / `I_GR` repetitiveness gate — my recommendation, and it supersedes (b) and (c).**
Record: Chanfreau **`planning/690`**. Measured on 110,219 Station C rows:

| class | frac SUS>32 |
|---|---:|
| real spliceosomal introns | **0.049** |
| genome background | 0.054 |
| short/mid novel `admit_candidate` | 0.049 / 0.019 |
| **long (>10 kb) novel `admit_candidate`** | **0.875** |

18× separation, and what Station C *rightly* admits is indistinguishable from real introns.
**Mechanistic close:** `PoolGateConfig.overhang_cap = 60`, but 72 % of long-admitted termini need
>64 bases to be unique — **the scoring window is smaller than the uniqueness requirement**, so `q`
measures information that cannot disambiguate. Proposed rule, parameter-free:
*an overhang cannot justify a placement at a locus whose SUS exceeds the overhang length used to
score it.*

**(b) Length pre-gate on the verdict.** Works, needs an organism constant. Nuclear max annotated
intron is **1,002 bp**; all long ones are mitochondrial (self-splicing, different machinery).
Measured cost: 99 rows, 294 reads, nothing else changes.

**(c) Consult `flagged` on BOTH verdict branches** + fix the docstring contradiction that caused it.
`station_c.py:383-393` reads `flagged` only on the non-canonical branch, so a 111 kb LTR-to-LTR
artifact with `repeat_flag=long_terminal_repeat` AND `selfhom_flag=1` was **admitted with both flags
lit**. And `station_c.py` §Repeat-context flags promises unconditional demotion while its own verdict
table restricts it to non-canonical — the code follows the table, and that contradiction is how the
hole survived review.

---

## 4. Open engineering, ranked

1. **2–10 kb residual band.** The 10 kb guard misses it (~15 % of the >2 kb population). **Cannot**
   be closed by lowering the bound — minimap2's clean max N-op is 5,064 bp. Needs a
   minimap2-disagreement discriminator.
2. **`--junction-pool-max-intron-len` is unreachable in production** — exposed only on
   `rectify consensus`, absent from `align`/`run-all`. "The guard did not fail, it was never
   reachable."
3. **D/N motif decision** — `dev/DN_THRESHOLD_ARTEFACT_20260812.md`. Both placements costed;
   reporting-side is *structurally* incapable of removing an annotated intron (Station C skips them
   before writing rows).
4. **`correct` does not propagate `Xa`/`Xc`/`Xn`/`Xt`** into `corrected_3ends.tsv`. If you add it:
   **SEPARATE columns, never a pre-computed concordance class** — the collapse is a viewer
   presentation decision and must stay reversible (rbrowse's constraint, and they are right).
5. 674 needs `h_data` raised (OOM'd at 22.469 G against a 24 G ceiling); 668b exit-1 root cause
   still unconfirmed.

---

## 5. The concordance taxonomy — AGREED with rbrowse, adopt it, do not diverge

rectify emits **no** three-class concordance. `Xa`=selected aligner, `Xc`=**selection confidence**
(a different axis — do NOT map it to concordance), `Xn`=`n_aligners_agree` (a COUNT, not a class),
`Xt`=tied aligners. Agreed definition, to be matched pipeline-side:

    class 0  identical  CIGARs byte-identical across every aligner in Xt
    class 1  agree      5' ends equal (PRE-rescue) AND 3' ends within the same PER-READ
                        ambiguity window AND junction MULTISETS equal under
                        same_junction(max_shift=60); CIGARs may differ
    class 2  disagree   anything else
    class 3  unknown    fewer than two arms — no comparison possible; draw UNKNOWN

**Three traps embedded there, each of which would have produced a large structured error:**
- Use `overhang_informativeness.same_junction` / `canonicalize_junction`, **never coordinate
  equality** — ambiguity-shifted junctions cluster at repeat-adjacent sites, so naive comparison
  inflates "disagreement" *non-uniformly by locus*.
- Compare 5′ ends **pre-rescue**. `five_prime_rescued` marks reads the consensus stage moved —
  261k in `wt_rep1` alone. Comparing rescued-vs-raw scores the pipeline working as disagreement.
- 3′ ends use the **per-read** `ambiguity_min`/`ambiguity_max`/`ambiguity_range`, not a constant.
- Class 3 is **~20 % of reads** and concentrates where other aligners failed. Not an edge case.

Status: **NOT BUILT.** rbrowse is holding pending Kevin; they will send the measured class
distribution before shipping. `682_drs1m` is the **only** source of concordance data that exists —
tell rbrowse before it moves or is re-run.

---

## 6. Two hazards worth carrying, both learned the hard way

**A numeric tag written by two paths is invisible when the ranges overlap.** My impossible-intron
guard wrote `Xn`, which already meant `n_aligners_agree`, and ran *after* it. rbrowse put it best:
**`Xn:i:3` is simultaneously a plausible agreement count and a plausible offending length, so a
contaminated bundle draws confident concordance classes with no symptom.** Fixed in `0128840`
(renamed `Xi`), with a regression test. `682_drs1m` is **verified clean** (only 1/2/3 across 600k
reads) despite being under the `ea4401e` umbrella — its consensus predates the guard. **Do not
order a needless re-tag of it.**

**A guard that references a name not in scope fails only when the guard fires** — so the NEGATIVE
test catches it, never the happy path. `os` in `station_c.py` and `logger` in `analyze_command.py`
both raised `NameError` instead of working. Related: my resolver ceiling test passed **vacuously**
until an anti-vacuity assert was added, and again when a tandem-repeat fixture met master's period
gate. **Write the negative case, and assert the guard's precondition was actually reached.**

**And: tRNA introns contaminate ANY annotation-derived control in this genome.** They broke a
positive control at 0.758 against a 0.049 background (`planning/690`) and skewed the intron-length
analysis (58 of 60 sub-40 bp annotated introns are tRNA). Split them out by default.

---

## 7. How to work with the peers

They are good, and three of my four hypotheses on the junction question were refuted by their
measurements — correctly. The pattern that produced every real finding overnight: **chase the
absurd number into the source rather than explaining it away.** rbrowse checked a distribution
instead of spot-checking a read; the cDNA session went to measure the cost of its own proposal and
found an n=5 that shouldn't exist; DRS caught its own false alarm as a tree/runtime skew rather
than reporting it as a regression.

Verify their claims yourself before acting — they expect it and it has caught real errors in both
directions, including mine.
