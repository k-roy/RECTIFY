# Discovery canonical-snap TIEBREAKER — C3 loose-end CLOSED (2026-06-30)

Closes the open question C3 deferred to the Discovery facet (see `dev/C3_DESIGN.md`
§"The `junction_score` path … is closed structurally"): **does the
`canonical_count`/`_n_annotated` TIEBREAKER in
`select.py::select_best_alignment` (the `junction_score` arbiter, `tiebreak='rectify'`)
ever actually decide a snap-vs-truth contest — and if so, does reweighting it
reduce canonical-snap FDR without harming annotated-junction recall?**

C3's structural spot-check PREDICTED the tiebreaker is unreachable: a snap is ≥1bp
off → induces a junction-proximity mismatch → `_count_junction_proximity_errors`
penalizes it on the PRIMARY `junction_score` → snap loses before the tiebreaker
fires. This probe tests that prediction empirically by driving the **actual shipped
`select_best_alignment(tiebreak='rectify')`** (pulled from `drs-validation-rebuild`)
on a minimal 2-member family per JUNCTION_DISCOVERY read.

## Method
`scripts/benchmark/discovery_tiebreak_probe.py` (mirrors `c3_junction_headroom.py`):
- corpus = `controlled.generate_corpus`; members per read = **mm2** (real
  `minimap2 -ax splice`, the snap) + **truth** (true_cigar, orthogonal placer,
  position-exact by construction). Built into BAMs, turned into `AlignmentInfo`
  via the **real `extract_alignment_info`**, then fed to the **real
  `select_best_alignment(tiebreak='rectify')`**. junction_score, ties, and the
  winner are all read off the shipped code path (not reconstructed).
- `annotated_junctions` reconstructed from truth rows with `klass==ANNOTATED`.
- **Validation gate (advisor):** on the ANNOTATED stratum, the truth member must
  yield `n_annotated>=1` — else the annotated set is mis-keyed and `_n_annotated`
  silently collapses to 0. **Gate PASSED 30/30.**
- 2 members is the construction MOST favorable to the snap winning (3'-agreement
  is degenerate → `n_annotated`/`canonical_count` are forced to decide).

## Result (reps=120; `discovery_tiebreak_probe_result.txt`, evidence `discovery_tiebreak_evidence.tsv`)
| stratum | n | mm2snap | snap&tru | tie/score | tb→snap | win=snap |
|---|---|---|---|---|---|---|
| ANNOTATED | 30 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| JUNCTION_AMB | 30 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| JUNCTION_DISCOVERY | 480 | 0.460 | 0.460 | 0.023 | **0.009** | 0.009 |
| NIC | 30 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

**C3's prediction is PARTIALLY REFUTED.** The tiebreaker DOES decide a snap-vs-truth
contest, but rarely: of 221 snapped discovery reads, **5 tie on the PRIMARY
`junction_score`** and **2 are decided toward the snap** (≈0.9% of snapped reads,
≈0.4% of discovery reads). `win=snap == tb→snap` — the snap NEVER wins on primary
score; only via the tiebreaker.

### Why the snap ties on primary score (the gap C3 missed)
On the tie reads BOTH members have `junction_proximity_errors = 0.00`. Precise
mechanism (verified against ref + scorer on the 2 harmful reads `jd_non_novel_r082`,
`_r105`):
1. The snap is an **exact-match INSERTION**, not a flanking substitution: minimap2
   finds a shifted placement whose flanking exon bases all match (`=` ops) by
   exploiting a 2bp donor repeat (`ref[200:202]==ref[400:402]`), parking a 1–3bp
   insertion at the intron boundary. E.g. r082 truth `200M200N200M` (intron
   [200,400)) vs snap `202=201N1I197=` (intron [202,403); donor +2 / acceptor +3).
2. `_count_junction_proximity_errors` inspects only the **exon side** of each
   boundary (`[istart-5,istart)` and `[iend,iend+5)`) and attributes a post-`N`
   insertion to `prev_rp = intron_end-1` — the intron's LAST base, ~1bp INSIDE the
   intron, just outside the exon-side acceptor window. So the snap's inserted bases
   are invisible to the primary score → tie.
The miss is therefore a **one-sided (exon-only) window + intron-base insertion
attribution**, NOT a ">5bp distance" effect (the snap error is ~1bp from the
junction). Had minimap2 placed the insertion one base into exon-2 it WOULD have
been counted — the gap is representation-dependent. This is the "NM≥1 vs
proximity_errors=0" residual, and it lives in `scoring.py`, not the tiebreaker.
**Rate caveat:** measured on error-free sim reads (which maximize exact-match `=`
snaps); on real ~5–10% ONT error the rate will differ — treat ~0.4–0.9% as a
mechanism existence-proof, not a real-DRS estimate.

### Which cell, and which tuple element
- `jd_non_anno` (truth non-canonical but ANNOTATED): tiebreaker correctly picks
  **truth** via `n_annotated` (1 > 0). HARMLESS — n_annotated is the lever.
- `jd_non_novel` (truth non-canonical AND novel): `n_annotated` ties at 0;
  `canonical_count` (snap=1) or — when both are 0 — pure **insertion order** picks
  the SNAP. This is the only harmful cell.

## Reweight ablation — SUBTRACTION is inert; an ADDITIVE edit-distance key HELPS
Per-read winner recomputed under variant tuples, both member orders (fitness=truth).
**Corrected after adversarial panel (iii):** the first ablation tested only
SUBTRACTIVE variants (subsets of the shipped tuple) and wrongly concluded "no
truth-favoring signal exists." It does — the snap pays its junction-shift as an
**insertion** (`mm2_ed=1` for r082, `=3` for r105 vs `truth_ed=0`), recorded by the
standard `NM` tag, which the shipped tuple OMITS.
| variant | tuple | win_exact | snap_harm (mm2-1st) | snap_harm (truth-1st) |
|---|---|---|---|---|
| V0 shipped | (agree, nann, canon) | 568 | 2 | 1 |
| V1 drop canon | (agree, nann) | 568 | 2 | 0 |
| V2 drop canon+nann | (agree,) | 565 | 5 | 0 |
| **V3 add −ED** | **(agree, nann, −ED, canon)** | **570** | **0** | **0** |
| **V4 −ED, no canon** | **(agree, nann, −ED)** | **570** | **0** | **0** |

- **SUBTRACTION is inert/harmful.** Dropping `canonical_count` (V1) leaves recall
  (568) and harm (2, mm2-first) unchanged — ties fall to insertion order (the harm
  count changes with member order → order-arbitrary). Dropping `n_annotated` (V2)
  HARMS recall (568→565); it is load-bearing (rescues `jd_non_anno`). Must be kept.
- **ADDING an alignment-edit-distance (−ED ≈ −NM) key BEFORE `canonical_count`
  (V3) HELPS, order-invariantly:** redirects BOTH harmful reads to truth, lifts
  position-exact recall 568→**570** (full on this corpus), `snap_harm=0` under BOTH
  member orders, and HARMS no other stratum (ANNOTATED/NIC/AMB recall unchanged;
  win_exact strictly increases). V4 shows `canonical_count` is then redundant for
  this residual, but keep it (it is harmless and helps determinism elsewhere).

## PRIMARY RECOMMENDATION — the surgical scoring.py attribution fix (dominates −ED)
The cleanest fix is NOT a tiebreaker change at all. In `_count_junction_proximity_errors`
a post-`N` insertion is attributed to `prev_rp = intron_end-1` (last intron base,
outside the exon-side window) — the blind spot. Change it to `intron_end` (first exon
base, in `[intron_end, intron_end+5)`):
```
# scoring.py _count_junction_proximity_errors, the N branch:
-        elif op == 3:  # N: intron skip
-            prev_rp = ref_pos + length - 1
+        elif op == 3:  # N: attribute a following insertion to the first EXON base
+            prev_rp = ref_pos + length            # intron_end, in the exon-side window
             ref_pos += length
```
Measured on this corpus (reps=120, 570 reads):
- Both harmful snaps gain `proxerr 0.00→1.00` (`202=201N1I197=`, `202=203N3I195=`) →
  junction_score 0→−1 → **snap loses on PRIMARY score → no tie → tiebreaker never
  fires.** Closes r082 AND r105.
- **0 / 570 truth members gain any proxerr** (zero recall cost on this corpus); 21
  snap members correctly gain a penalty.
- It NEVER touches `canonical_count` → no canonical-correct fence risk (the structural
  risk the −ED reorder carries; see below). This DOMINATES the −ED tiebreaker and is
  the recommended lever. (Earlier text wrongly conflated this with "widen the window,"
  which provably fails — the miss is mis-attribution, not width.)

## SECONDARY — an ADDITIVE −NM tiebreaker also helps (subtractive reweight does not)
Maps to the pre-committed option **"real residual snap-FDR; reweight helps"** — with
the correction that the helping reweight is **ADDITIVE** (insert `−edit_distance`),
NOT the subtractive "drop/demote canonical_count" the task framing suggested.
- Residual is small (2/221 snapped discovery reads ≈0.9%; ~0.4% of discovery reads),
  confined to the **non-canonical-novel** cell, and it ties on PRIMARY score because
  the snap's penalty (an insertion at the intron boundary) sits in an **intron-side
  blind spot** of `_count_junction_proximity_errors` (exon-side-only window; post-`N`
  insertion attributed to the last intron base). Verified the snap's cost IS visible
  in `NM` (1, 3) but invisible to the primary score even at window=15.
- **The scoring.py window-widen fix PROVABLY FAILS here** (panel iii, instrumented:
  `proxerr` stays 0.0 at windows 10 and 15) — the miss is structural (one-sided
  window + intron-base attribution), not a width tuning issue. So the cheap, in-scope
  lever is the tiebreaker `−NM`, not the scoring window.
- **Fence-safety of −ED is NOT fully earned by this corpus (advisor).** This corpus
  contains ONLY canonical=WRONG ties (the snap is the canonical member); it has ZERO
  canonical=CORRECT ties — the exact case `canonical_count` exists to protect and the
  case `−ED`-before-`canonical_count` would invert. Fixing r105 *requires* −ED to
  override `canonical_count` (r105 `mm2_canon=1`), so the fix and the risk are the SAME
  mechanism, and error-free sim can never produce the adversarial case. Also, on real
  ~5–10%-error ONT, whole-read `NM` is dominated by read-wide error and drowns the +1–3
  junction signal → prefer a junction-LOCAL edit count if −ED is pursued. This is why
  the scoring.py attribution fix above (local, primary-score, no canonical interaction)
  is preferred over −ED. What IS verified for −ED: the 4 pinned rectify tiebreaker tests
  (`tests/test_consensus_selection.py::test_tiebreaker_uses_annotated_junction_count`,
  `_uses_3prime_agreement`, `_default_is_rectify`, `_compass_vs_rectify_diverge`) are
  all decided by tuple elements 1–2 (3'-agreement, n_annotated); inserting `−ED` at
  position 3 leaves them invariant, and synthetic `AlignmentInfo`s default
  `edit_distance=0` → tie at `−ED` → fall through unchanged. `tests/test_gmap_fence_regression.py`
  pins a DIFFERENT module (`.score`/`n_canonical_junctions`/`_canonical_within_window`,
  the GMAP/COMPASS short-read scorer), untouched by a select.py tiebreak change.

## Prototype select.py diff (FLAGGED — not applied; validate before shipping)
Two-part change (productionizes V3):
1. `extract.py::AlignmentInfo`: add `edit_distance: int = 0`.
2. `extract.py::extract_alignment_info`: populate it, e.g.
   `edit_distance = read.get_tag('NM') if read.has_tag('NM') else _ref_edit_distance(read, genome)`
   (minimap2 runs with `--MD`, so `NM` is present on real BAMs).
3. `select.py::select_best_alignment` rectify `_tiebreak_key`:
   `return (_count_3prime_agreement(name), _n_annotated(name), -alignments[name].edit_distance, alignments[name].canonical_count)`
**Required validation before shipping** (NOT done here — thin pysam env lacks pytest;
M1-light): run `pytest tests/test_consensus_selection.py tests/test_junction_refiner.py
tests/test_validation_reads.py -m "not slow"`, the smoke gate, AND measure recall on
real multi-aligner DRS data (`wt_by4742_rep1`). The recall RISK: on real ~5–10%-error
ONT reads the truth-placing aligner carries baseline `NM` from read errors, so `−NM`
only discriminates by the snap's MARGINAL +insertion; a correct novel placement that
happens to carry higher whole-read `NM` than a competitor could lose. The 570/570 here
is the error-free-sim upper bound; the real-data recall cost is the gating measurement.

## Files
- `scripts/benchmark/discovery_tiebreak_probe.py` (+ `_result.txt`, `_evidence.tsv`)
- M1-light, LLR-free, fitness=truth. NO product code changed. Smoke gate GREEN.
- NOTE: probe pulls the canonical `rectify/core/consensus/` + the benchmark infra
  (`scripts/benchmark/`, `rectify/core/benchmark/`) into the worktree UNCOMMITTED
  (from `drs-validation-rebuild` / a317) purely to run; do not commit those.
