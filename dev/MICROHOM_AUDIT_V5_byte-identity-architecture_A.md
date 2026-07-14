# MICROHOM AUDIT V5 — Task "byte-identity-architecture" — Auditor A

**Role:** Adversarial auditor A (Opus-Max), independent. Verifying RECTIFY's positional-distinctiveness
CLOSE for the motif-blind junction re-placer's microhomology-drift guard.

**Scope of THIS task:**
(a) BYTE-IDENTITY: guard+cap+positional stack byte-identical at default (all four kwargs = 0.0).
    Independent BAM diff vs pre-close parent (sequential raw-order SHA256 + parallel position-sorted).
    pytest -m "not slow" green. Confirm `_positional_signal`/`_semiglobal_ed` never called at default.
(b) ARCHITECTURE: is a 4th stacked guard parameter the right design vs fixing `_score_junction`'s
    free-k soft-clip so delta_improve itself discriminates? Is a SECOND alignment (edit distance) of
    the same rescue redundant? Does the gate interact cleanly with the cap (spare-then-veto ordering)?
    Complexity/maintainability vs benefit on a default-off feature.

**Working dir:** /Users/kevinroy/work/rectify/.claude/worktrees/agent-a25a2c1e784ad37dc (READ-ONLY on repo)
**Python:** /Users/kevinroy/miniconda3/bin/python
**Parent commit (pre-close):** d1fd08d is the last commit BEFORE cd3de46 (the CLOSE). Actually need to
  verify: cd3de46 = "positional-distinctiveness gate closes...". Its parent chain determines the byte-diff base.

---

## CHECKPOINT LOG

### CP0 — record created
- Scratch: /private/tmp/claude-501/-Users-kevinroy-work-rectify/798e5203-2185-4c27-87b4-f4e50c56c3ba/scratchpad/audit_v5/byte-identity-architecture_A
- Commits confirmed: a97ff6d (donor invest) <- cd3de46 (CLOSE) <- 4440c2c (panel seed2) <- c0b5ae5 <- d1fd08d
- BYTE-DIFF PARENT for the CLOSE = **4440c2c** (`git log --format='%H %P'` of cd3de46). junction_refiner.py
  was touched ONLY by cd3de46 in that window (c0b5ae5/4440c2c are docs/tests only). So the code-relevant
  pre-close parent state of junction_refiner.py == the file at 4440c2c == the file at d1fd08d.
- Next: read junction_refiner.py CLOSE code + wiring.

### CP1 — static byte-identity analysis (STRUCTURAL PROOF)
The CLOSE diff to junction_refiner.py is PURELY ADDITIVE + guarded:
1. Adds `_semiglobal_ed` (l.580) and `_positional_signal` (l.600) — new top-level funcs, no existing code
   references them at import time.
2. Threads `drift_positional_gate: float = 0.0` as a DEFAULTED kwarg through 4 signatures
   (refine_read_junctions l.655, _run_sequential l.1366, _run_parallel l.1469→worker-state l.1525,
   refine_bam_junctions l.1701) + forwards at 3 internal call sites. All defaults 0.0.
3. The ONLY behavioural change is the veto block (l.915-930). Pre-close (parent) that block was:
       if best_score_cmp > incumbent_score - veto_margin:
           moves = False; profile.inc('move_margin_vetoes')
   Post-close it is wrapped:
       if best_score_cmp > incumbent_score - veto_margin:
           spared = False
           if drift_positional_gate > 0.0:        # <-- GATE. default 0.0 -> block SKIPPED
               psig = _positional_signal(...)
               if psig is not None and psig >= drift_positional_gate: spared = True; inc(...)
           if not spared:                          # default: spared always False -> ALWAYS taken
               moves = False; profile.inc('move_margin_vetoes')
   At gate==0.0 the inner block is dead, `spared` is always False, so the net effect is IDENTICAL to the
   parent (`moves=False` + same profile counter). **Structurally byte-identical at default.**

GREP PROOF (call-site census):
- `_positional_signal` production call sites: ONLY junction_refiner.py:922 (inside `if gate>0.0`). Others
  are tests + its own docstring.
- `_semiglobal_ed` production call sites: ONLY junction_refiner.py:629 (inside _positional_signal). Others
  are tests + dev/discovery_loss_panel.py (a standalone panel script, NOT imported by rectify).
- `drift_positional_gate` set to non-zero ANYWHERE in rectify/ : NONE. Only tests + dev panel pass >0.
- CLI wiring (correct_command.py:746 refine_bam_junctions call): passes NONE of hp_drift_margin /
  microhom_drift_margin / drift_near_tie_cap / drift_positional_gate. No config/argparse key exists for
  any of them (grep over rectify/ excluding junction_refiner.py = EMPTY). => guard genuinely UNWIRED.
- 3rd refine_read_junctions call (l.1961, an analysis/measurement fn) passes NO drift kwargs either.

Next: (a) monkeypatch trip-wire to PROVE _positional_signal/_semiglobal_ed never called at default on a
real refiner run; (b) pytest -m "not slow"; (c) real BAM diff parent-vs-HEAD.

### CP2 — TRIP-WIRE HARNESS GREEN (never-called-at-default PROVEN empirically)
Harness: scratch/tripwire_harness.py (monkeypatches JR._positional_signal + JR._semiglobal_ed to RAISE).
Scenario = the repo's own veto-band microhomology case, reconstructed independently (read from cryptic JE
placed at incumbent NE; unit ACGTAC + 1-mismatch copy ACGTAG => microhomology 5/6; distinctive exon2 tail).
ALL 6 checks PASS (scratch/tripwire_results.json):
- pure-default run (NO drift kwargs, == correct_command): new funcs called 0x. [veto block not even entered]
- WORST CASE drift-on/gate-off (microhom_drift_margin=8.0, drift_positional_gate=0.0): veto block IS
  entered (eff_margin>0), wrapper code executes, yet new funcs called 0x. acc stayed NE-region (veto fired).
- drift+cap on / gate off (microhom_drift_margin=8, drift_near_tie_cap=2, gate=0): new funcs called 0x.
- POSITIVE CONTROL: gate=1.0 -> _positional_signal called 1x, _semiglobal_ed 2x (2 EDs per signal). Proves
  the trip-wire would FIRE if the code path were live — the 0-counts above are real, not a dead harness.
- OUTPUT MATRIX (50 cells = 5 geometries x from_je{T,F} x 5 drift-setting combos): gate=0.0 identical to
  OMITTING the kwarg in every cell. mismatches=[].
- parent-equivalence (analytic): HEAD veto predicate at gate=0.0 == parent predicate over 1331 (best,inc,
  margin) triples. [NOTE: this is a re-implementation sanity check, NOT the real BAM diff — that's CP3.]

Next: CP3 real BAM diff (refine_bam_junctions default kwargs) parent@4440c2c vs HEAD, sequential raw-order
SHA256 + parallel position-sorted SHA256; CP4 pytest; CP5 architecture critique.

### CP3 — REAL BAM DIFF GREEN (parent@4440c2c vs HEAD, sequential + parallel — IDENTICAL)
**IMPORTANT PROVENANCE FINDING first:** the audited state is the WORKING TREE, not a clean commit.
- `hp_penalty.py` has uncommitted WIP (+91 lines) that adds `_precompute_refcol_ins_costs`. That symbol
  is defined in NO commit (cd3de46/a97ff6d(HEAD)/4440c2c all grep `def ...`=0) but IS imported by BOTH
  committed AND working-tree `junction_scoring.py` (count=3). => **HEAD-clean is NOT importable**; the
  code only runs with the hp_penalty WIP. (Not a CLOSE concern, but a repo-hygiene note: HEAD as committed
  can't `import rectify.core.splice.junction_refiner`.) Confirmed the ONLY code WIP is hp_penalty.py;
  junction_refiner.py is CLEAN vs HEAD (empty `git diff HEAD` -> the CLOSE code == committed cd3de46).
- Correct byte-diff BASE therefore = {live WIP hp_penalty} + {pre-close junction_refiner @4440c2c}. Built
  in scratch/parent_wt (git worktree @4440c2c, overlaid live hp_penalty.py, cleared stale pycache). Verified
  isolation: hp_penalty.py IDENTICAL between roots; junction_scoring.py IDENTICAL (0 diff lines, 4440c2c==HEAD);
  junction_refiner.py differs ONLY by the CLOSE (positional/semiglobal/spared/gate lines).

Harness: scratch/bam_diff_harness.py (spawn-safe; forces mp start_method='fork' because the repo's parallel
path is Linux-fork-only — module-global _WORKER_POOL_STATE isn't inherited under macOS 'spawn' -> KeyError
'header'; fork matches production Linux semantics). Input = 12-read synthetic BAM hitting the veto/drift path
(6 from-cryptic-JE + 6 from-incumbent-NE, 6 geometries). refine_bam_junctions at DEFAULT kwargs, sort_and_index
=False, digests via `samtools view`.

RESULT (digest_parent.json vs digest_head.json — bamdiff_compare.json):
  seq_raworder_sha256  PARENT == HEAD : e81d03a8478841187e51ce710dc47f307fdd7b6abbb1255700b87c5e96c602d9
  seq_possorted_sha256 PARENT == HEAD : e81d03a8...  (same)
  par_possorted_sha256 PARENT == HEAD : e81d03a8...  (same; n_workers=4)
  bytes 1902 both; stats IDENTICAL: total 12, n_op_reads 12, refined 6, unchanged 6, errors 0.
=> **Byte-identical at default across sequential (raw-order) AND parallel (position-sorted).**

**CORRECTION (important, self-caught via advisor):** at PURE DEFAULT the changed veto block is BYPASSED
ENTIRELY — line 915 gates the whole block on `eff_margin > 0.0`, and at default eff_margin=hold_margin=0.0.
Profiled proof: default run has microhom_drift_flagged=0, near_tie_cap_applied=0, move_margin_vetoes=0. So
the 6 refinements come from the plain `if moves:` emit (l.933), NOT the veto/spare logic. The default BAM
diff is a valid identity check but proves identity where the CLOSE's changed code is DEAD. => needs CP3b
(force the block live). [Earlier draft of this CP wrongly said "the changed veto block executes" at default
— that was FALSE; corrected here.]

### CP3b — VETO-FIRING BAM DIFF (the CLOSE's changed lines EXECUTED, parent vs HEAD IDENTICAL)
To exercise the EXACT restructured inner block on real code, pass `microhom_drift_margin=8.0` (NO cap) to
BOTH roots (parent supports it — predates the CLOSE via 3a716aa) => veto_margin=8 => the 6 near-tie moves
are VETOED. Profiled (HEAD, gate=0.0): microhom_drift_flagged=6, move_margin_vetoes=6, replacements_emitted=0,
positional_gate_spared=0. This drives:
  - parent: `if best>inc-margin: moves=False; inc('move_margin_vetoes')`
  - HEAD  : `if best>inc-margin: spared=False; [gate branch dead]; if not spared: moves=False; inc(...)`
Harness scratch/bam_diff_veto_harness.py. RESULT (veto_parent.json vs veto_head.json — bamdiff_veto_compare.json):
  refined/unchanged: parent 0/12, head 0/12 (all moves vetoed — veto FIRED).
  seq_raworder_sha256  PARENT == HEAD : 033c0976637fb021e92eac1f67cee5900f6f2b578a297a763ec4676610460a1a
  seq_possorted_sha256 PARENT == HEAD : 033c0976...  (same)
  par_possorted_sha256 PARENT == HEAD : 033c0976...  (same; n_workers=4)
  (digest differs from the default e81d03a8... precisely because the veto suppressed the moves — and is
   IDENTICAL between parent and HEAD.)
=> **The restructured `spared`/`if not spared` veto block is byte-identical to the parent's plain `moves=
   False` veto, ON LIVE CODE, IN THE REGIME WHERE THE CHANGED LINES RUN.** Combined with CP3 (default, block
   dead) this covers both regimes: block-bypassed AND block-entered-and-firing. The re-implementation in the
   trip-wire harness (part 3) is now SUPERSEDED by real diffs in both regimes.

Next: CP4 pytest -m "not slow"; CP5 architecture critique.

### CP4 — PYTEST
- Guard-relevant SUBSET (all 31 microhom-drift-guard tests + splice_junction + splice_summary +
  terminal_exon_refiner): **191 passed, 4 skipped, 1 xfailed, 9.69s**. GREEN.
- Full `pytest -m "not slow"` running to file (bo04nsivn -> pytest_A_done.txt); result appended at CP4b.
  [NOTE: first attempt reported exit 0 but stdout was progress-truncated by the bg pipe at
  test_validation_reads.py — that file spawns real `rectify.cli correct` subprocesses (slow-but-not-slow-
  marked), not a hang. Re-running cleanly.]
- IMPORTANT ENV NOTE for reproducibility: the suite/harness needs OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
  on macOS + mp start_method fork for the parallel refiner path; on Linux (production) this is moot.

### CP4b — FULL `pytest -m "not slow"` RESULT (authoritative RC): 1670 passed, 1 ERROR (CLOSE-INDEPENDENT)
Authoritative run bo04nsivn (`pytest > file 2>&1; echo RC=$?` — NO pipe, real exit code): pytest_A_done.txt
tail = **`1670 passed, 39 skipped, 4 deselected, 1 xfailed, 15 warnings, 1 error in 879.52s`; RC=1**.
[Correcting my earlier note: the FIRST run's "exit 0" was `pytest | tail -30` with no pipefail => that was
tail's rc, NOT pytest's. Discarded as non-evidence, per adversarial hygiene.]

The 1 error is **NOT a regression and NOT caused by the CLOSE** — proven:
- Error is at FIXTURE SETUP of `tests/test_restore_polya_from_parquet.py::test_restore_cat3_plus_2`:
  `AssertionError: Validation-bundle trim metadata TSV not found: .../scripts/validation_data/
  rebuild_2026_05/trimmed/validation_reads_polya_trim_metadata.tsv`. A MISSING DATA FILE in this worktree
  (the dir exists but `trimmed/` is empty; the path is UNTRACKED in git). Environmental provisioning gap.
- That test references junction_refiner / positional / drift / microhom / refine **0 times** (grep) — it is
  a polyA-restoration test, wholly unrelated to the audited guard.
- **The IDENTICAL error occurs on the PARENT (pre-close) worktree** (ran that single test against parent_wt:
  `1 error`, same AssertionError at line 117). => pre-existing, CLOSE-independent, would fail the same way
  with or without the CLOSE.
CONCLUSION: the audited change causes NO test regression. All guard-relevant tests pass; the sole error is a
missing-fixture-file present identically pre/post-close. "pytest green modulo one pre-existing missing-data
error unrelated to the change." Does NOT affect the byte-identity verdict (carried by the BAM diffs).

### CP5 — ARCHITECTURE CRITIQUE (the substance; advisory — does NOT block byte-identity)

**FINDING A5.1 — the CLOSE re-implements, with one deliberate change, an alignment the scorer ALREADY
computes. The redundancy is real; the fix arguably belongs IN the scorer.**
- `_score_junction` (junction_scoring.py:983) scores a candidate by `min_k (t1(k)+t2(k))` where
  t1(k)=`_score_hp_anchored(rescue[k:], genome[intron_end:intron_end+BUF])` (forward, left-anchored at
  the acceptor) and t2(k)=reverse-anchored prefix score. rescue=`query[q_split:q_split+30]`.
- `_positional_signal` (junction_refiner.py:600) computes `_semiglobal_ed(rescue, genome[ne:ne+34]) -
  _semiglobal_ed(rescue, genome[new_je:new_je+34])`, rescue=`q[q_split:q_split+28]`. `_semiglobal_ed`
  is a plain Levenshtein DP, left-anchored at ref[0], free ref suffix.
- => The CLOSE re-runs the SAME forward exon2-body alignment t1 already does, for BOTH the incumbent (ne)
  and the moved (new_je) acceptor, differing in exactly TWO ways: (1) NO free-k prefix split (hard-anchored
  at ref[0]); (2) plain edit distance vs HP-aware penalty. **Difference (1) is the entire load-bearing
  mechanism** — the commit itself says so: the free-k split "removes the scorer's soft-clip escape that
  hid the discriminating exon2 bases from delta_improve." Difference (2) (dropping HP-awareness) is an
  unremarked SIDE-EFFECT that makes the two signals live on different penalty models.
- ARCHITECTURAL CONSEQUENCE: the scorer's free-k loop is a SOFT-CLIP ESCAPE. For the WRONG (incumbent)
  junction the scorer can choose a k that aligns rescue[k:] cleanly to the exon2 body and pays only t2 on
  the ignored prefix — so `score(incumbent)` ~ `score(moved)`, compressing `delta_improve =
  incumbent_score - best_score_cmp` into the near-tie band where drift and real cryptics overlap (exactly
  the panel's "delta overlap band"). The CLOSE recovers discrimination by re-aligning HARD-anchored. But
  that means the RIGHT place for the fix is the scorer: a delta_improve computed from (or additionally
  gated by) a k=0 / anchored term would discriminate ON ITS OWN AXIS, making BOTH the 4th parameter AND
  the second alignment unnecessary. The current design bolts a corrective 2nd metric onto the OUTSIDE of a
  scorer whose known blind spot it leaves in place. => the 2nd O(m·n) alignment IS redundant with what the
  scorer computes, modulo the anchoring choice the scorer could itself expose (e.g. return (min_k score,
  k0_anchored score) or a "distinctiveness" term). This is the strongest architecture challenge.

**FINDING A5.2 — 4th stacked guard PARAMETER: growing a config surface to paper over a scorer deficiency.**
The move-gate now carries FOUR tunables on the SAME score axis: hold_margin (blunt), hp_drift_margin
(HP-context), microhom_drift_margin (microhomology-context), drift_near_tie_cap (read-evidence ceiling),
PLUS drift_positional_gate (read-evidence spare). The commit's own honest framing (l.899-905) admits the
cap "BOUNDS but does NOT add a new discriminating signal" — i.e. three of these are read-BLIND genomic-
context margins whose only effect is to make a move HARDER, and the 4th/5th are read-evidence patches to
undo the over-veto the first three cause. That is a control-theory smell: each layer corrects the previous
layer's imprecision rather than fixing the source. A scorer that discriminated real-vs-fab on delta_improve
would collapse this stack. Maintainability cost: 5 interacting thresholds, each threaded through 4 function
signatures + a worker-state dict, all needing joint tuning + the COMPASS gate — for a DEFAULT-OFF feature.

**FINDING A5.3 — gate/cap ordering is CLEAN but only because they don't co-fire at default (they can't
even be sanity-checked live at default).** At default (all 0.0) the cap and gate DO NOT INTERACT — both
branches are dead. When ENABLED the ordering is: compute veto_margin (cap-adjusted) -> if best>inc-margin,
FIRST try to SPARE via positional gate, THEN veto. This "spare-then-veto" is logically coherent: the cap
sets HOW WIDE the veto band is (by capping eff_margin), the gate decides WHETHER a read in that band is
spared. They compose without contradiction (cap can only shrink the band; gate can only rescue within it).
NO BUG found in the ordering. BUT: because both are off by default, the interaction is UNEXERCISED in any
default run — the only coverage is synthetic unit tests + the dev panel. On real data (COMPASS) the joint
behavior is unvalidated. (Verified empirically: my trip-wire's drift+cap+gate-off cell and the positive-
control gate-on cell both return acc=176 — the cap alone and the gate+cap agree there, so no ordering
inversion in that case.)

**FINDING A5.4 — window mismatch scorer(30) vs signal(28), and HP-awareness dropped, silently.** The
scorer uses rescue length 30 and ref buffer max(L+20,50); the signal uses 28 and W+6=34. Minor, but it
means the two metrics are NOT measuring the same span, and the signal is NOT HP-aware while the scorer is
(when a penalty_table is loaded — the DRS default). So on a homopolymer-adjacent acceptor the signal's
plain edit distance and the scorer's HP-penalized score can DISAGREE about which junction the read favors.
Not a byte-identity issue (default off), but a correctness concern for the ENABLED path that the docs do
not flag: the "positional signal" is a DIFFERENT, cruder scorer, not a hard-anchored view of the SAME one.

**FINDING A5.5 — the SHARPEST form of the redundancy: `_semiglobal_ed` is a FLAT-COST (HP-BLIND) CLONE of
a primitive the scorer ALREADY HAS (`_score_hp_anchored`).** `_score_hp_anchored` (junction_scoring.py:766)
is documented verbatim as "Left-anchored ... left end of ref fixed (no free prefix) and the right suffix of
ref free" — that is EXACTLY `_semiglobal_ed`'s contract (hard-anchored at ref[0], free ref suffix), except
`_score_hp_anchored` is HP-AWARE (uses the penalty_table). `_score_junction` only fails to discriminate
because it WRAPS this primitive in the `min_k (t1(k)+t2(k))` free-k loop (the query-side soft-clip escape).
=> The CLOSE did not need a new alignment at all: `_score_hp_anchored(rescue, genome[ne:...])` and
`_score_hp_anchored(rescue, genome[new_je:...])` at k=0 (no free-k loop) yield the identical, HP-AWARE
distinctiveness signal reusing the SAME penalty model — with ZERO new DP code. Instead the CLOSE wrote a
whole new flat-cost Levenshtein DP (`_semiglobal_ed`) that is both REDUNDANT with an existing primitive AND
CRUDER than it (HP-blind, so it disagrees with the scorer on homopolymer-adjacent acceptors — the DRS
default). This is the concrete, minimal-diff form of "the fix belongs in the scorer": expose the k=0 anchored
score from `_score_junction` (or a distinctiveness = anchored - free_k term) and delete both the new function
and the parameter.

**Net architecture verdict:** For a DEFAULT-OFF, UNWIRED, synthetic-only-validated feature, the byte-
identity safety is sound (proven above). But the design — a 4th/5th stacked margin + a second, cruder,
partially-redundant alignment bolted outside a scorer whose free-k soft-clip blind spot is left intact —
is a BAND-AID, not the deep fix the commit claims ("the deep fix the audit prescribed"). The deep fix is
to make delta_improve itself discriminate (expose/penalize the anchored term in _score_junction). The
"fault CLOSED" claim is OVERCLAIMED: it rests entirely on synthetic panels (seeds 1&2), the 98-99%
balanced-accuracy separation and the ~0.4%/4.3% operating point are NOT verified on real data, and the
task's own stated enablement gate (COMPASS) is untested. Recommend: do not present as "closed"; treat as
a synthetic-validated PROTOTYPE pending (i) a scorer-level alternative evaluation and (ii) COMPASS.

---

## FINAL VERDICT — Auditor A, task "byte-identity-architecture"

**(a) BYTE-IDENTITY: CONFIRMED CLEAR.** The guard+cap+positional stack is byte-identical at the default
(all four kwargs 0.0). Evidence, independent of the repo's own tests:
- Static: the CLOSE diff to junction_refiner.py is purely additive + gated; the only behavioural change is
  wrapped in `if drift_positional_gate > 0.0:` (dead at 0.0) and `if not spared:` (always taken at 0.0).
- Trip-wire (6/6 checks): `_positional_signal` and `_semiglobal_ed` are NEVER called at default — even with
  the drift margin + cap engaged (veto block entered) — and ARE called when gate>0 (positive control fires).
- Call-site census: those functions' ONLY production call site is inside the gate branch; NO code path in
  rectify/ passes a non-zero gate; the CLI (correct_command.py) passes NONE of the four drift kwargs and no
  config/argparse key exists for any of them -> the guard is genuinely UNWIRED.
- Real BAM diff, parent (pre-close junction_refiner @4440c2c + live WIP hp_penalty) vs HEAD, DEFAULT kwargs:
  seq raw-order + seq pos-sorted + parallel(4) pos-sorted ALL identical (e81d03a8...). [block bypassed here]
- Real BAM diff, VETO-FIRING (microhom_drift_margin=8, no cap) — the changed inner block EXECUTES (veto fires,
  refined 0/12) — parent vs HEAD ALL identical (033c0976...), sequential + parallel. [block live here]
- pytest -m "not slow": **1670 passed, 39 skipped, 4 deselected, 1 xfailed, 1 error** (RC=1). The single
  error is a PRE-EXISTING missing-fixture-file (test_restore_polya_from_parquet::test_restore_cat3_plus_2 —
  a polyA test referencing the guard 0 times) that fails IDENTICALLY on the pre-close parent. It is NOT a
  regression and NOT caused by the CLOSE (CP4b proof). All guard-relevant tests pass (subset: 191 passed).
  => "green modulo one pre-existing, CLOSE-independent, missing-data error." The task's "pytest green"
  requirement is met for the audited change (no change-induced failure); flagged honestly rather than
  claimed as a clean 0-error green.

**(b) ARCHITECTURE: FAULT FOUND (design, advisory — does not block the byte-identity verdict).**
Strongest challenge (unrebuttable, from the commit's OWN rationale): the CLOSE's mechanism is "removing the
scorer's free-k soft-clip escape" — a property of `_score_junction`. So the discrimination BELONGS IN THE
SCORER: `_positional_signal` re-runs the same forward exon2-body alignment `_score_junction`'s t1 already
computes, for BOTH junctions, differing ONLY by hard-anchoring (no free-k) + dropping HP-awareness. The
correct fix is to have `_score_junction` return the anchored/distinctiveness term ALONGSIDE its free-k score
(the free-k escape IS load-bearing for terminal-noise tolerance, so it must stay), making delta_improve
itself discriminate — which renders BOTH the 4th/5th stacked parameter AND the second O(m·n) alignment
unnecessary. Second challenge: "the fault CLOSED" is OVERCLAIMED — synthetic-only (seeds 1&2), the 98-99%
separation + ~0.4%/4.3% operating point are unverified on real data, the feature is unwired, and the task's
own stated enablement gate (COMPASS) is untested. No BUG in the gate/cap ordering (spare-then-veto composes
cleanly; verified they don't co-fire at default). Minor: window mismatch (scorer 30 vs signal 28) + HP-
awareness silently dropped make the "positional signal" a DIFFERENT, cruder scorer, not a hard-anchored view
of the same one — a correctness concern for the ENABLED path the docs don't flag.

**VERDICT = CLEAR** on the audited safety property (byte-identity at default): the change is default-off,
byte-identical, unwired, and inert — there is nothing to block; you cannot HOLD an inert change. **fault_found
= true** on the architecture: it is a band-aid beside a known scorer deficiency, and the "closed" framing
overclaims a synthetic-only, unwired result. Recommendation: safe to keep as default-off; do NOT enable or
call the fault "closed" until (i) a scorer-level alternative (anchored delta_improve term) is evaluated and
(ii) COMPASS real-data validation runs.

## FILES (this audit)
- Durable record: dev/MICROHOM_AUDIT_V5_byte-identity-architecture_A.md (this file)
- Scratch: /private/tmp/claude-501/.../audit_v5/byte-identity-architecture_A/
  - tripwire_harness.py + tripwire_results.json (never-called-at-default proof)
  - bam_diff_harness.py + build_input_bam.py + digest_{parent,head}.json + bamdiff_compare.json (default diff)
  - bam_diff_veto_harness.py + veto_{parent,head}.json + bamdiff_veto_compare.json (veto-firing diff)
  - parent_wt/ (git worktree @4440c2c + overlaid live hp_penalty.py — the runnable pre-close base)
  - pytest_A_done.txt (full suite), pytest_full.txt (subset)
