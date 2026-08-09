# MICROHOM AUDIT V5 — Auditor A — task "independent-remeasure"

**Role:** adversarial auditor (Opus-Max), READ-ONLY, own harness, own verdict. One of two independent auditors
on the same task; NO coordination. FIND THE FAULT.

**Under audit:** the POSITIONAL-DISTINCTIVENESS CLOSE for RECTIFY's motif-blind junction re-placer's
microhomology-drift guard. Specifically:
- `_positional_signal(genome, q, q_split, ne, new_je, W=28)` in `rectify/core/splice/junction_refiner.py`
  (lines 600-629): `ed(rescue, genome[ne:]) - ed(rescue, genome[new_je:])`, `rescue = q[q_split:q_split+W]`,
  hard-anchored edit distance (`_semiglobal_ed`, lines 580-597: free ref SUFFIX, NO free-prefix split).
- Wired as `drift_positional_gate` (default 0.0 = OFF = byte-identical). A drift-flagged would-be-veto is
  SPARED when `_positional_signal >= gate` (lines 915-930).
- CLAIMED (dev/DISCOVERY_LOSS_PANEL_RESULT.md, CP4/CP5): WIRED m3/cap2/gate1 → **~0.4% discovery loss /
  ~4.3% fab-residual**; ed_signal separates the delta overlap band [0.5,1.5] at **98-99% balanced accuracy**.

**My task:** INDEPENDENT discovery-loss RE-MEASUREMENT. Do NOT reuse `dev/discovery_loss_panel.py`. Build my
OWN panel, DIFFERENT construction + error model. Questions:
1. Does ~0.4% / ~4.3% at m3/cap2/gate1 hold on MY panel?
2. Is the original construction BIASED (U+U'+TAIL over-representing easy separation; longer repeats k>10;
   imperfect tails; higher error 10-15%; HP-adjacent boundaries; multi-candidate pools not just 2)?
3. At what microhomology/error regime does the ed signal's separation DEGRADE below usefulness? Report the
   frontier.
VERDICT: is "CLOSED" robust or construction-dependent?

**Env:** /Users/kevinroy/miniconda3/bin/python, repo root importable. Scratch:
/private/tmp/claude-501/-Users-kevinroy-work-rectify/798e5203-2185-4c27-87b4-f4e50c56c3ba/scratchpad/audit_v5/independent-remeasure_A

---

## CHECKPOINT 0 — orientation (code read, plan set)

**Code verified (not trusted from summary):**
- `_semiglobal_ed(query, ref)` (line 580): DP, `prev = list(range(n+1))` init (row0 = j deletions),
  `cur[0]=i` (col0 = i insertions), answer = `min(prev)` over last row = free ref SUFFIX, hard anchor at
  ref[0]. Confirmed: aligns ALL of query to a PREFIX of ref, free ref suffix, NO free query/ref prefix.
- `_positional_signal` (line 600): `rescue = q[q_split:q_split+W]` (W=28); `ref_inc = genome[ne:ne+W+6]`,
  `ref_mov = genome[new_je:new_je+W+6]`; returns `ed(rescue,ref_inc) - ed(rescue,ref_mov)`. Returns None if
  `new_je==ne` or rescue empty. So >0 ⇒ rescue matches MOVED acceptor better ⇒ genuine read evidence.
- Wiring (lines 915-930): veto fires when `best_score_cmp > incumbent_score - veto_margin`; if
  `drift_positional_gate>0` and `_positional_signal >= gate`, SPARE (moves stays True). ACCEPTOR-only (uses
  ne, new_je; donor/both-boundary fall through to margin/cap — no spare).
- `_move_microhomology` (line 511): min over MOVED boundaries of `_frac_match(shifted k-mer, k-mer at drift
  distance)`. Acceptor: `frac(seq[lo:hi], seq[hi:hi+k])`. Off-genome-edge → 0.0 (does not flag).

**Original panel's construction (the thing I must NOT reuse, and must challenge for bias):**
- genome = LPAD + exon1(60) + canonical intron GT..AG(90) + EXON2; EXON2 = U(k) + U'(=U w/ mm mismatch) +
  TAIL(random 30) + rand(40). ns=donor, ne=start of U (incumbent), je=ne+k=start of U' (cryptic).
- cry read: exon2 = genome[je:je+E2] = U' + part(TAIL). fab read: exon2 = genome[ne:ne+k+E2] = U + U' + tail.
- ONT error: 6% sub, 3% indel (indel = 1.5% del + 1.5% ins). Only 2 candidates in pool: [(ns,ne),(ns,je)].
- delta_improve = score_junction(inc) - score_junction(cry). mh = _move_microhomology.

**BIAS HYPOTHESES to test with a DIFFERENT construction:**
- H1 (TAIL-offset artifact): the discrimination signal is the RANDOM TAIL appearing at a k-bp offset between
  inc-exon2 (U+U'+TAIL at offset k) and cry-exon2 (U'+TAIL at offset 0). A random 30-40bp tail is MAXIMALLY
  distinguishing (ed jump ~k+). Real cryptics may have far less distinctive downstream context (low-complexity,
  HP, another repeat) ⇒ ed signal collapses. TEST: make the post-repeat context LOW-complexity / HP / a second
  repeat, not iid-random.
- H2 (perfect 2-candidate pool): only inc vs cry offered. Real refiner has many candidates. TEST: multi-candidate
  pool incl. decoys at je±1, je±2 that also sit in microhom.
- H3 (repeat length): k in {4,6,8,10} only. TEST k up to 20+ (longer direct repeats → longer ambiguity → ed on
  W=28 window may not reach past the repeat).
- H4 (error rate): only 6%/3%. TEST 10-15% (real ONT DRS worst-case), which floods the W=28 rescue with errors
  → ed(rescue, either ref) both large, DIFFERENCE shrinks to noise.
- H5 (imperfect/variable tail take, HP-adjacent boundaries).

**MY panel design (independent):** genome model = exon1 + intron + [direct-repeat block] + [downstream context
of controllable complexity]. Incumbent = repeat start, cryptic = repeat start + k. Post-repeat context KIND ∈
{random, homopolymer, dinuc-STR, second-repeat}. Sweep k up to 22, error ∈ {6/3, 10/5, 15/6}, imperfect repeat
(mm), multi-candidate pool. Measure guard-OFF discovery/fab, then m3/cap2/gate1 close (via the REAL refiner) +
the ed_signal separation. Compare frontier to the claim.

**STATUS:** code read, plan set, no compute yet. NEXT: write harness to scratch, smoke-test on trivial case,
then run the panel.

---

## CHECKPOINT 1 — REPRODUCTION CORNER MATCHES (harness validated)

Harness: `scratchpad/.../panel_A.py` (build), `analyze_A.py` (measure). DIFFERENT construction from the
original (repeat-block + downstream-context model, drives REAL `_positional_signal` W=28 fn + REAL refiner),
but the `mode=repro` corner (2 copies U+U', random 30bp tail, 6/3 error, 2-cand pool, motif_blind,
boundary_error_window=0, donor fixed) SHOULD recover the claim. It DOES:

```
repro COMBINED (seeds 1&2, 2160 reads):
  cry-at-risk (moved-OFF, mh>=0.5): 1080; delta med=7.0; esig med=8.0
  fab-drift guard-OFF: 2.1%; delta med=1.0; esig med=-4.0
  CAP-ALONE (m3/cap2)      disc-loss 0.6%   fab-residual 0.0%
  WIRED CLOSE (m3/cap2/g1) disc-loss 0.0%   fab-residual 4.3%   <- CLAIM: 0.4%/4.3%  ✓ MATCH
  SEP delta[0.5,1.5]: cry n=2 fab n=23 | bal-acc 98%            <- CLAIM: 98-99%     ✓ MATCH
```

**⇒ HARNESS VALIDATED.** My independent build recovers 0.0%/4.3% + 98% band-separation on the easy regime.
So any DIVERGENT frontier number below is a real construction-dependence, not a harness artifact.

**ALREADY VISIBLE (the tell):** in the easy regime cry-at-risk have delta med=7 / esig med=8 — they almost
NEVER land in the [0.5,1.5] band (only 2 of 1080). The "98% band separation" is carried by where the 23 FAB
reads land (esig med=-4), with essentially NO real cryptics in the band to lose. i.e. the easy construction
puts real cryptics FAR from the near-tie band. The frontier question = does a realistic construction pull
real cryptics INTO the band with esig that no longer separates? That is exactly the microsatellite/window
attack. NEXT: run array (microsatellite), longk, context modes.

---

## ★★★★ CHECKPOINT 2 — THE FAULT: the W=28 fixed window creates a microsatellite BLIND SPOT

Ran array/longk/ctx (mode panels) then a FOCUSED microsatellite panel (`panel_micro.py`, 2×40/cell,
2880 reads, perfect tandem arrays of period∈{3,4,6}, copies∈{2,3,4,6,8,12}). Data:
`micro_s1.csv`, `micro_s2.csv`. Result (COMBINED, WIRED close via REAL refiner):

```
                            claim (easy)      MY microsatellite panel
  WIRED m3/cap2/gate1:      0.4% / 4.3%       disc-loss 7.2% / fab-residual 9.1%   <- BOTH worse
  band[0.5,1.5] bal-acc:    98-99%            73%                                  <- separation collapses
```

**STRATIFIED BY `span` = (copies-1)*period (how far the identical tandem array extends past ne):**
```
   span   cryN   esig(cry)med   close-loss
    <=21   ...      2..4           0.0%        <- array shorter than W=28: signal clean, close holds
     28     62      0.0           69.4%        <- array reaches W: signal BLINDED -> massive loss
     30     19      0.0           94.7%
     33      2      0.0          100.0%
     42/44   .      0.0          100.0%
     66      2      0.0          100.0%
```
Stratified by `distinguishable` (first inc/cry divergence < W=28?): distinguishable=0 → esig med=0.0,
**close-loss 79.2%** (vs distinguishable=1 → 0.0%). Clean threshold at span≈W=28.

**★ THE SMOKING GUN (hand case, period=6 copies=6, span=30):** a read GENUINELY from the cryptic je carries
exon2 = `GGCTTT×5 + TAG TCGTGGGATGATC` — the 30bp array THEN a distinctive tail at read-offset 30.
Over the FULL exon2: `ed(read, genome[je:]) = 0` (perfect cryptic match), `ed(read, genome[ne:]) = 6`
(incumbent mismatch). **The read carries DECISIVE evidence for je.** But `_positional_signal` hard-codes
`rescue = q[q_split : q_split+28]` and `ref = genome[...: ...+W+6=34]`, so the window TRUNCATES at offset 28,
BEFORE the distinctive tail at offset 30 → `_positional_signal` returns **0** → 0 < gate(1) → the real
cryptic is NOT spared → VETOED. Verified live: `_positional_signal(g,q,q_split,ne,je) == 0` on the clean read.

**⇒ THIS IS A REAL, NON-AMBIGUOUS FAULT, NOT FUNDAMENTAL SINGLE-READ AMBIGUITY.** The information that
distinguishes the real cryptic from a fabrication EXISTS in the read (ed 0 vs 6 full-length); the fixed W=28
window DISCARDS it whenever the microhomology repeat span reaches ~28 bp. The close's own design doc frames
"hard-anchored ed (removing the free-k soft-clip escape)" as the fix — but it silently traded the soft-clip
escape for a FIXED-WINDOW truncation escape. Microsatellites/tandem arrays ≥28bp (biologically ubiquitous:
(CAG)n, poly-pyrimidine tracts, (AT)n) are precisely where this bites, and are exactly the regime the
original U+U'+TAIL panel (repeat span = k ≤ 10, always << 28) STRUCTURALLY EXCLUDES.

**BIAS CONFIRMED (question 2):** the original panel's construction IS biased — it caps the microhomology
repeat length at k≤10, so the identical-window regime that blinds the W=28 signal is never sampled. The
"98% balanced separation / 0.4% loss" is a property of short-repeat microhomology only.

**FRONTIER (question 3):** the ed signal's separation is useful iff the first inc/cry divergence falls
INSIDE the W=28 window (equivalently, tandem-array span < ~28 bp). At span ≳ 28 the signal → 0 and
discovery loss → 70-100%. Threshold ≈ 28 bp of contiguous microhomology, sharp.

**STATUS:** decisive fault found + hand-verified + persisted. NEXT: (a) confirm the mode=array aggregate
(disc-loss 1.1% there) only looks mild because most array cells have span<28 — i.e. verify the dilution;
(b) sanity that fab-residual worsening (9.1%) is also real (some fabs get esig>0 in the ambiguous window and
are wrongly spared); (c) then verdict.

---

## ★★★★★ CHECKPOINT 3 — FRONTIER LOCKED + END-TO-END CONFIRMED + DEFAULT-OFF context

**(a) End-to-end through the REAL refiner (period=6, copies=6, span=30, clean read):** mh=1.0,
`_positional_signal=0`, guard-OFF moves ne→je (discovery), but WIRED CLOSE → `reps: []` (held = LOST).
Profile confirms the veto path executed. NOT an overlay artifact — the shipped code vetoes the real cryptic.

**(b) Frontier is SHARP at span≈28 and ROBUST to error rate** (period=4, copies sweep, 2 error models):
```
   span     close-loss(6/3)   close-loss(12/6)     close-fab(6/3)
     8          0.0%              0.0%                 0%
    16          0.0%              2.0%                50%   <- fab CONVERSE: fabrications wrongly SPARED
    24          4.0%              6.1%                33%
    28         55.3%             61.0%                 —    <- BLIND SPOT: signal->0, massive discovery loss
    32        100.0%            100.0%                 —
```
Higher error (12/6) does NOT rescue the blind spot (advisor was right: the signal is a DIFFERENCE of two eds
vs the SAME corrupted rescue → error is common-mode & cancels; the blinding is pure window-truncation
GEOMETRY). Higher error makes the OVERALL number worse (disc-loss 13%→19%) by pulling more cryptics toward
the near-tie band, and band[0.5,1.5] balanced-acc collapses to **58-63%** (≈ coin flip) vs claimed 98%.

**(c) FAB CONVERSE is real:** at intermediate span (16) close-fab-residual = 50% — inside the partially-
ambiguous window, error-driven FABRICATIONS score `esig>0` and are WRONGLY SPARED. So the close both (i)
LOSES real discovery (span≥28) and (ii) ADMITS fabrication (intermediate span). It is not a safe monotone
guard in the microhomology regime.

**Dilution note:** the aggregate `mode=array` panel showed only 1.1% close-loss — because most of its cells
(copies 3/5/8, period 4/6/8/10) have span < 28 where the signal is fine. The fault is CONCENTRATED and
INVISIBLE to any aggregate that under-samples span≥28. The original panel (k≤10) never samples span>10 at all.

**DEFAULT-OFF / UNWIRED (verified in code):** `drift_positional_gate` (and `microhom_drift_margin`,
`drift_near_tie_cap`) default 0.0 in all 4 refiner fns (refine_read_junctions / _run_sequential /
_run_parallel / refine_bam_junctions) and `correct_command.py` passes NO drift kwargs → production is
byte-identical (guard off). So the fault is LATENT: it harms only if a user OPTS IN to the guard at the
documented operating point (m3/cap2/gate1). This bounds the blast radius but does NOT rescue the CLAIM.

---

## ⇒ VERDICT (Auditor A, independent) — FAULT FOUND; the CLOSE is CONSTRUCTION-DEPENDENT

**Q1 (does 0.4%/4.3% hold on my panel?):** On the EASY corner (matching the original construction) YES — my
independent harness reproduces 0.0%/4.3% + 98% band-separation, validating the harness. On a MICROSATELLITE
construction (tandem-array microhomology) it does NOT: disc-loss **7.2%** / fab-residual **9.1%** aggregate,
degrading to **55-100% discovery loss** for arrays whose span ≥ ~28 bp.

**Q2 (is the original panel biased?):** YES, decisively. `EXON2 = U + U' + TAIL(random)` caps the contiguous
microhomology at the repeat length k ≤ 10 and appends a MAXIMALLY-distinctive random tail immediately after —
placing the discriminating divergence at offset ≤ ~k+mm, always deep INSIDE the W=28 window. The regime that
BREAKS the signal (contiguous microhomology ≥ 28 bp, so the divergence is pushed OUTSIDE the window) is
structurally excluded. The "98% separation / 0.4% loss" is a property of short-repeat microhomology, not a
general result.

**Q3 (frontier where the ed signal degrades):** The signal is useful iff the first incumbent/cryptic exon2
divergence falls INSIDE the hard-coded W=28 window in `_positional_signal` (equivalently, contiguous
microhomology / tandem-array span < ~28 bp). At span ≳ 28 the two reference windows `genome[ne:ne+34]` and
`genome[new_je:new_je+34]` become IDENTICAL → `_positional_signal ≡ 0` → 0 < gate(1) → the real cryptic is
vetoed even though the read carries DECISIVE full-length evidence (hand case: full-read ed 0 vs 6). Threshold
≈ 28 bp of contiguous microhomology; sharp; error-rate-independent.

**ROOT CAUSE:** the close swaps the scorer's "free-k soft-clip escape" for a FIXED-WINDOW (W=28) truncation
escape. Info that distinguishes real cryptic from fabrication EXISTS in the read but is DISCARDED by the fixed
window whenever the microhomology span reaches the window length — a biologically ubiquitous regime ((CAG)n,
poly-pyrimidine tracts, (AT)n microsatellites, tandem repeats). W is not caller-tunable from the guard's
public knobs; it is hard-coded at 28 in `_positional_signal`.

**Is "CLOSED" robust or construction-dependent?** → **CONSTRUCTION-DEPENDENT.** The close holds only in the
read-distinguishable regime (divergence inside W=28), which the U+U'+TAIL panel exclusively samples, and
degrades to worse-than-useless (up to 100% discovery loss AND admits fabrication at intermediate span) exactly
where single-copy-slide microhomology is longest. The claim "the fault CLOSED" overstates a
regime-limited result.

**SEVERITY / SCOPE (fair to the code):** the guard is DEFAULT-OFF + UNWIRED, so production is unaffected today;
the fault is latent, triggered only by opt-in at the documented operating point. The right disposition is NOT
"ship blocked" (nothing ships on by default) but: the CLAIM in DISCOVERY_LOSS_PANEL_RESULT.md CP4/CP5 must be
DOWNSCOPED to "closes for contiguous microhomology < ~28 bp" and the microsatellite blind spot documented as
an OPEN limitation before any enablement. A fix (windowing to the first divergence, or extending W past the
detected repeat span, or falling through to margin/cap when both windows are identical rather than treating
sig=0 as "no evidence for the move") is required before the guard is enabled on real data with repeat content.

**HOLD** (verdict axis = the CLAIM as written): the positional close does NOT robustly close the read-blind
fault; it is construction-dependent with a sharp, realistic, hand-verified failure regime.

**Reproduce:** `panel_A.py --mode {repro,array,longk,ctx}`; `panel_micro.py --periods P --copies C
[--sub --indel]`; `analyze_A.py <csv...> --strata span,distinguishable`. All in
`scratchpad/audit_v5/independent-remeasure_A/`. Key CSVs: repro_s{1,2}, micro_s{1,2}, micro_err{06,12}.
