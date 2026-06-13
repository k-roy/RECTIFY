# TODO (fresh agent): stop the walkback lone-A regression structurally

**Status:** scoped, not started · **Drafted:** 2026-06-13 · **Owner:** unassigned
**Read first:** `AGENT_FIXES.md` entries [2026-06-04] + [2026-06-12]; this spec is self-contained otherwise.

---

## 0. TL;DR for a context-less agent

The DRS 3′-end walkback has a bug that has been "definitively fixed" 2–3 times and keeps coming
back: some reads end up with `corrected_3prime` parked **on a genomic stop-base (A)** instead of being
walked back to the first non-A (the true cleavage anchor). The clean core algorithm is correct; the
bug lives in the **guarded** production variant, where accreted early-exit/skip guards suppress
legitimate walkback for edge cases. Each past fix patched one guard → a different guard reopens it.

**Your job:** make the recurrence *structurally impossible* (not another point patch) + add a CI
regression metric so it can never silently come back. Do **not** just delete the guards — they reject
real artifacts. The fix is to separate "skip the walk" from "restrict where the walk may anchor."

---

## 1. The two functions (verified line numbers, `rectify/core/correct/walkback.py`, HEAD 2026-06-13)

- **`walkback_3prime`** (lines 43–186) — the **clean core**. Scans from the 3′ end inward, sets
  `corrected = rp` at the first position where `read_base == ref_base and read_base != stop_base`
  (line 181–182). **Always anchors on a non-stop base. Never buggy.** Not used as the DRS production path.
- **`walkback_3prime_guarded`** (lines 319+) — the **DRS production path**. Same core scan, plus:
  - **`early_exit_homopolymer_check`** (lines 437–468, default ON, `min_homopolymer_len=4`) — returns
    `None` (NO correction) if the genome around the 3′ end lacks a `stop_base×4` homopolymer. **This is
    the recurring culprit.** It SKIPS the walk wholesale. The [2026-06-04] fix added a `_base_at_3p !=
    stop_base` bypass (line 462–467) so lone-A termini still get walked — a point patch on this one guard.
  - **Guard 1 — large-deletion pre-scan** (restricts scan range around aligner-inserted ≥5bp deletions).
  - **Guard 2 — intron-boundary** (real N-ops clip the scan; artifact N-ops crossed; via
    `artifact_n_ref_starts`).
  - **Guard 3 — tail-context false-stop** (`tail_context_k=4`; rejects candidate anchors inside the poly-A tail).

Guards 1–3 **restrict where the walk may anchor** (legitimate artifact rejection — keep them). The
`early_exit_homopolymer_check` **decides whether to walk at all** (a suppressor — the problem child).

## 2. Evidence it recurs and that it matters (don't re-derive)

- Sumner human chr5 DRS, May-26 outputs (pre-[2026-06-04] fix): pooled **7.09%** of corrected ends on a
  gene-strand A; signature **+strand 105:1, 98% no-movement, 79% `atract_ambiguity`** (= the suppressed-walk
  reads). Re-corrected at HEAD (with the [2026-06-04] bypass) → **0.08%**. So the current patch *works* —
  but it's the 3rd such patch (git log: `30d2280` remove Case-2, `77ced6e` bridging guard, `1b1db38`
  lone-A bypass, `a1728eb` cat1 regression). The structural cause (guard-soup) is untouched.
- yeast is ~1.6% (lighter); human A-rich 3′UTRs + v5.2.0 basecaller make it ~6× worse — so human is the
  sensitive test bed.
- **Coordinate trap that has bitten analysis (not the bug itself):** `corrected_3prime` is
  **0-BASED-INCLUSIVE** (`# 0-based inclusive` at walkback.py:142,214; AGENT_FIXES 1059–1066). Measure
  on-A with 0-based genome fetch. A 1-based read produced a phantom "77% on A / over-walk" false alarm.

## 3. Goal / invariant to enforce

**Post-condition invariant:** after DRS walkback, `corrected_3prime` is on a stop-base (A) only when the
read genuinely has no inward non-stop anchor in the poly-A region (vanishingly rare). No guard may leave
the end parked on an A for a read that has a valid non-A anchor inward.

## 4. Two implementation options (recommend evaluating both; B is the safety net)

**Option A — guard-logic inversion (the principled fix).**
Reframe every guard as a constraint on the *anchor position*, never as a reason to *not walk*:
- Demote `early_exit_homopolymer_check` to a **pure performance hint** that can short-circuit ONLY when
  it provably cannot change the result (i.e., the 3′ base is non-stop AND already matches genome — the
  same condition as the core's terminal gate at line 168–170). Otherwise always run the scan.
- Keep Guards 1–3 as anchor-range restrictions (they reject real artifacts — large-del bridges, real
  introns, tail false-stops). They must never cause a *stop-base* anchor; if the only in-range anchor
  is a stop-base, that's a signal the restriction was wrong for this read.

**Option B — post-condition fallback (minimal, robust safety net; can ship alongside A).**
After `walkback_3prime_guarded` returns, if the result lands on a stop-base, **re-run the clean core
`walkback_3prime`** and prefer its (non-A) anchor. This guarantees the §3 invariant regardless of which
guard misfires, with minimal blast radius. Cheap insurance even if Option A is done well.

Recommended: implement **B as a guaranteed invariant** (it makes the whole class of regression
impossible), and do **A** to remove the now-redundant suppressor logic so the code stops accreting guards.

## 5. CI regression metric (REQUIRED — this is how it stops coming back)

Add a test that fails if the regression returns:
- On a fixed DRS fixture (extend `tests/test_walkback_readvsref.py`, or a small bundled DRS BAM),
  assert **corrected-on-stop-base rate ≤ ~2%** AND **no large +/- strand skew** in that rate (the skew is
  the canary — the bug is +strand-specific). Pick the threshold from the clean-core baseline on the fixture.
- Keep the existing **byte-identity** checks (guarded vs `find_polya_boundary` on Cat1–9) green — the
  refactor must not change correct results, only stop the stop-base parking.

## 6. Guardrails (do NOT regress these)

- **Byte-identical on the bundled Cat1–9 validation reads** (`tests/test_walkback_readvsref.py`) — the
  guards exist to pass these + reject mapPacBio GAA-bridging / large-del / real-intron artifacts. Don't
  delete artifact rejection.
- **No splice-motif gating** anywhere (lab principle; see memory `no-motifs-unbiased-discovery`).
- Walkback is **protocol-agnostic** (DRS both strands, cDNA, QSrev share this core). Changes must keep
  QSrev's `right_side_bridging_guard` path and cDNA behavior intact. Run the full walkback test set.
- `corrected_3prime` stays **0-based-inclusive** (don't "fix" the convention).

## 7. How to verify on real data (reuse, don't rebuild)

- Sherlock Sumner verification scripts: `<drive>/collaborations/sumner_lab/apa_3prime/{onA_check.py,
  onA_decomp.py,onA_all.py}` (0-based, correct). Re-correct a sample and confirm on-A ≈ HEAD's 0.08–0.27%.
- Re-correction submit template: `/scratch/users/kevinroy/sumner_lab/scripts/20_recorrect_HEAD.sh`
  (owners + AVX-512 constraint + idempotent skip-check, per CLAUDE.md). Outputs to `correct_HEAD_*/`.

## 8. Deliverables checklist

- [ ] Option B invariant implemented (guarded never emits a stop-base anchor; fallback to core).
- [ ] Option A: `early_exit_homopolymer_check` demoted to provably-safe short-circuit (or removed).
- [ ] CI on-A-rate + strand-skew regression test added and green.
- [ ] Cat1–9 byte-identity + full walkback test set still green.
- [ ] Re-verify on one Sumner sample (on-A ≤ ~0.3%); note result in AGENT_FIXES.
- [ ] CHANGELOG + AGENT_FIXES updated; this TODO marked done.

---

## ADDENDUM (2026-06-13, Kevin) — target is **100% non-A**; SUPERSEDES the "≤2%" framing in §3/§5

The body above (drafted by an agent) hedged toward "on-A ≤ ~2%" citing PolyASite's cross-tissue
observation that ~2% of CPA sites sit on a genomic A. **That biological caveat does NOT apply here —
RECTIFY's DRS policy overrides it.** Authoritative design policy:

- **DRS captures the poly-A tail by definition.** Assume **every** read 3′ end carries an (A)₁₀₊ tail.
  A DRS read's 3′ end IS the cleavage site + poly-A; DRS reads truncate at the **5′** end, never the 3′.
  So there is **no** "read legitimately ends mid-transcript on a genomic A" case to protect (delete that
  category from §1/§3 reasoning).
- **Poly-A not visible in the read ⇒ basecalling error** (Dorado mis-/under-called the homopolymer). The
  true 3′ end is still the CPA, so **still walk back.**
- **RECTIFY policy = left-shift every 3′ end to the first non-A, uniformly.** Therefore the TARGET is
  **100% of corrected DRS 3′ ends on a non-A.** A nonzero on-A fraction (incl. the current 0.08%) is an
  **inconsistency/bug to eliminate**, not a feature. Partial correction is sloppy and wrong by policy.
- Genuine genomic-A CPAs are **deliberately left-shifted too** — because (a) an A-CPA is indistinguishable
  from a tail-A aligned over it, (b) DRS semantics make the poly-A assumption universal, (c) consistency.
  **Do NOT add exceptions to "preserve" A-CPAs.**

**Implications (override §3–§5 accordingly):**
- The `early_exit_homopolymer_check`'s genomic-A-tract proxy is **invalid** under this policy: most CPAs
  are non-A and the poly-A may be basecall-truncated, so "genome has no A×4 near the 3′ end" is NOT a
  reason to skip the walk. It must **never** short-circuit a read whose 3′ end is on a stop-base.
- **Option B (post-condition) is a HARD GUARANTEE, not a safety net:** no DRS 3′ end may remain on a
  stop-base when any inward non-stop anchor exists in the read's aligned span. If none exists within
  `max_scan_depth`, **extend the scan** — an anchor must exist inward unless the entire aligned read is a
  stop-base run (enumerate + log those as the ONLY permitted residual; expected ≈ 0).
- **CI metric (revised, override §5):** assert corrected-on-stop-base rate **≈ 0%** (target 0 on the
  fixture; alarm above ~0.1%) and enumerate every residual with its reason. Keep the +/- strand-skew canary.
- **Acceptance (revised §8):** Sumner re-corrected on-A should be **~0%**, not 0.08% — decompose the
  current 0.08% (`onA_decomp.py` on `correct_HEAD_20260612/`: guard-misfire vs scan-depth-exceeded vs
  all-A-read) and drive it to zero, don't accept it.
