# 683 — G-aware 5′ clip arbiter: EXPLORED AND NOT NEEDED

**Status: PARKED, deliberately. Do not land. 2026-08-12, Kevin's call.**
Author of the work: session `cdna-trim-fix` (unit 683). Note filed by `rectify-realigner`.

## What it was

A guard for the PCB114 bridge-G: after `pretrim_consensus` strips the SSP + 27-nt UMI + GGG
bridge, a residual 1-bp `G` soft clip is left on a large share of molecules, and
`splice_aware_5prime.py` reads 5′ clips with **no minimum length** while engaging within
`DEFAULT_JUNCTION_PROXIMITY_BP = 10` of an intron edge. Yeast 3′SS is `AG` with the exon
starting immediately after, so a spurious terminal `G` sits exactly where it could tip an
acceptor placement by 1 bp. The arbiter added extra scrutiny for that case.

The reasoning was sound and the implementation is correct and well-gated. It is parked because
it was **measured**, not because it was doubted.

## Why it is not needed — the measurement

**0 of 23,663 suspect bridge-G clips fall within 10 bp of a 3′SS** (unit 683 CP2). The bridge G
sits at the **TSS** — the 5′ end of the transcript — while the 5′ rescue arbitrates **3′ splice
sites**. The two are at opposite ends of the molecule, so the hazard the guard defends against
does not arise in this chemistry. **The guard is correct and fires on nothing.**

⚠️ **Two numbers here are easy to conflate, and one supersedes the other in this context:**

| number | what it measures | status |
|---|---|---|
| ~34 % of molecules | 1-bp G **prevalence** after trimming | true, and not the relevant one |
| **0 / 23,663** | how often such a clip **reaches the rescue arbiter** | **the decisive figure** |

The prevalence figure was stated before the reach figure was measured, and an earlier relay of
this work (including mine, to Kevin) described the exposure as "frequent, not rare" on the
strength of prevalence alone. That framing is withdrawn.

## Where the work lives

Nothing was lost — it was extracted deliberately, not discarded:

- `~/work/UCLA/Chanfreau_Lab/planning/683_arbiter.patch` — reapply with `git apply`
- `~/work/UCLA/Chanfreau_Lab/planning/683_test_bridge_g_arbiter.py`
- `~/work/UCLA/Chanfreau_Lab/planning/683_g_aware_5p_clip_arbiter.md` — design + measurements

(It was extracted from the shared working tree because two sessions were editing one checkout.)

## When to revisit

Only if the premise changes — i.e. if a chemistry or protocol puts the bridge G, or any
systematic 1-bp 5′ artefact, at the **3′** end of the molecule where the rescue arbitrates. Then
re-measure the 0/23,663 figure first; the patch becomes relevant only if that number moves.

## The generalizable lesson

A guard that is correct, well-gated, and fires on nothing is still not worth landing: it adds a
code path and a maintenance surface for zero measured benefit. The discipline that caught this
was measuring **reach** (does the condition actually arrive at the code that cares?) rather than
**prevalence** (does the condition exist?). Prevalence is the easier number and the misleading
one.
