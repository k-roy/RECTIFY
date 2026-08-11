# 644i — the Station-C repeat-context flag, three ways: annotation is free, self-homology earns its keep, variant-multiplicity is REJECTED

**Date:** 2026-08-11. **Context:** Kevin's 644h correction — the copy-misplacement blindness
of the overhang-quality score is handled by a direct repeat flag (rDNA, CDS tandem arrays,
CUP1/ENA-class families, Ty/LTR, subtelomeres), not a per-candidate uniqueness computation.
**Flag semantics: DEMOTE to require-orthogonal-evidence, never discard.** Tool:
`scripts/benchmark/644i_stationc_repeat_flag.py` (runs locally against the 644h census +
bundled GFF; optional `--selfhom-paf`). Frame: non-canonical track, phase-1 survivors
(644h baseline q≥80 → 11/14 Gould gold, 1,537 junk).

## The three candidate flags, measured

| flag (composed with q≥80) | junk passing (of 1,537) | gold passing (of 11) | verdict |
|---|---:|---:|---|
| annotation (rRNA/Ty/LTR/tRNA/telomere/X/Y′ ± 500 bp + rDNA region) | 1,318 (−14%) | 11 (−0) | **ADOPT — free**: zero gold cost at every threshold |
| + genome self-homology track | **971 (−37%)** | 9 (−2) | **ADOPT** — the CDS-tandem-array leg |
| + variant-multiplicity (≥5 jitter-collapsed variants/locus) | 110 (−93%) | 1 (−91%!) | **REJECT as repeat proxy** |

- **Self-homology track**: `minimap2 -DP -k19 -w19 -m200 G G` → 3,334 segments, 1.2 Mb
  (9.9% of R64) after ±200 bp merge; seconds of one-time compute per genome. Demotes 2–3
  gold — all singletons (chrXVI:795133, chrIV:306804, chrXIII:551950) that genuinely sit in
  self-homologous sequence, i.e. exactly the reads whose placement needs orthogonal
  evidence. Canonical track cost: ZERO (11/11 gold + 3/3 shortlist junk all unflagged).
  Sensitivity is tunable (lower `-m`, `asm20` presets) — unexplored headroom.
- **Variant multiplicity REJECTED**: the number of distinct (jitter-collapsed) non-canonical
  variants at a locus measures *mapPacBio disagreement density*, which is high at repeat
  spray AND at genuinely contested/novel loci alike — at every operating point it demotes
  most of the Gould gold (11/14 at the locus-chain version; 8–13/14 at per-candidate overlap
  versions with K=4–6). The 644h "jitter-twin votes" observation survives, but as evidence
  FOR a locus under two-sided enumeration, not as a repeat flag.

## Station-C phase-2 admission (updated spec, all terms measured)

```
candidate (non-canonical track, phase-1 survivor)
  → q = short-exon-side overhang quality (644h)
  → flags = annotation ∪ self-homology (644i)
  unflagged AND q≥θ            → admit on (q × recurrence)     [9/11 Gould gold at θ=80]
  flagged OR q<θ               → DEMOTE: admit only with orthogonal evidence
                                  (short-read/COMPASS, cross-sample recurrence,
                                   mm2-side distress co-localization — 644d machinery)
```

Residual after both flags at q≥80: 971 junk vs 9 gold (≈0.9%) — dispersed, unique-sequence
junk. The remaining discriminators are the orthogonal-evidence classes (644f §consequence),
of which mm2-side distress co-localization (644d census) is the cheapest unmeasured probe.

**Files:** tool `scripts/benchmark/644i_stationc_repeat_flag.py`; self-homology recipe above
(PAF not committed — regenerate in seconds); builds on `PHASE2_OVERHANG_644H_20260811.md`
(644h) and the 644f/644g census docs.
