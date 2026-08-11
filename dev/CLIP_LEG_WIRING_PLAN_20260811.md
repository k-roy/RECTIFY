# Clip-leg wiring plan — triage terminal-clip + 5′-clip legs (2026-08-11)

Branch: `feat/triage-clip-legs` (from `master` @ `5df199e`).
Implements landscape doc `dev/REALIGNER_LANDSCAPE_AND_PATH_20260809.md` §4 step 2:
wire the triage clip legs — 5′ clips → the Cat3/rescue path, terminal (3′) clips →
the resolver machinery. One refusal discipline across both. Everything OFF by
default (`TriagePolicy.clip_legs_enable = False`); default behavior byte-identical.

## 1. API contract (verified by reading, this worktree)

### Triage layer — `rectify/core/consensus/triage.py`
- `TriagePolicy` dataclass: `max_junction_proximal_errors=1.0, max_clip_5p=30,
  max_clip_3p=30, triage_unannotated_junctions=True, junction_window_bp=5`.
- `TriageResult`: `read_id, label, reasons, junction_proximal_errors, clip_5p,
  clip_3p, n_junctions, n_unannotated`; property `junction_leg`.
  Reasons include `REASON_CLIP_5P='softclip_5p'`, `REASON_CLIP_3P='softclip_3p'`.
  `clip_5p`/`clip_3p` are STRAND-AWARE (from `extract.get_softclip_lengths`:
  5′ = left clip for `+`, right clip for `-`).
- `reentry_accept(original, candidate, genome, penalty_table=None) -> bool` —
  the STRICT arbiter: accepts only when `_cigar_hp_edit_distance(candidate)` <
  `original - 1e-9` AND `(refname, ref_start, cigartuples)` signature differs.
  `_cigar_hp_edit_distance` (corrected_consensus.py:59) is CIGAR+genome based
  (no MD/NM needed; S/H cost 1.0/base, N free) — so a genuine clip→exon/N
  resolution strictly improves it.
- `triage_realign_bam(input_bam, output_bam, genome, annotated_junctions,
  policy=None, pool_bams=None, penalty_table=None, realign=True,
  refine_kwargs=None, sort_and_index=True) -> (rows, stats)`.
  Junction leg: classify → build_junction_pool → refine_bam_junctions
  (motif_blind) on triaged junction-leg reads → per-read `reentry_accept` →
  final BAM with accepted replacements swapped in.

### Resolver — `rectify/core/align/overhang_resolver.py`
- `resolve_read(read, genome, index, cfg, stats) -> bool` — per-read entry;
  MUTATES the read in place (CIGAR/pos rewrite, XJ tag, MD/NM dropped);
  handles BOTH terminal clips, then (if `cfg.arb_enable`, default True) the v2
  junction re-arbitration. Chrom resolution: `standardize_chrom_name` then
  `CHROM_TO_GENOME` fallback.
- `ResolverConfig`: `alpha=0.01, max_intron=5000, min_clip=8,
  max_clip_match=200, max_edit_frac=0.2, min_margin=1.0, arb_enable=True, …`.
- `ResolverStats`: `refused_low_info`, `refused_repeat`, `resolved`,
  `candidates_evaluated`, `no_candidates`, `rejected_edit`,
  `rejected_ambiguous`, `extra{}` …
- Refusal is INTERNAL to `resolve_clip`: `is_repeat_expansion` then
  `assess_overhang(clip_used, alpha, max_window=max_intron)`; refused clips are
  never candidate-searched (641 T1).

### Shared gate — `rectify/core/splice/overhang_informativeness.py`
- `assess_overhang(seq, alpha=0.01, max_window=None, hp_discount=0.5)
  -> OverhangAssessment(length, i_eff_bits, w_max_bp, refused, alpha)`;
  bumps module `COUNTERS['assessed'/'refused']`. `reset_counters()` for tests.
- `rescue_3ss_truncation` already carries this gate INTERNALLY but DARK by
  default (`RECTIFY_OVERHANG_INFO_GATE=1` env enables; splice_aware_5prime.py
  lines 1436–1452 — refused ⇒ `rescue_seq=''`, structural paths stay live).

### Splice-site index — `rectify/core/splice/splice_site_index.py`
- `SpliceSiteIndex.build(genome_dict)` / `load_or_build(genome_path, genome)`;
  `sites_in(chrom, kind, lo, hi)` binary-search range query.

### Cat3 rescue path — `rectify/core/splice/splice_aware_5prime.py` +
### application machinery `rectify/core/bam/read_edits.py` / `bam_writer.py`
- `rescue_3ss_truncation(read, genome, candidate_junctions, strand=None,
  max_edit_frac=0.2, junction_proximity_bp=…, scan_bp=50, terminal_peel=True, …)
  -> dict` — does NOT mutate the read (try/finally restore). Candidate
  junctions are `(chrom, intron_start, intron_end)` 3-tuples with
  STANDARDIZED chrom names. Internally: skip-region bypass, repeat-expansion
  short-circuit on the 5′ clip, reanchor pre-pass, sequence rescue, terminal
  peel, intronic snap, proximity.
  Result keys: `rescued, rescue_type ('softclip'|'mpb_mismatch'|
  'intronic_snap'|'proximity'|'none'), five_prime_corrected,
  rescued_junction, edit_distance, query_bp`, plus (on sequence rescue)
  `five_prime_exon_cigar`, `five_prime_upstream_trim`, `reanchor_clip_len`.
- BAM materialization (how bam_writer applies it, `apply_corrected_edits_to_read`
  lines 266–307): (1) `_apply_reanchor_from_clip_len(read, reanchor_clip_len)`
  (read_edits.py:1719) when reanchor fired; (2)
  `extend_read_5prime_for_junction_rescue(read, five_prime_position,
  soft_clip_len, strand, exon_cigar_str, upstream_trim) -> bool` (in-place;
  converts the 5′ S op to exon ops + N; requires a live 5′ S op, else returns
  False); (3) the intronic reroute/softclip variants for in-intron 5′ ends.
- Correct-stage invocation (`bam_processor.py:491`):
  `_rescue_3ss(read, genome, _ss_junctions, strand)` with `_ss_junctions` =
  annotated ∪ read-derived ∪ nearby pool junctions.

### Bundled data (test harness, matches tests/test_triage.py)
- `rectify/data/validation/aligners/validation_reads.minimap2.bam` (refs
  `chrI…` style), genome `rectify/data/genomes/saccharomyces_cerevisiae/
  S288C…fsa.gz` — `load_genome` returns `chrI`-keyed dict (verified),
  annotation gff.gz via `load_annotated_junctions` → `(chrom, s, e, strand)`
  with standardized chroms.

## 2. Wiring design

### Config gate
- `TriagePolicy.clip_legs_enable: bool = False` (same dataclass style). With
  the flag False, `triage_realign_bam` takes the exact pre-change code path:
  no index build, no resolver import side effects on output, BAM byte-identical.

### Routing (one leg per failure mode, exclusive by clip end)
- `TriageResult.clip5_leg` property: `REASON_CLIP_5P in reasons` → Cat3 leg.
- `TriageResult.clip3_leg` property: `REASON_CLIP_3P in reasons` → resolver leg.
- A read carrying both reasons goes through BOTH legs, each on its own clip
  end — this is what enforces "every overhang assessed exactly once": the
  resolver is restricted to the 3′-end genomic side, the Cat3 leg owns the 5′
  clip.

### Terminal-clip (3′) leg
1. `resolve_read` gains an optional `sides: str = 'LR'` kwarg (deviation from
   brief, see §5): gates the LEFT/RIGHT clip blocks. Default `'LR'` preserves
   behavior for all existing callers. Triage passes the genomic side of the 3′
   clip: `'R'` for `+` strand reads, `'L'` for `-`.
2. Triage builds/receives ONE shared `SpliceSiteIndex` (new kwarg
   `splice_site_index=None`, built from `genome` on demand) and ONE
   `ResolverConfig` (new kwarg `resolver_config=None`; default
   `ResolverConfig(arb_enable=False)` — arbitration is the junction leg's /
   Station-C's territory, the clip leg stays clip-scoped).
3. Per routed read: copy the incumbent (`AlignedSegment.from_dict`), run
   `resolve_read(copy, genome, index, cfg, resolver_stats, sides=side)`.
   - Refusals (Δ of `refused_low_info + refused_repeat`) → `clip3_refused`,
     read passes through unchanged — first-class outcome.
   - `resolve_read` returned True → a PROPOSAL (`clip3_proposed`) → enters the
     existing arbiter: `reentry_accept(incumbent, copy, genome, penalty_table)`
     decides (`clip3_accepted`). The resolver proposes, triage disposes.
4. §4b seam: `_two_sided_candidate_hook(read, index, cfg)` — documented no-op
   returning `[]`, called after the resolver in the terminal-clip leg;
   post-landing resolver-library work (two-sided A/B enumeration) plugs in
   there. NOT implemented here by design.

### 5′-clip leg (Cat3 rescue)
1. Leg-level gate (the ONE refusal discipline): extract the 5′ clip
   (strand-aware), take its junction-proximal ≤`max_clip_match` bases
   (mirroring `resolve_clip`), then `is_repeat_expansion` +
   `assess_overhang(clip_used, alpha=cfg.alpha, max_window=cfg.max_intron)`.
   Refused → `clip5_refused`, passthrough, rescue NEVER invoked (the refused
   overhang is never sequence-searched). Note: `rescue_3ss_truncation`'s own
   internal gate stays env-dark (default off), so the overhang is assessed
   exactly once, in the leg. (If a user exports RECTIFY_OVERHANG_INFO_GATE=1
   the rescue re-assesses its imperfection-derived rescue_seq — a different
   sequence; documented, not fought.)
2. Rescue invocation mirrors bam_processor: candidate set = annotated
   3-tuples (chrom-standardized). `rescue_3ss_truncation(copy, genome,
   cand_set, strand)`.
3. Materialization reuses the bam_writer machinery, not a reimplementation:
   on `rescued` + sequence rescue (`rescue_type` in {'softclip',
   'mpb_mismatch'}): `_apply_reanchor_from_clip_len` if `reanchor_clip_len>0`,
   then `extend_read_5prime_for_junction_rescue(copy, five_prime_corrected,
   soft_clip_len=reanchor_clip_len or clip5_len, strand, exon_cigar_str,
   upstream_trim)`. Extension returning False (e.g. in-intron geometry the
   extension can't express) ⇒ no proposal. `intronic_snap`/`proximity`
   results are TSV-level in the correct stage; the triage leg emits no BAM
   proposal for them (MVP scope; counted as non-proposals).
4. Proposal → same arbiter: `reentry_accept` (`clip5_proposed` /
   `clip5_accepted`).

### Ordering & bookkeeping in `triage_realign_bam`
- Clip legs run AFTER the junction leg; the incumbent for a clip-leg read is
  the junction-leg-accepted record when one exists, else the original —
  proposals always compete against the current best under the same hp_ed gate.
- New stats keys (unconditional, zero when disabled): `clip5_leg, clip3_leg,
  clip5_refused, clip3_refused, clip5_proposed, clip3_proposed,
  clip5_accepted, clip3_accepted`. Existing `realigned`/`accepted` stats stay
  junction-leg-only (preserves the existing smoke inequality); the per-read
  rows' `realigned`/`accepted` fields cover both legs (accepted clip-leg reads
  are also marked realigned).
- Row schema unchanged (no new columns) — flag-off rows byte-identical.

## 3. Test plan (`tests/test_triage_clip_legs.py`)

Synthetic genome (test_overhang_resolver style): ~4 kb random seq, planted
`+` intron [1200,1500) (GT..AG), chrom `chrI`; annotated = {(chrI,1200,1500,'+')}.

(a) **Default-off identity**: TriagePolicy().clip_legs_enable is False; run
    `triage_realign_bam` on the bundled minimap2 validation BAM with
    `policy=None` and with `clip_legs_enable=False`, `sort_and_index=False`;
    SAM text byte-identical, rows equal, all clip stats zero, and
    `COUNTERS['assessed']` untouched (proves no leg/gate code ran).
    Dev-side (not in test): pre/post-change md5 of the default-path output.
(b) **Refusal first-class**: 3′ 40×A clip (poly(A)) and 5′ (AAG)-repeat/low-
    complexity clip on synthetic reads, legs ENABLED → `clip3_refused`/
    `clip5_refused` increment, zero candidates evaluated, output records
    byte-equal to input, nothing accepted.
(c) **Terminal-leg proposal → adjudication**: + strand read, 60M aligned
    exon-1 + 40 bp right clip = exon-2 head → resolver proposes M/N/M,
    arbiter accepts (`clip3_proposed==1, clip3_accepted==1`), output read has
    the N op + XJ tag; also assert the STRICT gate by feeding a garbage-clip
    read that yields no proposal.
(d) **5′-leg gate behavior**: + strand read starting at intron_end with a
    32 bp 5′ clip = exon-1 tail → gate assesses ONCE
    (`COUNTERS['assessed']==1` for that run), rescue proposes, arbiter
    accepts (`clip5_accepted==1`), output has exon+N structure; refused case
    covered in (b).

Suites: `python -m pytest tests/test_triage*.py tests/test_overhang_resolver.py -q`
then full gate `pytest -m "not slow" -q` (pipefail; only acceptable
pre-existing failure: test_consensus_tag_restoration.py::test_process_batch_…).

## 4. Commit sequence
1. this plan doc; 2. `resolve_read(sides=…)`; 3. triage clip legs + config
gate; 4. tests; 5. gate run + status appends here.

## 5. Deviations from the brief
- `resolve_read` gains a `sides='LR'` kwarg (default-preserving). Reason: the
  brief demands both "hand terminal-clip reads to resolve_read" AND "every
  overhang assessed exactly once" — but stock `resolve_read` assesses BOTH
  clips, and 5′ clips belong to the Cat3 leg. The kwarg is the minimal way to
  satisfy both clauses.
- Branch created from `5df199e` exactly as the brief states; actual `master`
  tip is `2701f3b` (one docs-only handoff commit ahead). Noted, not rebased.

## 6. Status checkpoints
- [2026-08-11] Plan written; contract verified by reading all listed modules.
