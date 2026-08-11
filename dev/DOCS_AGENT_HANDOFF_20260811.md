# Documentation Agent handoff — README / ARCHITECTURE / RTD refresh after the Re-aligner landing

**Written 2026-08-11 by the rectify-realigner agent. Status: READY TO EXECUTE (pass 1).**
Kevin's decision on timing: **two passes.** Pass 1 (now): document what is LANDED and PUBLIC —
`origin/master` jumped `255a06d → 03bc44f` (~400 commits) on 2026-08-11, so README/RTD are
materially stale about behavior users can already install. Pass 2 (later): the panel-narrative
and Station-C docs, AFTER Station C + phase 2 are built and the `[[617]]` leave-one-out settles
mapPacBio's panel status — writing those now would mean rewriting them.

## Build facts

- Docs = **mkdocs** (`mkdocs.yml`) + **ReadTheDocs** (`.readthedocs.yaml`). Verify with
  `mkdocs build --strict` before committing; treat warnings as errors.
- README.md at repo root. `docs/ARCHITECTURE.md`, `docs/index.md`, `docs/quickstart*.md`,
  `docs/user_guide/commands/*.md`, `docs/ALIGNER_RECOMMENDATIONS.md`.
- **Do NOT push to origin** — commit locally on a `docs/refresh-20260811` branch and hand
  back; master pushes are Kevin-gated.
- Known suite red (pre-existing, not yours): `test_consensus_tag_restoration` — ignore;
  docs-only diffs don't gate on it.

## Kevin's own edit list — reconcile FIRST

`Rectify_readme_changes_KR_proposed.txt` (repo root of the MAIN checkout
`~/work/rectify`, untracked) — Kevin's proposed README changes. Incorporate:
- Species-agnostic positioning (yeast → human; error profiles similar).
- Add branding to README (see `dev/branding/` canonical-asset decision — if unresolved, ask).
- DRS 5′-end caveat wording (degradation OR incomplete pore transit — corrects alignment
  artifacts, does not guarantee TSS recovery; cDNA distinguishes full-length via adapters).
- Walkback wording (pA-tail basecall errors → non-genomic tail fragments aligning to genome).
- The chimeric-reconstruction "explore how well this works" item is a SCIENCE task, not docs —
  leave it; flag to the realigner agent.
- Also in that file: the v1.0.0 checklist (CHANGELOG, version bump, doc consolidation P7
  git-mv batch, README_KR_edits.md reconciliation) — reconcile if those files exist.

## Pass 1 scope — landed, user-facing, currently undocumented or wrong

1. **Overhang resolver** (align stage): `core/align/overhang_resolver.py`, registered as the
   `overhang_resolver` arm in `multi_aligner.py` (see `AlignerConfig`, default OFF). Document:
   what it does (information-bounded soft-clip splice resolution on minimap2 output; refusal
   first-class), when to enable, measured value (mm2 769→852 gold junctions on 4.57 M upf1Δ
   DRS reads at 334 junk, ~10 CPU-h), and that it is NOT yet default-panel. Source of truth:
   `~/work/UCLA/Chanfreau_Lab/planning/644h_realigner_integration.md` §1–2 +
   `docs/`-side module docstrings. Human-genome use: NOT yet recommended (T6 decoy
   mitigation pending) — say so.
2. **`RECTIFY_SKIP_REGIONS`** env var (+ `yeast-rdna` shorthand): user-facing perf knob,
   measured 2.6× on chrXII, byte-identical non-rDNA output. Document in the HPC/perf page.
3. **`rectify triage` CLI**: consensus-triage layer — document as **EXPERIMENTAL**, one
   paragraph + `--help` dump; do not present as a pipeline stage yet.
4. **ONT-cDNA Path A default** (`3d893ef`): UMI-collapse-to-molecules is now the run-all
   default for ONT cDNA; both entry points. Update `quickstart_cdna.md` +
   `user_guide/commands/correct_cdna_overview.md`.
5. **Containment-first gene attribution default** (`7648725`): new DEFAULT in analyze;
   `region_class` per cluster, transcript model, ncRNA atlases. Update analyze/output docs
   (`output_formats.md` — new columns!).
6. **Browser-pack / QC final stage** (`5792b22`): run-all's fail-soft final stage. One
   section in the run command doc.
7. **run-all manifest fixes** (`70e9664`) — behavior now matches docs; check the multi-sample
   page (`docs/architecture/multi_sample_pipeline.md`) still describes reality.
8. ARCHITECTURE.md: add the align-stage resolver arm + consensus-stage triage layer to the
   stage diagram; the three-station framing gets ONE forward-looking paragraph (Station A
   landed; B = triage, experimental; C = in design) — no more until pass 2.

## Explicitly OUT of scope (pass 2, do not document now)

- Station C / phase-2 scorer / repeat flags (644h/644i are dev-only findings).
- Two-sided enumeration design; mapPacBio scout/probe role; any panel-recommendation change
  (all pending `[[617]]` LOO).
- `RECTIFY_OVERHANG_INFO_GATE` / `RECTIFY_OVERHANG_INFO_COUNTS` (dark/dev knobs).
- The 644-series numbers as marketing claims — cite only the resolver ladder fact in (1).

## Verification checklist before handing back

- `mkdocs build --strict` clean; README renders on GitHub (check tables/anchors).
- Every documented flag/command verified against `--help` output or the code (verify-don't-
  confabulate: run `rectify triage --help`, grep the AlignerConfig defaults).
- No doc references to pre-public 2.x/3.x versions (house rule, see CLAUDE.md).
- Diff is docs-only (`git diff --stat` shows no .py).
