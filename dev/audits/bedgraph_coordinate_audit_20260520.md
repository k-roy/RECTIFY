# Bedgraph coordinate audit — `corrected_3prime` 0-based vs BED half-open

**Date opened:** 2026-05-20
**Date closed:** 2026-05-20
**Triggered by:** `analyses/cross_modality_trt_20260519/MANUSCRIPT_ANCHORED_TRT.md` §4.4 — CST6 cross-modality check exposed a 1-bp left shift in `rectify/core/analyze/bedgraph.py::generate_bedgraphs`.
**Status:** **AUDIT COMPLETE.** Three emitters fixed (`analyze/bedgraph.py`, `analyze/manifest.py`, `scripts/generate_bedgraph_from_polished.py`); all other bedgraph writers under `rectify/core/` verified to already follow the 0-based half-open convention. All fixes uncommitted on M1 working tree at close-of-audit.

**Scope-expansion note (2026-05-20, post-close):** initial audit was scoped to `rectify/core/` plus `rectify/data/`. A parallel multi-agent audit run by the user expanded the search into `scripts/` and adjacent coordinate-convention code; one bedgraph-scope miss (`scripts/generate_bedgraph_from_polished.py:113`) and five non-bedgraph coordinate-convention findings were surfaced. The bedgraph miss is now fixed in this audit's scope; the five adjacent findings are tracked separately in `dev/BUGS_TO_FIX.md` (NEW-077 through NEW-081) since they are not bedgraph emitters. **Audit recipe updated below** to grep `scripts/` and project-side helpers in future audits.

## The convention this audit defends

Throughout `rectify`, the `corrected_3prime` column on the corrected-positions TSV is **0-based-inclusive** — the 0-based reference coordinate of the last aligned (and walked-back, where applicable) base on the genome. This convention is established in pysam and inherited end-to-end:

- `read.reference_end - 1` for `is_reverse=True` (3' end is the rightmost aligned base, 0-based)
- `read.reference_start` for `is_reverse=False` (3' end is the leftmost aligned base, 0-based)

Documented inline at multiple call sites with the comment `# 0-based inclusive`:
- `rectify/core/correct/walkback.py:142, 214, 471`
- `rectify/core/correct/indel_corrector.py:1661, 2027`
- `rectify/core/correct/protocols/quantseq_rev.py:72, 77`

**Therefore** any bedgraph emitter that turns a `corrected_3prime` value into a BED row must write `start = pos, end = pos + 1` — the half-open interval `[pos, pos+1)`. Any emitter that writes `start = pos - 1, end = pos` is shifting the position 1 bp to the LEFT of the true coordinate.

## Bug signature

```python
# CORRECT — pos is 0-based; BED is 0-based half-open
start = int(pos)
end   = int(pos) + 1

# BUGGY — treats pos as 1-based when it isn't
start = int(pos) - 1
end   = int(pos)
```

Inverse pattern (legitimate, do NOT "fix"): if `pos` was loaded from a 1-based source (e.g. GFF, SAM POS column, a TSV that explicitly documents 1-based), then `start = pos - 1` is correct.

## Audit targets — files that emit BED / bedgraph rows

The audit must check every site that writes `chrom\t<start>\t<end>` to a `.bedgraph`, `.bed`, or `.bw`-source file. Files identified at audit open:

| File | Status | Notes |
| --- | --- | --- |
| `rectify/core/analyze/bedgraph.py` | **FIXED** | `analyze` per-condition bedgraph (in-memory mode). Lines 99-100 corrected. 3 regression tests in `tests/test_analyze.py::TestBedgraphCoordinates`. See AGENT_FIXES.md `[2026-05-20]` entry. |
| `rectify/core/analyze/manifest.py` | **FIXED** | Manifest-mode equivalent of `generate_bedgraphs` (streams per-sample TSVs). Line 598 had the same off-by-1 spelling against the same 0-based `corrected_position`/`corrected_3prime` column. Corrected. No standalone unit test — relies on the documented invariant; the `analyze/bedgraph.py` tests pin the convention this file inherits. |
| `rectify/core/bam/bedgraph_writers.py` | VERIFIED CORRECT | `write_netseq_assigned_bedgraph` and `write_corrected_reads_bedgraph` both delegate the final row emission to `netseq/netseq_output.py::write_bedgraph` (line 48 / 144 imports). Convention propagates correctly. |
| `rectify/core/netseq/netseq_output.py` | VERIFIED CORRECT | `write_bedgraph` line 122 emits `f"{chrom}\t{pos}\t{pos + 1}\t..."` with an explicit comment `# BedGraph: 0-based, half-open coordinates`. `write_bigwig` line 191 uses `[p[0] + 1 for p in positions]` — also correct. |
| `rectify/core/commands/export_command.py` | VERIFIED CORRECT | `write_bedgraph` line 88 emits `f"{chrom}\t{pos}\t{pos+1}\t..."` with comment `# bedGraph: chrom, start (0-based), end, value`. `write_bigwig` line 128 uses `chrom_df['position'].values + 1` for ends — also correct. |
| `rectify/core/commands/analyze_command.py` | VERIFIED (delegate) | Delegates to `analyze/bedgraph.py::generate_bedgraphs` at line 163. No direct emission. Inherits the fix. |
| `rectify/core/commands/consensus_command.py` | VERIFIED (delegate) | `_generate_bedgraphs` at line 156 imports `write_bigwig, write_bedgraph` from `export_command` and dispatches. No direct emission. |
| `rectify/core/commands/correct_command.py` | VERIFIED (delegate) | `write_netseq_assigned_bedgraph` (line 1101), `write_corrected_reads_bedgraph` (line 1125), and a direct call to `netseq_output.write_bedgraph` (line 1140 import, line 1165 call). All paths route through the verified-correct emitters. |
| `rectify/core/commands/run_command.py` | NOT AN EMITTER | Sole `bedgraph` mention is a docstring at line 561 — not a writer. |
| `rectify/core/commands/split_command.py` | NOT AN EMITTER | Sole `bedgraph` mention is `--no-bedgraph` flag inside a template string at line 1242 — not a writer. |
| `rectify/data/validation/generate_igv_html.py` | NOT AN EMITTER | References existing bedgraph paths in HTML output (lines 408-410). Does not write bedgraphs. |
| `scripts/generate_bedgraph_from_polished.py` | **FIXED** (scope-expansion post-close) | Standalone CLI utility (not the production analyze path). Line 113 had the same `start = pos - 1; end = pos` against `corrected_3prime` (which `detect_position_column` at line 36 prefers over the legacy `position` column). Comment on lines 111-112 incorrectly asserted "Position is 1-based". Corrected with an extended comment recording the pre-fix state. Caught by the multi-agent audit run; missed by the initial audit because the recipe didn't grep `scripts/`. |

(Audit source: `grep -rEln "bedgraph|to_bedgraph|write_bedgraph" rectify/ --include="*.py"` 2026-05-20; full row-emission survey: `grep -nE "f\.write.*chrom.*\\t|\bpos[+-]?1?\b.*\\t" rectify/`. **Expanded post-close** to include `scripts/` and `analyses/` helper modules.)

## Verification recipe — per file

For each target:

1. Open the file and grep for every `start` / `end` assignment that writes to a BED row.
2. For each, trace the variable holding the position back to its source:
   - If it's `corrected_3prime`, `corrected_pos`, or anything derived from `reference_end - 1` / `reference_start` → **0-based** → must be `start = pos, end = pos + 1`.
   - If it's loaded from `pos`, `position`, or a 1-based source → 1-based → `start = pos - 1, end = pos` is correct.
   - If it's ambiguous, add a clarifying comment and a regression test that fixes the semantic.
3. Add or extend a unit test under `tests/` that constructs a tiny input with a known coordinate and asserts the emitted BED row matches BED 0-based half-open semantics.
4. Update this audit doc's table: mark FIXED + commit ref, or VERIFIED (no fix needed) with the line(s) inspected.

## Downstream blast radius

`rectify analyze` per-condition bedgraphs from every prior production run are 1 bp shifted. Specifically:

- Every set1/set2/set3 DRS `rectify analyze` output's `bedgraph/` subdir.
- Han 2023 wbfix bedgraphs in `han2023_analyze_wbfix_20260512/analyze/bedgraph/` (the run that triggered this audit).
- Any PCR-cDNA `rectify analyze` outputs.

**What downstream products use them and need to be re-checked:**

- IGV per-condition tracks loaded for visual inspection. The shift is sub-pixel at typical zooms but visible at the per-base zoom. Documents like figure panels that show the bedgraph at single-base resolution should be re-rendered after the fix lands and bedgraphs are regenerated.
- Any cross-track analysis that compared the `rectify analyze` per-condition bedgraph against an externally generated 3'-end track (bedtools, deepTools, etc.). The He et al manuscript-anchored TRT analysis was the trigger; others may exist.

**What is NOT affected:**

- Clustering, shift analysis, and per-position attribution all read directly from the `corrected_reads.tsv` / `corrected_3ends.tsv` position columns, not from the bedgraphs. The high-confidence TRT set in `analyses/cross_modality_trt_20260519/derived/table_s1_trt_genes_v2.tsv` and the v0.9.2 reprocess classification do not change with this fix.

## Action items

1. ☑ Audit each "NOT AUDITED" file above (2026-05-20, this session).
2. ☑ Add regression tests covering the convention (`TestBedgraphCoordinates` in `tests/test_analyze.py`).
3. ☑ For confirmed buggy emitters, apply the fix and update this doc + AGENT_FIXES.md (2 emitters fixed: `analyze/bedgraph.py`, `analyze/manifest.py`).
4. ☐ **Commit fixes + tests on `drs-validation-rebuild`, push, pull on H2 + Sherlock.** Pending — uncommitted on M1 at audit close.
5. ☐ **Decide whether to regenerate prior per-condition bedgraphs in production output dirs.** Recommendation: do not regenerate retroactively — the shift is sub-pixel at typical IGV zooms and the per-position attribution work that drives publication figures reads from `corrected_reads.tsv` directly. Document the shift in run-PROVENANCE for affected outputs and regenerate naturally on the next `rectify analyze` invocation.
6. ☐ Close out by adding a CHANGELOG entry once the commit lands.

## Related references

- Bug report + first fix: `AGENT_FIXES.md` entry `[2026-05-20] rectify analyze per-condition bedgraph: 1-bp left shift`
- Trigger analysis: `analyses/cross_modality_trt_20260519/MANUSCRIPT_ANCHORED_TRT.md` §4.4
- Open-bugs index entry: `dev/BUGS_TO_FIX.md` NEW-076
