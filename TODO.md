# RECTIFY TODO

Known gaps, planned improvements, and deferred work.
Add items here rather than inline `# TODO` comments where possible.

---

## Alignment / Consensus

### 5' Soft-Clip Rescue — Sequence-Based Matching
**File:** `rectify/core/consensus.py` — `score_alignment()`
**Priority:** High

The current 5' rescue is a blind soft-clip length penalty (`score -= five_prime_softclip * 2`).
This causes aligners that *missed* a real intron (and therefore start after it with 0 soft clip)
to beat aligners that *found* the intron but have unrescued 5' bases.

Observed example: read SRR32518284.448567 at SNC1 (chrI:87,446-87,878).
- minimap2 + gapmm2: found the intron, reference start = 87,329, 25bp 5' soft clip → score penalised
- mapPacBio: missed the intron, reference start = 87,445, 0bp soft clip → wins incorrectly

**Intended implementation:**
1. When a read has a 5' soft clip, take those clipped bases (from `read.query_sequence`)
2. Find annotated junctions within range of the read's 5' alignment start (using the
   junction BED already loaded during consensus)
3. Fetch the genomic sequence at the end of the upstream exon (same length as the clip)
4. Compute edit distance between soft-clipped query bases and exon-end reference sequence
5. If edit distance ≤ threshold (~20% of clip length), the clip is "explained" — reward
   the aligner rather than penalising it

`edit_distance()` already exists in `rectify/core/spikein_filter.py` and can be imported.

---

## Statistics / Observability

### Tie-Break Rate per Sample
**File:** `rectify/core/consensus.py` — now tracks `stats['tied_score']`
**Status:** Counter added (2026-03-29). Not yet surfaced in HTML report or corrected_3ends_stats.tsv.

### Aligner Combo Breakdown per Sample
**File:** `rectify/core/consensus.py` — now tracks `stats['by_aligner_combo']`
**Status:** Logged at INFO level (2026-03-29). Not yet surfaced in HTML report or stats TSV.

---

## Performance

### mapPacBio Index Caching
**File:** `rectify/core/multi_aligner.py` — `run_map_pacbio()`
**Status:** `nodisk` removed, `path=bbmap_index/` added (2026-03-29). Index pre-built for
bundled S. cerevisiae genome. New genomes will build on first run and cache thereafter.

---

## Analysis / Downstream

### Bedgraph and Genomic Distribution in Manifest Mode
**File:** `rectify/core/analyze_command.py` — `_run_analyze_manifest()`
**Priority:** Medium

Manifest mode currently skips bedgraph generation and genomic distribution analysis
because both require per-read alignment coordinates (`alignment_start`, `alignment_end`),
which are not stored in the position index. Two options:
1. During Pass 2, fall back to streaming the full TSV for samples that need these columns
2. Add a `rectify export --bedgraph --manifest` subcommand (cleaner — bedgraphs are
   large and users may not always want them)

### Expose Aligner Stats in HTML Report
**Files:** `rectify/core/consensus.py`, HTML report template
**Priority:** Low

`stats['tied_score']` and `stats['by_aligner_combo']` are tracked but not shown in
the HTML report or `corrected_3ends_stats.tsv`. Surfacing these would help diagnose
aligner quality per sample.

---

## Known Non-Issues

- `minratio=0.4` in `run_map_pacbio()` cmd — duplicates BBTools default, harmless.
- Combined TSV (`corrected_3ends_combined.tsv`) in `combined/` is now obsolete for the
  standard pipeline (manifest mode replaces it) but harmless to leave in place.
