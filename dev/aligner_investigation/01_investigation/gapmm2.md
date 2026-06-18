# gapmm2 — Source-Level Investigation

**Tool:** `gapmm2` ("gapped alignment using minimap2 — align transcripts to genome")
**Author:** Jon Palmer (`nextgenusfs`) · **Repo:** https://github.com/nextgenusfs/gapmm2 · **PyPI/Bioconda:** `gapmm2`
**Role in RECTIFY:** Tier‑1 aligner in the 5‑aligner ensemble; invoked as `gapmm2 -t <n> -i 5000 -o out.paf genome reads`; PAF (with `cs` tags) is converted to BAM. Win rate ≈ **0.8%**.

> **FACT vs INFERENCE.** Statements grounded in RECTIFY source (read directly) or quoted from the public gapmm2 README/source are marked **[FACT]**. Statements reasoned from behavior, RECTIFY's notes, or minimap2/edlib internals are marked **[INFER]**. The public gapmm2 source is relatively sparse and has shifted across versions (`master`/`main` branch, v25.4.x wheels), so some algorithmic detail is reconstructed from README excerpts + RECTIFY's defensive code rather than a line-by-line read of the installed wheel.

---

## 1. Relationship to minimap2 / mappy

**[FACT]** gapmm2 is a thin Python wrapper around minimap2. It does **not** implement its own seed‑chain‑align core; it delegates the primary spliced alignment to minimap2 and then performs a *post‑hoc terminal refinement* pass using `edlib`. Its declared dependencies are **mappy** (the Python minimap2 binding), **edlib** (fast banded/global pairwise aligner), and **natsort**.

**[FACT]** Two invocation paths exist depending on version:
- Early versions call mappy in‑process: `mp.Aligner(reference, preset="splice", n_threads=threads)`.
- **As of v25.4.13** gapmm2 "will run `minimap2` directly instead of mappy if `minimap2` is installed," because of **segmentation faults observed with the mappy Python bindings.** The direct command observed is:
  ```
  minimap2 -x splice --cs=short -t <threads> reference query
  ```

**[FACT]** The preset is `splice` (i.e. minimap2's spliced long‑read RNA preset) and `cs` short‑form tags are requested. The `splice` preset implies minimap2's `-k15 -w5 ... --splice` spliced‑alignment parameters (canonical GT‑AG modeling). **[INFER]** gapmm2 does **not** expose minimap2's `-uf` forward‑only flag, junction‑BED hints, or `-G`; the public CLI surfaces only a small set of flags (below). The intron bound it cares about is its **own** `-i/--max-intron`, which gapmm2 uses to size the *terminal refinement search window*, not (necessarily) as minimap2's `-G`. **[INFER]** Because the splice preset's default max‑intron applies to the minimap2 pass, the `-i` value primarily governs how far past the alignment edge gapmm2 will look for a missing terminal exon.

### Why minimap2 needs help at all
**[FACT]** minimap2's chaining can **drop very short terminal exons** entirely: a seed in a 13 bp 5′ exon is too short to anchor a chain across the intron, so minimap2 soft‑clips/omits it. gapmm2's README documents exactly this — a 543 bp transcript whose 13 bp 5′ exon (`409609–409621`) was missing from minimap2's output and recovered by gapmm2:
- minimap2: `[(409320,409554),(409090,409255),(408904,409032)]`
- gapmm2: `[(409609,409621),(409320,409554),(409090,409255),(408904,409032)]`

---

## 2. CLI options and defaults

**[FACT]** (from README / `--help`):

| Flag | Long | Default | Meaning |
|------|------|--------:|---------|
| `reference` | — | (req) | reference genome FASTA |
| `query` | — | (req) | transcripts/reads FASTA or FASTQ |
| `-o` | `--out` | stdout | output path |
| `-f` | `--out-format` | `paf` | `paf` or `gff3` |
| `-t` | `--threads` | **3** | minimap2 threads |
| `-m` | `--min-mapq` | **1** | min mapping quality; alignments below are dropped (`low-mapq` counter) |
| `-i` | `--max-intron` | **500** | "controls terminal search space" for refinement |
| `-d` | `--debug` | False | debug to stderr |

**[FACT — RECTIFY]** RECTIFY overrides the intron default with **`-i 5000`** (`run_gapmm2`, `multi_aligner.py:1116`, comment "Max intron size (yeast splice range)"). This widens the terminal search window from 500 → 5000 bp so gapmm2 can recover terminal exons across longer yeast introns; the gapmm2 default of 500 would miss them. **[FACT — RECTIFY]** RECTIFY deliberately does **not** pass `-m`: the `'-m 1'` it once passed (mislabeled "Mode 1: direct RNA") is actually gapmm2's `--min-mapq`, and some 25.4.5 wheels shipped that argparse arg without `type=int`, so `-m 1` made `min_mapq` the string `"1"` and raised `TypeError` on `h.mapq < min_mapq`. Default int `1` is what is wanted, so RECTIFY omits the flag (`multi_aligner.py:1110‑1115`).

**API:** **[FACT]** `aligner(reference, query, out_fmt=...)` returns a stats dict `{'n', 'low-mapq', 'refine-left', 'refine-right', ...}`. Helper `cs2coords()` parses a `cs`/CIGAR into genomic exon/intron coordinates.

---

## 3. Terminal‑Gap Refinement Algorithm

**[FACT]** The core idea: after minimap2 aligns a read, gapmm2 inspects each primary hit for **unaligned query bases at the termini** and tries to re‑place them as a terminal exon with `edlib`.

**Detection.** **[FACT]** For a hit `h` over read `seq`:
- `if h.q_st > 0:` → there are unaligned **5′/left** query bases → attempt **left** refinement.
- `if len(seq) > h.q_en:` → unaligned **3′/right** query bases → attempt **right** refinement.

**Realignment with edlib.** **[FACT]** Strand‑specific helpers (`left_plus_best_align`, `right_plus_best_align`, `left_minus_best_align`, …) carve a reference window near the alignment edge and run a **semi‑global** edlib alignment of the orphaned query segment into it:
```python
# plus strand, left side
if refstart - maxlen < 0: r_st = 0
else:                     r_st = refstart - maxlen
ref = refseq.seq(contig, r_st, refstart + offset)      # offset default ~6 bp
query = queryseq[0 : querystart + slide]
align = edlib.align(query, ref, mode="HW", k=0, task="path",
                    additionalEqualities=degenNuc)
```
Key points **[FACT]**:
- **`mode="HW"`** = "infix"/semi‑global: the *query* must be fully consumed but it can land anywhere in the reference window (free end gaps on the reference). This is exactly right for "where does this short orphan exon sit upstream of the intron."
- **`maxlen` = `-i/--max-intron`** bounds the reference window (`refstart - maxlen … refstart + offset`), i.e. the intron length gapmm2 is willing to span. This is why RECTIFY raises it to 5000.
- **`k=0`** lets edlib find the optimal edit distance with no band cap (within the window); `task="path"` returns a CIGAR.
- **`additionalEqualities=degenNuc`** makes IUPAC degenerate bases (N, R, Y…) count as matches — tolerant to ambiguity codes.
- A **`slide`** parameter scans candidate splice‑donor/acceptor offsets; helper scans (`find_all_splice_CT`, `filter4_splice_AC`, and plus‑strand GT/AG analogues) bias the chosen intron boundary toward **canonical splice dinucleotides**.

**Stitching back into `cs`.** **[FACT]** The refined terminal exon's edlib CIGAR is converted to `cs` ops and an explicit intron token is synthesized with the canonical donor/acceptor and length:
```python
# plus strand left:
new_cs = "cs:Z::{}~gt{}ag{}".format(align["cigar"].strip("="), intron_len, new_cs_str)
# minus strand left:
new_cs = "cs:Z:{}~ct{}ac:{}".format(new_cs_str, intron_len, align["cigar"].strip("="))
```
So gapmm2 emits a `~gtNNNNag` (plus) or `~ctNNNNac` (minus) splice token joining the recovered terminal exon to the rest of the alignment, and rewrites the PAF record's `cs`, query start/end, and target coordinates. **[FACT]** A refinement is counted only when the `cs` actually changed: `if paf[-1] != cs: stats["refine-left"] += 1`.

**Yield (README benchmark).** **[FACT]** On 6,926 alignments: **409 left‑refined, 63 right‑refined** (≈ 6.8% touched, heavily 5′‑biased). **[INFER]** The left bias reflects that 5′ terminal exons are most often the short ones dropped by chaining.

**Sequence lookup.** **[FACT]** gapmm2 retrieves the read sequence for refinement via the mappy/minimap2 query index: `query_idx.seq(paf[0])` keyed on the PAF query name. **[INFER]** This is the single most fragile assumption for RECTIFY (see §5): the lookup returns `None` if the name is duplicated in the indexed FASTQ.

---

## 4. PAF + `cs` Output and BAM Reconstruction (RECTIFY side)

**[FACT]** gapmm2 emits **PAF only** (or GFF3) — never SAM/BAM. PAF lines are `qname qlen qstart qend strand tname tlen tstart tend nmatch alnlen mapq tags…`. **PAF carries no SEQ/QUAL** — only coordinates and `cs`. Therefore RECTIFY must inject sequences when building a BAM.

**`cs` long‑form grammar** (parsed by `_cs_long_to_cigar`, `multi_aligner.py:837`): **[FACT]**
- `:n` → `n` matches → CIGAR `nM`
- `*xy` → 1 bp substitution → CIGAR `1X` (mismatch)
- `+seq` → insertion (query only) → `|seq|I`
- `-seq` → deletion (ref only) → `|seq|D`
- `~nnNNNnn` → splice/intron; the digit run is the intron length → `NNNN N` (splice skip)

**[FACT]** Because `cs` covers only the aligned block, RECTIFY recomputes terminal **soft‑clips** from PAF coordinates: `left_clip = query_start`, `right_clip = query_len - query_end`, placed as `S` ops on the strand‑appropriate end (5′ = left for `+`, right for `−`).

**`_paf_to_bam` reconstruction** (`multi_aligner.py:894`): **[FACT]**
1. Build a header from the genome `.fai` (`SO:queryname`, one `@SQ` per contig, `@PG gapmm2`).
2. For each PAF line: extract `cs:Z`, `NM:i`, `tp:A`. Records with **no `cs` tag are skipped**.
3. Non‑primary records (`tp:A` ≠ `P`, i.e. secondary `S` / inversion `I`) get FLAG `0x100` so only one primary per read survives downstream filtering.
4. Validate position: `target_start < 0` or unknown contig → **skip** (gapmm2 can emit out‑of‑bounds artifacts).
5. **Sequence injection** from the **original** FASTQ via `_load_fastq_sequences` (`{name:(seq,qual)}`). Plus strand → forward sequence; minus strand → `_reverse_complement(seq)` and reversed qualities (pysam expects alignment orientation).
6. Name‑sort the BAM (`pysam.sort -n`) so consensus selection can stream across aligners; then `_apply_calmd_eq` recomputes MD/`=`/`X`.

---

## 5. RECTIFY Integration Quirks (verbatim‑accurate)

All of the following are RECTIFY‑specific hardening, grounded in `multi_aligner.py` + CLAUDE.md.

### 5.1 Cleaned, deduplicated FASTQ is mandatory (`_clean_fastq`, `run_gapmm2`)
**[FACT]** `run_gapmm2` always writes a temp cleaned FASTQ before calling gapmm2 because of two distinct issues:

**Issue 1 — DRS auxiliary‑tag headers.** `samtools fastq -T pt` (DRS poly‑A trimming) embeds a tab‑separated tag in the header: `@UUID\tpt:i:25`. minimap2 truncates the name at the first whitespace → PAF query name is the bare UUID; mappy strips tabs similarly. RECTIFY strips the tag anyway (`header[1:].split()[0]`) "for robustness."

**Issue 2 — Duplicate UUIDs with empty sequences (root cause of a hard crash).** DRS‑trimmed FASTQs from Dorado contain ~1–2% reads with **duplicate UUIDs**: one entry has an empty (placeholder) sequence, one has the real sequence. **When mappy indexes a FASTQ containing duplicate read names, `seq(name)` returns `None` for *both* entries** regardless of which has real sequence. This makes gapmm2's refinement loop call `len(None)` → **`TypeError` at `align.py:883`**, crashing gapmm2 after processing only the reads *before the first duplicate* — deterministically ~22k of 695k. Diagnostic signature: the output PAF has exactly N lines where N = aligned reads before the first duplicate UUID.

**The fix in `_clean_fastq` (`multi_aligner.py:51`)** writes a FASTQ that:
1. strips DRS tags from names (`@UUID\tpt:i:N` → `@UUID`),
2. skips reads with empty sequence (`if not seq.rstrip(): continue`),
3. skips any UUID seen more than once (`seen_uuids` set).

**[FACT]** `run_gapmm2` passes this cleaned FASTQ to gapmm2, but passes the **original** `reads_path` (tags intact) to `_paf_to_bam` for sequence injection (`_load_fastq_sequences` re‑strips the tags so keys match the bare‑UUID PAF names).

### 5.2 `cs`‑overrun bug — query length mismatch (~0.02% of reads)
**[FACT]** gapmm2 "occasionally emits `cs` tags that over‑consume 1–4 query bases past `query_end`." In `_paf_to_bam`, RECTIFY recomputes the CIGAR's query length (`sum` of `M/I/S/X/=` ops) and **if `cigar_qlen != len(fwd_seq)` the record is skipped** (counted in `n_skipped`) rather than attempting CIGAR surgery — pysam would reject the malformed record otherwise. **[INFER]** This is the cs‑path analogue of the affine‑extension edge case where edlib's recovered terminal exon CIGAR slightly over‑runs the true query span.

### 5.3 `-m` / `--min-mapq` TypeError on some wheels
**[FACT]** See §2: RECTIFY omits `-m` entirely because some 25.4.5 wheels lack `type=int` on the argparse arg, turning `-m 1` into a string and raising `TypeError` on the mapq comparison. The default (int `1`) is the desired behavior.

### 5.4 Out‑of‑bounds / negative target positions
**[FACT]** `_paf_to_bam` validates `target_start >= 0` and that the contig exists, skipping artifacts. **[INFER]** These arise when edlib places a recovered terminal exon at a window edge and the reconstructed coordinate underflows.

### 5.5 Subprocess hygiene
**[FACT]** Run under `ALIGNER_TIMEOUT` (6 h). Empty/zero‑byte PAF → `RuntimeError`. Temp cleaned FASTQ deleted in `finally`. Non‑zero return code surfaces `stderr`.

### 5.6 Scoring footing in consensus (`consensus/scoring.py`)
**[FACT]** RECTIFY's consensus treats gapmm2 as a soft‑clipping aligner: `five_prime_softclip_seq` from "minimap2/gapmm2" is used as the rescue sequence (scoring.py:134, 706). It is also penalized for **junction‑proximity errors** — aligners (minimap2, gapmm2) that "accumulate mismatches/indels within a few bp of each junction" lose score relative to clean‑junction aligners like mapPacBio (scoring.py:745). **[INFER]** This penalty is part of why gapmm2 rarely wins (§6).

---

## 6. Strengths

- **[FACT]** **Recovers short terminal exons minimap2 drops** — the entire reason the tool exists; demonstrated on a 13 bp 5′ exon. Directly relevant to RECTIFY's mission (accurate terminal‑exon/3′‑end placement).
- **[FACT]** **Canonical‑junction‑aware terminal placement.** The edlib refinement biases the synthesized intron toward `gt…ag` / `ct…ac` donor/acceptor pairs, so recovered terminal exons tend to land at real splice sites.
- **[INFER]** **`HW` (infix) edlib + degenerate‑base equalities** gives tolerant, optimal placement of the orphan segment within a bounded window — good on noisy nanopore ends and ambiguity codes.
- **[FACT]** **Inherits minimap2's splice quality** for the bulk of each read, adding refinement only at the ~7% of termini that need it — cheap marginal cost.
- **[INFER]** **Independent failure modes** from the other ensemble members (it is the only one whose terminal exons come from edlib, not chaining), so it occasionally supplies the unique correct 3′ end no other aligner found — justifying its inclusion despite low win rate.

## 7. Weaknesses / Failure Modes

- **[FACT]** **PAF‑only output.** No SEQ/QUAL, no SAM/BAM → RECTIFY must inject sequences and synthesize soft‑clips/CIGAR from `cs` + coordinates (`_paf_to_bam`). Extra fragility surface.
- **[FACT]** **Duplicate‑UUID / empty‑seq crash.** mappy's `seq()` returns `None` for *all* copies of a duplicated name → `TypeError`, deterministic partial PAF. Mitigated only by RECTIFY's mandatory `_clean_fastq` dedup pass.
- **[FACT]** **`cs`‑overrun (~0.02%).** Over‑consumes 1–4 query bases past `query_end`; those records are dropped (lost reads), not repaired.
- **[FACT]** **mappy segfaults** → upstream switched to spawning `minimap2` directly (v25.4.13). **[INFER]** Implies the Python‑binding path is unstable on some inputs/builds.
- **[FACT]** **Argparse `--min-mapq` type bug** on some 25.4.x wheels (string vs int) → `TypeError`.
- **[FACT]** **Out‑of‑bounds/negative coordinates** can be emitted and must be filtered.
- **[FACT/INFER]** **Rarely wins in the ensemble (~0.8%).** Causes: (a) its body alignment *is* minimap2, so where minimap2 is right gapmm2 ≈ minimap2 and a dedicated aligner (deSALT, mapPacBio) usually scores ≥ as well; (b) it only changes ~7% of termini, so for 93% of reads it offers nothing distinctive; (c) RECTIFY's junction‑proximity penalty docks minimap2/gapmm2 for near‑junction mismatches, so even on refined termini a cleaner‑junction aligner (deSALT 78.9%, mapPacBio 18.2%) wins the read. gapmm2's value is the rare read where its edlib terminal recovery is uniquely correct — high precision, low frequency, low marginal value vs deSALT.
- **[INFER]** **Default `-i 500` is too small for yeast introns** — without RECTIFY's `-i 5000` override gapmm2 would miss exactly the terminal exons it is meant to recover.

---

## 8. Source References

**gapmm2 (upstream):**
- Repo: https://github.com/nextgenusfs/gapmm2 — README (worked 13 bp exon example; 6,926‑aln benchmark: 409 left / 63 right; v25.4.13 minimap2‑direct note; defaults `-i 500 -t 3 -m 1`).
- `gapmm2/align.py` — `aligner()` loop (`h.q_st > 0` / `len(seq) > h.q_en`), `left_plus_best_align` / `left_minus_best_align` edlib refinement (`mode="HW", k=0, task="path", additionalEqualities=degenNuc`, `maxlen`‑bounded window), `cs` reconstruction (`~gt{len}ag`, `~ct{len}ac`), `query_idx.seq(paf[0])`, `cs2coords()`, refine counters.
- Dependencies: mappy, edlib, natsort. PyPI/Bioconda package `gapmm2`.

**RECTIFY (read directly):**
- `rectify/core/align/multi_aligner.py`
  - `run_gapmm2` (~L1052): cleaned FASTQ, `-t -i 5000 -o`, `-m` omission rationale, timeout, PAF→BAM call, `_apply_calmd_eq`.
  - `_paf_to_bam` (~L894): header build, `cs`/`NM`/`tp` parse, primary/secondary FLAG, position validation, sequence injection, `cigar_qlen` mismatch skip, name‑sort.
  - `_clean_fastq` (L51): dedup + empty‑seq skip + tag strip.
  - `_load_fastq_sequences` (L88): bare‑UUID keys for injection.
  - `_cs_long_to_cigar` (L837): long‑form `cs` → CIGAR, soft‑clip derivation.
  - `_reverse_complement` (L47), `ALIGNER_TIMEOUT` (L119).
- `rectify/core/consensus/scoring.py` (L134, L706, L745): gapmm2 as soft‑clip aligner; junction‑proximity penalty.
- `CLAUDE.md` — "gapmm2: cleaned FASTQ required (DRS auxiliary tags + duplicate UUIDs)"; "gapmm2 PAF → BAM: sequence injection required"; win‑rate table (gapmm2 0.8%).
