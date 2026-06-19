# uLTRA — Source-Level Technical Investigation

> Verified against `origin/drs-validation-rebuild` @ 366c885 (2026-06-19).

**Tool:** uLTRA — annotation-guided long-read RNA splice aligner (Sahlin & Mäkinen, *Bioinformatics* 2021).
**Repo:** `github.com/ksahlin/ultra` (GPL-3.0).
**Role in RECTIFY:** Tier-2, annotation-**REQUIRED** aligner. Invoked from `rectify/core/align/multi_aligner.py::run_ultra` (~line 1265) as:

```
uLTRA pipeline --ont --disable_infer --t <n> --prefix <p> <ref.fa> <annotation.gtf> <reads.fastq> <outdir>
```

uLTRA wins ~2% of reads in RECTIFY's consensus. This report explains, at the source level, why: its annotation-guided model produces *very precise canonical junctions* on annotated isoforms but *structurally cannot win* reads whose true structure diverges from the GTF.

Throughout, **[FACT]** = verified from source/paper; **[INFERENCE]** = reasoned from evidence; **[RECTIFY]** = from the RECTIFY wrapper.

---

## Annotation Segmentation Model

uLTRA's core idea is to **stop treating the genome as a flat string** and instead index the annotation as a partitioned set of exon segments per gene, then chain a read across those segments.

**[FACT] Indexing entry point.** The `pipeline`/`index` subcommands call `create_augmented_gene.create_graph_from_exon_parts()`, which parses the GTF with **gffutils** and builds the segmentation. `--disable_infer` (used by RECTIFY) tells gffutils to skip *inferring* gene/transcript features from exon lines — it requires `gene` and `transcript` features to already exist in the GTF, which "can speed up the indexing considerably." **[RECTIFY]** this is why `run_ultra` normalizes SGD GTFs (`_normalize_gtf_for_ultra`: `mRNA`→`transcript`, derives `exon` lines, forces `transcript_id != gene_id`) before invoking uLTRA — `--disable_infer` would otherwise crash on a GTF lacking those features.

**[FACT] "Parts" — the first-level decomposition.** Exons are clustered into non-overlapping genomic regions called **parts**. The construction loop walks sorted exon records per chromosome: `if exon.start - 1 > active_stop + 20` it opens a *new part*; otherwise the exon extends `active_stop` (overlapping/proximal exons within 20 bp merge into one part). Key structures:
- `part_count_to_choord` — `(start, stop)` coordinates per part.
- `part_to_active_gene` — part → gene id.
- `part_to_canonical_pos` — `(chr_id, part_counter)` → set of exon boundary positions.
- `pos_to_exon_ids` — exon ids indexed by position and boundary type (start/end).

**[FACT] "Segments" — the second-level decomposition.** `get_canonical_segments()` cuts each part at every annotated exon boundary into **canonical segments**, producing:
- `parts_to_segments`, `segment_to_gene`, `segment_id_to_choordinates`, `segment_to_ref`.
- Each segment keyed as a packed byte string: `array("L", [chr_id, p1, p2]).tobytes()`.
- Segments respect a `min_segment_size` (tiny gaps merged).

So a gene becomes a linear tiling of disjoint segments; any annotated exon is a *union of consecutive segments*, and an isoform is a *subset path* over segments. This is the "exon-segment graph."

**[FACT] Splice graph / junction tables.** Per chromosome, transcripts are decomposed into junctions: `splices_to_transcripts`, `transcripts_to_splices`, plus annotated-junction lookup tables `all_splice_pairs_annotations` and `all_splice_sites_annotations`. These are the `valid_introns_sites` consulted later in `align.py`.

**[FACT] Flanks.** Inter-part / intergenic buffer regions are stored as `flank_ids` (same byte encoding as segments) so a read overhanging an annotated region's edge still has reference to align to.

**[FACT] Reference materialization.** `get_sequences_from_choordinates()` extracts the genomic sequence for each segment/flank (unpacking the byte coordinates), producing the **augmented reference** (`refs_sequences.fa` in the output dir) that seeding runs against. The whole index is pickled (`dill`) into `database.db`. **[RECTIFY]** `run_ultra` treats a non-empty `refs_sequences.fa` as a valid, reusable cache across chunks, and rebuilds only if it is empty/missing.

---

## Seeding (NAMs / MEMs / strobemers)

**[FACT] Two seeding backends exist; the default has changed across versions.**

1. **MEMs (original, paper default).** First pass uses **Maximal Exact Matches** between the read and the *transcriptome/segment* sequences. The paper's rationale: MEMs have *variable length*, so "provide more information on the relevance of the hit than fixed-length seeds." The paper used **slaMEM** to compute MEMs.
2. **NAMs via strobemers (current default).** Modern uLTRA delegates seeding to **namfinder** (a separate binary, a relative of `StrobeMap`/strobealign). namfinder finds **NAMs** (Non-overlapping Approximate Matches) using **strobemers** — randomized multi-part k-mers robust to indels/mismatches. namfinder is ">10× faster" and more memory-efficient than the old MEM finder, and produces smaller intermediate files. The legacy MEM path is still reachable via `--use_NAM_seeds`/`--nams` toggles; a newer `--nams` collinear-chaining mode reproduces minimap2-style chaining for mapping accuracy.

**[FACT][RECTIFY] namfinder is a hard runtime dependency in RECTIFY.** `run_ultra` checks `shutil.which('namfinder')` and, if absent, prepends a vendored `namfinder` binary's directory to `PATH` via `_extra_env`; otherwise it raises `FileNotFoundError`. This confirms RECTIFY runs the **NAM/strobemer** seeding path, not slaMEM.

**[FACT] `--ont` preset and seeding params.** `--ont` sets `min_acc=0.6`, strobe size `s=9`, and the minimap2-fallback k-mer `mm2_ksize=14` (vs `--isoseq`: `min_acc=0.8`, `s=10`). `min_acc` is the minimum accuracy threshold for a MAM to be considered. A `--thinning` (0–2) parameter controls seed reduction. **[INFERENCE]** the lower `min_acc=0.6` for ONT reflects the higher per-base error rate of nanopore vs IsoSeq.

**[FACT] Role of `--disable_infer` in seeding.** None directly — `--disable_infer` is purely an *indexing* speedup (skips gffutils feature inference). It does not change the seeding algorithm. **[RECTIFY]** it is mandatory here only because the normalized yeast GTF already supplies `gene`/`transcript`/`exon` features.

---

## Colinear Chaining over Exon Segments

This is uLTRA's algorithmic heart — a **two-pass collinear chaining** scheme. **[FACT]** implemented in `modules/colinear_solver.py`, which explicitly cites **"Algorithm 15.1 in *Genome-Scale Algorithmic Design*, Mäkinen et al."**

### Pass 1 — MEM/NAM chaining (read coverage)

**[FACT]** Functions `read_coverage(mems, max_intron)` (O(n²)) and `n_logn_read_coverage(mems)` (O(n log n)). The fast variant builds two **range-max segment trees** (`RMaxQST`, "T" and "I" trees, leaves padded to a power of 2) to answer range-maximum queries in O(log n) instead of a linear predecessor scan. The DP at each seed `v` (interval `[c..d]` on the read) computes:
- `C_a` = best chain ending in a **disjoint** predecessor: `(v.d - v.c + 1) + max_case_a`.
- `C_b` = best chain ending in an **overlapping** predecessor: `max_case_b + (v.d - prev_end)` (only the non-overlapping tail of `v` is credited).
- `C[j+1] = max(C_a, C_b)`, with backpointers for traceback.

The objective is **maximum read coverage** by a colinear chain of seeds — i.e., explain as much of the read as possible with monotonically increasing read/reference coordinates, subject to `max_intron`. This selects the **candidate gene(s)** the read belongs to.

### Pass 2 — MAM chaining with overlaps and gaps

**[FACT]** This is the "novel chaining problem" the paper highlights. Functions `read_coverage_mam_score(mams, overlap_threshold)` and `n_logn_read_coverage_mams(mams, overlap_threshold)` chain **MAMs (Maximal Approximate Matches)** rather than exact matches. Distinctive properties verified in source:
- **Approximate matches:** MAMs are scored, not just exact; computed with **edlib** (see Base-Level Alignment).
- **Overlap cost:** overlapping MAMs are *allowed* with a fractional penalty `0.1 * (gap_length)`; overlaps beyond `overlap_threshold` are not considered. Genomic-coordinate consistency is checked.
- **Gap cost:** gaps between MAMs are penalized in the objective.

**[FACT] The decisive trick — inject all annotated exons.** Pass 2 adds **all annotated exons of the candidate gene(s), including exons that had no MEM/NAM seed at all.** This is what lets uLTRA align reads across *very small exons* (the paper's headline **8-nt exon** example, missed by minimap2 and deSALT): the segment is in the chaining DP as a candidate tile even though no seed landed in it. The chain then "tiles" exon segments onto the read.

**[INFERENCE]** Complexity: both passes are O(n log n) in the seed/MAM count via the segment-tree acceleration (O(n²) fallbacks exist). The chaining-with-overlaps formulation is the paper's novel contribution over standard sparse dynamic-programming chaining (which forbids overlaps).

---

## Novel vs Annotated Junctions

**[FACT] Junction classification at assembly time.** `align.py::find_exons()` walks the chosen MAM chain and decides, between consecutive MAMs `m1`,`m2`, whether the reference gap is an intron or a deletion:
- `if (m1.y, m2.x) in valid_introns_sites` → **annotated intron** (exact known donor/acceptor) → emit an `N` CIGAR op at the canonical coordinates.
- `elif m2.x - m1.y > 10` → large unannotated gap → treated as a **novel intron** (`N`).
- else → small gap → **deletion** (`D`).

**[FACT] SQANTI-style read categories.** `classify_alignment2.py` assigns each read one of `FSM` (Full Splice Match — all junctions annotated), `ISM/NIC_known`, `NIC_novel`, `NNC` (Novel Not in Catalog). FSM reads are privileged downstream: in `align.py`, the accuracy filter `if alignment_score < 2*alignment_threshold*len(read_seq)` rejects low-scoring alignments **unless the read is FSM**.

**[FACT] Why annotation-guidance gives precise canonical junctions.** Because annotated junctions are snapped to the exact `valid_introns_sites` coordinates from the GTF, uLTRA's `N`-op boundaries on FSM/known reads are *exactly* the annotated donor/acceptor — no aligner-noise wobble. This is uLTRA's strength and the reason RECTIFY keeps it as a Tier-2 annotation-required source.

**[FACT] The minimap2 safety net.** uLTRA does **not** force every read into the annotation. `prefilter_genomic_reads`/the mm2 QC step (controlled by `--disable_mm2`, k-mer `mm2_ksize=14` under `--ont`) classifies a read as **genomic** if more than `--genomic_frac` (default 10%) of its aligned length falls outside uLTRA-indexed regions; such reads default to the **minimap2** alignment, and uLTRA vs minimap2 results are merged by score. **[INFERENCE]** for an intergenic / fully-novel read, uLTRA effectively *is* minimap2 — so it cannot offer a better junction than minimap2 there, capping its win rate on novel structure.

---

## Base-Level Alignment

**[FACT] Two engines (`modules/help_functions.py`).**
- **parasail** (semi-global, vectorized): `parasail.sg_trace_scan_16(read, ref, opening_penalty, gap_ext, matrix)` where `matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)`. **Defaults: match=+2, mismatch=−2, gap_open=3, gap_extend=1.** A local variant `parasail_local` uses `sw_trace_scan_16` (Smith-Waterman). This produces the final per-exon CIGAR.
- **edlib**: `edlib.align(read, ref, task="path", mode=...)`. Used (a) to *score MAMs* (approximate-match edit distance) in Pass 2, and (b) as the alignment engine when `len(ref) > 20000 or len(read) > 20000` with `mode="HW"` (infix/semi-global) for memory efficiency. Default mode otherwise `"NW"` (global).

**[FACT] Assembly into a SAM record.** `align.py::get_exact_alignment()` builds the per-segment reference (`created_ref_seq`) for the chosen exon path, aligns the read to it (parasail, or edlib for very long refs), stitches the exon CIGARs with `N` ops at junctions (annotated or novel, per `find_exons()`), and `sam_output.main()` emits the SAM record. Mapping quality is set from competing alignments; `coverage = covered / len(read_seq)`.

**[FACT] Poly-A handling.** `help_functions.remove_read_polyA_ends()` trims terminal homopolymer `A`/`T` runs within a window of `min(len(seq)//2, 100)` bp before alignment (runs over `threshold_len` collapsed to `to_len`). **[INFERENCE]** uLTRA therefore *strips* the poly-A tail rather than modeling it — it does not attempt precise CPA/3'-end placement; its 3' end is wherever the trimmed read's genomic alignment ends. This matters for RECTIFY, whose entire purpose is 3'-end correction: uLTRA contributes a *clean spliced backbone with canonical junctions*, not a refined 3' terminus.

**[RECTIFY] Post-processing.** `run_ultra` collects `<outdir>/ultra.sam`, converts to coordinate-sorted BAM (`samtools view -bS | samtools sort`), then runs `_apply_calmd_eq` (regenerates MD/`=`/`X` ops against the genome). uLTRA cannot read gzipped FASTA/GTF, so RECTIFY decompresses to a temp dir first.

---

## Strengths

1. **[FACT] Small / micro-exon recovery.** Injecting *all* annotated exons of the candidate gene into Pass-2 chaining (even seedless ones) recovers reads spanning exons as short as 8 nt that seed-and-extend aligners (minimap2, deSALT) drop. This is uLTRA's signature advantage.
2. **[FACT] Exact canonical junctions.** Annotated junctions are snapped to GTF coordinates → zero junction wobble on FSM/known reads. **[RECTIFY/INFERENCE]** this is precisely where uLTRA *wins its ~2%* in RECTIFY: reads on well-annotated multi-exon isoforms where its junction boundaries beat the noisier de-novo aligners.
3. **[FACT] Principled overlap/gap chaining** with O(n log n) segment-tree DP — handles indel-rich ONT reads at the seed level via NAMs/strobemers.
4. **[FACT] Graceful novel-read fallback** to minimap2 (no silent forcing into wrong annotation).

## Weaknesses

1. **[FACT] Annotation-bound.** Quality is capped by GTF completeness. For yeast (compact, well-annotated, few isoforms) the upside is small — most reads are simple, so uLTRA rarely *out-classes* the de-novo aligners, explaining the **low ~2% overall win rate**. **[INFERENCE]** uLTRA shines on complex multi-isoform genomes (human/mouse); on *S. cerevisiae* its niche is narrow.
2. **[INFERENCE] Novel-isoform / intergenic blind spot.** Reads >10% outside indexed regions are handed to minimap2, so uLTRA can never beat minimap2 there — by construction it loses those reads in consensus.
3. **[FACT] Poly-A stripped, no CPA modeling.** 3'-end is just the trimmed alignment end; uLTRA offers RECTIFY no 3'-end precision (the thing RECTIFY corrects). Its value is the junction backbone, not the terminus.
4. **[FACT][RECTIFY] Operationally heavy/brittle.** Requires GTF (not GFF) with `gene`/`transcript`/`exon` features and `transcript_id != gene_id`; cannot read gzip; needs a separate `namfinder` binary and a minimap2 install; builds a large pickled `database.db` + `refs_sequences.fa` index. RECTIFY adds `_normalize_gtf_for_ultra`, decompression, namfinder vendoring, and stale-cache detection to keep it running.
5. **[INFERENCE] Memory/runtime.** Per-gene MAM chaining plus an all-exon Pass-2 candidate set and pickled index inflate memory vs a pure seed-chain aligner; the >20 kb edlib switch in `get_exact_alignment()` is an explicit memory guard.

---

## Source / Paper References

- **Paper:** Kristoffer Sahlin, Veli Mäkinen. *Accurate spliced alignment of long RNA sequencing reads.* **Bioinformatics 37(24):4643–4651, 2021.** DOI: 10.1093/bioinformatics/btab540. PMID 34302453. Open access: PMC8665758. Preprint: bioRxiv 2020.09.02.279208.
- **Repo:** `github.com/ksahlin/ultra` (GPL-3.0). Key modules referenced:
  - `uLTRA` (main): subcommands `index`/`align`/`pipeline`; `--ont` (`min_acc=0.6, s=9, mm2_ksize=14`), `--isoseq` (`min_acc=0.8, s=10`), `--disable_infer`, `--disable_mm2`, `--genomic_frac` (0.10), `--use_NAM_seeds`/`--nams`, `--thinning`, `--t`, `--prefix`, `--index`.
  - `modules/create_augmented_gene.py` — `create_graph_from_exon_parts`, `get_canonical_segments` (parts/segments/flanks).
  - `modules/colinear_solver.py` — `read_coverage`, `n_logn_read_coverage`, `read_coverage_mam_score`, `n_logn_read_coverage_mams` (Algorithm 15.1, Mäkinen et al.; RMaxQST segment trees).
  - `modules/classify_alignment2.py` — FSM/ISM/NIC/NNC categories.
  - `modules/align.py` — `find_exons`, `get_exact_alignment`; `valid_introns_sites`.
  - `modules/help_functions.py` — `parasail_alignment` (match 2 / mismatch −2 / gap_open 3 / gap_ext 1, `sg_trace_scan_16`), `parasail_local` (`sw_trace_scan_16`), `edlib_alignment` (HW/NW), `remove_read_polyA_ends`.
- **Dependencies:** namfinder (NAMs/strobemers, supersedes slaMEM MEMs), minimap2 (genomic fallback), parasail, edlib, gffutils, intervaltree, dill, pysam.
- **RECTIFY:** `rectify/core/align/multi_aligner.py::run_ultra` (~L1265), `_normalize_gtf_for_ultra`, `_get_vendored_binary('namfinder')`, `_apply_calmd_eq`.

*Note on FACT vs INFERENCE:* Algorithm names, function names, CLI parameters, parasail/edlib numerics, the two-pass structure, and the all-exon Pass-2 injection are **verified from the GitHub source and the paper**. Complexity attributions to specific passes, the ~2%-win rationale, and yeast-specific niche claims are **inference** consistent with the source and RECTIFY's measured win rates.
