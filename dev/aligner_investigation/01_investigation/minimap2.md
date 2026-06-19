# minimap2 — Source-Level Investigation (Long-Read Splice Aligner)

> **Build note:** code-level claims here were verified against `master`; see `../CORRECTIONS_vs_DRS_BUILD.md` for re-verification vs `origin/drs-validation-rebuild` (flag sets CONFIRMED — minimap2 set + `-y`; penalty tables present & auto-resolved on the build).

**Scope:** minimap2 as used by RECTIFY for ONT direct RNA-seq (DRS) / cDNA splice
alignment in *S. cerevisiae*. Grounded in the RECTIFY wrapper
(`rectify/core/align/multi_aligner.py::run_minimap2`), the Li 2018 *Bioinformatics*
paper and its LaTeX source (`tex/minimap2.tex`), the 2021 "New strategies" paper,
the man page (`minimap2.1`), and the splice DP source (`ksw2_exts2_sse.c`,
`align.c`, `options.c`).

**Convention used below:** *FACT* = directly from source/paper/man page;
*INFERENCE* = my reasoning, clearly flagged. Numeric defaults are FACT from
`options.c::mm_idxopt_init` / `mm_mapopt_init` / `mm_set_opt` unless noted.

---

## 1. Algorithm Pipeline (seed → chain → align)

minimap2 is a **seed–chain–align** aligner. Three stages:

### 1.1 Minimizer seeding

*FACT.* A **(w,k)-minimizer** is the smallest k-mer (by a hash of its 2-bit
encoding) within every window of `w` consecutive k-mers. minimap2 collects
minimizers of the reference and stores them in a hash table: **key = hash of the
minimizer, value = list of (position, strand) of its occurrences** (`tex/minimap2.tex`:
"It collects minimizers of the reference sequences and indexes them in a hash
table, with the key being the hash of a minimizer and the value being a list of
locations of the minimizer copies").

Default index parameters (`mm_idxopt_init`): `k=15`, `w=10`, `bucket_bits=14`
(the hash table is bucketed; the low 14 bits of the hash select a bucket).
**Preset `splice` overrides to `k=15, w=5`** (`mm_set_opt`). RECTIFY then forces
`k=14` on the command line (see §3).

- **w=5 (splice preset):** denser sampling than genomic `w=10` → on average a
  minimizer every ~3 bp (`(w+1)/2`), raising sensitivity for short terminal exons.
- **HPC minimizers (`MM_I_HPC`, `-H`):** homopolymer-compressed seeds, set only by
  `map-pb`. *FACT (paper/cookbook):* HPC helps PacBio CLR but **hurts Nanopore**,
  which is why ONT presets and RECTIFY do not use it.

**Occurrence filtering of repetitive minimizers.** *FACT.* To bound the
quadratic chaining cost, minimap2 ignores minimizers that are too frequent. The
man-page/`mm_mapopt` knobs:
- `-f FLOAT` (`mid_occ_frac`, default `2e-4`): minimizers in the top fraction `f`
  of occurrence counts are flagged as repetitive and filtered before chaining.
- `--mid-occ` / `min_mid_occ`..`max_mid_occ`: the frequency cutoff derived from `f`
  is clamped into `[min_mid_occ, max_mid_occ]` (HiFi preset raises these to
  `50..500`). For yeast (12 Mb, low repeat content) this filter is rarely binding.

### 1.2 Chaining (SDP)

*FACT (paper §2.2 / `tex/minimap2.tex`, `lchain.c::mm_chain_dp`).* Anchors are
collinear seed matches `(x_i, y_i)` (reference pos, query pos), sorted by `x`. The
chain score of anchor `i` is the **sparse-DP recurrence**:

```
f(i) = max{ max_{i>j>=1} [ f(j) + α(j,i) − β(j,i) ],  w_i }
```

- `α(j,i) = min{ min(y_i−y_j, x_i−x_j), w_i }` — number of matching bases added
  (capped at the seed length `w_i`).
- `β(j,i) = γ_c( (y_i−y_j) − (x_i−x_j) )` — gap penalty as a function of the
  difference in reference vs query offset (collinearity deviation).
- Gap cost (genomic):
  ```
  γ_c(l) = 0.01 · w̄ · |l| + 0.5·log2|l|   (l ≠ 0);   0   (l = 0)
  ```
  where `w̄` is the average seed length. The `log2` term lets one chain absorb a
  long gap cheaply (intron-friendly).

**Inner-loop heuristic.** *FACT.* The naive recurrence is O(n²). minimap2 starts
the predecessor search at `i−1` and stops after up to **h = 50** iterations
without improvement; predecessors farther than the max-gap `G` get infinite cost.
This makes chaining near-linear in practice.

**Backtracking.** *FACT.* `P(i)` stores the best predecessor; chains are recovered
by following `P(i)` from high-`f` anchors until `P(i)=0` or an already-used anchor,
marking anchors used. Sub-optimal/secondary chains are kept (subject to
`--secondary` and a chain-score threshold `min_chain_score = 40`).

**2021 update (RMQ chaining).** *FACT (Li 2021, "New strategies").* The newer
release adds a **range-minimum-query (RMQ)-based chaining** path that allows much
longer gaps/joins efficiently (relevant for long introns and structural events),
replacing/augmenting the original bounded "range chaining." For yeast splice
alignment (introns < 1 kb) this is a minor effect.

### 1.3 Base-level alignment (DP)

*FACT (paper §2.3; `align.c::mm_align1`; `ksw2_ext*2_sse.c`).* Between adjacent
anchors in a chain, minimap2 fills the gaps with **KSW2**, a banded global DP
using the **Suzuki–Kasahara difference-based ("diagonal") formulation** with a
**two-piece affine gap cost**:

```
γ_a(l) = min{ q + |l|·e ,  q̃ + |l|·ẽ }
```

(`q` = short-gap open, `e` = short-gap extend; `q̃`,`ẽ` = long-gap open/extend.)
The two-piece form lets a single long deletion/intron be penalized by the cheaper
linear piece. KSW2 is **SIMD-vectorized (SSE)**: `ksw_extz2_sse` (no splice),
`ksw_extd2_sse` (two-piece affine, the genomic default), and **`ksw_exts2_sse`**
(splice DP — the one used here).

**Z-drop heuristic.** *FACT.* To avoid forcing alignment through unrelated
sequence between anchors, the alignment is broken where the score drops too fast:

```
break if ∃ (i',j'),(i,j) with i'<i, j'<j :  S(i',j') − S(i,j) > Z + e·|(i−i') − (j−j')|
```

Default `-z 400,200` (genomic); **splice preset uses `zdrop=200, zdrop_inv=100`**.
Unlike BLAST X-drop, Z-drop tolerates a single long gap (so introns survive).

**CIGAR / MD.** *FACT.* Backtracking emits the CIGAR. Splice gaps are emitted as
**`N` ops** (`KSW_CIGAR_N_SKIP`); short reference gaps as `D`. The MD tag (and
`=`/`X`) are produced by post-processing against the reference when `--MD` is set.

---

## 2. Splice Model Internals

*FACT (paper "Spliced alignment"; `ksw2_exts2_sse.c`; `tex/minimap2.tex`).*

### 2.1 Intron as a long gap with a separate DP state

`ksw_exts2_sse` adds a **dedicated splice gap state** (variable `x2`, initialized
to `−q2`) alongside the normal gap states (`x,y,u,v`, initialized to `−(q+e)`).
The intron state opens at cost `q2` (the **large** second gap-open, `q2=32` in the
splice preset) and extends at `e2=0` — i.e. **once an intron opens, extending it
across thousands of bases costs nothing**. This is exactly the "long-gap" piece of
the two-piece affine cost specialized for introns. The state transition
(`z = max(z, a, b, a2a)`, where `a2a` is the splice pathway) lets the DP choose
"open an intron" vs "normal deletion" per cell; the boundary length implied by
`q2` vs `q+|l|·e` decides when an `N` op beats a `D` op.

*INFERENCE:* This is why minimap2 (unlike a plain deletion model) will emit a long
`N` rather than a long `D` once the reference gap exceeds the break-even length —
and why setting `-G` (max intron) constrains where introns may be placed.

### 2.2 Donor/acceptor signal scoring (GT-AG)

*FACT (`tex/minimap2.tex`).* Position-specific donor `d(i)` and acceptor `a(i)`
penalties are added at candidate intron boundaries (precomputed `donor[]` /
`acceptor[]` arrays in `ksw_exts2_sse`):

- Donor `d(i)`: `0` if `T[i+1..i+3]` is `GTA` or `GTG`; `p/2` if `GTC` or `GTT`;
  `p` otherwise.
- Acceptor `a(i)`: `0` if `T[i−2..i]` is `CAG` or `TAG`; `p/2` if `AAG` or `GAG`;
  `p` otherwise.

So canonical **GT…AG** with the favorable flanking base costs nothing; the
non-canonical penalty `p` is the `noncan`/`-C` value (**splice: `C=9`**;
`splice:hq`: `C=5`). In the DP, the canonical dinucleotide checks are literal:
donor `target[t+1]==G(2) && target[t+2]==T(3)`, acceptor `target[t−1]==A(0) &&
target[t]==G(2)`.

### 2.3 `--splice-flank` (and why RECTIFY disables it)

*FACT.* The **flanking model** is the "extra base" beyond GT/AG: a real donor's
4th base tends to be `A/G` (the `GTA/GTG` → 0 vs `GTC/GTT` → p/2 split above), and
the acceptor's preceding base tends to be `C/T`. minimap2 docs: this trend is
91–92% in human/mouse and "evolutionarily conservative, all the way to
*S. cerevisiae*"; `--splice-flank=yes` (default in `splice`, flag
`MM_F_SPLICE_FLANK`) "generally leads to higher junction accuracy by several
percent." **`--splice-flank=no` makes the model score only GT..AG, ignoring the
extra base.**

*INFERENCE (why RECTIFY sets `--splice-flank=no`):* The flanking bonus biases the
**exact intron-boundary placement** toward the statistically-favored flanking base.
RECTIFY's downstream junction refinement (Module 2H, `junction_refiner.py`) and its
3'-end correction are **sequence-first** and use an *empirical/annotation-driven*
junction model; a built-in aligner prior that nudges boundaries by ±1 bp toward
`GTA/GTG` is an *uncontrolled* second opinion that can shift the reported 3' end /
junction off the true site. Disabling it makes minimap2's junction placement
depend only on the raw GT..AG signal + the read sequence, leaving boundary
refinement to RECTIFY. (The CLAUDE.md note "important for 3' end accuracy"
corroborates this intent.) This matches the documented use of `--splice-flank=no`
for organisms/controls (e.g. SIRV) that don't honor the flanking trend.

### 2.4 Strand semantics: `-uf` / `-ub` / `-un`

*FACT.* `-u` controls how the transcript strand is found from the splice signal.
By default (`-x splice` sets `-ub`) minimap2 aligns each chain **twice** — once
assuming `GT…AG`, once assuming the reverse-complement `CT…AC` — and keeps the
higher-scoring strand (flags `MM_F_SPLICE_FOR | MM_F_SPLICE_REV`). `-uf` =
**consider the forward transcript strand only** (`GT…AG` only); `-ub` = both;
`-un` = don't use the splice signal to determine strand. ONT direct-RNA reads are
**stranded** (sequenced 3'→5' as RNA, reported as the sense strand), so `-uf` is
correct and ~halves splice DP work. *FACT:* the official ONT DRS recipe is exactly
`minimap2 -ax splice -uf -k14 ref.fa drs.fq`.

### 2.5 `--junc-bed` / `--junc-bonus` (how annotation priors enter the DP)

*FACT (`ksw2_exts2_sse.c`, man page).* A junction BED (RECTIFY builds it from the
GFF: `annotation.junc.bed`, via `utils/junction_bed.py`) marks known donor/acceptor
sites. In the DP, junction-array entries adjust the per-position donor/acceptor
penalty:
- **Bonus mode:** positions matching a known splice site get an **additive
  `junc_bonus`** (default and RECTIFY value **9**) to the score, biasing the DP to
  place the intron at the annotated boundary.
- **Score mode (`KSW_EZ_SPLICE_SCORE`):** `junc[]` can carry per-site donor scores
  that *override* the canonical penalties:
  `donor[t] += junc[t+1]==0xff ? −junc_pen : (junc[t+1]>>1) − KSW_SPSC_OFFSET`
  (`junc_pen=5` in the splice preset). This is the hook for externally-predicted
  junction scores.

*FACT (CLAUDE.md / wrapper docstring):* the annotation is used to **improve** the
alignment but scoring stays **annotation-blind enough to still find novel
junctions** — `junc-bonus` is a tie-breaker bonus, not a hard constraint. RECTIFY
relies on this (novel CPA-proximal junctions must still be discoverable).

### 2.6 Canonical vs non-canonical & `-C`

*FACT.* Non-canonical junctions are not forbidden — they pay the `noncan`/`-C`
penalty (`9` for `splice`, `5` for `splice:hq`/Iso-Seq). Lowering `-C` helps
low-error data align genuine non-canonical sites; RECTIFY keeps the `splice`
default (`C=9`, never overridden in the wrapper).

---

## 3. Exact RECTIFY Invocation — Flag-by-Flag Analysis

From `multi_aligner.py::run_minimap2` (the live command; note the wrapper also
appends `-y`, not shown in CLAUDE.md):

```
minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD -y \
         -t <n> --junc-bed annotation.junc.bed --junc-bonus 9 genome reads
```

| Flag | Source default | RECTIFY value | Effect / rationale |
|------|----------------|---------------|--------------------|
| `-a` | — | set | Output **SAM** (then piped to `samtools sort -n`). |
| `-x splice` | — | set | Activates `ksw_exts2_sse` splice DP + preset `k=15,w=5, A=1,B=2, O=2,32, E=1,0, C=9, z=200,100, bw=200k, max_gap=2000, max_gap_ref=200000, junc_bonus=9, junc_pen=5`. **Match `A=1`/mismatch `B=2` are low** (noisy-read scoring). |
| `-u f` | `splice`→`b` | `f` | **Forward transcript strand only.** Correct for stranded DRS; skips the `CT…AC` pass. |
| `-k14` | splice `k=15` | `14` | **Smaller k → more, denser minimizers → higher sensitivity** on ~5–10% error DRS reads (better seeding of short first/last exons). Costs specificity, mitigated by yeast's small low-repeat genome. This is the documented ONT-DRS value. |
| `-G 5000` | splice `200000` | `5000` | **Max intron / max ref gap = 5 kb.** Yeast introns are < ~1 kb, so 5 kb is safely permissive while preventing spurious cross-gene "long-join" introns. Also lowers the chaining `G` cutoff. |
| `--splice-flank=no` | `yes` (splice) | `no` | Disables the GT**A/G**..**C/T**AG flanking bonus. **Hands intron-boundary/3'-end placement to RECTIFY** (see §2.3). |
| `--secondary=no` | `yes` | `no` | One primary alignment per read only. RECTIFY's consensus/correct logic expects a single record per UUID per aligner (the wrapper later strips secondaries for other aligners too). |
| `--MD` | off | set | Emit **MD tag** → required by RECTIFY's indel-artifact correction (Module 2C) and identity computation. Wrapper additionally runs `samtools calmd -e` to convert `M`→`=`/`X`. |
| `-y` | off | set | Copy FASTQ comment fields into BAM aux tags (carries DRS poly-A/pt metadata when present). |
| `-t <n>` | 3 | threads | Thread count. |
| `--junc-bed annotation.junc.bed` | none | set | Known yeast junctions from the GFF → additive bias toward annotated boundaries. |
| `--junc-bonus 9` | 9 | 9 | Bonus added at annotated sites (equals the splice default; explicit for clarity). |

**Not set (defaults inherited):** `B=2`, `O=2,32`, `E=1,0`, `C=9`, `z=200,100`,
`w=5`, `-f 2e-4`, `--end-bonus 0`. *INFERENCE:* `--end-bonus 0` (no reward for
extending the alignment to the read ends) is relevant to RECTIFY's mission —
minimap2 has **no incentive to push the 3' alignment end into a poly-A/homopolymer
tail**, so it soft-clips ambiguous 3' bases, contributing to the imprecise
3'-end behavior RECTIFY corrects (see §6).

**Pipeline note (wrapper):** SAM is streamed to `samtools sort -n` (name-sorted,
no coordinate index) so RECTIFY's cross-aligner consensus can stream-merge; a
background thread drains minimap2 stderr to avoid a 64 KB pipe-buffer deadlock.

---

## 4. Data Structures & Complexity

- **Index:** open hash table, key = 64-bit minimizer hash, value = position list;
  `bucket_bits=14`. Build is O(reference length); memory ≈ proportional to number
  of minimizers (~`2/(w+1)` of genome length). Yeast index is trivially small.
- **Seeding:** O(read length); minimizer extraction is a sliding-window minimum.
- **Chaining:** worst case O(n²) over `n` anchors; bounded to ~O(n·h) by the
  `h=50` predecessor cap and the `G` gap cutoff. RMQ chaining (2021) gives
  near-linear long-gap chaining.
- **Base DP (KSW2):** banded; O(band · alignment length) per inter-anchor gap;
  SIMD throughput ~8–16 cells/instruction. Z-drop and the band keep introns from
  blowing up the matrix.
- **Splice DP state:** the extra `x2` intron state adds a small constant factor
  over `ksw_extd2_sse`.

---

## 5. Strengths

- **Speed + good sensitivity** at ~5–10% ONT error, with `w=5`/`k=14` tuned for
  short terminal exons.
- **Principled splice model:** GT-AG donor/acceptor scoring, separate intron DP
  state with near-free intron extension, two-piece affine gaps, Z-drop that
  tolerates single long gaps.
- **Annotation-assisted but novel-junction-capable:** `--junc-bed`/`--junc-bonus`
  bias without hard constraint.
- **Stranded handling (`-uf`)** matches DRS biology and halves splice work.
- **De-facto standard** (the canonical ONT DRS recipe is minimap2's own), so its
  alignments are a strong baseline in RECTIFY's multi-aligner panel.

---

## 6. Weaknesses (3'-end / homopolymer / intron-boundary precision)

This is precisely the failure surface RECTIFY exists to repair (CLAUDE.md: in
correct-first win rates minimap2 wins only **0.1%**).

- **Imprecise 3' ends at poly-A / homopolymers.** *FACT + INFERENCE.* With
  `--end-bonus 0` and noisy reads, minimap2 has no incentive to extend the 3'
  alignment cleanly through an A-run; it tends to **soft-clip the poly-A tail and
  trailing homopolymer bases**, and the reported `reference_end` (the CPA on `+`)
  drifts. Published nanopore work notes "the 3' end of a read may not represent
  the mature 3' end" and reads ending "in the middle of internal exons." RECTIFY's
  Modules 2B/2E/2G/AG and `find_polya_boundary` exist to walk these back.
- **Homopolymer indel noise.** *FACT.* ONT mis-calls homopolymer length; minimap2's
  DP distributes the resulting indels somewhat arbitrarily within the run (the
  low `A=1,B=2` scoring + affine gaps don't pin a unique boundary). This is the
  staircase/ambiguity problem RECTIFY's HP-aware edit distance addresses.
- **Intron-boundary jitter (±1–few bp).** *FACT/INFERENCE.* When exon-end +
  intron-start sequence resembles intron-end + exon2-start (common at HP-rich
  junctions), the splice DP can place the `N` op a few bp off the true site;
  `--junc-bonus` helps only at annotated sites. RECTIFY's Module 2H re-scores
  candidate junctions sequence-first to fix this.
- **Local-optimum N-vs-D choice.** *INFERENCE.* The `q2` break-even can emit a `D`
  where a short `N` is biologically correct (or vice versa), especially near read
  ends. RECTIFY's chimeric/junction logic reclassifies short ref gaps.
- **`--secondary=no` discards alternatives** that might have placed the 3' end
  better — by design RECTIFY recovers cross-aligner consensus instead.

*Net:* minimap2 is a fast, accurate **chain/splice baseline** but is **not
optimized for single-base 3'-end (CPA) precision** in homopolymer/poly-A context —
the exact gap RECTIFY's correct-first pipeline fills.

---

## 7. PacBio vs ONT Presets (`options.c::mm_set_opt`)

| Preset | k | w | HPC | A | B | O (q,q2) | E (e,e2) | C | z | Use |
|--------|---|---|-----|---|---|----------|----------|---|---|-----|
| **splice** (RECTIFY) | 15→**14** | 5 | no | 1 | 2 | 2,32 | 1,0 | 9 | 200,100 | Noisy long RNA (ONT cDNA/DRS) |
| **splice:hq** | 15 | 5 | no | 1 | **4** | **6,24** | 1,0 | **5** | 200,100 | Accurate RNA (Iso-Seq, HiFi cDNA): higher mismatch/gap-open penalties, cheaper non-canonical |
| **map-ont** (genomic default) | 15 | 10 | no | 2 | 4 | 4,24 | 2,1 | n/a | 400,200 | Genomic ONT (not spliced) |
| **map-pb** | 19 | 10 | **yes** | 2 | 4 | 4,24 | 2,1 | n/a | 400,200 | PacBio CLR: HPC minimizers tame homopolymer error |
| **map-hifi / lr:hq** | 19 | 19 | no | 1 | 4 | 6,26 | 2,1 | n/a | (s=200) | PacBio HiFi: long k/w (low error), strict scoring |

Key contrasts relevant here:
- **HPC** is PacBio-only; it *hurts* ONT, so DRS uses ordinary minimizers + small
  `k`.
- **splice vs splice:hq** differ mainly in scoring stringency (`B`, `O`, `C`):
  HiFi/Iso-Seq can afford harsher mismatch/gap penalties and a lower
  non-canonical cost; noisy ONT cannot. RECTIFY correctly uses `splice` (noisy),
  not `splice:hq`.
- All splice presets share `e2=0` (free intron extension) and the GT-AG model;
  only the genomic `map-*` presets lack the splice state entirely.

---

## 8. Key Source References

| Component | File / function | Paper section |
|-----------|-----------------|---------------|
| Index / minimizers | `index.c`, `sketch.c`; hash table w/ `bucket_bits` | §2.1 "Indexing" |
| Defaults / presets | `options.c::mm_idxopt_init`, `mm_mapopt_init`, `mm_set_opt` | — |
| Chaining SDP | `lchain.c::mm_chain_dp` (recurrence `f(i)`, `α`,`β`,`γ_c`, `h=50`, `P(i)`) | §2.2 "Chaining" |
| RMQ / long-gap chaining | `lchain.c` (2021) | Li 2021 "New strategies" |
| Per-read driver | `map.c`, `align.c::mm_align1` | §2.3 |
| Genomic base DP | `ksw2_extd2_sse.c::ksw_extd2_sse` (two-piece affine `γ_a`) | §2.3 |
| **Splice base DP** | **`ksw2_exts2_sse.c::ksw_exts2_sse`** (intron state `x2`, `donor[]`/`acceptor[]`, `junc_bonus`, `junc_pen`, `noncan`, `KSW_SPSC_OFFSET`, `KSW_CIGAR_N_SKIP`) | "Spliced alignment" |
| Z-drop | `align.c` / KSW2 (`zdrop`, `zdrop_inv`) | §2.3 "Z-drop" |
| Splice signal model | `tex/minimap2.tex` donor `d(i)` / acceptor `a(i)` tables | "Spliced alignment" |
| Strand / `-uf` | `MM_F_SPLICE_FOR`/`MM_F_SPLICE_REV` in `options.c`/`map.c` | "Spliced alignment" |
| `--junc-bed`/`--junc-bonus` | `mm_idx_bed_read` → `junc[]` in `ksw_exts2_sse.c` | man page |
| Equations (verbatim) | `tex/minimap2.tex` | whole |

**Primary sources consulted:** Li 2018, *Minimap2: pairwise alignment for
nucleotide sequences*, Bioinformatics 34(18):3094 (and arXiv:1708.01492);
Li 2021, *New strategies to improve minimap2 alignment accuracy*, Bioinformatics
37(23):4572; minimap2 man page (`minimap2.1`); minimap2 source
(`options.c`, `lchain.c`, `align.c`, `ksw2_exts2_sse.c`, `tex/minimap2.tex`);
RECTIFY `rectify/core/align/multi_aligner.py::run_minimap2`.

*Caveat:* a few numeric defaults were read from the man page / `options.c` via
web mirrors (arXiv/PMC/Oxford direct fetches were 403-blocked); the splice-preset
scoring (`A=1,B=2,O=2,32,E=1,0,C=9,z=200,100,w=5,junc_bonus=9,junc_pen=5`) is
internally consistent across the man-page preset string and `mm_set_opt`.
