# mapPacBio (BBMap): parameters, ONT-DRS suitability, long-read fragility, chunking

`mapPacBio.sh` is BBMap's long-read mode (BBTools, Brian Bushnell). RECTIFY runs
it via `run_map_pacbio()` in `multi_aligner.py`. It is the only non-minimap2-based
aligner in the panel, so it adds algorithmic diversity to the consensus — but it
is **PacBio-tuned and a weak fit for noisy ONT (especially direct-RNA) spliced
alignment** (see "ONT suitability" below). Read this before trusting its intron
calls or adding it to an ONT panel.

---

## Parameter semantics — `intronlen` and `maxindel` (read this first)

Both are easy to misconfigure, and misconfiguration silently produces almost **no
spliced (`N`) introns**. Verbatim from the installed `bbmap.sh --help`:

> `intronlen=999999999` — *"Set to a lower number like 10 to change 'D' to 'N' in
> cigar strings for deletions of at least that length. This is used by Cufflinks;
> 'N' implies an intron while 'D' implies a deletion."*
>
> `maxindel=16000` — *"Don't look for indels longer than this. Lower is faster.
> Set to >=100k for RNAseq with long introns like mammals."*

- **`intronlen` is a MINIMUM D→N relabel threshold, not a maximum.** It does not
  change mapping — only whether a found gap is written `N` (intron) or `D`
  (deletion). At `intronlen=500000`, only deletions ≥500 kb become `N`; every
  human intron (mostly 1–50 kb) stays `D`. **Set it LOW: `intronlen=10`–`20`.**
- **`maxindel` caps the largest gap BBMap will SEARCH for.** At the default
  16000, mapPacBio cannot span introns >16 kb — it soft-clips/fragments the read
  instead. For mammalian RNA, **set `maxindel=200000`** (Bushnell's explicit
  guidance). `strictmaxindel=f` (default) means larger indels may still be found
  opportunistically, but the search is not directed beyond `maxindel`.

These compound: with `intronlen` too high AND `maxindel` at default, you get
near-zero `N` ops. **Empirical (human chr5 ONT DRS, 2026-05-25):** mapPacBio
emitted only ~1,046 introns vs 200k–418k for minimap2/uLTRA/deSALT on the same
reads — the signature of this misconfiguration.

> ⚠️ The 2026-05-24 "intronlen fix" (`intronlen=50` → `intronlen=max_intron`) was
> based on a misreading of the docs and made labeling worse. For yeast (introns
> <1 kb, default `maxindel` adequate) the symptom is masked. See AGENT_FIXES.md.

**Correct RECTIFY invocation** (`run_map_pacbio`, multi_aligner.py ~L723):

```
intronlen=10            # relabel any deletion >=10 bp as an intron (N)
maxindel=200000         # actually search across mammalian-scale introns
```

Decouple these from the genome-wide `--max-intron` that legitimately feeds
minimap2's `-G` and gapmm2's terminal window — they are different parameters.

**Fix validated (SMA_GSB2394 human chr5 ONT DRS, full 61k reads via rectify,
2026-05-25):** fixed params recovered **108,840 introns vs 158** with the broken
params — a **689×** difference. The broken params were catastrophic for splicing;
the fix works. (multi_aligner.py corrected 2026-05-25; original Sherlock copy
backed up as `multi_aligner.py.bak_20260525`.)

---

## ONT direct-RNA suitability

BBMap's error model assumes PacBio-style short indels, which does not match noisy
ONT (especially direct-RNA) reads. In the Křižanović 2018 long-read splice-aligner
benchmark, BBMap correctly spanned exon–exon junctions in only **~26.8 %** of ONT
reads (GMAP 87.1 %), aligned a smaller fraction of reads overall, and produced
more mispositioned alignments. The community standard for ONT spliced alignment is
**minimap2 `-ax splice -uf -k14`** (`-uf` forces forward transcript strand for
stranded DRS; `-k14` for noisy reads); uLTRA and deSALT are the strong
splice-aware alternates.

**Measured on this data (SMA_GSB2394, fixed params, vs minimap2/uLTRA/deSALT,
2026-05-25):** even correctly configured, mapPacBio is the **noisy outlier**.
**21.7 %** of its introns (23,590) are unique — found by no other aligner for that
read (>5 bp) — and **97.7 %** of those uniques (23,037) are also >5 bp from any
annotated GENCODE junction ("novel"). It contributes only ~553 unique *annotated*
junctions — the same as gapmm2 — among ~23k low-quality ones. Note "novel" here
conflates two cases a ±5 bp check can't separate: **fully spurious** introns
(invented where no junction exists) and **mispositioned-real** introns (a real
junction placed >5 bp off — BBMap's known ONT imprecision). Both are bad for
precise junction / 3′-end calling, but it is NOT proven they are all fake. The
robust signal is the **disagreement rate**: 21.7 % unique vs the panel, against
gapmm2's 1.8 % — mapPacBio is the outlier regardless. (To split spurious from
mispositioned, re-classify the novel-unique introns at a wider, e.g. ±50 bp,
annotation tolerance.)

**Consensus-level proof — and an HP-ED metric vulnerability (SMA_GSB2394,
2026-05-26).** The artifact-count alone isn't the verdict (the consensus is
supposed to discard losing alignments). So we ran the 4-aligner consensus with
fixed mapPacBio and looked at what it actually *selects*: mapPacBio **sole-wins
38,126 reads (39 %)** — won, with no other aligner agreeing (`Xn==1`). But on
exactly those reads, the three aligners it beat had **97–98 %** annotated
junctions while mapPacBio's winning junctions were only **77 %** annotated
(18,290 novel). Its unique contribution there: **503 real vs 17,067 novel
(1:34)**. So **HP-ED selected the strictly worse alignment.** Mechanism: HP-ED
scores N-ops as *free* and penalizes soft-clips ~1/base, so mapPacBio's
`maxindel=200000` habit of converting soft-clips into (spurious) introns makes it
win regardless of junction correctness. The consensus *selects* mapPacBio's
artifacts rather than filtering them. (Counterintuitively, the param fix makes
mapPacBio MORE harmful: broken-params it won mostly *unspliced* reads — no
junctions to get wrong; fixed-params it wins *spliced* reads with bad junctions.)

> ⚠️ **Latent HP-ED vulnerability (independent of mapPacBio):** any aligner that
> over-converts soft-clips to N-ops can game the consensus score, since N is free
> and soft-clips are penalized. Worth a metric rebalance (e.g. a penalty for
> unannotated/non-canonical N-ops). See [[feedback_hp_edit_distance_semantics]].

**Recommendation:** do **not** include mapPacBio in an ONT-DRS splice panel — it
*degrades* the consensus, replacing 97–98 %-annotated alignments with 77 % ones on
the reads it wins. Use minimap2 + uLTRA + deSALT. Dropped from the human-DRS panel
2026-05-26. Keep mapPacBio for **yeast**, where BBMap's short-indel error model
fits and it has earned its place. (Sources: Křižanović 2018 PMC6192213; uLTRA
PMC8665758.)

---

## Long-read fragility (fastareadlen, read splitting, assertions)

mapPacBio asserts/crashes on long reads unless handled:

- **`fastareadlen`**: the `mapPacBio.sh` wrapper defaults to 6000; reads >6019 bp
  trigger `AssertionError`. RECTIFY sets `fastareadlen=100000` **and** splits
  reads >`MAX_MPB_READ_LENGTH` (~6 kb) via `split_long_reads()` /
  `stitch_split_bam()` (mpb_split_reads.py) before alignment.
- **`printOutputStats` count assertion** (e.g. `198+137+0+135+0+0 = 470 != 471`,
  `align2.BBMapPacBio`): BBMap crashes reconciling read counts. RECTIFY tolerates
  the related `realign_new` AssertionError ("BBMap spliced greflen overflow",
  upstream BBTools issue #19) by **accepting the partial SAM** and logging
  `"accepting partial SAM — N read(s) lost"`. A **direct** `mapPacBio.sh` call
  that bypasses RECTIFY's read-splitting + partial-SAM tolerance will hit these
  crashes on long ONT reads — do not benchmark mapPacBio outside `run_map_pacbio`.
- **QNAME sanitization** (`_sanitize_mpb_fastq`, `trd=t`): ONT Dorado headers
  embed run metadata; mapPacBio copies the full header into the SAM QNAME, which
  violates the 254-char SAM limit and makes `samtools view` exit 1
  ("query name too long"). RECTIFY strips space- and tab-separated comments.

---

## 6 h `ALIGNER_TIMEOUT` and chunking

`ALIGNER_TIMEOUT = 21600 s` (multi_aligner.py:232) caps each aligner subprocess.
Large samples (≳400k reads) exceed it on a single mapPacBio pass. Use the built-in
interleaved chunking:

```bash
# one task per chunk (1/N of reads, even read-length distribution):
rectify align <reads> --aligners mapPacBio --mapPacBio-chunks 8 --mapPacBio-chunk-idx K ...
# merge (NO --chunk-idx) once all N chunk BAMs exist:
rectify align <reads> --aligners mapPacBio --mapPacBio-chunks 8 ...
```

> ⚠️ **chunk-idx exit-1 bug:** in chunk-idx mode the chunk alignment succeeds and
> writes `{prefix}.mapPacBio.chunk_K_of_N.bam`, but the task **exits 1** —
> `align_command._run_one_aligner` validates the standard `mapPacBio.bam` path,
> not the chunk path `run_map_pacbio` actually returned. The chunk BAMs are
> valid. **Do not put `--dependency=afterok` on the chunk array** (it will never
> be satisfied); gate the merge on all N chunk BAMs existing instead. Merge mode
> is unaffected (it writes the standard path). See AGENT_FIXES.md 2026-05-25.

Chunk-local `RN:i` tags collide across chunks (each chunk reindexes its own
FASTQ); they are corrected downstream by the consensus step's QNAME→RN
reinjection, so merged-then-consensus output is correct.

---

## `pt:i` tag embedded in the read NAME

mapPacBio embeds the poly-A `pt:i:N` tag into the QNAME (not as a SAM aux tag).
The separator depends on processing stage:

- **Direct output (pre-sort):** space-separated — `UUID pt:i:25`.
- **After `samtools sort`/`merge`:** the space becomes `_` (BAM forbids spaces in
  QNAME) — `UUID_pt:i:25`.

```python
if " pt:i:" in name: clean = name.split(" pt:i:")[0]   # pre-sort
if "_pt:i:" in name: clean = name.split("_pt:i:")[0]   # sorted/merged
```

---

## Verifying mapPacBio works

```bash
rectify align --reads <reads.fastq.gz> --genome <genome.fa> \
    --aligners mapPacBio --max-intron 200000 \
    --output /tmp/mpb_smoke --threads 8
samtools flagstat /tmp/mpb_smoke/*.mapPacBio.bam
# sanity: spliced reads should produce N ops, not just D:
samtools view /tmp/mpb_smoke/*.mapPacBio.bam | awk '$6 ~ /N/' | wc -l
```

If the `N`-op count is near zero on a mammalian RNA sample, suspect the
`intronlen`/`maxindel` misconfiguration above.

---

## Failure modes quick-reference

| Symptom | Cause | Fix |
|---------|-------|-----|
| ~0 `N` introns on mammalian RNA | `intronlen` too high (D not relabeled to N) | `intronlen=10`–`20` |
| Few junctions; long introns soft-clipped | `maxindel` at default 16000 | `maxindel=200000` |
| `AssertionError` on reads >6019 bp | `fastareadlen` too small / no read splitting | RECTIFY patches `fastareadlen=100000` + splits >6 kb reads; run via `run_map_pacbio` |
| `printOutputStats ... X != Y` crash | BBMap stats / `realign_new` overflow on long reads | RECTIFY accepts partial SAM; don't invoke `mapPacBio.sh` directly |
| `samtools view: query name too long`, exit 1 | Full ONT header copied to QNAME | `trd=t` / `_sanitize_mpb_fastq` |
| chunk task exits 1 but chunk BAM exists | chunk-idx path-validation bug | ignore exit code; gate merge on chunk-BAM presence (no `afterok`) |
| Subprocess killed at 21600 s | `ALIGNER_TIMEOUT` on a large sample | chunk with `--mapPacBio-chunks N` |
| Mispositioned / missing junctions on ONT | BBMap PacBio error model ≠ ONT | use minimap2/uLTRA/deSALT; consider dropping mapPacBio from ONT panels |

---

## Primary-alignment & duplicate handling

`rectify align` passes **no secondary-suppression flag** to mapPacBio — its command
is `ref/in/out/threads/path/fastareadlen/intronlen/maxindel/minratio/-Xmx`
(`multi_aligner.py:739-757`), with no `secondary=` argument and no post-hoc dedup.
BBMap's default secondary-output behavior is not determinable from rectify source.
This is safe for 3′-end counting because `rectify correct` skips
`is_secondary`/`is_supplementary` records regardless of producing aligner.

> **Empirical (2026-05-29): mapPacBio emits MULTIPLE records all flagged
> PRIMARY for a multi-mapping read** — it does *not* set the secondary
> (`0x100`) bit on the extra hits. Observed on yeast upf1Δ DRS: one read had
> three primary records (`flag=0` chrI+, `flag=0` chrIX+, `flag=16` chrXI−),
> another had two (chrII− no-intron, chrVII− with-intron). Consequences:
> - **Selecting "the" mapPacBio alignment by first-non-secondary record is
>   unsafe** — it lands on an arbitrary locus. Disambiguate by expected
>   chromosome / best score (e.g. the validation-bundle builder
>   `scripts/validation_data/upf1d_2026_05/build_upf1d_validation.py` picks the
>   primary on the intended chrom, preferring the N-op-bearing record).
> - **The `is_secondary` skip does NOT filter these** (they aren't flagged
>   secondary), so a mapPacBio multi-mapper reaches `rectify correct` as several
>   primaries → it can be processed/counted more than once. This is the same
>   class of hazard as the duplicate-primary case noted below, but produced by
>   the aligner itself rather than a doubled input FASTQ. Another reliability
>   mark against mapPacBio (cf. its junk 5′ force-extension that raw-ED
>   under-penalizes — see ALIGNER_RECOMMENDATIONS / the panel-drop rationale).

Note that mapPacBio reindexes per chunk, so chunk-local `RN:i` tags collide across
chunks (corrected downstream by the consensus QNAME→RN reinjection — see "6 h
`ALIGNER_TIMEOUT` and chunking" above); that is a tag-collision issue, not record
duplication.

The unguarded hazard for *all* aligners is duplicate **primary** records in an
external BAM (e.g. a doubled input FASTQ), which `rectify correct` does **not**
dedup → 2× double-counted 3′ ends. See the canonical writeup and cross-aligner
table in [minimap2.md](minimap2.md#-duplicate-primary-alignments--2-double-counted-3-ends-external-bam-hazard).
